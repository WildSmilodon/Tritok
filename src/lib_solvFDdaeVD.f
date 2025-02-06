
C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsB_DC_fill4(mesh,smatH,smatG,smatB,beta,diff,diffGrad,qdiff,sysm,rhsm,
     &                      velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)
C
C     $: Iz pravokotnih matrik g, h in b naredi CRS sistemsko in rhs, verzija brez laplace
C        pri kateri je hitrost v sistemski
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      REAL(8) smatH(mesh%nicnsp,mesh%npoc),smatG(mesh%nicnsp,mesh%npofc)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc),beta
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) diff(mesh%nnodes)
      REAL(8) qdiff(mesh%nq)
      REAL(8) diffGrad(mesh%nnodes,3)
      REAL(8) val
      TYPE(matrix) :: sysm,rhsm

      INTEGER wqkode(mesh%nq)
      INTEGER wkode(mesh%nnodes)
      INTEGER nunk, nb

      INTEGER ic,it,row,col,ii,jj,nuq,j,is,ir

      wqkode=mesh%DCqkode
      wkode=mesh%DCkode
      nunk=mesh%Dcnunk
      nb=mesh%Dcnb

      is=0
      ir=0
      DO ic=1,mesh%nicell ! zanka po celicah
        DO it=1,mesh%nsp ! po izvornih tockah znotraj celice
          row=(ic-1)*mesh%nsp+it ! vrstica v sistemski matriki
C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            val=diff(j)*smatH(row,col)+smatB(row,col)*beta
     &               +(
C                    advection
     &               +smatAbdx(row,col)*velocity(j,1)
     &               +smatAbdy(row,col)*velocity(j,2)
     &               +smatAbdz(row,col)*velocity(j,3)
     &                 -smatDx(row,col)*(velocity(j,1)+diffGrad(j,1))
     &                 -smatDy(row,col)*(velocity(j,2)+diffGrad(j,2))
     &                 -smatDz(row,col)*(velocity(j,3)+diffGrad(j,3))
     &                )

            nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              is=is+1
              sysm%v(sysm%d(is))=val
            ELSE ! znana vrednost - rhs
              ir=ir+1
              rhsm%v(rhsm%d(ir))=-val  ! minus zato, ker gre na drugo stran enacbe
            END IF
          END DO
C
C         nato G matrika (q)
C
          DO ii=1,mesh%nside
            DO jj=1,mesh%npof
              col=(ii-1)*mesh%npof+jj
              j=mesh%ibf(ic,ii,jj) ! stolpec za G
              nuq=wqkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                is=is+1
                sysm%v(sysm%d(is))=-smatG(row,col)*qDiff(j)
              ELSE ! znana vrednost - rhs
                ir=ir+1
                rhsm%v(rhsm%d(ir))=smatG(row,col)*qDiff(j) ! minus zato, ker gre na drugo stran enacbe
              END IF
            END DO
          END DO
        END DO
      END DO

      END





C -----------------------------------------------------------------------------
C
C       Solve diffusion advection equation using finite difference discretization
C       of time derivative where diffusivity varies in space
C       overdeterminet system of equations
C
      SUBROUTINE SolveFDdaeVD4(env,io,inp,cpu,mesh,sysm,rhsm,precv,
     &                      u,qu,velocity,rhsv,smatB,
     &                      smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nit,
     &                      diff,diffGrad)
C
C     $: resi H*u=G*q+B*rhsv-advekcija
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(matrix) :: sysm,rhsm
      TYPE(IOtype) io
      TYPE(inputtype) inp
      TYPE(CPUtype) cpu
      TYPE(penv) :: env

      INTEGER i,ic,it,j,nit,ierr
      INTEGER node

      REAL(8) rhsv(mesh%nnodes)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc)
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc)

      REAL(8) precv(mesh%DCnunk)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) u(mesh%nnodes)
      REAL(8) qu(mesh%nq)
      REAL(8) diff(mesh%nnodes)
      REAL(8) diffGrad(mesh%nnodes,3)


      REAL(8), ALLOCATABLE :: x(:), b(:), r(:)
      REAL(4) cptime,rcpu

      ALLOCATE (x(mesh%DCnunk),b(mesh%DCnb),r(mesh%nicnsp))

      rcpu=cptime(0.)
C
C     Set up rhs vector and initial approximation of the unknown vector
C
      DO i=1,mesh%nnodes
        IF (mesh%DCkode(i).LT.0) THEN
          x(ABS(mesh%Dckode(i)))=u(i)
        ELSE
          b(mesh%DCkode(i))=u(i)
        END IF
      END DO
      DO i=1,mesh%nq
        IF (mesh%Dcqkode(i).LT.0) THEN
          x(ABS(mesh%Dcqkode(i)))=qu(i)
        ELSE
          b(mesh%DCqkode(i))=qu(i)
        END IF
      END DO
C
C     r = rhs * b
C
      CALL CRSxV(rhsm,b,mesh%DCnb,r)
c
c     r = r - B * rhsv
c
      DO j=1,mesh%npoc
        i=0
        DO ic=1,mesh%nicell
          DO it=1,51 ! 51 izvornih tock v vsaki celici
            i=i+1  ! vrstica v matriki
              node=mesh%idc(ic,j)
              r(i)=r(i)
C                    time derivative + sources
     &               -smatB(i,j)*rhsv(node)
c     &               -u(node)*(
C                    advection
c     &               +smatAbdx(i,j)*velocity(node,1)
c     &               +smatAbdy(i,j)*velocity(node,2)
c     &               +smatAbdz(i,j)*velocity(node,3)
c     &                 -smatDx(i,j)*(velocity(node,1)+diffGrad(node,1))
c     &                 -smatDy(i,j)*(velocity(node,2)+diffGrad(node,2))
c     &                 -smatDz(i,j)*(velocity(node,3)+diffGrad(node,3))
c     &                )
          END DO
        END DO
      END DO
      cpu%time(21)=cpu%time(21)+cptime(rcpu)
      rcpu=cptime(0.)
C
C     solve overdetermined system of linear equations
C
      CALL slvlsqr2(sysm%neq,mesh%DCnunk,sysm%nnz,inp%lsqs_maxit,inp%dlse(8),nit,ierr,
     &              precv,sysm%i,sysm%j,sysm%v,r,x)


      cpu%time(22)=cpu%time(22)+cptime(rcpu)

      IF (ierr.NE.0) CALL WarnErr(env,io,inp,4,"SolveFDdaeVD","NAPAKA V SOVLERJU",ierr)
C
C     Copy solution vector to q and u using under-relaxation
c     under-relaxation should be 1.0 for linear problems
C
      DO i=1,mesh%nnodes
        IF (mesh%DCkode(i).LT.0) THEN
          u(i)=(1.0D0-inp%urDC)*u(i)+inp%urDC*x(ABS(mesh%DCkode(i)))
        END IF
      END DO
      DO i=1,mesh%nq
        IF (mesh%DCqkode(i).LT.0) THEN
          qu(i)=(1.0D0-inp%urDC)*qu(i)+inp%urDC*x(ABS(mesh%DCqkode(i)))
        END IF
      END DO


      DEALLOCATE (x,b,r)

      END






C -----------------------------------------------------------------------------      
C
C       Solve diffusion advection equation using finite difference discretization
C       of time derivative where diffusivity varies in space
C       overdeterminet system of equations
C
      SUBROUTINE SolveFDdaeVD2sq(env,io,inp,cpu,mesh,sysm,rhsm,prec,
     &                      u,qu,velocity,rhsv,smatB,
     &                      smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nit,
     &                      diff,diffGrad)
C
C     $: resi H*u=G*q+B*rhsv-advekcija
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(matrix) :: sysm,rhsm,prec
      TYPE(IOtype) io
      TYPE(inputtype) inp
      TYPE(CPUtype) cpu
      TYPE(penv) :: env

      INTEGER i,ic,j,nit,ierr
      INTEGER node,inode,isqr,isuma
      REAL(8) predznak

      REAL(8) rhsv(mesh%nnodes)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc)
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc)

      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) u(mesh%nnodes)
      REAL(8) qu(mesh%nq)
      REAL(8) diff(mesh%nnodes)
      REAL(8) diffGrad(mesh%nnodes,3)


      REAL(8), ALLOCATABLE :: x(:), b(:), r(:)
      REAL(4) cptime,rcpu,xcpu

      ALLOCATE (x(mesh%DCnunk),b(mesh%DCnb),r(mesh%DCnunk))

      rcpu=cptime(0.)
C
C     Set up rhs vector and initial approximation of the unknown vector
C
      DO i=1,mesh%nnodes
        IF (mesh%DCkode(i).LT.0) THEN
          x(ABS(mesh%Dckode(i)))=u(i)
        ELSE
          b(mesh%DCkode(i))=u(i)
        END IF
      END DO
      DO i=1,mesh%nq
        IF (mesh%Dcqkode(i).LT.0) THEN
          x(ABS(mesh%Dcqkode(i)))=qu(i)
        ELSE
          b(mesh%DCqkode(i))=qu(i)
        END IF
      END DO
C
C     r = rhs * b
C
      CALL CRSxV(rhsm,b,mesh%DCnb,r)
c
c     r = r - B * rhsv
c
      DO isqr=1,sysm%neq !nunk ! zanka po vrsticah kvadratne matrike
        inode=ABS(mesh%DCeql(isqr))  ! vozlisce v U ali Q mrezi, kateremu pripada ta enacba
        IF (mesh%DCeql(isqr).GT.0) THEN ! to pomeni, da je neznanka funkcija
c         dodajamo prispevke enack, ki se sestejejo
          DO isuma=1,mesh%sqUlistNO(inode)
            i=mesh%sqUlist(inode,isuma)
            ic=mesh%sqUlistIC(inode,isuma)
            DO j=1,mesh%npoc
              node=mesh%idc(ic,j)
              r(isqr)=r(isqr)
C                    time derivative + sources
     &               -smatB(i,j)*rhsv(node)
     &               -u(node)*(
C                    advection
     &               +smatAbdx(i,j)*velocity(node,1)
     &               +smatAbdy(i,j)*velocity(node,2)
     &               +smatAbdz(i,j)*velocity(node,3)
     &                 -smatDx(i,j)*(velocity(node,1)+diffGrad(node,1))
     &                 -smatDy(i,j)*(velocity(node,2)+diffGrad(node,2))
     &                 -smatDz(i,j)*(velocity(node,3)+diffGrad(node,3))
     &                )
            END DO
          END DO
        ELSE ! neznanka je fluks (CE JE NEZNANKA FLUKS ODSTEVAMO DRUGO ENACBO !!!!)
c       dodajamo prispevke enack,  ki se ODSTEJEJO
          DO isuma=1,mesh%sqQlistNO(inode) ! to grem maksimalno do dva
            i=mesh%sqQlist(inode,isuma)
            ic=mesh%sqQlistIC(inode,isuma)
            IF (isuma.EQ.2) THEN
              predznak=-1.0D0 ! drugo enacbo odstejemo od prve
            ELSE
              predznak=+1.0D0
            END IF
            DO j=1,mesh%npoc
              node=mesh%idc(ic,j)
              r(isqr)=r(isqr)+predznak*(
C                    time derivative + sources
     &               -smatB(i,j)*rhsv(node)
     &               -u(node)*(
C                    advection
     &               +smatAbdx(i,j)*velocity(node,1)
     &               +smatAbdy(i,j)*velocity(node,2)
     &               +smatAbdz(i,j)*velocity(node,3)
     &                 -smatDx(i,j)*(velocity(node,1)+diffGrad(node,1))
     &                 -smatDy(i,j)*(velocity(node,2)+diffGrad(node,2))
     &                 -smatDz(i,j)*(velocity(node,3)+diffGrad(node,3))
     &                ))
            END DO
          END DO
        END IF
      END DO

      cpu%time(21)=cpu%time(21)+cptime(rcpu)
      rcpu=cptime(0.)
C
C     solve system of linear equations
C
      CALL SolvSLE(inp%sqrs_type,inp%sqrs_prec,2,inp%sqrs_maxit,5,inp%dlse(8),sysm%neq,sysm%nnz,
     &                   prec,sysm,r,x,nit,xcpu,ierr)

      cpu%time(22)=cpu%time(22)+cptime(rcpu)

      IF (ierr.NE.0) CALL WarnErr(env,io,inp,4,"SolveFDdaeVD2sq","NAPAKA V SOVLERJU",ierr)
C
C     Copy solution vector to q and u using under-relaxation
c     under-relaxation should be 1.0 for linear problems
C
      DO i=1,mesh%nnodes
        IF (mesh%DCkode(i).LT.0) THEN
          u(i)=(1.0D0-inp%urDC)*u(i)+inp%urDC*x(ABS(mesh%DCkode(i)))
        END IF
      END DO
      DO i=1,mesh%nq
        IF (mesh%DCqkode(i).LT.0) THEN
          qu(i)=x(ABS(mesh%DCqkode(i)))
c          qu(i)=(1.0D0-inp%urDC)*qu(i)+inp%urDC*x(ABS(mesh%DCqkode(i)))
        END IF
      END DO


      DEALLOCATE (x,b,r)

      END



C -----------------------------------------------------------------------------
C
C       Solve diffusion advection equation using finite difference discretization
C       of time derivative where diffusivity varies in space
C       overdeterminet system of equations
C
      SUBROUTINE SolveFDdaeVD2(env,io,inp,cpu,mesh,sysm,rhsm,precv,
     &                      u,qu,velocity,rhsv,smatB,
     &                      smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nit,
     &                      diff,diffGrad)
C
C     $: resi H*u=G*q+B*rhsv-advekcija
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(matrix) :: sysm,rhsm
      TYPE(IOtype) io
      TYPE(inputtype) inp
      TYPE(CPUtype) cpu
      TYPE(penv) :: env

      INTEGER i,ic,it,j,nit,ierr
      INTEGER node

      REAL(8) rhsv(mesh%nnodes)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc)
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc)

      REAL(8) precv(mesh%DCnunk)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) u(mesh%nnodes)
      REAL(8) qu(mesh%nq)
      REAL(8) diff(mesh%nnodes)
      REAL(8) diffGrad(mesh%nnodes,3)


      REAL(8), ALLOCATABLE :: x(:), b(:), r(:)
      REAL(4) cptime,rcpu

      ALLOCATE (x(mesh%DCnunk),b(mesh%DCnb),r(mesh%nicnsp))

      rcpu=cptime(0.)
C
C     Set up rhs vector and initial approximation of the unknown vector
C
      DO i=1,mesh%nnodes
        IF (mesh%DCkode(i).LT.0) THEN
          x(ABS(mesh%Dckode(i)))=u(i)
        ELSE
          b(mesh%DCkode(i))=u(i)
        END IF
      END DO
      DO i=1,mesh%nq
        IF (mesh%Dcqkode(i).LT.0) THEN
          x(ABS(mesh%Dcqkode(i)))=qu(i)
        ELSE
          b(mesh%DCqkode(i))=qu(i)
        END IF
      END DO
C
C     r = rhs * b
C
      CALL CRSxV(rhsm,b,mesh%DCnb,r)
c
c     r = r - B * rhsv
c
      DO j=1,mesh%npoc
        i=0
        DO ic=1,mesh%nicell
          DO it=1,51 ! 51 izvornih tock v vsaki celici
            i=i+1  ! vrstica v matriki
              node=mesh%idc(ic,j)
              r(i)=r(i)
C                    time derivative + sources
     &               -smatB(i,j)*rhsv(node)
     &               -u(node)*(
C                    advection
     &               +smatAbdx(i,j)*velocity(node,1)
     &               +smatAbdy(i,j)*velocity(node,2)
     &               +smatAbdz(i,j)*velocity(node,3)
     &                 -smatDx(i,j)*(velocity(node,1)+diffGrad(node,1))
     &                 -smatDy(i,j)*(velocity(node,2)+diffGrad(node,2))
     &                 -smatDz(i,j)*(velocity(node,3)+diffGrad(node,3))
     &                )
          END DO
        END DO
      END DO
      cpu%time(21)=cpu%time(21)+cptime(rcpu)
      rcpu=cptime(0.)
C
C     solve overdetermined system of linear equations
C
      CALL slvlsqr2(sysm%neq,mesh%DCnunk,sysm%nnz,inp%lsqs_maxit,inp%dlse(8),nit,ierr,
     &              precv,sysm%i,sysm%j,sysm%v,r,x)


      cpu%time(22)=cpu%time(22)+cptime(rcpu)

      IF (ierr.NE.0) CALL WarnErr(env,io,inp,4,"SolveFDdaeVD","NAPAKA V SOVLERJU",ierr)
C
C     Copy solution vector to q and u using under-relaxation
c     under-relaxation should be 1.0 for linear problems
C
      DO i=1,mesh%nnodes
        IF (mesh%DCkode(i).LT.0) THEN
          u(i)=(1.0D0-inp%urDC)*u(i)+inp%urDC*x(ABS(mesh%DCkode(i)))
        END IF
      END DO
      DO i=1,mesh%nq
        IF (mesh%DCqkode(i).LT.0) THEN
          qu(i)=(1.0D0-inp%urDC)*qu(i)+inp%urDC*x(ABS(mesh%DCqkode(i)))
        END IF
      END DO


      DEALLOCATE (x,b,r)

      END




C -----------------------------------------------------------------------------
C
C       Solve diffusion advection equation using finite difference discretization
C       of time derivative where diffusivity varies in space
C
C
      SUBROUTINE SolveFDdaeVD(env,io,inp,cpu,mesh,sysm,rhsm,precv,
     &                      u,qu,velocity,rhsv,smatB,
     &                      smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nit,
     &                      diff,diffGrad,diffLap)
C
C     $: resi H*u=G*q+B*rhsv-advekcija
C
C -----------------------------------------------------------------------------
      USE inc_types   

      TYPE(meshType) :: mesh
      TYPE(matrix) :: sysm,rhsm
      TYPE(IOtype) io
      TYPE(inputtype) inp
      TYPE(CPUtype) cpu     
      TYPE(penv) :: env                     
      
      INTEGER i,ic,it,j,nit,ierr
      INTEGER node
      
      REAL(8) rhsv(mesh%nnodes)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc)
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z   
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc) 
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc) 
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc) 
               
      REAL(8) precv(mesh%DCnunk)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) u(mesh%nnodes)
      REAL(8) qu(mesh%nq)   
      REAL(8) diff(mesh%nnodes)
      REAL(8) diffGrad(mesh%nnodes,3)
      REAL(8) diffLap(mesh%nnodes)
  

      REAL(8), ALLOCATABLE :: x(:), b(:), r(:)
      REAL(4) cptime,rcpu      
      
      ALLOCATE (x(mesh%DCnunk),b(mesh%DCnb),r(mesh%nicnsp))

      rcpu=cptime(0.)
C
C     Set up rhs vector and initial approximation of the unknown vector
C      
      DO i=1,mesh%nnodes
        IF (mesh%DCkode(i).LT.0) THEN
          x(ABS(mesh%Dckode(i)))=u(i)
        ELSE
          b(mesh%DCkode(i))=u(i)                
        END IF
      END DO      
      DO i=1,mesh%nq
        IF (mesh%Dcqkode(i).LT.0) THEN
          x(ABS(mesh%Dcqkode(i)))=qu(i)
        ELSE
          b(mesh%DCqkode(i))=qu(i)        
        END IF
      END DO
C
C     r = rhs * b
C      
      CALL CRSxV(rhsm,b,mesh%DCnb,r)      
c
c     r = r - B * rhsv
c
      DO j=1,mesh%npoc
        i=0
        DO ic=1,mesh%nicell
          DO it=1,51 ! 51 izvornih tock v vsaki celici
            i=i+1  ! vrstica v matriki      
              node=mesh%idc(ic,j)
              r(i)=r(i)
C                    time derivative + sources              
     &               -smatB(i,j)*rhsv(node)/diff(node)
     &               -u(node)/diff(node)*(
C                    advection 
     &               +smatAbdx(i,j)*(velocity(node,1)-diffGrad(node,1))
     &               +smatAbdy(i,j)*(velocity(node,2)-diffGrad(node,2))
     &               +smatAbdz(i,j)*(velocity(node,3)-diffGrad(node,3))
     &                 -smatDx(i,j)*(velocity(node,1)-diffGrad(node,1))
     &                 -smatDy(i,j)*(velocity(node,2)-diffGrad(node,2))
     &                 -smatDz(i,j)*(velocity(node,3)-diffGrad(node,3))
c                    variable diffusion
     &               +smatB(i,j)/diff(node)*(
c                    (v-grad) sklarno grad     
     &               +(  (velocity(node,1)-diffGrad(node,1))*diffGrad(node,1)
     &                  +(velocity(node,2)-diffGrad(node,2))*diffGrad(node,2)
     &                  +(velocity(node,3)-diffGrad(node,3))*diffGrad(node,3)
     &                ))     
c                    lapDiff
     &               +smatB(i,j)*diffLap(node)     
     &                )     
          END DO
        END DO
      END DO
      cpu%time(21)=cpu%time(21)+cptime(rcpu)
      rcpu=cptime(0.)
     
C
C     solve overdetermined system of linear equations
C      
      CALL slvlsqr2(sysm%neq,mesh%DCnunk,sysm%nnz,inp%lsqs_maxit,inp%dlse(8),nit,ierr,
     &              precv,sysm%i,sysm%j,sysm%v,r,x)
     
      cpu%time(22)=cpu%time(22)+cptime(rcpu)

      IF (ierr.NE.0) CALL WarnErr(env,io,inp,4,"SolveFDdaeVD","NAPAKA V SOVLERJU",ierr)
C
C     Copy solution vector to q and u using under-relaxation
c     under-relaxation should be 1.0 for linear problems
C      
      DO i=1,mesh%nnodes
        IF (mesh%DCkode(i).LT.0) THEN
          u(i)=(1.0D0-inp%urDC)*u(i)+inp%urDC*x(ABS(mesh%DCkode(i)))
        END IF
      END DO      
      DO i=1,mesh%nq
        IF (mesh%DCqkode(i).LT.0) THEN
          qu(i)=(1.0D0-inp%urDC)*qu(i)+inp%urDC*x(ABS(mesh%DCqkode(i)))
        END IF
      END DO


      DEALLOCATE (x,b,r)
      
      END




C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsB_DC_fill2sq(mesh,smatH,smatG,smatB,beta,diff,qdiff,sysm,rhsm,
     &                      velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)
C
C     $: Iz pravokotnih matrik g, h in b naredi CRS sistemsko in rhs, verzija brez laplace
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      REAL(8) smatH(mesh%nicnsp,mesh%npoc),smatG(mesh%nicnsp,mesh%npofc)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc),beta
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) diff(mesh%nnodes)
      REAL(8) qdiff(mesh%nq)
      REAL(8) val
      TYPE(matrix) :: sysm,rhsm

      INTEGER wqkode(mesh%nq)
      INTEGER wkode(mesh%nnodes)
      INTEGER nunk, nb

      INTEGER ic,row,col,ii,jj,nuq,j

      INTEGER, ALLOCATABLE :: eql(:)
      REAL(8), ALLOCATABLE :: FRsys(:),FRrhs(:)
      REAL(8) predznak
      INTEGER inode,isqr,isuma,i,kr,ks



      wqkode=mesh%DCqkode
      wkode=mesh%DCkode
      nunk=mesh%Dcnunk
      nb=mesh%Dcnb


      ALLOCATE (eql(nunk))
      ALLOCATE (FRsys(nunk),FRrhs(nb))
      eql=mesh%DCeql


      kr=0
      ks=0
      DO isqr=1,nunk ! zanka po vrsticah kvadratne matrike
        FRsys=0.0D0
        FRrhs=0.0D0
        inode=ABS(eql(isqr))  ! vozlisce v U ali Q mrezi, kateremu pripada ta enacba
        IF (eql(isqr).GT.0) THEN ! to pomeni, da je neznanka funkcija
c         dodajamo prispevke enack, ki se sestejejo
          DO isuma=1,mesh%sqUlistNO(inode)
            row=mesh%sqUlist(inode,isuma)
            ic=mesh%sqUlistIC(inode,isuma)

C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            val=diff(j)*smatH(row,col)+smatB(row,col)*beta
            nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              FRsys(ABS(nuq))=FRsys(ABS(nuq))+val
            ELSE ! znana vrednost - rhs
              FRrhs(nuq)=FRrhs(nuq)-val  ! minus zato, ker gre na drugo stran enacbe
            END IF
          END DO
C
C         nato G matrika (q)
C
          DO ii=1,mesh%nside
            DO jj=1,mesh%npof
              col=(ii-1)*mesh%npof+jj
              j=mesh%ibf(ic,ii,jj) ! stolpec za G
              nuq=wqkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                FRsys(ABS(nuq))=FRsys(ABS(nuq))-smatG(row,col)*qDiff(j)
              ELSE ! znana vrednost - rhs
                FRrhs(nuq)=FRrhs(nuq)+smatG(row,col)*qDiff(j)
              END IF
            END DO
          END DO


          END DO
        ELSE ! neznanka je fluks (CE JE NEZNANKA FLUKS ODSTEVAMO DRUGO ENACBO !!!!)
c         dodajamo prispevke enack, ki se ODSTEJEJO
          DO isuma=1,mesh%sqQlistNO(inode) ! to grem maksimalno do dva
            row=mesh%sqQlist(inode,isuma)
            ic=mesh%sqQlistIC(inode,isuma)
            IF (isuma.EQ.2) THEN
              predznak=-1.0D0 ! drugo enacbo odstejemo od prve
            ELSE
              predznak=+1.0D0
            END IF
C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            val=diff(j)*smatH(row,col)+smatB(row,col)*beta
            nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              FRsys(ABS(nuq))=FRsys(ABS(nuq))+val*predznak
            ELSE ! znana vrednost - rhs
              FRrhs(nuq)=FRrhs(nuq)-val*predznak  ! minus zato, ker gre na drugo stran enacbe
            END IF
          END DO
C
C         nato G matrika (q)
C
          DO ii=1,mesh%nside
            DO jj=1,mesh%npof
              col=(ii-1)*mesh%npof+jj
              j=mesh%ibf(ic,ii,jj) ! stolpec za G
              nuq=wqkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                FRsys(ABS(nuq))=FRsys(ABS(nuq))-smatG(row,col)*qDiff(j)*predznak
              ELSE ! znana vrednost - rhs
                FRrhs(nuq)=FRrhs(nuq)+smatG(row,col)*qDiff(j)*predznak
              END IF
            END DO
          END DO


          END DO
        END IF
C
C         predelam celo vrsto v CRS
C
C         sys
          DO i=sysm%i(isqr),sysm%i(isqr+1)-1
            ks=ks+1
            sysm%v(ks)=FRsys(sysm%j(ks))
          END DO
C         rhs
          DO i=rhsm%i(isqr),rhsm%i(isqr+1)-1
            kr=kr+1
            rhsm%v(kr)=FRrhs(rhsm%j(kr))
          END DO

      END DO

      DEALLOCATE (eql)

      END







C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsB_DC_fill(mesh,smatH,smatG,smatB,beta,diff,sysm,rhsm,
     &                      velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)      
C
C     $: Iz pravokotnih matrik g, h in b naredi CRS sistemsko in rhs
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      REAL(8) smatH(mesh%nicnsp,mesh%npoc),smatG(mesh%nicnsp,mesh%npofc)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc),beta
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z   
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc) 
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc) 
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc)   
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) diff(mesh%nnodes)
      REAL(8) val   
      TYPE(matrix) :: sysm,rhsm
      
      INTEGER wqkode(mesh%nq)
      INTEGER wkode(mesh%nnodes)
      INTEGER nunk, nb
            
      INTEGER ic,it,row,col,ii,jj,nuq,j,is,ir
      
      wqkode=mesh%DCqkode
      wkode=mesh%DCkode
      nunk=mesh%Dcnunk
      nb=mesh%Dcnb

      is=0
      ir=0
      DO ic=1,mesh%nicell ! zanka po celicah
        DO it=1,mesh%nsp ! po izvornih tockah znotraj celice
          row=(ic-1)*mesh%nsp+it ! vrstica v sistemski matriki
C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            val=smatH(row,col)+smatB(row,col)*beta/diff(j)         
            nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              is=is+1
              sysm%v(sysm%d(is))=val
            ELSE ! znana vrednost - rhs
              ir=ir+1
              rhsm%v(rhsm%d(ir))=-val  ! minus zato, ker gre na drugo stran enacbe
            END IF
          END DO            
C
C         nato G matrika (q)
C
          DO ii=1,mesh%nside
            DO jj=1,mesh%npof
              col=(ii-1)*mesh%npof+jj
              j=mesh%ibf(ic,ii,jj) ! stolpec za G
              nuq=wqkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                is=is+1
                sysm%v(sysm%d(is))=-smatG(row,col)
              ELSE ! znana vrednost - rhs
                ir=ir+1
                rhsm%v(rhsm%d(ir))=smatG(row,col) ! minus zato, ker gre na drugo stran enacbe
              END IF
            END DO
          END DO
        END DO
      END DO
      
      END    


C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsB_DC_fill2(mesh,smatH,smatG,smatB,beta,diff,qdiff,sysm,rhsm,
     &                      velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)
C
C     $: Iz pravokotnih matrik g, h in b naredi CRS sistemsko in rhs, verzija brez laplace
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      REAL(8) smatH(mesh%nicnsp,mesh%npoc),smatG(mesh%nicnsp,mesh%npofc)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc),beta
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) diff(mesh%nnodes)
      REAL(8) qdiff(mesh%nq)
      REAL(8) val
      TYPE(matrix) :: sysm,rhsm

      INTEGER wqkode(mesh%nq)
      INTEGER wkode(mesh%nnodes)
      INTEGER nunk, nb

      INTEGER ic,it,row,col,ii,jj,nuq,j,is,ir

      wqkode=mesh%DCqkode
      wkode=mesh%DCkode
      nunk=mesh%Dcnunk
      nb=mesh%Dcnb

      is=0
      ir=0
      DO ic=1,mesh%nicell ! zanka po celicah
        DO it=1,mesh%nsp ! po izvornih tockah znotraj celice
          row=(ic-1)*mesh%nsp+it ! vrstica v sistemski matriki
C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            val=diff(j)*smatH(row,col)+smatB(row,col)*beta
            nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              is=is+1
              sysm%v(sysm%d(is))=val
            ELSE ! znana vrednost - rhs
              ir=ir+1
              rhsm%v(rhsm%d(ir))=-val  ! minus zato, ker gre na drugo stran enacbe
            END IF
          END DO
C
C         nato G matrika (q)
C
          DO ii=1,mesh%nside
            DO jj=1,mesh%npof
              col=(ii-1)*mesh%npof+jj
              j=mesh%ibf(ic,ii,jj) ! stolpec za G
              nuq=wqkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                is=is+1
                sysm%v(sysm%d(is))=-smatG(row,col)*qDiff(j)
              ELSE ! znana vrednost - rhs
                ir=ir+1
                rhsm%v(rhsm%d(ir))=smatG(row,col)*qDiff(j) ! minus zato, ker gre na drugo stran enacbe
              END IF
            END DO
          END DO
        END DO
      END DO

      END




C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsB_SQsetup(eqn,mesh)
C
C     $: Priprava za kvadratno matriko
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh

      INTEGER wqkode(mesh%nq)
      INTEGER wkode(mesh%nnodes)
      INTEGER nunk, nb, eqn
      INTEGER, ALLOCATABLE :: eql(:)

      INTEGER i,j

      IF (eqn.LE.3) THEN
        wqkode=mesh%wqkode(:,eqn)
        wkode=mesh%wkode(:,eqn)
        nunk=mesh%nunk(eqn)
        nb=mesh%nb(eqn)
      ELSE IF (eqn.LE.4) THEN
        wqkode=mesh%Tqkode
        wkode=mesh%Tkode
        nunk=mesh%Tnunk
        nb=mesh%Tnb
      ELSE
        wqkode=mesh%DCqkode
        wkode=mesh%DCkode
        nunk=mesh%DCnunk
        nb=mesh%DCnb
      END IF

      ALLOCATE (eql(nunk))


c     v EQL(i) je stevilka vozlisca, kateremu pripada i-ta enacba
c     zanka po vrstica v matriki
      DO i=1,nunk
        DO j=1,mesh%nnodes
          IF (wkode(j).EQ.-i) THEN
            eql(i)=j
          END IF
        END DO
        DO j=1,mesh%nq
          IF (wqkode(j).EQ.-i) THEN
            eql(i)=-j  ! pri fluksih z minus da vem da so fluks izvorne tocke
          END IF
        END DO
      END DO


      IF (eqn.LE.3) THEN
        IF (eqn.EQ.1) ALLOCATE (mesh%Weql(nunk,3)) ! najprej klici z 1 !!!!!!!!!!!!!!!!!
        mesh%Weql(:,eqn)=eql
      ELSE IF (eqn.LE.4) THEN
        ALLOCATE (mesh%Teql(nunk))
        mesh%Teql=eql
      ELSE
        ALLOCATE (mesh%DCeql(nunk))
        mesh%DCeql=eql
      END IF

      DEALLOCATE (eql)

      END



C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsB_SQcount(eqn,mesh,sysm,rhsm)
C
C     $: Iz pravokotnih matrik g, h in b naredi CRS sistemsko in rhs
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      INTEGER, ALLOCATABLE :: FRsys(:),FRrhs(:)
      INTEGER, ALLOCATABLE :: eql(:)
      TYPE(matrix) :: sysm,rhsm

      INTEGER wqkode(mesh%nq)
      INTEGER wkode(mesh%nnodes)
      INTEGER nunk, nb, eqn

      INTEGER ic,row,col,ii,jj,nuq,j,i,ks,kr
      INTEGER inode,isqr,isuma

      IF (eqn.LE.3) THEN
        wqkode=mesh%wqkode(:,eqn)
        wkode=mesh%wkode(:,eqn)
        nunk=mesh%nunk(eqn)
        ALLOCATE (eql(nunk))
        eql=mesh%Weql(:,eqn)
        nb=mesh%nb(eqn)
      ELSE IF (eqn.LE.4) THEN
        wqkode=mesh%Tqkode
        wkode=mesh%Tkode
        nunk=mesh%Tnunk
        nb=mesh%Tnb
        ALLOCATE (eql(nunk))
        eql=mesh%Teql
      ELSE
        wqkode=mesh%DCqkode
        wkode=mesh%DCkode
        nunk=mesh%DCnunk
        ALLOCATE (eql(nunk))
        eql=mesh%DCeql
        nb=mesh%DCnb
      END IF

C
C     stevilo enacb
C
      sysm%neq=nunk
      rhsm%neq=nunk
      ALLOCATE (FRsys(nunk),FRrhs(nb))
C
C     stevilo nenicelnih clenov
C
      sysm%nnz=0
      rhsm%nnz=0
      DO isqr=1,nunk ! zanka po vrsticah kvadratne matrike
        inode=ABS(eql(isqr))  ! vozlisce v U ali Q mrezi, kateremu pripada ta enacba
        FRsys=0
        FRrhs=0
        IF (eql(isqr).GT.0) THEN ! to pomeni, da je neznanka funkcija
c         dodajamo prispevke enack, ki se sestejejo
          DO isuma=1,mesh%sqUlistNO(inode)
            row=mesh%sqUlist(inode,isuma)
            ic=mesh%sqUlistIC(inode,isuma)

C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              FRsys(ABS(nuq))=FRsys(ABS(nuq))+1
            ELSE ! znana vrednost - rhs
              FRrhs(nuq)=FRrhs(nuq)+1
            END IF
          END DO
C
C         nato G matrika (q)
C
          DO ii=1,mesh%nside
            DO jj=1,mesh%npof
              col=(ii-1)*mesh%npof+jj
              j=mesh%ibf(ic,ii,jj) ! stolpec za G
              nuq=wqkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                FRsys(ABS(nuq))=FRsys(ABS(nuq))+1
              ELSE ! znana vrednost - rhs
                FRrhs(nuq)=FRrhs(nuq)+1
              END IF
            END DO
          END DO


          END DO
        ELSE ! neznanka je fluks
c         dodajamo prispevke enack, ki se sestejejo
          DO isuma=1,mesh%sqQlistNO(inode)
            row=mesh%sqQlist(inode,isuma)
            ic=mesh%sqQlistIC(inode,isuma)

C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              FRsys(ABS(nuq))=FRsys(ABS(nuq))+1
            ELSE ! znana vrednost - rhs
              FRrhs(nuq)=FRrhs(nuq)+1
            END IF
          END DO
C
C         nato G matrika (q)
C
          DO ii=1,mesh%nside
            DO jj=1,mesh%npof
              col=(ii-1)*mesh%npof+jj
              j=mesh%ibf(ic,ii,jj) ! stolpec za G
              nuq=wqkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                FRsys(ABS(nuq))=FRsys(ABS(nuq))+1
              ELSE ! znana vrednost - rhs
                FRrhs(nuq)=FRrhs(nuq)+1
              END IF
            END DO
          END DO


          END DO
        END IF

c       prestejem nnz
        DO i=1,nunk
          IF (FRsys(i).NE.0) sysm%nnz=sysm%nnz+1
        END DO
        DO i=1,nb
          IF (FRrhs(i).NE.0) rhsm%nnz=rhsm%nnz+1
        END DO

      END DO
C
C     Alokacija pomnilnika
C
      ALLOCATE(sysm%v(sysm%nnz),sysm%i(sysm%neq+1),sysm%j(sysm%nnz),sysm%d(sysm%neq))
      ALLOCATE(rhsm%v(rhsm%nnz),rhsm%i(rhsm%neq+1),rhsm%j(rhsm%nnz),rhsm%d(rhsm%neq))
C
C     Polnjene matrik (napolnim polno vrsti, nato v CRS (da mi ni treba sortirat)
C
      ks=0
      kr=0
      DO isqr=1,nunk ! zanka po vrsticah kvadratne matrike
        inode=ABS(eql(isqr))  ! vozlisce v U ali Q mrezi, kateremu pripada ta enacba
        FRsys=0
        FRrhs=0
        IF (eql(isqr).GT.0) THEN ! to pomeni, da je neznanka funkcija
c         dodajamo prispevke enack, ki se sestejejo
          DO isuma=1,mesh%sqUlistNO(inode)
            row=mesh%sqUlist(inode,isuma)
            ic=mesh%sqUlistIC(inode,isuma)

C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              FRsys(ABS(nuq))=FRsys(ABS(nuq))+1
            ELSE ! znana vrednost - rhs
              FRrhs(nuq)=FRrhs(nuq)+1
            END IF
          END DO
C
C         nato G matrika (q)
C
          DO ii=1,mesh%nside
            DO jj=1,mesh%npof
              col=(ii-1)*mesh%npof+jj
              j=mesh%ibf(ic,ii,jj) ! stolpec za G
              nuq=wqkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                FRsys(ABS(nuq))=FRsys(ABS(nuq))+1
              ELSE ! znana vrednost - rhs
                FRrhs(nuq)=FRrhs(nuq)+1
              END IF
            END DO
          END DO


          END DO
        ELSE ! neznanka je fluks
c         dodajamo prispevke enack, ki se sestejejo
          DO isuma=1,mesh%sqQlistNO(inode)
            row=mesh%sqQlist(inode,isuma)
            ic=mesh%sqQlistIC(inode,isuma)

C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              FRsys(ABS(nuq))=FRsys(ABS(nuq))+1
            ELSE ! znana vrednost - rhs
              FRrhs(nuq)=FRrhs(nuq)+1
            END IF
          END DO
C
C         nato G matrika (q)
C
          DO ii=1,mesh%nside
            DO jj=1,mesh%npof
              col=(ii-1)*mesh%npof+jj
              j=mesh%ibf(ic,ii,jj) ! stolpec za G
              nuq=wqkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                FRsys(ABS(nuq))=FRsys(ABS(nuq))+1
              ELSE ! znana vrednost - rhs
                FRrhs(nuq)=FRrhs(nuq)+1
              END IF
            END DO
          END DO


          END DO
        END IF

C
C         predelam v CRS
C
C         zacetek vrste
          sysm%i(isqr)=ks+1
          rhsm%i(isqr)=kr+1
C         sys
          DO i=1,nunk
            IF (FRsys(i).NE.0) THEN
              ks=ks+1
              sysm%j(ks)=i
              IF (i.EQ.isqr) sysm%d(isqr)=ks ! postavimo diagonalo
            END IF
          END DO
C         rhs
          DO i=1,nb
            IF (FRrhs(i).NE.0) THEN
              kr=kr+1
c              rhsm%d(FRrhs(i))=kr
              rhsm%j(kr)=i
            END IF
          END DO

      END DO
c     zadnji clen, zacetek neobstojece vrste
      sysm%i(nunk+1)=ks+1
      rhsm%i(nunk+1)=kr+1

      DEALLOCATE (FRsys,FRrhs)

      END







C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsB_count(eqn,mesh,sysm,rhsm)
C
C     $: Iz pravokotnih matrik g, h in b naredi CRS sistemsko in rhs
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      INTEGER, ALLOCATABLE :: FRsys(:),FRrhs(:)
      TYPE(matrix) :: sysm,rhsm
      
      INTEGER wqkode(mesh%nq)
      INTEGER wkode(mesh%nnodes)
      INTEGER nunk, nb, eqn
            
      INTEGER ic,it,row,col,ii,jj,nuq,j,i,ks,kr,is,ir
      
      IF (eqn.LE.3) THEN
        wqkode=mesh%wqkode(:,eqn)
        wkode=mesh%wkode(:,eqn)
        nunk=mesh%nunk(eqn)
        nb=mesh%nb(eqn)
      ELSE IF (eqn.LE.4) THEN
        wqkode=mesh%Tqkode
        wkode=mesh%Tkode
        nunk=mesh%Tnunk
        nb=mesh%Tnb
      ELSE IF (eqn.LE.5) THEN
        wqkode=mesh%DCqkode
        wkode=mesh%DCkode
        nunk=mesh%DCnunk
        nb=mesh%DCnb 
      ELSE
        wqkode=mesh%mHqkode
        wkode=mesh%mHkode
        nunk=mesh%mHnunk
        nb=mesh%mHnb       
      END IF        
C      
C     stevilo enacb
C
      sysm%neq=mesh%nicnsp      
      rhsm%neq=mesh%nicnsp
C
C     stevilo nenicelnih clenov
C      
      sysm%nnz=0
      rhsm%nnz=0
      DO ic=1,mesh%nicell ! zanka po celicah
        DO it=1,mesh%nsp ! po izvornih tockah znotraj celice
          row=(ic-1)*mesh%nsp+it ! vrstica v sistemski matriki        
C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
              j=mesh%idc(ic,col) ! stolpec za H
              nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                sysm%nnz=sysm%nnz+1
              ELSE ! znana vrednost - rhs
                rhsm%nnz=rhsm%nnz+1
              END IF
          END DO            
C
C         nato G matrika (q)
C
          DO ii=1,mesh%nside
            DO jj=1,mesh%npof
              col=(ii-1)*mesh%npof+jj
              j=mesh%ibf(ic,ii,jj) ! stolpec za G
              nuq=wqkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                sysm%nnz=sysm%nnz+1
              ELSE ! znana vrednost - rhs
                rhsm%nnz=rhsm%nnz+1
              END IF
            END DO
          END DO
C                    
        END DO
      END DO
C
C     Alokacija pomnilnika, v sysm%d bom dal preslikavo za hitro polnjenje matrike
C
      ALLOCATE(sysm%v(sysm%nnz),sysm%i(sysm%neq+1),sysm%j(sysm%nnz),sysm%d(sysm%nnz))
      ALLOCATE(rhsm%v(rhsm%nnz),rhsm%i(rhsm%neq+1),rhsm%j(rhsm%nnz),rhsm%d(rhsm%nnz))
C
C     Polnjene matrik (napolnim polno vrsti, nato v CRS (da mi ni treba sortirat)
C         
      ALLOCATE (FRsys(nunk),FRrhs(nb))
      ks=0
      kr=0
      is=0
      ir=0
      DO ic=1,mesh%nicell ! zanka po celicah
        DO it=1,mesh%nsp ! po izvornih tockah znotraj celice
          row=(ic-1)*mesh%nsp+it ! vrstica v sistemski matriki
          FRsys=0
          FRrhs=0
C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              is=is+1
              FRsys(ABS(nuq))=is
            ELSE ! znana vrednost - rhs
              ir=ir+1
              FRrhs(nuq)=ir
            END IF
          END DO            
C
C         nato G matrika (q)
C
          DO ii=1,mesh%nside
            DO jj=1,mesh%npof
              col=(ii-1)*mesh%npof+jj
              j=mesh%ibf(ic,ii,jj) ! stolpec za G
              nuq=wqkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                is=is+1
                FRsys(ABS(nuq))=is
              ELSE ! znana vrednost - rhs
                ir=ir+1
                FRrhs(nuq)=ir
              END IF
            END DO
          END DO
C               
C         predelam v CRS
C     
C         zacetek vrste
          sysm%i(row)=ks+1
          rhsm%i(row)=kr+1          
C         sys
          DO i=1,nunk
            IF (FRsys(i).NE.0) THEN
              ks=ks+1
              sysm%d(FRsys(i))=ks
              sysm%j(ks)=i
            END IF
          END DO
C         rhs
          DO i=1,nb
            IF (FRrhs(i).NE.0) THEN
              kr=kr+1
              rhsm%d(FRrhs(i))=kr
              rhsm%j(kr)=i
            END IF
          END DO
c          
        END DO
      END DO     
c     zadnji clen, zacetek neobstojece vrste
      sysm%i(row+1)=ks+1
      rhsm%i(row+1)=kr+1          
       
      
      DEALLOCATE (FRsys,FRrhs)
      
      END




C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsBvDiff(mesh,smatH,smatG,smatB,beta,diff,sysm,rhsm)
C
C     $: Iz pravokotnih matrik g, h in b naredi CRS sistemsko in rhs
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      REAL(8) smatH(mesh%nicnsp,mesh%npoc),smatG(mesh%nicnsp,mesh%npofc)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc),beta,diff(mesh%nnodes)
      REAL(8), ALLOCATABLE :: FRsys(:),FRrhs(:)
      REAL(8) val
      TYPE(matrix) :: sysm,rhsm
      
      INTEGER wqkode(mesh%nq)
      INTEGER wkode(mesh%nnodes)
      INTEGER nunk, nb
            
      INTEGER ic,it,row,col,ii,jj,nuq,j,i,ks,kr
           
      wqkode=mesh%DCqkode
      wkode=mesh%DCkode
      nunk=mesh%DCnunk
      nb=mesh%DCnb
C      
C     stevilo enacb
C
      sysm%neq=mesh%nicnsp      
      rhsm%neq=mesh%nicnsp
C
C     stevilo nenicelnih clenov
C      
      sysm%nnz=0
      rhsm%nnz=0
      DO ic=1,mesh%nicell ! zanka po celicah
        DO it=1,mesh%nsp ! po izvornih tockah znotraj celice
          row=(ic-1)*mesh%nsp+it ! vrstica v sistemski matriki        
C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            val=smatH(row,col)+smatB(row,col)*beta/diff(j)
            IF (val.NE.0.0D0) THEN
              nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                sysm%nnz=sysm%nnz+1
              ELSE ! znana vrednost - rhs
                rhsm%nnz=rhsm%nnz+1
              END IF
            END IF
          END DO            
C
C         nato G matrika (q)
C
          DO ii=1,mesh%nside
            DO jj=1,mesh%npof
              col=(ii-1)*mesh%npof+jj
              IF (smatG(row,col).NE.0.0D0) THEN
                j=mesh%ibf(ic,ii,jj) ! stolpec za G
                nuq=wqkode(j) ! stevilka v vektorju neznank oziroma znank
                IF (nuq.LT.0) THEN ! neznaka - sys
                  sysm%nnz=sysm%nnz+1
                ELSE ! znana vrednost - rhs
                  rhsm%nnz=rhsm%nnz+1
                END IF
              END IF
            END DO
          END DO
C                    
        END DO
      END DO
C
C     Alokacija pomnilnika
C
      ALLOCATE(sysm%v(sysm%nnz),sysm%i(sysm%neq+1),sysm%j(sysm%nnz))
      ALLOCATE(rhsm%v(rhsm%nnz),rhsm%i(rhsm%neq+1),rhsm%j(rhsm%nnz))   
C
C     Polnjene matrik (napolnim polno vrsti, nato v CRS (da mi ni treba sortirat)
C         
      ALLOCATE (FRsys(nunk),FRrhs(nb))
      ks=0
      kr=0
      DO ic=1,mesh%nicell ! zanka po celicah
        DO it=1,mesh%nsp ! po izvornih tockah znotraj celice
          row=(ic-1)*mesh%nsp+it ! vrstica v sistemski matriki
          FRsys=0.0D0
          FRrhs=0.0D0
C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              FRsys(ABS(nuq))=smatH(row,col)+smatB(row,col)*beta/diff(j)
            ELSE ! znana vrednost - rhs
              FRrhs(nuq)=-smatH(row,col)-smatB(row,col)*beta/diff(j)  ! minus zato, ker gre na drugo stran enacbe
            END IF
          END DO            
C
C         nato G matrika (q)
C
          DO ii=1,mesh%nside
            DO jj=1,mesh%npof
              col=(ii-1)*mesh%npof+jj
              j=mesh%ibf(ic,ii,jj) ! stolpec za G
              nuq=wqkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                FRsys(ABS(nuq))=-smatG(row,col)
              ELSE ! znana vrednost - rhs
                FRrhs(nuq)=smatG(row,col) ! minus zato, ker gre na drugo stran enacbe
              END IF
            END DO
          END DO
C               
C         predelam v CRS
C     
C         zacetek vrste
          sysm%i(row)=ks+1
          rhsm%i(row)=kr+1          
C         sys
          DO i=1,nunk
            IF (FRsys(i).NE.0.0D0) THEN
              ks=ks+1
              sysm%v(ks)=FRsys(i)
              sysm%j(ks)=i
            END IF
          END DO
C         rhs
          DO i=1,nb
            IF (FRrhs(i).NE.0.0D0) THEN
              kr=kr+1
              rhsm%v(kr)=FRrhs(i)
              rhsm%j(kr)=i
            END IF
          END DO
c          
        END DO
      END DO     
c     zadnji clen, zacetek neobstojece vrste
      sysm%i(row+1)=ks+1
      rhsm%i(row+1)=kr+1          
       
      
      DEALLOCATE (FRsys,FRrhs)
      
      END
      
C -----------------------------------------------------------------------------      
      SUBROUTINE CalDiffGrad(mesh,DCdiff,DCdiffGrad,DCdiffLap)
c
c     Calculate diffusivity gradients
c      
C -----------------------------------------------------------------------------
      USE inc_types   

      TYPE(meshType) :: mesh      
      
      REAL(8) DCdiffGrad(mesh%nnodes,3) 
      REAL(8) DCdiffLap(mesh%nnodes)
      REAL(8) DCdiff(mesh%nnodes)       
      
      REAL(8), ALLOCATABLE :: tmp(:,:)
      
      
      CALL setGrad(mesh,DCdiff,DCdiffGrad)
      
      ALLOCATE (tmp(mesh%nnodes,3))
      
      CALL setGrad(mesh,DCdiffGrad(:,1),tmp)  ! grad (d(DCdiff)/dx))
      DCdiffLap=tmp(:,1)

      CALL setGrad(mesh,DCdiffGrad(:,2),tmp)  ! grad (d(DCdiff)/dy))
      DCdiffLap=DCdiffLap+tmp(:,2)

      CALL setGrad(mesh,DCdiffGrad(:,3),tmp)  ! grad (d(DCdiff)/dz))
      DCdiffLap=DCdiffLap+tmp(:,3)
      
      
      DEALLOCATE (tmp)
      END
      
      
C -----------------------------------------------------------------------------
      SUBROUTINE SetUpDiffSources(mesh,inp,DCdiff,DCrhsv,DCtime,DCsource)
C
C     Set up concentration diffusivity and sources
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(InputType) inp
      TYPE(meshType) :: mesh
      REAL(8) DCrhsv(mesh%nnodes)    ! body force (time derivatives + sources)
      REAL(8) DCtime(mesh%nnodes)    ! body force (time derivatives)
      REAL(8) DCdiff(mesh%nnodes)    ! difuzivnost
      REAL(8) DCsource(mesh%nnodes)    ! sources
      INTEGER i
      REAL(8) S

      DO i=1,mesh%nnodes
         DCdiff(i)=1.0D0/(inp%ReN*inp%ScN)
c              DCdiff(i)=mesh%x(i,1)  ! diff = x (mesh med 1 in 2, doktorat Tabela 5.1) u(1.5)=0.5849625
c              DCdiff(i)=mesh%x(i,1)**2  ! diff = x^2 (mesh med 1 in 2, doktorat Tabela 5.1) u(1.5)=0.6666666
      END DO


C     set up sources DCrhsv(i)=DCtime(i)+S(i)
      DO i=1,mesh%nnodes
        DCrhsv(i)=DCtime(i)-DCsource(i)
      END DO

      END


C -----------------------------------------------------------------------------
      SUBROUTINE SetUpDCtest(env,io,mesh,inp,DCsource,anal,time,dc,qdc)
C
C     Set up concentration diffusivity and sources
C
C -----------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(InputType) inp
      TYPE(meshType) :: mesh
      TYPE(penv) :: env
      TYPE(IOtype)    :: io
      REAL(8) x,y,z,time,u,f,pi,dvapi
      INTEGER i,ii,dn,j
      REAL(8) DCsource(mesh%nnodes)    ! sources
      REAL(8) anal(mesh%nnodes)    ! analytic

      REAL(8) dc(mesh%nnodes)
      REAL(8) qdc(mesh%nq)

      REAL(8), ALLOCATABLE :: dudx(:),dudy(:),dudz(:)

      pi=4.0D0*ATAN(1.0D0)
      dvapi = 2.0D0 * pi


      IF (inp%iTest.EQ.0) THEN
        DCsource=0.0D0
        RETURN
      END IF

      ALLOCATE (dudx(mesh%nq),dudy(mesh%nq),dudz(mesh%nq))
      dudx=0.0D0 ! za tiste, kjer ne podam
      dudy=0.0D0
      dudz=0.0D0

c     add additional sources
      DO i=1,mesh%nnodes
         x=mesh%x(i,1)
         y=mesh%x(i,2)
         z=mesh%x(i,3)
         u=dc(i)

         IF (inp%itest.EQ.1) THEN ! 1D (difcon-Dirichlet-1D.bic)
           DCsource(i)=-2.0D0
           anal(i)=x**2.0d0+sin(x)*sin(time)

         ELSE IF (inp%itest.EQ.10) THEN ! 3D (difcon-Dirichlet-3D.bic)
            f = 1.0D0 - time - x - y - z + (time + x + y + z)**3 
            DCsource(i)=  f  + u - u**3
            anal(i)=x+y+z+time

          ELSE IF (inp%itest.EQ.11) THEN ! 1D (difcon-Dirichlet-1D.bic)
            f = exp(-3.0D0*time)*sin(pi*x)*(exp(2.0D0*time)*(-2.0D0+pi**2)+ (sin(pi*x))**2  )
            DCsource(i)=  f  + u - u**3
            anal(i)=exp(-time)*sin(pi*x)

          ELSE IF (inp%itest.EQ.12) THEN ! 1D (difcon-Dirichlet-1D.bic)
            f = exp(-3.0D0*time)*sin(dvapi*x)*(2.0D0*exp(2.0D0*time)*(-1.0D0+2.0D0*pi**2)+ (sin(dvapi*x))**2  )
            DCsource(i)=  f  + u - u**3
            anal(i)=exp(-time)*sin(dvapi*x)


          ELSE IF (inp%itest.EQ.13) THEN ! 2D (allenCahn-test13.bic)
            f = exp(-3.0D0*time)*sin(pi*x)*sin(pi*y) *
     &       (2.0D0*exp(2.0D0*time)*(-1.0D0+pi**2)+ (sin(pi*x))**2 * (sin(pi*y))**2 )
            DCsource(i)=  f  + u - u**3
            anal(i)=exp(-time)*sin(pi*x)*sin(pi*y)


          ELSE IF (inp%itest.EQ.14) THEN ! 2D (allenCahn-test13.bic)
            f = exp(-3.0D0*time)*sin(dvapi*x)*sin(dvapi*y) *
     &       (2.0D0*exp(2.0D0*time)*(-1.0D0+4.0D0*pi**2)+ (sin(pi*x))**2 * (sin(pi*y))**2 )
            DCsource(i)=  f  + u - u**3
            anal(i)=exp(-time)*sin(dvapi*x)*sin(dvapi*y)


          ELSE IF (inp%itest.EQ.15) THEN ! 3D (allenCahn-test16.bic)
            f = exp(-3.0D0*time)*sin(pi*x)*sin(pi*y)*sin(pi*z) *
     &       (exp(2.0D0*time)*(-2.0D0+3.0D0*pi**2)+ (sin(pi*x))**2 * (sin(pi*y))**2 * (sin(pi*z))**2)
            DCsource(i)=  f  + u - u**3
            anal(i)=exp(-time)*sin(pi*x)*sin(pi*y)*sin(pi*z)


          ELSE IF (inp%itest.EQ.16) THEN ! 3D (allenCahn-test16.bic)
            f = exp(-3.0D0*time)*sin(dvapi*x)*sin(dvapi*y)*sin(dvapi*z) *
     &       (2.0D0*exp(2.0D0*time)*(-1.0D0+6.0D0*pi**2)+ (sin(dvapi*x))**2 * (sin(dvapi*y))**2 * (sin(dvapi*z))**2)
            DCsource(i)=  f  + u - u**3
            anal(i)=exp(-time)*sin(dvapi*x)*sin(dvapi*y)*sin(dvapi*z)




          ELSE IF (inp%itest.EQ.6) THEN ! 3D (difcon-Dirichlet-3D.bic)
           DCsource(i)=(x*y*z*time)**2-dc(i)**2
           anal(i)=x*y*z*time
         ELSE IF (inp%itest.EQ.3) THEN ! 2D (difcon-Dirichlet-2D.bic)
           DCsource(i)=0.0D0
           anal(i)=exp(x+y-sqrt(2.0D0)*time)
         ELSE IF (inp%itest.EQ.4) THEN ! 2D (difcon-Dirichlet-2D.bic)
           DCsource(i)=(x*y*time)**2-dc(i)**2
           anal(i)=x*y*time
         ELSE IF (inp%itest.EQ.2) THEN ! 1D (difcon-Dirichlet-1D.bic)
           DCsource(i)=dc(i)+dc(i)**2-x*time-(x*time)**2
           anal(i)=x*time
         ELSE IF (inp%itest.EQ.5) THEN ! 3D (difcon-Dirichlet-3D.bic)
           DCsource(i)=sin(x)+sin(y)
           anal(i)=sin(x)+sin(y)+sin(z)*sin(time)
         ELSE IF (inp%itest.NE.0) THEN
           CALL WarnErr(env,io,inp,2,'SetUpDCtest','Unknown test!',inp%itest)
         END IF
      END DO


c     Consider fluxes
      DO i=1,mesh%nbelem
        DO j=1,mesh%npof       
          dn=mesh%ibcf(i,j)
          x=mesh%xq(dn,1)
          y=mesh%xq(dn,2)
          z=mesh%xq(dn,3)
        
          IF (inp%itest.EQ.13) THEN ! 2D (allenCahn-test13.bic)
            dudx(dn)=pi*exp(-time)*sin(pi*y)*cos(pi*x)
          END IF            
          IF (inp%itest.EQ.14) THEN ! 2D (allenCahn-test13.bic)
              dudx(dn)=dvapi*exp(-time)*sin(dvapi*y)*cos(dvapi*x)
          END IF
          IF (inp%itest.EQ.15) THEN ! 3D (allenCahn-test16-mixed.bic)
            dudx(dn)=pi*exp(-time)*cos(pi*x)*sin(pi*y)*sin(pi*z)
            dudy(dn)=pi*exp(-time)*sin(pi*x)*cos(pi*y)*sin(pi*z)
            dudz(dn)=pi*exp(-time)*sin(pi*x)*sin(pi*y)*cos(pi*z)
          END IF            
          IF (inp%itest.EQ.16) THEN !  3D (allenCahn-test16-mixed.bic)
              dudx(dn)=dvapi*exp(-time)*cos(dvapi*x)*sin(dvapi*y)*sin(dvapi*z)
              dudy(dn)=dvapi*exp(-time)*sin(dvapi*x)*cos(dvapi*y)*sin(dvapi*z)
              dudz(dn)=dvapi*exp(-time)*sin(dvapi*x)*sin(dvapi*y)*cos(dvapi*z)
          END IF



        END DO
      END DO

C
C     Set up boundary conditions
C

C     Dirichlet, kode=0
      DO i=1,mesh%nbnodes
        dn=mesh%gbn(i)
        IF (mesh%DCkode(dn).GT.0) THEN ! to je znana vrednost, mesh%iwn(i) je stevilka stene
          dc(dn)=anal(dn)
        END IF
      END DO

C     Neumann, kode=2
      ii=0
      DO i=1,mesh%nbelem
        DO j=1,mesh%npof
          ii=ii+1
          dn=mesh%ibcf(i,j)
          IF (mesh%DCqkode(dn).GT.0) THEN ! to je znana vrednost, mesh%bewn(i) je stevilka stene
            IF (mesh%bewn(i).EQ.1) THEN ! stena pri z=1
              qdc(dn)=-dudz(dn)
            END IF
            IF (mesh%bewn(i).EQ.6) THEN ! stena pri z=2
              qdc(dn)=dudz(dn)
            END IF
            IF (mesh%bewn(i).EQ.5) THEN ! stena pri x=1
              qdc(dn)=-dudx(dn)
            END IF
            IF (mesh%bewn(i).EQ.3) THEN ! stena pri x=2
              qdc(dn)=dudx(dn)
            END IF
            IF (mesh%bewn(i).EQ.2) THEN ! stena pri y=1
              qdc(dn)=-dudy(dn)
            END IF
            IF (mesh%bewn(i).EQ.4) THEN ! stena pri y=2
              qdc(dn)=dudy(dn)
            END IF
c            fmQdc(ii)=qdc(dn)
          END IF
        END DO
      END DO

      DEALLOCATE (dudx,dudy,dudz)


      END


C -----------------------------------------------------------------------------
      SUBROUTINE DCinitial(env,io,mesh,inp,dc,ptsdc,pptsdc)
C
C     $: Set initial conditions
C
C -----------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(meshType) :: mesh
      TYPE(InputType) inp
      TYPE(IOtype)    :: io
      TYPE(penv) :: env
      REAL(8) dc(mesh%nnodes)
      REAL(8) ptsdc(mesh%nnodes)
      REAL(8) pptsdc(mesh%nnodes)

      INTEGER i
      REAL(8) x,y,z,pi,dvapi

      pi=4.0D0*ATAN(1.0D0)
      dvapi = 2.0D0 * pi

c     init
      DO i=1,mesh%nnodes
         x=mesh%x(i,1)
         y=mesh%x(i,2)
         z=mesh%x(i,3)

         IF (inp%itest.EQ.1) THEN ! 1D (difcon-Dirichlet-1D.bic)
           dc(i)=x**2.0D0
           ptsdc(i)=dc(i)-inp%tstep*sin(x)
           pptsdc(i)=0.0D0

          ELSE IF (inp%itest.EQ.10) THEN ! 3D (difcon-Dirichlet-3D.bic)
            dc(i)=x+y+z
            ptsdc(i)=dc(i)-inp%tstep
            pptsdc(i)=0.0D0

          ELSE IF (inp%itest.EQ.11) THEN ! 3D (difcon-Dirichlet-1D.bic)
            dc(i)=sin(pi*x)
            ptsdc(i)=dc(i) + inp%tstep * dc(i)
            pptsdc(i)=0.0D0

          ELSE IF (inp%itest.EQ.12) THEN ! 3D (difcon-Dirichlet-1D.bic)
            dc(i)=sin(dvapi*x)
            ptsdc(i)=dc(i) + inp%tstep * dc(i)
            pptsdc(i)=0.0D0


          ELSE IF (inp%itest.EQ.13) THEN ! 3D (allenCahn-test13.bic)
            dc(i)=sin(pi*x)*sin(pi*y)
            ptsdc(i)=dc(i) + inp%tstep * dc(i)
            pptsdc(i)=0.0D0


          ELSE IF (inp%itest.EQ.14) THEN ! 3D (allenCahn-test13.bic)
            dc(i)=sin(dvapi*x)*sin(dvapi*y)
            ptsdc(i)=dc(i) + inp%tstep * dc(i)
            pptsdc(i)=0.0D0


          ELSE IF (inp%itest.EQ.15) THEN ! 3D (allenCahn-test16.bic)
            dc(i)=sin(pi*x)*sin(pi*y)*sin(pi*z)
            ptsdc(i)=dc(i) + inp%tstep * dc(i)
            pptsdc(i)=0.0D0


          ELSE IF (inp%itest.EQ.16) THEN ! 3D (allenCahn-test16.bic)
            dc(i)=sin(dvapi*x)*sin(dvapi*y)*sin(dvapi*z)
            ptsdc(i)=dc(i) + inp%tstep * dc(i)
            pptsdc(i)=0.0D0


          ELSE IF (inp%itest.EQ.6) THEN ! 3D (difcon-Dirichlet-3D.bic)
           dc(i)=0.0D0
           ptsdc(i)=-inp%tstep*x*y*z
           pptsdc(i)=0.0D0
         ELSE IF (inp%itest.EQ.3) THEN ! 2D (difcon-Dirichlet-2D.bic)
           dc(i)=exp(x+y)
           ptsdc(i)=(1.0D0+inp%tstep*sqrt(2.0D0))*exp(x+y)
           pptsdc(i)=0.0D0
         ELSE IF (inp%itest.EQ.4) THEN ! 2D (difcon-Dirichlet-2D.bic)
           dc(i)=0.0D0
           ptsdc(i)=-inp%tstep*x*y
           pptsdc(i)=0.0D0
         ELSE IF (inp%itest.EQ.2) THEN ! 1D (difcon-Dirichlet-1D.bic)
           dc(i)=0.0D0
           ptsdc(i)=-inp%tstep*x
           pptsdc(i)=0.0D0
         ELSE IF (inp%itest.EQ.5) THEN ! 3D (difcon-Dirichlet-3D.bic)
           dc(i)=sin(x)+sin(y)
           ptsdc(i)=dc(i)-inp%tstep*sin(z)
           pptsdc(i)=0.0D0
         ELSE IF (inp%itest.NE.0) THEN
           CALL WarnErr(env,io,inp,2,'DCinitial','Unknown test!',inp%itest)
         END IF
      END DO

      END

C -----------------------------------------------------------------------------
      SUBROUTINE CopyTestFields(mesh,dc,anal,temp,vorticity)
C
C     $: copy fields for export and postprocessing
C
C -----------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(meshType) :: mesh

      REAL(8) dc(mesh%nnodes)
      REAL(8) anal(mesh%nnodes)
      REAL(8) temp(mesh%nnodes)
      REAL(8) vorticity(mesh%nnodes,3)

      INTEGER i

      DO i=1,mesh%nnodes
        temp(i)=anal(i)
        vorticity(i,1)=ABS(dc(i)-anal(i))
        IF (anal(i).NE.0.0D0) THEN
          vorticity(i,2)=ABS(dc(i)-anal(i))/ABS(anal(i))
        ELSE
          vorticity(i,2)=ABS(dc(i)-anal(i))
        END IF
      END DO

      END




C -----------------------------------------------------------------------------
      SUBROUTINE CalRMS(io,iTest,mesh,sol,anal,time)
C
C     $: Calculate RMS
C
C -----------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(meshType) :: mesh
      TYPE(IOtype)    :: io
      INTEGER iTest,i,n,lun,iDC,n1,n2,n3,TestType
      REAL(8) sol(mesh%nnodes),time
      REAL(8) anal(mesh%nnodes),vs,vsa,rms


      vs=0.0D0
      vsa=0.0D0
      n=0

      DO i=1,mesh%nnodes
        IF (mesh%lbn(i).EQ.0) THEN ! obmocni
          n=n+1
          vs=vs+(sol(i)-anal(i))**2
          vsa=vsa+(anal(i))**2
        END IF
      END DO

      rms=SQRT(vs/vsa)

      WRITE (io%dct,*) time,rms


      END

      
