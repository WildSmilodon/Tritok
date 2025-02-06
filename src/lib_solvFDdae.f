C -----------------------------------------------------------------------------      
C
C       Solve diffusion advection equation using finite difference discretization
C       of time derivative
C
C
      SUBROUTINE SolveFDdaeSQ(env,io,inp,cpu,mesh,sysm,rhsm,prec,
     &                      u,qu,velocity,rhsv,smatB,
     &                      smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nit,enazalfa)
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

c      REAL(8) precv(mesh%Tnunk)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) u(mesh%nnodes)
      REAL(8) qu(mesh%nq)
      REAL(8) enazalfa

      REAL(8), ALLOCATABLE :: x(:), b(:), r(:)
      REAL(4) cptime,rcpu,xcpu

      ALLOCATE (x(mesh%Tnunk),b(mesh%Tnb),r(mesh%Tnunk))

      rcpu=cptime(0.)
C
C     Set up rhs vector and initial approximation of the unknown vector
C
      DO i=1,mesh%nnodes
        IF (mesh%Tkode(i).LT.0) THEN
          x(ABS(mesh%Tkode(i)))=u(i)
        ELSE
          b(mesh%Tkode(i))=u(i)
        END IF
      END DO
      DO i=1,mesh%nq
        IF (mesh%Tqkode(i).LT.0) THEN
          x(ABS(mesh%Tqkode(i)))=qu(i)
        ELSE
          b(mesh%Tqkode(i))=qu(i)
        END IF
      END DO
C
C     r = rhs * b
C
      CALL CRSxV(rhsm,b,mesh%Tnb,r)
c
c     r = r - B * rhsv
c
      DO isqr=1,sysm%neq !nunk ! zanka po vrsticah kvadratne matrike
        inode=ABS(mesh%Teql(isqr))  ! vozlisce v U ali Q mrezi, kateremu pripada ta enacba
        IF (mesh%Teql(isqr).GT.0) THEN ! to pomeni, da je neznanka funkcija
c         dodajamo prispevke enack, ki se sestejejo
          DO isuma=1,mesh%sqUlistNO(inode)
            i=mesh%sqUlist(inode,isuma)
            ic=mesh%sqUlistIC(inode,isuma)
            DO j=1,mesh%npoc
              node=mesh%idc(ic,j)
              r(isqr)=r(isqr)
C                    time derivative
     &               -smatB(i,j)*rhsv(node)*enazalfa
     &               -enazalfa*(
C                    advection and vortex twisting and stretching
     &               +smatAbdx(i,j)*velocity(node,1)*u(node)
     &               +smatAbdy(i,j)*velocity(node,2)*u(node)
     &               +smatAbdz(i,j)*velocity(node,3)*u(node)
     &                 -smatDx(i,j)*velocity(node,1)*u(node)
     &                 -smatDy(i,j)*velocity(node,2)*u(node)
     &                 -smatDz(i,j)*velocity(node,3)*u(node)
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
C                    time derivative
     &               -smatB(i,j)*rhsv(node)*enazalfa
     &               -enazalfa*(
C                    advection and vortex twisting and stretching
     &               +smatAbdx(i,j)*velocity(node,1)*u(node)
     &               +smatAbdy(i,j)*velocity(node,2)*u(node)
     &               +smatAbdz(i,j)*velocity(node,3)*u(node)
     &                 -smatDx(i,j)*velocity(node,1)*u(node)
     &                 -smatDy(i,j)*velocity(node,2)*u(node)
     &                 -smatDz(i,j)*velocity(node,3)*u(node)
     &                ))
            END DO
          END DO
        END IF
      END DO




      cpu%time(19)=cpu%time(19)+cptime(rcpu)
      rcpu=cptime(0.)
C
C     solve system of linear equations
C
      CALL SolvSLE(inp%sqrs_type,inp%sqrs_prec,2,inp%sqrs_maxit,5,inp%dlse(7),sysm%neq,sysm%nnz,
     &                   prec,sysm,r,x,nit,xcpu,ierr)

      cpu%time(20)=cpu%time(20)+cptime(rcpu)

      IF (ierr.NE.0) CALL WarnErr(env,io,inp,4,"SolveFDwTE","NAPAKA V SOVLERJU",ierr)
C
C     Copy solution vector to q and u using under-relaxation
c     under-relaxation should be 1.0 for linear problems
C
      DO i=1,mesh%nnodes
        IF (mesh%Tkode(i).LT.0) THEN
          u(i)=(1.0D0-inp%urDT)*u(i)+inp%urDT*x(ABS(mesh%Tkode(i)))
        END IF
      END DO
      DO i=1,mesh%nq
        IF (mesh%Tqkode(i).LT.0) THEN
          qu(i)=x(ABS(mesh%Tqkode(i)))
        END IF
      END DO


      DEALLOCATE (x,b,r)

      END





C -----------------------------------------------------------------------------
C
C       Solve diffusion advection equation using finite difference discretization
C       of time derivative
C
C
      SUBROUTINE SolveFDdae(env,io,inp,cpu,mesh,sysm,rhsm,precv,
     &                      u,qu,velocity,rhsv,smatB,
     &                      smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nit,enazalfa)
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
               
      REAL(8) precv(mesh%Tnunk)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) u(mesh%nnodes)
      REAL(8) qu(mesh%nq)   
      REAL(8) enazalfa   

      REAL(8), ALLOCATABLE :: x(:), b(:), r(:)
      REAL(4) cptime,rcpu      
      
      ALLOCATE (x(mesh%Tnunk),b(mesh%Tnb),r(mesh%nicnsp))

      rcpu=cptime(0.)
C
C     Set up rhs vector and initial approximation of the unknown vector
C      
      DO i=1,mesh%nnodes
        IF (mesh%Tkode(i).LT.0) THEN
          x(ABS(mesh%Tkode(i)))=u(i)
        ELSE
          b(mesh%Tkode(i))=u(i)                
        END IF
      END DO      
      DO i=1,mesh%nq
        IF (mesh%Tqkode(i).LT.0) THEN
          x(ABS(mesh%Tqkode(i)))=qu(i)
        ELSE
          b(mesh%Tqkode(i))=qu(i)        
        END IF
      END DO

C
C     r = rhs * b
C      
      CALL CRSxV(rhsm,b,mesh%Tnb,r)

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
C                    time derivative              
     &               -smatB(i,j)*rhsv(node)*enazalfa
     &               -enazalfa*(
C                    advection and vortex twisting and stretching     
     &               +smatAbdx(i,j)*velocity(node,1)*u(node)
     &               +smatAbdy(i,j)*velocity(node,2)*u(node)
     &               +smatAbdz(i,j)*velocity(node,3)*u(node)
     &                 -smatDx(i,j)*velocity(node,1)*u(node)
     &                 -smatDy(i,j)*velocity(node,2)*u(node)
     &                 -smatDz(i,j)*velocity(node,3)*u(node)
     &                )
          END DO
        END DO
      END DO
      cpu%time(19)=cpu%time(19)+cptime(rcpu)
      rcpu=cptime(0.)
C
C     solve overdetermined system of linear equations
C      
      CALL slvlsqr2(sysm%neq,mesh%Tnunk,sysm%nnz,inp%lsqs_maxit,inp%dlse(7),nit,ierr,
     &              precv,sysm%i,sysm%j,sysm%v,r,x)
     
      cpu%time(20)=cpu%time(20)+cptime(rcpu)

      IF (ierr.NE.0) CALL WarnErr(env,io,inp,4,"SolveFDwTE","NAPAKA V SOVLERJU",ierr)
C
C     Copy solution vector to q and u using under-relaxation
c     under-relaxation should be 1.0 for linear problems
C      
      DO i=1,mesh%nnodes
        IF (mesh%Tkode(i).LT.0) THEN
          u(i)=(1.0D0-inp%urDT)*u(i)+inp%urDT*x(ABS(mesh%Tkode(i)))
        END IF
      END DO      
      DO i=1,mesh%nq
        IF (mesh%Tqkode(i).LT.0) THEN
c          qu(i)=(1.0D0-inp%urDT)*qu(i)+inp%urDT*x(ABS(mesh%Tqkode(i)))
          qu(i)=x(ABS(mesh%Tqkode(i)))
        END IF
      END DO

      DEALLOCATE (x,b,r)
      
      END
      
      

C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsB_Temp_fill(mesh,smatH,smatG,smatB,beta,ren,sysm,rhsm,
     &                      velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)      
C
C     $: Iz pravokotnih matrik g, h in b naredi CRS sistemsko in rhs
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      REAL(8) smatH(mesh%nicnsp,mesh%npoc),smatG(mesh%nicnsp,mesh%npofc)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc),beta,ren
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z   
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc) 
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc) 
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc)   
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) val
      TYPE(matrix) :: sysm,rhsm
      
      INTEGER wqkode(mesh%nq)
      INTEGER wkode(mesh%nnodes)
      INTEGER nunk, nb
            
      INTEGER ic,it,row,col,ii,jj,nuq,j,is,ir
      
      wqkode=mesh%Tqkode
      wkode=mesh%Tkode
      nunk=mesh%Tnunk
      nb=mesh%Tnb

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
            val=smatH(row,col)+smatB(row,col)*beta*ren          
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
C
C       Solve diffusion advection equation using finite difference discretization
C       of time derivative
C
C
      SUBROUTINE SolveFDdaeRCPL(env,io,inp,cpu,mesh,sysm,rhsm,precv,
     &                      u,qu,velocity,timev,sources,smatB,
     &                      smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nit,repr,
     &                      rho,cp,lambda,gradL)
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

      REAL(8) sources(mesh%nnodes)
      REAL(8) timev(mesh%nnodes)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc)
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc)

      REAL(8) rho(mesh%nnodes)
      REAL(8) cp(mesh%nnodes)
      REAL(8) lambda(mesh%nnodes)
      REAL(8) gradL(mesh%nnodes,3)
      REAL(8) velocity(mesh%nnodes,3)

      REAL(8) precv(mesh%Tnunk)
      REAL(8) u(mesh%nnodes)
      REAL(8) qu(mesh%nq)
      REAL(8) repr

      REAL(8), ALLOCATABLE :: x(:), b(:), r(:)
      REAL(4) cptime,rcpu

      ALLOCATE (x(mesh%Tnunk),b(mesh%Tnb),r(mesh%nicnsp))

      rcpu=cptime(0.)
C
C     Set up rhs vector and initial approximation of the unknown vector
C
      DO i=1,mesh%nnodes
        IF (mesh%Tkode(i).LT.0) THEN
          x(ABS(mesh%Tkode(i)))=u(i)
        ELSE
          b(mesh%Tkode(i))=u(i)
        END IF
      END DO
      DO i=1,mesh%nq
        IF (mesh%Tqkode(i).LT.0) THEN
          x(ABS(mesh%Tqkode(i)))=qu(i)
        ELSE
          b(mesh%Tqkode(i))=qu(i)
        END IF
      END DO

C
C     r = rhs * b
C
      CALL CRSxV(rhsm,b,mesh%Tnb,r)

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
C                    time derivative and sources
     &               -smatB(i,j)*(timev(node)*repr*rho(node)*cp(node)-sources(node))
c     &               -u(node)*(
C                    advection
c     &               +smatAbdx(i,j)*velocity(node,1)*rho(node)*cp(node)
c     &               +smatAbdy(i,j)*velocity(node,2)*rho(node)*cp(node)
c     &               +smatAbdz(i,j)*velocity(node,3)*rho(node)*cp(node)
c     &                 -smatDx(i,j)*(velocity(node,1)*rho(node)*cp(node)+gradL(node,1))
c     &                 -smatDy(i,j)*(velocity(node,2)*rho(node)*cp(node)+gradL(node,2))
c     &                 -smatDz(i,j)*(velocity(node,3)*rho(node)*cp(node)+gradL(node,3))
c     &                )
          END DO
        END DO
      END DO
      cpu%time(19)=cpu%time(19)+cptime(rcpu)
      rcpu=cptime(0.)
C
C     solve overdetermined system of linear equations
C
      CALL slvlsqr2(sysm%neq,mesh%Tnunk,sysm%nnz,inp%lsqs_maxit,inp%dlse(7),nit,ierr,
     &              precv,sysm%i,sysm%j,sysm%v,r,x)

      cpu%time(20)=cpu%time(20)+cptime(rcpu)

      IF (ierr.NE.0) CALL WarnErr(env,io,inp,4,"SolveFDwTE","NAPAKA V SOVLERJU",ierr)
C
C     Copy solution vector to q and u using under-relaxation
c     under-relaxation should be 1.0 for linear problems
C
      DO i=1,mesh%nnodes
        IF (mesh%Tkode(i).LT.0) THEN
          u(i)=(1.0D0-inp%urDT)*u(i)+inp%urDT*x(ABS(mesh%Tkode(i)))
        END IF
      END DO
      DO i=1,mesh%nq
        IF (mesh%Tqkode(i).LT.0) THEN
c          qu(i)=(1.0D0-inp%urDT)*qu(i)+inp%urDT*x(ABS(mesh%Tqkode(i)))
          qu(i)=x(ABS(mesh%Tqkode(i)))
        END IF
      END DO

      DEALLOCATE (x,b,r)

      END



C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsB_TempRCPL_fill(mesh,smatH,smatG,smatB,beta,repr,sysm,rhsm,
     &                      smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,
     &                      rho,cp,lambda,qlambda,dRhoCpdt,velocity,gradL)
C
C     $: Iz pravokotnih matrik g, h in b naredi CRS sistemsko in rhs
C        za primer, kjer so rho, cp in lambda odvisni od kraja
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      REAL(8) smatH(mesh%nicnsp,mesh%npoc),smatG(mesh%nicnsp,mesh%npofc)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc),beta,repr
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc)
      REAL(8) rho(mesh%nnodes)
      REAL(8) cp(mesh%nnodes)
      REAL(8) lambda(mesh%nnodes)
      REAL(8) qlambda(mesh%nq)
      REAL(8) dRhoCpdt(mesh%nnodes)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) gradL(mesh%nnodes,3)

      REAL(8) val
      TYPE(matrix) :: sysm,rhsm

      INTEGER wqkode(mesh%nq)
      INTEGER wkode(mesh%nnodes)
      INTEGER nunk, nb

      INTEGER ic,it,row,col,ii,jj,nuq,j,is,ir

      wqkode=mesh%Tqkode
      wkode=mesh%Tkode
      nunk=mesh%Tnunk
      nb=mesh%Tnb

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
            val=lambda(j)*smatH(row,col)
     &                   +smatB(row,col)*repr*( beta*rho(j)*cp(j)+dRhoCpdt(j) )
     &               +(
C                    advection
     &               +smatAbdx(row,col)*velocity(j,1)*rho(j)*cp(j)*repr
     &               +smatAbdy(row,col)*velocity(j,2)*rho(j)*cp(j)*repr
     &               +smatAbdz(row,col)*velocity(j,3)*rho(j)*cp(j)*repr
     &                 -smatDx(row,col)*(velocity(j,1)*rho(j)*cp(j)*repr+gradL(j,1))
     &                 -smatDy(row,col)*(velocity(j,2)*rho(j)*cp(j)*repr+gradL(j,2))
     &                 -smatDz(row,col)*(velocity(j,3)*rho(j)*cp(j)*repr+gradL(j,3))
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
                sysm%v(sysm%d(is))=-smatG(row,col)*qlambda(j)
              ELSE ! znana vrednost - rhs
                ir=ir+1
                rhsm%v(rhsm%d(ir))=smatG(row,col)*qlambda(j) ! minus zato, ker gre na drugo stran enacbe
              END IF
            END DO
          END DO
        END DO
      END DO

      END







C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsB_Temp_SQfill(mesh,smatH,smatG,smatB,beta,ren,sysm,rhsm,
     &                      velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)
C
C     $: Iz pravokotnih matrik g, h in b naredi CRS sistemsko in rhs
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      REAL(8) smatH(mesh%nicnsp,mesh%npoc),smatG(mesh%nicnsp,mesh%npofc)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc),beta,ren
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc)
      REAL(8) velocity(mesh%nnodes,3)
      INTEGER, ALLOCATABLE :: eql(:)
      REAL(8), ALLOCATABLE :: FRsys(:),FRrhs(:)
      REAL(8) val,predznak
      TYPE(matrix) :: sysm,rhsm

      INTEGER wqkode(mesh%nq)
      INTEGER wkode(mesh%nnodes)
      INTEGER nunk, nb

      INTEGER ic,row,col,ii,jj,nuq,j
      INTEGER inode,isqr,isuma,i,kr,ks

      wqkode=mesh%Tqkode
      wkode=mesh%Tkode
      nunk=mesh%Tnunk
      nb=mesh%Tnb
      ALLOCATE (eql(nunk))
      ALLOCATE (FRsys(nunk),FRrhs(nb))
      eql=mesh%Teql


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
            val=smatH(row,col)+smatB(row,col)*beta*ren
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
                FRsys(ABS(nuq))=FRsys(ABS(nuq))-smatG(row,col)
              ELSE ! znana vrednost - rhs
                FRrhs(nuq)=FRrhs(nuq)+smatG(row,col)
              END IF
            END DO
          END DO


          END DO
        ELSE ! neznanka je fluks (CE JE NEZNANKA FLUKS ODSTEVAMO DRUGO ENACBO !!!!)
c         dodajamo prispevke enack, ki se ODSTEJEJO
          DO isuma=1,mesh%sqQlistNO(inode) ! to grem maksimalno do dva
            row=mesh%sqQlist(inode,isuma)
            ic=mesh%sqQlistIC(inode,isuma)

c      v idnode imaš številko vozlišča v Q mreži,
c      na podlagi tega najdi predznak iz flpm ???????

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
            val=smatH(row,col)+smatB(row,col)*beta*ren
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
                FRsys(ABS(nuq))=FRsys(ABS(nuq))-smatG(row,col)*predznak
              ELSE ! znana vrednost - rhs
                FRrhs(nuq)=FRrhs(nuq)+smatG(row,col)*predznak
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
      SUBROUTINE TimeDerRhoCp(mesh,dRhoCpdt,prhocp,rho,cp,tstep)
C
C     $: calculate time derivative of rho*cp
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh

      REAL(8) rho(mesh%nnodes)
      REAL(8) cp(mesh%nnodes)
      REAL(8) dRhoCpdt(mesh%nnodes)
      REAL(8) prhocp(mesh%nnodes)
      REAL(8) tstep

      INTEGER i

      DO i=1,mesh%nnodes
        drhocpdt(i)=( rho(i)*cp(i)-prhocp(i) ) / tstep
        prhocp(i)  =  rho(i)*cp(i)
      END DO

      END


C -----------------------------------------------------------------------------
      SUBROUTINE SetUpAnal(env,io,inp,mesh,Anal,t,model)
C
C     $: sets up analytical solution for model
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(inputtype) :: inp
      TYPE(IOtype)    :: io
      TYPE(penv) :: env

      REAL(8) Anal(mesh%nnodes)
      REAL(8) x,y,z,t

      INTEGER model,i

      IF (model.EQ.11) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          anal(i)=x
        END DO
        RETURN ! dont change BC
      ELSE IF (model.EQ.14) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          anal(i)=t + x + y + z
        END DO
      ELSE IF (model.EQ.15) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          anal(i)=t + x + y + z
        END DO
      ELSE IF (model.EQ.16) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          anal(i)=t + x + y + z
        END DO
      ELSE IF (model.EQ.17) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          anal(i)=t**2.0D0 + x**2.0D0 + y**2.0D0 + z**2.0D0
        END DO
      ELSE IF (model.EQ.18.OR.model.EQ.19) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          anal(i)=t + x + y + z
        END DO
      ELSE IF (model.EQ.20.OR.model.EQ.22.OR.model.EQ.23.OR.model.EQ.24
     &     .OR.model.EQ.25.OR.model.EQ.26.OR.model.EQ.27) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          anal(i)=t*x*y*z
        END DO
      ELSE IF (model.EQ.21) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          anal(i)=t**2.0D0+x**2.0D0+y**2.0D0+z**2.0D0
        END DO
      ELSE
        CALL WarnErr(env,io,inp,4,"SetUpAnal","wrong material and BC model",model)
      END IF

      END



C -----------------------------------------------------------------------------
      SUBROUTINE CalTempRMS(env,io,inp,mesh,cnt,Temp,t,model)
C
C     $: sets up analytical solution for model
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(inputtype) :: inp
      TYPE(countType) :: cnt
      TYPE(IOtype)    :: io
      TYPE(penv) :: env

      REAL(8) Temp(mesh%nnodes)
      REAL(8), ALLOCATABLE :: Anal(:)
      REAL(8) x,y,z,t

      REAL(8) v1,v2,rms

      INTEGER model,i,lun

C
C     First run, set up files
C
      lun=io%rms
      IF (cnt%ftime.EQ.1) THEN
        OPEN (lun,FILE=TRIM(io%rms_name),STATUS='UNKNOWN')
        WRITE (lun,'(A)') 'VARIABLES = "Time", "rms"'
        WRITE (lun,'(A,I3,A)') 'ZONE T = "model=',model, '"'
      END IF
c
C     Prepare analytical solution based on model
C
      ALLOCATE (Anal(mesh%nnodes))
      CALL SetUpAnal(env,io,inp,mesh,Anal,t,model)
C
C     Calculate RMS
C
      v1=0.0D0
      v2=0.0D0

      DO i=1,mesh%nnodes
        IF (mesh%lbn(i).EQ.0) THEN ! only domain nodes
          v1=v1+(anal(i)-temp(i))**2.0D0
          v2=v2+anal(i)**2.0D0
c          print *,anal(i),temp(i)
        END IF
      END DO

      rms = SQRT(v1/v2)

      WRITE (lun,*) t,rms

      DEALLOCATE (Anal)

      END

C -----------------------------------------------------------------------------
      SUBROUTINE SetUpRhoCpLambdaAndBC(env,io,inp,mesh,gp,rho,cp,lambda,Temp,qTemp,velocity,model,sources,t,init,PartVolFrac)
C
C     $: sets up variable material properties, sources and boundary conditions
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(inputtype) :: inp
      TYPE(IOtype)    :: io
      TYPE(penv) :: env
      TYPE(GaussType) :: gp

      REAL(8) rho(mesh%nnodes)
      REAL(8) cp(mesh%nnodes)
      REAL(8) lambda(mesh%nnodes)
      REAL(8) x,y,z,t

      REAL(8) temp(mesh%nnodes)
      REAL(8) PartVolFrac(mesh%nnodes)
      REAL(8) sources(mesh%nnodes)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) qtemp(mesh%nq)

      REAL(8), ALLOCATABLE :: dudx(:),dudy(:),dudz(:),anal(:)
      REAL(8) pi,omega

      INTEGER model,i,j,ii,dn,init

      ALLOCATE (dudx(mesh%nq),dudy(mesh%nq),dudz(mesh%nq),anal(mesh%nnodes))
      dudx=0.0D0 ! za tiste, kjer ne podam
      dudy=0.0D0
      dudz=0.0D0
      sources=0.0D0
      pi=4.0D0*ATAN(1.0D0)
C
C     Set up analytical solution based on t and model
C
C     just for test models
      IF (model.GT.10) CALL SetUpAnal(env,io,inp,mesh,Anal,t,model)

      IF (model.EQ.1) THEN
        DO i=1,mesh%nnodes
          rho(i)   = 1.0D0
          cp(i)    = 1.0D0
          lambda(i)= 1.0D0
        END DO
        RETURN ! dont change BC
      ELSE IF (model.EQ.2) THEN
        CALL NanofluidMaterialProperties(inp,mesh,gp,rho,cp,lambda,PartVolFrac)
        RETURN ! dont change BC
      ELSE IF (model.EQ.14) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          rho(i)   = 1.0D0
          cp(i)    = y*t
          lambda(i)= 1.0D0 + x + y + z
          velocity(i,1)=x
          velocity(i,2)=0.0D0
          velocity(i,3)=0.0D0
          sources(i)=-3.0D0 + t*y + t*x*y + y*(t + x + y + z) + t*y*(t + x + y + z) ! sources
        END DO
        DO i=1,mesh%nq !  iks=mesh%xq(i,1),  ips=mesh%xq(i,2),  ze=mesh%xq(i,3)
          dudx(i)=1.0D0
          dudy(i)=1.0D0
          dudz(i)=1.0D0
        END DO
C       Initial conditions
        IF (init.EQ.1) THEN
          DO i=1,mesh%nnodes
            Temp(i)=anal(i)
          END DO
        END IF
      ELSE IF (model.EQ.15) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          rho(i)   = 1.0D0
          cp(i)    = y*t
          lambda(i)= 1.0D0 + x + y + z
          velocity(i,1)=x
          velocity(i,2)=-y
          velocity(i,3)=0.0D0
          sources(i)=-3.0D0 + t*y + t*x*y - t*y**2 + y*(t + x + y + z) - t*y*(t + x + y + z) ! sources
        END DO
        DO i=1,mesh%nq !  iks=mesh%xq(i,1),  ips=mesh%xq(i,2),  ze=mesh%xq(i,3)
          dudx(i)=1.0D0
          dudy(i)=1.0D0
          dudz(i)=1.0D0
        END DO
C       Initial conditions
        IF (init.EQ.1) THEN
          DO i=1,mesh%nnodes
            Temp(i)=anal(i)
          END DO
        END IF
      ELSE IF (model.EQ.16) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          rho(i)   = 1.0D0
          cp(i)    = 1.0D0
          lambda(i)= 3.0D0 + Sin(x*y*z)
          velocity(i,1)=x
          velocity(i,2)=-y
          velocity(i,3)=0.0D0
          sources(i)=1.0D0 + x - y - (y*z + x*(y + z))*Cos(x*y*z) ! sources
        END DO
        DO i=1,mesh%nq !  iks=mesh%xq(i,1),  ips=mesh%xq(i,2),  ze=mesh%xq(i,3)
          dudx(i)=1.0D0
          dudy(i)=1.0D0
          dudz(i)=1.0D0
        END DO
C       Initial conditions
        IF (init.EQ.1) THEN
          DO i=1,mesh%nnodes
            Temp(i)=anal(i)
          END DO
        END IF
      ELSE IF (model.EQ.17) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          rho(i)   = 1.0D0
          cp(i)    = 1.0D0
          lambda(i)= 1.0D0 + x**2.0D0 + y**2.0D0 + z**2.0D0
          velocity(i,1)=x
          velocity(i,2)=-y
          velocity(i,3)=0.0D0
          sources(i)=2.0D0* (-3.0D0 + t - 4.0D0*x**2.0D0 - 6.0D0*y**2.0D0 - 5.0D0*z**2.0D0 )    ! sources
        END DO
        DO i=1,mesh%nq !  iks=mesh%xq(i,1),  ips=mesh%xq(i,2),  ze=mesh%xq(i,3)
          x=mesh%xq(i,1)
          y=mesh%xq(i,2)
          z=mesh%xq(i,3)
          dudx(i)=2.0D0*x
          dudy(i)=2.0D0*y
          dudz(i)=2.0D0*z
        END DO
C       Initial conditions
        IF (init.EQ.1) THEN
          DO i=1,mesh%nnodes
            Temp(i)=anal(i)
          END DO
        END IF
      ELSE IF (model.EQ.18) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          rho(i)   = 1.0D0
          cp(i)    = 1.0D0
          lambda(i)= 1.0D0
          velocity(i,1)=x
          velocity(i,2)=-2.0D0*y
          velocity(i,3)=z
          sources(i)=1.0D0+x-2.0D0*y+z ! sources
        END DO
        DO i=1,mesh%nq !  iks=mesh%xq(i,1),  ips=mesh%xq(i,2),  ze=mesh%xq(i,3)
          x=mesh%xq(i,1)
          y=mesh%xq(i,2)
          z=mesh%xq(i,3)
          dudx(i)=1.0D0
          dudy(i)=1.0D0
          dudz(i)=1.0D0
        END DO
C       Initial conditions
        IF (init.EQ.1) THEN
          DO i=1,mesh%nnodes
            Temp(i)=anal(i)
          END DO
        END IF
      ELSE IF (model.EQ.19) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          rho(i)   = 1.0D0
          cp(i)    = 1.0D0+x*y
          lambda(i)= 1.0D0+z*y
          velocity(i,1)=x
          velocity(i,2)=-2.0D0*y
          velocity(i,3)=z
          sources(i)=  1.0D0 - 3.0D0* y + x* (1.0D0 + y - t* y - 3.0D0* y**2.0D0) ! sources
        END DO
        DO i=1,mesh%nq !  iks=mesh%xq(i,1),  ips=mesh%xq(i,2),  ze=mesh%xq(i,3)
          x=mesh%xq(i,1)
          y=mesh%xq(i,2)
          z=mesh%xq(i,3)
          dudx(i)=1.0D0
          dudy(i)=1.0D0
          dudz(i)=1.0D0
        END DO
C       Initial conditions
        IF (init.EQ.1) THEN
          DO i=1,mesh%nnodes
            Temp(i)=anal(i)
          END DO
        END IF
      ELSE IF (model.EQ.20) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          rho(i)   = 1.0D0
          cp(i)    = 1.0D0+x*y
          lambda(i)= 1.0D0+z*y
          velocity(i,1)=x
          velocity(i,2)=-2.0D0*y
          velocity(i,3)=z
          sources(i)=   x*(y*(1.0D0 + x*y)*z - t*(z**2.0D0 + y**2.0D0*(1.0D0 + x*z)))   ! sources
        END DO
        DO i=1,mesh%nq !  iks=mesh%xq(i,1),  ips=mesh%xq(i,2),  ze=mesh%xq(i,3)
          x=mesh%xq(i,1)
          y=mesh%xq(i,2)
          z=mesh%xq(i,3)
          dudx(i)=t*y*z
          dudy(i)=t*x*z
          dudz(i)=t*x*y
        END DO
C       Initial conditions
        IF (init.EQ.1) THEN
          DO i=1,mesh%nnodes
            Temp(i)=anal(i)
          END DO
        END IF

      ELSE IF (model.EQ.21) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          rho(i)   = 1.0D0
          cp(i)    = 1.0D0
          lambda(i)= 1.0D0
          velocity(i,1)=x
          velocity(i,2)=-2.0D0*y
          velocity(i,3)=z
          sources(i)= 2.0D0* (-3.0D0 + t + x**2.0D0 - 2.0D0* y**2.0D0 + z**2.0D0);  ! sources
        END DO
        DO i=1,mesh%nq !  iks=mesh%xq(i,1),  ips=mesh%xq(i,2),  ze=mesh%xq(i,3)
          x=mesh%xq(i,1)
          y=mesh%xq(i,2)
          z=mesh%xq(i,3)
          dudx(i)=2.0D0*x
          dudy(i)=2.0D0*y
          dudz(i)=2.0D0*z
        END DO
C       Initial conditions
        IF (init.EQ.1) THEN
          DO i=1,mesh%nnodes
            Temp(i)=anal(i)
          END DO
        END IF

      ELSE IF (model.EQ.22) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          rho(i)   = 1.0D0
          cp(i)    = 1.0D0
          lambda(i)= 2.0D0  + Cos(omega*pi*x*z*y)
          velocity(i,1)=x
          velocity(i,2)=-2.0D0*y
          velocity(i,3)=z
          omega=1.0D0
          sources(i)= x*y*z + pi*t*(y**2.0D0*z**2.0D0 + x**2.0D0*(y**2.0D0 + z**2.0D0))*omega*sin(pi*x*y*z*omega)  ! sources
        END DO
        DO i=1,mesh%nq !  iks=mesh%xq(i,1),  ips=mesh%xq(i,2),  ze=mesh%xq(i,3)
          x=mesh%xq(i,1)
          y=mesh%xq(i,2)
          z=mesh%xq(i,3)
          dudx(i)=t*y*z
          dudy(i)=t*x*z
          dudz(i)=t*x*y
        END DO
C       Initial conditions
        IF (init.EQ.1) THEN
          DO i=1,mesh%nnodes
            Temp(i)=anal(i)
          END DO
        END IF


      ELSE IF (model.EQ.23) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          rho(i)   = 1.0D0
          cp(i)    = 1.0D0
          lambda(i)= 2.0D0  + Cos(omega*pi*x*z*y)
          velocity(i,1)=x
          velocity(i,2)=-2.0D0*y
          velocity(i,3)=z
          omega=2.0D0
          sources(i)= x*y*z + pi*t*(y**2.0D0*z**2.0D0 + x**2.0D0*(y**2.0D0 + z**2.0D0))*omega*sin(pi*x*y*z*omega)  ! sources
        END DO
        DO i=1,mesh%nq !  iks=mesh%xq(i,1),  ips=mesh%xq(i,2),  ze=mesh%xq(i,3)
          x=mesh%xq(i,1)
          y=mesh%xq(i,2)
          z=mesh%xq(i,3)
          dudx(i)=t*y*z
          dudy(i)=t*x*z
          dudz(i)=t*x*y
        END DO
C       Initial conditions
        IF (init.EQ.1) THEN
          DO i=1,mesh%nnodes
            Temp(i)=anal(i)
          END DO
        END IF


      ELSE IF (model.EQ.24) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          rho(i)   = 1.0D0
          cp(i)    = 1.0D0
          lambda(i)= 2.0D0  + Cos(omega*pi*x*z*y)
          velocity(i,1)=x
          velocity(i,2)=-2.0D0*y
          velocity(i,3)=z
          omega=4.0D0
          sources(i)= x*y*z + pi*t*(y**2.0D0*z**2.0D0 + x**2.0D0*(y**2.0D0 + z**2.0D0))*omega*sin(pi*x*y*z*omega)  ! sources
        END DO
        DO i=1,mesh%nq !  iks=mesh%xq(i,1),  ips=mesh%xq(i,2),  ze=mesh%xq(i,3)
          x=mesh%xq(i,1)
          y=mesh%xq(i,2)
          z=mesh%xq(i,3)
          dudx(i)=t*y*z
          dudy(i)=t*x*z
          dudz(i)=t*x*y
        END DO
C       Initial conditions
        IF (init.EQ.1) THEN
          DO i=1,mesh%nnodes
            Temp(i)=anal(i)
          END DO
        END IF


      ELSE IF (model.EQ.25) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          omega=1.0D0
          rho(i)   = 1.0D0
          cp(i)    = 1.0D0  + x**omega
          lambda(i)= 1.0D0  + z**omega
          velocity(i,1)=x
          velocity(i,2)=-2.0D0*y
          velocity(i,3)=z
          sources(i)=x*y*(-t*omega*z**(omega-1.0D0) +  z*(1.0D0 + x**omega*(1.0D0 + t*omega)))   ! sources
        END DO
        DO i=1,mesh%nq !  iks=mesh%xq(i,1),  ips=mesh%xq(i,2),  ze=mesh%xq(i,3)
          x=mesh%xq(i,1)
          y=mesh%xq(i,2)
          z=mesh%xq(i,3)
          dudx(i)=t*y*z
          dudy(i)=t*x*z
          dudz(i)=t*x*y
        END DO
C       Initial conditions
        IF (init.EQ.1) THEN
          DO i=1,mesh%nnodes
            Temp(i)=anal(i)
          END DO
        END IF



      ELSE IF (model.EQ.26) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          omega=2.0D0
          rho(i)   = 1.0D0
          cp(i)    = 1.0D0  + x**omega
          lambda(i)= 1.0D0  + z**omega
          velocity(i,1)=x
          velocity(i,2)=-2.0D0*y
          velocity(i,3)=z
          sources(i)=x*y*(-t*omega*z**(omega-1.0D0) +  z*(1.0D0 + x**omega*(1.0D0 + t*omega)))   ! sources
        END DO
        DO i=1,mesh%nq !  iks=mesh%xq(i,1),  ips=mesh%xq(i,2),  ze=mesh%xq(i,3)
          x=mesh%xq(i,1)
          y=mesh%xq(i,2)
          z=mesh%xq(i,3)
          dudx(i)=t*y*z
          dudy(i)=t*x*z
          dudz(i)=t*x*y
        END DO
C       Initial conditions
        IF (init.EQ.1) THEN
          DO i=1,mesh%nnodes
            Temp(i)=anal(i)
          END DO
        END IF



      ELSE IF (model.EQ.27) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
          omega=4.0D0
          rho(i)   = 1.0D0
          cp(i)    = 1.0D0  + x**omega
          lambda(i)= 1.0D0  + z**omega
          velocity(i,1)=x
          velocity(i,2)=-2.0D0*y
          velocity(i,3)=z
          sources(i)=x*y*(-t*omega*z**(omega-1.0D0) +  z*(1.0D0 + x**omega*(1.0D0 + t*omega)))   ! sources
        END DO
        DO i=1,mesh%nq !  iks=mesh%xq(i,1),  ips=mesh%xq(i,2),  ze=mesh%xq(i,3)
          x=mesh%xq(i,1)
          y=mesh%xq(i,2)
          z=mesh%xq(i,3)
          dudx(i)=t*y*z
          dudy(i)=t*x*z
          dudz(i)=t*x*y
        END DO
C       Initial conditions
        IF (init.EQ.1) THEN
          DO i=1,mesh%nnodes
            Temp(i)=anal(i)
          END DO
        END IF


      ELSE
        CALL WarnErr(env,io,inp,4,"SetUpRhoCpLambdaAndBC","wrong material and BC model",model)
      END IF




C     Change boundary conditions

C     Dirichlet, kode=0
      DO i=1,mesh%nbnodes
        dn=mesh%gbn(i)
        IF (mesh%Tkode(dn).GT.0) THEN ! to je znana vrednost, mesh%iwn(i) je stevilka stene
          Temp(dn)=anal(dn)
        END IF
      END DO

C     Neumann, kode=2
      ii=0
      DO i=1,mesh%nbelem
        DO j=1,mesh%npof
          ii=ii+1
          dn=mesh%ibcf(i,j)
          IF (mesh%Tqkode(dn).GT.0) THEN ! to je znana vrednost, mesh%bewn(i) je stevilka stene
            IF (mesh%bewn(i).EQ.1) THEN ! stena pri z=1
              qTemp(dn)=-dudz(dn)
            END IF
            IF (mesh%bewn(i).EQ.6) THEN ! stena pri z=2
              qTemp(dn)=dudz(dn)
            END IF
            IF (mesh%bewn(i).EQ.5) THEN ! stena pri x=1
              qTemp(dn)=-dudx(dn)
            END IF
            IF (mesh%bewn(i).EQ.3) THEN ! stena pri x=2
              qTemp(dn)=dudx(dn)
            END IF
            IF (mesh%bewn(i).EQ.2) THEN ! stena pri y=1
              qTemp(dn)=-dudy(dn)
            END IF
            IF (mesh%bewn(i).EQ.4) THEN ! stena pri y=2
              qTemp(dn)=dudy(dn)
            END IF
          END IF
        END DO
      END DO


      DEALLOCATE (dudx,dudy,dudz,anal)

      END



C -----------------------------------------------------------------------------
      SUBROUTINE RhoCpsetGrad(mesh,rho,cp,grad)
C
C     $: grad(rho*cp)
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh

      REAL(8) rho(mesh%nnodes)
      REAL(8) cp(mesh%nnodes)
      REAL(8) grad(mesh%nnodes,3)

      REAL(8), ALLOCATABLE :: rhocp(:)

      INTEGER i

      ALLOCATE (rhocp(mesh%nnodes))

      DO i=1,mesh%nnodes
        rhocp(i)=rho(i)*cp(i)
      END DO

      CALL  setGrad(mesh,rhocp,grad)

      DEALLOCATE (rhocp)

      END


C -----------------------------------------------------------------------------
      SUBROUTINE SetUpFvectorDivF(env,io,inp,mesh,rho,cp,lambda,velocity,RePr,f,Trhsv)
C
C     $: f = Re*Pr*rho*cp/lambda - (grad lambda) / lambda
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(inputtype) :: inp
      TYPE(IOtype)    :: io
      TYPE(penv) :: env

      REAL(8) rho(mesh%nnodes)
      REAL(8) cp(mesh%nnodes)
      REAL(8) lambda(mesh%nnodes)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) RePr
      REAL(8), ALLOCATABLE :: grad(:,:)

      REAL(8) f(mesh%nnodes,3)
      REAL(8) Trhsv(mesh%nnodes)

      INTEGER i,j

      ALLOCATE (grad(mesh%nnodes,3))

C     Calculate gradient of thermal conductivity
      CALL setGrad(mesh,lambda,grad)
C     Sum up f vector
      DO i=1,mesh%nnodes
        DO j=1,3
          f(i,j)=RePr*rho(i)*cp(i)/lambda(i)*velocity(i,j)-grad(i,j)/lambda(i)
        END DO
      END DO

C     Calculate divergence of f and add it to sources
      CALL setGrad(mesh,f,grad)
      DO i=1,mesh%nnodes
        Trhsv(i)=Trhsv(i) - (grad(i,1)+grad(i,2)+grad(i,3)) / (RePr*rho(i)*cp(i)/lambda(i))
c       delim zato, ker kasneje nazaj mnozim, da se lahko sesteva skupaj s casovnim odovodom
      END DO

      DEALLOCATE (grad)

      END


