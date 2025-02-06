C -----------------------------------------------------------------------------      
C
C       Solve vorticity transport equation using finite difference discretization
C       of time derivative (DOLOCEN SISTEM ENACB)
C
C
      SUBROUTINE SolveFDwTEsq(eqn,env,io,inp,cpu,mesh,sysm,rhsm,prec,
     &                      vorticity,qvorticity,temp,dc,velocity,kbf,rhsv,smatB,
     &                      smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nit,nanoA,nanoB)
C
C     $: resi H*u=G*q+B*rhsv-advekcija-vortextwisting
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
      INTEGER node
      INTEGER eqn  ! equation number 1=w_x unknown, 2=y, 3=z
      INTEGER inode,isqr,isuma
      REAL(8) predznak

      REAL(8) rhsv(mesh%nnodes)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc)
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc)

c      REAL(8) precv(mesh%nunk(eqn))
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) kbf(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3)
      REAL(8) qvorticity(mesh%nq,3)
      REAL(8) temp(mesh%nnodes)
      REAL(8) dc(mesh%nnodes)
      REAL(8) nanoA,nanoB

      REAL(8), ALLOCATABLE :: x(:), b(:), r(:)
      REAL(4) cptime,rcpu,xcpu

      ALLOCATE (x(mesh%nunk(eqn)),b(mesh%nb(eqn)),r(mesh%nunk(eqn)))

      rcpu=cptime(0.)
C
C     Set up rhs vector and initial approximation of the unknown vector
C
      DO i=1,mesh%nnodes
        IF (mesh%wkode(i,eqn).LT.0) THEN
          x(ABS(mesh%wkode(i,eqn)))=vorticity(i,eqn)
        ELSE
          b(mesh%wkode(i,eqn))=vorticity(i,eqn)
        END IF
      END DO
      DO i=1,mesh%nq
        IF (mesh%wqkode(i,eqn).LT.0) THEN
          x(ABS(mesh%wqkode(i,eqn)))=qvorticity(i,eqn)
        ELSE
          b(mesh%wqkode(i,eqn))=qvorticity(i,eqn)
        END IF
      END DO
C
C     r = rhs * b
C
      CALL CRSxV(rhsm,b,mesh%nb(eqn),r)
c
c     r = r - B * rhsv
c
      DO isqr=1,sysm%neq !nunk ! zanka po vrsticah kvadratne matrike
        inode=ABS(mesh%Weql(isqr,eqn))  ! vozlisce v U ali Q mrezi, kateremu pripada ta enacba
        IF (mesh%Weql(isqr,eqn).GT.0) THEN ! to pomeni, da je neznanka funkcija
c         dodajamo prispevke enack, ki se sestejejo
          DO isuma=1,mesh%sqUlistNO(inode)
            i=mesh%sqUlist(inode,isuma)
            ic=mesh%sqUlistIC(inode,isuma)
            DO j=1,mesh%npoc
              node=mesh%idc(ic,j)
              r(isqr)=r(isqr)
C                    time derivative
     &               -smatB(i,j)*rhsv(node)*inp%ren/nanoA
     &               -inp%ren/nanoA*(
C                    advection and vortex twisting and stretching
     &               +smatAbdx(i,j)*(velocity(node,1)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,1))
     &               +smatAbdy(i,j)*(velocity(node,2)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,2))
     &               +smatAbdz(i,j)*(velocity(node,3)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,3))
     &                 -smatDx(i,j)*(velocity(node,1)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,1))
     &                 -smatDy(i,j)*(velocity(node,2)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,2))
     &                 -smatDz(i,j)*(velocity(node,3)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,3))
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
     &               -smatB(i,j)*rhsv(node)*inp%ren/nanoA
     &               -inp%ren/nanoA*(
C                    advection and vortex twisting and stretching
     &               +smatAbdx(i,j)*(velocity(node,1)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,1))
     &               +smatAbdy(i,j)*(velocity(node,2)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,2))
     &               +smatAbdz(i,j)*(velocity(node,3)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,3))
     &                 -smatDx(i,j)*(velocity(node,1)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,1))
     &                 -smatDy(i,j)*(velocity(node,2)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,2))
     &                 -smatDz(i,j)*(velocity(node,3)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,3))
     &                ))
            END DO
          END DO
        END IF
      END DO

c
c     Bouyancy
c
      IF (inp%iDT.GT.0) THEN
        CALL BouyancyRHSsq(eqn,mesh,inp,temp,dc,kbf,r,
     &                   smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nanoA,nanoB)
      END IF

      cpu%time((eqn-1)*6+5)=cpu%time((eqn-1)*6+5)+cptime(rcpu)
      rcpu=cptime(0.)

C
C     solve system of linear equations
C
      CALL SolvSLE(inp%sqrs_type,inp%sqrs_prec,2,inp%sqrs_maxit,5,inp%dlse(3+eqn),sysm%neq,sysm%nnz,
     &                   prec,sysm,r,x,nit,xcpu,ierr)

      cpu%time((eqn-1)*6+6)=cpu%time((eqn-1)*6+6)+cptime(rcpu)

      IF (ierr.NE.0) CALL WarnErr(env,io,inp,4,"SolveFDwTE","NAPAKA V SOVLERJU",ierr)
C
C     Copy solution vector to q and u using under-relaxation
c     under-relaxation should be 1.0 for linear problems
C

      IF (inp%iDWUR.GT.0) THEN ! proti steni je ur povecana
        DO i=1,mesh%nnodes
          IF (mesh%wkode(i,eqn).LT.0) THEN
            vorticity(i,eqn)=(1.0D0-inp%urDw(eqn)*mesh%dwur(i))*vorticity(i,eqn)
     &                             +inp%urDw(eqn)*mesh%dwur(i)*x(ABS(mesh%wkode(i,eqn)))
          END IF
        END DO
      ELSE ! normalno, ur ves cas enaka
        DO i=1,mesh%nnodes
          IF (mesh%wkode(i,eqn).LT.0) THEN
            vorticity(i,eqn)=(1.0D0-inp%urDw(eqn))*vorticity(i,eqn)+inp%urDw(eqn)*x(ABS(mesh%wkode(i,eqn)))
          END IF
        END DO
      END IF

c     fluksi niso podrelaksirani, ker jih nikjer ne rabim.
c     sluzijo samo kot zacetni priblizek za novo iteracijo in zmanjsujejo stevilo iteracij
      DO i=1,mesh%nq
        IF (mesh%wqkode(i,eqn).LT.0) THEN
c          qvorticity(i,eqn)=(1.0D0-inp%urDw(eqn))*qvorticity(i,eqn)+inp%urDw(eqn)*x(ABS(mesh%wqkode(i,eqn)))
          qvorticity(i,eqn)=x(ABS(mesh%wqkode(i,eqn)))
        END IF
      END DO

      DEALLOCATE (x,b,r)

      END



C -----------------------------------------------------------------------------
C
C       Solve vorticity transport equation using finite difference discretization
C       of time derivative
C
C
      SUBROUTINE SolveFDwTE(eqn,env,io,inp,cpu,mesh,sysm,rhsm,precv,
     &                      vorticity,qvorticity,temp,dc,velocity,kbf,rhsv,smatB,
     &                      smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nit,nanoA,nanoB)
C
C     $: resi H*u=G*q+B*rhsv-advekcija-vortextwisting
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
      INTEGER eqn  ! equation number 1=w_x unknown, 2=y, 3=z
      
      REAL(8) rhsv(mesh%nnodes)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc)
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z   
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc) 
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc) 
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc) 
               
      REAL(8) precv(mesh%nunk(eqn))
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) kbf(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3)
      REAL(8) qvorticity(mesh%nq,3)      
      REAL(8) temp(mesh%nnodes)
      REAL(8) dc(mesh%nnodes)
      REAL(8) nanoA,nanoB

      REAL(8), ALLOCATABLE :: x(:), b(:), r(:)
      REAL(4) cptime,rcpu      
      
      ALLOCATE (x(mesh%nunk(eqn)),b(mesh%nb(eqn)),r(mesh%nicnsp))

      rcpu=cptime(0.)
C
C     Set up rhs vector and initial approximation of the unknown vector
C      
      DO i=1,mesh%nnodes
        IF (mesh%wkode(i,eqn).LT.0) THEN
          x(ABS(mesh%wkode(i,eqn)))=vorticity(i,eqn)
        ELSE
          b(mesh%wkode(i,eqn))=vorticity(i,eqn)                
        END IF
      END DO      
      DO i=1,mesh%nq
        IF (mesh%wqkode(i,eqn).LT.0) THEN
          x(ABS(mesh%wqkode(i,eqn)))=qvorticity(i,eqn)
        ELSE
          b(mesh%wqkode(i,eqn))=qvorticity(i,eqn)        
        END IF
      END DO
C
C     r = rhs * b
C      
      CALL CRSxV(rhsm,b,mesh%nb(eqn),r)      
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
     &               -smatB(i,j)*rhsv(node)*inp%ren/nanoA
     &               -inp%ren/nanoA*(
C                    advection and vortex twisting and stretching     
     &               +smatAbdx(i,j)*(velocity(node,1)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,1))
     &               +smatAbdy(i,j)*(velocity(node,2)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,2))
     &               +smatAbdz(i,j)*(velocity(node,3)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,3))
     &                 -smatDx(i,j)*(velocity(node,1)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,1))
     &                 -smatDy(i,j)*(velocity(node,2)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,2))
     &                 -smatDz(i,j)*(velocity(node,3)*vorticity(node,eqn)-velocity(node,eqn)*vorticity(node,3))
     &                )
          END DO
        END DO
      END DO
c
c     Bouyancy
c      
      IF (inp%iDT.GT.0) THEN
        CALL BouyancyRHS(eqn,mesh,inp,temp,dc,kbf,r,
     &                   smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nanoA,nanoB)
      END IF
      
      cpu%time((eqn-1)*6+5)=cpu%time((eqn-1)*6+5)+cptime(rcpu)
      rcpu=cptime(0.)
C
C     solve overdetermined system of linear equations
C      
      CALL slvlsqr2(sysm%neq,mesh%nunk(eqn),sysm%nnz,inp%lsqs_maxit,inp%dlse(3+eqn),nit,ierr,
     &              precv,sysm%i,sysm%j,sysm%v,r,x)
     
      cpu%time((eqn-1)*6+6)=cpu%time((eqn-1)*6+6)+cptime(rcpu)

      IF (ierr.NE.0) CALL WarnErr(env,io,inp,4,"SolveFDwTE","NAPAKA V SOVLERJU",ierr)
C
C     Copy solution vector to q and u using under-relaxation
c     under-relaxation should be 1.0 for linear problems
C      

      IF (inp%iDWUR.GT.0) THEN ! proti steni je ur povecana
        DO i=1,mesh%nnodes
          IF (mesh%wkode(i,eqn).LT.0) THEN
            vorticity(i,eqn)=(1.0D0-inp%urDw(eqn)*mesh%dwur(i))*vorticity(i,eqn)
     &                             +inp%urDw(eqn)*mesh%dwur(i)*x(ABS(mesh%wkode(i,eqn)))
          END IF
        END DO
      ELSE ! normalno, ur ves cas enaka
        DO i=1,mesh%nnodes
          IF (mesh%wkode(i,eqn).LT.0) THEN
            vorticity(i,eqn)=(1.0D0-inp%urDw(eqn))*vorticity(i,eqn)+inp%urDw(eqn)*x(ABS(mesh%wkode(i,eqn)))
          END IF
        END DO
      END IF

c     fluksi niso podrelaksirani, ker jih nikjer ne rabim.
c     sluzijo samo kot zacetni priblizek za novo iteracijo in zmanjsujejo stevilo iteracij
      DO i=1,mesh%nq
        IF (mesh%wqkode(i,eqn).LT.0) THEN
c          qvorticity(i,eqn)=(1.0D0-inp%urDw(eqn))*qvorticity(i,eqn)+inp%urDw(eqn)*x(ABS(mesh%wqkode(i,eqn)))
          qvorticity(i,eqn)=x(ABS(mesh%wqkode(i,eqn)))
        END IF
      END DO

      DEALLOCATE (x,b,r)
      
      END


C -----------------------------------------------------------------------------
      SUBROUTINE BouyancyRHSsq(eqn,mesh,inp,temp,dc,kbf,r,
     &                       smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nanoA,nanoB)
C
C     $: sets up bouyancy on right hand side
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(inputtype) inp

      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc)

      REAL(8) temp(mesh%nnodes)
      REAL(8) dc(mesh%nnodes)
      REAL(8) kbf(mesh%nnodes,3)
      REAL(8) fac,facm,nanoA,nanoB

      INTEGER eqn
      INTEGER i,j,ic,node

      REAL(8) r(mesh%nunk(eqn))
      INTEGER inode,isqr,isuma
      REAL(8) predznak

      fac=inp%ran/inp%prn/inp%ren*nanoB/nanoA
      facm=inp%ram/inp%prn/inp%ren*nanoB/nanoA


      IF (eqn.EQ.1) THEN
c       *** x ***
c      DO j=1,mesh%npoc
        DO isqr=1,mesh%nunk(eqn) !sysm%neq !nunk ! zanka po vrsticah kvadratne matrike
          inode=ABS(mesh%Weql(isqr,eqn))  ! vozlisce v U ali Q mrezi, kateremu pripada ta enacba
          IF (mesh%Weql(isqr,eqn).GT.0) THEN ! to pomeni, da je neznanka funkcija
c           dodajamo prispevke enack, ki se sestejejo
            DO isuma=1,mesh%sqUlistNO(inode)
              i=mesh%sqUlist(inode,isuma)
              ic=mesh%sqUlistIC(inode,isuma)
              DO j=1,mesh%npoc
              node=mesh%idc(ic,j)
              r(isqr)=r(isqr)-fac*(
     &               +inp%gz*(smatAbdy(i,j)-smatDy(i,j))*( temp(node) - dc(node) / inp%Rro )
     &               -inp%gy*(smatAbdz(i,j)-smatDz(i,j))*( temp(node) - dc(node) / inp%Rro ) )
     &               +facm*(
     &               +(smatAbdy(i,j)-smatDy(i,j))*kbf(node,3)
     &               -(smatAbdz(i,j)-smatDz(i,j))*kbf(node,2))
            END DO
            END DO
          ELSE ! neznanka je fluks (CE JE NEZNANKA FLUKS ODSTEVAMO DRUGO ENACBO !!!!)
c         dodajamo prispevke enack,  ki se ODSTEJEJO
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
              r(isqr)=r(isqr)+predznak*(-fac*(
     &               +inp%gz*(smatAbdy(i,j)-smatDy(i,j))*( temp(node) - dc(node) / inp%Rro )
     &               -inp%gy*(smatAbdz(i,j)-smatDz(i,j))*( temp(node) - dc(node) / inp%Rro ) )
     &               +facm*(
     &               +(smatAbdy(i,j)-smatDy(i,j))*kbf(node,3)
     &               -(smatAbdz(i,j)-smatDz(i,j))*kbf(node,2)))
            END DO
            END DO
          END IF
        END DO
c      END DO
      ELSE IF (eqn.EQ.2) THEN
c       *** y ***
c      DO j=1,mesh%npoc
        DO isqr=1,mesh%nunk(eqn) !sysm%neq !nunk ! zanka po vrsticah kvadratne matrike
          inode=ABS(mesh%Weql(isqr,eqn))  ! vozlisce v U ali Q mrezi, kateremu pripada ta enacba
          IF (mesh%Weql(isqr,eqn).GT.0) THEN ! to pomeni, da je neznanka funkcija
c           dodajamo prispevke enack, ki se sestejejo
            DO isuma=1,mesh%sqUlistNO(inode)
              i=mesh%sqUlist(inode,isuma)
              ic=mesh%sqUlistIC(inode,isuma)
              DO j=1,mesh%npoc
              node=mesh%idc(ic,j)
              r(isqr)=r(isqr)-fac*(
     &               -inp%gz*(smatAbdx(i,j)-smatDx(i,j))*( temp(node) - dc(node) / inp%Rro )
     &               +inp%gx*(smatAbdz(i,j)-smatDz(i,j))*( temp(node) - dc(node) / inp%Rro ) )
     &               +facm*(
     &               +(smatAbdx(i,j)-smatDx(i,j))*kbf(node,3)
     &               -(smatAbdz(i,j)-smatDz(i,j))*kbf(node,1))

            END DO
            END DO
          ELSE ! neznanka je fluks (CE JE NEZNANKA FLUKS ODSTEVAMO DRUGO ENACBO !!!!)
c         dodajamo prispevke enack,  ki se ODSTEJEJO
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
              r(isqr)=r(isqr)+predznak*(-fac*(
     &               -inp%gz*(smatAbdx(i,j)-smatDx(i,j))*( temp(node) - dc(node) / inp%Rro )
     &               +inp%gx*(smatAbdz(i,j)-smatDz(i,j))*( temp(node) - dc(node) / inp%Rro ) )
     &               +facm*(
     &               +(smatAbdx(i,j)-smatDx(i,j))*kbf(node,3)
     &               -(smatAbdz(i,j)-smatDz(i,j))*kbf(node,1)))
            END DO
            END DO
          END IF
        END DO
c      END DO
      ELSE
c       *** z ***
c      DO j=1,mesh%npoc
        DO isqr=1,mesh%nunk(eqn) !sysm%neq !nunk ! zanka po vrsticah kvadratne matrike
          inode=ABS(mesh%Weql(isqr,eqn))  ! vozlisce v U ali Q mrezi, kateremu pripada ta enacba
          IF (mesh%Weql(isqr,eqn).GT.0) THEN ! to pomeni, da je neznanka funkcija
c           dodajamo prispevke enack, ki se sestejejo
            DO isuma=1,mesh%sqUlistNO(inode)
              i=mesh%sqUlist(inode,isuma)
              ic=mesh%sqUlistIC(inode,isuma)
              DO j=1,mesh%npoc
              node=mesh%idc(ic,j)
              r(isqr)=r(isqr)-fac*(
     &               +inp%gy*(smatAbdx(i,j)-smatDx(i,j))*( temp(node) - dc(node) / inp%Rro )
     &               -inp%gx*(smatAbdy(i,j)-smatDy(i,j))*( temp(node) - dc(node) / inp%Rro ) )
     &               +facm*(
     &               +(smatAbdx(i,j)-smatDx(i,j))*kbf(node,2)
     &               -(smatAbdy(i,j)-smatDy(i,j))*kbf(node,1))
            END DO
            END DO
          ELSE ! neznanka je fluks (CE JE NEZNANKA FLUKS ODSTEVAMO DRUGO ENACBO !!!!)
c         dodajamo prispevke enack,  ki se ODSTEJEJO
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
              r(isqr)=r(isqr)+predznak*(-fac*(
     &               +inp%gy*(smatAbdx(i,j)-smatDx(i,j))*( temp(node) - dc(node) / inp%Rro )
     &               -inp%gx*(smatAbdy(i,j)-smatDy(i,j))*( temp(node) - dc(node) / inp%Rro ) )
     &               +facm*(
     &               +(smatAbdx(i,j)-smatDx(i,j))*kbf(node,2)
     &               -(smatAbdy(i,j)-smatDy(i,j))*kbf(node,1)))
            END DO
            END DO
          END IF
        END DO
c      END DO
      END IF

      END



      
C -----------------------------------------------------------------------------
      SUBROUTINE BouyancyRHS(eqn,mesh,inp,temp,dc,kbf,r,
     &                       smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nanoA,nanoB)
C
C     $: sets up bouyancy on right hand side
C
C -----------------------------------------------------------------------------
      USE inc_types   

      TYPE(meshType) :: mesh
      TYPE(inputtype) inp

      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z   
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc) 
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc) 
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc) 
      
      REAL(8) temp(mesh%nnodes)
      REAL(8) dc(mesh%nnodes)
      REAL(8) kbf(mesh%nnodes,3)
      REAL(8) r(mesh%nicnsp)
      REAL(8) fac,facm,nanoA,nanoB
      
      INTEGER eqn
      INTEGER i,j,it,ic,node

      fac=inp%ran/inp%prn/inp%ren*nanoB/nanoA
      facm=inp%ram/inp%prn/inp%ren*nanoB/nanoA

      IF (eqn.EQ.1) THEN
c       *** x ***
        DO j=1,mesh%npoc
          i=0
          DO ic=1,mesh%nicell
            DO it=1,51 ! 51 izvornih tock v vsaki celici
              i=i+1  ! vrstica v matriki      
                node=mesh%idc(ic,j)
                r(i)=r(i)-fac*(
     &               +inp%gz*(smatAbdy(i,j)-smatDy(i,j))*( temp(node) - dc(node) / inp%Rro )
     &               -inp%gy*(smatAbdz(i,j)-smatDz(i,j))*( temp(node) - dc(node) / inp%Rro ) )
     &               +facm*(
     &               +(smatAbdy(i,j)-smatDy(i,j))*kbf(node,3)
     &               -(smatAbdz(i,j)-smatDz(i,j))*kbf(node,2))
            END DO
          END DO
        END DO
      ELSE IF (eqn.EQ.2) THEN
c       *** y ***   
        DO j=1,mesh%npoc
          i=0
          DO ic=1,mesh%nicell
            DO it=1,51 ! 51 izvornih tock v vsaki celici
              i=i+1  ! vrstica v matriki      
                node=mesh%idc(ic,j)
                r(i)=r(i)-fac*(
     &               -inp%gz*(smatAbdx(i,j)-smatDx(i,j))*( temp(node) - dc(node) / inp%Rro )
     &               +inp%gx*(smatAbdz(i,j)-smatDz(i,j))*( temp(node) - dc(node) / inp%Rro ) )
     &               +facm*(
     &               +(smatAbdx(i,j)-smatDx(i,j))*kbf(node,3)
     &               -(smatAbdz(i,j)-smatDz(i,j))*kbf(node,1))
            END DO
          END DO
        END DO         
      ELSE
c       *** z ***   
        DO j=1,mesh%npoc
          i=0
          DO ic=1,mesh%nicell
            DO it=1,51 ! 51 izvornih tock v vsaki celici
              i=i+1  ! vrstica v matriki      
                node=mesh%idc(ic,j)
                r(i)=r(i)-fac*(
     &               +inp%gy*(smatAbdx(i,j)-smatDx(i,j))*( temp(node) - dc(node) / inp%Rro )
     &               -inp%gx*(smatAbdy(i,j)-smatDy(i,j))*( temp(node) - dc(node) / inp%Rro ) )
     &               +facm*(
     &               +(smatAbdx(i,j)-smatDx(i,j))*kbf(node,2)
     &               -(smatAbdy(i,j)-smatDy(i,j))*kbf(node,1))
            END DO
          END DO
        END DO         
      END IF
      
      END
      
      

C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsB_w_fill(eqn,mesh,smatH,smatG,smatB,beta,ren,sysm,rhsm,
     &           velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)      
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
      REAL(8) velocity(mesh%nnodes,3),val
      TYPE(matrix) :: sysm,rhsm
      
      INTEGER wqkode(mesh%nq)
      INTEGER wkode(mesh%nnodes)
      INTEGER nunk, nb, eqn
            
      INTEGER ic,it,row,col,ii,jj,nuq,j,is,ir
      
      wqkode=mesh%wqkode(:,eqn)
      wkode=mesh%wkode(:,eqn)
      nunk=mesh%nunk(eqn)
      nb=mesh%nb(eqn)

      is=0
      ir=0
      IF (eqn.EQ.1) THEN
        DO ic=1,mesh%nicell ! zanka po celicah
          DO it=1,mesh%nsp ! po izvornih tockah znotraj celice
            row=(ic-1)*mesh%nsp+it ! vrstica v sistemski matriki
C
C           najprej H matrika (u)
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
C           nato G matrika (q)
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
      ELSE IF (eqn.EQ.2) THEN
        DO ic=1,mesh%nicell ! zanka po celicah
          DO it=1,mesh%nsp ! po izvornih tockah znotraj celice
            row=(ic-1)*mesh%nsp+it ! vrstica v sistemski matriki
C
C           najprej H matrika (u)
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
C           nato G matrika (q)
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
      ELSE IF (eqn.EQ.3) THEN
        DO ic=1,mesh%nicell ! zanka po celicah
          DO it=1,mesh%nsp ! po izvornih tockah znotraj celice
            row=(ic-1)*mesh%nsp+it ! vrstica v sistemski matriki
C
C           najprej H matrika (u)
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
C           nato G matrika (q)
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
      END IF

      
      END      








C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsB_w_SQfill(eqn,mesh,smatH,smatG,smatB,beta,ren,sysm,rhsm,
     &           velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)
C
C     $: Iz pravokotnih matrik g, h in b naredi KVADRATNO CRS sistemsko in rhs
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
      REAL(8) velocity(mesh%nnodes,3),val
      TYPE(matrix) :: sysm,rhsm
      INTEGER, ALLOCATABLE :: eql(:)
      REAL(8), ALLOCATABLE :: FRsys(:),FRrhs(:)
      REAL(8) predznak

      INTEGER wqkode(mesh%nq)
      INTEGER wkode(mesh%nnodes)
      INTEGER nunk, nb, eqn

      INTEGER ic,row,col,ii,jj,nuq,j
      INTEGER i,inode,isqr,isuma,kr,ks

      wqkode=mesh%wqkode(:,eqn)
      wkode=mesh%wkode(:,eqn)
      nunk=mesh%nunk(eqn)
      nb=mesh%nb(eqn)
      ALLOCATE (eql(nunk))
      ALLOCATE (FRsys(nunk),FRrhs(nb))
      eql=mesh%Weql(:,eqn)



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
