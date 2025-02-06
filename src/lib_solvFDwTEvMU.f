
C -----------------------------------------------------------------------------
      SUBROUTINE BouyancyRHS_vMU(eqn,mesh,inp,temp,dc,kbf,r,
     &                       smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nanoA,nanoB,DvaEpsGVisc,cvte)
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

      REAL(8) DvaEpsGVisc(mesh%nnodes,3) !  2\epika\cdot\nabla\mu
      REAL(8) cvte(mesh%nnodes)
      REAL(8) temp(mesh%nnodes)
      REAL(8) dc(mesh%nnodes)
      REAL(8) kbf(mesh%nnodes,3)
      REAL(8) r(mesh%nicnsp)
      REAL(8) fac,facm,facv,nanoA,nanoB

      INTEGER eqn
      INTEGER i,j,it,ic,node

      fac=inp%ran/inp%prn/inp%ren*nanoB/nanoA
      facv=nanoB/nanoA
      facm=inp%ram/inp%prn/inp%ren*nanoB/nanoA

      IF (eqn.EQ.1) THEN
c       *** x ***
        DO j=1,mesh%npoc
          i=0
          DO ic=1,mesh%nicell
            DO it=1,51 ! 51 izvornih tock v vsaki celici
              i=i+1  ! vrstica v matriki
                node=mesh%idc(ic,j)
                r(i)=r(i)
     &               -cvte(node)*fac*(  ! bouyancy
     &               +inp%gz*(smatAbdy(i,j)-smatDy(i,j))*( temp(node) - dc(node) / inp%Rro )
     &               -inp%gy*(smatAbdz(i,j)-smatDz(i,j))*( temp(node) - dc(node) / inp%Rro ) )
     &               +facm*(  ! Kelvin body force
     &               +(smatAbdy(i,j)-smatDy(i,j))*kbf(node,3)
     &               -(smatAbdz(i,j)-smatDz(i,j))*kbf(node,2))
     &               +facv*(  ! Variable viscosity
     &               +(smatAbdy(i,j)-smatDy(i,j))*DvaEpsGVisc(node,3)
     &               -(smatAbdz(i,j)-smatDz(i,j))*DvaEpsGVisc(node,2))
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
                r(i)=r(i)
     &               -cvte(node)*fac*(  ! bouyancy
     &               -inp%gz*(smatAbdx(i,j)-smatDx(i,j))*( temp(node) - dc(node) / inp%Rro )
     &               +inp%gx*(smatAbdz(i,j)-smatDz(i,j))*( temp(node) - dc(node) / inp%Rro ) )
     &               +facm*( ! Kelvin body force
     &               +(smatAbdx(i,j)-smatDx(i,j))*kbf(node,3)
     &               -(smatAbdz(i,j)-smatDz(i,j))*kbf(node,1))
     &               +facv*( ! Variable viscosity
     &               +(smatAbdx(i,j)-smatDx(i,j))*DvaEpsGVisc(node,3)
     &               -(smatAbdz(i,j)-smatDz(i,j))*DvaEpsGVisc(node,1))
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
                r(i)=r(i)
     &               -cvte(node)*fac*(  ! bouyancy
     &               +inp%gy*(smatAbdx(i,j)-smatDx(i,j))*( temp(node) - dc(node) / inp%Rro )
     &               -inp%gx*(smatAbdy(i,j)-smatDy(i,j))*( temp(node) - dc(node) / inp%Rro ) )
     &               +facm*(! Kelvin body force
     &               +(smatAbdx(i,j)-smatDx(i,j))*kbf(node,2)
     &               -(smatAbdy(i,j)-smatDy(i,j))*kbf(node,1))
     &               +facv*( ! Variable viscosity
     &               +(smatAbdx(i,j)-smatDx(i,j))*DvaEpsGVisc(node,2)
     &               -(smatAbdy(i,j)-smatDy(i,j))*DvaEpsGVisc(node,1))
            END DO
          END DO
        END DO
      END IF

      END



C -----------------------------------------------------------------------------
C
C       Calculate values that depend on gradient of viscosity
C
      SUBROUTINE ViscGradTerms(mesh,velocity,vorticity,GradVisc,DvaEpsGVisc,CuOmCrGrVi)
C
C -----------------------------------------------------------------------------
      USE inc_types

      INTEGER i,j

      TYPE(meshType) :: mesh
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3)
      REAL(8) GradVisc(mesh%nnodes,3)

      REAL(8) DvaEpsGVisc(mesh%nnodes,3) !  2\epika\cdot\nabla\mu
      REAL(8) CuOmCrGrVi(mesh%nnodes,3)  !  \nabla\time\omega\times\nabla\mu

      REAL(8), ALLOCATABLE :: GradVx(:,:)   ! gradient of velocity component vx
      REAL(8), ALLOCATABLE :: GradVy(:,:)   ! gradient of velocity component vy
      REAL(8), ALLOCATABLE :: GradVz(:,:)   ! gradient of velocity component vz
      REAL(8), ALLOCATABLE :: CurlOfVort(:,:) ! curl of vorticity \nabla\times\omega
      REAL(8), ALLOCATABLE :: tmp(:,:) ! tmp vector

C
C     Allocate temporary variables
C
      ALLOCATE( GradVx(mesh%nnodes,3) )
      ALLOCATE( GradVy(mesh%nnodes,3) )
      ALLOCATE( GradVz(mesh%nnodes,3) )
      ALLOCATE( CurlOfVort(mesh%nnodes,3) )
      ALLOCATE( tmp(mesh%nnodes,3) )
C
C     Calculate
C

C     velocity gradient
      CALL setGrad(mesh,Velocity(:,1),GradVx)
      CALL setGrad(mesh,Velocity(:,2),GradVy)
      CALL setGrad(mesh,Velocity(:,3),GradVz)

c     2\epika\cdot\nabla\mu

      CALL SetDotProduct(mesh,GradVisc,GradVx,DvaEpsGVisc(:,1))
      CALL SetDotProduct(mesh,GradVisc,GradVy,DvaEpsGVisc(:,2))
      CALL SetDotProduct(mesh,GradVisc,GradVz,DvaEpsGVisc(:,3))

      CALL SetCrossProduct(mesh,tmp,GradVisc,vorticity)

      DO i=1,mesh%nnodes
        DO j=1,3
          DvaEpsGVisc(i,j)=2.0D0*DvaEpsGVisc(i,j)+tmp(i,j)
        END DO
      END DO

c     calculate curl of vorticity
      CALL setCurl(mesh,Vorticity,CurlOfVort)
c     \nabla\time\omega\times\nabla\mu
      CALL setCrossProduct(mesh,CuOmCrGrVi,CurlOfVort,GradVisc)

      DEALLOCATE(CurlOfVort,GradVx,GradVy,GradVz,tmp)

      END

C -----------------------------------------------------------------------------
C
C       Solve vorticity transport equation using finite difference discretization
C       of time derivative, variable material properties
C
C
      SUBROUTINE SolveFDwTE_vMU(eqn,env,io,inp,cpu,mesh,sysm,rhsm,precv,
     &                      vorticity,qvorticity,temp,dc,velocity,kbf,rhsvi,smatB,
     &                      smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nit,nanoA,nanoB,
     &                      DvaEpsGVisc,CuOmCrGrVi,cvte)
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

      REAL(8) rhsvi(mesh%nnodes)
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
      REAL(8) cvte(mesh%nnodes)
      REAL(8) nanoA,nanoB

      REAL(8) DvaEpsGVisc(mesh%nnodes,3) !  2\epika\cdot\nabla\mu
      REAL(8) CuOmCrGrVi(mesh%nnodes,3)  !  \nabla\time\omega\times\nabla\mu

      REAL(8), ALLOCATABLE :: x(:), b(:), r(:), rhsv(:)
      REAL(4) cptime,rcpu
C
C     Start time measurement
C
      rcpu=cptime(0.)
C
C     Add variable material properties contribution to rhs vector
C
      ALLOCATE (rhsv(mesh%nnodes))
      DO i=1,mesh%nnodes
        rhsv(i)=rhsvi(i)-CuOmCrGrVi(i,eqn)/inp%ren*nanoA
      END DO
C
C     Allocate variables
C
      ALLOCATE (x(mesh%nunk(eqn)),b(mesh%nb(eqn)),r(mesh%nicnsp))
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
c     Bouyancy, kelvin body force, variable viscosity
c
      CALL BouyancyRHS_vMU(eqn,mesh,inp,temp,dc,kbf,r,
     &                   smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nanoA,nanoB,DvaEpsGVisc,cvte)

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

      DEALLOCATE (x,b,r,rhsv)

      END





C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsB_w_fill_vMU(eqn,mesh,smatH,smatG,smatB,beta,ren,sysm,rhsm,
     &           velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,visc,qvisc)
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

      REAL(8) visc(mesh%nnodes)
      REAL(8) qvisc(mesh%nq)

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
              val=visc(j)*smatH(row,col)+smatB(row,col)*beta*ren
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
                  sysm%v(sysm%d(is))=-smatG(row,col)*qvisc(j)
                ELSE ! znana vrednost - rhs
                  ir=ir+1
                  rhsm%v(rhsm%d(ir))=smatG(row,col)*qvisc(j) ! minus zato, ker gre na drugo stran enacbe
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
              val=visc(j)*smatH(row,col)+smatB(row,col)*beta*ren
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
                  sysm%v(sysm%d(is))=-smatG(row,col)*qvisc(j)
                ELSE ! znana vrednost - rhs
                  ir=ir+1
                  rhsm%v(rhsm%d(ir))=smatG(row,col)*qvisc(j) ! minus zato, ker gre na drugo stran enacbe
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
              val=visc(j)*smatH(row,col)+smatB(row,col)*beta*ren
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
                  sysm%v(sysm%d(is))=-smatG(row,col)*qvisc(j)
                ELSE ! znana vrednost - rhs
                  ir=ir+1
                  rhsm%v(rhsm%d(ir))=smatG(row,col)*qvisc(j) ! minus zato, ker gre na drugo stran enacbe
                END IF
              END DO
            END DO
          END DO
        END DO
      END IF

      END





C -----------------------------------------------------------------------------
      SUBROUTINE SetUpViscCvte(env,io,inp,mesh,viscosity,cvte,temp)
C
C     $: sets up variable material properties
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(inputtype) :: inp
      TYPE(IOtype)    :: io
      TYPE(penv) :: env

      REAL(8) z
      INTEGER i

      REAL(8) viscosity(mesh%nnodes)
      REAL(8) cvte(mesh%nnodes)
      REAL(8) temp(mesh%nnodes)

      IF (inp%veqm.EQ.1) THEN
        DO i=1,mesh%nnodes
          viscosity(i)= 1.0D0
          cvte(i)     = 1.0D0
        END DO

      ELSE IF (inp%veqm.EQ.11) THEN
        DO i=1,mesh%nnodes
          z=mesh%x(i,3)
          viscosity(i)= 1.0D0+z**2.0D0
          cvte(i)     = -12.0D0*(1.0D0-3*z)*z
          inp%gx=1.0D0 ! ker test, uporabim ctve kot source
          inp%gy=0.0D0 ! glej Mathematico!
          inp%gz=0.0D0
          temp(i)=-1.0D0 ! ker je vzgonski clen z minusom (glej 2.10 v doktoratu)
        END DO

      ELSE IF (inp%veqm.EQ.12) THEN
        DO i=1,mesh%nnodes
          z=mesh%x(i,3)
          viscosity(i)= 1.0D0+z
          cvte(i)     = 24.0D0*z-6.0D0
          inp%gx=1.0D0 ! ker test, uporabim ctve kot source
          inp%gy=0.0D0 ! glej Mathematico!
          inp%gz=0.0D0
          temp(i)=-1.0D0 ! ker je vzgonski clen z minusom (glej 2.10 v doktoratu)
        END DO


      ELSE
        CALL WarnErr(env,io,inp,4,"SetUpViscCvte","wrong material model",inp%veqm)
      END IF

      END





C -----------------------------------------------------------------------------
      SUBROUTINE SetCurl(mesh,a,b)
C
C     $: b = curl(a)
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh

      REAL(8) a(mesh%nnodes,3)
      REAL(8) b(mesh%nnodes,3)

      REAL(8), ALLOCATABLE :: gradAx(:,:)
      REAL(8), ALLOCATABLE :: gradAy(:,:)
      REAL(8), ALLOCATABLE :: gradAz(:,:)

      INTEGER i

      ALLOCATE (gradAx(mesh%nnodes,3))
      ALLOCATE (gradAy(mesh%nnodes,3))
      ALLOCATE (gradAz(mesh%nnodes,3))

      CALL setGrad(mesh,a(:,1),GradAx)
      CALL setGrad(mesh,a(:,2),GradAy)
      CALL setGrad(mesh,a(:,3),GradAz)

      DO i=1,mesh%nnodes
        b(i,1)= + gradAz(i,2) - gradAy(i,3)
        b(i,2)= - gradAz(i,1) - gradAx(i,3)
        b(i,3)= + gradAy(i,1) - gradAx(i,2)
      END DO

      DEALLOCATE (gradAx,gradAy,gradAz)

      END


