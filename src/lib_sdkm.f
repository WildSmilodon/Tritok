C
C
C               Single Domain Kinematics - lib_sdkm.f
C
C
c
c       List of subroutines
c
c       sdkm_solve
c       sdkm_rhs
c       sdkm_initguess
c       sdkm_distunk
c       sdKm_sysm
c       FMATsdkm
c       INTEBc9sdkm
c       sdKm_corInt
c       CheckSDIntFile
c       WriteSDIntDisk
c       ReadSDIntDisk


C -----------------------------------------------------------------------------
      SUBROUTINE sdkm_solve(eqn,env,inp,cpu,mesh,lnits,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &                      sdA,sdPivot,velocity,vorticity)
C
C     $: Solve EQNth single domain kinematics equation
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      TYPE(inputtype) inp 
      TYPE(CPUtype) cpu     
      TYPE(penv) :: env         
      
      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion
      REAL(8) sdA(env%zac:env%kon,mesh%nbnodes)  ! system matrix
      REAL(8) sdPivot(env%zac:env%kon)           ! system matrix      

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral 
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni
      
      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDy(env%zac:env%kon,mesh%nnodes) ! to je za drugi integral 
      REAL(8) sdDz(env%zac:env%kon,mesh%nnodes) ! na desni 
      
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3)  
      INTEGER lslvt,lpret,lprep,lmaxit,lstopt,lnits
      REAL(4) lslveps
      
      REAL(8), ALLOCATABLE :: rhsv(:),xsol(:)
      
      INTEGER eqn,ierr,i
      REAL(4) cptime,rcpu,rrcpu
      
      ALLOCATE (rhsv(env%zac:env%kon),xsol(mesh%nbnodes))
      
      rcpu=cptime(0.)

      IF (env%nproc.EQ.1) THEN
        lslvt=0
        lpret=2
        lprep=2
        lmaxit=500
        lstopt=5
        lslveps=1.0E-15
        ierr=0
      ELSE
        lslvt=2
        lpret=1
        lprep=2
        lmaxit=inp%pkms_maxit
        lstopt=5
        lslveps=inp%pkms_eps
        ierr=0
      END IF

c
c     set up initial guess for solver
c      
c      CALL sdkm_initguess(mesh,vorticity(:,eqn),xsol)
      xsol=0.0D0
C
C     Set up r.h.s.
C
      CALL sdkm_rhs(eqn,env,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &              velocity,vorticity,rhsv)      
      cpu%time((eqn-1)*6+1)=cpu%time((eqn-1)*6+1)+cptime(rcpu)
      rcpu=cptime(0.)
c
c     solve
c      
      IF (env%nproc.EQ.1) THEN      
        CALL SolvEQNfm(lslvt,lpret,lprep,lmaxit,lstopt,lslveps,
     &              mesh%nbnodes,sdPivot,sdA,rhsv,xsol,lnits,rrcpu,ierr)
      ELSE
        CALL pSolvEQNfm(env,lslvt,lpret,lprep,lmaxit,lstopt,lslveps,
     &              mesh%nbnodes,env%nmn,sdPivot,sdA,rhsv,xsol,lnits,rrcpu,ierr)
      END IF     
      cpu%time((eqn-1)*6+2)=cpu%time((eqn-1)*6+2)+cptime(rcpu)
c
c     copy solution to vorticity vector
c
      CALL sdkm_distunk(eqn,inp,mesh,vorticity(:,eqn),xsol)  ! tu bi morda podrelaksiral!!      
      
      DEALLOCATE (rhsv,xsol)

C     Employ periodic boundary condition
c      IF (mesh%iPeri.GT.0) THEN
c        DO i=1,mesh%nPeri
c          vorticity(mesh%pwFT(i,2),eqn)=vorticity(mesh%pwFT(i,1),eqn)
c        END DO
c      END IF
      
      END



C -----------------------------------------------------------------------------
      SUBROUTINE sdkm_solveBCN(eqn,env,inp,cpu,mesh,lnits,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &                      sdA,sdPivot,velocity,vorticity)
C
C     $: Solve EQNth single domain kinematics equation
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(inputtype) inp
      TYPE(CPUtype) cpu
      TYPE(penv) :: env

      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion
      REAL(8) sdA(env%zac:env%kon,mesh%nbnodes)  ! system matrix
      REAL(8) sdPivot(env%zac:env%kon)           ! system matrix

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni

      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDy(env%zac:env%kon,mesh%nnodes) ! to je za drugi integral
      REAL(8) sdDz(env%zac:env%kon,mesh%nnodes) ! na desni

      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3)
      INTEGER lslvt,lpret,lprep,lmaxit,lstopt,lnits
      REAL(4) lslveps

      REAL(8), ALLOCATABLE :: rhsv(:),xsol(:)

      INTEGER eqn,ierr,i
      REAL(4) cptime,rcpu,rrcpu

      ALLOCATE (rhsv(env%zac:env%kon),xsol(mesh%nbnodes))

      rcpu=cptime(0.)

      IF (env%nproc.EQ.1) THEN
        lslvt=0
        lpret=2
        lprep=2
        lmaxit=500
        lstopt=5
        lslveps=1.0E-15
        ierr=0
      ELSE
        lslvt=2
        lpret=1
        lprep=2
        lmaxit=inp%pkms_maxit
        lstopt=5
        lslveps=inp%pkms_eps
        ierr=0
      END IF

c
c     set up initial guess for solver
c
c      CALL sdkm_initguess(mesh,vorticity(:,eqn),xsol)
      xsol=0.0D0
C
C     Set up r.h.s.
C
      CALL sdkm_rhsBCN(eqn,env,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &              velocity,vorticity,rhsv)
      cpu%time((eqn-1)*6+1)=cpu%time((eqn-1)*6+1)+cptime(rcpu)
      rcpu=cptime(0.)
c
c     solve
c
      IF (env%nproc.EQ.1) THEN
        CALL SolvEQNfm(lslvt,lpret,lprep,lmaxit,lstopt,lslveps,
     &              mesh%nbnodes,sdPivot,sdA,rhsv,xsol,lnits,rrcpu,ierr)
      ELSE
        CALL pSolvEQNfm(env,lslvt,lpret,lprep,lmaxit,lstopt,lslveps,
     &              mesh%nbnodes,env%nmn,sdPivot,sdA,rhsv,xsol,lnits,rrcpu,ierr)
      END IF
      cpu%time((eqn-1)*6+2)=cpu%time((eqn-1)*6+2)+cptime(rcpu)
c
c     copy solution to vorticity vector
c
      CALL sdkm_distunk(eqn,inp,mesh,vorticity(:,eqn),xsol)  ! tu bi morda podrelaksiral!!

      DEALLOCATE (rhsv,xsol)

C     Employ periodic boundary condition
c      IF (mesh%iPeri.GT.0) THEN
c        DO i=1,mesh%nPeri
c          vorticity(mesh%pwFT(i,2),eqn)=vorticity(mesh%pwFT(i,1),eqn)
c        END DO
c      END IF

      END


C -----------------------------------------------------------------------------
      SUBROUTINE testw(mesh,Vorticity)

C
C     $: preverimo koliko je normalna komponenta vrtincnosti
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh

      REAL(8) vorticity(mesh%nnodes,3)

      INTEGER i,j,k

      DO i=1,mesh%nbnodes
        j=mesh%gbn(i)
        print *,i,mesh%iwn(i),
     & mesh%RotMat(i)%v(1,1)*Vorticity(j,1)+mesh%RotMat(i)%v(1,2)*Vorticity(j,2)+mesh%RotMat(i)%v(1,3)*Vorticity(j,3)
      END DO


      END

C -----------------------------------------------------------------------------
      SUBROUTINE sdkm_RotSolve(env,inp,cpu,mesh,lnits,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &                      sdA,sdPivot,velocity,vorticity)
C
C     $: Solve EQNth single domain kinematics equation
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(inputtype) inp
      TYPE(CPUtype) cpu
      TYPE(penv) :: env

      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion
      REAL(8) sdA(env%zac:env%kon,mesh%nbnodes)  ! system matrix
      REAL(8) sdPivot(env%zac:env%kon)           ! system matrix

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni

      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDy(env%zac:env%kon,mesh%nnodes) ! to je za drugi integral
      REAL(8) sdDz(env%zac:env%kon,mesh%nnodes) ! na desni

      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3)
      INTEGER lslvt,lpret,lprep,lmaxit,lstopt,lnits
      REAL(4) lslveps

      REAL(8), ALLOCATABLE :: rhsv1(:),xsol1(:)
      REAL(8), ALLOCATABLE :: rhsv2(:),xsol2(:),xsolNormal(:)

      INTEGER ierr,i
      REAL(4) cptime,rcpu,rrcpu

      ALLOCATE (rhsv1(env%zac:env%kon),xsol1(mesh%nbnodes))
      ALLOCATE (rhsv2(env%zac:env%kon),xsol2(mesh%nbnodes))
      ALLOCATE (xsolNormal(mesh%nbnodes))

      rcpu=cptime(0.)

      IF (env%nproc.EQ.1) THEN
        lslvt=0
        lpret=2
        lprep=2
        lmaxit=500
        lstopt=5
        lslveps=1.0E-15
        ierr=0
      ELSE
        lslvt=2
        lpret=1
        lprep=2
        lmaxit=inp%pkms_maxit
        lstopt=5
        lslveps=inp%pkms_eps
        ierr=0
      END IF

c
c     set up initial guess for solver
c
      xsol1=0.0D0
      xsol2=0.0D0
C
C     Set up r.h.s.
C
      print *,"ww0"
      CALL sdkm_RotRhs(env,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &              velocity,vorticity,rhsv1,rhsv2,xsolNormal)
      cpu%time((1-1)*6+1)=cpu%time((1-1)*6+1)+cptime(rcpu)
      rcpu=cptime(0.)
c
c     solve
c
      print *,"ww"
      IF (env%nproc.EQ.1) THEN
        CALL SolvEQNfm(lslvt,lpret,lprep,lmaxit,lstopt,lslveps,
     &              mesh%nbnodes,sdPivot,sdA,rhsv1,xsol1,lnits,rrcpu,ierr)
      ELSE
        CALL pSolvEQNfm(env,lslvt,lpret,lprep,lmaxit,lstopt,lslveps,
     &              mesh%nbnodes,env%nmn,sdPivot,sdA,rhsv1,xsol1,lnits,rrcpu,ierr)
      END IF
      print *,"ww1"
      IF (env%nproc.EQ.1) THEN
        CALL SolvEQNfm(lslvt,lpret,lprep,lmaxit,lstopt,lslveps,
     &              mesh%nbnodes,sdPivot,sdA,rhsv2,xsol2,lnits,rrcpu,ierr)
      ELSE
        CALL pSolvEQNfm(env,lslvt,lpret,lprep,lmaxit,lstopt,lslveps,
     &              mesh%nbnodes,env%nmn,sdPivot,sdA,rhsv2,xsol2,lnits,rrcpu,ierr)
      END IF
      cpu%time((1-1)*6+2)=cpu%time((1-1)*6+2)+cptime(rcpu)
c
c     copy solution to vorticity (or velocity) vector
c
      CALL sdkm_Rotdistunk(inp,mesh,vorticity,velocity,xsol1,xsol2,xsolNormal)

      DEALLOCATE (rhsv1,xsol1)
      DEALLOCATE (rhsv2,xsol2,xsolNormal)

      END


C -----------------------------------------------------------------------------
      SUBROUTINE sdkm_rhs(eqn,env,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &                    velocity,vorticity,dfi)
C
C     $: Form r.h.s. for EQNth single domain kinematics equation
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      TYPE(penv) :: env

      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion
 
      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral 
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni
      
      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDy(env%zac:env%kon,mesh%nnodes) ! to je za drugi integral 
      REAL(8) sdDz(env%zac:env%kon,mesh%nnodes) ! na desni 
      
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3) 
      REAL(8) dfi(env%zac:env%kon) 
      
      INTEGER eqn,i,j,k
C
C____ INITIALIZE ARRAY :
C      
      dfi=0.00D0

      IF (eqn.EQ.1) THEN ! x eqaution
C
C____ LOOP OVER COLUMNS :
C
      DO k=1,mesh%nnodes
        j=mesh%lbn(k)
C
C....   CONSIDER BOUNDARY VORTICITY CONTRIBUTION :
C
        IF (j.NE.0) THEN                                  ! BOUNDARY
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdH(i,j)*(mesh%ny(i)*velocity(k,3)-mesh%nz(i)*velocity(k,2))
c     &                       +sdH(i,j)*mesh%ny(i)*velocity(k,3)
c     &                       -sdH(i,j)*mesh%nz(i)*velocity(k,2)
     &                       -sdHtx(i,j)*(mesh%ny(i)*velocity(k,2)+mesh%nz(i)*velocity(k,3))
c     &                       -sdHtx(i,j)*mesh%ny(i)*velocity(k,2)
c     &                       -sdHtx(i,j)*mesh%nz(i)*velocity(k,3)     
     &                       +velocity(k,1)*(sdHty(i,j)*mesh%ny(i)+sdHtz(i,j)*mesh%nz(i))
c     &                       +sdHty(i,j)*mesh%ny(i)*velocity(k,1)
c     &                       +sdHtz(i,j)*mesh%nz(i)*velocity(k,1)
     &                       +sdDx(i,k)*(
     &                                  +mesh%nx(i)*vorticity(k,1)
     &                                  +mesh%ny(i)*vorticity(k,2)
     &                                  +mesh%nz(i)*vorticity(k,3))          
          END DO !loop over rows
C
C....   CONSIDER DOMAIN VORTICITY CONTRIBUTION :
C
        ELSE                                              ! DOMAIN
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdDx(i,k)*(mesh%ny(i)*vorticity(k,2)+mesh%nz(i)*vorticity(k,3))
     &                       -vorticity(k,1)*(sdDy(i,k)*mesh%ny(i)+sdDz(i,k)*mesh%nz(i))
c     &                       +sdDx(i,k)*mesh%ny(i)*vorticity(k,2)
c     &                       +sdDx(i,k)*mesh%nz(i)*vorticity(k,3)
c     &                       -sdDy(i,k)*mesh%ny(i)*vorticity(k,1)
c     &                       -sdDz(i,k)*mesh%nz(i)*vorticity(k,1)
          END DO !loop over rows
        END IF !boundary or domain node
      END DO !loop over columns
      ELSE IF (eqn.EQ.2) THEN ! y eqaution
C
C____ LOOP OVER COLUMNS :
C
      DO k=1,mesh%nnodes
        j=mesh%lbn(k)
C
C....   CONSIDER BOUNDARY VORTICITY CONTRIBUTION :
C
        IF (j.NE.0) THEN                                  ! BOUNDARY
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdH(i,j)*(mesh%nz(i)*velocity(k,1)-mesh%nx(i)*velocity(k,3))
c     &                       +sdH(i,j)*mesh%nz(i)*velocity(k,1)
c     &                       -sdH(i,j)*mesh%nx(i)*velocity(k,3)
     &                       -sdHty(i,j)*(mesh%nx(i)*velocity(k,1)+mesh%nz(i)*velocity(k,3))
c     &                       -sdHty(i,j)*mesh%nx(i)*velocity(k,1)
c     &                       -sdHty(i,j)*mesh%nz(i)*velocity(k,3)     
     &                       +velocity(k,2)*(sdHtx(i,j)*mesh%nx(i)+sdHtz(i,j)*mesh%nz(i))
c     &                       +sdHtx(i,j)*mesh%nx(i)*velocity(k,2)
c     &                       +sdHtz(i,j)*mesh%nz(i)*velocity(k,2)
     &                       +sdDy(i,k)*(
     &                                   mesh%nx(i)*vorticity(k,1)
     &                                  +mesh%ny(i)*vorticity(k,2)
     &                                  +mesh%nz(i)*vorticity(k,3))          
          END DO !loop over rows
C
C....   CONSIDER DOMAIN VORTICITY CONTRIBUTION :
C
        ELSE                                              ! DOMAIN
          DO i=env%zac,env%kon      
                dfi(i)=dfi(i)
     &                       +sdDy(i,k)*(mesh%nx(i)*vorticity(k,1)+mesh%nz(i)*vorticity(k,3))
c     &                       +sdDy(i,k)*mesh%nx(i)*vorticity(k,1)
c     &                       +sdDy(i,k)*mesh%nz(i)*vorticity(k,3)
     &                       -vorticity(k,2)*(sdDx(i,k)*mesh%nx(i)+sdDz(i,k)*mesh%nz(i))
c     &                       -sdDx(i,k)*mesh%nx(i)*vorticity(k,2)
c     &                       -sdDz(i,k)*mesh%nz(i)*vorticity(k,2)
          END DO !loop over rows
        END IF !boundary or domain node
      END DO !loop over columns

      ELSE IF (eqn.EQ.3) THEN ! z eqaution
C
C____ LOOP OVER COLUMNS :
C
      DO k=1,mesh%nnodes
        j=mesh%lbn(k)
C
C....   CONSIDER BOUNDARY VORTICITY CONTRIBUTION :
C
        IF (j.NE.0) THEN                                  ! BOUNDARY
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdH(i,j)*(mesh%nx(i)*velocity(k,2)-mesh%ny(i)*velocity(k,1))
c     &                       +sdH(i,j)*mesh%nx(i)*velocity(k,2)
c     &                       -sdH(i,j)*mesh%ny(i)*velocity(k,1)
     &                       -sdHtz(i,j)*(mesh%nx(i)*velocity(k,1)+mesh%ny(i)*velocity(k,2))
c     &                       -sdHtz(i,j)*mesh%nx(i)*velocity(k,1)
c     &                       -sdHtz(i,j)*mesh%ny(i)*velocity(k,2)     
     &                       +velocity(k,3)*(sdHtx(i,j)*mesh%nx(i)+sdHty(i,j)*mesh%ny(i))
c     &                       +sdHtx(i,j)*mesh%nx(i)*velocity(k,3)
c     &                       +sdHty(i,j)*mesh%ny(i)*velocity(k,3)
     &                       +sdDz(i,k)*(
     &                                   mesh%nx(i)*vorticity(k,1)
     &                                  +mesh%ny(i)*vorticity(k,2)
     &                                  +mesh%nz(i)*vorticity(k,3))          
          END DO !loop over rows
C
C....   CONSIDER DOMAIN VORTICITY CONTRIBUTION :
C
        ELSE                                              ! DOMAIN
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdDz(i,k)*(mesh%nx(i)*vorticity(k,1)+mesh%ny(i)*vorticity(k,2))
c     &                       +sdDz(i,k)*mesh%nx(i)*vorticity(k,1)
c     &                       +sdDz(i,k)*mesh%ny(i)*vorticity(k,2)
     &                       -vorticity(k,3)*(sdDx(i,k)*mesh%nx(i)+sdDy(i,k)*mesh%ny(i))
c     &                       -sdDx(i,k)*mesh%nx(i)*vorticity(k,3)
c     &                       -sdDy(i,k)*mesh%ny(i)*vorticity(k,3)
          END DO !loop over rows
        END IF !boundary or domain node
      END DO !loop over columns

      END IF

C     resitev za bug; Ce ne resujem enacbe, dam na desno vrtincnost.
      DO i=env%zac,env%kon
        IF (mesh%kode(i,eqn).EQ.1.OR.mesh%kode(i,eqn).EQ.2) THEN   ! NE RESUJEM TE ENACBE        
          dfi(i)=vorticity(mesh%gbn(i),eqn)  
        END IF
      END DO

           
      END



C -----------------------------------------------------------------------------
      SUBROUTINE sdkm_rhsBCN(eqn,env,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &                    velocity,vorticity,dfi)
C
C     $: Form r.h.s. for EQNth single domain kinematics equation
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(penv) :: env

      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni

      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDy(env%zac:env%kon,mesh%nnodes) ! to je za drugi integral
      REAL(8) sdDz(env%zac:env%kon,mesh%nnodes) ! na desni

      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3)
      REAL(8) dfi(env%zac:env%kon)

      INTEGER eqn,i,j,k
C
C____ INITIALIZE ARRAY :
C
      dfi=0.00D0

      IF (eqn.EQ.1) THEN ! x eqaution
C
C____ LOOP OVER COLUMNS :
C
      DO k=1,mesh%nnodes
        j=mesh%lbn(k)
C
C....   CONSIDER BOUNDARY VORTICITY CONTRIBUTION :
C
        IF (j.NE.0) THEN                                  ! BOUNDARY
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdH(i,j)*(mesh%ny(i)*velocity(k,3)-mesh%nz(i)*velocity(k,2))
c     &                       +sdH(i,j)*mesh%ny(i)*velocity(k,3)
c     &                       -sdH(i,j)*mesh%nz(i)*velocity(k,2)
     &                       -sdHtx(i,j)*(mesh%ny(i)*velocity(k,2)+mesh%nz(i)*velocity(k,3))
c     &                       -sdHtx(i,j)*mesh%ny(i)*velocity(k,2)
c     &                       -sdHtx(i,j)*mesh%nz(i)*velocity(k,3)
     &                       +velocity(k,1)*(sdHty(i,j)*mesh%ny(i)+sdHtz(i,j)*mesh%nz(i))
c     &                       +sdHty(i,j)*mesh%ny(i)*velocity(k,1)
c     &                       +sdHtz(i,j)*mesh%nz(i)*velocity(k,1)
     &                       +sdDx(i,k)*(
     &                                  +mesh%nx(i)*vorticity(k,1)
     &                                  +mesh%ny(i)*vorticity(k,2)
     &                                  +mesh%nz(i)*vorticity(k,3))
          END DO !loop over rows
C
C....   CONSIDER DOMAIN VORTICITY CONTRIBUTION :
C
        ELSE                                              ! DOMAIN
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdDx(i,k)*(mesh%ny(i)*vorticity(k,2)+mesh%nz(i)*vorticity(k,3))
     &                       -vorticity(k,1)*(sdDy(i,k)*mesh%ny(i)+sdDz(i,k)*mesh%nz(i))
c     &                       +sdDx(i,k)*mesh%ny(i)*vorticity(k,2)
c     &                       +sdDx(i,k)*mesh%nz(i)*vorticity(k,3)
c     &                       -sdDy(i,k)*mesh%ny(i)*vorticity(k,1)
c     &                       -sdDz(i,k)*mesh%nz(i)*vorticity(k,1)
          END DO !loop over rows
        END IF !boundary or domain node
      END DO !loop over columns
      ELSE IF (eqn.EQ.2) THEN ! y eqaution
C
C____ LOOP OVER COLUMNS :
C
      DO k=1,mesh%nnodes
        j=mesh%lbn(k)
C
C....   CONSIDER BOUNDARY VORTICITY CONTRIBUTION :
C
        IF (j.NE.0) THEN                                  ! BOUNDARY
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdH(i,j)*(mesh%nz(i)*velocity(k,1)-mesh%nx(i)*velocity(k,3))
c     &                       +sdH(i,j)*mesh%nz(i)*velocity(k,1)
c     &                       -sdH(i,j)*mesh%nx(i)*velocity(k,3)
     &                       -sdHty(i,j)*(mesh%nx(i)*velocity(k,1)+mesh%nz(i)*velocity(k,3))
c     &                       -sdHty(i,j)*mesh%nx(i)*velocity(k,1)
c     &                       -sdHty(i,j)*mesh%nz(i)*velocity(k,3)
     &                       +velocity(k,2)*(sdHtx(i,j)*mesh%nx(i)+sdHtz(i,j)*mesh%nz(i))
c     &                       +sdHtx(i,j)*mesh%nx(i)*velocity(k,2)
c     &                       +sdHtz(i,j)*mesh%nz(i)*velocity(k,2)
     &                       +sdDy(i,k)*(
     &                                   mesh%nx(i)*vorticity(k,1)
     &                                  +mesh%ny(i)*vorticity(k,2)
     &                                  +mesh%nz(i)*vorticity(k,3))
          END DO !loop over rows
C
C....   CONSIDER DOMAIN VORTICITY CONTRIBUTION :
C
        ELSE                                              ! DOMAIN
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdDy(i,k)*(mesh%nx(i)*vorticity(k,1)+mesh%nz(i)*vorticity(k,3))
c     &                       +sdDy(i,k)*mesh%nx(i)*vorticity(k,1)
c     &                       +sdDy(i,k)*mesh%nz(i)*vorticity(k,3)
     &                       -vorticity(k,2)*(sdDx(i,k)*mesh%nx(i)+sdDz(i,k)*mesh%nz(i))
c     &                       -sdDx(i,k)*mesh%nx(i)*vorticity(k,2)
c     &                       -sdDz(i,k)*mesh%nz(i)*vorticity(k,2)
          END DO !loop over rows
        END IF !boundary or domain node
      END DO !loop over columns

      ELSE IF (eqn.EQ.3) THEN ! z eqaution
C
C____ LOOP OVER COLUMNS :
C
      DO k=1,mesh%nnodes
        j=mesh%lbn(k)
C
C....   CONSIDER BOUNDARY VORTICITY CONTRIBUTION :
C
        IF (j.NE.0) THEN                                  ! BOUNDARY
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdH(i,j)*(mesh%nx(i)*velocity(k,2)-mesh%ny(i)*velocity(k,1))
c     &                       +sdH(i,j)*mesh%nx(i)*velocity(k,2)
c     &                       -sdH(i,j)*mesh%ny(i)*velocity(k,1)
     &                       -sdHtz(i,j)*(mesh%nx(i)*velocity(k,1)+mesh%ny(i)*velocity(k,2))
c     &                       -sdHtz(i,j)*mesh%nx(i)*velocity(k,1)
c     &                       -sdHtz(i,j)*mesh%ny(i)*velocity(k,2)
     &                       +velocity(k,3)*(sdHtx(i,j)*mesh%nx(i)+sdHty(i,j)*mesh%ny(i))
c     &                       +sdHtx(i,j)*mesh%nx(i)*velocity(k,3)
c     &                       +sdHty(i,j)*mesh%ny(i)*velocity(k,3)
     &                       +sdDz(i,k)*(
     &                                   mesh%nx(i)*vorticity(k,1)
     &                                  +mesh%ny(i)*vorticity(k,2)
     &                                  +mesh%nz(i)*vorticity(k,3))
          END DO !loop over rows
C
C....   CONSIDER DOMAIN VORTICITY CONTRIBUTION :
C
        ELSE                                              ! DOMAIN
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdDz(i,k)*(mesh%nx(i)*vorticity(k,1)+mesh%ny(i)*vorticity(k,2))
c     &                       +sdDz(i,k)*mesh%nx(i)*vorticity(k,1)
c     &                       +sdDz(i,k)*mesh%ny(i)*vorticity(k,2)
     &                       -vorticity(k,3)*(sdDx(i,k)*mesh%nx(i)+sdDy(i,k)*mesh%ny(i))
c     &                       -sdDx(i,k)*mesh%nx(i)*vorticity(k,3)
c     &                       -sdDy(i,k)*mesh%ny(i)*vorticity(k,3)
          END DO !loop over rows
        END IF !boundary or domain node
      END DO !loop over columns

      END IF

C     resitev za bug; Ce ne resujem enacbe, dam na desno vrtincnost.
      DO i=env%zac,env%kon
        IF (mesh%kode(i,eqn).EQ.1.OR.mesh%kode(i,eqn).EQ.2) THEN   ! NE RESUJEM TE ENACBE
          dfi(i)=vorticity(mesh%gbn(i),eqn)
        ELSE IF (mesh%bcnl(i,1).EQ.eqn) THEN ! enacba eqn je normalna, jo racunam posredno
          dfi(i)=(mesh%wnwall(mesh%iwn(i)) !to je vrednost normalne vrtincnosti na robu
     &           -mesh%bcn(i,mesh%bcnl(i,2))*vorticity(mesh%gbn(i),mesh%bcnl(i,2))
     &           -mesh%bcn(i,mesh%bcnl(i,3))*vorticity(mesh%gbn(i),mesh%bcnl(i,3))
     &          )/mesh%bcn(i,mesh%bcnl(i,1))
        END IF
      END DO


      END


C -----------------------------------------------------------------------------
      SUBROUTINE sdkm_RotRhs(env,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &                    OriVelocity,OriVorticity,dfi1,dfi2,dfi3)

C
C     $: Form r.h.s. for EQNth single domain kinematics equation
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(penv) :: env

      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni

      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDy(env%zac:env%kon,mesh%nnodes) ! to je za drugi integral
      REAL(8) sdDz(env%zac:env%kon,mesh%nnodes) ! na desni

      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3)
      REAL(8) OriVelocity(mesh%nnodes,3)
      REAL(8) OriVorticity(mesh%nnodes,3)
      REAL(8) dfi1(env%zac:env%kon)
      REAL(8) dfi2(env%zac:env%kon)
      REAL(8) dfi3(env%zac:env%kon)

      INTEGER i,j,k
C
C____ INITIALIZE ARRAY :
C
      dfi1=0.00D0
      dfi2=0.00D0
C
C____ CREATE RIGHT HAND SIDE VECTOR
C
      DO i=env%zac,env%kon ! LOOP OVER ROWS
C       Rotate veloctiy and vorticity fields
        CALL RotateVecByR(OriVelocity,velocity,mesh%RotMat(i)%v,mesh%nnodes)
        CALL RotateVecByR(OriVorticity,vorticity,mesh%RotMat(i)%v,mesh%nnodes)

c          print *, mesh%RotMat(i)%v

c        do ii=1,mesh%nbnodes
c          if (mesh%bkmc(ii).Eq.1) then
c          print *,mesh%iwn(ii),mesh%bkmc(ii)
c          print *,velocity(mesh%gbn(ii),1),velocity(mesh%gbn(ii),2),velocity(mesh%gbn(ii),3)
c          print *,vorticity(mesh%gbn(ii),1),vorticity(mesh%gbn(ii),2),vorticity(mesh%gbn(ii),3)
c          end if
c        end do
c        stop

        IF (mesh%bkmc(i).EQ.1) THEN ! TANGENTIAL BOUNDARY VORTICTIES ARE UNKNOWN

          DO k=1,mesh%nnodes ! LOOP OVER COLUMNS
            j=mesh%lbn(k)
C
C  ....     CONSIDER BOUNDARY VORTICITY CONTRIBUTION :
C
            IF (j.NE.0) THEN                                  ! BOUNDARY

                dfi1(i)=dfi1(i)
     &                       -sdH(i,j)  *velocity(k,3)
     &                       -sdHty(i,j)*velocity(k,1)
     &                       +sdHtx(i,j)*velocity(k,2)
     &                       +sdDy(i,k) *vorticity(k,1)
c     &                       -sdDx(i,k) *vorticity(k,2) ! samo za test to brisi

                dfi2(i)=dfi2(i)
     &                       +sdH(i,j)  *velocity(k,2)
     &                       -sdHtz(i,j)*velocity(k,1)
     &                       +sdHtx(i,j)*velocity(k,3)
     &                       +sdDz(i,k) *vorticity(k,1)
c     &                       -sdDx(i,k) *vorticity(k,3) ! samo za test to brisi
C
C....       CONSIDER DOMAIN VORTICITY CONTRIBUTION :
C
            ELSE                                              ! DOMAIN
                dfi1(i)=dfi1(i)
     &                       +sdDy(i,k)*vorticity(k,1)-sdDx(i,k)*vorticity(k,2)

                dfi2(i)=dfi2(i)
     &                       +sdDz(i,k)*vorticity(k,1)-sdDx(i,k)*vorticity(k,3)

            END IF !boundary or domain node
          END DO !loop over columns
          print *,"dfi1",i,dfi1(i),dfi2(i)
        ELSE IF (mesh%bkmc(i).EQ.2) THEN ! TANGENTIAL BOUNDARY VELOCITIES ARE UNKNOWN
           print *,"kujku"
        ELSE ! THIS EQUATION IS NOT SOLVED
          print *,"eger"
          dfi1(i)=vorticity(mesh%gbn(i),2)
          dfi2(i)=vorticity(mesh%gbn(i),3)
        END IF
        dfi3(i)=vorticity(mesh%gbn(i),1) ! normal vorticity component is always known
      END DO !loop over rows


      END



C -----------------------------------------------------------------------------
      SUBROUTINE sdkm_solve8(eqn,env,inp,cpu,mesh,lnits,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,sdDxB,sdDyB,sdDzB,
     &                      sdA,sdPivot,velocity,vorticity)
C
C     $: Solve EQNth single domain kinematics equation
C        (linear interpolation fo domain vorticity)
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(inputtype) inp
      TYPE(CPUtype) cpu
      TYPE(penv) :: env

      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion
      REAL(8) sdA(env%zac:env%kon,mesh%nbnodes)  ! system matrix
      REAL(8) sdPivot(env%zac:env%kon)           ! system matrix

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni

      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes8) ! doktorat, enacba (4.7)
      REAL(8) sdDy(env%zac:env%kon,mesh%nnodes8) ! to je za drugi integral
      REAL(8) sdDz(env%zac:env%kon,mesh%nnodes8) ! na desni (linearno)

      REAL(8) sdDxB(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDyB(env%zac:env%kon,mesh%nbnodes) ! to je za drugi integral
      REAL(8) sdDzB(env%zac:env%kon,mesh%nbnodes) ! na desni (samo robni del, kvadratna)


      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3)
      INTEGER lslvt,lpret,lprep,lmaxit,lstopt,lnits
      REAL(4) lslveps

      REAL(8), ALLOCATABLE :: rhsv(:),xsol(:)

      INTEGER eqn,ierr
      REAL(4) cptime,rcpu,rrcpu

      ALLOCATE (rhsv(env%zac:env%kon),xsol(mesh%nbnodes))

      rcpu=cptime(0.)

      IF (env%nproc.EQ.1) THEN
        lslvt=0
        lpret=2
        lprep=2
        lmaxit=500
        lstopt=5
        lslveps=1.0E-15
        ierr=0
      ELSE
        lslvt=2
        lpret=1
        lprep=2
        lmaxit=inp%pkms_maxit
        lstopt=5
        lslveps=inp%pkms_eps
        ierr=0
      END IF

c
c     set up initial guess for solver
c
c      CALL sdkm_initguess(mesh,vorticity(:,eqn),xsol)
      xsol=0.0D0
C
C     Set up r.h.s.
C
      CALL sdkm_rhs8(eqn,env,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,sdDxB,sdDyB,sdDzB,
     &              velocity,vorticity,rhsv)
      cpu%time((eqn-1)*6+1)=cpu%time((eqn-1)*6+1)+cptime(rcpu)
      rcpu=cptime(0.)

c      print *,rhsv


c
c     solve
c
      IF (env%nproc.EQ.1) THEN
        CALL SolvEQNfm(lslvt,lpret,lprep,lmaxit,lstopt,lslveps,
     &              mesh%nbnodes,sdPivot,sdA,rhsv,xsol,lnits,rrcpu,ierr)
      ELSE
        CALL pSolvEQNfm(env,lslvt,lpret,lprep,lmaxit,lstopt,lslveps,
     &              mesh%nbnodes,env%nmn,sdPivot,sdA,rhsv,xsol,lnits,rrcpu,ierr)
      END IF
      cpu%time((eqn-1)*6+2)=cpu%time((eqn-1)*6+2)+cptime(rcpu)
c
c     copy solution to vorticity vector
c
      CALL sdkm_distunk(eqn,inp,mesh,vorticity(:,eqn),xsol)  ! tu bi morda podrelaksiral!!

      DEALLOCATE (rhsv,xsol)

      END


C -----------------------------------------------------------------------------
      SUBROUTINE sdkm_rhs8(eqn,env,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,sdDxB,sdDyB,sdDzB,
     &                    velocity,vorticity,dfi)
C
C     $: Form r.h.s. for EQNth single domain kinematics equation
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(penv) :: env

      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni

      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes8) ! doktorat, enacba (4.7)
      REAL(8) sdDy(env%zac:env%kon,mesh%nnodes8) ! to je za drugi integral
      REAL(8) sdDz(env%zac:env%kon,mesh%nnodes8) ! na desni

      REAL(8) sdDxB(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDyB(env%zac:env%kon,mesh%nbnodes) ! to je za drugi integral
      REAL(8) sdDzB(env%zac:env%kon,mesh%nbnodes) ! na desni


      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3)
      REAL(8) dfi(env%zac:env%kon)

      INTEGER eqn,i,j,k,v
C
C____ INITIALIZE ARRAY :
C
      dfi=0.00D0

      IF (eqn.EQ.1) THEN ! x eqaution
C
C____ LOOP OVER COLUMNS : for matrices with quadratic interpolation
C
      DO k=1,mesh%nnodes
        j=mesh%lbn(k)
        IF (j.NE.0) THEN                                  ! BOUNDARY
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdH(i,j)*(mesh%ny(i)*velocity(k,3)-mesh%nz(i)*velocity(k,2))
     &                       -sdHtx(i,j)*(mesh%ny(i)*velocity(k,2)+mesh%nz(i)*velocity(k,3))
     &                       +velocity(k,1)*(sdHty(i,j)*mesh%ny(i)+sdHtz(i,j)*mesh%nz(i))
     &                       +sdDxB(i,j)*(
     &                                  +mesh%nx(i)*vorticity(k,1)
     &                                  +mesh%ny(i)*vorticity(k,2)
     &                                  +mesh%nz(i)*vorticity(k,3))
          END DO !loop over rows
        END IF !boundary or domain node
      END DO !loop over columns
C
C____ LOOP OVER COLUMNS : for matrices with linear interpolation
C
      DO k=1,mesh%nnodes8
        j=mesh%lbn8(k)
        v=mesh%p827(k)
        IF (j.EQ.0) THEN                                   ! DOMAIN
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdDx(i,k)*(mesh%ny(i)*vorticity(v,2)+mesh%nz(i)*vorticity(v,3))
     &                       -vorticity(v,1)*(sdDy(i,k)*mesh%ny(i)+sdDz(i,k)*mesh%nz(i))
          END DO !loop over rows
        END IF !boundary or domain node
      END DO !loop over columns

      ELSE IF (eqn.EQ.2) THEN ! y eqaution
C
C____ LOOP OVER COLUMNS : : for matrices with quadratic interpolation
C
      DO k=1,mesh%nnodes
        j=mesh%lbn(k)
        IF (j.NE.0) THEN                                  ! BOUNDARY
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdH(i,j)*(mesh%nz(i)*velocity(k,1)-mesh%nx(i)*velocity(k,3))
     &                       -sdHty(i,j)*(mesh%nx(i)*velocity(k,1)+mesh%nz(i)*velocity(k,3))
     &                       +velocity(k,2)*(sdHtx(i,j)*mesh%nx(i)+sdHtz(i,j)*mesh%nz(i))
     &                       +sdDyB(i,j)*(
     &                                   mesh%nx(i)*vorticity(k,1)
     &                                  +mesh%ny(i)*vorticity(k,2)
     &                                  +mesh%nz(i)*vorticity(k,3))
          END DO !loop over rows
        END IF !boundary or domain node
      END DO !loop over columns
C
C____ LOOP OVER COLUMNS : for matrices with linear interpolation
C
      DO k=1,mesh%nnodes8
        j=mesh%lbn8(k)
        v=mesh%p827(k)
        IF (j.EQ.0) THEN                                  ! DOMAIN
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdDy(i,k)*(mesh%nx(i)*vorticity(v,1)+mesh%nz(i)*vorticity(v,3))
     &                       -vorticity(v,2)*(sdDx(i,k)*mesh%nx(i)+sdDz(i,k)*mesh%nz(i))
          END DO !loop over rows
        END IF !boundary or domain node
      END DO !loop over columns


      ELSE IF (eqn.EQ.3) THEN ! z eqaution
C
C____ LOOP OVER COLUMNS : for matrices with quadratic interpolation
C
      DO k=1,mesh%nnodes
        j=mesh%lbn(k)
        IF (j.NE.0) THEN                                  ! BOUNDARY
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdH(i,j)*(mesh%nx(i)*velocity(k,2)-mesh%ny(i)*velocity(k,1))
     &                       -sdHtz(i,j)*(mesh%nx(i)*velocity(k,1)+mesh%ny(i)*velocity(k,2))
     &                       +velocity(k,3)*(sdHtx(i,j)*mesh%nx(i)+sdHty(i,j)*mesh%ny(i))
     &                       +sdDzB(i,j)*(
     &                                   mesh%nx(i)*vorticity(k,1)
     &                                  +mesh%ny(i)*vorticity(k,2)
     &                                  +mesh%nz(i)*vorticity(k,3))
          END DO !loop over rows
        END IF !boundary or domain node
      END DO !loop over columns
C
C____ LOOP OVER COLUMNS : : for matrices with linear interpolation
C
      DO k=1,mesh%nnodes8
        j=mesh%lbn8(k)
        v=mesh%p827(k)
        IF (j.EQ.0) THEN                                  ! DOMAIN
          DO i=env%zac,env%kon
                dfi(i)=dfi(i)
     &                       +sdDz(i,k)*(mesh%nx(i)*vorticity(v,1)+mesh%ny(i)*vorticity(v,2))
     &                       -vorticity(v,3)*(sdDx(i,k)*mesh%nx(i)+sdDy(i,k)*mesh%ny(i))
          END DO !loop over rows
        END IF !boundary or domain node
      END DO !loop over columns


      END IF

C     resitev za bug; Ce ne resujem enacbe, dam na desno vrtincnost.
      DO i=env%zac,env%kon
        IF (mesh%kode(i,eqn).EQ.1.OR.mesh%kode(i,eqn).EQ.2) THEN   ! NE RESUJEM TE ENACBE
          dfi(i)=vorticity(mesh%gbn(i),eqn)
        END IF
      END DO


      END



C -----------------------------------------------------------------------------
      SUBROUTINE sdkm_initguess(mesh,vorticity,xsol)
C
C     $: Copy boundary vorticities to initial guess vector
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      REAL(8) xsol(mesh%nbnodes),vorticity(mesh%nnodes)
      INTEGER i,j
      
      DO i=1,mesh%nbnodes
        j=mesh%gbn(i)
        xsol(i)=vorticity(j)
      END DO
      
      END
C -----------------------------------------------------------------------------
      SUBROUTINE sdkm_distunk(eqn,inp,mesh,vorticity,xsol)
C
C     $: Copy solution vector to boundary vorticities
C
C -----------------------------------------------------------------------------
      USE inc_types 
      
      TYPE(inputtype) inp            
      TYPE(meshType) :: mesh
      REAL(8) xsol(mesh%nbnodes),vorticity(mesh%nnodes)
      INTEGER i,j,eqn
      
      DO i=1,mesh%nbnodes
        j=mesh%gbn(i)
        IF (mesh%kode(i,eqn).EQ.0.OR.mesh%kode(i,eqn).EQ.3) THEN
          vorticity(j)=vorticity(j)*(1.0D0-inp%urBw(eqn))+inp%urBw(eqn)*xsol(i)
        END IF
      END DO
      
      END

C -----------------------------------------------------------------------------
      SUBROUTINE sdkm_Rotdistunk(inp,mesh,vorticity,velocity,xsol1,xsol2,xsolNormal)
C
C     $: Copy solution vector to boundary vorticities
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(inputtype) inp
      TYPE(meshType) :: mesh
      REAL(8) xsol1(mesh%nbnodes),xsol2(mesh%nbnodes),xsolNormal(mesh%nbnodes)
      REAL(8) vorticity(mesh%nnodes,3),velocity(mesh%nnodes,3),v(3),vr(3)
      INTEGER i,j,eqn

      DO i=1,mesh%nbnodes
        j=mesh%gbn(i)
c        IF (mesh%bkmc(i).EQ.1) THEN
          v(1)=xsolNormal(i)
          v(2)=xsol1(i)
          v(3)=xsol2(i)
          vr=MATMUL(mesh%RotMatTransp(i)%v,v)
          DO eqn=1,3
            vorticity(j,eqn)=vorticity(j,eqn)*(1.0D0-inp%urBw(eqn))+inp%urBw(eqn)*vr(eqn)
          END DO
c        ELSE
c          Print *,"erere"
c        END IF
      END DO

      END



C -----------------------------------------------------------------------------
      SUBROUTINE sdKm_sysmBCN(eqn,env,mesh,sdA,sdDx,sdDy,sdDz)
C
C     $: Form single domain kinematics system matrix, based on BC normals
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(penv) :: env

      REAL(8) sdA(env%zac:env%kon,mesh%nbnodes)  ! system matrix

      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes)
      REAL(8) sdDy(env%zac:env%kon,mesh%nnodes)
      REAL(8) sdDz(env%zac:env%kon,mesh%nnodes)

      INTEGER i,j,dn,eqn

      DO j=1,mesh%nbnodes
        dn=mesh%gbn(j)
        DO i=env%zac,env%kon !1,mesh%nbnodes
          IF ( (mesh%kode(i,eqn).EQ.0.OR.mesh%kode(i,eqn).EQ.3).AND.(mesh%bcnl(i,1).NE.eqn) ) THEN   ! Tangential (kode=0 or 3)
            sdA(i,j)=mesh%nx(i)*sdDx(i,dn)
     &              +mesh%ny(i)*sdDy(i,dn)
     &              +mesh%nz(i)*sdDz(i,dn)
          ELSE                              ! do nothing (kode=1 or 2)
            IF (i.NE.j) THEN
              sdA(i,j)=0.0D0
            ELSE
              sdA(i,j)=1.0D0
            END IF
          END IF
        END DO
      END DO

      END




C -----------------------------------------------------------------------------
      SUBROUTINE sdKm_sysm(eqn,env,mesh,sdA,sdDx,sdDy,sdDz)
C
C     $: Form single domain kinematics system matrix
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      TYPE(penv) :: env
      
      REAL(8) sdA(env%zac:env%kon,mesh%nbnodes)  ! system matrix
      
      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes) 
      REAL(8) sdDy(env%zac:env%kon,mesh%nnodes) 
      REAL(8) sdDz(env%zac:env%kon,mesh%nnodes) 
      
      INTEGER i,j,dn,eqn

      DO j=1,mesh%nbnodes
        dn=mesh%gbn(j)
        DO i=env%zac,env%kon !1,mesh%nbnodes
          IF (mesh%kode(i,eqn).EQ.0.OR.mesh%kode(i,eqn).EQ.3) THEN   ! Tangential (kode=0 or 3)
            sdA(i,j)=mesh%nx(i)*sdDx(i,dn)
     &              +mesh%ny(i)*sdDy(i,dn)
     &              +mesh%nz(i)*sdDz(i,dn)     
          ELSE                              ! do nothing (kode=1 or 2)
            IF (i.NE.j) THEN
              sdA(i,j)=0.0D0
            ELSE
              sdA(i,j)=1.0D0
            END IF
          END IF
        END DO
      END DO
     
      END



C -----------------------------------------------------------------------------
      SUBROUTINE sdKm_RotSysm(env,mesh,sdA,sdDx,sdH)
C
C     $: Form single domain kinematics system matrix
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(penv) :: env

      REAL(8) sdA(env%zac:env%kon,mesh%nbnodes)  ! system matrix
      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes)
      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)

      INTEGER i,j,dn,eqn

      DO j=1,mesh%nbnodes
        dn=mesh%gbn(j)
        DO i=env%zac,env%kon !1,mesh%nbnodes
          IF (mesh%bkmc(i).EQ.1) THEN  ! tangential boundary vorticity is unknown
            sdA(i,j)=sdDx(i,dn)
          ELSE IF (mesh%bkmc(i).EQ.2) THEN  ! tangential boundary velocity is unknown
            sdA(i,j)=sdH(i,dn)
          ELSE                              ! do nothing
            IF (i.NE.j) THEN
              sdA(i,j)=0.0D0
            ELSE
              sdA(i,j)=1.0D0
            END IF
          END IF
        END DO
      END DO

      END


C -----------------------------------------------------------------------------
      SUBROUTINE sdKm_sysm8(eqn,env,mesh,sdA,sdDx,sdDy,sdDz)
C
C     $: Form single domain kinematics system matrix (sdD are boundary only)
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(penv) :: env

      REAL(8) sdA(env%zac:env%kon,mesh%nbnodes)  ! system matrix

      REAL(8) sdDx(env%zac:env%kon,mesh%nbnodes)
      REAL(8) sdDy(env%zac:env%kon,mesh%nbnodes)
      REAL(8) sdDz(env%zac:env%kon,mesh%nbnodes)

      INTEGER i,j,dn,eqn

      DO j=1,mesh%nbnodes
        dn=j !mesh%gbn(j)
        DO i=env%zac,env%kon !1,mesh%nbnodes
          IF (mesh%kode(i,eqn).EQ.0.OR.mesh%kode(i,eqn).EQ.3) THEN   ! Tangential (kode=0 or 3)
            sdA(i,j)=mesh%nx(i)*sdDx(i,dn)
     &              +mesh%ny(i)*sdDy(i,dn)
     &              +mesh%nz(i)*sdDz(i,dn)
          ELSE                              ! do nothing (kode=1 or 2)
            IF (i.NE.j) THEN
              sdA(i,j)=0.0D0
            ELSE
              sdA(i,j)=1.0D0
            END IF
          END IF
        END DO
      END DO

      END

C -----------------------------------------------------------------------------
      SUBROUTINE FMATsdkm(env,mesh,gauss,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &                ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)
C
C     $: Form Matrices,  kinematics eqaution (H, Htx, Hty, Htz, Dx, Dy, Dz)
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(penv) :: env
      
      INTEGER ic,je,dn,isrc,isip,i,j
      REAL(8) xp,yp,zp,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
      REAL(8) x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8

      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral 
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni
      
      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDy(env%zac:env%kon,mesh%nnodes) ! to je za drugi integral 
      REAL(8) sdDz(env%zac:env%kon,mesh%nnodes) ! na desni       
      REAL(8) ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz ! ocena natancnosti integralov
      REAL(8) ontHtxDz,ontHtyDx,ontHtzDy

      REAL(8) vsota,vsota1,vsota2,vsota3,hv,c
      REAL(8) htzkratx,htxkraty,htykratz

      REAL(8), ALLOCATABLE :: he(:),be(:),dxe(:),dye(:),dze(:)
      REAL(8), ALLOCATABLE :: htx(:),hty(:),htz(:)

      ALLOCATE (he(mesh%npob),be(mesh%npoc))
      ALLOCATE (dxe(mesh%npoc),dye(mesh%npoc),dze(mesh%npoc))
      ALLOCATE (htx(mesh%npob),hty(mesh%npob),htz(mesh%npob))

      sdH=0.0D0                    
      sdHtx=0.0D0
      sdHty=0.0D0
      sdHtz=0.0D0                       
      sdDx=0.0D0
      sdDy=0.0D0
      sdDz=0.0D0     
      
      ontH=0.0D0
      ontHtx=0.0D0
      ontHty=0.0D0
      ontHtz=0.0D0
      ontDx=0.0D0
      ontDy=0.0D0
      ontDz=0.0D0
      ontHtxDz=0.0D0
      ontHtyDx=0.0D0
      ontHtzDy=0.0D0
C
C     Zanka po vseh robnih vozliscih = izvornih tockah
C    
      DO i=env%zac,env%kon !1,mesh%nbnodes 
        dn=mesh%gbn(i)
        xp=mesh%x(dn,1)  ! izvorna tocka x
        yp=mesh%x(dn,2)
        zp=mesh%x(dn,3)  
C 
C       zanka po robnih elementih  - INTEGRACIJA PO ROBU
C
        DO je=1,mesh%nbelem
          x1=mesh%x(mesh%ibc(je,1),1)
          y1=mesh%x(mesh%ibc(je,1),2)
          z1=mesh%x(mesh%ibc(je,1),3)
          x2=mesh%x(mesh%ibc(je,3),1)
          y2=mesh%x(mesh%ibc(je,3),2)
          z2=mesh%x(mesh%ibc(je,3),3)
          x3=mesh%x(mesh%ibc(je,5),1)
          y3=mesh%x(mesh%ibc(je,5),2)
          z3=mesh%x(mesh%ibc(je,5),3)
          x4=mesh%x(mesh%ibc(je,7),1)
          y4=mesh%x(mesh%ibc(je,7),2)
          z4=mesh%x(mesh%ibc(je,7),3)                   
C         SET SOURCE POINT :
          isrc=0
          DO isip=1,mesh%npob
            IF (dn .EQ. mesh%ibc(je,isip)) isrc=isip
          END DO        

C         integracija po robnih elementih celice            
          CALL INTEBc9sdkm(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                     isrc,mesh%npob,he,htx,hty,htz)            
c         zlozimo po matrikah            
          DO isip=1,mesh%npob
            j=mesh%lbn(mesh%ibc(je,isip)) ! stolpec
c          if (isrc.ne.0.and.j.eq.i) print *,"w",htx(isip),isip
            sdH(i,j)=sdH(i,j)+he(isip)             
            sdHtx(i,j)=sdHtx(i,j)-Htx(isip)
            sdHty(i,j)=sdHty(i,j)-Hty(isip)
            sdHtz(i,j)=sdHtz(i,j)-Htz(isip)                                          
          END DO      
        END DO ! po robnih elementih
        
c       izracunam prosti koeficient c in ga dam na diagonalo
c       predpostavim togi premik telesa u=1, q=0          
        c=0.0D0
        DO j=1,mesh%nbnodes        
          c=c-sdH(i,j)         
        END DO
        ontH=max(ontH,c)
        sdH(i,i)=sdH(i,i)+c
c       preverim se ht matrike, vsota vseh clenov mora biti 0
c test s couettom ali z rigid body movment je isti - da isti popravek singularnega integrala
        vsota=0.0D0
        vsota1=0.0D0
        DO j=1,mesh%nbnodes         
          vsota=vsota+sdHtx(i,j)         
          vsota1=vsota1+sdHtx(i,j)*mesh%x(mesh%gbn(j),1)
        END DO
c        print *,dn,vsota,vsota1/mesh%x(mesh%gbn(i),1)
        ontHtx=max(ontHtx,ABS(vsota))
        sdhtx(i,i)=sdhtx(i,i)-vsota  ! popravimo singularni integral
        vsota=0.0D0
        vsota1=0.0D0
        DO j=1,mesh%nbnodes         
          vsota=vsota+sdHty(i,j)         
          vsota1=vsota1+sdHty(i,j)*mesh%x(mesh%gbn(j),2)
        END DO
        ontHty=max(ontHty,ABS(vsota))
        sdhty(i,i)=sdhty(i,i)-vsota  ! popravimo singularni integral
        vsota=0.0D0
        vsota1=0.0D0
        DO j=1,mesh%nbnodes         
           vsota=vsota+sdHtz(i,j)         
           vsota1=vsota1+sdHtz(i,j)*mesh%x(mesh%gbn(j),3)
        END DO
        ontHtz=max(ontHtz,ABS(vsota))          
        sdhtz(i,i)=sdhtz(i,i)-vsota  ! popravimo singularni integral        
C 
C       zanka po notranjih celicah  - INTEGRACIJA PO OBMOCJU
C
        DO ic=1,mesh%nicell
          x1=mesh%x(mesh%idc(ic,6),1)
          y1=mesh%x(mesh%idc(ic,6),2)
          z1=mesh%x(mesh%idc(ic,6),3)          
          x2=mesh%x(mesh%idc(ic,7),1)
          y2=mesh%x(mesh%idc(ic,7),2)
          z2=mesh%x(mesh%idc(ic,7),3)          
          x3=mesh%x(mesh%idc(ic,2),1)
          y3=mesh%x(mesh%idc(ic,2),2)
          z3=mesh%x(mesh%idc(ic,2),3)          
          x4=mesh%x(mesh%idc(ic,3),1)
          y4=mesh%x(mesh%idc(ic,3),2)
          z4=mesh%x(mesh%idc(ic,3),3)          
          x5=mesh%x(mesh%idc(ic,5),1)
          y5=mesh%x(mesh%idc(ic,5),2)
          z5=mesh%x(mesh%idc(ic,5),3)          
          x6=mesh%x(mesh%idc(ic,8),1)
          y6=mesh%x(mesh%idc(ic,8),2)
          z6=mesh%x(mesh%idc(ic,8),3)          
          x7=mesh%x(mesh%idc(ic,1),1)
          y7=mesh%x(mesh%idc(ic,1),2)
          z7=mesh%x(mesh%idc(ic,1),3)          
          x8=mesh%x(mesh%idc(ic,4),1)
          y8=mesh%x(mesh%idc(ic,4),2)
          z8=mesh%x(mesh%idc(ic,4),3)          
C         SET SOURCE POINT :                 
          isrc=0
          DO isip=1,mesh%npoc
            IF (dn .EQ. mesh%idc(ic,isip)) isrc=isip
          END DO
          CALL INTEDc27dcLap(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                  x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,be,dxe,dye,dze,isrc,mesh%npoc)
c         zlozimo vrstico po matrikah      
          DO isip=1,mesh%npoc  
            j=mesh%idc(ic,isip)                  
            sdDx(i,j)=sdDx(i,j)+dxe(isip)
            sdDy(i,j)=sdDy(i,j)+dye(isip)
            sdDz(i,j)=sdDz(i,j)+dze(isip)
          END DO 
        END DO ! po notranjih celicah
                  
c       preverimo natancnost Dx,Dy,Dz integralov
c       v=(z,0,0), w=(0,1,0)
        hv=0.0D0
        vsota=0.0D0
        DO isip=1,mesh%nbnodes
          j=mesh%gbn(isip)
          hv=hv+sdH(i,isip)*mesh%x(j,3)                  
        END DO
        DO j=1,mesh%nnodes
          vsota=vsota+sdDz(i,j)
        END DO
        ontDz=max(ontDz,ABS(hv-vsota))
        sdDz(i,mesh%gbn(i))=sdDz(i,mesh%gbn(i))+hv-vsota  ! popravimo singularni integral
c       v=(0,x,0), w=(0,0,1)
        hv=0.0D0
        vsota=0.0D0
        DO isip=1,mesh%nbnodes
          j=mesh%gbn(isip)
          hv=hv+sdH(i,isip)*mesh%x(j,1)                  
        END DO
        DO j=1,mesh%nnodes
          vsota=vsota+sdDx(i,j)
        END DO
        ontDx=max(ontDx,ABS(hv-vsota))
        sdDx(i,mesh%gbn(i))=sdDx(i,mesh%gbn(i))+hv-vsota  ! popravimo singularni integral
c       v=(0,0,y), w=(1,0,0)
        hv=0.0D0
        vsota=0.0D0
        DO isip=1,mesh%nbnodes
          j=mesh%gbn(isip)
          hv=hv+sdH(i,isip)*mesh%x(j,2)                  
        END DO
        DO j=1,mesh%nnodes
          vsota=vsota+sdDy(i,j)
        END DO
        ontDy=max(ontDy,ABS(hv-vsota))                  
        sdDy(i,mesh%gbn(i))=sdDy(i,mesh%gbn(i))+hv-vsota  ! popravimo singularni integral

c       singularni popravljeni - se enkrat preverimo
        htxkraty=0.0D0
        htykratz=0.0D0      
        htzkratx=0.0D0          
        vsota1=0.0D0
        vsota2=0.0D0
        vsota3=0.0D0                
        DO j=1,mesh%nnodes         
          vsota1=vsota1+sdDz(i,j)
          vsota2=vsota2+sdDx(i,j)      
          vsota3=vsota3+sdDy(i,j)           
        END DO
        DO j=1,mesh%nbnodes         
          htxkraty=htxkraty+sdHtx(i,j)*mesh%x(mesh%gbn(j),2)
          htykratz=htykratz+sdHty(i,j)*mesh%x(mesh%gbn(j),3)
          htzkratx=htzkratx+sdHtz(i,j)*mesh%x(mesh%gbn(j),1)
        END DO
        ontHtzDy=max(ontHtzDy,ABS(htzkratx+vsota3))
        ontHtyDx=max(ontHtyDx,ABS(htykratz+vsota2))
        ontHtxDz=max(ontHtxDz,ABS(htxkraty+vsota1))
        
      END DO ! po robnih vozliscih = izvornih tockah
      ontH=ABS(ontH-0.5D0)  ! ker vem, da je najvecji C = 1/2

      END


C -----------------------------------------------------------------------------
      SUBROUTINE FMATsdkmRot(env,mesh,gauss,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &                ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)
C
C     $: Form Matrices,  kinematics eqaution (H, Htx, Hty, Htz, Dx, Dy, Dz)
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(penv) :: env

      INTEGER ic,je,dn,isrc,isip,i,j
      REAL(8) xp,yp,zp,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
      REAL(8) x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8

      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni

      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDy(env%zac:env%kon,mesh%nnodes) ! to je za drugi integral
      REAL(8) sdDz(env%zac:env%kon,mesh%nnodes) ! na desni
      REAL(8) ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz ! ocena natancnosti integralov
      REAL(8) ontHtxDz,ontHtyDx,ontHtzDy

      REAL(8) vsota,vsota1,vsota2,vsota3,hv,c
      REAL(8) htzkratx,htxkraty,htykratz
      REAL(8) tmp(3),t2(3),dp,omega(3),velo(3)

      REAL(8), ALLOCATABLE :: he(:),be(:),dxe(:),dye(:),dze(:)
      REAL(8), ALLOCATABLE :: htx(:),hty(:),htz(:)
      REAL(8), ALLOCATABLE :: orix(:,:)

      ALLOCATE (he(mesh%npob),be(mesh%npoc))
      ALLOCATE (dxe(mesh%npoc),dye(mesh%npoc),dze(mesh%npoc))
      ALLOCATE (htx(mesh%npob),hty(mesh%npob),htz(mesh%npob))

      ALLOCATE (orix(mesh%nnodes,mesh%npx))

      sdH=0.0D0
      sdHtx=0.0D0
      sdHty=0.0D0
      sdHtz=0.0D0
      sdDx=0.0D0
      sdDy=0.0D0
      sdDz=0.0D0

      ontH=0.0D0
      ontHtx=0.0D0
      ontHty=0.0D0
      ontHtz=0.0D0
      ontDx=0.0D0
      ontDy=0.0D0
      ontDz=0.0D0
      ontHtxDz=0.0D0
      ontHtyDx=0.0D0
      ontHtzDy=0.0D0
C
C     Zanka po vseh robnih vozliscih = izvornih tockah
C
      DO i=env%zac,env%kon !1,mesh%nbnodes
C
C       Rotate mesh node locations by rotation matrix R(i),
C       so that for each source point, we have
C       v=(v_normal,vtang1,v_tang2) and
C       omega = (omega_normal,omega_tang1, omega_tang2)
C
        orix=mesh%x
        CALL RotateByR(mesh%x,mesh%RotMat(i)%v,mesh%nnodes)
C
C       Set source point
C
        dn=mesh%gbn(i)
        xp=mesh%x(dn,1)  ! izvorna tocka x
        yp=mesh%x(dn,2)
        zp=mesh%x(dn,3)
C
C       zanka po robnih elementih  - INTEGRACIJA PO ROBU
C
        DO je=1,mesh%nbelem
          x1=mesh%x(mesh%ibc(je,1),1)
          y1=mesh%x(mesh%ibc(je,1),2)
          z1=mesh%x(mesh%ibc(je,1),3)
          x2=mesh%x(mesh%ibc(je,3),1)
          y2=mesh%x(mesh%ibc(je,3),2)
          z2=mesh%x(mesh%ibc(je,3),3)
          x3=mesh%x(mesh%ibc(je,5),1)
          y3=mesh%x(mesh%ibc(je,5),2)
          z3=mesh%x(mesh%ibc(je,5),3)
          x4=mesh%x(mesh%ibc(je,7),1)
          y4=mesh%x(mesh%ibc(je,7),2)
          z4=mesh%x(mesh%ibc(je,7),3)
C         SET SOURCE POINT :
          isrc=0
          DO isip=1,mesh%npob
            IF (dn .EQ. mesh%ibc(je,isip)) isrc=isip
          END DO

C         integracija po robnih elementih celice
          CALL INTEBc9sdkm(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                     isrc,mesh%npob,he,htx,hty,htz)

c         zlozimo po matrikah
          DO isip=1,mesh%npob
            j=mesh%lbn(mesh%ibc(je,isip)) ! stolpec
c          if (isrc.ne.0.and.j.eq.i) print *,"w",htx(isip),isip
            sdH(i,j)=sdH(i,j)+he(isip)
            sdHtx(i,j)=sdHtx(i,j)-Htx(isip)
            sdHty(i,j)=sdHty(i,j)-Hty(isip)
            sdHtz(i,j)=sdHtz(i,j)-Htz(isip)
          END DO
        END DO ! po robnih elementih

c       izracunam prosti koeficient c in ga dam na diagonalo
c       predpostavim togi premik telesa u=1, q=0
        c=0.0D0
        DO j=1,mesh%nbnodes
          c=c-sdH(i,j)
        END DO
        ontH=max(ontH,c)
        sdH(i,i)=sdH(i,i)+c
c       preverim se ht matrike, vsota vseh clenov mora biti 0
c test s couettom ali z rigid body movment je isti - da isti popravek singularnega integrala
        vsota=0.0D0
        vsota1=0.0D0
        DO j=1,mesh%nbnodes
          vsota=vsota+sdHtx(i,j)
          vsota1=vsota1+sdHtx(i,j)*mesh%x(mesh%gbn(j),1)
        END DO
c        print *,dn,vsota,vsota1/mesh%x(mesh%gbn(i),1)
        ontHtx=max(ontHtx,ABS(vsota))
        sdhtx(i,i)=sdhtx(i,i)-vsota  ! popravimo singularni integral
        vsota=0.0D0
        vsota1=0.0D0
        DO j=1,mesh%nbnodes
          vsota=vsota+sdHty(i,j)
          vsota1=vsota1+sdHty(i,j)*mesh%x(mesh%gbn(j),2)
        END DO
        ontHty=max(ontHty,ABS(vsota))
        sdhty(i,i)=sdhty(i,i)-vsota  ! popravimo singularni integral
        vsota=0.0D0
        vsota1=0.0D0
        DO j=1,mesh%nbnodes
           vsota=vsota+sdHtz(i,j)
           vsota1=vsota1+sdHtz(i,j)*mesh%x(mesh%gbn(j),3)
        END DO
        ontHtz=max(ontHtz,ABS(vsota))
        sdhtz(i,i)=sdhtz(i,i)-vsota  ! popravimo singularni integral
C
C       zanka po notranjih celicah  - INTEGRACIJA PO OBMOCJU
C
        DO ic=1,mesh%nicell
          x1=mesh%x(mesh%idc(ic,6),1)
          y1=mesh%x(mesh%idc(ic,6),2)
          z1=mesh%x(mesh%idc(ic,6),3)
          x2=mesh%x(mesh%idc(ic,7),1)
          y2=mesh%x(mesh%idc(ic,7),2)
          z2=mesh%x(mesh%idc(ic,7),3)
          x3=mesh%x(mesh%idc(ic,2),1)
          y3=mesh%x(mesh%idc(ic,2),2)
          z3=mesh%x(mesh%idc(ic,2),3)
          x4=mesh%x(mesh%idc(ic,3),1)
          y4=mesh%x(mesh%idc(ic,3),2)
          z4=mesh%x(mesh%idc(ic,3),3)
          x5=mesh%x(mesh%idc(ic,5),1)
          y5=mesh%x(mesh%idc(ic,5),2)
          z5=mesh%x(mesh%idc(ic,5),3)
          x6=mesh%x(mesh%idc(ic,8),1)
          y6=mesh%x(mesh%idc(ic,8),2)
          z6=mesh%x(mesh%idc(ic,8),3)
          x7=mesh%x(mesh%idc(ic,1),1)
          y7=mesh%x(mesh%idc(ic,1),2)
          z7=mesh%x(mesh%idc(ic,1),3)
          x8=mesh%x(mesh%idc(ic,4),1)
          y8=mesh%x(mesh%idc(ic,4),2)
          z8=mesh%x(mesh%idc(ic,4),3)
C         SET SOURCE POINT :
          isrc=0
          DO isip=1,mesh%npoc
            IF (dn .EQ. mesh%idc(ic,isip)) isrc=isip
          END DO
          CALL INTEDc27dcLap(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                  x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,be,dxe,dye,dze,isrc,mesh%npoc)

c         zlozimo vrstico po matrikah
          DO isip=1,mesh%npoc
            j=mesh%idc(ic,isip)
            sdDx(i,j)=sdDx(i,j)+dxe(isip)
            sdDy(i,j)=sdDy(i,j)+dye(isip)
            sdDz(i,j)=sdDz(i,j)+dze(isip)
          END DO
        END DO ! po notranjih celicah

c       preverimo natancnost Dx,Dy,Dz integralov
c       v=(z,0,0), w=(0,1,0)
        hv=0.0D0
        vsota=0.0D0
        DO isip=1,mesh%nbnodes
          j=mesh%gbn(isip)
          hv=hv+sdH(i,isip)*mesh%x(j,3)
        END DO
        DO j=1,mesh%nnodes
          vsota=vsota+sdDz(i,j)
        END DO
        ontDz=max(ontDz,ABS(hv-vsota))
        sdDz(i,mesh%gbn(i))=sdDz(i,mesh%gbn(i))+hv-vsota  ! popravimo singularni integral
c       v=(0,x,0), w=(0,0,1)
        hv=0.0D0
        vsota=0.0D0
        DO isip=1,mesh%nbnodes
          j=mesh%gbn(isip)
          hv=hv+sdH(i,isip)*mesh%x(j,1)
        END DO
        DO j=1,mesh%nnodes
          vsota=vsota+sdDx(i,j)
        END DO
        ontDx=max(ontDx,ABS(hv-vsota))
        sdDx(i,mesh%gbn(i))=sdDx(i,mesh%gbn(i))+hv-vsota  ! popravimo singularni integral
c       v=(0,0,y), w=(1,0,0)
        hv=0.0D0
        vsota=0.0D0
        DO isip=1,mesh%nbnodes
          j=mesh%gbn(isip)
          hv=hv+sdH(i,isip)*mesh%x(j,2)
        END DO
        DO j=1,mesh%nnodes
          vsota=vsota+sdDy(i,j)
        END DO
        ontDy=max(ontDy,ABS(hv-vsota))
        sdDy(i,mesh%gbn(i))=sdDy(i,mesh%gbn(i))+hv-vsota  ! popravimo singularni integral

c       singularni popravljeni - se enkrat preverimo
        htxkraty=0.0D0
        htykratz=0.0D0
        htzkratx=0.0D0
        vsota1=0.0D0
        vsota2=0.0D0
        vsota3=0.0D0
        DO j=1,mesh%nnodes
          vsota1=vsota1+sdDz(i,j)
          vsota2=vsota2+sdDx(i,j)
          vsota3=vsota3+sdDy(i,j)
        END DO
        DO j=1,mesh%nbnodes
          htxkraty=htxkraty+sdHtx(i,j)*mesh%x(mesh%gbn(j),2)
          htykratz=htykratz+sdHty(i,j)*mesh%x(mesh%gbn(j),3)
          htzkratx=htzkratx+sdHtz(i,j)*mesh%x(mesh%gbn(j),1)
        END DO
        ontHtzDy=max(ontHtzDy,ABS(htzkratx+vsota3))
        ontHtyDx=max(ontHtyDx,ABS(htykratz+vsota2))
        ontHtxDz=max(ontHtxDz,ABS(htxkraty+vsota1))

c       Rotate mesh back
        mesh%x=orix

      END DO ! po robnih vozliscih = izvornih tockah
      ontH=ABS(ontH-0.5D0)  ! ker vem, da je najvecji C = 1/2

      END



C -----------------------------------------------------------------------------
      SUBROUTINE FMATsdkm8(env,mesh,gauss,sdH,sdHtx,sdHty,sdHtz,
     &                sdDx8,sdDy8,sdDz8,  ! domain part linear
     &                sdDxB,sdDyB,sdDzB,  ! boundary part quadratic
     &                ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)
C
C     $: Form Matrices,  kinematics eqaution (H, Htx, Hty, Htz, Dx, Dy, Dz)
c     v Dx,Dy in Dz je interpolacija linerana
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(penv) :: env

      INTEGER ic,je,dn,isrc,isip,i,j
      REAL(8) xp,yp,zp,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
      REAL(8) x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8

      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni

      REAL(8) sdDx8(env%zac:env%kon,mesh%nnodes8) ! doktorat, enacba (4.7)
      REAL(8) sdDy8(env%zac:env%kon,mesh%nnodes8) ! to je za drugi integral
      REAL(8) sdDz8(env%zac:env%kon,mesh%nnodes8) ! na desni

      REAL(8) sdDxB(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDyB(env%zac:env%kon,mesh%nbnodes) ! to je za drugi integral
      REAL(8) sdDzB(env%zac:env%kon,mesh%nbnodes) ! na desni SAMO ROBNE

      REAL(8) ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz ! ocena natancnosti integralov
      REAL(8) ontHtxDz,ontHtyDx,ontHtzDy

      REAL(8) vsota,vsota1,vsota2,vsota3,hv,c
      REAL(8) vsotaX,vsotaY,vsotaZ
      REAL(8) htzkratx,htxkraty,htykratz

      REAL(8), ALLOCATABLE :: he(:)
      REAL(8), ALLOCATABLE :: be8(:),dxe8(:),dye8(:),dze8(:)
      REAL(8), ALLOCATABLE :: beB(:),dxeB(:),dyeB(:),dzeB(:)
      REAL(8), ALLOCATABLE :: htx(:),hty(:),htz(:)

      ALLOCATE (he(mesh%npob))
      ALLOCATE (htx(mesh%npob),hty(mesh%npob),htz(mesh%npob))
      ALLOCATE (be8(mesh%npoc8),dxe8(mesh%npoc8),dye8(mesh%npoc8),dze8(mesh%npoc8))
      ALLOCATE (beB(mesh%npoc),dxeB(mesh%npoc),dyeB(mesh%npoc),dzeB(mesh%npoc))

      sdH=0.0D0
      sdHtx=0.0D0
      sdHty=0.0D0
      sdHtz=0.0D0

      sdDx8=0.0D0
      sdDy8=0.0D0
      sdDz8=0.0D0

      sdDxB=0.0D0
      sdDyB=0.0D0
      sdDzB=0.0D0


      ontH=0.0D0
      ontHtx=0.0D0
      ontHty=0.0D0
      ontHtz=0.0D0
      ontDx=0.0D0
      ontDy=0.0D0
      ontDz=0.0D0
      ontHtxDz=0.0D0
      ontHtyDx=0.0D0
      ontHtzDy=0.0D0
C
C     Zanka po vseh robnih vozliscih = izvornih tockah
C
      DO i=env%zac,env%kon !1,mesh%nbnodes
        dn=mesh%gbn(i)
        xp=mesh%x(dn,1)  ! izvorna tocka x
        yp=mesh%x(dn,2)
        zp=mesh%x(dn,3)
C
C       zanka po robnih elementih  - INTEGRACIJA PO ROBU
C
        DO je=1,mesh%nbelem
          x1=mesh%x(mesh%ibc(je,1),1)
          y1=mesh%x(mesh%ibc(je,1),2)
          z1=mesh%x(mesh%ibc(je,1),3)
          x2=mesh%x(mesh%ibc(je,3),1)
          y2=mesh%x(mesh%ibc(je,3),2)
          z2=mesh%x(mesh%ibc(je,3),3)
          x3=mesh%x(mesh%ibc(je,5),1)
          y3=mesh%x(mesh%ibc(je,5),2)
          z3=mesh%x(mesh%ibc(je,5),3)
          x4=mesh%x(mesh%ibc(je,7),1)
          y4=mesh%x(mesh%ibc(je,7),2)
          z4=mesh%x(mesh%ibc(je,7),3)
C         SET SOURCE POINT :
          isrc=0
          DO isip=1,mesh%npob
            IF (dn .EQ. mesh%ibc(je,isip)) isrc=isip
          END DO

C         integracija po robnih elementih celice
          CALL INTEBc9sdkm(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                     isrc,mesh%npob,he,htx,hty,htz)
c         zlozimo po matrikah
          DO isip=1,mesh%npob
            j=mesh%lbn(mesh%ibc(je,isip)) ! stolpec
c          if (isrc.ne.0.and.j.eq.i) print *,"w",htx(isip),isip
            sdH(i,j)=sdH(i,j)+he(isip)
            sdHtx(i,j)=sdHtx(i,j)-Htx(isip)
            sdHty(i,j)=sdHty(i,j)-Hty(isip)
            sdHtz(i,j)=sdHtz(i,j)-Htz(isip)
          END DO
        END DO ! po robnih elementih

c       izracunam prosti koeficient c in ga dam na diagonalo
c       predpostavim togi premik telesa u=1, q=0
        c=0.0D0
        DO j=1,mesh%nbnodes
          c=c-sdH(i,j)
        END DO
        ontH=max(ontH,c)
        sdH(i,i)=sdH(i,i)+c
c       preverim se ht matrike, vsota vseh clenov mora biti 0
c test s couettom ali z rigid body movment je isti - da isti popravek singularnega integrala
        vsota=0.0D0
        vsota1=0.0D0
        DO j=1,mesh%nbnodes
          vsota=vsota+sdHtx(i,j)
          vsota1=vsota1+sdHtx(i,j)*mesh%x(mesh%gbn(j),1)
        END DO
c        print *,dn,vsota,vsota1/mesh%x(mesh%gbn(i),1)
        ontHtx=max(ontHtx,ABS(vsota))
        sdhtx(i,i)=sdhtx(i,i)-vsota  ! popravimo singularni integral
        vsota=0.0D0
        vsota1=0.0D0
        DO j=1,mesh%nbnodes
          vsota=vsota+sdHty(i,j)
          vsota1=vsota1+sdHty(i,j)*mesh%x(mesh%gbn(j),2)
        END DO
        ontHty=max(ontHty,ABS(vsota))
        sdhty(i,i)=sdhty(i,i)-vsota  ! popravimo singularni integral
        vsota=0.0D0
        vsota1=0.0D0
        DO j=1,mesh%nbnodes
           vsota=vsota+sdHtz(i,j)
           vsota1=vsota1+sdHtz(i,j)*mesh%x(mesh%gbn(j),3)
        END DO
        ontHtz=max(ontHtz,ABS(vsota))
        sdhtz(i,i)=sdhtz(i,i)-vsota  ! popravimo singularni integral
C
C       zanka po notranjih celicah  - INTEGRACIJA PO OBMOCJU
C
        vsotaX=0.0D0
        vsotaY=0.0D0
        vsotaZ=0.0D0
        vsota1=0.0D0
        vsota2=0.0D0
        vsota3=0.0D0
        DO ic=1,mesh%nicell
          x1=mesh%x(mesh%idc(ic,6),1)
          y1=mesh%x(mesh%idc(ic,6),2)
          z1=mesh%x(mesh%idc(ic,6),3)
          x2=mesh%x(mesh%idc(ic,7),1)
          y2=mesh%x(mesh%idc(ic,7),2)
          z2=mesh%x(mesh%idc(ic,7),3)
          x3=mesh%x(mesh%idc(ic,2),1)
          y3=mesh%x(mesh%idc(ic,2),2)
          z3=mesh%x(mesh%idc(ic,2),3)
          x4=mesh%x(mesh%idc(ic,3),1)
          y4=mesh%x(mesh%idc(ic,3),2)
          z4=mesh%x(mesh%idc(ic,3),3)
          x5=mesh%x(mesh%idc(ic,5),1)
          y5=mesh%x(mesh%idc(ic,5),2)
          z5=mesh%x(mesh%idc(ic,5),3)
          x6=mesh%x(mesh%idc(ic,8),1)
          y6=mesh%x(mesh%idc(ic,8),2)
          z6=mesh%x(mesh%idc(ic,8),3)
          x7=mesh%x(mesh%idc(ic,1),1)
          y7=mesh%x(mesh%idc(ic,1),2)
          z7=mesh%x(mesh%idc(ic,1),3)
          x8=mesh%x(mesh%idc(ic,4),1)
          y8=mesh%x(mesh%idc(ic,4),2)
          z8=mesh%x(mesh%idc(ic,4),3)
C         SET SOURCE POINT :
          isrc=0
          DO isip=1,mesh%npoc
            IF (dn .EQ. mesh%idc(ic,isip)) isrc=isip
          END DO
c         quadratic part  (integriram VSE celice, notranjega ne shranim, rabim samo za izracun singularnega integrala)
c          IF (mesh%bcells(ic).EQ.1) THEN ! ta celica ima robna vozlicsa
            CALL INTEDc27dcLap(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                  x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,beB,dxeB,dyeB,dzeB,isrc,mesh%npoc)
c           zlozimo vrstico po matrikah
            DO isip=1,mesh%npoc
              j=mesh%idc(ic,isip)
              vsotaX=vsotaX+dxeB(isip)
              vsotaY=vsotaY+dyeB(isip)
              vsotaZ=vsotaZ+dzeB(isip)
c             je na robu?
              IF (mesh%lbn(j).NE.0) THEN
                sdDxB(i,mesh%lbn(j))=sdDxB(i,mesh%lbn(j))+dxeB(isip)
                sdDyB(i,mesh%lbn(j))=sdDyB(i,mesh%lbn(j))+dyeB(isip)
                sdDzB(i,mesh%lbn(j))=sdDzB(i,mesh%lbn(j))+dzeB(isip)
              END IF
            END DO
c          END IF
c         linear part  TU BI LAHKO ROBNIH NE SHRANJEVAL !!!
          CALL INTEDc8dcLap(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                  x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,be8,dxe8,dye8,dze8,isrc,mesh%npoc8)
c         zlozimo vrstico po matrikah
          DO isip=1,mesh%npoc8
            j=mesh%p278(mesh%idc8(ic,isip))  !!!!!!!!!!!! TO POMENI DA RABIM NOV IDC !!!!!!!!!!!!!
            vsota1=vsota1+dze8(isip)
            vsota2=vsota2+dxe8(isip)
            vsota3=vsota3+dye8(isip)
            IF (j.NE.0) THEN ! da ne dobim robnih
              sdDx8(i,j)=sdDx8(i,j)+dxe8(isip)
              sdDy8(i,j)=sdDy8(i,j)+dye8(isip)
              sdDz8(i,j)=sdDz8(i,j)+dze8(isip)
            END IF
          END DO
        END DO ! po notranjih celicah

c       QUADRATIC PART
c       preverimo natancnost Dx,Dy,Dz integralov
c       v=(z,0,0), w=(0,1,0)
        hv=0.0D0
        vsota=0.0D0
        DO isip=1,mesh%nbnodes
          j=mesh%gbn(isip)
          hv=hv+sdH(i,isip)*mesh%x(j,3)
        END DO
        DO j=1,mesh%nnodes8
          vsota=vsota+sdDz8(i,j)  ! ta vsota je identicno enaka ce jo izracunam s kvadratno ali s linerano interpolacijo
        END DO
        ontDz=max(ontDz,ABS(hv-vsotaZ))
        sdDzB(i,i)=sdDzB(i,i)+hv-vsotaZ  ! popravimo singularni integral
c       v=(0,x,0), w=(0,0,1)
        hv=0.0D0
        vsota=0.0D0
        DO isip=1,mesh%nbnodes
          j=mesh%gbn(isip)
          hv=hv+sdH(i,isip)*mesh%x(j,1)
        END DO
        DO j=1,mesh%nnodes8
          vsota=vsota+sdDx8(i,j)
        END DO
        ontDx=max(ontDx,ABS(hv-vsotaX))
        sdDxB(i,i)=sdDxB(i,i)+hv-vsotaX  ! popravimo singularni integral
c       v=(0,0,y), w=(1,0,0)
        hv=0.0D0
        vsota=0.0D0
        DO isip=1,mesh%nbnodes
          j=mesh%gbn(isip)
          hv=hv+sdH(i,isip)*mesh%x(j,2)
        END DO
        DO j=1,mesh%nnodes8
          vsota=vsota+sdDy8(i,j)
        END DO
        ontDy=max(ontDy,ABS(hv-vsotaY))
        sdDyB(i,i)=sdDyB(i,i)+hv-vsotaY  ! popravimo singularni integral

c       singularni popravljeni - se enkrat preverimo
        htxkraty=0.0D0
        htykratz=0.0D0
        htzkratx=0.0D0
        DO j=1,mesh%nbnodes
          htxkraty=htxkraty+sdHtx(i,j)*mesh%x(mesh%gbn(j),2)
          htykratz=htykratz+sdHty(i,j)*mesh%x(mesh%gbn(j),3)
          htzkratx=htzkratx+sdHtz(i,j)*mesh%x(mesh%gbn(j),1)
        END DO
        ontHtzDy=max(ontHtzDy,ABS(htzkratx+vsota3))
        ontHtyDx=max(ontHtyDx,ABS(htykratz+vsota2))
        ontHtxDz=max(ontHtxDz,ABS(htxkraty+vsota1))

      END DO ! po robnih vozliscih = izvornih tockah
      ontH=ABS(ontH-0.5D0)  ! ker vem, da je najvecji C = 1/2

      END

C -----------------------------------------------------------------------------            
      SUBROUTINE INTEBc9sdkm(gp,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                       isrc,nsipex,he,heyz,hezx,hexy)
C
C     $: Integracija Robni 9 tockovni element s 4 tockovno geometrijo
C        za single domain enacbo kinematike
C
C -----------------------------------------------------------------------------      
      USE inc_types 
      TYPE(gausstype) :: gp
      
      INTEGER i,j,k,isrc,isip,nsipex,ng1s,ng2s,ng1,ng2
c
      REAL(8) xp,yp,zp,
     &        x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &        xc0,yc0,zc0,xet,yet,zet,xks,yks,zks,
     &        xx1,yy1,
     &        rx,ry,rz,ra,ira2,
     &        ro,roth,th,ajac,
     &        pi,pi2,
     &        eti,etj,eta1m,eta2m,eta1p,eta2p,
     &        anx,any,anz,anx1,any1,anz1,
     &        fung,funx,funy,funz,funh,
     &        dgamma
      REAL(8) he(nsipex)
      REAL(8) hexy(nsipex),heyz(nsipex),hezx(nsipex)
      REAL(8) fih(9),fig4(4)      
      REAL(8) al(4),fii(4),th0(4),th1(4),ksi(13),eta(13)
c
c      integral divison
      REAL(8) a,b,c,d,dex,gii,gij
c      REAL(8) d1,d2,d3,d4,d13max,d24max,minedge
      INTEGER idivXi,ndivXi,idivEt,ndivEt,ising

c
      DATA ksi / 0.0D00, 0.5D00, 1.0D00, 1.0D00, 1.0D00, 0.5D00, 0.0D00, 0.0D00, 0.5D00, 0.125D0, 0.875D0, 0.875D0, 0.125D0/
      DATA eta / 0.0D00, 0.0D00, 0.0D00, 0.5D00, 1.0D00, 1.0D00, 1.0D00, 0.5D00, 0.5D00, 0.125D0, 0.125D0, 0.875D0, 0.875D0/
C
C*** SET NUMBER OF INTEGRATION POINTS
C
c     singular      
      ng1s=gp%ng1(gp%kmBs)
      ng2s=gp%ng2(gp%kmBs)
c     regular      
      ng1=gp%ng1(gp%kmBr)
      ng2=gp%ng2(gp%kmBr)
C
C***  9 NODE CONTINUOUS BOUNDARY ELEMENT
C
      PI=2.0D0*ASIN(1.0D0)
      PI2=2.0D0*PI
c
      he=0.00D00
      hexy=0.00D00
      heyz=0.00D00
      hezx=0.00D00

C
C     Integral razdelimo tako, da integriramo po kvdadratkih
C
c      d1=SQRT((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
c      d2=SQRT((x2-x3)**2+(y2-y3)**2+(z2-z3)**2)
c      d3=SQRT((x3-x4)**2+(y3-y4)**2+(z3-z4)**2)
c      d4=SQRT((x4-x1)**2+(y4-y1)**2+(z4-z1)**2)   
c      d13max=max(d1,d3)
c      d24max=max(d2,d4)      

      ndivXi=1 !max(2,INT(d13max/minedge+0.5D0))
      ndivEt=1 !max(2,INT(d24max/minedge+0.5D0))        
c
c     glavna integracijska zanka
c      
c      a=-1.0D0
c      b=+1.0D0
c      c=-1.0D0
c      d=+1.0D0                  
c      dex=1.0D0      
      DO idivXi=1,ndivXi
        DO idivEt=1,ndivEt
          a=-1.0D0+(idivXi-1)*2.0D0/ndivXi
          b=-1.0D0+(idivXi)*2.0D0/ndivXi          
          c=-1.0D0+(idivEt-1)*2.0D0/ndivEt
          d=-1.0D0+(idivEt)*2.0D0/ndivEt
          dex=0.25D0*(b-a)*(d-c)        
          IF (isrc.NE.0) THEN
            XX1=ksi(isrc)  ! med 0 in 1
            XX1=-1.0D0+2.0D0*XX1  ! med -1 in 1
            XX1=(XX1-a)/(b-a)  ! med 0 in 1 v intervalu a,b
      
            YY1=eta(isrc)  ! med 0 in 1
            YY1=-1.0D0+2.0D0*YY1  ! med -1 in 1
            YY1=(YY1-c)/(d-c)  ! med 0 in 1 v intervalu c,d
          
            IF (XX1.GE.0.0D0.AND.XX1.LE.1.0D0.AND.YY1.GE.0.0D0.AND.YY1.LE.1.0D0) THEN
              ising=1
            ELSE
              ising=0
            END IF
          ELSE
            ising=0
          ENDIF         
c
C*** SINGULAR INTEGRALS
c
c      IF (isrc.NE.0) THEN
      IF (ising.NE.0) THEN      
c        XX1=ksi(isrc)  ! to stran ce imas razdelitev na podobmocja
c        YY1=eta(isrc)
c
        DO K=1,4
c
          IF(K.EQ.1) THEN
            IF (1.0D0-XX1.EQ.0.0D0) GOTO 1000
            AL(K)=SQRT(YY1**2+(1.0D0-XX1)**2)
            FII(K)=ACOS(YY1/AL(K))
            TH0(K)=1.5D0*PI+ACOS(YY1/AL(K))
            TH1(K)=PI2+ATAN((1.0D0-YY1)/(1.0D0-XX1))
c
          ELSE IF(K.EQ.2) THEN
            IF (1.0D0-YY1.EQ.0.0D0) GOTO 1000
            AL(K)=SQRT((1.0D0-YY1)**2+(1.0D0-XX1)**2)
            FII(K)=ACOS((1.0D0-XX1)/AL(K))
            TH0(K)=ASIN((1.0D0-YY1)/AL(K))
            TH1(K)=0.5D0*PI+ATAN(XX1/(1.0D0-YY1))
c
          ELSE IF(K.EQ.3) THEN
            IF (XX1.EQ.0.0D0) GOTO 1000
            AL(K)=SQRT(XX1**2+(1.0D0-YY1)**2)
            FII(K)=ACOS((1.0D0-YY1)/AL(K))
            TH0(K)=0.5D0*PI+ACOS((1.0D0-YY1)/AL(K))
            TH1(K)=PI+ATAN(YY1/XX1)
c
          ELSE
            IF (YY1.EQ.0.0D0) GOTO 1000
            AL(K)=SQRT(YY1**2+XX1**2)
            FII(K)=ACOS(XX1/AL(K))
            TH0(K)=PI+ACOS(XX1/AL(K))
            TH1(K)=1.5D0*PI+ATAN((1.0D0-XX1)/YY1)
          END IF
C
C***      GAUSS INTEGRATION 
C
          DO I=ng1s,ng2s
            DO J=ng1s,ng2s
c
              TH=(TH1(K)-TH0(K))*gp%GI(I)/2.0D0+(TH1(K)+TH0(K))/2.0D0
              ROTH=AL(K)*SIN(FII(K))/SIN(TH-TH0(K)+FII(K))
              RO=ROTH*gp%GI(J)/2.0D0+ROTH/2.0D0
c
              ETI=2.0D0*(XX1+RO*COS(TH))-1.0D0
              ETJ=2.0D0*(YY1+RO*SIN(TH))-1.0D0
c
              ETI=a+0.5D0*(ETI+1.0D0)*(b-a)     
              ETJ=c+0.5D0*(ETJ+1.0D0)*(d-c)                   
c
              ETA1M=1.D0-ETI
              ETA1P=1.D0+ETI
              ETA2M=1.D0-ETJ
              ETA2P=1.D0+ETJ
 
              FIG4(1)=0.25D0*ETA1M*ETA2M
              FIG4(2)=0.25D0*ETA1P*ETA2M
              FIG4(3)=0.25D0*ETA1P*ETA2P
              FIG4(4)=0.25D0*ETA1M*ETA2P
c
              XC0=FIG4(1)*X1+FIG4(2)*X2+FIG4(3)*X3+FIG4(4)*X4
              YC0=FIG4(1)*Y1+FIG4(2)*Y2+FIG4(3)*Y3+FIG4(4)*Y4
              ZC0=FIG4(1)*Z1+FIG4(2)*Z2+FIG4(3)*Z3+FIG4(4)*Z4
c
              XKS=0.25D0*(-ETA2M*X1+ETA2M*X2+ETA2P*X3-ETA2P*X4)
              YKS=0.25D0*(-ETA2M*Y1+ETA2M*Y2+ETA2P*Y3-ETA2P*Y4)
              ZKS=0.25D0*(-ETA2M*Z1+ETA2M*Z2+ETA2P*Z3-ETA2P*Z4)
C
              XET=0.25D0*(-ETA1M*X1-ETA1P*X2+ETA1P*X3+ETA1M*X4)
              YET=0.25D0*(-ETA1M*Y1-ETA1P*Y2+ETA1P*Y3+ETA1M*Y4)
              ZET=0.25D0*(-ETA1M*Z1-ETA1P*Z2+ETA1P*Z3+ETA1M*Z4)
c
              ANX=YKS*ZET-YET*ZKS
              ANY=XET*ZKS-XKS*ZET
              ANZ=XKS*YET-XET*YKS
c
              AJAC=1.0D0/SQRT(ANX**2+ANY**2+ANZ**2)
c
              ANX1=ANX*AJAC
              ANY1=ANY*AJAC
              ANZ1=ANZ*AJAC
c
              AJAC=dex*(TH1(K)-TH0(K))*ROTH*RO*gp%OME(I)*gp%OME(J)/AJAC
c
              RX=XP-XC0
              RY=YP-YC0
              RZ=ZP-ZC0
              RA=SQRT(RX*RX+RY*RY+RZ*RZ)
c             IF (RA.LT.1.0D-12) RA=1.0D-12
              IRA2=1.0D00/(RA*RA)
C
C***          ELLIPTIC LAPLACE FUNDAMENTAL SOLUTION
C 
              FUNG=0.25D0/(PI*RA)
              FUNX=(XP-XC0)*IRA2*FUNG
              FUNY=(YP-YC0)*IRA2*FUNG
              FUNZ=(ZP-ZC0)*IRA2*FUNG
              FUNH=FUNX*ANX1+FUNY*ANY1+FUNZ*ANZ1

c             calculate H continous shape functions
              CALL cshape9(fih,ETI,ETJ,nsipex)
              DO isip=1,nsipex
                DGAMMA=AJAC*FIH(isip)
                he(isip)=he(isip)+FUNH*DGAMMA                            
                hexy(isip)=hexy(isip)-(FUNX*ANY1-FUNY*ANX1)*DGAMMA
                heyz(isip)=heyz(isip)-(FUNY*ANZ1-FUNZ*ANY1)*DGAMMA
                hezx(isip)=hezx(isip)-(FUNZ*ANX1-FUNX*ANZ1)*DGAMMA
              END DO
            END DO
          END DO
 1000   END DO
c      
      ELSE
C
C*** REGULAR INTEGRALS
C
        DO i=ng1,ng2
          gii=a+0.5D0*(gp%GI(I)+1.0D0)*(b-a)
          ETA1M=1.0D0-gii
          ETA1P=1.0D0+gii
          DO j=ng1,ng2
            gij=c+0.5D0*(gp%GI(J)+1.0D0)*(d-c)                   
            ETA2M=1.0D0-gij
            ETA2P=1.0D0+gij
C
            FIG4(1)=0.25D0*ETA1M*ETA2M
            FIG4(2)=0.25D0*ETA1P*ETA2M
            FIG4(3)=0.25D0*ETA1P*ETA2P
            FIG4(4)=0.25D0*ETA1M*ETA2P
c
            XC0=FIG4(1)*X1+FIG4(2)*X2+FIG4(3)*X3+FIG4(4)*X4
            YC0=FIG4(1)*Y1+FIG4(2)*Y2+FIG4(3)*Y3+FIG4(4)*Y4
            ZC0=FIG4(1)*Z1+FIG4(2)*Z2+FIG4(3)*Z3+FIG4(4)*Z4
c
            XKS=0.25D0*(-ETA2M*X1+ETA2M*X2+ETA2P*X3-ETA2P*X4)
            YKS=0.25D0*(-ETA2M*Y1+ETA2M*Y2+ETA2P*Y3-ETA2P*Y4)
            ZKS=0.25D0*(-ETA2M*Z1+ETA2M*Z2+ETA2P*Z3-ETA2P*Z4)
C
            XET=0.25D0*(-ETA1M*X1-ETA1P*X2+ETA1P*X3+ETA1M*X4)
            YET=0.25D0*(-ETA1M*Y1-ETA1P*Y2+ETA1P*Y3+ETA1M*Y4)
            ZET=0.25D0*(-ETA1M*Z1-ETA1P*Z2+ETA1P*Z3+ETA1M*Z4)
c
            ANX=YKS*ZET-YET*ZKS
            ANY=XET*ZKS-XKS*ZET
            ANZ=XKS*YET-XET*YKS
c
            AJAC=1.0D0/SQRT(ANX**2+ANY**2+ANZ**2)
c
            ANX1=ANX*AJAC
            ANY1=ANY*AJAC
            ANZ1=ANZ*AJAC
c
            AJAC=dex*gp%OME(I)*gp%OME(J)/AJAC
c
            RX=XP-XC0
            RY=YP-YC0
            RZ=ZP-ZC0
            RA=SQRT(RX*RX+RY*RY+RZ*RZ)
c           IF (RA.LT.1.0D-12) RA=1.0D-12
            IRA2=1.0D00/(RA*RA)
C
C***        ELLIPTIC LAPLACE FUNDAMENTAL SOLUTION
C     
            FUNG=0.25D0/(PI*RA)
            FUNX=(XP-XC0)*IRA2*FUNG
            FUNY=(YP-YC0)*IRA2*FUNG
            FUNZ=(ZP-ZC0)*IRA2*FUNG
            FUNH=FUNX*ANX1+FUNY*ANY1+FUNZ*ANZ1
c
c           calculate H continous shape functions
            CALL cshape9(fih,gii,gij,nsipex)            
            DO isip=1,nsipex
              DGAMMA=AJAC*FIh(isip)
              he(isip)=he(isip)+FUNH*DGAMMA
              hexy(isip)=hexy(isip)-(FUNX*ANY1-FUNY*ANX1)*DGAMMA
              heyz(isip)=heyz(isip)-(FUNY*ANZ1-FUNZ*ANY1)*DGAMMA
              hezx(isip)=hezx(isip)-(FUNZ*ANX1-FUNX*ANZ1)*DGAMMA
            END DO           
          END DO
        END DO
      END IF
      
        END DO
      END DO        
      
c            
      
      END


C -----------------------------------------------------------------------------
      SUBROUTINE sdKm_corInt(env,io,inp,mesh,gauss,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz)
C
C     $: Compute or Read
C        Matrices, single domain kinematics eqaution
C        (H, Htx, Hty, Htz, Dx, Dy, Dz)
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(IOtype) io
      TYPE(inputtype) inp   
      TYPE(penv) :: env           
      
      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral 
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni
      
      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDy(env%zac:env%kon,mesh%nnodes) ! to je za drugi integral 
      REAL(8) sdDz(env%zac:env%kon,mesh%nnodes) ! na desni       
      REAL(8) ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz ! ocena natancnosti integralov
      REAL(8) ontHtxDz,ontHtyDx,ontHtzDy


      INTEGER iok
      
C     Check if integrals exist on file
      CALL CheckSDIntFile(io,inp%INTversion,mesh%nnodes,mesh%nbnodes,env%nproc,iok)      

      IF (iok.EQ.1) THEN
c       Read integrals from disk 
        CALL WarnErr(env,io,inp,0,"corSDin_kmSteady","Reading single domain integrals!",0)   
        CALL ReadSDIntDisk(env,io,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &                      ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)       
c       communicate errors (in case some of the procesors calculated integrals)
        CALL CommErr(env,ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)

      ELSE
c       Calcualte integrals      
        CALL WarnErr(env,io,inp,0,"corSDin_kmSteady","Calculating single domain integrals!",0)     
        CALL FMATsdkm(env,mesh,gauss,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &                ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)
c       communicate errors
        CALL CommErr(env,ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)
          
        CALL WriteSDIntDisk(env,io,inp,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &                      ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)

      END IF
      
      IF (env%myproc.EQ.1) THEN
        WRITE(io%l,'(A)') "Ocena natancnosti single domain integralov!"
        WRITE(io%l,'(A,G15.10)') "H  = ",ontH 
        WRITE(io%l,'(A,G15.10)') "Htx= ",ontHtx 
        WRITE(io%l,'(A,G15.10)') "Hty= ",ontHty
        WRITE(io%l,'(A,G15.10)') "Htz= ",ontHtz
        WRITE(io%l,'(A,G15.10)') "Dx = ",ontDx 
        WRITE(io%l,'(A,G15.10)') "Dy = ",ontDy
        WRITE(io%l,'(A,G15.10)') "Dz = ",ontDz
        WRITE(io%l,'(A)') "Po popravljanju singularnih integralov:"
        WRITE(io%l,'(A,G15.10)') "HtxDz = ",ontHtxDz
        WRITE(io%l,'(A,G15.10)') "HtyDx = ",ontHtyDx
        WRITE(io%l,'(A,G15.10)') "HtzDy = ",ontHtzDy
        WRITE(io%l,'(A)') ""
      END IF

      END

C -----------------------------------------------------------------------------
      SUBROUTINE sdKm_corIntRot(env,io,inp,mesh,gauss,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz)
C
C     $: Compute or Read
C        Matrices, single domain kinematics eqaution
C        (H, Htx, Hty, Htz, Dx, Dy, Dz)
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(IOtype) io
      TYPE(inputtype) inp
      TYPE(penv) :: env

      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni

      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDy(env%zac:env%kon,mesh%nnodes) ! to je za drugi integral
      REAL(8) sdDz(env%zac:env%kon,mesh%nnodes) ! na desni
      REAL(8) ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz ! ocena natancnosti integralov
      REAL(8) ontHtxDz,ontHtyDx,ontHtzDy


      INTEGER iok

C     Check if integrals exist on file
      CALL CheckSDIntFile(io,inp%INTversion,mesh%nnodes,mesh%nbnodes,env%nproc,iok)

      IF (iok.EQ.1) THEN
c       Read integrals from disk
        CALL WarnErr(env,io,inp,0,"corSDin_kmSteady","Reading single domain integrals!",0)
        CALL ReadSDIntDisk(env,io,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &                      ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)
c       communicate errors (in case some of the procesors calculated integrals)
        CALL CommErr(env,ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)

      ELSE
c       Calcualte integrals
        CALL WarnErr(env,io,inp,0,"corSDin_kmSteady","Calculating single domain integrals!",0)
        CALL FMATsdkmRot(env,mesh,gauss,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &                ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)
c       communicate errors
        CALL CommErr(env,ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)

        CALL WriteSDIntDisk(env,io,inp,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &                      ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)

      END IF

      IF (env%myproc.EQ.1) THEN
        WRITE(io%l,'(A)') "Ocena natancnosti single domain integralov!"
        WRITE(io%l,'(A,G15.10)') "H  = ",ontH
        WRITE(io%l,'(A,G15.10)') "Htx= ",ontHtx
        WRITE(io%l,'(A,G15.10)') "Hty= ",ontHty
        WRITE(io%l,'(A,G15.10)') "Htz= ",ontHtz
        WRITE(io%l,'(A,G15.10)') "Dx = ",ontDx
        WRITE(io%l,'(A,G15.10)') "Dy = ",ontDy
        WRITE(io%l,'(A,G15.10)') "Dz = ",ontDz
        WRITE(io%l,'(A)') "Po popravljanju singularnih integralov:"
        WRITE(io%l,'(A,G15.10)') "HtxDz = ",ontHtxDz
        WRITE(io%l,'(A,G15.10)') "HtyDx = ",ontHtyDx
        WRITE(io%l,'(A,G15.10)') "HtzDy = ",ontHtzDy
        WRITE(io%l,'(A)') ""
      END IF

      END
      


C -----------------------------------------------------------------------------
      SUBROUTINE sdKm_corInt8(env,io,inp,mesh,gauss,sdH,sdHtx,sdHty,sdHtz,sdDx8,sdDy8,sdDz8,sdDxB,sdDyB,sdDzB)
C
C     $: Compute or Read
C        Matrices, single domain kinematics eqaution
C        (H, Htx, Hty, Htz, Dx, Dy, Dz)
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(IOtype) io
      TYPE(inputtype) inp
      TYPE(penv) :: env

      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni

      REAL(8) sdDx8(env%zac:env%kon,mesh%nnodes8) ! doktorat, enacba (4.7)
      REAL(8) sdDy8(env%zac:env%kon,mesh%nnodes8) ! to je za drugi integral
      REAL(8) sdDz8(env%zac:env%kon,mesh%nnodes8) ! na desni

      REAL(8) sdDxB(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDyB(env%zac:env%kon,mesh%nbnodes) ! to je za drugi integral
      REAL(8) sdDzB(env%zac:env%kon,mesh%nbnodes) ! na desni

      REAL(8) ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz ! ocena natancnosti integralov
      REAL(8) ontHtxDz,ontHtyDx,ontHtzDy


      INTEGER iok

C     Check if integrals exist on file
      CALL CheckSDIntFile8(io,inp%INTversion,mesh%nnodes,mesh%nbnodes,env%nproc,iok)

      IF (iok.EQ.1) THEN
c       Read integrals from disk
        CALL WarnErr(env,io,inp,0,"corSDin_kmSteady","Reading single domain integrals!",0)
        CALL ReadSDIntDisk8(env,io,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx8,sdDy8,sdDz8,sdDxB,sdDyB,sdDzB,
     &                      ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)
c       communicate errors (in case some of the procesors calculated integrals)
        CALL CommErr(env,ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)

      ELSE
c       Calcualte integrals
        CALL WarnErr(env,io,inp,0,"corSDin_kmSteady","Calculating single domain integrals!",0)
        CALL FMATsdkm8(env,mesh,gauss,sdH,sdHtx,sdHty,sdHtz,sdDx8,sdDy8,sdDz8,sdDxB,sdDyB,sdDzB,
     &                ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)
c       communicate errors
        CALL CommErr(env,ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)

        CALL WriteSDIntDisk8(env,io,inp,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx8,sdDy8,sdDz8,sdDxB,sdDyB,sdDzB,
     &                      ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)

      END IF

      IF (env%myproc.EQ.1) THEN
        WRITE(io%l,'(A)') "Ocena natancnosti single domain integralov!"
        WRITE(io%l,'(A,G15.10)') "H  = ",ontH
        WRITE(io%l,'(A,G15.10)') "Htx= ",ontHtx
        WRITE(io%l,'(A,G15.10)') "Hty= ",ontHty
        WRITE(io%l,'(A,G15.10)') "Htz= ",ontHtz
        WRITE(io%l,'(A,G15.10)') "Dx = ",ontDx
        WRITE(io%l,'(A,G15.10)') "Dy = ",ontDy
        WRITE(io%l,'(A,G15.10)') "Dz = ",ontDz
        WRITE(io%l,'(A)') "Po popravljanju singularnih integralov:"
        WRITE(io%l,'(A,G15.10)') "HtxDz = ",ontHtxDz
        WRITE(io%l,'(A,G15.10)') "HtyDx = ",ontHtyDx
        WRITE(io%l,'(A,G15.10)') "HtzDy = ",ontHtzDy
        WRITE(io%l,'(A)') ""
      END IF

      END

C -----------------------------------------------------------------------------
      SUBROUTINE CheckSDIntFile(io,a,b,c,d,iok)
C -----------------------------------------------------------------------------      
      USE inc_types
      TYPE(IOtype) io 
      INTEGER a,b,c,d
      INTEGER aa,bb,cc,dd,iok

      iok=0
      OPEN (io%ikm,FILE=TRIM(io%ikm_name),FORM='UNFORMATTED',STATUS='OLD',ERR=10)
      READ(io%ikm) aa
      READ(io%ikm) bb,cc,dd
            
      IF (a.EQ.aa.AND.b.EQ.bb.AND.c.EQ.cc.AND.d.EQ.dd) iok=1    
      CLOSE(io%ikm)
      
10    RETURN
      END    

C -----------------------------------------------------------------------------
      SUBROUTINE CheckSDIntFile8(io,a,b,c,d,iok)
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(IOtype) io
      INTEGER a,b,c,d
      INTEGER aa,bb,cc,dd,iok

      iok=0
      OPEN (io%ikmL,FILE=TRIM(io%ikmL_name),FORM='UNFORMATTED',STATUS='OLD',ERR=10)
      READ(io%ikmL) aa
      READ(io%ikmL) bb,cc,dd

      IF (a.EQ.aa.AND.b.EQ.bb.AND.c.EQ.cc.AND.d.EQ.dd) iok=1
      CLOSE(io%ikmL)

10    RETURN
      END
      
C -----------------------------------------------------------------------------
      SUBROUTINE WriteSDIntDisk(env,io,inp,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &                          ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)       
C
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      TYPE(IOtype) io
      TYPE(inputtype) inp        
      TYPE(penv) :: env
      
      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral 
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni
      
      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDy(env%zac:env%kon,mesh%nnodes) ! to je za drugi integral 
      REAL(8) sdDz(env%zac:env%kon,mesh%nnodes) ! na desni       
      REAL(8) ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy ! ocena natancnosti integralov     

      OPEN (io%ikm,FILE=TRIM(io%ikm_name),FORM='UNFORMATTED',STATUS='UNKNOWN')      
      WRITE(io%ikm) inp%INTversion
      WRITE(io%ikm) mesh%nnodes,mesh%nbnodes,env%nproc
     
      CALL WrSMat(sdH,env%nmn,mesh%nbnodes,io%ikm)
      CALL WrSMat(sdHtx,env%nmn,mesh%nbnodes,io%ikm)
      CALL WrSMat(sdHty,env%nmn,mesh%nbnodes,io%ikm)
      CALL WrSMat(sdHtz,env%nmn,mesh%nbnodes,io%ikm)
      CALL WrSMat(sdDx,env%nmn,mesh%nnodes,io%ikm)
      CALL WrSMat(sdDy,env%nmn,mesh%nnodes,io%ikm)
      CALL WrSMat(sdDz,env%nmn,mesh%nnodes,io%ikm)  
      WRITE(io%ikm) ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy                                              
      
      CLOSE(io%ikm)          
      
      END
C -----------------------------------------------------------------------------
      SUBROUTINE ReadSDIntDisk(env,io,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,
     &                         ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)       
C
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      TYPE(IOtype) io
      TYPE(penv) :: env
      INTEGER a,b,c,d      
      
      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral 
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni
      
      REAL(8) sdDx(env%zac:env%kon,mesh%nnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDy(env%zac:env%kon,mesh%nnodes) ! to je za drugi integral 
      REAL(8) sdDz(env%zac:env%kon,mesh%nnodes) ! na desni       
      REAL(8) ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy ! ocena natancnosti integralov     

      OPEN (io%ikm,FILE=TRIM(io%ikm_name),FORM='UNFORMATTED',STATUS='OLD')      
      READ(io%ikm) a
      READ(io%ikm) b,c,d
     
      CALL RdSMat(sdH,env%nmn,mesh%nbnodes,io%ikm)
      CALL RdSMat(sdHtx,env%nmn,mesh%nbnodes,io%ikm)
      CALL RdSMat(sdHty,env%nmn,mesh%nbnodes,io%ikm)
      CALL RdSMat(sdHtz,env%nmn,mesh%nbnodes,io%ikm)
      CALL RdSMat(sdDx,env%nmn,mesh%nnodes,io%ikm)
      CALL RdSMat(sdDy,env%nmn,mesh%nnodes,io%ikm)
      CALL RdSMat(sdDz,env%nmn,mesh%nnodes,io%ikm)  
      READ(io%ikm) ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy                                              
      
      CLOSE(io%ikm)          
      
      END


C -----------------------------------------------------------------------------
      SUBROUTINE WriteSDIntDisk8(env,io,inp,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx8,sdDy8,sdDz8,sdDxB,sdDyB,sdDzB,
     &                          ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)
C
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(IOtype) io
      TYPE(inputtype) inp
      TYPE(penv) :: env

      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni

      REAL(8) sdDx8(env%zac:env%kon,mesh%nnodes8) ! doktorat, enacba (4.7)
      REAL(8) sdDy8(env%zac:env%kon,mesh%nnodes8) ! to je za drugi integral
      REAL(8) sdDz8(env%zac:env%kon,mesh%nnodes8) ! na desni

      REAL(8) sdDxB(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDyB(env%zac:env%kon,mesh%nbnodes) ! to je za drugi integral
      REAL(8) sdDzB(env%zac:env%kon,mesh%nbnodes) ! na desni

      REAL(8) ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy ! ocena natancnosti integralov

      OPEN (io%ikmL,FILE=TRIM(io%ikmL_name),FORM='UNFORMATTED',STATUS='UNKNOWN')
      WRITE(io%ikmL) inp%INTversion
      WRITE(io%ikmL) mesh%nnodes,mesh%nbnodes,env%nproc

      CALL WrSMat(sdH,env%nmn,mesh%nbnodes,io%ikmL)
      CALL WrSMat(sdHtx,env%nmn,mesh%nbnodes,io%ikmL)
      CALL WrSMat(sdHty,env%nmn,mesh%nbnodes,io%ikmL)
      CALL WrSMat(sdHtz,env%nmn,mesh%nbnodes,io%ikmL)
      CALL WrSMat(sdDx8,env%nmn,mesh%nnodes8,io%ikmL)
      CALL WrSMat(sdDy8,env%nmn,mesh%nnodes8,io%ikmL)
      CALL WrSMat(sdDz8,env%nmn,mesh%nnodes8,io%ikmL)
      CALL WrSMat(sdDxB,env%nmn,mesh%nbnodes,io%ikmL)
      CALL WrSMat(sdDyB,env%nmn,mesh%nbnodes,io%ikmL)
      CALL WrSMat(sdDzB,env%nmn,mesh%nbnodes,io%ikmL)

      WRITE(io%ikmL) ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy

      CLOSE(io%ikmL)

      END
C -----------------------------------------------------------------------------
      SUBROUTINE ReadSDIntDisk8(env,io,mesh,sdH,sdHtx,sdHty,sdHtz,sdDx8,sdDy8,sdDz8,sdDxB,sdDyB,sdDzB,
     &                         ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy)
C
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(IOtype) io
      TYPE(penv) :: env
      INTEGER a,b,c,d

      REAL(8) sdH(env%zac:env%kon,mesh%nbnodes)  ! diffusion

      REAL(8) sdHtx(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdHty(env%zac:env%kon,mesh%nbnodes) ! to je za prvi integral
      REAL(8) sdHtz(env%zac:env%kon,mesh%nbnodes) ! na desni

      REAL(8) sdDx8(env%zac:env%kon,mesh%nnodes8) ! doktorat, enacba (4.7)
      REAL(8) sdDy8(env%zac:env%kon,mesh%nnodes8) ! to je za drugi integral
      REAL(8) sdDz8(env%zac:env%kon,mesh%nnodes8) ! na desni

      REAL(8) sdDxB(env%zac:env%kon,mesh%nbnodes) ! doktorat, enacba (4.7)
      REAL(8) sdDyB(env%zac:env%kon,mesh%nbnodes) ! to je za drugi integral
      REAL(8) sdDzB(env%zac:env%kon,mesh%nbnodes) ! na desni

      REAL(8) ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy ! ocena natancnosti integralov

      OPEN (io%ikmL,FILE=TRIM(io%ikmL_name),FORM='UNFORMATTED',STATUS='OLD')
      READ(io%ikmL) a
      READ(io%ikmL) b,c,d

      CALL RdSMat(sdH,env%nmn,mesh%nbnodes,io%ikmL)
      CALL RdSMat(sdHtx,env%nmn,mesh%nbnodes,io%ikmL)
      CALL RdSMat(sdHty,env%nmn,mesh%nbnodes,io%ikmL)
      CALL RdSMat(sdHtz,env%nmn,mesh%nbnodes,io%ikmL)
      CALL RdSMat(sdDx8,env%nmn,mesh%nnodes8,io%ikmL)
      CALL RdSMat(sdDy8,env%nmn,mesh%nnodes8,io%ikmL)
      CALL RdSMat(sdDz8,env%nmn,mesh%nnodes8,io%ikmL)
      CALL RdSMat(sdDxB,env%nmn,mesh%nbnodes,io%ikmL)
      CALL RdSMat(sdDyB,env%nmn,mesh%nbnodes,io%ikmL)
      CALL RdSMat(sdDzB,env%nmn,mesh%nbnodes,io%ikmL)
      READ(io%ikmL) ontH,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtxDz,ontHtyDx,ontHtzDy

      CLOSE(io%ikmL)

      END
