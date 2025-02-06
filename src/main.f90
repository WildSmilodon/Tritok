! Tritok, a 3D BEM based fluid flow solver
! Copyright (C) 2008 - 2025  Jure Ravnik
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

!
!     Start of TriTok CFD code
!
      PROGRAM Tritok
!
!      
!      
!           3D BEM based fluid flow solver
!       
!           written by 
!
!           - Jure Ravnik, jure.ravnik@um.si
!
!           University of Maribor
!           Faculty of Mechanical Engineering
!           Smetanova 17, SI-2000, Maribor, Slovenia
!
!     kode = 0  - resuje kinematiko za robne vrtincnosti, predpise hitrost iz BIC
!     kode = 1  - NE resuje kinematike za robne vrtincnosti, predpise hitrost iz BIC
!     kode = 2  - NE resuje kinematike za robne vrtincnosti, open B.C. za hitrost
!     kode = 3  - resuje kinematiko za robne vrtincnosti, open B.C. za hitrost
!
!     EQNS iBw iDv iDw iDT iDC imH
!
!     iBw = 1  - single domain resitev enacbe kinematike za robne vrtincnosti
!     iBw = 2  - single domain resitev enacbe kinematike za robne vrtincnosti
!                z linearno aproksimacijo za obmocni clen (DELA SLABO)
!     iBw = 3  - single domain resitev enacbe kinematike za robne vrtincnosti
!                s pomocjo rotacijskih matrik
!     iDv = 1  - least squares resitev predolocenega sistema enacb za subdomain kinematiko
!                za obmocne hitrosti
!     iDv = 2  - resitev dolocenega sistema enacb za subdomain kinematiko
!                za obmocne hitrosti
!     iDw = 1  - least squares resitev predolocenega sistema enacb za subdomain kinetiko
!                za obmocne vrtincnosti
!     iDw = 2  - resitev dolocenega sistema enacb za subdomain kinetiko za obmocne vrtincnosti
!     iDw = 3  - least squares resitev predolocenega sistema enacb za subdomain kinetiko
!                za obmocne vrtincnosti kjer sta viskoznost in beta spremenljivi
!     iDT = 1  - least squares resitev predolocenega sistema enacb za subdomain kinetiko
!                za obmocne temperature
!     iDT = 2  - resitev dolocenega sistema enacb za subdomain kinetiko za obmocne temperature
!     iDT = 3  - resitev predolocenega sistema enacb za subdomain kinetiko za obmocne temperature
!                z spremeljivimi snovnimi lastnostmi
!     iDC = 1  - moj (laplace) least squares resitev predolocenega sistema enacb za subdomain
!                kinetiko za DC enacbo s spremenljivo difuzivnostjo
!     iDC = 2  - parametrix least squares resitev predolocenega sistema enacb za subdomain
!                kinetiko za DC enacbo s spremenljivo difuzivnostjo
!     iDC = 3  - parametrix resitev dolocenega sistema enacb za subdomain
!                kinetiko za DC enacbo s spremenljivo difuzivnostjo
!     iDC = 4  - parametrix least squares resitev predolocenega sistema enacb za subdomain
!                kinetiko za DC enacbo s spremenljivo difuzivnostjo, hitrost je v sistemski,
!                rabi samo eno itreracijo v linearnih primerih.
!     imH = 1  - Hemlholtz equation (lap u - mu^2 u = f)
!     imH = 2  - (subdomains) diffusion convection with variable diff solved as Hemlholtz equation
!     imH = 3  - (single) diffusion convection with variable diff solved as Hemlholtz equation
!     imH = 4  - (single+matmul) diffusion convection with variable diff solved as Hemlholtz equation
!     imH = 5  - (single+wawt) diffusion convection with variable diff solved as Hemlholtz equation
!     imH = 6  - (single+aca) diffusion convection with variable diff solved as Hemlholtz equation
!     imH = 7  - (single+Hmat+(no ali WT ali aca)) diffusion convection with variable diff solved as Hemlholtz equation
!
!
      USE inc_types
      IMPLICIT NONE

      TYPE(meshType)  :: mesh
      TYPE(gausstype) :: gauss
      TYPE(IOtype)    :: io
      TYPE(inputtype) :: inp      
      TYPE(KMpointer) :: kmp
      TYPE(CPUtype)   :: cpu
      TYPE(countType) :: cnt    
      TYPE(penv) :: env
      TYPE(particleType), ALLOCATABLE :: part(:) ! seznam delcev
      
      TYPE(matrix) :: VTEsysmX,VTErhsmX,VTEprecX ! vorticity transport eqaution sysm in rhsm X
      TYPE(matrix) :: VTEsysmY,VTErhsmY,VTEprecY
      TYPE(matrix) :: VTEsysmZ,VTErhsmZ,VTEprecZ
      TYPE(matrix) :: Tsysm,Trhsm,Tprec      ! energy equation
      TYPE(matrix) :: DCsysm,DCrhsm,DCprec      ! general diff-conv. equation
      TYPE(matrix) :: mHsysm,mHrhsm,mHprec      ! modified Helmholtz equation
      TYPE(matrix) :: KMsysm,KMrhsm,KMprec     ! kinematics - for domain velocities
!
!     Macro integrals matrices
!
      REAL(8), ALLOCATABLE :: smatH(:,:), smatG(:,:)
      REAL(8), ALLOCATABLE :: smatB(:,:) ! source term (and FD time term)
      REAL(8), ALLOCATABLE :: smatAbdX(:,:),smatAbdY(:,:),smatAbdZ(:,:) ! advection boundary-domain
      REAL(8), ALLOCATABLE :: smatHtx(:,:),smatHty(:,:),smatHtz(:,:) ! kinematics, Ht
      REAL(8), ALLOCATABLE :: smatDx(:,:),smatDy(:,:),smatDz(:,:) ! kinematics, D

!     modified Helmholtz equation
      REAL(8), ALLOCATABLE :: mHsmatH(:,:), mHsmatG(:,:)
      REAL(8), ALLOCATABLE :: mHsmatB(:,:) ! source term 
      REAL(8), ALLOCATABLE :: mHprecv(:) ! precondition vector
      REAL(8), ALLOCATABLE :: mHsmatAbdX(:,:),mHsmatAbdY(:,:),mHsmatAbdZ(:,:) ! advection boundary-domain
      REAL(8), ALLOCATABLE :: mHsmatHtx(:,:),mHsmatHty(:,:),mHsmatHtz(:,:) ! kinematics, Ht
      REAL(8), ALLOCATABLE :: mHsmatDx(:,:),mHsmatDy(:,:),mHsmatDz(:,:) ! kinematics, D
!     system, rsh matrices, preconditioner
      REAL(8), ALLOCATABLE :: fmmHA(:,:),fmmHPivot(:)  ! fm sistemska matrika
      REAL(8), ALLOCATABLE :: fmmHRHS(:,:)  ! fm matrika na desni strani
!
!     Single domain kinematics matrices
!
      REAL(8), ALLOCATABLE :: sdH(:,:)
      REAL(8), ALLOCATABLE :: sdHtx(:,:),sdHty(:,:),sdHtz(:,:) ! kinematics, Ht
      REAL(8), ALLOCATABLE :: sdDx(:,:),sdDy(:,:),sdDz(:,:) ! kinematics, D
      REAL(8), ALLOCATABLE :: sdDxB(:,:),sdDyB(:,:),sdDzB(:,:) ! kinematics, D
      REAL(8), ALLOCATABLE :: sdAx(:,:), sdPivotx(:) ! system matrix
      REAL(8), ALLOCATABLE :: sdAy(:,:), sdPivoty(:) ! system matrix
      REAL(8), ALLOCATABLE :: sdAz(:,:), sdPivotz(:) ! system matrix            

      REAL(8), ALLOCATABLE :: KMprecv(:) ! preconditioning vectors
      REAL(8), ALLOCATABLE :: VTEprecvX(:),VTEprecvY(:),VTEprecvZ(:) ! preconditioning vectors      
      REAL(8), ALLOCATABLE :: Tprecv(:), DCprecv(:)      
!
!     Field functions
!     
      REAL(8), ALLOCATABLE :: vorticity(:,:)
      REAL(8), ALLOCATABLE :: qvorticity(:,:)
      REAL(8), ALLOCATABLE :: ptsvorticity(:,:)  ! previous time step
      REAL(8), ALLOCATABLE :: pptsvorticity(:,:) ! pred previous time step               
      REAL(8), ALLOCATABLE :: pnlsvorticity(:,:) ! previous nonlinear iteration   
      REAL(8), ALLOCATABLE :: rhsv(:,:)  ! vorticity right hand side vector 
      
      REAL(8), ALLOCATABLE :: velocity(:,:),pnlsvelocity(:,:)  

      REAL(8), ALLOCATABLE :: temp(:)
      REAL(8), ALLOCATABLE :: qtemp(:)
      REAL(8), ALLOCATABLE :: ptstemp(:)  ! previous time step
      REAL(8), ALLOCATABLE :: pptstemp(:) ! pred previous time step               
      REAL(8), ALLOCATABLE :: pnlstemp(:) ! previous nonlinear iteration   
      REAL(8), ALLOCATABLE :: Trhsv(:)  ! temperature right hand side vector
      REAL(8), ALLOCATABLE :: Ttime(:)  ! temperature source due to previous time step

      REAL(8), ALLOCATABLE :: dc(:)    ! polje za splosno D-C enacbo, s spremenljivo diff.
      REAL(8), ALLOCATABLE :: qdc(:)
      REAL(8), ALLOCATABLE :: ptsdc(:)  ! previous time step
      REAL(8), ALLOCATABLE :: pptsdc(:) ! pred previous time step               
      REAL(8), ALLOCATABLE :: pnlsdc(:) ! previous nonlinear iteration   
      REAL(8), ALLOCATABLE :: DCrhsv(:)  ! DC right hand side vector
      REAL(8), ALLOCATABLE :: DCtime(:)  ! DC time derivative right hand side vector
      REAL(8), ALLOCATABLE :: DCdiff(:)  ! DC difuzivnost
      REAL(8), ALLOCATABLE :: DCqDiff(:)  ! DC difuzivnost v macro tockah
      REAL(8), ALLOCATABLE :: DCdiffGrad(:,:)  ! gradient DC difuzivnost   
      REAL(8), ALLOCATABLE :: DCdiffLap(:)  ! laplace DC difuzivnost
      REAL(8), ALLOCATABLE :: DCanal(:) ! analytical solution, for rms calculation
      REAL(8), ALLOCATABLE :: DCsource(:) ! source

!     modified Helmholtz equation
      REAL(8), ALLOCATABLE :: mHu(:) 
      REAL(8), ALLOCATABLE :: mHq(:),fmqmH(:)
      REAL(8), ALLOCATABLE :: mHrhsV(:) 
      REAL(8), ALLOCATABLE :: mHtime(:)  ! mH time derivative right hand side vector
      REAL(8), ALLOCATABLE :: mHdiff(:)  ! mH difuzivnost
      REAL(8), ALLOCATABLE :: mHqDiff(:)  ! mH difuzivnost v macro tockah
      REAL(8), ALLOCATABLE :: mHdiffGrad(:,:)  ! gradient mH difuzivnost
      REAL(8), ALLOCATABLE :: ptsmHu(:)  ! previous time step
      REAL(8), ALLOCATABLE :: pptsmHu(:) ! pred previous time step               
      REAL(8), ALLOCATABLE :: pnlsmHu(:) ! previous nonlinear iteration
      REAL(8), ALLOCATABLE :: mHanal(:) ! analytical solution, for rms calculation
!
!     Single domain kinetics integrals matrices (for mod Helmholtz fundamental solution)
!
!     U source points
      REAL(8), ALLOCATABLE :: UmatH(:,:), UmatG(:,:)
      REAL(8), ALLOCATABLE :: UmatB(:,:),UmatBp(:,:) ! source term (and FD time term)
      REAL(8), ALLOCATABLE :: UmatAbdX(:,:),UmatAbdY(:,:),UmatAbdZ(:,:) ! advection boundary-domain
      REAL(8), ALLOCATABLE :: UmatDx(:,:),UmatDy(:,:),UmatDz(:,:) ! kinematics, D
!     Q source points
      REAL(8), ALLOCATABLE :: QmatH(:,:), QmatG(:,:)
      REAL(8), ALLOCATABLE :: QmatB(:,:),QmatBp(:,:) ! source term (and FD time term)
      REAL(8), ALLOCATABLE :: QmatAbdX(:,:),QmatAbdY(:,:),QmatAbdZ(:,:) ! advection boundary-domain
      REAL(8), ALLOCATABLE :: QmatDx(:,:),QmatDy(:,:),QmatDz(:,:) ! kinematics, D

!     wavelet compressed version
      TYPE(crsMatrixType) :: cUmatB,cUmatBp,cUmatAbdX,cUmatAbdY,cUmatAbdZ,cUmatDx,cUmatDy,cUmatDz
      TYPE(crsMatrixType) :: cQmatB,cQmatBp,cQmatAbdX,cQmatAbdY,cQmatAbdZ,cQmatDx,cQmatDy,cQmatDz
      REAL(8) WTnnz,WTn

!     aca compressed version
      TYPE(acaMatrixType) :: aUmatB,aUmatBp,aUmatAbdX,aUmatAbdY,aUmatAbdZ,aUmatDx,aUmatDy,aUmatDz
      TYPE(acaMatrixType) :: aQmatB,aQmatBp,aQmatAbdX,aQmatAbdY,aQmatAbdZ,aQmatDx,aQmatDy,aQmatDz

!     H matrix
      TYPE(HmatrixType) :: hUmatB,hUmatBp,hUmatAbdX,hUmatAbdY,hUmatAbdZ,hUmatDx,hUmatDy,hUmatDz
      TYPE(HmatrixType) :: hQmatB,hQmatBp,hQmatAbdX,hQmatAbdY,hQmatAbdZ,hQmatDx,hQmatDy,hQmatDz

      REAL(8), ALLOCATABLE :: tmp(:)  ! zacasno, pomozno polje
!     nanofluids
      REAL(8) nanoA, nanoB, nanoC
      REAL(8), ALLOCATABLE :: PartVolFrac(:) ! particle volume fraction
!     ferrofluids
      REAL(8), ALLOCATABLE :: mfs(:,:)  ! magnetic field strength \vec H
      REAL(8), ALLOCATABLE :: hshg(:,:)  ! \vec H \cdot \vec H, H skalarno H
      REAL(8), ALLOCATABLE :: kbf(:,:)  ! Kelvin body force
      REAL(8), ALLOCATABLE :: tempGrad(:,:)  ! Temperature gradient
!     energy equation with variable material propersites
      REAL(8), ALLOCATABLE :: Density(:)  ! nondimensional density
      REAL(8), ALLOCATABLE :: HeatCapacity(:)  ! nondimensional heat capacity
      REAL(8), ALLOCATABLE :: ThermalConductivity(:)  ! nondimensional thermal conductivity
      REAL(8), ALLOCATABLE :: qThermalConductivity(:)  ! nondimensional thermal conductivity in flux points
      REAL(8), ALLOCATABLE :: gradTC(:,:)     ! gradient of nondimensional thermal conductivity
      REAL(8), ALLOCATABLE :: dRhoCpdt(:)  ! d(rho*cp)/dt
      REAL(8), ALLOCATABLE :: prhocp(:)    ! rho*cp in previous time step
!     vorticity transport equation with variable material propersites
      REAL(8), ALLOCATABLE :: Viscosity(:)  ! nondimensional mu=mu^dim / mu^inReynolds
      REAL(8), ALLOCATABLE :: qViscosity(:)  ! nondimensional Viscosity  in flux points
      REAL(8), ALLOCATABLE :: GradVisc(:,:)   ! gradient of nondimensional Viscosity
      REAL(8), ALLOCATABLE :: cvte(:)  ! nondimensional beta=beta^dim / beta^inRayleigh
      REAL(8), ALLOCATABLE :: DvaEpsGVisc(:,:) !  2\epika\cdot\nabla\mu
      REAL(8), ALLOCATABLE :: CuOmCrGrVi(:,:)  !  \nabla\time\omega\times\nabla\mu
!     Functions
      REAL GetFMmem

      REAL cptime
      REAL rcpu

      INTEGER i,j,itim,istop
      INTEGER eqn,ierr
!
!____ Init measuring cpu time
!
      cpu%t0=cptime(0.)
      cpu%t00=cpu%t0
      CALL InitCpu(cpu,itim)
! 
!____ Version data 
!
      CALL DatumInUra(inp%StartTime)
      inp%IDname='TriTok'
      inp%IDversion='8.0'
      inp%IDdate='Februar 2025'
      inp%INTversion=1
      inp%RSTversion=4
      inp%cTreeFileVers=2
!
!____ Init random generator
!
      CALL RANDOM_SEED
!
!____ MPI init
!
      CALL par_init(env)
!
!____ Input / output files
!
      CALL SetUpIO(env,io,inp)
! 
!     Read input file
!
      CALL ReadInputFile(env,mesh,io,inp,gauss)
!
!     Set up Gauss quadrature positions and weights
!
      CALL GaussPosWeights(gauss)
!
!     Prebere mrezo 
!
      CALL ReadMeshGid3D(env,io,inp,mesh,gauss)
      CALL VmesniCas(cpu,itim)
!
!     Generira tocke xq in povezave ibf za nezvezni del mreze
!
      CALL GenMakroMesh(mesh)         
      
      IF (env%myproc.EQ.1) THEN
      WRITE(io%l,'(A,I8)') "nnodes      =",mesh%nnodes  
      WRITE(io%l,'(A,I8)') "nq          =",mesh%nq
      WRITE(io%l,'(A,I8)') "vozlis!     =",mesh%nq+mesh%nnodes 
      WRITE(io%l,'(A)') ""
      WRITE(io%l,'(A,I8)') "nbnodes     =",mesh%nbnodes
      WRITE(io%l,'(A,I8)') "flux na robu=",mesh%nbelem*4
      WRITE(io%l,'(A,I8)') "robnih      =",mesh%nbelem*4+mesh%nbnodes
      WRITE(io%l,'(A)') ""
      END IF
      CALL VmesniCas(cpu,itim)

! 
!     Allocate field functions
!
      ALLOCATE (velocity(mesh%nnodes,3))
      ALLOCATE (vorticity(mesh%nnodes,3))
      ALLOCATE (qvorticity(mesh%nq,3))        
      ALLOCATE (pnlsvelocity(mesh%nnodes,3)) ! previous non linear iteration
      ALLOCATE (ptsvorticity(mesh%nnodes,3))  ! previous time step
      ALLOCATE (pptsvorticity(mesh%nnodes,3)) ! pred previous time step               
      ALLOCATE (pnlsvorticity(mesh%nnodes,3)) ! previous nonlinear iteration
      ALLOCATE (rhsv(mesh%nnodes,3))   ! body force (includes time derivatives)          

      ALLOCATE (temp(mesh%nnodes))
      ALLOCATE (qtemp(mesh%nq))        
      ALLOCATE (ptstemp(mesh%nnodes))  ! previous time step
      ALLOCATE (pptstemp(mesh%nnodes)) ! pred previous time step               
      ALLOCATE (pnlstemp(mesh%nnodes)) ! previous nonlinear iteration
      ALLOCATE (Trhsv(mesh%nnodes))    ! body force (total)
      ALLOCATE (Ttime(mesh%nnodes))    ! body force (only time derivatives)

      ALLOCATE (dc(mesh%nnodes))
      ALLOCATE (qdc(mesh%nq))        
      ALLOCATE (ptsdc(mesh%nnodes))  ! previous time step
      ALLOCATE (pptsdc(mesh%nnodes)) ! pred previous time step               
      ALLOCATE (pnlsdc(mesh%nnodes)) ! previous nonlinear iteration
      ALLOCATE (DCrhsv(mesh%nnodes))    ! body force (time derivatives + sources)
      ALLOCATE (DCtime(mesh%nnodes))    ! body force (time derivatives)
      ALLOCATE (DCdiff(mesh%nnodes))    ! difuzivnost
      ALLOCATE (DCqDiff(mesh%nq))    ! difuzivnost v macro tockah
      ALLOCATE (DCdiffGrad(mesh%nnodes,3))    ! gradient difuzivnosti 
      ALLOCATE (DCdiffLap(mesh%nnodes))    ! nabla kvadrat difuzivnosti
      ALLOCATE (DCanal(mesh%nnodes))    ! analytic solution
      ALLOCATE (DCsource(mesh%nnodes))    ! source


      ALLOCATE (mHu(mesh%nnodes)) ! modified Helmholtz
      ALLOCATE (mHq(mesh%nq))     ! equation
      ALLOCATE (fmqmH(mesh%nbnpof))     ! full matrix boundary flux
      ALLOCATE (mHrhsV(mesh%nnodes)) ! modified Helmholtz right hand side source
      ALLOCATE (mHtime(mesh%nnodes))    ! body force (time derivatives)
      ALLOCATE (mHdiff(mesh%nnodes))    ! difuzivnost
      ALLOCATE (mHqDiff(mesh%nq))    ! difuzivnost v macro tockah
      ALLOCATE (mHdiffGrad(mesh%nnodes,3))    ! gradient difuzivnosti
      ALLOCATE (ptsmHu(mesh%nnodes))  ! previous time step
      ALLOCATE (pptsmHu(mesh%nnodes)) ! pred previous time step               
      ALLOCATE (pnlsmHu(mesh%nnodes)) ! previous nonlinear iteration
      ALLOCATE (mHanal(mesh%nnodes)) ! previous nonlinear iteration


      ALLOCATE (tmp(mesh%nnodes))    ! zacasno pomozno polje

! 
!     Prebere BiC datoteko
!
      CALL ReadBic(env,io,inp,mesh,velocity,vorticity,qvorticity,temp,qtemp,dc,qdc,mHu,mhq,fmqmH)
      IF (inp%iacha.NE.0) CALL Channel3Danal(mesh,velocity)
      IF (env%myproc.EQ.1) CALL OutputInitial(io,inp,mesh,velocity,vorticity,temp,dc)

!
!     Set finite difference time scheme approximation
!
      CALL SetTimeScheme(env,io,inp)
!
!     V bic je koncna razdelitev vozlisc po stenah
!     za vsako robno vozlisce najdem najblizje vozlisce, ki
!     ni na isti steni. To potrebujem za dv/dn=0 robni pogoj
!
      CALL FillKobc(mesh)
      CALL VmesniCas(cpu,itim)
!
!     Divide boundary nodes across processors
!                          
      CALL SetEnvGlob(env,mesh)
      CALL DivideNodes(env,mesh)
!
!     Estimate requred RAM per processor
!
      IF (env%myproc.EQ.1) CALL EstRAM(env,mesh,io,inp)
      IF (inp%iestr.EQ.1) THEN
        CALL WarnErr(env,io,inp,0,'Tritok','Estimated RAM only, exiting!',0)
        CALL StopProgram(env,io,inp,0)
      END IF
!
!     Nanofluid
!
      IF (inp%inano.EQ.1) THEN
        CALL NanofluidProp(inp,io,env,nanoA,nanoB,nanoC)
        IF (inp%iDT.EQ.3) nanoC=1.0D0 ! ker racunam z variable properties
      ELSE
        nanoA=1.0D0
        nanoB=1.0D0
        nanoC=1.0D0
      END IF
!
!     Ferrofluid
!
      ALLOCATE (kbf(mesh%nnodes,3)) ! Kelvin body force
      IF (inp%iFerro.GT.0) THEN
        WRITE (io%l,'(A)') "Simulation of ferrofluid"
        WRITE (io%l,'(A,G18.9)') "Magnetic Rayleigh number =",inp%Ram
        WRITE (io%l,'(A,G18.9)') "Beta * Delta T =",inp%bdt
        WRITE (io%l,'(A,G18.9/)') "hi_0 =",inp%hi0

        ALLOCATE (mfs(mesh%nnodes,3)) ! magnetic field strength \vec H
        ALLOCATE (hshg(mesh%nnodes,3)) ! grad \vec H \cdot \vec H, gradient (H skalarno H)
        ALLOCATE (tempGrad(mesh%nnodes,3)) ! Temperature gradient
        CALL SetUpMagneticField(mesh,inp,mfs,hshg)
      ELSE
        kbf=0.0D0
        inp%ram=0.0D0
      END IF
!
!     Single Domain kinematika
!

!     Robne matrike
      IF (inp%iBw.NE.0) THEN
        ALLOCATE (sdH(env%zac:env%kon,mesh%nbnodes))  ! diffusion
        ALLOCATE (sdHtx(env%zac:env%kon,mesh%nbnodes)) ! doktorat, enacba (4.7)
        ALLOCATE (sdHty(env%zac:env%kon,mesh%nbnodes)) ! to je za prvi integral
        ALLOCATE (sdHtz(env%zac:env%kon,mesh%nbnodes)) ! na desni

!     Kvadratna interpolacija vseh funkcij
      IF (inp%iBw.EQ.1) THEN
!       Sistemske matrike in predpogojevalci
        ALLOCATE (sdAx(env%zac:env%kon,mesh%nbnodes))
        ALLOCATE (sdAy(env%zac:env%kon,mesh%nbnodes))
        ALLOCATE (sdAz(env%zac:env%kon,mesh%nbnodes))
        ALLOCATE (sdPivotx(env%zac:env%kon))
        ALLOCATE (sdPivoty(env%zac:env%kon))
        ALLOCATE (sdPivotz(env%zac:env%kon))
!       Obmocne matrike
        ALLOCATE (sdDx(env%zac:env%kon,mesh%nnodes)) ! doktorat, enacba (4.7)
        ALLOCATE (sdDy(env%zac:env%kon,mesh%nnodes)) ! to je za drugi integral
        ALLOCATE (sdDz(env%zac:env%kon,mesh%nnodes)) ! na desni
!
!       Compute or read single domain kinematics integrals
!
        CALL sdKm_corInt(env,io,inp,mesh,gauss,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz)
!
!       Form single domain kinematics system matrices and preconditioners
!
        CALL sdKm_sysm(1,env,mesh,sdAx,sdDx,sdDy,sdDz)
        CALL sdKm_sysm(2,env,mesh,sdAy,sdDx,sdDy,sdDz)
        CALL sdKm_sysm(3,env,mesh,sdAz,sdDx,sdDy,sdDz)
!        CALL CalculateRotationMatrices(mesh) ! test

!     Omega na desni strani linearno interpolirana
      ELSE IF (inp%iBw.EQ.2) THEN ! omega interpolirana linearno
!       Sistemske matrike in predpogojevalci
        ALLOCATE (sdAx(env%zac:env%kon,mesh%nbnodes))
        ALLOCATE (sdAy(env%zac:env%kon,mesh%nbnodes))
        ALLOCATE (sdAz(env%zac:env%kon,mesh%nbnodes))
        ALLOCATE (sdPivotx(env%zac:env%kon))
        ALLOCATE (sdPivoty(env%zac:env%kon))
        ALLOCATE (sdPivotz(env%zac:env%kon))
!       Obmocni integrali izracuznani z linearno interpolacijo za omego
        ALLOCATE (sdDx(env%zac:env%kon,mesh%nnodes8)) ! doktorat, enacba (4.7)
        ALLOCATE (sdDy(env%zac:env%kon,mesh%nnodes8)) ! to je za drugi integral
        ALLOCATE (sdDz(env%zac:env%kon,mesh%nnodes8)) ! na desni
!       robni del (kvadratna interpolacija)
        ALLOCATE (sdDxB(env%zac:env%kon,mesh%nbnodes)) ! doktorat, enacba (4.7)
        ALLOCATE (sdDyB(env%zac:env%kon,mesh%nbnodes)) ! to je za drugi integral
        ALLOCATE (sdDzB(env%zac:env%kon,mesh%nbnodes)) ! na desni

        CALL sdKm_corInt8(env,io,inp,mesh,gauss,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,sdDxB,sdDyB,sdDzB)
!
!       Form single domain kinematics system matrices
!
        CALL sdKm_sysm8(1,env,mesh,sdAx,sdDxB,sdDyB,sdDzB)
        CALL sdKm_sysm8(2,env,mesh,sdAy,sdDxB,sdDyB,sdDzB)
        CALL sdKm_sysm8(3,env,mesh,sdAz,sdDxB,sdDyB,sdDzB)

!     We recognise that normal vorticity component is zero at the wall,
!     so we solve only for tangential components
      ELSE IF (inp%iBw.EQ.3) THEN
!
!       For each source point, we find the largest component of the normal
!       This vorticity component, will be calculated by use of known w_normal
!
        CALL AnalyseSourcePointNormals(mesh)

!       Sistemske matrike in predpogojevalci
        ALLOCATE (sdAx(env%zac:env%kon,mesh%nbnodes))
        ALLOCATE (sdAy(env%zac:env%kon,mesh%nbnodes))
        ALLOCATE (sdAz(env%zac:env%kon,mesh%nbnodes))
        ALLOCATE (sdPivotx(env%zac:env%kon))
        ALLOCATE (sdPivoty(env%zac:env%kon))
        ALLOCATE (sdPivotz(env%zac:env%kon))
!       Obmocne matrike
        ALLOCATE (sdDx(env%zac:env%kon,mesh%nnodes)) ! doktorat, enacba (4.7)
        ALLOCATE (sdDy(env%zac:env%kon,mesh%nnodes)) ! to je za drugi integral
        ALLOCATE (sdDz(env%zac:env%kon,mesh%nnodes)) ! na desni
!
!       Compute or read single domain kinematics integrals
!
        CALL sdKm_corInt(env,io,inp,mesh,gauss,sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz)
!
!       Form single domain kinematics system matrices and preconditioners
!
        CALL sdKm_sysmBCN(1,env,mesh,sdAx,sdDx,sdDy,sdDz)
        CALL sdKm_sysmBCN(2,env,mesh,sdAy,sdDx,sdDy,sdDz)
        CALL sdKm_sysmBCN(3,env,mesh,sdAz,sdDx,sdDy,sdDz)

      END IF
!
!     Preconditioners
!

        IF (env%nproc.EQ.1) THEN
!         1 proc:  0 2 2 500 5 1.E-6
          CALL FormPREfm(2,mesh%nbnodes,sdPivotx,sdAx,rcpu,ierr)  ! pret=2 - LU faktorizacija
          CALL FormPREfm(2,mesh%nbnodes,sdPivoty,sdAy,rcpu,ierr)  ! pret=2 - LU faktorizacija
          CALL FormPREfm(2,mesh%nbnodes,sdPivotz,sdAz,rcpu,ierr)  ! pret=2 - LU faktorizacija
        ELSE
!          vec proc:  2 1 2 500 5 1.E-6
          CALL pFormPREfm(env,1,mesh%nbnodes,env%nmn,sdPivotx,sdAx,rcpu,ierr) ! pret=1 - diagonal
          CALL pFormPREfm(env,1,mesh%nbnodes,env%nmn,sdPivoty,sdAy,rcpu,ierr) ! pret=1 - diagonal
          CALL pFormPREfm(env,1,mesh%nbnodes,env%nmn,sdPivotz,sdAz,rcpu,ierr) ! pret=1 - diagonal
        END IF

      IF (env%myproc.EQ.1) THEN
        WRITE (io%l,'(A)') "single domain kinematika za robne vrtincnosti"
        WRITE (io%l,'(A,I7)')  "stevilo enacb   =",mesh%nbnodes
        WRITE (io%l,'(A,I7/)') "stevilo neznank =",mesh%nbnodes
        IF (inp%iBw.EQ.2) THEN
          WRITE (io%l,'(A)') "Uporabljam linearno aproksimacijo za vrtincnost na desni!"
          WRITE (io%l,'(A,I7/)')  "stevilo linearnih vozlisc   =",mesh%nnodes8
        END IF
      END IF
      CALL VmesniCas(cpu,itim)

      END IF


!
!     modified Helmholtz -  Set mu value
!
      IF (inp%imH.GT.0) THEN

!       set up problem
        IF (inp%imH.EQ.1) THEN
          DO i=1,mesh%nnodes
            IF (inp%mHf.EQ.1.OR.inp%mHf.EQ.2) THEN
              mHdiff(i)= inp%beta / ( inp%mHmu * inp%mHmu )
            ELSE IF (inp%mHf.EQ.3) THEN
              mHdiff(i)= inp%beta / ( (1+mesh%x(i,1)*mesh%x(i,2)*mesh%x(i,3))**2 )
            ELSE
              CALL WarnErr(env,io,inp,4,"TriTok","Wrong mod.Helm test case!",inp%mHf)
              STOP
            END IF
          END DO
        ELSE
          CALL modHelm_SetProblem(inp,mesh,mHdiff,mHrhsV,mHtime,mHu,mHq,fmqmH,cnt%rtime,velocity,mHanal)
        END IF

!       calculate difusivitiy gradients
        CALL setGrad(mesh,mHdiff,mHdiffGrad)
!       interpolate diffusivity to macro mesh
        CALL Int2MacroNodes(mesh,mHdiff,mHqDiff)

!       single domain solutions
        IF (inp%imH.GE.3) THEN

          IF (env%myproc.EQ.1) THEN
            WRITE (io%l,'(A)') "single domain kinetika za modificirano Helmholtzovo enacbo"
            WRITE (io%l,'(A,F13.3,A)') "full matrix memory estimate =",GetFMmem(mesh)," Gb."
            CALL FLUSH(io%l)
          END IF


!
!         Single domain KINETIKA za modified Helmholtz
!
          ALLOCATE (UmatH(mesh%nnodes,mesh%nbnodes))   ! grad u* skalarno normala po gamma
          ALLOCATE (UmatG(mesh%nnodes,mesh%nbnpof))  ! u* po gamma
          ALLOCATE (UmatB(mesh%nnodes,mesh%nnodes))   ! u* po omega
          ALLOCATE (UmatBp(mesh%nnodes,mesh%nnodes))   ! u* po omega za MATMUL pri var diff
          ALLOCATE (UmatAbdx(mesh%nnodes,mesh%nbnodes)) ! u* krat nx po gamma
          ALLOCATE (UmatAbdy(mesh%nnodes,mesh%nbnodes)) ! u* krat ny po gamma
          ALLOCATE (UmatAbdz(mesh%nnodes,mesh%nbnodes)) ! u* krat nz po gamma
          ALLOCATE (UmatDx(mesh%nnodes,mesh%nnodes)) ! grad u*_x po omega
          ALLOCATE (UmatDy(mesh%nnodes,mesh%nnodes)) ! grad u*_y po omega
          ALLOCATE (UmatDz(mesh%nnodes,mesh%nnodes)) ! grad u*_z po omega

          ALLOCATE (QmatH(mesh%nbnpof,mesh%nbnodes))   ! grad u* skalarno normala po gamma
          ALLOCATE (QmatG(mesh%nbnpof,mesh%nbnpof))  ! u* po gamma
          ALLOCATE (QmatB(mesh%nbnpof,mesh%nnodes))   ! u* po omega
          ALLOCATE (QmatBp(mesh%nbnpof,mesh%nnodes))   ! u* po omega ! u* po omega za MATMUL pri var diff
          ALLOCATE (QmatAbdx(mesh%nbnpof,mesh%nbnodes)) ! u* krat nx po gamma
          ALLOCATE (QmatAbdy(mesh%nbnpof,mesh%nbnodes)) ! u* krat ny po gamma
          ALLOCATE (QmatAbdz(mesh%nbnpof,mesh%nbnodes)) ! u* krat nz po gamma
          ALLOCATE (QmatDx(mesh%nbnpof,mesh%nnodes)) ! grad u*_x po omega
          ALLOCATE (QmatDy(mesh%nbnpof,mesh%nnodes)) ! grad u*_y po omega
          ALLOCATE (QmatDz(mesh%nbnpof,mesh%nnodes)) ! grad u*_z po omega

          CALL corMin_SDkin_mHfs(env,io,inp,mesh,gauss,mHdiff,mHqDiff, &
                     UmatH,UmatG,UmatB,UmatBp,UmatAbdx,UmatAbdy,UmatAbdz,UmatDx,UmatDy,UmatDz, &
                     QmatH,QmatG,QmatB,QmatBp,QmatAbdx,QmatAbdy,QmatAbdz,QmatDx,QmatDy,QmatDz)

!         mesh%mHnb     ! ta pove koliko je znanih u ali q
!         mesh%mHfmnunk ! ta pove koliko je neznanih u ali q
          ALLOCATE (fmmHA(mesh%mHfmnunk,mesh%mHfmnunk))
          ALLOCATE (fmmHPivot(mesh%mHfmnunk))
          ALLOCATE (fmmHRHS(mesh%mHfmnunk,mesh%mHnb))

!         set up system matrix
          CALL fmmHSysM(mesh,fmmHA,fmmHRHS,mHdiff,mHqDiff,UmatH,UmatG,QmatH,QmatG)
          DEALLOCATE (UmatH,UmatG,QmatH,QmatG)
!         calculate preconditioner
          CALL FormPREfm(2,mesh%mHfmnunk,fmmHPivot,fmmHA,rcpu,ierr)  ! pret=2 - LU faktorizacija

!
!         wavelet compression
!
          IF (inp%imH.EQ.5) THEN

            WRITE(io%l,'(A,G10.5)') "wavelet compression, kappa =",inp%WT_kappa
            WRITE(io%l,'(A)') "matrix,MB,comp,spar,norm"
            WTnnz=0.0D0
            WTn=0.0D0
            CALL WaveletCompSmat(io,"UmatB   ",UmatB,mesh%nnodes,mesh%nnodes,cUmatB,WTnnz,WTn,inp%WT_kappa)
            DEALLOCATE(UmatB)
            CALL WaveletCompSmat(io,"UmatBp  ",UmatBp,mesh%nnodes,mesh%nnodes,cUmatBp,WTnnz,WTn,inp%WT_kappa)   ! u* po omega
            DEALLOCATE(UmatBp)
            CALL WaveletCompSmat(io,"UmatAbdx",UmatAbdx,mesh%nnodes,mesh%nbnodes,cUmatAbdx,WTnnz,WTn,inp%WT_kappa) ! u* krat nx po gamma
            DEALLOCATE(UmatAbdx)
            CALL WaveletCompSmat(io,"UmatAbdy",UmatAbdy,mesh%nnodes,mesh%nbnodes,cUmatAbdy,WTnnz,WTn,inp%WT_kappa) ! u* krat ny po gamma
            DEALLOCATE(UmatAbdy)
            CALL WaveletCompSmat(io,"UmatAbdz",UmatAbdz,mesh%nnodes,mesh%nbnodes,cUmatAbdz,WTnnz,WTn,inp%WT_kappa) ! u* krat nz po gamma
            DEALLOCATE(UmatAbdz)
            CALL WaveletCompSmat(io,"UmatDx  ",UmatDx,mesh%nnodes,mesh%nnodes,cUmatDx,WTnnz,WTn,inp%WT_kappa) ! grad u*_x po omega
            DEALLOCATE(UmatDx)
            CALL WaveletCompSmat(io,"UmatDy  ",UmatDy,mesh%nnodes,mesh%nnodes,cUmatDy,WTnnz,WTn,inp%WT_kappa) ! grad u*_y po omega
            DEALLOCATE(UmatDy)
            CALL WaveletCompSmat(io,"UmatDz  ",UmatDz,mesh%nnodes,mesh%nnodes,cUmatDz,WTnnz,WTn,inp%WT_kappa) ! grad u*_z po omega
            DEALLOCATE(UmatDz)

            CALL WaveletCompSmat(io,"QmatB   ",QmatB,mesh%nbnpof,mesh%nnodes,cQmatB,WTnnz,WTn,inp%WT_kappa)   ! u* po omega
            DEALLOCATE(QmatB)
            CALL WaveletCompSmat(io,"QmatBp  ",QmatBp,mesh%nbnpof,mesh%nnodes,cQmatBp,WTnnz,WTn,inp%WT_kappa)   ! u* po omega
            DEALLOCATE(QmatBp)
            CALL WaveletCompSmat(io,"QmatAbdx",QmatAbdx,mesh%nbnpof,mesh%nbnodes,cQmatAbdx,WTnnz,WTn,inp%WT_kappa) ! u* krat nx po gamma
            DEALLOCATE(QmatAbdx)
            CALL WaveletCompSmat(io,"QmatAbdy",QmatAbdy,mesh%nbnpof,mesh%nbnodes,cQmatAbdy,WTnnz,WTn,inp%WT_kappa) ! u* krat ny po gamma
            DEALLOCATE(QmatAbdy)
            CALL WaveletCompSmat(io,"QmatAbdz",QmatAbdz,mesh%nbnpof,mesh%nbnodes,cQmatAbdz,WTnnz,WTn,inp%WT_kappa) ! u* krat nz po gamma
            DEALLOCATE(QmatAbdz)
            CALL WaveletCompSmat(io,"QmatDx  ",QmatDx,mesh%nbnpof,mesh%nnodes,cQmatDx,WTnnz,WTn,inp%WT_kappa) ! grad u*_x po omega
            DEALLOCATE(QmatDx)
            CALL WaveletCompSmat(io,"QmatDy  ",QmatDy,mesh%nbnpof,mesh%nnodes,cQmatDy,WTnnz,WTn,inp%WT_kappa) ! grad u*_y po omega
            DEALLOCATE(QmatDy)
            CALL WaveletCompSmat(io,"QmatDz  ",QmatDz,mesh%nbnpof,mesh%nnodes,cQmatDz,WTnnz,WTn,inp%WT_kappa) ! grad u*_z po omega
            DEALLOCATE(QmatDz)

            WRITE(io%l,'(A,G20.10)') "TOTAL matrix element ratio (compressed / full) =  ",WTnnz/WTn
            FLUSH(io%l)


!            CALL WaveletComp(env,io,inp,mesh,
!     &               UmatB,UmatBp,UmatAbdx,UmatAbdy,UmatAbdz,UmatDx,UmatDy,UmatDz,
!     &               QmatB,QmatBp,QmatAbdx,QmatAbdy,QmatAbdz,QmatDx,QmatDy,QmatDz,
!     &               cUmatB,cUmatBp,cUmatAbdx,cUmatAbdy,cUmatAbdz,cUmatDx,cUmatDy,cUmatDz,
!     &               cQmatB,cQmatBp,cQmatAbdx,cQmatAbdy,cQmatAbdz,cQmatDx,cQmatDy,cQmatDz)

!            DEALLOCATE (UmatB,UmatBp,UmatAbdx,UmatAbdy,UmatAbdz,UmatDx,UmatDy,UmatDz,
!     &                  QmatB,QmatBp,QmatAbdx,QmatAbdy,QmatAbdz,QmatDx,QmatDy,QmatDz)


          END IF

!
!         ACA compression
!

          IF (inp%imH.EQ.6) THEN

            WRITE(io%l,'(A,G10.5)') "ACA compression, eps =",inp%ACA_eps
            WRITE(io%l,'(A)') "matrix,MB,comp,norm"
            WTnnz=0.0D0
            WTn=0.0D0

            CALL AcaCompSmat(io,"UmatB   ",UmatB,mesh%nnodes,mesh%nnodes,aUmatB,WTnnz,WTn,inp%ACA_eps,inp%ACA_type)
            DEALLOCATE(UmatB)
            CALL AcaCompSmat(io,"UmatBp  ",UmatBp,mesh%nnodes,mesh%nnodes,aUmatBp,WTnnz,WTn,inp%ACA_eps,inp%ACA_type)   ! u* po omega
            DEALLOCATE(UmatBp)
            CALL AcaCompSmat(io,"UmatAbdx",UmatAbdx,mesh%nnodes,mesh%nbnodes,aUmatAbdx,WTnnz,WTn,inp%ACA_eps,inp%ACA_type) ! u* krat nx po gamma
            DEALLOCATE(UmatAbdx)
            CALL AcaCompSmat(io,"UmatAbdy",UmatAbdy,mesh%nnodes,mesh%nbnodes,aUmatAbdy,WTnnz,WTn,inp%ACA_eps,inp%ACA_type) ! u* krat ny po gamma
            DEALLOCATE(UmatAbdy)
            CALL AcaCompSmat(io,"UmatAbdz",UmatAbdz,mesh%nnodes,mesh%nbnodes,aUmatAbdz,WTnnz,WTn,inp%ACA_eps,inp%ACA_type) ! u* krat nz po gamma
            DEALLOCATE(UmatAbdz)
            CALL AcaCompSmat(io,"UmatDx  ",UmatDx,mesh%nnodes,mesh%nnodes,aUmatDx,WTnnz,WTn,inp%ACA_eps,inp%ACA_type) ! grad u*_x po omega
            DEALLOCATE(UmatDx)
            CALL AcaCompSmat(io,"UmatDy  ",UmatDy,mesh%nnodes,mesh%nnodes,aUmatDy,WTnnz,WTn,inp%ACA_eps,inp%ACA_type) ! grad u*_y po omega
            DEALLOCATE(UmatDy)
            CALL AcaCompSmat(io,"UmatDz  ",UmatDz,mesh%nnodes,mesh%nnodes,aUmatDz,WTnnz,WTn,inp%ACA_eps,inp%ACA_type) ! grad u*_z po omega
            DEALLOCATE(UmatDz)

            CALL AcaCompSmat(io,"QmatB   ",QmatB,mesh%nbnpof,mesh%nnodes,aQmatB,WTnnz,WTn,inp%ACA_eps,inp%ACA_type)   ! u* po omega
            DEALLOCATE(QmatB)
            CALL AcaCompSmat(io,"QmatBp  ",QmatBp,mesh%nbnpof,mesh%nnodes,aQmatBp,WTnnz,WTn,inp%ACA_eps,inp%ACA_type)   ! u* po omega
            DEALLOCATE(QmatBp)
            CALL AcaCompSmat(io,"QmatAbdx",QmatAbdx,mesh%nbnpof,mesh%nbnodes,aQmatAbdx,WTnnz,WTn,inp%ACA_eps,inp%ACA_type) ! u* krat nx po gamma
            DEALLOCATE(QmatAbdx)
            CALL AcaCompSmat(io,"QmatAbdy",QmatAbdy,mesh%nbnpof,mesh%nbnodes,aQmatAbdy,WTnnz,WTn,inp%ACA_eps,inp%ACA_type) ! u* krat ny po gamma
            DEALLOCATE(QmatAbdy)
            CALL AcaCompSmat(io,"QmatAbdz",QmatAbdz,mesh%nbnpof,mesh%nbnodes,aQmatAbdz,WTnnz,WTn,inp%ACA_eps,inp%ACA_type) ! u* krat nz po gamma
            DEALLOCATE(QmatAbdz)
            CALL AcaCompSmat(io,"QmatDx  ",QmatDx,mesh%nbnpof,mesh%nnodes,aQmatDx,WTnnz,WTn,inp%ACA_eps,inp%ACA_type) ! grad u*_x po omega
            DEALLOCATE(QmatDx)
            CALL AcaCompSmat(io,"QmatDy  ",QmatDy,mesh%nbnpof,mesh%nnodes,aQmatDy,WTnnz,WTn,inp%ACA_eps,inp%ACA_type) ! grad u*_y po omega
            DEALLOCATE(QmatDy)
            CALL AcaCompSmat(io,"QmatDz  ",QmatDz,mesh%nbnpof,mesh%nnodes,aQmatDz,WTnnz,WTn,inp%ACA_eps,inp%ACA_type) ! grad u*_z po omega
            DEALLOCATE(QmatDz)



            WRITE(io%l,'(A,G20.10)') "TOTAL matrix element ratio (compressed / full) =  ",WTnnz/WTn
            FLUSH(io%l)

          END IF

!
!         Cluster tree + compression
!

          IF (inp%imH.EQ.7) THEN

!           Prepare H matrix structure based on admisibility criterion
            CALL SetUpHmatrices(env,io,inp,mesh, &
                    hUmatB,hUmatBp,hUmatAbdx,hUmatAbdy,hUmatAbdz,hUmatDx,hUmatDy,hUmatDz, &
                    hQmatB,hQmatBp,hQmatAbdx,hQmatAbdy,hQmatAbdz,hQmatDx,hQmatDy,hQmatDz)


            IF (inp%ct_type.EQ.2) THEN
              WRITE(io%l,'(A,G10.5)') "wavelet compression, kappa =",inp%WT_kappa
              WRITE(io%l,'(A)') "matrix,MB,comp,norm,nMparts"
            ELSE IF (inp%ct_type.EQ.3) THEN
              WRITE(io%l,'(A,G10.5)') "ACA compression, eps =",inp%ACA_eps
              WRITE(io%l,'(A)') "matrix,MB,comp,norm,nMparts"
            ELSE
              WRITE(io%l,'(A)') "No cluster tree compression!"
            END IF
            WTnnz=0.0D0
            WTn=0.0D0

            CALL cTree_FillMatrices(io,inp,"UmatB   ",UmatB,mesh%nnodes,mesh%nnodes,hUmatB,WTnnz,WTn,inp%ct_type)
            DEALLOCATE(UmatB)
            CALL cTree_FillMatrices(io,inp,"UmatBp  ",UmatBp,mesh%nnodes,mesh%nnodes,hUmatBp,WTnnz,WTn,inp%ct_type)
            DEALLOCATE(UmatBp)
            CALL cTree_FillMatrices(io,inp,"UmatAbdx",UmatAbdx,mesh%nnodes,mesh%nbnodes,hUmatAbdx,WTnnz,WTn,inp%ct_type)
            DEALLOCATE(UmatAbdx)
            CALL cTree_FillMatrices(io,inp,"UmatAbdy",UmatAbdy,mesh%nnodes,mesh%nbnodes,hUmatAbdy,WTnnz,WTn,inp%ct_type)
            DEALLOCATE(UmatAbdy)
            CALL cTree_FillMatrices(io,inp,"UmatAbdz",UmatAbdz,mesh%nnodes,mesh%nbnodes,hUmatAbdz,WTnnz,WTn,inp%ct_type)
            DEALLOCATE(UmatAbdz)
            CALL cTree_FillMatrices(io,inp,"UmatDx  ",UmatDx,mesh%nnodes,mesh%nnodes,hUmatDx,WTnnz,WTn,inp%ct_type)
            DEALLOCATE(UmatDx)
            CALL cTree_FillMatrices(io,inp,"UmatDy  ",UmatDy,mesh%nnodes,mesh%nnodes,hUmatDy,WTnnz,WTn,inp%ct_type)
            DEALLOCATE(UmatDy)
            CALL cTree_FillMatrices(io,inp,"UmatDz  ",UmatDz,mesh%nnodes,mesh%nnodes,hUmatDz,WTnnz,WTn,inp%ct_type)
            DEALLOCATE(UmatDz)
            CALL cTree_FillMatrices(io,inp,"QmatB   ",QmatB,mesh%nbnpof,mesh%nnodes,hQmatB,WTnnz,WTn,inp%ct_type)
            DEALLOCATE(QmatB)
            CALL cTree_FillMatrices(io,inp,"QmatBp  ",QmatBp,mesh%nbnpof,mesh%nnodes,hQmatBp,WTnnz,WTn,inp%ct_type)
            DEALLOCATE(QmatBp)
            CALL cTree_FillMatrices(io,inp,"QmatAbdx",QmatAbdx,mesh%nbnpof,mesh%nbnodes,hQmatAbdx,WTnnz,WTn,inp%ct_type)
            DEALLOCATE(QmatAbdx)
            CALL cTree_FillMatrices(io,inp,"QmatAbdy",QmatAbdy,mesh%nbnpof,mesh%nbnodes,hQmatAbdy,WTnnz,WTn,inp%ct_type)
            DEALLOCATE(QmatAbdy)
            CALL cTree_FillMatrices(io,inp,"QmatAbdz",QmatAbdz,mesh%nbnpof,mesh%nbnodes,hQmatAbdz,WTnnz,WTn,inp%ct_type)
            DEALLOCATE(QmatAbdz)
            CALL cTree_FillMatrices(io,inp,"QmatDx  ",QmatDx,mesh%nbnpof,mesh%nnodes,hQmatDx,WTnnz,WTn,inp%ct_type)
            DEALLOCATE(QmatDx)
            CALL cTree_FillMatrices(io,inp,"QmatDy  ",QmatDy,mesh%nbnpof,mesh%nnodes,hQmatDy,WTnnz,WTn,inp%ct_type)
            DEALLOCATE(QmatDy)
            CALL cTree_FillMatrices(io,inp,"QmatDz  ",QmatDz,mesh%nbnpof,mesh%nnodes,hQmatDz,WTnnz,WTn,inp%ct_type)
            DEALLOCATE(QmatDz)

            WRITE(io%l,'(A,G20.10)') "TOTAL matrix element ratio (compressed / full) =  ",WTnnz/WTn
            FLUSH(io%l)

          END IF

        END IF
      END IF

      CALL VmesniCas(cpu,itim)        

!
!     Macro
!

!
!     Determine kinematics knowns and unknowns
!
      CALL DetKMunk3(mesh,kmp)
       
!
!     Integracija, shranim v pravokotne matrike
!      
      ALLOCATE (smatH(mesh%nicnsp,mesh%npoc),smatG(mesh%nicnsp,mesh%npofc))
!
!     Kinetics matrices
!      
      ALLOCATE (smatB(mesh%nicnsp,mesh%npoc))   ! source term
      ALLOCATE (smatAbdx(mesh%nicnsp,mesh%npoc)) ! advection - boundary  part X
      ALLOCATE (smatAbdy(mesh%nicnsp,mesh%npoc)) ! advection - boundary  part Y
      ALLOCATE (smatAbdz(mesh%nicnsp,mesh%npoc)) ! advection - boundary  part Z
!
!     kinematics matrices
!            
      ALLOCATE (smatHtx(mesh%nicnpoc,mesh%npoc)) ! doktorat, enacba (4.7)
      ALLOCATE (smatHty(mesh%nicnpoc,mesh%npoc)) ! to je za prvi integral 
      ALLOCATE (smatHtz(mesh%nicnpoc,mesh%npoc)) ! na desni
      
      ALLOCATE (smatDx(mesh%nicnsp,mesh%npoc)) ! doktorat, enacba (4.7)
      ALLOCATE (smatDy(mesh%nicnsp,mesh%npoc)) ! to je za drugi integral 
      ALLOCATE (smatDz(mesh%nicnsp,mesh%npoc)) ! na desni            
!
!     Integration
!
      IF ( (inp%iDV.NE.0) .OR. (inp%iDW.NE.0) .OR. (inp%iDT.NE.0) .OR. (inp%iDC.NE.0) ) THEN
        CALL corMin_kmdcSteady(env,io,inp,mesh,gauss,smatH,smatG,smatB, &
                                smatAbdx,smatAbdy,smatAbdz, &
                                smatHtx,smatHty,smatHtz,smatDx,smatDy,smatDz)          
      END IF

!
!     modified Helmholzt
!      
      IF (inp%imH.EQ.1.OR.inp%imH.EQ.2) THEN
        ALLOCATE (mHsmatH(mesh%nicnsp,mesh%npoc))  ! H matrix
        ALLOCATE (mHsmatG(mesh%nicnsp,mesh%npofc)) ! G matrix
        ALLOCATE (mHsmatB(mesh%nicnsp,mesh%npoc))  ! source term
        ALLOCATE (mHsmatAbdx(mesh%nicnsp,mesh%npoc)) ! advection - boundary  part X
        ALLOCATE (mHsmatAbdy(mesh%nicnsp,mesh%npoc)) ! advection - boundary  part Y
        ALLOCATE (mHsmatAbdz(mesh%nicnsp,mesh%npoc)) ! advection - boundary  part Z
!
!       kinematics matrices
!
        ALLOCATE (mHsmatHtx(mesh%nicnpoc,mesh%npoc)) ! doktorat, enacba (4.7)
        ALLOCATE (mHsmatHty(mesh%nicnpoc,mesh%npoc)) ! to je za prvi integral
        ALLOCATE (mHsmatHtz(mesh%nicnpoc,mesh%npoc)) ! na desni

        ALLOCATE (mHsmatDx(mesh%nicnsp,mesh%npoc)) ! doktorat, enacba (4.7)
        ALLOCATE (mHsmatDy(mesh%nicnsp,mesh%npoc)) ! to je za drugi integral
        ALLOCATE (mHsmatDz(mesh%nicnsp,mesh%npoc)) ! na desni
!
!       Integration
!
        CALL corMin_modHelm(env,io,inp,mesh,gauss,mHdiff,mHqDiff,mHsmatH,mHsmatG,mHsmatB, &
                                mHsmatAbdx,mHsmatAbdy,mHsmatAbdz, &
                                mHsmatHtx,mHsmatHty,mHsmatHtz,mHsmatDx,mHsmatDy,mHsmatDz)

      END IF

      CALL VmesniCas(cpu,itim)
!
!     tvorim CRS sysm in rhsm za kinetiko (notranje vrtincnosti)
!
      IF (inp%iDw.EQ.1) THEN  ! predolocen sistem enacb
        CALL sMat2crsSysRhsB_count(1,mesh,VTEsysmX,VTErhsmX)
        CALL sMat2crsSysRhsB_count(2,mesh,VTEsysmY,VTErhsmY)
        CALL sMat2crsSysRhsB_count(3,mesh,VTEsysmZ,VTErhsmZ)
!       Set up preconditioning vector (for diagonal preconditioning)
        ALLOCATE(VTEprecvX(mesh%nunk(1)))
        ALLOCATE(VTEprecvY(mesh%nunk(2)))
        ALLOCATE(VTEprecvZ(mesh%nunk(3)))            

        CALL sMat2crsSysRhsB_w_fill(1,mesh,smatH,smatG,smatB,inp%beta,inp%ren/nanoA,VTEsysmX,VTErhsmX, &
                                      velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)
        CALL sMat2crsSysRhsB_w_fill(2,mesh,smatH,smatG,smatB,inp%beta,inp%ren/nanoA,VTEsysmY,VTErhsmY, &
                                      velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)
        CALL sMat2crsSysRhsB_w_fill(3,mesh,smatH,smatG,smatB,inp%beta,inp%ren/nanoA,VTEsysmZ,VTErhsmZ, &
                                      velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)

!       calculate preconditioner        
        CALL prelsqr2(VTEsysmX%neq,mesh%nunk(1),VTEsysmX%nnz,VTEprecvX,VTEsysmX%i,VTEsysmX%j,VTEsysmX%v)      
        CALL prelsqr2(VTEsysmY%neq,mesh%nunk(2),VTEsysmY%nnz,VTEprecvY,VTEsysmY%i,VTEsysmY%j,VTEsysmY%v)      
        CALL prelsqr2(VTEsysmZ%neq,mesh%nunk(3),VTEsysmZ%nnz,VTEprecvZ,VTEsysmZ%i,VTEsysmZ%j,VTEsysmZ%v)



      ELSE IF (inp%iDw.EQ.2) THEN  ! dolocen sistem enacb
        CALL sMat2crsSysRhsB_SQsetup(1,mesh)
        CALL sMat2crsSysRhsB_SQsetup(2,mesh)
        CALL sMat2crsSysRhsB_SQsetup(3,mesh)

        CALL sMat2crsSysRhsB_SQcount(1,mesh,VTEsysmX,VTErhsmX)
        CALL sMat2crsSysRhsB_SQcount(2,mesh,VTEsysmY,VTErhsmY)
        CALL sMat2crsSysRhsB_SQcount(3,mesh,VTEsysmZ,VTErhsmZ)

        CALL sMat2crsSysRhsB_w_SQfill(1,mesh,smatH,smatG,smatB,inp%beta,inp%ren/nanoA,VTEsysmX,VTErhsmX, &
                                      velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)
        CALL sMat2crsSysRhsB_w_SQfill(2,mesh,smatH,smatG,smatB,inp%beta,inp%ren/nanoA,VTEsysmY,VTErhsmY, &
                                      velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)
        CALL sMat2crsSysRhsB_w_SQfill(3,mesh,smatH,smatG,smatB,inp%beta,inp%ren/nanoA,VTEsysmZ,VTErhsmZ, &
                                      velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)

!       Set up preconditioning matrix : 2 - ilu, 1 - diag, 0 - none
        CALL CopyMatrix(VTEsysmX,VTEprecX)
        CALL FormPRMjr(inp%sqrs_prec,VTEprecX,VTEsysmX)

        CALL CopyMatrix(VTEsysmY,VTEprecY)
        CALL FormPRMjr(inp%sqrs_prec,VTEprecY,VTEsysmY)

        CALL CopyMatrix(VTEsysmZ,VTEprecZ)
        CALL FormPRMjr(inp%sqrs_prec,VTEprecZ,VTEsysmZ)

      ELSE IF (inp%iDw.EQ.3) THEN  ! variable viscosity and beta
!       Set up systems of equations
        CALL sMat2crsSysRhsB_count(1,mesh,VTEsysmX,VTErhsmX)
        CALL sMat2crsSysRhsB_count(2,mesh,VTEsysmY,VTErhsmY)
        CALL sMat2crsSysRhsB_count(3,mesh,VTEsysmZ,VTErhsmZ)
!       Set up preconditioning vector (for diagonal preconditioning)
        ALLOCATE(VTEprecvX(mesh%nunk(1)))
        ALLOCATE(VTEprecvY(mesh%nunk(2)))
        ALLOCATE(VTEprecvZ(mesh%nunk(3)))
!       Set up variable material properties
        ALLOCATE (Viscosity(mesh%nnodes))  ! nondimensional mu=mu^dim / mu^inReynolds
        ALLOCATE (qViscosity(mesh%nq))  ! nondimensional mu=mu^dim / mu^inReynolds (v flux tockah)
        ALLOCATE (GradVisc(mesh%nnodes,3)) ! gradient of viscosity
        ALLOCATE (cvte(mesh%nnodes)) ! nondimensional beta=beta^dim / beta^inRayleigh
        ALLOCATE( DvaEpsGVisc(mesh%nnodes,3)) !  2\epika\cdot\nabla\mu
        ALLOCATE( CuOmCrGrVi(mesh%nnodes,3) ) !  \nabla\time\omega\times\nabla\mu

        CALL SetUpViscCvte(env,io,inp,mesh,viscosity,cvte,temp)

      ELSE IF (inp%iDw.NE.0) THEN
        CALL WarnErr(env,io,inp,4,"Tritok","Wrong iDw code!",inp%iDw)
      END IF 
      IF (inp%iDw.GT.0) THEN
        IF (env%myproc.EQ.1) THEN
          WRITE (io%l,'(A)') "makro kinetika za notranje vrtincnosti"
          WRITE (io%l,'(A,I7)')  "stevilo enacb   X =",VTEsysmX%neq
          WRITE (io%l,'(A,I7)') "stevilo neznank X =",mesh%nunk(1)
          WRITE (io%l,'(A,I7)')  "stevilo enacb   Y =",VTEsysmY%neq
          WRITE (io%l,'(A,I7)') "stevilo neznank Y =",mesh%nunk(2)
          WRITE (io%l,'(A,I7)')  "stevilo enacb   Z =",VTEsysmZ%neq
          WRITE (io%l,'(A,I7/)') "stevilo neznank Z =",mesh%nunk(3)
        END IF
      END IF
!
!     Set up particles
!
      IF (inp%iPart.EQ.1) THEN
        CALL CoR_NeighbourList(env,mesh,io,inp)
        CALL CalNodeVolume(mesh)
        ALLOCATE (part(inp%npart))
        CALL SetUpParticles(inp,part,mesh)
        ALLOCATE (PartVolFrac(mesh%nnodes))
        CALL CalParticleVolumeFraction(env,io,inp,gauss,part,mesh,PartVolFrac)
        IF (inp%iBrown.EQ.1) CALL SetUpBrown(env,io,inp)
        IF (inp%iThph.EQ.1)  CALL SetUpThermophoresis(env,io,inp)
      END IF


!
!     tvorim CRS sysm in rhsm za kinetiko (notranje temperature)
!
      IF (inp%iDT.EQ.1) THEN ! predolocen sistem enacb
        CALL sMat2crsSysRhsB_count(4,mesh,Tsysm,Trhsm)
!       Set up preconditioning vector (for diagonal preconditioning)
        ALLOCATE(Tprecv(mesh%Tnunk))
        CALL sMat2crsSysRhsB_Temp_fill(mesh,smatH,smatG,smatB,inp%beta,inp%ren*inp%prn/nanoC,Tsysm,Trhsm, &
                            velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)              
!       preconditioner            
        CALL prelsqr2(Tsysm%neq,mesh%Tnunk,Tsysm%nnz,Tprecv,Tsysm%i,Tsysm%j,Tsysm%v)
      ELSE IF (inp%iDT.EQ.2) THEN ! dolocen sistem enacb
        CALL sMat2crsSysRhsB_SQsetup(4,mesh)
        CALL sMat2crsSysRhsB_SQcount(4,mesh,Tsysm,Trhsm)
        CALL sMat2crsSysRhsB_Temp_SQfill(mesh,smatH,smatG,smatB,inp%beta,inp%ren*inp%prn/nanoC,Tsysm,Trhsm, &
                            velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)

!       Set up preconditioning matrix : 2 - ilu, 1 - diag, 0 - none
        CALL CopyMatrix(Tsysm,Tprec)
        CALL FormPRMjr(inp%sqrs_prec,Tprec,Tsysm)
      ELSE IF (inp%iDT.EQ.3) THEN ! spremenljive lastnosti tekocine
!       Set up variable material properties
        ALLOCATE (Density(mesh%nnodes))  ! nondimensional density
        ALLOCATE (HeatCapacity(mesh%nnodes))  ! nondimensional heat capacity
        ALLOCATE (ThermalConductivity(mesh%nnodes))  ! nondimensional thermal conductivity
        ALLOCATE (qThermalConductivity(mesh%nq))  ! nondimensional thermal conductivity
        ALLOCATE (gradTC(mesh%nnodes,3))
        ALLOCATE (dRhoCpdt(mesh%nnodes)) ! d(rho*cp)/dt
        ALLOCATE (prhocp(mesh%nnodes))   ! rho*cp in previous time step

!       Set up  initial conditions
        CALL SetUpRhoCpLambdaAndBC(env,io,inp,mesh,gauss,Density,HeatCapacity, &
                                   ThermalConductivity,Temp,qTemp,velocity,inp%eeqm,Trhsv,cnt%rtime,1,PartVolFrac)

!       Count and set up crs matrices
        CALL sMat2crsSysRhsB_count(4,mesh,Tsysm,Trhsm)
!       Set up preconditioning vector (for diagonal preconditioning)
        ALLOCATE(Tprecv(mesh%Tnunk))

      ELSE IF (inp%iDT.NE.0) THEN
        CALL WarnErr(env,io,inp,4,"Tritok","Wrong iDT code!",inp%iDT)
      END IF
      IF (inp%iDT.GT.0) THEN
        IF (env%myproc.EQ.1) THEN
          WRITE (io%l,'(A)') "makro kinetika za notranje temperature"
          WRITE (io%l,'(A,I7)') "stevilo enacb   =",Tsysm%neq
          WRITE (io%l,'(A,I7/)') "stevilo neznank =",mesh%Tnunk
        END IF
      END IF
!
!     tvorim CRS sysm in rhsm za kinetiko (general diff-conv equation)
!
      IF (inp%iDC.EQ.1.OR.inp%iDC.EQ.2.OR.inp%iDC.EQ.4) THEN ! predolocen sistem enacb
        CALL sMat2crsSysRhsB_count(5,mesh,DCsysm,DCrhsm)
!       Set up preconditioning vector (for diagonal preconditioning)
        ALLOCATE(DCprecv(mesh%DCnunk))
      ELSE IF (inp%iDC.EQ.3) THEN ! dolocen sistem enacb
        CALL sMat2crsSysRhsB_SQsetup(5,mesh)
        CALL sMat2crsSysRhsB_SQcount(5,mesh,DCsysm,DCrhsm)
        CALL CopyMatrix(DCsysm,DCprec)
      ELSE IF (inp%iDC.NE.0) THEN
        CALL WarnErr(env,io,inp,4,"Tritok","Wrong iDC code!",inp%iDC)
      END IF

      IF (inp%iDC.GT.0) THEN
        IF (env%myproc.EQ.1) THEN
          WRITE (io%l,'(A)') "makro kinetika za splosno diff-conv enacbo"
          WRITE (io%l,'(A,I7)') "stevilo enacb   =",DCsysm%neq
          WRITE (io%l,'(A,I7/)') "stevilo neznank =",mesh%DCnunk
        END IF
      END IF

!
!     tvorim CRS sysm in rhsm za kinetiko (modified Helmholtz equation)
!
      IF (inp%imH.EQ.1) THEN ! predolocen sistem enacb
        CALL sMat2crsSysRhsB_count(6,mesh,mHsysm,mHrhsm)
!       Set up preconditioning vector (for diagonal preconditioning)
        ALLOCATE(mHprecv(mesh%mHnunk))
!       Fill system and rhs matrices
        CALL sMat2crsSysRhsB_modHelm_fill(mesh,mHsmatH,mHsmatG,mHsysm,mHrhsm)
!       preconditioner            
        CALL prelsqr2(mHsysm%neq,mesh%mHnunk,mHsysm%nnz,mHprecv,mHsysm%i,mHsysm%j,mHsysm%v)
      ELSE IF (inp%imH.EQ.2) THEN ! predolocen sistem enacb
        CALL sMat2crsSysRhsB_count(6,mesh,mHsysm,mHrhsm)
!       Set up preconditioning vector (for diagonal preconditioning)
        ALLOCATE(mHprecv(mesh%mHnunk))
      ELSE IF (inp%imH.NE.0.AND.inp%imH.NE.3.AND.inp%imH.NE.4.AND.inp%imH.NE.5.AND.inp%imH.NE.6.AND.inp%imH.NE.7) THEN
        CALL WarnErr(env,io,inp,4,"Tritok","Wrong imH code!",inp%imH)
      END IF

      IF (inp%imH.EQ.1.OR.inp%imH.EQ.2) THEN
        IF (env%myproc.EQ.1) THEN
          WRITE (io%l,'(A)') "makro kinetika za modificirano Helmholtzovo enacbo"
          WRITE (io%l,'(A,I7)') "stevilo enacb   =",mHsysm%neq
          WRITE (io%l,'(A,I7/)') "stevilo neznank =",mesh%mHnunk
        END IF
      END IF


!
!     tvorim CRS sysm in rhsm za makro kinematiko za notranje hitrosti
!     ker so robni pogoji vedno predpisane hitrosti imam samo eno sistemsko matriko za vse tri enacbe
!     
      IF (inp%iDv.EQ.1) THEN
        CALL sMat2crsSysRhsKM3(mesh,smatH,KMsysm,KMrhsm,kmp)
!       Set up preconditioning vector (for diagonal preconditioning)
        ALLOCATE(KMprecv(kmp%nsmcol))
        CALL prelsqr2(KMsysm%neq,kmp%nsmcol,KMsysm%nnz,KMprecv,KMsysm%i,KMsysm%j,KMsysm%v)      
      ELSE IF (inp%iDv.EQ.2) THEN
        CALL sMat2crsSysRhsKM3sq(mesh,smatH,KMsysm,KMrhsm,kmp)
!       Set up preconditioning matrix : 2 - ilu, 1 - diag, 0 - none
        CALL CopyMatrix(KMsysm,KMprec)
        CALL FormPRMjr(inp%sqrs_prec,KMprec,KMsysm)
      ELSE IF (inp%iDv.NE.0) THEN
        CALL WarnErr(env,io,inp,4,"Tritok","Wrong iDv code!",inp%iDv)
      END IF
      IF (inp%iDv.GT.0) THEN
        IF (env%myproc.EQ.1) THEN
          WRITE (io%l,'(A)') "makro kinematika za notranje hitrosti"
          WRITE (io%l,'(A,I7)')  "stevilo enacb   =",KMsysm%neq
          WRITE (io%l,'(A,I7/)') "stevilo neznank =",kmp%nsmcol
        END IF
      END IF

!     Dont need this matrices anymore
      IF (inp%iDC.EQ.0.AND.inp%iDT.NE.3.AND.inp%iDw.NE.3) THEN
        DEALLOCATE (smatG)
        DEALLOCATE (smatH)
      END IF

      CALL VmesniCas(cpu,itim)
      
!      ----------------------------------------------------------------------------
!                                 PRE TIME LOOP INIT
!      ----------------------------------------------------------------------------      
      cnt%glnnlit=0
      cnt%nit=0
      cnt%tstep=0
      cnt%trtstep=0
      cnt%rtime=0.0D0
      pptsvorticity=vorticity ! u pred previous = u previous
      ptsvorticity=vorticity ! u previous time step = u      
      pnlsvorticity=vorticity
      pnlsvelocity=velocity

      cnt%ftime=1
      IF (inp%iDT.GT.0) THEN
        pptstemp=temp ! u pred previous = u previous
        ptstemp=temp ! u previous time step = u      
        pnlstemp=temp
        IF (inp%iDT.EQ.3) THEN
          DO i=1,mesh%nnodes
            prhocp(i)=Density(i)*HeatCapacity(i)
          END DO
        END IF
      END IF      

      IF (inp%iDC.GT.0) THEN
        pptsdc=dc ! u pred previous = u previous
        ptsdc=dc ! u previous time step = u      
        pnlsdc=dc
!       set up initial sources
        DCsource=0.0D0
!       Set time derivative initial condition
        IF (inp%iTest.NE.0) THEN
          CALL DCinitial(env,io,mesh,inp,dc,ptsdc,pptsdc)
        END IF

      END IF   

      IF (inp%imH.GT.0) THEN
        IF (inp%imH.EQ.2.OR.inp%imH.EQ.3.OR.inp%imH.EQ.4.OR.inp%imH.EQ.5.OR.inp%imH.EQ.6.OR.inp%imH.EQ.7) THEN
          CALL modHelm_SetProblem(inp,mesh,mHdiff,mHrhsV,mHtime,mHu,mHq,fmqmH,cnt%rtime,velocity,mHanal)
        END IF
        pptsmHu=mHu ! u pred previous = u previous
        ptsmHu=mHu ! u previous time step = u      
        pnlsmHu=mHu
      END IF

!
!     Read restart file
!          
      IF (inp%rst.EQ.1) CALL RrstF(env,mesh,io,inp,cnt,velocity,vorticity, &
                   ptsvorticity,pptsvorticity,qvorticity,pnlsvorticity,pnlsvelocity, &
                   temp,ptstemp,pptstemp,qtemp,pnlstemp, &
                   dc,ptsdc,pptsdc,qdc,pnlsdc)


!      ----------------------------------------------------------------------------
!                                 START TIME LOOP
!      ----------------------------------------------------------------------------
1618  CONTINUE
        cnt%trtstep=cnt%trtstep+1
        cnt%rtime=cnt%rtime+inp%tstep
        cnt%tstep=cnt%tstep+1
        cnt%tsnnlit=0  ! number of nonlinear iteration within time step=0
        
        inp%dlse=inp%lsqs_eps ! start with full solver accuracy
                  

        pptsvorticity=ptsvorticity ! u pred previous = u previous
        ptsvorticity=vorticity     ! u previous time step = u       
        DO j=1,3
          DO i=1,mesh%nnodes
!           od parcialnega odvoda po casu
            rhsv(i,j)=-ptsvorticity(i,j)*inp%beta2+pptsvorticity(i,j)*inp%beta3  
          END DO
        END DO

!       \p T / \p t  - temperature equation
        IF (inp%iDT.GT.0) THEN
          pptstemp=ptstemp ! u pred previous = u previous
          ptstemp=temp     ! u previous time step = u       
          DO i=1,mesh%nnodes
            Ttime(i)=-ptsTemp(i)*inp%beta2+pptsTemp(i)*inp%beta3
          END DO

          IF (inp%iDT.EQ.3) THEN
!           Calculate time derivative of rho*cp
            CALL TimeDerRhoCp(mesh,dRhoCpdt,prhocp,Density,HeatCapacity,inp%tstep)
          END IF

        END IF

!       general D-C equation
        IF (inp%iDC.GT.0) THEN         
          pptsdc=ptsdc ! u pred previous = u previous
          ptsdc=dc     ! u previous time step = u       
          DO i=1,mesh%nnodes
            DCtime(i)=-ptsdc(i)*inp%beta2+pptsdc(i)*inp%beta3
          END DO
        END IF

!       modified Helmholtz equation
        IF (inp%imH.GT.0) THEN         
          pptsmHu=ptsmHu ! u pred previous = u previous
          ptsmHu=mHu     ! u previous time step = u       
          DO i=1,mesh%nnodes
            mHtime(i)=-ptsmHu(i)*inp%beta2+pptsmHu(i)*inp%beta3
          END DO
        END IF


!      ----------------------------------------------------------------------------
!                                 START NONLINEAR LOOP
!      ----------------------------------------------------------------------------
2718    CONTINUE
          cnt%tsnnlit=cnt%tsnnlit+1  ! number of nonlinear iterations within time step
          cnt%glnnlit=cnt%glnnlit+1  ! number of all nonlinear iterations
!      ----------------------------------------------------------------------------
!                                 BOUNDARY VORTICITIES
!      ----------------------------------------------------------------------------

          IF (inp%iBw.EQ.1) THEN
          eqn=1
          CALL sdkm_solve(eqn,env,inp,cpu,mesh,cnt%nit(eqn),sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz, &
                          sdAx,sdPivotx,velocity,vorticity)

          eqn=2
          CALL sdkm_solve(eqn,env,inp,cpu,mesh,cnt%nit(eqn),sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz, &
                          sdAy,sdPivoty,velocity,vorticity)

          eqn=3
          CALL sdkm_solve(eqn,env,inp,cpu,mesh,cnt%nit(eqn),sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz, &
                          sdAz,sdPivotz,velocity,vorticity)


          ELSE IF (inp%iBw.EQ.2) THEN
          eqn=1
          CALL sdkm_solve8(eqn,env,inp,cpu,mesh,cnt%nit(eqn),sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,sdDxB,sdDyB,sdDzB, &
                          sdAx,sdPivotx,velocity,vorticity)

          eqn=2
          CALL sdkm_solve8(eqn,env,inp,cpu,mesh,cnt%nit(eqn),sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,sdDxB,sdDyB,sdDzB, &
                          sdAy,sdPivoty,velocity,vorticity)

          eqn=3
          CALL sdkm_solve8(eqn,env,inp,cpu,mesh,cnt%nit(eqn),sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz,sdDxB,sdDyB,sdDzB, &
                          sdAz,sdPivotz,velocity,vorticity)

          ELSE IF (inp%iBw.EQ.3) THEN

          eqn=1
          CALL sdkm_solveBCN(eqn,env,inp,cpu,mesh,cnt%nit(eqn),sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz, &
                          sdAx,sdPivotx,velocity,vorticity)

          eqn=2
          CALL sdkm_solveBCN(eqn,env,inp,cpu,mesh,cnt%nit(eqn),sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz, &
                          sdAy,sdPivoty,velocity,vorticity)

          eqn=3
          CALL sdkm_solveBCN(eqn,env,inp,cpu,mesh,cnt%nit(eqn),sdH,sdHtx,sdHty,sdHtz,sdDx,sdDy,sdDz, &
                          sdAz,sdPivotz,velocity,vorticity)

          END IF

!      ----------------------------------------------------------------------------
!                                 DOMAIN VELOCITIES
!      ----------------------------------------------------------------------------
          IF (inp%iDv.EQ.1) THEN
            eqn=1
            CALL solv_kmMacro3(eqn,inp,cpu,mesh,cnt%nit(3+eqn),KMsysm,KMrhsm,KMprecv,kmp, &
                             smatHtx,smatHty,smatHtz, &
                             smatDx,smatDy,smatDz,velocity,vorticity)

            eqn=2
            CALL solv_kmMacro3(eqn,inp,cpu,mesh,cnt%nit(3+eqn),KMsysm,KMrhsm,KMprecv,kmp, &
                             smatHtx,smatHty,smatHtz, &
                             smatDx,smatDy,smatDz,velocity,vorticity)

            eqn=3
            CALL solv_kmMacro3(eqn,inp,cpu,mesh,cnt%nit(3+eqn),KMsysm,KMrhsm,KMprecv,kmp, &
                             smatHtx,smatHty,smatHtz, &
                             smatDx,smatDy,smatDz,velocity,vorticity)
          ELSE IF (inp%iDv.EQ.2) THEN
            eqn=1
            CALL solv_kmMacro3sq(eqn,inp,cpu,mesh,cnt%nit(3+eqn),KMsysm,KMrhsm,KMprec,kmp, &
                             smatHtx,smatHty,smatHtz, &
                             smatDx,smatDy,smatDz,velocity,vorticity)             

            eqn=2
            CALL solv_kmMacro3sq(eqn,inp,cpu,mesh,cnt%nit(3+eqn),KMsysm,KMrhsm,KMprec,kmp, &
                             smatHtx,smatHty,smatHtz, &
                             smatDx,smatDy,smatDz,velocity,vorticity)             

            eqn=3
            CALL solv_kmMacro3sq(eqn,inp,cpu,mesh,cnt%nit(3+eqn),KMsysm,KMrhsm,KMprec,kmp, &
                             smatHtx,smatHty,smatHtz, &
                             smatDx,smatDy,smatDz,velocity,vorticity)    
          END IF
!      ----------------------------------------------------------------------------
!                                   DOMAIN TEMPERATURES
!      ----------------------------------------------------------------------------
          IF (inp%iDT.EQ.1) THEN
            CALL SolveFDdae(env,io,inp,cpu,mesh,Tsysm,Trhsm,Tprecv,temp,qtemp,velocity,Ttime, &
                            smatB,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz, &
                            cnt%nit(10),inp%Ren*inp%prn/nanoC)
          ELSE IF (inp%iDT.EQ.2) THEN
            CALL SolveFDdaeSQ(env,io,inp,cpu,mesh,Tsysm,Trhsm,Tprec,temp,qtemp,velocity,Ttime, &
                            smatB,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz, &
                            cnt%nit(10),inp%Ren*inp%prn/nanoC)
          ELSE IF (inp%iDT.EQ.3) THEN
!           set up density, heat capacity, thermalconductivity and  Energy equation boundary conditions
            CALL SetUpRhoCpLambdaAndBC(env,io,inp,mesh,gauss,Density,HeatCapacity, &
                                       ThermalConductivity,Temp,qTemp,velocity,inp%eeqm,Trhsv,cnt%rtime,0,PartVolFrac)
!           interpolate diffusivity to macro mesh
            CALL Int2MacroNodes(mesh,ThermalConductivity,qThermalConductivity)
!           Calculate gradient of thermal conductivity
            CALL setGrad(mesh,ThermalConductivity,gradTC)
!           system matrix
            CALL sMat2crsSysRhsB_TempRCPL_fill(mesh,smatH,smatG,smatB,inp%beta,inp%ren*inp%prn/nanoC,Tsysm,Trhsm, &
                            smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz, &
                            Density,HeatCapacity,ThermalConductivity,qThermalConductivity,dRhoCpdt,velocity,gradTC)
!           preconditioner
            CALL prelsqr2(Tsysm%neq,mesh%Tnunk,Tsysm%nnz,Tprecv,Tsysm%i,Tsysm%j,Tsysm%v)
!           solve
            CALL SolveFDdaeRCPL(env,io,inp,cpu,mesh,Tsysm,Trhsm,Tprecv,temp,qtemp,velocity,Ttime,Trhsv, &
                                smatB,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz, &
                                cnt%nit(10),inp%Ren*inp%prn/nanoC, &
                                Density,HeatCapacity,ThermalConductivity,gradTC)

          END IF
!      ----------------------------------------------------------------------------
!                                   MODIFIED HELMHOLZ EQUATION
!      ----------------------------------------------------------------------------

          IF (inp%imH.GT.0) THEN

            IF (inp%imH.EQ.1) THEN ! resi nabla^2 u - mu^2 u = f
!             set up rhs
              CALL modHelm_SetUpRHS(inp,mesh,mHrhsV,mHtime,mHu,mHdiff)
              CALL SolveModHelm(env,io,inp,cpu,mesh,mHdiff,mHqDiff,mHsysm,mHrhsm,mHprecv,mHu,mHq,mHrhsV, &
                            mHsmatB,cnt%nit(12))
              CALL modHelm_calRMS(inp,gauss,mesh,mHu)

            ELSE IF (inp%imH.EQ.2) THEN ! resi D-C enacbo
!             set up diffusicity, rhs, velocity, etc
              CALL modHelm_SetProblem(inp,mesh,mHdiff,mHrhsV,mHtime,mHu,mHq,fmqmH,cnt%rtime,velocity,mHanal)
!             calculate difusivitiy gradients
              CALL setGrad(mesh,mHdiff,mHdiffGrad)
!             interpolate diffusivity to macro mesh
              CALL Int2MacroNodes(mesh,mHdiff,mHqDiff)
!             set up system matrix
              CALL sMat2crsSysRhsB_modHelm_fill2(mesh,mHsmatH,mHsmatG,mHdiff,mHqDiff,mHsysm,mHrhsm)
!             preconditioner
              CALL prelsqr2(mHsysm%neq,mesh%mHnunk,mHsysm%nnz,mHprecv,mHsysm%i,mHsysm%j,mHsysm%v)
!             solve
              CALL SolveModHelmDC(env,io,inp,cpu,mesh,mHsysm,mHrhsm,mHprecv,mHu,mHq,velocity,mHrhsV, &
                            mHsmatB,mHsmatAbdx,mHsmatAbdy,mHsmatAbdz,mHsmatDx,mHsmatDy,mHsmatDz, &
                            cnt%nit(12),mHdiff,mHqDiff,mHdiffGrad)

            ELSE IF (inp%imH.EQ.3) THEN ! resi D-C enacbo z full matrix single domain
!             set up diffusicity, rhs, velocity, etc
              CALL modHelm_SetProblem(inp,mesh,mHdiff,mHrhsV,mHtime,mHu,mHq,fmqmH,cnt%rtime,velocity,mHanal)
!             solve
              CALL fmSolveModHelmDC(mesh,cpu,inp,env,io,fmmHA,fmmHRHS,fmmHPivot,mHu,fmqmH,velocity,mHrhsV, &
                     UmatB,UmatAbdx,UmatAbdy,UmatAbdz,UmatDx,UmatDy,UmatDz, &
                     QmatB,QmatAbdx,QmatAbdy,QmatAbdz,QmatDx,QmatDy,QmatDz, &
                     cnt%nit(12),mHdiff,mHqDiff,mHdiffGrad)

            ELSE IF (inp%imH.EQ.4) THEN ! resi D-C enacbo z full matrix single domain
!                                         MATMUL verzija
!             set up diffusicity, rhs, velocity, etc
              CALL modHelm_SetProblem(inp,mesh,mHdiff,mHrhsV,mHtime,mHu,mHq,fmqmH,cnt%rtime,velocity,mHanal)
!             solve
              CALL fmSolveModHelmDCmm(mesh,cpu,inp,env,io,fmmHA,fmmHRHS,fmmHPivot,mHu,fmqmH,velocity,mHrhsV, &
                     UmatB,UmatBp,UmatAbdx,UmatAbdy,UmatAbdz,UmatDx,UmatDy,UmatDz, &
                     QmatB,QmatBp,QmatAbdx,QmatAbdy,QmatAbdz,QmatDx,QmatDy,QmatDz, &
                     cnt%nit(12),mHdiff,mHqDiff,mHdiffGrad)

            ELSE IF (inp%imH.EQ.5) THEN ! resi D-C enacbo z full matrix single domain
!                                         Wavelet verzija
!             set up diffusicity, rhs, velocity, etc
              CALL modHelm_SetProblem(inp,mesh,mHdiff,mHrhsV,mHtime,mHu,mHq,fmqmH,cnt%rtime,velocity,mHanal)
!             solve
              CALL fmSolveModHelmDCwawt(mesh,cpu,inp,env,io,fmmHA,fmmHRHS,fmmHPivot,mHu,fmqmH,velocity,mHrhsV, &
                     cUmatB,cUmatBp,cUmatAbdx,cUmatAbdy,cUmatAbdz,cUmatDx,cUmatDy,cUmatDz, &
                     cQmatB,cQmatBp,cQmatAbdx,cQmatAbdy,cQmatAbdz,cQmatDx,cQmatDy,cQmatDz, &
                     cnt%nit(12),mHdiff,mHqDiff,mHdiffGrad)

            ELSE IF (inp%imH.EQ.6) THEN ! resi D-C enacbo z full matrix single domain
!                                         ACA verzija
!             set up diffusicity, rhs, velocity, etc
              CALL modHelm_SetProblem(inp,mesh,mHdiff,mHrhsV,mHtime,mHu,mHq,fmqmH,cnt%rtime,velocity,mHanal)
!             solve
              CALL fmSolveModHelmDCaca(mesh,cpu,inp,env,io,fmmHA,fmmHRHS,fmmHPivot,mHu,fmqmH,velocity,mHrhsV, &
                     aUmatB,aUmatBp,aUmatAbdx,aUmatAbdy,aUmatAbdz,aUmatDx,aUmatDy,aUmatDz, &
                     aQmatB,aQmatBp,aQmatAbdx,aQmatAbdy,aQmatAbdz,aQmatDx,aQmatDy,aQmatDz, &
                     cnt%nit(12),mHdiff,mHqDiff,mHdiffGrad)


            ELSE IF (inp%imH.EQ.7) THEN ! resi D-C enacbo z full matrix single domain
!                                         H matrix version
!             set up diffusicity, rhs, velocity, etc
              CALL modHelm_SetProblem(inp,mesh,mHdiff,mHrhsV,mHtime,mHu,mHq,fmqmH,cnt%rtime,velocity,mHanal)
!             solve
              CALL fmSolveModHelmDChmat(mesh,cpu,inp,env,io,fmmHA,fmmHRHS,fmmHPivot,mHu,fmqmH,velocity,mHrhsV, &
                     hUmatB,hUmatBp,hUmatAbdx,hUmatAbdy,hUmatAbdz,hUmatDx,hUmatDy,hUmatDz, &
                     hQmatB,hQmatBp,hQmatAbdx,hQmatAbdy,hQmatAbdz,hQmatDx,hQmatDy,hQmatDz, &
                     cnt%nit(12),mHdiff,mHqDiff,mHdiffGrad)

            END IF

            IF (inp%imH.GE.3) THEN
                CALL FluxFM2macro(mesh,fmqmH,mHq)
            END IF


          END IF

!      ----------------------------------------------------------------------------
!                                   KELVIN BODY FORCE
!      ----------------------------------------------------------------------------
          IF (inp%iFerro.GT.0) THEN
!           calculate temperature gradient
            CALL setGrad(mesh,temp,tempGrad)
            CALL CalKelvinBodyForce(mesh,inp,temp,mfs,hshg,kbf,tempGrad)
          END IF

!      ----------------------------------------------------------------------------
!                                   DOMAIN VORTICITIES
!      ----------------------------------------------------------------------------

          IF (inp%iDw.EQ.1) THEN
            eqn=1
            CALL SolveFDwTE(eqn,env,io,inp,cpu,mesh,VTEsysmX,VTErhsmX,VTEprecvX, &
                          vorticity,qvorticity,temp,dc,velocity,kbf,rhsv(:,eqn), &
                          smatB,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,cnt%nit(6+eqn),nanoA,nanoB)

            eqn=2
            CALL SolveFDwTE(eqn,env,io,inp,cpu,mesh,VTEsysmY,VTErhsmY,VTEprecvY, &
                          vorticity,qvorticity,temp,dc,velocity,kbf,rhsv(:,eqn), &
                          smatB,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,cnt%nit(6+eqn),nanoA,nanoB)

            eqn=3
            CALL SolveFDwTE(eqn,env,io,inp,cpu,mesh,VTEsysmZ,VTErhsmZ,VTEprecvZ, &
                        vorticity,qvorticity,temp,dc,velocity,kbf,rhsv(:,eqn), &
                        smatB,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,cnt%nit(6+eqn),nanoA,nanoB)


          ELSE IF (inp%iDw.EQ.2) THEN
            eqn=1
            CALL SolveFDwTEsq(eqn,env,io,inp,cpu,mesh,VTEsysmX,VTErhsmX,VTEprecX, &
                          vorticity,qvorticity,temp,dc,velocity,kbf,rhsv(:,eqn), &
                          smatB,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,cnt%nit(6+eqn),nanoA,nanoB)
            eqn=2
            CALL SolveFDwTEsq(eqn,env,io,inp,cpu,mesh,VTEsysmY,VTErhsmY,VTEprecY, &
                          vorticity,qvorticity,temp,dc,velocity,kbf,rhsv(:,eqn), &
                          smatB,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,cnt%nit(6+eqn),nanoA,nanoB)
            eqn=3
            CALL SolveFDwTEsq(eqn,env,io,inp,cpu,mesh,VTEsysmZ,VTErhsmZ,VTEprecZ, &
                          vorticity,qvorticity,temp,dc,velocity,kbf,rhsv(:,eqn), &
                          smatB,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,cnt%nit(6+eqn),nanoA,nanoB)

          ELSE IF (inp%iDw.EQ.3) THEN
!           use model to set up viscosity and cvte
            CALL SetUpViscCvte(env,io,inp,mesh,viscosity,cvte,temp)
!           interpolate Viscosity to macro mesh
            CALL Int2MacroNodes(mesh,Viscosity,qViscosity)
!           Calculate gradient of Viscosity
            CALL setGrad(mesh,Viscosity,GradVisc)
!           set up system matrices

            CALL sMat2crsSysRhsB_w_fill_vMU(1,mesh,smatH,smatG,smatB,inp%beta,inp%ren/nanoA,VTEsysmX,VTErhsmX, &
                                      velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,Viscosity,qViscosity)
            CALL sMat2crsSysRhsB_w_fill_vMU(2,mesh,smatH,smatG,smatB,inp%beta,inp%ren/nanoA,VTEsysmY,VTErhsmY, &
                                      velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,Viscosity,qViscosity)
            CALL sMat2crsSysRhsB_w_fill_vMU(3,mesh,smatH,smatG,smatB,inp%beta,inp%ren/nanoA,VTEsysmZ,VTErhsmZ, &
                                      velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,Viscosity,qViscosity)

!           calculate preconditioner
            CALL prelsqr2(VTEsysmX%neq,mesh%nunk(1),VTEsysmX%nnz,VTEprecvX,VTEsysmX%i,VTEsysmX%j,VTEsysmX%v)
            CALL prelsqr2(VTEsysmY%neq,mesh%nunk(2),VTEsysmY%nnz,VTEprecvY,VTEsysmY%i,VTEsysmY%j,VTEsysmY%v)
            CALL prelsqr2(VTEsysmZ%neq,mesh%nunk(3),VTEsysmZ%nnz,VTEprecvZ,VTEsysmZ%i,VTEsysmZ%j,VTEsysmZ%v)
!           calculate terms that involve viscosity gradient
            CALL ViscGradTerms(mesh,velocity,vorticity,GradVisc,DvaEpsGVisc,CuOmCrGrVi)
!
!           solve
!
            eqn=1
            CALL SolveFDwTE_vMU(eqn,env,io,inp,cpu,mesh,VTEsysmX,VTErhsmX,VTEprecvX, &
                          vorticity,qvorticity,temp,dc,velocity,kbf,rhsv(:,eqn), &
                          smatB,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,cnt%nit(6+eqn),nanoA,nanoB, &
                          DvaEpsGVisc,CuOmCrGrVi,cvte)

            eqn=2
            CALL SolveFDwTE_vMU(eqn,env,io,inp,cpu,mesh,VTEsysmY,VTErhsmY,VTEprecvY, &
                          vorticity,qvorticity,temp,dc,velocity,kbf,rhsv(:,eqn), &
                          smatB,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,cnt%nit(6+eqn),nanoA,nanoB, &
                          DvaEpsGVisc,CuOmCrGrVi,cvte)

            eqn=3
            CALL SolveFDwTE_vMU(eqn,env,io,inp,cpu,mesh,VTEsysmZ,VTErhsmZ,VTEprecvZ, &
                          vorticity,qvorticity,temp,dc,velocity,kbf,rhsv(:,eqn), &
                          smatB,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,cnt%nit(6+eqn),nanoA,nanoB, &
                          DvaEpsGVisc,CuOmCrGrVi,cvte)


          END IF 

!      ----------------------------------------------------------------------------
!                                   GENERAL DIFF-CONV EQUATION
!      ----------------------------------------------------------------------------
          IF (inp%iDC.GT.0) THEN
!           set up sources
            IF (inp%itest.GT.0) CALL SetUpDCtest(env,io,mesh,inp,DCsource,DCanal,cnt%rtime,dc,qdc)
!           set up diffusivity and rhs vector
            CALL SetUpDiffSources(mesh,inp,DCdiff,DCrhsv,DCtime,DCsource)

            IF (inp%iDC.EQ.1) THEN  ! verzija, ki rabi laplace
!             calculate difusivitiy gradients
              CALL CalDiffGrad(mesh,DCdiff,DCdiffGrad,DCdiffLap)
!             set up system matrix
              CALL sMat2crsSysRhsB_DC_fill(mesh,smatH,smatG,smatB,inp%beta,DCdiff,DCsysm,DCrhsm, &
                                     velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)              
!             preconditioner
              CALL prelsqr2(Dcsysm%neq,mesh%DCnunk,DCsysm%nnz,Dcprecv,DCsysm%i,Dcsysm%j,DCsysm%v)
!             solve
              CALL SolveFDdaeVD(env,io,inp,cpu,mesh,DCsysm,DCrhsm,DCprecv,dc,qdc,velocity,DCrhsv, &
                            smatB,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz, &
                            cnt%nit(11),DCdiff,DCdiffGrad,DCdiffLap)
            ELSE IF (inp%iDC.EQ.2) THEN ! verzija ki ne rabi Laplace
!             calculate difusivitiy gradients
              CALL setGrad(mesh,DCdiff,DCdiffGrad)
!             interpolate diffusivity to macro mesh
              CALL Int2MacroNodes(mesh,DCdiff,DCqDiff)
!             set up system matrix
              CALL sMat2crsSysRhsB_DC_fill2(mesh,smatH,smatG,smatB,inp%beta,DCdiff,DCqDiff,DCsysm,DCrhsm, &
                                     velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)
!             preconditioner
              CALL prelsqr2(Dcsysm%neq,mesh%DCnunk,DCsysm%nnz,Dcprecv,DCsysm%i,Dcsysm%j,DCsysm%v)
!             solve
              CALL SolveFDdaeVD2(env,io,inp,cpu,mesh,DCsysm,DCrhsm,DCprecv,dc,qdc,velocity,DCrhsv, &
                            smatB,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz, &
                            cnt%nit(11),DCdiff,DCdiffGrad)

            ELSE IF (inp%iDC.EQ.3) THEN ! verzija ki ne rabi Laplace, DOLOCEN SISTEM ENACB

!             calculate difusivitiy gradients
              CALL setGrad(mesh,DCdiff,DCdiffGrad)
!             interpolate diffusivity to macro mesh
              CALL Int2MacroNodes(mesh,DCdiff,DCqDiff)
!             set up system matrix
              CALL sMat2crsSysRhsB_DC_fill2sq(mesh,smatH,smatG,smatB,inp%beta,DCdiff,DCqDiff,DCsysm,DCrhsm, &
                                     velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)
!             preconditioner
              CALL dcopy (DCprec%nnz, DCsysm%v, 1, DCprec%v,1)
              CALL FormPRMjr(inp%sqrs_prec,DCprec,DCsysm)
!             solve
              CALL SolveFDdaeVD2sq(env,io,inp,cpu,mesh,DCsysm,DCrhsm,DCprec,dc,qdc,velocity,DCrhsv, &
                            smatB,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz, &
                            cnt%nit(11),DCdiff,DCdiffGrad)

            ELSE IF (inp%iDC.EQ.4) THEN ! verzija ki ne rabi Laplace, IN JE ZA LINEARNE PROBLEME
!             calculate difusivitiy gradients
              CALL setGrad(mesh,DCdiff,DCdiffGrad)
!             interpolate diffusivity to macro mesh
              CALL Int2MacroNodes(mesh,DCdiff,DCqDiff)
!             set up system matrix
              CALL sMat2crsSysRhsB_DC_fill4(mesh,smatH,smatG,smatB,inp%beta,DCdiff,DCdiffGrad,DCqDiff,DCsysm,DCrhsm, &
                                    velocity,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz)
!             preconditioner
              CALL prelsqr2(Dcsysm%neq,mesh%DCnunk,DCsysm%nnz,Dcprecv,DCsysm%i,Dcsysm%j,DCsysm%v)
!             solve
              CALL SolveFDdaeVD4(env,io,inp,cpu,mesh,DCsysm,DCrhsm,DCprecv,dc,qdc,velocity,DCrhsv, &
                            smatB,smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz, &
                            cnt%nit(11),DCdiff,DCdiffGrad)
            END IF
          END IF

!      ----------------------------------------------------------------------------
!                                   CHECK CONVERGENCE
!      ----------------------------------------------------------------------------
          CALL CheckConvergencevw(env,io,inp,mesh,cnt,velocity,pnlsvelocity,vorticity,pnlsvorticity, &
                                 temp,pnlstemp,dc,pnlsdc,mHu,pnlsmHu)
     
          CALL CheckControlFiles(env,io,inp,istop) 
          
          IF (inp%iDLSE.EQ.1.AND.cnt%tsnnlit.GT.1) THEN  ! prva dva naredi s 1E-7, nato dvigne
            IF (cnt%tsnnlit.EQ.2) inp%dlse=inp%dlse_epsth
            CALL DLSE(env,io,cnt,inp)
          END IF
!      ----------------------------------------------------------------------------
!                                 END NONLINEAR LOOP
!      ----------------------------------------------------------------------------
        IF (cnt%tsnnlit.LT.inp%nnlit &
                                      .AND.(cnt%nlierr(1).GT.inp%eps &
                                        .OR.cnt%nlierr(2).GT.inp%eps &
                                        .OR.cnt%nlierr(3).GT.inp%eps &
                                        .OR.cnt%nlierr(4).GT.inp%eps &
                                        .OR.cnt%nlierr(5).GT.inp%eps &
                                        .OR.cnt%nlierr(6).GT.inp%eps &
                                        .OR.cnt%nlierr(7).GT.inp%eps &     
                                        .OR.cnt%nlierr(8).GT.inp%eps &   
                                        .OR.cnt%nlierr(9).GT.inp%eps) &     
                                      .AND.istop.EQ.0) GOTO 2718

!      ----------------------------------------------------------------------------
!                                 MOVE PARTICLES
!      ----------------------------------------------------------------------------
      IF (inp%iPart.GT.0) THEN
        CALL MoveParticles(inp,part,mesh,cpu,velocity,vorticity,Temp)
        CALL CalParticleVolumeFraction(env,io,inp,gauss,part,mesh,PartVolFrac)
      END IF


        IF (env%myproc.EQ.1) THEN
          IF ((MOD(cnt%trtstep,inp%iwrit).EQ.0).OR.istop.NE.0) THEN
!           ce testiranje, daj anal -> temp, abs,rel error -> vorticity
            IF (inp%itest.GT.0) CALL CopyTestFields(mesh,dc,DCanal,temp,vorticity)
!           izpise tri.walls.dat
            IF (inp%iRawl.GT.0) CALL ExportWallToTecplot(io,inp,cnt,mesh,velocity,vorticity,temp,dc,qtemp,qdc)
!           tri.rms
            IF (inp%eeqm.GT.10.AND.inp%iDT.EQ.3) CALL CalTempRMS(env,io,inp,mesh,cnt,Temp,cnt%rtime,inp%eeqm)
!           izpise delce
            IF (inp%iPart.GT.0) CALL OutputParticles(io,inp,cnt,part)
!           izpise tri.materials.dat
            IF (inp%iDt.EQ.3) CALL OutputTDMaterialProperties(io,inp,cnt,mesh,Density,HeatCapacity,ThermalConductivity,PartVolFrac)
!           izpise tri.dat
            CALL OutputTDResultsTecplotvw(io,inp,cnt,mesh,velocity,vorticity,temp,dc,mHu)
!           izpise rezultate po stenah (tri.twfl.dat)
            CALL ExportWallFlux(io%twfl,mesh,cnt,qtemp,qvorticity,mHq)
!           izpise tri.vtu
            IF (inp%iPara.GT.0) CALL OutputTDResultsParaview(io,inp,cnt,mesh,velocity,vorticity,temp,dc,mHu)
            IF (inp%iParaVTK.GT.0) CALL OutputTDResultsParaviewVTK(io,inp,cnt,mesh,velocity,vorticity,temp,dc,mHu)
          END IF
!           zapise RMS napako v log datoteko
            IF (inp%itest.GT.0) CALL CalRMS(io,inp%iTest,mesh,dc,DCanal,cnt%rtime)
        END IF

        CALL CheckTIMEConvergence(env,io,mesh,cnt,vorticity,ptsvorticity)
        
        IF (env%myproc.EQ.1) THEN
          IF (inp%iDT.EQ.1.OR.inp%iDT.EQ.2) THEN
            CALL WallFlux(io%tfl,mesh,cnt,gauss,qtemp," qT")
            CALL WallFuncInt(io%tfl,mesh,cnt,gauss,temp,"  T")
            CALL WallLocalFlux(io%lnu,mesh,cnt,gauss,qtemp)
          ELSE IF (inp%iDT.EQ.3) THEN
            CALL WallFlux(io%tfl,mesh,cnt,gauss,qtemp," qT")
            CALL WallFluxLambda(io%tfl,mesh,cnt,gauss,qtemp,"qTl",ThermalConductivity)
          END IF

          CALL MassFlux(io%mass,mesh,cnt,gauss,velocity)

          IF (inp%iDw.GT.0) THEN
            CALL WallFlux(io%wfl,mesh,cnt,gauss,qvorticity(:,1),"qwx")
            CALL WallFlux(io%wfl,mesh,cnt,gauss,qvorticity(:,2),"qwy")
            CALL WallFlux(io%wfl,mesh,cnt,gauss,qvorticity(:,3),"qwz")

            CALL WallFuncInt(io%wfl,mesh,cnt,gauss,vorticity(:,1)," wx")
            CALL WallFuncInt(io%wfl,mesh,cnt,gauss,vorticity(:,2)," wy")
            CALL WallFuncInt(io%wfl,mesh,cnt,gauss,vorticity(:,3)," wz")

            tmp=abs(vorticity(:,1))
            CALL WallFuncInt(io%wfl,mesh,cnt,gauss,tmp,"Awx")
            tmp=abs(vorticity(:,2))
            CALL WallFuncInt(io%wfl,mesh,cnt,gauss,tmp,"Awy")
            tmp=abs(vorticity(:,3))
            CALL WallFuncInt(io%wfl,mesh,cnt,gauss,tmp,"Awz")
          END IF

          IF (inp%imH.GT.0) THEN
            CALL WallFlux(io%wfl,mesh,cnt,gauss,mHq,"mHq") ! zapisem kar v worticity flux file
          END IF

          IF (inp%iDC.GT.0) THEN
            CALL WallFlux(io%wfl,mesh,cnt,gauss,qdc,"qdc") ! zapisem kar v worticity flux file
          END IF

!         Write restart file
          CALL WrstF(mesh,io,inp,cnt,velocity,vorticity, &
                  ptsvorticity,pptsvorticity,qvorticity,pnlsvorticity,pnlsvelocity, &
                  temp,ptstemp,pptstemp,qtemp,pnlstemp, &
                  dc,ptsdc,pptsdc,qdc,pnlsdc)

!         modified Helmholz, test cases, chck RMS
          IF (inp%imH.EQ.2.OR.inp%imH.EQ.3.OR.inp%imH.EQ.4.OR.inp%imH.EQ.5.OR.inp%imH.EQ.6.OR.inp%imH.EQ.7) THEN ! resi D-C enacbo
            CALL modHelm_calRMS2(io,inp,gauss,mesh,mHu,mHanal,cnt%rtime)
          END IF

         END IF

!      ----------------------------------------------------------------------------
!                                 END TIME LOOP
!      ----------------------------------------------------------------------------
      IF (cnt%trtstep.LT.inp%nstep.AND.istop.EQ.0) GOTO 1618               
!
!     Output data for stochastic postprocessing
!
      IF (inp%estd.EQ.1) CALL OutputStochasticResults(io,inp,cnt,gauss,mesh,velocity,vorticity,temp,qtemp,qvorticity)


      CALL VmesniCas(cpu,itim)
      cpu%time(itim+1)=cptime(cpu%t00)        
!
!     Write cpu information
!
      IF (env%myproc.EQ.1) CALL WriteCpu(io,cpu)        
!
!     Finish
!      
      CALL StopProgram(env,io,inp,0)    
      END
