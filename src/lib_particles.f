C
C     Particle Tracking Library
C
C -----------------------------------------------------------------------------
      SUBROUTINE MoveParticles(inp,part,mesh,cpu,velocity,vorticity,temperature)
C
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType)  :: mesh
      TYPE(CPUtype)   :: cpu
      TYPE(InputType) inp
      TYPE(particleType) part(inp%npart)
      TYPE(fluidAtPartType) :: fap ! tekocina na mestu delca

      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3)
      REAL(8) temperature(mesh%nnodes)
      REAL(4) cptime,rcpu
      INTEGER i,ierr,j
      REAL(8) rnd1,rnd2,rnd3 ! rnd number -1 to 1

      REAL(8) tstep
      REAL(8), ALLOCATABLE :: gradT(:,:)

      rcpu=cptime(0.)

      ALLOCATE (gradT(mesh%nnodes,3))
C     for thermophoresis caluclate temperature gradient
      IF (inp%iThph.GT.0) THEN
        CALL setGrad(mesh,temperature,gradT)
      ELSE
        gradT=0.0D0
      END IF

c     particle time step is smaller than fluid time step
      tstep=inp%tstep / DBLE(inp%PartTSR)

c     loop over particle time
      DO j=1,inp%PartTSR
        DO i=1,inp%npart
          IF (part(i)%active) THEN
            CALL GetFFFap(mesh,velocity,vorticity,temperature,gradT,fap,part(i),ierr)
            IF (ierr.EQ.0) THEN

C             massles particles velocity = fluid velocity + thermophoresis
              part(i)%vx=fap%vx+inp%thphMul*fap%gradT(1)/(fap%T+inp%thphT0/inp%thphDeltaT)
              part(i)%vy=fap%vy+inp%thphMul*fap%gradT(2)/(fap%T+inp%thphT0/inp%thphDeltaT)
              part(i)%vz=fap%vz+inp%thphMul*fap%gradT(3)/(fap%T+inp%thphT0/inp%thphDeltaT)

c             generate random fluctuation
              CALL rndXiEtaZeta(rnd1,rnd2,rnd3)

c             move
              part(i)%x=part(i)%x+part(i)%vx*tstep+rnd1*inp%BrownT0mul*SQRT(tstep)
              part(i)%y=part(i)%y+part(i)%vy*tstep+rnd2*inp%BrownT0mul*SQRT(tstep)
              part(i)%z=part(i)%z+part(i)%vz*tstep+rnd3*inp%BrownT0mul*SQRT(tstep)

c             preverimo ce so se zaleteli v steno
c             in apliciramo ustrezne robne pogoje
              CALL CavityBC(part(i))
            ELSE
              part(i)%active=.FALSE.
            END IF
          END IF
        END DO
      END DO

      DEALLOCATE (gradT)

      cpu%time(25)=cpu%time(25)+cptime(rcpu)

      END





C -----------------------------------------------------------------------------
      SUBROUTINE SetUpThermophoresis(env,io,inp)
C
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(InputType) inp
      TYPE(IOtype)    :: io
      TYPE(penv) :: env

      REAL(8) g0,DeltaT
      g0=9.81D0

      inp%thphMul=-inp%thphA * (inp%thphPartDiam / (2.0D0*1.0D-9) )**(-inp%thphB)*inp%PrN

      inp%thphDeltaT=inp%RaN*inp%PrN*inp%nanoF_k**2 /
     &   ( g0* inp%nanoF_beta* inp%thphDomainSize**3 * inp%nanoF_rho**2 *  inp%nanoF_cp**2 )


      IF (env%myproc.EQ.1) THEN
        WRITE (io%l,'(/A)') "Thermophoresis motion of particles"
        WRITE (io%l,'(A,G18.9)') "A      ",inp%thphA
        WRITE (io%l,'(A,G18.9)') "B      ",inp%thphB
        WRITE (io%l,'(A,G18.9)') "Particle diameter [m]",inp%thphPartDiam
        WRITE (io%l,'(A,G18.9)') "T0 [K] ",inp%thphT0
        WRITE (io%l,'(A,G18.9)') "Domain size [m]",inp%thphDomainSize
        WRITE (io%l,'(A,G18.9)') "Thermophoresis delta T [K]   =",inp%thphDeltaT
        WRITE (io%l,'(A,G18.9/)') "Thermophoresis multiplyer =",inp%thphMul
      END IF

      END


C -----------------------------------------------------------------------------
      SUBROUTINE SetUpBrown(env,io,inp)
C
C     Calculate particle volume fraction
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(InputType) inp
      TYPE(IOtype)    :: io
      TYPE(penv) :: env


      REAL(8) g0,deltat,kB,pi,visc,diff,ScN

      kB=1.380648D-23
      g0=9.81D0
      pi=4.0D0*ATAN(1.0D0)

      DeltaT=inp%RaN*inp%PrN*inp%nanoF_k**2 /
     &   ( g0* inp%nanoF_beta* inp%BrownDomainSize**3 * inp%nanoF_rho**2 *  inp%nanoF_cp**2 )

      visc = inp%PrN * inp%nanoF_k / ( inp%nanoF_rho *  inp%nanoF_cp )

      diff = kB * inp%BrownT0 / ( 3.0D0 * pi *  inp%nanoF_rho * visc * inp%BrownPartDiam)

      ScN = visc / diff

      inp%BrownT0mul=SQRT(6.0D0*inp%PrN/ScN)

      IF (env%myproc.EQ.1) THEN
        WRITE (io%l,'(/A)') "Brownian motion of particles"
        WRITE (io%l,'(A,G18.9)') "T0 [K] ",inp%BrownT0
        WRITE (io%l,'(A,G18.9)') "Domain size [m]",inp%BrownDomainSize
        WRITE (io%l,'(A,G18.9)') "Particle diameter [m]",inp%BrownPartDiam
        WRITE (io%l,'(A,G18.9)') "Delta T [K] =",deltaT
        WRITE (io%l,'(A,G18.9)') "Viscosity [m2/s] =",visc
        WRITE (io%l,'(A,G18.9)') "T0 Particle Diff. [m2/s] =",diff
        WRITE (io%l,'(A,G18.9)') "T0 Particle Schmidt number [-] =",ScN
        WRITE (io%l,'(A,G18.9/)') "T0 Brown multiplayer =",inp%BrownT0mul
      END IF

      END

C -----------------------------------------------------------------------------
      SUBROUTINE CalParticleVolumeFraction(env,io,inp,gp,part,mesh,PartVolFrac)
C
C     Calculate particle volume fraction
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType)  :: mesh
      TYPE(InputType) inp
      TYPE(particleType) part(inp%npart)
      TYPE(IOtype)    :: io
      TYPE(penv) :: env
      TYPE(gausstype) :: gp


      REAL(8) PartVolFrac(mesh%nnodes)
      REAL(8) fi0,ficrtica,sigma
      REAL(8) Ni(mesh%npoc) ! shape functions
      REAL(8) sumni,dist,mindist,PVFdomint
      REAL(8), ALLOCATABLE :: NpartCell(:),strho(:)

      INTEGER i,j,jj,node,k



      ficrtica=inp%nanoVsh * inp%nanovf
      fi0=inp%nanovf - ficrtica
      sigma=inp%nanoSigma

      PartVolFrac=0.0D0

      IF (inp%nanoMod.EQ.1) THEN

c       distribute particle to mesh
        DO i=1,inp%npart
          CALL GetPARTlc(mesh,part(i),Ni)
c         sum ni
          sumni=0.0D0
          DO j=1,mesh%npoc
            sumni=sumni+Ni(j)**sigma
          END DO
c         distribute
          DO j=1,mesh%npoc
            node = mesh%idc(part(i)%cell,j)
            PartVolFrac(node)=PartVolFrac(node)+Ni(j)**sigma/sumni
          END DO
        END DO

c       local particle density
        DO i=1,mesh%nnodes
          PartVolFrac(i)=PartVolFrac(i)/mesh%nodeVol(i)
        END DO

c       domain integral of local particle density
        CALL DomainIntegral(mesh,gp,PartVolFrac,PVFdomint)

c       final change
        DO i=1,mesh%nnodes
          PartVolFrac(i)=fi0+ficrtica*PartVolFrac(i)/PVFdomint*mesh%MeshVolume
        END DO

c       final change
c        DO i=1,mesh%nnodes
c          PartVolFrac(i)=fi0+ficrtica*DBLE(mesh%nnodes)/DBLE(inp%npart)*PartVolFrac(i)
c        END DO

      ELSE IF (inp%nanoMod.EQ.2) THEN
C       closest node
        DO i=1,inp%npart
          CALL GetPARTlc(mesh,part(i),Ni)

          mindist=1.0D10
c         distribute
          DO j=1,mesh%npoc
            node = mesh%idc(part(i)%cell,j)
            dist=SQRT( (part(i)%x - mesh%x(node,1))**2.0D0 +
     &                 (part(i)%y - mesh%x(node,2))**2.0D0 +
     &                 (part(i)%z - mesh%x(node,3))**2.0D0 )
            IF (dist.LT.mindist) THEN
              mindist=dist
              jj=j
            END IF
          END DO
          node = mesh%idc(part(i)%cell,jj)
          PartVolFrac(node)=PartVolFrac(node)+1.0D0 ! contribution to only one node
        END DO

c       local particle density
c        DO i=1,mesh%nnodes
c          PartVolFrac(i)=PartVolFrac(i)/mesh%nodeVol(i)
c        END DO

c       domain integral of local particle density
        CALL DomainIntegral(mesh,gp,PartVolFrac,PVFdomint)

c       final change
        DO i=1,mesh%nnodes
          PartVolFrac(i)=fi0+ficrtica*PartVolFrac(i)/PVFdomint*mesh%MeshVolume
        END DO


      ELSE IF (inp%nanoMod.EQ.3) Then

C       Prestejemo delce v celicah okoli nodea
        DO i=1,inp%npart
          j=part(i)%cell
          DO k=1,mesh%npoc
            PartVolFrac(mesh%idc(j,k))=PartVolFrac(mesh%idc(j,k))+1.0D0
          END DO
        END DO

c       spremenimo stdel/m3 v volumski delez
c       tako da lokalno konc. delcev delimo s povprečno koncentracijo delcev
c       (Nlok/Vlok) / (Ntot/Vtot)
        DO i=1,mesh%nnodes
          PartVolFrac(i)=fi0+ficrtica*PartVolFrac(i)/mesh%nodeVol(i)*mesh%MeshVolume/DBLE(inp%npart)
        END DO


      ELSE
        CALL WarnErr(env,io,inp,4,"CalParticleVolumeFraction","Wrong model!",0)
      END IF


c     check
c      print *,PartVolFrac
c      CALL DomainIntegral(mesh,gp,PartVolFrac,PVFdomint)
c      print *,PVFdomint/mesh%MeshVolume

c      sumni=0.0D0
c      DO i=1,mesh%nnodes
c        sumni=sumni+PartVolFrac(i)
c        print *,PartVolFrac(i)
c      END DO
c      print *,sumni/DBLE(mesh%nnodes)
c      print *,mesh%MeshVolume,ficrtica,fi0
c      stop



      END

C -----------------------------------------------------------------------------
      SUBROUTINE CalNodeVolume(mesh)
C
C     Calculate volume around each node
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType)  :: mesh

      INTEGER i,j,k
      LOGICAL je

      ALLOCATE (mesh%nodeVol(mesh%nnodes))

      mesh%nodeVolTot=0.0D0
      mesh%nodeVol=0.0D0

      DO i=1,mesh%nnodes
        DO j=1,mesh%nicell
          je=.FALSE.
          DO k=1,mesh%npoc
            IF (i.EQ.mesh%idc(j,k)) je=.TRUE.
          END DO
          IF (je) mesh%nodeVol(i)=mesh%nodeVol(i)+mesh%CellVolume(j)
        END DO
        mesh%nodeVolTot=mesh%nodeVolTot+mesh%nodeVol(i)
      END DO

      END



C -----------------------------------------------------------------------------
      SUBROUTINE CavityBC(part)
C
C     $: Elasticen odboj od vseh sten
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(particleType) :: part

c     spodaj
      IF (part%z.LT.0.0D0) THEN
        part%z=-part%z
        part%vz=-part%vz
      END IF

c     zgoraj
      IF (part%z.GT.1.0D0) THEN
        part%z=2.0D0-part%z
        part%vz=-part%vz
      END IF

c     levo
      IF (part%x.LT.0.0D0) THEN
        part%x=-part%x
        part%vx=-part%vx
      END IF

c     desno
      IF (part%x.GT.1.0D0) THEN
        part%x=2.0D0-part%x
        part%vx=-part%vx
      END IF

c     spredaj
      IF (part%y.LT.0.0D0) THEN
        part%y=-part%y
        part%vy=-part%vy
      END IF

c     zadaj
      IF (part%y.GT.1.0D0) THEN
        part%y=2.0D0-part%y
        part%vy=-part%vy
      END IF

      END


C -----------------------------------------------------------------------------
      SUBROUTINE GetPARTlc(mesh,part,fi)
C
C     $: get shape function based on location of the particle in cell
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType)  :: mesh
      TYPE(particleType) :: part ! delec

      REAL(8) fi(mesh%npoc) ! shape functions

      REAL(8) x,y,z,xi,eta,zeta
      INTEGER ierr,ipart,cell,i,j,node

      x=part%x
      y=part%y
      z=part%z
c     init
      xi=0.0D0
      eta=0.0D0
      zeta=0.0D0


c     naprej pogledamo ce je ostal v isti celici
c     is the particle located in the same mesh cell?
      cell=part%cell
      CALL tcfp_cross(mesh,cell,x,y,z,ierr)
      IF (ierr.EQ.0) GOTO 40

c     nato pogledamo v sosednje celice
c     it the particle in neighbouring cells?
      DO i=1,mesh%neil(part%cell,27)
        cell=mesh%neil(part%cell,i)
        CALL tcfp_cross(mesh,cell,x,y,z,ierr)
        IF (ierr.EQ.0) GOTO 40
      END DO

c     nazadnje gremo cez vse (do tega ne bi smelo nikoli priti)
c     oziroma to se zgodi takrat, ki je delec izven obmocja
c     finally, we check each cell
      DO cell=1,mesh%nicell
        CALL tcfp_cross(mesh,cell,x,y,z,ierr)
        IF (ierr.EQ.0) GOTO 40
      END DO

C     Nisem ga nasel, torej ga ni nikjer v mrezi -> sel v zid??
c     Particle was not found in any mesh cell
      RETURN

40    CONTINUE
c     naprej pogledamo ce je v celici, ki smo jo nasli zgoraj
c     calculate local particle coordinates (xi,eta,zeta) in mesh cell
      CALL tcfp(mesh,cell,xi,eta,zeta,x,y,z,ierr)
      IF (ierr.EQ.0) GOTO 100

C     Nisem ga nasel, torej ga ni nikjer v mrezi -> sel v zid??
c     Paticle was not found
      RETURN

100   CONTINUE
c     delec je bil najden v celici "cell" s koordiantami xi,eta,zeta
c     particle was found in CELL having local coordinates XI,ETA,ZETA
      part%cell=cell
c     calculate shape functions
      CALL cshape27(xi,eta,zeta,fi,mesh%npoc)

      END


C -----------------------------------------------------------------------------
      SUBROUTINE GetFFFap(mesh,velocity,vorticity,temperature,gradT,fap,part,ierr)
C
C     $: get fluid properties at the location of the ipart-th particle
C        ierr=0 particle was found, fluid properties are in fap
C        ierr=1 particles was NOT found anywhere in the mesh
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType)  :: mesh
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3)
      REAL(8) temperature(mesh%nnodes)
      REAL(8) gradT(mesh%nnodes,3)
      TYPE(fluidAtPartType) :: fap ! tekocina na mestu delca

      TYPE(particleType) :: part ! delec


      REAL(8) x,y,z,xi,eta,zeta
      REAL(8) fi(mesh%npoc) ! shape functions
      INTEGER ierr,ipart,cell,i,j,node

      x=part%x
      y=part%y
      z=part%z
c     init flow field at the location of the particle
      fap%vx=0.0D0
      fap%vy=0.0D0
      fap%vz=0.0D0
      fap%wx=0.0D0
      fap%wy=0.0D0
      fap%wz=0.0D0
      fap%T =0.0D0
      fap%gradVx=0.0D0
      fap%gradVy=0.0D0
      fap%gradVz=0.0D0
      fap%dvxdt=0.0D0
      fap%dvydt=0.0D0
      fap%dvzdt=0.0D0
      fap%gradT(1)=0.0D0
      fap%gradT(2)=0.0D0
      fap%gradT(3)=0.0D0
c     init
      xi=0.0D0
      eta=0.0D0
      zeta=0.0D0


c     naprej pogledamo ce je ostal v isti celici
c     is the particle located in the same mesh cell?
      cell=part%cell
      CALL tcfp_cross(mesh,cell,x,y,z,ierr)
      IF (ierr.EQ.0) GOTO 40

c     nato pogledamo v sosednje celice
c     it the particle in neighbouring cells?
      DO i=1,mesh%neil(part%cell,27)
        cell=mesh%neil(part%cell,i)
        CALL tcfp_cross(mesh,cell,x,y,z,ierr)
        IF (ierr.EQ.0) GOTO 40
      END DO

c     nazadnje gremo cez vse (do tega ne bi smelo nikoli priti)
c     oziroma to se zgodi takrat, ki je delec izven obmocja
c     finally, we check each cell
      DO cell=1,mesh%nicell
        CALL tcfp_cross(mesh,cell,x,y,z,ierr)
        IF (ierr.EQ.0) GOTO 40
      END DO

C     Nisem ga nasel, torej ga ni nikjer v mrezi -> sel v zid??
c     Particle was not found in any mesh cell
      RETURN

40    CONTINUE
c     remember the cell number in which the particle was found
      part%cell=cell

c     naprej pogledamo ce je v celici, ki smo jo nasli zgoraj
c     calculate local particle coordinates (xi,eta,zeta) in mesh cell
      cell=part%cell
      CALL tcfp(mesh,cell,xi,eta,zeta,x,y,z,ierr)
      IF (ierr.EQ.0) GOTO 100

C     Nisem ga nasel, torej ga ni nikjer v mrezi -> sel v zid??
c     Paticle was not found
      RETURN

100   CONTINUE
c     delec je bil najden v celici "cell" s koordiantami xi,eta,zeta
c     particle was found in CELL having local coordinates XI,ETA,ZETA

c     calculate shape functions
      CALL cshape27(xi,eta,zeta,fi,mesh%npoc)

c     interpoliram lastnosi tekocine na mestu delca v fap
c     interpolate fluid parameters
      DO i=1,mesh%npoc
        node=mesh%idc(cell,i)
        fap%vx=fap%vx+velocity(node,1)*fi(i)
        fap%vy=fap%vy+velocity(node,2)*fi(i)
        fap%vz=fap%vz+velocity(node,3)*fi(i)
        fap%wx=fap%wx+vorticity(node,1)*fi(i)
        fap%wy=fap%wy+vorticity(node,2)*fi(i)
        fap%wz=fap%wz+vorticity(node,3)*fi(i)
        fap%T =fap%T+ temperature(node)*fi(i)
        fap%gradT(1)=fap%gradT(1)+gradT(node,1)*fi(i)
        fap%gradT(2)=fap%gradT(2)+gradT(node,2)*fi(i)
        fap%gradT(3)=fap%gradT(3)+gradT(node,3)*fi(i)
      END DO



      END


C -----------------------------------------------------------------------------
      SUBROUTINE SetUpParticles(inp,part,mesh)
C
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType)  :: mesh
      TYPE(InputType) inp
      TYPE(particleType) part(inp%npart)

      REAL(8) xi,eta,zeta,x,y,z
      INTEGER rndcell,i

c      integer bla(mesh%nicell)
c      real(8) v

      DO i=1,inp%npart
        CALL rndCellWei(mesh%CellVolume,rndcell,mesh%nicell) ! izbere nakljucno celico
        CALL rndXiEtaZeta(xi,eta,zeta)   ! izbere nakljucno lokacijo znotraj celice
        CALL ppos(mesh,rndcell,xi,eta,zeta,x,y,z) ! izracuna lokacijo delca

        part(i)%id=i
        part(i)%x=x
        part(i)%y=y
        part(i)%z=z
        part(i)%cell=rndcell
        part(i)%active=.TRUE.
        part(i)%vx=0.0D0
        part(i)%vy=0.0D0
        part(i)%vz=0.0D0

      END DO

c      bla=0
c      do i=1,inp%npart
c        bla(part(i)%cell)=bla(part(i)%cell)+1
c      end do

c      v=0.0D0
c      do i=1,mesh%nicell
c        print *,DBLE(bla(i))/mesh%cellvolume(i)
c        v=v+DBLE(bla(i))/mesh%cellvolume(i)
c     end do

c      print *,"ee",1.0D0-v/DBLE(mesh%nicell)/DBLE(inp%npart)

c      stop

      END



C -----------------------------------------------------------------------------

      SUBROUTINE OutputParticles(io,inp,cnt,part)
C
C     $: Outputs function
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(IOtype) io
      TYPE(inputtype) inp
      TYPE(particleType) part(inp%npart)
      TYPE(countType) cnt


      INTEGER itp,nic,i
      CHARACTER*100 zonetitle,pfn
      CHARACTER*255 vrstica

      itp=io%pres

      Write(zonetitle,'(I5,A1,G12.7)') cnt%tstep,'-',cnt%rtime
c      zonetitle="u(x,y,z)"
C     First time step - requires header
      IF (cnt%ftime.EQ.1) THEN
        OPEN (itp,FILE=TRIM(io%pres_name),STATUS='UNKNOWN')

        WRITE (itp,'(A)') '# |---------------------------|'
        WRITE (itp,'(A)') '# |        T R I T O K        |'
        WRITE (itp,'(A)') '# |         PARTICLES         |'
        WRITE (itp,'(A)') '# |---------------------------|'

        WRITE (itp,'(A)') 'VARIABLES = "X", "Y", "Z","vx","vy","vz"'
        CALL GETARG(0,pfn)
        WRITE (itp,'(7A)') 'DATASETAUXDATA code = "',trim(inp%IDname),' v',trim(inp%IDversion),', ',trim(inp%IDdate),'"'
        WRITE (itp,'(3A)') 'DATASETAUXDATA program = "',trim(pfn),'"'
      END IF

      WRITE (itp,'(A,A,A)') 'ZONE T="',trim(zonetitle),'"'
      CALL ZoneAuxData(itp,cnt%tstep,cnt%rtime)

C     Output results
      DO  i=1,inp%npart
        IF (part(i)%active) THEN
          WRITE (vrstica,'(6G18.9)') part(i)%x,part(i)%y,part(i)%z,part(i)%vx,part(i)%vy,part(i)%vz
          CALL sqblnk(itp,vrstica)
        END IF
      END DO

      END




C -----------------------------------------------------------------------------
      SUBROUTINE CoR_NeighbourList(env,mesh,io,inp)
C
C     $: Preberi ali izracunan na novo seznam sosedov
C        Read from file (.neil) or Calculate a list of neighbouring cells
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType)  :: mesh
      TYPE(IOtype)    :: io
      TYPE(inputtype) :: inp
      TYPE(penv) :: env

      INTEGER iok,dummy

c     maksimalno stevilo sosedov ene celice je 26
c     v 27-ti stolpec shranim stevilo sostedov
      ALLOCATE (mesh%neil(mesh%nicell,27))

      CALL CheckNeilFile(io,mesh%nnodes,mesh%nbnodes,mesh%nicell,iok)

      IF (iok.EQ.1) THEN
c       Read neil from disk
        IF (env%myproc.EQ.1) CALL WarnErr(env,io,inp,0,"CoR_NeighbourList","Reading neighbour list!",0)
        OPEN (io%neil,FILE=TRIM(io%neil_name),FORM='UNFORMATTED',STATUS='OLD')
        READ(io%neil) dummy,dummy,dummy
        CALL rdIMat(mesh%neil,mesh%nicell,27,io%neil)
        CLOSE(io%neil)
      ELSE
c       Calculate neil
        IF (env%myproc.EQ.1) CALL WarnErr(env,io,inp,0,"CoR_NeighbourList","Making neighbour list!",0)
        CALL MakeNeighList(mesh)
c       write neil to disk
        IF (env%myproc.EQ.1) THEN
          OPEN (io%neil,FILE=TRIM(io%neil_name),FORM='UNFORMATTED',STATUS='UNKNOWN')
          WRITE(io%neil) mesh%nnodes,mesh%nbnodes,mesh%nicell
          CALL WrIMat(mesh%neil,mesh%nicell,27,io%neil)
          CLOSE(io%neil)
        END IF
      END IF



      END

C -----------------------------------------------------------------------------
      SUBROUTINE CheckNeilFile(io,a,b,c,iok)
C
C     $: Preveri ce je Neil na disku za pravo mrezo
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(IOtype) io

      INTEGER a,b,c
      INTEGER aa,bb,cc,iok

      iok=0
      OPEN(io%neil,FILE=TRIM(io%neil_name),FORM='UNFORMATTED',STATUS='OLD',ERR=10)
      READ(io%neil) aa,bb,cc

      IF (a.EQ.aa.AND.b.EQ.bb.AND.c.EQ.cc) iok=1
      CLOSE(io%neil)

10    RETURN
      END



C -----------------------------------------------------------------------------
      SUBROUTINE MakeNeighList(mesh)
C
C     $: naredi seznam sosedov vsake celice
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh

      INTEGER i,j,k,l
      LOGICAL jeSosed,imamNaSeznamu


      mesh%neil=0

      DO i=1,mesh%nicell
        DO j=1,mesh%nicell
          IF (i.NE.j) THEN
c           celica j je soseda celice i, ce je katerikoli node v obeh celicah
            jeSosed=.FALSE.
            DO k=1,mesh%npoc
              DO l=1,mesh%npoc
                IF (mesh%idc(i,k).EQ.mesh%idc(j,l)) THEN
                  jeSosed=.TRUE.
                  GOTO 10
                END IF
              END DO  !l
            END DO  !k
10          CONTINUE
c            print *,i,j,k,l,jesosed
            IF (jeSosed) THEN
c             poiscem, ce ga ze imam na seznamu
              imamNaSeznamu=.FALSE.
              DO l=1,mesh%neil(i,27)
                IF (mesh%neil(i,l).EQ.j) THEN
                  imamNaSeznamu=.TRUE.
                  GOTO 20
                END IF
              END DO
20            CONTINUE
              IF (imamNaSeznamu.EQV..FALSE.) THEN
                mesh%neil(i,27)=mesh%neil(i,27)+1
                mesh%neil(i,mesh%neil(i,27))=j
              END IF
            END IF
          END IF
        END DO  !j
      END DO !i


      END


C______________________________________________________________________C
C______________________________________________________________________C
      SUBROUTINE WrIMat(mat,nrow,ncol,io)
C        __    _      ___
C        Write Single Matrix
C______________________________________________________________________C
C______________________________________________________________________C
      INTEGER nrow,ncol,io,j
      INTEGER mat(nrow,ncol)

      DO j=1,ncol
        CALL wrIvec(io,nrow,mat(1,j))
      END DO

      END
C______________________________________________________________________C
C______________________________________________________________________C
      SUBROUTINE RdIMat(mat,nrow,ncol,io)
C        _  _ _      ___
C        Read Single Matrix
C______________________________________________________________________C
C______________________________________________________________________C
      INTEGER nrow,ncol,io,j
      INTEGER mat(nrow,ncol)

      DO j=1,ncol
        CALL rdIvec(io,nrow,mat(1,j))
      END DO

      END

C----------------------------------------------------------------------C
c----------------------------------------------------------------------c
c                                                                      c
      SUBROUTINE rdIvec(ifr,nnx,vec)
c                                                                      c
c----------------------------------------------------------------------c
C----------------------------------------------------------------------C
c.......................................................................
c..                                                                   ..
c..   REad VECtor using MAXSIZE chunks                                ..
c..   --   ---                                                        ..
c.......................................................................
      INTEGER ifr,nnx,maxsize,nblock,i,j,k
      INTEGER  vec(nnx)
      PARAMETER (maxsize=8192)

      nblock=INT(nnx/maxsize)
      DO j=1,nblock
        k=maxsize*(j-1)
        READ(ifr) (vec(i),i=k+1,k+maxsize)
      END DO
      k=maxsize*nblock
      READ(ifr) (vec(i),i=k+1,nnx)
      RETURN
      END
C----------------------------------------------------------------------C
c----------------------------------------------------------------------c
c                                                                      c
      SUBROUTINE wrIvec(ifr,nnx,vec)
c                                                                      c
c----------------------------------------------------------------------c
C----------------------------------------------------------------------C
c.......................................................................
c..                                                                   ..
c..   WRite VECtor using MAXSIZE chunks                               ..
c..   --    ---                                                       ..
c.......................................................................
      INTEGER ifr,nnx,maxsize,nblock,i,j,k
      INTEGER  vec(nnx)
      PARAMETER (maxsize=8192)
c....
      nblock=INT(nnx/maxsize)
      DO j=1,nblock
        k=maxsize*(j-1)
        WRITE(ifr) (vec(i),i=k+1,k+maxsize)
      END DO
      k=maxsize*nblock
      WRITE(ifr) (vec(i),i=k+1,nnx)
      RETURN
      END



C -----------------------------------------------------------------------------
      SUBROUTINE rndCellWei(utezi,stevilo,mxst)
C
C     $: vrne nakljucno naravno stevilo med 1 in maxst
C        v zacetku programa daj CALL RANDOM_SEED
C        uposteva utezi
C        vsota utezi mora biti 1
C
C -----------------------------------------------------------------------------
      INTEGER stevilo,mxst,i
      REAL(8) x,utezi(mxst),vsota

      CALL RANDOM_NUMBER(x)

      vsota=0.0D0

      DO i=1,mxst
        vsota=vsota+utezi(i)
        IF (vsota.GE.x) THEN
          stevilo=i
          EXIT
        END IF
      END DO



      IF (stevilo.GT.mxst.OR.stevilo.LT.1) PRINT *,"ERROR"

      END



C -----------------------------------------------------------------------------
      SUBROUTINE rndIno(stevilo,mxst)
C
C     $: vrne nakljucno naravno stevilo med 1 in maxst
C        v zacetku programa daj CALL RANDOM_SEED
C
C -----------------------------------------------------------------------------
      REAL(8) x
      INTEGER stevilo,mxst

      CALL RANDOM_NUMBER(x)
      stevilo=INT(x*mxst)+1

      END

C -----------------------------------------------------------------------------
      SUBROUTINE rndXiEtaZeta(xi,eta,zeta)
C
C     $: vrne 3 nakljucna stevila med -1 in 1
C        v zacetku programa daj CALL RANDOM_SEED
C
C -----------------------------------------------------------------------------
      REAL(8) xi,eta,zeta,x

      CALL RANDOM_NUMBER(x)
      xi=-1.0D0+2.0D0*x

      CALL RANDOM_NUMBER(x)
      eta=-1.0D0+2.0D0*x

      CALL RANDOM_NUMBER(x)
      zeta=-1.0D0+2.0D0*x


      END


C -----------------------------------------------------------------------------
      SUBROUTINE ppos(mesh,ic,xi,eta,zeta,x,y,z)
C
C     $: izracuna x,y,z iz stevilke celice in xi,eta,zeta
C        calculate x,y,z from xi,eta,zeta
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType)  :: mesh
      REAL(8) xi,eta,zeta,x,y,z
      INTEGER ic

      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8) x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8


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

      CALL lks2kks(xi,eta,zeta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                         x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x,y,z)


      END


C -----------------------------------------------------------------------------
      SUBROUTINE tcfp_cross(mesh,ic,x,y,z,ierr)
C
C     $: preveri ali je v tej celici delec
C        check if particle is in a mesh cell (ic), using cross products
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType)  :: mesh
      REAL(8) x,y,z
      INTEGER ic,ierr

      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8) x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8


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


      CALL cipic(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                   x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,ierr)


      END



C -----------------------------------------------------------------------------
      SUBROUTINE tcfp(mesh,ic,xi,eta,zeta,x,y,z,ierr)
C
C     $: preveri ali je v tej celici delec
C        find local coordinates of the particle (xi,eta,zeta)
C        if particle is not in cell "ic", ierr not zero
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType)  :: mesh
      REAL(8) xi,eta,zeta,x,y,z
      INTEGER ic,ierr

      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8) x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8


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


      CALL kks2lks(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                   x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,xi,eta,zeta,ierr)

      END
