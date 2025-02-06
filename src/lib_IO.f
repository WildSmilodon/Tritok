C -----------------------------------------------------------------------------
      SUBROUTINE ReadInputFile(env,mesh,io,inp,gauss)
C
C     $: reads bem.inp
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(meshType) :: mesh
      TYPE(InputType) inp
      TYPE(IOtype) io
      TYPE(gausstype) :: gauss
      TYPE(penv) :: env
      
      INTEGER lun
      CHARACTER KeyWord*4,OneLine*255,dummy*64,tmp*255

      lun=io%inp
      OPEN (lun,FILE=io%inp_name,ERR=10,STATUS='OLD') !,SHARED)

C
C***    Set up default settings:
c
      gauss%kmBs=8 ! number of gauss points for KineMatics Boundary Singular integrals (8=48)
      gauss%kmBr=8 ! numbee of gauss points for KineMatics Boundary Regular integrals (8=48)
      inp%copy=0  ! do not copy
      inp%nnlit=1000
      inp%ur=0.01D0
      inp%eps=1.0D-7
      inp%ReN=1.0D0
      inp%RaN=1.0D0
      inp%PrN=1.0D0
      inp%PeN=-1.0D0
      inp%ScN=1.0D0   ! Pe=Re*Sc
      inp%Rro=1.0D20  ! Density stability ratio
      inp%gx=0.0D0
      inp%gy=0.0D0
      inp%gz=-1.0D0

c     least squares solver      
      inp%lsqs_maxit=5000
      inp%lsqs_eps=1.0E-6
      inp%idlse=0
      inp%dlse_r=100.0
      inp%dlse_epsth=1.0E-3

c     square solver
      inp%sqrs_type=2 ! rbicstab
      inp%sqrs_prec=2 ! ilu
      inp%sqrs_maxit=500
      inp%sqrs_eps=1.0E-6

c     parallel kinematics solver
      inp%pkms_maxit=500
      inp%pkms_eps=1.0E-6

c     underrelaxation defaults
      inp%urdv(1)=0.1D0
      inp%urdv(2)=0.1D0
      inp%urdv(3)=0.1D0
      inp%urdw(1)=0.1D0
      inp%urdw(2)=0.1D0
      inp%urdw(3)=0.1D0
      inp%urbw(1)=0.1D0
      inp%urbw(2)=0.1D0
      inp%urbw(3)=0.1D0      
      inp%urDT=0.1D0
      inp%urDC=0.1D0          
      inp%urmH=1.0D0      

      inp%inano=0 ! not a nanofluid
      inp%iDWUR=0 ! under-relaxation larger at the walls

      inp%imixc=0 ! not a mixture of forced and natural convection problem

      inp%iestr=0 ! do not only estimate RAM

c     not a ferrofluid !!
      inp%iFerro=0
      inp%Ram=0.0D0
      inp%bdt=1.0D0
      inp%hi0=0.0D0

c     izvoz v paraview (vtu,vtk)
      inp%iPara=0
      inp%iParaVTK=0

c     izvoz rezultatov po stenah
      inp%iRawl=1


c     energy equation model (for iDT=3)
      inp%eeqm=0

c     particle tracking
      inp%iPart=0
      inp%nPart=0
      inp%PartTSR=1

c     nano
      inp%nanoVsh=0.1D0
      inp%nanoSigma=1.0D0
      inp%nanoMod=1

c     Browninan motion
      inp%iBrown=0  ! 1/0
      inp%BrownT0=293.0D0
      inp%BrownDomainSize=1.0D0
      inp%BrownPartDiam=50D-9
      inp%BrownT0mul=0.0D0 ! must be zero !!!

C       Thermophoresis
      inp%iThph=0 ! 1/0
      inp%thphA=1264.0 D0
      inp%thphB=1.417D0
      inp%thphPartDiam=50.0D-9
      inp%thphMul=0.0D0 ! must be zero !!!
      inp%thphT0=293.0D0
      inp%thphDomainSize=1.0D0

C     Mesh transromation
      inp%iMTrot=0
      inp%iMTstr=0
      inp%iMTtra=0
      inp%MTrot=0.0D0
      inp%MTstr=0.0D0
      inp%MTtra=0.0D0

C     Vorticity equation material model
      inp%veqm=0

C     Internal integral divisions
      gauss%iiDiv=1

c       modified Helmholz equation solver
c       nabla^2 u - mu^2 u = f
      inp%mHmu = 1.0D0
      inp%mHf = 1

c     wawelet compression
      inp%WT_kappa = 1.0D-5
c     ACA compression
      inp%ACA_eps = 0.5
      inp%ACA_type = 1
      inp%ctree=""

C     cluster trees
      inp%ct_type=1 ! no compression
      inp%ct_eta=1.0D0 ! admisiblity criterion

c     test DC equation
      inp%iTest=0

c     time scheme
      inp%TimeScheme=1

c     export stohastic data
      inp%estd=0

c     analytic channel flow
      inp%iacha=0

      CALL rOneTL(lun,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')
        READ (OneLine,*) KeyWord
C
C***    GEO AND BC FILE LOCATION :
c
        IF (KeyWord.EQ.'MDIR') THEN
          READ(OneLine,'(A5,A)') dummy,mesh%path

        ELSE IF (KeyWord.EQ.'LGEO') THEN
          READ(OneLine,'(A5,A)') dummy,mesh%filename

        ELSE IF (KeyWord.EQ.'LBIC') THEN
          READ(OneLine,'(A5,A)') dummy,inp%lbic

        ELSE IF (KeyWord.EQ.'SCRO') THEN
          READ(OneLine,*) dummy,inp%scro
          
        ELSE IF (KeyWord.EQ.'TIME') THEN
          READ(OneLine,*) dummy,inp%tstep,inp%nstep,inp%iwrit
          IF (inp%iwrit.GT.inp%nstep) inp%iwrit=inp%nstep          

        ELSE IF (KeyWord.EQ.'TIMS') THEN
          READ(OneLine,*) dummy,inp%TimeScheme

        ELSE IF (KeyWord.EQ.'PROB') THEN
          READ(OneLine,*) dummy,inp%prob

        ELSE IF (KeyWord.EQ.'ESTD') THEN
          READ(OneLine,*) dummy,inp%estd
          
        ELSE IF (KeyWord.EQ.'REST') THEN
          READ(OneLine,*) dummy,inp%rst          

        ELSE IF (KeyWord.EQ.'PART') THEN
          READ(OneLine,*) dummy,inp%iPart,inp%nPart,inp%PartTSR

        ELSE IF (KeyWord.EQ.'LSQS') THEN
          READ(OneLine,*) dummy,inp%lsqs_maxit,inp%lsqs_eps

        ELSE IF (KeyWord.EQ.'SQRS') THEN
          READ(OneLine,*) dummy,inp%sqrs_type,inp%sqrs_prec,inp%sqrs_maxit,inp%sqrs_eps

        ELSE IF (KeyWord.EQ.'MTRO') THEN
          READ(OneLine,*) dummy,inp%iMTrot,inp%MTrot(1),inp%MTrot(2),inp%MTrot(3)
        ELSE IF (KeyWord.EQ.'MTST') THEN
          READ(OneLine,*) dummy,inp%iMTstr,inp%MTstr(1),inp%MTstr(2),inp%MTstr(3)
        ELSE IF (KeyWord.EQ.'MTTR') THEN
          READ(OneLine,*) dummy,inp%iMTtra,inp%MTtra(1),inp%MTtra(2),inp%MTtra(3)

        ELSE IF (KeyWord.EQ.'DLSE') THEN
          READ(OneLine,*) dummy,inp%iDLSE,inp%dlse_r,inp%dlse_epsth          

        ELSE IF (KeyWord.EQ.'DWUR') THEN
          READ(OneLine,*) dummy,inp%iDWUR,inp%dwur_max,inp%dwur_a

        ELSE IF (KeyWord.EQ.'PKMS') THEN ! parallel kinematics solver 212
          READ(OneLine,*) dummy,inp%pkms_maxit,inp%pkms_eps

        ELSE IF (KeyWord.EQ.'EQNS') THEN 
          READ(OneLine,*) dummy,inp%iBw,inp%iDv,inp%iDw,inp%iDT,inp%iDC,inp%imH

        ELSE IF (KeyWord.EQ.'GRAV') THEN 
          READ(OneLine,*) dummy,inp%gx,inp%gy,inp%gz

        ELSE IF (KeyWord.EQ.'REYN') THEN 
          READ(OneLine,*) dummy,inp%ReN          

        ELSE IF (KeyWord.EQ.'PECL') THEN
          READ(OneLine,*) dummy,inp%PeN

        ELSE IF (KeyWord.EQ.'MHEL') THEN
          READ(OneLine,*) dummy,inp%mHmu,inp%mhf ! mu and f

        ELSE IF (KeyWord.EQ.'SRRO') THEN
          READ(OneLine,*) dummy,inp%ScN,inp%Rro

        ELSE IF (KeyWord.EQ.'RAPR') THEN 
          READ(OneLine,*) dummy,inp%RaN,inp%PrN     

        ELSE IF (KeyWord.EQ.'MIXC') THEN
          READ(OneLine,*) dummy,inp%imixc

        ELSE IF (KeyWord.EQ.'ACHA') THEN
          READ(OneLine,*) dummy,inp%iacha

        ELSE IF (KeyWord.EQ.'ESTR') THEN
          READ(OneLine,*) dummy,inp%iestr

        ELSE IF (KeyWord.EQ.'TEST') THEN
          READ(OneLine,*) dummy,inp%iTest ! test DC equation

        ELSE IF (KeyWord.EQ.'FERR') THEN
          READ(OneLine,*) dummy,inp%iFerro,inp%Ram,inp%bdt,inp%hi0

        ELSE IF (KeyWord.EQ.'URBW') THEN 
          READ(OneLine,*) dummy,inp%urbw(1),inp%urbw(2),inp%urbw(3)

        ELSE IF (KeyWord.EQ.'URDV') THEN 
          READ(OneLine,*) dummy,inp%urdv(1),inp%urdv(2),inp%urdv(3)
          
        ELSE IF (KeyWord.EQ.'URDW') THEN 
          READ(OneLine,*) dummy,inp%urdw(1),inp%urdw(2),inp%urdw(3)

        ELSE IF (KeyWord.EQ.'URDT') THEN 
          READ(OneLine,*) dummy,inp%urDT

        ELSE IF (KeyWord.EQ.'URMH') THEN 
          READ(OneLine,*) dummy,inp%urmh

        ELSE IF (KeyWord.EQ.'URDC') THEN 
          READ(OneLine,*) dummy,inp%urDC

        ELSE IF (KeyWord.EQ.'PARA') THEN
          READ(OneLine,*) dummy,inp%iPara,inp%iParaVTK

        ELSE IF (KeyWord.EQ.'WAWT') THEN
          READ(OneLine,*) dummy,inp%WT_kappa

        ELSE IF (KeyWord.EQ.'CTRT') THEN
          READ(OneLine,*) dummy,inp%ct_type,inp%ct_eta

        ELSE IF (KeyWord.EQ.'ACAC') THEN
          READ(OneLine,*) dummy,inp%ACA_type,inp%ACA_eps

        ELSE IF (KeyWord.EQ.'RAWL') THEN
          READ(OneLine,*) dummy,inp%iRawl

        ELSE IF (KeyWord.EQ.'EEQM') THEN ! energy equation model
          READ(OneLine,*) dummy,inp%eeqm

        ELSE IF (KeyWord.EQ.'VEQM') THEN ! vorticity equation material model
          READ(OneLine,*) dummy,inp%veqm

        ELSE IF (KeyWord.EQ.'NANO') THEN
          READ(OneLine,*) dummy,inp%inano,inp%nanovf,
     &                                    inp%nanoF_cp,inp%nanoF_rho,inp%nanoF_k,inp%nanoF_beta,
     &                                    inp%nanoS_cp,inp%nanoS_rho,inp%nanoS_k,inp%nanoS_beta

        ELSE IF (KeyWord.EQ.'NANP') THEN
          READ(OneLine,*) dummy,inp%nanoMod,inp%nanoVsh,inp%nanoSigma

        ELSE IF (KeyWord.EQ.'BROW') THEN
          READ(OneLine,*) dummy,inp%iBrown,inp%BrownT0,inp%BrownDomainSize,inp%BrownPartDiam

        ELSE IF (KeyWord.EQ.'THPH') THEN
          READ(OneLine,*) dummy,inp%iThph,inp%thphA,inp%thphB,inp%thphT0,inp%thphDomainSize,inp%thphPartDiam

        ELSE IF (KeyWord.EQ.'RELA') THEN ! max number of nonlinear itera, underrelax, epsilon
          READ(OneLine,*) dummy,inp%nnlit,inp%ur,inp%eps

        ELSE IF (KeyWord.EQ.'COPY') THEN  ! kopira integrale (vse celice enake)
          READ(OneLine,*) dummy,inp%copy

        ELSE IF (KeyWord.EQ.'CTRE') THEN  ! cluster tree filename
          READ(OneLine,'(A5,A)') dummy,inp%ctree
C
C***    NUMBER OF GAUSS POINTS FOR INTEGRATION :
c
        ELSE IF (KeyWord.EQ.'GAUS') THEN
          READ(OneLine,*) dummy,gauss%kmBr,gauss%kmBs,gauss%kmDr,gauss%kmDs
        ELSE IF (KeyWord.EQ.'IIDI') THEN
          READ(OneLine,*) dummy,gauss%iiDiv
        END IF
        CALL rOneTL(lun,OneLine)
      END DO

      CLOSE (lun)
      IF (inp%iDLSE.EQ.1) THEN
        inp%dlse=inp%dlse_epsth
      ELSE
        inp%dlse=inp%lsqs_eps
      END IF

      mesh%fullname=Trim(mesh%path)//Trim(mesh%filename)  

      inp%mdir=mesh%path
      inp%lgeo=mesh%filename
      inp%lgeo_withpath=mesh%fullname
      inp%lbic_withpath=Trim(mesh%path)//Trim(inp%lbic)
      io%ctree_name=TRIM(mesh%path)//TRIM(inp%ctree)

C     nondimensionlization such that we suppose Re*Pr=1 (only natural convection problem)
      IF (inp%iDT.GT.0.AND.inp%imixc.EQ.0) THEN
        inp%reN=1.0D0/inp%PrN
      END IF

      IF (inp%PeN.GT.0.0D0) THEN
        inp%ScN=inp%PeN/inp%ReN
      END IF
     
      RETURN
      
10    continue ! error when opening input file
      WRITE (tmp,'(A,A)') "could not open : ",io%inp_name
      CALL WarnErr(env,io,inp,2,"ReadInputFile",trim(tmp),0)
      
      END






C______________________________________________________________________C
      SUBROUTINE ExportWallToTecplot(io,inp,cnt,mesh,v,w,T,dc,qT,qdc)
C
C     Naredi ploskovno mrezo po stenah in izpise rezultate
C
C______________________________________________________________________C
      USE inc_types

      TYPE(IOtype) io
      TYPE(meshType) mesh
      TYPE(inputtype) inp
      TYPE(countType) cnt

      REAL(8) v(mesh%nnodes,3)
      REAL(8) w(mesh%nnodes,3)
      REAL(8) T(mesh%nnodes)
      REAL(8) dc(mesh%nnodes)

      REAL(8) qT(mesh%nq)
      REAL(8) qdc(mesh%nq)

      REAL(8), ALLOCATABLE :: intqTemp(:,:),tmp(:)
      REAL(8), ALLOCATABLE :: intqdc(:,:)

      CHARACTER*255 vrstica,zonetitle

      INTEGER, ALLOCATABLE :: imon(:)

      INTEGER walln,nwbe,nwn,i,lun,j,ii,jj,node

      WRITE(zonetitle,'(I5,A1,G12.7)') cnt%tstep,'-',cnt%rtime

C     Interpoliram flukse na mrezo za funkcijo
      ALLOCATE (intqTemp(mesh%nnodes,mesh%nofw),tmp(mesh%nnodes))
      ALLOCATE (intqdc(mesh%nnodes,mesh%nofw))
      DO i=1,mesh%nofw
        CALL InterWallFluxtoFun(mesh,qT,tmp,i)
        DO j=1,mesh%nnodes
          intqTemp(j,i)=tmp(j)
        END DO
        CALL InterWallFluxtoFun(mesh,qdc,tmp,i)
        DO j=1,mesh%nnodes
          intqdc(j,i)=tmp(j)
        END DO
      END DO


      lun=io%rawl

      ALLOCATE (imon(mesh%nnodes))

C      Izberem zid
      DO walln=1,mesh%nofw

C     Prestejem stevilo robnih elementov v streni
      nwbe=0
      imon=0
      DO i=1,mesh%nbelem
        IF (mesh%bewn(i).EQ.walln) THEN ! sem na pravi steni
          nwbe=nwbe+1
          DO j=1,mesh%npob
            imon(mesh%ibc(i,j))=1
          END DO
        END IF
      END DO
C     In stevilo vozlisc v steni
      nwn=0
      DO i=1,mesh%nnodes
        IF (imon(i).EQ.1) nwn=nwn+1
      END DO

      IF (cnt%ftime.EQ.1) THEN
        WRITE (lun,'(A,I3,A1,A,A1,A,I7,A,I7,A)') 'ZONE T="wall=',walln,' ',TRIM(zonetitle),'"',
     & ', DATAPACKING=POINT, N=',nwn,' ,E= ',nwbe*4,'  ZONETYPE=FEQUADRILATERAL'
      ELSE
        WRITE (lun,'(A,I3,A1,A,A1,A,I7,A,I7,A,A,I2,A,I2)')
     & 'ZONE T="wall=',walln,' ',TRIM(zonetitle),'"',
     & ', DATAPACKING=POINT, N=',nwn,' ,E= ',nwbe*4,'  ZONETYPE=FEQUADRILATERAL',
     & ', VARSHARELIST=([1,2,3,4]=',walln,'), CONNECTIVITYSHAREZONE = ',walln
      END IF
      WRITE (lun,'(A,I2,A)') 'AUXDATA wall="',walln,'"'
      CALL ZoneAuxData(lun,cnt%tstep,cnt%rtime)

C     naredimo nov seznam vozlisc
      imon=0
      j=0

      DO ii=1,mesh%nbelem
        IF (mesh%bewn(ii).EQ.walln) THEN ! sem na pravi steni
          DO jj=1,mesh%npob
            node=mesh%ibc(ii,jj)
            IF (imon(node).EQ.0) THEN ! samo nove
              j=j+1
              imon(node)=j  ! obratno
              IF (cnt%ftime.EQ.1) THEN
                WRITE (vrstica,'(3G18.9,I8,10G18.9)')
     &          mesh%x(node,1),mesh%x(node,2),mesh%x(node,3),node,v(node,1),v(node,2),v(node,3),
     &          w(node,1),w(node,2),w(node,3),T(node),dc(node),intqTemp(node,walln),intqdc(node,walln)
              ELSE
                WRITE (vrstica,'(10G18.9)')
     &          v(node,1),v(node,2),v(node,3),
     &          w(node,1),w(node,2),w(node,3),T(node),dc(node),intqTemp(node,walln),intqdc(node,walln)
              END IF
              CALL sqblnk(lun,vrstica)
            END IF
          END DO
        END IF
      END DO


c     conectivity
      IF (cnt%ftime.EQ.1) THEN
        DO i=1,mesh%nbelem
          IF (mesh%bewn(i).EQ.walln) THEN ! sem na pravi steni
            WRITE(vrstica,*) imon(mesh%ibc(i,1)),imon(mesh%ibc(i,2)),imon(mesh%ibc(i,9)),imon(mesh%ibc(i,8))
            CALL sqblnk(lun,vrstica)
            WRITE(vrstica,*) imon(mesh%ibc(i,2)),imon(mesh%ibc(i,3)),imon(mesh%ibc(i,4)),imon(mesh%ibc(i,9))
            CALL sqblnk(lun,vrstica)
            WRITE(vrstica,*) imon(mesh%ibc(i,8)),imon(mesh%ibc(i,9)),imon(mesh%ibc(i,6)),imon(mesh%ibc(i,7))
            CALL sqblnk(lun,vrstica)
            WRITE(vrstica,*) imon(mesh%ibc(i,9)),imon(mesh%ibc(i,4)),imon(mesh%ibc(i,5)),imon(mesh%ibc(i,6))
            CALL sqblnk(lun,vrstica)
          END IF
        END DO
      END IF

c     konec zanke po stenah
      END DO

C     konec
      DEALLOCATE (imon,intqTemp,tmp,intqdc)

      END





C -----------------------------------------------------------------------------
      SUBROUTINE RrstF(env,mesh,io,inp,cnt,velocity,vorticity,
     &             ptsvorticity,pptsvorticity,qvorticity,pnlsvorticity,pnlsvelocity,
     &             temp,ptstemp,pptstemp,qtemp,pnlstemp,
     &             dc,ptsdc,pptsdc,qdc,pnlsdc)
C
C     $: reads restart file
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(IOtype) io
      TYPE(countType) cnt
      TYPE(InputType) :: inp
      TYPE(meshType) :: mesh
      TYPE(penv) :: env      
      
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) pnlsvelocity(mesh%nnodes,3)      
      REAL(8) vorticity(mesh%nnodes,3)      
      REAL(8) qvorticity(mesh%nq,3)      
      REAL(8) ptsvorticity(mesh%nnodes,3)                  
      REAL(8) pnlsvorticity(mesh%nnodes,3)      
      REAL(8) pptsvorticity(mesh%nnodes,3)  

      REAL(8) temp(mesh%nnodes)      
      REAL(8) qtemp(mesh%nq)      
      REAL(8) ptstemp(mesh%nnodes)                  
      REAL(8) pnlstemp(mesh%nnodes)      
      REAL(8) pptstemp(mesh%nnodes)                  

      REAL(8) dc(mesh%nnodes)
      REAL(8) qdc(mesh%nq)
      REAL(8) ptsdc(mesh%nnodes)
      REAL(8) pnlsdc(mesh%nnodes)
      REAL(8) pptsdc(mesh%nnodes)

      INTEGER version
      INTEGER j  

      OPEN (io%rst,FILE=io%rst_name,FORM='UNFORMATTED',STATUS='OLD',ERR=10)
      READ(io%rst) version
           
      IF (version.LE.inp%RSTversion.AND.version.GE.3) THEN
        CALL WarnErr(env,io,inp,0,"RrstF","Reading restart file ...",0)
        READ(io%rst) cnt%tstep,cnt%glnnlit
        READ(io%rst) cnt%rtime
C
C____   FIELD FUNCTIONS :
C     
      DO j=1,3
        CALL rdvec(io%rst,mesh%nnodes,velocity(1,j))
      END DO
      DO j=1,3
        CALL rdvec(io%rst,mesh%nnodes,pnlsvelocity(1,j))        
      END DO
      DO j=1,3      
        CALL rdvec(io%rst,mesh%nnodes,vorticity(1,j))
      END DO
      DO j=1,3      
        CALL rdvec(io%rst,mesh%nnodes,ptsvorticity(1,j))
      END DO
      DO j=1,3      
        CALL rdvec(io%rst,mesh%nnodes,pptsvorticity(1,j))
      END DO
      DO j=1,3      
        CALL rdvec(io%rst,mesh%nnodes,pnlsvorticity(1,j))          
      END DO
      DO j=1,3
        CALL rdvec(io%rst,mesh%nq,qvorticity(1,j))                              
      END DO
      CALL rdvec(io%rst,mesh%nnodes,temp)
      CALL rdvec(io%rst,mesh%nnodes,ptstemp)
      CALL rdvec(io%rst,mesh%nnodes,pptstemp)
      CALL rdvec(io%rst,mesh%nnodes,pnlstemp)          
      CALL rdvec(io%rst,mesh%nq,qtemp)       

      IF (version.GE.4) THEN ! v verziji restarta 4 dodal DC v restart datoteko
        CALL rdvec(io%rst,mesh%nnodes,dc)
        CALL rdvec(io%rst,mesh%nnodes,ptsdc)
        CALL rdvec(io%rst,mesh%nnodes,pptsdc)
        CALL rdvec(io%rst,mesh%nnodes,pnlsdc)
        CALL rdvec(io%rst,mesh%nq,qdc)
      END IF


        CALL WarnErr(env,io,inp,0,"RrstF","Finished reading restart file.",0)
      ELSE ! wrong restart file version
        CALL WarnErr(env,io,inp,3,"RrstF","Wrong restart file version",version)
      END IF

      CLOSE (io%rst)  
      RETURN

10    CONTINUE  ! Could not open the restart file        
      CALL WarnErr(env,io,inp,1,"RrstF","Could not open the restart file",0)
        
      END


C -----------------------------------------------------------------------------
      SUBROUTINE WrstF(mesh,io,inp,cnt,velocity,vorticity,
     &             ptsvorticity,pptsvorticity,qvorticity,pnlsvorticity,pnlsvelocity,
     &             temp,ptstemp,pptstemp,qtemp,pnlstemp,
     &             dc,ptsdc,pptsdc,qdc,pnlsdc)
C
C     $: writes restart file
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(IOtype) io
      TYPE(countType) cnt
      TYPE(InputType) :: inp
      TYPE(meshType) :: mesh
      
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) pnlsvelocity(mesh%nnodes,3)      
      REAL(8) vorticity(mesh%nnodes,3)      
      REAL(8) qvorticity(mesh%nq,3)      
      REAL(8) ptsvorticity(mesh%nnodes,3)                  
      REAL(8) pnlsvorticity(mesh%nnodes,3)      
      REAL(8) pptsvorticity(mesh%nnodes,3)            

      REAL(8) temp(mesh%nnodes)      
      REAL(8) qtemp(mesh%nq)      
      REAL(8) ptstemp(mesh%nnodes)                  
      REAL(8) pnlstemp(mesh%nnodes)      
      REAL(8) pptstemp(mesh%nnodes)       

      REAL(8) dc(mesh%nnodes)
      REAL(8) qdc(mesh%nq)
      REAL(8) ptsdc(mesh%nnodes)
      REAL(8) pnlsdc(mesh%nnodes)
      REAL(8) pptsdc(mesh%nnodes)
            
      INTEGER j  

      OPEN (io%rst,FILE=io%rst_name,FORM='UNFORMATTED',STATUS='UNKNOWN')
      WRITE(io%rst) inp%RSTversion
      WRITE(io%rst) cnt%tstep,cnt%glnnlit
      WRITE(io%rst) cnt%rtime
C
C____ FIELD FUNCTIONS :
C     
      DO j=1,3
        CALL wrvec(io%rst,mesh%nnodes,velocity(1,j))
      END DO
      DO j=1,3
        CALL wrvec(io%rst,mesh%nnodes,pnlsvelocity(1,j))        
      END DO
      DO j=1,3      
        CALL wrvec(io%rst,mesh%nnodes,vorticity(1,j))
      END DO
      DO j=1,3      
        CALL wrvec(io%rst,mesh%nnodes,ptsvorticity(1,j))
      END DO
      DO j=1,3      
        CALL wrvec(io%rst,mesh%nnodes,pptsvorticity(1,j))
      END DO
      DO j=1,3      
        CALL wrvec(io%rst,mesh%nnodes,pnlsvorticity(1,j))          
      END DO
      DO j=1,3
        CALL wrvec(io%rst,mesh%nq,qvorticity(1,j))                              
      END DO
      CALL wrvec(io%rst,mesh%nnodes,temp)
      CALL wrvec(io%rst,mesh%nnodes,ptstemp)
      CALL wrvec(io%rst,mesh%nnodes,pptstemp)
      CALL wrvec(io%rst,mesh%nnodes,pnlstemp)          
      CALL wrvec(io%rst,mesh%nq,qtemp)                              

      CALL wrvec(io%rst,mesh%nnodes,dc)
      CALL wrvec(io%rst,mesh%nnodes,ptsdc)
      CALL wrvec(io%rst,mesh%nnodes,pptsdc)
      CALL wrvec(io%rst,mesh%nnodes,pnlsdc)
      CALL wrvec(io%rst,mesh%nq,qdc)

      CLOSE (io%rst)  
        
      END


C -----------------------------------------------------------------------------

      SUBROUTINE AddLeadingZeros(st,niz)

C
C     $: doda leading zeros
C
C -----------------------------------------------------------------------------
      INTEGER st
      CHARACTER(6) niz

      IF (st.LT.10) THEN
        WRITE (niz,'(A,I1)') "00000",st
      ELSE IF (st.LT.100) THEN
        WRITE (niz,'(A,I2)') "0000",st
      ELSE IF (st.LT.1000) THEN
        WRITE (niz,'(A,I3)') "000",st
      ELSE IF (st.LT.10000) THEN
        WRITE (niz,'(A,I4)') "00",st
      ELSE IF (st.LT.100000) THEN
        WRITE (niz,'(A,I5)') "0",st
      ELSE IF (st.LT.1000000) THEN
        WRITE (niz,'(I6)') st
      END IF

      END


C -----------------------------------------------------------------------------

      SUBROUTINE OutputTDResultsParaview(io,inp,cnt,mesh,v,w,T,dc,mH)

C
C     $: Outputs function
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(IOtype) io
      TYPE(meshType) mesh
      TYPE(inputtype) inp
      TYPE(countType) cnt

      REAL(8) v(mesh%nnodes,3)
      REAL(8) w(mesh%nnodes,3)
      REAL(8) T(mesh%nnodes)
      REAL(8) dc(mesh%nnodes)
      REAL(8) mH(mesh%nnodes)

      INTEGER itp,nic,i
      CHARACTER*6 cifra
      CHARACTER*255 vrstica


      itp=96

      CALL AddLeadingZeros(cnt%tstep,cifra)
      WRITE (vrstica,'(A,A6,A)') "tri_",trim(cifra),".vtu"
      OPEN (itp,FILE=trim(vrstica),STATUS='UNKNOWN')

      nic=8*mesh%nicell



      WRITE (itp,'(A)') '<?xml version="1.0"?>'
      WRITE (itp,'(A)')
     &'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
      WRITE (itp,'(A)') '<UnstructuredGrid>'
      WRITE (itp,'(A,I10,A,I10,A)') '<Piece NumberOfPoints="',mesh%NNODES,'" NumberOfCells="',nic,'">'

      WRITE (itp,'(A)') '<PointData Vectors="velocity,vorticity" Scalars="temperature,dc">'

C      HITROST
      WRITE (itp,'(A)') '<DataArray type="Float32" Name="velocity" NumberOfComponents="3" format="ascii">'
      DO  i=1,mesh%nnodes
        WRITE (vrstica,'(3G18.9)') v(i,1),v(i,2),v(i,3)
        CALL sqblnk(itp,vrstica)
      END DO
      WRITE (itp,'(A)') '</DataArray>'

C      VRTINCNOST
      WRITE (itp,'(A)') '<DataArray type="Float32" Name="vorticity" NumberOfComponents="3" format="ascii">'
      DO  i=1,mesh%nnodes
        WRITE (vrstica,'(3G18.9)') w(i,1),w(i,2),w(i,3)
        CALL sqblnk(itp,vrstica)
      END DO
      WRITE (itp,'(A)') '</DataArray>'


C      TEMPERATURA
      WRITE (itp,'(A)') '<DataArray type="Float32" Name="temperature" format="ascii">'
      DO  i=1,mesh%nnodes
        WRITE (itp,*) T(i)
      END DO
      WRITE (itp,'(A)') '</DataArray>'


C      DC
      WRITE (itp,'(A)') '<DataArray type="Float32" Name="dc" format="ascii">'
      DO  i=1,mesh%nnodes
        WRITE (itp,*) dc(i)
      END DO
      WRITE (itp,'(A)') '</DataArray>'

C      mH
      WRITE (itp,'(A)') '<DataArray type="Float32" Name="modHelm" format="ascii">'
      DO  i=1,mesh%nnodes
        WRITE (itp,*) mh(i)
      END DO
      WRITE (itp,'(A)') '</DataArray>'

      WRITE (itp,'(A)') '</PointData>'

C      KOORDINATE VOZLISC
      WRITE (itp,'(A)') '<Points>'
      WRITE (itp,'(A)') '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'

      DO  i=1,mesh%nnodes
        WRITE (vrstica,'(3G18.9)') mesh%x(i,1),mesh%x(i,2),mesh%x(i,3)
        CALL sqblnk(itp,vrstica)
      END DO

      WRITE (itp,'(A)') '</DataArray>'
      WRITE (itp,'(A)') '</Points>'

      WRITE (itp,'(A)') '<Cells>'
        WRITE (itp,'(A)') '<DataArray type="Int32" Name="connectivity" format="ascii">'
        CALL TECwIDCpara(itp,mesh)
        WRITE (itp,'(A)') '</DataArray>'
        WRITE (itp,'(A)') '<DataArray type="Int32" Name="offsets" format="ascii">'
          DO i=1,nic
            WRITE (itp,*) i*8
          END DO
        WRITE (itp,'(A)') '</DataArray>'
        WRITE (itp,'(A)') '<DataArray type="UInt8" Name="types" format="ascii">'
          DO i=1,nic
            WRITE (itp,*) "12" ! ker je 12 heksaeder
          END DO
        WRITE (itp,'(A)') '</DataArray>'
      WRITE (itp,'(A)') '</Cells>'


      WRITE (itp,'(A)') '</Piece>'
      WRITE (itp,'(A)') '</UnstructuredGrid>'
      WRITE (itp,'(A)') '</VTKFile>'

      CLOSE (itp)


      END


C -----------------------------------------------------------------------------

      SUBROUTINE OutputTDResultsParaviewVTK(io,inp,cnt,mesh,v,w,T,dc,mH)

C
C     $: Outputs results in VTK format
C
C -----------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(IOtype) io
      TYPE(meshType) mesh
      TYPE(inputtype) inp
      TYPE(countType) cnt

      REAL(8) v(mesh%nnodes,3)
      REAL(8) w(mesh%nnodes,3)
      REAL(8) T(mesh%nnodes)
      REAL(8) dc(mesh%nnodes)
      REAL(8) mH(mesh%nnodes)

      INTEGER itp,nic,i
      CHARACTER*6 cifra
      CHARACTER*255 vrstica


      itp=96

      CALL AddLeadingZeros(cnt%tstep,cifra)
      WRITE (vrstica,'(A,A6,A)') "tri_",trim(cifra),".vtk"
      OPEN (itp,FILE=trim(vrstica),STATUS='UNKNOWN')

      nic=8*mesh%nicell


      WRITE (itp,'(A)') '# vtk DataFile Version 2.0'
      WRITE (itp,'(A)') 'Tritok VTK results'
      WRITE (itp,'(A)') 'ASCII'
      WRITE (itp,'(A)') 'DATASET UNSTRUCTURED_GRID'
C
C     Export points
C
      WRITE (itp,'(A,I10,A)') 'POINTS',mesh%nnodes," float"
      DO  i=1,mesh%nnodes
        WRITE (vrstica,'(3G18.9)') mesh%x(i,1),mesh%x(i,2),mesh%x(i,3)
        CALL sqblnk(itp,vrstica)
      END DO
C
C     Export cell conectivity
C
      WRITE (itp,'(A,I10,I10)') 'CELLS ',nic,nic*9 ! (8 za heksa + 1 za st)
      CALL VTKwIDC(itp,mesh)
C
C     Export cell types
C
      WRITE (itp,'(A,I10)') 'CELL_TYPES ',nic
      DO i=1,nic
        WRITE (itp,'(A)') "12"  ! ker je 12 heksaeder
      END DO
C
C     Export field data
C
      WRITE (itp,'(A,I10)') 'POINT_DATA ',mesh%nnodes
      WRITE (itp,'(A,I10)') 'FIELD attributes 5'

      WRITE (itp,'(A,I10,A)') 'U 3 ',mesh%nnodes," float"
      DO i=1,mesh%nnodes
        WRITE (vrstica,*) v(i,1),v(i,2),v(i,3)
        CALL sqblnk(itp,vrstica)
      END DO

      WRITE (itp,'(A,I10,A)') 'Vort 3 ',mesh%nnodes," float"
      DO i=1,mesh%nnodes
        WRITE (vrstica,*) w(i,1),w(i,2),w(i,3)
        CALL sqblnk(itp,vrstica)
      END DO

      WRITE (itp,'(A,I10,A)') 'T 1 ',mesh%nnodes," float"
      DO i=1,mesh%nnodes
        WRITE (vrstica,*) T(i)
        CALL sqblnk(itp,vrstica)
      END DO


      WRITE (itp,'(A,I10,A)') 'DC 1 ',mesh%nnodes," float"
      DO i=1,mesh%nnodes
        WRITE (vrstica,*) dc(i)+0.3
        CALL sqblnk(itp,vrstica)
      END DO

      WRITE (itp,'(A,I10,A)') 'mH 1 ',mesh%nnodes," float"
      DO i=1,mesh%nnodes
        WRITE (vrstica,*) mH(i)
        CALL sqblnk(itp,vrstica)
      END DO


      CLOSE (itp)


      END



C -----------------------------------------------------------------------------
      
      SUBROUTINE OutputTDResultsTecplotvw(io,inp,cnt,mesh,v,w,T,dc,mHu)
C
C     $: Outputs function
C
C -----------------------------------------------------------------------------        
      USE inc_types
      TYPE(IOtype) io
      TYPE(meshType) mesh
      TYPE(inputtype) inp       
      TYPE(countType) cnt  
            
      REAL(8) v(mesh%nnodes,3)
      REAL(8) w(mesh%nnodes,3)      
      REAL(8) T(mesh%nnodes)
      REAL(8) dc(mesh%nnodes)        
      REAL(8) mHu(mesh%nnodes)
 
      INTEGER itp,nic,i
      CHARACTER*100 zonetitle,pfn
      CHARACTER*255 vrstica
      
      itp=io%res

      nic=8*mesh%nicell
      Write(zonetitle,'(I5,A1,G12.7)') cnt%tstep,'-',cnt%rtime
c      zonetitle="u(x,y,z)"
C     First time step - requires header
      IF (cnt%ftime.EQ.1) THEN
        cnt%ftime=0
        OPEN (itp,FILE=TRIM(io%res_name),STATUS='UNKNOWN')

        WRITE (itp,'(A)') '# |---------------------------|' 
        WRITE (itp,'(A)') '# |                           |' 
        WRITE (itp,'(A)') '# |        T R I T O K        |' 
        WRITE (itp,'(A)') '# |                           |'
        WRITE (itp,'(A)') '# |---------------------------|'  
        
        WRITE (itp,'(A)') 'VARIABLES = "X", "Y", "Z", "node","vx","vy","vz","wx","wy","wz","T","dc","mH"'
        CALL GETARG(0,pfn)
        WRITE (itp,'(7A)') 'DATASETAUXDATA code = "',trim(inp%IDname),' v',trim(inp%IDversion),', ',trim(inp%IDdate),'"'
        WRITE (itp,'(3A)') 'DATASETAUXDATA program = "',trim(pfn),'"'
        CALL CopyINP2AUX(io)

        WRITE (itp,'(A,A,A,I8,A,I8)') 
     &      'ZONE T="',trim(zonetitle),'",F=FEPOINT, ET=BRICK,N= ', mesh%NNODES, ',E=', nic
        CALL ZoneAuxData(itp,cnt%tstep,cnt%rtime)
C       Output results        
        DO  i=1,mesh%nnodes
          WRITE (vrstica,'(3G18.9,I8,9G18.9)') 
     &        mesh%x(i,1),mesh%x(i,2),mesh%x(i,3),i,v(i,1),v(i,2),v(i,3),w(i,1),w(i,2),w(i,3),T(i),dc(i),mHu(i)
          CALL sqblnk(itp,vrstica)           
        END DO
        CALL TECwIDC(itp,mesh)
        
      ELSE ! not the first entry - append
        WRITE (itp,'(A8,A,A,I8,A,I8,A)') 
     &      'ZONE T="',trim(zonetitle),'",F=FEPOINT, ET=BRICK,N= ', mesh%NNODES, ',E=', nic,
     &      ', VARSHARELIST=([1,2,3,4]=1), CONNECTIVITYSHAREZONE = 1'
        CALL ZoneAuxData(itp,cnt%tstep,cnt%rtime)      

C       Output results        
        DO  i=1,mesh%nnodes
          WRITE (vrstica,'(9G18.9)') v(i,1),v(i,2),v(i,3),w(i,1),w(i,2),w(i,3),T(i),dc(i),mHu(i)
          CALL sqblnk(itp,vrstica)           
        END DO        
      END IF
      
      END   


C -----------------------------------------------------------------------------

      SUBROUTINE OutputTDMaterialProperties(io,inp,cnt,mesh,rho,cp,lambda,pvf)
C
C     $: Outputs function
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(IOtype) io
      TYPE(meshType) mesh
      TYPE(inputtype) inp
      TYPE(countType) cnt

      REAL(8) lambda(mesh%nnodes)
      REAL(8) rho(mesh%nnodes)
      REAL(8) cp(mesh%nnodes)
      REAL(8) pvf(mesh%nnodes)

      INTEGER itp,nic,i
      CHARACTER*100 zonetitle,pfn
      CHARACTER*255 vrstica

      itp=io%Mres

      nic=8*mesh%nicell
      Write(zonetitle,'(I5,A1,G12.7)') cnt%tstep,'-',cnt%rtime
c      zonetitle="u(x,y,z)"
C     First time step - requires header
      IF (cnt%ftime.EQ.1) THEN
        OPEN (itp,FILE=TRIM(io%Mres_name),STATUS='UNKNOWN')

        WRITE (itp,'(A)') '# |---------------------------|'
        WRITE (itp,'(A)') '# |                           |'
        WRITE (itp,'(A)') '# |        T R I T O K        |'
        WRITE (itp,'(A)') '# |         MATERIALS         |'
        WRITE (itp,'(A)') '# |---------------------------|'

        WRITE (itp,'(A)') 'VARIABLES = "X", "Y", "Z", "node","rho","cp","lambda","pvf"'
        CALL GETARG(0,pfn)
        WRITE (itp,'(7A)') 'DATASETAUXDATA code = "',trim(inp%IDname),' v',trim(inp%IDversion),', ',trim(inp%IDdate),'"'
        WRITE (itp,'(3A)') 'DATASETAUXDATA program = "',trim(pfn),'"'

        WRITE (itp,'(A,A,A,I8,A,I8)')
     &      'ZONE T="',trim(zonetitle),'",F=FEPOINT, ET=BRICK,N= ', mesh%NNODES, ',E=', nic
        CALL ZoneAuxData(itp,cnt%tstep,cnt%rtime)
C       Output results
        DO  i=1,mesh%nnodes
          WRITE (vrstica,'(3G18.9,I8,4G18.9)')
     &        mesh%x(i,1),mesh%x(i,2),mesh%x(i,3),i,rho(i),cp(i),lambda(i),pvf(i)
          CALL sqblnk(itp,vrstica)
        END DO
        CALL TECwIDC(itp,mesh)

      ELSE ! not the first entry - append
        WRITE (itp,'(A8,A,A,I8,A,I8,A)')
     &      'ZONE T="',trim(zonetitle),'",F=FEPOINT, ET=BRICK,N= ', mesh%NNODES, ',E=', nic,
     &      ', VARSHARELIST=([1,2,3,4]=1), CONNECTIVITYSHAREZONE = 1'
        CALL ZoneAuxData(itp,cnt%tstep,cnt%rtime)

C       Output results
        DO  i=1,mesh%nnodes
          WRITE (vrstica,'(4G18.9)') rho(i),cp(i),lambda(i),pvf(i)
          CALL sqblnk(itp,vrstica)
        END DO
      END IF

      END



C -----------------------------------------------------------------------------

      SUBROUTINE OutputInitial(io,inp,mesh,v,w,T,dc)
C
C     $: Outputs functions as defined in tri.bic
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(IOtype) io
      TYPE(meshType) mesh
      TYPE(inputtype) inp

      REAL(8) v(mesh%nnodes,3)
      REAL(8) w(mesh%nnodes,3)
      REAL(8) T(mesh%nnodes)
      REAL(8) dc(mesh%nnodes)
      REAL(8), ALLOCATABLE :: nx(:),ny(:),nz(:)
      REAL(8), ALLOCATABLE :: bcnx(:),bcny(:),bcnz(:)

      INTEGER itp,nic,i,dn
      CHARACTER*100 zonetitle,pfn
      CHARACTER*355 vrstica

      itp=io%resBIC

      ALLOCATE (nx(mesh%nnodes))
      ALLOCATE (ny(mesh%nnodes))
      ALLOCATE (nz(mesh%nnodes))
      ALLOCATE (bcnx(mesh%nnodes))
      ALLOCATE (bcny(mesh%nnodes))
      ALLOCATE (bcnz(mesh%nnodes))
      nx=0.0D0
      ny=0.0D0
      nz=0.0D0
      bcnx=0.0D0
      bcny=0.0D0
      bcnz=0.0D0


      DO i=1,mesh%nbnodes
        dn=mesh%gbn(i)
        nx(dn)=mesh%nx(i)
        ny(dn)=mesh%ny(i)
        nz(dn)=mesh%nz(i)
        bcnx(dn)=mesh%bcn(i,1)
        bcny(dn)=mesh%bcn(i,2)
        bcnz(dn)=mesh%bcn(i,3)
      END DO


      nic=8*mesh%nicell
      Write(zonetitle,'(A)') "Initial"
        OPEN (itp,FILE=TRIM(io%resBIC_name),STATUS='UNKNOWN')

        WRITE (itp,'(A)') '# |---------------------------|'
        WRITE (itp,'(A)') '# |                           |'
        WRITE (itp,'(A)') '# |        T R I T O K        |'
        WRITE (itp,'(A)') '# |   INITIAL TRI.BIC DATA    |'
        WRITE (itp,'(A)') '# |---------------------------|'

        WRITE (itp,'(A,A)') 'VARIABLES = "X", "Y", "Z", "node","vx","vy","vz","wx","wy","wz","T"',
     &                      ',"dc","nx","ny","nz","bcnx","bcny","bcnz"'
        CALL GETARG(0,pfn)
        WRITE (itp,'(7A)') 'DATASETAUXDATA code = "',trim(inp%IDname),' v',trim(inp%IDversion),', ',trim(inp%IDdate),'"'
        WRITE (itp,'(3A)') 'DATASETAUXDATA program = "',trim(pfn),'"'

        WRITE (itp,'(A,A,A,I8,A,I8)')
     &      'ZONE T="',trim(zonetitle),'",F=FEPOINT, ET=BRICK,N= ', mesh%NNODES, ',E=', nic
C       Output results
        DO  i=1,mesh%nnodes
          WRITE (vrstica,'(3G18.10,I8,16G18.10)')
     &        mesh%x(i,1),mesh%x(i,2),mesh%x(i,3),i,v(i,1),v(i,2),v(i,3),w(i,1),w(i,2),w(i,3),T(i),dc(i),
     &        nx(i),ny(i),nz(i),bcnx(i),bcny(i),bcnz(i)
          CALL sqblnk(itp,vrstica)
        END DO
        CALL TECwIDC(itp,mesh)

      CLOSE(itp)

      DEALLOCATE (nx,ny,nz)
      DEALLOCATE (bcnx,bcny,bcnz)

      END


C______________________________________________________________________C
      SUBROUTINE TECwIDC(lun,mesh)
C
C     Tecplot write cell conectivity                                   C
C______________________________________________________________________C
      USE inc_types
      TYPE(meshType) mesh
C....
      CHARACTER*200 vrstica
      INTEGER i,lun
      
      INTEGER a1(8),a2(8),a3(8),a4(8),a5(8),a6(8),a7(8),a8(8)
      
      DATA a1/6,18,23,14,17,26,27,22/
      DATA a2/18,7,15,23,26,19,24,27/
      DATA a3/14,23,10,2,22,27,21,9/
      DATA a4/23,15,3,10,27,24,11,21/
      DATA a5/17,26,27,22,5,20,25,13/
      DATA a6/26,19,24,27,20,8,16,25/
      DATA a7/22,27,21,9,13,25,12,1/
      DATA a8/27,24,11,21,25,16,4,12/
C
C     Write cell connectivity
C
      IF (mesh%npoc.EQ.9) THEN
        DO i=1,mesh%nicell
          WRITE(vrstica,*) mesh%idc(i,1),mesh%idc(i,2),mesh%idc(i,9),mesh%idc(i,8)
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,2),mesh%idc(i,3),mesh%idc(i,4),mesh%idc(i,9)
          CALL sqblnk(lun,vrstica)          
          WRITE(vrstica,*) mesh%idc(i,8),mesh%idc(i,9),mesh%idc(i,6),mesh%idc(i,7)
          CALL sqblnk(lun,vrstica)          
          WRITE(vrstica,*) mesh%idc(i,9),mesh%idc(i,4),mesh%idc(i,5),mesh%idc(i,6)
          CALL sqblnk(lun,vrstica)          
        END DO
      ELSE IF (mesh%npoc.EQ.4) THEN
        DO i=1,mesh%nicell
          WRITE(vrstica,*) mesh%idc(i,1),mesh%idc(i,2),mesh%idc(i,3),mesh%idc(i,4)
          CALL sqblnk(lun,vrstica)          
        END DO
      ELSE IF (mesh%npoc.EQ.27) THEN
        DO i=1,mesh%nicell
          WRITE(vrstica,*) mesh%idc(i,a1(1)),mesh%idc(i,a1(2)),mesh%idc(i,a1(3)),mesh%idc(i,a1(4))
     &                    ,mesh%idc(i,a1(5)),mesh%idc(i,a1(6)),mesh%idc(i,a1(7)),mesh%idc(i,a1(8))
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,a2(1)),mesh%idc(i,a2(2)),mesh%idc(i,a2(3)),mesh%idc(i,a2(4))
     &                    ,mesh%idc(i,a2(5)),mesh%idc(i,a2(6)),mesh%idc(i,a2(7)),mesh%idc(i,a2(8))
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,a3(1)),mesh%idc(i,a3(2)),mesh%idc(i,a3(3)),mesh%idc(i,a3(4))
     &                    ,mesh%idc(i,a3(5)),mesh%idc(i,a3(6)),mesh%idc(i,a3(7)),mesh%idc(i,a3(8))
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,a4(1)),mesh%idc(i,a4(2)),mesh%idc(i,a4(3)),mesh%idc(i,a4(4))
     &                    ,mesh%idc(i,a4(5)),mesh%idc(i,a4(6)),mesh%idc(i,a4(7)),mesh%idc(i,a4(8))
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,a5(1)),mesh%idc(i,a5(2)),mesh%idc(i,a5(3)),mesh%idc(i,a5(4))
     &                    ,mesh%idc(i,a5(5)),mesh%idc(i,a5(6)),mesh%idc(i,a5(7)),mesh%idc(i,a5(8))
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,a6(1)),mesh%idc(i,a6(2)),mesh%idc(i,a6(3)),mesh%idc(i,a6(4))
     &                    ,mesh%idc(i,a6(5)),mesh%idc(i,a6(6)),mesh%idc(i,a6(7)),mesh%idc(i,a6(8))
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,a7(1)),mesh%idc(i,a7(2)),mesh%idc(i,a7(3)),mesh%idc(i,a7(4))
     &                    ,mesh%idc(i,a7(5)),mesh%idc(i,a7(6)),mesh%idc(i,a7(7)),mesh%idc(i,a7(8))
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,a8(1)),mesh%idc(i,a8(2)),mesh%idc(i,a8(3)),mesh%idc(i,a8(4))
     &                    ,mesh%idc(i,a8(5)),mesh%idc(i,a8(6)),mesh%idc(i,a8(7)),mesh%idc(i,a8(8))
          CALL sqblnk(lun,vrstica)

        END DO

      END IF

      END       
      

C______________________________________________________________________C
      SUBROUTINE TECwIDCpara(lun,mesh)
C
C     Tecplot write cell conectivity                                   C
C______________________________________________________________________C
      USE inc_types
      TYPE(meshType) mesh
C....
      CHARACTER*200 vrstica
      INTEGER i,lun

      INTEGER a1(8),a2(8),a3(8),a4(8),a5(8),a6(8),a7(8),a8(8)

      DATA a1/6,18,23,14,17,26,27,22/
      DATA a2/18,7,15,23,26,19,24,27/
      DATA a3/14,23,10,2,22,27,21,9/
      DATA a4/23,15,3,10,27,24,11,21/
      DATA a5/17,26,27,22,5,20,25,13/
      DATA a6/26,19,24,27,20,8,16,25/
      DATA a7/22,27,21,9,13,25,12,1/
      DATA a8/27,24,11,21,25,16,4,12/
C
C     Write cell connectivity
C
      IF (mesh%npoc.EQ.9) THEN
        DO i=1,mesh%nicell
          WRITE(vrstica,*) mesh%idc(i,1)-1,mesh%idc(i,2)-1,mesh%idc(i,9)-1,mesh%idc(i,8)-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,2)-1,mesh%idc(i,3)-1,mesh%idc(i,4)-1,mesh%idc(i,9)-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,8)-1,mesh%idc(i,9)-1,mesh%idc(i,6)-1,mesh%idc(i,7)-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,9)-1,mesh%idc(i,4)-1,mesh%idc(i,5)-1,mesh%idc(i,6)-1
          CALL sqblnk(lun,vrstica)
        END DO
      ELSE IF (mesh%npoc.EQ.4) THEN
        DO i=1,mesh%nicell
          WRITE(vrstica,*) mesh%idc(i,1)-1,mesh%idc(i,2)-1,mesh%idc(i,3)-1,mesh%idc(i,4)-1
          CALL sqblnk(lun,vrstica)
        END DO
      ELSE IF (mesh%npoc.EQ.27) THEN
        DO i=1,mesh%nicell
          WRITE(vrstica,*) mesh%idc(i,a1(1))-1,mesh%idc(i,a1(2))-1,mesh%idc(i,a1(3))-1,mesh%idc(i,a1(4))-1
     &                    ,mesh%idc(i,a1(5))-1,mesh%idc(i,a1(6))-1,mesh%idc(i,a1(7))-1,mesh%idc(i,a1(8))-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,a2(1))-1,mesh%idc(i,a2(2))-1,mesh%idc(i,a2(3))-1,mesh%idc(i,a2(4))-1
     &                    ,mesh%idc(i,a2(5))-1,mesh%idc(i,a2(6))-1,mesh%idc(i,a2(7))-1,mesh%idc(i,a2(8))-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,a3(1))-1,mesh%idc(i,a3(2))-1,mesh%idc(i,a3(3))-1,mesh%idc(i,a3(4))-1
     &                    ,mesh%idc(i,a3(5))-1,mesh%idc(i,a3(6))-1,mesh%idc(i,a3(7))-1,mesh%idc(i,a3(8))-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,a4(1))-1,mesh%idc(i,a4(2))-1,mesh%idc(i,a4(3))-1,mesh%idc(i,a4(4))-1
     &                    ,mesh%idc(i,a4(5))-1,mesh%idc(i,a4(6))-1,mesh%idc(i,a4(7))-1,mesh%idc(i,a4(8))-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,a5(1))-1,mesh%idc(i,a5(2))-1,mesh%idc(i,a5(3))-1,mesh%idc(i,a5(4))-1
     &                    ,mesh%idc(i,a5(5))-1,mesh%idc(i,a5(6))-1,mesh%idc(i,a5(7))-1,mesh%idc(i,a5(8))-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,a6(1))-1,mesh%idc(i,a6(2))-1,mesh%idc(i,a6(3))-1,mesh%idc(i,a6(4))-1
     &                    ,mesh%idc(i,a6(5))-1,mesh%idc(i,a6(6))-1,mesh%idc(i,a6(7))-1,mesh%idc(i,a6(8))-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,a7(1))-1,mesh%idc(i,a7(2))-1,mesh%idc(i,a7(3))-1,mesh%idc(i,a7(4))-1
     &                    ,mesh%idc(i,a7(5))-1,mesh%idc(i,a7(6))-1,mesh%idc(i,a7(7))-1,mesh%idc(i,a7(8))-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) mesh%idc(i,a8(1))-1,mesh%idc(i,a8(2))-1,mesh%idc(i,a8(3))-1,mesh%idc(i,a8(4))-1
     &                    ,mesh%idc(i,a8(5))-1,mesh%idc(i,a8(6))-1,mesh%idc(i,a8(7))-1,mesh%idc(i,a8(8))-1
          CALL sqblnk(lun,vrstica)

        END DO

      END IF

      END


C______________________________________________________________________C
      SUBROUTINE VTKwIDC(lun,mesh)
C
C     Tecplot write cell conectivity                                   C
C______________________________________________________________________C
      USE inc_types
      TYPE(meshType) mesh
C....
      CHARACTER*200 vrstica
      INTEGER i,lun

      INTEGER a1(8),a2(8),a3(8),a4(8),a5(8),a6(8),a7(8),a8(8)

      DATA a1/6,18,23,14,17,26,27,22/
      DATA a2/18,7,15,23,26,19,24,27/
      DATA a3/14,23,10,2,22,27,21,9/
      DATA a4/23,15,3,10,27,24,11,21/
      DATA a5/17,26,27,22,5,20,25,13/
      DATA a6/26,19,24,27,20,8,16,25/
      DATA a7/22,27,21,9,13,25,12,1/
      DATA a8/27,24,11,21,25,16,4,12/
C
C     Write cell connectivity
C
      IF (mesh%npoc.EQ.27) THEN
        DO i=1,mesh%nicell
          WRITE(vrstica,*) "8 ",mesh%idc(i,a1(1))-1,mesh%idc(i,a1(2))-1,mesh%idc(i,a1(3))-1,mesh%idc(i,a1(4))-1
     &                    ,mesh%idc(i,a1(5))-1,mesh%idc(i,a1(6))-1,mesh%idc(i,a1(7))-1,mesh%idc(i,a1(8))-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) "8 ",mesh%idc(i,a2(1))-1,mesh%idc(i,a2(2))-1,mesh%idc(i,a2(3))-1,mesh%idc(i,a2(4))-1
     &                    ,mesh%idc(i,a2(5))-1,mesh%idc(i,a2(6))-1,mesh%idc(i,a2(7))-1,mesh%idc(i,a2(8))-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) "8 ",mesh%idc(i,a3(1))-1,mesh%idc(i,a3(2))-1,mesh%idc(i,a3(3))-1,mesh%idc(i,a3(4))-1
     &                    ,mesh%idc(i,a3(5))-1,mesh%idc(i,a3(6))-1,mesh%idc(i,a3(7))-1,mesh%idc(i,a3(8))-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) "8 ",mesh%idc(i,a4(1))-1,mesh%idc(i,a4(2))-1,mesh%idc(i,a4(3))-1,mesh%idc(i,a4(4))-1
     &                    ,mesh%idc(i,a4(5))-1,mesh%idc(i,a4(6))-1,mesh%idc(i,a4(7))-1,mesh%idc(i,a4(8))-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) "8 ",mesh%idc(i,a5(1))-1,mesh%idc(i,a5(2))-1,mesh%idc(i,a5(3))-1,mesh%idc(i,a5(4))-1
     &                    ,mesh%idc(i,a5(5))-1,mesh%idc(i,a5(6))-1,mesh%idc(i,a5(7))-1,mesh%idc(i,a5(8))-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) "8 ",mesh%idc(i,a6(1))-1,mesh%idc(i,a6(2))-1,mesh%idc(i,a6(3))-1,mesh%idc(i,a6(4))-1
     &                    ,mesh%idc(i,a6(5))-1,mesh%idc(i,a6(6))-1,mesh%idc(i,a6(7))-1,mesh%idc(i,a6(8))-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) "8 ",mesh%idc(i,a7(1))-1,mesh%idc(i,a7(2))-1,mesh%idc(i,a7(3))-1,mesh%idc(i,a7(4))-1
     &                    ,mesh%idc(i,a7(5))-1,mesh%idc(i,a7(6))-1,mesh%idc(i,a7(7))-1,mesh%idc(i,a7(8))-1
          CALL sqblnk(lun,vrstica)
          WRITE(vrstica,*) "8 ",mesh%idc(i,a8(1))-1,mesh%idc(i,a8(2))-1,mesh%idc(i,a8(3))-1,mesh%idc(i,a8(4))-1
     &                    ,mesh%idc(i,a8(5))-1,mesh%idc(i,a8(6))-1,mesh%idc(i,a8(7))-1,mesh%idc(i,a8(8))-1
          CALL sqblnk(lun,vrstica)

        END DO
      ELSE
        Print *,"error in VTKwIDC"
        STOP
      END IF

      END

C______________________________________________________________________C  

      SUBROUTINE WarnErr(env,io,inp,WorE,routine,text,ierr)
C                       
C     $: Writes warning or error to error file
C______________________________________________________________________C
      USE inc_types
      TYPE(IOtype) io
      TYPE(inputtype) inp
      TYPE(penv) :: env      
      
      CHARACTER*(*) routine,text
      CHARACTER*255 msgtype
      INTEGER WorE,ierr

C
C       Message, warning or error
C
      IF (env%myproc.EQ.1) THEN

      IF (WorE.EQ.0) THEN
      WRITE(msgtype,'(A7,A4,A,A4,A)') 
     &      "Message"," :: ",trim(routine)," :: ",trim(text)
      ELSE IF (WorE.EQ.5) THEN
       WRITE(msgtype,'(A7,A4,A,A4,A,I8)') 
     &      "Message"," :: ",trim(routine)," :: ",trim(text),ierr
      ELSE IF (WorE.EQ.1) THEN
       WRITE(msgtype,'(A7,A4,A,A4,A)') 
     &      "Warning"," :: ",trim(routine)," :: ",trim(text)
      ELSE IF (WorE.EQ.2) THEN
       WRITE(msgtype,'(A7,A4,A,A4,A)') 
     &      "Error  "," :: ",trim(routine)," :: ",trim(text)
      ELSE IF (WorE.EQ.3) THEN
       WRITE(msgtype,'(A7,A4,A,A4,A,I8)') 
     &      "Warning"," :: ",trim(routine)," :: ",trim(text),ierr
      ELSE IF (WorE.EQ.4) THEN
          WRITE(msgtype,'(A7,A4,A,A4,A,I8)') 
     &      "Error  "," :: ",trim(routine)," :: ",trim(text),ierr
       ELSE
         msgtype="Unkonwn"
        END IF
C
C       write to error file
C   
        WRITE(io%errwarn,'(A)') trim(msgtype)
        CALL FLUSH(io%errwarn)
C
C       stop if error ocured
C     
      END IF

      IF (WorE.EQ.2.OR.WorE.EQ.4) THEN
        CALL StopProgram(env,io,inp,1)
      END IF

      END

C -----------------------------------------------------------------------------
      
      SUBROUTINE SetUpIO(env,io,inp)
C
C     $: sets up input and output file names and numbers
C
C -----------------------------------------------------------------------------       
      USE inc_types
      TYPE(IOtype) io
      TYPE(inputtype) inp
      TYPE(penv) env
      
      CHARACTER(50) hostn      

C
C     Files only processor number 1 opens :
C
      IF (env%myproc.EQ.1) THEN
c     log file      
      io%l=11
      io%l_name='tri.log'
      OPEN (io%l,FILE=io%l_name,STATUS='UNKNOWN')
      WRITE (io%l,'(/6A/)') trim(inp%IDname)," ",trim(inp%IDversion),", ",trim(inp%IDdate),"."      
      WRITE (io%l,'(A)') inp%StartTime
      CALL hostnm(hostn)
      WRITE (io%l,'(A,A)') "Running on: ",hostn
      WRITE (io%l,'(A,I5,A/)') "Using ",env%nproc," processors."

c     error/warning file      
      io%errwarn=13
      io%errwarn_name='tri.err'
      OPEN (io%errwarn,FILE=io%errwarn_name,STATUS='UNKNOWN')

c     Status file - convergence
      io%sta=20
      io%sta_name='tri.sta'
      OPEN (io%sta,FILE=io%sta_name,STATUS='UNKNOWN')
      WRITE(io%sta,'(5A)') trim(inp%IDname),' v',trim(inp%IDversion),', ',trim(inp%IDdate)
      WRITE(io%sta,'(A,8(10X,A2))') " tstep  tsite      vx","vy","vz","wx","wy","wz"," T","DC","mH"

c     Time Status file - convergence
      io%tim=25
      io%tim_name='tri.tim'
      OPEN (io%tim,FILE=io%tim_name,STATUS='UNKNOWN')
      WRITE(io%tim,'(5A)') trim(inp%IDname),' v',trim(inp%IDversion),', ',trim(inp%IDdate)
      WRITE(io%tim,'(A,3(10X,A2))') " tstep  tsite","wx","wy","wz"

      
c     Iteration of solves
      io%ite=19
      io%ite_name='tri.ite'
      OPEN (io%ite,FILE=io%ite_name,STATUS='UNKNOWN')
      WRITE(io%ite,'(5A)') trim(inp%IDname),' v',trim(inp%IDversion),', ',trim(inp%IDdate)
      WRITE(io%ite,'(A,12(2X,A3))') " tstep  tsite","Bwx","Bwy","Bwz","Dvx","Dvy","Dvz","Dwx","Dwy","Dwz","  T"," DC"," mH"

c     epsilon of solves
      io%sol=28
      io%sol_name='tri.sol'
      OPEN (io%sol,FILE=io%sol_name,STATUS='UNKNOWN')
      WRITE(io%sol,'(5A)') trim(inp%IDname),' v',trim(inp%IDversion),', ',trim(inp%IDdate)
      WRITE(io%sol,'(A,12(6X,A3))') " tstep  tsite","Dvx","Dvy","Dvz","Dwx","Dwy","Dwz","  T"," DC"," mH"


c     Temperature flux through wall
      io%tfl=26
      io%tfl_name='tri.tfl'
      OPEN (io%tfl,FILE=io%tfl_name,STATUS='UNKNOWN')
      WRITE(io%tfl,'(5A)') trim(inp%IDname),' v',trim(inp%IDversion),', ',trim(inp%IDdate)
      WRITE(io%tfl,'(A)') " tstep , temperature flux through walls"

c     Temperature flux through wall export
      io%twfl=29
      io%twfl_name='tri.twfl.dat'
      OPEN (io%twfl,FILE=io%twfl_name,STATUS='UNKNOWN')
      WRITE (io%twfl,'(A)') '# |---------------------------|'
      WRITE (io%twfl,'(A)') '# |                           |'
      WRITE (io%twfl,'(A)') '# |    TEMPERATURE  F L U X   |'
      WRITE (io%twfl,'(A)') '# |                           |'
      WRITE (io%twfl,'(A)') '# |---------------------------|'
      WRITE (io%twfl,'(A)') 'VARIABLES = "X", "Y", "Z","qTemp","qwx","qwy","qwz","qmH"'


c     Vorticity flux through wall
      io%wfl=27
      io%wfl_name='tri.wfl'
      OPEN (io%wfl,FILE=io%wfl_name,STATUS='UNKNOWN')
      WRITE(io%wfl,'(5A)') trim(inp%IDname),' v',trim(inp%IDversion),', ',trim(inp%IDdate)
      WRITE(io%wfl,'(A)') " tstep , vorticity flux through walls"

c     Vorticity flux through wall
      io%mass=30
      io%mass_name='tri.mass'
      OPEN (io%mass,FILE=io%mass_name,STATUS='UNKNOWN')
      WRITE(io%mass,'(5A)') trim(inp%IDname),' v',trim(inp%IDversion),', ',trim(inp%IDdate)
      WRITE(io%mass,'(A)') " tstep , mass flux through walls"

c     Temperature flux through wall export
      io%rawl=36
      io%rawl_name='tri.walls.dat'
      OPEN (io%rawl,FILE=io%rawl_name,STATUS='UNKNOWN')
      WRITE (io%rawl,'(A)') '# |---------------------------|'
      WRITE (io%rawl,'(A)') '# |        T R I T O K        |'
      WRITE (io%rawl,'(A)') '# |    results along walls    |'
      WRITE (io%rawl,'(A)') '# |---------------------------|'
      WRITE (io%rawl,'(A)') 'VARIABLES = "X", "Y", "Z", "node","vx","vy","vz","wx","wy","wz","T","dc","qT","qdc"'

c     Local nusselt number
      io%lnu=35
      io%lnu_name='tri.lnu.dat'
      OPEN (io%lnu,FILE=io%lnu_name,STATUS='UNKNOWN')
c      WRITE(io%lnu,'(5A)') trim(inp%IDname),' v',trim(inp%IDversion),', ',trim(inp%IDdate)
      WRITE(io%lnu,'(A)') "VARIABLES = X,Y,Z,FluxInt,Area,Nu"


c     modified Helmholtz RMS file
      io%mhr=37
      io%mhr_name='tri.mH-rms.dat'
      OPEN (io%mHr,FILE=io%mHr_name,STATUS='UNKNOWN')
      WRITE(io%mhr,'(A)') "VARIABLES = time,rms"

c     DC test RMS file
      io%dct=38
      io%dct_name='tri.dc-rms.dat'
      OPEN (io%dct,FILE=io%dct_name,STATUS='UNKNOWN')
      WRITE(io%dct,'(A)') "VARIABLES = time,rms"

c     Stohastic export
      io%std=43
      io%std_name='tri.std.dat'

      END IF


      
c     input file      
      io%inp=12
      io%inp_name='tri.inp'
      
      
c     mesh file
      io%mesh=14      

c     BiC file
      io%bic=16

c     Tecplot results file
      io%res=18
      io%res_name='tri.dat'

c     Particle Tecplot results file
      io%pres=40
      io%pres_name='tri.part.dat'

c     Materials Tecplot results file
      io%mres=41
      io%mres_name='tri.materials.dat'


c     Tecplot results file fo initial data (from tri.bic)
      io%resBIC=32
      io%resBIC_name='tri.initial.dat'


c     flux results file
      io%resf=31
      io%resf_name='tri.q.dat'

c     macro integrals file      
      io%imi=24
      io%imi_name='tri.macro.int'

c     mod. Helmholtz macro integrals file
      io%imh=40
      io%imh_name='tri.modHelm.macro.int'

c     mod. Helmholtz full integrals file
      io%imhf=41
      io%imhf_name='tri.modHelm.full.int'

c     ctree file
      io%ctree=42

c     integrals file      
      io%ikm=22
      IF (env%myproc.LT.10) THEN
        WRITE(io%ikm_name,'(A,I1)') 'tri.sdkm.int.',env%myproc
      ELSE IF (env%myproc.LT.100) THEN
        WRITE(io%ikm_name,'(A,I2)') 'tri.sdkm.int.',env%myproc
      ELSE IF (env%myproc.LT.1000) THEN
        WRITE(io%ikm_name,'(A,I3)') 'tri.sdkm.int.',env%myproc
      END IF

c     integrals file
      io%ikmL=37
      IF (env%myproc.LT.10) THEN
        WRITE(io%ikmL_name,'(A,I1)') 'tri.sdkm-L.int.',env%myproc
      ELSE IF (env%myproc.LT.100) THEN
        WRITE(io%ikmL_name,'(A,I2)') 'tri.sdkm-L.int.',env%myproc
      ELSE IF (env%myproc.LT.1000) THEN
        WRITE(io%ikmL_name,'(A,I3)') 'tri.sdkm-L.int.',env%myproc
      END IF

c     restart file 
      io%rst=21
      io%rst_name='tri.rst'

c     rms difference file
      io%rms=38
      io%rms_name='tri.rms.dat'

c     neil file
      io%neil=39
      io%neil_name='tri.neil'



      RETURN
      

c     integration error file      
      io%interr=15
      io%interr_name='tri.interr.dat'
      OPEN (io%interr,FILE=io%interr_name,STATUS='UNKNOWN')


c     initial flow conditions file
      io%init=17
      io%init_name='tri.initial.dat'

      

c     mass conservation file
      io%mas=23
      io%mas_name='tri.mass'
      OPEN (io%mas,FILE=io%mas_name,STATUS='UNKNOWN')
      WRITE(io%mas,'(A)') " tstep            diff            tot_in          tot_out walls in walls out"  
     
      END       
      
C *********************************************************************
C **                                                                 **
C ** Squeeze multiple blanks into single blank                       **
C **                                                                 **
C *********************************************************************
C---------------------------------------------------------------------C
      SUBROUTINE sqblnk(lun,line)
C---------------------------------------------------------------------C
      CHARACTER*(*) line
      LOGICAL flag
      INTEGER i,ii,j,lnblnk,lun

      j=1
      flag=.FALSE.
      DO i=1,LNBLNK(line)+1
        IF (line(i:i).NE.' ' .AND. .NOT.flag) THEN
          flag=.TRUE.
          ii=i
        ELSE IF (line(i:i).EQ.' ' .AND. flag) THEN
          flag=.FALSE.
          line(j:j+i-ii)=line(ii:i)
          j=j+i-ii+1
        END IF
      END DO
      WRITE(lun,'(A)') line(1:j-2)
      RETURN
      END
      INTEGER FUNCTION LNBLNK (string)
C
C     LNBLNK returns the index of the last non-blank character in string
C
      CHARACTER*(*) string
      INTEGER i

      DO i=LEN(string),1,-1
        IF (string(i:i).NE.' ') THEN
          lnblnk=i
          RETURN
        END IF
      END DO
      lnblnk=0
      RETURN
      END 
            
C______________________________________________________________________C

      SUBROUTINE CopyINP2AUX(io)
C
C     Tecplot Copy bem.inp to Dataset auxillay data in tecplot results C
C______________________________________________________________________C
      USE inc_types
      TYPE(IOtype) io   
C....
      CHARACTER*255 OneLine
      INTEGER i,j

      i=0
      OPEN (io%inp,FILE=io%inp_name,ERR=10,STATUS='OLD') !,SHARED)
      CALL rOneTL(io%inp,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')
c       Tecplot ima probleme z znakom \ -> ga zamenjamo s /
        DO j=1,len_trim(OneLine)
          IF (ichar(OneLine(j:j)).EQ.92) OneLine(j:j)='/'
        END DO
        i=i+1
        IF (i.LT.10) WRITE (io%res,'(A,I1,A,A,A)') 'DATASETAUXDATA INP0',i,' = "',trim(OneLine),'"' 
        IF (i.GE.10) WRITE (io%res,'(A,I2,A,A,A)') 'DATASETAUXDATA INP',i,' = "',trim(OneLine),'"' 
        CALL rOneTL(io%inp,OneLine)
      END DO

      CLOSE(io%inp)

10    RETURN
      END


C______________________________________________________________________C

      SUBROUTINE CopyINP2FILE(io,lun)
C
C     copy tri.inp to file as a remark
C______________________________________________________________________C
      USE inc_types
      IMPLICIT NONE
      TYPE(IOtype) io
      INTEGER lun
C....
      CHARACTER*255 OneLine

      OPEN (io%inp,FILE=io%inp_name,ERR=10,STATUS='OLD') !,SHARED)
      CALL rOneTL(io%inp,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')
        WRITE (lun,'(A,A)') '# ',trim(OneLine)
        CALL rOneTL(io%inp,OneLine)
      END DO

      CLOSE(io%inp)

10    RETURN
      END




      
C______________________________________________________________________C
      SUBROUTINE WriteCpu(io,cpu)
C
C     Write cpu time information to log file
C
      USE inc_types
      TYPE(IOtype) io   
      TYPE(CPUtype) cpu
      INTEGER i
      
      WRITE (io%l,'(//A)') '   CPU time information [s]'
      WRITE (io%l,'(A)') '|-------------------------------|-------------|'
      do i=1,cpu%ncpu
        WRITE (io%l,'(A1,1x,A30,A1,F12.3,1x,A1)') '|',cpu%desc(i),'|',cpu%time(i),'|'
      END DO 
      WRITE (io%l,'(A)') '|-------------------------------|-------------|'

      END      

C -----------------------------------------------------------------------------

      SUBROUTINE OutputStochasticResults(io,inp,cnt,gauss,mesh,v,w,T,qT,qw)
C
C     $: Outputs function
C
C -----------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(IOtype) io
      TYPE(meshType) mesh
      TYPE(gausstype) :: gauss
      TYPE(inputtype) inp
      TYPE(countType) cnt

      REAL(8) v(mesh%nnodes,3)
      REAL(8) w(mesh%nnodes,3)
      REAL(8) T(mesh%nnodes)
      REAL(8) qT(mesh%nq)
      REAL(8) qw(mesh%nq,3)

      REAL(8), ALLOCATABLE :: flux(:)

      INTEGER lun,i,j

      lun=io%std

      ALLOCATE (flux(mesh%nofw))

      OPEN (lun,FILE=io%std_name,STATUS='UNKNOWN')
      WRITE (lun,'(A)') '# |---------------------------------------------|'
      WRITE (lun,'(A)') '# |                 T R I T O K                 |'
      WRITE (lun,'(A)') '# |    results for stochastic postprocessing    |'
      WRITE (lun,'(A)') '# |---------------------------------------------|'
      WRITE (lun,'(6A)')'# ',trim(inp%IDname),' v',trim(inp%IDversion),', ',trim(inp%IDdate)
      CALL CopyINP2FILE(io,lun)
C
C     Temperature flux through walls
C
      WRITE (lun,'(A)') '#'
      WRITE (lun,'(A)') '#  Temperature flux through walls'
      WRITE (lun,'(A)') '#'
      CALL GetWallFlux(mesh,cnt,gauss,qT,flux)
      WRITE (lun,*) mesh%nofw
      DO i=1,mesh%nofw
        WRITE (lun,*) flux(i)
      END DO
C
C     Mass flux through walls
C
      WRITE (lun,'(A)') '#'
      WRITE (lun,'(A)') '#  Mass flux through walls'
      WRITE (lun,'(A)') '#'
      CALL GetMassFlux(mesh,cnt,gauss,v,flux)
      WRITE (lun,*) mesh%nofw
      DO i=1,mesh%nofw
        WRITE (lun,*) flux(i)
      END DO
C
C     Vorticity flux through walls
C
      WRITE (lun,'(A)') '#'
      WRITE (lun,'(A)') '#  Vorticity(x,y,z) flux through walls'
      WRITE (lun,'(A)') '#'
      DO j=1,3
        CALL GetWallFlux(mesh,cnt,gauss,qw(:,j),flux)
        WRITE (lun,*) mesh%nofw
        DO i=1,mesh%nofw
          WRITE (lun,*) flux(i)
        END DO
      END DO
C
C     Velocity nodal values (x,y,z)
C
      WRITE (lun,'(A)') '#'
      WRITE (lun,'(A)') '#  Velocity nodal values (x,y,z)'
      WRITE (lun,'(A)') '#'
      DO j=1,3
        WRITE (lun,*) mesh%nnodes
        DO i=1,mesh%nnodes
          WRITE (lun,*) v(i,j)
        END DO
      END DO
C
C     Vorticity nodal values (x,y,z)
C
      WRITE (lun,'(A)') '#'
      WRITE (lun,'(A)') '#  Vorticity nodal values (x,y,z)'
      WRITE (lun,'(A)') '#'
      DO j=1,3
        WRITE (lun,*) mesh%nnodes
        DO i=1,mesh%nnodes
          WRITE (lun,*) w(i,j)
        END DO
      END DO
C
C     Temperature nodal values
C
      WRITE (lun,'(A)') '#'
      WRITE (lun,'(A)') '#  Temperature nodal values'
      WRITE (lun,'(A)') '#'
      WRITE (lun,*) mesh%nnodes
      DO i=1,mesh%nnodes
        WRITE (lun,*) T(i)
      END DO

      CLOSE(lun)

      DEALLOCATE (flux)

      END
            
