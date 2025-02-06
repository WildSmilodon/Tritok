      MODULE inc_types

!
! -----------------------------------------------------------------------------------------
!      
      INCLUDE 'mpif.h'
      TYPE penv
        INTEGER          :: comm        ! communicator number
        INTEGER          :: nproc       ! number of processes
        INTEGER          :: mpr         ! my process number (zero based)
        INTEGER          :: myproc      ! my process number (one based)
!
!***  Global data
!
        INTEGER          :: nbnodes     ! number of boundary nodes
        INTEGER          :: nnodes      ! number of total nodes
        INTEGER          :: nbelem      ! number of boundary elements
        INTEGER          :: nicell      ! number of internal cells
        INTEGER          :: npob        ! nodal points of boundary element
        INTEGER          :: npoc        ! nodal points of cell
        INTEGER          :: npx         ! number of dimensions
!
!***  Local data
!
        INTEGER          :: zac         ! my first node
        INTEGER          :: kon         ! my last node
        INTEGER          :: nmn         ! number of my nodes
        INTEGER, POINTER :: imac(:)     ! beginning of my cells
        INTEGER, POINTER :: nmac(:)     ! number of cells
        INTEGER, POINTER :: inod(:)     ! beginning of my nodes
        INTEGER, POINTER :: nnod(:)     ! number of nodes
        INTEGER, POINTER :: imat(:)     ! beginning of my rows in matrix 
        INTEGER, POINTER :: nmat(:)     ! number of rows
        
        INTEGER          :: ncp         ! 
        INTEGER          :: ncex
        INTEGER, POINTER :: ice(:)
        INTEGER, POINTER :: jce(:)
        INTEGER, POINTER :: prop(:)
        INTEGER, POINTER :: beop(:)
        INTEGER, POINTER :: bein(:)
        INTEGER          :: CLst        ! number of elements in CLlist
        INTEGER          :: CLstkom     ! number of communication pairs
        INTEGER, POINTER :: CLlist(:)   ! equations to be communicationed
        INTEGER, POINTER :: CLp2(:)     ! processor 2
        INTEGER, POINTER :: CLzac1(:)   ! start of list in CLlist - sent by processor 1 -> 2
        INTEGER, POINTER :: CLkon1(:)   ! end of list in CLlist
        INTEGER, POINTER :: CLste1(:)   ! number of items to be communicated = CLkon-Clzac+1
        INTEGER, POINTER :: CLzac2(:)   ! start of list in CLlist - sent by processor 2 -> 1
        INTEGER, POINTER :: CLkon2(:)   ! end of list in CLlist
        INTEGER, POINTER :: CLste2(:)   ! number of items to be communicated = CLkon-Clzac+1
      END TYPE
!
! -----------------------------------------------------------------------------------------
!      
      TYPE pslv
        INTEGER    :: LDA                 ! neq
        INTEGER    :: N                   ! neq
        INTEGER    :: LOCLEN              ! peq
        INTEGER    :: BLKSZ
        INTEGER    :: BASIS
        INTEGER    :: NPROCS              ! npr
        INTEGER    :: PROCID              ! mpr
        INTEGER    :: PRECONT
        INTEGER    :: STOPT
        INTEGER    :: MAXIT
        INTEGER    :: ITNO
        INTEGER    :: STATUS
        INTEGER    :: STEPERR
        REAL(8)    :: EPSILON
        REAL(8)    :: EXITNORM
        TYPE(penv) :: ENV
      END TYPE pslv
!
! -----------------------------------------------------------------------------------------
!      
      TYPE meshType
        CHARACTER*255 ::  FileName
        CHARACTER*50  ::  Version
        CHARACTER*255 ::  Path
        CHARACTER*255 ::  Title        
        CHARACTER*255 ::  FullName            ! path & file name
        REAL(8) :: id                         ! unique ID
        INTEGER :: nbelem,nicell,nnodes,nbnodes,npob,npoc,npx !,nx,ny,simx,simy,nofc,
        INTEGER :: nofw
        INTEGER :: nside,npof  ! stevilo stranic = 6, stevilo nodeov na nezvezni stranici =4
        INTEGER :: nq      ! stevilo nezveznih nodeov (vseh skupaj)
        INTEGER :: npofc   !  stevilo nezveznih nodeov v eni celici = nside*npof=24
        INTEGER :: nicnpoc ! nicell * npoc 27
        INTEGER :: nsp     !=mesh%npoc+mesh%npofc   ! 27+24=51
        INTEGER :: nicnsp  !=nicell * nsp
        INTEGER :: nbnpof  ! nbelem * npof
        
        INTEGER, POINTER :: ibc(:,:)  ! nbelem,npob
        INTEGER, POINTER :: idc(:,:)  ! nicell,npoc
        INTEGER, POINTER :: ibf(:,:,:)! nicell,nside,npof
        INTEGER, POINTER :: ibcf(:,:) ! nbelem,npof        
        REAL(8), POINTER :: x(:,:)    ! nnodes,npx
        REAL(8), POINTER :: xq(:,:)   ! nq,npx
        REAL(8), POINTER :: CellVolume(:)   ! nicell, mesh element volume
        REAL(8), POINTER :: nodeVol(:) ! nnodes, volume around a single node
        REAL(8) nodeVolTot
        REAL(8) MeshVolume  ! total volume of the mesh

        INTEGER, POINTER :: ukode(:) ! doloci ali je vozlice znano ali neznano
        INTEGER, POINTER :: qkode(:) ! na katero po vrsti je v vektrju

!       makro vorticity transport equation        
        INTEGER, POINTER :: wkode(:,:) ! doloci ali je vozlice znano ali neznano
        INTEGER, POINTER :: wqkode(:,:) ! na katero po vrsti je v vektrju
        INTEGER, POINTER :: nunk(:),nb(:) ! stevilo neznank in znanih vrednosti

!       makro energy transport equation        
        INTEGER, POINTER :: Tkode(:) ! doloci ali je vozlice znano ali neznano
        INTEGER, POINTER :: Tqkode(:) ! na katero po vrsti je v vektrju
        INTEGER Tnunk,Tnb ! stevilo neznank in znanih vrednosti

!       makro general diff-conv transport equation        
        INTEGER, POINTER :: DCkode(:) ! doloci ali je vozlice znano ali neznano
        INTEGER, POINTER :: DCqkode(:) ! na katero po vrsti je v vektrju
        INTEGER DCnunk,DCnb ! stevilo neznank in znanih vrednosti

!       makro modified Helmholtz  equation        
        INTEGER, POINTER :: mHkode(:) ! doloci ali je vozlice znano ali neznano
        INTEGER, POINTER :: mHqkode(:) ! na katero po vrsti je v vektrju
        INTEGER mHnunk,mHnb ! stevilo neznank in znanih vrednosti
!       za single domain
        INTEGER mHfmnunk
        INTEGER, POINTER :: mHqfmkode(:) ! na katero po vrsti je v vektrju
        
        INTEGER, POINTER :: flpm(:,:) ! flukc plus minus, katera G integral je s plusom kateri z munisim


        INTEGER, POINTER :: gbn(:)    ! boundary node list (nbnodes)
        
        INTEGER, POINTER :: lbn(:)    ! inverse boundary node list (nnodes), 0=obmocno vozlisce

        REAL(8), POINTER :: nx(:),ny(:),nz(:) ! noramala na tocko na robu (nbnodes)
        REAL(8), POINTER :: bcn(:,:) !x(:),bcny(:),bcnz(:) ! noramala na tocko na robu (nbnodes,3) glede na robne pogoje
        REAL(8), POINTER :: dwur(:)  ! dynamic vorticity under-relaxation

        INTEGER, POINTER :: bewn(:)   ! boundary element wall number (nbelem)
        INTEGER, POINTER :: iwn(:)    ! boundary node wall number (nbnodes)

        
!        INTEGER :: nbnodm1,nnmnb,nxpny2        
!        REAL(8), POINTER :: area(:),aleng(:),un(:,:),ut(:,:)
!        REAL(8), POINTER :: nneX(:,:),nneY(:,:) ! normala na robni element celice = (nicell,4)
!        REAL(8), POINTER :: dre(:,:)            ! dolzina robnega elementa celice = (nicell,4)        

!        INTEGER, POINTER :: iwc(:,:)        
         INTEGER, POINTER :: kode(:,:)                
         INTEGER, POINTER :: kobc(:) 
!        REAL(8), POINTER :: kobcKsi(:),kobcEta(:)                    
!        INTEGER, POINTER :: icn(:)        
!        INTEGER, POINTER :: isw(:)        
!        REAL(8), POINTER :: intsize(:)         ! makro mesh interpolation, total size of cells node belongs to
       
!        INTEGER, POINTER :: dnl(:)    ! domain node list (nnodes-nbnodes), 
!        INTEGER, POINTER :: idnl(:)   ! inverse domain node list (nnodes), 0=robno vozlisce

!        Linearna interpolacija omege za desno stran kinematike
         INTEGER :: nnodes8  ! stevilo linearnih vozlisc (samo notranje, brez robnih)
         INTEGER :: npoc8    ! 8
         INTEGER, POINTER :: idc8(:,:)  ! nicell,npoc8
         INTEGER, POINTER :: p278(:)    ! nnodes, kaze na stolpec v linearni matriki
         INTEGER, POINTER :: p827(:)    ! nnodes8, kaze iz linearne mreze nazaj v kvadratno
         INTEGER, POINTER :: lbn8(:)    ! kateri izmed vozlics so na robu
!         INTEGER, POINTER :: bcells(:)  ! ali ima celica kateri node na robu


!       za periodicne robne pogoje
        INTEGER iPeri ! stikalo
        INTEGER pwFrom,pwTo ! kopiram iz From na To, stevilki zidov
        INTEGER nPeri ! stevilo periodicnih nodeov
        INTEGER, POINTER :: pwFT(:,:) ! preslikava - stevilka noda na streni pwFrom, stevilka noda na steni pwTo


!       za kvadratni sistem enacb
      INTEGER, POINTER :: sqUlist(:,:) ! nnodes,8 zato, ker imam lahko najvec osem sosedov
      INTEGER, POINTER :: sqUlistKM(:,:) ! nnodes,8 zato, ker imam lahko najvec osem sosedov (za kinematiko)
      INTEGER, POINTER :: sqUlistIC(:,:) ! nnodes,8 zato, ker imam lahko najvec osem sosedov
      INTEGER, POINTER :: sqUlistNO(:) ! nnodes stevilo enacb za izvorno tocko

      INTEGER, POINTER :: sqQlist(:,:) ! nq,2 zato, ker imam lahko najvec dva soseda
      INTEGER, POINTER :: sqQlistIC(:,:) ! nq,2 zato, ker imam lahko najvec osem sosedov
      INTEGER, POINTER :: sqQlistNO(:) ! nq stevilo enacb za izvorno tocko

      INTEGER, POINTER :: Teql(:)
      INTEGER, POINTER :: Weql(:,:)
      INTEGER, POINTER :: DCeql(:)
      INTEGER, POINTER :: KMeql(:)

!      Rotation matrices
       TYPE(MATRIX3X3), POINTER :: RotMat(:)
       TYPE(MATRIX3X3), POINTER :: RotMatTransp(:)

!      novi bic - .NBC
       INTEGER, POINTER :: bkmc(:) ! boundary kinematics code (1-vort, 2-vel)


!     BCNL - boundary conditions normal list
!     1st value = 1=x,2=y,3=x component pointing towards normal direction
!     2nd and 3rd value = the other two components
      INTEGER, POINTER :: bcnl(:,:)
!     normal vorticity componenta value at the walls
      REAL(8), POINTER :: wnwall(:)


        INTEGER, POINTER  :: neil(:,:) ! list of cell neighbours

      END TYPE meshType 

!
! -----------------------------------------------------------------------------------------
!            

      TYPE KMpointer
        INTEGER :: nsmcol,nrhscol
        INTEGER, POINTER :: vx(:),vy(:),vz(:)
        INTEGER, POINTER :: wx(:),wy(:),wz(:)
        INTEGER, POINTER :: u(:)
        INTEGER, POINTER :: bunt(:),bunn(:)
      END TYPE KMpointer
!
! -----------------------------------------------------------------------------------------
!            
      TYPE InputType      
        CHARACTER*100 IDname,IDversion,IDdate
        CHARACTER*39 StartTime,EndTime      
        CHARACTER*255 mdir,lgeo,lbic,lgeo_withpath,lbic_withpath
        CHARACTER*255 inp_title, bic_title,ctree
        INTEGER prob,copy
!        INTEGER nffx
!        INTEGER ikm,iknw
        INTEGER RSTversion,cTreeFileVers
        INTEGER INTversion
        INTEGER iBw,iDv,iDw,iDT,iDC,imH ! switches EQNS
        REAL(8) ReN ! Reynolds
        REAL(8) RaN ! Rayleigh
        REAL(8) PrN ! Prandtl     
        REAL(8) PeN ! Peclet
        REAL(8) ScN ! Schmidt
        REAL(8) Rro ! Density stability ratio (double diffusion)
        REAL(8) gx,gy,gz ! gravity
        REAL(8) tstep
        INTEGER nstep,iwrit
!        INTEGER nit
        REAL(8) ur,eps
        INTEGER nnlit
!        INTEGER slvt(4),pret(4),prep(4),maxit(4),stopt(4)
!        REAL(4) slveps(4) 

        INTEGER lsqs_maxit,pkms_maxit
        REAL(4) lsqs_eps,pkms_eps
        INTEGER sqrs_type,sqrs_prec,sqrs_maxit
        REAL(4) sqrs_eps
        INTEGER iDLSE
        REAL(4) dlse(9)  ! dynamic least equares epsilon
        REAL(4) dlse_epsth ! dynamic least equares epsilon maximal value
        REAL(4) dlse_r ! dynamic least equares epsilon ratio

        REAL(8) beta,beta2,beta3
!        INTEGER iTimeSH                         ! time scheme type =3   
        INTEGER rst                             ! =1 - read restart file   
        INTEGER scro                            ! output data on screen
        REAL(8) urBw(3)  ! under relaxation Boundary vorticities (kinematics)
        REAL(8) urDw(3)  ! under relaxation Domain vorticities (kinetics)
        REAL(8) urDv(3)  ! under relaxation Domain velocities (kinematics)    
        REAL(8) urDT     ! under relaxation Domain temperatures (kinetics)
        REAL(8) urDC     ! under relaxation general DC equation (kinetics)
        REAL(8) urmh     ! under relaxation modified Helmholtz equation

        INTEGER iDWUR
        REAL(8) dwur_max,dwur_a

!       nanofluids
        INTEGER inano
        REAL(8) nanovf
        REAL(8) nanoF_cp, nanoF_rho, nanoF_k, nanoF_beta
        REAL(8) nanoS_cp, nanoS_rho, nanoS_k, nanoS_beta
        REAL(8) nanoVsh, nanoSigma ! share of fi0, potenca v modelu za razdelitev delcev
        INTEGER nanoMod

!       mixture of forced and natural convection
        INTEGER imixc

!       only estimate RAM
        INTEGER iEstr

!       mesh tranformations
        INTEGER iMTrot,iMTstr,iMTtra
        REAL(8) MTrot(3),MTstr(3),MTtra(3)

!       ferrofluids
        INTEGER iFerro
        REAL(8) Ram  ! magnetic Rayleigh number
        REAL(8) bdt  ! beta * Delta T
        REAL(8) hi0  ! magnetic susceptibility (differential)

!       modified Helmholz equation solver
!       nabla^2 u - mu^2 u = f
        REAL(8) mhmu ! mu
        INTEGER mHf  ! f  (type, forumals in code)
 
!       izvoz rezultatov v paraview (*.vtu)
        INTEGER iPara
!       izvoz rezultatov v paraview (*.vtk)
        INTEGER iParaVTK

!       izvoz rezultatov po stenah (tri.walls.dat)
        INTEGER iRawl

!       energy equation model
        INTEGER EEQM

!       vorticity equation model
        INTEGER VEQM


!       particle tracking, part time step ratio
        INTEGER iPart, nPart,PartTSR

!       Brownian motion
        INTEGER iBrown  ! 1/0
        REAL(8) BrownT0 ! T'=(T-T0)/Delta T
        REAL(8) BrownDomainSize ! [m]
        REAL(8) BrownPartDiam ! Browninan particle diameter [m]
        REAL(8) BrownT0mul ! SQRT(6*Pr/Sc)

!       Thermophoresis
        INTEGER iThph
        REAL(8) thphA,thphB
        REAL(8) thphPartDiam ! Browninan particle diameter [m]
        REAL(8) thphMul ! -A*(d/2r0)^-B*Pr
        REAL(8) thphT0 ! T'=(T-T0)/Delta T
        REAL(8) thphDeltaT ! ghddcgh
        REAL(8) thphDomainSize ! [m]

!       Wavelet compression
        REAL(8) WT_kappa

!       ACA compression
        REAL(8) ACA_eps
        INTEGER ACA_type

!       Cluster trees
        INTEGER ct_type ! 1 = no compression, 2 = WT, 3=ACA
        REAL(8) ct_eta  ! admissiblity criterion

!       test DC equation
        INTEGER iTest

!       time scheme
        INTEGER TimeScheme

!       export stocastih data
        INTEGER estd

!       prescribe analytic channel flow
        INTEGER iacha

      END TYPE       
!
! -----------------------------------------------------------------------------------------
!            
      TYPE GaussType
         INTEGER :: kmBs,kmBr,kmDs,kmDr    ! gre od 1 do 8 = 48
         INTEGER :: iiDiv                  ! integral integral divisions
         REAL(8), POINTER :: gi(:),ome(:)
         INTEGER, POINTER :: ng1(:),ng2(:)
      END TYPE   
!
! -----------------------------------------------------------------------------------------
!            
      TYPE IOtype     
        INTEGER l,inp,errwarn,mesh,interr,bic,init,res,ite,sta,rst,ikm,mas,resf,imi,tim
        INTEGER tfl,wfl,sol,twfl,mass,ikmL,resBIC,lnu,rawl,rms,neil,pres,mres
        INTEGER mHr,imh,imhf,ctree,dct,std
        CHARACTER*30 l_name,inp_name,errwarn_name,interr_name,init_name,res_name,ite_name
        CHARACTER*30 sta_name,rst_name,ikm_name,mas_name,resf_name,imi_name,tim_name
        CHARACTER*30 tfl_name, wfl_name, sol_name, twfl_name, mass_name, ikmL_name
        CHARACTER*30 resBIC_name,lnu_name,rawl_name,rms_name,neil_name,pres_name,mres_name
        CHARACTER*30 mHr_name,imh_name,imhf_name
        CHARACTER*255 ctree_name,dct_name,std_name
      END TYPE   
  
!
! -----------------------------------------------------------------------------------------
!            
      TYPE MATRIX
        INTEGER :: neq,nnz
        INTEGER, POINTER :: i(:),j(:),d(:),ij(:,:)
        REAL(8), POINTER :: v(:)
      END TYPE MATRIX
!
! -----------------------------------------------------------------------------------------
!
      TYPE MATRIX3X3
        REAL(8) :: v(3,3)
      END TYPE MATRIX3X3

!
! -----------------------------------------------------------------------------------------
!            
      TYPE countType
        INTEGER :: tstep         ! global time step number
        INTEGER :: tsnnlit       ! number of nonlinear iterations within a timestep
        INTEGER :: glnnlit       ! global number of nonlinear iterations
        INTEGER :: trtstep       ! number of time steps in this run
        REAL(8) :: rtime         ! real time
        INTEGER :: nit(12)        ! solver iterations Bw, Dv, Dw, DT, DC,mH
        REAL(8) :: errw          ! nonlinear iteration error
        REAL(8) :: nlierr(9)
        INTEGER :: ftime
      END TYPE countType
!
! -----------------------------------------------------------------------------------------
!            
      TYPE CPUtype
        INTEGER ncpu     
        REAL(4) t0,t00
        REAL(4), POINTER :: time(:) 
        CHARACTER*30, POINTER :: desc(:) 
      END TYPE  
!
! -----------------------------------------------------------------------------------------
!

      TYPE particleType
        INTEGER :: id       ! identifikacijska stevilka delca, particle ID number
        REAL(8) :: x,y,z    ! lokacija, particle location
        REAL(8) :: vx,vy,vz ! hitrost delca, particle velocity
        LOGICAL :: active   ! ali je v toku ali ga je odneslo ven ali se je posedel, TRUE=in the domain,
!                                                      FALSE=out of the domain or stuck at the wall
        INTEGER :: cell     ! celica v kateri je, mesh cell number, where particle is located
      END TYPE
!
! -----------------------------------------------------------------------------------------
!
      TYPE fluidAtPartType
        REAL(8) vx,vy,vz  ! hitrost tekocine / fluid velocity
        REAL(8) wx,wy,wz  ! vrtincnost / vorticity
        REAL(8) T,gradT(3)  ! temperature and gradient
        REAL(8) gradVx(3),gradVy(3),gradVz(3) ! velocity gradients
        REAL(8) dvxdt,dvydt,dvzdt ! velocity time derivatives
        REAL(8) nhhx,nhhy,nhhz ! nabla (\vec H \cdot \vec H)
      END TYPE
!
! -----------------------------------------------------------------------------------------
!
      TYPE crsMatrixType
        INTEGER nrow        ! stevilo vrstic
        INTEGER nz          ! stevilo nicel v polni matriki
        INTEGER nnz         ! stevilo nenicelnih elementov
        INTEGER, POINTER :: zvr(:),c(:) ! zvr(nrow+1), c(nnz)
        REAL(8), POINTER :: v(:)  ! v(nnz) vrednosti
      END TYPE

      TYPE MatrixKolazType
        INTEGER ulr,ulc
        INTEGER nrow,ncol
        REAL(8) meanabs     ! povprecna vrednost absolutnih vrednosti vseh elementov v matriki
        REAL(8) meja        ! meja za zanemarjanje
        INTEGER nz          ! stevilo nicel v polni matriki

        INTEGER :: nnz      ! stevilo nenicelnih elementov
        INTEGER, POINTER :: zvr(:),c(:) ! zvr(nrow+1), c(nnz)
        REAL(8), POINTER :: v(:)  ! v(nnz) vrednosti
      END TYPE
      TYPE VectorKolazType
        INTEGER length
        INTEGER start
      END TYPE
!
! -----------------------------------------------------------------------------------------
!
      TYPE acaMatrixType
        INTEGER nrow,ncol,rank  ! stevilo vrstic
        REAL(8), POINTER :: a(:,:),b(:,:) ! a(nrow,rank) b(rank,ncol)
      END TYPE


!
! -----------------------------------------------------------------------------------------
!
      TYPE fullMatrixType
        INTEGER nrow,ncol  ! stevilo vrstic
        REAL(8), POINTER :: v(:,:) ! v(nrow,ncol)
      END TYPE
!
! -----------------------------------------------------------------------------------------
!
      TYPE HmatrixType

        INTEGER Nparts    ! stevilo kosckov na katero je razdeljena matrika
        INTEGER nRow,nCol ! size of the original matrix
        TYPE(HmatrixEntryType), POINTER :: H(:)

      END TYPE


!
! -----------------------------------------------------------------------------------------
!
      TYPE HmatrixEntryType
!                      ! za vsak koscek
        INTEGER  t1v   ! veja iz prvega drevesa
        INTEGER  t2v   ! veja iz drugega drevesa
        INTEGER  nRow  ! st vozlisc v prvem drevesu
        INTEGER  nCol  ! st vozlisc v drugem drevesu
        INTEGER, POINTER :: rNodes(:) ! list of nodes in rows
        INTEGER, POINTER :: cNodes(:) ! list of nodes in cols
        INTEGER  admiss! ali je koscek admissible ali ne
        INTEGER  cType  ! 1=full, 2=wavelet, 3=ACA
        TYPE(acaMatrixType) acaM ! aca compressed matrix
        TYPE(crsMatrixType) crsM ! wavelet compressed matrix
        TYPE(fullMatrixType) fmat ! full (inadmissible) matrix (nRow,nCol)

      END TYPE
!
! -----------------------------------------------------------------------------------------
!
      TYPE HtreeType

        INTEGER Nvej    ! stevilo vej
        INTEGER maxLev  ! najglobji level v drevesu
        TYPE(vejaType), POINTER :: veja(:) ! tree

      END TYPE
!
! -----------------------------------------------------------------------------------------
!
      TYPE vejaType

        INTEGER lev ! level
        INTEGER parent ! parent ID
        INTEGER child(2) ! chile ID's
        INTEGER nchild  ! number of children
        REAL(8) diameter !  estimate of size of cluster of nodes

        INTEGER nelems ! number of mesh cells in this branch
        INTEGER, POINTER :: elems(:) ! list of mesh elements in this branch
        INTEGER nnodes
        INTEGER, POINTER :: nodes(:) ! list of nodes in elems
        INTEGER nq ! stevilo q izvornih tock na robnih elementih
        INTEGER, POINTER :: qnodes(:) ! list of q nodes in elemens

!       data for combined trees
        INTEGER admiss
        INTEGER t1v,t2v  ! veje prvega in drugega drevesa v kombianciji
        REAL(8) dist  ! min distance between t1v and t2v clusters


      END TYPE

      END MODULE
