C -----------------------------------------------------------------------------
      SUBROUTINE SetWallValues(mesh,wall,bctype,a,b,c,d,u)

C
C     $: Sets values to wall nodes
C
C -----------------------------------------------------------------------------
      USE inc_types
      
      TYPE(meshType) :: mesh
      INTEGER wall,bctype,dn
      REAL(8) a,b,c,d
      REAL(8) u(mesh%nnodes)
      
      INTEGER i
      
      IF (bcType.EQ.0) THEN   ! fixed value  u=a
        DO i=1,mesh%nbnodes
          dn=mesh%gbn(i)
          IF (mesh%iwn(i).EQ.wall) THEN
              u(dn)=a
          END IF
        END DO
      ELSE IF (bcType.EQ.1) THEN   ! linear u=a*x+b*y+c*z+d        
        DO i=1,mesh%nbnodes
          dn=mesh%gbn(i)
          IF (mesh%iwn(i).EQ.wall) THEN
            u(dn)=a*mesh%x(dn,1)+b*mesh%x(dn,2)+c*mesh%x(dn,3)+d
          END IF
        END DO     
      ELSE IF (bcType.EQ.2) THEN   ! quadratic po z u=a*z**2+b*z+c+d
        DO i=1,mesh%nbnodes
          dn=mesh%gbn(i)
          IF (mesh%iwn(i).EQ.wall) THEN
            u(dn)=a*mesh%x(dn,3)**2+b*mesh%x(dn,3)+c+d
          END IF
        END DO
      ELSE IF (bcType.EQ.3) THEN   ! developed laminar flow profile (za vx v smeri x)
        DO i=1,mesh%nbnodes
          dn=mesh%gbn(i)
          IF (mesh%iwn(i).EQ.wall) THEN
            CALL Channel3Danal1point(mesh%x(dn,2),mesh%x(dn,3),u(dn),a,b,c,d)
          END IF
        END DO
      END IF
        
      END

C -----------------------------------------------------------------------------
      SUBROUTINE ReadBic(env,io,inp,mesh,velocity,vorticity,qvorticity,temp,qtemp,dc,qdc,mHu,mhq,fmqmH)
C
C     $: Prebere initial and boundary conditions file
C
C -----------------------------------------------------------------------------
      USE inc_types
      
      TYPE(meshType) :: mesh
      TYPE(IOtype) :: io
      TYPE(InputType) inp
      TYPE(penv) :: env
        
      INTEGER lun,i,j,k,kode
      CHARACTER KeyWord*4,OneLine*255,tmp*255 

      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3)
      REAL(8) qvorticity(mesh%nq,3)          
      REAL(8) temp(mesh%nnodes)
      REAL(8) qtemp(mesh%nq)    
      REAL(8) dc(mesh%nnodes)
      REAL(8) qdc(mesh%nq)    
      REAL(8) mHu(mesh%nnodes)
      REAL(8) mHq(mesh%nq)  
      REAL(8) fmqmH(mesh%nbnpof)
            
      INTEGER wall1,wall2,wall3
      INTEGER ff,wall,bctype,smer,dn,ii
      REAL(8) value,value1,value2,value3,value4
     
      lun=io%bic

c     iwn = bounday node wall number
c     na zacetku jih razdelim kar po bewn, nato iz bic preberm in popravim 
c     robove in vogale
      ALLOCATE (mesh%iwn(mesh%nbnodes))
      DO i=1,mesh%nbelem
        DO j=1,mesh%npob
          k=mesh%lbn(mesh%ibc(i,j))
          mesh%iwn(k)=mesh%bewn(i)
        END DO
      END DO      
  
      ALLOCATE (mesh%kode(mesh%nbnodes,3))
      
C
C     Boundary condtition matrices for vorticity transport equation
C      
      ALLOCATE (mesh%wqkode(mesh%nq,3))      
      ALLOCATE (mesh%wkode(mesh%nnodes,3))      
      ALLOCATE (mesh%nunk(3))
      ALLOCATE (mesh%nb(3))
C
C     Boundary condtition matrices for temperature equation
C  
      ALLOCATE (mesh%Tqkode(mesh%nq))      
      ALLOCATE (mesh%Tkode(mesh%nnodes))      
C
C     Boundary condtition matrices for genreal D-C equation
C  
      ALLOCATE (mesh%DCqkode(mesh%nq))      
      ALLOCATE (mesh%DCkode(mesh%nnodes))
C
C     Boundary condtition matrices for modified Helmholtz equation
C  
      ALLOCATE (mesh%mHqkode(mesh%nq))      
      ALLOCATE (mesh%mHkode(mesh%nnodes))
      ALLOCATE (mesh%mHqfmkode(mesh%nbnpof)) ! za full matrix
  
c     na zacektu postavim vse q in u za neznanke v vorticity transport enacbi
      mesh%wkode=0
      mesh%wqkode=0      
      mesh%Tkode=0
      mesh%Tqkode=0      
      mesh%DCkode=0
      mesh%DCqkode=0   
      mesh%mHkode=0
      mesh%mHqkode=0
      mesh%mHqfmkode=0
C     stevilo neznank
      mesh%nunk=0
      mesh%Tnunk=0
      mesh%DCnunk=0
      mesh%mHnunk=0      
C     stevilo znanih
      mesh%nb=0               
      mesh%Tnb=0
      mesh%DCnb=0
      mesh%mHnb=0                     

      velocity=0.0D0
      vorticity=0.0D0      
      qvorticity=0.0D0
      temp=0.0D0
      qtemp=0.0D0
      dc=0.0D0
      qdc=0.0D0
      mHu=0.0D0
      mHq=0.0D0

C     periodicni robni pogoji
      mesh%iPeri=0

      
      OPEN (lun,FILE=trim(inp%lbic_withpath),ERR=10,STATUS='OLD')


      CALL rOneTL(lun,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')
        READ (OneLine,*) KeyWord
C
C____    Edges - distribute nodes on edges among walls
C
        IF (KeyWord.EQ.'EDGE') THEN
            CALL rOneTL(lun,OneLine)
            DO WHILE (OneLine(1:3).NE.'END'.AND.OneLine(1:3).NE.'EOF')
              READ(Oneline,*) wall1,wall2,wall
              CALL bewniwnEdge(mesh,wall1,wall2,wall)
              CALL rOneTL(lun,OneLine)
            END DO             
C
C____    Corners - distribute nodes on corners among walls
C
        ELSE IF (KeyWord.EQ.'CORN') THEN
            CALL rOneTL(lun,OneLine)
            DO WHILE (OneLine(1:3).NE.'END'.AND.OneLine(1:3).NE.'EOF')
              READ(Oneline,*) wall1,wall2,wall3,wall
              CALL bewniwnCorner(mesh,wall1,wall2,wall3,wall)
              CALL rOneTL(lun,OneLine)
            END DO      
C
C____    Periodic BC
C
        ELSE IF (KeyWord.EQ.'PERI') THEN
            CALL rOneTL(lun,OneLine)
            DO WHILE (OneLine(1:3).NE.'END'.AND.OneLine(1:3).NE.'EOF')
              READ(Oneline,*) mesh%pwFrom,mesh%pwTo
              mesh%iPeri=1
              CALL FindPeriodicNodes(mesh,env,io,inp)
              CALL rOneTL(lun,OneLine)
            END DO
C
C____    Initial conditions
C
        ELSE IF (KeyWord.EQ.'INIT') THEN
            CALL rOneTL(lun,OneLine)
            DO WHILE (OneLine(1:3).NE.'END'.AND.OneLine(1:3).NE.'EOF')
              READ(Oneline,*) ff,value
              IF (ff.EQ.10) THEN
                temp(:)=value
              ELSE IF (ff.EQ.11) THEN
                dc(:)=value
              ELSE IF (ff.EQ.12) THEN
                mHu(:)=value
              ELSE IF (ff.GE.7) THEN
                vorticity(:,ff-6)=value
              ELSE IF (ff.GE.4) THEN
                vorticity(:,ff-3)=value
              ELSE
                velocity(:,ff)=value
              END IF
              CALL rOneTL(lun,OneLine)
            END DO             
C
C____    Wall normal component of vorticity
C
        ELSE IF (KeyWord.EQ.'N-VO') THEN
            CALL rOneTL(lun,OneLine)
            DO WHILE (OneLine(1:3).NE.'END'.AND.OneLine(1:3).NE.'EOF')
              READ(Oneline,*) ff,value
              IF (ff.GT.0.AND.ff.LE.mesh%nofw) THEN
                mesh%wnwall(ff)=value
              ELSE
                CALL WarnErr(env,io,inp,1,'ReadBiC','Unknown wall!',0)
              END IF
              CALL rOneTL(lun,OneLine)
            END DO
C
C____    Boundary conditions for signle domain kinematics
C           kode = 0 - calculate w value
C                = 1 - ignore, do nothing
C
C
        ELSE IF (KeyWord.EQ.'X-KI'.OR.KeyWord.EQ.'Y-KI'.OR.KeyWord.EQ.'Z-KI') THEN
            IF (KeyWord.EQ.'X-KI') THEN
              ff=1
            ELSE IF (KeyWord.EQ.'Y-KI') THEN
              ff=2
            ELSE 
              ff=3
            END IF
            CALL rOneTL(lun,OneLine)
            DO WHILE (OneLine(1:3).NE.'END'.AND.OneLine(1:3).NE.'EOF')
              READ(Oneline,*) wall,kode,bctype
              IF (bctype.EQ.0) THEN
                CALL rOneTL(lun,OneLine)
                READ (Oneline,*) value1
              ELSE IF (bctype.EQ.1.or.bctype.EQ.2.or.bctype.EQ.3) THEN
                CALL rOneTL(lun,OneLine)
                READ (Oneline,*) value1,value2,value3,value4
              END IF
              CALL SetWallValues(mesh,wall,bctype,value1,value2,value3,value4,velocity(:,ff))
c             set kode
              DO i=1,mesh%nbnodes
                IF (mesh%iwn(i).EQ.wall) THEN
                  mesh%kode(i,ff)=kode
                END IF
              END DO              
              CALL rOneTL(lun,OneLine)
            END DO             

C
C____    Boundary conditions for vorticity transport equation
C
        ELSE IF (KeyWord.EQ.'X-VO'.OR.KeyWord.EQ.'Y-VO'.OR.KeyWord.EQ.'Z-VO') THEN
            IF (KeyWord.EQ.'X-VO') THEN
              smer=1
            ELSE IF (KeyWord.EQ.'Y-VO') THEN
              smer=2
            ELSE 
              smer=3
            END IF        
            CALL rOneTL(lun,OneLine)
            DO WHILE (OneLine(1:3).NE.'END'.AND.OneLine(1:3).NE.'EOF')        
              READ(Oneline,*) ff,wall,value        
              IF (ff.EQ.1) THEN ! u
                DO i=1,mesh%nbnodes
                  IF (mesh%iwn(i).EQ.wall) THEN
                  dn=mesh%gbn(i)
c                DO i=1,mesh%nbelem
c                  IF (mesh%bewn(i).EQ.wall) THEN
c                    DO j=1,mesh%npob
c                      vorticity(mesh%ibc(i,j),smer)=value    ! nastavim robno vrednost (ki jo potem povozi kinematika)
c                      mesh%nb(smer)=mesh%nb(smer)+1
c                      mesh%wkode(mesh%ibc(i,j),smer)=mesh%nb(smer)    ! nastavim kodo  = znanka
                      vorticity(dn,smer)=value    ! nastavim robno vrednost (ki jo potem povozi kinematika)
                      mesh%nb(smer)=mesh%nb(smer)+1
                      mesh%wkode(dn,smer)=mesh%nb(smer)    ! nastavim kodo  = znanka
c                    END DO
                  END IF
                END DO                
              ELSE IF (ff.EQ.2) THEN ! q
                DO i=1,mesh%nbelem
                  IF (mesh%bewn(i).EQ.wall) THEN
                    DO j=1,mesh%npof
                      qvorticity(mesh%ibcf(i,j),smer)=value    ! nastavim robno vrednost
                      mesh%nb(smer)=mesh%nb(smer)+1
                      mesh%wqkode(mesh%ibcf(i,j),smer)=mesh%nb(smer)    ! nastavim kodo  = znanka
                    END DO
                  END IF
                END DO                                
              END IF
              CALL rOneTL(lun,OneLine)
            END DO             
C
C____    Boundary conditions for energy transport equation
C
        ELSE IF (KeyWord.EQ.'TEMP') THEN
            CALL rOneTL(lun,OneLine)
            DO WHILE (OneLine(1:3).NE.'END'.AND.OneLine(1:3).NE.'EOF')        
              READ(Oneline,*) ff,wall,value        
              IF (ff.EQ.1) THEN ! u
                DO i=1,mesh%nbnodes
                  IF (mesh%iwn(i).EQ.wall) THEN
                    dn=mesh%gbn(i)
                    temp(dn)=value    ! nastavim robno vrednost
                    mesh%Tnb=mesh%Tnb+1
                    mesh%Tkode(dn)=mesh%Tnb    ! nastavim kodo  = znanka
                  END IF
                END DO                
              ELSE IF (ff.EQ.2) THEN ! q
                DO i=1,mesh%nbelem
                  IF (mesh%bewn(i).EQ.wall) THEN
                    DO j=1,mesh%npof
                      qtemp(mesh%ibcf(i,j))=value    ! nastavim robno vrednost
                      mesh%Tnb=mesh%Tnb+1
                      mesh%Tqkode(mesh%ibcf(i,j))=mesh%Tnb    ! nastavim kodo  = znanka
                    END DO
                  END IF
                END DO                                
              END IF
              CALL rOneTL(lun,OneLine)
            END DO     
C
C____    Boundary conditions for general D-C transport equation
C
        ELSE IF (KeyWord.EQ.'GEDC') THEN
            CALL rOneTL(lun,OneLine)
            DO WHILE (OneLine(1:3).NE.'END'.AND.OneLine(1:3).NE.'EOF')        
              READ(Oneline,*) ff,wall,value        
              IF (ff.EQ.1) THEN ! u
                DO i=1,mesh%nbnodes
                  IF (mesh%iwn(i).EQ.wall) THEN
                    dn=mesh%gbn(i)
                    dc(dn)=value    ! nastavim robno vrednost
                    mesh%DCnb=mesh%DCnb+1
                    mesh%DCkode(dn)=mesh%DCnb    ! nastavim kodo  = znanka
                  END IF
                END DO                
              ELSE IF (ff.EQ.2) THEN ! q
                DO i=1,mesh%nbelem
                  IF (mesh%bewn(i).EQ.wall) THEN
                    DO j=1,mesh%npof
                      qdc(mesh%ibcf(i,j))=value    ! nastavim robno vrednost
                      mesh%DCnb=mesh%DCnb+1
                      mesh%DCqkode(mesh%ibcf(i,j))=mesh%DCnb    ! nastavim kodo  = znanka
                    END DO
                  END IF
                END DO                                
              END IF
              CALL rOneTL(lun,OneLine)
            END DO                       
C
C____    Boundary conditions for modified Helmholtz equation
C
        ELSE IF (KeyWord.EQ.'MHEL') THEN
            CALL rOneTL(lun,OneLine)
            DO WHILE (OneLine(1:3).NE.'END'.AND.OneLine(1:3).NE.'EOF')        
              READ(Oneline,*) ff,wall,value        
              IF (ff.EQ.1) THEN ! u
                DO i=1,mesh%nbnodes
                  IF (mesh%iwn(i).EQ.wall) THEN
                    dn=mesh%gbn(i)
                    mHu(dn)=value    ! nastavim robno vrednost
                    mesh%mHnb=mesh%mHnb+1
                    mesh%mHkode(dn)=mesh%mHnb    ! nastavim kodo  = znanka
                  END IF
                END DO                
              ELSE IF (ff.EQ.2) THEN ! q
                ii=0
                DO i=1,mesh%nbelem
                  DO j=1,mesh%npof
                    ii=ii+1
                    IF (mesh%bewn(i).EQ.wall) THEN
                      mHq(mesh%ibcf(i,j))=value    ! nastavim robno vrednost
                      mesh%mHnb=mesh%mHnb+1
                      mesh%mHqkode(mesh%ibcf(i,j))=mesh%mHnb    ! nastavim kodo  = znanka
c                     za polno matriko (fluski samo po zunanjem robu)
                      fmQmH(ii)=value
                      mesh%mHqfmkode(ii)=mesh%mHnb    ! nastavim kodo  = znanka
                    END IF
                  END DO
                END DO                                
              END IF
              CALL rOneTL(lun,OneLine)
            END DO  
C
C____ KEYWORD NOT FOUND :
C
        ELSE
          CALL WarnErr(env,io,inp,1,'ReadBiC','KeyWord not found!',0)
        END IF
        CALL rOneTL(lun,OneLine)
      END DO        
           
      CLOSE(lun)
c
c     vorticity transport enacba
c     Neznanke prestejemo in razporedimo po vektroju neznank
c     oznacimo z negativnim predznakom v ukode in qkode
c       
      DO smer=1,3
        DO i=1,mesh%nnodes
          IF (mesh%wkode(i,smer).EQ.0) THEN
            mesh%nunk(smer)=mesh%nunk(smer)+1
            mesh%wkode(i,smer)=-mesh%nunk(smer)
          END IF
        END DO
        DO i=1,mesh%nq
          IF (mesh%wqkode(i,smer).EQ.0) THEN
            mesh%nunk(smer)=mesh%nunk(smer)+1
            mesh%wqkode(i,smer)=-mesh%nunk(smer)
          END IF
        END DO
      END DO

c
c     energy transport enacba
c     Neznanke prestejemo in razporedimo po vektroju neznank
c     oznacimo z negativnim predznakom v ukode in qkode
c       
      DO i=1,mesh%nnodes
        IF (mesh%Tkode(i).EQ.0) THEN
          mesh%Tnunk=mesh%Tnunk+1
          mesh%Tkode(i)=-mesh%Tnunk
        END IF
      END DO
      DO i=1,mesh%nq
        IF (mesh%Tqkode(i).EQ.0) THEN
          mesh%Tnunk=mesh%Tnunk+1
          mesh%Tqkode(i)=-mesh%Tnunk
        END IF
      END DO

c
c     genreal D-C transport enacba
c     Neznanke prestejemo in razporedimo po vektroju neznank
c     oznacimo z negativnim predznakom v ukode in qkode
c       
      DO i=1,mesh%nnodes
        IF (mesh%DCkode(i).EQ.0) THEN
          mesh%DCnunk=mesh%DCnunk+1
          mesh%DCkode(i)=-mesh%DCnunk
        END IF
      END DO
      DO i=1,mesh%nq
        IF (mesh%DCqkode(i).EQ.0) THEN
          mesh%DCnunk=mesh%DCnunk+1
          mesh%DCqkode(i)=-mesh%DCnunk
        END IF
      END DO
c
c     modified Helmholtz enacba
c     Neznanke prestejemo in razporedimo po vektroju neznank
c     oznacimo z negativnim predznakom v ukode in qkode
c       
      DO i=1,mesh%nnodes
        IF (mesh%mHkode(i).EQ.0) THEN
          mesh%mHnunk=mesh%mHnunk+1
          mesh%mHkode(i)=-mesh%mHnunk
        END IF
      END DO

      mesh%mHfmnunk=mesh%mHnunk

c     za single domain
      DO i=1,mesh%nbnpof
        IF (mesh%mHqfmkode(i).EQ.0) THEN
          mesh%mHfmnunk=mesh%mHfmnunk+1
          mesh%mHqfmkode(i)=-mesh%mHfmnunk
        END IF
      END DO

      DO i=1,mesh%nq
        IF (mesh%mHqkode(i).EQ.0) THEN
          mesh%mHnunk=mesh%mHnunk+1
          mesh%mHqkode(i)=-mesh%mHnunk
        END IF
      END DO
C
C     Creates source point normals based on boundary condtions
C
      CALL BCSourcePointNormals(mesh)
      
      RETURN

10    continue ! error when opening BiC file

c      WRITE (tmp,'(A,A)') "could not open BiC file: ",trim(inp%lbic_withpath)
      WRITE (tmp,'(A)') "could not open BiC file!"      
      CALL WarnErr(env,io,inp,2,"ReadBic",trim(tmp),0)

      END


C -----------------------------------------------------------------------------
      SUBROUTINE bewniwnEdge(mesh,wall1,wall2,wall)
C
C     $: za node, ki so v wall1 in hkrati v wall2 postavi v iwn vrednost wall
C
C -----------------------------------------------------------------------------
      USE inc_types
      
      TYPE(meshType) :: mesh
      INTEGER wall1,wall2,wall
      INTEGER i,j,dn,node,ok1,ok2

      DO node=1,mesh%nbnodes
        dn=mesh%gbn(node)
        ok1=0
        ok2=0
        DO i=1,mesh%nbelem
          DO j=1,mesh%npob
            IF (dn.EQ.mesh%ibc(i,j)) THEN
              IF (mesh%bewn(i).EQ.wall1) THEN
                ok1=1
              ELSE IF (mesh%bewn(i).EQ.wall2) THEN
                ok2=1
              END IF
            END IF
          END DO
        END DO    
        IF (ok1.EQ.1.AND.ok2.EQ.1) THEN ! node je na obeh ploskvah
c          print *,dn,mesh%iwn(node),wall
          mesh%iwn(node)=wall
        END IF
      END DO

      END 

C -----------------------------------------------------------------------------
      SUBROUTINE bewniwnCorner(mesh,wall1,wall2,wall3,wall)
C
C     $: za node, ki so v wall1 in hkrati v wall2 in wall3 postavi v iwn vrednost wall
C
C -----------------------------------------------------------------------------
      USE inc_types
      
      TYPE(meshType) :: mesh
      INTEGER wall1,wall2,wall3,wall
      INTEGER i,j,dn,node,ok1,ok2,ok3

      DO node=1,mesh%nbnodes
        dn=mesh%gbn(node)
        ok1=0
        ok2=0
        ok3=0
        DO i=1,mesh%nbelem
          DO j=1,mesh%npob
            IF (dn.EQ.mesh%ibc(i,j)) THEN
              IF (mesh%bewn(i).EQ.wall1) THEN
                ok1=1
              ELSE IF (mesh%bewn(i).EQ.wall2) THEN
                ok2=1
              ELSE IF (mesh%bewn(i).EQ.wall3) THEN
                ok3=1
              END IF
            END IF
          END DO
        END DO    
        IF (ok1.EQ.1.AND.ok2.EQ.1.AND.ok3.EQ.1) THEN ! node je vogalu treh ploskev
c          print *,dn,mesh%iwn(node),wall
          mesh%iwn(node)=wall
        END IF
      END DO

      END 


C -----------------------------------------------------------------------------
      SUBROUTINE FindPeriodicNodes(mesh,env,io,inp)
C
C     $: poisce vozlisca, ki se bodo kopirala za periodicne B.C.
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(IOtype) :: io
      TYPE(InputType) inp
      TYPE(penv) :: env

      INTEGER ifrom,ito,i,j,k,n
      REAL(8) eps

      eps=1.0D-10

C     najprej prestejemo stevilo vozlisc
      ifrom=0
      ito=0
      DO i=1,mesh%nbnodes
        IF (mesh%iwn(i).EQ.mesh%pwFrom) ifrom=ifrom+1
        IF (mesh%iwn(i).EQ.mesh%pwTo) ito=ito+1
      END DO

      IF (ifrom.NE.ito) THEN ! ERROR
        CALL WarnErr(env,io,inp,2,"FindPeriodicNodes","Periodic walls do not match!",0)
      END IF

      mesh%nPeri=ito
c     alociramo
      ALLOCATE (mesh%pwFT(mesh%nPeri,2))
      mesh%pwFT=-1 ! za test
      ito=0
c     iscemo
      DO i=1,mesh%nbnodes
        IF (mesh%iwn(i).EQ.mesh%pwFrom) THEN
          ito=ito+1
          mesh%pwFT(ito,1)=mesh%gbn(i)
            DO j=1,mesh%nbnodes
              IF (mesh%iwn(j).EQ.mesh%pwTo) THEN
                 n=0
c                dve od treh koordinat morata biti enaki, to bi bilo bolje s kaksnim smernim vektorjem
                 DO k=1,3
                   IF (abs(mesh%x(mesh%gbn(i),k)-mesh%x(mesh%gbn(j),k)).LT.eps) n=n+1
                 END DO
                 IF (n.EQ.2) mesh%pwFT(ito,2)=mesh%gbn(j)
              END IF
            END DO
        END IF
      END DO

      DO i=1,mesh%nPeri
        DO j=1,2
          IF (mesh%pwFT(i,j).LT.1) CALL WarnErr(env,io,inp,2,"FindPeriodicNodes","Periodic walls do not match!",0)
        END DO
c        print *,mesh%pwFT(i,1),mesh%pwFT(i,2)
      END DO

      END
