C -----------------------------------------------------------------------------
      SUBROUTINE RotateByR(x,R,n)
C
C     $: rotates n vectors by matrix R
C
C -----------------------------------------------------------------------------
      INTEGER i,n
      REAL(8) x(n,3),R(3,3),v(3),vr(3)

      DO i=1,n
        v(1)=x(i,1)
        v(2)=x(i,2)
        v(3)=x(i,3)
        vr=MATMUL(R,v)
        x(i,1)=vr(1)
        x(i,2)=vr(2)
        x(i,3)=vr(3)
      END DO


      END

C -----------------------------------------------------------------------------
      SUBROUTINE RotateVecByR(x,y,R,n)
C
C     $: rotates n vectors by matrix R, stores results in y
C
C -----------------------------------------------------------------------------
      INTEGER i,n
      REAL(8) x(n,3),y(n,3),R(3,3),v(3),vr(3)

      DO i=1,n
        v(1)=x(i,1)
        v(2)=x(i,2)
        v(3)=x(i,3)
        vr=MATMUL(R,v)
        y(i,1)=vr(1)
        y(i,2)=vr(2)
        y(i,3)=vr(3)
      END DO

      END



C -----------------------------------------------------------------------------
      SUBROUTINE AnalyseSourcePointNormals(mesh)
C
C     $: for each source point normal, calculate a rotation matrix
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh

      REAL(8) n(3)
      INTEGER i,nmax

C     BCNL - boundary conditions normal list
C     1st value = 1=x,2=y,3=x component pointing towards normal direction
C     2nd and 3rd value = the other two components
      ALLOCATE (mesh%bcnl(mesh%nbnodes,3))

      DO i=1,mesh%nbnodes
        n(1)=ABS(mesh%bcn(i,1))
        n(2)=ABS(mesh%bcn(i,2))
        n(3)=ABS(mesh%bcn(i,3))
C       sort
        nmax=-1
        IF (n(1).GE.n(2).AND.n(1).GE.n(3)) THEN
          nmax=1
          mesh%bcnl(i,1)=nmax
          mesh%bcnl(i,2)=2
          mesh%bcnl(i,3)=3
        ELSE IF (n(2).GE.n(1).AND.n(2).GE.n(3)) THEN
          nmax=2
          mesh%bcnl(i,1)=nmax
          mesh%bcnl(i,2)=1
          mesh%bcnl(i,3)=3
        ELSE IF (n(3).GE.n(1).AND.n(3).GE.n(2)) THEN
          nmax=3
          mesh%bcnl(i,1)=nmax
          mesh%bcnl(i,2)=1
          mesh%bcnl(i,3)=2
        END IF
        IF (nmax.LT.0) THEN
          PRINT *,"ERROR IN AnalyseSourcePointNormals"
          STOP
        END IF
      END DO

      END

C -----------------------------------------------------------------------------
      SUBROUTINE CalculateRotationMatrices(mesh)
C
C     $: for each source point normal, calculate a rotation matrix
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh

      REAL(8) n(3)
      INTEGER i,j,k

      ALLOCATE (mesh%RotMat(mesh%nbnodes))
      ALLOCATE (mesh%RotMatTransp(mesh%nbnodes))

      DO i=1,mesh%nbnodes
        n(1)=mesh%bcn(i,1)
        n(2)=mesh%bcn(i,2)
        n(3)=mesh%bcn(i,3)
        CALL GetRotationMatrix(mesh%RotMat(i)%v,mesh%RotMatTransp(i)%v,n)

        print *,i
        print *,n(1),n(2),n(3)
        print *,mesh%RotMat(i)%v(1,1),mesh%RotMat(i)%v(1,2),mesh%RotMat(i)%v(1,3)
c        print *,n
c        print *,MATMUL(mesh%RotMat(i)%v,n)
c        print *,MATMUL(mesh%RotMatTransp(i)%v,MATMUL(mesh%RotMat(i)%v,n))
c        do j=1,3
c          do k=1,3
c             print *,j,k,mesh%RotMat(i)%v(j,k),mesh%RotMatTransp(i)%v(j,k)
c             end do
c             end do
c        print *,mesh%RotMat(i)%v
c        stop
      END DO

      END


C -----------------------------------------------------------------------------
      SUBROUTINE CheckLinearGeometry(mesh)
C
C     $: imamo 27 tockovno mrezo, vendar je geometrija samo 8 tockovna
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh

      REAL(8) lc1(27),lc2(27),lc3(27) ! lokalne koordinate tock v mrezi
      REAL(8) fi(8) ! interpolacijske funcije za 8 tockovno geometrijo
      REAL(8) x(8),y(8),z(8),xx,yy,zz,nap
      INTEGER idc8(8),ic,i,j

      idc8(1)=6
      idc8(2)=7
      idc8(3)=2
      idc8(4)=3
      idc8(5)=5
      idc8(6)=8
      idc8(7)=1
      idc8(8)=4
      nap=0.0D0

      CALL SetLC(lc1,lc2,lc3)

      DO ic=1,mesh%nicell
        DO i=1,8
          x(i)=mesh%x(mesh%idc(ic,idc8(i)),1)
          y(i)=mesh%x(mesh%idc(ic,idc8(i)),2)
          z(i)=mesh%x(mesh%idc(ic,idc8(i)),3)
        END DO

C       Izracunam, kje bi morale biti tocke v mrezi, ce bi geometrija bila 8 tockovna
        DO i=9,27
          CALL cshape8(lc1(i),lc2(i),lc3(i),fi,8)
          xx=0.0D0
          yy=0.0D0
          zz=0.0D0
          DO j=1,8
            xx=xx+fi(j)*x(j)
            yy=yy+fi(j)*y(j)
            zz=zz+fi(j)*z(j)
          END DO

c        print *,xx,yy,zz
c        print *,mesh%x(mesh%idc(ic,i),1),mesh%x(mesh%idc(ic,i),2),mesh%x(mesh%idc(ic,i),3)
c        print *,xx-mesh%x(mesh%idc(ic,i),1),yy-mesh%x(mesh%idc(ic,i),2),zz-mesh%x(mesh%idc(ic,i),3)
          nap=nap+abs(xx-mesh%x(mesh%idc(ic,i),1))+abs(yy-mesh%x(mesh%idc(ic,i),2))+abs(zz-mesh%x(mesh%idc(ic,i),3))

C         Popravim mrezo
          mesh%x(mesh%idc(ic,i),1)=xx
          mesh%x(mesh%idc(ic,i),2)=yy
          mesh%x(mesh%idc(ic,i),3)=zz

        END DO
c      print *,ic, nap

      END DO
c      print *,nap

      END





C -----------------------------------------------------------------------------
      SUBROUTINE RotateMesh(mesh)
C
C     $: zavrti mrezo okoli x osi
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh

      REAL(8) novy,novz,theta
      INTEGER i

      theta=0.5D0

      DO i=1,mesh%nnodes
c        theta = 1.0D0*mesh%x(i,3)
        novy=cos(theta)*mesh%x(i,2)-sin(theta)*mesh%x(i,3)
        novz=sin(theta)*mesh%x(i,2)+cos(theta)*mesh%x(i,3)
        mesh%x(i,2)=novy
        mesh%x(i,3)=novz
      END DO

      END



C -----------------------------------------------------------------------------
      SUBROUTINE FillKobc(mesh)
C
C     $: za vsako robno vozlisce najdem najblizje vozlisce, ki
C        ni na isti steni in ima normalno razdaljo vecjo od nic. 
C        To potrebujem za dv/dn=0 robni pogoj
C
C -----------------------------------------------------------------------------
      USE inc_types
      
      TYPE(meshType) :: mesh
      
      INTEGER i,j,dn,ok
      REAL(8) dist,mindist,eps,ndist
      
      ALLOCATE (mesh%kobc(mesh%nbnodes))
      mesh%kobc=-1
      eps=1.0D-5 ! za ne dobim vozlisc v isti ravnini
           
      DO i=1,mesh%nbnodes
        dn=mesh%gbn(i)
        mindist=1.0D20
        DO j=1,mesh%nnodes
          ok=0
          IF (mesh%lbn(j).EQ.0) THEN
            ok=1
          ELSE IF (mesh%iwn(i).NE.mesh%iwn(mesh%lbn(j))) THEN 
             ndist=SQRT( (mesh%nx(i)*(mesh%x(dn,1)-mesh%x(j,1)))**2       
     &                  +(mesh%ny(i)*(mesh%x(dn,2)-mesh%x(j,2)))**2
     &                  +(mesh%nz(i)*(mesh%x(dn,3)-mesh%x(j,3)))**2 )
c             if (dn.eq.6)              print *,ndist,j
            IF (ndist.GT.eps) THEN
              ok=1
            END IF
          END IF
          IF (ok.EQ.1) THEN
C           vozlisce j je ali v obmocju ali na drugi steni kot i
            dist=SQRT((mesh%x(dn,1)-mesh%x(j,1))**2+(mesh%x(dn,2)-mesh%x(j,2))**2+(mesh%x(dn,3)-mesh%x(j,3))**2)
            IF (dist.LT.mindist) THEN
              mindist=dist
              mesh%kobc(i)=j
            END IF
          END IF
        END DO
c        print *,dn,mesh%kobc(i)
      END DO
      
      END

C -----------------------------------------------------------------------------
      SUBROUTINE ReadMeshGid3D(env,io,inp,mesh,gp)
C
C     $: reads mesh
C
C -----------------------------------------------------------------------------
      USE inc_types
      
      TYPE(meshType) :: mesh
      TYPE(IOtype) :: io
      TYPE(InputType) inp
      TYPE(penv) :: env      
      TYPE(GaussType) :: gp
        
      INTEGER lun,i,j,k,kk,dn,be
      INTEGER, ALLOCATABLE :: itmp(:)
      CHARACTER KeyWord*4,OneLine*255,tmp*255 
      
      lun=io%mesh
      
      OPEN (lun,FILE=trim(mesh%fullname),ERR=10,STATUS='OLD') !,SHARED)

c
C***  READ VERSION NUMBER OR TITLE OF GEO INPUT FILE :
c
      CALL rOneTL(lun,OneLine)
      READ (OneLine,*) KeyWord
      mesh%version=TRIM(OneLine(6:len_trim(OneLine)))

      CALL rOneTL(lun,OneLine)
      mesh%title=OneLine(1:len_trim(OneLine))
c
C***  READ MESH SIZE :
c
      CALL rOneTL(lun,OneLine)
      READ(OneLine,*) mesh%nbelem,mesh%nicell,mesh%nnodes,mesh%nbnodes,mesh%npob,mesh%npoc,mesh%npx
      IF (mesh%npx.NE.3) THEN
        CALL WarnErr(env,io,inp,4,"ReadMeshGid3D","only npx=3 supported!",0)
      ELSE
        READ(OneLine,*) mesh%nbelem,mesh%nicell,mesh%nnodes,mesh%nbnodes,mesh%npob,mesh%npoc,mesh%npx
c     &                  ,mesh%nx,mesh%ny,mesh%simx,mesh%simy,mesh%nofc,mesh%nofw
      END IF
c      mesh%nbnodm1=mesh%nbnodes-1
c      mesh%nnmnb=mesh%nnodes-mesh%nbnodes
c      mesh%nxpny2=(mesh%nx+mesh%ny)/2

      ALLOCATE (mesh%ibc(mesh%nbelem,mesh%npob),mesh%idc(mesh%nicell,mesh%npoc))
      ALLOCATE (mesh%x(mesh%nnodes,mesh%npx))
      ALLOCATE (mesh%bewn(mesh%nbelem))
      ALLOCATE (mesh%gbn(mesh%nbnodes),mesh%lbn(mesh%nnodes))
      ALLOCATE (mesh%CellVolume(mesh%nicell))

      
c      ALLOCATE (mesh%area(mesh%nicell),mesh%aleng(mesh%nbelem))
c      ALLOCATE (mesh%un(mesh%nbnodes,mesh%npx),mesh%ut(mesh%nbnodes,mesh%npx))

c      ALLOCATE (mesh%iwn(mesh%nbnodes),mesh%iwc(mesh%nofw,mesh%npob))
c      ALLOCATE (mesh%icn(mesh%nofw))
c      ALLOCATE (mesh%dnl(mesh%nnmnb),mesh%idnl(mesh%nnodes))  
c      ALLOCATE (mesh%isw(mesh%nofw))    
c      ALLOCATE (mesh%kobc(mesh%nbnodes))            
c
C____ READ IN MESH POINT COORDINATES:
c
      mesh%id=0.0D00
      DO i=1,mesh%nnodes
        CALL rOneTL(lun,OneLine)
        READ (OneLine,*) k,(mesh%x(k,j),j=1,mesh%npx)
        DO j=1,mesh%npx
          mesh%id=mesh%id+SQRT(ABS(mesh%x(k,j)))
        END DO        
      END DO
C
C____ READ BOUNDARY ELEMENT CONNECTIVITY  :
C
      DO i=1,mesh%nbelem
        CALL rOneTL(lun,OneLine)
        READ (OneLine,*) k,(mesh%ibc(k,j),j=1,mesh%npob)
      END DO                  
C
C____ READ IN INTERNAL CELL CONNECTIVITY  :
C
      DO i=1,mesh%nicell
        CALL rOneTL(lun,OneLine)
        READ (OneLine,*) k,(mesh%idc(k,j),j=1,mesh%npoc)
      END DO
C
C____ FORM BOUNDARY NODES VECTOR
C       
      k=0
      mesh%gbn=0
      mesh%lbn=0
      DO i=1,mesh%nbelem
        DO j=1,mesh%npob
          dn=mesh%ibc(i,j)
          DO kk=1,k
            IF (mesh%gbn(kk).EQ.dn) GOTO 100
          END DO
          k=k+1
          mesh%lbn(dn)=k
          mesh%gbn(k)=dn
100      END DO
      END DO
C
C____ READ BOUNDARY ELEMENT WALL NUMBER :
C
      DO i=1,mesh%nbelem
        CALL rOneTL(lun,OneLine)
        READ (OneLine,*) be,mesh%bewn(be)
      END DO
C
C____ NUMBER OF CELLS TIMES NUMBER OF POINTS IN A CELL
C          
      mesh%nicnpoc=mesh%nicell*mesh%npoc

      CLOSE (lun)

C
C     Transform mesh (rotation, stretching, translation)
C
      IF (inp%iMTstr.EQ.1) CALL MeshStretch(inp,mesh)
      IF (inp%iMTrot.EQ.1) CALL MeshRotate(inp,mesh)
      IF (inp%iMTtra.EQ.1) CALL MeshTranslate(inp,mesh)
C
C     We have linear geometry, check if mesh ok, or correct
C
      CALL CheckLinearGeometry(mesh)

C
C     Calculate normal in every boundary point
C
      CALL SourcePointNormals(mesh)
C
C     find number of walls
C
      mesh%nofw=MAXVAL(mesh%bewn)
C
C     normal component of vorticity at the wall
C
      ALLOCATE (mesh%wnwall(mesh%nofw))
      mesh%wnwall=0.0D0
C
C     Calculate element volumes
C
      CALL CalculateMeshElementVolume(mesh,gp)
c
c     Calculate dynamic omega underrelaxation
c
      IF (inp%iDWUR.GT.0) THEN
        ALLOCATE (mesh%dwur(mesh%nnodes))
        CALL CalculateDWUR(mesh,inp)
      END IF
C
C     Podatki, ki jih rabim, da v kinematiki na desni linearno interpoliram omego
c     iBw=2
C
      mesh%npoc8=8
      ALLOCATE (mesh%idc8(mesh%nicell,mesh%npoc8))
      DO i=1,mesh%nicell
        mesh%idc8(i,1)=mesh%idc(i,6)
        mesh%idc8(i,2)=mesh%idc(i,7)
        mesh%idc8(i,3)=mesh%idc(i,2)
        mesh%idc8(i,4)=mesh%idc(i,3)
        mesh%idc8(i,5)=mesh%idc(i,5)
        mesh%idc8(i,6)=mesh%idc(i,8)
        mesh%idc8(i,7)=mesh%idc(i,1)
        mesh%idc8(i,8)=mesh%idc(i,4)
      END DO
      ALLOCATE (itmp(mesh%nnodes))
      ALLOCATE (mesh%p278(mesh%nnodes))
      mesh%p278=0
      itmp=0
      DO i=1,mesh%nicell
        DO j=1,mesh%npoc8
c         odstranim tiste, ki so na robu. Ali je vozlisce mesh%idc8(i,j) na robu?
          IF (mesh%lbn(mesh%idc8(i,j)).EQ.0) THEN
            itmp(mesh%idc8(i,j))=1
          END IF
        END DO
      END DO
      mesh%nnodes8=0
      j=0
      DO i=1,mesh%nnodes
        mesh%nnodes8=mesh%nnodes8+itmp(i)
        IF (itmp(i).NE.0) THEN
          j=j+1
          mesh%p278(i)=j
        END IF
      END DO
      DEALLOCATE(itmp)

      ALLOCATE (mesh%p827(mesh%nnodes8))
      j=0
      DO i=1,mesh%nnodes
        IF (mesh%p278(i).GT.0) THEN
          j=j+1
          mesh%p827(j)=i
        END IF
      END DO

      ALLOCATE (mesh%lbn8(mesh%nnodes8))
      mesh%lbn8=0
      DO i=1,mesh%nnodes
        IF (mesh%p278(i).GT.0) THEN ! torej je i-vozlisce v linearni mrezi
          IF (mesh%lbn(i).NE.0) THEN ! je tudi na robu
            mesh%lbn8(mesh%p278(i))=1
          END IF
        END IF
      END DO

c      ALLOCATE (mesh%bcells(mesh%nicell)) ! ali ima celica vozlisca na robu
c      mesh%bcells=0
c      DO i=1,mesh%nicell
c        k=0
c        DO j=1,mesh%npoc
c          IF (mesh%lbn(mesh%idc(i,j)).NE.0) k=1
c        END DO
c        IF (k.EQ.1) mesh%bcells(i)=1
c      END DO

      RETURN
      
10    continue ! error when opening mesh
      WRITE (tmp,'(A,A)') "could not open mesh: ",trim(mesh%fullname)
      CALL WarnErr(env,io,inp,2,"ReadMeshGid3D",trim(tmp),0)

      END




C -----------------------------------------------------------------------------
      SUBROUTINE MeshStretch(inp,mesh)
C
C     $: stretches mesh
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(InputType) inp
      INTEGER i,j

      DO i=1,mesh%nnodes
        DO j=1,mesh%npx
          mesh%x(i,j)=mesh%x(i,j)*inp%MTstr(j)
        END DO
      END DO

      END

C -----------------------------------------------------------------------------
      SUBROUTINE MeshTranslate(inp,mesh)
C
C     $: translates mesh
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(InputType) inp
      INTEGER i,j

      DO i=1,mesh%nnodes
        DO j=1,mesh%npx
          mesh%x(i,j)=mesh%x(i,j)+inp%MTtra(j)
        END DO
      END DO

      END



C -----------------------------------------------------------------------------
      SUBROUTINE MeshRotate(inp,mesh)
C
C     $: rotates mesh so that (1,0,0) points to inp%MTrot
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(InputType) inp
      INTEGER i,j

      REAL(8) R(3,3),RT(3,3),v(3),vr(3)

      CALL NormVector(inp%MTrot)
      CALL GetRotationMatrix(R,RT,inp%MTrot)

      DO i=1,mesh%nnodes
        v(1)=mesh%x(i,1)
        v(2)=mesh%x(i,2)
        v(3)=mesh%x(i,3)
        vr=MATMUL(RT,v)
        mesh%x(i,1)=vr(1)
        mesh%x(i,2)=vr(2)
        mesh%x(i,3)=vr(3)
      END DO

      END


C -----------------------------------------------------------------------------
      SUBROUTINE CalculateDWUR(mesh,inp)
C
C     Calculates dynamic omega underrelaxation
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(InputType) inp

      TYPE(meshType) :: mesh
      REAL(8) minim,maxim,dist
      REAL(8) b,c
      INTEGER i,j,k

      DO i=1,mesh%nnodes
        IF (mesh%lbn(i).EQ.0) THEN ! tocka je noter
          minim=1.0D10
          DO k=1,mesh%nbnodes
            j=mesh%gbn(k)
            dist=SQRT( (mesh%x(i,1)-mesh%x(j,1))**2+
     &                 (mesh%x(i,2)-mesh%x(j,2))**2+
     &                 (mesh%x(i,3)-mesh%x(j,3))**2)
            IF (dist.LT.minim) THEN
              minim=dist
              mesh%dwur(i)=minim
            END IF
          END DO
        ELSE
          mesh%dwur(i)=0.0D0
        END IF
      END DO

c     apply underrelaxation based on wall distance

      IF (inp%iDWUR.EQ.1) THEN
        DO i=1,mesh%nnodes
          IF (mesh%lbn(i).EQ.0) THEN ! tocka je noter
            IF (mesh%dwur(i).LT.inp%dwur_a) THEN
              mesh%dwur(i)=inp%dwur_max
            ELSE
              mesh%dwur(i)=1.0D0
            END IF
          END IF
        END DO

      ELSE IF (inp%iDWUR.EQ.2) THEN
c     mormmiram za max razdaljo, ki naj bo ena, ostale manj.
        maxim=-1.0D0
        DO i=1,mesh%nnodes
          IF (mesh%dwur(i).GT.maxim) maxim=mesh%dwur(i)
        END DO
        DO i=1,mesh%nnodes
          mesh%dwur(i)=mesh%dwur(i)/maxim
        END DO
        b = (1.0D0 - inp%dwur_max)/(Exp(-inp%dwur_a) - 1.0D0)
        c = inp%dwur_max - b
        DO i=1,mesh%nnodes
          IF (mesh%lbn(i).EQ.0) THEN ! tocka je noter
            mesh%dwur(i)= c + b*EXP(-inp%dwur_a*mesh%dwur(i))
          END IF
        END DO
      END IF

!      print *,b,c
!      j=0
!      do i=1,mesh%nnodes
!      if (mesh%dwur(i).gt.1.01) then
!        j=j+1 !print *,mesh%dwur(i)
!        print *,mesh%x(i,1),mesh%x(i,2),mesh%x(i,3)
!        endif
!      end do
!      print *,j
!      stop

      END

C -----------------------------------------------------------------------------
      SUBROUTINE SourcePointNormals(mesh)
C
C     Calculates normals in source (boundary) points
C 
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh      

      REAL(8) x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
      REAL(8) eti,etj,eta1m,eta2m,eta1p,eta2p
      REAL(8) xks,yks,zks,xet,yet,zet
      REAL(8), ALLOCATABLE :: anx(:),any(:),anz(:)
      REAL(8) ajac
      
      INTEGER je,i,j,k,dn,n

c     normale na robne elemente      
      ALLOCATE (anx(mesh%nbelem))
      ALLOCATE (any(mesh%nbelem))
      ALLOCATE (anz(mesh%nbelem))      

C     normalo racunam v tocki xi=0,eta=0
      ETI=0.0D0
      ETJ=0.0D0
      ETA1M=1.D0-ETI
      ETA1P=1.D0+ETI
      ETA2M=1.D0-ETJ
      ETA2P=1.D0+ETJ

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
        
        XKS=0.25D0*(-ETA2M*X1+ETA2M*X2+ETA2P*X3-ETA2P*X4)
        YKS=0.25D0*(-ETA2M*Y1+ETA2M*Y2+ETA2P*Y3-ETA2P*Y4)
        ZKS=0.25D0*(-ETA2M*Z1+ETA2M*Z2+ETA2P*Z3-ETA2P*Z4)
C
        XET=0.25D0*(-ETA1M*X1-ETA1P*X2+ETA1P*X3+ETA1M*X4)
        YET=0.25D0*(-ETA1M*Y1-ETA1P*Y2+ETA1P*Y3+ETA1M*Y4)
        ZET=0.25D0*(-ETA1M*Z1-ETA1P*Z2+ETA1P*Z3+ETA1M*Z4)
c
        ANX(je)=YKS*ZET-YET*ZKS
        ANY(je)=XET*ZKS-XKS*ZET
        ANZ(je)=XKS*YET-XET*YKS
c
        AJAC=1.0D0/SQRT(ANX(je)**2+ANY(je)**2+ANZ(je)**2)
c
        ANX(je)=ANX(je)*AJAC
        ANY(je)=ANY(je)*AJAC
        ANZ(je)=ANZ(je)*AJAC

      END DO
      
      ALLOCATE (mesh%nx(mesh%nbnodes))
      ALLOCATE (mesh%ny(mesh%nbnodes))
      ALLOCATE (mesh%nz(mesh%nbnodes))    
      mesh%nx=0.0D0        
      mesh%ny=0.0D0
      mesh%nz=0.0D0      
      
      DO i=1,mesh%nbnodes
        dn=mesh%gbn(i)
        n=0
        DO j=1,mesh%nbelem
          DO k=1,mesh%npob
            IF (mesh%ibc(j,k).EQ.dn) THEN
              n=n+1  ! stevilo robnih elementov na katerih je tocka i
              mesh%nx(i)=mesh%nx(i)+anx(j)
              mesh%ny(i)=mesh%ny(i)+any(j)
              mesh%nz(i)=mesh%nz(i)+anz(j)                            
            END IF
          END DO
        END DO
c        mesh%nx(i)=mesh%nx(i)/DBLE(n)  ! povprecna smer, ne glede na velikost ploskev
c        mesh%ny(i)=mesh%ny(i)/DBLE(n)
c        mesh%nz(i)=mesh%nz(i)/DBLE(n)
c        mesh%nx(i)=1.0D0
c        mesh%ny(i)=1.0D0
c        mesh%nz(i)=1.0D0

        AJAC=1.0D0/SQRT(mesh%nx(i)**2+mesh%ny(i)**2+mesh%nz(i)**2)
        mesh%nx(i)=mesh%nx(i)*AJAC  ! povprecna smer, ne glede na velikost ploskev
        mesh%ny(i)=mesh%ny(i)*AJAC
        mesh%nz(i)=mesh%nz(i)*AJAC
c        print *,i,n
c        print *,mesh%x(dn,1),mesh%x(dn,2),mesh%x(dn,3)
c        print *,mesh%nx(i),mesh%ny(i),mesh%nz(i)

      END DO
      
      DEALLOCATE (anx,any,anz)
      
      END



C -----------------------------------------------------------------------------
      SUBROUTINE BCSourcePointNormals(mesh)
C
C     Calculates normals in source (boundary) points, takes boundary conditions into accout
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh

      REAL(8) x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
      REAL(8) eti,etj,eta1m,eta2m,eta1p,eta2p
      REAL(8) xks,yks,zks,xet,yet,zet
      REAL(8), ALLOCATABLE :: anx(:),any(:),anz(:)
      REAL(8) ajac

      INTEGER je,i,j,k,dn,n

c     normale na robne elemente
      ALLOCATE (anx(mesh%nbelem))
      ALLOCATE (any(mesh%nbelem))
      ALLOCATE (anz(mesh%nbelem))

C     normalo racunam v tocki xi=0,eta=0
      ETI=0.0D0
      ETJ=0.0D0
      ETA1M=1.D0-ETI
      ETA1P=1.D0+ETI
      ETA2M=1.D0-ETJ
      ETA2P=1.D0+ETJ

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

        XKS=0.25D0*(-ETA2M*X1+ETA2M*X2+ETA2P*X3-ETA2P*X4)
        YKS=0.25D0*(-ETA2M*Y1+ETA2M*Y2+ETA2P*Y3-ETA2P*Y4)
        ZKS=0.25D0*(-ETA2M*Z1+ETA2M*Z2+ETA2P*Z3-ETA2P*Z4)
C
        XET=0.25D0*(-ETA1M*X1-ETA1P*X2+ETA1P*X3+ETA1M*X4)
        YET=0.25D0*(-ETA1M*Y1-ETA1P*Y2+ETA1P*Y3+ETA1M*Y4)
        ZET=0.25D0*(-ETA1M*Z1-ETA1P*Z2+ETA1P*Z3+ETA1M*Z4)
c
        ANX(je)=YKS*ZET-YET*ZKS
        ANY(je)=XET*ZKS-XKS*ZET
        ANZ(je)=XKS*YET-XET*YKS
c
        AJAC=1.0D0/SQRT(ANX(je)**2+ANY(je)**2+ANZ(je)**2)
c
        ANX(je)=ANX(je)*AJAC
        ANY(je)=ANY(je)*AJAC
        ANZ(je)=ANZ(je)*AJAC

      END DO

c        INTEGER, POINTER :: bewn(:)   ! boundary element wall number (nbelem)
c        INTEGER, POINTER :: iwn(:)    ! boundary node wall number (nbnodes)

      ALLOCATE (mesh%bcn(mesh%nbnodes,3))
      mesh%bcn=0.0D0

      DO i=1,mesh%nbnodes
        dn=mesh%gbn(i)
        DO j=1,mesh%nbelem
          DO k=1,mesh%npob
            IF (mesh%ibc(j,k).EQ.dn.AND.mesh%bewn(j).EQ.mesh%iwn(i)) THEN
c              mesh%bcn(i,1)=anx(j)
c              mesh%bcn(i,2)=any(j)
c              mesh%bcn(i,3)=anz(j)
              mesh%bcn(i,1)=mesh%bcn(i,1)+anx(j)
              mesh%bcn(i,2)=mesh%bcn(i,2)+any(j)
              mesh%bcn(i,3)=mesh%bcn(i,3)+anz(j)
            END IF
          END DO
        END DO

C       Ker sem sesteval vektorje, na novo normiram
        AJAC=1.0D0/SQRT(mesh%bcn(i,1)**2+mesh%bcn(i,2)**2+mesh%bcn(i,3)**2)
        mesh%bcn(i,1)=mesh%bcn(i,1)*AJAC  ! povprecna smer, ne glede na velikost ploskev
        mesh%bcn(i,2)=mesh%bcn(i,2)*AJAC
        mesh%bcn(i,3)=mesh%bcn(i,3)*AJAC

      END DO

      DEALLOCATE (anx,any,anz)

      END

      
C -----------------------------------------------------------------------------
      SUBROUTINE GenMakroMesh(mesh)
C
C     $: Iz mreze naredi mrezo za makro
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      
      INTEGER i,j,k,l,maxq,nq,ii,jj
      REAL(8), ALLOCATABLE :: tmpq(:,:)
      

      mesh%nside=6 ! imamo heksaedre
      mesh%npof=4  ! imamo linearno 4 tockovno ploskev za flukse
      mesh%npofc=mesh%nside*mesh%npof
      
c     alociram ibf = conectivity za makro nezvezne flukse
      ALLOCATE(mesh%ibf(mesh%nicell,mesh%nside,mesh%npof))
      
c     Najprej ugotovimo koliko tock za q-je potrebujem
c     na vsaki ploskvi vsake celice imam 4 tocke za linearno interpolacijo q-jev.      
c     maksimalno setevilo q jev      
      maxq=mesh%nside*mesh%npof*mesh%nicell
      ALLOCATE (tmpq(maxq,3))
      nq=0
      
      DO i=1,mesh%nicell
        CALL AddSide(mesh,i,tmpq,maxq,nq,1,6,2,3,7) ! Spodnja ploskev
        CALL AddSide(mesh,i,tmpq,maxq,nq,2,6,7,8,5) ! spredaj ploskev
        CALL AddSide(mesh,i,tmpq,maxq,nq,3,7,3,4,8) ! desno ploskev
        CALL AddSide(mesh,i,tmpq,maxq,nq,4,3,2,1,4) ! zadaj ploskev
        CALL AddSide(mesh,i,tmpq,maxq,nq,5,2,6,5,1) ! levo ploskev
        CALL AddSide(mesh,i,tmpq,maxq,nq,6,5,8,4,1) ! zgoraj ploskev     
      END DO

c     skopiramo lokacije nodeov iz zacasnega v trajno polje
      mesh%nq=nq
      ALLOCATE (mesh%xq(mesh%nq,mesh%npx))
      DO i=1,mesh%nq
        DO j=1,mesh%npx
          mesh%xq(i,j)=tmpq(i,j)
        END DO
      END DO
      DEALLOCATE (tmpq)
      
c      ALLOCATE (mesh%qkode(mesh%nq))      
      ALLOCATE (mesh%ibcf(mesh%nbelem,mesh%npof))      

C
C     Sestavimo se ibcf = boundary element conectivity za flukse
C
      DO i=1,mesh%nbelem
        CALL AddSideBound(mesh,i)
      END DO
C
C     fluksi skozi stranico, ki si jo delita dva elementa sta samo
C     nasprotno predznacena. Naredim polje (nicell,nside) ko pove
C     kateri bo predznacen s plus kateri z minus
C      
      ALLOCATE (mesh%flpm(mesh%nicell,mesh%nside))
      mesh%flpm=0
      DO i=1,mesh%nicell
        DO j=1,mesh%nside
          IF (mesh%flpm(i,j).EQ.0) THEN ! tega se nimam
c           poiscem tisto celico ii in ploskev jj, ki se ujema s to
            CALL FindAdjCell(mesh,i,j,ii,jj)   
            mesh%flpm(i,j)=1
            IF (ii.NE.-1.AND.jj.NE.-1) THEN ! nisem na robu
              mesh%flpm(ii,jj)=-1
            END IF   
          END IF
        END DO
      END DO
C
C____ NUMBER OF CELLS TIMES NUMBER OF SOURCE POINTS IN A CELL
C          
      mesh%nsp=mesh%npoc+mesh%npofc   ! 27+24=51
      mesh%nicnsp=mesh%nicell*mesh%nsp   ! = 51*nicell=stevilo enacb   
      mesh%nbnpof=mesh%nbelem*mesh%npof  ! za polne matrike, seznam fluksov na robu


C
C     za kvadratne sisteme enacb
C
      ALLOCATE (mesh%sqUlist(mesh%nnodes,8)) ! 8 zato, ker imam lahko najvec osem sosedov
      ALLOCATE (mesh%sqUlistKM(mesh%nnodes,8)) ! 8 zato, ker imam lahko najvec osem sosedov
      ALLOCATE (mesh%sqUlistIC(mesh%nnodes,8)) ! 8 zato, ker imam lahko najvec osem sosedov
      ALLOCATE (mesh%sqUlistNO(mesh%nnodes)) ! stevilo enacb za izvorno tocko

      ALLOCATE (mesh%sqQlist(mesh%nq,2)) ! 2 zato, ker imam lahko najvec dva soseda
      ALLOCATE (mesh%sqQlistIC(mesh%nq,2)) ! 2 zato, ker imam lahko najvec osem sosedov
      ALLOCATE (mesh%sqQlistNO(mesh%nq)) ! stevilo enacb za izvorno tocko

      mesh%sqUlistNO=0
      DO j=1,mesh%nicell
        DO k=1,mesh%npoc
          i=mesh%idc(j,k)
          mesh%sqUlistNO(i)=mesh%sqUlistNO(i)+1
          mesh%sqUlist(i,mesh%sqUlistNO(i))=(j-1)*mesh%nsp+k ! vrstica v matriki integralov
          mesh%sqUlistKM(i,mesh%sqUlistNO(i))=(j-1)*mesh%npoc+k ! vrstica v matriki integralov za htx,hty,htz
          mesh%sqUlistIC(i,mesh%sqUlistNO(i))=j
        END DO
      END DO

      mesh%sqQlistNO=0
      DO j=1,mesh%nicell
        DO k=1,mesh%nside
          DO l=1,mesh%npof
            i=mesh%ibf(j,k,l)
            mesh%sqQlistNO(i)=mesh%sqQlistNO(i)+1
            mesh%sqQlist(i,mesh%sqQlistNO(i))=(j-1)*mesh%nsp+mesh%npoc+(k-1)*mesh%npof+l ! vrstica
            mesh%sqQlistIC(i,mesh%sqQlistNO(i))=j
          END DO
        END DO
      END DO
      
      END
C -----------------------------------------------------------------------------        
      SUBROUTINE FindAdjCell(mesh,c1,s1,c2,s2)
C
C     $: Iz celice c1 in strani s1 najde c2,s2, ki se ujemata
C
C -----------------------------------------------------------------------------
      USE inc_types 
      TYPE(meshType) :: mesh    
      INTEGER c1,s1,c2,s2
      INTEGER i,j,k,nistih
      c2=-1
      s2=-1 ! ce je c1,s1 na robu, potem c2,s2 ne bom nasel.
      
      DO i=1,mesh%nicell
        IF (i.NE.c1) THEN
          DO j=1,mesh%nside
            nistih=0
            DO k=1,mesh%npof
              IF (mesh%ibf(i,j,k).EQ.mesh%ibf(c1,s1,1).OR.mesh%ibf(i,j,k).EQ.mesh%ibf(c1,s1,2).OR.
     &            mesh%ibf(i,j,k).EQ.mesh%ibf(c1,s1,3).OR.mesh%ibf(i,j,k).EQ.mesh%ibf(c1,s1,4)) THEN
                nistih=nistih+1
              END IF           
            END DO
            IF (nistih.EQ.mesh%npof) THEN ! nasel sorodnika
              c2=i
              s2=j
              EXIT
            END IF
          END DO
        END IF
      END DO      
      
      END
C -----------------------------------------------------------------------------        
      SUBROUTINE AddSideBound(mesh,i)
C
C     $: Izracuna node za flux za eno boundary stranico 
C
C -----------------------------------------------------------------------------
      USE inc_types 
      TYPE(meshType) :: mesh     
      INTEGER i,n1,n2,n3,n4,jf,j     
      REAL(8) qx(4),qy(4),qz(4)
      LOGICAL imam      
      
      n1=1
      n2=3
      n3=5
      n4=7

c     izracunam lokacije tock za flux      
      CALL zv4nez4(mesh%x(mesh%ibc(i,n1),1),mesh%x(mesh%ibc(i,n1),2),mesh%x(mesh%ibc(i,n1),3),
     &               mesh%x(mesh%ibc(i,n2),1),mesh%x(mesh%ibc(i,n2),2),mesh%x(mesh%ibc(i,n2),3),
     &               mesh%x(mesh%ibc(i,n3),1),mesh%x(mesh%ibc(i,n3),2),mesh%x(mesh%ibc(i,n3),3),
     &               mesh%x(mesh%ibc(i,n4),1),mesh%x(mesh%ibc(i,n4),2),mesh%x(mesh%ibc(i,n4),3),
     &               qx(1),qy(1),qz(1),
     &               qx(2),qy(2),qz(2),
     &               qx(3),qy(3),qz(3),
     &               qx(4),qy(4),qz(4))   
     
c     preverim ce te tocke ze imam c seznamu
      DO j=1,4
        CALL ImamVsez(mesh%xq,mesh%nq,mesh%nq,qx(j),qy(j),qz(j),imam,jf)
        IF (imam.EQV..FALSE.) THEN  ! jure EQ EQV
          print *,"ERROR!!"
        ELSE 
          mesh%ibcf(i,j)=jf
        END IF
      END DO
             
      END

C -----------------------------------------------------------------------------        
      SUBROUTINE AddSide(mesh,i,tmpq,maxq,nq,ns,n1,n2,n3,n4)
C
C     $: Izracuna node za flux za eno stranico 
C
C -----------------------------------------------------------------------------
      USE inc_types 
      TYPE(meshType) :: mesh     
      INTEGER maxq,i,nq,n1,n2,n3,n4,ns,jf,j
      REAL(8) tmpq(maxq,3)       
      REAL(8) qx(4),qy(4),qz(4)
      LOGICAL imam      

c     izracunam lokacije tock za flux      
      CALL zv4nez4(mesh%x(mesh%idc(i,n1),1),mesh%x(mesh%idc(i,n1),2),mesh%x(mesh%idc(i,n1),3),
     &               mesh%x(mesh%idc(i,n2),1),mesh%x(mesh%idc(i,n2),2),mesh%x(mesh%idc(i,n2),3),
     &               mesh%x(mesh%idc(i,n3),1),mesh%x(mesh%idc(i,n3),2),mesh%x(mesh%idc(i,n3),3),
     &               mesh%x(mesh%idc(i,n4),1),mesh%x(mesh%idc(i,n4),2),mesh%x(mesh%idc(i,n4),3),
     &               qx(1),qy(1),qz(1),
     &               qx(2),qy(2),qz(2),
     &               qx(3),qy(3),qz(3),
     &               qx(4),qy(4),qz(4))   
     
c     preverim ce te tocke ze imam v seznamu
      DO j=1,4
        CALL ImamVsez(tmpq,maxq,nq,qx(j),qy(j),qz(j),imam,jf)
        IF (imam.EQV..FALSE.) THEN  ! jure EQ EQV
          nq=nq+1
          tmpq(nq,1)=qx(j)
          tmpq(nq,2)=qy(j)
          tmpq(nq,3)=qz(j)
          mesh%ibf(i,ns,j)=nq
        ELSE 
          mesh%ibf(i,ns,j)=jf
        END IF
      END DO
             
      END
     
C -----------------------------------------------------------------------------
      SUBROUTINE ImamVsez(tmpq,maxq,nq,qx,qy,qz,imam,zapst)
C
C     $: Ugotovi, ali v seznamu ze imam izbrano tocko
C
C -----------------------------------------------------------------------------
      INTEGER maxq,zapst,i,nq
      REAL(8) tmpq(maxq,3)
      REAL(8) qx,qy,qz
      LOGICAL imam
      
      imam=.FALSE.
      
      DO i=1,nq
        IF (ABS(tmpq(i,1)-qx).LT.1.0D-10.AND.
     &      ABS(tmpq(i,2)-qy).LT.1.0D-10.AND.
     &      ABS(tmpq(i,3)-qz).LT.1.0D-10) THEN
          imam=.TRUE.
          zapst=i
          EXIT
        END IF
      END DO
      
      END
     

