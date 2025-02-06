C ------------------------------------------------------------------------------------
      SUBROUTINE cTree_FillMatrices(io,inp,mName,fMat,nrow,ncol,hMat,nnz,n,cType)

C
C     Prebere drevesa in priprave H matrike
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(IOtype) io
      TYPE(inputtype) inp
c      TYPE(HtreeType) t1,t2
      TYPE(HmatrixType) :: hMat

      INTEGER nrow,ncol,i,ir,ic,r,c
      INTEGER cType ! compression type
      REAL(8) fMat(nrow,ncol)
      CHARACTER(8) mName
      REAL(8) MB,comp,norm,nnz,n,eps,nCOMe,nFULe,spar


      REAL(8), ALLOCATABLE :: vec(:),rez1(:),rez2(:)


      ALLOCATE ( vec(ncol),rez1(nrow),rez2(nrow) )
c     set up random vector for accuracy test
      DO i=1,ncol
        CALL RANDOM_NUMBER(vec(i))
      END DO
c     full matrix x vec results
      rez1=MATMUL(fMat,vec)

C     set compression type for all matrix parts (inadmissible are always full)
      DO i=1,hMat%Nparts
        IF (hMat%H(i)%admiss.EQ.1) THEN
          hMat%H(i)%cType = cType
        ELSE
          hMat%H(i)%cType = 1
        END IF
      END DO

      nCOMe = 0.0D0
C     Copy and/or compress matrix
      DO i=1,hMat%Nparts
C       first, no compression, just copy matrix to H format
        hMat%H(i)%fmat%nrow = hMat%H(i)%nRow  ! number of rows in matrix part
        hMat%H(i)%fmat%ncol = hMat%H(i)%nCol  ! number of cols in matrix part

        ALLOCATE (hMat%H(i)%fmat%v(hMat%H(i)%fmat%nrow,hMat%H(i)%fmat%ncol))
c       copy
        DO ir=1,hMat%H(i)%fmat%nrow
          r = hMat%H(i)%Rnodes(ir)
          DO ic=1,hMat%H(i)%fmat%ncol
            c = hMat%H(i)%Cnodes(ic)  ! LBN is preapplied
            hMat%H(i)%fmat%v(ir,ic) = fmat(r,c)
          END DO
        END DO

        IF (hMat%H(i)%cType.EQ.2) THEN ! wavelet compression
          CALL Che_WcomMatrix(hMat%H(i)%fmat%v,hMat%H(i)%fmat%nrow,hMat%H(i)%fmat%ncol,
     &                        hMat%H(i)%crsM,inp%WT_kappa,MB,comp,spar,norm)
          nCOMe=nCOMe+spar*DBLE(hMat%H(i)%fmat%nrow)*DBLE(hMat%H(i)%fmat%ncol)
          DEALLOCATE (hMat%H(i)%fmat%v)
        ELSE IF (hMat%H(i)%cType.EQ.3) THEN ! ACA compression
           CALL AcaDo(hMat%H(i)%fmat%v,hMat%H(i)%fmat%nrow,hMat%H(i)%fmat%ncol,hMat%H(i)%acaM,inp%ACA_eps,norm,inp%ACA_type)
           nCOMe = nCOMe + DBLE(hMat%H(i)%acaM%nrow)*DBLE(hMat%H(i)%acaM%rank)+DBLE(hMat%H(i)%acaM%rank)*DBLE(hMat%H(i)%acaM%ncol)
        ELSE ! no compression
          nCOMe = nCOMe + DBLE(hMat%H(i)%fmat%nrow)*DBLE(hMat%H(i)%fmat%ncol)
        END IF

      END DO

      nFULe = DBLE(nrow)*DBLE(ncol)

      comp = nCOMe / nFULe
      MB = nCOMe * 8.0 / 1024.0 / 1024.0

c     verify with random vector multiplication
      CALL cTree_MxV(hMat,vec,rez2)
c     get norm

      CALL getRMSnorm(rez1,rez2,nrow,norm)

C     output in log file
      WRITE(io%l,'(A8,3G20.10,I8)') mName,MB,comp,norm,hmat%nparts
C     total compression characteristics
      nnz=nnz+nCOMe
      n=n+nFULe


c      WRITE(*,'(A8,4G20.10)') mName,MB,comp,norm


      END


C ------------------------------------------------------------------------------------
      SUBROUTINE cTree_MxV(hMat,vec,rez)

C
C     Mnozenje H matrike z vektorjem
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HmatrixType) :: hMat

      INTEGER i,ir,ic,r,c
      REAL(8) vec(Hmat%ncol)
      REAL(8) rez(Hmat%nrow)

      REAL(8), ALLOCATABLE :: vv(:),rr(:)

C     set result to zero
      DO i=1,Hmat%nrow
        rez(i)=0.0D0
      END DO

C     multiply
      DO i=1,hMat%Nparts

c          DO ir=1,hMat%H(i)%fmat%nrow
c             r = hMat%H(i)%Rnodes(ir)
c            DO ic=1,hMat%H(i)%fmat%ncol
c               c = hMat%H(i)%Cnodes(ic) ! lbn is preapplied
c              rez(r) = rez(r) + hMat%H(i)%fmat%v(ir,ic) * vec(c)
c            END DO
c          END DO

C       allocate temp vectors for matmul
        ALLOCATE (vv(hMat%H(i)%fmat%ncol))
        ALLOCATE (rr(hMat%H(i)%fmat%nrow))
c       prepare column vektor
        DO ic=1,hMat%H(i)%fmat%ncol
          c = hMat%H(i)%Cnodes(ic) ! lbn is preapplied
          vv(ic) = vec(c)
        END DO

C       perform matrix time vector multiplication
        IF (hMat%H(i)%cType.EQ.1) THEN ! no compression
          rr= MATMUL(hMat%H(i)%fmat%v,vv)
        ELSE IF (hMat%H(i)%cType.EQ.2) THEN ! wavelet compression
          CALL Che_WTAWx(hMat%H(i)%crsM,vv,hMat%H(i)%fmat%ncol,rr,hMat%H(i)%fmat%nrow)
        ELSE ! ACA compression
          CALL aca_abxv(hMat%H(i)%acaM,vv,rr)
        END IF

C       copy solution back
        DO ir=1,hMat%H(i)%fmat%nrow
          r = hMat%H(i)%Rnodes(ir)
          rez(r)=rez(r)+ rr(ir)
c          print *,i,ir,rr(ir)
        END DO
C       deallocate temp vectors
        DEALLOCATE (vv,rr)

      END DO


      END

C ------------------------------------------------------------------------------------
      SUBROUTINE SetUpHmatrices(env,io,inp,mesh,
     &               hUmatB,hUmatBp,hUmatAbdx,hUmatAbdy,hUmatAbdz,hUmatDx,hUmatDy,hUmatDz,
     &               hQmatB,hQmatBp,hQmatAbdx,hQmatAbdy,hQmatAbdz,hQmatDx,hQmatDy,hQmatDz)

C
C     Prebere drevesa in priprave H matrike
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(IOtype) io
      TYPE(penv) :: env
      TYPE(inputtype) inp
      TYPE(HtreeType) bt,dt,bdt,ddt,bbt,dbt
      TYPE(meshType) :: mesh

      TYPE(HmatrixType) :: hUmatB,hUmatBp,hUmatAbdX,hUmatAbdY,hUmatAbdZ,hUmatDx,hUmatDy,hUmatDz
      TYPE(HmatrixType) :: hQmatB,hQmatBp,hQmatAbdX,hQmatAbdY,hQmatAbdZ,hQmatDx,hQmatDy,hQmatDz


      INTEGER n


C     Read file and prepare booundary and domain trees
      CALL ReadcTreeFile(env,io,inp,mesh,bt,dt)

C     Create tree combinations
C     Boundary - Domain tree
      CALL cTree_CombineTrees(mesh,bdt,bt,dt)
C     Boundary - Boundary tree
      CALL cTree_CombineTrees(mesh,bbt,bt,bt)
C     Domain - Domain tree
      CALL cTree_CombineTrees(mesh,ddt,dt,dt)
C     Domain - Boundary tree
      CALL cTree_CombineTrees(mesh,dbt,dt,bt)


C     Verify admisibility

      CALL cTree_Admis(bdt,bt,dt,1,inp%ct_eta)
      CALL cTree_Admis(ddt,dt,dt,1,inp%ct_eta)
      CALL cTree_Admis(bbt,bt,bt,1,inp%ct_eta)
      CALL cTree_Admis(dbt,dt,bt,1,inp%ct_eta)

C     Draw matrix
      CALL cTree_drawMatrix(mesh,dbt,dt,bt,"BD")
      CALL cTree_drawMatrix(mesh,ddt,dt,dt,"DD")
      CALL cTree_drawMatrix(mesh,bbt,bt,bt,"BB")
      CALL cTree_drawMatrix(mesh,dbt,dt,bt,"DB")


C     Create H matrices
      CALL cTree_CreateHmatrix(mesh,hUmatB,ddt,dt,dt,mesh%nnodes,mesh%nnodes,1,1,0)
      CALL cTree_CreateHmatrix(mesh,hUmatBp,ddt,dt,dt,mesh%nnodes,mesh%nnodes,1,1,0)

      CALL cTree_CreateHmatrix(mesh,hUmatAbdx,dbt,dt,bt,mesh%nnodes,mesh%nbnodes,1,1,1)
      CALL cTree_CreateHmatrix(mesh,hUmatAbdy,dbt,dt,bt,mesh%nnodes,mesh%nbnodes,1,1,1)
      CALL cTree_CreateHmatrix(mesh,hUmatAbdz,dbt,dt,bt,mesh%nnodes,mesh%nbnodes,1,1,1)

      CALL cTree_CreateHmatrix(mesh,hUmatDx,ddt,dt,dt,mesh%nnodes,mesh%nnodes,1,1,0)
      CALL cTree_CreateHmatrix(mesh,hUmatDy,ddt,dt,dt,mesh%nnodes,mesh%nnodes,1,1,0)
      CALL cTree_CreateHmatrix(mesh,hUmatDz,ddt,dt,dt,mesh%nnodes,mesh%nnodes,1,1,0)

      CALL cTree_CreateHmatrix(mesh,hQmatB,bdt,bt,dt,mesh%nbnpof,mesh%nnodes,2,1,0)
      CALL cTree_CreateHmatrix(mesh,hQmatBp,bdt,bt,dt,mesh%nbnpof,mesh%nnodes,2,1,0)

      CALL cTree_CreateHmatrix(mesh,hQmatAbdx,bbt,bt,bt,mesh%nbnpof,mesh%nbnodes,2,1,1)
      CALL cTree_CreateHmatrix(mesh,hQmatAbdy,bbt,bt,bt,mesh%nbnpof,mesh%nbnodes,2,1,1)
      CALL cTree_CreateHmatrix(mesh,hQmatAbdz,bbt,bt,bt,mesh%nbnpof,mesh%nbnodes,2,1,1)

      CALL cTree_CreateHmatrix(mesh,hQmatDx,bdt,bt,dt,mesh%nbnpof,mesh%nnodes,2,1,0)
      CALL cTree_CreateHmatrix(mesh,hQmatDy,bdt,bt,dt,mesh%nbnpof,mesh%nnodes,2,1,0)
      CALL cTree_CreateHmatrix(mesh,hQmatDz,bdt,bt,dt,mesh%nbnpof,mesh%nnodes,2,1,0)

      END

C ------------------------------------------------------------------------------------
      SUBROUTINE cTree_CreateHmatrix(mesh,Hmat,t,t1,t2,nr,nc,t1UQ,t2UQ,iLBN)
C
C     Izdela H matriko iz dreves
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t,t1,t2
      TYPE(HmatrixType) Hmat
      TYPE(meshType) :: mesh
      INTEGER n,nr,nc
      INTEGER t1UQ,t2UQ ! U or Q nodes (1=U, 2=Q)
      INTEGER iLBN ! =1 appy LBN

C     REmeber size of the whole matrix
      hMat%nRow = nr
      hMat%nCol = nc

C     count number of matrix parts
      n=0
      CALL cTree_CountMatrixParts(t,t1,t2,1,n)

C     allocate matrix parts
      Hmat%Nparts = n
      ALLOCATE(Hmat%H(Hmat%Nparts))

C     set up matrix parts
      n=0
      CALL cTree_SetUpMatrix(mesh,Hmat,t,t1,t2,1,n,t1UQ,t2UQ,iLBN)

      END



C ------------------------------------------------------------------------------------
      SUBROUTINE ReadcTreeFile(env,io,inp,mesh,bt,dt)
C
C     Prebere drevesa in seznam na disk
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(penv) :: env
      TYPE(IOtype) io
      TYPE(inputtype) inp
      TYPE(HtreeType) bt,dt
      TYPE(meshType) :: mesh

      INTEGER nme    ! stevilo vseh vej
      INTEGER ctfv,i,j

      OPEN(UNIT=io%ctree,FILE=io%ctree_name,ERR=100,FORM='UNFORMATTED',STATUS="OLD")

c     version
      READ (io%ctree) ctfv

      IF (inp%cTreeFileVers.EQ.ctfv) THEN
c       array sizes
        READ (io%ctree) dt%nvej
        READ (io%ctree) bt%nvej
        READ (io%ctree) nme
      ELSE
        GOTO 101
      END IF

C     Read domain tree
      ALLOCATE ( dt%veja(dt%nvej) )
      CALL ReadSingleTree(dt,io%ctree)
C     Read boundary tree
      ALLOCATE ( bt%veja(bt%nvej) )
      CALL ReadSingleTree(bt,io%ctree)

      CLOSE (io%ctree)


      WRITE(io%l,'(A)') "Cluster trees"
      WRITE(io%l,'(A,I8)') "dNvej      =",bt%Nvej
      WRITE(io%l,'(A,I8)') "bNvej      =",dt%Nvej


C     Create node lists
      CALL cTree_NodeList(mesh%idc,mesh%npoc,mesh%nicell,mesh%nnodes,dt,1)
      CALL cTree_NodeList(mesh%ibc,mesh%npob,mesh%nbelem,mesh%nnodes,bt,1)
C     flux nodes
      CALL cTree_fNodeList(mesh%ibcf,mesh%npof,mesh%nbelem,bt)
C     get depth of trees
      CALL cTree_GetMaxLev(dt)
      CALL cTree_GetMaxLev(bt)
C     Estimate sizes of tree branches
      CALL cTree_EstClusterSize(mesh,bt)
      CALL cTree_EstClusterSize(mesh,dt)

      RETURN

100   CALL WarnErr(env,io,inp,2,"ReadcTreeSizesDisk","Can't open .ctree file",0)
101   CALL WarnErr(env,io,inp,2,"ReadcTreeSizesDisk","Wrong cTreeFileVers",0)

      END





C ------------------------------------------------------------------------------------
      SUBROUTINE cTree_drawMatrix(mesh,t,t1,t2,text)
C
C     $: zapise strukturo matrike po elementih na disk
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t,t1,t2
      TYPE(meshType) :: mesh

      INTEGER lun,i
      REAL(8) xtl,ytl
      CHARACTER(255) filename
      CHARACTER(2) text

      lun=99
      xtl=0.0D0
      ytl=0.0D0

c     Izris strukture matrike v tecplot
      WRITE (filename,'(A,A,A2,A)') trim(mesh%filename),".BlockStruc.",text,".lay"
      OPEN(UNIT=lun,FILE=filename,STATUS="UNKNOWN")
      WRITE (lun,'(A)') "#!MC 1100"
      WRITE (lun,'(A)') "$!FRAMELAYOUT"
      WRITE (lun,'(A)') "  SHOWBORDER = NO"

      CALL cTree_RecDraw(t,t1,t2,1,lun,xtl,ytl)

      CLOSE (lun)

      END


C ------------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE cTree_RecDraw(t,t1,t2,veja,lun,xtl,ytl)
C
C     $: bruhh
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t,t1,t2
      INTEGER veja,i,t1v,t2v

      INTEGER lun,c1,c2,me
      REAL(8) xtl,ytl,lxtl,lytl


      t1v = t%veja(veja)%t1v
      t2v = t%veja(veja)%t2v

C     ali sem admissible ali pa list potem risem
      IF ( t%veja(veja)%admiss .EQ. 1 . OR. t%veja(veja)%nchild .EQ. 0 ) THEN

        WRITE (lun,'(A,F15.5,A,F15.5,A)') "$!ATTACHGEOM GEOMTYPE = RECTANGLE ANCHORPOS { X = ",xtl," Y = ",ytl," }"

C       inadmissible so RDECI = RDECE ne kompresiram
        IF (t%veja(veja)%admiss.EQ.0) THEN
          WRITE (lun,'(A)') "ISFILLED = YES"
          WRITE (lun,'(A)') "FILLCOLOR = RED"
        END IF
        WRITE (lun,'(A)') "RAWDATA"
        WRITE (lun,'(2F15.5)') DBLE(t2%veja(t2v)%nelems),-DBLE(t1%veja(t1v)%nelems)


      ELSE  ! moram iti globje
        DO i=1,t%veja(veja)%nchild
          lxtl = xtl
          lytl = ytl

          IF (i.EQ.2) THEN
              c1 = t1%veja( t%veja(t%veja(veja)%child(1))%t1v )%nelems
              c2 = t1%veja( t%veja(t%veja(veja)%child(2))%t1v )%nelems
              me = t1%veja( t%veja(veja)%t1v )%nelems

              IF (c1+c2.EQ.me) THEN
                lytl = ytl - c1
              ELSE
                c1 = t2%veja( t%veja(t%veja(veja)%child(1))%t2v )%nelems
                lxtl = xtl + c1
              END IF

          END IF
          CALL cTree_RecDraw(t,t1,t2,t%veja(veja)%child(i),lun,lxtl,lytl)
        END DO
      END IF

      END



C ------------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE cTree_SetUpMatrix(mesh,Hmat,t,t1,t2,veja,n,t1UQ,t2UQ,iLBN)
C
C     $: bruhh
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t,t1,t2
      TYPE(HmatrixType) Hmat
      TYPE(meshType) :: mesh
      INTEGER veja,i,n,t1UQ,t2UQ,iLBN


C     ali sem admissible ali pa list
      IF ( t%veja(veja)%admiss .EQ. 1 . OR. t%veja(veja)%nchild .EQ. 0 ) THEN

        n = n + 1

        Hmat%H(n)%t1v = t%veja(veja)%t1v
        Hmat%H(n)%t2v = t%veja(veja)%t2v

        IF (t1UQ.EQ.1) THEN ! U nodes
          Hmat%H(n)%nRow = t1%veja(Hmat%H(n)%t1v)%nnodes
          ALLOCATE (Hmat%H(n)%Rnodes(Hmat%H(n)%nRow))
          DO i=1,Hmat%H(n)%nRow
            Hmat%H(n)%Rnodes(i)=t1%veja(Hmat%H(n)%t1v)%nodes(i)
          END DO
        ELSE  ! Q nodes
          Hmat%H(n)%nRow = t1%veja(Hmat%H(n)%t1v)%nq
          ALLOCATE (Hmat%H(n)%Rnodes(Hmat%H(n)%nRow))
          DO i=1,Hmat%H(n)%nRow
            Hmat%H(n)%Rnodes(i)=t1%veja(Hmat%H(n)%t1v)%qnodes(i)
          END DO
        END IF

        IF (t2UQ.EQ.1) THEN ! U nodes
          Hmat%H(n)%nCol = t2%veja(Hmat%H(n)%t2v)%nnodes
          ALLOCATE (Hmat%H(n)%Cnodes(Hmat%H(n)%nCol))
          DO i=1,Hmat%H(n)%nCol
            IF (iLBN.EQ.0) THEN
              Hmat%H(n)%Cnodes(i)=t2%veja(Hmat%H(n)%t2v)%nodes(i)
            ELSE
              Hmat%H(n)%Cnodes(i)=mesh%lbn(t2%veja(Hmat%H(n)%t2v)%nodes(i))
            END IF
          END DO
        ELSE  ! Q nodes
          Hmat%H(n)%nCol = t2%veja(Hmat%H(n)%t2v)%nq
          ALLOCATE (Hmat%H(n)%Cnodes(Hmat%H(n)%nCol))
          DO i=1,Hmat%H(n)%nCol
            Hmat%H(n)%Cnodes(i)=t2%veja(Hmat%H(n)%t2v)%qnodes(i)
          END DO
        END IF

        Hmat%H(n)%admiss = t%veja(veja)%admiss

      ELSE  ! moram iti globje
        DO i=1,t%veja(veja)%nchild
          CALL cTree_SetUpMatrix(mesh,Hmat,t,t1,t2,t%veja(veja)%child(i),n,t1UQ,t2UQ,iLBN)
        END DO
      END IF

      END


C ------------------------------------------------------------------------------------
      SUBROUTINE cTree_CombineTrees(mesh,t,t1,t2)
C
C     combine trees
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t,t1,t2
      TYPE(meshType) :: mesh
      INTEGER tv

C     Predpostavljam, da sta t1 in t2 enako globoka !!
      t%maxLev = 2 * t1%maxLev - 1

C     Count branches in combined tree
      t%nVej = 0
      CALL cTree_CountCreateTreeCombo(t,t1,t2,1,t%nVej)

C     Allocate memory for branches
      ALLOCATE (t%veja(t%nVej))

C     Set up branch combination (lev in lev+1)
      tv=0
      CALL cTree_CreateTreeCombo(t,t1,t2,1,tv)

C     Set up parents
      CALL cTree_SetUpParents(t,t1,t2)
C     Set up children
      CALL cTree_SetUpChildren(t)
C     Set up admissible and inadmissible blocks
      CALL cTree_GetDistBetCl(mesh,t,t1,t2,1)

      END

C ------------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE cTree_CountMatrixParts(t,t1,t2,veja,n)
C
C     $: bruhh
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t,t1,t2
      INTEGER veja,i,t1v,t2v,n


C     ali sem admissible ali pa list
      IF ( t%veja(veja)%admiss .EQ. 1 . OR. t%veja(veja)%nchild .EQ. 0 ) THEN

        t1v = t%veja(veja)%t1v
        t2v = t%veja(veja)%t2v
        n = n + 1

      ELSE  ! moram iti globje
        DO i=1,t%veja(veja)%nchild
          CALL cTree_CountMatrixParts(t,t1,t2,t%veja(veja)%child(i),n)
        END DO
      END IF

      END



C ------------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE cTree_Admis(t,t1,t2,veja,eta)
C
C     $: Naredi addmisibility kriterij
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t,t1,t2
      INTEGER veja,i
      REAL(8) d1,d2,d,eta

      d1 = t1%veja(t%veja(veja)%t1v)%diameter
      d2 = t2%veja(t%veja(veja)%t2v)%diameter
      d  = t%veja(veja)%dist

      IF ( MIN(d1,d2) .LT. eta * d ) THEN
        t%veja(veja)%admiss = 1
      ELSE
        t%veja(veja)%admiss = 0
      END IF

c      print *,veja,t%veja(veja)%admiss,d1,d2,d

      DO i=1,t%veja(veja)%nchild
        CALL cTree_Admis(t,t1,t2,t%veja(veja)%child(i),eta)
      END DO

      END



C ------------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE cTree_GetDistBetCl(mesh,t,t1,t2,veja)
C
C     $: bruhh
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t,t1,t2
      TYPE(meshType) :: mesh
      INTEGER veja,i,j,t1v,t2v,n1,n2

      REAL(8) d1,d2,d

      d1 = t1%veja(t%veja(veja)%t1v)%diameter
      d2 = t2%veja(t%veja(veja)%t2v)%diameter

      t1v = t%veja(veja)%t1v
      t2v = t%veja(veja)%t2v

      t%veja(veja)%dist=1.0E10

c     Calculate distance between clusters
      DO i=1,t1%veja(t1v)%nnodes
        n1 = t1%veja(t1v)%nodes(i)
        DO j=1,t2%veja(t2v)%nnodes
          n2 = t2%veja(t2v)%nodes(j)

            d = SQRT( (mesh%x(n1,1)-mesh%x(n2,1))**2
     &               +(mesh%x(n1,2)-mesh%x(n2,2))**2
     &               +(mesh%x(n1,3)-mesh%x(n2,3))**2)
            t%veja(veja)%dist = MIN(d,t%veja(veja)%dist)
        END DO
      END DO


c      print *,veja,d1,d2,d

      DO i=1,t%veja(veja)%nchild
        CALL cTree_GetDistBetCl(mesh,t,t1,t2,t%veja(veja)%child(i))
      END DO

      END


C ------------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE cTree_print(t,veja,tv)
C
C     $: bruhh
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t
      INTEGER tv,veja,i

      tv=tv+1

      DO i=1,t%veja(veja)%nchild
        CALL cTree_print(t,t%veja(veja)%child(i),tv)
      END DO

      END


C ------------------------------------------------------------------------------------
      SUBROUTINE cTree_SetUpChildren(t)
C
C     na podlagi znanih starsev pogrunta otroke
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE


      TYPE(HtreeType) t
      INTEGER veja,p

c     Set number of children to zero
      DO veja = 1,t%nVej
        t%veja(veja)%nchild = 0
      END DO

      DO veja = 2,t%nVej ! prva veja nima starsev
        p = t%veja(veja)%parent
        t%veja(p)%nchild = t%veja(p)%nchild + 1
        t%veja(p)%child(t%veja(p)%nchild)=veja
      END DO

      END



C ------------------------------------------------------------------------------------
      SUBROUTINE cTree_SetUpParents(t,t1,t2)
C
C     pogrunta starse vej v drevesu
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t,t1,t2
      INTEGER veja,j

c      print *,t%maxLev,t%nVej

      DO veja = 1,t%nVej
c        if (veja<5) print *,"veja",veja
        t%veja(veja)%parent=0
        DO j = 1,t%nVej
          IF ( t%veja(veja)%lev .EQ. t%veja(j)%lev + 1 ) THEN
c            if (veja<5) print *,"j",t%veja(veja)%t1v,t%veja(j)%t1v
c            if (veja<5) print *,"j",t%veja(veja)%t2v,t%veja(j)%t1v
c            if (veja<5) print *,"j",t%veja(j)%t1v,t%veja(j)%t2v,t%veja(veja)%t1v,t%veja(veja)%t2v
c            IF ( ( t%veja(veja)%t1v .EQ.  t%veja(j)%t1v ) .OR. ( t%veja(veja)%t2v .EQ.  t%veja(j)%t2v ) )   THEN
c              t%veja(veja)%parent=j
c              if (veja<6) print *,"parent",t%veja(veja)%parent
c            END IF
c            if (veja==5) print *,t%veja(veja)%t1v,t%veja(j)%t1v,"ENA"
c           if (veja==5) print *,t2%veja(t%veja(j)%t2v)%child(1),t%veja(veja)%t2v
c           if (veja==5) print *,t2%veja(t%veja(j)%t2v)%child(2),t%veja(veja)%t2v
            IF ( ( t%veja(veja)%t1v .EQ.  t%veja(j)%t1v ) . AND . (
     &         ( t2%veja(t%veja(j)%t2v)%child(1) .EQ. t%veja(veja)%t2v ) . OR .
     &         ( t2%veja(t%veja(j)%t2v)%child(2) .EQ. t%veja(veja)%t2v ) ) ) THEN
                t%veja(veja)%parent=j
c                if (veja<6) print *,"parent",t%veja(veja)%parent
            END IF

            IF ( ( t%veja(veja)%t2v .EQ.  t%veja(j)%t2v ) . AND . (
     &         ( t1%veja(t%veja(j)%t1v)%child(1) .EQ. t%veja(veja)%t1v ) . OR .
     &         ( t1%veja(t%veja(j)%t1v)%child(2) .EQ. t%veja(veja)%t1v ) ) ) THEN
                t%veja(veja)%parent=j
c                if (veja<6) print *,"parent",t%veja(veja)%parent
            END IF



          END IF
        END DO
c        if (veja<100) print *,veja,t%veja(veja)%parent,t%veja(veja)%t1v,t%veja(veja)%t2v
c        print *,t%veja(veja)%parent
      END DO

      END




C ------------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE cTree_CreateTreeCombo(t,t1,t2,t1v,tv)
C
C     $: bruhh
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t,t1,t2
      INTEGER i,t1v,t2v,tv,mojLvl

C     poisci vse v t2 na istem levelu
      mojLvl =  t1%veja(t1v)%lev
      DO t2v=1,t2%Nvej
        IF (t2%veja(t2v)%lev.EQ.mojLvl) THEN
         tv=tv+1
         t%veja(tv)%lev= 2 * mojLvl - 1
         t%veja(tv)%t1v=t1v
         t%veja(tv)%t2v=t2v
c         print *,tv,t1v,t2v
        END IF
        IF (t2%veja(t2v)%lev.EQ.mojLvl+1) THEN
         tv=tv+1
         t%veja(tv)%lev= 2 * mojLvl
         t%veja(tv)%t1v=t1v
         t%veja(tv)%t2v=t2v
c         print *,tv,t1v,t2v
        END IF
      END DO

      DO i=1,t1%veja(t1v)%nchild
        CALL cTree_CreateTreeCombo(t,t1,t2,t1%veja(t1v)%child(i),tv)
      END DO

      END

C ------------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE cTree_CountCreateTreeCombo(t,t1,t2,t1v,tv)
C
C     $: bruhh
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t,t1,t2
      INTEGER i,t1v,t2v,tv,mojLvl

C     poisci vse v t2 na istem levelu
      mojLvl =  t1%veja(t1v)%lev
      DO t2v=1,t2%Nvej
        IF (t2%veja(t2v)%lev.EQ.mojLvl) THEN
         tv=tv+1
        END IF
        IF (t2%veja(t2v)%lev.EQ.mojLvl+1) THEN
         tv=tv+1
        END IF
      END DO

      DO i=1,t1%veja(t1v)%nchild
        CALL cTree_CountCreateTreeCombo(t,t1,t2,t1%veja(t1v)%child(i),tv)
      END DO

      END




C ------------------------------------------------------------------------------------
      SUBROUTINE cTree_GetMaxLev(t)
C
C     globina drevesa
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t
      INTEGER veja

      t%maxLev=0
      DO veja=1,t%Nvej
        t%maxLev = MAX(t%maxLev,t%veja(veja)%lev)
      END DO

      END



C ------------------------------------------------------------------------------------
      SUBROUTINE cTree_EstClusterSize(mesh,t)
C
C     za velikost clustra vzame največjo razdaljo med node v clustru
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t
      TYPE(meshType) :: mesh
      INTEGER veja,i,j,n1,n2
      REAL(8) d

C     zanka po vejah drevesa
      DO veja = 1,t%nVej
        t%veja(veja)%diameter = 0.0D0
        DO i=1,t%veja(veja)%nnodes
          n1 = t%veja(veja)%nodes(i)
          DO j=1,t%veja(veja)%nnodes
            n2 = t%veja(veja)%nodes(j)
            d = SQRT( (mesh%x(n1,1)-mesh%x(n2,1))**2
     &               +(mesh%x(n1,2)-mesh%x(n2,2))**2
     &               +(mesh%x(n1,3)-mesh%x(n2,3))**2)
            t%veja(veja)%diameter = MAX(d,t%veja(veja)%diameter)
          END DO
        END DO
c        print *,veja,t%veja(veja)%diameter
      END DO

      END

C ------------------------------------------------------------------------------------
      SUBROUTINE cTree_fNodeList(ibcf,npof,nbelem,t)
C
C     $: Izdela sezname vozlisc, tako da se ne podvajajo med obmocji
C     pri fluksih ni problema, ker se ne podvajajo med elementi
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t
      INTEGER npof,nbelem,nnodes
      INTEGER ibcf(nbelem,npof)

      INTEGER i,j,k,ele,veja,l,inode


C     zanka po vejah drevesa
      DO veja = 1,t%nVej
C       za vsako vejo vem koliko nodeov bo
        t%veja(veja)%nq = t%veja(veja)%nelems * npof
        ALLOCATE ( t%veja(veja)%qnodes(t%veja(veja)%nq) )
C       naredim seznam vseh vozlisc na elementih
        l=0
        DO j=1,t%veja(veja)%nelems
          ele=t%veja(veja)%elems(j)
          DO k=1,npof
            l=l+1
            t%veja(veja)%qnodes(l)=(ele-1)*npof+k !ibcf(ele,k)
          END DO
        END DO

      END DO



      END



C ------------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE cTree_NodeList(idc,npoc,nelem,nnodes,t,veja)
C
C     $: Izdela sezname vozlisc, tako da se ne podvajajo med obmocji
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t
      INTEGER npoc,nnodes,nelem
      INTEGER idc(nelem,npoc)

      INTEGER i,j,k,ele,veja,l,inode

      INTEGER, ALLOCATABLE :: inodes(:)
      INTEGER, ALLOCATABLE :: veknn(:)

      LOGICAL imastars,imaistilev

c     za vsako tocko iz mojega robnega elementa preverim:
c        - ali jo ima moj stars
c        - ali jo ima katera od vej na istem levelu

      ALLOCATE (veknn(nnodes))
      ALLOCATE (inodes(nnodes))

c     naredim seznam vseh robnih vozlisc v tej veji
      veknn=0
      DO j=1,t%veja(veja)%nelems
        ele=t%veja(veja)%elems(j)
        DO k=1,npoc
          veknn(idc(ele,k))=1
        END DO
      END DO

c     odlocim se katera so v tej veji in katera v drugih
      inodes=0
      DO inode=1,nnodes
        IF (veknn(inode).EQ.1) THEN
c          ali ima inode moj stars
          IF (veja.EQ.1) THEN
            imastars=.TRUE.
          ELSE
            imastars=.FALSE.
            DO k=1,t%veja(t%veja(veja)%parent)%nnodes
              IF (inode.EQ.t%veja(t%veja(veja)%parent)%nodes(k)) THEN
                imastars=.TRUE.
              END IF
            END DO
          END IF
          IF (imastars.EQV..TRUE.) THEN
c           ali ima nekdo na istem levelu
            imaistilev=.FALSE.
            DO l=1,t%nvej
              IF (l.NE.veja.AND.t%veja(l)%lev.EQ.t%veja(veja)%lev) THEN
                DO k=1,t%veja(l)%nnodes
                  IF (inode.EQ.t%veja(l)%nodes(k)) THEN
                    imaistilev=.TRUE.
                    EXIT
                  END IF
                END DO
              END IF
            END DO
            IF (imaistilev.EQV..FALSE.) THEN ! ta je moj !!!
              t%veja(veja)%nnodes=t%veja(veja)%nnodes+1
              inodes(t%veja(veja)%nnodes)=inode
            END IF
          END IF
        END IF
      END DO

c     naredim seznam tistih, ki so v moji veji
      ALLOCATE (t%veja(veja)%nodes(t%veja(veja)%nnodes))
      DO k=1,t%veja(veja)%nnodes
        t%veja(veja)%nodes(k)=inodes(k)
      END DO

c      print *,veja,t%veja(veja)%lev,t%veja(veja)%nnodes

      DEALLOCATE(inodes,veknn)

      DO i=1,t%veja(veja)%nchild
        CALL cTree_NodeList(idc,npoc,nelem,nnodes,t,t%veja(veja)%child(i))
      END DO

      END


C ------------------------------------------------------------------------------------
      SUBROUTINE ReadSingleTree(t,lun)
C
C     Prebere drevesa in seznam na disk
C
C ------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(HtreeType) t
      INTEGER lun,i

      DO i=1,t%nvej
        READ(lun) t%veja(i)%nelems,t%veja(i)%lev,t%veja(i)%parent,t%veja(i)%nchild,
     &            t%veja(i)%child(1),t%veja(i)%child(2),t%veja(i)%nnodes
        ALLOCATE (t%veja(i)%elems(t%veja(i)%nelems))
        ALLOCATE (t%veja(i)%nodes(t%veja(i)%nnodes))
        CALL rdivec(lun,t%veja(i)%nelems,t%veja(i)%elems)
        CALL rdivec(lun,t%veja(i)%nnodes,t%veja(i)%nodes)
        t%veja(i)%nnodes=0 ! to moram narediti se enkrat, ker so zraven robni, ki jih ne rabim
        DEALLOCATE (t%veja(i)%nodes)
      END DO

      END
