C     ******************************************************************
C     *                                                                *
C     *                          lib_ACA.f                              *
C     *                                                                *
C     ******************************************************************
C     *                                                                *
C     *                 Authors: Jure Ravnik, Jan Tibaut               *
C     *                                                                *
C     ******************************************************************
C     *                                                                *
C     *                     Version:  20.2.2018                        *
C     *                                                                *
C     ******************************************************************
C
C     needs module with
C
C      TYPE acaMatrixType
C        INTEGER nrow,ncol,rank  ! stevilo vrstic, stoplcev, rank
C        REAL(8), POINTER :: a(:,:),b(:,:) ! a(nrow,rank) b(rank,ncol)
C      END TYPE


C     ------------------------------------------------------------------
      SUBROUTINE aca_RankCrossAprox(FullMat,nrow,ncol,aM,eps,norm)
C
C     Izdela a,b = low rank approximation of fullmat
C
C     ------------------------------------------------------------------
      USE inc_types

      TYPE(acaMatrixType) :: acaMat,aM
      REAL(8) eps,norm,s,s2

      INTEGER nrow,ncol,rank
      REAL(8) FullMat(nrow,ncol)

      INTEGER lrow,lcol,irank
      REAL(8) delta
      REAL(8), ALLOCATABLE :: vec(:),rez1(:),rez2(:)

      LOGICAL done

      done = .FALSE.
      rank = MIN(nrow,ncol)

      ALLOCATE ( vec(ncol),rez1(nrow),rez2(nrow) )
c     set up random vector for accuracy test
      DO i=1,ncol
        CALL RANDOM_NUMBER(vec(i))
      END DO
c     full matrix x vec results
      rez1=MATMUL(FullMat,vec)


      acaMat%nrow=nrow
      acaMat%ncol=ncol
      acaMat%rank=rank
      ALLOCATE (acaMat%a(acaMat%nrow,acaMat%rank))
      ALLOCATE (acaMat%b(acaMat%rank,acaMat%ncol))
      acaMat%a=0.0D0
      acaMat%b=0.0D0


      DO irank=1,rank
        CALL FindLargestValue(FullMat,nrow,ncol,lrow,lcol,delta)
c        IF (delta.EQ.0.0D0) Print *,"nula"
        CALL CopyCR2ab(FullMat,nrow,ncol,lrow,lcol,delta,acaMat%a,acaMat%b,rank,irank)
        CALL SubsApp(FullMat,nrow,ncol,acaMat%a,acaMat%b,rank,irank)

c       verify with random vector multiplication
        CALL aca_abxv(acaMat,vec,rez2)
c       get norm
        CALL getRMSnorm(rez1,rez2,nrow,norm)

c       if norm ok, exit
        IF (norm.LT.eps) THEN
            aM%nrow=nrow
            aM%ncol=ncol
            aM%rank=irank
            ALLOCATE (aM%a(aM%nrow,aM%rank))
            ALLOCATE (aM%b(aM%rank,aM%ncol))
            DO j=1,aM%rank
                DO i=1,aM%nrow
                    aM%a(i,j)=acaMat%a(i,j)
                END DO
            END DO
            DO j=1,aM%ncol
                DO i=1,aM%rank
                    aM%b(i,j)=acaMat%b(i,j)
                END DO
            END DO
            DEALLOCATE (acaMat%a,acaMat%b)
            done = .TRUE.
           EXIT
        END IF
      END DO

      IF (done.EQV..FALSE.) THEN ! something went wrong, copy everything
            print *,"EE",rank
            aM%nrow=nrow
            aM%ncol=ncol
            aM%rank=rank
            ALLOCATE (aM%a(aM%nrow,aM%rank))
            ALLOCATE (aM%b(aM%rank,aM%ncol))
            DO j=1,aM%rank
                DO i=1,aM%nrow
                    aM%a(i,j)=acaMat%a(i,j)
c                    print *,acaMat%a(i,j)
                END DO
            END DO
            DO j=1,aM%ncol
                DO i=1,aM%rank
                    aM%b(i,j)=acaMat%b(i,j)
                END DO
            END DO
            DEALLOCATE (acaMat%a,acaMat%b)
      END IF


      DEALLOCATE (vec,rez1,rez2)

      END


C     ------------------------------------------------------------------
      SUBROUTINE aca_CrossAprox(FullMat,nrow,ncol,acaMat,eps,norm)
C
C     Izdela a,b = low rank approximation of fullmat
C
C     ------------------------------------------------------------------
      USE inc_types

      TYPE(acaMatrixType) :: acaMat
      REAL(8) eps,norm,s,s2

      INTEGER nrow,ncol,rank
      REAL(8) FullMat(nrow,ncol)

      INTEGER lrow,lcol,irank
      REAL(8) delta
      REAL(8), ALLOCATABLE :: vec(:),rez1(:),rez2(:)

      rank = INT(eps*DBLE(MIN(nrow,ncol)))



      ALLOCATE ( vec(ncol),rez1(nrow),rez2(nrow) )
c     set up random vector for accuracy test
      DO i=1,ncol
        CALL RANDOM_NUMBER(vec(i))
      END DO
c     full matrix x vec results
      rez1=MATMUL(FullMat,vec)


      acaMat%nrow=nrow
      acaMat%ncol=ncol
      acaMat%rank=rank
      ALLOCATE (acaMat%a(acaMat%nrow,acaMat%rank))
      ALLOCATE (acaMat%b(acaMat%rank,acaMat%ncol))


      DO irank=1,rank
        CALL FindLargestValue(FullMat,nrow,ncol,lrow,lcol,delta)
        CALL CopyCR2ab(FullMat,nrow,ncol,lrow,lcol,delta,acaMat%a,acaMat%b,rank,irank)
        CALL SubsApp(FullMat,nrow,ncol,acaMat%a,acaMat%b,rank,irank)
      END DO

c     verify with random vector multiplication
      CALL aca_abxv(acaMat,vec,rez2)
c     get norm
      CALL getRMSnorm(rez1,rez2,nrow,norm)

c      print *,nrow,ncol,rank,norm

      DEALLOCATE (vec,rez1,rez2)

      END


C     ------------------------------------------------------------------
      SUBROUTINE acaj_CrossAprox(FullMat,nrow,ncol,acaMat,eps,norm)
C
C     Izdela a,b = low rank approximation of fullmat, Frobenius norm, Jan version
C
C     ------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(acaMatrixType) :: acaMat
      REAL(8) eps,norm,s,s2
      INTEGER i,j

      INTEGER nrow,ncol,rank
      REAL(8) FullMat(nrow,ncol)

      INTEGER lrow,lcol,irank
      REAL(8) delta
      REAL(8), ALLOCATABLE :: vec(:),rez1(:),rez2(:)


      ALLOCATE ( vec(ncol),rez1(nrow),rez2(nrow) )
c     set up random vector for accuracy test
      DO i=1,ncol
        CALL RANDOM_NUMBER(vec(i))
      END DO
c     full matrix x vec results
      rez1=MATMUL(FullMat,vec)

c     do ACA
      CALL ACA_AB(acaMat,FullMat,nrow,ncol,EPS)

c     verify with random vector multiplication
      CALL aca_abxv(acaMat,vec,rez2)
c     get norm
      CALL getRMSnorm(rez1,rez2,nrow,norm)

      print *,nrow,ncol,acaMat%rank,norm

      DEALLOCATE (vec,rez1,rez2)

      END


C     ------------------------------------------------------------------
      SUBROUTINE aca_abxv(acaMat,vekX,vekB)
C
C     Mnozenje aca low rank app. z vektorjem
C
C     ------------------------------------------------------------------
      USE inc_types

      TYPE(acaMatrixType) :: acaMat

      REAL(8) vekX(acaMat%ncol)
      REAL(8) vekB(acaMat%nrow)

      REAL(8), ALLOCATABLE :: tmpv(:)

      ALLOCATE (tmpv(acaMat%rank))

      tmpv=MATMUL(acaMat%b,vekX)
      vekB=MATMUL(acaMat%a,tmpv)

      DEALLOCATE (tmpv)

      END

C     ------------------------------------------------------------------
      SUBROUTINE aca_abxvSLOW(acaMat,vekX,vekB)
C
C     Mnozenje aca low rank app. z vektorjem
C
C     ------------------------------------------------------------------
      USE inc_types

      TYPE(acaMatrixType) :: acaMat

      INTEGER irank,i,j

      REAL(8) vekX(acaMat%ncol)
      REAL(8) vekB(acaMat%nrow)

      REAL(8), ALLOCATABLE :: tmpv(:)

      ALLOCATE (tmpv(acaMat%rank))

      tmpv=0.0D0
      DO irank=1,acaMat%rank
        DO j=1,acaMat%ncol
          tmpv(irank)=tmpv(irank)+acaMat%b(irank,j)*vekX(j)
        END DO
      END DO

      vekB=0.0D0
      DO i=1,acaMat%nrow
        DO irank=1,acaMat%rank
          vekB(i)=vekB(i)+acaMat%a(i,irank)*tmpv(irank)
        END DO
      END DO

      DEALLOCATE (tmpv)

      END



C     ------------------------------------------------------------------
      SUBROUTINE SubsApp(FullMat,nrow,ncol,a,b,rank,irank)
C
C     Od matrike odsteje aproximacijo
C
C     ------------------------------------------------------------------
      INTEGER nrow,ncol,rank
      REAL(8) FullMat(nrow,ncol)

      INTEGER irank,i,j

      REAL(8) a(nrow,rank)
      REAL(8) b(rank,ncol)

      DO j=1,ncol
        DO i=1,nrow
          FullMat(i,j)=FullMat(i,j) - a(i,irank)*b(irank,j)
        END DO
      END DO

      END


C     ------------------------------------------------------------------
      SUBROUTINE CopyCR2ab(FullMat,nrow,ncol,lrow,lcol,delta,a,b,rank,irank)
C
C     Izpolni en stolpec v a in eno vrsto v b
C
C     ------------------------------------------------------------------
      INTEGER nrow,ncol,rank
      REAL(8) FullMat(nrow,ncol)

      INTEGER lrow,lcol,irank,i
      REAL(8) delta

      REAL(8) a(nrow,rank)
      REAL(8) b(rank,ncol)

      IF (delta.NE.0.0D0) THEN  ! če je cela matrika = 0
        DO i=1,nrow
          a(i,irank)=FullMat(i,lcol) / delta
        END DO
      ELSE
        DO i=1,nrow
          a(i,irank)=FullMat(i,lcol)
        END DO
      END IF

      DO i=1,ncol
        b(irank,i)=FullMat(lrow,i)
      END DO

      END

C     ------------------------------------------------------------------
      SUBROUTINE FindLargestValue(FullMat,nrow,ncol,lrow,lcol,delta)
C
C     Finds the larges element in FullMat in absolute sense
C
C     ------------------------------------------------------------------
      INTEGER nrow,ncol,lrow,lcol
      REAL(8) FullMat(nrow,ncol)
      REAL(8) delta,maxv

      INTEGER i,j

      maxv=-1.0D0

      DO j=1,ncol
        DO i=1,nrow
          IF (abs(FullMat(i,j)).GT.maxv) THEN
            maxv=abs(FullMat(i,j))
            lrow=i
            lcol=j
          END IF
        END DO
      END DO

      delta=FullMat(lrow,lcol)
      END



C---------------------------------------------------------------------------------------------------
      SUBROUTINE ACA_AB(acaMat,fMat,nrow,ncol,EPS)
C---------------------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE
      TYPE(acaMatrixType) :: acaMat

      INTEGER nrow,ncol
      REAL(8) fMat(nrow,ncol)
      REAL(8),DIMENSION(:,:),ALLOCATABLE::A,B
      INTEGER counter,X,DIFFR,DIFFC,RANKN
      INTEGER RANK
      REAL(8) EPS,NORM,SNORM

      RANK=MIN(nrow,ncol)

      ALLOCATE(A(nrow,RANK),B(RANK,ncol))
      X=0
      DIFFC=0
      DIFFR=0
      RANKN=0
      counter=RANK
      SNORM=0.0D0
      NORM=1.0D0
      CALL ADAPTIVE_CROSS_APPROXIMATION2(acaMat,counter,X,fMat,A,B,nrow,ncol,RANK,SNORM,NORM,EPS,RANKN,DIFFR,DIFFC)

      DEALLOCATE(A,B)

      END
C---------------------------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE ADAPTIVE_CROSS_APPROXIMATION2(acaMat,counter,X,R,A,B,M,N,RANK,SNORM,NORM,EPS,RANKN,DIFFR,DIFFC)
C---------------------------------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE
      REAL(8),DIMENSION(:),ALLOCATABLE::AVEC,BVEC
      INTEGER counter,M,N,RANK,STP,X,DIFFR,DIFFC,RANKN
      REAL(8) SNORM,SUM1,SUM2,SUM3,NORM,EPS
      REAL(8) R(M,N),A(M,RANK),B(RANK,N)
      TYPE(acaMatrixType) :: acaMat

      IF(counter.GE.1)THEN
        IF(NORM.GT.EPS)THEN
          ALLOCATE(AVEC(M),BVEC(N))
          CALL APPLAY_ACA(R,AVEC,BVEC,M,N)
          CALL EPS_ACA_NORM(A,B,AVEC,BVEC,SUM1,SUM2,SUM3,NORM,SNORM,M,N,RANK,counter)
          X=X+1
          RANKN=RANKN+1
        ELSE
          DIFFR=DIFFR+1
        END IF
        CALL ADAPTIVE_CROSS_APPROXIMATION2(acaMat,counter-1,X,R,A,B,M,N,RANK,SNORM,NORM,EPS,RANKN,DIFFR,DIFFC)
        IF(counter.GT.DIFFR.OR.DIFFR.EQ.RANK)THEN
          IF(DIFFC.EQ.0)THEN
            acaMat%nrow = M
            acaMat%ncol = N
            acaMat%rank = X
            ALLOCATE(acaMat%A(acaMat%nrow,acaMat%rank),acaMat%B(acaMat%rank,acaMat%ncol))
            DIFFC=X
          END IF
          CALL FILL_AB(acaMat%A,acaMat%B,AVEC,BVEC,M,N,DIFFC,X)
          X=X-1
          DEALLOCATE(AVEC,BVEC)
        END IF
      END IF

      END
C---------------------------------------------------------------------------------------------------
      SUBROUTINE EPS_ACA_NORM(A,B,AVEC,BVEC,SUM1,SUM2,SUM3,NORM,SNORM,M,N,RANK,counter)
C---------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER j,k
      INTEGER M,N,RANK,counter
      REAL(8) SUM1,SUM2,SUM3,NORM,SNORM
      REAL(8) AVEC(M),BVEC(N),A(M,RANK),B(RANK,N)
      SUM1=0.0D0
      SUM2=0.0D0
      SUM3=0.0D0
      DO j=1,M
       SUM1=SUM1+AVEC(j)**2
      END DO
      DO j=1,N
       SUM2=SUM2+BVEC(j)**2
      END DO
      DO j=1,M
       A(j,counter)=AVEC(j)
      END DO
      DO j=1,N
       B(counter,j)=BVEC(j)
      END DO
      SUM3=0.0D0
      DO k=counter,RANK
       DO j=1,M
        SUM3=SUM3+A(j,k)*AVEC(j)
       END DO
      END DO
      DO j=1,N
       DO k=counter,RANK
        SUM3=SUM3+B(k,j)*BVEC(j)
       END DO
      END DO
      SNORM=SNORM+SUM3+SUM1*SUM2
      NORM=SQRT(SUM1)*SQRT(SUM2)/SQRT(SNORM)
      END
C---------------------------------------------------------------------------------------------------
      SUBROUTINE APPLAY_ACA(R,AVEC,BVEC,M,N)
C---------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER i,j
      INTEGER iz,jz,M,N
      REAL(8) sigma,maxs
      REAL(8) R(M,N),AVEC(M),BVEC(N)
C     Poièe maksimalno vrednost v matriki in doloèi njegovo mesto
       maxs=-1.0D0
       DO i=1,M
        DO j=1,N
         IF(ABS(R(i,j)).GT.maxs) THEN
          iz=i
          jz=j
          maxs=ABS(R(i,j))
         END IF
        END DO
       END DO
C     Normira
      sigma=1.0D0/R(iz,jz)
C     doloèi vrednost matrike A
      DO i=1,M
       AVEC(i)=R(i,jz)
      END DO
C     doloèi vrednost matrike B
      DO j=1,N
       BVEC(j)=R(iz,j)*sigma
      END DO
C     Doloèi rezidualno matriko
      DO i=1,M
       DO j=1,N
        R(i,j)=R(i,j)-AVEC(i)*BVEC(j)
       END DO
      END DO
      END
C---------------------------------------------------------------------------------------------------
      SUBROUTINE FILL_AB(A,B,AVEC,BVEC,M,N,RANK,counter)
C---------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER i,j
      INTEGER M,N,RANK,counter,counter1,X
      REAL(8) A(M,RANK),B(RANK,N),AVEC(M),BVEC(N)

      DO j=1,M
       A(j,counter)=AVEC(j)
      END DO
      DO j=1,N
       B(counter,j)=BVEC(j)
      END DO
      END
C---------------------------------------------------------------------------------------------------


