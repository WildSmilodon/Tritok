C -----------------------------------------------------------------------------
      SUBROUTINE fmSolveModHelmDChmat(mesh,cpu,inp,env,io,A,RHS,Pivot,u,qu,velocity,rhsv,
     &               UmatB,UmatBp,UmatAbdx,UmatAbdy,UmatAbdz,UmatDx,UmatDy,UmatDz,
     &               QmatB,QmatBp,QmatAbdx,QmatAbdy,QmatAbdz,QmatDx,QmatDy,QmatDz,
     &               nit,diff,qdiff,diffGrad)
C
C     $: Slove, cTree version
C
C -----------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(meshType) :: mesh
      TYPE(CPUtype) cpu
      TYPE(inputtype) inp
      TYPE(IOtype) io
      TYPE(penv) :: env

      INTEGER i,ieq,iii,jjj,dn,node,row,col
      INTEGER nit


      REAL(8) rhsv(mesh%nnodes)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) u(mesh%nnodes)
      REAL(8) qu(mesh%nbnpof)
      REAL(8) diff(mesh%nnodes)
      REAL(8) diffGrad(mesh%nnodes,3)
      REAL(8) qDiff(mesh%nq)
      REAL(8) mu,mun

      TYPE(HmatrixType) :: UmatB,UmatBp,UmatAbdX,UmatAbdY,UmatAbdZ,UmatDx,UmatDy,UmatDz
      TYPE(HmatrixType) :: QmatB,QmatBp,QmatAbdX,QmatAbdY,QmatAbdZ,QmatDx,QmatDy,QmatDz

      REAL(8) A(mesh%mHfmnunk,mesh%mHfmnunk)
      REAL(8) RHS(mesh%mHfmnunk,mesh%mHnb)
      REAL(8) Pivot(mesh%mHfmnunk)

      REAL(8), ALLOCATABLE :: x(:), b(:), r(:)
      REAL(8), ALLOCATABLE :: rB(:),rDx(:),rDy(:),rDz(:),rAx(:),rAy(:),rAz(:),rBp(:)
      REAL(8), ALLOCATABLE :: qrB(:),qrDx(:),qrDy(:),qrDz(:),qrAx(:),qrAy(:),qrAz(:),qrBp(:)
      REAL(8), ALLOCATABLE :: vB(:),vDx(:),vDy(:),vDz(:),vAx(:),vAy(:),vAz(:),vBp(:)
      REAL(4) cptime,rcpu,rrcpu

      INTEGER lslvt,lpret,lprep,lmaxit,lstopt,lnits,ierr
      REAL(4) lslveps

        lslvt=0
        lpret=2
        lprep=2
        lmaxit=500
        lstopt=5
        lslveps=1.0E-15
        ierr=0

      ALLOCATE (x(mesh%mHfmnunk),b(mesh%mHnb),r(mesh%mHfmnunk))

      ALLOCATE (rB(mesh%nnodes),rBp(mesh%nnodes))
      ALLOCATE (rDx(mesh%nnodes),rDy(mesh%nnodes),rDz(mesh%nnodes))
      ALLOCATE (rAx(mesh%nnodes),rAy(mesh%nnodes),rAz(mesh%nnodes))

      ALLOCATE (qrB(mesh%nbnpof),qrBp(mesh%nbnpof))
      ALLOCATE (qrDx(mesh%nbnpof),qrDy(mesh%nbnpof),qrDz(mesh%nbnpof))
      ALLOCATE (qrAx(mesh%nbnpof),qrAy(mesh%nbnpof),qrAz(mesh%nbnpof))


      ALLOCATE (vB(mesh%nnodes),vBp(mesh%nnodes))
      ALLOCATE (vDx(mesh%nnodes),vDy(mesh%nnodes),vDz(mesh%nnodes))
      ALLOCATE (vAx(mesh%nbnodes),vAy(mesh%nbnodes),vAz(mesh%nbnodes))

      rcpu=cptime(0.)
C
C     Set up rhs vector and initial approximation of the unknown vector
C
      DO i=1,mesh%nnodes
        IF (mesh%mHkode(i).LT.0) THEN
          x(ABS(mesh%mHkode(i)))=u(i)
        ELSE
          b(mesh%mHkode(i))=u(i)
        END IF
      END DO
      DO i=1,mesh%nbnpof
        IF (mesh%mHqfmkode(i).LT.0) THEN
          x(ABS(mesh%mHqfmkode(i)))=qu(i)
        ELSE
          b(mesh%mHqfmkode(i))=qu(i)
        END IF
      END DO

C
C     r = rhs * b
C
      r=MATMUL(RHS,b)
c
c     r = r - B * rhsv
c

c     set up vectors for multiplictaion (base on boundary conditions)
C     obmocni integrali
      DO col=1,mesh%nnodes
          node=col
          vB(col) =-rhsv(node)
          vBp(col)=u(node)
          vDx(col)=u(node)*(velocity(node,1)+diffGrad(node,1))
          vDy(col)=u(node)*(velocity(node,2)+diffGrad(node,2))
          vDz(col)=u(node)*(velocity(node,3)+diffGrad(node,3))
      END DO
C       robni integrali
      DO col=1,mesh%nbnodes
          node=mesh%gbn(col)
          vAx(col)=-u(node)*velocity(node,1)
          vAy(col)=-u(node)*velocity(node,2)
          vAz(col)=-u(node)*velocity(node,3)
      END DO
c
c        u izvorne tocke
c
c     matmul
      CALL cTree_MxV(UmatB,vB,rB)
      CALL cTree_MxV(UmatBp,vBp,rBp)
      CALL cTree_MxV(UmatDx,vDx,rDx)
      CALL cTree_MxV(UmatDy,vDy,rDy)
      CALL cTree_MxV(UmatDz,vDz,rDz)
      CALL cTree_MxV(UmatAbdX,vAx,rAx)
      CALL cTree_MxV(UmatAbdY,vAy,rAy)
      CALL cTree_MxV(UmatAbdZ,vAz,rAz)

      ieq=0 ! stevilka enacbe
c     po u izvornih tockah
      DO row=1,mesh%nnodes
c        print *,row,rB(row)
        IF (mesh%mHkode(row).LT.0) THEN ! v tej izvorni tocki imam enacbo
          ieq=ieq+1
          r(ieq)=r(ieq)+rB(row)+rBp(row)+rDx(row)+rDy(row)+rDz(row)+rAx(row)+rAy(row)+rAz(row)
        END IF
      END DO
      DEALLOCATE (rB,rDx,rDy,rDz,rAx,rAy,rAz,rBp)
c
c        q izvorne tocke
c
c     matmul
      CALL cTree_MxV(QmatB,vB,qrB)
      CALL cTree_MxV(QmatBp,vBp,qrBp)
      CALL cTree_MxV(QmatDx,vDx,qrDx)
      CALL cTree_MxV(QmatDy,vDy,qrDy)
      CALL cTree_MxV(QmatDz,vDz,qrDz)
      CALL cTree_MxV(QmatAbdX,vAx,qrAx)
      CALL cTree_MxV(QmatAbdY,vAy,qrAy)
      CALL cTree_MxV(QmatAbdZ,vAz,qrAz)

c     po q izvornih tockah
C     Zanka po Q izvornih tockah
      DO iii=1,mesh%nbelem
        DO jjj=1,mesh%npof
          row=(iii-1)*mesh%npof+jjj
          IF (mesh%mHqfmkode(row).LT.0) THEN ! v tej izvorni tocki imam enacbo
            ieq=ieq+1
            r(ieq)=r(ieq)+qrB(row)+qrBp(row)+qrDx(row)+qrDy(row)+qrDz(row)+qrAx(row)+qrAy(row)+qrAz(row)
          END IF
        END DO
      END DO
      DEALLOCATE (qrB,qrDx,qrDy,qrDz,qrAx,qrAy,qrAz,qrBp)
      DEALLOCATE (vB,vDx,vDy,vDz,vAx,vAy,vAz,vBp)


      cpu%time(23)=cpu%time(23)+cptime(rcpu)
      rcpu=cptime(0.)

C
C     solve determined system of linear equations
C
      CALL SolvEQNfm(lslvt,lpret,lprep,lmaxit,lstopt,lslveps,
     &              mesh%mHfmnunk,Pivot,A,r,x,lnits,rrcpu,ierr)
      nit=1 ! ker je LU dekompozicija
      cpu%time(24)=cpu%time(24)+cptime(rcpu)
      IF (ierr.NE.0) CALL WarnErr(env,io,inp,4,"fmSolveModHelmDChmat","NAPAKA V SOVLERJU",ierr)

C
C     Copy solution vector to q and u using under-relaxation
c     under-relaxation should be 1.0 for linear problems
C
      DO i=1,mesh%nnodes
        IF (mesh%mHkode(i).LT.0) THEN
          u(i)=(1.0D0-inp%urmH)*u(i)+inp%urmH*x(ABS(mesh%mHkode(i)))
        END IF
      END DO
      DO i=1,mesh%nbnpof
        IF (mesh%mHqfmkode(i).LT.0) THEN
          qu(i)=(1.0D0-inp%urmH)*qu(i)+inp%urmH*x(ABS(mesh%mHqfmkode(i)))
        END IF
      END DO

      DEALLOCATE (x,b,r)


      END




C -----------------------------------------------------------------------------
      SUBROUTINE fmSolveModHelmDCaca(mesh,cpu,inp,env,io,A,RHS,Pivot,u,qu,velocity,rhsv,
     &               aUmatB,aUmatBp,aUmatAbdx,aUmatAbdy,aUmatAbdz,aUmatDx,aUmatDy,aUmatDz,
     &               aQmatB,aQmatBp,aQmatAbdx,aQmatAbdy,aQmatAbdz,aQmatDx,aQmatDy,aQmatDz,
     &               nit,diff,qdiff,diffGrad)
C
C     $: Slove, aca version
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(CPUtype) cpu
      TYPE(inputtype) inp
      TYPE(IOtype) io
      TYPE(penv) :: env

      INTEGER i,ieq,iii,jjj,dn,node,row,col
      INTEGER nit


      REAL(8) rhsv(mesh%nnodes)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) u(mesh%nnodes)
      REAL(8) qu(mesh%nbnpof)
      REAL(8) diff(mesh%nnodes)
      REAL(8) diffGrad(mesh%nnodes,3)
      REAL(8) qDiff(mesh%nq)
      REAL(8) mu,mun

      TYPE(acaMatrixType) :: aUmatB,aUmatBp,aUmatAbdX,aUmatAbdY,aUmatAbdZ,aUmatDx,aUmatDy,aUmatDz
      TYPE(acaMatrixType) :: aQmatB,aQmatBp,aQmatAbdX,aQmatAbdY,aQmatAbdZ,aQmatDx,aQmatDy,aQmatDz

      REAL(8) A(mesh%mHfmnunk,mesh%mHfmnunk)
      REAL(8) RHS(mesh%mHfmnunk,mesh%mHnb)
      REAL(8) Pivot(mesh%mHfmnunk)

      REAL(8), ALLOCATABLE :: x(:), b(:), r(:)
      REAL(8), ALLOCATABLE :: rB(:),rDx(:),rDy(:),rDz(:),rAx(:),rAy(:),rAz(:),rBp(:)
      REAL(8), ALLOCATABLE :: qrB(:),qrDx(:),qrDy(:),qrDz(:),qrAx(:),qrAy(:),qrAz(:),qrBp(:)
      REAL(8), ALLOCATABLE :: vB(:),vDx(:),vDy(:),vDz(:),vAx(:),vAy(:),vAz(:),vBp(:)
      REAL(4) cptime,rcpu,rrcpu

      INTEGER lslvt,lpret,lprep,lmaxit,lstopt,lnits,ierr
      REAL(4) lslveps

        lslvt=0
        lpret=2
        lprep=2
        lmaxit=500
        lstopt=5
        lslveps=1.0E-15
        ierr=0

      ALLOCATE (x(mesh%mHfmnunk),b(mesh%mHnb),r(mesh%mHfmnunk))

      ALLOCATE (rB(mesh%nnodes),rBp(mesh%nnodes))
      ALLOCATE (rDx(mesh%nnodes),rDy(mesh%nnodes),rDz(mesh%nnodes))
      ALLOCATE (rAx(mesh%nnodes),rAy(mesh%nnodes),rAz(mesh%nnodes))

      ALLOCATE (qrB(mesh%nbnpof),qrBp(mesh%nbnpof))
      ALLOCATE (qrDx(mesh%nbnpof),qrDy(mesh%nbnpof),qrDz(mesh%nbnpof))
      ALLOCATE (qrAx(mesh%nbnpof),qrAy(mesh%nbnpof),qrAz(mesh%nbnpof))


      ALLOCATE (vB(mesh%nnodes),vBp(mesh%nnodes))
      ALLOCATE (vDx(mesh%nnodes),vDy(mesh%nnodes),vDz(mesh%nnodes))
      ALLOCATE (vAx(mesh%nbnodes),vAy(mesh%nbnodes),vAz(mesh%nbnodes))

      rcpu=cptime(0.)
C
C     Set up rhs vector and initial approximation of the unknown vector
C
      DO i=1,mesh%nnodes
        IF (mesh%mHkode(i).LT.0) THEN
          x(ABS(mesh%mHkode(i)))=u(i)
        ELSE
          b(mesh%mHkode(i))=u(i)
        END IF
      END DO
      DO i=1,mesh%nbnpof
        IF (mesh%mHqfmkode(i).LT.0) THEN
          x(ABS(mesh%mHqfmkode(i)))=qu(i)
        ELSE
          b(mesh%mHqfmkode(i))=qu(i)
        END IF
      END DO

C
C     r = rhs * b
C
      r=MATMUL(RHS,b)
c
c     r = r - B * rhsv
c

c     set up vectors for multiplictaion (base on boundary conditions)
C     obmocni integrali
      DO col=1,mesh%nnodes
          node=col
          vB(col) =-rhsv(node)
          vBp(col)=u(node)
          vDx(col)=u(node)*(velocity(node,1)+diffGrad(node,1))
          vDy(col)=u(node)*(velocity(node,2)+diffGrad(node,2))
          vDz(col)=u(node)*(velocity(node,3)+diffGrad(node,3))
      END DO
C       robni integrali
      DO col=1,mesh%nbnodes
          node=mesh%gbn(col)
          vAx(col)=-u(node)*velocity(node,1)
          vAy(col)=-u(node)*velocity(node,2)
          vAz(col)=-u(node)*velocity(node,3)
      END DO
c
c        u izvorne tocke
c
c     matmul,
      CALL aca_abxv(aUmatB,vB,rB)
      CALL aca_abxv(aUmatBp,vBp,rBp)
      CALL aca_abxv(aUmatDx,vDx,rDx)
      CALL aca_abxv(aUmatDy,vDy,rDy)
      CALL aca_abxv(aUmatDz,vDz,rDz)
      CALL aca_abxv(aUmatAbdX,vAx,rAx)
      CALL aca_abxv(aUmatAbdY,vAy,rAy)
      CALL aca_abxv(aUmatAbdZ,vAz,rAz)

      ieq=0 ! stevilka enacbe
c     po u izvornih tockah
      DO row=1,mesh%nnodes
c        print *,row,rB(row)
        IF (mesh%mHkode(row).LT.0) THEN ! v tej izvorni tocki imam enacbo
          ieq=ieq+1
          r(ieq)=r(ieq)+rB(row)+rBp(row)+rDx(row)+rDy(row)+rDz(row)+rAx(row)+rAy(row)+rAz(row)
        END IF
      END DO
      DEALLOCATE (rB,rDx,rDy,rDz,rAx,rAy,rAz,rBp)
c
c        q izvorne tocke
c
c     matmul,
      CALL aca_abxv(aQmatB,vB,qrB)
      CALL aca_abxv(aQmatBp,vBp,qrBp)
      CALL aca_abxv(aQmatDx,vDx,qrDx)
      CALL aca_abxv(aQmatDy,vDy,qrDy)
      CALL aca_abxv(aQmatDz,vDz,qrDz)
      CALL aca_abxv(aQmatAbdX,vAx,qrAx)
      CALL aca_abxv(aQmatAbdY,vAy,qrAy)
      CALL aca_abxv(aQmatAbdZ,vAz,qrAz)

c     po q izvornih tockah
C     Zanka po Q izvornih tockah
      DO iii=1,mesh%nbelem
        DO jjj=1,mesh%npof
          row=(iii-1)*mesh%npof+jjj
          IF (mesh%mHqfmkode(row).LT.0) THEN ! v tej izvorni tocki imam enacbo
            ieq=ieq+1
            r(ieq)=r(ieq)+qrB(row)+qrBp(row)+qrDx(row)+qrDy(row)+qrDz(row)+qrAx(row)+qrAy(row)+qrAz(row)
          END IF
        END DO
      END DO
      DEALLOCATE (qrB,qrDx,qrDy,qrDz,qrAx,qrAy,qrAz,qrBp)
      DEALLOCATE (vB,vDx,vDy,vDz,vAx,vAy,vAz,vBp)


      cpu%time(23)=cpu%time(23)+cptime(rcpu)
      rcpu=cptime(0.)
C
C     solve determined system of linear equations
C
      CALL SolvEQNfm(lslvt,lpret,lprep,lmaxit,lstopt,lslveps,
     &              mesh%mHfmnunk,Pivot,A,r,x,lnits,rrcpu,ierr)
      nit=1 ! ker je LU dekompozicija
      cpu%time(24)=cpu%time(24)+cptime(rcpu)
      IF (ierr.NE.0) CALL WarnErr(env,io,inp,4,"fmSolveModHelmDCaca","NAPAKA V SOVLERJU",ierr)

C
C     Copy solution vector to q and u using under-relaxation
c     under-relaxation should be 1.0 for linear problems
C
      DO i=1,mesh%nnodes
        IF (mesh%mHkode(i).LT.0) THEN
          u(i)=(1.0D0-inp%urmH)*u(i)+inp%urmH*x(ABS(mesh%mHkode(i)))
        END IF
      END DO
      DO i=1,mesh%nbnpof
        IF (mesh%mHqfmkode(i).LT.0) THEN
          qu(i)=(1.0D0-inp%urmH)*qu(i)+inp%urmH*x(ABS(mesh%mHqfmkode(i)))
        END IF
      END DO

      DEALLOCATE (x,b,r)


      END




C -----------------------------------------------------------------------------
      SUBROUTINE fmSolveModHelmDCwawt(mesh,cpu,inp,env,io,A,RHS,Pivot,u,qu,velocity,rhsv,
     &               cUmatB,cUmatBp,cUmatAbdx,cUmatAbdy,cUmatAbdz,cUmatDx,cUmatDy,cUmatDz,
     &               cQmatB,cQmatBp,cQmatAbdx,cQmatAbdy,cQmatAbdz,cQmatDx,cQmatDy,cQmatDz,
     &               nit,diff,qdiff,diffGrad)
C
C     $: Slove, wawt version
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(CPUtype) cpu
      TYPE(inputtype) inp
      TYPE(IOtype) io
      TYPE(penv) :: env

      INTEGER i,ieq,iii,jjj,dn,node,row,col
      INTEGER nit


      REAL(8) rhsv(mesh%nnodes)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) u(mesh%nnodes)
      REAL(8) qu(mesh%nbnpof)
      REAL(8) diff(mesh%nnodes)
      REAL(8) diffGrad(mesh%nnodes,3)
      REAL(8) qDiff(mesh%nq)
      REAL(8) mu,mun

      TYPE(crsMatrixType) :: cUmatB,cUmatBp,cUmatAbdX,cUmatAbdY,cUmatAbdZ,cUmatDx,cUmatDy,cUmatDz
      TYPE(crsMatrixType) :: cQmatB,cQmatBp,cQmatAbdX,cQmatAbdY,cQmatAbdZ,cQmatDx,cQmatDy,cQmatDz

      REAL(8) A(mesh%mHfmnunk,mesh%mHfmnunk)
      REAL(8) RHS(mesh%mHfmnunk,mesh%mHnb)
      REAL(8) Pivot(mesh%mHfmnunk)

      REAL(8), ALLOCATABLE :: x(:), b(:), r(:)
      REAL(8), ALLOCATABLE :: rB(:),rDx(:),rDy(:),rDz(:),rAx(:),rAy(:),rAz(:),rBp(:)
      REAL(8), ALLOCATABLE :: qrB(:),qrDx(:),qrDy(:),qrDz(:),qrAx(:),qrAy(:),qrAz(:),qrBp(:)
      REAL(8), ALLOCATABLE :: vB(:),vDx(:),vDy(:),vDz(:),vAx(:),vAy(:),vAz(:),vBp(:)
      REAL(4) cptime,rcpu,rrcpu

      INTEGER lslvt,lpret,lprep,lmaxit,lstopt,lnits,ierr
      REAL(4) lslveps

        lslvt=0
        lpret=2
        lprep=2
        lmaxit=500
        lstopt=5
        lslveps=1.0E-15
        ierr=0

      ALLOCATE (x(mesh%mHfmnunk),b(mesh%mHnb),r(mesh%mHfmnunk))

      ALLOCATE (rB(mesh%nnodes),rBp(mesh%nnodes))
      ALLOCATE (rDx(mesh%nnodes),rDy(mesh%nnodes),rDz(mesh%nnodes))
      ALLOCATE (rAx(mesh%nnodes),rAy(mesh%nnodes),rAz(mesh%nnodes))

      ALLOCATE (qrB(mesh%nbnpof),qrBp(mesh%nbnpof))
      ALLOCATE (qrDx(mesh%nbnpof),qrDy(mesh%nbnpof),qrDz(mesh%nbnpof))
      ALLOCATE (qrAx(mesh%nbnpof),qrAy(mesh%nbnpof),qrAz(mesh%nbnpof))


      ALLOCATE (vB(mesh%nnodes),vBp(mesh%nnodes))
      ALLOCATE (vDx(mesh%nnodes),vDy(mesh%nnodes),vDz(mesh%nnodes))
      ALLOCATE (vAx(mesh%nbnodes),vAy(mesh%nbnodes),vAz(mesh%nbnodes))

      rcpu=cptime(0.)
C
C     Set up rhs vector and initial approximation of the unknown vector
C
      DO i=1,mesh%nnodes
        IF (mesh%mHkode(i).LT.0) THEN
          x(ABS(mesh%mHkode(i)))=u(i)
        ELSE
          b(mesh%mHkode(i))=u(i)
        END IF
      END DO
      DO i=1,mesh%nbnpof
        IF (mesh%mHqfmkode(i).LT.0) THEN
          x(ABS(mesh%mHqfmkode(i)))=qu(i)
        ELSE
          b(mesh%mHqfmkode(i))=qu(i)
        END IF
      END DO

C
C     r = rhs * b
C
      r=MATMUL(RHS,b)
c
c     r = r - B * rhsv
c

c     set up vectors for multiplictaion (base on boundary conditions)
C     obmocni integrali
      DO col=1,mesh%nnodes
          node=col
          vB(col) =-rhsv(node)
          vBp(col)=u(node)
          vDx(col)=u(node)*(velocity(node,1)+diffGrad(node,1))
          vDy(col)=u(node)*(velocity(node,2)+diffGrad(node,2))
          vDz(col)=u(node)*(velocity(node,3)+diffGrad(node,3))
      END DO
C       robni integrali
      DO col=1,mesh%nbnodes
          node=mesh%gbn(col)
          vAx(col)=-u(node)*velocity(node,1)
          vAy(col)=-u(node)*velocity(node,2)
          vAz(col)=-u(node)*velocity(node,3)
      END DO
c
c        u izvorne tocke
c
c     matmul, r = W^T(A*Wx)
      CALL Che_WTAWx(cUmatB,vB,mesh%nnodes,rB,mesh%nnodes)
      CALL Che_WTAWx(cUmatBp,vBp,mesh%nnodes,rBp,mesh%nnodes)
      CALL Che_WTAWx(cUmatDx,vDx,mesh%nnodes,rDx,mesh%nnodes)
      CALL Che_WTAWx(cUmatDy,vDy,mesh%nnodes,rDy,mesh%nnodes)
      CALL Che_WTAWx(cUmatDz,vDz,mesh%nnodes,rDz,mesh%nnodes)
      CALL Che_WTAWx(cUmatAbdX,vAx,mesh%nbnodes,rAx,mesh%nnodes)
      CALL Che_WTAWx(cUmatAbdY,vAy,mesh%nbnodes,rAy,mesh%nnodes)
      CALL Che_WTAWx(cUmatAbdZ,vAz,mesh%nbnodes,rAz,mesh%nnodes)

      ieq=0 ! stevilka enacbe
c     po u izvornih tockah
      DO row=1,mesh%nnodes
c        print *,row,rB(row)
        IF (mesh%mHkode(row).LT.0) THEN ! v tej izvorni tocki imam enacbo
          ieq=ieq+1
          r(ieq)=r(ieq)+rB(row)+rBp(row)+rDx(row)+rDy(row)+rDz(row)+rAx(row)+rAy(row)+rAz(row)
        END IF
      END DO
      DEALLOCATE (rB,rDx,rDy,rDz,rAx,rAy,rAz,rBp)
c
c        q izvorne tocke
c
c     matmul, r = W^T(A*Wx)
      CALL Che_WTAWx(cQmatB,vB,mesh%nnodes,qrB,mesh%nbnpof)
      CALL Che_WTAWx(cQmatBp,vBp,mesh%nnodes,qrBp,mesh%nbnpof)
      CALL Che_WTAWx(cQmatDx,vDx,mesh%nnodes,qrDx,mesh%nbnpof)
      CALL Che_WTAWx(cQmatDy,vDy,mesh%nnodes,qrDy,mesh%nbnpof)
      CALL Che_WTAWx(cQmatDz,vDz,mesh%nnodes,qrDz,mesh%nbnpof)
      CALL Che_WTAWx(cQmatAbdX,vAx,mesh%nbnodes,qrAx,mesh%nbnpof)
      CALL Che_WTAWx(cQmatAbdY,vAy,mesh%nbnodes,qrAy,mesh%nbnpof)
      CALL Che_WTAWx(cQmatAbdZ,vAz,mesh%nbnodes,qrAz,mesh%nbnpof)

c     po q izvornih tockah
C     Zanka po Q izvornih tockah
      DO iii=1,mesh%nbelem
        DO jjj=1,mesh%npof
          row=(iii-1)*mesh%npof+jjj
          IF (mesh%mHqfmkode(row).LT.0) THEN ! v tej izvorni tocki imam enacbo
            ieq=ieq+1
            r(ieq)=r(ieq)+qrB(row)+qrBp(row)+qrDx(row)+qrDy(row)+qrDz(row)+qrAx(row)+qrAy(row)+qrAz(row)
          END IF
        END DO
      END DO
      DEALLOCATE (qrB,qrDx,qrDy,qrDz,qrAx,qrAy,qrAz,qrBp)
      DEALLOCATE (vB,vDx,vDy,vDz,vAx,vAy,vAz,vBp)


      cpu%time(23)=cpu%time(23)+cptime(rcpu)
      rcpu=cptime(0.)
C
C     solve determined system of linear equations
C
      CALL SolvEQNfm(lslvt,lpret,lprep,lmaxit,lstopt,lslveps,
     &              mesh%mHfmnunk,Pivot,A,r,x,lnits,rrcpu,ierr)
      nit=1 ! ker je LU dekompozicija
      cpu%time(24)=cpu%time(24)+cptime(rcpu)
      IF (ierr.NE.0) CALL WarnErr(env,io,inp,4,"fmSolveModHelmDCwawt","NAPAKA V SOVLERJU",ierr)

C
C     Copy solution vector to q and u using under-relaxation
c     under-relaxation should be 1.0 for linear problems
C
      DO i=1,mesh%nnodes
        IF (mesh%mHkode(i).LT.0) THEN
          u(i)=(1.0D0-inp%urmH)*u(i)+inp%urmH*x(ABS(mesh%mHkode(i)))
        END IF
      END DO
      DO i=1,mesh%nbnpof
        IF (mesh%mHqfmkode(i).LT.0) THEN
          qu(i)=(1.0D0-inp%urmH)*qu(i)+inp%urmH*x(ABS(mesh%mHqfmkode(i)))
        END IF
      END DO

      DEALLOCATE (x,b,r)


      END


C -----------------------------------------------------------------------------
      SUBROUTINE fmSolveModHelmDCmm(mesh,cpu,inp,env,io,A,RHS,Pivot,u,qu,velocity,rhsv,
     &               UmatB,UmatBp,UmatAbdx,UmatAbdy,UmatAbdz,UmatDx,UmatDy,UmatDz,
     &               QmatB,QmatBp,QmatAbdx,QmatAbdy,QmatAbdz,QmatDx,QmatDy,QmatDz,
     &               nit,diff,qdiff,diffGrad)
C
C     $: Slove, matmul version
C

C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(CPUtype) cpu
      TYPE(inputtype) inp
      TYPE(IOtype) io
      TYPE(penv) :: env

      INTEGER i,ieq,iii,jjj,dn,node,row,col
      INTEGER nit


      REAL(8) rhsv(mesh%nnodes)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) u(mesh%nnodes)
      REAL(8) qu(mesh%nbnpof)
      REAL(8) diff(mesh%nnodes)
      REAL(8) diffGrad(mesh%nnodes,3)
      REAL(8) qDiff(mesh%nq)
      REAL(8) mu,mun

      REAL(8) UmatB(mesh%nnodes,mesh%nnodes)   ! u* po omega
      REAL(8) UmatBp(mesh%nnodes,mesh%nnodes)   ! u* po omega
      REAL(8) UmatAbdx(mesh%nnodes,mesh%nbnodes) ! u* krat nx po gamma
      REAL(8) UmatAbdy(mesh%nnodes,mesh%nbnodes) ! u* krat ny po gamma
      REAL(8) UmatAbdz(mesh%nnodes,mesh%nbnodes) ! u* krat nz po gamma
      REAL(8) UmatDx(mesh%nnodes,mesh%nnodes) ! grad u*_x po omega
      REAL(8) UmatDy(mesh%nnodes,mesh%nnodes) ! grad u*_y po omega
      REAL(8) UmatDz(mesh%nnodes,mesh%nnodes) ! grad u*_z po omega

      REAL(8) QmatB(mesh%nbnpof,mesh%nnodes)   ! u* po omega
      REAL(8) QmatBp(mesh%nbnpof,mesh%nnodes)   ! u* po omega
      REAL(8) QmatAbdx(mesh%nbnpof,mesh%nbnodes) ! u* krat nx po gamma
      REAL(8) QmatAbdy(mesh%nbnpof,mesh%nbnodes) ! u* krat ny po gamma
      REAL(8) QmatAbdz(mesh%nbnpof,mesh%nbnodes) ! u* krat nz po gamma
      REAL(8) QmatDx(mesh%nbnpof,mesh%nnodes) ! grad u*_x po omega
      REAL(8) QmatDy(mesh%nbnpof,mesh%nnodes) ! grad u*_y po omega
      REAL(8) QmatDz(mesh%nbnpof,mesh%nnodes) ! grad u*_z po omega

      REAL(8) A(mesh%mHfmnunk,mesh%mHfmnunk)
      REAL(8) RHS(mesh%mHfmnunk,mesh%mHnb)
      REAL(8) Pivot(mesh%mHfmnunk)

      REAL(8), ALLOCATABLE :: x(:), b(:), r(:)
      REAL(8), ALLOCATABLE :: rB(:),rDx(:),rDy(:),rDz(:),rAx(:),rAy(:),rAz(:),rBp(:)
      REAL(8), ALLOCATABLE :: qrB(:),qrDx(:),qrDy(:),qrDz(:),qrAx(:),qrAy(:),qrAz(:),qrBp(:)
      REAL(8), ALLOCATABLE :: vB(:),vDx(:),vDy(:),vDz(:),vAx(:),vAy(:),vAz(:),vBp(:)
      REAL(4) cptime,rcpu,rrcpu

      INTEGER lslvt,lpret,lprep,lmaxit,lstopt,lnits,ierr
      REAL(4) lslveps

        lslvt=0
        lpret=2
        lprep=2
        lmaxit=500
        lstopt=5
        lslveps=1.0E-15
        ierr=0

      ALLOCATE (x(mesh%mHfmnunk),b(mesh%mHnb),r(mesh%mHfmnunk))

      ALLOCATE (rB(mesh%nnodes),rBp(mesh%nnodes))
      ALLOCATE (rDx(mesh%nnodes),rDy(mesh%nnodes),rDz(mesh%nnodes))
      ALLOCATE (rAx(mesh%nnodes),rAy(mesh%nnodes),rAz(mesh%nnodes))

      ALLOCATE (qrB(mesh%nbnpof),qrBp(mesh%nbnpof))
      ALLOCATE (qrDx(mesh%nbnpof),qrDy(mesh%nbnpof),qrDz(mesh%nbnpof))
      ALLOCATE (qrAx(mesh%nbnpof),qrAy(mesh%nbnpof),qrAz(mesh%nbnpof))

      ALLOCATE (vB(mesh%nnodes),vBp(mesh%nnodes))
      ALLOCATE (vDx(mesh%nnodes),vDy(mesh%nnodes),vDz(mesh%nnodes))
      ALLOCATE (vAx(mesh%nbnodes),vAy(mesh%nbnodes),vAz(mesh%nbnodes))

      rcpu=cptime(0.)
C
C     Set up rhs vector and initial approximation of the unknown vector
C
      DO i=1,mesh%nnodes
        IF (mesh%mHkode(i).LT.0) THEN
          x(ABS(mesh%mHkode(i)))=u(i)
        ELSE
          b(mesh%mHkode(i))=u(i)
        END IF
      END DO
      DO i=1,mesh%nbnpof
        IF (mesh%mHqfmkode(i).LT.0) THEN
          x(ABS(mesh%mHqfmkode(i)))=qu(i)
        ELSE
          b(mesh%mHqfmkode(i))=qu(i)
        END IF
      END DO

C
C     r = rhs * b
C
      r=MATMUL(RHS,b)
c
c     r = r - B * rhsv
c

c     set up vectors for multiplictaion (base on boundary conditions)
C     obmocni integrali
      DO col=1,mesh%nnodes
          node=col
          vB(col) =-rhsv(node)
          vBp(col)=u(node)
          vDx(col)=u(node)*(velocity(node,1)+diffGrad(node,1))
          vDy(col)=u(node)*(velocity(node,2)+diffGrad(node,2))
          vDz(col)=u(node)*(velocity(node,3)+diffGrad(node,3))
      END DO
C       robni integrali
      DO col=1,mesh%nbnodes
          node=mesh%gbn(col)
          vAx(col)=-u(node)*velocity(node,1)
          vAy(col)=-u(node)*velocity(node,2)
          vAz(col)=-u(node)*velocity(node,3)
      END DO
c
c        u izvorne tocke
c
      rB =MATMUL(UmatB,vB)
      rBp=MATMUL(UmatBp,vBp)
      rDx=MATMUL(UmatDx,vDx)
      rDy=MATMUL(UmatDy,vDy)
      rDz=MATMUL(UmatDz,vDz)
      rAx=MATMUL(UmatAbdx,vAx)
      rAy=MATMUL(UmatAbdy,vAy)
      rAz=MATMUL(UmatAbdz,vAz)

      ieq=0 ! stevilka enacbe
c     po u izvornih tockah
      DO row=1,mesh%nnodes
        IF (mesh%mHkode(row).LT.0) THEN ! v tej izvorni tocki imam enacbo
          ieq=ieq+1
          r(ieq)=r(ieq)+rB(row)+rBp(row)+rDx(row)+rDy(row)+rDz(row)+rAx(row)+rAy(row)+rAz(row)
        END IF
      END DO
c
c        q izvorne tocke
c
      qrB =MATMUL(QmatB,vB)
      qrBp=MATMUL(QmatBp,vBp)
      qrDx=MATMUL(QmatDx,vDx)
      qrDy=MATMUL(QmatDy,vDy)
      qrDz=MATMUL(QmatDz,vDz)
      qrAx=MATMUL(QmatAbdx,vAx)
      qrAy=MATMUL(QmatAbdy,vAy)
      qrAz=MATMUL(QmatAbdz,vAz)

c     po q izvornih tockah
C     Zanka po Q izvornih tockah
      DO iii=1,mesh%nbelem
        DO jjj=1,mesh%npof
          row=(iii-1)*mesh%npof+jjj
          IF (mesh%mHqfmkode(row).LT.0) THEN ! v tej izvorni tocki imam enacbo
            ieq=ieq+1
            r(ieq)=r(ieq)+qrB(row)+qrBp(row)+qrDx(row)+qrDy(row)+qrDz(row)+qrAx(row)+qrAy(row)+qrAz(row)
          END IF
        END DO
      END DO

      cpu%time(23)=cpu%time(23)+cptime(rcpu)
      rcpu=cptime(0.)
C
C     solve determined system of linear equations
      CALL SolvEQNfm(lslvt,lpret,lprep,lmaxit,lstopt,lslveps,
     &              mesh%mHfmnunk,Pivot,A,r,x,lnits,rrcpu,ierr)
      nit=1 ! ker je LU dekompozicija
      cpu%time(24)=cpu%time(24)+cptime(rcpu)
      IF (ierr.NE.0) CALL WarnErr(env,io,inp,4,"fmSolveModHelmDC","NAPAKA V SOVLERJU",ierr)

C
C     Copy solution vector to q and u using under-relaxation
c     under-relaxation should be 1.0 for linear problems
C
      DO i=1,mesh%nnodes
        IF (mesh%mHkode(i).LT.0) THEN
          u(i)=(1.0D0-inp%urmH)*u(i)+inp%urmH*x(ABS(mesh%mHkode(i)))
        END IF
      END DO
      DO i=1,mesh%nbnpof
        IF (mesh%mHqfmkode(i).LT.0) THEN
          qu(i)=(1.0D0-inp%urmH)*qu(i)+inp%urmH*x(ABS(mesh%mHqfmkode(i)))
        END IF
      END DO

      DEALLOCATE (x,b,r)
      DEALLOCATE (rB,rDx,rDy,rDz,rAx,rAy,rAz,rBp)
      DEALLOCATE (qrB,qrDx,qrDy,qrDz,qrAx,qrAy,qrAz,qrBp)
      DEALLOCATE (vB,vDx,vDy,vDz,vAx,vAy,vAz,vBp)

      END

C -----------------------------------------------------------------------------
      SUBROUTINE MakeBprime(mesh,inp,UmatB,UmatBp,QmatB,QmatBp,diff,qdiff)
C
C     $: creates Bp matrices
C

C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(inputtype) inp

      INTEGER i,iii,jjj,dn,node,row,col
      INTEGER nit

      REAL(8) diff(mesh%nnodes)
      REAL(8) qDiff(mesh%nq)
      REAL(8) mu,mun

      REAL(8) UmatB(mesh%nnodes,mesh%nnodes)   ! u* po omega
      REAL(8) QmatB(mesh%nbnpof,mesh%nnodes)   ! u* po omega

      REAL(8) UmatBp(mesh%nnodes,mesh%nnodes)   ! u* po omega
      REAL(8) QmatBp(mesh%nbnpof,mesh%nnodes)   ! u* po omega

c
c        u izvorne tocke
c

c     zaradi tega, ker je mu za vsako vrsto drugacen
c     po u izvornih tockah
      DO row=1,mesh%nnodes
        mu=SQRT(inp%beta/diff(row)) ! mu v izvorni tocki (ta je sel v integral)
        DO col=1,mesh%nnodes
           node=col
           mun=SQRT(inp%beta/diff(node)) ! mu v obmocni tocki
           UmatBp(row,col)=UmatB(row,col)*inp%beta*( (mu/mun)**2 - 1.0D0)
        END DO
      END DO
c
c        q izvorne tocke
c
      DO iii=1,mesh%nbelem
        DO jjj=1,mesh%npof
          row=(iii-1)*mesh%npof+jjj
          dn=mesh%ibcf(iii,jjj)
          mu=SQRT(inp%beta/qDiff(dn)) ! mu v izvorni tocki (ta je sel v integral)
          DO col=1,mesh%nnodes
            node=col
            mun=SQRT(inp%beta/diff(node)) ! mu v obmocni tocki
            QmatBp(row,col)=QmatB(row,col)*inp%beta*( (mu/mun)**2 - 1.0D0)
          END DO
        END DO
      END DO

      END



C -----------------------------------------------------------------------------
      SUBROUTINE fmSolveModHelmDC(mesh,cpu,inp,env,io,A,RHS,Pivot,u,qu,velocity,rhsv,
     &               UmatB,UmatAbdx,UmatAbdy,UmatAbdz,UmatDx,UmatDy,UmatDz,
     &               QmatB,QmatAbdx,QmatAbdy,QmatAbdz,QmatDx,QmatDy,QmatDz,
     &               nit,diff,qdiff,diffGrad)
C
C     $: Set up system and rhs matrices
C

C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(CPUtype) cpu
      TYPE(inputtype) inp
      TYPE(IOtype) io
      TYPE(penv) :: env

      INTEGER i,ieq,iii,jjj,dn,node,row,col
      INTEGER nit


      REAL(8) rhsv(mesh%nnodes)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) u(mesh%nnodes)
      REAL(8) qu(mesh%nbnpof)
      REAL(8) diff(mesh%nnodes)
      REAL(8) diffGrad(mesh%nnodes,3)
      REAL(8) qDiff(mesh%nq)
      REAL(8) mu,mun

      REAL(8) UmatB(mesh%nnodes,mesh%nnodes)   ! u* po omega
      REAL(8) UmatAbdx(mesh%nnodes,mesh%nbnodes) ! u* krat nx po gamma
      REAL(8) UmatAbdy(mesh%nnodes,mesh%nbnodes) ! u* krat ny po gamma
      REAL(8) UmatAbdz(mesh%nnodes,mesh%nbnodes) ! u* krat nz po gamma
      REAL(8) UmatDx(mesh%nnodes,mesh%nnodes) ! grad u*_x po omega
      REAL(8) UmatDy(mesh%nnodes,mesh%nnodes) ! grad u*_y po omega
      REAL(8) UmatDz(mesh%nnodes,mesh%nnodes) ! grad u*_z po omega

      REAL(8) QmatB(mesh%nbnpof,mesh%nnodes)   ! u* po omega
      REAL(8) QmatAbdx(mesh%nbnpof,mesh%nbnodes) ! u* krat nx po gamma
      REAL(8) QmatAbdy(mesh%nbnpof,mesh%nbnodes) ! u* krat ny po gamma
      REAL(8) QmatAbdz(mesh%nbnpof,mesh%nbnodes) ! u* krat nz po gamma
      REAL(8) QmatDx(mesh%nbnpof,mesh%nnodes) ! grad u*_x po omega
      REAL(8) QmatDy(mesh%nbnpof,mesh%nnodes) ! grad u*_y po omega
      REAL(8) QmatDz(mesh%nbnpof,mesh%nnodes) ! grad u*_z po omega

      REAL(8) A(mesh%mHfmnunk,mesh%mHfmnunk)
      REAL(8) RHS(mesh%mHfmnunk,mesh%mHnb)
      REAL(8) Pivot(mesh%mHfmnunk)

      REAL(8), ALLOCATABLE :: x(:), b(:), r(:)
      REAL(4) cptime,rcpu,rrcpu

      INTEGER lslvt,lpret,lprep,lmaxit,lstopt,lnits,ierr
      REAL(4) lslveps

        lslvt=0
        lpret=2
        lprep=2
        lmaxit=500
        lstopt=5
        lslveps=1.0E-15
        ierr=0

      ALLOCATE (x(mesh%mHfmnunk),b(mesh%mHnb),r(mesh%mHfmnunk))

      rcpu=cptime(0.)
C
C     Set up rhs vector and initial approximation of the unknown vector
C
      DO i=1,mesh%nnodes
        IF (mesh%mHkode(i).LT.0) THEN
          x(ABS(mesh%mHkode(i)))=u(i)
        ELSE
          b(mesh%mHkode(i))=u(i)
        END IF
      END DO
      DO i=1,mesh%nbnpof
        IF (mesh%mHqfmkode(i).LT.0) THEN
          x(ABS(mesh%mHqfmkode(i)))=qu(i)
        ELSE
          b(mesh%mHqfmkode(i))=qu(i)
        END IF
      END DO

C
C     r = rhs * b
C
      r=MATMUL(RHS,b)
c
c     r = r - B * rhsv
c
      ieq=0 ! stevilka enacbe
c     po u izvornih tockah
      DO row=1,mesh%nnodes
        IF (mesh%mHkode(row).LT.0) THEN ! v tej izvorni tocki imam enacbo
          ieq=ieq+1
          mu=SQRT(inp%beta/diff(row)) ! mu v izvorni tocki (ta je sel v integral)

C         obmocni integrali
          DO col=1,mesh%nnodes
              node=col
              mun=SQRT(inp%beta/diff(node)) ! mu v obmocni tocki
              r(ieq)=r(ieq)
C                    time derivative + sources
     &               -UmatB(row,col)*rhsv(node)
     &               -u(node)*(
     &                 -UmatDx(row,col)*(velocity(node,1)+diffGrad(node,1))
     &                 -UmatDy(row,col)*(velocity(node,2)+diffGrad(node,2))
     &                 -UmatDz(row,col)*(velocity(node,3)+diffGrad(node,3))
     &                  -UmatB(row,col)*inp%beta*( (mu/mun)**2 - 1.0D0)
     &                )
          END DO
C         robni integrali
          DO col=1,mesh%nbnodes
              node=mesh%gbn(col)
              r(ieq)=r(ieq)
     &               -u(node)*(
C                    advection
     &               +UmatAbdx(row,col)*velocity(node,1)
     &               +UmatAbdy(row,col)*velocity(node,2)
     &               +UmatAbdz(row,col)*velocity(node,3)
     &                )
           END DO

        END IF
      END DO

c     po q izvornih tockah
C     Zanka po Q izvornih tockah
      DO iii=1,mesh%nbelem
        DO jjj=1,mesh%npof
          row=(iii-1)*mesh%npof+jjj
          dn=mesh%ibcf(iii,jjj)
          mu=SQRT(inp%beta/qDiff(dn)) ! mu v izvorni tocki (ta je sel v integral)
          IF (mesh%mHqfmkode(row).LT.0) THEN ! v tej izvorni tocki imam enacbo
            ieq=ieq+1
C           obmocni integrali
            DO col=1,mesh%nnodes
              node=col
              mun=SQRT(inp%beta/diff(node)) ! mu v obmocni tocki
              r(ieq)=r(ieq)
C                    time derivative + sources
     &               -QmatB(row,col)*rhsv(node)
     &               -u(node)*(
C                    advection
     &                 -QmatDx(row,col)*(velocity(node,1)+diffGrad(node,1))
     &                 -QmatDy(row,col)*(velocity(node,2)+diffGrad(node,2))
     &                 -QmatDz(row,col)*(velocity(node,3)+diffGrad(node,3))
     &                  -QmatB(row,col)*inp%beta*( (mu/mun)**2 - 1.0D0)
     &                )
            END DO
C           robni integrali
            DO col=1,mesh%nbnodes
              node=mesh%gbn(col)
              r(ieq)=r(ieq)
     &               -u(node)*(
C                    advection
     &               +QmatAbdx(row,col)*velocity(node,1)
     &               +QmatAbdy(row,col)*velocity(node,2)
     &               +QmatAbdz(row,col)*velocity(node,3)
     &                )
            END DO
          END IF
        END DO
      END DO

      cpu%time(23)=cpu%time(23)+cptime(rcpu)
      rcpu=cptime(0.)
C
C     solve determined system of linear equations
      CALL SolvEQNfm(lslvt,lpret,lprep,lmaxit,lstopt,lslveps,
     &              mesh%mHfmnunk,Pivot,A,r,x,lnits,rrcpu,ierr)
      nit=1 ! ker je LU dekompozicija
      cpu%time(24)=cpu%time(24)+cptime(rcpu)
      IF (ierr.NE.0) CALL WarnErr(env,io,inp,4,"fmSolveModHelmDC","NAPAKA V SOVLERJU",ierr)



C
C     Copy solution vector to q and u using under-relaxation
c     under-relaxation should be 1.0 for linear problems
C
      DO i=1,mesh%nnodes
        IF (mesh%mHkode(i).LT.0) THEN
          u(i)=(1.0D0-inp%urmH)*u(i)+inp%urmH*x(ABS(mesh%mHkode(i)))
        END IF
      END DO
      DO i=1,mesh%nbnpof
        IF (mesh%mHqfmkode(i).LT.0) THEN
          qu(i)=(1.0D0-inp%urmH)*qu(i)+inp%urmH*x(ABS(mesh%mHqfmkode(i)))
        END IF
      END DO

      DEALLOCATE (x,b,r)

      END




C -----------------------------------------------------------------------------
      SUBROUTINE fmmHSysM(mesh,fmA,fmRHS,diff,qdiff,UmatH,UmatG,QmatH,QmatG)
C
C     $: Set up system and rhs matrices
C

C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh

      INTEGER ieq,ii,jj,j,nuq,row,col,iii,jjj,dn
      REAL(8) val

      REAL(8) diff(mesh%nnodes)
      REAL(8) qdiff(mesh%nq)

      REAL(8) UmatH(mesh%nnodes,mesh%nbnodes)   ! grad u* skalarno normala po gamma
      REAL(8) UmatG(mesh%nnodes,mesh%nbnpof)  ! u* po gamma

      REAL(8) QmatH(mesh%nbnpof,mesh%nbnodes)   ! grad u* skalarno normala po gamma
      REAL(8) QmatG(mesh%nbnpof,mesh%nbnpof)  ! u* po gamma

      REAL(8) fmA(mesh%mHfmnunk,mesh%mHfmnunk)
      REAL(8) fmRHS(mesh%mHfmnunk,mesh%mHnb)

      fmA=-99999.999D0
      fmRHS=-99999.999D0
c     po u izvornih tockah
      ieq=0 ! stevilka enacbe
      DO row=1,mesh%nnodes
        IF (mesh%mHkode(row).LT.0) THEN ! v tej izvorni tocki imam enacbo
          ieq=ieq+1
C
C         najprej H matrika (u)
C
            DO j=1,mesh%nnodes ! zanka po stolpcih pravokotne matrike
              val=0.0D0
              IF (mesh%lbn(j).NE.0) THEN ! je robna
                col=mesh%lbn(j)
                val=val+diff(j)*UmatH(row,col)
              END IF

              nuq=mesh%mHkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                fmA(ieq,ABS(nuq))=val
              ELSE ! znana vrednost - rhs
                fmRHS(ieq,nuq)=-val  ! minus zato, ker gre na drugo stran enacbe
              END IF
            END DO
C
C           nato G matrika (q)
C
            DO ii=1,mesh%nbelem
              DO jj=1,mesh%npof
                j=mesh%ibcf(ii,jj)
                col=(ii-1)*mesh%npof+jj
                nuq=mesh%mHqfmkode(col) ! stevilka v vektorju neznank oziroma znank
                IF (nuq.LT.0) THEN ! neznaka - sys
                  fmA(ieq,ABS(nuq))=-UmatG(row,col)*qDiff(j) ! v bistvu bi rabil samo qDiff po robu
                ELSE ! znana vrednost - rhs
                  fmRHS(ieq,nuq)=UmatG(row,col)*qDiff(j) ! minus zato, ker gre na drugo stran enacbe
                END IF
              END DO
            END DO

          IF (mesh%lbn(row).EQ.0) THEN ! ta tocka je v obmocju, dodati moram c(ksi) na diagonalo
              fmA(ieq,ieq)=fmA(ieq,ieq)+1.0D0*diff(row) ! v obmocju je c=1
          END IF

        END IF
      END DO
C     Zanka po Q izvornih tockah
      DO iii=1,mesh%nbelem
        DO jjj=1,mesh%npof
          row=(iii-1)*mesh%npof+jjj
          dn=mesh%ibcf(iii,jjj)
          IF (mesh%mHqfmkode(row).LT.0) THEN ! v tej izvorni tocki imam enacbo
            ieq=ieq+1
C
C           najprej H matrika (u)
C
            DO j=1,mesh%nnodes ! zanka po stolpcih pravokotne matrike
              val=0.0D0
              IF (mesh%lbn(j).NE.0) THEN ! je robna
                col=mesh%lbn(j)
                val=diff(j)*QmatH(row,col)
              END IF

              nuq=mesh%mHkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                fmA(ieq,ABS(nuq))=val
              ELSE ! znana vrednost - rhs
                fmRHS(ieq,nuq)=-val  ! minus zato, ker gre na drugo stran enacbe
              END IF
            END DO
C
C           nato G matrika (q)
C
            DO ii=1,mesh%nbelem
              DO jj=1,mesh%npof
                j=mesh%ibcf(ii,jj)
                col=(ii-1)*mesh%npof+jj
                nuq=mesh%mHqfmkode(col) ! stevilka v vektorju neznank oziroma znank
                IF (nuq.LT.0) THEN ! neznaka - sys
                  fmA(ieq,ABS(nuq))=-QmatG(row,col)*qDiff(j) ! v bistvu bi rabil samo qDiff po robu
                ELSE ! znana vrednost - rhs
                  fmRHS(ieq,nuq)=QmatG(row,col)*qDiff(j) ! minus zato, ker gre na drugo stran enacbe
                END IF
              END DO
            END DO
          END IF
        END DO
      END DO

C     mesh%mHnb     ! ta pove koliko je znanih u ali q
c     mesh%mHfmnunk ! ta pove koliko je neznanih u ali q
      END




C -----------------------------------------------------------------------------
      SUBROUTINE FMATsdmHfs(env,io,inp,mesh,gauss,mHdiff,mHqDiff,
     &               UmatH,UmatG,UmatB,UmatAbdx,UmatAbdy,UmatAbdz,UmatDx,UmatDy,UmatDz,
     &               QmatH,QmatG,QmatB,QmatAbdx,QmatAbdy,QmatAbdz,QmatDx,QmatDy,QmatDz,
     &                          ontH,ontG,ontDx,ontDy,ontDz,mumax)
C
C     $: Form Matrices, unsteady diffusion advection equation and kinematics??
C        (H, G, Ab, Adx, Ady, Adz, B, Htx, Hty, Htz, Dx, Dy, Dz)
C

C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(IOtype) io
      TYPE(inputtype) inp
      TYPE(penv) :: env

      INTEGER i,ic,it,je,isrc,isip,j,itf,iti
      REAL(8) xp,yp,zp,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
      REAL(8) x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8
      REAL(8), ALLOCATABLE :: ge(:),he(:),be(:),dxe(:),dye(:),dze(:),abx(:),aby(:),abz(:)
      REAL(8), ALLOCATABLE :: heyz(:),hezx(:),hexy(:) ! teh ne rabim, samo da je ista rutina

c     Difuzivnost za osnovno resitev
      REAL(8) mHdiff(mesh%nnodes),mHqDiff(mesh%nq)

      REAL(8) UmatH(mesh%nnodes,mesh%nbnodes)   ! grad u* skalarno normala po gamma
      REAL(8) UmatG(mesh%nnodes,mesh%nbnpof)  ! u* po gamma
      REAL(8) UmatB(mesh%nnodes,mesh%nnodes)   ! u* po omega
      REAL(8) UmatAbdx(mesh%nnodes,mesh%nbnodes) ! u* krat nx po gamma
      REAL(8) UmatAbdy(mesh%nnodes,mesh%nbnodes) ! u* krat ny po gamma
      REAL(8) UmatAbdz(mesh%nnodes,mesh%nbnodes) ! u* krat nz po gamma
      REAL(8) UmatDx(mesh%nnodes,mesh%nnodes) ! grad u*_x po omega
      REAL(8) UmatDy(mesh%nnodes,mesh%nnodes) ! grad u*_y po omega
      REAL(8) UmatDz(mesh%nnodes,mesh%nnodes) ! grad u*_z po omega

      REAL(8) QmatH(mesh%nbnpof,mesh%nbnodes)   ! grad u* skalarno normala po gamma
      REAL(8) QmatG(mesh%nbnpof,mesh%nbnpof)  ! u* po gamma
      REAL(8) QmatB(mesh%nbnpof,mesh%nnodes)   ! u* po omega
      REAL(8) QmatAbdx(mesh%nbnpof,mesh%nbnodes) ! u* krat nx po gamma
      REAL(8) QmatAbdy(mesh%nbnpof,mesh%nbnodes) ! u* krat ny po gamma
      REAL(8) QmatAbdz(mesh%nbnpof,mesh%nbnodes) ! u* krat nz po gamma
      REAL(8) QmatDx(mesh%nbnpof,mesh%nnodes) ! grad u*_x po omega
      REAL(8) QmatDy(mesh%nbnpof,mesh%nnodes) ! grad u*_y po omega
      REAL(8) QmatDz(mesh%nbnpof,mesh%nnodes) ! grad u*_z po omega

      REAL(8) c,minedge
      REAL(8) ontH,ontG,ontDx,ontDy,ontDz ! ocena natancnosti integralov
      REAL(8) xi,eta,fi(mesh%npob)
      REAL(8) vsota
      REAL(8) mu,mumax,vH,vB,vD
      REAL(4) cptime,rcpu
      INTEGER rrcpu

      ontG=0.0D0
      ontH=0.0D0
      ontDx=0.0D0
      ontDy=0.0D0
      ontDz=0.0D0
      mumax=0.0D0


c     ker izracunam integrale pomnozene z int. funkcijami naenkrat rabim
      ALLOCATE (he(mesh%npob),ge(mesh%npof),be(mesh%npoc))
      ALLOCATE (abx(mesh%npob),aby(mesh%npob),abz(mesh%npob))
      ALLOCATE (heyz(mesh%npob),hezx(mesh%npob),hexy(mesh%npob)) ! teh ne rabim
      ALLOCATE (dxe(mesh%npoc),dye(mesh%npoc),dze(mesh%npoc))


      UmatH=0.0D0
      UmatG=0.0D0
      UmatB=0.0D0
      UmatAbdx=0.0D0
      UmatAbdy=0.0D0
      UmatAbdz=0.0D0
      UmatDx=0.0D0
      UmatDy=0.0D0
      UmatDz=0.0D0

      QmatH=0.0D0
      QmatG=0.0D0
      QmatB=0.0D0
      QmatAbdx=0.0D0
      QmatAbdy=0.0D0
      QmatAbdz=0.0D0
      QmatDx=0.0D0
      QmatDy=0.0D0
      QmatDz=0.0D0


c     Estimate time for integration
      rcpu=cptime(0.)


      CALL WarnErr(env,io,inp,5,"FMATsdmHfs","single-domain, mHfundSol, U source points",mesh%nnodes)

C     Zanka po U izvornih tockah
      DO it=1,mesh%nnodes

c       mu v izvorni tocki
        mu=SQRT(inp%beta/mHdiff(it))
        mumax=MAX(mu,mumax)

c       lokacija izvorne tocke
        xp=mesh%x(it,1)
        yp=mesh%x(it,2)
        zp=mesh%x(it,3)
C
C       INTEGRACIJA PO ROBU
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
C           SET SOURCE POINT :
            isrc=0
            DO isip=1,mesh%npob
              IF (it .EQ. mesh%ibc(je,isip)) isrc=isip
            END DO
C           poisce najkrajsi rob
            CALL FindMinEdgeOneElem(mesh,minedge,je)
C           integracija po robnih elementih celice
            CALL INTEBc9modHelm(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                   isrc,mesh%npob,mesh%npof,he,ge,abx,aby,abz,
     &                   heyz,hezx,hexy,minedge,mu)
c           zlozimo H-je po matrikah
            DO isip=1,mesh%npob
              j=mesh%lbn(mesh%ibc(je,isip)) ! stolpec
              UmatH(it,j)=UmatH(it,j)+he(isip)
              UmatAbdx(it,j)=UmatAbdx(it,j)+Abx(isip)
              UmatAbdy(it,j)=UmatAbdy(it,j)+Aby(isip)
              UmatAbdz(it,j)=UmatAbdz(it,j)+Abz(isip)
            END DO
c           zlozim G-je po matrikah
            DO isip=1,mesh%npof
c              j=mesh%ibcf(je,isip) ! stolpec
c             ko zlagam v pravokotno matriko, jih dam kar po vrsti
              j=(je-1)*mesh%npof+isip
              UmatG(it,j)=UmatG(it,j)+ge(isip)
            END DO
        END DO ! po robnih elementih
C
C       INTEGRACIJA PO OBMOCJU
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
C
C....     SET SOURCE POINT :
C
          isrc=0
          DO isip=1,mesh%npoc
            IF (it .EQ. mesh%idc(ic,isip)) isrc=isip
          END DO

          CALL INTEDc27modHelm(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                  x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,be,dxe,dye,dze,isrc,mesh%npoc,mu)
c         zlozimo vrstico po matrikah
          DO isip=1,mesh%npoc
            j=mesh%idc(ic,isip)
            UmatB(it,j)=UmatB(it,j)+be(isip)
            UmatDx(it,j)=UmatDx(it,j)+dxe(isip)
            UmatDy(it,j)=UmatDy(it,j)+dye(isip)
            UmatDz(it,j)=UmatDz(it,j)+dze(isip)
          END DO
        END DO ! po celicah

C       ****************************************************************
C                    Correct for singular integrals
C       ****************************************************************
        mu=SQRT(inp%beta/mHdiff(it))
c       izracunam prosti koeficient c in ga dam na diagonalo
c       predpostavim togi premik telesa u=1, q=0, f=-mu^2
        vH=0.0D0
        DO j=1,mesh%nbnodes
          vH=Vh+UmatH(it,j)
        END DO
        vB=0.0D0
        DO j=1,mesh%nnodes
          vB=vB+UmatB(it,j)
        END DO
        c = - vH + mu*mu*vB

        IF (mesh%lbn(it).GT.0) THEN ! je robno
          UmatH(it,mesh%lbn(it))=UmatH(it,mesh%lbn(it))+c
          IF (it.EQ.1) ontH=max(ontH,ABS(0.125D0-c)) ! prvo je vogal (nisem ziher)
        END IF

c       preverimo natancnost integralov
c       v=(z,0,0), w=(0,1,0)
        vD=0.0D0
        vH=0.0D0
        vB=0.0D0
        DO j=1,mesh%nbnodes
          vH=vH+UmatH(it,j)*mesh%x(mesh%gbn(j),3)
        END DO
        DO j=1,mesh%nnodes
          vB=vB+UmatB(it,j)*mesh%x(j,3)
          vD=vD+UmatDz(it,j)
        END DO
        vsota= Vh-mu*mu*Vb-vD
        IF (mesh%lbn(it).EQ.0) vsota = vsota + mesh%x(it,3) ! je obmocnom expl. dodam c*u (1*z)
        ontDz=max(ontDz,ABS(vsota))
        UmatDz(it,it)=UmatDz(it,it)+vsota  ! popravimo singularni integral  ??XX??

c       v=(0,x,0), w=(0,0,1)
        vD=0.0D0
        vH=0.0D0
        vB=0.0D0
        DO j=1,mesh%nbnodes
          vH=vH+UmatH(it,j)*mesh%x(mesh%gbn(j),1)
        END DO
        DO j=1,mesh%nnodes
          vB=vB+UmatB(it,j)*mesh%x(j,1)
          vD=vD+UmatDx(it,j)
        END DO
        vsota= Vh-mu*mu*Vb-vD
        IF (mesh%lbn(it).EQ.0) vsota = vsota + mesh%x(it,1) ! je obmocnom expl. dodam c*u (1*z)
        ontDx=max(ontDx,ABS(vsota))
        UmatDx(it,it)=UmatDx(it,it)+vsota  ! popravimo singularni integral  ??XX??

c       v=(0,0,y), w=(1,0,0)
        vD=0.0D0
        vH=0.0D0
        vB=0.0D0
        DO j=1,mesh%nbnodes
          vH=vH+UmatH(it,j)*mesh%x(mesh%gbn(j),2)
        END DO
        DO j=1,mesh%nnodes
          vB=vB+UmatB(it,j)*mesh%x(j,2)
          vD=vD+UmatDy(it,j)
        END DO
        vsota= Vh-mu*mu*Vb-vD
        IF (mesh%lbn(it).EQ.0) vsota = vsota + mesh%x(it,2) ! je obmocnom expl. dodam c*u (1*z)
        ontDy=max(ontDy,ABS(vsota))
        UmatDy(it,it)=UmatDy(it,it)+vsota  ! popravimo singularni integral  ??XX??

c         preverimo natancnost integralov
c         Predpostavim, da je mu=0 ali nesk. nisem ziher ampak dela
c         SMER X
c          vsota=0.0D0
c          DO j=1,mesh%nnodes
c            vsota=vsota-UmatDx(it,j)
c          END DO
c          DO j=1,mesh%nbnodes
c            vsota=vsota+UmatAbdx(it,j) ! to bi moralo biti nic
c          END DO
c          print *,vsota
c           if (it.eq.10) print *,vsota
c          print *,it,vsota,UmatDx(it,1),xp,yp,zp
c          ontDx=max(ontDx,ABS(vsota))
c          UmatDx(it,it)=UmatDx(it,it)+vsota
c         SMER Y
c          vsota=0.0D0
c          DO j=1,mesh%nnodes
c            vsota=vsota-UmatDy(it,j)
c          END DO
c          DO j=1,mesh%nbnodes
c            vsota=vsota+UmatAbdy(it,j)! to bi moralo biti nic
c          END DO
c          ontDy=max(ontDy,ABS(vsota))
c          UmatDy(it,it)=UmatDy(it,it)+vsota
c          print *,vsota
c         SMER Z
c          vsota=0.0D0
c          DO j=1,mesh%nnodes
c            vsota=vsota-UmatDz(it,j)
c          END DO
c          DO j=1,mesh%nbnodes
c            vsota=vsota+UmatAbdz(it,j)  ! to bi moralo biti nic
c          END DO
c          print *,vsota,"ww"
c          ontDz=max(ontDz,ABS(vsota))
c          UmatDz(it,it)=UmatDz(it,it)+vsota

c       Estimate time needed for integration
        IF (it.EQ.10) THEN
          rrcpu = INT(cptime(rcpu)*0.1*(mesh%nnodes+mesh%nbelem*mesh%npof))
          CALL WarnErr(env,io,inp,5,"FMATsdmHfs","single-domain, mHfundSol, integration time [s]",rrcpu)
        END IF


      END DO ! po izvornih tockah u


      CALL WarnErr(env,io,inp,5,"FMATsdmHfs","single-domain, mHfundSol, Q source points",mesh%nbelem*mesh%npof)

C     Zanka po Q izvornih tockah
      DO iti=1,mesh%nbelem
        DO itf=1,mesh%npof
          it=mesh%ibcf(iti,itf)
c         mu v izvorni tocki
          mu=SQRT(inp%beta/mHqDiff(it))
          mumax=MAX(mu,mumax)
c         lokacija izvorne tocke
          xp=mesh%xq(it,1)
          yp=mesh%xq(it,2)
          zp=mesh%xq(it,3)
c         vrstica v matriki
          i=(iti-1)*mesh%npof+itf
C
C         INTEGRACIJA PO ROBU
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
C           SET SOURCE POINT :
            isrc=0
            DO isip=1,mesh%npof
              IF (it .EQ. mesh%ibcf(je,isip)) isrc=isip+9 ! ker je tako v INTEBc9dc
            END DO
C           poisce najkrajsi rob
            CALL FindMinEdgeOneElem(mesh,minedge,je)
C           integracija po robnih elementih celice
            CALL INTEBc9modHelm(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                   isrc,mesh%npob,mesh%npof,he,ge,abx,aby,abz,
     &                   heyz,hezx,hexy,minedge,mu)
c           zlozimo H-je po matrikah
            DO isip=1,mesh%npob
              j=mesh%lbn(mesh%ibc(je,isip)) ! stolpec
              QmatH(i,j)=QmatH(i,j)+he(isip)
              QmatAbdx(i,j)=QmatAbdx(i,j)+Abx(isip)
              QmatAbdy(i,j)=QmatAbdy(i,j)+Aby(isip)
              QmatAbdz(i,j)=QmatAbdz(i,j)+Abz(isip)
            END DO
c           zlozim G-je po matrikah
            DO isip=1,mesh%npof
c              j=mesh%ibcf(je,isip) ! stolpec
c             ko zlagam v pravokotno matriko, jih dam kar po vrsti
              j=(je-1)*mesh%npof+isip
              QmatG(i,j)=QmatG(i,j)+ge(isip)
            END DO
          END DO ! po robnih elementih
C
C       INTEGRACIJA PO OBMOCJU
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
C
C....     SET SOURCE POINT :
C
          isrc=0
          DO je=1,mesh%nside
            DO isip=1,mesh%npof
              IF (it .EQ. mesh%ibf(ic,je,isip)) isrc=27+(je-1)*mesh%npof+isip
            END DO
          END DO

          CALL INTEDc27modHelm(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                  x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,be,dxe,dye,dze,isrc,mesh%npoc,mu)
c         zlozimo vrstico po matrikah
          DO isip=1,mesh%npoc
            j=mesh%idc(ic,isip)
            QmatB(i,j)=QmatB(i,j)+be(isip)
            QmatDx(i,j)=QmatDx(i,j)+dxe(isip)
            QmatDy(i,j)=QmatDy(i,j)+dye(isip)
            QmatDz(i,j)=QmatDz(i,j)+dze(isip)
          END DO
        END DO ! po celicah


c       izracunam prosti koeficient c in ga dam na diagonalo
        vH=0.0D0
        vB=0.0D0
        DO j=1,mesh%nbnodes
          vH=vH+QmatH(i,j)
        END DO
        DO j=1,mesh%nnodes
          vB=vB+QmatB(i,j)
        END DO
        c = - vH + mu*mu*vB
c       ker je fluks nezvezen, je izvorna tocka vedno na ploskvi, zato
c       vem, da je c=1/2. To uporabim za oceno natancnosti integracije
        ontG=max(ontG,ABS(0.5D0-c))

c         ker nimam neznanega u-ja v izvorni tocki za fluks, moram prispevek c
c         zinterpolirati na ostale
c         ugotovimo xi in eta za izvorno tocko
          CALL xieta34(xi,eta,itf)
c         izracunamo interpolacijske funkcije 9 tockovne zvezne v teh xi in eta
          CALL cshape9(fi,xi,eta,mesh%npob)
c         dodamo del C ja k vsaki izmed devetih tock
          DO isip=1,mesh%npob
            j=mesh%lbn(mesh%ibc(iti,isip)) ! stolpec
            QmatH(i,j)=QmatH(i,j)+c*fi(isip)
          END DO


        END DO
      END DO ! po izvornih tockah q


      DEALLOCATE (ge,he,be,abx,aby,abz,dxe,dye,dze,heyz,hezx,hexy)
      END



C -----------------------------------------------------------------------------
      SUBROUTINE corMin_SDkin_mHfs(env,io,inp,mesh,gauss,mHdiff,mHqDiff,
     &               UmatH,UmatG,UmatB,UmatBp,UmatAbdx,UmatAbdy,UmatAbdz,UmatDx,UmatDy,UmatDz,
     &               QmatH,QmatG,QmatB,QmatBp,QmatAbdx,QmatAbdy,QmatAbdz,QmatDx,QmatDy,QmatDz)
C
C     $: Compute or Read Matrices, modified Helmholtz fundamental solution
C        Single domain
C        (H, G, Adx, Ady, Adz, B, Dx, Dy, Dz)
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(IOtype) io
      TYPE(inputtype) inp
      TYPE(penv) :: env

      REAL(8) mHdiff(mesh%nnodes),mHqDiff(mesh%nq)

      REAL(8) UmatH(mesh%nnodes,mesh%nbnodes)   ! grad u* skalarno normala po gamma
      REAL(8) UmatG(mesh%nnodes,mesh%nbnpof)  ! u* po gamma
      REAL(8) UmatB(mesh%nnodes,mesh%nnodes)   ! u* po omega
      REAL(8) UmatBp(mesh%nnodes,mesh%nnodes)   ! u* po omega
      REAL(8) UmatAbdx(mesh%nnodes,mesh%nbnodes) ! u* krat nx po gamma
      REAL(8) UmatAbdy(mesh%nnodes,mesh%nbnodes) ! u* krat ny po gamma
      REAL(8) UmatAbdz(mesh%nnodes,mesh%nbnodes) ! u* krat nz po gamma
      REAL(8) UmatDx(mesh%nnodes,mesh%nnodes) ! grad u*_x po omega
      REAL(8) UmatDy(mesh%nnodes,mesh%nnodes) ! grad u*_y po omega
      REAL(8) UmatDz(mesh%nnodes,mesh%nnodes) ! grad u*_z po omega

      REAL(8) QmatH(mesh%nbnpof,mesh%nbnodes)   ! grad u* skalarno normala po gamma
      REAL(8) QmatG(mesh%nbnpof,mesh%nbnpof)  ! u* po gamma
      REAL(8) QmatB(mesh%nbnpof,mesh%nnodes)   ! u* po omega
      REAL(8) QmatBp(mesh%nbnpof,mesh%nnodes)   ! u* po omega
      REAL(8) QmatAbdx(mesh%nbnpof,mesh%nbnodes) ! u* krat nx po gamma
      REAL(8) QmatAbdy(mesh%nbnpof,mesh%nbnodes) ! u* krat ny po gamma
      REAL(8) QmatAbdz(mesh%nbnpof,mesh%nbnodes) ! u* krat nz po gamma
      REAL(8) QmatDx(mesh%nbnpof,mesh%nnodes) ! grad u*_x po omega
      REAL(8) QmatDy(mesh%nbnpof,mesh%nnodes) ! grad u*_y po omega
      REAL(8) QmatDz(mesh%nbnpof,mesh%nnodes) ! grad u*_z po omega

      REAL(8) ontH,ontG,ontDx,ontDy,ontDz ! ocena natancnosti integralov
      REAL(8) mumax,s

      INTEGER iok

C     sum of diff as parameter for integral integritiy
      s = 0.0D0
      DO i=1,mesh%nnodes
        s = s + mHdiff(i)
      END DO

C     Check if integrals exist on file
      CALL CheckmHfullIntFile(io,inp%INTversion,mesh%nnodes,mesh%nbnodes,inp%beta,s,iok)

      IF (iok.EQ.1) THEN
c       Read integrals from disk
        CALL WarnErr(env,io,inp,0,"corMin_SDkin_mHfs","Reading single-domain mH integrals!",0)
        CALL ReadmHfullIntDisk(io,inp,mesh,
     &               UmatH,UmatG,UmatB,UmatBp,UmatAbdx,UmatAbdy,UmatAbdz,UmatDx,UmatDy,UmatDz,
     &               QmatH,QmatG,QmatB,QmatBp,QmatAbdx,QmatAbdy,QmatAbdz,QmatDx,QmatDy,QmatDz,
     &                             ontH,ontG,ontDx,ontDy,ontDz,mumax)

      ELSE
c       Calcualte integrals
        CALL WarnErr(env,io,inp,0,"corMin_SDkin_mHfs","Calculating single-domain mH integrals!",0)

        CALL FMATsdmHfs(env,io,inp,mesh,gauss,mHdiff,mHqDiff,
     &               UmatH,UmatG,UmatB,UmatAbdx,UmatAbdy,UmatAbdz,UmatDx,UmatDy,UmatDz,
     &               QmatH,QmatG,QmatB,QmatAbdx,QmatAbdy,QmatAbdz,QmatDx,QmatDy,QmatDz,
     &                          ontH,ontG,ontDx,ontDy,ontDz,mumax)
C       integral matrices needed for MATMUL version due to variable diff.
        CALL MakeBprime(mesh,inp,UmatB,UmatBp,QmatB,QmatBp,mHdiff,mhQdiff)


        IF (env%myproc.EQ.1) THEN ! to ni paralelizirano, samo en zapise
          CALL WritemHfullIntDisk(io,inp,mesh,
     &               UmatH,UmatG,UmatB,UmatBp,UmatAbdx,UmatAbdy,UmatAbdz,UmatDx,UmatDy,UmatDz,
     &               QmatH,QmatG,QmatB,QmatBp,QmatAbdx,QmatAbdy,QmatAbdz,QmatDx,QmatDy,QmatDz,
     &                             ontH,ontG,ontDx,ontDy,ontDz,mumax,inp%beta,s)
        END IF
      END IF

      IF (env%myproc.EQ.1) THEN
        WRITE(io%l,'(A)') "Ocena natancnosti mH single domain integralov!"
        WRITE(io%l,'(A,G15.10)') "H (u i.t.)= ",ontH
        WRITE(io%l,'(A,G15.10)') "H (q i.t.)= ",ontG
        WRITE(io%l,'(A,G15.10)') "Dx= ",ontDx
        WRITE(io%l,'(A,G15.10)') "Dy= ",ontDy
        WRITE(io%l,'(A,G15.10)') "Dz= ",ontDz
        WRITE(io%l,'(A,G15.10)') "max. mu= ",mumax
        WRITE(io%l,'(A)') ""
      END IF

c      print *,ontH,ontG,ontDx,ontDy,ontDz,mumax

      END


C -----------------------------------------------------------------------------
      SUBROUTINE ReadmHfullIntDisk(io,inp,mesh,
     &               UmatH,UmatG,UmatB,UmatBp,UmatAbdx,UmatAbdy,UmatAbdz,UmatDx,UmatDy,UmatDz,
     &               QmatH,QmatG,QmatB,QmatBp,QmatAbdx,QmatAbdy,QmatAbdz,QmatDx,QmatDy,QmatDz,
     &                             ontH,ontG,ontDx,ontDy,ontDz,mumax)
C
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(IOtype) io
      TYPE(inputtype) inp

      REAL(8) UmatH(mesh%nnodes,mesh%nbnodes)   ! grad u* skalarno normala po gamma
      REAL(8) UmatG(mesh%nnodes,mesh%nbnpof)  ! u* po gamma
      REAL(8) UmatB(mesh%nnodes,mesh%nnodes)   ! u* po omega
      REAL(8) UmatBp(mesh%nnodes,mesh%nnodes)   ! u* po omega
      REAL(8) UmatAbdx(mesh%nnodes,mesh%nbnodes) ! u* krat nx po gamma
      REAL(8) UmatAbdy(mesh%nnodes,mesh%nbnodes) ! u* krat ny po gamma
      REAL(8) UmatAbdz(mesh%nnodes,mesh%nbnodes) ! u* krat nz po gamma
      REAL(8) UmatDx(mesh%nnodes,mesh%nnodes) ! grad u*_x po omega
      REAL(8) UmatDy(mesh%nnodes,mesh%nnodes) ! grad u*_y po omega
      REAL(8) UmatDz(mesh%nnodes,mesh%nnodes) ! grad u*_z po omega

      REAL(8) QmatH(mesh%nbnpof,mesh%nbnodes)   ! grad u* skalarno normala po gamma
      REAL(8) QmatG(mesh%nbnpof,mesh%nbnpof)  ! u* po gamma
      REAL(8) QmatB(mesh%nbnpof,mesh%nnodes)   ! u* po omega
      REAL(8) QmatBp(mesh%nbnpof,mesh%nnodes)   ! u* po omega
      REAL(8) QmatAbdx(mesh%nbnpof,mesh%nbnodes) ! u* krat nx po gamma
      REAL(8) QmatAbdy(mesh%nbnpof,mesh%nbnodes) ! u* krat ny po gamma
      REAL(8) QmatAbdz(mesh%nbnpof,mesh%nbnodes) ! u* krat nz po gamma
      REAL(8) QmatDx(mesh%nbnpof,mesh%nnodes) ! grad u*_x po omega
      REAL(8) QmatDy(mesh%nbnpof,mesh%nnodes) ! grad u*_y po omega
      REAL(8) QmatDz(mesh%nbnpof,mesh%nnodes) ! grad u*_z po omega

      REAL(8) ontH,ontG,ontDx,ontDy,ontDz ! ocena natancnosti integralov
      REAL(8) mumax
      INTEGER a,b

      OPEN (io%imhf,FILE=TRIM(io%imhf_name),FORM='UNFORMATTED',STATUS='OLD')
      READ(io%imhf) a
      READ(io%imhf) a,b
      READ(io%imhf) mumax
      READ(io%imhf) mumax
      READ(io%imhf) mumax


      CALL RdSMat(UmatH,mesh%nnodes,mesh%nbnodes,io%imhf)   ! grad u* skalarno normala po gamma
      CALL RdSMat(UmatG,mesh%nnodes,mesh%nbnpof,io%imhf)  ! u* po gamma
      CALL RdSMat(UmatB,mesh%nnodes,mesh%nnodes,io%imhf)   ! u* po omega
      CALL RdSMat(UmatBp,mesh%nnodes,mesh%nnodes,io%imhf)   ! u* po omega
      CALL RdSMat(UmatAbdx,mesh%nnodes,mesh%nbnodes,io%imhf) ! u* krat nx po gamma
      CALL RdSMat(UmatAbdy,mesh%nnodes,mesh%nbnodes,io%imhf) ! u* krat ny po gamma
      CALL RdSMat(UmatAbdz,mesh%nnodes,mesh%nbnodes,io%imhf) ! u* krat nz po gamma
      CALL RdSMat(UmatDx,mesh%nnodes,mesh%nnodes,io%imhf) ! grad u*_x po omega
      CALL RdSMat(UmatDy,mesh%nnodes,mesh%nnodes,io%imhf) ! grad u*_y po omega
      CALL RdSMat(UmatDz,mesh%nnodes,mesh%nnodes,io%imhf) ! grad u*_z po omega

      CALL RdSMat(QmatH,mesh%nbnpof,mesh%nbnodes,io%imhf)   ! grad u* skalarno normala po gamma
      CALL RdSMat(QmatG,mesh%nbnpof,mesh%nbnpof,io%imhf)  ! u* po gamma
      CALL RdSMat(QmatB,mesh%nbnpof,mesh%nnodes,io%imhf)   ! u* po omega
      CALL RdSMat(QmatBp,mesh%nbnpof,mesh%nnodes,io%imhf)   ! u* po omega
      CALL RdSMat(QmatAbdx,mesh%nbnpof,mesh%nbnodes,io%imhf) ! u* krat nx po gamma
      CALL RdSMat(QmatAbdy,mesh%nbnpof,mesh%nbnodes,io%imhf) ! u* krat ny po gamma
      CALL RdSMat(QmatAbdz,mesh%nbnpof,mesh%nbnodes,io%imhf) ! u* krat nz po gamma
      CALL RdSMat(QmatDx,mesh%nbnpof,mesh%nnodes,io%imhf) ! grad u*_x po omega
      CALL RdSMat(QmatDy,mesh%nbnpof,mesh%nnodes,io%imhf) ! grad u*_y po omega
      CALL RdSMat(QmatDz,mesh%nbnpof,mesh%nnodes,io%imhf) ! grad u*_z po omega

      READ(io%imhf) ontH,ontG,ontDx,ontDy,ontDz

      CLOSE(io%imhf)

      END

C -----------------------------------------------------------------------------
      SUBROUTINE WritemHfullIntDisk(io,inp,mesh,
     &               UmatH,UmatG,UmatB,UmatBp,UmatAbdx,UmatAbdy,UmatAbdz,UmatDx,UmatDy,UmatDz,
     &               QmatH,QmatG,QmatB,QmatBp,QmatAbdx,QmatAbdy,QmatAbdz,QmatDx,QmatDy,QmatDz,
     &                             ontH,ontG,ontDx,ontDy,ontDz,mumax,b,s)
C
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(IOtype) io
      TYPE(inputtype) inp

      REAL(8) UmatH(mesh%nnodes,mesh%nbnodes)   ! grad u* skalarno normala po gamma
      REAL(8) UmatG(mesh%nnodes,mesh%nbnpof)  ! u* po gamma
      REAL(8) UmatB(mesh%nnodes,mesh%nnodes)   ! u* po omega
      REAL(8) UmatBp(mesh%nnodes,mesh%nnodes)   ! u* po omega
      REAL(8) UmatAbdx(mesh%nnodes,mesh%nbnodes) ! u* krat nx po gamma
      REAL(8) UmatAbdy(mesh%nnodes,mesh%nbnodes) ! u* krat ny po gamma
      REAL(8) UmatAbdz(mesh%nnodes,mesh%nbnodes) ! u* krat nz po gamma
      REAL(8) UmatDx(mesh%nnodes,mesh%nnodes) ! grad u*_x po omega
      REAL(8) UmatDy(mesh%nnodes,mesh%nnodes) ! grad u*_y po omega
      REAL(8) UmatDz(mesh%nnodes,mesh%nnodes) ! grad u*_z po omega

      REAL(8) QmatH(mesh%nbnpof,mesh%nbnodes)   ! grad u* skalarno normala po gamma
      REAL(8) QmatG(mesh%nbnpof,mesh%nbnpof)  ! u* po gamma
      REAL(8) QmatB(mesh%nbnpof,mesh%nnodes)   ! u* po omega
      REAL(8) QmatBp(mesh%nbnpof,mesh%nnodes)   ! u* po omega
      REAL(8) QmatAbdx(mesh%nbnpof,mesh%nbnodes) ! u* krat nx po gamma
      REAL(8) QmatAbdy(mesh%nbnpof,mesh%nbnodes) ! u* krat ny po gamma
      REAL(8) QmatAbdz(mesh%nbnpof,mesh%nbnodes) ! u* krat nz po gamma
      REAL(8) QmatDx(mesh%nbnpof,mesh%nnodes) ! grad u*_x po omega
      REAL(8) QmatDy(mesh%nbnpof,mesh%nnodes) ! grad u*_y po omega
      REAL(8) QmatDz(mesh%nbnpof,mesh%nnodes) ! grad u*_z po omega

      REAL(8) ontH,ontG,ontDx,ontDy,ontDz ! ocena natancnosti integralov
      REAL(8) mumax,b,s

      OPEN (io%imhf,FILE=TRIM(io%imhf_name),FORM='UNFORMATTED',STATUS='UNKNOWN')
      WRITE(io%imhf) inp%INTversion
      WRITE(io%imhf) mesh%nnodes,mesh%nbnodes
      WRITE(io%imhf) b
      WRITE(io%imhf) s
      WRITE(io%imhf) mumax


      CALL WrSMat(UmatH,mesh%nnodes,mesh%nbnodes,io%imhf)   ! grad u* skalarno normala po gamma
      CALL WrSMat(UmatG,mesh%nnodes,mesh%nbnpof,io%imhf)  ! u* po gamma
      CALL WrSMat(UmatB,mesh%nnodes,mesh%nnodes,io%imhf)   ! u* po omega
      CALL WrSMat(UmatBp,mesh%nnodes,mesh%nnodes,io%imhf)   ! u* po omega
      CALL WrSMat(UmatAbdx,mesh%nnodes,mesh%nbnodes,io%imhf) ! u* krat nx po gamma
      CALL WrSMat(UmatAbdy,mesh%nnodes,mesh%nbnodes,io%imhf) ! u* krat ny po gamma
      CALL WrSMat(UmatAbdz,mesh%nnodes,mesh%nbnodes,io%imhf) ! u* krat nz po gamma
      CALL WrSMat(UmatDx,mesh%nnodes,mesh%nnodes,io%imhf) ! grad u*_x po omega
      CALL WrSMat(UmatDy,mesh%nnodes,mesh%nnodes,io%imhf) ! grad u*_y po omega
      CALL WrSMat(UmatDz,mesh%nnodes,mesh%nnodes,io%imhf) ! grad u*_z po omega

      CALL WrSMat(QmatH,mesh%nbnpof,mesh%nbnodes,io%imhf)   ! grad u* skalarno normala po gamma
      CALL WrSMat(QmatG,mesh%nbnpof,mesh%nbnpof,io%imhf)  ! u* po gamma
      CALL WrSMat(QmatB,mesh%nbnpof,mesh%nnodes,io%imhf)   ! u* po omega
      CALL WrSMat(QmatBp,mesh%nbnpof,mesh%nnodes,io%imhf)   ! u* po omega
      CALL WrSMat(QmatAbdx,mesh%nbnpof,mesh%nbnodes,io%imhf) ! u* krat nx po gamma
      CALL WrSMat(QmatAbdy,mesh%nbnpof,mesh%nbnodes,io%imhf) ! u* krat ny po gamma
      CALL WrSMat(QmatAbdz,mesh%nbnpof,mesh%nbnodes,io%imhf) ! u* krat nz po gamma
      CALL WrSMat(QmatDx,mesh%nbnpof,mesh%nnodes,io%imhf) ! grad u*_x po omega
      CALL WrSMat(QmatDy,mesh%nbnpof,mesh%nnodes,io%imhf) ! grad u*_y po omega
      CALL WrSMat(QmatDz,mesh%nbnpof,mesh%nnodes,io%imhf) ! grad u*_z po omega

      WRITE(io%imhf) ontH,ontG,ontDx,ontDy,ontDz

      CLOSE(io%imhf)

      END

C______________________________________________________________________C
      SUBROUTINE CheckmHfullIntFile(io,a,b,c,d,e,iok)
      USE inc_types
      TYPE(IOtype) io
      INTEGER a,b,c
      INTEGER aa,bb,cc,iok
      REAL(8) d,e,dd,ee

      iok=0
      OPEN (io%imhf,FILE=TRIM(io%imhf_name),FORM='UNFORMATTED',STATUS='OLD',ERR=10)
      READ(io%imhf) aa
      READ(io%imhf) bb,cc
      READ(io%imhf) dd
      READ(io%imhf) ee

      IF (a.EQ.aa.AND.b.EQ.bb.AND.c.EQ.cc.AND.d.EQ.dd.AND.e.EQ.ee) iok=1
      CLOSE(io%imhf)

10    RETURN
      END



C -----------------------------------------------------------------------------
      SUBROUTINE FindMinEdgeOneElem(mesh,minedge,i)
C
C     $: Poisce najkrajsi rob elementa
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      REAL(8) minedge,elen
      INTEGER i

      minedge=1.0D20

        elen=SQRT((mesh%x(mesh%ibc(i,1),1)-mesh%x(mesh%ibc(i,3),1))**2+
     &            (mesh%x(mesh%ibc(i,1),2)-mesh%x(mesh%ibc(i,3),2))**2+
     &            (mesh%x(mesh%ibc(i,1),3)-mesh%x(mesh%ibc(i,3),3))**2)
        IF (elen.LT.minedge) minedge=elen

        elen=SQRT((mesh%x(mesh%ibc(i,3),1)-mesh%x(mesh%ibc(i,5),1))**2+
     &            (mesh%x(mesh%ibc(i,3),2)-mesh%x(mesh%ibc(i,5),2))**2+
     &            (mesh%x(mesh%ibc(i,3),3)-mesh%x(mesh%ibc(i,5),3))**2)
        IF (elen.LT.minedge) minedge=elen

        elen=SQRT((mesh%x(mesh%ibc(i,5),1)-mesh%x(mesh%ibc(i,7),1))**2+
     &            (mesh%x(mesh%ibc(i,5),2)-mesh%x(mesh%ibc(i,7),2))**2+
     &            (mesh%x(mesh%ibc(i,5),3)-mesh%x(mesh%ibc(i,7),3))**2)
        IF (elen.LT.minedge) minedge=elen

        elen=SQRT((mesh%x(mesh%ibc(i,7),1)-mesh%x(mesh%ibc(i,1),1))**2+
     &            (mesh%x(mesh%ibc(i,7),2)-mesh%x(mesh%ibc(i,1),2))**2+
     &            (mesh%x(mesh%ibc(i,7),3)-mesh%x(mesh%ibc(i,1),3))**2)
        IF (elen.LT.minedge) minedge=elen

      END

C     ------------------------------------------------------------------
      REAL FUNCTION GetFMmem(mesh)

C     ------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      REAL nnodes,nbnodes,nbnpof,mHfmnunk,mHnb

      nnodes=REAL(mesh%nnodes)
      nbnodes=REAL(mesh%nbnodes)
      nbnpof=REAL(mesh%nbnpof)
      mHfmnunk=REAL(mesh%mHfmnunk)
      mHnb=REAL(mesh%mHnb)

      GetFMmem = ( 5.0*nnodes*nnodes + nnodes*nbnpof + 4.0*nnodes*nbnodes
     &         +   5.0*nbnpof*nnodes + nbnpof*nbnpof + 4.0*nbnpof*nbnodes
     &         +   mHfmnunk*mHfmnunk + mHfmnunk*mHnb )
     &         * 8.0 / 1024.0 / 1024.0 / 1024.0

      END


C -----------------------------------------------------------------------------
      SUBROUTINE WaveletCompSmat(io,mName,fMat,nrow,ncol,crsMat,nnz,n,kappa)
C
C     $: Compresses matrices and produces CRS version of matrices
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(IOtype) io
      REAL(8) MB,comp,spar,norm,nnz,n,kappa
      INTEGER nrow,ncol
      CHARACTER(8) mName

      REAL(8) fMat(nrow,ncol)       ! full matrix
      TYPE(crsMatrixType) :: crsMat ! compressed matrix

      CALL Che_WcomMatrix(fMat,nrow,ncol,crsMat,kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A8,4G20.10)') mName,MB,comp,spar,norm
      nnz=nnz+spar*DBLE(nrow)*DBLE(ncol)
      n=n+DBLE(nrow)*DBLE(ncol)

      END


C -----------------------------------------------------------------------------
      SUBROUTINE AcaDo(fMat,nrow,ncol,acaMat,eps,norm,type)
C
C     $: Compresses matrices and produces CRS version of matrices
C
C -----------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      REAL(8) norm,eps
      INTEGER nrow,ncol,type

      REAL(8) fMat(nrow,ncol)       ! full matrix
      TYPE(acaMatrixType) :: acaMat ! compressed matrix

      IF (type.EQ.1) THEN ! fiksed rank
        CALL aca_CrossAprox(fMat,nrow,ncol,acaMat,eps,norm)
      ELSE IF (type.EQ.2) THEN ! Frobenuis norm rank
        CALL acaj_CrossAprox(fMat,nrow,ncol,acaMat,eps,norm)
      ELSE ! set rank based on vector multiplication
        CALL aca_RankCrossAprox(fMat,nrow,ncol,acaMat,eps,norm)
      END IF

      END


C -----------------------------------------------------------------------------
      SUBROUTINE AcaCompSmat(io,mName,fMat,nrow,ncol,acaMat,nnz,n,eps,type)
C
C     $: Compresses matrices and produces CRS version of matrices
C
C -----------------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(IOtype) io
      REAL(8) MB,comp,norm,nnz,n,eps,nACAe,nFULe
      INTEGER nrow,ncol,type
      INTEGER nz,i,j
      CHARACTER(8) mName

      REAL(8) fMat(nrow,ncol)       ! full matrix
      TYPE(acaMatrixType) :: acaMat ! compressed matrix

      CALL AcaDo(fMat,nrow,ncol,acaMat,eps,norm,type)

c      IF (type.EQ.1) THEN ! fiksed rank
c        CALL aca_CrossAprox(fMat,nrow,ncol,acaMat,eps,norm)
c      ELSE IF (type.EQ.2) THEN ! Frobenuis norm rank
c        CALL acaj_CrossAprox(fMat,nrow,ncol,acaMat,eps,norm)
c      ELSE ! set rank based on vector multiplication
c        CALL aca_RankCrossAprox(fMat,nrow,ncol,acaMat,eps,norm)
c      END IF

      nACAe = DBLE(acaMat%nrow)*DBLE(acaMat%rank)+DBLE(acaMat%rank)*DBLE(acaMat%ncol)
      nFULe = DBLE(acaMat%nrow)*DBLE(acaMat%ncol)

      comp = nACAe / nFULe
      MB = nACAe * 8.0 / 1024.0 / 1024

      WRITE(io%l,'(A8,4G20.10)') mName,MB,comp,norm
      nnz=nnz+nACAe
      n=n+nFULe

C     check number of zero elements (v CRS se splača če ničel več kot 1/3 vseh)
c            nz=0
c            DO j=1,acaMat%rank
c                DO i=1,acaMat%nrow
c                    IF (acaMat%a(i,j).EQ.0.0D0) nz=nz+1
c                END DO
c            END DO
c            print *,mName,nz,acaMat%rank*acaMat%nrow,DBLE(nz)/DBLE(acaMat%rank*acaMat%ncol)
c            nz=0
c            DO j=1,acaMat%ncol
c                DO i=1,acaMat%rank
c                    IF (acaMat%b(i,j).EQ.0.0D0) nz=nz+1
c                END DO
c            END DO
c            print *,mName,nz,acaMat%rank*acaMat%ncol,DBLE(nz)/DBLE(acaMat%rank*acaMat%ncol)

      END



C -----------------------------------------------------------------------------
      SUBROUTINE WaveletComp(env,io,inp,mesh,
     &               UmatB,UmatBp,UmatAbdx,UmatAbdy,UmatAbdz,UmatDx,UmatDy,UmatDz,
     &               QmatB,QmatBp,QmatAbdx,QmatAbdy,QmatAbdz,QmatDx,QmatDy,QmatDz,
     &               cUmatB,cUmatBp,cUmatAbdx,cUmatAbdy,cUmatAbdz,cUmatDx,cUmatDy,cUmatDz,
     &               cQmatB,cQmatBp,cQmatAbdx,cQmatAbdy,cQmatAbdz,cQmatDx,cQmatDy,cQmatDz)

C
C     $: Compresses matrices and produces CRS version of matrices
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(IOtype) io
      TYPE(inputtype) inp
      TYPE(penv) :: env

      REAL(8) MB,comp,spar,norm,nnz,n

      REAL(8) UmatB(mesh%nnodes,mesh%nnodes)   ! u* po omega
      REAL(8) UmatBp(mesh%nnodes,mesh%nnodes)   ! u* po omega
      REAL(8) UmatAbdx(mesh%nnodes,mesh%nbnodes) ! u* krat nx po gamma
      REAL(8) UmatAbdy(mesh%nnodes,mesh%nbnodes) ! u* krat ny po gamma
      REAL(8) UmatAbdz(mesh%nnodes,mesh%nbnodes) ! u* krat nz po gamma
      REAL(8) UmatDx(mesh%nnodes,mesh%nnodes) ! grad u*_x po omega
      REAL(8) UmatDy(mesh%nnodes,mesh%nnodes) ! grad u*_y po omega
      REAL(8) UmatDz(mesh%nnodes,mesh%nnodes) ! grad u*_z po omega

      REAL(8) QmatB(mesh%nbnpof,mesh%nnodes)   ! u* po omega
      REAL(8) QmatBp(mesh%nbnpof,mesh%nnodes)   ! u* po omega
      REAL(8) QmatAbdx(mesh%nbnpof,mesh%nbnodes) ! u* krat nx po gamma
      REAL(8) QmatAbdy(mesh%nbnpof,mesh%nbnodes) ! u* krat ny po gamma
      REAL(8) QmatAbdz(mesh%nbnpof,mesh%nbnodes) ! u* krat nz po gamma
      REAL(8) QmatDx(mesh%nbnpof,mesh%nnodes) ! grad u*_x po omega
      REAL(8) QmatDy(mesh%nbnpof,mesh%nnodes) ! grad u*_y po omega
      REAL(8) QmatDz(mesh%nbnpof,mesh%nnodes) ! grad u*_z po omega

      TYPE(crsMatrixType) :: cUmatB,cUmatBp,cUmatAbdX,cUmatAbdY,cUmatAbdZ,cUmatDx,cUmatDy,cUmatDz
      TYPE(crsMatrixType) :: cQmatB,cQmatBp,cQmatAbdX,cQmatAbdY,cQmatAbdZ,cQmatDx,cQmatDy,cQmatDz

      WRITE(io%l,'(A,G10.5)') "wavelet compression, kappa =",inp%WT_kappa
      WRITE(io%l,'(A)') "matrix,MB,comp,spar,norm"

      nnz=0.0D0
      nz=0.0D0

      CALL Che_WcomMatrix(UmatB,mesh%nnodes,mesh%nnodes,cUmatB,inp%WT_kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A,4G20.10)') "UmatB   ",MB,comp,spar,norm
      nnz=nnz+spar*DBLE(mesh%nnodes)*DBLE(mesh%nnodes)
      n=n+DBLE(mesh%nnodes)*DBLE(mesh%nnodes)
      CALL Che_WcomMatrix(UmatBp,mesh%nnodes,mesh%nnodes,cUmatBp,inp%WT_kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A,4G20.10)') "UmatBp  ",MB,comp,spar,norm
      nnz=nnz+spar*DBLE(mesh%nnodes)*DBLE(mesh%nnodes)
      n=n+DBLE(mesh%nnodes)*DBLE(mesh%nnodes)

      CALL Che_WcomMatrix(UmatAbdx,mesh%nnodes,mesh%nbnodes,cUmatAbdx,inp%WT_kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A,4G20.10)') "UmatAbdx",MB,comp,spar,norm
      nnz=nnz+spar*DBLE(mesh%nnodes)*DBLE(mesh%nbnodes)
      n=n+DBLE(mesh%nnodes)*DBLE(mesh%nbnodes)
      CALL Che_WcomMatrix(UmatAbdy,mesh%nnodes,mesh%nbnodes,cUmatAbdy,inp%WT_kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A,4G20.10)') "UmatAbdy",MB,comp,spar,norm
      nnz=nnz+spar*DBLE(mesh%nnodes)*DBLE(mesh%nbnodes)
      n=n+DBLE(mesh%nnodes)*DBLE(mesh%nbnodes)
      CALL Che_WcomMatrix(UmatAbdz,mesh%nnodes,mesh%nbnodes,cUmatAbdz,inp%WT_kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A,4G20.10)') "UmatAbdz",MB,comp,spar,norm
      nnz=nnz+spar*DBLE(mesh%nnodes)*DBLE(mesh%nbnodes)
      n=n+DBLE(mesh%nnodes)*DBLE(mesh%nbnodes)

      CALL Che_WcomMatrix(UmatDx,mesh%nnodes,mesh%nnodes,cUmatDx,inp%WT_kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A,4G20.10)') "UmatDx  ",MB,comp,spar,norm
      nnz=nnz+spar*DBLE(mesh%nnodes)*DBLE(mesh%nnodes)
      n=n+DBLE(mesh%nnodes)*DBLE(mesh%nnodes)
      CALL Che_WcomMatrix(UmatDy,mesh%nnodes,mesh%nnodes,cUmatDy,inp%WT_kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A,4G20.10)') "UmatDy  ",MB,comp,spar,norm
      nnz=nnz+spar*DBLE(mesh%nnodes)*DBLE(mesh%nnodes)
      n=n+DBLE(mesh%nnodes)*DBLE(mesh%nnodes)
      CALL Che_WcomMatrix(UmatDz,mesh%nnodes,mesh%nnodes,cUmatDz,inp%WT_kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A,4G20.10)') "UmatDx  ",MB,comp,spar,norm
      nnz=nnz+spar*DBLE(mesh%nnodes)*DBLE(mesh%nnodes)
      n=n+DBLE(mesh%nnodes)*DBLE(mesh%nnodes)

      CALL Che_WcomMatrix(QmatB,mesh%nbnpof,mesh%nnodes,cQmatB,inp%WT_kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A,4G20.10)') "QmatB   ",MB,comp,spar,norm
      nnz=nnz+spar*DBLE(mesh%nbnpof)*DBLE(mesh%nnodes)
      n=n+DBLE(mesh%nbnpof)*DBLE(mesh%nnodes)
      CALL Che_WcomMatrix(QmatBp,mesh%nbnpof,mesh%nnodes,cQmatBp,inp%WT_kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A,4G20.10)') "QmatBp  ",MB,comp,spar,norm
      nnz=nnz+spar*DBLE(mesh%nbnpof)*DBLE(mesh%nnodes)
      n=n+DBLE(mesh%nbnpof)*DBLE(mesh%nnodes)


      CALL Che_WcomMatrix(QmatAbdx,mesh%nbnpof,mesh%nbnodes,cQmatAbdx,inp%WT_kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A,4G20.10)') "QmatAbdx",MB,comp,spar,norm
      nnz=nnz+spar*DBLE(mesh%nbnpof)*DBLE(mesh%nbnodes)
      n=n+DBLE(mesh%nbnpof)*DBLE(mesh%nbnodes)
      CALL Che_WcomMatrix(QmatAbdy,mesh%nbnpof,mesh%nbnodes,cQmatAbdy,inp%WT_kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A,4G20.10)') "QmatAbdy",MB,comp,spar,norm
      nnz=nnz+spar*DBLE(mesh%nbnpof)*DBLE(mesh%nbnodes)
      n=n+DBLE(mesh%nbnpof)*DBLE(mesh%nbnodes)
      CALL Che_WcomMatrix(QmatAbdz,mesh%nbnpof,mesh%nbnodes,cQmatAbdz,inp%WT_kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A,4G20.10)') "QmatAbdz",MB,comp,spar,norm
      nnz=nnz+spar*DBLE(mesh%nbnpof)*DBLE(mesh%nbnodes)
      n=n+DBLE(mesh%nbnpof)*DBLE(mesh%nbnodes)

      CALL Che_WcomMatrix(QmatDx,mesh%nbnpof,mesh%nnodes,cQmatDx,inp%WT_kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A,4G20.10)') "QmatDx  ",MB,comp,spar,norm
      nnz=nnz+spar*DBLE(mesh%nbnpof)*DBLE(mesh%nnodes)
      n=n+DBLE(mesh%nbnpof)*DBLE(mesh%nnodes)
      CALL Che_WcomMatrix(QmatDy,mesh%nbnpof,mesh%nnodes,cQmatDy,inp%WT_kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A,4G20.10)') "QmatDy  ",MB,comp,spar,norm
      nnz=nnz+spar*DBLE(mesh%nbnpof)*DBLE(mesh%nnodes)
      n=n+DBLE(mesh%nbnpof)*DBLE(mesh%nnodes)
      CALL Che_WcomMatrix(QmatDz,mesh%nbnpof,mesh%nnodes,cQmatDz,inp%WT_kappa,MB,comp,spar,norm)
      WRITE(io%l,'(A,4G20.10)') "QmatDx  ",MB,comp,spar,norm
      nnz=nnz+spar*DBLE(mesh%nbnpof)*DBLE(mesh%nnodes)
      n=n+DBLE(mesh%nbnpof)*DBLE(mesh%nnodes)

      WRITE(io%l,'(A,G20.10)') "TOTAL matrix element ratio (compressed / full) =  ",nnz/n

      FLUSH(io%l)

      END
