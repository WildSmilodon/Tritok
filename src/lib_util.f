C     ------------------------------------------------------------------
      SUBROUTINE SetTimeScheme(env,io,inp)
C
C     Set finite difference time scheme approximation
C
C     ------------------------------------------------------------------
      USE inc_types
      IMPLICIT NONE

      TYPE(inputtype) :: inp
      TYPE(IOtype)    :: io
      TYPE(penv) :: env

      IF (inp%TimeScheme.EQ.1) THEN   ! du/dt
C
C       first order derivative over time
C
        inp%beta =+1.50D00/inp%tstep
        inp%beta2=+2.00D00/inp%tstep
        inp%beta3=+0.50D00/inp%tstep
      ELSE IF (inp%TimeScheme.EQ.2) THEN  ! d^2u/dt^2
C
C       second order derivative over time
C
        inp%beta =+1.00D00/inp%tstep**2.0D0
        inp%beta2=+2.00D00/inp%tstep**2.0D0
        inp%beta3=+1.00D00/inp%tstep**2.0D0
      ELSE
        CALL WarnErr(env,io,inp,2,'SetTimeScheme','Unknown time scheme!',0)
      END IF

      END

C     ------------------------------------------------------------------
      SUBROUTINE getRMSnorm(vec1,vec2,n,norm)
C
C     Izracuna RMS normo med vec1 in vec2 relativno na vec1
C
C     ------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n,i
      REAL(8) s,s2,vec1(n),vec2(n),norm
      s=0.0D0
      s2=0.0D0

      DO i=1,n
        s=s+(vec1(i)-vec2(i))**2
        s2=s2+vec1(i)**2
      ENDDO
      IF (s2>0.0D0) THEN
        norm = SQRT(s/s2)
      ELSE
        norm = SQRT(s)
      END IF

      END


C______________________________________________________________________C
      SUBROUTINE InterWallFluxtoFun(mesh,q,u,walln)
C
C     Interpolira iz fluks tock v funkcijske tocke, samo za steno WALLN
C
C______________________________________________________________________C
      USE inc_types

      TYPE(meshType) mesh

      REAL(8) q(mesh%nq)
      REAL(8) u(mesh%nnodes)

      REAL(8) xi(9),eta(9),fi(4)
      INTEGER, ALLOCATABLE :: ndata(:)

      INTEGER i,j,k,node,walln

       xi(1)=-1.0D0
      eta(1)=-1.0D0

       xi(2)= 0.0D0
      eta(2)=-1.0D0

       xi(3)=+1.0D0
      eta(3)=-1.0D0

       xi(4)=+1.0D0
      eta(4)= 0.0D0

       xi(5)=+1.0D0
      eta(5)=+1.0D0

       xi(6)= 0.0D0
      eta(6)=+1.0D0

       xi(7)=-1.0D0
      eta(7)=+1.0D0

       xi(8)=-1.0D0
      eta(8)= 0.0D0

       xi(9)= 0.0D0
      eta(9)= 0.0D0



      ALLOCATE (ndata(mesh%nnodes))
      ndata=0
      u=0.0D0

      DO i=1,mesh%nbelem
        IF (mesh%bewn(i).EQ.walln) THEN ! sem na pravi steni
          DO j=1,mesh%npob
            node=mesh%ibc(i,j)
            ndata(node)=ndata(node)+1
            CALL dl34shape4(fi,xi(j),eta(j),4)
            DO k=1,4
              u(node)=u(node)+fi(k)*q(mesh%ibcf(i,k))
            END DO
          END DO
        END IF
      END DO

c     izracunam povprecje
      DO i=1,mesh%nnodes
        IF (ndata(i).GT.0) u(i)=u(i)/DBLE(ndata(i))
      END DO

      DEALLOCATE (ndata)

      END






C ********************************************************************
      SUBROUTINE GetRotationMatrix(R,RT,n1)
C     Creates a rotation matrix R, so that
C     vector n1 is rotated into e1=(1,0,0)
C     e1=R*n1
C     n1=RT*e1
      REAL(8) e1(3),e2(3),n1(3),n2(3),tmp(3),dp
      REAL(8) e3(3),e4(3),n3(3),n4(3)
      REAL(8) Re(3,3),RnT(3,3),R(3,3),RT(3,3)
C
      e1(1)=1.0D0
      e1(2)=0.0D0
      e1(3)=0.0D0
C
      e2(1)=0.0D0
      e2(2)=1.0D0
      e2(3)=0.0D0

C     Create a vector n2 perpendicular to vector n1
C     tmp is a vector, which is not colilear with n1
      tmp(1)=n1(2)+3.345345
      tmp(2)=n1(3)-7.345344
      tmp(3)=n1(1)+6.243234

      CALL CrossProduct(n2,n1,tmp)
      CALL NormVector(n2)
C     Dot product of perpendicular vectors should be zero
      CALL DotProduct(n1,n2,dp)
      IF (ABS(dp).GT.1.0D-14) THEN
        PRINT *,"ERROR IN GetRotationMatrix"
        STOP
      END IF
C     Use TRIAD algorithm http://en.wikipedia.org/wiki/Triad_Method
      CALL CrossProduct(n3,n1,n2)
      CALL NormVector(n3)
      CALL CrossProduct(n4,n1,n3)
      CALL NormVector(n4)
      CALL CrossProduct(e3,e1,e2)
      CALL NormVector(e3)
      CALL CrossProduct(e4,e1,e3)
      CALL NormVector(e4)

      RnT(1,1)=n1(1)
      RnT(1,2)=n1(2)
      RnT(1,3)=n1(3)

      RnT(2,1)=n3(1)
      RnT(2,2)=n3(2)
      RnT(2,3)=n3(3)

      RnT(3,1)=n4(1)
      RnT(3,2)=n4(2)
      RnT(3,3)=n4(3)

      Re(1,1)=e1(1)
      Re(2,1)=e1(2)
      Re(3,1)=e1(3)

      Re(1,2)=e3(1)
      Re(2,2)=e3(2)
      Re(3,2)=e3(3)

      Re(1,3)=e4(1)
      Re(2,3)=e4(2)
      Re(3,3)=e4(3)

C     Rotational matrix
      R=MATMUL(Re,RnT)

C     Transpose rotation matrix
      RT(1,1)=R(1,1)
      RT(1,2)=R(2,1)
      RT(1,3)=R(3,1)
      RT(2,1)=R(1,2)
      RT(2,2)=R(2,2)
      RT(2,3)=R(3,2)
      RT(3,1)=R(1,3)
      RT(3,2)=R(2,3)
      RT(3,3)=R(3,3)

c      print *,R
c      print *,RT

      END


C ********************************************************************
      SUBROUTINE SetCrossProduct(mesh,a,b,c)
C     a=b x c
      USE inc_types
      TYPE(meshType) :: mesh
      INTEGER i
      REAL(8) a(mesh%nnodes,3),b(mesh%nnodes,3),c(mesh%nnodes,3)

      DO i=1,mesh%nnodes
        a(i,1)= b(i,2)*c(i,3)-b(i,3)*c(i,2)
        a(i,2)=-b(i,1)*c(i,3)+b(i,3)*c(i,1)
        a(i,3)= b(i,1)*c(i,2)-b(i,2)*c(i,1)
      END DO

      END

C ********************************************************************
      SUBROUTINE CrossProduct(a,b,c)
C     a=b x c
      REAL(8) a(3),b(3),c(3)

      a(1)= b(2)*c(3)-b(3)*c(2)
      a(2)=-b(1)*c(3)+b(3)*c(1)
      a(3)= b(1)*c(2)-b(2)*c(1)

      END

C ********************************************************************
      SUBROUTINE DotProduct(a,b,dp)
C     dp = a \cdot b
      REAL(8) a(3),b(3),dp

      dp=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

      END

C ********************************************************************
      SUBROUTINE SetDotProduct(mesh,a,b,dp)
C     dp = a \cdot b
      USE inc_types
      TYPE(meshType) :: mesh
      INTEGER i
      REAL(8) a(mesh%nnodes,3),b(mesh%nnodes,3),dp(mesh%nnodes)

      DO i=1,mesh%nnodes
        dp(i)=a(i,1)*b(i,1)+a(i,2)*b(i,2)+a(i,3)*b(i,3)
      END DO

      END



C ********************************************************************
      SUBROUTINE NormVector(a)
C     a=a / norm(a)
      REAL(8) a(3),norm
      INTEGER i

      DO i=1,3
        norm=norm+a(i)**2
      END DO
      norm=SQRT(norm)
      DO i=1,3
        a(i)=a(i)/norm
      END DO

      END


C______________________________________________________________________C
      SUBROUTINE Channel3Danal1point(yy,zz,velocity,Ly,Lz,Cy,Cz)
      REAL(8) velocity

      REAL(8) Ly,Lz,up,mu,pi,Cy,Cz
      INTEGER n,nn,i
      REAL(8) val,suma,sumao,dpdx
      REAL(8) y,z,h,yy,zz

      up=1.0D0
      mu=1.0D0
      PI=2.*ASIN(1.0D0)
      suma=0.0D0
      do n=1,100
        nn=2*n-1
        suma=suma+Tanh(DBLE(nn)*Pi*ly*0.5D0/lz)/DBLE(nn)**5
      end do
      dpdx=12.0D0*mu*up/lz**2/(1.0D0-192.0D0/Pi**5*Lz/Ly*suma)

c v tri.bic podaj (na primer)
c 5 1 3
c #Ly, Lz, Cy, Cz
c 4 16 2 0
c Za kvadrat s stranico 1 in centrom 0.5,0.5
c 1 1 0.5 0.5

C
C     formula je za kanal, kjer je 0,0 v sredini, Lz/2,Ly/2 pa na vogalu
C     mreza naj bo v ravnini y-z med 0,0 in Ly,Lz - v x smeri pa je poljubna
C
c      y=yy-0.5D0*Ly
c      z=zz-0.5D0*Lz

C     Cy in Cz je center pravokotnika, Ly in Lz sta širini
      y=yy-Cy
      z=zz-Cz


      suma=0.0D0
      n=0

10    CONTINUE
        n=n+1
        nn=2*n-1
        sumao=suma
        suma=suma+(-1.0D0)**((nn-1)/2)/nn**3*COS(nn*Pi*z/Lz)*(1.0D0-COSH(nn*Pi*y/Lz)/COSH(0.5D0*nn*Pi*Ly/Lz))
      IF (ABS(suma-sumao).GT.1.0D-6) GOTO 10

      velocity=4.0D0*Lz**2/mu/Pi**3*dpdx*suma


      END



C______________________________________________________________________C
      SUBROUTINE Channel3Danal(mesh,velocity)
C
c  @ARTICLE{che04,
c  author = {C. S. Chen},
c  title = {Numerical method for predicting three-dimensional steady
compressible
c    flow in long microchannels},
c  journal = {J. Micromech. Microeng.},
c  year = {2004},
c  volume = {14},
c  pages = {1091--1100},
C
C

      USE inc_types

      TYPE(meshType) :: mesh
      REAL(8) velocity(mesh%nnodes,3)

      REAL(8) Ly,Lz,up,mu,pi
      INTEGER n,nn,i
      REAL(8) val,suma,sumao,dpdx
      REAL(8) y,z,h

      Ly=1.0D0
      Lz=1.0D0
      up=1.0D0
      mu=1.0D0
      PI=2.*ASIN(1.0D0)
      velocity=0.0D0

      suma=0.0D0
      do n=1,100
        nn=2*n-1
        suma=suma+Tanh(DBLE(nn)*Pi*ly*0.5D0/lz)/DBLE(nn)**5
      end do
      dpdx=12.0D0*mu*up/lz**2/(1.0D0-192.0D0/Pi**5*Lz/Ly*suma)

C
C     loop through the mesh
C
      DO i=1,mesh%nnodes
C
C     formula je za kanal, kjer je 0,0 v sredini, Lz/2,Ly/2 pa na vogalu
C     mreza naj bo v ravnini y-z med 0,0 in Ly,Lz - v x smeri pa je poljubna
C
      y=mesh%x(i,2)-0.5D0*Ly
      z=mesh%x(i,3)-0.5D0*Lz

      suma=0.0D0
      n=0

10    CONTINUE
        n=n+1
        nn=2*n-1
        sumao=suma
        suma=suma+(-1.0D0)**((nn-1)/2)/nn**3*COS(nn*Pi*z/Lz)*(1.0D0-COSH(nn*Pi*y/Lz)/COSH(0.5D0*nn*Pi*Ly/Lz))
      IF (ABS(suma-sumao).GT.1.0D-6) GOTO 10

      velocity(i,1)=4.0D0*Lz**2/mu/Pi**3*dpdx*suma

      END DO

      END



C______________________________________________________________________C
      SUBROUTINE CopyMatrix(a,q)
C
C     $: Skopira strukturo vrednosti a matrike v q
C        uporablja blas rutine, saj je prepros q%v=a%v za visoke nnz
C        sesuje s segmentation fault.
C
C______________________________________________________________________C
      USE inc_types

      TYPE(matrix) :: a,q

      q%nnz=a%nnz
      q%neq=a%neq

      ALLOCATE(q%v(q%nnz))
      ALLOCATE(q%i(q%neq+1))
      ALLOCATE(q%j(q%nnz))
      ALLOCATE(q%d(q%neq))

      CALL ICOPY (q%neq+1, a%i, 1, q%i, 1)
      CALL ICOPY (q%nnz, a%j, 1, q%j, 1)
      CALL ICOPY (q%neq, a%d, 1, q%d, 1)
      CALL dcopy (q%nnz, a%v, 1, q%v,1)

      END

C______________________________________________________________________C
      SUBROUTINE DLSE(env,io,cnt,inp)
C
C     $: Nastavi epsilon za least squares solver
C
C______________________________________________________________________C
      USE inc_types   

      TYPE(penv) :: env          
      TYPE(InputType) inp
      TYPE(countType) cnt
      TYPE(IOtype) :: io      
      INTEGER i
      REAl(4) oldeps

c     eps za lsqr bo inp%dlse_r-krat manjsi od konvergence

      DO i=1,8
        oldeps=inp%dlse(i)
        IF (cnt%nlierr(i)/inp%dlse_r.LT.inp%dlse_epsth.AND.cnt%nlierr(i).NE.0.0D0) THEN
          inp%dlse(i)=cnt%nlierr(i)/inp%dlse_r
        ELSE
          inp%dlse(i)=inp%dlse_epsth
        END IF
        IF (oldeps.LT.inp%dlse(i))       inp%dlse(i)=oldeps       ! ne dovolim povecevanja
        IF (inp%dlse(i).LT.inp%lsqs_eps) inp%dlse(i)=inp%lsqs_eps ! ne pod spodnjo mejo
      END DO

      IF (env%myproc.EQ.1) THEN   
        WRITE (io%sol,'(I6,1x,I6,2x,G7.2,2x,G7.2,2x,G7.2,2x,G7.2,2x,G7.2,2x,G7.2,2x,G7.2,2x,G7.2,2x,G7.2)') 
     &         cnt%tstep,cnt%tsnnlit,(inp%dlse(i),i=1,9)
      END IF
      
      END


C______________________________________________________________________C
      SUBROUTINE WallFlux(lun,mesh,cnt,gauss,q,id)
C
C     $: Integrira q po stenah
C
C______________________________________________________________________C
      USE inc_types   
      
      TYPE(meshType) :: mesh      
      TYPE(gausstype) :: gauss 
      TYPE(countType) cnt            

      REAL(8) q(mesh%nq)
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8), ALLOCATABLE :: ge(:),flux(:),area(:),geb(:)
      CHARACTER(3) id
      CHARACTER(20) FMT
      
      INTEGER lun,i,j
      
      ALLOCATE (ge(mesh%npof),flux(mesh%nofw),area(mesh%nofw),geb(mesh%npob))
      
      flux=0.0D0
      area=0.0D0
      
      DO i=1,mesh%nbelem
          x1=mesh%x(mesh%ibc(i,1),1)
          y1=mesh%x(mesh%ibc(i,1),2)
          z1=mesh%x(mesh%ibc(i,1),3)
          x2=mesh%x(mesh%ibc(i,3),1)
          y2=mesh%x(mesh%ibc(i,3),2)
          z2=mesh%x(mesh%ibc(i,3),3)
          x3=mesh%x(mesh%ibc(i,5),1)
          y3=mesh%x(mesh%ibc(i,5),2)
          z3=mesh%x(mesh%ibc(i,5),3)
          x4=mesh%x(mesh%ibc(i,7),1)
          y4=mesh%x(mesh%ibc(i,7),2)
          z4=mesh%x(mesh%ibc(i,7),3)                   
          
          CALL INTEBqFlux(gauss,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,mesh%npof,ge)

          DO j=1,mesh%npof
            flux(mesh%bewn(i))=flux(mesh%bewn(i))+ge(j)*q(mesh%ibcf(i,j))
          END DO

          CALL INTEBfunc(gauss,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,mesh%npob,geb)

          DO j=1,mesh%npob
            area(mesh%bewn(i))=area(mesh%bewn(i))+geb(j)
          END DO

      END DO

      WRITE(FMT,'("(A3,1x,I6,", I0, "G16.8)")') mesh%nofw
      WRITE (lun,FMT) id,cnt%tstep,(flux(i),i=1,mesh%nofw)
      WRITE (lun,FMT) "ARE",cnt%tstep,(area(i),i=1,mesh%nofw)
      WRITE (lun,FMT) "q/A",cnt%tstep,(flux(i)/area(i),i=1,mesh%nofw)
      
      DEALLOCATE (ge,flux,geb,area)
      
      END


C______________________________________________________________________C
      SUBROUTINE GetWallFlux(mesh,cnt,gauss,q,flux)
C
C     $: Integrira q po stenah, vrne flux [W/m2]
C
C______________________________________________________________________C
      USE inc_types
      IMPLICIT NONE

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(countType) cnt

      REAL(8) q(mesh%nq),flux(mesh%nofw)
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8), ALLOCATABLE :: ge(:),area(:),geb(:)
      CHARACTER(3) id
      CHARACTER(20) FMT

      INTEGER i,j

      ALLOCATE (ge(mesh%npof),area(mesh%nofw),geb(mesh%npob))

      flux=0.0D0
      area=0.0D0

      DO i=1,mesh%nbelem
          x1=mesh%x(mesh%ibc(i,1),1)
          y1=mesh%x(mesh%ibc(i,1),2)
          z1=mesh%x(mesh%ibc(i,1),3)
          x2=mesh%x(mesh%ibc(i,3),1)
          y2=mesh%x(mesh%ibc(i,3),2)
          z2=mesh%x(mesh%ibc(i,3),3)
          x3=mesh%x(mesh%ibc(i,5),1)
          y3=mesh%x(mesh%ibc(i,5),2)
          z3=mesh%x(mesh%ibc(i,5),3)
          x4=mesh%x(mesh%ibc(i,7),1)
          y4=mesh%x(mesh%ibc(i,7),2)
          z4=mesh%x(mesh%ibc(i,7),3)

          CALL INTEBqFlux(gauss,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,mesh%npof,ge)

          DO j=1,mesh%npof
            flux(mesh%bewn(i))=flux(mesh%bewn(i))+ge(j)*q(mesh%ibcf(i,j))
          END DO

          CALL INTEBfunc(gauss,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,mesh%npob,geb)

          DO j=1,mesh%npob
            area(mesh%bewn(i))=area(mesh%bewn(i))+geb(j)
          END DO

      END DO

      DO i=1,mesh%nofw
        flux(i)=flux(i)/area(i)
      END DO

      DEALLOCATE (ge,geb,area)

      END


C______________________________________________________________________C
      SUBROUTINE WallFluxLambda(lun,mesh,cnt,gauss,q,id,lambda)
C
C     $: Integrira q po stenah
C
C______________________________________________________________________C
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(countType) cnt

      REAL(8) q(mesh%nq),lambda(mesh%nnodes),qlambda
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8), ALLOCATABLE :: ge(:),flux(:),ze(:),xi(:),eta(:)
      CHARACTER(3) id
      CHARACTER(20) FMT

      INTEGER lun,i,j,jj

      ALLOCATE (ge(mesh%npof),flux(mesh%nofw),ze(mesh%npob))
      ALLOCATE (xi(mesh%npof),eta(mesh%npof))

      flux=0.0D0

      CALL SetLC2Dnzv(xi,eta)

      DO i=1,mesh%nbelem
          x1=mesh%x(mesh%ibc(i,1),1)
          y1=mesh%x(mesh%ibc(i,1),2)
          z1=mesh%x(mesh%ibc(i,1),3)
          x2=mesh%x(mesh%ibc(i,3),1)
          y2=mesh%x(mesh%ibc(i,3),2)
          z2=mesh%x(mesh%ibc(i,3),3)
          x3=mesh%x(mesh%ibc(i,5),1)
          y3=mesh%x(mesh%ibc(i,5),2)
          z3=mesh%x(mesh%ibc(i,5),3)
          x4=mesh%x(mesh%ibc(i,7),1)
          y4=mesh%x(mesh%ibc(i,7),2)
          z4=mesh%x(mesh%ibc(i,7),3)

          CALL INTEBqFlux(gauss,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,mesh%npof,ge)


          DO j=1,mesh%npof ! po fluks tockah
            qlambda=0.0D0
c           dolocimo xi in eta glede na j
            CALL cshape9(ze,xi(j),eta(j),mesh%npob)
c           preinterpolira iz zveznih v nezvezne tocke
            DO jj=1,mesh%npob ! po zveznih tockah
              qlambda=qlambda+ze(jj)*lambda(mesh%ibc(i,jj))
            END DO
c           integriramo
            flux(mesh%bewn(i))=flux(mesh%bewn(i))+ge(j)*q(mesh%ibcf(i,j))*qlambda
c            print *,qlambda
          END DO
      END DO

      WRITE(FMT,'("(A3,1x,I6,", I0, "G16.8)")') mesh%nofw
      WRITE (lun,FMT) id,cnt%tstep,(flux(i),i=1,mesh%nofw)


      DEALLOCATE (ze,ge,flux,xi,eta)

      END



C______________________________________________________________________C
      SUBROUTINE WallLocalFlux(lun,mesh,cnt,gauss,q)
C
C     $: Integrira q po stenah
C
C______________________________________________________________________C
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(countType) cnt

      REAL(8) q(mesh%nq)
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8), ALLOCATABLE :: ge(:),flux(:)
      REAL(8), ALLOCATABLE :: geb(:),area(:)
      CHARACTER(3) id
      CHARACTER(20) FMT

      INTEGER lun,i,j,node

      ALLOCATE (ge(mesh%npof),flux(mesh%nbelem),geb(mesh%npob),area(mesh%nbelem))

      flux=0.0D0
      area=0.0D0

      DO i=1,mesh%nbelem
          x1=mesh%x(mesh%ibc(i,1),1)
          y1=mesh%x(mesh%ibc(i,1),2)
          z1=mesh%x(mesh%ibc(i,1),3)
          x2=mesh%x(mesh%ibc(i,3),1)
          y2=mesh%x(mesh%ibc(i,3),2)
          z2=mesh%x(mesh%ibc(i,3),3)
          x3=mesh%x(mesh%ibc(i,5),1)
          y3=mesh%x(mesh%ibc(i,5),2)
          z3=mesh%x(mesh%ibc(i,5),3)
          x4=mesh%x(mesh%ibc(i,7),1)
          y4=mesh%x(mesh%ibc(i,7),2)
          z4=mesh%x(mesh%ibc(i,7),3)

          CALL INTEBqFlux(gauss,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,mesh%npof,ge)

          DO j=1,mesh%npof
            flux(i)=flux(i)+ge(j)*q(mesh%ibcf(i,j))
          END DO

          CALL INTEBfunc(gauss,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,mesh%npob,geb)

          DO j=1,mesh%npob
            area(i)=area(i)+geb(j)
          END DO

      END DO


      DO i=1,mesh%nofw
        WRITE (lun,'(A,I3,A,I6,A1)') 'ZONE T="wall=',i," tstep=",cnt%tstep,'"'
        DO j=1,mesh%nbelem
          IF (mesh%bewn(j).EQ.i) THEN
            node=mesh%ibc(j,9)
            WRITE (lun,'(6G16.8)') mesh%x(node,1),mesh%x(node,2),mesh%x(node,3),flux(j),area(j),flux(j)/area(j)
          END IF
        END DO
      END DO

      DEALLOCATE (ge,flux,geb,area)


      END



C______________________________________________________________________C
      SUBROUTINE WallFuncInt(lun,mesh,cnt,gauss,f,id)
C
C     $: Integrira funkcijo po stenah
C
C______________________________________________________________________C
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(countType) cnt

      REAL(8) f(mesh%nnodes)
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8), ALLOCATABLE :: ge(:),func(:)
      CHARACTER(3) id
      CHARACTER(20) FMT

      INTEGER lun,i,j

      ALLOCATE (ge(mesh%npob),func(mesh%nofw))

      func=0.0D0

      DO i=1,mesh%nbelem
          x1=mesh%x(mesh%ibc(i,1),1)
          y1=mesh%x(mesh%ibc(i,1),2)
          z1=mesh%x(mesh%ibc(i,1),3)
          x2=mesh%x(mesh%ibc(i,3),1)
          y2=mesh%x(mesh%ibc(i,3),2)
          z2=mesh%x(mesh%ibc(i,3),3)
          x3=mesh%x(mesh%ibc(i,5),1)
          y3=mesh%x(mesh%ibc(i,5),2)
          z3=mesh%x(mesh%ibc(i,5),3)
          x4=mesh%x(mesh%ibc(i,7),1)
          y4=mesh%x(mesh%ibc(i,7),2)
          z4=mesh%x(mesh%ibc(i,7),3)

          CALL INTEBfunc(gauss,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,mesh%npob,ge)

          DO j=1,mesh%npob
            func(mesh%bewn(i))=func(mesh%bewn(i))+ge(j)*f(mesh%ibc(i,j))
          END DO
      END DO
      WRITE(FMT,'("(A3,1x,I6,", I0, "G16.8)")') mesh%nofw
      WRITE (lun,FMT) id,cnt%tstep,(func(i),i=1,mesh%nofw)

      DEALLOCATE (ge,func)

      END


C______________________________________________________________________C
      SUBROUTINE MassFlux(lun,mesh,cnt,gauss,velocity)
C
C     $: Integrira masni pretok po stenah
C
C______________________________________________________________________C
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(countType) cnt

      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8), ALLOCATABLE :: ge(:),func(:)
      CHARACTER*20 FMT

      INTEGER lun,i,j,k,kk

      ALLOCATE (ge(mesh%npob),func(mesh%nofw))

      func=0.0D0

      DO i=1,mesh%nbelem
          x1=mesh%x(mesh%ibc(i,1),1)
          y1=mesh%x(mesh%ibc(i,1),2)
          z1=mesh%x(mesh%ibc(i,1),3)
          x2=mesh%x(mesh%ibc(i,3),1)
          y2=mesh%x(mesh%ibc(i,3),2)
          z2=mesh%x(mesh%ibc(i,3),3)
          x3=mesh%x(mesh%ibc(i,5),1)
          y3=mesh%x(mesh%ibc(i,5),2)
          z3=mesh%x(mesh%ibc(i,5),3)
          x4=mesh%x(mesh%ibc(i,7),1)
          y4=mesh%x(mesh%ibc(i,7),2)
          z4=mesh%x(mesh%ibc(i,7),3)

          CALL INTEBfunc(gauss,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,mesh%npob,ge)

c         normalo imam za tocke na robu, rabim pa za robne elemente,
c         vzamem tisto iz sredinske tocke
          kk=mesh%lbn(mesh%ibc(i,9))

          DO j=1,mesh%npob
            k=mesh%ibc(i,j)
            func(mesh%bewn(i))=func(mesh%bewn(i))+ge(j)*
     &      (velocity(k,1)*mesh%nx(kk)+velocity(k,2)*mesh%ny(kk)+velocity(k,3)*mesh%nz(kk))
          END DO
      END DO

      WRITE(FMT,'("(I6,", I0, "G16.8)")') mesh%nofw
      WRITE (lun,FMT) cnt%tstep,(func(i),i=1,mesh%nofw)

      DEALLOCATE (ge,func)

      END


C______________________________________________________________________C
      SUBROUTINE GetMassFlux(mesh,cnt,gauss,velocity,flux)
C
C     $: Integrira masni pretok po stenah
C
C______________________________________________________________________C
      USE inc_types
      IMPLICIT NONE

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(countType) cnt

      REAL(8) velocity(mesh%nnodes,3),flux(mesh%nofw)
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8), ALLOCATABLE :: ge(:)
      CHARACTER*20 FMT

      INTEGER lun,i,j,k,kk

      ALLOCATE (ge(mesh%npob))

      flux=0.0D0

      DO i=1,mesh%nbelem
          x1=mesh%x(mesh%ibc(i,1),1)
          y1=mesh%x(mesh%ibc(i,1),2)
          z1=mesh%x(mesh%ibc(i,1),3)
          x2=mesh%x(mesh%ibc(i,3),1)
          y2=mesh%x(mesh%ibc(i,3),2)
          z2=mesh%x(mesh%ibc(i,3),3)
          x3=mesh%x(mesh%ibc(i,5),1)
          y3=mesh%x(mesh%ibc(i,5),2)
          z3=mesh%x(mesh%ibc(i,5),3)
          x4=mesh%x(mesh%ibc(i,7),1)
          y4=mesh%x(mesh%ibc(i,7),2)
          z4=mesh%x(mesh%ibc(i,7),3)

          CALL INTEBfunc(gauss,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,mesh%npob,ge)

c         normalo imam za tocke na robu, rabim pa za robne elemente,
c         vzamem tisto iz sredinske tocke
          kk=mesh%lbn(mesh%ibc(i,9))

          DO j=1,mesh%npob
            k=mesh%ibc(i,j)
            flux(mesh%bewn(i))=flux(mesh%bewn(i))+ge(j)*
     &      (velocity(k,1)*mesh%nx(kk)+velocity(k,2)*mesh%ny(kk)+velocity(k,3)*mesh%nz(kk))
          END DO
      END DO

      DEALLOCATE (ge)

      END



C______________________________________________________________________C
      SUBROUTINE ExportWallFlux(lun,mesh,cnt,q,qw,mHq)
C
C     $: Izpise q po stenah
C
C______________________________________________________________________C
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(countType) cnt

      REAL(8) q(mesh%nq),qw(mesh%nq,3),mHq(mesh%nq)
      CHARACTER*255 vrstica

      INTEGER lun,i,j,wall,node

      DO wall=1,mesh%nofw
        WRITE (lun,'(A,I2,A,I6,A1)') 'ZONE T="wall=',wall,' tstep=',cnt%tstep,'"'
        DO i=1,mesh%nbelem
          DO j=1,mesh%npof
            node=mesh%ibcf(i,j)
            IF (mesh%bewn(i).EQ.wall) THEN
              WRITE (vrstica,'(3G18.9,5G18.9)')
     &               mesh%xq(node,1),mesh%xq(node,2),mesh%xq(node,3),q(node),qw(node,1),qw(node,2),qw(node,3),mHq(node)
              CALL sqblnk(lun,vrstica)
            END IF
          END DO
        END DO
      END DO

      END



C______________________________________________________________________C
      SUBROUTINE FluxFM2macro(mesh,fmq,mq)
C
C     $: Skopira flux iz FM v macro
C
C______________________________________________________________________C
      USE inc_types

      TYPE(meshType) :: mesh
      REAL(8) mq(mesh%nq)
      REAL(8) fmq(mesh%nbnpof)

      INTEGER i,j,ii,dn

C     Neumann
      ii=0
      DO i=1,mesh%nbelem
        DO j=1,mesh%npof
          ii=ii+1
          dn=mesh%ibcf(i,j)
          mq(dn)=fmq(ii)
        END DO
      END DO

      END

C______________________________________________________________________C
      SUBROUTINE CheckControlFiles(env,io,inp,istop)
C
C     Get user requests from system
C
      USE inc_types   
      
      TYPE(InputType) :: inp
      TYPE(IOtype) :: io
      TYPE(penv) :: env      

      INTEGER istop,ierr
      LOGICAL lexist

      istop=0

      IF (env%myproc.EQ.1) THEN        
        INQUIRE (FILE='stop',EXIST=lexist)
        IF (lexist) THEN
          CALL WarnErr(env,io,inp,0,'CheckControlFiles','User requested writing of restart file.',0)
          CALL WarnErr(env,io,inp,0,'CheckControlFiles','User requested end of program.',0)
C         Delete the stop file
C          CALL SYSTEM("rm stop")
          OPEN  (UNIT=5, FILE="stop", STATUS="OLD")   ! Current directory
          CLOSE (UNIT=5, STATUS="DELETE")
          istop=1
        END IF      
      END IF

c     root (myproc=1) poslje istop vsem ostalim
      CALL MPI_BCAST(istop,1,MPI_INTEGER,0,env%comm,ierr)

      END



C______________________________________________________________________C
      SUBROUTINE ZoneAuxData(itp,tstep,rtime)
C
C     Tecplot write auxillay data for each zone                        C
C______________________________________________________________________C
      INTEGER itp,tstep,fnblc
      REAL*8 rtime
      CHARACTER*30 tmp

C v tecplotu naredi txt objekt in vanj napisi : t=&(AUXZONE[ACTIVEOFFSET=1]:rtime)s
C to potem izpise itime za vsako zono svoj.

      WRITE (tmp,*) tstep
      WRITE (itp,'(A,A,A)') 'AUXDATA tstep="',tmp(fnblc(tmp):len_trim(tmp)),'"'
      WRITE (tmp,'(G15.3)') rtime
      WRITE (itp,'(A,A,A)') 'AUXDATA rtime="',tmp(fnblc(tmp):len_trim(tmp)),'"'

      RETURN
      END
C______________________________________________________________________C
      INTEGER FUNCTION fnblc(niz)
C
C      returns the number of the first non blank character in a string
C______________________________________________________________________C
      CHARACTER*(*) niz
      INTEGER i
   
      DO i=1,Len_trim(niz)
        IF (niz(i:i).NE.' ') THEN
          fnblc=i
          RETURN
        END IF
      END DO

c     string containes only blanks
      fnblc=1
C
      RETURN
      END

C______________________________________________________________________C
      SUBROUTINE CheckConvergencevw(env,io,inp,mesh,cnt,v,vpnls,w,wpnls,T,Tpnls,dc,dcpnls,mHu,mHupnls)
C  
C     Calculates residual between nonlinear iterations
C______________________________________________________________________C
      USE inc_types
      TYPE(InputType) inp  
      TYPE(meshType) :: mesh  
      TYPE(IOtype) io      
      TYPE(countType) cnt  
      TYPE(penv) :: env    
      
      REAL(8) v(mesh%nnodes,3)
      REAL(8) vpnls(mesh%nnodes,3)      
      REAL(8) w(mesh%nnodes,3)
      REAL(8) wpnls(mesh%nnodes,3)            
      REAL(8) T(mesh%nnodes)
      REAL(8) Tpnls(mesh%nnodes)            
      REAL(8) DC(mesh%nnodes)
      REAL(8) DCpnls(mesh%nnodes)
      REAL(8) mHu(mesh%nnodes)
      REAL(8) mHupnls(mesh%nnodes)            
      REAL(8) dif(9),rndif(9),rnorm(9)
      REAL(8) maxnorm,eps
      
      DATA eps/1.0D-15/

      INTEGER i,j
      rndif=0.0D0
      rnorm=0.0D0
      cnt%nlierr=9.999999999D00
      
      DO i=1,mesh%nnodes
        dif(1)=v(i,1)-vpnls(i,1)
        dif(2)=v(i,2)-vpnls(i,2)
        dif(3)=v(i,3)-vpnls(i,3)
        dif(4)=w(i,1)-wpnls(i,1)
        dif(5)=w(i,2)-wpnls(i,2)
        dif(6)=w(i,3)-wpnls(i,3)                                        
        dif(7)=T(i)-Tpnls(i)
        dif(8)=DC(i)-DCpnls(i) 
        dif(9)=mHu(i)-mHupnls(i)                      

        DO j=1,9
          rndif(j)=rndif(j)+dif(j)**2
        END DO
        rnorm(1)=rnorm(1)+v(i,1)**2        
        rnorm(2)=rnorm(2)+v(i,2)**2        
        rnorm(3)=rnorm(3)+v(i,3)**2        
        rnorm(4)=rnorm(4)+w(i,1)**2        
        rnorm(5)=rnorm(5)+w(i,2)**2        
        rnorm(6)=rnorm(6)+w(i,3)**2                       
        rnorm(7)=rnorm(7)+T(i)**2
        rnorm(8)=rnorm(8)+DC(i)**2   
        rnorm(9)=rnorm(9)+mHu(i)**2                               

        vpnls(i,1)=v(i,1)
        vpnls(i,2)=v(i,2)
        vpnls(i,3)=v(i,3)                    
        wpnls(i,1)=w(i,1)
        wpnls(i,2)=w(i,2)
        wpnls(i,3)=w(i,3)                    
        Tpnls(i)=T(i)
        DCpnls(i)=DC(i)     
        mHupnls(i)=mHu(i)                            
      END DO

      maxnorm=MAXVAL(rnorm(1:3))
      IF (maxnorm.GT.eps) cnt%nlierr(1:3)=SQRT(rndif(1:3)/maxnorm)

      maxnorm=MAXVAL(rnorm(4:6))
      IF (maxnorm.GT.eps) cnt%nlierr(4:6)=SQRT(rndif(4:6)/maxnorm)

      maxnorm=rnorm(7)
      IF (maxnorm.GT.eps) cnt%nlierr(7)=SQRT(rndif(7)/maxnorm)

      maxnorm=rnorm(8)
      IF (maxnorm.GT.eps) cnt%nlierr(8)=SQRT(rndif(8)/maxnorm)

      maxnorm=rnorm(9)
      IF (maxnorm.GT.eps) cnt%nlierr(9)=SQRT(rndif(9)/maxnorm)
      
      DO i=1,9
        IF (rndif(i).EQ.0.AND.rnorm(i).EQ.0.0D0) cnt%nlierr(i)=0.0D0
        IF (cnt%nlierr(i).GT.9.999999999D00) cnt%nlierr(i)=9.999999999D00
      END DO
      
      IF (env%myproc.EQ.1) THEN
      
      IF (inp%scro.EQ.1) THEN
        WRITE (*,'(I6,1x,I6,2x,G10.5,2x,G10.5,2x,G10.5,2x,G10.5,2x,G10.5,2x,G10.5,2x,G10.5,2x,G10.5,2x,G10.5)') 
     &         cnt%tstep,cnt%tsnnlit,(cnt%nlierr(j),j=1,9)
      END IF

      WRITE (io%sta,'(I6,1x,I6,2x,G10.5,2x,G10.5,2x,G10.5,2x,G10.5,2x,G10.5,2x,G10.5,2x,G10.5,2x,G10.5,2x,G10.5)') 
     &         cnt%tstep,cnt%tsnnlit,(cnt%nlierr(j),j=1,9)
      CALL FLUSH(io%sta)

      WRITE (io%ite,'(I6,1x,I6,12I5)') 
     &         cnt%tstep,cnt%tsnnlit,(cnt%nit(j),j=1,12)

      END IF
      
      END

C______________________________________________________________________C
      SUBROUTINE CheckTIMEConvergence(env,io,mesh,cnt,w,wpnls)
C  
C     Calculates residual between nonlinear iterations
C______________________________________________________________________C
      USE inc_types
      TYPE(meshType) :: mesh  
      TYPE(IOtype) io      
      TYPE(countType) cnt      
      TYPE(penv) env
      
      REAL(8) w(mesh%nnodes,3)
      REAL(8) wpnls(mesh%nnodes,3)            
      REAL(8) dif(3),rndif(3),rnorm(3),erru(3)  
      REAL(8) maxnorm,eps
      
      DATA eps/1.0D-3/

      INTEGER i,j
      rndif=0.0D0
      rnorm=0.0D0
      erru=9.999999999D00
      
      DO i=1,mesh%nnodes
        dif(1)=w(i,1)-wpnls(i,1)
        dif(2)=w(i,2)-wpnls(i,2)
        dif(3)=w(i,3)-wpnls(i,3)                                        
        
        DO j=1,3
          rndif(j)=rndif(j)+dif(j)**2
        END DO
        rnorm(1)=rnorm(1)+w(i,1)**2        
        rnorm(2)=rnorm(2)+w(i,2)**2        
        rnorm(3)=rnorm(3)+w(i,3)**2                       

      END DO

      maxnorm=MAXVAL(rnorm(1:3))
      IF (maxnorm.GT.eps) erru(1:3)=SQRT(rndif(1:3)/maxnorm)
      continue
  
      IF (env%myproc.EQ.1) THEN
        WRITE (io%tim,'(I6,1x,I6,2x,G10.5,2x,G10.5,2x,G10.5)') 
     &         cnt%tstep,cnt%tsnnlit,(erru(j),j=1,3)
      END IF
      
      END



C______________________________________________________________________C
C______________________________________________________________________C
      SUBROUTINE rOneTL(lun,OneLine)
C     _    ___ _    _ 
C     Read One Text Line 
C
C______________________________________________________________________C
C     Returns the first nonempty text line in file LUN, which does not
C     include the # character. If end of file is encoutered, it returns EOF
      CHARACTER*(*) OneLine
      INTEGER lun,i
  
10    READ(lun,'(A)',END=20) OneLine  

C     Check if line is empty
      IF (len_trim(OneLine).EQ.0) GOTO 10

C     Check if line contains # character
      DO i=1,len_trim(OneLine)
        IF (OneLine(i:i).EQ.'#') GOTO 10
      ENDDO

      RETURN

20    OneLine='EOF'
      END     

C -----------------------------------------------------------------------------
      SUBROUTINE DatumInUra(cas)
C
C     $: writes date and time to cas
C
C -----------------------------------------------------------------------------            
      CHARACTER*50 D,T,Z
      INTEGER Value(8)
      CHARACTER*39 cas     
      
      CALL DATE_AND_TIME(D,T,Z,Value)
      WRITE (cas,33) 
     &    'Date and time = ',Value(3),'.',Value(2),'.',Value(1),' ',
     &    Value(5),':',Value(6),':',Value(7),'.',Value(8)
33    FORMAT (A16,I2,A1,I2,A1,I4,A1,I2,A1,I2,A1,I2,A1,I3)

      END      
      
C______________________________________________________________________C
      SUBROUTINE InitCpu(cpu,itim)
C
C     Init names for cpu
C
      USE inc_types   
      TYPE(CPUtype) cpu      
      INTEGER i,itim
      
      cpu%ncpu=34
      
      ALLOCATE (cpu%time(cpu%ncpu))
      ALLOCATE (cpu%desc(cpu%ncpu))      

      i=1
      cpu%desc(i)="x Bw r.h.s"
      i=i+1
      cpu%desc(i)="x Bw solve"      
      i=i+1
      cpu%desc(i)="x Dv r.h.s"
      i=i+1      
      cpu%desc(i)="x Dv solve"      
      i=i+1      
      cpu%desc(i)="x Dw r.h.s"
      i=i+1      
      cpu%desc(i)="x Dw solve"      
      i=i+1            
      cpu%desc(i)="y Bw r.h.s"
      i=i+1
      cpu%desc(i)="y Bw solve"      
      i=i+1
      cpu%desc(i)="y Dv r.h.s"
      i=i+1      
      cpu%desc(i)="y Dv solve"      
      i=i+1      
      cpu%desc(i)="y Dw r.h.s"
      i=i+1      
      cpu%desc(i)="y Dw solve"      
      i=i+1            
      cpu%desc(i)="z Bw r.h.s"
      i=i+1
      cpu%desc(i)="z Bw solve"      
      i=i+1
      cpu%desc(i)="z Dv r.h.s"
      i=i+1      
      cpu%desc(i)="z Dv solve"      
      i=i+1      
      cpu%desc(i)="z Dw r.h.s"
      i=i+1      
      cpu%desc(i)="z Dw solve"      
      i=i+1      
      cpu%desc(i)="DT r.h.s"
      i=i+1      
      cpu%desc(i)="DT solve"    
      i=i+1      
      cpu%desc(i)="DC r.h.s"
      i=i+1      
      cpu%desc(i)="DC solve" 
      i=i+1
      cpu%desc(i)="mH r.h.s"
      i=i+1      
      cpu%desc(i)="mH solve" 
      i=i+1
      cpu%desc(i)="particles"
      itim=i
      i=i+1      
      cpu%desc(i)="init, read input files"
      i=i+1      
      cpu%desc(i)="generate macro mesh"
      i=i+1      
      cpu%desc(i)="read and generate BiC"            
      i=i+1      
      cpu%desc(i)="integration single domain"      
      i=i+1      
      cpu%desc(i)="form sys, rhs matrix single domain"      
      i=i+1      
      cpu%desc(i)="integration makro domain"            
      i=i+1      
      cpu%desc(i)="form makro sys, rhs matrix"
      i=i+1      
      cpu%desc(i)="solve"       
      i=i+1      
      cpu%desc(i)="TOTAL"      
      
      cpu%time=0.          
      
      END             
      
C______________________________________________________________________C
      SUBROUTINE VmesniCas(cpu,itim)
C
C     Remember cpu time
C
      USE inc_types   
      TYPE(CPUtype) cpu     
      REAL cptime
      INTEGER itim 
      
      itim=itim+1
      cpu%time(itim)=cptime(cpu%t0)
      cpu%t0=cptime(0.)

      END      
C -----------------------------------------------------------------------------      
      SUBROUTINE CRSxV(mat,vek,DolVek,rez)
C
C     $: mnozi CRS matriko z vektorjem
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(matrix) :: mat
      INTEGER DolVek
      INTEGER i,j
      
      REAL(8) vek(DolVek)
      REAL(8) rez(mat%neq)

      DO i=1,mat%neq
        rez(i)=0.0D00
        DO j=mat%i(i),mat%i(i+1)-1
          rez(i)=rez(i)+mat%v(j)*vek(mat%j(j))
        END DO
      END DO
      
      END



c-------------------------------------------------------------------C
      SUBROUTINE Int2MacroNodes(mesh,diff,QDiff)
c....
      USE inc_types
      TYPE(meshType) :: mesh
c....
      REAL(8) xsi(mesh%npofc),eta(mesh%npofc),ceta(mesh%npofc)
      INTEGER ic,ii,jj,ll,kk
      REAL(8) fi(mesh%npoc)
      REAL(8) qdiff(mesh%nq)
      REAL(8) diff(mesh%nnodes)
c....
      xsi(1)=-0.750000000000000
      xsi(2)=-0.750000000000000
      xsi(3)=0.750000000000000
      xsi(4)=0.750000000000000
      xsi(5)=-0.750000000000000
      xsi(6)=0.750000000000000
      xsi(7)=0.750000000000000
      xsi(8)=-0.750000000000000
      xsi(9)=1.00000000000000
      xsi(10)=1.00000000000000
      xsi(11)=1.00000000000000
      xsi(12)=1.00000000000000
      xsi(13)=0.750000000000000
      xsi(14)=-0.750000000000000
      xsi(15)=-0.750000000000000
      xsi(16)=0.750000000000000
      xsi(17)=-1.00000000000000
      xsi(18)=-1.00000000000000
      xsi(19)=-1.00000000000000
      xsi(20)=-1.00000000000000
      xsi(21)=-0.750000000000000
      xsi(22)=0.750000000000000
      xsi(23)=0.750000000000000
      xsi(24)=-0.750000000000000
c....
      eta(1)=-0.750000000000000
      eta(2)=0.750000000000000
      eta(3)=0.750000000000000
      eta(4)=-0.750000000000000
      eta(5)=-1.00000000000000
      eta(6)=-1.00000000000000
      eta(7)=-1.00000000000000
      eta(8)=-1.00000000000000
      eta(9)=-0.750000000000000
      eta(10)=0.750000000000000
      eta(11)=0.750000000000000
      eta(12)=-0.750000000000000
      eta(13)=1.00000000000000
      eta(14)=1.00000000000000
      eta(15)=1.00000000000000
      eta(16)=1.00000000000000
      eta(17)=0.750000000000000
      eta(18)=-0.750000000000000
      eta(19)=-0.750000000000000
      eta(20)=0.750000000000000
      eta(21)=-0.750000000000000
      eta(22)=-0.750000000000000
      eta(23)=0.750000000000000
      eta(24)=0.750000000000000
c....
      ceta(1)=-1.00000000000000
      ceta(2)=-1.00000000000000
      ceta(3)=-1.00000000000000
      ceta(4)=-1.00000000000000
      ceta(5)=-0.750000000000000
      ceta(6)=-0.750000000000000
      ceta(7)=0.750000000000000
      ceta(8)=0.750000000000000
      ceta(9)=-0.750000000000000
      ceta(10)=-0.750000000000000
      ceta(11)=0.750000000000000
      ceta(12)=0.750000000000000
      ceta(13)=-0.750000000000000
      ceta(14)=-0.750000000000000
      ceta(15)=0.750000000000000
      ceta(16)=0.750000000000000
      ceta(17)=-0.750000000000000
      ceta(18)=-0.750000000000000
      ceta(19)=0.750000000000000
      ceta(20)=0.750000000000000
      ceta(21)=1.00000000000000
      ceta(22)=1.00000000000000
      ceta(23)=1.00000000000000
      ceta(24)=1.00000000000000
c
c....
      DO ic=1,mesh%nicell
        DO ii=1,mesh%nside
          DO jj=1,mesh%npof
            ll=(ii-1)*mesh%npof+jj
            CALL cshape27(xsi(ll),eta(ll),ceta(ll),fi,mesh%npoc)
            qdiff(mesh%ibf(ic,ii,jj))=0.0D0
            DO kk=1,27
              qdiff(mesh%ibf(ic,ii,jj))=qdiff(mesh%ibf(ic,ii,jj))+fi(kk)*diff(mesh%idc(ic,kk))
            END DO
          END DO
        END DO
      END DO


!      do ii=1,mesh%nq
!      print*,qdiff(ii)-mesh%xq(ii,1)
!      end do
!      stop

      END            
