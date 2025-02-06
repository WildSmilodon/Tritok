c
c
c     lib_3Dpt  - 3D particle tracking
c
c     spisal Jure Ravnik
c
c     List of  Routines
c                      - kks2lks  ! find xi,eta,zeta from x,y,z
c                      - lks2kks  ! find x,y,z from xi,eta,zeta
c                      - cipic    ! check if particle in cell
c
c     Version : 5. May 2008
c

C ------------------------------------------------------------------------------------    
      SUBROUTINE cipic(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                   x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,ierr)
C
C     $: izracuna normale na ploskve elementa in racuna skalarne 
C        produkte z vektorji do delca. Ce so vsi skalarni produkti
C        negativni, je delec v celici
C        heksaeder, 8 ogljisc
C
C        vrne ierr=0, ce je delec v elementu, in ierr=1 ce ga ni
C      
C ------------------------------------------------------------------------------------    
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8) x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x,y,z
      REAL(8) v1x,v1y,v1z,v2x,v2y,v2z,nx,ny,nz,r
      INTEGER ierr

      ierr=1


c     spredaj
      v1x=x2-x1
      v1y=y2-y1
      v1z=z2-z1

      v2x=x5-x1
      v2y=y5-y1
      v2z=z5-z1

      CALL ncrossp(v1x,v1y,v1z,v2x,v2y,v2z,nx,ny,nz)
      v2x=x-x1
      v2y=y-y1
      v2z=z-z1
      CALL dotp(nx,ny,nz,v2x,v2y,v2z,r)
      IF (r.GT.0.0D0) RETURN
        

c     desno
      v1x=x4-x2
      v1y=y4-y2
      v1z=z4-z2

      v2x=x6-x2
      v2y=y6-y2
      v2z=z6-z2

      CALL ncrossp(v1x,v1y,v1z,v2x,v2y,v2z,nx,ny,nz)
      v2x=x-x2
      v2y=y-y2
      v2z=z-z2
      CALL dotp(nx,ny,nz,v2x,v2y,v2z,r)
      IF (r.GT.0.0D0) RETURN
      
c     zadaj
      v1x=x3-x4
      v1y=y3-y4
      v1z=z3-z4

      v2x=x8-x4
      v2y=y8-y4
      v2z=z8-z4

      CALL ncrossp(v1x,v1y,v1z,v2x,v2y,v2z,nx,ny,nz)
      v2x=x-x4
      v2y=y-y4
      v2z=z-z4
      CALL dotp(nx,ny,nz,v2x,v2y,v2z,r)
      IF (r.GT.0.0D0) RETURN

c     levo
      v1x=x1-x3
      v1y=y1-y3
      v1z=z1-z3

      v2x=x7-x3
      v2y=y7-y3
      v2z=z7-z3

      CALL ncrossp(v1x,v1y,v1z,v2x,v2y,v2z,nx,ny,nz)
      v2x=x-x3
      v2y=y-y3
      v2z=z-z3
      CALL dotp(nx,ny,nz,v2x,v2y,v2z,r)
      IF (r.GT.0.0D0) RETURN

c     spodaj
      v1x=x2-x4
      v1y=y2-y4
      v1z=z2-z4

      v2x=x3-x4
      v2y=y3-y4
      v2z=z3-z4

      CALL ncrossp(v1x,v1y,v1z,v2x,v2y,v2z,nx,ny,nz)
      v2x=x-x4
      v2y=y-y4
      v2z=z-z4
      CALL dotp(nx,ny,nz,v2x,v2y,v2z,r)
      IF (r.GT.0.0D0) RETURN

c     zgoraj
      v1x=x8-x6
      v1y=y8-y6
      v1z=z8-z6

      v2x=x5-x6
      v2y=y5-y6
      v2z=z5-z6

      CALL ncrossp(v1x,v1y,v1z,v2x,v2y,v2z,nx,ny,nz)
      v2x=x-x6
      v2y=y-y6
      v2z=z-z6
      CALL dotp(nx,ny,nz,v2x,v2y,v2z,r)
      IF (r.GT.0.0D0) RETURN

c     delec je v elementu            
      ierr=0

      END
      

C -----------------------------------------------------------------------------
      SUBROUTINE ncrossp(v1x,v1y,v1z,v2x,v2y,v2z,rx,ry,rz)
C
C     $: NORMIRAN vektorski produkt
C
C -----------------------------------------------------------------------------            
      REAL(8) v1x,v1y,v1z,v2x,v2y,v2z,rx,ry,rz,ajac
      
      rx=v1y*v2z-v1z*v2y
      ry=-v1x*v2z+v1z*v2x
      rz=v1x*v2y-v1y*v2x

      ajac=SQRT(rx**2+ry**2+rz**2)
      rx=rx/ajac
      ry=ry/ajac
      rz=rz/ajac

       
      END

C -----------------------------------------------------------------------------
      SUBROUTINE dotp(v1x,v1y,v1z,v2x,v2y,v2z,r)
C
C     $: NORMIRAN vektorski produkt
C
C -----------------------------------------------------------------------------            
      REAL(8) v1x,v1y,v1z,v2x,v2y,v2z,r
      
      r=v1x*v2x+v1y*v2y+v1z*v2z
       
      END      
      
C ------------------------------------------------------------------------------------    
      SUBROUTINE kks2lks(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                   x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,xi,eta,zeta,ierr)

C
C     $: na podlagi x,y,z najde xi,eta,zeta znotraj 8ih vogalov
C        glej NR stran 372 poglavje 9-6 Newton-Rapson method
C
C        ierr=0, tocka (x,y,z) je v heksaedru, ki ga doloca 8 vozlisc
C        ierr=1, tocka (x,y,z) NI v heksaedru, ki ga doloca 8 vozlisc
C      
C ------------------------------------------------------------------------------------    
      REAL(8) xi,eta,zeta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8) x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x,y,z
      REAL(8) xx(3),eps
      INTEGER maxit,nit,ierr

c     zacetni priblizek
      xx=0.0D0
c     toleranca 
      eps=1.0D-10 !-6 !-13
c     najvecje stevilo iteracij
      maxit=1000 !10000
            
      CALL mnewt(maxit,xx,3,eps,eps,x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                         x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,nit)

c     zapomnimo rezultate, preverimo ce skonvergiralo !!
      IF (nit.GE.maxit) THEN
        WRITE (*,*) "Dosezeno maximalno stevilo iteracij, kks2lks",nit
      END IF

      xi=xx(1)
      eta=xx(2)
      zeta=xx(3)       

c     preverimo, ce je noter ali zunaj
      IF ( xi.GE.-1.0D0.AND.  xi.LE.1.0D0.AND.
     &    eta.GE.-1.0D0.AND. eta.LE.1.0D0.AND.
     &   zeta.GE.-1.0D0.AND.zeta.LE.1.0D0) THEN
        ierr=0
      ELSE 
        ierr=1 
      END IF

      END
      
      
C ------------------------------------------------------------------------------------    
      SUBROUTINE kkk2lks_userfun(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                           x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,xi,eta,zeta,F,J)

C
C     $: na podlagi x,y,z in priblizka xi,eta,zeta izracuna f, ki ga
C       minimiziramo in odvode jacobijeve 
C      
C ------------------------------------------------------------------------------------    
      REAL(8) xi,eta,zeta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8) x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x,y,z

      REAL(8) FIG(8),eta1m,eta1p,eta2m,eta2p,eta3m,eta3p
      
      REAL(8) F(3),J(3,3)

      ETA1M=1.0D0-xi
      ETA1P=1.0D0+xi   
      ETA2M=1.0D0-eta
      ETA2P=1.0D0+eta   
      ETA3M=1.0D0-zeta   
      ETA3P=1.0D0+zeta  

      FIG(1)=0.125D0*ETA1M*ETA2M*ETA3M
      FIG(2)=0.125D0*ETA1P*ETA2M*ETA3M
      FIG(3)=0.125D0*ETA1M*ETA2P*ETA3M
      FIG(4)=0.125D0*ETA1P*ETA2P*ETA3M
      FIG(5)=0.125D0*ETA1M*ETA2M*ETA3P
      FIG(6)=0.125D0*ETA1P*ETA2M*ETA3P
      FIG(7)=0.125D0*ETA1M*ETA2P*ETA3P
      FIG(8)=0.125D0*ETA1P*ETA2P*ETA3P

C     funkcije, ki jih minimiziramo
      F(1)=FIG(1)*X1+FIG(2)*X2+FIG(3)*X3+FIG(4)*X4+FIG(5)*X5+FIG(6)*X6+FIG(7)*X7+FIG(8)*X8-X
      F(2)=FIG(1)*Y1+FIG(2)*Y2+FIG(3)*Y3+FIG(4)*Y4+FIG(5)*Y5+FIG(6)*Y6+FIG(7)*Y7+FIG(8)*Y8-Y
      F(3)=FIG(1)*Z1+FIG(2)*Z2+FIG(3)*Z3+FIG(4)*Z4+FIG(5)*Z5+FIG(6)*Z6+FIG(7)*Z7+FIG(8)*Z8-Z

C     Jacobijeva matrika J_ij=\p F_i / \p x_j      
      J(1,1)=-ETA2M*ETA3M*X1+ETA2M*ETA3M*X2-ETA2P*ETA3M*X3+ETA2P*ETA3M*X4
     &       -ETA2M*ETA3P*X5+ETA2M*ETA3P*X6-ETA2P*ETA3P*X7+ETA2P*ETA3P*X8

      J(2,1)=-ETA2M*ETA3M*Y1+ETA2M*ETA3M*Y2-ETA2P*ETA3M*Y3+ETA2P*ETA3M*Y4
     &       -ETA2M*ETA3P*Y5+ETA2M*ETA3P*Y6-ETA2P*ETA3P*Y7+ETA2P*ETA3P*Y8

      J(3,1)=-ETA2M*ETA3M*Z1+ETA2M*ETA3M*Z2-ETA2P*ETA3M*Z3+ETA2P*ETA3M*Z4
     &       -ETA2M*ETA3P*Z5+ETA2M*ETA3P*Z6-ETA2P*ETA3P*Z7+ETA2P*ETA3P*Z8
c     
      J(1,2)=-ETA1M*ETA3M*X1-ETA1P*ETA3M*X2+ETA1M*ETA3M*X3+ETA1P*ETA3M*X4
     &       -ETA1M*ETA3P*X5-ETA1P*ETA3P*X6+ETA1M*ETA3P*X7+ETA1P*ETA3P*X8

      J(2,2)=-ETA1M*ETA3M*Y1-ETA1P*ETA3M*Y2+ETA1M*ETA3M*Y3+ETA1P*ETA3M*Y4
     &       -ETA1M*ETA3P*Y5-ETA1P*ETA3P*Y6+ETA1M*ETA3P*Y7+ETA1P*ETA3P*Y8

      J(3,2)=-ETA1M*ETA3M*Z1-ETA1P*ETA3M*Z2+ETA1M*ETA3M*Z3+ETA1P*ETA3M*Z4
     &       -ETA1M*ETA3P*Z5-ETA1P*ETA3P*Z6+ETA1M*ETA3P*Z7+ETA1P*ETA3P*Z8
c
      J(1,3)=-ETA1M*ETA2M*X1-ETA1P*ETA2M*X2-ETA1M*ETA2P*X3-ETA1P*ETA2P*X4
     &       +ETA1M*ETA2M*X5+ETA1P*ETA2M*X6+ETA1M*ETA2P*X7+ETA1P*ETA2P*X8

      J(2,3)=-ETA1M*ETA2M*Y1-ETA1P*ETA2M*Y2-ETA1M*ETA2P*Y3-ETA1P*ETA2P*Y4
     &       +ETA1M*ETA2M*Y5+ETA1P*ETA2M*Y6+ETA1M*ETA2P*Y7+ETA1P*ETA2P*Y8

      J(3,3)=-ETA1M*ETA2M*Z1-ETA1P*ETA2M*Z2-ETA1M*ETA2P*Z3-ETA1P*ETA2P*Z4
     &       +ETA1M*ETA2M*Z5+ETA1P*ETA2M*Z6+ETA1M*ETA2P*Z7+ETA1P*ETA2P*Z8      
      
      END

C ------------------------------------------------------------------------------------
      SUBROUTINE lks2kks(xi,eta,zeta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                               x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x,y,z)
C
C     transformacija iz lokalnega koor. sistema (xi,eta,zeta) v kartezijevega (x,y,z)
C     za heksaeder dolocen z 8 vogali
C    
C ------------------------------------------------------------------------------------          
      REAL(8) xi,eta,zeta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8) x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x,y,z
      REAL(8) FIG(8),eta1m,eta1p,eta2m,eta2p,eta3m,eta3p

      ETA1M=1.0D0-xi
      ETA1P=1.0D0+xi   
      ETA2M=1.0D0-eta
      ETA2P=1.0D0+eta   
      ETA3M=1.0D0-zeta   
      ETA3P=1.0D0+zeta  

      FIG(1)=0.125D0*ETA1M*ETA2M*ETA3M
      FIG(2)=0.125D0*ETA1P*ETA2M*ETA3M
      FIG(3)=0.125D0*ETA1M*ETA2P*ETA3M
      FIG(4)=0.125D0*ETA1P*ETA2P*ETA3M
      FIG(5)=0.125D0*ETA1M*ETA2M*ETA3P
      FIG(6)=0.125D0*ETA1P*ETA2M*ETA3P
      FIG(7)=0.125D0*ETA1M*ETA2P*ETA3P
      FIG(8)=0.125D0*ETA1P*ETA2P*ETA3P

      X=FIG(1)*X1+FIG(2)*X2+FIG(3)*X3+FIG(4)*X4+FIG(5)*X5+FIG(6)*X6+FIG(7)*X7+FIG(8)*X8
      Y=FIG(1)*Y1+FIG(2)*Y2+FIG(3)*Y3+FIG(4)*Y4+FIG(5)*Y5+FIG(6)*Y6+FIG(7)*Y7+FIG(8)*Y8
      Z=FIG(1)*Z1+FIG(2)*Z2+FIG(3)*Z3+FIG(4)*Z4+FIG(5)*Z5+FIG(6)*Z6+FIG(7)*Z7+FIG(8)*Z8
      
      END
     


      
C ------------------------------------------------------------------------------------        
      SUBROUTINE mnewt(ntrial,xx,n,tolx,tolf,x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                         x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,nit)
      INTEGER n,ntrial
      REAL(8) tolf,tolx,xx(n)
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8) x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x,y,z      
CU    USES lubksb,ludcmp,usrfun
      INTEGER i,k,indx(n),nit
      REAL(8) d,errf,errx,fjac(n,n),fvec(n),p(n)
      do 14  k=1,ntrial
        nit=k
        CALL kkk2lks_userfun(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,xx(1),xx(2),xx(3),Fvec,Fjac)
        errf=0.
        do 11 i=1,n
          errf=errf+abs(fvec(i))
11      continue
        if(errf.le.tolf)return
        do 12 i=1,n
          p(i)=-fvec(i)
12      continue
        call ludcmp(fjac,n,n,indx,d)
        call lubksb(fjac,n,n,indx,p)
        errx=0.
        do 13 i=1,n
          errx=errx+abs(p(i))
          xx(i)=xx(i)+p(i)
13      continue
        if(errx.le.tolx)return
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *5ji6.)+.  
C ------------------------------------------------------------------------------------        
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL(8) a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL(8) sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *5ji6.)+.
C ------------------------------------------------------------------------------------      
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL(8) d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL(8) aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) then
          print *,'singular matrix in ludcmp'
          stop
        end if
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *5ji6.)+.
