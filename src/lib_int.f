
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE INTEDc27dcLap
     &(
     &    gp,xp,yp,zp,
     &    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &    x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &    be,dxe,dye,dze,isrc,nsipmx
     &)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  INTEgration over Domain cell                                    **
c **  ----             -                                              **
c **  internal cell geometry : LINEAR ( 8-node )                      **
c **  internal cell function interpolation : QUADRATIC ( 27-node )    **
c **  fundamental solution   : ELLIPTIC LAPLACE                       **
c **                                                                  **
c **********************************************************************
      USE inc_types 
      TYPE(gausstype) :: gp
c
      INTEGER i,j,k,l,n,p,nsipmx,isrc,ng1s,ng2s,ng1,ng2
c      
      REAL(8) xp,yp,zp,
     &        x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &        x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &        xc0,yc0,zc0,
     &        xx1,yy1,zz1,
     &        pi,pi2,
     &        aa0,gama,ffi0,ffi1,ffi,ro,ro0,ro1,th,th0,th1,aa1,rr1,eps,
     &        eti,etj,etk,
     &        eta1m,eta1p,eta2m,eta2p,eta3m,eta3p,
     &        aj11,aj12,aj13,aj21,aj22,aj23,aj31,aj32,aj33,
     &        ajac,rx,ry,rz,ra,ira2,
     &        fung,funx,funy,funz,
     &        ksi,eta,zet,domega
c
      REAL(8) fig27(nsipmx),fig(8)
      REAL(8) dxe(nsipmx),dye(nsipmx),dze(nsipmx)
      REAL(8) be(nsipmx)
      DIMENSION ksi(51),eta(51),zet(51)

      REAL(8) a,b,c,d,e,f,dex,gii,gij,gik
      INTEGER idivXi,ndivXi,idivEt,ndivEt,idivZt,ndivZt,ising

c
      DATA ksi /0.0D+0,0.0D+0,1.0D+0,1.0D+0,0.0D+0,0.0D+0,1.0D+0,1.0D+0,0.0D+0,
     &          5.0D-1,1.0D+0,5.0D-1,0.0D+0,0.0D+0,1.0D+0,1.0D+0,0.0D+0,5.0D-1,
     &          1.0D+0,5.0D-1,5.0D-1,0.0D+0,5.0D-1,1.0D+0,5.0D-1,5.0D-1,5.0D-1,
     &          0.125D0,0.125D0,0.875D0,0.875D0,0.125D0,0.875D0,0.875D0,0.125D0,
     &          1.000D0,1.000D0,1.000D0,1.000D0,0.875D0,0.125D0,0.125D0,0.875D0,
     &          0.000D0,0.000D0,0.000D0,0.000D0,0.125D0,0.875D0,0.875D0,0.125D0/
      DATA eta /1.0D+0,1.0D+0,1.0D+0,1.0D+0,0.0D+0,0.0D+0,0.0D+0,0.0D+0,1.0D+0,
     &          1.0D+0,1.0D+0,1.0D+0,5.0D-1,5.0D-1,5.0D-1,5.0D-1,0.0D+0,0.0D+0,
     &          0.0D+0,0.0D+0,1.0D+0,5.0D-1,5.0D-1,5.0D-1,5.0D-1,0.0D+0,5.0D-1,
     &          0.125D0,0.875D0,0.875D0,0.125D0,0.000D0,0.000D0,0.000D0,0.000D0,
     &          0.125D0,0.875D0,0.875D0,0.125D0,1.000D0,1.000D0,1.000D0,1.000D0,
     &          0.875D0,0.125D0,0.125D0,0.875D0,0.125D0,0.125D0,0.875D0,0.875D0/
      DATA zet /1.0D+0,0.0D+0,0.0D+0,1.0D+0,1.0D+0,0.0D+0,0.0D+0,1.0D+0,5.0D-1,
     &          0.0D+0,5.0D-1,1.0D+0,1.0D+0,0.0D+0,0.0D+0,1.0D+0,5.0D-1,0.0D+0,
     &          5.0D-1,1.0D+0,5.0D-1,5.0D-1,0.0D+0,5.0D-1,1.0D+0,5.0D-1,5.0D-1,
     &          0.000D0,0.000D0,0.000D0,0.000D0,0.125D0,0.125D0,0.875D0,0.875D0,
     &          0.125D0,0.125D0,0.875D0,0.875D0,0.125D0,0.125D0,0.875D0,0.875D0,
     &          0.125D0,0.125D0,0.875D0,0.875D0,1.000D0,1.000D0,1.000D0,1.000D0/          
C
C*** SET NUMBER OF INTEGRATION POINTS
C
c     singular      
      ng1s=gp%ng1(gp%kmDs)
      ng2s=gp%ng2(gp%kmDs)
c     regular      
      ng1=gp%ng1(gp%kmDr)
      ng2=gp%ng2(gp%kmDr)
C
C***  27 NODE CONTINUOUS DOMAIN CELL (8 NODE GEOMETRY)
C
      PI=2.*ASIN(1.0D0)
      PI2=PI*2.0D0
c
      dxe=0.00D00
      dye=0.00D00
      dze=0.00D00
      be=0.00D00

C     razdelitev na podobmocja
      IF (isrc.NE.0) THEN
        ndivXi=gp%iiDiv
        ndivEt=gp%iiDiv
        ndivZt=gp%iiDiv
      ELSE
        ndivXi=1
        ndivEt=1
        ndivZt=1
      END IF

      DO idivXi=1,ndivXi
        DO idivEt=1,ndivEt
          DO idivZt=1,ndivZt
          a=-1.0D0+(idivXi-1)*2.0D0/ndivXi
          b=-1.0D0+(idivXi)*2.0D0/ndivXi          
          c=-1.0D0+(idivEt-1)*2.0D0/ndivEt
          d=-1.0D0+(idivEt)*2.0D0/ndivEt
          e=-1.0D0+(idivZt-1)*2.0D0/ndivZt
          f=-1.0D0+(idivZt)*2.0D0/ndivZt
          dex=0.125D0*(b-a)*(d-c)*(f-e) 

          IF (isrc.NE.0) THEN
            XX1=ksi(isrc)  ! med 0 in 1
            XX1=-1.0D0+2.0D0*XX1  ! med -1 in 1
            XX1=(XX1-a)/(b-a)  ! med 0 in 1 v intervalu a,b
      
            YY1=eta(isrc)  ! med 0 in 1
            YY1=-1.0D0+2.0D0*YY1  ! med -1 in 1
            YY1=(YY1-c)/(d-c)  ! med 0 in 1 v intervalu c,d

            ZZ1=zet(isrc)  ! med 0 in 1
            ZZ1=-1.0D0+2.0D0*ZZ1  ! med -1 in 1
            ZZ1=(ZZ1-e)/(f-e)  ! med 0 in 1 v intervalu e,f
          
            IF (XX1.GE.0.0D0.AND.XX1.LE.1.0D0.AND.YY1.GE.0.0D0.AND.YY1.LE.1.0D0.AND.ZZ1.GE.0.0D0.AND.ZZ1.LE.1.0D0) THEN
              ising=1
            ELSE
              ising=0
            END IF
          ELSE
            ising=0
          ENDIF   
c
C*** SINGULAR INTEGRALS
c
      IF (ising.NE.0) THEN
c        XX1=ksi(isrc)
c        YY1=eta(isrc)
c        ZZ1=zet(isrc)
c
        DO L=1,4
c
          IF (L.EQ.1) THEN
            IF (1.0D0-XX1.EQ.0.0D0) GOTO 1000
            aA0=SQRT(YY1**2+(1.0D0-XX1)**2)
            gama=ACOS(YY1/aA0)
            ffi0=1.5D0*PI+ACOS(YY1/aa0)
            ffi1=PI2+ATAN((1.0D0-YY1)/(1.0D0-XX1))
          ELSE IF (L.EQ.2) THEN
            IF (1.0D0-YY1.EQ.0.0D0) GOTO 1000
            aA0=SQRT((1.0D0-YY1)**2+(1.0D0-XX1)**2)
            gama=ACOS((1.0D0-XX1)/aA0)
            ffi0=ASIN((1.0D0-YY1)/aA0)
            ffi1=0.5D0*PI+ATAN(XX1/(1.0D0-YY1))
          ELSE IF (L.EQ.3) THEN
            IF (XX1.EQ.0.0D0) GOTO 1000
            aA0=SQRT(XX1**2+(1.0D0-YY1)**2)
            gama=ACOS((1.0D0-YY1)/aA0)
            ffi0=0.5D0*PI+ACOS((1.0D0-YY1)/aA0)
            ffi1=PI+ATAN(YY1/XX1)
          ELSE
            IF (YY1.EQ.0.0D0) GOTO 1000
            aA0=SQRT(YY1**2+XX1**2)
            gama=ACOS(XX1/aA0)
            ffi0=PI+ACOS(XX1/aA0)
            ffi1=1.5D0*PI+ATAN((1.0D0-XX1)/YY1)
          END IF
c        
          DO K=ng1s,ng2s
            ffi=(ffi1-ffi0)*gp%GI(k)/2.0D0+(ffi1+ffi0)/2.0D0
            RO0=Aa0*SIN(gama)/SIN(ffi-ffi0+gama)
c
            DO N=1,3
              IF (N.eq.1) THEN
                aa1=1.0D0-zz1
                eps=pi/2.0D0
                th0=0.0D0
                th1=ACOS(aa1/SQRT(aa1**2+ro0**2))
              ELSE IF (N.eq.2) THEN
                aa1=SQRT(ro0**2+(1.0D0-zz1)**2)
                eps=ACOS((1.0D0-zz1)/aa1)
                th0=eps
                rr1=MIN(1.0D00,ro0/SQRT(zz1**2+ro0**2))
                th1=0.5D0*pi+ACOS(rr1)
              ELSE
                aa1=SQRT(ro0**2+zz1**2)
                rr1=MIN(1.0D00,ro0/aa1)
                eps=ACOS(rr1)
                th0=0.5D0*pi+eps
                th1=pi
              END IF
c
              DO I=ng1s,ng2s
                th=(th1-th0)*gp%GI(I)/2.0D0+(th1+th0)/2.0D0
                RO1=Aa1*SIN(eps)/SIN(th-th0+eps)
c
                DO J=ng1s,ng2s
                    ro=(ro1*gp%gi(j)+ro1)/2.0D0
c
                  eti=2.0D0*(xx1+ro*SIN(th)*COS(ffi))-1.0D0
                  etj=2.0D0*(yy1+ro*SIN(th)*SIN(ffi))-1.0D0
                  etk=2.0D0*(zz1+ro*COS(th))-1.0D0

                  ETI=a+0.5D0*(ETI+1.0D0)*(b-a)     
                  ETJ=c+0.5D0*(ETJ+1.0D0)*(d-c) 
                  ETK=e+0.5D0*(ETK+1.0D0)*(f-e) 

                 ETA1M=1.0D0-ETI   
                 ETA1P=1.0D0+ETI   
                 ETA2M=1.0D0-ETJ   
                 ETA2P=1.0D0+ETJ   
                 ETA3M=1.0D0-ETK   
                 ETA3P=1.0D0+ETK   

                  FIG(1)=0.125D0*ETA1M*ETA2M*ETA3M
                  FIG(2)=0.125D0*ETA1P*ETA2M*ETA3M
                  FIG(3)=0.125D0*ETA1M*ETA2P*ETA3M
                  FIG(4)=0.125D0*ETA1P*ETA2P*ETA3M
                  FIG(5)=0.125D0*ETA1M*ETA2M*ETA3P
                  FIG(6)=0.125D0*ETA1P*ETA2M*ETA3P
                  FIG(7)=0.125D0*ETA1M*ETA2P*ETA3P
                  FIG(8)=0.125D0*ETA1P*ETA2P*ETA3P

                  XC0=FIG(1)*X1+FIG(2)*X2+FIG(3)*X3+FIG(4)*X4
     &               +FIG(5)*X5+FIG(6)*X6+FIG(7)*X7+FIG(8)*X8
                  YC0=FIG(1)*Y1+FIG(2)*Y2+FIG(3)*Y3+FIG(4)*Y4
     &               +FIG(5)*Y5+FIG(6)*Y6+FIG(7)*Y7+FIG(8)*Y8
                  ZC0=FIG(1)*Z1+FIG(2)*Z2+FIG(3)*Z3+FIG(4)*Z4
     &               +FIG(5)*Z5+FIG(6)*Z6+FIG(7)*Z7+FIG(8)*Z8
C
C***  Jacobian derivatives
C
              AJ11=-ETA2M*ETA3M*X1+ETA2M*ETA3M*X2-ETA2P*ETA3M*X3+ETA2P*ETA3M*X4
     &                 -ETA2M*ETA3P*X5+ETA2M*ETA3P*X6-ETA2P*ETA3P*X7+ETA2P*ETA3P*X8

                 AJ12=-ETA2M*ETA3M*Y1+ETA2M*ETA3M*Y2-ETA2P*ETA3M*Y3+ETA2P*ETA3M*Y4
     &                 -ETA2M*ETA3P*Y5+ETA2M*ETA3P*Y6-ETA2P*ETA3P*Y7+ETA2P*ETA3P*Y8

                 AJ13=-ETA2M*ETA3M*Z1+ETA2M*ETA3M*Z2-ETA2P*ETA3M*Z3+ETA2P*ETA3M*Z4
     &                 -ETA2M*ETA3P*Z5+ETA2M*ETA3P*Z6-ETA2P*ETA3P*Z7+ETA2P*ETA3P*Z8
c     
                 AJ21=-ETA1M*ETA3M*X1-ETA1P*ETA3M*X2+ETA1M*ETA3M*X3+ETA1P*ETA3M*X4
     &                 -ETA1M*ETA3P*X5-ETA1P*ETA3P*X6+ETA1M*ETA3P*X7+ETA1P*ETA3P*X8

                 AJ22=-ETA1M*ETA3M*Y1-ETA1P*ETA3M*Y2+ETA1M*ETA3M*Y3+ETA1P*ETA3M*Y4
     &                 -ETA1M*ETA3P*Y5-ETA1P*ETA3P*Y6+ETA1M*ETA3P*Y7+ETA1P*ETA3P*Y8

                 AJ23=-ETA1M*ETA3M*Z1-ETA1P*ETA3M*Z2+ETA1M*ETA3M*Z3+ETA1P*ETA3M*Z4
     &                 -ETA1M*ETA3P*Z5-ETA1P*ETA3P*Z6+ETA1M*ETA3P*Z7+ETA1P*ETA3P*Z8
c
                AJ31=-ETA1M*ETA2M*X1-ETA1P*ETA2M*X2-ETA1M*ETA2P*X3-ETA1P*ETA2P*X4
     &                 +ETA1M*ETA2M*X5+ETA1P*ETA2M*X6+ETA1M*ETA2P*X7+ETA1P*ETA2P*X8

                 AJ32=-ETA1M*ETA2M*Y1-ETA1P*ETA2M*Y2-ETA1M*ETA2P*Y3-ETA1P*ETA2P*Y4
     &                 +ETA1M*ETA2M*Y5+ETA1P*ETA2M*Y6+ETA1M*ETA2P*Y7+ETA1P*ETA2P*Y8

                 AJ33=-ETA1M*ETA2M*Z1-ETA1P*ETA2M*Z2-ETA1M*ETA2P*Z3-ETA1P*ETA2P*Z4
     &                 +ETA1M*ETA2M*Z5+ETA1P*ETA2M*Z6+ETA1M*ETA2P*Z7+ETA1P*ETA2P*Z8
c
                  AJAC=(Aj11*Aj22*Aj33+Aj12*Aj23*Aj31+Aj13*Aj21*Aj32-
     &                  Aj11*Aj23*Aj32-Aj12*Aj21*Aj33-Aj13*Aj22*Aj31)
c
                  AJAC=1.953125D-3*AJAC*gp%OME(I)*gp%OME(J)*gp%OME(K)*(ffi1-ffi0)*(th1-th0)*ro1*ro**2*SIN(th)
                  AJAC=dex*AJAC

                  RX=XP-XC0
                  RY=YP-YC0
                  RZ=ZP-ZC0
                  RA=SQRT(RX*RX+RY*RY+RZ*RZ)
                  IF (RA.EQ.0.0D00) GO TO 800
                  IF (RA.LT.1.0D-12) RA=1.0D-12
                  IRA2=1.0D00/(RA*RA)
C
C***              ELLIPTIC LAPLACE FUNDAMENTAL SOLUTION
C
                  FUNG=0.25D0/(PI*RA)
                  FUNX=(XP-XC0)*IRA2*FUNG
                  FUNY=(YP-YC0)*IRA2*FUNG
                  FUNZ=(ZP-ZC0)*IRA2*FUNG
c
                  CALL cshape27(eti,etj,etk,fig27,nsipmx)
c
                  DO P=1,nsipmx
                    DOMEGA=AJAC*FIG27(P)
                    dxe(p)=dxe(p)+FUNX*DOMEGA
                    dye(p)=dye(p)+FUNY*DOMEGA
                    dze(p)=dze(p)+FUNZ*DOMEGA
                    be(p) =be(p) +FUNG*DOMEGA
                  END DO      
c
  800           END DO !J
              END DO !I
            END DO !N
          END DO !K
 1000   END DO !L
      ELSE
        DO K=ng1,ng2
          gik=e+0.5D0*(gp%GI(K)+1.0D0)*(f-e)    
          ETA3M=1.0D0-gik
          ETA3P=1.0D0+gik
          DO J=ng1,ng2
            gij=c+0.5D0*(gp%GI(J)+1.0D0)*(d-c)   
            ETA2M=1.0D0-gij
            ETA2P=1.0D0+gij
            DO I=ng1,ng2
              gii=a+0.5D0*(gp%GI(I)+1.0D0)*(b-a)
              ETA1M=1.0D0-gii
              ETA1P=1.0D0+gii

              FIG(1)=0.125D0*ETA1M*ETA2M*ETA3M
              FIG(2)=0.125D0*ETA1P*ETA2M*ETA3M
              FIG(3)=0.125D0*ETA1M*ETA2P*ETA3M
              FIG(4)=0.125D0*ETA1P*ETA2P*ETA3M
              FIG(5)=0.125D0*ETA1M*ETA2M*ETA3P
              FIG(6)=0.125D0*ETA1P*ETA2M*ETA3P
              FIG(7)=0.125D0*ETA1M*ETA2P*ETA3P
              FIG(8)=0.125D0*ETA1P*ETA2P*ETA3P

              XC0=FIG(1)*X1+FIG(2)*X2+FIG(3)*X3+FIG(4)*X4
     &           +FIG(5)*X5+FIG(6)*X6+FIG(7)*X7+FIG(8)*X8
              YC0=FIG(1)*Y1+FIG(2)*Y2+FIG(3)*Y3+FIG(4)*Y4
     &           +FIG(5)*Y5+FIG(6)*Y6+FIG(7)*Y7+FIG(8)*Y8
              ZC0=FIG(1)*Z1+FIG(2)*Z2+FIG(3)*Z3+FIG(4)*Z4
     &           +FIG(5)*Z5+FIG(6)*Z6+FIG(7)*Z7+FIG(8)*Z8
C
C***  Jacobian derivatives
C
           AJ11=-ETA2M*ETA3M*X1+ETA2M*ETA3M*X2-ETA2P*ETA3M*X3+ETA2P*ETA3M*X4
     &             -ETA2M*ETA3P*X5+ETA2M*ETA3P*X6-ETA2P*ETA3P*X7+ETA2P*ETA3P*X8

         AJ12=-ETA2M*ETA3M*Y1+ETA2M*ETA3M*Y2-ETA2P*ETA3M*Y3+ETA2P*ETA3M*Y4
     &             -ETA2M*ETA3P*Y5+ETA2M*ETA3P*Y6-ETA2P*ETA3P*Y7+ETA2P*ETA3P*Y8

         AJ13=-ETA2M*ETA3M*Z1+ETA2M*ETA3M*Z2-ETA2P*ETA3M*Z3+ETA2P*ETA3M*Z4
     &             -ETA2M*ETA3P*Z5+ETA2M*ETA3P*Z6-ETA2P*ETA3P*Z7+ETA2P*ETA3P*Z8
c 
             AJ21=-ETA1M*ETA3M*X1-ETA1P*ETA3M*X2+ETA1M*ETA3M*X3+ETA1P*ETA3M*X4
     &             -ETA1M*ETA3P*X5-ETA1P*ETA3P*X6+ETA1M*ETA3P*X7+ETA1P*ETA3P*X8

         AJ22=-ETA1M*ETA3M*Y1-ETA1P*ETA3M*Y2+ETA1M*ETA3M*Y3+ETA1P*ETA3M*Y4
     &             -ETA1M*ETA3P*Y5-ETA1P*ETA3P*Y6+ETA1M*ETA3P*Y7+ETA1P*ETA3P*Y8

         AJ23=-ETA1M*ETA3M*Z1-ETA1P*ETA3M*Z2+ETA1M*ETA3M*Z3+ETA1P*ETA3M*Z4
     &             -ETA1M*ETA3P*Z5-ETA1P*ETA3P*Z6+ETA1M*ETA3P*Z7+ETA1P*ETA3P*Z8
c
              AJ31=-ETA1M*ETA2M*X1-ETA1P*ETA2M*X2-ETA1M*ETA2P*X3-ETA1P*ETA2P*X4
     &             +ETA1M*ETA2M*X5+ETA1P*ETA2M*X6+ETA1M*ETA2P*X7+ETA1P*ETA2P*X8

         AJ32=-ETA1M*ETA2M*Y1-ETA1P*ETA2M*Y2-ETA1M*ETA2P*Y3-ETA1P*ETA2P*Y4
     &             +ETA1M*ETA2M*Y5+ETA1P*ETA2M*Y6+ETA1M*ETA2P*Y7+ETA1P*ETA2P*Y8

         AJ33=-ETA1M*ETA2M*Z1-ETA1P*ETA2M*Z2-ETA1M*ETA2P*Z3-ETA1P*ETA2P*Z4
     &             +ETA1M*ETA2M*Z5+ETA1P*ETA2M*Z6+ETA1M*ETA2P*Z7+ETA1P*ETA2P*Z8
c
              AJAC=(Aj11*Aj22*Aj33+Aj12*Aj23*Aj31+Aj13*Aj21*Aj32-
     &              Aj11*Aj23*Aj32-Aj12*Aj21*Aj33-Aj13*Aj22*Aj31)
c
              AJAC=1.953125D-3*AJAC*gp%OME(I)*gp%OME(J)*gp%OME(K)*dex
c                  1/512=1.953125D-3
              RX=XP-XC0
              RY=YP-YC0
              RZ=ZP-ZC0
              RA=SQRT(RX*RX+RY*RY+RZ*RZ)
              IF (RA.LT.1.0D-12) RA=1.0D-12 
              IRA2=1.0D00/(RA*RA)
C
C***          ELLIPTIC LAPLACE FUNDAMENTAL SOLUTION
C
              FUNG=0.25D0/(PI*RA)
              FUNX=(XP-XC0)*IRA2*FUNG
              FUNY=(YP-YC0)*IRA2*FUNG
              FUNZ=(ZP-ZC0)*IRA2*FUNG
c              
C             interpolacijska funkcija, ki jo mnozim zraven
              CALL cshape27(gii,gij,gik,fig27,nsipmx)
c
              DO P=1,nsipmx
                DOMEGA=AJAC*FIG27(P)
                dxe(p)=dxe(p)+FUNX*DOMEGA
                dye(p)=dye(p)+FUNY*DOMEGA
                dze(p)=dze(p)+FUNZ*DOMEGA
                be(p) =be(p) +FUNG*DOMEGA
              END DO      
c
            END DO !I
          END DO !J
        END DO !K
      END IF

            END DO !gkI
          END DO !gkJ
        END DO !gkK

c
      RETURN
      END


C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE INTEDc27testIntFi
     &(
     &    gp,
     &    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &    x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &    IntFi
     &)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  INTEgration over Domain cell                                    **
c **  ----             -                                              **
c **  internal cell geometry : LINEAR ( 8-node )                      **
c **  internal cell function interpolation : QUADRATIC ( 27-node )    **
c **  fundamental solution   : ELLIPTIC LAPLACE                       **
c **                                                                  **
c **********************************************************************
      USE inc_types 
      TYPE(gausstype) :: gp
c
      INTEGER i,j,k,l,n,p,nsipmx,isrc,ng1s,ng2s,ng1,ng2
c      
      REAL(8) xp,yp,zp,
     &        x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &        x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &        xc0,yc0,zc0,
     &        xx1,yy1,zz1,
     &        pi,pi2,
     &        aa0,gama,ffi0,ffi1,ffi,ro,ro0,ro1,th,th0,th1,aa1,rr1,eps,
     &        eti,etj,etk,
     &        eta1m,eta1p,eta2m,eta2p,eta3m,eta3p,
     &        aj11,aj12,aj13,aj21,aj22,aj23,aj31,aj32,aj33,
     &        ajac,rx,ry,rz,ra,ira2,
     &        fung,funx,funy,funz,
     &        ksi,eta,zet,domega
c
      REAL(8) fig27(27),fig(8),IntFi(27)

      REAL(8) a,b,c,d,e,f,dex,gii,gij,gik
      INTEGER idivXi,ndivXi,idivEt,ndivEt,idivZt,ndivZt

C 
C*** SET NUMBER OF INTEGRATION POINTS
C
c     singular      
      ng1s=gp%ng1(gp%kmDs)
      ng2s=gp%ng2(gp%kmDs)
c     regular      
      ng1=gp%ng1(gp%kmDr)
      ng2=gp%ng2(gp%kmDr)
C
C***  27 NODE CONTINUOUS DOMAIN CELL (8 NODE GEOMETRY)
C

      PI=2.*ASIN(1.0D0)
      PI2=PI*2.0D0
c

      IntFi=0.00D00

C     razdelitev na podobmocja
      ndivXi=2
      ndivEt=3
      ndivZt=4


      DO idivXi=1,ndivXi
        DO idivEt=1,ndivEt
          DO idivZt=1,ndivZt
          a=-1.0D0+(idivXi-1)*2.0D0/ndivXi
          b=-1.0D0+(idivXi)*2.0D0/ndivXi          
          c=-1.0D0+(idivEt-1)*2.0D0/ndivEt
          d=-1.0D0+(idivEt)*2.0D0/ndivEt
          e=-1.0D0+(idivZt-1)*2.0D0/ndivZt
          f=-1.0D0+(idivZt)*2.0D0/ndivZt
          dex=0.125D0*(b-a)*(d-c)*(f-e) 

C     gauss point loops
        DO K=ng1,ng2
          gik=e+0.5D0*(gp%GI(K)+1.0D0)*(f-e)    
          ETA3M=1.0D0-gik
          ETA3P=1.0D0+gik
          DO J=ng1,ng2
            gij=c+0.5D0*(gp%GI(J)+1.0D0)*(d-c)   
            ETA2M=1.0D0-gij
            ETA2P=1.0D0+gij
            DO I=ng1,ng2
              gii=a+0.5D0*(gp%GI(I)+1.0D0)*(b-a)
              ETA1M=1.0D0-gii
              ETA1P=1.0D0+gii


              FIG(1)=0.125D0*ETA1M*ETA2M*ETA3M
              FIG(2)=0.125D0*ETA1P*ETA2M*ETA3M
              FIG(3)=0.125D0*ETA1M*ETA2P*ETA3M
              FIG(4)=0.125D0*ETA1P*ETA2P*ETA3M
              FIG(5)=0.125D0*ETA1M*ETA2M*ETA3P
              FIG(6)=0.125D0*ETA1P*ETA2M*ETA3P
              FIG(7)=0.125D0*ETA1M*ETA2P*ETA3P
              FIG(8)=0.125D0*ETA1P*ETA2P*ETA3P

C             tocka v (x,y,z) koordinatem sistemu
              XC0=FIG(1)*X1+FIG(2)*X2+FIG(3)*X3+FIG(4)*X4
     &           +FIG(5)*X5+FIG(6)*X6+FIG(7)*X7+FIG(8)*X8
              YC0=FIG(1)*Y1+FIG(2)*Y2+FIG(3)*Y3+FIG(4)*Y4
     &           +FIG(5)*Y5+FIG(6)*Y6+FIG(7)*Y7+FIG(8)*Y8
              ZC0=FIG(1)*Z1+FIG(2)*Z2+FIG(3)*Z3+FIG(4)*Z4
     &           +FIG(5)*Z5+FIG(6)*Z6+FIG(7)*Z7+FIG(8)*Z8
C
C***  Jacobian derivatives
C
           AJ11=-ETA2M*ETA3M*X1+ETA2M*ETA3M*X2-ETA2P*ETA3M*X3+ETA2P*ETA3M*X4
     &             -ETA2M*ETA3P*X5+ETA2M*ETA3P*X6-ETA2P*ETA3P*X7+ETA2P*ETA3P*X8

         AJ12=-ETA2M*ETA3M*Y1+ETA2M*ETA3M*Y2-ETA2P*ETA3M*Y3+ETA2P*ETA3M*Y4
     &             -ETA2M*ETA3P*Y5+ETA2M*ETA3P*Y6-ETA2P*ETA3P*Y7+ETA2P*ETA3P*Y8

         AJ13=-ETA2M*ETA3M*Z1+ETA2M*ETA3M*Z2-ETA2P*ETA3M*Z3+ETA2P*ETA3M*Z4
     &             -ETA2M*ETA3P*Z5+ETA2M*ETA3P*Z6-ETA2P*ETA3P*Z7+ETA2P*ETA3P*Z8
c 
             AJ21=-ETA1M*ETA3M*X1-ETA1P*ETA3M*X2+ETA1M*ETA3M*X3+ETA1P*ETA3M*X4
     &             -ETA1M*ETA3P*X5-ETA1P*ETA3P*X6+ETA1M*ETA3P*X7+ETA1P*ETA3P*X8

         AJ22=-ETA1M*ETA3M*Y1-ETA1P*ETA3M*Y2+ETA1M*ETA3M*Y3+ETA1P*ETA3M*Y4
     &             -ETA1M*ETA3P*Y5-ETA1P*ETA3P*Y6+ETA1M*ETA3P*Y7+ETA1P*ETA3P*Y8

         AJ23=-ETA1M*ETA3M*Z1-ETA1P*ETA3M*Z2+ETA1M*ETA3M*Z3+ETA1P*ETA3M*Z4
     &             -ETA1M*ETA3P*Z5-ETA1P*ETA3P*Z6+ETA1M*ETA3P*Z7+ETA1P*ETA3P*Z8
c
              AJ31=-ETA1M*ETA2M*X1-ETA1P*ETA2M*X2-ETA1M*ETA2P*X3-ETA1P*ETA2P*X4
     &             +ETA1M*ETA2M*X5+ETA1P*ETA2M*X6+ETA1M*ETA2P*X7+ETA1P*ETA2P*X8

         AJ32=-ETA1M*ETA2M*Y1-ETA1P*ETA2M*Y2-ETA1M*ETA2P*Y3-ETA1P*ETA2P*Y4
     &             +ETA1M*ETA2M*Y5+ETA1P*ETA2M*Y6+ETA1M*ETA2P*Y7+ETA1P*ETA2P*Y8

         AJ33=-ETA1M*ETA2M*Z1-ETA1P*ETA2M*Z2-ETA1M*ETA2P*Z3-ETA1P*ETA2P*Z4
     &             +ETA1M*ETA2M*Z5+ETA1P*ETA2M*Z6+ETA1M*ETA2P*Z7+ETA1P*ETA2P*Z8
c
              AJAC=(Aj11*Aj22*Aj33+Aj12*Aj23*Aj31+Aj13*Aj21*Aj32-
     &              Aj11*Aj23*Aj32-Aj12*Aj21*Aj33-Aj13*Aj22*Aj31)
c
              AJAC=1.953125D-3*AJAC*gp%OME(I)*gp%OME(J)*gp%OME(K)*dex
c                  1/512=1.953125D-3


              FUNG = XC0*YC0*ZC0 ! funkcija, ki jo integrairam * int F = x*y*z

C             interpolacijska funkcija, ki jo mnozim zraven
              CALL cshape27(gii,gij,gik,fig27,27)

c
              DO P=1,27
                IntFi(P) = IntFi(P) + AJAC*FIG27(P)*FUNG
              END DO      
c
            END DO !I
          END DO !J
        END DO !K

            END DO !gkI
          END DO !gkJ
        END DO !gkK
c
 
      RETURN
      END





C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE INTEDc27testF
     &(
     &    gp,
     &    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &    x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &    rezultat
     &)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  INTEgration over Domain cell                                    **
c **  ----             -                                              **
c **  internal cell geometry : LINEAR ( 8-node )                      **
c **  internal cell function interpolation : QUADRATIC ( 27-node )    **
c **  fundamental solution   : ELLIPTIC LAPLACE                       **
c **                                                                  **
c **********************************************************************
      USE inc_types 
      TYPE(gausstype) :: gp
c
      INTEGER i,j,k,l,n,p,nsipmx,isrc,ng1s,ng2s,ng1,ng2
c      
      REAL(8) xp,yp,zp,
     &        x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &        x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &        xc0,yc0,zc0,
     &        xx1,yy1,zz1,
     &        pi,pi2,
     &        aa0,gama,ffi0,ffi1,ffi,ro,ro0,ro1,th,th0,th1,aa1,rr1,eps,
     &        eti,etj,etk,
     &        eta1m,eta1p,eta2m,eta2p,eta3m,eta3p,
     &        aj11,aj12,aj13,aj21,aj22,aj23,aj31,aj32,aj33,
     &        ajac,rx,ry,rz,ra,ira2,
     &        fung,funx,funy,funz,
     &        ksi,eta,zet,domega
c
      REAL(8) fig27(27),fig(8),rezultat

      REAL(8) a,b,c,d,e,f,dex,gii,gij,gik
      INTEGER idivXi,ndivXi,idivEt,ndivEt,idivZt,ndivZt

C 
C*** SET NUMBER OF INTEGRATION POINTS
C
c     singular      
      ng1s=gp%ng1(gp%kmDs)
      ng2s=gp%ng2(gp%kmDs)
c     regular      
      ng1=gp%ng1(gp%kmDr)
      ng2=gp%ng2(gp%kmDr)
C
C***  27 NODE CONTINUOUS DOMAIN CELL (8 NODE GEOMETRY)
C

      PI=2.*ASIN(1.0D0)
      PI2=PI*2.0D0
c

      rezultat=0.00D00

C     razdelitev na podobmocja
      ndivXi=2
      ndivEt=3
      ndivZt=4


      DO idivXi=1,ndivXi
        DO idivEt=1,ndivEt
          DO idivZt=1,ndivZt
          a=-1.0D0+(idivXi-1)*2.0D0/ndivXi
          b=-1.0D0+(idivXi)*2.0D0/ndivXi          
          c=-1.0D0+(idivEt-1)*2.0D0/ndivEt
          d=-1.0D0+(idivEt)*2.0D0/ndivEt
          e=-1.0D0+(idivZt-1)*2.0D0/ndivZt
          f=-1.0D0+(idivZt)*2.0D0/ndivZt
          dex=0.125D0*(b-a)*(d-c)*(f-e) 

C     gauss point loops
        DO K=ng1,ng2
          gik=e+0.5D0*(gp%GI(K)+1.0D0)*(f-e)    
          ETA3M=1.0D0-gik
          ETA3P=1.0D0+gik
          DO J=ng1,ng2
            gij=c+0.5D0*(gp%GI(J)+1.0D0)*(d-c)   
            ETA2M=1.0D0-gij
            ETA2P=1.0D0+gij
            DO I=ng1,ng2
              gii=a+0.5D0*(gp%GI(I)+1.0D0)*(b-a)
              ETA1M=1.0D0-gii
              ETA1P=1.0D0+gii


              FIG(1)=0.125D0*ETA1M*ETA2M*ETA3M
              FIG(2)=0.125D0*ETA1P*ETA2M*ETA3M
              FIG(3)=0.125D0*ETA1M*ETA2P*ETA3M
              FIG(4)=0.125D0*ETA1P*ETA2P*ETA3M
              FIG(5)=0.125D0*ETA1M*ETA2M*ETA3P
              FIG(6)=0.125D0*ETA1P*ETA2M*ETA3P
              FIG(7)=0.125D0*ETA1M*ETA2P*ETA3P
              FIG(8)=0.125D0*ETA1P*ETA2P*ETA3P

C             tocka v (x,y,z) koordinatem sistemu
              XC0=FIG(1)*X1+FIG(2)*X2+FIG(3)*X3+FIG(4)*X4
     &           +FIG(5)*X5+FIG(6)*X6+FIG(7)*X7+FIG(8)*X8
              YC0=FIG(1)*Y1+FIG(2)*Y2+FIG(3)*Y3+FIG(4)*Y4
     &           +FIG(5)*Y5+FIG(6)*Y6+FIG(7)*Y7+FIG(8)*Y8
              ZC0=FIG(1)*Z1+FIG(2)*Z2+FIG(3)*Z3+FIG(4)*Z4
     &           +FIG(5)*Z5+FIG(6)*Z6+FIG(7)*Z7+FIG(8)*Z8
C
C***  Jacobian derivatives
C
           AJ11=-ETA2M*ETA3M*X1+ETA2M*ETA3M*X2-ETA2P*ETA3M*X3+ETA2P*ETA3M*X4
     &             -ETA2M*ETA3P*X5+ETA2M*ETA3P*X6-ETA2P*ETA3P*X7+ETA2P*ETA3P*X8

         AJ12=-ETA2M*ETA3M*Y1+ETA2M*ETA3M*Y2-ETA2P*ETA3M*Y3+ETA2P*ETA3M*Y4
     &             -ETA2M*ETA3P*Y5+ETA2M*ETA3P*Y6-ETA2P*ETA3P*Y7+ETA2P*ETA3P*Y8

         AJ13=-ETA2M*ETA3M*Z1+ETA2M*ETA3M*Z2-ETA2P*ETA3M*Z3+ETA2P*ETA3M*Z4
     &             -ETA2M*ETA3P*Z5+ETA2M*ETA3P*Z6-ETA2P*ETA3P*Z7+ETA2P*ETA3P*Z8
c 
             AJ21=-ETA1M*ETA3M*X1-ETA1P*ETA3M*X2+ETA1M*ETA3M*X3+ETA1P*ETA3M*X4
     &             -ETA1M*ETA3P*X5-ETA1P*ETA3P*X6+ETA1M*ETA3P*X7+ETA1P*ETA3P*X8

         AJ22=-ETA1M*ETA3M*Y1-ETA1P*ETA3M*Y2+ETA1M*ETA3M*Y3+ETA1P*ETA3M*Y4
     &             -ETA1M*ETA3P*Y5-ETA1P*ETA3P*Y6+ETA1M*ETA3P*Y7+ETA1P*ETA3P*Y8

         AJ23=-ETA1M*ETA3M*Z1-ETA1P*ETA3M*Z2+ETA1M*ETA3M*Z3+ETA1P*ETA3M*Z4
     &             -ETA1M*ETA3P*Z5-ETA1P*ETA3P*Z6+ETA1M*ETA3P*Z7+ETA1P*ETA3P*Z8
c
              AJ31=-ETA1M*ETA2M*X1-ETA1P*ETA2M*X2-ETA1M*ETA2P*X3-ETA1P*ETA2P*X4
     &             +ETA1M*ETA2M*X5+ETA1P*ETA2M*X6+ETA1M*ETA2P*X7+ETA1P*ETA2P*X8

         AJ32=-ETA1M*ETA2M*Y1-ETA1P*ETA2M*Y2-ETA1M*ETA2P*Y3-ETA1P*ETA2P*Y4
     &             +ETA1M*ETA2M*Y5+ETA1P*ETA2M*Y6+ETA1M*ETA2P*Y7+ETA1P*ETA2P*Y8

         AJ33=-ETA1M*ETA2M*Z1-ETA1P*ETA2M*Z2-ETA1M*ETA2P*Z3-ETA1P*ETA2P*Z4
     &             +ETA1M*ETA2M*Z5+ETA1P*ETA2M*Z6+ETA1M*ETA2P*Z7+ETA1P*ETA2P*Z8
c
              AJAC=(Aj11*Aj22*Aj33+Aj12*Aj23*Aj31+Aj13*Aj21*Aj32-
     &              Aj11*Aj23*Aj32-Aj12*Aj21*Aj33-Aj13*Aj22*Aj31)
c
              AJAC=1.953125D-3*AJAC*gp%OME(I)*gp%OME(J)*gp%OME(K)*dex
c                  1/512=1.953125D-3

              FUNG = SIN(100.0D0*XC0*YC0*ZC0) ! funkcija, ki jo integrairam
c              
	      rezultat = rezultat + FUNG*AJAC
C

c
            END DO !I
          END DO !J
        END DO !K

            END DO !gkI
          END DO !gkJ
        END DO !gkK
c
      RETURN
      END



C----------------------------------------------------------------------C
c
      SUBROUTINE Test3DIntegration(mesh,gp)
C
C     $: Testiram 3D integral
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(meshType) :: mesh
      TYPE(GaussType) :: gp
      INTEGER ic
      REAL(8) x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
      REAL(8) x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8
      REAL(8) rezultat,vsota

      REAL(8), ALLOCATABLE :: T(:),IntFi(:)

	vsota = 0.0D0
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


        CALL INTEDc27testF(gp,
     &    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &    x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &    rezultat)
	vsota=vsota+rezultat

      END DO
  	print *,"vsota",vsota
C             // Mathematica    0.130174132171184615 - ni do konca skonvergirala
C             // 101010         0.13017412657731764
C             // 444            0.13017412657586619
C             // 333            0.13017412637393846
C             // 222            0.13017420582202685
C             // 111            0.13045507871710810


C     sedaj verzija integrala iterpolirane funkcije
      ALLOCATE (T(mesh%nnodes))
      DO ic=1,mesh%nnodes
        T(ic)=1.0D0+mesh%x(ic,1)**2+mesh%x(ic,2)**2+mesh%x(ic,3)**2
      END DO

      ALLOCATE (IntFi(27))
c     integral interpolacijskih funckij
	vsota = 0.0D0
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


        CALL INTEDc27testIntFi(gp,
     &    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &    x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &    IntFi)
	
        DO i=1,27
          vsota=vsota+IntFi(i)*T(mesh%idc(ic,i))
        END DO
        
      END DO
  	print *,"vsota",vsota
C       integral (1+x^2+y^2+z^2)*(x*y*z) od -1,1,1-1,-1.1
C       analitino 0.3125



      DEALLOCATE(T,IntFi)

      STOP

      END SUBROUTINE


C----------------------------------------------------------------------C
c
      SUBROUTINE CalculateMeshElementVolume(mesh,gp)
C
C     $: Z integracijo determinante Jacobijeve matrike izracuna volumen elementa
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(meshType) :: mesh
      TYPE(GaussType) :: gp
      INTEGER ic
      REAL(8) x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
      REAL(8) x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8

c     total mesh volume
      mesh%MeshVolume=0.0D0

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

c        CALL IntCellVolume(gp,
        CALL IntCellVolume8(gp,
     &    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &    x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &    mesh%CellVolume(ic))

C       total volume of the mesh
        mesh%MeshVolume=mesh%MeshVolume+mesh%CellVolume(ic)

      END DO

      END SUBROUTINE
c
C----------------------------------------------------------------------C
c
      SUBROUTINE IntCellVolume(gp,
     &    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &    x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &    CellVolume)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  INTEgration over Domain cell                                    **
c **  ----             -                                              **
c **  internal cell geometry : LINEAR ( 8-node )                      **
c **  internal cell function interpolation : QUADRATIC ( 27-node )    **
c **  integrira determinanto jacobijeve matrike = izracuna volumen    **
c **                                                                  **
c **********************************************************************
      USE inc_types
      TYPE(gausstype) :: gp
c
      INTEGER i,j,k,p,ng1,ng2
c
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &        x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &        eta1m,eta1p,eta2m,eta2p,eta3m,eta3p,
     &        aj11,aj12,aj13,aj21,aj22,aj23,aj31,aj32,aj33,
     &        ajac

      REAL(8) CellVolume
      REAL(8) fig27(27)
C
C*** SET NUMBER OF INTEGRATION POINTS
      ng1=gp%ng1(2)  ! 2 pomeni pomoje dovolj gaussovnih tock
      ng2=gp%ng2(2)

      CellVolume=0.0D0
C
C***  27 NODE CONTINUOUS DOMAIN CELL (8 NODE GEOMETRY)
C
      DO K=ng1,ng2
        DO J=ng1,ng2
          DO I=ng1,ng2
            ETA1M=1.0D0-gp%GI(I)
            ETA1P=1.0D0+gp%GI(I)
            ETA2M=1.0D0-gp%GI(J)
            ETA2P=1.0D0+gp%GI(J)
            ETA3M=1.0D0-gp%GI(K)
            ETA3P=1.0D0+gp%GI(K)
C
C***  Jacobian derivatives
C
            AJ11=-ETA2M*ETA3M*X1+ETA2M*ETA3M*X2-ETA2P*ETA3M*X3+ETA2P*ETA3M*X4
     &           -ETA2M*ETA3P*X5+ETA2M*ETA3P*X6-ETA2P*ETA3P*X7+ETA2P*ETA3P*X8
            AJ12=-ETA2M*ETA3M*Y1+ETA2M*ETA3M*Y2-ETA2P*ETA3M*Y3+ETA2P*ETA3M*Y4
     &           -ETA2M*ETA3P*Y5+ETA2M*ETA3P*Y6-ETA2P*ETA3P*Y7+ETA2P*ETA3P*Y8
            AJ13=-ETA2M*ETA3M*Z1+ETA2M*ETA3M*Z2-ETA2P*ETA3M*Z3+ETA2P*ETA3M*Z4
     &           -ETA2M*ETA3P*Z5+ETA2M*ETA3P*Z6-ETA2P*ETA3P*Z7+ETA2P*ETA3P*Z8
            AJ21=-ETA1M*ETA3M*X1-ETA1P*ETA3M*X2+ETA1M*ETA3M*X3+ETA1P*ETA3M*X4
     &           -ETA1M*ETA3P*X5-ETA1P*ETA3P*X6+ETA1M*ETA3P*X7+ETA1P*ETA3P*X8
            AJ22=-ETA1M*ETA3M*Y1-ETA1P*ETA3M*Y2+ETA1M*ETA3M*Y3+ETA1P*ETA3M*Y4
     &           -ETA1M*ETA3P*Y5-ETA1P*ETA3P*Y6+ETA1M*ETA3P*Y7+ETA1P*ETA3P*Y8
            AJ23=-ETA1M*ETA3M*Z1-ETA1P*ETA3M*Z2+ETA1M*ETA3M*Z3+ETA1P*ETA3M*Z4
     &           -ETA1M*ETA3P*Z5-ETA1P*ETA3P*Z6+ETA1M*ETA3P*Z7+ETA1P*ETA3P*Z8
            AJ31=-ETA1M*ETA2M*X1-ETA1P*ETA2M*X2-ETA1M*ETA2P*X3-ETA1P*ETA2P*X4
     &           +ETA1M*ETA2M*X5+ETA1P*ETA2M*X6+ETA1M*ETA2P*X7+ETA1P*ETA2P*X8
            AJ32=-ETA1M*ETA2M*Y1-ETA1P*ETA2M*Y2-ETA1M*ETA2P*Y3-ETA1P*ETA2P*Y4
     &           +ETA1M*ETA2M*Y5+ETA1P*ETA2M*Y6+ETA1M*ETA2P*Y7+ETA1P*ETA2P*Y8
            AJ33=-ETA1M*ETA2M*Z1-ETA1P*ETA2M*Z2-ETA1M*ETA2P*Z3-ETA1P*ETA2P*Z4
     &           +ETA1M*ETA2M*Z5+ETA1P*ETA2M*Z6+ETA1M*ETA2P*Z7+ETA1P*ETA2P*Z8

            AJAC=(Aj11*Aj22*Aj33+Aj12*Aj23*Aj31+Aj13*Aj21*Aj32-
     &            Aj11*Aj23*Aj32-Aj12*Aj21*Aj33-Aj13*Aj22*Aj31)
            AJAC=1.953125D-3*AJAC*gp%OME(I)*gp%OME(J)*gp%OME(K)
c                  1/512=1.953125D-3

            CALL cshape27(gp%GI(I),gp%GI(J),gp%GI(K),fig27,27)
c
            DO P=1,27
              CellVolume=CellVolume+AJAC*FIG27(P)
            END DO
c
          END DO !I
        END DO !J
      END DO !K
c
      END

c
C----------------------------------------------------------------------C
c
      SUBROUTINE IntCellVolume8(gp,
     &    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &    x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &    CellVolume)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  INTEgration over Domain cell                                    **
c **  ----             -                                              **
c **  internal cell geometry : LINEAR ( 8-node )                      **
c **  integrira determinanto jacobijeve matrike = izracuna volumen    **
c **                                                                  **
c **********************************************************************
      USE inc_types
      TYPE(gausstype) :: gp
c
      INTEGER i,j,k,p,ng1,ng2
c
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &        x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &        eta1m,eta1p,eta2m,eta2p,eta3m,eta3p,
     &        aj11,aj12,aj13,aj21,aj22,aj23,aj31,aj32,aj33,
     &        ajac

      REAL(8) CellVolume
      REAL(8) fig8(8)
C
C*** SET NUMBER OF INTEGRATION POINTS
      ng1=gp%ng1(2)  ! 2 pomeni pomoje dovolj gaussovnih tock
      ng2=gp%ng2(2)

      CellVolume=0.0D0
C
C***  8 NODE CONTINUOUS DOMAIN CELL (8 NODE GEOMETRY)
C
      DO K=ng1,ng2
        DO J=ng1,ng2
          DO I=ng1,ng2
            ETA1M=1.0D0-gp%GI(I)
            ETA1P=1.0D0+gp%GI(I)
            ETA2M=1.0D0-gp%GI(J)
            ETA2P=1.0D0+gp%GI(J)
            ETA3M=1.0D0-gp%GI(K)
            ETA3P=1.0D0+gp%GI(K)
C
C***  Jacobian derivatives
C
            AJ11=-ETA2M*ETA3M*X1+ETA2M*ETA3M*X2-ETA2P*ETA3M*X3+ETA2P*ETA3M*X4
     &           -ETA2M*ETA3P*X5+ETA2M*ETA3P*X6-ETA2P*ETA3P*X7+ETA2P*ETA3P*X8
            AJ12=-ETA2M*ETA3M*Y1+ETA2M*ETA3M*Y2-ETA2P*ETA3M*Y3+ETA2P*ETA3M*Y4
     &           -ETA2M*ETA3P*Y5+ETA2M*ETA3P*Y6-ETA2P*ETA3P*Y7+ETA2P*ETA3P*Y8
            AJ13=-ETA2M*ETA3M*Z1+ETA2M*ETA3M*Z2-ETA2P*ETA3M*Z3+ETA2P*ETA3M*Z4
     &           -ETA2M*ETA3P*Z5+ETA2M*ETA3P*Z6-ETA2P*ETA3P*Z7+ETA2P*ETA3P*Z8
            AJ21=-ETA1M*ETA3M*X1-ETA1P*ETA3M*X2+ETA1M*ETA3M*X3+ETA1P*ETA3M*X4
     &           -ETA1M*ETA3P*X5-ETA1P*ETA3P*X6+ETA1M*ETA3P*X7+ETA1P*ETA3P*X8
            AJ22=-ETA1M*ETA3M*Y1-ETA1P*ETA3M*Y2+ETA1M*ETA3M*Y3+ETA1P*ETA3M*Y4
     &           -ETA1M*ETA3P*Y5-ETA1P*ETA3P*Y6+ETA1M*ETA3P*Y7+ETA1P*ETA3P*Y8
            AJ23=-ETA1M*ETA3M*Z1-ETA1P*ETA3M*Z2+ETA1M*ETA3M*Z3+ETA1P*ETA3M*Z4
     &           -ETA1M*ETA3P*Z5-ETA1P*ETA3P*Z6+ETA1M*ETA3P*Z7+ETA1P*ETA3P*Z8
            AJ31=-ETA1M*ETA2M*X1-ETA1P*ETA2M*X2-ETA1M*ETA2P*X3-ETA1P*ETA2P*X4
     &           +ETA1M*ETA2M*X5+ETA1P*ETA2M*X6+ETA1M*ETA2P*X7+ETA1P*ETA2P*X8
            AJ32=-ETA1M*ETA2M*Y1-ETA1P*ETA2M*Y2-ETA1M*ETA2P*Y3-ETA1P*ETA2P*Y4
     &           +ETA1M*ETA2M*Y5+ETA1P*ETA2M*Y6+ETA1M*ETA2P*Y7+ETA1P*ETA2P*Y8
            AJ33=-ETA1M*ETA2M*Z1-ETA1P*ETA2M*Z2-ETA1M*ETA2P*Z3-ETA1P*ETA2P*Z4
     &           +ETA1M*ETA2M*Z5+ETA1P*ETA2M*Z6+ETA1M*ETA2P*Z7+ETA1P*ETA2P*Z8

            AJAC=(Aj11*Aj22*Aj33+Aj12*Aj23*Aj31+Aj13*Aj21*Aj32-
     &            Aj11*Aj23*Aj32-Aj12*Aj21*Aj33-Aj13*Aj22*Aj31)
            AJAC=1.953125D-3*AJAC*gp%OME(I)*gp%OME(J)*gp%OME(K)
c                  1/512=1.953125D-3

            CALL cshape8(gp%GI(I),gp%GI(J),gp%GI(K),fig8,8)
c
            DO P=1,8
              CellVolume=CellVolume+AJAC*FIG8(P)
            END DO
c
          END DO !I
        END DO !J
      END DO !K
c
      END

C -----------------------------------------------------------------------------            
      SUBROUTINE INTEBqFlux(gp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,nsipexg,ge)
C
C     $: 2D integracija elementa s 4 tockovno geometrijo
C        v katerm je nezvezno podana funkcija
C        integrale izracunam za vse interpolacijske funkcije
C
C -----------------------------------------------------------------------------      
      USE inc_types 
      TYPE(gausstype) :: gp
      
      INTEGER i,j,isip,nsipexg,ng1,ng2
c
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &        xc0,yc0,zc0,xet,yet,zet,xks,yks,zks,
     &        ajac,
     &        eta1m,eta2m,eta1p,eta2p,
     &        anx,any,anz
      REAL(8) ge(nsipexg)
      REAL(8) fig(4),fig4(4)      
c
c      integral divison
      REAL(8) a,b,c,d,dex,gii,gij
c
C
C*** SET NUMBER OF INTEGRATION POINTS
C

c     regular      
      ng1=gp%ng1(4)
      ng2=gp%ng2(4)
c
      ge=0.00D00
      
      a=-1.0D0 !+(idivXi-1)*2.0D0/ndivXi
      b=1.0D0 !-1.0D0+(idivXi)*2.0D0/ndivXi          
      c=-1.0D0 !+(idivEt-1)*2.0D0/ndivEt
      d=1.0D0 !-1.0D0+(idivEt)*2.0D0/ndivEt
      dex=0.25D0*(b-a)*(d-c)        
C
        DO i=ng1,ng2
          gii=a+0.5D0*(gp%GI(I)+1.0D0)*(b-a)
          ETA1M=1.0D0-gii
          ETA1P=1.0D0+gii
          DO j=ng1,ng2
            gij=c+0.5D0*(gp%GI(J)+1.0D0)*(d-c)                   
            ETA2M=1.0D0-gij
            ETA2P=1.0D0+gij
C
            FIG4(1)=0.25D0*ETA1M*ETA2M
            FIG4(2)=0.25D0*ETA1P*ETA2M
            FIG4(3)=0.25D0*ETA1P*ETA2P
            FIG4(4)=0.25D0*ETA1M*ETA2P
c
            XC0=FIG4(1)*X1+FIG4(2)*X2+FIG4(3)*X3+FIG4(4)*X4
            YC0=FIG4(1)*Y1+FIG4(2)*Y2+FIG4(3)*Y3+FIG4(4)*Y4
            ZC0=FIG4(1)*Z1+FIG4(2)*Z2+FIG4(3)*Z3+FIG4(4)*Z4
c
            XKS=0.25D0*(-ETA2M*X1+ETA2M*X2+ETA2P*X3-ETA2P*X4)
            YKS=0.25D0*(-ETA2M*Y1+ETA2M*Y2+ETA2P*Y3-ETA2P*Y4)
            ZKS=0.25D0*(-ETA2M*Z1+ETA2M*Z2+ETA2P*Z3-ETA2P*Z4)
C
            XET=0.25D0*(-ETA1M*X1-ETA1P*X2+ETA1P*X3+ETA1M*X4)
            YET=0.25D0*(-ETA1M*Y1-ETA1P*Y2+ETA1P*Y3+ETA1M*Y4)
            ZET=0.25D0*(-ETA1M*Z1-ETA1P*Z2+ETA1P*Z3+ETA1M*Z4)
c
            ANX=YKS*ZET-YET*ZKS
            ANY=XET*ZKS-XKS*ZET
            ANZ=XKS*YET-XET*YKS
c
            AJAC=dex*gp%OME(I)*gp%OME(J)*SQRT(ANX**2+ANY**2+ANZ**2)

c           calculate G dis-continous shape functions
            CALL dl34shape4(fig,gii,gij,nsipexg)
            DO isip=1,nsipexg
              ge(isip)=ge(isip)+AJAC*FIG(isip)
            END DO
            
          END DO
        END DO
      
      END


C -----------------------------------------------------------------------------
      SUBROUTINE INTEBfunc(gp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,nsipexg,ge)
C
C     $: 2D integracija elementa s 4 tockovno geometrijo
C        v katerem je zvezno podana funkcija
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(gausstype) :: gp

      INTEGER i,j,isip,nsipexg,ng1,ng2
c
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &        xc0,yc0,zc0,xet,yet,zet,xks,yks,zks,
     &        ajac,
     &        eta1m,eta2m,eta1p,eta2p,
     &        anx,any,anz
      REAL(8) ge(nsipexg)
      REAL(8) fig(9),fig4(4)
c
c      integral divison
      REAL(8) a,b,c,d,dex,gii,gij
c
C
C*** SET NUMBER OF INTEGRATION POINTS
C

c     regular
      ng1=gp%ng1(4)
      ng2=gp%ng2(4)
c
      ge=0.00D00

      a=-1.0D0 !+(idivXi-1)*2.0D0/ndivXi
      b=1.0D0 !-1.0D0+(idivXi)*2.0D0/ndivXi
      c=-1.0D0 !+(idivEt-1)*2.0D0/ndivEt
      d=1.0D0 !-1.0D0+(idivEt)*2.0D0/ndivEt
      dex=0.25D0*(b-a)*(d-c)
C
        DO i=ng1,ng2
          gii=a+0.5D0*(gp%GI(I)+1.0D0)*(b-a)
          ETA1M=1.0D0-gii
          ETA1P=1.0D0+gii
          DO j=ng1,ng2
            gij=c+0.5D0*(gp%GI(J)+1.0D0)*(d-c)
            ETA2M=1.0D0-gij
            ETA2P=1.0D0+gij
C
            FIG4(1)=0.25D0*ETA1M*ETA2M
            FIG4(2)=0.25D0*ETA1P*ETA2M
            FIG4(3)=0.25D0*ETA1P*ETA2P
            FIG4(4)=0.25D0*ETA1M*ETA2P
c
            XC0=FIG4(1)*X1+FIG4(2)*X2+FIG4(3)*X3+FIG4(4)*X4
            YC0=FIG4(1)*Y1+FIG4(2)*Y2+FIG4(3)*Y3+FIG4(4)*Y4
            ZC0=FIG4(1)*Z1+FIG4(2)*Z2+FIG4(3)*Z3+FIG4(4)*Z4
c
            XKS=0.25D0*(-ETA2M*X1+ETA2M*X2+ETA2P*X3-ETA2P*X4)
            YKS=0.25D0*(-ETA2M*Y1+ETA2M*Y2+ETA2P*Y3-ETA2P*Y4)
            ZKS=0.25D0*(-ETA2M*Z1+ETA2M*Z2+ETA2P*Z3-ETA2P*Z4)
C
            XET=0.25D0*(-ETA1M*X1-ETA1P*X2+ETA1P*X3+ETA1M*X4)
            YET=0.25D0*(-ETA1M*Y1-ETA1P*Y2+ETA1P*Y3+ETA1M*Y4)
            ZET=0.25D0*(-ETA1M*Z1-ETA1P*Z2+ETA1P*Z3+ETA1M*Z4)
c
            ANX=YKS*ZET-YET*ZKS
            ANY=XET*ZKS-XKS*ZET
            ANZ=XKS*YET-XET*YKS
c
            AJAC=dex*gp%OME(I)*gp%OME(J)*SQRT(ANX**2+ANY**2+ANZ**2)

c           calculate shape functions
            CALL cshape9(fig,gii,gij,nsipexg)
            DO isip=1,nsipexg
              ge(isip)=ge(isip)+AJAC*FIG(isip)
            END DO

          END DO
        END DO

      END

C -----------------------------------------------------------------------------
      SUBROUTINE corMin_kmdcSteady(env,io,inp,mesh,gauss,smatH,smatG,smatB,
     &                             smatAbdx,smatAbdy,smatAbdz,
     &                             smatHtx,smatHty,smatHtz,
     &                             smatDx,smatDy,smatDz)
C
C     $: Compute or Read Matrices, Steady diffusion advection equation and kinematics eqaution
C        (H, G, Ab, Adx, Ady, Adz, B, Htx, Hty, Htz, Dx, Dy, Dz)
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(IOtype) io
      TYPE(inputtype) inp 
      TYPE(penv) :: env             
      
      REAL(8) smatH(mesh%nicnsp,mesh%npoc),smatG(mesh%nicnsp,mesh%npofc) ! diffusion
      REAL(8) smatB(mesh%nicnsp,mesh%npoc)  ! source    
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary - domain part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary - domain part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary - domain part Z
C
C     kinematics matrices
C      
      REAL(8) smatHtx(mesh%nicnpoc,mesh%npoc) ! doktorat, enacba (4.7)
      REAL(8) smatHty(mesh%nicnpoc,mesh%npoc) ! to je za prvi integral 
      REAL(8) smatHtz(mesh%nicnpoc,mesh%npoc) ! na desni
      
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc) ! doktorat, enacba (4.7)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc) ! to je za drugi integral 
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc) ! na desni       
      REAL(8) ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz ! ocena natancnosti integralov
      REAL(8) ontHtyDx,ontHtzDy,ontHtxDz

      INTEGER iok
      
C     Check if integrals exist on file
      CALL CheckMacroIntFile(io,inp%INTversion,mesh%nnodes,mesh%nbnodes,iok)

      IF (iok.EQ.1) THEN
c       Read integrals from disk 
        CALL WarnErr(env,io,inp,0,"corMin_kmdcSteady","Reading macro integrals!",0)   
        CALL ReadMacroIntDisk(io,mesh,smatH,smatG,smatB,
     &                             smatAbdx,smatAbdy,smatAbdz,
     &                             smatHtx,smatHty,smatHtz,
     &                             smatDx,smatDy,smatDz,
     &                             ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtyDx,ontHtzDy,ontHtxDz)

      ELSE
c       Calcualte integrals      
        CALL WarnErr(env,io,inp,0,"corMin_kmdcSteady","Calculating macro integrals!",0)    
         
        CALL FMATkmdcSteady(env,io,inp,mesh,gauss,smatH,smatG,smatB,
     &                          smatAbdx,smatAbdy,smatAbdz,
     &                          smatHtx,smatHty,smatHtz,
     &                          smatDx,smatDy,smatDz,
     &                          ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtyDx,ontHtzDy,ontHtxDz)

        IF (env%myproc.EQ.1) THEN ! to ni paralelizirano, samo en zapise
          CALL WriteMacroIntDisk(io,inp,mesh,smatH,smatG,smatB,
     &                             smatAbdx,smatAbdy,smatAbdz,
     &                             smatHtx,smatHty,smatHtz,
     &                             smatDx,smatDy,smatDz,
     &                             ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtyDx,ontHtzDy,ontHtxDz)     
        END IF
      END IF

      IF (env%myproc.EQ.1) THEN
        WRITE(io%l,'(A)') "Ocena natancnosti makro integralov!"
        WRITE(io%l,'(A,G15.10)') "H (u i.t.)= ",ontH 
        WRITE(io%l,'(A,G15.10)') "H (q i.t.)= ",ontG       
        WRITE(io%l,'(A,G15.10)') "Htx= ",ontHtx 
        WRITE(io%l,'(A,G15.10)') "Hty= ",ontHty
        WRITE(io%l,'(A,G15.10)') "Htz= ",ontHtz
        WRITE(io%l,'(A,G15.10)') "Dx= ",ontDx 
        WRITE(io%l,'(A,G15.10)') "Dy= ",ontDy
        WRITE(io%l,'(A,G15.10)') "Dz= ",ontDz
        WRITE(io%l,'(A)') "Po popravljanju singularnih integralov:"
        WRITE(io%l,'(A,G15.10)') "HtxDz = ",ontHtxDz
        WRITE(io%l,'(A,G15.10)') "HtyDx = ",ontHtyDx
        WRITE(io%l,'(A,G15.10)') "HtzDy = ",ontHtzDy
        WRITE(io%l,'(A)') ""
      END IF      

      END
C -----------------------------------------------------------------------------
      SUBROUTINE WriteMacroIntDisk(io,inp,mesh,smatH,smatG,smatB,
     &                             smatAbdx,smatAbdy,smatAbdz,
     &                             smatHtx,smatHty,smatHtz,
     &                             smatDx,smatDy,smatDz,
     &                             ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtyDx,ontHtzDy,ontHtxDz)       
C
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      TYPE(IOtype) io
      TYPE(inputtype) inp        
      
      REAL(8) smatH(mesh%nicnsp,mesh%npoc),smatG(mesh%nicnsp,mesh%npofc) ! diffusion
      REAL(8) smatB(mesh%nicnsp,mesh%npoc)  ! source    
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary - domain part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary - domain part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary - domain part Z  
      REAL(8) smatHtx(mesh%nicnpoc,mesh%npoc) ! doktorat, enacba (4.7)
      REAL(8) smatHty(mesh%nicnpoc,mesh%npoc) ! to je za prvi integral 
      REAL(8) smatHtz(mesh%nicnpoc,mesh%npoc) ! na desni 
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc) ! doktorat, enacba (4.7)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc) ! to je za drugi integral 
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc) ! na desni      
      REAL(8) ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz ! ocena natancnosti integralov   
      REAL(8) ontHtyDx,ontHtzDy,ontHtxDz         

      OPEN (io%imi,FILE=TRIM(io%imi_name),FORM='UNFORMATTED',STATUS='UNKNOWN')      
      WRITE(io%imi) inp%INTversion
      WRITE(io%imi) mesh%nnodes,mesh%nbnodes
     
      CALL WrSMat(smatH,mesh%nicnsp,mesh%npoc,io%imi)
      CALL WrSMat(smatG,mesh%nicnsp,mesh%npofc,io%imi)
      CALL WrSMat(smatB,mesh%nicnsp,mesh%npoc,io%imi)      
      CALL WrSMat(smatAbdX,mesh%nicnsp,mesh%npoc,io%imi)
      CALL WrSMat(smatAbdY,mesh%nicnsp,mesh%npoc,io%imi)
      CALL WrSMat(smatAbdZ,mesh%nicnsp,mesh%npoc,io%imi)
      CALL WrSMat(smatHtx,mesh%nicnpoc,mesh%npoc,io%imi)
      CALL WrSMat(smatHty,mesh%nicnpoc,mesh%npoc,io%imi)
      CALL WrSMat(smatHtz,mesh%nicnpoc,mesh%npoc,io%imi)
      CALL WrSMat(smatDx,mesh%nicnsp,mesh%npoc,io%imi)
      CALL WrSMat(smatDy,mesh%nicnsp,mesh%npoc,io%imi)
      CALL WrSMat(smatDz,mesh%nicnsp,mesh%npoc,io%imi)  
      WRITE(io%imi) ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtyDx,ontHtzDy,ontHtxDz                                              
      

      CLOSE(io%imi)          
      
      END

C -----------------------------------------------------------------------------
      SUBROUTINE ReadMacroIntDisk(io,mesh,smatH,smatG,smatB,
     &                             smatAbdx,smatAbdy,smatAbdz,
     &                             smatHtx,smatHty,smatHtz,
     &                             smatDx,smatDy,smatDz,
     &                             ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtyDx,ontHtzDy,ontHtxDz)  
C
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      TYPE(IOtype) io
      INTEGER a,b,c    
      
      REAL(8) smatH(mesh%nicnsp,mesh%npoc),smatG(mesh%nicnsp,mesh%npofc) ! diffusion
      REAL(8) smatB(mesh%nicnsp,mesh%npoc)  ! source    
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary - domain part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary - domain part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary - domain part Z  
      REAL(8) smatHtx(mesh%nicnpoc,mesh%npoc) ! doktorat, enacba (4.7)
      REAL(8) smatHty(mesh%nicnpoc,mesh%npoc) ! to je za prvi integral 
      REAL(8) smatHtz(mesh%nicnpoc,mesh%npoc) ! na desni 
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc) ! doktorat, enacba (4.7)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc) ! to je za drugi integral 
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc) ! na desni      
      REAL(8) ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz ! ocena natancnosti integralov 
      REAL(8) ontHtyDx,ontHtzDy,ontHtxDz           

      OPEN (io%imi,FILE=TRIM(io%imi_name),FORM='UNFORMATTED',STATUS='OLD')      
      READ(io%imi) a
      READ(io%imi) b,c
      
      CALL RdSMat(smatH,mesh%nicnsp,mesh%npoc,io%imi)
      CALL RdSMat(smatG,mesh%nicnsp,mesh%npofc,io%imi)
      CALL RdSMat(smatB,mesh%nicnsp,mesh%npoc,io%imi)       
      CALL RdSMat(smatAbdX,mesh%nicnsp,mesh%npoc,io%imi)
      CALL RdSMat(smatAbdY,mesh%nicnsp,mesh%npoc,io%imi)
      CALL RdSMat(smatAbdZ,mesh%nicnsp,mesh%npoc,io%imi)
      CALL RdSMat(smatHtx,mesh%nicnpoc,mesh%npoc,io%imi)
      CALL RdSMat(smatHty,mesh%nicnpoc,mesh%npoc,io%imi)
      CALL RdSMat(smatHtz,mesh%nicnpoc,mesh%npoc,io%imi)
      CALL RdSMat(smatDx,mesh%nicnsp,mesh%npoc,io%imi)
      CALL RdSMat(smatDy,mesh%nicnsp,mesh%npoc,io%imi)
      CALL RdSMat(smatDz,mesh%nicnsp,mesh%npoc,io%imi)                                                
      READ (io%imi) ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtyDx,ontHtzDy,ontHtxDz  

      CLOSE(io%imi)          
      
      END      
      
C______________________________________________________________________C
      SUBROUTINE CheckMacroIntFile(io,a,b,c,iok)
      USE inc_types
      TYPE(IOtype) io 
      INTEGER a,b,c
      INTEGER aa,bb,cc,iok

      iok=0
      OPEN (io%imi,FILE=TRIM(io%imi_name),FORM='UNFORMATTED',STATUS='OLD',ERR=10)
      READ(io%imi) aa
      READ(io%imi) bb,cc
            
      IF (a.EQ.aa.AND.b.EQ.bb.AND.c.EQ.cc) iok=1    
      CLOSE(io%imi)
      
10    RETURN
      END
            
C -----------------------------------------------------------------------------            
      SUBROUTINE INTEBc9kmdclap(gp,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                   isrc,nsipex,nsipexg,he,ge,adx,ady,adz,
     &                   heyz,hezx,hexy,minedge)
C
C     $: Integracija Robni 9 tockovni element s 4 tockovno geometrijo
C        za steady state difuzivno konvektivno enacbo + enacbo kinematike
C
C -----------------------------------------------------------------------------      
      USE inc_types 
      TYPE(gausstype) :: gp
      
      INTEGER i,j,k,isrc,isip,nsipex,nsipexg,ng1s,ng2s,ng1,ng2
c
      REAL(8) xp,yp,zp,
     &        x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &        xc0,yc0,zc0,xet,yet,zet,xks,yks,zks,
     &        xx1,yy1,
     &        rx,ry,rz,ra,ira2,
     &        ro,roth,th,ajac,
     &        pi,pi2,
     &        eti,etj,eta1m,eta2m,eta1p,eta2p,
     &        anx,any,anz,anx1,any1,anz1,
     &        fung,funx,funy,funz,funh,
     &        dgamma
      REAL(8) ge(nsipexg),he(nsipex),adx(nsipex),ady(nsipex),adz(nsipex)
      REAL(8) hexy(nsipex),heyz(nsipex),hezx(nsipex)
c      REAL(8) fig(nsipexg),fih(nsipex),fig4(4)
      REAL(8) fig(4),fih(9),fig4(4)      
      REAL(8) al(4),fii(4),th0(4),th1(4),ksi(13),eta(13)
c
c      integral divison
      REAL(8) a,b,c,d,dex,gii,gij,d1,d2,d3,d4,d13max,d24max,minedge
      INTEGER idivXi,ndivXi,idivEt,ndivEt,ising

c
      DATA ksi / 0.0D00, 0.5D00, 1.0D00, 1.0D00, 1.0D00, 0.5D00, 0.0D00, 0.0D00, 0.5D00, 0.125D0, 0.875D0, 0.875D0, 0.125D0/
      DATA eta / 0.0D00, 0.0D00, 0.0D00, 0.5D00, 1.0D00, 1.0D00, 1.0D00, 0.5D00, 0.5D00, 0.125D0, 0.125D0, 0.875D0, 0.875D0/
C
C*** SET NUMBER OF INTEGRATION POINTS
C
c     singular      
      ng1s=gp%ng1(gp%kmBs)
      ng2s=gp%ng2(gp%kmBs)
c     regular      
      ng1=gp%ng1(gp%kmBr)
      ng2=gp%ng2(gp%kmBr)
C
C***  9 NODE CONTINUOUS BOUNDARY ELEMENT
C
      PI=2.0D0*ASIN(1.0D0)
      PI2=2.0D0*PI
c
      he=0.00D00
      ge=0.00D00
      adx=0.0D0
      ady=0.0D0
      adz=0.0D0            
      hexy=0.00D00
      heyz=0.00D00
      hezx=0.00D00

C
C     Integral razdelimo tako, da integriramo po kvdadratkih
C
      d1=SQRT((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
      d2=SQRT((x2-x3)**2+(y2-y3)**2+(z2-z3)**2)
      d3=SQRT((x3-x4)**2+(y3-y4)**2+(z3-z4)**2)
      d4=SQRT((x4-x1)**2+(y4-y1)**2+(z4-z1)**2)   
      d13max=max(d1,d3)
      d24max=max(d2,d4)      

      ndivXi=max(gp%iidiv,INT(d13max/minedge+0.5D0))
      ndivEt=max(gp%iidiv,INT(d24max/minedge+0.5D0))
c
c     glavna integracijska zanka
c      
      DO idivXi=1,ndivXi
        DO idivEt=1,ndivEt
          a=-1.0D0+(idivXi-1)*2.0D0/ndivXi
          b=-1.0D0+(idivXi)*2.0D0/ndivXi          
          c=-1.0D0+(idivEt-1)*2.0D0/ndivEt
          d=-1.0D0+(idivEt)*2.0D0/ndivEt
          dex=0.25D0*(b-a)*(d-c)        
c          print *,b-a,d-c
          IF (isrc.NE.0) THEN
            XX1=ksi(isrc)  ! med 0 in 1
            XX1=-1.0D0+2.0D0*XX1  ! med -1 in 1
            XX1=(XX1-a)/(b-a)  ! med 0 in 1 v intervalu a,b
      
            YY1=eta(isrc)  ! med 0 in 1
            YY1=-1.0D0+2.0D0*YY1  ! med -1 in 1
            YY1=(YY1-c)/(d-c)  ! med 0 in 1 v intervalu c,d
          
            IF (XX1.GE.0.0D0.AND.XX1.LE.1.0D0.AND.YY1.GE.0.0D0.AND.YY1.LE.1.0D0) THEN
              ising=1
            ELSE
              ising=0
            END IF
          ELSE
            ising=0
          ENDIF         
c
C*** SINGULAR INTEGRALS
c
      IF (ising.NE.0) THEN
c        XX1=ksi(isrc)
c        YY1=eta(isrc)
c
        DO K=1,4
c
          IF(K.EQ.1) THEN
            IF (1.0D0-XX1.EQ.0.0D0) GOTO 1000
            AL(K)=SQRT(YY1**2+(1.0D0-XX1)**2)
            FII(K)=ACOS(YY1/AL(K))
            TH0(K)=1.5D0*PI+ACOS(YY1/AL(K))
            TH1(K)=PI2+ATAN((1.0D0-YY1)/(1.0D0-XX1))
c
          ELSE IF(K.EQ.2) THEN
            IF (1.0D0-YY1.EQ.0.0D0) GOTO 1000
            AL(K)=SQRT((1.0D0-YY1)**2+(1.0D0-XX1)**2)
            FII(K)=ACOS((1.0D0-XX1)/AL(K))
            TH0(K)=ASIN((1.0D0-YY1)/AL(K))
            TH1(K)=0.5D0*PI+ATAN(XX1/(1.0D0-YY1))
c
          ELSE IF(K.EQ.3) THEN
            IF (XX1.EQ.0.0D0) GOTO 1000
            AL(K)=SQRT(XX1**2+(1.0D0-YY1)**2)
            FII(K)=ACOS((1.0D0-YY1)/AL(K))
            TH0(K)=0.5D0*PI+ACOS((1.0D0-YY1)/AL(K))
            TH1(K)=PI+ATAN(YY1/XX1)
c
          ELSE
            IF (YY1.EQ.0.0D0) GOTO 1000
            AL(K)=SQRT(YY1**2+XX1**2)
            FII(K)=ACOS(XX1/AL(K))
            TH0(K)=PI+ACOS(XX1/AL(K))
            TH1(K)=1.5D0*PI+ATAN((1.0D0-XX1)/YY1)
          END IF
C
C***      GAUSS INTEGRATION (48 points)
C
          DO I=ng1s,ng2s
            DO J=ng1s,ng2s
c
              TH=(TH1(K)-TH0(K))*gp%GI(I)/2.0D0+(TH1(K)+TH0(K))/2.0D0
              ROTH=AL(K)*SIN(FII(K))/SIN(TH-TH0(K)+FII(K))
              RO=ROTH*gp%GI(J)/2.0D0+ROTH/2.0D0
c
              ETI=2.0D0*(XX1+RO*COS(TH))-1.0D0
              ETJ=2.0D0*(YY1+RO*SIN(TH))-1.0D0
c
              ETI=a+0.5D0*(ETI+1.0D0)*(b-a)     
              ETJ=c+0.5D0*(ETJ+1.0D0)*(d-c)                   
c
              ETA1M=1.D0-ETI
              ETA1P=1.D0+ETI
              ETA2M=1.D0-ETJ
              ETA2P=1.D0+ETJ
 
              FIG4(1)=0.25D0*ETA1M*ETA2M
              FIG4(2)=0.25D0*ETA1P*ETA2M
              FIG4(3)=0.25D0*ETA1P*ETA2P
              FIG4(4)=0.25D0*ETA1M*ETA2P
c
              XC0=FIG4(1)*X1+FIG4(2)*X2+FIG4(3)*X3+FIG4(4)*X4
              YC0=FIG4(1)*Y1+FIG4(2)*Y2+FIG4(3)*Y3+FIG4(4)*Y4
              ZC0=FIG4(1)*Z1+FIG4(2)*Z2+FIG4(3)*Z3+FIG4(4)*Z4
c
              XKS=0.25D0*(-ETA2M*X1+ETA2M*X2+ETA2P*X3-ETA2P*X4)
              YKS=0.25D0*(-ETA2M*Y1+ETA2M*Y2+ETA2P*Y3-ETA2P*Y4)
              ZKS=0.25D0*(-ETA2M*Z1+ETA2M*Z2+ETA2P*Z3-ETA2P*Z4)
C
              XET=0.25D0*(-ETA1M*X1-ETA1P*X2+ETA1P*X3+ETA1M*X4)
              YET=0.25D0*(-ETA1M*Y1-ETA1P*Y2+ETA1P*Y3+ETA1M*Y4)
              ZET=0.25D0*(-ETA1M*Z1-ETA1P*Z2+ETA1P*Z3+ETA1M*Z4)
c
              ANX=YKS*ZET-YET*ZKS
              ANY=XET*ZKS-XKS*ZET
              ANZ=XKS*YET-XET*YKS
c
              AJAC=1.0D0/SQRT(ANX**2+ANY**2+ANZ**2)
c
              ANX1=ANX*AJAC
              ANY1=ANY*AJAC
              ANZ1=ANZ*AJAC
c
              AJAC=dex*(TH1(K)-TH0(K))*ROTH*RO*gp%OME(I)*gp%OME(J)/AJAC
c
              RX=XP-XC0
              RY=YP-YC0
              RZ=ZP-ZC0
              RA=SQRT(RX*RX+RY*RY+RZ*RZ)
c             IF (RA.LT.1.0D-12) RA=1.0D-12
              IRA2=1.0D00/(RA*RA)
C
C***          ELLIPTIC LAPLACE FUNDAMENTAL SOLUTION
C 
              FUNG=0.25D0/(PI*RA)
              FUNX=(XP-XC0)*IRA2*FUNG
              FUNY=(YP-YC0)*IRA2*FUNG
              FUNZ=(ZP-ZC0)*IRA2*FUNG
              FUNH=FUNX*ANX1+FUNY*ANY1+FUNZ*ANZ1

c             calculate H continous shape functions
              CALL cshape9(fih,ETI,ETJ,nsipex)
              DO isip=1,nsipex
                DGAMMA=AJAC*FIH(isip)
                he(isip)=he(isip)+FUNH*DGAMMA
                adx(isip)=adx(isip)+FUNG*ANX1*DGAMMA
                ady(isip)=ady(isip)+FUNG*ANY1*DGAMMA
                adz(isip)=adz(isip)+FUNG*ANZ1*DGAMMA                                
                hexy(isip)=hexy(isip)-(FUNX*ANY1-FUNY*ANX1)*DGAMMA
                heyz(isip)=heyz(isip)-(FUNY*ANZ1-FUNZ*ANY1)*DGAMMA
                hezx(isip)=hezx(isip)-(FUNZ*ANX1-FUNX*ANZ1)*DGAMMA
              END DO

c             calculate G dis-continous shape functions
              CALL dl34shape4(fig,ETI,ETJ,nsipexg)
              DO isip=1,nsipexg
                DGAMMA=AJAC*FIG(isip)
                ge(isip)=ge(isip)+FUNG*DGAMMA
              END DO

            END DO
          END DO
 1000   END DO
c      
      ELSE
C
C*** REGULAR INTEGRALS
C
        DO i=ng1,ng2
          gii=a+0.5D0*(gp%GI(I)+1.0D0)*(b-a)
          ETA1M=1.0D0-gii
          ETA1P=1.0D0+gii
          DO j=ng1,ng2
            gij=c+0.5D0*(gp%GI(J)+1.0D0)*(d-c)                   
            ETA2M=1.0D0-gij
            ETA2P=1.0D0+gij
C
            FIG4(1)=0.25D0*ETA1M*ETA2M
            FIG4(2)=0.25D0*ETA1P*ETA2M
            FIG4(3)=0.25D0*ETA1P*ETA2P
            FIG4(4)=0.25D0*ETA1M*ETA2P
c
            XC0=FIG4(1)*X1+FIG4(2)*X2+FIG4(3)*X3+FIG4(4)*X4
            YC0=FIG4(1)*Y1+FIG4(2)*Y2+FIG4(3)*Y3+FIG4(4)*Y4
            ZC0=FIG4(1)*Z1+FIG4(2)*Z2+FIG4(3)*Z3+FIG4(4)*Z4
c
            XKS=0.25D0*(-ETA2M*X1+ETA2M*X2+ETA2P*X3-ETA2P*X4)
            YKS=0.25D0*(-ETA2M*Y1+ETA2M*Y2+ETA2P*Y3-ETA2P*Y4)
            ZKS=0.25D0*(-ETA2M*Z1+ETA2M*Z2+ETA2P*Z3-ETA2P*Z4)
C
            XET=0.25D0*(-ETA1M*X1-ETA1P*X2+ETA1P*X3+ETA1M*X4)
            YET=0.25D0*(-ETA1M*Y1-ETA1P*Y2+ETA1P*Y3+ETA1M*Y4)
            ZET=0.25D0*(-ETA1M*Z1-ETA1P*Z2+ETA1P*Z3+ETA1M*Z4)
c
            ANX=YKS*ZET-YET*ZKS
            ANY=XET*ZKS-XKS*ZET
            ANZ=XKS*YET-XET*YKS
c
            AJAC=1.0D0/SQRT(ANX**2+ANY**2+ANZ**2)
c
            ANX1=ANX*AJAC
            ANY1=ANY*AJAC
            ANZ1=ANZ*AJAC
c
            AJAC=dex*gp%OME(I)*gp%OME(J)/AJAC
c
            RX=XP-XC0
            RY=YP-YC0
            RZ=ZP-ZC0
            RA=SQRT(RX*RX+RY*RY+RZ*RZ)
c           IF (RA.LT.1.0D-12) RA=1.0D-12
            IRA2=1.0D00/(RA*RA)
C
C***        ELLIPTIC LAPLACE FUNDAMENTAL SOLUTION
C     
            FUNG=0.25D0/(PI*RA)
            FUNX=(XP-XC0)*IRA2*FUNG
            FUNY=(YP-YC0)*IRA2*FUNG
            FUNZ=(ZP-ZC0)*IRA2*FUNG
            FUNH=FUNX*ANX1+FUNY*ANY1+FUNZ*ANZ1
c
c           calculate H continous shape functions
            CALL cshape9(fih,gii,gij,nsipex)            
            DO isip=1,nsipex
              DGAMMA=AJAC*FIh(isip)
              he(isip)=he(isip)+FUNH*DGAMMA
              adx(isip)=adx(isip)+FUNG*ANX1*DGAMMA
              ady(isip)=ady(isip)+FUNG*ANY1*DGAMMA
              adz(isip)=adz(isip)+FUNG*ANZ1*DGAMMA                
              hexy(isip)=hexy(isip)-(FUNX*ANY1-FUNY*ANX1)*DGAMMA
              heyz(isip)=heyz(isip)-(FUNY*ANZ1-FUNZ*ANY1)*DGAMMA
              hezx(isip)=hezx(isip)-(FUNZ*ANX1-FUNX*ANZ1)*DGAMMA
            END DO

c             calculate G dis-continous shape functions
            CALL dl34shape4(fig,gii,gij,nsipexg)
            DO isip=1,nsipexg
              DGAMMA=AJAC*FIG(isip)
              ge(isip)=ge(isip)+FUNG*DGAMMA
            END DO
            
          END DO
        END DO
        CONTINUE
      END IF
      
        END DO
      END DO        
      
c            
      
      END

C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE INTEDc27dcLapOld
     &(
     &    gp,xp,yp,zp,
     &    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &    x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &    be,dxe,dye,dze,isrc,nsipmx
     &)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  INTEgration over Domain cell                                    **
c **  ----             -                                              **
c **  internal cell geometry : LINEAR ( 8-node )                      **
c **  internal cell function interpolation : QUADRATIC ( 27-node )    **
c **  fundamental solution   : ELLIPTIC LAPLACE                       **
c **                                                                  **
c **********************************************************************
      USE inc_types 
      TYPE(gausstype) :: gp
c
      INTEGER i,j,k,l,n,p,nsipmx,isrc,ng1s,ng2s,ng1,ng2
c      
      REAL(8) xp,yp,zp,
     &        x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &        x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &        xc0,yc0,zc0,
     &        xx1,yy1,zz1,
     &        pi,pi2,
     &        aa0,gama,ffi0,ffi1,ffi,ro,ro0,ro1,th,th0,th1,aa1,rr1,eps,
     &        eti,etj,etk,
     &        eta1m,eta1p,eta2m,eta2p,eta3m,eta3p,
     &        aj11,aj12,aj13,aj21,aj22,aj23,aj31,aj32,aj33,
     &        ajac,rx,ry,rz,ra,ira2,
     &        fung,funx,funy,funz,
     &        ksi,eta,zet,domega
c
      REAL(8) fig27(nsipmx),fig(8)
      REAL(8) dxe(nsipmx),dye(nsipmx),dze(nsipmx)
      REAL(8) be(nsipmx)
      DIMENSION ksi(51),eta(51),zet(51)
c
      DATA ksi /0.0D+0,0.0D+0,1.0D+0,1.0D+0,0.0D+0,0.0D+0,1.0D+0,1.0D+0,0.0D+0,
     &          5.0D-1,1.0D+0,5.0D-1,0.0D+0,0.0D+0,1.0D+0,1.0D+0,0.0D+0,5.0D-1,
     &          1.0D+0,5.0D-1,5.0D-1,0.0D+0,5.0D-1,1.0D+0,5.0D-1,5.0D-1,5.0D-1,
     &          0.125D0,0.125D0,0.875D0,0.875D0,0.125D0,0.875D0,0.875D0,0.125D0,
     &          1.000D0,1.000D0,1.000D0,1.000D0,0.875D0,0.125D0,0.125D0,0.875D0,
     &          0.000D0,0.000D0,0.000D0,0.000D0,0.125D0,0.875D0,0.875D0,0.125D0/
      DATA eta /1.0D+0,1.0D+0,1.0D+0,1.0D+0,0.0D+0,0.0D+0,0.0D+0,0.0D+0,1.0D+0,
     &          1.0D+0,1.0D+0,1.0D+0,5.0D-1,5.0D-1,5.0D-1,5.0D-1,0.0D+0,0.0D+0,
     &          0.0D+0,0.0D+0,1.0D+0,5.0D-1,5.0D-1,5.0D-1,5.0D-1,0.0D+0,5.0D-1,
     &          0.125D0,0.875D0,0.875D0,0.125D0,0.000D0,0.000D0,0.000D0,0.000D0,
     &          0.125D0,0.875D0,0.875D0,0.125D0,1.000D0,1.000D0,1.000D0,1.000D0,
     &          0.875D0,0.125D0,0.125D0,0.875D0,0.125D0,0.125D0,0.875D0,0.875D0/
      DATA zet /1.0D+0,0.0D+0,0.0D+0,1.0D+0,1.0D+0,0.0D+0,0.0D+0,1.0D+0,5.0D-1,
     &          0.0D+0,5.0D-1,1.0D+0,1.0D+0,0.0D+0,0.0D+0,1.0D+0,5.0D-1,0.0D+0,
     &          5.0D-1,1.0D+0,5.0D-1,5.0D-1,0.0D+0,5.0D-1,1.0D+0,5.0D-1,5.0D-1,
     &          0.000D0,0.000D0,0.000D0,0.000D0,0.125D0,0.125D0,0.875D0,0.875D0,
     &          0.125D0,0.125D0,0.875D0,0.875D0,0.125D0,0.125D0,0.875D0,0.875D0,
     &          0.125D0,0.125D0,0.875D0,0.875D0,1.000D0,1.000D0,1.000D0,1.000D0/          
C
C*** SET NUMBER OF INTEGRATION POINTS
C
c     singular      
      ng1s=gp%ng1(gp%kmDs)
      ng2s=gp%ng2(gp%kmDs)
c     regular      
      ng1=gp%ng1(gp%kmDr)
      ng2=gp%ng2(gp%kmDr)
C
C***  27 NODE CONTINUOUS DOMAIN CELL (8 NODE GEOMETRY)
C
      PI=2.*ASIN(1.0D0)
      PI2=PI*2.0D0
c
      dxe=0.00D00
      dye=0.00D00
      dze=0.00D00
      be=0.00D00
c
C*** SINGULAR INTEGRALS
c
      IF (isrc.NE.0) THEN
        XX1=ksi(isrc)
        YY1=eta(isrc)
        ZZ1=zet(isrc)
c
        DO L=1,4
c
          IF (L.EQ.1) THEN
            IF (1.0D0-XX1.EQ.0.0D0) GOTO 1000
            aA0=SQRT(YY1**2+(1.0D0-XX1)**2)
            gama=ACOS(YY1/aA0)
            ffi0=1.5D0*PI+ACOS(YY1/aa0)
            ffi1=PI2+ATAN((1.0D0-YY1)/(1.0D0-XX1))
          ELSE IF (L.EQ.2) THEN
            IF (1.0D0-YY1.EQ.0.0D0) GOTO 1000
            aA0=SQRT((1.0D0-YY1)**2+(1.0D0-XX1)**2)
            gama=ACOS((1.0D0-XX1)/aA0)
            ffi0=ASIN((1.0D0-YY1)/aA0)
            ffi1=0.5D0*PI+ATAN(XX1/(1.0D0-YY1))
          ELSE IF (L.EQ.3) THEN
            IF (XX1.EQ.0.0D0) GOTO 1000
            aA0=SQRT(XX1**2+(1.0D0-YY1)**2)
            gama=ACOS((1.0D0-YY1)/aA0)
            ffi0=0.5D0*PI+ACOS((1.0D0-YY1)/aA0)
            ffi1=PI+ATAN(YY1/XX1)
          ELSE
            IF (YY1.EQ.0.0D0) GOTO 1000
            aA0=SQRT(YY1**2+XX1**2)
            gama=ACOS(XX1/aA0)
            ffi0=PI+ACOS(XX1/aA0)
            ffi1=1.5D0*PI+ATAN((1.0D0-XX1)/YY1)
          END IF
c        
          DO K=ng1s,ng2s
            ffi=(ffi1-ffi0)*gp%GI(k)/2.0D0+(ffi1+ffi0)/2.0D0
            RO0=Aa0*SIN(gama)/SIN(ffi-ffi0+gama)
c
            DO N=1,3
              IF (N.eq.1) THEN
                aa1=1.0D0-zz1
                eps=pi/2.0D0
                th0=0.0D0
                th1=ACOS(aa1/SQRT(aa1**2+ro0**2))
              ELSE IF (N.eq.2) THEN
                aa1=SQRT(ro0**2+(1.0D0-zz1)**2)
                eps=ACOS((1.0D0-zz1)/aa1)
                th0=eps
                rr1=MIN(1.0D00,ro0/SQRT(zz1**2+ro0**2))
                th1=0.5D0*pi+ACOS(rr1)
              ELSE
                aa1=SQRT(ro0**2+zz1**2)
                rr1=MIN(1.0D00,ro0/aa1)
                eps=ACOS(rr1)
                th0=0.5D0*pi+eps
                th1=pi
              END IF
c
              DO I=ng1s,ng2s
                th=(th1-th0)*gp%GI(I)/2.0D0+(th1+th0)/2.0D0
                RO1=Aa1*SIN(eps)/SIN(th-th0+eps)
c
                DO J=ng1s,ng2s
                    ro=(ro1*gp%gi(j)+ro1)/2.0D0
c
                  eti=2.0D0*(xx1+ro*SIN(th)*COS(ffi))-1.0D0
                  etj=2.0D0*(yy1+ro*SIN(th)*SIN(ffi))-1.0D0
                  etk=2.0D0*(zz1+ro*COS(th))-1.0D0

                 ETA1M=1.0D0-ETI   
                 ETA1P=1.0D0+ETI   
                 ETA2M=1.0D0-ETJ   
                 ETA2P=1.0D0+ETJ   
                 ETA3M=1.0D0-ETK   
                 ETA3P=1.0D0+ETK   

                  FIG(1)=0.125D0*ETA1M*ETA2M*ETA3M
                  FIG(2)=0.125D0*ETA1P*ETA2M*ETA3M
                  FIG(3)=0.125D0*ETA1M*ETA2P*ETA3M
                  FIG(4)=0.125D0*ETA1P*ETA2P*ETA3M
                  FIG(5)=0.125D0*ETA1M*ETA2M*ETA3P
                  FIG(6)=0.125D0*ETA1P*ETA2M*ETA3P
                  FIG(7)=0.125D0*ETA1M*ETA2P*ETA3P
                  FIG(8)=0.125D0*ETA1P*ETA2P*ETA3P

                  XC0=FIG(1)*X1+FIG(2)*X2+FIG(3)*X3+FIG(4)*X4
     &               +FIG(5)*X5+FIG(6)*X6+FIG(7)*X7+FIG(8)*X8
                  YC0=FIG(1)*Y1+FIG(2)*Y2+FIG(3)*Y3+FIG(4)*Y4
     &               +FIG(5)*Y5+FIG(6)*Y6+FIG(7)*Y7+FIG(8)*Y8
                  ZC0=FIG(1)*Z1+FIG(2)*Z2+FIG(3)*Z3+FIG(4)*Z4
     &               +FIG(5)*Z5+FIG(6)*Z6+FIG(7)*Z7+FIG(8)*Z8
C
C***  Jacobian derivatives
C
              AJ11=-ETA2M*ETA3M*X1+ETA2M*ETA3M*X2-ETA2P*ETA3M*X3+ETA2P*ETA3M*X4
     &                 -ETA2M*ETA3P*X5+ETA2M*ETA3P*X6-ETA2P*ETA3P*X7+ETA2P*ETA3P*X8

                 AJ12=-ETA2M*ETA3M*Y1+ETA2M*ETA3M*Y2-ETA2P*ETA3M*Y3+ETA2P*ETA3M*Y4
     &                 -ETA2M*ETA3P*Y5+ETA2M*ETA3P*Y6-ETA2P*ETA3P*Y7+ETA2P*ETA3P*Y8

                 AJ13=-ETA2M*ETA3M*Z1+ETA2M*ETA3M*Z2-ETA2P*ETA3M*Z3+ETA2P*ETA3M*Z4
     &                 -ETA2M*ETA3P*Z5+ETA2M*ETA3P*Z6-ETA2P*ETA3P*Z7+ETA2P*ETA3P*Z8
c     
                 AJ21=-ETA1M*ETA3M*X1-ETA1P*ETA3M*X2+ETA1M*ETA3M*X3+ETA1P*ETA3M*X4
     &                 -ETA1M*ETA3P*X5-ETA1P*ETA3P*X6+ETA1M*ETA3P*X7+ETA1P*ETA3P*X8

                 AJ22=-ETA1M*ETA3M*Y1-ETA1P*ETA3M*Y2+ETA1M*ETA3M*Y3+ETA1P*ETA3M*Y4
     &                 -ETA1M*ETA3P*Y5-ETA1P*ETA3P*Y6+ETA1M*ETA3P*Y7+ETA1P*ETA3P*Y8

                 AJ23=-ETA1M*ETA3M*Z1-ETA1P*ETA3M*Z2+ETA1M*ETA3M*Z3+ETA1P*ETA3M*Z4
     &                 -ETA1M*ETA3P*Z5-ETA1P*ETA3P*Z6+ETA1M*ETA3P*Z7+ETA1P*ETA3P*Z8
c
                AJ31=-ETA1M*ETA2M*X1-ETA1P*ETA2M*X2-ETA1M*ETA2P*X3-ETA1P*ETA2P*X4
     &                 +ETA1M*ETA2M*X5+ETA1P*ETA2M*X6+ETA1M*ETA2P*X7+ETA1P*ETA2P*X8

                 AJ32=-ETA1M*ETA2M*Y1-ETA1P*ETA2M*Y2-ETA1M*ETA2P*Y3-ETA1P*ETA2P*Y4
     &                 +ETA1M*ETA2M*Y5+ETA1P*ETA2M*Y6+ETA1M*ETA2P*Y7+ETA1P*ETA2P*Y8

                 AJ33=-ETA1M*ETA2M*Z1-ETA1P*ETA2M*Z2-ETA1M*ETA2P*Z3-ETA1P*ETA2P*Z4
     &                 +ETA1M*ETA2M*Z5+ETA1P*ETA2M*Z6+ETA1M*ETA2P*Z7+ETA1P*ETA2P*Z8
c
                  AJAC=(Aj11*Aj22*Aj33+Aj12*Aj23*Aj31+Aj13*Aj21*Aj32-
     &                  Aj11*Aj23*Aj32-Aj12*Aj21*Aj33-Aj13*Aj22*Aj31)
c
                  AJAC=1.953125D-3*AJAC*gp%OME(I)*gp%OME(J)*gp%OME(K)*(ffi1-ffi0)*(th1-th0)*ro1*ro**2*SIN(th)

                  RX=XP-XC0
                  RY=YP-YC0
                  RZ=ZP-ZC0
                  RA=SQRT(RX*RX+RY*RY+RZ*RZ)
                  IF (RA.EQ.0.0D00) GO TO 800
                  IF (RA.LT.1.0D-12) RA=1.0D-12
                  IRA2=1.0D00/(RA*RA)
C
C***              ELLIPTIC LAPLACE FUNDAMENTAL SOLUTION
C
                  FUNG=0.25D0/(PI*RA)
                  FUNX=(XP-XC0)*IRA2*FUNG
                  FUNY=(YP-YC0)*IRA2*FUNG
                  FUNZ=(ZP-ZC0)*IRA2*FUNG
c
                  CALL cshape27(eti,etj,etk,fig27,nsipmx)
c
                  DO P=1,nsipmx
                    DOMEGA=AJAC*FIG27(P)
                    dxe(p)=dxe(p)+FUNX*DOMEGA
                    dye(p)=dye(p)+FUNY*DOMEGA
                    dze(p)=dze(p)+FUNZ*DOMEGA
                    be(p) =be(p) +FUNG*DOMEGA
                  END DO      
c
  800           END DO !J
              END DO !I
            END DO !N
          END DO !K
 1000   END DO !L
      ELSE
        DO K=ng1,ng2
          DO J=ng1,ng2
            DO I=ng1,ng2

            ETA1M=1.0D0-gp%GI(I)
            ETA1P=1.0D0+gp%GI(I)
            ETA2M=1.0D0-gp%GI(J)
            ETA2P=1.0D0+gp%GI(J)
            ETA3M=1.0D0-gp%GI(K)
            ETA3P=1.0D0+gp%GI(K)

              FIG(1)=0.125D0*ETA1M*ETA2M*ETA3M
              FIG(2)=0.125D0*ETA1P*ETA2M*ETA3M
              FIG(3)=0.125D0*ETA1M*ETA2P*ETA3M
              FIG(4)=0.125D0*ETA1P*ETA2P*ETA3M
              FIG(5)=0.125D0*ETA1M*ETA2M*ETA3P
              FIG(6)=0.125D0*ETA1P*ETA2M*ETA3P
              FIG(7)=0.125D0*ETA1M*ETA2P*ETA3P
              FIG(8)=0.125D0*ETA1P*ETA2P*ETA3P

              XC0=FIG(1)*X1+FIG(2)*X2+FIG(3)*X3+FIG(4)*X4
     &           +FIG(5)*X5+FIG(6)*X6+FIG(7)*X7+FIG(8)*X8
              YC0=FIG(1)*Y1+FIG(2)*Y2+FIG(3)*Y3+FIG(4)*Y4
     &           +FIG(5)*Y5+FIG(6)*Y6+FIG(7)*Y7+FIG(8)*Y8
              ZC0=FIG(1)*Z1+FIG(2)*Z2+FIG(3)*Z3+FIG(4)*Z4
     &           +FIG(5)*Z5+FIG(6)*Z6+FIG(7)*Z7+FIG(8)*Z8
C
C***  Jacobian derivatives
C
           AJ11=-ETA2M*ETA3M*X1+ETA2M*ETA3M*X2-ETA2P*ETA3M*X3+ETA2P*ETA3M*X4
     &             -ETA2M*ETA3P*X5+ETA2M*ETA3P*X6-ETA2P*ETA3P*X7+ETA2P*ETA3P*X8

         AJ12=-ETA2M*ETA3M*Y1+ETA2M*ETA3M*Y2-ETA2P*ETA3M*Y3+ETA2P*ETA3M*Y4
     &             -ETA2M*ETA3P*Y5+ETA2M*ETA3P*Y6-ETA2P*ETA3P*Y7+ETA2P*ETA3P*Y8

         AJ13=-ETA2M*ETA3M*Z1+ETA2M*ETA3M*Z2-ETA2P*ETA3M*Z3+ETA2P*ETA3M*Z4
     &             -ETA2M*ETA3P*Z5+ETA2M*ETA3P*Z6-ETA2P*ETA3P*Z7+ETA2P*ETA3P*Z8
c 
             AJ21=-ETA1M*ETA3M*X1-ETA1P*ETA3M*X2+ETA1M*ETA3M*X3+ETA1P*ETA3M*X4
     &             -ETA1M*ETA3P*X5-ETA1P*ETA3P*X6+ETA1M*ETA3P*X7+ETA1P*ETA3P*X8

         AJ22=-ETA1M*ETA3M*Y1-ETA1P*ETA3M*Y2+ETA1M*ETA3M*Y3+ETA1P*ETA3M*Y4
     &             -ETA1M*ETA3P*Y5-ETA1P*ETA3P*Y6+ETA1M*ETA3P*Y7+ETA1P*ETA3P*Y8

         AJ23=-ETA1M*ETA3M*Z1-ETA1P*ETA3M*Z2+ETA1M*ETA3M*Z3+ETA1P*ETA3M*Z4
     &             -ETA1M*ETA3P*Z5-ETA1P*ETA3P*Z6+ETA1M*ETA3P*Z7+ETA1P*ETA3P*Z8
c
              AJ31=-ETA1M*ETA2M*X1-ETA1P*ETA2M*X2-ETA1M*ETA2P*X3-ETA1P*ETA2P*X4
     &             +ETA1M*ETA2M*X5+ETA1P*ETA2M*X6+ETA1M*ETA2P*X7+ETA1P*ETA2P*X8

         AJ32=-ETA1M*ETA2M*Y1-ETA1P*ETA2M*Y2-ETA1M*ETA2P*Y3-ETA1P*ETA2P*Y4
     &             +ETA1M*ETA2M*Y5+ETA1P*ETA2M*Y6+ETA1M*ETA2P*Y7+ETA1P*ETA2P*Y8

         AJ33=-ETA1M*ETA2M*Z1-ETA1P*ETA2M*Z2-ETA1M*ETA2P*Z3-ETA1P*ETA2P*Z4
     &             +ETA1M*ETA2M*Z5+ETA1P*ETA2M*Z6+ETA1M*ETA2P*Z7+ETA1P*ETA2P*Z8
c
              AJAC=(Aj11*Aj22*Aj33+Aj12*Aj23*Aj31+Aj13*Aj21*Aj32-
     &              Aj11*Aj23*Aj32-Aj12*Aj21*Aj33-Aj13*Aj22*Aj31)
c
              AJAC=1.953125D-3*AJAC*gp%OME(I)*gp%OME(J)*gp%OME(K)
c                  1/512=1.953125D-3
              RX=XP-XC0
              RY=YP-YC0
              RZ=ZP-ZC0
              RA=SQRT(RX*RX+RY*RY+RZ*RZ)
              IF (RA.LT.1.0D-12) RA=1.0D-12 
              IRA2=1.0D00/(RA*RA)
C
C***          ELLIPTIC LAPLACE FUNDAMENTAL SOLUTION
C
              FUNG=0.25D0/(PI*RA)
              FUNX=(XP-XC0)*IRA2*FUNG
              FUNY=(YP-YC0)*IRA2*FUNG
              FUNZ=(ZP-ZC0)*IRA2*FUNG
c              
              CALL cshape27(gp%GI(I),gp%GI(J),gp%GI(K),fig27,nsipmx)
c
              DO P=1,nsipmx
                DOMEGA=AJAC*FIG27(P)
                dxe(p)=dxe(p)+FUNX*DOMEGA
                dye(p)=dye(p)+FUNY*DOMEGA
                dze(p)=dze(p)+FUNZ*DOMEGA
                be(p) =be(p) +FUNG*DOMEGA
              END DO      
c
            END DO !I
          END DO !J
        END DO !K
      END IF
c
      RETURN
      END



C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE INTEDc8dcLap
     &(
     &    gp,xp,yp,zp,
     &    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &    x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &    be,dxe,dye,dze,isrc,nsipmx
     &)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  INTEgration over Domain cell                                    **
c **  ----             -                                              **
c **  internal cell geometry : LINEAR ( 8-node )                      **
c **  internal cell function interpolation : LINEAR ( 8-node )        **
c **  fundamental solution   : ELLIPTIC LAPLACE                       **
c **                                                                  **
c **********************************************************************
      USE inc_types
      TYPE(gausstype) :: gp
c
      INTEGER i,j,k,l,n,p,nsipmx,isrc,ng1s,ng2s,ng1,ng2
c
      REAL(8) xp,yp,zp,
     &        x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &        x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &        xc0,yc0,zc0,
     &        xx1,yy1,zz1,
     &        pi,pi2,
     &        aa0,gama,ffi0,ffi1,ffi,ro,ro0,ro1,th,th0,th1,aa1,rr1,eps,
     &        eti,etj,etk,
     &        eta1m,eta1p,eta2m,eta2p,eta3m,eta3p,
     &        aj11,aj12,aj13,aj21,aj22,aj23,aj31,aj32,aj33,
     &        ajac,rx,ry,rz,ra,ira2,
     &        fung,funx,funy,funz,
     &        ksi,eta,zet,domega

c     nsipmx=8
      REAL(8) fig(nsipmx)
      REAL(8) dxe(nsipmx),dye(nsipmx),dze(nsipmx)
      REAL(8) be(nsipmx)
      DIMENSION ksi(51),eta(51),zet(51)
c
      DATA ksi /0.0D+0,0.0D+0,1.0D+0,1.0D+0,0.0D+0,0.0D+0,1.0D+0,1.0D+0,0.0D+0,
     &          5.0D-1,1.0D+0,5.0D-1,0.0D+0,0.0D+0,1.0D+0,1.0D+0,0.0D+0,5.0D-1,
     &          1.0D+0,5.0D-1,5.0D-1,0.0D+0,5.0D-1,1.0D+0,5.0D-1,5.0D-1,5.0D-1,
     &          0.125D0,0.125D0,0.875D0,0.875D0,0.125D0,0.875D0,0.875D0,0.125D0,
     &          1.000D0,1.000D0,1.000D0,1.000D0,0.875D0,0.125D0,0.125D0,0.875D0,
     &          0.000D0,0.000D0,0.000D0,0.000D0,0.125D0,0.875D0,0.875D0,0.125D0/
      DATA eta /1.0D+0,1.0D+0,1.0D+0,1.0D+0,0.0D+0,0.0D+0,0.0D+0,0.0D+0,1.0D+0,
     &          1.0D+0,1.0D+0,1.0D+0,5.0D-1,5.0D-1,5.0D-1,5.0D-1,0.0D+0,0.0D+0,
     &          0.0D+0,0.0D+0,1.0D+0,5.0D-1,5.0D-1,5.0D-1,5.0D-1,0.0D+0,5.0D-1,
     &          0.125D0,0.875D0,0.875D0,0.125D0,0.000D0,0.000D0,0.000D0,0.000D0,
     &          0.125D0,0.875D0,0.875D0,0.125D0,1.000D0,1.000D0,1.000D0,1.000D0,
     &          0.875D0,0.125D0,0.125D0,0.875D0,0.125D0,0.125D0,0.875D0,0.875D0/
      DATA zet /1.0D+0,0.0D+0,0.0D+0,1.0D+0,1.0D+0,0.0D+0,0.0D+0,1.0D+0,5.0D-1,
     &          0.0D+0,5.0D-1,1.0D+0,1.0D+0,0.0D+0,0.0D+0,1.0D+0,5.0D-1,0.0D+0,
     &          5.0D-1,1.0D+0,5.0D-1,5.0D-1,0.0D+0,5.0D-1,1.0D+0,5.0D-1,5.0D-1,
     &          0.000D0,0.000D0,0.000D0,0.000D0,0.125D0,0.125D0,0.875D0,0.875D0,
     &          0.125D0,0.125D0,0.875D0,0.875D0,0.125D0,0.125D0,0.875D0,0.875D0,
     &          0.125D0,0.125D0,0.875D0,0.875D0,1.000D0,1.000D0,1.000D0,1.000D0/
C
C*** SET NUMBER OF INTEGRATION POINTS
C
c     singular
      ng1s=gp%ng1(gp%kmDs)
      ng2s=gp%ng2(gp%kmDs)
c     regular
      ng1=gp%ng1(gp%kmDr)
      ng2=gp%ng2(gp%kmDr)
C
C***  8 NODE CONTINUOUS DOMAIN CELL (8 NODE GEOMETRY)
C
      PI=2.*ASIN(1.0D0)
      PI2=PI*2.0D0
c
      dxe=0.00D00
      dye=0.00D00
      dze=0.00D00
      be=0.00D00
c
C*** SINGULAR INTEGRALS
c
      IF (isrc.NE.0) THEN
        XX1=ksi(isrc)
        YY1=eta(isrc)
        ZZ1=zet(isrc)
c
        DO L=1,4
c
          IF (L.EQ.1) THEN
            IF (1.0D0-XX1.EQ.0.0D0) GOTO 1000
            aA0=SQRT(YY1**2+(1.0D0-XX1)**2)
            gama=ACOS(YY1/aA0)
            ffi0=1.5D0*PI+ACOS(YY1/aa0)
            ffi1=PI2+ATAN((1.0D0-YY1)/(1.0D0-XX1))
          ELSE IF (L.EQ.2) THEN
            IF (1.0D0-YY1.EQ.0.0D0) GOTO 1000
            aA0=SQRT((1.0D0-YY1)**2+(1.0D0-XX1)**2)
            gama=ACOS((1.0D0-XX1)/aA0)
            ffi0=ASIN((1.0D0-YY1)/aA0)
            ffi1=0.5D0*PI+ATAN(XX1/(1.0D0-YY1))
          ELSE IF (L.EQ.3) THEN
            IF (XX1.EQ.0.0D0) GOTO 1000
            aA0=SQRT(XX1**2+(1.0D0-YY1)**2)
            gama=ACOS((1.0D0-YY1)/aA0)
            ffi0=0.5D0*PI+ACOS((1.0D0-YY1)/aA0)
            ffi1=PI+ATAN(YY1/XX1)
          ELSE
            IF (YY1.EQ.0.0D0) GOTO 1000
            aA0=SQRT(YY1**2+XX1**2)
            gama=ACOS(XX1/aA0)
            ffi0=PI+ACOS(XX1/aA0)
            ffi1=1.5D0*PI+ATAN((1.0D0-XX1)/YY1)
          END IF
c
          DO K=ng1s,ng2s
            ffi=(ffi1-ffi0)*gp%GI(k)/2.0D0+(ffi1+ffi0)/2.0D0
            RO0=Aa0*SIN(gama)/SIN(ffi-ffi0+gama)
c
            DO N=1,3
              IF (N.eq.1) THEN
                aa1=1.0D0-zz1
                eps=pi/2.0D0
                th0=0.0D0
                th1=ACOS(aa1/SQRT(aa1**2+ro0**2))
              ELSE IF (N.eq.2) THEN
                aa1=SQRT(ro0**2+(1.0D0-zz1)**2)
                eps=ACOS((1.0D0-zz1)/aa1)
                th0=eps
                rr1=MIN(1.0D00,ro0/SQRT(zz1**2+ro0**2))
                th1=0.5D0*pi+ACOS(rr1)
              ELSE
                aa1=SQRT(ro0**2+zz1**2)
                rr1=MIN(1.0D00,ro0/aa1)
                eps=ACOS(rr1)
                th0=0.5D0*pi+eps
                th1=pi
              END IF
c
              DO I=ng1s,ng2s
                th=(th1-th0)*gp%GI(I)/2.0D0+(th1+th0)/2.0D0
                RO1=Aa1*SIN(eps)/SIN(th-th0+eps)
c
                DO J=ng1s,ng2s
                    ro=(ro1*gp%gi(j)+ro1)/2.0D0
c
                  eti=2.0D0*(xx1+ro*SIN(th)*COS(ffi))-1.0D0
                  etj=2.0D0*(yy1+ro*SIN(th)*SIN(ffi))-1.0D0
                  etk=2.0D0*(zz1+ro*COS(th))-1.0D0

                 ETA1M=1.0D0-ETI
                 ETA1P=1.0D0+ETI
                 ETA2M=1.0D0-ETJ
                 ETA2P=1.0D0+ETJ
                 ETA3M=1.0D0-ETK
                 ETA3P=1.0D0+ETK

                  FIG(1)=0.125D0*ETA1M*ETA2M*ETA3M
                  FIG(2)=0.125D0*ETA1P*ETA2M*ETA3M
                  FIG(3)=0.125D0*ETA1M*ETA2P*ETA3M
                  FIG(4)=0.125D0*ETA1P*ETA2P*ETA3M
                  FIG(5)=0.125D0*ETA1M*ETA2M*ETA3P
                  FIG(6)=0.125D0*ETA1P*ETA2M*ETA3P
                  FIG(7)=0.125D0*ETA1M*ETA2P*ETA3P
                  FIG(8)=0.125D0*ETA1P*ETA2P*ETA3P

                  XC0=FIG(1)*X1+FIG(2)*X2+FIG(3)*X3+FIG(4)*X4
     &               +FIG(5)*X5+FIG(6)*X6+FIG(7)*X7+FIG(8)*X8
                  YC0=FIG(1)*Y1+FIG(2)*Y2+FIG(3)*Y3+FIG(4)*Y4
     &               +FIG(5)*Y5+FIG(6)*Y6+FIG(7)*Y7+FIG(8)*Y8
                  ZC0=FIG(1)*Z1+FIG(2)*Z2+FIG(3)*Z3+FIG(4)*Z4
     &               +FIG(5)*Z5+FIG(6)*Z6+FIG(7)*Z7+FIG(8)*Z8
C
C***  Jacobian derivatives
C
              AJ11=-ETA2M*ETA3M*X1+ETA2M*ETA3M*X2-ETA2P*ETA3M*X3+ETA2P*ETA3M*X4
     &                 -ETA2M*ETA3P*X5+ETA2M*ETA3P*X6-ETA2P*ETA3P*X7+ETA2P*ETA3P*X8

                 AJ12=-ETA2M*ETA3M*Y1+ETA2M*ETA3M*Y2-ETA2P*ETA3M*Y3+ETA2P*ETA3M*Y4
     &                 -ETA2M*ETA3P*Y5+ETA2M*ETA3P*Y6-ETA2P*ETA3P*Y7+ETA2P*ETA3P*Y8

                 AJ13=-ETA2M*ETA3M*Z1+ETA2M*ETA3M*Z2-ETA2P*ETA3M*Z3+ETA2P*ETA3M*Z4
     &                 -ETA2M*ETA3P*Z5+ETA2M*ETA3P*Z6-ETA2P*ETA3P*Z7+ETA2P*ETA3P*Z8
c
                 AJ21=-ETA1M*ETA3M*X1-ETA1P*ETA3M*X2+ETA1M*ETA3M*X3+ETA1P*ETA3M*X4
     &                 -ETA1M*ETA3P*X5-ETA1P*ETA3P*X6+ETA1M*ETA3P*X7+ETA1P*ETA3P*X8

                 AJ22=-ETA1M*ETA3M*Y1-ETA1P*ETA3M*Y2+ETA1M*ETA3M*Y3+ETA1P*ETA3M*Y4
     &                 -ETA1M*ETA3P*Y5-ETA1P*ETA3P*Y6+ETA1M*ETA3P*Y7+ETA1P*ETA3P*Y8

                 AJ23=-ETA1M*ETA3M*Z1-ETA1P*ETA3M*Z2+ETA1M*ETA3M*Z3+ETA1P*ETA3M*Z4
     &                 -ETA1M*ETA3P*Z5-ETA1P*ETA3P*Z6+ETA1M*ETA3P*Z7+ETA1P*ETA3P*Z8
c
                AJ31=-ETA1M*ETA2M*X1-ETA1P*ETA2M*X2-ETA1M*ETA2P*X3-ETA1P*ETA2P*X4
     &                 +ETA1M*ETA2M*X5+ETA1P*ETA2M*X6+ETA1M*ETA2P*X7+ETA1P*ETA2P*X8

                 AJ32=-ETA1M*ETA2M*Y1-ETA1P*ETA2M*Y2-ETA1M*ETA2P*Y3-ETA1P*ETA2P*Y4
     &                 +ETA1M*ETA2M*Y5+ETA1P*ETA2M*Y6+ETA1M*ETA2P*Y7+ETA1P*ETA2P*Y8

                 AJ33=-ETA1M*ETA2M*Z1-ETA1P*ETA2M*Z2-ETA1M*ETA2P*Z3-ETA1P*ETA2P*Z4
     &                 +ETA1M*ETA2M*Z5+ETA1P*ETA2M*Z6+ETA1M*ETA2P*Z7+ETA1P*ETA2P*Z8
c
                  AJAC=(Aj11*Aj22*Aj33+Aj12*Aj23*Aj31+Aj13*Aj21*Aj32-
     &                  Aj11*Aj23*Aj32-Aj12*Aj21*Aj33-Aj13*Aj22*Aj31)
c
                  AJAC=1.953125D-3*AJAC*gp%OME(I)*gp%OME(J)*gp%OME(K)*(ffi1-ffi0)*(th1-th0)*ro1*ro**2*SIN(th)

                  RX=XP-XC0
                  RY=YP-YC0
                  RZ=ZP-ZC0
                  RA=SQRT(RX*RX+RY*RY+RZ*RZ)
                  IF (RA.EQ.0.0D00) GO TO 800
                  IF (RA.LT.1.0D-12) RA=1.0D-12
                  IRA2=1.0D00/(RA*RA)
C
C***              ELLIPTIC LAPLACE FUNDAMENTAL SOLUTION
C
                  FUNG=0.25D0/(PI*RA)
                  FUNX=(XP-XC0)*IRA2*FUNG
                  FUNY=(YP-YC0)*IRA2*FUNG
                  FUNZ=(ZP-ZC0)*IRA2*FUNG
c
c                  CALL cshape27(eti,etj,etk,fig27,nsipmx)
c
                  DO P=1,nsipmx
                    DOMEGA=AJAC*FIG(P)
                    dxe(p)=dxe(p)+FUNX*DOMEGA
                    dye(p)=dye(p)+FUNY*DOMEGA
                    dze(p)=dze(p)+FUNZ*DOMEGA
                    be(p) =be(p) +FUNG*DOMEGA
                  END DO
c
  800           END DO !J
              END DO !I
            END DO !N
          END DO !K
 1000   END DO !L
      ELSE
        DO K=ng1,ng2
          DO J=ng1,ng2
            DO I=ng1,ng2

            ETA1M=1.0D0-gp%GI(I)
            ETA1P=1.0D0+gp%GI(I)
            ETA2M=1.0D0-gp%GI(J)
            ETA2P=1.0D0+gp%GI(J)
            ETA3M=1.0D0-gp%GI(K)
            ETA3P=1.0D0+gp%GI(K)

              FIG(1)=0.125D0*ETA1M*ETA2M*ETA3M
              FIG(2)=0.125D0*ETA1P*ETA2M*ETA3M
              FIG(3)=0.125D0*ETA1M*ETA2P*ETA3M
              FIG(4)=0.125D0*ETA1P*ETA2P*ETA3M
              FIG(5)=0.125D0*ETA1M*ETA2M*ETA3P
              FIG(6)=0.125D0*ETA1P*ETA2M*ETA3P
              FIG(7)=0.125D0*ETA1M*ETA2P*ETA3P
              FIG(8)=0.125D0*ETA1P*ETA2P*ETA3P

              XC0=FIG(1)*X1+FIG(2)*X2+FIG(3)*X3+FIG(4)*X4
     &           +FIG(5)*X5+FIG(6)*X6+FIG(7)*X7+FIG(8)*X8
              YC0=FIG(1)*Y1+FIG(2)*Y2+FIG(3)*Y3+FIG(4)*Y4
     &           +FIG(5)*Y5+FIG(6)*Y6+FIG(7)*Y7+FIG(8)*Y8
              ZC0=FIG(1)*Z1+FIG(2)*Z2+FIG(3)*Z3+FIG(4)*Z4
     &           +FIG(5)*Z5+FIG(6)*Z6+FIG(7)*Z7+FIG(8)*Z8
C
C***  Jacobian derivatives
C
           AJ11=-ETA2M*ETA3M*X1+ETA2M*ETA3M*X2-ETA2P*ETA3M*X3+ETA2P*ETA3M*X4
     &             -ETA2M*ETA3P*X5+ETA2M*ETA3P*X6-ETA2P*ETA3P*X7+ETA2P*ETA3P*X8

         AJ12=-ETA2M*ETA3M*Y1+ETA2M*ETA3M*Y2-ETA2P*ETA3M*Y3+ETA2P*ETA3M*Y4
     &             -ETA2M*ETA3P*Y5+ETA2M*ETA3P*Y6-ETA2P*ETA3P*Y7+ETA2P*ETA3P*Y8

         AJ13=-ETA2M*ETA3M*Z1+ETA2M*ETA3M*Z2-ETA2P*ETA3M*Z3+ETA2P*ETA3M*Z4
     &             -ETA2M*ETA3P*Z5+ETA2M*ETA3P*Z6-ETA2P*ETA3P*Z7+ETA2P*ETA3P*Z8
c
             AJ21=-ETA1M*ETA3M*X1-ETA1P*ETA3M*X2+ETA1M*ETA3M*X3+ETA1P*ETA3M*X4
     &             -ETA1M*ETA3P*X5-ETA1P*ETA3P*X6+ETA1M*ETA3P*X7+ETA1P*ETA3P*X8

         AJ22=-ETA1M*ETA3M*Y1-ETA1P*ETA3M*Y2+ETA1M*ETA3M*Y3+ETA1P*ETA3M*Y4
     &             -ETA1M*ETA3P*Y5-ETA1P*ETA3P*Y6+ETA1M*ETA3P*Y7+ETA1P*ETA3P*Y8

         AJ23=-ETA1M*ETA3M*Z1-ETA1P*ETA3M*Z2+ETA1M*ETA3M*Z3+ETA1P*ETA3M*Z4
     &             -ETA1M*ETA3P*Z5-ETA1P*ETA3P*Z6+ETA1M*ETA3P*Z7+ETA1P*ETA3P*Z8
c
              AJ31=-ETA1M*ETA2M*X1-ETA1P*ETA2M*X2-ETA1M*ETA2P*X3-ETA1P*ETA2P*X4
     &             +ETA1M*ETA2M*X5+ETA1P*ETA2M*X6+ETA1M*ETA2P*X7+ETA1P*ETA2P*X8

         AJ32=-ETA1M*ETA2M*Y1-ETA1P*ETA2M*Y2-ETA1M*ETA2P*Y3-ETA1P*ETA2P*Y4
     &             +ETA1M*ETA2M*Y5+ETA1P*ETA2M*Y6+ETA1M*ETA2P*Y7+ETA1P*ETA2P*Y8

         AJ33=-ETA1M*ETA2M*Z1-ETA1P*ETA2M*Z2-ETA1M*ETA2P*Z3-ETA1P*ETA2P*Z4
     &             +ETA1M*ETA2M*Z5+ETA1P*ETA2M*Z6+ETA1M*ETA2P*Z7+ETA1P*ETA2P*Z8
c
              AJAC=(Aj11*Aj22*Aj33+Aj12*Aj23*Aj31+Aj13*Aj21*Aj32-
     &              Aj11*Aj23*Aj32-Aj12*Aj21*Aj33-Aj13*Aj22*Aj31)
c
              AJAC=1.953125D-3*AJAC*gp%OME(I)*gp%OME(J)*gp%OME(K)
c                  1/512=1.953125D-3
              RX=XP-XC0
              RY=YP-YC0
              RZ=ZP-ZC0
              RA=SQRT(RX*RX+RY*RY+RZ*RZ)
              IF (RA.LT.1.0D-12) RA=1.0D-12
              IRA2=1.0D00/(RA*RA)
C
C***          ELLIPTIC LAPLACE FUNDAMENTAL SOLUTION
C
              FUNG=0.25D0/(PI*RA)
              FUNX=(XP-XC0)*IRA2*FUNG
              FUNY=(YP-YC0)*IRA2*FUNG
              FUNZ=(ZP-ZC0)*IRA2*FUNG
c
c              CALL cshape27(gp%GI(I),gp%GI(J),gp%GI(K),fig27,nsipmx)
c
              DO P=1,nsipmx
                DOMEGA=AJAC*FIG(P)
                dxe(p)=dxe(p)+FUNX*DOMEGA
                dye(p)=dye(p)+FUNY*DOMEGA
                dze(p)=dze(p)+FUNZ*DOMEGA
                be(p) =be(p) +FUNG*DOMEGA
              END DO
c
            END DO !I
          END DO !J
        END DO !K
      END IF
c
      RETURN
      END



C -----------------------------------------------------------------------------            
      SUBROUTINE FindMinEdge(mesh,ibc,minedge)      
C
C     $: Poisce najkrajsi rob elementa
C
C -----------------------------------------------------------------------------      
      USE inc_types   

      TYPE(meshType) :: mesh
      REAL(8) minedge,elen
      INTEGER ibc(mesh%nside,mesh%npob)
      INTEGER i

      minedge=1.0D20
 
      DO i=1,mesh%nside
        elen=SQRT((mesh%x(ibc(i,1),1)-mesh%x(ibc(i,3),1))**2+
     &            (mesh%x(ibc(i,1),2)-mesh%x(ibc(i,3),2))**2+
     &            (mesh%x(ibc(i,1),3)-mesh%x(ibc(i,3),3))**2)
        IF (elen.LT.minedge) minedge=elen

        elen=SQRT((mesh%x(ibc(i,3),1)-mesh%x(ibc(i,5),1))**2+
     &            (mesh%x(ibc(i,3),2)-mesh%x(ibc(i,5),2))**2+
     &            (mesh%x(ibc(i,3),3)-mesh%x(ibc(i,5),3))**2)
        IF (elen.LT.minedge) minedge=elen

        elen=SQRT((mesh%x(ibc(i,5),1)-mesh%x(ibc(i,7),1))**2+
     &            (mesh%x(ibc(i,5),2)-mesh%x(ibc(i,7),2))**2+
     &            (mesh%x(ibc(i,5),3)-mesh%x(ibc(i,7),3))**2)
        IF (elen.LT.minedge) minedge=elen

        elen=SQRT((mesh%x(ibc(i,7),1)-mesh%x(ibc(i,1),1))**2+
     &            (mesh%x(ibc(i,7),2)-mesh%x(ibc(i,1),2))**2+
     &            (mesh%x(ibc(i,7),3)-mesh%x(ibc(i,1),3))**2)
        IF (elen.LT.minedge) minedge=elen
      END DO
        
      END

C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsB(eqn,mesh,smatH,smatG,smatB,beta,ren,sysm,rhsm)
C
C     $: Iz pravokotnih matrik g, h in b naredi CRS sistemsko in rhs
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      REAL(8) smatH(mesh%nicnsp,mesh%npoc),smatG(mesh%nicnsp,mesh%npofc)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc),beta,ren
      REAL(8), ALLOCATABLE :: FRsys(:),FRrhs(:)
      TYPE(matrix) :: sysm,rhsm
      
      INTEGER wqkode(mesh%nq)
      INTEGER wkode(mesh%nnodes)
      INTEGER nunk, nb, eqn
            
      INTEGER ic,it,row,col,ii,jj,nuq,j,i,ks,kr
           
      IF (eqn.LE.3) THEN
        wqkode=mesh%wqkode(:,eqn)
        wkode=mesh%wkode(:,eqn)
        nunk=mesh%nunk(eqn)
        nb=mesh%nb(eqn)
      ELSE IF (eqn.EQ.4) THEN
        wqkode=mesh%Tqkode
        wkode=mesh%Tkode
        nunk=mesh%Tnunk
        nb=mesh%Tnb
      ELSE
        wqkode=mesh%DCqkode
        wkode=mesh%DCkode
        nunk=mesh%DCnunk
        nb=mesh%DCnb
      END IF        
C      
C     stevilo enacb
C
      sysm%neq=mesh%nicnsp      
      rhsm%neq=mesh%nicnsp
C
C     stevilo nenicelnih clenov
C      
      sysm%nnz=0
      rhsm%nnz=0
      DO ic=1,mesh%nicell ! zanka po celicah
        DO it=1,mesh%nsp ! po izvornih tockah znotraj celice
          row=(ic-1)*mesh%nsp+it ! vrstica v sistemski matriki        
C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            IF (smatH(row,col)+smatB(row,col)*beta*ren.NE.0.0D0) THEN
              j=mesh%idc(ic,col) ! stolpec za H
              nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                sysm%nnz=sysm%nnz+1
              ELSE ! znana vrednost - rhs
                rhsm%nnz=rhsm%nnz+1
              END IF
            END IF
          END DO            
C
C         nato G matrika (q)
C
          DO ii=1,mesh%nside
            DO jj=1,mesh%npof
              col=(ii-1)*mesh%npof+jj
              IF (smatG(row,col).NE.0.0D0) THEN
                j=mesh%ibf(ic,ii,jj) ! stolpec za G
                nuq=wqkode(j) ! stevilka v vektorju neznank oziroma znank
                IF (nuq.LT.0) THEN ! neznaka - sys
                  sysm%nnz=sysm%nnz+1
                ELSE ! znana vrednost - rhs
                  rhsm%nnz=rhsm%nnz+1
                END IF
              END IF
            END DO
          END DO
C                    
        END DO
      END DO
C
C     Alokacija pomnilnika
C
      ALLOCATE(sysm%v(sysm%nnz),sysm%i(sysm%neq+1),sysm%j(sysm%nnz))
      ALLOCATE(rhsm%v(rhsm%nnz),rhsm%i(rhsm%neq+1),rhsm%j(rhsm%nnz))   
C
C     Polnjene matrik (napolnim polno vrsti, nato v CRS (da mi ni treba sortirat)
C         
      ALLOCATE (FRsys(nunk),FRrhs(nb))
      ks=0
      kr=0
      DO ic=1,mesh%nicell ! zanka po celicah
        DO it=1,mesh%nsp ! po izvornih tockah znotraj celice
          row=(ic-1)*mesh%nsp+it ! vrstica v sistemski matriki
          FRsys=0.0D0
          FRrhs=0.0D0
C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              FRsys(ABS(nuq))=smatH(row,col)+smatB(row,col)*beta*ren
            ELSE ! znana vrednost - rhs
              FRrhs(nuq)=-smatH(row,col)-smatB(row,col)*beta*ren  ! minus zato, ker gre na drugo stran enacbe
            END IF
          END DO            
C
C         nato G matrika (q)
C
          DO ii=1,mesh%nside
            DO jj=1,mesh%npof
              col=(ii-1)*mesh%npof+jj
              j=mesh%ibf(ic,ii,jj) ! stolpec za G
              nuq=wqkode(j) ! stevilka v vektorju neznank oziroma znank
              IF (nuq.LT.0) THEN ! neznaka - sys
                FRsys(ABS(nuq))=-smatG(row,col)
              ELSE ! znana vrednost - rhs
                FRrhs(nuq)=smatG(row,col) ! minus zato, ker gre na drugo stran enacbe
              END IF
            END DO
          END DO
C               
C         predelam v CRS
C     
C         zacetek vrste
          sysm%i(row)=ks+1
          rhsm%i(row)=kr+1          
C         sys
          DO i=1,nunk
            IF (FRsys(i).NE.0.0D0) THEN
              ks=ks+1
              sysm%v(ks)=FRsys(i)
              sysm%j(ks)=i
            END IF
          END DO
C         rhs
          DO i=1,nb
            IF (FRrhs(i).NE.0.0D0) THEN
              kr=kr+1
              rhsm%v(kr)=FRrhs(i)
              rhsm%j(kr)=i
            END IF
          END DO
c          
        END DO
      END DO     
c     zadnji clen, zacetek neobstojece vrste
      sysm%i(row+1)=ks+1
      rhsm%i(row+1)=kr+1          
       
      
      DEALLOCATE (FRsys,FRrhs)
      
      END
      
      
C -----------------------------------------------------------------------------
      SUBROUTINE FindColInIdc(mesh,ic,col,cinidc)
C
C     $: Poisce zaporedno stevilko col v idc-ju
C
C -----------------------------------------------------------------------------
      USE inc_types       
      TYPE(meshType) :: mesh
      INTEGER ic,col,cinidc,i    
      
      cinidc=-1
      DO i=1,mesh%npoc
        IF (col.EQ.mesh%idc(ic,i)) THEN
          cinidc=i
          EXIT
        END IF
      END DO
      
      END
      
      
      
      
C______________________________________________________________________________
      SUBROUTINE GaussPosWeights(gaus)
C
C     Defines positions and weights for gauss quadrature
C     Bronstein 730, 441
C      
      USE inc_types
      TYPE(gausstype) :: gaus
      
      ALLOCATE(gaus%gi(220),gaus%ome(220))
      ALLOCATE(gaus%ng1(9),gaus%ng2(9))
      
C---  N =   2  ( ng1 =  1 , ng2 =  2 )
      gaus%ng1(1)=1
      gaus%ng2(1)=2
      gaus%gi(1)   = -0.5773502691896258D0
      gaus%gi(2)   =  0.5773502691896258D0
      gaus%ome(1)  =  1.0000000000000000D0
      gaus%ome(2)  =  1.0000000000000000D0
      
C---  N =   4 ( ng1 =  3 , ng2 =  6 )
      gaus%ng1(2)=3
      gaus%ng2(2)=6
      gaus%gi(3)   = -0.8611363115940526D0
      gaus%gi(4)   = -0.3399810435848563D0
      gaus%gi(5)   =  0.3399810435848563D0
      gaus%gi(6)   =  0.8611363115940526D0
      gaus%ome(3)  =  0.3478548451374539D0
      gaus%ome(4)  =  0.6521451548625462D0
      gaus%ome(5)  =  0.6521451548625462D0
      gaus%ome(6)  =  0.3478548451374539D0
      
C---  N =   6 ( ng1 =  7 , ng2 = 12 )
      gaus%ng1(3)=7
      gaus%ng2(3)=12
      gaus%gi(7)   = -0.9324695142031520D0
      gaus%gi(8)   = -0.6612093864662645D0
      gaus%gi(9)   = -0.2386191860831969D0
      gaus%gi(10)  =  0.2386191860831971D0
      gaus%gi(11)  =  0.6612093864662646D0
      gaus%gi(12)  =  0.9324695142031521D0
      gaus%ome(7)  =  0.1713244923791704D0
      gaus%ome(8)  =  0.3607615730481386D0
      gaus%ome(9)  =  0.4679139345726910D0
      gaus%ome(10) =  0.4679139345726910D0
      gaus%ome(11) =  0.3607615730481386D0
      gaus%ome(12) =  0.1713244923791703D0
      
C---  N =   8 ( ng1 =  13 , ng2 = 20 )
      gaus%ng1(4)=13
      gaus%ng2(4)=20
      gaus%gi(13)  = -0.9602898564975362D0
      gaus%gi(14)  = -0.7966664774136267D0
      gaus%gi(15)  = -0.5255324099163290D0
      gaus%gi(16)  = -0.1834346424956499D0
      gaus%gi(17)  =  0.1834346424956499D0
      gaus%gi(18)  =  0.5255324099163290D0
      gaus%gi(19)  =  0.7966664774136267D0
      gaus%gi(20)  =  0.9602898564975362D0
      gaus%ome(13) =  0.1012285362903762D0
      gaus%ome(14) =  0.2223810344533745D0
      gaus%ome(15) =  0.3137066458778873D0
      gaus%ome(16) =  0.3626837833783620D0
      gaus%ome(17) =  0.3626837833783620D0
      gaus%ome(18) =  0.3137066458778873D0
      gaus%ome(19) =  0.2223810344533745D0
      gaus%ome(20) =  0.1012285362903764D0
      
C---  N =  16 ( ng1 =  21 , ng2 = 36 )
      gaus%ng1(5)=21
      gaus%ng2(5)=36
      gaus%gi(21) = -0.9894009349916499D0
      gaus%gi(22) = -0.9445750230732326D0
      gaus%gi(23) = -0.8656312023878318D0
      gaus%gi(24) = -0.7554044083550030D0
      gaus%gi(25) = -0.6178762444026437D0
      gaus%gi(26) = -0.4580167776572274D0
      gaus%gi(27) = -0.2816035507792589D0
      gaus%gi(28) = -0.0950125098376375D0
      gaus%gi(29) =  0.0950125098376374D0
      gaus%gi(30) =  0.2816035507792590D0
      gaus%gi(31) =  0.4580167776572274D0
      gaus%gi(32) =  0.6178762444026437D0
      gaus%gi(33) =  0.7554044083550031D0
      gaus%gi(34) =  0.8656312023878318D0
      gaus%gi(35) =  0.9445750230732326D0
      gaus%gi(36) =  0.9894009349916500D0
      gaus%ome(21)=  0.0271524594117541D0
      gaus%ome(22)=  0.0622535239386479D0
      gaus%ome(23)=  0.0951585116824928D0
      gaus%ome(24)=  0.1246289712555339D0
      gaus%ome(25)=  0.1495959888165767D0
      gaus%ome(26)=  0.1691565193950025D0
      gaus%ome(27)=  0.1826034150449236D0
      gaus%ome(28)=  0.1894506104550685D0
      gaus%ome(29)=  0.1894506104550685D0
      gaus%ome(30)=  0.1826034150449236D0
      gaus%ome(31)=  0.1691565193950025D0
      gaus%ome(32)=  0.1495959888165767D0
      gaus%ome(33)=  0.1246289712555338D0
      gaus%ome(34)=  0.0951585116824927D0
      gaus%ome(35)=  0.0622535239386479D0
      gaus%ome(36)=  0.0271524594117539D0
      
C---  N =  32( ng1 =  37 , ng2 = 68 )
      gaus%ng1(6)=37
      gaus%ng2(6)=68
      gaus%gi(37) = -0.9972638618494815D0
      gaus%gi(38) = -0.9856115115452683D0
      gaus%gi(39) = -0.9647622555875064D0
      gaus%gi(40) = -0.9349060759377397D0
      gaus%gi(41) = -0.8963211557660521D0
      gaus%gi(42) = -0.8493676137325700D0
      gaus%gi(43) = -0.7944837959679424D0
      gaus%gi(44) = -0.7321821187402897D0
      gaus%gi(45) = -0.6630442669302152D0
      gaus%gi(46) = -0.5877157572407623D0
      gaus%gi(47) = -0.5068999089322294D0
      gaus%gi(48) = -0.4213512761306354D0
      gaus%gi(49) = -0.3318686022821277D0
      gaus%gi(50) = -0.2392873622521372D0
      gaus%gi(51) = -0.1444719615827965D0
      gaus%gi(52) = -0.0483076656877383D0
      gaus%gi(53) =  0.0483076656877384D0
      gaus%gi(54) =  0.1444719615827965D0
      gaus%gi(55) =  0.2392873622521372D0
      gaus%gi(56) =  0.3318686022821278D0
      gaus%gi(57) =  0.4213512761306355D0
      gaus%gi(58) =  0.5068999089322294D0
      gaus%gi(59) =  0.5877157572407624D0
      gaus%gi(60) =  0.6630442669302153D0
      gaus%gi(61) =  0.7321821187402897D0
      gaus%gi(62) =  0.7944837959679424D0
      gaus%gi(63) =  0.8493676137325699D0
      gaus%gi(64) =  0.8963211557660521D0
      gaus%gi(65) =  0.9349060759377397D0
      gaus%gi(66) =  0.9647622555875065D0
      gaus%gi(67) =  0.9856115115452684D0
      gaus%gi(68) =  0.9972638618494816D0
      gaus%ome(37)=  0.0070186100094701D0
      gaus%ome(38)=  0.0162743947309057D0
      gaus%ome(39)=  0.0253920653092621D0
      gaus%ome(40)=  0.0342738629130214D0
      gaus%ome(41)=  0.0428358980222267D0
      gaus%ome(42)=  0.0509980592623762D0
      gaus%ome(43)=  0.0586840934785356D0
      gaus%ome(44)=  0.0658222227763618D0
      gaus%ome(45)=  0.0723457941088485D0
      gaus%ome(46)=  0.0781938957870703D0
      gaus%ome(47)=  0.0833119242269468D0
      gaus%ome(48)=  0.0876520930044038D0
      gaus%ome(49)=  0.0911738786957639D0
      gaus%ome(50)=  0.0938443990808046D0
      gaus%ome(51)=  0.0956387200792749D0
      gaus%ome(52)=  0.0965400885147278D0
      gaus%ome(53)=  0.0965400885147278D0
      gaus%ome(54)=  0.0956387200792749D0
      gaus%ome(55)=  0.0938443990808046D0
      gaus%ome(56)=  0.0911738786957639D0
      gaus%ome(57)=  0.0876520930044038D0
      gaus%ome(58)=  0.0833119242269468D0
      gaus%ome(59)=  0.0781938957870703D0
      gaus%ome(60)=  0.0723457941088485D0
      gaus%ome(61)=  0.0658222227763619D0
      gaus%ome(62)=  0.0586840934785355D0
      gaus%ome(63)=  0.0509980592623762D0
      gaus%ome(64)=  0.0428358980222267D0
      gaus%ome(65)=  0.0342738629130214D0
      gaus%ome(66)=  0.0253920653092620D0
      gaus%ome(67)=  0.0162743947309056D0
      gaus%ome(68)=  0.0070186100094699D0

C---  N =  40 ( ng1 =  69 , ng2 = 108 )
      gaus%ng1(7)=69
      gaus%ng2(7)=108
      gaus%gi(69)= -0.998237709710559200350D00
      gaus%gi(70)= -0.990726238699457006453D00
      gaus%gi(71)= -0.977259949983774262663D00
      gaus%gi(72)= -0.957916819213791655805D00
      gaus%gi(73)= -0.932812808278676533361D00
      gaus%gi(74)= -0.902098806968874296728D00
      gaus%gi(75)= -0.865959503212259503821D00
      gaus%gi(76)= -0.824612230833311663196D00
      gaus%gi(77)= -0.778305651426519387695D00
      gaus%gi(78)= -0.727318255189927103281D00
      gaus%gi(79)= -0.671956684614179548379D00
      gaus%gi(80)= -0.612553889667980237953D00
      gaus%gi(81)= -0.549467125095128202076D00
      gaus%gi(82)= -0.483075801686178712909D00
      gaus%gi(83)= -0.413779204371605001525D00
      gaus%gi(84)= -0.341994090825758473007D00
      gaus%gi(85)= -0.268152185007253681141D00
      gaus%gi(86)= -0.192697580701371099716D00
      gaus%gi(87)= -0.116084070675255208483D00
      gaus%gi(88)= -0.038772417506050821933D00

      gaus%gi(89)=  0.038772417506050821933D00
      gaus%gi(90)=  0.116084070675255208483D00
      gaus%gi(91)=  0.192697580701371099716D00
      gaus%gi(92)=  0.268152185007253681141D00
      gaus%gi(93)=  0.341994090825758473007D00
      gaus%gi(94)=  0.413779204371605001525D00
      gaus%gi(95)=  0.483075801686178712909D00
      gaus%gi(96)=  0.549467125095128202076D00
      gaus%gi(97)=  0.612553889667980237953D00
      gaus%gi(98)=  0.671956684614179548379D00
      gaus%gi(99)=  0.727318255189927103281D00
      gaus%gi(100)= 0.778305651426519387695D00
      gaus%gi(101)= 0.824612230833311663196D00
      gaus%gi(102)= 0.865959503212259503821D00
      gaus%gi(103)= 0.902098806968874296728D00
      gaus%gi(104)= 0.932812808278676533361D00
      gaus%gi(105)= 0.957916819213791655805D00
      gaus%gi(106)= 0.977259949983774262663D00
      gaus%gi(107)= 0.990726238699457006453D00
      gaus%gi(108)= 0.998237709710559200350D00

      gaus%ome(69)= 0.004521277098533191258D00
      gaus%ome(70)= 0.010498284531152813615D00
      gaus%ome(71)= 0.016421058381907888713D00
      gaus%ome(72)= 0.022245849194166957262D00
      gaus%ome(73)= 0.027937006980023401098D00
      gaus%ome(74)= 0.033460195282547847393D00
      gaus%ome(75)= 0.038782167974472017640D00
      gaus%ome(76)= 0.043870908185673271992D00
      gaus%ome(77)= 0.048695807635072232061D00
      gaus%ome(78)= 0.053227846983936824355D00
      gaus%ome(79)= 0.057439769099391551367D00
      gaus%ome(80)= 0.061306242492928939167D00
      gaus%ome(81)= 0.064804013456601038075D00
      gaus%ome(82)= 0.067912045815233903826D00
      gaus%ome(83)= 0.070611647391286779695D00
      gaus%ome(84)= 0.072886582395804059061D00
      gaus%ome(85)= 0.074723169057968264200D00
      gaus%ome(86)= 0.076110361900626242372D00
      gaus%ome(87)= 0.077039818164247965588D00
      gaus%ome(88)= 0.077505947978424811264D00

      gaus%ome(89)= 0.077505947978424811264D00
      gaus%ome(90)= 0.077039818164247965588D00
      gaus%ome(91)= 0.076110361900626242372D00
      gaus%ome(92)= 0.074723169057968264200D00
      gaus%ome(93)= 0.072886582395804059061D00
      gaus%ome(94)= 0.070611647391286779695D00
      gaus%ome(95)= 0.067912045815233903826D00
      gaus%ome(96)= 0.064804013456601038075D00
      gaus%ome(97)= 0.061306242492928939167D00
      gaus%ome(98)= 0.057439769099391551367D00
      gaus%ome(99)= 0.053227846983936824355D00
      gaus%ome(100)=0.048695807635072232061D00
      gaus%ome(101)=0.043870908185673271992D00
      gaus%ome(102)=0.038782167974472017640D00
      gaus%ome(103)=0.033460195282547847393D00
      gaus%ome(104)=0.027937006980023401098D00
      gaus%ome(105)=0.022245849194166957262D00
      gaus%ome(106)=0.016421058381907888713D00
      gaus%ome(107)=0.010498284531152813615D00
      gaus%ome(108)=0.004521277098533191258D00

C---  N =  48 ( ng1 = 109 , ng2 = 156 )
      gaus%ng1(8)=109
      gaus%ng2(8)=156
      gaus%gi(109)=-0.998771007252426118601D00
      gaus%gi(110)=-0.993530172266350757548D00
      gaus%gi(111)=-0.984124583722826857745D00
      gaus%gi(112)=-0.970591592546247250461D00
      gaus%gi(113)=-0.952987703160430860723D00
      gaus%gi(114)=-0.931386690706554333114D00
      gaus%gi(115)=-0.905879136715569672822D00
      gaus%gi(116)=-0.876572020274247885906D00
      gaus%gi(117)=-0.843588261624393530711D00
      gaus%gi(118)=-0.807066204029442627083D00
      gaus%gi(119)=-0.767159032515740339254D00
      gaus%gi(120)=-0.724034130923814654674D00
      gaus%gi(121)=-0.677872379632663905212D00
      gaus%gi(122)=-0.628867396776513623995D00
      gaus%gi(123)=-0.577224726083972703818D00
      gaus%gi(124)=-0.523160974722233033678D00
      gaus%gi(125)=-0.466902904750958404545D00
      gaus%gi(126)=-0.408686481990716729916D00
      gaus%gi(127)=-0.348755886292160738160D00
      gaus%gi(128)=-0.287362487355455576736D00
      gaus%gi(129)=-0.224763790394689061225D00
      gaus%gi(130)=-0.161222356068891718056D00
      gaus%gi(131)=-0.097004699209462698930D00
      gaus%gi(132)=-0.032380170962869362033D00

      gaus%gi(133)= 0.032380170962869362033D00
      gaus%gi(134)= 0.097004699209462698930D00
      gaus%gi(135)= 0.161222356068891718056D00
      gaus%gi(136)= 0.224763790394689061225D00
      gaus%gi(137)= 0.287362487355455576736D00
      gaus%gi(138)= 0.348755886292160738160D00
      gaus%gi(139)= 0.408686481990716729916D00
      gaus%gi(140)= 0.466902904750958404545D00
      gaus%gi(141)= 0.523160974722233033678D00
      gaus%gi(142)= 0.577224726083972703818D00
      gaus%gi(143)= 0.628867396776513623995D00
      gaus%gi(144)= 0.677872379632663905212D00
      gaus%gi(145)= 0.724034130923814654674D00
      gaus%gi(146)= 0.767159032515740339254D00
      gaus%gi(147)= 0.807066204029442627083D00
      gaus%gi(148)= 0.843588261624393530711D00
      gaus%gi(149)= 0.876572020274247885906D00
      gaus%gi(150)= 0.905879136715569672822D00
      gaus%gi(151)= 0.931386690706554333114D00
      gaus%gi(152)= 0.952987703160430860723D00
      gaus%gi(153)= 0.970591592546247250461D00
      gaus%gi(154)= 0.984124583722826857745D00
      gaus%gi(155)= 0.993530172266350757548D00
      gaus%gi(156)= 0.998771007252426118601D00

      gaus%ome(109)= 0.003153346052305838633D00
      gaus%ome(110)= 0.007327553901276262102D00
      gaus%ome(111)= 0.011477234579234539490D00
      gaus%ome(112)= 0.015579315722943848728D00
      gaus%ome(113)= 0.019616160457355527814D00
      gaus%ome(114)= 0.023570760839324379141D00
      gaus%ome(115)= 0.027426509708356948200D00
      gaus%ome(116)= 0.031167227832798088902D00
      gaus%ome(117)= 0.034777222564770438893D00
      gaus%ome(118)= 0.038241351065830706317D00
      gaus%ome(119)= 0.041545082943464749214D00
      gaus%ome(120)= 0.044674560856694280419D00
      gaus%ome(121)= 0.047616658492490474826D00
      gaus%ome(122)= 0.050359035553854474958D00
      gaus%ome(123)= 0.052890189485193667096D00
      gaus%ome(124)= 0.055199503699984162868D00
      gaus%ome(125)= 0.057277292100403215705D00
      gaus%ome(126)= 0.059114839698395635746D00
      gaus%ome(127)= 0.060704439165893880053D00
      gaus%ome(128)= 0.062039423159892663904D00
      gaus%ome(129)= 0.063114192286254025657D00
      gaus%ome(130)= 0.063924238584648186624D00
      gaus%ome(131)= 0.064466164435950082207D00
      gaus%ome(132)= 0.064737696812683922503D00

      gaus%ome(133)= 0.064737696812683922503D00
      gaus%ome(134)= 0.064466164435950082207D00
      gaus%ome(135)= 0.063924238584648186624D00
      gaus%ome(136)= 0.063114192286254025657D00
      gaus%ome(137)= 0.062039423159892663904D00
      gaus%ome(138)= 0.060704439165893880053D00
      gaus%ome(139)= 0.059114839698395635746D00
      gaus%ome(140)= 0.057277292100403215705D00
      gaus%ome(141)= 0.055199503699984162868D00
      gaus%ome(142)= 0.052890189485193667096D00
      gaus%ome(143)= 0.050359035553854474958D00
      gaus%ome(144)= 0.047616658492490474826D00
      gaus%ome(145)= 0.044674560856694280419D00
      gaus%ome(146)= 0.041545082943464749214D00
      gaus%ome(147)= 0.038241351065830706317D00
      gaus%ome(148)= 0.034777222564770438893D00
      gaus%ome(149)= 0.031167227832798088902D00
      gaus%ome(150)= 0.027426509708356948200D00
      gaus%ome(151)= 0.023570760839324379141D00
      gaus%ome(152)= 0.019616160457355527814D00
      gaus%ome(153)= 0.015579315722943848728D00
      gaus%ome(154)= 0.011477234579234539490D00
      gaus%ome(155)= 0.007327553901276262102D00
      gaus%ome(156)= 0.003153346052305838633D00
c
C---  N =  64 ( ng1 = 157 , ng2 = 220 )
      gaus%ng1(9)=157
      gaus%ng2(9)=220

      gaus%gi(157)=-0.9993050417357721394569D0
      gaus%gi(158)=-0.9963401167719552793469D0
      gaus%gi(159)=-0.9910133714767443207394D0
      gaus%gi(160)=-0.9833362538846259569313D0
      gaus%gi(161)=-0.9733268277899109637419D0
      gaus%gi(162)=-0.9610087996520537189186D0
      gaus%gi(163)=-0.9464113748584028160625D0
      gaus%gi(164)=-0.9295691721319395758215D0
      gaus%gi(165)=-0.9105221370785028057564D0
      gaus%gi(166)=-0.8893154459951141058534D0
      gaus%gi(167)=-0.8659993981540928197608D0
      gaus%gi(168)=-0.840629296252580362752D0
      gaus%gi(169)=-0.813265315122797559742D0
      gaus%gi(170)=-0.7839723589433414076102D0
      gaus%gi(171)=-0.752819907260531896612D0
      gaus%gi(172)=-0.7198818501716108268489D0
      gaus%gi(173)=-0.6852363130542332425636D0
      gaus%gi(174)=-0.648965471254657339858D0
      gaus%gi(175)=-0.6111553551723932502489D0
      gaus%gi(176)=-0.5718956462026340342839D0
      gaus%gi(177)=-0.531279464019894545658D0
      gaus%gi(178)=-0.4894031457070529574785D0
      gaus%gi(179)=-0.446366017253464087985D0
      gaus%gi(180)=-0.4022701579639916036958D0
      gaus%gi(181)=-0.3572201583376681159504D0
      gaus%gi(182)=-0.3113228719902109561575D0
      gaus%gi(183)=-0.264687162208767416374D0
      gaus%gi(184)=-0.2174236437400070841497D0
      gaus%gi(185)=-0.1696444204239928180373D0
      gaus%gi(186)=-0.1214628192961205544704D0
      gaus%gi(187)=-0.0729931217877990394495D0
      gaus%gi(188)=-0.02435029266342443250896D0
      gaus%gi(189)=0.024350292663424432509D0
      gaus%gi(190)=0.0729931217877990394495D0
      gaus%gi(191)=0.1214628192961205544704D0
      gaus%gi(192)=0.1696444204239928180373D0
      gaus%gi(193)=0.2174236437400070841497D0
      gaus%gi(194)=0.264687162208767416374D0
      gaus%gi(195)=0.311322871990210956158D0
      gaus%gi(196)=0.3572201583376681159504D0
      gaus%gi(197)=0.4022701579639916036958D0
      gaus%gi(198)=0.446366017253464087985D0
      gaus%gi(199)=0.489403145707052957479D0
      gaus%gi(200)=0.531279464019894545658D0
      gaus%gi(201)=0.5718956462026340342839D0
      gaus%gi(202)=0.611155355172393250249D0
      gaus%gi(203)=0.6489654712546573398578D0
      gaus%gi(204)=0.6852363130542332425636D0
      gaus%gi(205)=0.7198818501716108268489D0
      gaus%gi(206)=0.7528199072605318966119D0
      gaus%gi(207)=0.7839723589433414076102D0
      gaus%gi(208)=0.8132653151227975597419D0
      gaus%gi(209)=0.8406292962525803627517D0
      gaus%gi(210)=0.8659993981540928197608D0
      gaus%gi(211)=0.8893154459951141058534D0
      gaus%gi(212)=0.9105221370785028057564D0
      gaus%gi(213)=0.9295691721319395758215D0
      gaus%gi(214)=0.9464113748584028160625D0
      gaus%gi(215)=0.9610087996520537189186D0
      gaus%gi(216)=0.9733268277899109637419D0
      gaus%gi(217)=0.9833362538846259569313D0
      gaus%gi(218)=0.9910133714767443207394D0
      gaus%gi(219)=0.9963401167719552793469D0
      gaus%gi(220)=0.9993050417357721394569D0

      gaus%ome(157)= 0.0017832807216964329473D0
      gaus%ome(158)= 0.0041470332605624676353D0
      gaus%ome(159)= 0.006504457968978362856D0
      gaus%ome(160)= 0.008846759826363947723D0
      gaus%ome(161)= 0.0111681394601311288186D0
      gaus%ome(162)= 0.0134630478967186425981D0
      gaus%ome(163)= 0.015726030476024719322D0
      gaus%ome(164)= 0.01795171577569734308505D0
      gaus%ome(165)= 0.020134823153530209372D0
      gaus%ome(166)= 0.0222701738083832541593D0
      gaus%ome(167)= 0.0243527025687108733382D0
      gaus%ome(168)= 0.026377469715054658672D0
      gaus%ome(169)= 0.0283396726142594832275D0
      gaus%ome(170)= 0.03023465707240247886797D0
      gaus%ome(171)= 0.0320579283548515535855D0
      gaus%ome(172)= 0.03380516183714160939157D0
      gaus%ome(173)= 0.0354722132568823838107D0
      gaus%ome(174)= 0.0370551285402400460404D0
      gaus%ome(175)= 0.038550153178615629129D0
      gaus%ome(176)= 0.0399537411327203413867D0
      gaus%ome(177)= 0.04126256324262352861D0
      gaus%ome(178)= 0.04247351512365358900734D0
      gaus%ome(179)= 0.043583724529323453377D0
      gaus%ome(180)= 0.04459055816375656306D0
      gaus%ome(181)= 0.0454916279274181444798D0
      gaus%ome(182)= 0.046284796581314417296D0
      gaus%ome(183)= 0.04696818281621001732533D0
      gaus%ome(184)= 0.0475401657148303086623D0
      gaus%ome(185)= 0.047999388596458307728D0
      gaus%ome(186)= 0.0483447622348029571698D0
      gaus%ome(187)= 0.0485754674415034269348D0
      gaus%ome(188)= 0.0486909570091397203834D0
      gaus%ome(189)= 0.0486909570091397203834D0
      gaus%ome(190)= 0.048575467441503426935D0
      gaus%ome(191)= 0.04834476223480295717D0
      gaus%ome(192)= 0.0479993885964583077281D0
      gaus%ome(193)= 0.0475401657148303086623D0
      gaus%ome(194)= 0.046968182816210017325D0
      gaus%ome(195)= 0.046284796581314417296D0
      gaus%ome(196)= 0.04549162792741814448D0
      gaus%ome(197)= 0.0445905581637565630601D0
      gaus%ome(198)= 0.043583724529323453377D0
      gaus%ome(199)= 0.042473515123653589007D0
      gaus%ome(200)= 0.0412625632426235286102D0
      gaus%ome(201)= 0.0399537411327203413867D0
      gaus%ome(202)= 0.038550153178615629129D0
      gaus%ome(203)= 0.0370551285402400460404D0
      gaus%ome(204)= 0.03547221325688238381069D0
      gaus%ome(205)= 0.033805161837141609392D0
      gaus%ome(206)= 0.032057928354851553585D0
      gaus%ome(207)= 0.030234657072402478868D0
      gaus%ome(208)= 0.0283396726142594832275D0
      gaus%ome(209)= 0.026377469715054658672D0
      gaus%ome(210)= 0.024352702568710873338D0
      gaus%ome(211)= 0.0222701738083832541593D0
      gaus%ome(212)= 0.0201348231535302093723D0
      gaus%ome(213)= 0.017951715775697343085D0
      gaus%ome(214)= 0.015726030476024719322D0
      gaus%ome(215)= 0.0134630478967186425981D0
      gaus%ome(216)= 0.0111681394601311288186D0
      gaus%ome(217)= 0.008846759826363947723D0
      gaus%ome(218)= 0.006504457968978362856D0
      gaus%ome(219)= 0.0041470332605624676353D0
      gaus%ome(220)= 0.0017832807216964329473D0



      END
      
C______________________________________________________________________C
C______________________________________________________________________C
      SUBROUTINE WrSMat(mat,nrow,ncol,io)
C        __    _      ___
C        Write Single Matrix
C______________________________________________________________________C
C______________________________________________________________________C
      INTEGER nrow,ncol,io,j
      REAL*8 mat(nrow,ncol)

      DO j=1,ncol
        CALL wrvec(io,nrow,mat(1,j))
      END DO

      END
C______________________________________________________________________C
C______________________________________________________________________C
      SUBROUTINE RdSMat(mat,nrow,ncol,io)
C        _  _ _      ___
C        Read Single Matrix
C______________________________________________________________________C
C______________________________________________________________________C
      INTEGER nrow,ncol,io,j
      REAL*8 mat(nrow,ncol)

      DO j=1,ncol
        CALL rdvec(io,nrow,mat(1,j))
      END DO

      END

C----------------------------------------------------------------------C
c----------------------------------------------------------------------c
c                                                                      c
      SUBROUTINE rdvec(ifr,nnx,vec)
c                                                                      c
c----------------------------------------------------------------------c
C----------------------------------------------------------------------C
c.......................................................................
c..                                                                   ..
c..   REad VECtor using MAXSIZE chunks                                ..
c..   --   ---                                                        ..
c.......................................................................
      INTEGER ifr,nnx,maxsize,nblock,i,j,k
      REAL*8  vec(nnx)
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
      SUBROUTINE wrvec(ifr,nnx,vec)
c                                                                      c
c----------------------------------------------------------------------c
C----------------------------------------------------------------------C
c.......................................................................
c..                                                                   ..
c..   WRite VECtor using MAXSIZE chunks                               ..
c..   --    ---                                                       ..
c.......................................................................
      INTEGER ifr,nnx,maxsize,nblock,i,j,k
      REAL*8  vec(nnx)
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
      


C----------------------------------------------------------------------C
c
      SUBROUTINE DomainIntegral(mesh,gp,f,res)
C
C     $: Integrira f po celotnem volumnu
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(meshType) :: mesh
      TYPE(GaussType) :: gp
      INTEGER ic
      REAL(8) x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
      REAL(8) x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8
      REAL(8) f(mesh%nnodes),res,rescell,vol,resvol

      res=0.0D0


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

        CALL IntCellVolumeFunc(mesh,gp,
     &    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &    x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &    f,ic,rescell)

        res=res+rescell

      END DO

      END SUBROUTINE
c
C----------------------------------------------------------------------C
c
      SUBROUTINE IntCellVolumeFunc(mesh,gp,
     &    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &    x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &    f,ic,res)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  INTEgration over Domain cell                                    **
c **  ----             -                                              **
c **  internal cell geometry : LINEAR ( 8-node )                      **
c **  internal cell function interpolation : QUADRATIC ( 27-node )    **
c **  integrira determinanto jacobijeve matrike = izracuna volumen    **
c **                                                                  **
c **********************************************************************
      USE inc_types
      TYPE(gausstype) :: gp
      TYPE(meshType) :: mesh
      REAL(8) f(mesh%nnodes)
      INTEGER ic
c
      INTEGER i,j,k,p,ng1,ng2
c
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &        x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &        eta1m,eta1p,eta2m,eta2p,eta3m,eta3p,
     &        aj11,aj12,aj13,aj21,aj22,aj23,aj31,aj32,aj33,
     &        ajac

      REAL(8) res,vol
      REAL(8) fig27(27)
C
C*** SET NUMBER OF INTEGRATION POINTS
      ng1=gp%ng1(2)  ! 2 pomeni pomoje dovolj gaussovnih tock
      ng2=gp%ng2(2)

      res=0.0D0
      vol=0.0D0
C
C***  27 NODE CONTINUOUS DOMAIN CELL (8 NODE GEOMETRY)
C
      DO K=ng1,ng2
        DO J=ng1,ng2
          DO I=ng1,ng2
            ETA1M=1.0D0-gp%GI(I)
            ETA1P=1.0D0+gp%GI(I)
            ETA2M=1.0D0-gp%GI(J)
            ETA2P=1.0D0+gp%GI(J)
            ETA3M=1.0D0-gp%GI(K)
            ETA3P=1.0D0+gp%GI(K)
C
C***  Jacobian derivatives
C
            AJ11=-ETA2M*ETA3M*X1+ETA2M*ETA3M*X2-ETA2P*ETA3M*X3+ETA2P*ETA3M*X4
     &           -ETA2M*ETA3P*X5+ETA2M*ETA3P*X6-ETA2P*ETA3P*X7+ETA2P*ETA3P*X8
            AJ12=-ETA2M*ETA3M*Y1+ETA2M*ETA3M*Y2-ETA2P*ETA3M*Y3+ETA2P*ETA3M*Y4
     &           -ETA2M*ETA3P*Y5+ETA2M*ETA3P*Y6-ETA2P*ETA3P*Y7+ETA2P*ETA3P*Y8
            AJ13=-ETA2M*ETA3M*Z1+ETA2M*ETA3M*Z2-ETA2P*ETA3M*Z3+ETA2P*ETA3M*Z4
     &           -ETA2M*ETA3P*Z5+ETA2M*ETA3P*Z6-ETA2P*ETA3P*Z7+ETA2P*ETA3P*Z8
            AJ21=-ETA1M*ETA3M*X1-ETA1P*ETA3M*X2+ETA1M*ETA3M*X3+ETA1P*ETA3M*X4
     &           -ETA1M*ETA3P*X5-ETA1P*ETA3P*X6+ETA1M*ETA3P*X7+ETA1P*ETA3P*X8
            AJ22=-ETA1M*ETA3M*Y1-ETA1P*ETA3M*Y2+ETA1M*ETA3M*Y3+ETA1P*ETA3M*Y4
     &           -ETA1M*ETA3P*Y5-ETA1P*ETA3P*Y6+ETA1M*ETA3P*Y7+ETA1P*ETA3P*Y8
            AJ23=-ETA1M*ETA3M*Z1-ETA1P*ETA3M*Z2+ETA1M*ETA3M*Z3+ETA1P*ETA3M*Z4
     &           -ETA1M*ETA3P*Z5-ETA1P*ETA3P*Z6+ETA1M*ETA3P*Z7+ETA1P*ETA3P*Z8
            AJ31=-ETA1M*ETA2M*X1-ETA1P*ETA2M*X2-ETA1M*ETA2P*X3-ETA1P*ETA2P*X4
     &           +ETA1M*ETA2M*X5+ETA1P*ETA2M*X6+ETA1M*ETA2P*X7+ETA1P*ETA2P*X8
            AJ32=-ETA1M*ETA2M*Y1-ETA1P*ETA2M*Y2-ETA1M*ETA2P*Y3-ETA1P*ETA2P*Y4
     &           +ETA1M*ETA2M*Y5+ETA1P*ETA2M*Y6+ETA1M*ETA2P*Y7+ETA1P*ETA2P*Y8
            AJ33=-ETA1M*ETA2M*Z1-ETA1P*ETA2M*Z2-ETA1M*ETA2P*Z3-ETA1P*ETA2P*Z4
     &           +ETA1M*ETA2M*Z5+ETA1P*ETA2M*Z6+ETA1M*ETA2P*Z7+ETA1P*ETA2P*Z8

            AJAC=(Aj11*Aj22*Aj33+Aj12*Aj23*Aj31+Aj13*Aj21*Aj32-
     &            Aj11*Aj23*Aj32-Aj12*Aj21*Aj33-Aj13*Aj22*Aj31)
            AJAC=1.953125D-3*AJAC*gp%OME(I)*gp%OME(J)*gp%OME(K)
c                  1/512=1.953125D-3

            CALL cshape27(gp%GI(I),gp%GI(J),gp%GI(K),fig27,27)
c
            DO P=1,27
              res=res+AJAC*FIG27(P)*f(mesh%idc(ic,P))
            END DO
c
          END DO !I
        END DO !J
      END DO !K
c
      END

c




C----------------------------------------------------------------------C
c
      SUBROUTINE TestIfifj(mesh)
C
C     $: Testiram integracijo z intepolacijo dveh funkcij pod integralom
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(meshType) :: mesh

      INTEGER ic
      REAL(8) total,vsota

      REAL(8), ALLOCATABLE :: temp(:),Ififj(:,:),uz(:)

      INTEGER inode, jnode, i, j

C     Set up integrals of interpolation functions fi*fj
      ALLOCATE (Ififj(27,27))
      CALL setIfiji27(Ififj)

C     Set up fake temperature data, for test only
      ALLOCATE(temp(mesh%nnodes))
      DO i=1,mesh%nnodes
        temp(i) = 1.0D0 + mesh%x(i,1)**2  + mesh%x(i,2)**2  + mesh%x(i,3)**2   ! T = 1 + x^2 + y^2 + z^2
      END DO

C     tu bo prisla zanka po izvornih tockah

c     za izbrano izvorno tocko zgeneriram uzvezdica
      ALLOCATE(uz(mesh%nnodes))
      DO i=1,mesh%nnodes
        uz(i) = 1.0D0 + mesh%x(i,1)  + mesh%x(i,2)  + mesh%x(i,3)   ! uz = 1 + x + y + z
      END DO

c     rezultat integrala int (temp*uz) dOmega po celi mrezi
      total=0.0D0

      DO ic=1,mesh%nicell ! po vseh celicah v mrezi
        DO i=1,mesh%npoc   ! po interpolacijskih funkcijah Temperature
          inode = mesh%idc(ic,i)
          T = temp(inode) ! temperatura, s katero se integral mnozi

          vsota = 0.0D0
          DO j=1,mesh%npoc   ! po interpolacijskih funkcijah uz
            jnode = mesh%idc(ic,j)
c           mnozim matrika Ififj krat vektor uz iz vozlisc v celici
            vsota = vsota + Ififj(i,j)*uz(jnode)
          END DO
c         skalarni produkt vektor T * (Ififj.uz) * Jacobi ( 8 ker je volumen (-1,1)^3 = 8
          total = total + vsota * T * mesh%CellVolume(ic) / 8.0D0
        END DO
      END DO

      print *,"tot ",total ! analiticno = 21/4=5.25

      DEALLOCATE (Ififj,temp,uz)
      STOP

      END SUBROUTINE
c

