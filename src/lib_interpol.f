C----------------------------------------------------------------------C
      SUBROUTINE SetLC(lc1,lc2,lc3)
      
      REAL(8) lc1(27),lc2(27),lc3(27)
      INTEGER i


      lc1(1)=0.0D+0
      lc1(2)=0.0D+0
      lc1(3)=1.0D+0
      lc1(4)=1.0D+0
      lc1(5)=0.0D+0
      lc1(6)=0.0D+0
      lc1(7)=1.0D+0
      lc1(8)=1.0D+0
      lc1(9)=0.0D+0
      lc1(10)=5.0D-1
      lc1(11)=1.0D+0
      lc1(12)=5.0D-1
      lc1(13)=0.0D+0
      lc1(14)=0.0D+0
      lc1(15)=1.0D+0
      lc1(16)=1.0D+0
      lc1(17)=0.0D+0
      lc1(18)=5.0D-1
      lc1(19)=1.0D+0
      lc1(20)=5.0D-1
      lc1(21)=5.0D-1
      lc1(22)=0.0D+0
      lc1(23)=5.0D-1
      lc1(24)=1.0D+0
      lc1(25)=5.0D-1
      lc1(26)=5.0D-1
      lc1(27)=5.0D-1
     
      lc2(1)=1.0D+0
      lc2(2)=1.0D+0
      lc2(3)=1.0D+0
      lc2(4)=1.0D+0
      lc2(5)=0.0D+0
      lc2(6)=0.0D+0
      lc2(7)=0.0D+0
      lc2(8)=0.0D+0
      lc2(9)=1.0D+0
      lc2(10)=1.0D+0
      lc2(11)=1.0D+0
      lc2(12)=1.0D+0
      lc2(13)=5.0D-1
      lc2(14)=5.0D-1
      lc2(15)=5.0D-1
      lc2(16)=5.0D-1
      lc2(17)=0.0D+0
      lc2(18)=0.0D+0
      lc2(19)=0.0D+0
      lc2(20)=0.0D+0
      lc2(21)=1.0D+0
      lc2(22)=5.0D-1
      lc2(23)=5.0D-1
      lc2(24)=5.0D-1
      lc2(25)=5.0D-1
      lc2(26)=0.0D+0
      lc2(27)=5.0D-1
     
      lc3(1)=1.0D+0
      lc3(2)=0.0D+0
      lc3(3)=0.0D+0
      lc3(4)=1.0D+0
      lc3(5)=1.0D+0
      lc3(6)=0.0D+0
      lc3(7)=0.0D+0
      lc3(8)=1.0D+0
      lc3(9)=5.0D-1
      lc3(10)=0.0D+0
      lc3(11)=5.0D-1
      lc3(12)=1.0D+0
      lc3(13)=1.0D+0
      lc3(14)=0.0D+0
      lc3(15)=0.0D+0
      lc3(16)=1.0D+0
      lc3(17)=5.0D-1
      lc3(18)=0.0D+0
      lc3(19)=5.0D-1
      lc3(20)=1.0D+0
      lc3(21)=5.0D-1
      lc3(22)=5.0D-1
      lc3(23)=0.0D+0
      lc3(24)=5.0D-1
      lc3(25)=1.0D+0
      lc3(26)=5.0D-1
      lc3(27)=5.0D-1  


c     ker se mi ni dalo na novo natipkat     
      DO i=1,27
        lc1(i)=lc1(i)*2.0D0-1.0D0            
        lc2(i)=lc2(i)*2.0D0-1.0D0      
        lc3(i)=lc3(i)*2.0D0-1.0D0              
      END DO
 
      END

C----------------------------------------------------------------------C
      SUBROUTINE SetLC2Dnzv(lc1,lc2)
      REAL(8) lc1(4),lc2(4),tc

      tc=3.0D0/4.0D0

      lc1(1)=-tc
      lc2(1)=-tc

      lc1(2)=+tc
      lc2(2)=-tc

      lc1(3)=+tc
      lc2(3)=+tc

      lc1(4)=-tc
      lc2(4)=+tc

      END

C----------------------------------------------------------------------C

      SUBROUTINE setGRAD(mesh,u,gradu)
C
C     Calculate field function derivaties
C      
      USE inc_types   

      TYPE(meshType) :: mesh
      INTEGER l,ic,it,dn,i

      REAL(8) dfidxi(mesh%npoc),dfideta(mesh%npoc),dfidzeta(mesh%npoc)
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &        x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &        dxidx,detadx,dzetadx,
     &        dxidy,detady,dzetady,
     &        dxidz,detadz,dzetadz,ajac,
     &        dfildx,dfildy,dfildz    
     
      REAL(8) lc1(27),lc2(27),lc3(27)
      REAL(8) u(mesh%nnodes),gradu(mesh%nnodes,3)
      REAL(8), ALLOCATABLE :: SumVol(:)
      
      CALL SetLC(lc1,lc2,lc3)

c     init
      gradu=0.0D0
      ALLOCATE (SumVol(mesh%nnodes)) ! vsota volumnom elementov katerim pripada to vozlisce
      sumVol=0.0D0
            
      DO ic=1,mesh%nicell  ! po vseh celicah v mrezi
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
        DO it=1,mesh%npoc ! po vseh vozliscih znotraj ene celice
          dn=mesh%idc(ic,it)
          SumVol(dn)=SumVol(dn)+mesh%cellVolume(ic)  ! izracunam volumen vseh elemntov, katerim pripada dn
c         lokalne koordinate it-tega vozlica so v lc1,lc2,lc3
c         rabim odvode vseh interpolacijskih funkcij v tej tocki
          CALL cshaped27(lc1(it),lc2(it),lc3(it),dfidxi,dfideta,dfidzeta,mesh%npoc)
c         rabim d(xi)/dx, itd
          CALL detJac3Dlcd(lc1(it),lc2(it),lc3(it),
     &        x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &        x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &        dxidx,detadx,dzetadx,
     &        dxidy,detady,dzetady,
     &        dxidz,detadz,dzetadz,ajac)

          DO l=1,mesh%npoc ! vsota za izracun odvoda v dn-tem vozliscu
            dfildx=dfidxi(l)*dxidx+dfideta(l)*detadx+dfidzeta(l)*dzetadx
            dfildy=dfidxi(l)*dxidy+dfideta(l)*detady+dfidzeta(l)*dzetady
            dfildz=dfidxi(l)*dxidz+dfideta(l)*detadz+dfidzeta(l)*dzetadz

            gradu(dn,1)=gradu(dn,1)+u(mesh%idc(ic,l))*dfildx*mesh%CellVolume(ic)
            gradu(dn,2)=gradu(dn,2)+u(mesh%idc(ic,l))*dfildy*mesh%CellVolume(ic)
            gradu(dn,3)=gradu(dn,3)+u(mesh%idc(ic,l))*dfildz*mesh%CellVolume(ic)
          END DO  
        END DO         
      END DO

c     Dokoncno izracunam utezeno povprecje
      DO i=1,mesh%nnodes
        gradu(i,1)=gradu(i,1)/sumVol(i)
        gradu(i,2)=gradu(i,2)/sumVol(i)
        gradu(i,3)=gradu(i,3)/sumVol(i)
      END DO

      DEALLOCATE (SumVol)

      END   



C----------------------------------------------------------------------C

      SUBROUTINE detJac3Dlcd(xi,eta,zeta,
     &        x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &        x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &        dxidx,detadx,dzetadx,
     &        dxidy,detady,dzetady,
     &        dxidz,detadz,dzetadz,ajac)     

C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Calculate Jacobian Determinant and derivatives of local coors, 3D  **
c **  -    -----                                 -                    **
c **********************************************************************
      REAL(8) eta1m,eta1p,eta2m,eta2p,eta3m,eta3p,
     &        xi,eta,zeta,
     &        x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &        x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &        aj11,aj12,aj13,aj21,aj22,aj23,aj31,aj32,aj33,
     &        ajac,
     &        dxidx,detadx,dzetadx,
     &        dxidy,detady,dzetady,
     &        dxidz,detadz,dzetadz       
C
C
C
      ETA1M=1.0D0-xi   
      ETA1P=1.0D0+xi
      ETA2M=1.0D0-eta
      ETA2P=1.0D0+eta
      ETA3M=1.0D0-zeta
      ETA3P=1.0D0+zeta
C
C***  Jacobian derivatives  AJ11 = d(x)/d(xi), Aj33= d(z)/d(zeta)
C
      AJ11=0.125D0*(-ETA2M*ETA3M*X1+ETA2M*ETA3M*X2-ETA2P*ETA3M*X3+ETA2P*ETA3M*X4
     &              -ETA2M*ETA3P*X5+ETA2M*ETA3P*X6-ETA2P*ETA3P*X7+ETA2P*ETA3P*X8)

      AJ12=0.125D0*(-ETA2M*ETA3M*Y1+ETA2M*ETA3M*Y2-ETA2P*ETA3M*Y3+ETA2P*ETA3M*Y4
     &              -ETA2M*ETA3P*Y5+ETA2M*ETA3P*Y6-ETA2P*ETA3P*Y7+ETA2P*ETA3P*Y8)

      AJ13=0.125D0*(-ETA2M*ETA3M*Z1+ETA2M*ETA3M*Z2-ETA2P*ETA3M*Z3+ETA2P*ETA3M*Z4
     &              -ETA2M*ETA3P*Z5+ETA2M*ETA3P*Z6-ETA2P*ETA3P*Z7+ETA2P*ETA3P*Z8)
c     
      AJ21=0.125D0*(-ETA1M*ETA3M*X1-ETA1P*ETA3M*X2+ETA1M*ETA3M*X3+ETA1P*ETA3M*X4
     &              -ETA1M*ETA3P*X5-ETA1P*ETA3P*X6+ETA1M*ETA3P*X7+ETA1P*ETA3P*X8)

      AJ22=0.125D0*(-ETA1M*ETA3M*Y1-ETA1P*ETA3M*Y2+ETA1M*ETA3M*Y3+ETA1P*ETA3M*Y4
     &              -ETA1M*ETA3P*Y5-ETA1P*ETA3P*Y6+ETA1M*ETA3P*Y7+ETA1P*ETA3P*Y8)

      AJ23=0.125D0*(-ETA1M*ETA3M*Z1-ETA1P*ETA3M*Z2+ETA1M*ETA3M*Z3+ETA1P*ETA3M*Z4
     &              -ETA1M*ETA3P*Z5-ETA1P*ETA3P*Z6+ETA1M*ETA3P*Z7+ETA1P*ETA3P*Z8)
c
      AJ31=0.125D0*(-ETA1M*ETA2M*X1-ETA1P*ETA2M*X2-ETA1M*ETA2P*X3-ETA1P*ETA2P*X4
     &              +ETA1M*ETA2M*X5+ETA1P*ETA2M*X6+ETA1M*ETA2P*X7+ETA1P*ETA2P*X8)

      AJ32=0.125D0*(-ETA1M*ETA2M*Y1-ETA1P*ETA2M*Y2-ETA1M*ETA2P*Y3-ETA1P*ETA2P*Y4
     &              +ETA1M*ETA2M*Y5+ETA1P*ETA2M*Y6+ETA1M*ETA2P*Y7+ETA1P*ETA2P*Y8)

      AJ33=0.125D0*(-ETA1M*ETA2M*Z1-ETA1P*ETA2M*Z2-ETA1M*ETA2P*Z3-ETA1P*ETA2P*Z4
     &              +ETA1M*ETA2M*Z5+ETA1P*ETA2M*Z6+ETA1M*ETA2P*Z7+ETA1P*ETA2P*Z8)
c
c           Determinanata Jacobija
c
      AJAC=(Aj11*Aj22*Aj33+Aj12*Aj23*Aj31+Aj13*Aj21*Aj32-
     &      Aj11*Aj23*Aj32-Aj12*Aj21*Aj33-Aj13*Aj22*Aj31)      
c
c           Racunam d(xi)/d(x)
c           dobim z obratom jacobijeve matrike
c
c     prvi stolpec
        dxidx  =(aj33*aj22-aj23*aj32)/ajac
        dxidy =(aj23*aj31-aj21*aj33)/ajac
        dxidz=(aj21*aj32-aj22*aj31)/ajac
c     drugi stolpec           
        detadx  =(aj13*aj32-aj33*aj12)/ajac
        detady =(aj33*aj11-aj13*aj31)/ajac
        detadz=(aj12*aj31-aj11*aj32)/ajac
c     tretji stolpec           
        dzetadx  =(aj23*aj12-aj13*aj22)/ajac
        dzetady =(aj21*aj13-aj23*aj11)/ajac
        dzetadz=(aj11*aj22-aj21*aj12)/ajac

c     preverim, da je produt matrik enak identiteti
c      Print *,aj11*dxidx+aj12*detadx+aj13*dzetadx
c      Print *,aj11*dxidy+aj12*detady+aj13*dzetady
c      Print *,aj11*dxidz+aj12*detadz+aj13*dzetadz

c      Print *,aj21*dxidx+aj22*detadx+aj23*dzetadx
c      Print *,aj21*dxidy+aj22*detady+aj23*dzetady
c      Print *,aj21*dxidz+aj22*detadz+aj23*dzetadz

c      Print *,aj31*dxidx+aj32*detadx+aj33*dzetadx
c      Print *,aj31*dxidy+aj32*detady+aj33*dzetady
c      Print *,aj31*dxidz+aj32*detadz+aj33*dzetadz
c      stop


      END

C----------------------------------------------------------------------C

      SUBROUTINE cshape8(xi,eta,zeta,fi,mpom)

C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **   Cell SHAPE functions for continuous linear 8-node cell         **
c **   -    -----                                 -                   **
c **********************************************************************
      INTEGER mpom
      REAL(8) xi,eta,zeta,fi(mpom)
      REAL(8) eta1m,eta1p,eta2m,eta2p,eta3m,eta3p

      ETA1M=1.0D0-xi
      ETA1P=1.0D0+xi
      ETA2M=1.0D0-eta
      ETA2P=1.0D0+eta
      ETA3M=1.0D0-zeta
      ETA3P=1.0D0+zeta

      FI(1)=0.125D0*ETA1M*ETA2M*ETA3M
      FI(2)=0.125D0*ETA1P*ETA2M*ETA3M
      FI(3)=0.125D0*ETA1M*ETA2P*ETA3M
      FI(4)=0.125D0*ETA1P*ETA2P*ETA3M
      FI(5)=0.125D0*ETA1M*ETA2M*ETA3P
      FI(6)=0.125D0*ETA1P*ETA2M*ETA3P
      FI(7)=0.125D0*ETA1M*ETA2P*ETA3P
      FI(8)=0.125D0*ETA1P*ETA2P*ETA3P

      END


C----------------------------------------------------------------------C

      SUBROUTINE cshape27(xi,eta,zeta,fi,mpom)

C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Cell SHAPE functions for continuous quadr. 27-node cell         **
c **  -    -----                                 -                    **
c **********************************************************************
      INTEGER mpom
      REAL(8) xi,eta,zeta,fi
      DIMENSION fi(mpom)
      
      REAL(8) ena8,ena4,ena2,epz,emz,epe,eme,emx,epx
      
      ena8=1.0D0/8.0D0
      ena4=1.0D0/4.0D0
      ena2=1.0D0/2.0D0
      epz=1.0D0+zeta
      emz=1.0D0-zeta
      epe=1.0D0+eta
      eme=1.0D0-eta
      epx=1.0D0+xi
      emx=1.0D0-xi   


      fi( 1)=-ena8*zeta*epz*eta*epe*emx*xi  
      fi( 2)= ena8* emz* zeta*  eta*  epe* emx* xi  
      fi( 3)=-ena8* emz* zeta*  eta*  epe* xi*  epx 
      fi( 4)= ena8* zeta*  epz* eta*  epe* xi*  epx 
      fi( 5)= ena8* zeta*  epz* eme* eta*  emx* xi  
      fi( 6)=-ena8* emz* zeta*  eme* eta*  emx* xi  
      fi( 7)= ena8* emz* zeta*  eme* eta*  xi*  epx 
      fi( 8)=-ena8* zeta*  epz* eme* eta*  xi*  epx 
      fi( 9)=-ena4* emz* epz* eta*  epe* emx* xi  
      fi(10)=-ena4* emz* zeta*  eta*  epe* emx* epx 
      fi(11)= ena4* emz* epz* eta*  epe* xi*  epx 
      fi(12)= ena4* zeta*  epz* eta*  epe* emx* epx 
      fi(13)=-ena4* zeta*  epz* eme* epe* emx* xi  
      fi(14)= ena4* emz* zeta*  eme* epe* emx* xi  
      fi(15)=-ena4* emz* zeta*  eme* epe* xi * epx 
      fi(16)= ena4* zeta*  epz* eme* epe* xi * epx 
      fi(17)= ena4* emz* epz* eme* eta*  emx* xi  
      fi(18)= ena4* emz* zeta*  eme* eta*  emx* epx 
      fi(19)=-ena4* emz* epz* eme* eta*  xi * epx 
      fi(20)=-ena4* zeta*  epz* eme* eta*  emx* epx 
      fi(21)= ena2* emz* epz* eta*  epe* emx* epx 
      fi(22)=-ena2* emz* epz* eme* epe* emx* xi  
      fi(23)=-ena2* emz* zeta*  eme* epe* emx* epx 
      fi(24)= ena2* emz* epz* eme* epe* xi*  epx 
      fi(25)= ena2* zeta*  epz* eme* epe* emx* epx 
      fi(26)=-ena2* emz* epz* eme* eta*  emx* epx 
      fi(27)= emz* epz* eme* epe* emx* epx 
     
      END
C----------------------------------------------------------------------C

      SUBROUTINE cshaped27(xi,eta,zeta,dfidxi,dfideta,dfidzeta,mpom)

C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Cell SHAPE functions for continuous quadr. 27-node cell         **
c **  -    -----                                 -                    **
c **********************************************************************
      INTEGER mpom
      REAL(8) xi,eta,zeta
      REAL(8) dfidxi(mpom),dfideta(mpom),dfidzeta(mpom)
      REAL(8) ena8,ena4,ena2,epz,emz,epe,eme,emx,epx
      REAL(8) e2p1,e2m1,x2p1,x2m1,z2p1,z2m1
      REAL(8) en2m1,xn2m1,zn2m1
      
      ena8=1.0D0/8.0D0
      ena4=1.0D0/4.0D0
      ena2=1.0D0/2.0D0
      epz=1.0D0+zeta
      emz=1.0D0-zeta
      epe=1.0D0+eta
      eme=1.0D0-eta
      epx=1.0D0+xi
      emx=1.0D0-xi
      
      e2p1=2.0D0*eta+1.0D0
      e2m1=2.0D0*eta-1.0D0      
      x2p1=2.0D0*xi+1.0D0
      x2m1=2.0D0*xi-1.0D0
      z2p1=2.0D0*zeta+1.0D0
      z2m1=2.0D0*zeta-1.0D0      
      
      en2m1=eta*eta-1.0D0
      xn2m1=xi*xi-1.0D0
      zn2m1=zeta*zeta-1.0D0            
c
c     d(fi)/d(xi)
c
      dfidxi( 1)=ena8 *zeta  *epz *eta  *epe *x2m1  
      dfidxi( 2)=- ena8 *emz *zeta  *eta  *epe *x2m1  
      dfidxi( 3)=- ena8 *emz *zeta  *eta  *epe *x2p1  
      dfidxi( 4)= ena8 *zeta  *epz *eta  *epe *x2p1  
      dfidxi( 5)=- ena8 *zeta  *epz *eme *eta  *x2m1  
      dfidxi( 6)= ena8 *emz *zeta  *eme *eta  *x2m1  
      dfidxi( 7)= ena8 *emz *zeta  *eme *eta  *x2p1  
      dfidxi( 8)=- ena8 *zeta  *epz *eme *eta  *x2p1  
      dfidxi( 9)= ena4 *emz *epz *eta  *epe *x2m1  
      dfidxi(10)= ena2 *emz *zeta  *eta  *epe *xi   
      dfidxi(11)= ena4 *emz *epz *eta  *epe *x2p1  
      dfidxi(12)= -ena2 *zeta  *epz *eta  *epe *xi   
      dfidxi(13)= ena4 *zeta  *epz *eme *epe *x2m1  
      dfidxi(14)= -ena4 *emz *zeta  *eme *epe *x2m1  
      dfidxi(15)= -ena4 *emz *zeta  *eme *epe *x2p1  
      dfidxi(16)= ena4 *zeta  *epz *eme *epe *x2p1  
      dfidxi(17)= -ena4 *emz *epz *eme *eta  *x2m1  
      dfidxi(18)= -ena2 *emz *zeta  *eme *eta  *xi   
      dfidxi(19)= -ena4 *emz *epz *eme *eta  *x2p1  
      dfidxi(20)= ena2 *zeta  *epz *eme *eta  *xi   
      dfidxi(21)= zn2m1 *eta  *epe *xi   
      dfidxi(22)= ena2 *zn2m1 *en2m1 *x2m1  
      dfidxi(23)= -emz *zeta  *en2m1 *xi   
      dfidxi(24)= ena2 *emz *epz *eme *epe *x2p1  
      dfidxi(25)= zeta  *epz *en2m1 *xi   
      dfidxi(26)= -zn2m1 *eme *eta  *xi   
      dfidxi(27)= -2.0D0 *zn2m1 *en2m1 *xi             
c
c     d(fi)/d(eta)
c
      dfideta( 1)= -ena8 *zeta  *epz *e2p1 *emx *xi   
      dfideta( 2)= ena8 *emz *zeta  *e2p1 *emx *xi   
      dfideta( 3)= -ena8 *emz *zeta  *e2p1 *xi  *epx  
      dfideta( 4)= ena8 *zeta  *epz *e2p1 *xi  *epx  
      dfideta( 5)= -ena8 *zeta  *epz *e2m1 *emx *xi   
      dfideta( 6)= ena8 *emz *zeta  *e2m1 *emx *xi   
      dfideta( 7)= -ena8 *emz *zeta  *e2m1 *xi  *epx  
      dfideta( 8)= ena8 *zeta  *epz *e2m1 *xi  *epx  
      dfideta( 9)= -ena4 *emz *epz *e2p1 *emx *xi   
      dfideta(10)= -ena4 *emz *zeta  *e2p1 *emx *epx  
      dfideta(11)= ena4 *emz *epz *e2p1 *xi  *epx  
      dfideta(12)= ena4 *zeta  *epz *e2p1 *emx *epx  
      dfideta(13)= ena2 *zeta  *epz *eta  *emx *xi   
      dfideta(14)= -ena2 *emz *zeta  *eta  *emx *xi   
      dfideta(15)= ena2 *emz *zeta  *eta  *xi  *epx  
      dfideta(16)= -ena2 *zeta  *epz *eta  *xi  *epx  
      dfideta(17)= -ena4 *emz *epz *e2m1 *emx *xi   
      dfideta(18)= -ena4 *emz *zeta  *e2m1 *emx *epx  
      dfideta(19)= ena4 *emz *epz *e2m1 *xi  *epx  
      dfideta(20)= ena4 *zeta  *epz *e2m1 *emx *epx  
      dfideta(21)= ena2 *emz *epz *e2p1 *emx *epx  
      dfideta(22)= -zn2m1 *eta  *emx *xi   
      dfideta(23)= -emz *zeta  *eta  *xn2m1  
      dfideta(24)= zn2m1 *eta  *xi  *epx  
      dfideta(25)= zeta  *epz *eta  *xn2m1  
      dfideta(26)= ena2 *zn2m1 *e2m1 *xn2m1  
      dfideta(27)= -2.0D0 *zn2m1 *eta  *xn2m1
c
c     d(fi)/d(zeta)
c
      dfidzeta( 1)= -ena8 *z2p1 *eta  *epe *emx *xi   
      dfidzeta( 2)= -ena8 *z2m1 *eta  *epe *emx *xi   
      dfidzeta( 3)= ena8 *z2m1 *eta  *epe *xi  *epx  
      dfidzeta( 4)= ena8 *z2p1 *eta  *epe *xi  *epx  
      dfidzeta( 5)= ena8 *z2p1 *eme *eta  *emx *xi   
      dfidzeta( 6)= ena8 *z2m1 *eme *eta *emx *xi   
      dfidzeta( 7)= -ena8 *z2m1 *eme *eta  *xi  *epx  
      dfidzeta( 8)= -ena8 *z2p1 *eme *eta  *xi  *epx  
      dfidzeta( 9)= ena2 *zeta  *eta  *epe *emx *xi   
      dfidzeta(10)= ena4 *z2m1 *eta  *epe *emx *epx  
      dfidzeta(11)= -ena2 *zeta  *eta  *epe *xi  *epx  
      dfidzeta(12)= ena4 *z2p1 *eta  *epe *emx *epx  
      dfidzeta(13)= -ena4 *z2p1 *eme *epe *emx *xi   
      dfidzeta(14)= -ena4 *z2m1 *eme *epe *emx *xi   
      dfidzeta(15)= ena4 *z2m1 *eme *epe *xi  *epx  
      dfidzeta(16)= ena4 *z2p1 *eme *epe *xi  *epx  
      dfidzeta(17)= -ena2 *zeta  *eme *eta  *emx *xi   
      dfidzeta(18)= -ena4 *z2m1 *eme *eta  *emx *epx  
      dfidzeta(19)= ena2 *zeta  *eme *eta  *xi  *epx  
      dfidzeta(20)= -ena4 *z2p1 *eme *eta  *emx *epx  
      dfidzeta(21)= zeta  *eta  *epe *xn2m1  
      dfidzeta(22)= -zeta  *en2m1 *emx *xi   
      dfidzeta(23)= ena2 *z2m1 *en2m1 *xn2m1  
      dfidzeta(24)= zeta  *en2m1 *xi  *epx  
      dfidzeta(25)= ena2 *z2p1 *eme *epe *emx *epx  
      dfidzeta(26)= -zeta  *eme *eta  *xn2m1  
      dfidzeta(27)= -2.0D0 *zeta  *en2m1 *xn2m1

      END

C----------------------------------------------------------------------C

      SUBROUTINE xieta34(xi,eta,n)

C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **            xi and eta for discontinuous linear 3/4  4-node       **
c **                                             -                    **
c **********************************************************************
      INTEGER n
      REAL(8) xi,eta,tc

      tc=3.0D0/4.0D0

      IF (n.EQ.1) THEN
        xi =-tc
        eta=-tc
      ELSE IF (n.EQ.2) THEN
        xi = tc
        eta=-tc
      ELSE IF (n.EQ.3) THEN
        xi = tc
        eta= tc
      ELSE IF (n.EQ.4) THEN
        xi =-tc
        eta= tc
      ELSE
        Print *,"ERROR, xieta34"
      END IF

      RETURN
      END
C----------------------------------------------------------------------C

      SUBROUTINE dl34shape4(fi,xsi,eta,mpom)

C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Cell SHAPE functions for discontinuous linear 3/4  4-node       **
c **                                             -                    **
c **********************************************************************
      INTEGER mpom
      REAL(8) xsi,eta,fi,sd,tc
      DIMENSION fi(mpom)

      sd=4.0D0/9.0D0
      tc=3.0D0/4.0D0

      fi(1)= sd*(xsi-tc)*(eta-tc)
      fi(2)=-sd*(xsi+tc)*(eta-tc)
      fi(3)= sd*(xsi+tc)*(eta+tc)
      fi(4)=-sd*(xsi-tc)*(eta+tc)

      RETURN
      END
C^L

C----------------------------------------------------------------------C

      SUBROUTINE cshape9(fi,xsi,eta,mpom)

C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Cell SHAPE functions for continuous quadr. 9-node cell    CQC9  **
c **  -    -----                                 -                    **
c **********************************************************************
      INTEGER mpom
      REAL(8) xsi,eta,fi
      DIMENSION fi(mpom)

      fi(1)=0.25D00*(xsi-xsi**2)*(eta-eta**2)
      fi(2)=0.50D00*(1.00D00-xsi**2)*(eta**2-eta)
      fi(3)=0.25D00*(xsi+xsi**2)*(eta**2-eta)
      fi(4)=0.50D00*(xsi+xsi**2)*(1.00D00-eta**2)
      fi(5)=0.25D00*(xsi+xsi**2)*(eta+eta**2)
      fi(6)=0.50D00*(1.00D00-xsi**2)*(eta**2+eta)
      fi(7)=0.25D00*(xsi**2-xsi)*(eta+eta**2)
      fi(8)=0.50D00*(xsi-xsi**2)*(eta**2-1.00D00)
      fi(9)=(1.00D00-xsi**2)*(1.00D00-eta**2)
      RETURN
      END
C^L     
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE cshder9(fk,fe,xsi,eta,mpom)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Cell SHape function DERivatives for CQC9                        **
c **  -    --             ---                -                        **
c **********************************************************************
      INTEGER mpom
      REAL(8) fk,fe,xsi,eta
      DIMENSION fk(mpom),fe(mpom)
c
C***  derivatives of shape functions w.r.t. local coordinate - xsi
c
      fk(1)=0.25*(1.0-2.0*xsi)*(eta-eta**2)
      fk(2)=-xsi*(eta**2-eta)
      fk(3)=0.25*(2.0*xsi+1.0)*(eta**2-eta)
      fk(4)=0.5*(1.0+2.0*xsi)*(1.0-eta**2)
      fk(5)=0.25*(1.0+2.0*xsi)*(eta+eta**2)
      fk(6)=-xsi*(eta**2+eta)
      fk(7)=0.25*(2.0*xsi-1.0)*(eta+eta**2)
      fk(8)=0.5*(1.0-2.0*xsi)*(eta**2-1.0)
      fk(9)=-2.0*xsi*(1.0-eta**2)
c
c*** derivatives of the shape functions w.r.t. local coordinate eta
c
      fe(1)=0.25*(xsi-xsi**2)*(1.0-2.0*eta)
      fe(2)=0.5*(1.0-xsi**2)*(2.0*eta-1.0)
      fe(3)=0.25*(xsi+xsi**2)*(2.0*eta-1.0)
      fe(4)=-eta*(xsi+xsi**2)
      fe(5)=0.25*(xsi+xsi**2)*(2.0*eta+1.0)
      fe(6)=0.5*(1.0-xsi**2)*(2.0*eta+1.0)
      fe(7)=0.25*(xsi**2-xsi)*(2.0*eta+1.0)
      fe(8)=eta*(xsi-xsi**2)
      fe(9)=-2.0*eta*(1.0-xsi**2)
      RETURN
      END      


C -----------------------------------------------------------------------------
      SUBROUTINE idc2ibc(icell,mesh,ibc)
C
C     $: Za eno celico zgradi ibc za robove
C
C -----------------------------------------------------------------------------
      USE inc_types 
      TYPE(meshType) :: mesh      
      INTEGER icell,i
      INTEGER ibc(mesh%nside,mesh%npob)
      
      INTEGER psp(9),pspr(9),pde(9),pza(9),ple(9),pzg(9)

c     rocno vtipkam vrstni red vozlisc v celici, tako da sestavim robne elemente      
c      DATA psp /1,8,7,6,5,4,3,2,9/
c      DATA pde /3,  4,  5,  14,  23,  22,  21,  12,  13/
c      DATA pza /5,  6,  7,  16,  25,  24,  23,  14,  15/
c      DATA ple /7,  8,  1,  10,  19,  26,  25,  16,  17/
c      DATA pspr /1,  2,  3,  12,  21,  20,  19,  10,  11/
c      DATA pzg /19,20,  21,  22,  23, 24,  25,  26,  27/
      
      DATA psp /6,14,2,10,3,15,7,18,23/
      DATA pspr/6,18,7,19,8,20,5,17,26/
      DATA pde /7,15,3,11,4,16,8,19,24/
      DATA pza /3,10,2,9,1,12,4,11,21/
      DATA ple /2,14,6,17,5,13,1,9,22/
      DATA pzg /5,20,8,16,4,12,1,13,25/
      
      
      DO i=1,mesh%npob
        ibc(1,i)=mesh%idc(icell,psp(i))
        ibc(2,i)=mesh%idc(icell,pspr(i))
        ibc(3,i)=mesh%idc(icell,pde(i))
        ibc(4,i)=mesh%idc(icell,pza(i))
        ibc(5,i)=mesh%idc(icell,ple(i))
        ibc(6,i)=mesh%idc(icell,pzg(i))    
      END DO
     
      END     


C -----------------------------------------------------------------------------
      SUBROUTINE zv4nez4(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                   qx1,qy1,qz1,qx2,qy2,qz2,qx3,qy3,qz3,qx4,qy4,qz4)
C
C     $: Prepacuna lokacije tock iz zvezne v nezvezno
C
C -----------------------------------------------------------------------------
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8) qx1,qy1,qz1,qx2,qy2,qz2,qx3,qy3,qz3,qx4,qy4,qz4
      REAL(8) mat(4,4),a,b,c
      
      a=49.0D0/64.0D0
      b=7.0D0/64.0D0
      c=1.0D0/64.0D0
      
      mat(1,1)=a
      mat(1,2)=b
      mat(1,3)=c
      mat(1,4)=b      

      mat(2,1)=b
      mat(2,2)=a
      mat(2,3)=b
      mat(2,4)=c      

      mat(3,1)=c
      mat(3,2)=b
      mat(3,3)=a
      mat(3,4)=b      

      mat(4,1)=b
      mat(4,2)=c
      mat(4,3)=b
      mat(4,4)=a      
      
      qx1=mat(1,1)*x1+mat(1,2)*x2+mat(1,3)*x3+mat(1,4)*x4
      qx2=mat(2,1)*x1+mat(2,2)*x2+mat(2,3)*x3+mat(2,4)*x4
      qx3=mat(3,1)*x1+mat(3,2)*x2+mat(3,3)*x3+mat(3,4)*x4
      qx4=mat(4,1)*x1+mat(4,2)*x2+mat(4,3)*x3+mat(4,4)*x4                  

      qy1=mat(1,1)*y1+mat(1,2)*y2+mat(1,3)*y3+mat(1,4)*y4
      qy2=mat(2,1)*y1+mat(2,2)*y2+mat(2,3)*y3+mat(2,4)*y4
      qy3=mat(3,1)*y1+mat(3,2)*y2+mat(3,3)*y3+mat(3,4)*y4
      qy4=mat(4,1)*y1+mat(4,2)*y2+mat(4,3)*y3+mat(4,4)*y4                  

      qz1=mat(1,1)*z1+mat(1,2)*z2+mat(1,3)*z3+mat(1,4)*z4
      qz2=mat(2,1)*z1+mat(2,2)*z2+mat(2,3)*z3+mat(2,4)*z4
      qz3=mat(3,1)*z1+mat(3,2)*z2+mat(3,3)*z3+mat(3,4)*z4
      qz4=mat(4,1)*z1+mat(4,2)*z2+mat(4,3)*z3+mat(4,4)*z4                  

      
      END


C -----------------------------------------------------------------------------
      SUBROUTINE setIfiji27(Ififj)
C
C     $: Z Mathematico izračunani integrali interpolacijskih funkcij
C        po enotskem elementu (-1,1)^3
C
C      = Integral ( fi(i) * fi(j) ) dOmega
C
C -----------------------------------------------------------------------------
      REAL(8) Ififj(27,27)

      Ififj(1,1) = 64.0D0 / 3375.0D0
      Ififj(1,2) = -16.0D0 / 3375.0D0
      Ififj(1,3) = 4.0D0 / 3375.0D0
      Ififj(1,4) = -16.0D0 / 3375.0D0
      Ififj(1,5) = -16.0D0 / 3375.0D0
      Ififj(1,6) = 4.0D0 / 3375.0D0
      Ififj(1,7) = -1.0D0 / 3375.0D0
      Ififj(1,8) = 4.0D0 / 3375.0D0
      Ififj(1,9) = 32.0D0 / 3375.0D0
      Ififj(1,10) = -8.0D0 / 3375.0D0
      Ififj(1,11) = -8.0D0 / 3375.0D0
      Ififj(1,12) = 32.0D0 / 3375.0D0
      Ififj(1,13) = 32.0D0 / 3375.0D0
      Ififj(1,14) = -8.0D0 / 3375.0D0
      Ififj(1,15) = 2.0D0 / 3375.0D0
      Ififj(1,16) = -8.0D0 / 3375.0D0
      Ififj(1,17) = -8.0D0 / 3375.0D0
      Ififj(1,18) = 2.0D0 / 3375.0D0
      Ififj(1,19) = 2.0D0 / 3375.0D0
      Ififj(1,20) = -8.0D0 / 3375.0D0
      Ififj(1,21) = 16.0D0 / 3375.0D0
      Ififj(1,22) = 16.0D0 / 3375.0D0
      Ififj(1,23) = -4.0D0 / 3375.0D0
      Ififj(1,24) = -4.0D0 / 3375.0D0
      Ififj(1,25) = 16.0D0 / 3375.0D0
      Ififj(1,26) = -4.0D0 / 3375.0D0
      Ififj(1,27) = 8.0D0 / 3375.0D0
      Ififj(2,1) = -16.0D0 / 3375.0D0
      Ififj(2,2) = 64.0D0 / 3375.0D0
      Ififj(2,3) = -16.0D0 / 3375.0D0
      Ififj(2,4) = 4.0D0 / 3375.0D0
      Ififj(2,5) = 4.0D0 / 3375.0D0
      Ififj(2,6) = -16.0D0 / 3375.0D0
      Ififj(2,7) = 4.0D0 / 3375.0D0
      Ififj(2,8) = -1.0D0 / 3375.0D0
      Ififj(2,9) = 32.0D0 / 3375.0D0
      Ififj(2,10) = 32.0D0 / 3375.0D0
      Ififj(2,11) = -8.0D0 / 3375.0D0
      Ififj(2,12) = -8.0D0 / 3375.0D0
      Ififj(2,13) = -8.0D0 / 3375.0D0
      Ififj(2,14) = 32.0D0 / 3375.0D0
      Ififj(2,15) = -8.0D0 / 3375.0D0
      Ififj(2,16) = 2.0D0 / 3375.0D0
      Ififj(2,17) = -8.0D0 / 3375.0D0
      Ififj(2,18) = -8.0D0 / 3375.0D0
      Ififj(2,19) = 2.0D0 / 3375.0D0
      Ififj(2,20) = 2.0D0 / 3375.0D0
      Ififj(2,21) = 16.0D0 / 3375.0D0
      Ififj(2,22) = 16.0D0 / 3375.0D0
      Ififj(2,23) = 16.0D0 / 3375.0D0
      Ififj(2,24) = -4.0D0 / 3375.0D0
      Ififj(2,25) = -4.0D0 / 3375.0D0
      Ififj(2,26) = -4.0D0 / 3375.0D0
      Ififj(2,27) = 8.0D0 / 3375.0D0
      Ififj(3,1) = 4.0D0 / 3375.0D0
      Ififj(3,2) = -16.0D0 / 3375.0D0
      Ififj(3,3) = 64.0D0 / 3375.0D0
      Ififj(3,4) = -16.0D0 / 3375.0D0
      Ififj(3,5) = -1.0D0 / 3375.0D0
      Ififj(3,6) = 4.0D0 / 3375.0D0
      Ififj(3,7) = -16.0D0 / 3375.0D0
      Ififj(3,8) = 4.0D0 / 3375.0D0
      Ififj(3,9) = -8.0D0 / 3375.0D0
      Ififj(3,10) = 32.0D0 / 3375.0D0
      Ififj(3,11) = 32.0D0 / 3375.0D0
      Ififj(3,12) = -8.0D0 / 3375.0D0
      Ififj(3,13) = 2.0D0 / 3375.0D0
      Ififj(3,14) = -8.0D0 / 3375.0D0
      Ififj(3,15) = 32.0D0 / 3375.0D0
      Ififj(3,16) = -8.0D0 / 3375.0D0
      Ififj(3,17) = 2.0D0 / 3375.0D0
      Ififj(3,18) = -8.0D0 / 3375.0D0
      Ififj(3,19) = -8.0D0 / 3375.0D0
      Ififj(3,20) = 2.0D0 / 3375.0D0
      Ififj(3,21) = 16.0D0 / 3375.0D0
      Ififj(3,22) = -4.0D0 / 3375.0D0
      Ififj(3,23) = 16.0D0 / 3375.0D0
      Ififj(3,24) = 16.0D0 / 3375.0D0
      Ififj(3,25) = -4.0D0 / 3375.0D0
      Ififj(3,26) = -4.0D0 / 3375.0D0
      Ififj(3,27) = 8.0D0 / 3375.0D0
      Ififj(4,1) = -16.0D0 / 3375.0D0
      Ififj(4,2) = 4.0D0 / 3375.0D0
      Ififj(4,3) = -16.0D0 / 3375.0D0
      Ififj(4,4) = 64.0D0 / 3375.0D0
      Ififj(4,5) = 4.0D0 / 3375.0D0
      Ififj(4,6) = -1.0D0 / 3375.0D0
      Ififj(4,7) = 4.0D0 / 3375.0D0
      Ififj(4,8) = -16.0D0 / 3375.0D0
      Ififj(4,9) = -8.0D0 / 3375.0D0
      Ififj(4,10) = -8.0D0 / 3375.0D0
      Ififj(4,11) = 32.0D0 / 3375.0D0
      Ififj(4,12) = 32.0D0 / 3375.0D0
      Ififj(4,13) = -8.0D0 / 3375.0D0
      Ififj(4,14) = 2.0D0 / 3375.0D0
      Ififj(4,15) = -8.0D0 / 3375.0D0
      Ififj(4,16) = 32.0D0 / 3375.0D0
      Ififj(4,17) = 2.0D0 / 3375.0D0
      Ififj(4,18) = 2.0D0 / 3375.0D0
      Ififj(4,19) = -8.0D0 / 3375.0D0
      Ififj(4,20) = -8.0D0 / 3375.0D0
      Ififj(4,21) = 16.0D0 / 3375.0D0
      Ififj(4,22) = -4.0D0 / 3375.0D0
      Ififj(4,23) = -4.0D0 / 3375.0D0
      Ififj(4,24) = 16.0D0 / 3375.0D0
      Ififj(4,25) = 16.0D0 / 3375.0D0
      Ififj(4,26) = -4.0D0 / 3375.0D0
      Ififj(4,27) = 8.0D0 / 3375.0D0
      Ififj(5,1) = -16.0D0 / 3375.0D0
      Ififj(5,2) = 4.0D0 / 3375.0D0
      Ififj(5,3) = -1.0D0 / 3375.0D0
      Ififj(5,4) = 4.0D0 / 3375.0D0
      Ififj(5,5) = 64.0D0 / 3375.0D0
      Ififj(5,6) = -16.0D0 / 3375.0D0
      Ififj(5,7) = 4.0D0 / 3375.0D0
      Ififj(5,8) = -16.0D0 / 3375.0D0
      Ififj(5,9) = -8.0D0 / 3375.0D0
      Ififj(5,10) = 2.0D0 / 3375.0D0
      Ififj(5,11) = 2.0D0 / 3375.0D0
      Ififj(5,12) = -8.0D0 / 3375.0D0
      Ififj(5,13) = 32.0D0 / 3375.0D0
      Ififj(5,14) = -8.0D0 / 3375.0D0
      Ififj(5,15) = 2.0D0 / 3375.0D0
      Ififj(5,16) = -8.0D0 / 3375.0D0
      Ififj(5,17) = 32.0D0 / 3375.0D0
      Ififj(5,18) = -8.0D0 / 3375.0D0
      Ififj(5,19) = -8.0D0 / 3375.0D0
      Ififj(5,20) = 32.0D0 / 3375.0D0
      Ififj(5,21) = -4.0D0 / 3375.0D0
      Ififj(5,22) = 16.0D0 / 3375.0D0
      Ififj(5,23) = -4.0D0 / 3375.0D0
      Ififj(5,24) = -4.0D0 / 3375.0D0
      Ififj(5,25) = 16.0D0 / 3375.0D0
      Ififj(5,26) = 16.0D0 / 3375.0D0
      Ififj(5,27) = 8.0D0 / 3375.0D0
      Ififj(6,1) = 4.0D0 / 3375.0D0
      Ififj(6,2) = -16.0D0 / 3375.0D0
      Ififj(6,3) = 4.0D0 / 3375.0D0
      Ififj(6,4) = -1.0D0 / 3375.0D0
      Ififj(6,5) = -16.0D0 / 3375.0D0
      Ififj(6,6) = 64.0D0 / 3375.0D0
      Ififj(6,7) = -16.0D0 / 3375.0D0
      Ififj(6,8) = 4.0D0 / 3375.0D0
      Ififj(6,9) = -8.0D0 / 3375.0D0
      Ififj(6,10) = -8.0D0 / 3375.0D0
      Ififj(6,11) = 2.0D0 / 3375.0D0
      Ififj(6,12) = 2.0D0 / 3375.0D0
      Ififj(6,13) = -8.0D0 / 3375.0D0
      Ififj(6,14) = 32.0D0 / 3375.0D0
      Ififj(6,15) = -8.0D0 / 3375.0D0
      Ififj(6,16) = 2.0D0 / 3375.0D0
      Ififj(6,17) = 32.0D0 / 3375.0D0
      Ififj(6,18) = 32.0D0 / 3375.0D0
      Ififj(6,19) = -8.0D0 / 3375.0D0
      Ififj(6,20) = -8.0D0 / 3375.0D0
      Ififj(6,21) = -4.0D0 / 3375.0D0
      Ififj(6,22) = 16.0D0 / 3375.0D0
      Ififj(6,23) = 16.0D0 / 3375.0D0
      Ififj(6,24) = -4.0D0 / 3375.0D0
      Ififj(6,25) = -4.0D0 / 3375.0D0
      Ififj(6,26) = 16.0D0 / 3375.0D0
      Ififj(6,27) = 8.0D0 / 3375.0D0
      Ififj(7,1) = -1.0D0 / 3375.0D0
      Ififj(7,2) = 4.0D0 / 3375.0D0
      Ififj(7,3) = -16.0D0 / 3375.0D0
      Ififj(7,4) = 4.0D0 / 3375.0D0
      Ififj(7,5) = 4.0D0 / 3375.0D0
      Ififj(7,6) = -16.0D0 / 3375.0D0
      Ififj(7,7) = 64.0D0 / 3375.0D0
      Ififj(7,8) = -16.0D0 / 3375.0D0
      Ififj(7,9) = 2.0D0 / 3375.0D0
      Ififj(7,10) = -8.0D0 / 3375.0D0
      Ififj(7,11) = -8.0D0 / 3375.0D0
      Ififj(7,12) = 2.0D0 / 3375.0D0
      Ififj(7,13) = 2.0D0 / 3375.0D0
      Ififj(7,14) = -8.0D0 / 3375.0D0
      Ififj(7,15) = 32.0D0 / 3375.0D0
      Ififj(7,16) = -8.0D0 / 3375.0D0
      Ififj(7,17) = -8.0D0 / 3375.0D0
      Ififj(7,18) = 32.0D0 / 3375.0D0
      Ififj(7,19) = 32.0D0 / 3375.0D0
      Ififj(7,20) = -8.0D0 / 3375.0D0
      Ififj(7,21) = -4.0D0 / 3375.0D0
      Ififj(7,22) = -4.0D0 / 3375.0D0
      Ififj(7,23) = 16.0D0 / 3375.0D0
      Ififj(7,24) = 16.0D0 / 3375.0D0
      Ififj(7,25) = -4.0D0 / 3375.0D0
      Ififj(7,26) = 16.0D0 / 3375.0D0
      Ififj(7,27) = 8.0D0 / 3375.0D0
      Ififj(8,1) = 4.0D0 / 3375.0D0
      Ififj(8,2) = -1.0D0 / 3375.0D0
      Ififj(8,3) = 4.0D0 / 3375.0D0
      Ififj(8,4) = -16.0D0 / 3375.0D0
      Ififj(8,5) = -16.0D0 / 3375.0D0
      Ififj(8,6) = 4.0D0 / 3375.0D0
      Ififj(8,7) = -16.0D0 / 3375.0D0
      Ififj(8,8) = 64.0D0 / 3375.0D0
      Ififj(8,9) = 2.0D0 / 3375.0D0
      Ififj(8,10) = 2.0D0 / 3375.0D0
      Ififj(8,11) = -8.0D0 / 3375.0D0
      Ififj(8,12) = -8.0D0 / 3375.0D0
      Ififj(8,13) = -8.0D0 / 3375.0D0
      Ififj(8,14) = 2.0D0 / 3375.0D0
      Ififj(8,15) = -8.0D0 / 3375.0D0
      Ififj(8,16) = 32.0D0 / 3375.0D0
      Ififj(8,17) = -8.0D0 / 3375.0D0
      Ififj(8,18) = -8.0D0 / 3375.0D0
      Ififj(8,19) = 32.0D0 / 3375.0D0
      Ififj(8,20) = 32.0D0 / 3375.0D0
      Ififj(8,21) = -4.0D0 / 3375.0D0
      Ififj(8,22) = -4.0D0 / 3375.0D0
      Ififj(8,23) = -4.0D0 / 3375.0D0
      Ififj(8,24) = 16.0D0 / 3375.0D0
      Ififj(8,25) = 16.0D0 / 3375.0D0
      Ififj(8,26) = 16.0D0 / 3375.0D0
      Ififj(8,27) = 8.0D0 / 3375.0D0
      Ififj(9,1) = 32.0D0 / 3375.0D0
      Ififj(9,2) = 32.0D0 / 3375.0D0
      Ififj(9,3) = -8.0D0 / 3375.0D0
      Ififj(9,4) = -8.0D0 / 3375.0D0
      Ififj(9,5) = -8.0D0 / 3375.0D0
      Ififj(9,6) = -8.0D0 / 3375.0D0
      Ififj(9,7) = 2.0D0 / 3375.0D0
      Ififj(9,8) = 2.0D0 / 3375.0D0
      Ififj(9,9) = 256.0D0 / 3375.0D0
      Ififj(9,10) = 16.0D0 / 3375.0D0
      Ififj(9,11) = -64.0D0 / 3375.0D0
      Ififj(9,12) = 16.0D0 / 3375.0D0
      Ififj(9,13) = 16.0D0 / 3375.0D0
      Ififj(9,14) = 16.0D0 / 3375.0D0
      Ififj(9,15) = -4.0D0 / 3375.0D0
      Ififj(9,16) = -4.0D0 / 3375.0D0
      Ififj(9,17) = -64.0D0 / 3375.0D0
      Ififj(9,18) = -4.0D0 / 3375.0D0
      Ififj(9,19) = 16.0D0 / 3375.0D0
      Ififj(9,20) = -4.0D0 / 3375.0D0
      Ififj(9,21) = 128.0D0 / 3375.0D0
      Ififj(9,22) = 128.0D0 / 3375.0D0
      Ififj(9,23) = 8.0D0 / 3375.0D0
      Ififj(9,24) = -32.0D0 / 3375.0D0
      Ififj(9,25) = 8.0D0 / 3375.0D0
      Ififj(9,26) = -32.0D0 / 3375.0D0
      Ififj(9,27) = 64.0D0 / 3375.0D0
      Ififj(10,1) = -8.0D0 / 3375.0D0
      Ififj(10,2) = 32.0D0 / 3375.0D0
      Ififj(10,3) = 32.0D0 / 3375.0D0
      Ififj(10,4) = -8.0D0 / 3375.0D0
      Ififj(10,5) = 2.0D0 / 3375.0D0
      Ififj(10,6) = -8.0D0 / 3375.0D0
      Ififj(10,7) = -8.0D0 / 3375.0D0
      Ififj(10,8) = 2.0D0 / 3375.0D0
      Ififj(10,9) = 16.0D0 / 3375.0D0
      Ififj(10,10) = 256.0D0 / 3375.0D0
      Ififj(10,11) = 16.0D0 / 3375.0D0
      Ififj(10,12) = -64.0D0 / 3375.0D0
      Ififj(10,13) = -4.0D0 / 3375.0D0
      Ififj(10,14) = 16.0D0 / 3375.0D0
      Ififj(10,15) = 16.0D0 / 3375.0D0
      Ififj(10,16) = -4.0D0 / 3375.0D0
      Ififj(10,17) = -4.0D0 / 3375.0D0
      Ififj(10,18) = -64.0D0 / 3375.0D0
      Ififj(10,19) = -4.0D0 / 3375.0D0
      Ififj(10,20) = 16.0D0 / 3375.0D0
      Ififj(10,21) = 128.0D0 / 3375.0D0
      Ififj(10,22) = 8.0D0 / 3375.0D0
      Ififj(10,23) = 128.0D0 / 3375.0D0
      Ififj(10,24) = 8.0D0 / 3375.0D0
      Ififj(10,25) = -32.0D0 / 3375.0D0
      Ififj(10,26) = -32.0D0 / 3375.0D0
      Ififj(10,27) = 64.0D0 / 3375.0D0
      Ififj(11,1) = -8.0D0 / 3375.0D0
      Ififj(11,2) = -8.0D0 / 3375.0D0
      Ififj(11,3) = 32.0D0 / 3375.0D0
      Ififj(11,4) = 32.0D0 / 3375.0D0
      Ififj(11,5) = 2.0D0 / 3375.0D0
      Ififj(11,6) = 2.0D0 / 3375.0D0
      Ififj(11,7) = -8.0D0 / 3375.0D0
      Ififj(11,8) = -8.0D0 / 3375.0D0
      Ififj(11,9) = -64.0D0 / 3375.0D0
      Ififj(11,10) = 16.0D0 / 3375.0D0
      Ififj(11,11) = 256.0D0 / 3375.0D0
      Ififj(11,12) = 16.0D0 / 3375.0D0
      Ififj(11,13) = -4.0D0 / 3375.0D0
      Ififj(11,14) = -4.0D0 / 3375.0D0
      Ififj(11,15) = 16.0D0 / 3375.0D0
      Ififj(11,16) = 16.0D0 / 3375.0D0
      Ififj(11,17) = 16.0D0 / 3375.0D0
      Ififj(11,18) = -4.0D0 / 3375.0D0
      Ififj(11,19) = -64.0D0 / 3375.0D0
      Ififj(11,20) = -4.0D0 / 3375.0D0
      Ififj(11,21) = 128.0D0 / 3375.0D0
      Ififj(11,22) = -32.0D0 / 3375.0D0
      Ififj(11,23) = 8.0D0 / 3375.0D0
      Ififj(11,24) = 128.0D0 / 3375.0D0
      Ififj(11,25) = 8.0D0 / 3375.0D0
      Ififj(11,26) = -32.0D0 / 3375.0D0
      Ififj(11,27) = 64.0D0 / 3375.0D0
      Ififj(12,1) = 32.0D0 / 3375.0D0
      Ififj(12,2) = -8.0D0 / 3375.0D0
      Ififj(12,3) = -8.0D0 / 3375.0D0
      Ififj(12,4) = 32.0D0 / 3375.0D0
      Ififj(12,5) = -8.0D0 / 3375.0D0
      Ififj(12,6) = 2.0D0 / 3375.0D0
      Ififj(12,7) = 2.0D0 / 3375.0D0
      Ififj(12,8) = -8.0D0 / 3375.0D0
      Ififj(12,9) = 16.0D0 / 3375.0D0
      Ififj(12,10) = -64.0D0 / 3375.0D0
      Ififj(12,11) = 16.0D0 / 3375.0D0
      Ififj(12,12) = 256.0D0 / 3375.0D0
      Ififj(12,13) = 16.0D0 / 3375.0D0
      Ififj(12,14) = -4.0D0 / 3375.0D0
      Ififj(12,15) = -4.0D0 / 3375.0D0
      Ififj(12,16) = 16.0D0 / 3375.0D0
      Ififj(12,17) = -4.0D0 / 3375.0D0
      Ififj(12,18) = 16.0D0 / 3375.0D0
      Ififj(12,19) = -4.0D0 / 3375.0D0
      Ififj(12,20) = -64.0D0 / 3375.0D0
      Ififj(12,21) = 128.0D0 / 3375.0D0
      Ififj(12,22) = 8.0D0 / 3375.0D0
      Ififj(12,23) = -32.0D0 / 3375.0D0
      Ififj(12,24) = 8.0D0 / 3375.0D0
      Ififj(12,25) = 128.0D0 / 3375.0D0
      Ififj(12,26) = -32.0D0 / 3375.0D0
      Ififj(12,27) = 64.0D0 / 3375.0D0
      Ififj(13,1) = 32.0D0 / 3375.0D0
      Ififj(13,2) = -8.0D0 / 3375.0D0
      Ififj(13,3) = 2.0D0 / 3375.0D0
      Ififj(13,4) = -8.0D0 / 3375.0D0
      Ififj(13,5) = 32.0D0 / 3375.0D0
      Ififj(13,6) = -8.0D0 / 3375.0D0
      Ififj(13,7) = 2.0D0 / 3375.0D0
      Ififj(13,8) = -8.0D0 / 3375.0D0
      Ififj(13,9) = 16.0D0 / 3375.0D0
      Ififj(13,10) = -4.0D0 / 3375.0D0
      Ififj(13,11) = -4.0D0 / 3375.0D0
      Ififj(13,12) = 16.0D0 / 3375.0D0
      Ififj(13,13) = 256.0D0 / 3375.0D0
      Ififj(13,14) = -64.0D0 / 3375.0D0
      Ififj(13,15) = 16.0D0 / 3375.0D0
      Ififj(13,16) = -64.0D0 / 3375.0D0
      Ififj(13,17) = 16.0D0 / 3375.0D0
      Ififj(13,18) = -4.0D0 / 3375.0D0
      Ififj(13,19) = -4.0D0 / 3375.0D0
      Ififj(13,20) = 16.0D0 / 3375.0D0
      Ififj(13,21) = 8.0D0 / 3375.0D0
      Ififj(13,22) = 128.0D0 / 3375.0D0
      Ififj(13,23) = -32.0D0 / 3375.0D0
      Ififj(13,24) = -32.0D0 / 3375.0D0
      Ififj(13,25) = 128.0D0 / 3375.0D0
      Ififj(13,26) = 8.0D0 / 3375.0D0
      Ififj(13,27) = 64.0D0 / 3375.0D0
      Ififj(14,1) = -8.0D0 / 3375.0D0
      Ififj(14,2) = 32.0D0 / 3375.0D0
      Ififj(14,3) = -8.0D0 / 3375.0D0
      Ififj(14,4) = 2.0D0 / 3375.0D0
      Ififj(14,5) = -8.0D0 / 3375.0D0
      Ififj(14,6) = 32.0D0 / 3375.0D0
      Ififj(14,7) = -8.0D0 / 3375.0D0
      Ififj(14,8) = 2.0D0 / 3375.0D0
      Ififj(14,9) = 16.0D0 / 3375.0D0
      Ififj(14,10) = 16.0D0 / 3375.0D0
      Ififj(14,11) = -4.0D0 / 3375.0D0
      Ififj(14,12) = -4.0D0 / 3375.0D0
      Ififj(14,13) = -64.0D0 / 3375.0D0
      Ififj(14,14) = 256.0D0 / 3375.0D0
      Ififj(14,15) = -64.0D0 / 3375.0D0
      Ififj(14,16) = 16.0D0 / 3375.0D0
      Ififj(14,17) = 16.0D0 / 3375.0D0
      Ififj(14,18) = 16.0D0 / 3375.0D0
      Ififj(14,19) = -4.0D0 / 3375.0D0
      Ififj(14,20) = -4.0D0 / 3375.0D0
      Ififj(14,21) = 8.0D0 / 3375.0D0
      Ififj(14,22) = 128.0D0 / 3375.0D0
      Ififj(14,23) = 128.0D0 / 3375.0D0
      Ififj(14,24) = -32.0D0 / 3375.0D0
      Ififj(14,25) = -32.0D0 / 3375.0D0
      Ififj(14,26) = 8.0D0 / 3375.0D0
      Ififj(14,27) = 64.0D0 / 3375.0D0
      Ififj(15,1) = 2.0D0 / 3375.0D0
      Ififj(15,2) = -8.0D0 / 3375.0D0
      Ififj(15,3) = 32.0D0 / 3375.0D0
      Ififj(15,4) = -8.0D0 / 3375.0D0
      Ififj(15,5) = 2.0D0 / 3375.0D0
      Ififj(15,6) = -8.0D0 / 3375.0D0
      Ififj(15,7) = 32.0D0 / 3375.0D0
      Ififj(15,8) = -8.0D0 / 3375.0D0
      Ififj(15,9) = -4.0D0 / 3375.0D0
      Ififj(15,10) = 16.0D0 / 3375.0D0
      Ififj(15,11) = 16.0D0 / 3375.0D0
      Ififj(15,12) = -4.0D0 / 3375.0D0
      Ififj(15,13) = 16.0D0 / 3375.0D0
      Ififj(15,14) = -64.0D0 / 3375.0D0
      Ififj(15,15) = 256.0D0 / 3375.0D0
      Ififj(15,16) = -64.0D0 / 3375.0D0
      Ififj(15,17) = -4.0D0 / 3375.0D0
      Ififj(15,18) = 16.0D0 / 3375.0D0
      Ififj(15,19) = 16.0D0 / 3375.0D0
      Ififj(15,20) = -4.0D0 / 3375.0D0
      Ififj(15,21) = 8.0D0 / 3375.0D0
      Ififj(15,22) = -32.0D0 / 3375.0D0
      Ififj(15,23) = 128.0D0 / 3375.0D0
      Ififj(15,24) = 128.0D0 / 3375.0D0
      Ififj(15,25) = -32.0D0 / 3375.0D0
      Ififj(15,26) = 8.0D0 / 3375.0D0
      Ififj(15,27) = 64.0D0 / 3375.0D0
      Ififj(16,1) = -8.0D0 / 3375.0D0
      Ififj(16,2) = 2.0D0 / 3375.0D0
      Ififj(16,3) = -8.0D0 / 3375.0D0
      Ififj(16,4) = 32.0D0 / 3375.0D0
      Ififj(16,5) = -8.0D0 / 3375.0D0
      Ififj(16,6) = 2.0D0 / 3375.0D0
      Ififj(16,7) = -8.0D0 / 3375.0D0
      Ififj(16,8) = 32.0D0 / 3375.0D0
      Ififj(16,9) = -4.0D0 / 3375.0D0
      Ififj(16,10) = -4.0D0 / 3375.0D0
      Ififj(16,11) = 16.0D0 / 3375.0D0
      Ififj(16,12) = 16.0D0 / 3375.0D0
      Ififj(16,13) = -64.0D0 / 3375.0D0
      Ififj(16,14) = 16.0D0 / 3375.0D0
      Ififj(16,15) = -64.0D0 / 3375.0D0
      Ififj(16,16) = 256.0D0 / 3375.0D0
      Ififj(16,17) = -4.0D0 / 3375.0D0
      Ififj(16,18) = -4.0D0 / 3375.0D0
      Ififj(16,19) = 16.0D0 / 3375.0D0
      Ififj(16,20) = 16.0D0 / 3375.0D0
      Ififj(16,21) = 8.0D0 / 3375.0D0
      Ififj(16,22) = -32.0D0 / 3375.0D0
      Ififj(16,23) = -32.0D0 / 3375.0D0
      Ififj(16,24) = 128.0D0 / 3375.0D0
      Ififj(16,25) = 128.0D0 / 3375.0D0
      Ififj(16,26) = 8.0D0 / 3375.0D0
      Ififj(16,27) = 64.0D0 / 3375.0D0
      Ififj(17,1) = -8.0D0 / 3375.0D0
      Ififj(17,2) = -8.0D0 / 3375.0D0
      Ififj(17,3) = 2.0D0 / 3375.0D0
      Ififj(17,4) = 2.0D0 / 3375.0D0
      Ififj(17,5) = 32.0D0 / 3375.0D0
      Ififj(17,6) = 32.0D0 / 3375.0D0
      Ififj(17,7) = -8.0D0 / 3375.0D0
      Ififj(17,8) = -8.0D0 / 3375.0D0
      Ififj(17,9) = -64.0D0 / 3375.0D0
      Ififj(17,10) = -4.0D0 / 3375.0D0
      Ififj(17,11) = 16.0D0 / 3375.0D0
      Ififj(17,12) = -4.0D0 / 3375.0D0
      Ififj(17,13) = 16.0D0 / 3375.0D0
      Ififj(17,14) = 16.0D0 / 3375.0D0
      Ififj(17,15) = -4.0D0 / 3375.0D0
      Ififj(17,16) = -4.0D0 / 3375.0D0
      Ififj(17,17) = 256.0D0 / 3375.0D0
      Ififj(17,18) = 16.0D0 / 3375.0D0
      Ififj(17,19) = -64.0D0 / 3375.0D0
      Ififj(17,20) = 16.0D0 / 3375.0D0
      Ififj(17,21) = -32.0D0 / 3375.0D0
      Ififj(17,22) = 128.0D0 / 3375.0D0
      Ififj(17,23) = 8.0D0 / 3375.0D0
      Ififj(17,24) = -32.0D0 / 3375.0D0
      Ififj(17,25) = 8.0D0 / 3375.0D0
      Ififj(17,26) = 128.0D0 / 3375.0D0
      Ififj(17,27) = 64.0D0 / 3375.0D0
      Ififj(18,1) = 2.0D0 / 3375.0D0
      Ififj(18,2) = -8.0D0 / 3375.0D0
      Ififj(18,3) = -8.0D0 / 3375.0D0
      Ififj(18,4) = 2.0D0 / 3375.0D0
      Ififj(18,5) = -8.0D0 / 3375.0D0
      Ififj(18,6) = 32.0D0 / 3375.0D0
      Ififj(18,7) = 32.0D0 / 3375.0D0
      Ififj(18,8) = -8.0D0 / 3375.0D0
      Ififj(18,9) = -4.0D0 / 3375.0D0
      Ififj(18,10) = -64.0D0 / 3375.0D0
      Ififj(18,11) = -4.0D0 / 3375.0D0
      Ififj(18,12) = 16.0D0 / 3375.0D0
      Ififj(18,13) = -4.0D0 / 3375.0D0
      Ififj(18,14) = 16.0D0 / 3375.0D0
      Ififj(18,15) = 16.0D0 / 3375.0D0
      Ififj(18,16) = -4.0D0 / 3375.0D0
      Ififj(18,17) = 16.0D0 / 3375.0D0
      Ififj(18,18) = 256.0D0 / 3375.0D0
      Ififj(18,19) = 16.0D0 / 3375.0D0
      Ififj(18,20) = -64.0D0 / 3375.0D0
      Ififj(18,21) = -32.0D0 / 3375.0D0
      Ififj(18,22) = 8.0D0 / 3375.0D0
      Ififj(18,23) = 128.0D0 / 3375.0D0
      Ififj(18,24) = 8.0D0 / 3375.0D0
      Ififj(18,25) = -32.0D0 / 3375.0D0
      Ififj(18,26) = 128.0D0 / 3375.0D0
      Ififj(18,27) = 64.0D0 / 3375.0D0
      Ififj(19,1) = 2.0D0 / 3375.0D0
      Ififj(19,2) = 2.0D0 / 3375.0D0
      Ififj(19,3) = -8.0D0 / 3375.0D0
      Ififj(19,4) = -8.0D0 / 3375.0D0
      Ififj(19,5) = -8.0D0 / 3375.0D0
      Ififj(19,6) = -8.0D0 / 3375.0D0
      Ififj(19,7) = 32.0D0 / 3375.0D0
      Ififj(19,8) = 32.0D0 / 3375.0D0
      Ififj(19,9) = 16.0D0 / 3375.0D0
      Ififj(19,10) = -4.0D0 / 3375.0D0
      Ififj(19,11) = -64.0D0 / 3375.0D0
      Ififj(19,12) = -4.0D0 / 3375.0D0
      Ififj(19,13) = -4.0D0 / 3375.0D0
      Ififj(19,14) = -4.0D0 / 3375.0D0
      Ififj(19,15) = 16.0D0 / 3375.0D0
      Ififj(19,16) = 16.0D0 / 3375.0D0
      Ififj(19,17) = -64.0D0 / 3375.0D0
      Ififj(19,18) = 16.0D0 / 3375.0D0
      Ififj(19,19) = 256.0D0 / 3375.0D0
      Ififj(19,20) = 16.0D0 / 3375.0D0
      Ififj(19,21) = -32.0D0 / 3375.0D0
      Ififj(19,22) = -32.0D0 / 3375.0D0
      Ififj(19,23) = 8.0D0 / 3375.0D0
      Ififj(19,24) = 128.0D0 / 3375.0D0
      Ififj(19,25) = 8.0D0 / 3375.0D0
      Ififj(19,26) = 128.0D0 / 3375.0D0
      Ififj(19,27) = 64.0D0 / 3375.0D0
      Ififj(20,1) = -8.0D0 / 3375.0D0
      Ififj(20,2) = 2.0D0 / 3375.0D0
      Ififj(20,3) = 2.0D0 / 3375.0D0
      Ififj(20,4) = -8.0D0 / 3375.0D0
      Ififj(20,5) = 32.0D0 / 3375.0D0
      Ififj(20,6) = -8.0D0 / 3375.0D0
      Ififj(20,7) = -8.0D0 / 3375.0D0
      Ififj(20,8) = 32.0D0 / 3375.0D0
      Ififj(20,9) = -4.0D0 / 3375.0D0
      Ififj(20,10) = 16.0D0 / 3375.0D0
      Ififj(20,11) = -4.0D0 / 3375.0D0
      Ififj(20,12) = -64.0D0 / 3375.0D0
      Ififj(20,13) = 16.0D0 / 3375.0D0
      Ififj(20,14) = -4.0D0 / 3375.0D0
      Ififj(20,15) = -4.0D0 / 3375.0D0
      Ififj(20,16) = 16.0D0 / 3375.0D0
      Ififj(20,17) = 16.0D0 / 3375.0D0
      Ififj(20,18) = -64.0D0 / 3375.0D0
      Ififj(20,19) = 16.0D0 / 3375.0D0
      Ififj(20,20) = 256.0D0 / 3375.0D0
      Ififj(20,21) = -32.0D0 / 3375.0D0
      Ififj(20,22) = 8.0D0 / 3375.0D0
      Ififj(20,23) = -32.0D0 / 3375.0D0
      Ififj(20,24) = 8.0D0 / 3375.0D0
      Ififj(20,25) = 128.0D0 / 3375.0D0
      Ififj(20,26) = 128.0D0 / 3375.0D0
      Ififj(20,27) = 64.0D0 / 3375.0D0
      Ififj(21,1) = 16.0D0 / 3375.0D0
      Ififj(21,2) = 16.0D0 / 3375.0D0
      Ififj(21,3) = 16.0D0 / 3375.0D0
      Ififj(21,4) = 16.0D0 / 3375.0D0
      Ififj(21,5) = -4.0D0 / 3375.0D0
      Ififj(21,6) = -4.0D0 / 3375.0D0
      Ififj(21,7) = -4.0D0 / 3375.0D0
      Ififj(21,8) = -4.0D0 / 3375.0D0
      Ififj(21,9) = 128.0D0 / 3375.0D0
      Ififj(21,10) = 128.0D0 / 3375.0D0
      Ififj(21,11) = 128.0D0 / 3375.0D0
      Ififj(21,12) = 128.0D0 / 3375.0D0
      Ififj(21,13) = 8.0D0 / 3375.0D0
      Ififj(21,14) = 8.0D0 / 3375.0D0
      Ififj(21,15) = 8.0D0 / 3375.0D0
      Ififj(21,16) = 8.0D0 / 3375.0D0
      Ififj(21,17) = -32.0D0 / 3375.0D0
      Ififj(21,18) = -32.0D0 / 3375.0D0
      Ififj(21,19) = -32.0D0 / 3375.0D0
      Ififj(21,20) = -32.0D0 / 3375.0D0
      Ififj(21,21) = 1024.0D0 / 3375.0D0
      Ififj(21,22) = 64.0D0 / 3375.0D0
      Ififj(21,23) = 64.0D0 / 3375.0D0
      Ififj(21,24) = 64.0D0 / 3375.0D0
      Ififj(21,25) = 64.0D0 / 3375.0D0
      Ififj(21,26) = -256.0D0 / 3375.0D0
      Ififj(21,27) = 512.0D0 / 3375.0D0
      Ififj(22,1) = 16.0D0 / 3375.0D0
      Ififj(22,2) = 16.0D0 / 3375.0D0
      Ififj(22,3) = -4.0D0 / 3375.0D0
      Ififj(22,4) = -4.0D0 / 3375.0D0
      Ififj(22,5) = 16.0D0 / 3375.0D0
      Ififj(22,6) = 16.0D0 / 3375.0D0
      Ififj(22,7) = -4.0D0 / 3375.0D0
      Ififj(22,8) = -4.0D0 / 3375.0D0
      Ififj(22,9) = 128.0D0 / 3375.0D0
      Ififj(22,10) = 8.0D0 / 3375.0D0
      Ififj(22,11) = -32.0D0 / 3375.0D0
      Ififj(22,12) = 8.0D0 / 3375.0D0
      Ififj(22,13) = 128.0D0 / 3375.0D0
      Ififj(22,14) = 128.0D0 / 3375.0D0
      Ififj(22,15) = -32.0D0 / 3375.0D0
      Ififj(22,16) = -32.0D0 / 3375.0D0
      Ififj(22,17) = 128.0D0 / 3375.0D0
      Ififj(22,18) = 8.0D0 / 3375.0D0
      Ififj(22,19) = -32.0D0 / 3375.0D0
      Ififj(22,20) = 8.0D0 / 3375.0D0
      Ififj(22,21) = 64.0D0 / 3375.0D0
      Ififj(22,22) = 1024.0D0 / 3375.0D0
      Ififj(22,23) = 64.0D0 / 3375.0D0
      Ififj(22,24) = -256.0D0 / 3375.0D0
      Ififj(22,25) = 64.0D0 / 3375.0D0
      Ififj(22,26) = 64.0D0 / 3375.0D0
      Ififj(22,27) = 512.0D0 / 3375.0D0
      Ififj(23,1) = -4.0D0 / 3375.0D0
      Ififj(23,2) = 16.0D0 / 3375.0D0
      Ififj(23,3) = 16.0D0 / 3375.0D0
      Ififj(23,4) = -4.0D0 / 3375.0D0
      Ififj(23,5) = -4.0D0 / 3375.0D0
      Ififj(23,6) = 16.0D0 / 3375.0D0
      Ififj(23,7) = 16.0D0 / 3375.0D0
      Ififj(23,8) = -4.0D0 / 3375.0D0
      Ififj(23,9) = 8.0D0 / 3375.0D0
      Ififj(23,10) = 128.0D0 / 3375.0D0
      Ififj(23,11) = 8.0D0 / 3375.0D0
      Ififj(23,12) = -32.0D0 / 3375.0D0
      Ififj(23,13) = -32.0D0 / 3375.0D0
      Ififj(23,14) = 128.0D0 / 3375.0D0
      Ififj(23,15) = 128.0D0 / 3375.0D0
      Ififj(23,16) = -32.0D0 / 3375.0D0
      Ififj(23,17) = 8.0D0 / 3375.0D0
      Ififj(23,18) = 128.0D0 / 3375.0D0
      Ififj(23,19) = 8.0D0 / 3375.0D0
      Ififj(23,20) = -32.0D0 / 3375.0D0
      Ififj(23,21) = 64.0D0 / 3375.0D0
      Ififj(23,22) = 64.0D0 / 3375.0D0
      Ififj(23,23) = 1024.0D0 / 3375.0D0
      Ififj(23,24) = 64.0D0 / 3375.0D0
      Ififj(23,25) = -256.0D0 / 3375.0D0
      Ififj(23,26) = 64.0D0 / 3375.0D0
      Ififj(23,27) = 512.0D0 / 3375.0D0
      Ififj(24,1) = -4.0D0 / 3375.0D0
      Ififj(24,2) = -4.0D0 / 3375.0D0
      Ififj(24,3) = 16.0D0 / 3375.0D0
      Ififj(24,4) = 16.0D0 / 3375.0D0
      Ififj(24,5) = -4.0D0 / 3375.0D0
      Ififj(24,6) = -4.0D0 / 3375.0D0
      Ififj(24,7) = 16.0D0 / 3375.0D0
      Ififj(24,8) = 16.0D0 / 3375.0D0
      Ififj(24,9) = -32.0D0 / 3375.0D0
      Ififj(24,10) = 8.0D0 / 3375.0D0
      Ififj(24,11) = 128.0D0 / 3375.0D0
      Ififj(24,12) = 8.0D0 / 3375.0D0
      Ififj(24,13) = -32.0D0 / 3375.0D0
      Ififj(24,14) = -32.0D0 / 3375.0D0
      Ififj(24,15) = 128.0D0 / 3375.0D0
      Ififj(24,16) = 128.0D0 / 3375.0D0
      Ififj(24,17) = -32.0D0 / 3375.0D0
      Ififj(24,18) = 8.0D0 / 3375.0D0
      Ififj(24,19) = 128.0D0 / 3375.0D0
      Ififj(24,20) = 8.0D0 / 3375.0D0
      Ififj(24,21) = 64.0D0 / 3375.0D0
      Ififj(24,22) = -256.0D0 / 3375.0D0
      Ififj(24,23) = 64.0D0 / 3375.0D0
      Ififj(24,24) = 1024.0D0 / 3375.0D0
      Ififj(24,25) = 64.0D0 / 3375.0D0
      Ififj(24,26) = 64.0D0 / 3375.0D0
      Ififj(24,27) = 512.0D0 / 3375.0D0
      Ififj(25,1) = 16.0D0 / 3375.0D0
      Ififj(25,2) = -4.0D0 / 3375.0D0
      Ififj(25,3) = -4.0D0 / 3375.0D0
      Ififj(25,4) = 16.0D0 / 3375.0D0
      Ififj(25,5) = 16.0D0 / 3375.0D0
      Ififj(25,6) = -4.0D0 / 3375.0D0
      Ififj(25,7) = -4.0D0 / 3375.0D0
      Ififj(25,8) = 16.0D0 / 3375.0D0
      Ififj(25,9) = 8.0D0 / 3375.0D0
      Ififj(25,10) = -32.0D0 / 3375.0D0
      Ififj(25,11) = 8.0D0 / 3375.0D0
      Ififj(25,12) = 128.0D0 / 3375.0D0
      Ififj(25,13) = 128.0D0 / 3375.0D0
      Ififj(25,14) = -32.0D0 / 3375.0D0
      Ififj(25,15) = -32.0D0 / 3375.0D0
      Ififj(25,16) = 128.0D0 / 3375.0D0
      Ififj(25,17) = 8.0D0 / 3375.0D0
      Ififj(25,18) = -32.0D0 / 3375.0D0
      Ififj(25,19) = 8.0D0 / 3375.0D0
      Ififj(25,20) = 128.0D0 / 3375.0D0
      Ififj(25,21) = 64.0D0 / 3375.0D0
      Ififj(25,22) = 64.0D0 / 3375.0D0
      Ififj(25,23) = -256.0D0 / 3375.0D0
      Ififj(25,24) = 64.0D0 / 3375.0D0
      Ififj(25,25) = 1024.0D0 / 3375.0D0
      Ififj(25,26) = 64.0D0 / 3375.0D0
      Ififj(25,27) = 512.0D0 / 3375.0D0
      Ififj(26,1) = -4.0D0 / 3375.0D0
      Ififj(26,2) = -4.0D0 / 3375.0D0
      Ififj(26,3) = -4.0D0 / 3375.0D0
      Ififj(26,4) = -4.0D0 / 3375.0D0
      Ififj(26,5) = 16.0D0 / 3375.0D0
      Ififj(26,6) = 16.0D0 / 3375.0D0
      Ififj(26,7) = 16.0D0 / 3375.0D0
      Ififj(26,8) = 16.0D0 / 3375.0D0
      Ififj(26,9) = -32.0D0 / 3375.0D0
      Ififj(26,10) = -32.0D0 / 3375.0D0
      Ififj(26,11) = -32.0D0 / 3375.0D0
      Ififj(26,12) = -32.0D0 / 3375.0D0
      Ififj(26,13) = 8.0D0 / 3375.0D0
      Ififj(26,14) = 8.0D0 / 3375.0D0
      Ififj(26,15) = 8.0D0 / 3375.0D0
      Ififj(26,16) = 8.0D0 / 3375.0D0
      Ififj(26,17) = 128.0D0 / 3375.0D0
      Ififj(26,18) = 128.0D0 / 3375.0D0
      Ififj(26,19) = 128.0D0 / 3375.0D0
      Ififj(26,20) = 128.0D0 / 3375.0D0
      Ififj(26,21) = -256.0D0 / 3375.0D0
      Ififj(26,22) = 64.0D0 / 3375.0D0
      Ififj(26,23) = 64.0D0 / 3375.0D0
      Ififj(26,24) = 64.0D0 / 3375.0D0
      Ififj(26,25) = 64.0D0 / 3375.0D0
      Ififj(26,26) = 1024.0D0 / 3375.0D0
      Ififj(26,27) = 512.0D0 / 3375.0D0
      Ififj(27,1) = 8.0D0 / 3375.0D0
      Ififj(27,2) = 8.0D0 / 3375.0D0
      Ififj(27,3) = 8.0D0 / 3375.0D0
      Ififj(27,4) = 8.0D0 / 3375.0D0
      Ififj(27,5) = 8.0D0 / 3375.0D0
      Ififj(27,6) = 8.0D0 / 3375.0D0
      Ififj(27,7) = 8.0D0 / 3375.0D0
      Ififj(27,8) = 8.0D0 / 3375.0D0
      Ififj(27,9) = 64.0D0 / 3375.0D0
      Ififj(27,10) = 64.0D0 / 3375.0D0
      Ififj(27,11) = 64.0D0 / 3375.0D0
      Ififj(27,12) = 64.0D0 / 3375.0D0
      Ififj(27,13) = 64.0D0 / 3375.0D0
      Ififj(27,14) = 64.0D0 / 3375.0D0
      Ififj(27,15) = 64.0D0 / 3375.0D0
      Ififj(27,16) = 64.0D0 / 3375.0D0
      Ififj(27,17) = 64.0D0 / 3375.0D0
      Ififj(27,18) = 64.0D0 / 3375.0D0
      Ififj(27,19) = 64.0D0 / 3375.0D0
      Ififj(27,20) = 64.0D0 / 3375.0D0
      Ififj(27,21) = 512.0D0 / 3375.0D0
      Ififj(27,22) = 512.0D0 / 3375.0D0
      Ififj(27,23) = 512.0D0 / 3375.0D0
      Ififj(27,24) = 512.0D0 / 3375.0D0
      Ififj(27,25) = 512.0D0 / 3375.0D0
      Ififj(27,26) = 512.0D0 / 3375.0D0
      Ififj(27,27) = 4096.0D0 / 3375.0D0

      END


