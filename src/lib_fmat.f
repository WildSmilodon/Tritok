
C -----------------------------------------------------------------------------
      SUBROUTINE FMATkmdcSteady(env,io,inp,mesh,gauss,smatH,smatG,smatB,
     &                          smatAbdx,smatAbdy,smatAbdz,
     &                          smatHtx,smatHty,smatHtz,
     &                          smatDx,smatDy,smatDz,
     &                          ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtyDx,ontHtzDy,ontHtxDz)
C
C     $: Form Matrices, Steady diffusion advection equation and kinematics eqaution
C        (H, G, Ab, Adx, Ady, Adz, B, Htx, Hty, Htz, Dx, Dy, Dz)
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(IOtype) io
      TYPE(inputtype) inp  
      TYPE(penv) :: env            
      
      INTEGER ic,it,je,dn,isrc,isip,i,j,cinidc,its,itf,ii
      REAL(8) xp,yp,zp,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
      REAL(8) x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8
      INTEGER, ALLOCATABLE :: ibc(:,:)
      REAL(8), ALLOCATABLE :: ge(:),he(:),be(:),dxe(:),dye(:),dze(:),abx(:),aby(:),abz(:)
      REAL(8), ALLOCATABLE :: htx(:),hty(:),htz(:)
      
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
  
      REAL(8) c,minedge
      REAL(8) ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz ! ocena natancnosti integralov
      REAL(8) ontHtyDx,ontHtzDy,ontHtxDz
      REAL(8) htykratz,htzkratx,htxkraty
      REAL(8) xi,eta,fi(mesh%npob)
      REAL(8) vsota,hv
      
      ontG=0.0D0
      ontH=0.0D0
      ontHtx=0.0D0
      ontHty=0.0D0
      ontHtz=0.0D0  
      ontDx=0.0D0
      ontDy=0.0D0
      ontDz=0.0D0  

C     ibc za makro element bom generiral sproti      
      ALLOCATE (ibc(mesh%nside,mesh%npob))
      
c     ker izracunam integrale pomnozene z int. funkcijami naenkrat rabim
      ALLOCATE (he(mesh%npob),ge(mesh%npof),be(mesh%npoc))
      ALLOCATE (abx(mesh%npob),aby(mesh%npob),abz(mesh%npob))
      ALLOCATE (dxe(mesh%npoc),dye(mesh%npoc),dze(mesh%npoc))
      ALLOCATE (htx(mesh%npob),hty(mesh%npob),htz(mesh%npob))     
      
      smatH=0.0D0
      smatG=0.0D0
      smatB=0.0D0
      smatAbdx=0.0D0
      smatAbdy=0.0D0
      smatAbdz=0.0D0                       
      smatHtx=0.0D0
      smatHty=0.0D0
      smatHtz=0.0D0                       
      smatDx=0.0D0
      smatDy=0.0D0
      smatDz=0.0D0   
      
      ontHtyDx=0.0D0  
      ontHtzDy=0.0D0
      ontHtxDz=0.0D0                  

           
C     Zanka po makro elementih (celicah)
      DO ic=1,mesh%nicell
        IF (ic.GT.1.AND.inp%copy.EQ.1) THEN
          it=0
          DO i=(ic-1)*mesh%nsp+1,ic*mesh%nsp
            it=it+1
            DO j=1,mesh%npoc ! H
              smatH(i,j)=smatH(it,j)                           
            END DO
            DO je=1,mesh%nside  ! G          
              DO isip=1,mesh%npof
                cinidc=(je-1)*mesh%npof+isip              
                smatG(i,cinidc)=smatG(it,cinidc)*DBLE(mesh%flpm(ic,je))*DBLE(mesh%flpm(1,je))
              END DO 
            END DO
            DO isip=1,mesh%npoc ! B
              smatB(i,isip)=smatB(it,isip)
              smatAbdx(i,isip)=smatAbdx(it,isip)
              smatAbdy(i,isip)=smatAbdy(it,isip)
              smatAbdz(i,isip)=smatAbdz(it,isip)
              smatDx(i,isip)=smatDx(it,isip)               
              smatDy(i,isip)=smatDy(it,isip)               
              smatDz(i,isip)=smatDz(it,isip)                                           
            END DO 
          END DO
          it=0
          DO i=(ic-1)*mesh%npoc+1,ic*mesh%npoc
            it=it+1
            DO isip=1,mesh%npoc ! Hti
              smatHtx(i,isip)=smatHtx(it,isip)               
              smatHty(i,isip)=smatHty(it,isip)               
              smatHtz(i,isip)=smatHtz(it,isip)                                           
            END DO 
          END DO
        ELSE
        CALL WarnErr(env,io,inp,5,"FMATdcSteady","integrate element",ic)
C
C       INTEGRACIJA PO ROBU
C
C       zgeneriram ibc za to celico
        CALL idc2ibc(ic,mesh,ibc)  
C       poisce najkrajsi rob
        CALL FindMinEdge(mesh,ibc,minedge)      
c        
C       Zanka po izvornih tockah v u
c
        DO it=1,mesh%npoc
          dn=mesh%idc(ic,it)
          xp=mesh%x(dn,1)
          yp=mesh%x(dn,2)
          zp=mesh%x(dn,3)
          i=(ic-1)*mesh%nsp+it ! vrstica
          ii=(ic-1)*mesh%npoc+it ! vrstica za htx,hty,htz
C         Zanka po robnih ploskvah makro elementa
          DO je=1,mesh%nside
            x1=mesh%x(ibc(je,1),1)
            y1=mesh%x(ibc(je,1),2)
            z1=mesh%x(ibc(je,1),3)
            x2=mesh%x(ibc(je,3),1)
            y2=mesh%x(ibc(je,3),2)
            z2=mesh%x(ibc(je,3),3)
            x3=mesh%x(ibc(je,5),1)
            y3=mesh%x(ibc(je,5),2)
            z3=mesh%x(ibc(je,5),3)
            x4=mesh%x(ibc(je,7),1)
            y4=mesh%x(ibc(je,7),2)
            z4=mesh%x(ibc(je,7),3)          
C           SET SOURCE POINT :
            isrc=0
            DO isip=1,mesh%npob
              IF (dn .EQ. ibc(je,isip)) isrc=isip
            END DO
C           integracija po robnih elementih celice            
            CALL INTEBc9kmdclap(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                          isrc,mesh%npob,mesh%npof,
     &                          he,ge,abx,aby,abz,htx,hty,htz,minedge)  
c           zlozimo H-je po matrikah            
            DO isip=1,mesh%npob
              j=ibc(je,isip) ! stolpec
c             poiscem kateri po vrsti je stolpec j v idcju celice (da lahko shranjuejm zgosceno)              
              CALL FindColInIdc(mesh,ic,j,cinidc)
c              matH(i,j)=matH(i,j)+he(isip)
              smatH(i,cinidc)=smatH(i,cinidc)+he(isip)             
              smatAbdx(i,cinidc)=smatAbdx(i,cinidc)+Abx(isip)              
              smatAbdy(i,cinidc)=smatAbdy(i,cinidc)+Aby(isip)              
              smatAbdz(i,cinidc)=smatAbdz(i,cinidc)+Abz(isip)
              smatHtx(ii,cinidc)=smatHtx(ii,cinidc)-Htx(isip)
              smatHty(ii,cinidc)=smatHty(ii,cinidc)-Hty(isip)
              smatHtz(ii,cinidc)=smatHtz(ii,cinidc)-Htz(isip)                                          
            END DO 
c           zlozim G-je po matrikah            
            DO isip=1,mesh%npof
              j=mesh%ibf(ic,je,isip) ! stolpec
c             ko zlagam v pravokotno matriko, jih dam kar po vrsti
              cinidc=(je-1)*mesh%npof+isip              
c             da nimam podvajanja fluxov na ploskvah, ki se dotikajo        
              smatG(i,cinidc)=smatG(i,cinidc)+ge(isip)*DBLE(mesh%flpm(ic,je))
            END DO                                  
          END DO
c         izracunam prosti koeficient c in ga dam na diagonalo
c         predpostavim togi premik telesa u=1, q=0          
          c=0.0D0
          DO j=1,mesh%npoc        
            c=c-smatH(i,j)         
          END DO
c         poiscem kateri po vrsti je stolpec j v idcju celice (da lahko shranjuejm zgosceno)              
          CALL FindColInIdc(mesh,ic,dn,cinidc)
          smatH(i,cinidc)=smatH(i,cinidc)+c
c         preverim se ht matrike, vsota vseh clenov mora biti 0
          vsota=0.0D0
          DO j=1,mesh%npoc        
            vsota=vsota+smatHtx(ii,j)
          END DO
          ontHtx=max(ontHtx,ABS(vsota))
          smathtx(ii,it)=smathtx(ii,it)-vsota  ! popravimo singularni integral ??XX??         
          vsota=0.0D0
          DO j=1,mesh%npoc        
            vsota=vsota+smatHty(ii,j)         
          END DO
          ontHty=max(ontHty,ABS(vsota))
          smathty(ii,it)=smathty(ii,it)-vsota  ! popravimo singularni integral          
          vsota=0.0D0
          DO j=1,mesh%npoc        
            vsota=vsota+smatHtz(ii,j)         
          END DO
          ontHtz=max(ontHtz,ABS(vsota))          
          smathtz(ii,it)=smathtz(ii,it)-vsota  ! popravimo singularni integral          
        END DO
        ontH=ontH+(ABS(1.0D0-c)) ! zadnja tocka je v sredini elementa
c
C       Zanka po izvornih tockah v q
c
        DO its=1,mesh%nside
          DO itf=1,mesh%npof
            dn=mesh%ibf(ic,its,itf)
            xp=mesh%xq(dn,1)
            yp=mesh%xq(dn,2)
            zp=mesh%xq(dn,3)
            i=(ic-1)*mesh%nsp+mesh%npoc+(its-1)*mesh%npof+itf ! vrstica
C           Zanka po robnih ploskvah makro elementa
            DO je=1,mesh%nside
              x1=mesh%x(ibc(je,1),1)
              y1=mesh%x(ibc(je,1),2)
              z1=mesh%x(ibc(je,1),3)
              x2=mesh%x(ibc(je,3),1)
              y2=mesh%x(ibc(je,3),2)
              z2=mesh%x(ibc(je,3),3)
              x3=mesh%x(ibc(je,5),1)
              y3=mesh%x(ibc(je,5),2)
              z3=mesh%x(ibc(je,5),3)
              x4=mesh%x(ibc(je,7),1)
              y4=mesh%x(ibc(je,7),2)
              z4=mesh%x(ibc(je,7),3)          
C             SET SOURCE POINT :
              IF (its.EQ.je) THEN ! integriram po tisti ploskvi v kateri je izvorna tocka
                isrc=itf+9  ! +9 zato, ker je tako v polju v rutini INTEBc9
              ELSE 
                isrc=0
              END IF
C             integracija po robnih elementih celice            
              CALL INTEBc9kmdclap(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                            isrc,mesh%npob,mesh%npof,
     &                            he,ge,abx,aby,abz,htx,hty,htz,minedge)  
c             zlozimo H-je po matrikah            
              DO isip=1,mesh%npob
                j=ibc(je,isip) ! stolpec
c               poiscem kateri po vrsti je stolpec j v idcju celice (da lahko shranjuejm zgosceno)              
                CALL FindColInIdc(mesh,ic,j,cinidc)
                smatH(i,cinidc)=smatH(i,cinidc)+he(isip)   
                smatAbdx(i,cinidc)=smatAbdx(i,cinidc)+Abx(isip)              
                smatAbdy(i,cinidc)=smatAbdy(i,cinidc)+Aby(isip)              
                smatAbdz(i,cinidc)=smatAbdz(i,cinidc)+Abz(isip) 

c                smatHtx(i,cinidc)=smatHtx(i,cinidc)-Htx(isip)
c                smatHty(i,cinidc)=smatHty(i,cinidc)-Hty(isip)
c                smatHtz(i,cinidc)=smatHtz(i,cinidc)-Htz(isip)                                          
                                             
              END DO 
c             zlozim G-je po matrikah            
              DO isip=1,mesh%npof
                j=mesh%ibf(ic,je,isip) ! stolpec
c               ko zlagam v pravokotno matriko, jih dam kar po vrsti
                cinidc=(je-1)*mesh%npof+isip              
c               da nimam podvajanja fluxov na ploskvah, ki se dotikajo        
                smatG(i,cinidc)=smatG(i,cinidc)+ge(isip)*DBLE(mesh%flpm(ic,je))
              END DO                                  
            END DO
c           izracunam prosti koeficient c 
c           predpostavim togi premik telesa u=1, q=0          
            c=0.0D0
            DO j=1,mesh%npoc        
              c=c-smatH(i,j)         
            END DO
c           ker je fluks nezvezen, je izvorna tocka vedno na ploskvi, zato 
c           vem, da je c=1/2. To uporabim za oceno natancnosti integracije
            ontG=ontG+ABS(0.5D0-c)
c           ker nimam neznanega u-ja v izvorni tocki za fluks, moram prispevek c
c           zinterpolirati na ostale 
c           ugotovimo xi in eta za izvorno tocko          
            CALL xieta34(xi,eta,itf)
c           izracunamo interpolacijske funkcije 9 tockovne zvezne v teh xi in eta
            CALL cshape9(fi,xi,eta,mesh%npob)            
c           dodamo del C ja k vsaki izmed devetih tock
            DO isip=1,mesh%npob
              j=ibc(its,isip) ! stolpec
c             poiscem kateri po vrsti je stolpec j v idcju celice (da lahko shranjuejm zgosceno)              
              CALL FindColInIdc(mesh,ic,j,cinidc)
c              smatH(i,cinidc)=smatH(i,cinidc)+0.5D0*fi(isip) ! to je slabse
              smatH(i,cinidc)=smatH(i,cinidc)+c*fi(isip)
            END DO   
c           ht matrike, ne vem zakaj vsota vseh clenov ni nic, ce je izvorna tocka v q                     
c           oziroma sume htz je nic na ploskvah 1 in 6, y 2 in 4, x pa 3 in 5
c           drugot pa 0.118643
c            vsota=0.0D0
c            DO j=1,mesh%npoc        
c              vsota=vsota+smatHtz(i,j)         
c            END DO
c            print *,vsota,"q",i
          END DO
        END DO
C
C       INTEGRACIJA PO OBMOCJU
C
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
c        
C       Zanka po izvornih tockah v u
c
        DO it=1,mesh%npoc
          dn=mesh%idc(ic,it)
          xp=mesh%x(dn,1)
          yp=mesh%x(dn,2)
          zp=mesh%x(dn,3)
          isrc=it ! source point ! VEDNO singularen integral !!!
          i=(ic-1)*mesh%nsp+it ! vrstica
          ii=(ic-1)*mesh%npoc+it ! vrstica za htx,hty,htz          
          CALL INTEDc27dcLap(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                  x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,be,dxe,dye,dze,isrc,mesh%npoc)
c         zlozimo vrstico po matrikah      
          DO isip=1,mesh%npoc
            smatB(i,isip)=be(isip)              
c            smatAbdx(i,isip)=smatAbdx(i,isip)-dxe(isip) ! to popravi, da po odstevalo sele ko bo imelo singularnega popravljanega
c            smatAbdy(i,isip)=smatAbdy(i,isip)-dye(isip)
c            smatAbdz(i,isip)=smatAbdz(i,isip)-dze(isip)                        
            smatDx(i,isip)=dxe(isip)
            smatDy(i,isip)=dye(isip)
            smatDz(i,isip)=dze(isip)
          END DO 
c         preverimo natancnost integralov
c         v=(z,0,0), w=(0,1,0)
          hv=0.0D0
          vsota=0.0D0
          DO j=1,mesh%npoc
            hv=hv+smatH(i,j)*mesh%x(mesh%idc(ic,j),3)                  
            vsota=vsota+smatDz(i,j)
          END DO
          ontDz=max(ontDz,ABS(hv-vsota))
          smatDz(i,it)=smatDz(i,it)+hv-vsota  ! popravimo singularni integral  ??XX??        
c         v=(0,x,0), w=(0,0,1)
          hv=0.0D0
          vsota=0.0D0
c          htykratz=0.0D0
          DO j=1,mesh%npoc
            hv=hv+smatH(i,j)*mesh%x(mesh%idc(ic,j),1)                  
            vsota=vsota+smatDx(i,j)
c            htykratz=htykratz+smatHty(ii,j)*mesh%x(mesh%idc(ic,j),3)
          END DO
c          ontHtyDx=max(ontHtyDx,ABS(htykratz+vsota))
          ontDx=max(ontDx,ABS(hv-vsota))
          smatDx(i,it)=smatDx(i,it)+hv-vsota  ! popravimo singularni integral  ??XX??
c         v=(0,0,y), w=(1,0,0)
          hv=0.0D0
          vsota=0.0D0
          DO j=1,mesh%npoc
            hv=hv+smatH(i,j)*mesh%x(mesh%idc(ic,j),2)                  
            vsota=vsota+smatDy(i,j)
          END DO
          ontDy=max(ontDy,ABS(hv-vsota))                  
          smatDy(i,it)=smatDy(i,it)+hv-vsota  ! popravimo singularni integral  ??XX??          
c       singularni popravljeni - se enkrat preverimo
          htxkraty=0.0D0
          htykratz=0.0D0      
          htzkratx=0.0D0          
          DO j=1,mesh%npoc            
            htxkraty=htxkraty+smatHtx(ii,j)*mesh%x(mesh%idc(ic,j),2)+smatDz(i,j)
            htykratz=htykratz+smatHty(ii,j)*mesh%x(mesh%idc(ic,j),3)+smatDx(i,j)     
            htzkratx=htzkratx+smatHtz(ii,j)*mesh%x(mesh%idc(ic,j),1)+smatDy(i,j)           
          END DO
          ontHtzDy=max(ontHtzDy,ABS(htzkratx))
          ontHtyDx=max(ontHtyDx,ABS(htykratz))
          ontHtxDz=max(ontHtxDz,ABS(htxkraty))
        END DO
c
C       Zanka po izvornih tockah v q
c
        DO its=1,mesh%nside
          DO itf=1,mesh%npof
            dn=mesh%ibf(ic,its,itf)
            xp=mesh%xq(dn,1)
            yp=mesh%xq(dn,2)
            zp=mesh%xq(dn,3)
            i=(ic-1)*mesh%nsp+mesh%npoc+(its-1)*mesh%npof+itf ! vrstica
            isrc=27+(its-1)*mesh%npof+itf
            CALL INTEDc27dcLap(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                  x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,be,dxe,dye,dze,isrc,mesh%npoc)
c           zlozimo vrstico po matrikah            
            DO isip=1,mesh%npoc
              smatB(i,isip)=be(isip)              
c            smatAbdx(i,isip)=smatAbdx(i,isip)-dxe(isip)
c            smatAbdy(i,isip)=smatAbdy(i,isip)-dye(isip)
c            smatAbdz(i,isip)=smatAbdz(i,isip)-dze(isip)
            smatDx(i,isip)=dxe(isip)
            smatDy(i,isip)=dye(isip)
            smatDz(i,isip)=dze(isip)                                                        
            END DO 
c         v=(z,0,0), w=(0,1,0) ! to spila !!
c          hv=0.0D0
c          vsota=0.0D0
c          DO j=1,mesh%npoc
c            hv=hv+smatH(i,j)*mesh%x(mesh%idc(ic,j),3)                  
c            vsota=vsota+smatDz(i,j)
c          END DO
c          print *,hv-vsota,"q"
          END DO
        END DO
        END IF
      END DO    
      IF (inp%copy.EQ.1) THEN
        ontG=ontG/mesh%npofc     
      ELSE
        ontG=ontG/mesh%nicell/mesh%npofc
        ontH=ontH/mesh%nicell
      END IF

      DEALLOCATE (ibc,ge,he,be,abx,aby,abz,dxe,dye,dze,htx,hty,htz)
      
      END

