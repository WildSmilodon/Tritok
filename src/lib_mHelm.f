C -----------------------------------------------------------------------------
      SUBROUTINE FMATmodHelm(env,io,inp,mesh,gauss,mHdiff,mHqDiff,smatH,smatG,smatB,
     &                          smatAbdx,smatAbdy,smatAbdz,
     &                          smatHtx,smatHty,smatHtz,
     &                          smatDx,smatDy,smatDz,
     &                          ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtyDx,ontHtzDy,ontHtxDz,mumax)
C
C     $: Form Matrices, modified Helmholtz eqaution
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

      REAL(8) mHdiff(mesh%nnodes),mHqDiff(mesh%nq),mu
      
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
      REAL(8) vsota,hv,vB,vH,mumax

c     ocene natancnosti integralov
      ontG=0.0D0
      ontH=0.0D0
      ontHtx=0.0D0
      ontHty=0.0D0
      ontHtz=0.0D0
      ontDx=0.0D0
      ontDy=0.0D0
      ontDz=0.0D0
      ontHtzDy=0.0D0
      ontHtyDx=0.0D0
      ontHtxDz=0.0D0
      mumax=0.D0


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
        CALL WarnErr(env,io,inp,5,"FMATmodHelm","integrate element",ic)
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
          mu=SQRT(inp%beta/mHdiff(dn))
          mumax=MAX(mu,mumax)
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
            CALL INTEBc9modHelm(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                          isrc,mesh%npob,mesh%npof,he,ge,abx,aby,abz,htx,hty,htz,minedge,mu)
c           zlozimo H-je po matrikah            
            DO isip=1,mesh%npob
              j=ibc(je,isip) ! stolpec
c             poiscem kateri po vrsti je stolpec j v idcju celice (da lahko shranjuejm zgosceno)              
              CALL FindColInIdc(mesh,ic,j,cinidc)
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
        END DO
c
C       Zanka po izvornih tockah v q
c
        DO its=1,mesh%nside
          DO itf=1,mesh%npof
            dn=mesh%ibf(ic,its,itf)
            mu=SQRT(inp%beta/mHqDiff(dn))
            mumax=MAX(mu,mumax)
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
              CALL INTEBc9modHelm(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                          isrc,mesh%npob,mesh%npof,he,ge,abx,aby,abz,htx,hty,htz,minedge,mu)
c             zlozimo H-je po matrikah            
              DO isip=1,mesh%npob
                j=ibc(je,isip) ! stolpec
c               poiscem kateri po vrsti je stolpec j v idcju celice (da lahko shranjuejm zgosceno)              
                CALL FindColInIdc(mesh,ic,j,cinidc)
                smatH(i,cinidc)=smatH(i,cinidc)+he(isip)
                smatAbdx(i,cinidc)=smatAbdx(i,cinidc)+Abx(isip)
                smatAbdy(i,cinidc)=smatAbdy(i,cinidc)+Aby(isip)
                smatAbdz(i,cinidc)=smatAbdz(i,cinidc)+Abz(isip)
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
          END DO
        END DO
C
C       INTEGRACIJA PO OBMOCJU **************************************************************************
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
          mu=SQRT(inp%beta/mHdiff(dn))
          mumax=MAX(mu,mumax)
          xp=mesh%x(dn,1)
          yp=mesh%x(dn,2)
          zp=mesh%x(dn,3)
          isrc=it ! source point ! VEDNO singularen integral !!!
          i=(ic-1)*mesh%nsp+it ! vrstica
          ii=(ic-1)*mesh%npoc+it ! vrstica za htx,hty,htz          
          CALL INTEDc27modHelm(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                  x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,be,dxe,dye,dze,isrc,mesh%npoc,mu)
c         zlozimo vrstico po matrikah      
          DO isip=1,mesh%npoc
            smatB(i,isip)=be(isip)
            smatDx(i,isip)=dxe(isip)
            smatDy(i,isip)=dye(isip)
            smatDz(i,isip)=dze(isip)
          END DO 
        END DO
c
C       Zanka po izvornih tockah v q
c
        DO its=1,mesh%nside
          DO itf=1,mesh%npof
            dn=mesh%ibf(ic,its,itf)
            mu=SQRT(inp%beta/mHqDiff(dn))
            mumax=MAX(mu,mumax)
            xp=mesh%xq(dn,1)
            yp=mesh%xq(dn,2)
            zp=mesh%xq(dn,3)
            i=(ic-1)*mesh%nsp+mesh%npoc+(its-1)*mesh%npof+itf ! vrstica
            isrc=27+(its-1)*mesh%npof+itf
            CALL INTEDc27modHelm(gauss,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                  x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,be,dxe,dye,dze,isrc,mesh%npoc,mu)
c           zlozimo vrstico po matrikah            
            DO isip=1,mesh%npoc
              smatB(i,isip)=be(isip)
              smatDx(i,isip)=dxe(isip)
              smatDy(i,isip)=dye(isip)
              smatDz(i,isip)=dze(isip)
            END DO 
          END DO
        END DO

C       ****************************************************************
C                    Correct for singular integrals
C       ****************************************************************
C
c        
C       Zanka po izvornih tockah v u
c
        DO it=1,mesh%npoc
          dn=mesh%idc(ic,it)
          mu=SQRT(inp%beta/mHdiff(dn))
          i=(ic-1)*mesh%nsp+it ! vrstica
          ii=(ic-1)*mesh%npoc+it ! vrstica za htx,hty,htz
c         izracunam prosti koeficient c in ga dam na diagonalo
c         predpostavim togi premik telesa u=1, q=0, f=-mu^2
          vH=0.0D0
          DO j=1,mesh%npoc
            vH=Vh+smatH(i,j)
          END DO
          vB=0.0D0
          DO j=1,mesh%npoc        
            vB=vB+smatB(i,j)
          END DO
          c = - vH + mu*mu*vB
c         poiscem kateri po vrsti je stolpec j v idcju celice (da lahko shranjuejm zgosceno)              
          CALL FindColInIdc(mesh,ic,dn,cinidc)
          smatH(i,cinidc)=smatH(i,cinidc)+c

c         preverim se ht matrike, vsota vseh clenov mora biti 0 ZAKAJ ??
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
            mu=SQRT(inp%beta/mHqDiff(dn))
            i=(ic-1)*mesh%nsp+mesh%npoc+(its-1)*mesh%npof+itf ! vrstica
c           izracunam prosti koeficient c 
c           predpostavim togi premik telesa u=1, q=0, f=-mu^2
            vH=0.0D0
            DO j=1,mesh%npoc
              vH=Vh+smatH(i,j)
            END DO
            vB=0.0D0
            DO j=1,mesh%npoc
              vB=vB+smatB(i,j)
            END DO
            c = - vH + mu*mu*vB
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
          END DO
        END DO
C       OBMOCNI *************
c
C       Zanka po izvornih tockah v u
c
        DO it=1,mesh%npoc
          dn=mesh%idc(ic,it)
          mu=SQRT(inp%beta/mHdiff(dn))
          i=(ic-1)*mesh%nsp+it ! vrstica
          ii=(ic-1)*mesh%npoc+it ! vrstica za htx,hty,htz

c         preverimo natancnost integralov
c         v=(z,0,0), w=(0,1,0)
          hv=0.0D0
          vsota=0.0D0
          DO j=1,mesh%npoc
            hv=hv+(smatH(i,j)-mu*mu*smatB(i,j) )*mesh%x(mesh%idc(ic,j),3)
            vsota=vsota+smatDz(i,j)
          END DO
          ontDz=max(ontDz,ABS(hv-vsota))
          smatDz(i,it)=smatDz(i,it)+hv-vsota  ! popravimo singularni integral  ??XX??
c         v=(0,x,0), w=(0,0,1)
          hv=0.0D0
          vsota=0.0D0
          DO j=1,mesh%npoc
            hv=hv+(smatH(i,j)-mu*mu*smatB(i,j))*mesh%x(mesh%idc(ic,j),1)
            vsota=vsota+smatDx(i,j)
          END DO
          ontDx=max(ontDx,ABS(hv-vsota))
          smatDx(i,it)=smatDx(i,it)+hv-vsota  ! popravimo singularni integral  ??XX??
c         v=(0,0,y), w=(1,0,0)
          hv=0.0D0
          vsota=0.0D0
          DO j=1,mesh%npoc
            hv=hv+(smatH(i,j)-mu*mu*smatB(i,j))*mesh%x(mesh%idc(ic,j),2)
            vsota=vsota+smatDy(i,j)
          END DO
          ontDy=max(ontDy,ABS(hv-vsota))
          smatDy(i,it)=smatDy(i,it)+hv-vsota  ! popravimo singularni integral  ??XX??

c         singularni popravljeni - se enkrat preverimo
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


        END IF  ! copy / no copy
      END DO  ! nicell   




      IF (inp%copy.EQ.1) THEN
        ontG=ontG/mesh%npofc     
      ELSE
        ontG=ontG/mesh%nicell/mesh%npofc
        ontH=ontH/mesh%nicell
      END IF

      DEALLOCATE (ibc,ge,he,be)

      END



C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsB_modHelm_fill2(mesh,smatH,smatG,diff,qdiff,sysm,rhsm)
C
C     $: Iz pravokotnih matrik g, h naredi CRS sistemsko in rhs
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      REAL(8) smatH(mesh%nicnsp,mesh%npoc),smatG(mesh%nicnsp,mesh%npofc)
      REAL(8) diff(mesh%nnodes)
      REAL(8) qdiff(mesh%nq)
      REAL(8) val
      TYPE(matrix) :: sysm,rhsm

      INTEGER wqkode(mesh%nq)
      INTEGER wkode(mesh%nnodes)
      INTEGER nunk, nb

      INTEGER ic,it,row,col,ii,jj,nuq,j,is,ir

      wqkode=mesh%mHqkode
      wkode=mesh%mHkode
      nunk=mesh%mHnunk
      nb=mesh%mHnb

      is=0
      ir=0
      DO ic=1,mesh%nicell ! zanka po celicah
        DO it=1,mesh%nsp ! po izvornih tockah znotraj celice
          row=(ic-1)*mesh%nsp+it ! vrstica v sistemski matriki
C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            val=diff(j)*smatH(row,col)
            nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              is=is+1
              sysm%v(sysm%d(is))=val
            ELSE ! znana vrednost - rhs
              ir=ir+1
              rhsm%v(rhsm%d(ir))=-val  ! minus zato, ker gre na drugo stran enacbe
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
                is=is+1
                sysm%v(sysm%d(is))=-smatG(row,col)*qDiff(j)
              ELSE ! znana vrednost - rhs
                ir=ir+1
                rhsm%v(rhsm%d(ir))=smatG(row,col)*qDiff(j) ! minus zato, ker gre na drugo stran enacbe
              END IF
            END DO
          END DO
        END DO
      END DO

      END









C -----------------------------------------------------------------------------
C
C       cal RMS for the modified Helmholtz equation test cases
C
      SUBROUTINE modHelm_calRMS(inp,gauss,mesh,mHu)
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(inputtype) inp
      TYPE(gausstype) :: gauss
      REAL(8), ALLOCATABLE ::  anal(:)
      REAL(8) mHu(mesh%nnodes)
      REAL(8) s2,s,rms

      INTEGER i

      ALLOCATE (anal(mesh%nnodes))

      IF (inp%mHf.EQ.1) THEN
        DO i=1,mesh%nnodes
          anal(i)=1.0D0-EXP(-inp%mHmu*mesh%x(i,1))
        END DO
      ELSE IF (inp%mHf.EQ.2.OR.inp%mHf.EQ.3) THEN
        DO i=1,mesh%nnodes
          anal(i)=mesh%x(i,1)**2+mesh%x(i,2)**2+mesh%x(i,3)**2
        END DO
      ELSE
        print *,"ERRRRR"
      END IF

      s=0.0D0
      s2=0.0D0

      DO i=1,mesh%nnodes
        IF (mesh%lbn(i).EQ.0) THEN
          s2=s2 + (anal(i)-mHu(i))**2
          s =s  + (anal(i))**2
c          print *,anal(i),mHu(i)
        END IF
      END DO

      rms = SQRT(s2/s)

      OPEN (68,FILE='tri.mH-rms.txt',STATUS='UNKNOWN',position="APPEND")
      WRITE (68,'(F13.3,I3,I3,G20.10,A20)') inp%mHmu,gauss%kmBr,gauss%iiDiv,rms,TRIM(mesh%FileName)
      CLOSE(68)

      DEALLOCATE (anal)


      END


C -----------------------------------------------------------------------------
C
C       cal RMS for the modified Helmholtz equation test cases
C
      SUBROUTINE modHelm_calRMS2(io,inp,gauss,mesh,mHu,anal,rtime)
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(IOtype)    :: io
      TYPE(meshType) :: mesh
      TYPE(inputtype) inp
      TYPE(gausstype) :: gauss
      REAL(8) mHu(mesh%nnodes)
      REAL(8) s2,s,rms,rtime
      REAL(8) anal(mesh%nnodes)

      INTEGER i

      s=0.0D0
      s2=0.0D0

      DO i=1,mesh%nnodes
        IF (mesh%lbn(i).EQ.0) THEN
          s2=s2 + (anal(i)-mHu(i))**2
          s =s  + (anal(i))**2
        END IF
      END DO

      rms = SQRT(s2/s)

      WRITE (io%mhr,'(2G20.10)') rtime,rms



      END



C -----------------------------------------------------------------------------
C
C       diffusivity of the modified Helmholtz equation
C
      SUBROUTINE modHelm_SetProblem(inp,mesh,diff,mHrhs,mHtime,mHu,mHq,fmqmH,rtime,velocity,anal)
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(inputtype) inp
      REAL(8) diff(mesh%nnodes)
      REAL(8) mHrhs(mesh%nnodes)
      REAL(8) mHtime(mesh%nnodes)
      REAL(8) mHu(mesh%nnodes)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) rtime
      REAL(8) x,y,z,t,f,nu
      REAL(8) pi

      REAL(8) anal(mesh%nnodes)
      REAL(8) mHq(mesh%nq)
      REAL(8) fmqmH(mesh%nbnpof)
      REAL(8), ALLOCATABLE :: dudx(:),dudy(:),dudz(:)

      INTEGER i,ii,dn,j

      ALLOCATE (dudx(mesh%nq),dudy(mesh%nq),dudz(mesh%nq))
      dudx=0.0D0 ! za tiste, kjer ne podam
      dudy=0.0D0
      dudz=0.0D0

      pi=4.0D0*ATAN(1.0D0)

      t=rtime
c     time derivatives are always on RHS

      IF (inp%mHf.EQ.6) THEN
        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
c         set diffusivity values
          diff(i)= 1 + x*y*z
c         set velocity values
          velocity(i,1)=x
          velocity(i,2)=y
          velocity(i,3)=-2.0D0*z
c         rhs sources
          f=-5.0d0 + 2.0d0*x**2 + 2.0d0*y**2 - 12.0d0*x*y*z -4.0d0* z**2
c         analytical solution
          anal(i) = x**2 + y**2 + z**2 + t
c         set
          mHrhs(i)=mHtime(i)-f
        END DO
        DO i=1,mesh%nq
          x=mesh%xq(i,1)
          y=mesh%xq(i,2)
          z=mesh%xq(i,3)
          dudx(i)=2.0D0*x
          dudy(i)=2.0d0*y
          dudz(i)=2.0D0*z
        END DO

      ELSE IF (inp%mHf.GE.7.AND.inp%mHf.LE.9) THEN
        IF (inp%mHf.EQ.7) THEN
            nu = 0.5D0*pi
        ELSE IF (inp%mHf.EQ.8) THEN
            nu = 1.5D0*pi
        ELSE IF (inp%mHf.EQ.9) THEN
            nu = 2.5D0*pi
        END IF

        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
c         set diffusivity values
          diff(i)= 2.0D0 + sin(nu*x*y*z)
c         set velocity values
          velocity(i,1)=x
          velocity(i,2)=y
          velocity(i,3)=-2.0D0*z
c         rhs sources
          f=-11.0D0+2.0D0*x**2+2.0D0*y**2-4.0D0*z**2-6.0D0*x*y*z*nu*Cos(x*y*z*nu)-6.0D0*Sin(x*y*z*nu)


c         analytical solution
          anal(i) = x**2 + y**2 + z**2 + t
c         set
          mHrhs(i)=mHtime(i)-f
        END DO
        DO i=1,mesh%nq
          x=mesh%xq(i,1)
          y=mesh%xq(i,2)
          z=mesh%xq(i,3)
          dudx(i)=2.0D0*x
          dudy(i)=2.0d0*y
          dudz(i)=2.0D0*z
        END DO

      ELSE IF (inp%mHf.GE.10.AND.inp%mHf.LE.12) THEN
        IF (inp%mHf.EQ.10) THEN
            nu = 2.0D0
        ELSE IF (inp%mHf.EQ.11) THEN
            nu = 5.0D0
        ELSE IF (inp%mHf.EQ.12) THEN
            nu = 10.0D0
        END IF

        DO i=1,mesh%nnodes
          x=mesh%x(i,1)
          y=mesh%x(i,2)
          z=mesh%x(i,3)
c         set diffusivity values
          diff(i)= 1.0D0 + x*y*z
c         set velocity values
          velocity(i,1)=nu
          velocity(i,2)=1.0D0
          velocity(i,3)=1.0D0
c         rhs sources
          IF (inp%mHf.EQ.10) THEN
            f=1+nu*(-1.0D0+nu-x*nu+x*y*z*nu)
          ELSE
            f=1+nu*(-1.0D0+nu-x*nu+x*y*z*nu)*x**(nu-2.0D0)
          END IF
c         analytical solution
          anal(i) = t-x**nu
c         set
          mHrhs(i)=mHtime(i)-f
        END DO
        DO i=1,mesh%nq
          x=mesh%xq(i,1)
          y=mesh%xq(i,2)
          z=mesh%xq(i,3)
          dudx(i)=-nu*x**(nu-1.0D0)
          dudy(i)=0.0D0
          dudz(i)=0.0D0
        END DO

      ELSE
        print *,"ERRRRR"
        STOP
      END IF

C     Dirichlet, kode=0
      DO i=1,mesh%nbnodes
        dn=mesh%gbn(i)
        IF (mesh%mHkode(dn).GT.0) THEN ! to je znana vrednost, mesh%iwn(i) je stevilka stene
          mHu(dn)=anal(dn)
        END IF
      END DO

C     Neumann, kode=2
      ii=0
      DO i=1,mesh%nbelem
        DO j=1,mesh%npof
          ii=ii+1
          dn=mesh%ibcf(i,j)
          IF (mesh%mHqkode(dn).GT.0) THEN ! to je znana vrednost, mesh%bewn(i) je stevilka stene
            IF (mesh%bewn(i).EQ.1) THEN ! stena pri z=1
              mHq(dn)=-dudz(dn)
            END IF
            IF (mesh%bewn(i).EQ.6) THEN ! stena pri z=2
              mHq(dn)=dudz(dn)
            END IF
            IF (mesh%bewn(i).EQ.5) THEN ! stena pri x=1
              mHq(dn)=-dudx(dn)
            END IF
            IF (mesh%bewn(i).EQ.3) THEN ! stena pri x=2
              mHq(dn)=dudx(dn)
            END IF
            IF (mesh%bewn(i).EQ.2) THEN ! stena pri y=1
              mHq(dn)=-dudy(dn)
            END IF
            IF (mesh%bewn(i).EQ.4) THEN ! stena pri y=2
              mHq(dn)=dudy(dn)
            END IF
            fmqmH(ii)=mHq(dn)

          END IF
        END DO
      END DO

      DEALLOCATE (dudx,dudy,dudz)


      END


C -----------------------------------------------------------------------------
C
C       diffusivity of the modified Helmholtz equation
C
      SUBROUTINE modHelm_SetProblemOld(inp,mesh,diff,mHrhs,mHtime,mHu,rtime,velocity,anal)
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(inputtype) inp
      REAL(8) diff(mesh%nnodes)
      REAL(8) mHrhs(mesh%nnodes)
      REAL(8) mHtime(mesh%nnodes)
      REAL(8) mHu(mesh%nnodes)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) rtime
      REAL(8) x,y,z,t,f,nu

      REAL(8) anal(mesh%nnodes)

      INTEGER i

      t=rtime
c     time derivatives are always on RHS
      DO i=1,mesh%nnodes
        mHrhs(i)=mHtime(i)
        x=mesh%x(i,1)
        y=mesh%x(i,2)
        z=mesh%x(i,3)

        IF (inp%mHf.EQ.1) THEN
c         set diffusivity values
          diff(i)=inp%beta  / ( inp%mHmu *inp%mHmu )
c         set velocity values
          velocity(i,1)=0.0D0
          velocity(i,2)=0.0D0
          velocity(i,3)=0.0D0
c         sources on rhs
          f = 1.0D0
c         set bc
          anal(i)=x+y+z+t

        ELSE IF (inp%mHf.EQ.2) THEN
c         set diffusivity values
          diff(i)=inp%beta  / ( inp%mHmu *inp%mHmu )
c         set velocity values
          velocity(i,1)=1.0D0
          velocity(i,2)=2.0D0
          velocity(i,3)=3.0D0
c         sources on rhs
          f=7.0D0
c         set bc
          anal(i) = x+y+z+t
        ELSE IF (inp%mHf.EQ.3) THEN
c         set diffusivity values
          diff(i)=  1.0D0
c         set velocity values
          velocity(i,1)=1.0D0
          velocity(i,2)=2.0D0
          velocity(i,3)=3.0D0
c         sources on rhs
          f=7.0D0
c         set bc
          anal(i)=x+y+z+t

        ELSE IF (inp%mHf.EQ.4) THEN
c         set diffusivity values
          diff(i)=  2.0D0
c         set velocity values
          velocity(i,1)=1.0D0
          velocity(i,2)=2.0D0
          velocity(i,3)=3.0D0
c         rhs sources
          f=-11.0D0 + 2.0D0*x + 4.0D0*y + 6.0D0*z
c         analytical solution
          anal(i) = x**2 + y**2 + z**2 + t

        ELSE IF (inp%mHf.EQ.5) THEN
c         set diffusivity values
          diff(i)=  2.0D0
c         set velocity values
          velocity(i,1)=x
          velocity(i,2)=y
          velocity(i,3)=-2.0D0*z
c         rhs sources
          f=-11.0D0 + 2.0D0*x**2 + 2.0D0*y**2 - 4.0D0*z**2
c         analytical solution
          anal(i) = x**2 + y**2 + z**2 + t

        ELSE IF (inp%mHf.EQ.6) THEN
c         set diffusivity values
          diff(i)= 1 + x*y*z
c         set velocity values
          velocity(i,1)=x
          velocity(i,2)=y
          velocity(i,3)=-2.0D0*z
c         rhs sources
          f=-5.0d0 + 2.0d0*x**2 + 2.0d0*y**2 - 12.0d0*x*y*z -4.0d0* z**2
c         analytical solution
          anal(i) = x**2 + y**2 + z**2 + t


        ELSE IF (inp%mHf.EQ.7) THEN
c         set diffusivity values
          nu = inp%mHmu
          diff(i)= 1.0D0
c         set velocity values
          velocity(i,1)=1.0D0
          velocity(i,2)=1.0D0
          velocity(i,3)=1.0D0
c         rhs sources
          f=(t + x + y + z)**(-2.0D0 +nu)*(3.0D0 + 4.0D0*t + 4.0D0*x + 4.0D0*y + 4.0D0*z - 3.0D0*nu)*nu
c         analytical solution
          anal(i) = (x + y + z + t)**nu

        ELSE
          print *,"ERRRRR"
          STOP
        END IF

c       set bc
        IF (mesh%lbn(i).NE.0) mHu(i)=anal(i)
c       sources on rhs
        mHrhs(i)=mHrhs(i)-f

      END DO

      END


C -----------------------------------------------------------------------------
C
C       R.H.S. of the modified Helmholtz equation
C
      SUBROUTINE modHelm_SetUpRHS(inp,mesh,mHrhs,mHtime,mHu,mHdiff)
C -----------------------------------------------------------------------------
      USE inc_types   

      TYPE(meshType) :: mesh
      TYPE(inputtype) inp
      REAL(8) mHrhs(mesh%nnodes)
      REAL(8) mHtime(mesh%nnodes)
      REAL(8) mHu(mesh%nnodes)
      REAL(8) mHdiff(mesh%nnodes)

      REAL(8) x,y,z,mu,anal

      INTEGER i

c     time derivatives are always on RHS
      DO i=1,mesh%nnodes
        mHrhs(i)=mHtime(i)
        x=mesh%x(i,1)
        y=mesh%x(i,2)
        z=mesh%x(i,3)
        mu=SQRT(inp%beta/mHdiff(i)) !

        IF (inp%mHf.EQ.1) THEN
            mHrhs(i)=-mu*mu
            anal=1.0D0-EXP(-mu*x)
        ELSE IF (inp%mHf.EQ.2) THEN
            mHrhs(i)=6.0D0-mu*mu*(x**2+y**2+z**2)
            anal=x**2+y**2+z**2
        ELSE IF (inp%mHf.EQ.3) THEN
            mHrhs(i)=6.0D0-mu*mu*(x**2+y**2+z**2)
            anal=x**2+y**2+z**2
        ELSE
          print *,"ERRRRR"
        END IF

        IF (mesh%lbn(i).NE.0) mHu(i)=anal

      END DO

      END
C -----------------------------------------------------------------------------
C
C       Solve modified Helmholtz equation, mu can vary with position
C
C
      SUBROUTINE SolveModHelm(env,io,inp,cpu,mesh,mHdiff,mHqDiff,sysm,rhsm,precv,
     &                      u,qu,rhsv,smatB,nit)
C
C     $: resi H*u=G*q+B*rhsv
C
C -----------------------------------------------------------------------------
      USE inc_types   

      TYPE(meshType) :: mesh
      TYPE(matrix) :: sysm,rhsm
      TYPE(IOtype) io
      TYPE(inputtype) inp
      TYPE(CPUtype) cpu     
      TYPE(penv) :: env                     
      
      INTEGER i,ic,it,j,nit,ierr
      INTEGER node,dn
      
      REAL(8) mHdiff(mesh%nnodes),mu,mun
      REAL(8) mHqDiff(mesh%nq)

      REAL(8) rhsv(mesh%nnodes)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc)
               
      REAL(8) precv(mesh%mHnunk)
      REAL(8) u(mesh%nnodes)
      REAL(8) qu(mesh%nq) 
  
      REAL(8), ALLOCATABLE :: x(:), b(:), r(:)
      REAL(4) cptime,rcpu      
      
      ALLOCATE (x(mesh%mHnunk),b(mesh%mHnb),r(mesh%nicnsp))

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
      DO i=1,mesh%nq
        IF (mesh%mHqkode(i).LT.0) THEN
          x(ABS(mesh%mHqkode(i)))=qu(i)
        ELSE
          b(mesh%mhqkode(i))=qu(i)        
        END IF
      END DO

C
C     r = rhs * b
C      
      CALL CRSxV(rhsm,b,mesh%mHnb,r)

c
c     r = r - B * rhsv
c
c      DO j=1,mesh%npoc
c        i=0
c        DO ic=1,mesh%nicell
c          DO it=1,51 ! 51 izvornih tock v vsaki celici
c            i=i+1  ! vrstica v matriki
c              node=mesh%idc(ic,j)
c              r(i)=r(i) - smatB(i,j)*rhsv(node)
c          END DO
c        END DO
c      END DO


c
c     r = r - B * rhsv
c
      DO j=1,mesh%npoc
        i=0
        DO ic=1,mesh%nicell

          node=mesh%idc(ic,j)
          mun=SQRT(inp%beta/mHdiff(node)) ! mu v obmocni tocki
c
C         Zanka po izvornih tockah v u
c
          DO it=1,mesh%npoc
            dn=mesh%idc(ic,it)
            mu=SQRT(inp%beta/mHdiff(dn)) ! mu v izvorni tocki (ta je sel v integral)
            i=i+1  ! vrstica v matriki
            r(i)=r(i) - smatB(i,j)*rhsv(node)
     &                - smatB(i,j)*(mun*mun-mu*mu)*u(node) ! zaradi variabilnega mu
          END DO
c
C         Zanka po izvornih tockah v q
c
          DO its=1,mesh%nside
            DO itf=1,mesh%npof
              dn=mesh%ibf(ic,its,itf)
              mu=SQRT(inp%beta/mHqDiff(dn)) ! mu v izvorni tocki (ta je sel v integral)
              i=i+1  ! vrstica v matriki
              r(i)=r(i) - smatB(i,j)*rhsv(node)
     &                  - smatB(i,j)*(mun*mun-mu*mu)*u(node) ! zaradi variabilnega mu
            END DO
          END DO
        END DO
      END DO




      cpu%time(23)=cpu%time(23)+cptime(rcpu)
      rcpu=cptime(0.)

C
C     solve overdetermined system of linear equations
C      
      CALL slvlsqr2(sysm%neq,mesh%mHnunk,sysm%nnz,inp%lsqs_maxit,inp%dlse(9),nit,ierr,
     &              precv,sysm%i,sysm%j,sysm%v,r,x)
     
      cpu%time(24)=cpu%time(24)+cptime(rcpu)

      IF (ierr.NE.0) CALL WarnErr(env,io,inp,4,"SolveModHelm","NAPAKA V SOVLERJU",ierr)
C
C     Copy solution vector to q and u using under-relaxation
c     under-relaxation should be 1.0 for linear problems
C      
      DO i=1,mesh%nnodes
        IF (mesh%mHkode(i).LT.0) THEN
          u(i)=(1.0D0-inp%urmH)*u(i)+inp%urmH*x(ABS(mesh%mHkode(i)))
        END IF
      END DO      
      DO i=1,mesh%nq
        IF (mesh%mHqkode(i).LT.0) THEN
          qu(i)=x(ABS(mesh%mHqkode(i)))
        END IF
      END DO

      DEALLOCATE (x,b,r)
      
      END
      
      

C -----------------------------------------------------------------------------
C
C       Solve diffusion advection equation using finite difference discretization
C       of time derivative where diffusivity varies in space
C       use modified Helmholtz fundamental solution and
C       overdeterminet system of equations
C
C -----------------------------------------------------------------------------
C
      SUBROUTINE SolveModHelmDC(env,io,inp,cpu,mesh,sysm,rhsm,precv,
     &                      u,qu,velocity,rhsv,smatB,
     &                      smatAbdx,smatAbdy,smatAbdz,smatDx,smatDy,smatDz,nit,
     &                      diff,qDiff,diffGrad)
C
C     $: resi H*u=G*q+B*rhsv-advekcija
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(matrix) :: sysm,rhsm
      TYPE(IOtype) io
      TYPE(inputtype) inp
      TYPE(CPUtype) cpu
      TYPE(penv) :: env

      INTEGER i,ic,it,j,nit,ierr
      INTEGER node,dn

      REAL(8) diff(mesh%nnodes),mu,mun
      REAL(8) qDiff(mesh%nq)
      REAL(8) diffGrad(mesh%nnodes,3)

      REAL(8) rhsv(mesh%nnodes)
      REAL(8) smatB(mesh%nicnsp,mesh%npoc)
      REAL(8) smatAbdX(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part X
      REAL(8) smatAbdY(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Y
      REAL(8) smatAbdZ(mesh%nicnsp,mesh%npoc) ! advection - boundary -  part Z
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc)
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc)

      REAL(8) precv(mesh%mHnunk)
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) u(mesh%nnodes)
      REAL(8) qu(mesh%nq)



      REAL(8), ALLOCATABLE :: x(:), b(:), r(:)
      REAL(4) cptime,rcpu

      ALLOCATE (x(mesh%mHnunk),b(mesh%mHnb),r(mesh%nicnsp))

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
      DO i=1,mesh%nq
        IF (mesh%mHqkode(i).LT.0) THEN
          x(ABS(mesh%mHqkode(i)))=qu(i)
        ELSE
          b(mesh%mhqkode(i))=qu(i)
        END IF
      END DO

C
C     r = rhs * b
C
      CALL CRSxV(rhsm,b,mesh%mHnb,r)
c
c     r = r - B * rhsv
c
      DO j=1,mesh%npoc
        i=0
        DO ic=1,mesh%nicell

          node=mesh%idc(ic,j) ! stolpec v matriki
          mun=SQRT(inp%beta/diff(node)) ! mu v obmocni tocki
c
C         Zanka po izvornih tockah v u
c
          DO it=1,mesh%npoc
            dn=mesh%idc(ic,it)
            mu=SQRT(inp%beta/diff(dn)) ! mu v izvorni tocki (ta je sel v integral)
            i=i+1  ! vrstica v matriki
            r(i)=r(i) - smatB(i,j)*rhsv(node)
     &               -u(node)*(
C                    advection
     &               +smatAbdx(i,j)*velocity(node,1)
     &               +smatAbdy(i,j)*velocity(node,2)
     &               +smatAbdz(i,j)*velocity(node,3)
     &                 -smatDx(i,j)*(velocity(node,1)+diffGrad(node,1))
     &                 -smatDy(i,j)*(velocity(node,2)+diffGrad(node,2))
     &                 -smatDz(i,j)*(velocity(node,3)+diffGrad(node,3))
     &                  -smatB(i,j)*inp%beta*( (mu/mun)**2 - 1.0D0)
     &                )
          END DO
c
C         Zanka po izvornih tockah v q
c
          DO its=1,mesh%nside
            DO itf=1,mesh%npof
              dn=mesh%ibf(ic,its,itf)
              mu=SQRT(inp%beta/qDiff(dn)) ! mu v izvorni tocki (ta je sel v integral)
              i=i+1  ! vrstica v matriki
              r(i)=r(i) - smatB(i,j)*rhsv(node)
     &               -u(node)*(
C                    advection
     &               +smatAbdx(i,j)*velocity(node,1)
     &               +smatAbdy(i,j)*velocity(node,2)
     &               +smatAbdz(i,j)*velocity(node,3)
     &                 -smatDx(i,j)*(velocity(node,1)+diffGrad(node,1))
     &                 -smatDy(i,j)*(velocity(node,2)+diffGrad(node,2))
     &                 -smatDz(i,j)*(velocity(node,3)+diffGrad(node,3))
     &                  -smatB(i,j)*inp%beta*( (mu/mun)**2 - 1.0D0)
     &                )
            END DO
          END DO
        END DO
      END DO




      cpu%time(23)=cpu%time(23)+cptime(rcpu)
      rcpu=cptime(0.)

C
C     solve overdetermined system of linear equations
C
      CALL slvlsqr2(sysm%neq,mesh%mHnunk,sysm%nnz,inp%lsqs_maxit,inp%dlse(9),nit,ierr,
     &              precv,sysm%i,sysm%j,sysm%v,r,x)

      cpu%time(24)=cpu%time(24)+cptime(rcpu)

      IF (ierr.NE.0) CALL WarnErr(env,io,inp,4,"SolveModHelmDC","NAPAKA V SOVLERJU",ierr)
C
C     Copy solution vector to q and u using under-relaxation
c     under-relaxation should be 1.0 for linear problems
C
      DO i=1,mesh%nnodes
        IF (mesh%mHkode(i).LT.0) THEN
          u(i)=(1.0D0-inp%urmH)*u(i)+inp%urmH*x(ABS(mesh%mHkode(i)))
        END IF
      END DO
      DO i=1,mesh%nq
        IF (mesh%mHqkode(i).LT.0) THEN
          qu(i)=x(ABS(mesh%mHqkode(i)))
        END IF
      END DO

      DEALLOCATE (x,b,r)

      END





C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsB_modHelm_fill(mesh,smatH,smatG,sysm,rhsm)      
C
C     $: Iz pravokotnih matrik g, h in b naredi CRS sistemsko in rhs
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      REAL(8) smatH(mesh%nicnsp,mesh%npoc),smatG(mesh%nicnsp,mesh%npofc)
      REAL(8) val
      TYPE(matrix) :: sysm,rhsm
      
      INTEGER wqkode(mesh%nq)
      INTEGER wkode(mesh%nnodes)
      INTEGER nunk, nb
            
      INTEGER ic,it,row,col,ii,jj,nuq,j,is,ir
      
      wqkode=mesh%mHqkode
      wkode=mesh%mHkode
      nunk=mesh%mHnunk
      nb=mesh%mHnb

      is=0
      ir=0
      DO ic=1,mesh%nicell ! zanka po celicah
        DO it=1,mesh%nsp ! po izvornih tockah znotraj celice
          row=(ic-1)*mesh%nsp+it ! vrstica v sistemski matriki
C
C         najprej H matrika (u)
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            val=smatH(row,col)    
            nuq=wkode(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              is=is+1
              sysm%v(sysm%d(is))=val
            ELSE ! znana vrednost - rhs
              ir=ir+1
              rhsm%v(rhsm%d(ir))=-val  ! minus zato, ker gre na drugo stran enacbe
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
                is=is+1
                sysm%v(sysm%d(is))=-smatG(row,col)
              ELSE ! znana vrednost - rhs
                ir=ir+1
                rhsm%v(rhsm%d(ir))=smatG(row,col) ! minus zato, ker gre na drugo stran enacbe
              END IF
            END DO
          END DO
        END DO
      END DO
      
      END      





C -----------------------------------------------------------------------------
      SUBROUTINE corMin_modHelm(env,io,inp,mesh,gauss,mHdiff,mHqDiff,smatH,smatG,smatB,
     &                             smatAbdx,smatAbdy,smatAbdz,
     &                             smatHtx,smatHty,smatHtz,
     &                             smatDx,smatDy,smatDz)
C
C     $: Compute or Read Matrices, modified Helmholtz eqaution
C        (H, G, Ab, Adx, Ady, Adz, B, Htx, Hty, Htz, Dx, Dy, Dz)
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      TYPE(gausstype) :: gauss
      TYPE(IOtype) io
      TYPE(inputtype) inp 
      TYPE(penv) :: env             
      
      REAL(8) mHdiff(mesh%nnodes),mHqDiff(mesh%nq)

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
      REAL(8) ontHtyDx,ontHtzDy,ontHtxDz,mumax,s

      INTEGER iok

C     sum of diff as parameter for integral integritiy
      s = 0.0D0
      DO i=1,mesh%nnodes
        s = s + mHdiff(i)
      END DO

     

c     Check if integrals exist on file
      CALL CheckmHMacroIntFile(io,inp%INTversion,mesh%nnodes,mesh%nbnodes,inp%beta,s,iok)

      IF (iok.EQ.1) THEN
c       Read integrals from disk 
        CALL WarnErr(env,io,inp,0,"corMin_modHelm","Reading modified Helmholtz integrals!",0)

        CALL ReadmHMacroIntDisk(io,mesh,smatH,smatG,smatB,
     &                             smatAbdx,smatAbdy,smatAbdz,
     &                             smatHtx,smatHty,smatHtz,
     &                             smatDx,smatDy,smatDz,
     &                             ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtyDx,ontHtzDy,ontHtxDz,mumax)

      ELSE
c       Calcualte integrals      
        CALL WarnErr(env,io,inp,0,"corMin_modHelm","Calculating modified Helmholtz integrals!",0)    

        CALL FMATmodHelm(env,io,inp,mesh,gauss,mHdiff,mHqDiff,smatH,smatG,smatB,
     &                          smatAbdx,smatAbdy,smatAbdz,
     &                          smatHtx,smatHty,smatHtz,
     &                          smatDx,smatDy,smatDz,
     &                          ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtyDx,ontHtzDy,ontHtxDz,mumax)


        IF (env%myproc.EQ.1) THEN ! to ni paralelizirano, samo en zapise
          CALL WritemHMacroIntDisk(io,inp,mesh,smatH,smatG,smatB,
     &                             smatAbdx,smatAbdy,smatAbdz,
     &                             smatHtx,smatHty,smatHtz,
     &                             smatDx,smatDy,smatDz,
     &                             ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtyDx,ontHtzDy,ontHtxDz,mumax,inp%beta,s)
        END IF
      END IF

      IF (env%myproc.EQ.1) THEN
        WRITE(io%l,'(A)') "Ocena natancnosti mod.Helmholtz integralov!"
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
        WRITE(io%l,'(A,G15.10)') "Max. mu value = ",mumax
        WRITE(io%l,'(A)') ""
      END IF      

      END


C______________________________________________________________________C
      SUBROUTINE CheckmHMacroIntFile(io,a,b,c,d,e,iok)
      USE inc_types
      TYPE(IOtype) io
      INTEGER a,b,c
      INTEGER aa,bb,cc,iok
      REAL(8) d,e,dd,ee

      iok=0
      OPEN (io%imh,FILE=TRIM(io%imh_name),FORM='UNFORMATTED',STATUS='OLD',ERR=10)
      READ(io%imh) aa
      READ(io%imh) bb,cc
      READ(io%imh) dd
      READ(io%imh) ee

      IF (a.EQ.aa.AND.b.EQ.bb.AND.c.EQ.cc.AND.d.EQ.dd.AND.e.EQ.ee) iok=1
      CLOSE(io%imh)

10    RETURN
      END

C -----------------------------------------------------------------------------
      SUBROUTINE WritemHMacroIntDisk(io,inp,mesh,smatH,smatG,smatB,
     &                             smatAbdx,smatAbdy,smatAbdz,
     &                             smatHtx,smatHty,smatHtz,
     &                             smatDx,smatDy,smatDz,
     &                             ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtyDx,ontHtzDy,ontHtxDz,mumax,b,s)
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
      REAL(8) ontHtyDx,ontHtzDy,ontHtxDz,mumax,b,s

      OPEN (io%imh,FILE=TRIM(io%imh_name),FORM='UNFORMATTED',STATUS='UNKNOWN')
      WRITE(io%imh) inp%INTversion
      WRITE(io%imh) mesh%nnodes,mesh%nbnodes
      WRITE(io%imh) b
      WRITE(io%imh) s
      WRITE(io%imh) mumax

      CALL WrSMat(smatH,mesh%nicnsp,mesh%npoc,io%imh)
      CALL WrSMat(smatG,mesh%nicnsp,mesh%npofc,io%imh)
      CALL WrSMat(smatB,mesh%nicnsp,mesh%npoc,io%imh)
      CALL WrSMat(smatAbdX,mesh%nicnsp,mesh%npoc,io%imh)
      CALL WrSMat(smatAbdY,mesh%nicnsp,mesh%npoc,io%imh)
      CALL WrSMat(smatAbdZ,mesh%nicnsp,mesh%npoc,io%imh)
      CALL WrSMat(smatHtx,mesh%nicnpoc,mesh%npoc,io%imh)
      CALL WrSMat(smatHty,mesh%nicnpoc,mesh%npoc,io%imh)
      CALL WrSMat(smatHtz,mesh%nicnpoc,mesh%npoc,io%imh)
      CALL WrSMat(smatDx,mesh%nicnsp,mesh%npoc,io%imh)
      CALL WrSMat(smatDy,mesh%nicnsp,mesh%npoc,io%imh)
      CALL WrSMat(smatDz,mesh%nicnsp,mesh%npoc,io%imh)
      WRITE(io%imh) ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtyDx,ontHtzDy,ontHtxDz


      CLOSE(io%imh)

      END



C -----------------------------------------------------------------------------
      SUBROUTINE ReadmHMacroIntDisk(io,mesh,smatH,smatG,smatB,
     &                             smatAbdx,smatAbdy,smatAbdz,
     &                             smatHtx,smatHty,smatHtz,
     &                             smatDx,smatDy,smatDz,
     &                             ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtyDx,ontHtzDy,ontHtxDz,mumax)
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
      REAL(8) ontHtyDx,ontHtzDy,ontHtxDz,mumax

      OPEN (io%imh,FILE=TRIM(io%imh_name),FORM='UNFORMATTED',STATUS='OLD')
      READ(io%imh) a
      READ(io%imh) b,c
      READ(io%imh) mumax
      READ(io%imh) mumax
      READ(io%imh) mumax

      CALL RdSMat(smatH,mesh%nicnsp,mesh%npoc,io%imh)
      CALL RdSMat(smatG,mesh%nicnsp,mesh%npofc,io%imh)
      CALL RdSMat(smatB,mesh%nicnsp,mesh%npoc,io%imh)
      CALL RdSMat(smatAbdX,mesh%nicnsp,mesh%npoc,io%imh)
      CALL RdSMat(smatAbdY,mesh%nicnsp,mesh%npoc,io%imh)
      CALL RdSMat(smatAbdZ,mesh%nicnsp,mesh%npoc,io%imh)
      CALL RdSMat(smatHtx,mesh%nicnpoc,mesh%npoc,io%imh)
      CALL RdSMat(smatHty,mesh%nicnpoc,mesh%npoc,io%imh)
      CALL RdSMat(smatHtz,mesh%nicnpoc,mesh%npoc,io%imh)
      CALL RdSMat(smatDx,mesh%nicnsp,mesh%npoc,io%imh)
      CALL RdSMat(smatDy,mesh%nicnsp,mesh%npoc,io%imh)
      CALL RdSMat(smatDz,mesh%nicnsp,mesh%npoc,io%imh)
      READ (io%imh) ontH,ontG,ontHtx,ontHty,ontHtz,ontDx,ontDy,ontDz,ontHtyDx,ontHtzDy,ontHtxDz

      CLOSE(io%imh)

      END


C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE INTEDc27modHelm
     &(
     &    gp,xp,yp,zp,
     &    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &    x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,
     &    be,dxe,dye,dze,isrc,nsipmx,mu
     &)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  INTEgration over Domain cell                                    **
c **  ----             -                                              **
c **  internal cell geometry : LINEAR ( 8-node )                      **
c **  internal cell function interpolation : QUADRATIC ( 27-node )    **
c **  fundamental solution   : MODIFIED HELMHOLTZ                     **
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

      REAL(8) a,b,c,d,e,f,dex,gii,gij,gik,mu
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
C***              MODIFIED HELMHOLTZ FUNDAMENTAL SOLUTION
C
                  FUNG=0.25D0/(PI*RA)*EXP(-mu*RA)
                  FUNX=(XP-XC0)*IRA2*FUNG*(1+mu*RA)
                  FUNY=(YP-YC0)*IRA2*FUNG*(1+mu*RA)
                  FUNZ=(ZP-ZC0)*IRA2*FUNG*(1+mu*RA)
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
C***              MODIFIED HELMHOLTZ FUNDAMENTAL SOLUTION
C
                  FUNG=0.25D0/(PI*RA)*EXP(-mu*RA)
                  FUNX=(XP-XC0)*IRA2*FUNG*(1+mu*RA)
                  FUNY=(YP-YC0)*IRA2*FUNG*(1+mu*RA)
                  FUNZ=(ZP-ZC0)*IRA2*FUNG*(1+mu*RA)
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




            
C -----------------------------------------------------------------------------            
      SUBROUTINE INTEBc9modHelm(gp,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                   isrc,nsipex,nsipexg,he,ge,adx,ady,adz,
     &                   heyz,hezx,hexy,minedge,mu)
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
      REAL(8) fig(4),fih(9),fig4(4)      
      REAL(8) al(4),fii(4),th0(4),th1(4),ksi(13),eta(13)
c
c      integral divison
      REAL(8) a,b,c,d,dex,gii,gij,d1,d2,d3,d4,d13max,d24max,minedge,mu
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
C***          MODIFIED HELMHOLTZ FUNDAMENTAL SOLUTION
C
              FUNG=0.25D0/(PI*RA)*EXP(-mu*RA)
              FUNX=(XP-XC0)*IRA2*FUNG*(1+mu*RA)
              FUNY=(YP-YC0)*IRA2*FUNG*(1+mu*RA)
              FUNZ=(ZP-ZC0)*IRA2*FUNG*(1+mu*RA)
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
C***        MODIFIED HELMHOLTZ FUNDAMENTAL SOLUTION
C
            FUNG=0.25D0/(PI*RA)*EXP(-mu*RA)
            FUNX=(XP-XC0)*IRA2*FUNG*(1+mu*RA)
            FUNY=(YP-YC0)*IRA2*FUNG*(1+mu*RA)
            FUNZ=(ZP-ZC0)*IRA2*FUNG*(1+mu*RA)
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


