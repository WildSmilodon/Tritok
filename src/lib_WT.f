C     ******************************************************************
C     *                                                                *
C     *                          lib_WT.f                              *
C     *                                                                *
C     ******************************************************************
C     *                                                                *
C     *                     Author: Jure Ravnik                        *
C     *                                                                *
C     ******************************************************************
C     *                                                                *
C     *                     Version:  15.2.2018                        *
C     *                                                                *
C     ******************************************************************


C     ******************************************************************
C     *                                                                *
C     *                    DAUBECHIES WAVELETS                         *
C     *                                                                *
C     * Rutina, ki postavi koeficiente h in g za izvajanje             *
C     * valcne transformacije with Daubechies wavelets                 *
C     * Reference:                                                     *
C     * Daubechies, I. "Ten Lectures on Wavelets"                      *
C     * Society for industrial and applied mathematics, 1992           *
C     *                                                                *
C     ******************************************************************

C       Usage:

C     --> FiltCoef(h,g,M)
c      INTEGER M                ! number of vanishing moments
c      REAL*8 h(2*M),g(2*M)     ! wavelet filter coefficients

C     ******************************************************************
C     *                                                                *
C     *                            CHE                                 *
C     *                                                                *
C     * Rutine za izvajanje valcne transformacije z Haar valcki        *
C     * na matrikah in vektorjih poljubnih dolzin, pri cemer           *
C     * uporabimo raztegovanje, da ravno pravo stevilo valcnih         *
C     * koeficientov postane nic                                       *
C     * Reference:                                                     *
C     * Ravnik, J.; Skerget, L. & Hribersek, M.                        *
C     * The wavelet transform for BEM computational fluid dynamics     *
C     * Eng. Anal. Bound. Elem., 2004, 28, 1303-1314                   *
C     *                                                                *
C     ******************************************************************

c       Uporaba :

c     --> Che_WcomMatrix(a,nVR,nST,crsMat,kappa,MB,comp,spar,norm)
c     --> Che_WTAWx(crsMat,vekX,ncol,vekB,nrow)  ! B = W^T ( crsMat * (W * x) )
c     --> Che_WaWT(a,nVR,nST)
c     --> Che_FormCRS(FullMat,nrow,ncol,crsMat)
c     --> Che_CRSmatMul(crsMat,vekX,ncol,vekB,nrow)
c     --> Che_WTxV(v,nvek)
c     --> Che_WxV(v,nvek)
c     --> Che_Zanemari(mat,nVR,nST,tresh,nz,nnz)        ! fiksna meja za vse razen scaling
c     --> Che_VTzanemari(mat,nrow,ncol,meja,fak,nz,nnz) ! spremenljiva meja fak**(J-(j+j')/2), razen scaling
c     --> Che_AbsPovp(mat,nVR,nST,povp)
c     --> Che_CRSdataSize(crsMat,nrow,ncol,MB,comp,spar)

c      Potrebujemo datoteko 'inc_types.f' z definiranim tipom
!      TYPE crsMatrixType
!        INTEGER nrow        ! stevilo vrstic
!        INTEGER nz          ! stevilo nicel v polni matriki
!        INTEGER nnz         ! stevilo nenicelnih elementov
!        INTEGER, POINTER :: zvr(:),c(:) ! zvr(nrow+1), c(nnz)
!        REAL(8), POINTER :: v(:)  ! v(nnz) vrednosti
!      END TYPE    

C     ******************************************************************
C     *                                                                *
C     *                          KOLAZ                                 *
C     *                                                                *
C     * Rutine za izvajanje valcne transformacije z Daubechie valcki   *
C     * na matrikah in vektorjih poljubnih dolzin, pri cemer           *
C     * matriko razrezemo v kolaz delov velikosti 2^n x 2^m            *
C     * Reference:                                                     *
C     * - the KOLAZ is still unpublished                               *
C     * - the 2^n fast wavelet transform using Daubechies wavelets     *
C     *   Beylkin, G.; Coifman, R. & Rokhlin, V.                       *
C     *   Fast wavelet transforms and numerical algorithms             *
C     *   Comm. Pure Appl. Math, 1991, 44, 141-183                     *
C     *                                                                *
C     ******************************************************************

c       Uporaba:
c
c       Spremenljivke

!      INTEGER     m    ! number of vanishing moments
!      REAL(8) kappa    ! meja = kappa * meanabs
!      REAL(8)   h,g    ! valcni koeficienti, polja dolzine (2*M)


c       Rutine, s katerimo naredimo kolaz strukturo:

c     -->  Kol_DivVec(n,nVecDiv,VecDiv)     
c     -->  Kol_CountDivVec(n,st)
c     -->  Kol_DivMat(nrow,ncol,fulr,fulc,ns,KolMat,st)
c     -->  Kol_CountDivMat(nrow,ncol,fulr,fulc,ns)                     

c       Izpis strukture na disk v tecplot formatu (optional)

c     -->  Kol_MatTecOut(nKolMat,KolMat)

c       Valcna tranfromacija matrike, zanemarjanje in zapis v CRS

c     -->  Kol_WaWT(FullMat,nrow,ncol,KolMat,nKolMat,M,h,g)   
c     -->  Kol_Set2Zero(FullMat,nrow,ncol,KolMat,nKolMat,kappa)       ! fiksna meja za vse razen scaling
c     -->  Kol_vtSet2Zero(FullMat,nrow,ncol,KolMat,nKolMat,kappa,fak) ! spremenljiva meja fak**(J-(j+j')/2), razen scaling
c     -->  Kol_FormCRS(FullMat,nrow,ncol,KolMat,nKolMat)

c       Izracun velikosti podatkovne strukture, kompresije, sparsity (optional)

c     -->  Kol_CRSdataSize(KolMat,nKolMat,nrow,ncol,MB,comp,spar)

c       Valcna in obratna valcna transformacija vektorja

c     -->  Kol_FWTxV(vec,length,KolVec,nKolVec,M,h,g)
c     -->  Kol_iFWTxV(vec,length,KolVec,nKolVec,M,h,g)      

c       Mnozenje kolaz matrike z vektorjem

c     -->  Kol_CRSmatMul(KolMat,nKolMat,vekX,ncol,vekB,nrow)

c      Potrebujemo datoteko 'inc_types.f' z definiranima tipoma
!      TYPE MatrixKolazType
!        INTEGER ulr,ulc
!        INTEGER nrow,ncol
!        REAL(8) meanabs     ! povprecna vrednost absolutnih vrednosti vseh elementov v matriki
!        REAL(8) meja        ! meja za zanemarjanje
!        INTEGER nz          ! stevilo nicel v polni matriki
!        
!        INTEGER :: nnz      ! stevilo nenicelnih elementov
!        INTEGER, POINTER :: zvr(:),c(:) ! zvr(nrow+1), c(nnz)
!        REAL(8), POINTER :: v(:)  ! v(nnz) vrednosti
!      END TYPE
!      
!      TYPE VectorKolazType
!        INTEGER length
!        INTEGER start
!      END TYPE


C     ------------------------------------------------------------------   
      SUBROUTINE Che_WcomMatrix(mat,nVR,nST,crsMat,kappa,MB,comp,spar,norm)
C
C     glede na kappa naredi WAWT -> CRS obliko
C
C     ------------------------------------------------------------------
      USE inc_types

      TYPE(crsMatrixType) :: crsMat
      INTEGER nVR,nST,nz,nnz
      REAL(8) kappa,povp,tresh
      REAL(8) mat(nVR,nST)
      REAL(8) MB,comp,spar
      REAL(8), ALLOCATABLE :: vec(:),rez1(:),rez2(:)
      REAL(8) s,s2,norm
      INTEGER i

      ALLOCATE ( vec(nSt),rez1(nVR),rez2(nVR) )
c     set up random vector for accuracy test
      DO i=1,nST
        CALL RANDOM_NUMBER(vec(i))
      END DO
c     full matrix x vec results
      rez1=MATMUL(mat,vec)

C     do WAWT (orig matrix destroyed !!)
      CALL Che_WaWT(mat,nVR,nST)

c     get average element size
      CALL Che_AbsPovp(mat,nVR,nST,povp)

c     zero out elements below threshold
      tresh = kappa * povp
      CALL Che_Zanemari(mat,nVR,nST,tresh,nz,nnz)  ! fiksna meja za vse razen scaling

c     create CRS representation of the matrix
      crsMat%nnz = nnz
      CALL Che_FormCRS(mat,nVR,nST,crsMat)

c     estimate CRS data size
      CALL Che_CRSdataSize(crsMat,nVR,nST,MB,comp,spar)

c     compressed CRS matrix x vec result
c     wavelet transform of vector
      CALL Che_WxV(vec,nSt)
c     matrix x vector
      CALL Che_CRSmatMul(crsMat,vec,nSt,rez2,nVR)
c     inverse wavelt transform of rez2
      CALL Che_WTxV(rez2,nVR)
c     get norm
      s=0.0D0
      s2=0.0D0
      DO i=1,nVR
        s=s+(rez1(i)-rez2(i))**2
        s2=s2+rez1(i)**2
      ENDDO
      norm = SQRT(s/s2)

      DEALLOCATE (vec,rez1,rez2)

      END


C     ------------------------------------------------------------------
      SUBROUTINE Che_VTzanemari(mat,nrow,ncol,meja,fak,nz,nnz)
C
C     Zanemarjam elemente, pri cemer mejo dvigujem s faktorjem
C     fak**(J-(j+j')/2)
C      
C     ------------------------------------------------------------------      
      INTEGER nrow,ncol,nz,nnz
      REAL(8) mat(nrow,ncol),meja,fak,tresh
      
      INTEGER kmaxR,kmaxC,maxk,eR,eC,sR,sC


      INTEGER i,j,ii,jj
      
      CALL DvaNaK(nrow,kmaxR)
      CALL DvaNaK(ncol,kmaxC)
      
      maxk=MAX(kmaxR,kmaxC)+1

c     prvi clen je scaling wavelet, ga ne bomo zanemarili
      nnz=nrow+ncol-1
      nz=0
      DO i=0,kmaxR
        IF (i.EQ.kmaxR) THEN ! zaradi moznosti, da matrika nima 2^(kmaxR+1) vrstic
          eR=nrow
          sR=2**i+1
        ELSE IF (i.EQ.0) THEN
          eR=2
          sR=2        
        ELSE
          sR=2**i+1 
          eR=2**(i+1)
        END IF
        DO j=0,kmaxC
          IF (j.EQ.kmaxC) THEN ! zaradi moznosti, da matrika nima 2^(kmaxC+1) stolpcev
            eC=ncol
            sC=2**j+1
          ELSE IF (j.EQ.0) THEN
            eC=2
            sC=2
          ELSE
            eC=2**(j+1)
            sC=2**j+1
          END IF
          tresh=meja*fak**(maxk-(i+j)/2.0D0)
          DO ii=sR,eR
            DO jj=sC,eC
              IF (ABS(mat(ii,jj)).LT.tresh) THEN
                mat(ii,jj)=0.0D0
                nz=nz+1
              ELSE
                nnz=nnz+1
              END IF
            END DO
          END DO
        END DO
      END DO

          
      END

C     ------------------------------------------------------------------      
      SUBROUTINE Che_calLastMp(mat,nrow,ncol,mp)
C
C     Izracuna povprecno vrednost absolutnih vrednosti
C     samo za del matrike, ki pripada valckom najvecjega reda
C      
C     ------------------------------------------------------------------      
      INTEGER nrow,ncol
      REAL(8) mat(nrow,ncol),mp
      
      INTEGER kmaxR,kmaxC,j,i
      
      CALL DvaNaK(nrow,kmaxR)
      CALL DvaNaK(ncol,kmaxC)

      mp=0.0D0
      DO i=2**kmaxR+1,nrow
        DO j=2**kmaxC+1,ncol
          mp=mp+ABS(mat(i,j))
        END DO
      END DO
      mp=mp/DBLE(nrow-2**kmaxR)/DBLE(ncol-2**kmaxC)
     
      END
C     ------------------------------------------------------------------
      SUBROUTINE Che_WTAWx(crsMat,X,ncol,vekB,nrow)
C
C     Mnozenje matrike z vektorjem, CRS verzija, vkljucuje W in WT
C     Rezultati vekB
C
C     ------------------------------------------------------------------
      USE inc_types

      TYPE(crsMatrixType) :: crsMat

      INTEGER nrow,ncol
      REAL(8) X(ncol),vekB(nrow)
      REAL(8), ALLOCATABLE :: vekX(:)

c     make a copy
      ALLOCATE (vekX(ncol))
      DO i=1,ncol
        vekX(i)=X(i)
      END DO

c     W*x
      CALL Che_WxV(vekX,ncol)
c     A*(W*x)
      CALL Che_CRSmatMul(crsMat,vekX,ncol,vekB,nrow)
C     W^T(A*(W*x))
      CALL Che_WTxV(vekB,nrow)

      DEALLOCATE (vekX)

      END


C     ------------------------------------------------------------------      
      SUBROUTINE Che_CRSmatMul(crsMat,vekX,ncol,vekB,nrow)
C
C     Mnozenje matrike z vektorjem, CRS verzija  
C      
C     ------------------------------------------------------------------
      USE inc_types 
      
      TYPE(crsMatrixType) :: crsMat
      
      INTEGER nrow,ncol
      REAL(8) vekX(ncol),vekB(nrow)
      
      INTEGER r,c
      

      vekB=0.0D0      
      DO r=1,crsMat%nrow
        DO c=crsMat%Zvr(r),crsMat%Zvr(r+1)-1
          vekB(r)=vekB(r)+crsMat%v(c)*vekX(crsMat%c(c))
        END DO
      END DO

      END

C     ------------------------------------------------------------------
      SUBROUTINE Che_FormCRS(FullMat,nrow,ncol,crsMat)
C
C     Rutina za izdelavo CRS zapisa matrik v kolazu  
C     predpostvim, da je v vsaki vrstici FullMat, vsaj en clen razlicen od nic
C      
C     ------------------------------------------------------------------
      USE inc_types 
      
      TYPE(crsMatrixType) :: crsMat
      INTEGER nrow,ncol
      REAL(8) FullMat(nrow,ncol)

      
      INTEGER i,j,st,stvr,ibu


      crsMat%nrow=nrow

      ALLOCATE (crsMat%v(crsMat%nnz))
      ALLOCATE (crsMat%c(crsMat%nnz))
      ALLOCATE (crsMat%zvr(crsMat%nrow+1))

      st=0
      stvr=0
      ibu=0

      crsMat%zvr=0
    
      DO i=1,nrow
        DO j=1,ncol
          IF (FullMat(i,j).NE.0.0D00) THEN
            st=st+1
            crsMat%v(st)=FullMat(i,j)
            crsMat%c(st)=j
            IF (i.ne.ibu) THEN
              stvr=stvr+1
              ibu=i
              crsMat%zvr(stvr)=st
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      crsMat%zvr(stvr+1)=st+1
      
      END

C     ------------------------------------------------------------------
      SUBROUTINE Che_CRSdataSize(crsMat,nrow,ncol,MB,comp,spar)
C
C     Izracuna stevilo podatkov v megabytih, Ki jo zavzema crsMat  
C     comp = datasizeCRS / datasizeFULL
C     spar = stevilo nenicelnih clenov / stevilo vseh elementov matrike
C      
C     ------------------------------------------------------------------
      USE inc_types 
      
      TYPE(crsMatrixType) :: crsMat

      INTEGER ni,nr,nrow,ncol
      REAL(8) MB,comp,spar

      nr=crsMat%nnz  ! stevilo real(8)
      ni=crsMat%nnz+crsMat%nrow+4 ! stevilo integerjev (4)
      
      MB  =(DBLE(ni)*4.0D0+DBLE(nr)*8.0D0)/1024.0D0/1024.0D0
      comp=(DBLE(ni)*0.5D0+DBLE(nr))/DBLE(nrow)/DBLE(ncol)
      spar=DBLE(nr)/DBLE(nrow)/DBLE(ncol)
      
      END



C     ------------------------------------------------------------------
      SUBROUTINE Che_WaWT(a,nVR,nST)
C
C     Naredi WaWt na matriki a
C      
C     ------------------------------------------------------------------
      INTEGER nVR,nST,ii,jj,kk,i
      INTEGER n,k,l,kmax,DvaNaKmaxmk,PodSP,PodZG,PodSR
      INTEGER podv
      REAL*8 a(nVR,nST)
      REAL*8 v,valcek,norma
      ALLOCATABLE v(:)

c   ***  Delamo SWRxA
c   kmax=int(log(real(nVR))/log(2.0)-1)+1
      CALL DvaNaK(nVR,kmax)
      n=2**(kmax+1)
      ALLOCATE (v(n))
      norma=1.0/Sqrt(Dfloat(n))

      DO kk=1,nST


      if (n.NE.nVR) then
c       podvajam
        podv=2*(nVR-2**kmax)+1

        DO ii=1,podv-1
        v(ii)=a(ii,kk)
        ENDDO

        ii=podv-1
        jj=ii
10      ii=ii+1
        jj=jj+1
        v(ii)=a(jj,kk)
        ii=ii+1
        v(ii)=a(jj,kk)
        if (ii.LT.n) goto 10
      else
        DO ii=1,n
        v(ii)=a(ii,kk)
        ENDDO
      endif

      DO ii=1,nVR
      a(ii,kk)=0  
c     Prva vsrtica je polna - scaling wavelet
      IF (ii.eq.1) THEN
        do i=1,n
        a(ii,kk)=a(ii,kk)+v(i)*norma
        enddo
      ELSE

c   Nato koeficiente pri ostalih valckih (=vrsticah) 
c   k=int(log(real(ii))/log(2.0)-1)+1
c   if (ii.EQ.2) k=0
      CALL DvaNaK(ii,k)
      l=ii-2**k
      DvaNaKmaxmk=2**(kmax-k)
      PodSP=(2*l-1-1)*DvaNaKmaxmk+1
      PodSR=(2*l-1)*DvaNaKmaxmk
      PodZG=(2*l)*DvaNaKmaxmk
      valcek=(2.0**(Dfloat(k)/2))*norma

    
            do i=PodSP,PodSR
                a(ii,kk)=a(ii,kk)+valcek*v(i)
            enddo
            
            do i=PodSR+1,PodZG
                a(ii,kk)=a(ii,kk)-valcek*v(i)
            enddo
    

      ENDIF
      ENDDO   

      ENDDO
      DEALLOCATE (v)
c   ***  Delamo AxSWR^-1
c   kmax=int(log(real(nST))/log(2.0)-1)+1
      CALL DvaNaK(nST,kmax)
      n=2**(kmax+1)
      ALLOCATE (v(n))
      norma=1.0/Sqrt(Dfloat(n))


      DO kk=1,nVR


c   najprej vt*R^-1 = raztegnem vrstivco s podvajanjem polovicnih clenov
      if (n.NE.nST) then
        podv=2*(nST-2**kmax)+1

        DO ii=1,podv-1
        v(ii)=a(kk,ii)
        ENDDO

        ii=podv-1
        jj=ii
100     ii=ii+1
        jj=jj+1
        v(ii)=a(kk,jj)
        ii=ii+1
        v(ii)=0.0D00
        if (ii.LT.n) goto 100
      else
        DO ii=1,n
        v(ii)=a(kk,ii)
        ENDDO
      endif

c nato podaljsana vrstica na Wt = to je enako kot W*V
      DO ii=1,nST
      a(kk,ii)=0  
c     Prva vsrtica je polna - scaling wavelet
      IF (ii.eq.1) THEN
        do i=1,n
        a(kk,ii)=a(kk,ii)+v(i)*norma
        enddo
      ELSE

c   Nato koeficiente pri ostalih valckih (=vrsticah) 
c   k=int(log(real(ii))/log(2.0)-1)+1
c   if (ii.EQ.2) k=0
      CALL DvaNaK(ii,k)
      l=ii-2**k
      DvaNaKmaxmk=2**(kmax-k)
      PodSP=(2*l-1-1)*DvaNaKmaxmk+1
      PodSR=(2*l-1)*DvaNaKmaxmk
      PodZG=(2*l)*DvaNaKmaxmk
      valcek=(2.0**(Dfloat(k)/2))*norma
    
            do i=PodSP,PodSR
                a(kk,ii)=a(kk,ii)+valcek*v(i)
            enddo
            
            do i=PodSR+1,PodZG
                a(kk,ii)=a(kk,ii)-valcek*v(i)
            enddo
    

      ENDIF
      ENDDO   


      ENDDO
      DEALLOCATE (v)

      END

C     ------------------------------------------------------------------
      SUBROUTINE Che_WxV(r,nvek)
C
c     Naredimo W * V, W je valcna matrika, ki jo naredim sproti, V je vektor
c     preverim, ce je stevilo elementov vektorja 2^n
c     ce imam elementov manj, je strategija taka, da spremenim
c     originalni vektor r tako, da bo transformiran imel same nicle
c     clene zacnem podvajati pri 2*(nvek-2^kmax)+1
c     To je ze SWR verzija, izpisano tisto kar naredi rutina mSWRxV
c      
C     ------------------------------------------------------------------    
      INTEGER n,i,ii,jj,k,l,kmax,DvaNaKmaxmk,PodSP,PodZG,PodSR
      INTEGER nvek,podv
      REAL*8 v,r(nvek),valcek,norma
      ALLOCATABLE v(:)

      CALL DvaNaK(nvek,kmax)
      n=2**(kmax+1)
      ALLOCATE (v(n))
      norma=1.0/Sqrt(Dfloat(n))   


      if (n.NE.nvek) then
c       podvajam
        podv=2*(nvek-2**kmax)+1

        DO ii=1,podv-1
        v(ii)=r(ii)
        ENDDO

        ii=podv-1
        jj=ii
10      ii=ii+1
        jj=jj+1
        v(ii)=r(jj)
        ii=ii+1
        v(ii)=r(jj)
        if (ii.LT.n) goto 10
      else
        v=r
      endif

      DO ii=1,nvek
      r(ii)=0 
c     Prva vsrtica je polna - scaling wavelet
      IF (ii.eq.1) THEN
        do i=1,n
        r(ii)=r(ii)+v(i)*norma
        enddo
      ELSE

c   Nato koeficiente pri ostalih valckih (=vrsticah) 
c   k=int(log(real(ii))/log(2.0)-1)+1
c   if (ii.EQ.2) k=0
      CALL DvaNaK(ii,k)
      l=ii-2**k
      DvaNaKmaxmk=2**(kmax-k)
      PodSP=(2*l-1-1)*DvaNaKmaxmk+1
      PodSR=(2*l-1)*DvaNaKmaxmk
      PodZG=(2*l)*DvaNaKmaxmk
      valcek=(2.0**(Dfloat(k)/2))*norma

    
            do i=PodSP,PodSR
                r(ii)=r(ii)+valcek*v(i)
            enddo
            
            do i=PodSR+1,PodZG
                r(ii)=r(ii)-valcek*v(i)
            enddo
    

      ENDIF
      ENDDO   

      END



C     ------------------------------------------------------------------
      SUBROUTINE Che_WTxV(v,nvek)
c      
c     Naredimo W^T * V, W^T je transponirana valcna matrika, ki jo naredim sproti, V je vektor
c     predpostavimo, da so vsi koeficienti nad nvek enaki nic (to smo dosegli z prvo WT)
c     resitev vsebuje podvovjene clene, ki jih ob koncu vrzemo ven
c     To je izpisna obratna SWR transformacija    
c      
C     ------------------------------------------------------------------      
      INTEGER n,ii,jj,k,l,kmax,DvaNaKmaxmk
      INTEGER nvek, podv, iks
      REAL*8 v(nvek),r,valcek,norma
      ALLOCATABLE r(:)



      CALL DvaNaK(nvek,kmax)
      n=2**(kmax+1)
      ALLOCATE (r(n))
      norma=1.0/Sqrt(Dfloat(n))


      DO ii=1,n
c Scaling wavelet = norma
      r(ii)=v(1)*norma
        DO k=0,kmax
            valcek=(2.0**(Dfloat(k)/2))*norma
            DvaNaKmaxmk=2**(kmax-k)
            iks=1+int((real(ii)-0.01)/real(DvaNaKmaxmk))
    
            if (mod(iks,2).eq.0) then
                l=iks/2
                jj=2**k+l
                if (jj.LE.nvek) r(ii)=r(ii)-valcek*v(jj)
            else
                l=(iks+1)/2
                jj=2**k+l
                if (jj.LE.nvek) r(ii)=r(ii)+valcek*v(jj)
            endif
        ENDDO
      ENDDO

      podv=2*(nvek-2**kmax)+1

      if (n.NE.nvek) then
c       Prepisemo nazaj, ker so v r podvojene
        DO ii=1,podv-1
        v(ii)=r(ii)
        ENDDO
        ii=podv-1
        jj=ii
10      ii=ii+2
        jj=jj+1
        v(jj)=r(ii)
        if (ii.LT.n) goto 10
      else
        v=r
      endif
      END

C     ------------------------------------------------------------------    
      SUBROUTINE Che_AbsPovp(mat,nVR,nST,povp)
c
c     Izracuna povprecno vrednost aboslutnih vrednosti vseh elementov      
c      
C     ------------------------------------------------------------------          
      INTEGER i,j,nVR,nST
      REAL*8 mat(nVR,nST),povp

      povp=0.0D0

      DO i=1,nST
        DO  j=1,nVR
            povp=povp+abs(mat(j,i))
        ENDDO
      ENDDO

      povp=povp/real(nVR)/real(nST)

      END


C     ------------------------------------------------------------------     
      SUBROUTINE Che_Zanemari(mat,nVR,nST,tresh,nz,nnz)
c
c     Postavi elementve na nic - poskrbeti moram, da vsaj en v vrsti ni nic      
c      
C     ------------------------------------------------------------------     
      INTEGER i,j,nVR,nST,nz,nnz
      REAL*8 mat(nVR,nST),tresh
      LOGICAL VsiNic
    
      nz=0
      nnz=0
      DO  i=2,nVR ! ne zanemarjam scaling (in mother) waveleta !!!!
        VsiNic=.TRUE.
        DO j=2,nST-1  ! ne zanemarjam scaling in mother waveleta !!!!
          IF (abs(mat(i,j)).LT.tresh) THEN
            mat(i,j)=0.0D00
            nz=nz+1
          ELSE
            VsiNic=.FALSE.
          END IF
        ENDDO
        j=nST
        IF (.NOT.(VsiNic)) THEN ! ga lahko zanemarim
          IF (abs(mat(i,j)).LT.tresh) THEN
            mat(i,j)=0.0D00
            nz=nz+1
          END IF
        ELSE   ! ne sme biti nic
          IF (mat(i,j).EQ.0.0D0) THEN
            mat(i,j)=1.0D-14
          END IF
        END IF
      ENDDO
      nnz=nVR*nST-nz
      END

C     ------------------------------------------------------------------     
      SUBROUTINE DvaNaK(n,k)
c      
c     ideja te rutine je v tem, da za 33 do 64 poda k=5, 65 - 128 - k=6, itd.
c
C     ------------------------------------------------------------------     
      INTEGER n,k,i


      IF (n.EQ.1) THEN
        k=0
        RETURN
      ENDIF

      DO i=0,100
        IF (2**i.GE.n) THEN
            k=i-1
            RETURN
        ENDIF
      ENDDO

      END

C     ------------------------------------------------------------------     
      SUBROUTINE DvaNaKp(n,k)
c      
c     2**k je prva potenca stevila 2 manjsa ali enaka n
c
C     ------------------------------------------------------------------     
      INTEGER n,k,i


      IF (n.EQ.1) THEN
        k=0
        RETURN
      ENDIF

      DO i=0,100
        IF (2**i.GT.n) THEN
            k=i-1
            RETURN
        ENDIF
      ENDDO

      END       

C     ------------------------------------------------------------------
      SUBROUTINE Kol_CRSdataSize(KolMat,nKolMat,nrow,ncol,MB,comp,spar)
C
C     Izracuna stevilo podatkov v megabytih, Ki jo zavzema KolMat  
C     comp = datasizeCRS / datasizeFULL
C     spar = stevilo nenicelnih clenov / stevilo vseh elementov matrike
C      
C     ------------------------------------------------------------------
      USE inc_types 
      
      INTEGER nKolMat
      TYPE(MatrixKolazType) :: KolMat(nKolMat)

      INTEGER ni,nr,i,nrow,ncol
      REAL(8) MB,comp,spar

      ni=0  ! stevilo integerjev (4)
      nr=0  ! stevilo real(8)

      DO i=1,nKolMat
        nr=nr+KolMat(i)%nnz+2
        ni=ni+KolMat(i)%nnz+KolMat(i)%nrow+7
      END DO
      
      MB  =(DBLE(ni)*4.0D0+DBLE(nr)*8.0D0)/1024.0D0/1024.0D0
      comp=(DBLE(ni)*0.5D0+DBLE(nr))/DBLE(nrow)/DBLE(ncol)
      spar=DBLE(nr-2*nKolMat)/DBLE(nrow)/DBLE(ncol)
      
      END


C     ------------------------------------------------------------------
      SUBROUTINE Kol_FormCRS(FullMat,nrow,ncol,KolMat,nKolMat)
C
C     Rutina za izdelavo CRS zapisa matrik v kolazu  
C     predpostvim, da je v vsaki vrstici FullMat, vsaj en clen razlicen od nic
C      
C     ------------------------------------------------------------------
      USE inc_types 
      
      INTEGER nKolMat
      TYPE(MatrixKolazType) :: KolMat(nKolMat)
      INTEGER nrow,ncol
      REAL(8) FullMat(nrow,ncol)

      
      INTEGER i,j,st,stvr,ibu,ii

      DO ii=1,nKolMat  ! zanka po kolazu matrik

        KolMat(ii)%nnz=Kolmat(ii)%nrow * Kolmat(ii)%ncol - Kolmat(ii)%nz 
        ALLOCATE (KolMat(ii)%v(KolMat(ii)%nnz))
        ALLOCATE (KolMat(ii)%c(KolMat(ii)%nnz))
        ALLOCATE (KolMat(ii)%zvr(KolMat(ii)%nrow+1))

        st=0
        stvr=0
        ibu=0

        KolMat(ii)%zvr=0
    
        DO i=1,KolMat(ii)%nrow
          DO j=1,KolMat(ii)%ncol
            IF (FullMat(KolMat(ii)%ulr-1+i,KolMat(ii)%ulc-1+j).NE.0.0D00) THEN
              st=st+1
              KolMat(ii)%v(st)=FullMat(KolMat(ii)%ulr-1+i,KolMat(ii)%ulc-1+j)
              KolMat(ii)%c(st)=j
              IF (i.ne.ibu) THEN
                stvr=stvr+1
                ibu=i
                KolMat(ii)%zvr(stvr)=st
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        KolMat(ii)%zvr(stvr+1)=st+1
      END DO
      
      END
      
C     ------------------------------------------------------------------      
      SUBROUTINE Kol_CRSmatMul(KolMat,nKolMat,vekX,ncol,vekB,nrow)
C
C     Mnozenje kolaz matrike z vektorjem, CRS verzija  
C      
C     ------------------------------------------------------------------
      USE inc_types 
      
      INTEGER nKolMat
      TYPE(MatrixKolazType) :: KolMat(nKolMat)
      
      INTEGER nrow,ncol
      REAL(8) vekX(ncol),vekB(nrow)
      
      INTEGER i,r,c,rr,cc
      

      vekB=0.0D0      
      DO i=1,nKolMat
        DO r=1,KolMat(i)%nrow
          rr=r+KolMat(i)%ulr-1
          DO c=KolMat(i)%Zvr(r),KolMat(i)%Zvr(r+1)-1
            cc=KolMat(i)%c(c)+KolMat(i)%ulc-1
            vekB(rr)=vekB(rr)+KolMat(i)%v(c)*vekX(cc)
          END DO
        END DO
      END DO


      END

C     ------------------------------------------------------------------
      SUBROUTINE Kol_vtSet2Zero(FullMat,nrow,ncol,KolMat,nKolMat,kappa,fak)
C
c     Naredi in-situ zanemarjanje clenov v WAWT matriki, variable threshold
C     C     fak**(J-(j+j')/2)
C      
C     ------------------------------------------------------------------
      USE inc_types 
      
      INTEGER nKolMat
      TYPE(MatrixKolazType) :: KolMat(nKolMat)
      INTEGER nrow,ncol
      REAL(8) FullMat(nrow,ncol),fak
      
      INTEGER maxk,kmaxR,kmaxC
      
      INTEGER i,j,r,c,eR,eC,sR,sC,ij
      REAL(8) kappa,tresh   

      DO ij=1,nKolMat  ! zanka po delih matrike
        KolMat(ij)%nz=0
        KolMat(ij)%meja = kappa * KolMat(ij)%meanabs
        
             
        CALL DvaNaK(KolMat(ij)%nrow,kmaxR)
        CALL DvaNaK(KolMat(ij)%ncol,kmaxC)

        maxk=MAX(kmaxR,kmaxC)+1

c     prvi clen je scaling wavelet, ga ne bomo zanemarili
c        nnz=nrow+ncol-1
        KolMat(ij)%nz=0
        DO i=0,kmaxR
          IF (i.EQ.0) THEN
            eR=KolMat(ij)%ulr+1
            sR=KolMat(ij)%ulr+1
            IF (KolMat(ij)%nrow.EQ.1) eR=sR-1 ! ce je samo 1 debelo, ni zanemarjanja
          ELSE
            sR=KolMat(ij)%ulr-1+2**i+1 
            eR=KolMat(ij)%ulr-1+2**(i+1)
          END IF
          DO j=0,kmaxC
            IF (j.EQ.0) THEN
              sC=KolMat(ij)%ulc+1
              eC=KolMat(ij)%ulc+1
              IF (KolMat(ij)%ncol.EQ.1) eC=sC-1
            ELSE
              eC=KolMat(ij)%ulc-1+2**(j+1)
              sC=KolMat(ij)%ulc-1+2**j+1
            END IF
            tresh=KolMat(ij)%meja*fak**(maxk-(i+j)/2.0D0)
            DO r=sR,eR
              DO c=sC,eC             
                IF (ABS(FullMat(r,c)).LT.tresh.OR.FullMat(r,c).EQ.0.0D0) THEN ! zanemari
                  FullMat(r,c)=0.0D0
                  KolMat(ij)%nz=KolMat(ij)%nz+1
                END IF 
              END DO     
            END DO
          END DO
        END DO
      END DO
      
      END


C     ------------------------------------------------------------------
      SUBROUTINE Kol_Set2Zero(FullMat,nrow,ncol,KolMat,nKolMat,kappa)
C
c     Naredi in-situ zanemarjanje clenov v WAWT matriki
c     scaling wavelet ne zanemarjam, za ostale fixed threshold
C      
C     ------------------------------------------------------------------
      USE inc_types 
      
      INTEGER nKolMat
      TYPE(MatrixKolazType) :: KolMat(nKolMat)
      INTEGER nrow,ncol
      REAL(8) FullMat(nrow,ncol)   
      
      INTEGER i,r,c
      REAL(8) kappa   

      DO i=1,nKolMat
        KolMat(i)%nz=0
        KolMat(i)%meja = kappa * KolMat(i)%meanabs
        DO r=KolMat(i)%ulr+1,KolMat(i)%ulr+KolMat(i)%nrow-1 !+1 da ne zanemarjam scaling waveleta
          DO c=KolMat(i)%ulc+1,KolMat(i)%ulc+KolMat(i)%ncol-1 
              IF (ABS(FullMat(r,c)).LT.KolMat(i)%meja.OR.FullMat(r,c).EQ.0.0D0) THEN ! zanemari
                FullMat(r,c)=0.0D0
                KolMat(i)%nz=KolMat(i)%nz+1
              END IF      
          END DO
        END DO
      END DO
      
      END

!C     ------------------------------------------------------------------
!      SUBROUTINE Kol_Set2Zero(FullMat,nrow,ncol,KolMat,nKolMat,kappa)
!C
!c     Naredi in-situ zanemarjanje clenov v WAWT matriki
!C      
!C     ------------------------------------------------------------------
!      USE inc_types 
!      
!      INTEGER nKolMat
!      TYPE(MatrixKolazType) :: KolMat(nKolMat)
!      INTEGER nrow,ncol
!      REAL(8) FullMat(nrow,ncol)   
!      
!      INTEGER i,r,c
!      REAL(8) kappa   
!      LOGICAL VsiNic
!
!      DO i=1,nKolMat
!        KolMat(i)%nz=0
!        KolMat(i)%meja = kappa * KolMat(i)%meanabs
!        DO r=KolMat(i)%ulr+1,KolMat(i)%ulr+KolMat(i)%nrow-1 !+1 da ne zanemarjam scaling waveleta
!          VsiNic=.TRUE.
!          DO c=KolMat(i)%ulc+1,KolMat(i)%ulc+KolMat(i)%ncol-2 ! ne smem kompletno cele vrste zanemariti.
!              IF (ABS(FullMat(r,c)).LT.KolMat(i)%meja.OR.FullMat(r,c).EQ.0.0D0) THEN ! zanemari
!                FullMat(r,c)=0.0D0
!                KolMat(i)%nz=KolMat(i)%nz+1
!              ELSE
!                VsiNic=.FALSE.
!              END IF      
!          END DO
!          c=KolMat(i)%ulc+KolMat(i)%ncol-1  ! zadnji, ce so bili vsi ostali 0, potem ta ne sme biti
!          IF (VsiNic.EQ..FALSE.) THEN ! niso vsi nic
!            IF (ABS(FullMat(r,c)).LT.KolMat(i)%meja.OR.FullMat(r,c).EQ.0.0D0) THEN ! zanemari
!              FullMat(r,c)=0.0D0
!              KolMat(i)%nz=KolMat(i)%nz+1
!            END IF
!          ELSE ! vsi so nic, torej zadnji ne sme biti
!            IF (FullMat(r,c).EQ.0.0D0) THEN ! ce je slucajno nic, ga umetno malo dvignem 
!              FullMat(r,c)=1.0D-14
!            END IF
!          END IF
!        END DO
!      END DO
!      
!      END


C     ------------------------------------------------------------------
      SUBROUTINE Kol_WaWT(FullMat,nrow,ncol,KolMat,nKolMat,M,h,g)
C
c     Naredi in-situ valcno transformacijo WaWT kolaz matrike  
C      
C     ------------------------------------------------------------------
      USE inc_types 
      
      INTEGER nKolMat
      TYPE(MatrixKolazType) :: KolMat(nKolMat)
      INTEGER nrow,ncol
      REAL(8) FullMat(nrow,ncol)      

      INTEGER M
      REAL(8) h(2*M),g(2*M)  ! wavelet coefficients

      INTEGER i
      
      DO i=1,nKolMat
        CALL isFWTxAcrhg(FullMat,nrow,ncol,KolMat(i)%meanabs,KolMat(i)%ulr,KolMat(i)%ulc,KolMat(i)%nrow,KolMat(i)%ncol,M,h,g)
      END DO

      END


C     ------------------------------------------------------------------
      SUBROUTINE Kol_FWTxV(vec,length,KolVec,nKolVec,M,h,g)
C
c     Naredi is-situ valcno transformacijo kolaz vektorja  
C      
C     ------------------------------------------------------------------
      USE inc_types 
      INTEGER nKolVec,length
      TYPE(VectorKolazType) :: KolVec(nKolVec)
      REAL(8) vec(length)

      INTEGER M
      REAL(8) h(2*M),g(2*M)  ! wavelet coefficients

      INTEGER i
      
      DO i=1,nKolVec
        CALL isFWTxVhg(vec,length,KolVec(i)%start,KolVec(i)%length,M,h,g)
      END DO

      END

C     ------------------------------------------------------------------
      SUBROUTINE Kol_iFWTxV(vec,length,KolVec,nKolVec,M,h,g)
C
c     Naredi inverzno valcno transformacijo kolaz vektorja  
C      
C     ------------------------------------------------------------------
      USE inc_types 
      INTEGER nKolVec,length
      TYPE(VectorKolazType) :: KolVec(nKolVec)
      REAL(8) vec(length)

      INTEGER M
      REAL(8) h(2*M),g(2*M)  ! wavelet coefficients

      INTEGER i
      
      DO i=1,nKolVec
        CALL isiFWTxVhg(vec,length,KolVec(i)%start,KolVec(i)%length,M,h,g)
      END DO

      END      



C     ------------------------------------------------------------------    
      SUBROUTINE Kol_MatTecOut(nKolMat,KolMat)
C
c     Izris strukture matrike v tecplot  
C      
C     ------------------------------------------------------------------
      USE inc_types 
      INTEGER nKolMat,i
      TYPE(MatrixKolazType) :: KolMat(nKolMat)

      OPEN(UNIT=12,FILE="KolMat.lay",STATUS="UNKNOWN")
      WRITE (12,'(A)') "#!MC 1100"
      
      DO i=1,nKolMat
        WRITE (12,'(A,F15.5,A,F15.5,A)') 
     &  "$!ATTACHGEOM GEOMTYPE = RECTANGLE ANCHORPOS { X = ",DBLE(KolMat(i)%ulc)," Y = ",-DBLE(KolMat(i)%ulr)," }"
        WRITE (12,'(A)') "RAWDATA"
        WRITE (12,'(2F15.5)') DBLE(KolMat(i)%ncol),-DBLE(KolMat(i)%nrow)
      END DO
      
      CLOSE (12) 

      END


C     ------------------------------------------------------------------    
      SUBROUTINE Kva_MatTecOut(nKolMat,KolMat)
C
c     Izris strukture matrike v tecplot  
C      
C     ------------------------------------------------------------------
      USE inc_types 
      INTEGER nKolMat,i
      TYPE(MatrixKolazType) :: KolMat(nKolMat)

      OPEN(UNIT=12,FILE="KvaMat.lay",STATUS="UNKNOWN")
      WRITE (12,'(A)') "#!MC 1100"
      
      DO i=1,nKolMat
        WRITE (12,'(A,F15.5,A,F15.5,A)') 
     &  "$!ATTACHGEOM GEOMTYPE = RECTANGLE ANCHORPOS { X = ",DBLE(KolMat(i)%ulc)," Y = ",-DBLE(KolMat(i)%ulr)," }"
        WRITE (12,'(A)') "RAWDATA"
        WRITE (12,'(2F15.5)') DBLE(KolMat(i)%ncol),-DBLE(KolMat(i)%nrow)
      END DO
      
      CLOSE (12) 

      END

      

C     ------------------------------------------------------------------    
      SUBROUTINE Kva_CountDivVec(n,nVecDiv,kvas)
C
C     Presteje seznam delov tipa 2^n, pri cemer najprej kvadratke
C      
C     ------------------------------------------------------------------
      USE inc_types 
      
      INTEGER n,k,nVecDiv,st
      INTEGER kvas ! velikost kvadratka
      INTEGER i,ostane
      
      nVecDiv=n/kvas
      ostane=n-kvas*nVecDiv
      
      CALL DvaNaKp(kvas,k)
      
      DO i=k-1,0,-1
        st=ostane/2**i
        IF (st.NE.0) THEN
          ostane=ostane-2**i
          nVecDiv=nVecDiv+1
        END IF
      END DO
      
      END


C     ------------------------------------------------------------------    
      SUBROUTINE Kva_DivVec(n,nVecDiv,VecDiv,kvas)
C
C     Naredi seznam delov tipa 2^n, pri cemer najprej kvadratke
C      
C     ------------------------------------------------------------------
      USE inc_types 
      
      INTEGER n,k,nVecDiv,st
      TYPE(VectorKolazType) :: vecDiv(nVecDiv)
      INTEGER kvas ! velikost kvadratka
      INTEGER stkvas,i,ostane
      
      stkvas=n/kvas
      
      DO i=1,stkvas
        vecDiv(i)%length=kvas
        vecDiv(i)%start=1+(i-1)*kvas
      END DO
      
      ostane=n-kvas*stkvas
      
      CALL DvaNaKp(kvas,k)
      
      DO i=k-1,0,-1
        st=ostane/2**i
        IF (st.NE.0) THEN
          ostane=ostane-2**i
          stkvas=stkvas+1
          vecDiv(stkvas)%length=2**i
          vecDiv(stkvas)%start=vecDiv(stkvas-1)%start+vecDiv(stkvas-1)%length      
        END IF
      END DO
      
      END
      

C     ------------------------------------------------------------------    
      SUBROUTINE Kol_DivVec(n,nVecDiv,VecDiv)
C
C     Naredi seznam delov tipa 2^n  
C      
C     ------------------------------------------------------------------
      USE inc_types 
      
      INTEGER n,k,nn,nVecDiv,st
      TYPE(VectorKolazType) :: vecDiv(nVecDiv)
      
      nn=n+1
      k=0
      st=0

10    CONTINUE
        nn=nn-2**k      
c        k=INT(log(dble(nn))/log(2.0D0))
        CALL DvaNaKp(nn,k)
        st=st+1
        vecDiv(st)%length=2**k
      
        IF (st.GT.1) THEN
          vecDiv(st)%start=vecDiv(st-1)%start+vecDiv(st-1)%length
        ELSE
          vecDiv(st)%start=1
        END IF      
      IF (2**k.NE.nn) GOTO 10

      END

C     ------------------------------------------------------------------    
      SUBROUTINE Kol_CountDivVec(n,st)
C
C     Naredi seznam delov tipa 2^n  
C      
C     ------------------------------------------------------------------
      INTEGER n,k,nn,st
      
      nn=n+1
      k=0
      st=0

10    CONTINUE
        nn=nn-2**k      
c        k=INT(log(dble(nn))/log(2.0D0))
        CALL DvaNaKp(nn,k)
        st=st+1
      IF (2**k.NE.nn) GOTO 10


      END


C     ------------------------------------------------------------------    
      SUBROUTINE Kva_DivMat(KolMat,nKolMat,RowDiv,nRowDiv,ColDiv,nColDiv)
C
C     Naredi seznam delov tipa 2^n x 2^m, najprej kvadratki  
C      
C     ------------------------------------------------------------------
      USE inc_types 

      INTEGER i,j,st
      INTEGER nkolMat,nrowdiv,ncoldiv
      
      TYPE(MatrixKolazType) :: KolMat(nKolMat)
      TYPE(VectorKolazType) :: rowDiv(nrowDiv)
      TYPE(VectorKolazType) :: colDiv(ncolDiv)

      st=0
      DO i=1,nColDiv
        DO j=1,nRowDiv
          st=st+1
          KolMat(st)%ulr=rowDiv(j)%start
          KolMat(st)%ulc=colDiv(i)%start
          KolMat(st)%nrow=rowDiv(j)%length
          KolMat(st)%ncol=colDiv(i)%length
        END DO
      END DO
      
      END


C     ------------------------------------------------------------------    
      RECURSIVE SUBROUTINE Kol_DivMat(nrow,ncol,fulr,fulc,ns,KolMat,st)
C
C     Naredi seznam delov tipa 2^n x 2^m  
C      
C     ------------------------------------------------------------------
      USE inc_types 

      INTEGER nrow,ncol,ns,st
      INTEGER fulr,fulc ! father upper left corner row,col
      INTEGER kr,kc,sinNrow,sinNcol,ulc,ulr
      
      TYPE(MatrixKolazType) :: KolMat(ns)

c     2**kr je prva potenca stevila 2 manjsa ali enaka nrow      
c      kr=INT(log(dble(nrow))/log(2.0D0))
c      kc=INT(log(dble(ncol))/log(2.0D0))
      CALL DvaNaKp(nrow,kr)
      CALL DvaNaKp(ncol,kc)      
      
c     prvi sin A      
      sinNrow=2**kr
      sinNcol=2**kc
      ulc=fulc
      ulr=fulr
      
      
      IF (sinNrow.NE.Nrow.AND.sinNcol.NE.ncol) THEN
c        print *,"1sin",sinNrow,sinNcol
        CALL Kol_DivMat(sinNrow,sinNcol,ulr,ulc,ns,KolMat,st)
        sinNrow=nrow-2**kr
        sinNcol=2**kc
        ulr=fulr+2**kr
        ulc=fulc
c        print *,"2sin",sinNrow,sinNcol
        CALL Kol_DivMat(sinNrow,sinNcol,ulr,ulc,ns,KolMat,st)
        sinNrow=2**kr
        sinNcol=ncol-2**kc
        ulr=fulr
        ulc=fulc+2**kc
c        print *,"3sin",sinNrow,sinNcol
        CALL Kol_DivMat(sinNrow,sinNcol,ulr,ulc,ns,KolMat,st)
        sinNrow=nrow-2**kr
        sinNcol=ncol-2**kc
        ulr=fulr+2**kr
        ulc=fulc+2**kc        
c        print *,"4sin",sinNrow,sinNcol
        CALL Kol_DivMat(sinNrow,sinNcol,ulr,ulc,ns,KolMat,st)

      ELSE IF (sinNrow.EQ.Nrow.AND.sinNcol.NE.ncol) THEN
c        print *,"1sin",sinNrow,sinNcol
        CALL Kol_DivMat(sinNrow,sinNcol,ulr,ulc,ns,KolMat,st)
        sinNrow=2**kr
        sinNcol=ncol-2**kc
        ulr=fulr
        ulc=fulc+2**kc        
c        print *,"2sin",sinNrow,sinNcol
        CALL Kol_DivMat(sinNrow,sinNcol,ulr,ulc,ns,KolMat,st)

      ELSE IF (sinNrow.NE.Nrow.AND.sinNcol.EQ.ncol) THEN
c        print *,"1sin",sinNrow,sinNcol
        CALL Kol_DivMat(sinNrow,sinNcol,ulr,ulc,ns,KolMat,st)
        sinNrow=nrow-2**kr
        sinNcol=2**kc
        ulr=fulr+2**kr
        ulc=fulc
c        print *,"2sin",sinNrow,sinNcol
        CALL Kol_DivMat(sinNrow,sinNcol,ulr,ulc,ns,KolMat,st)
      ELSE
        st=st+1
        KolMat(st)%ulr=fulr
        KolMat(st)%ulc=fulc
        KolMat(st)%nrow=nrow
        KolMat(st)%ncol=ncol
c        print *,"nimas",nrow,ncol,fulr,fulc,ns
      END IF
      END

      
      
C     ------------------------------------------------------------------    
      RECURSIVE SUBROUTINE Kol_CountDivMat(nrow,ncol,fulr,fulc,ns)
C
C     Presteje na koliko delov tipa 2^n x 2^m bo razdeljena matrika 
C      
C     ------------------------------------------------------------------
      INTEGER nrow,ncol,ns
      INTEGER fulr,fulc ! father upper left corner row,col
      INTEGER kr,kc,sinNrow,sinNcol,ulc,ulr

c     2**kr je prva potenca stevila 2 manjsa ali enaka nrow      
c      kr=INT(log(dble(nrow))/log(2.0D0))
c      kc=INT(log(dble(ncol))/log(2.0D0))
      CALL DvaNaKp(nrow,kr)
      CALL DvaNaKp(ncol,kc)      
      
c     prvi sin A      
      sinNrow=2**kr
      sinNcol=2**kc
      ulc=fulc
      ulr=fulr
      
      
      IF (sinNrow.NE.Nrow.AND.sinNcol.NE.ncol) THEN
c        print *,"1sin",sinNrow,sinNcol
        CALL Kol_CountDivMat(sinNrow,sinNcol,ulr,ulc,ns)
        sinNrow=nrow-2**kr
        sinNcol=2**kc
        ulr=fulr+2**kr
        ulc=fulc
c        print *,"2sin",sinNrow,sinNcol
        CALL Kol_CountDivMat(sinNrow,sinNcol,ulr,ulc,ns)
        sinNrow=2**kr
        sinNcol=ncol-2**kc
        ulr=fulr
        ulc=fulc+2**kc
c        print *,"3sin",sinNrow,sinNcol
        CALL Kol_CountDivMat(sinNrow,sinNcol,ulr,ulc,ns)
        sinNrow=nrow-2**kr
        sinNcol=ncol-2**kc
        ulr=fulr+2**kr
        ulc=fulc+2**kc        
c        print *,"4sin",sinNrow,sinNcol
        CALL Kol_CountDivMat(sinNrow,sinNcol,ulr,ulc,ns)

      ELSE IF (sinNrow.EQ.Nrow.AND.sinNcol.NE.ncol) THEN
c        print *,"1sin",sinNrow,sinNcol
        CALL Kol_CountDivMat(sinNrow,sinNcol,ulr,ulc,ns)
        sinNrow=2**kr
        sinNcol=ncol-2**kc
        ulr=fulr
        ulc=fulc+2**kc        
c        print *,"2sin",sinNrow,sinNcol
        CALL Kol_CountDivMat(sinNrow,sinNcol,ulr,ulc,ns)

      ELSE IF (sinNrow.NE.Nrow.AND.sinNcol.EQ.ncol) THEN
c        print *,"1sin",sinNrow,sinNcol
        CALL Kol_CountDivMat(sinNrow,sinNcol,ulr,ulc,ns)
        sinNrow=nrow-2**kr
        sinNcol=2**kc
        ulr=fulr+2**kr
        ulc=fulc
c        print *,"2sin",sinNrow,sinNcol
        CALL Kol_CountDivMat(sinNrow,sinNcol,ulr,ulc,ns)
      ELSE
        ns=ns+1
c        print *,"nimas",nrow,ncol,fulr,fulc,ns
      END IF
      END
      
      

c -----------------------------------------------------------------------------
c   Fast Wavelet Transform of a Matrix W*A*W^-1 by tranforming 
c   individual columns and rows
      SUBROUTINE isFWTxAcrhg(a,sizer,sizec,mean,ulr,ulc,nrow,ncol,M,h,g)
      INTEGER M,ii,jj,nrow,ncol,sizer,sizec,ulr,ulc
      REAL(8) a(sizer,sizec),mean
      REAL(8) h(2*M),g(2*M)
      REAL(8), ALLOCATABLE :: cr(:)
 
      ALLOCATE (cr(ncol))
      
c     najprej transformiramo vse stolpce
      DO ii=ulr,ulr+nrow-1
        DO jj=ulc,ulc+ncol-1
          cr(jj-ulc+1)=a(ii,jj)
        ENDDO
        CALL FWTxVhg(cr,ncol,M,h,g)
        DO jj=ulc,ulc+ncol-1
          a(ii,jj)=cr(jj-ulc+1)
        ENDDO
      ENDDO
      DEALLOCATE (cr)
      
      ALLOCATE (cr(nrow))      
      mean=0
c     nato se vrstice iz ze transformiranih stolpcev
      DO ii=ulc,ulc+ncol-1
        DO jj=ulr,ulr+nrow-1
          cr(jj-ulr+1)=a(jj,ii)
        ENDDO
        CALL FWTxVhg(cr,nrow,M,h,g)
        DO jj=ulr,ulr+nrow-1
          a(jj,ii)=cr(jj-ulr+1)
          mean=mean+ABS(a(jj,ii))
        ENDDO
      ENDDO
      DEALLOCATE (cr)
      
      mean=mean/DBLE(nrow)/DBLE(ncol)

      RETURN
      END

      
c -----------------------------------------------------------------------------
c   Fast Wavelet Transform, kjer ze imam h,g
      SUBROUTINE FWTxVhg(d,nvek,M,h,g)
c FWT transformacija vektorja
c M je stevilo vanishing moments. 
c M=1 imamo Haarove valcke
c M=2 imamo Daubechies D4 valcke
c M= in tako naprej, M je stevilo vanishing moments

      INTEGER n,j,k,maxk,M,nvek
      REAL*8 d(nvek),s(nvek),snov(nvek),h(2*M),g(2*M)

      
      n=int(log(real(nvek))/log(2.0)+0.49)

      s=d

      DO j=0,n-1
            maxk=2**(n-j-1)
            DO k=1,maxk
                  CALL Enk(h,g,M,s,nvek,k,d(2**(n-j-1)+k),snov(k),2**(n-j))
            ENDDO
            DO k=1,maxk
                  s(k)=snov(k)
            ENDDO
      ENDDO
c na prvo mesto v vektor d gre s(n,1), ker ga rabim za obratno transf. (d(2**(n-j-1)+k-1)=s(1))
      d(1)=s(1)

      RETURN
      END      


c -----------------------------------------------------------------------------
      SUBROUTINE Enk(h,g,M,s,nvek,k,dven,sven,Sts)
c      Ta podprogram mnozi za en "k".
      INTEGER nvek,ii,M,k,Sts,StsPeri,stevs
      REAL*8 g(2*M),h(2*M),s(nvek),dven,sven,tmp
      dven=0
      tmp=0

c StsPeri je stevilo s-jev, ki jih moram narediti periodicno
      StsPeri=(2*M+2*k-2)-Sts

      IF (StsPeri.LE.0) THEN
            DO ii=1,2*M
c            enacbi 3.13 in 3.14
            dven=dven+g(ii)*s(ii+2*k-2)
          tmp=tmp+h(ii)*s(ii+2*k-2)
            ENDDO
      ELSE
c            Najprej normalni del
            DO ii=1,2*M-StsPeri
                  dven=dven+g(ii)*s(ii+2*k-2)
                tmp=tmp+h(ii)*s(ii+2*k-2)
            ENDDO
c            Nato tiste, ki so prek - zacnem s(1) in tako naprej
c pazi, lahko se zgodi, da je potrebno periodo orbniti veckrat
c to je takrat, ko je StsPeri vecje od Sts
            DO ii=2*M-StsPeri+1,2*M
              stevs=ii-(2*M-StsPeri)
10                  IF (stevs.GT.sts) THEN
                        stevs=stevs-sts
                        goto 10
                  ENDIF
                  dven=dven+g(ii)*s(stevs)
                tmp=tmp+h(ii)*s(stevs)
            ENDDO
      ENDIF
      sven=tmp
      RETURN 
      END
      
c -----------------------------------------------------------------------------
c   Fast Wavelet Transform, kjer ze imam h,g
      SUBROUTINE isFWTxVhg(d,length,start,nvek,M,h,g)
c FWT transformacija vektorja, vendar je del vektorja d, dolzine nvek, zacenjsi s start
c
c     vektor d je dolg length
c     transformiram pa le njegov del : od start do start+nvek-1
c
c M je stevilo vanishing moments. 
c M=1 imamo Haarove valcke
c M=2 imamo Daubechies D4 valcke
c M= in tako naprej, M je stevilo vanishing moments

      INTEGER n,j,k,maxk,M,nvek,start,length
      REAL*8 d(length),s(nvek),snov(nvek),h(2*M),g(2*M)

      
      n=int(log(real(nvek))/log(2.0)+0.49)

      DO j=start,start+nvek-1
        s(j-start+1)=d(j)
      END DO

      DO j=0,n-1
            maxk=2**(n-j-1)
            DO k=1,maxk
                  CALL Enk(h,g,M,s,nvek,k,d(start-1+2**(n-j-1)+k),snov(k),2**(n-j))
            ENDDO
            DO k=1,maxk
                  s(k)=snov(k)
            ENDDO
      ENDDO
c na prvo mesto v vektor d gre s(n,1), ker ga rabim za obratno transf. (d(2**(n-j-1)+k-1)=s(1))
      d(start)=s(1)

      RETURN
      END


c -----------------------------------------------------------------------------
c   Inverse Fast Wavelet Transform
      SUBROUTINE isiFWTxVhg(v,length,start,nvek,M,h,g)
c     vektor v je dolg length
c     transformiram pa le njegov del : od start do start+nvek-1  
      INTEGER n,j,maxk,M,ii,nvek,length,start
      REAL*8 v(length),s(nvek),snov(nvek),h(2*M),g(2*M),s2ii,s2iim1

      
      n=int(log(real(nvek))/log(2.0)+0.49)

c M je stevilo vanishing moments. M=1 imamo Haarove valcke
      
      s(1)=v(start)
c To je pretvorba med d-jem iz clanka in vektorjem v, v katerega shranjujem
c      d(j,k)=v(2**(n-j+1)+k)
c     s(n,1)=v(1)

      DO j=n,1,-1
            maxk=2**(n-j)
            DO ii=1,maxk
                  CALL iEnk(h,g,M,s,v(start:start+nvek-1),nvek,j,n,ii,s2ii,s2iim1,maxk)
                  snov(2*ii)=s2ii
                  snov(2*ii-1)=s2iim1
            ENDDO
            DO ii=1,2*maxk
                  s(ii)=snov(ii)
            ENDDO
      ENDDO
 
      DO j=start,start+nvek-1
        v(j)=s(j-start+1)
      END DO

      RETURN
      END       
      
c -----------------------------------------------------------------------------
      SUBROUTINE iEnk(h,g,M,s,v,nvek,j,n,ii,s2ii,s2iim1,maxk)
      INTEGER nvek,ii,M,k,j,maxk,sts,n
      REAL*8 g(2*M),h(2*M),s(nvek),v(nvek),s2ii,s2iim1

      s2ii=0
      s2iim1=0

      DO k=1,M
c            enacbi 3.15
            IF (ii-k+1.GE.1) THEN
      s2ii=s2ii    +h(2*k)  *s(ii-k+1)+g(2*k)  *v(2**(n-j)+ii-k+1)
      s2iim1=s2iim1+h(2*k-1)*s(ii-k+1)+g(2*k-1)*v(2**(n-j)+ii-k+1)
            ELSE
c pazi mozno je, da je periodo potrebno ponoviti veckrat
c stevilo s-ov je maxk
c od maxk odstevamo toliko (ii-k+1)-jev, da je (ii-k+1) manjse od maxk
            sts=maxk+(ii-k+1)
10            IF (sts.LT.1) THEN
c            Print *,"sts=",sts
            sts=sts+maxk
            goto 10
            ENDIF

      s2ii=s2ii+
     &h(2*k)  *s(sts)+g(2*k)  *v(2**(n-j)+sts)
      s2iim1=s2iim1+
     &h(2*k-1)*s(sts)+g(2*k-1)*v(2**(n-j)+sts)
            ENDIF
      ENDDO

      RETURN
      END      
      
c -----------------------------------------------------------------------------
c * Coefficients for Daubechies wavelets 1-38
c * -----------------------------------------
c *
c * Computed by Kazuo Hatano, Aichi Institute of Technology.
c * ftp://phase.etl.go.jp/pub/phase/wavelet/index.html
c *
c * Compiled and verified by Olli Niemitalo.
c * http://www.biochem.oulu.fi/~oniemita/dsp/daub.h
c * 
c * Rewriten into FORTRAN by Jure Ravnik
c * jure.ravnik@uni-mb.si
c *
c * Last change : 2004-01-08
c *
c * Discrete Wavelet Transformation (DWT) breaks a signal down into
c * subbands distributed logarithimically in frequency, each sampled
c * at a rate that has a natural proportion to the frequencies in that
c * band. The traditional fourier transformation has no time domain
c * resolution at all, or when done using many short windows on a
c * longer data, equal resolution at all frequencies. The distribution
c * of samples in the time and frequency domain by DWT is of form:
c *
c * log f
c *  |XXXXXXXXXXXXXXXX  X = a sample
c *  |X X X X X X X X   f = frequency
c *  |X   X   X   X     t = time
c *  |X       X
c *  |X
c *   ----------------t
c *
c * Single
c * subband decomposition         and      reconstruction:
c *
c *     -> high -> decimate -------------> dilute -> high
c *    |   pass    by 2      high subband  by 2      pass \
c * in |                                                   + out
c *    |                                                  /  =in
c *     -> low  -> decimate -------------> dilute -> low
c *        pass    by 2      low subband   by 2      pass
c *
c * This creates two subbands from the input signal, both sampled at half
c * the original frequency. The filters approximate halfband FIR filters
c * and are determined by the choice of wavelet. Using Daubechies wavelets
c * (and most others), the data can be reconstructed to the exact original
c * even when the halfband filters are not perfect. Note that the amount
c * of information (samples) stays the same throughout the operation.
c *
c * Decimation by 2: ABCDEFGHIJKLMNOPQR -> ACEGIKMOQ
c *   Dilution by 2: ACEGIKMOQ -> A0C0E0G0I0K0M0O0Q0
c *
c * To get the logarithmic resolution in frequency, the low subband is
c * re-transformed, and again, the low subband from this transformation
c * gets the same treatment etc.
c *
c * Decomposition:
c *
c *     -> high -> decimate --------------------------------> subband0
c *    |   pass    by 2
c * in |                 -> high -> decimate ---------------> subband1
c *    |                |   pass    by 2
c *     -> low -> decim |                 -> high -> decim -> subband2
c *        pass   by 2  |                |   pass    by 2
c *                      -> low -> decim |
c *                         pass   by 2  |   .   down to what suffices
c *                                       -> .    or if periodic data,
c *                                          .     until short of data
c * Reconstruction:
c *
c * subband0 -----------------------------------> dilute -> high
c *                                               by 2      pass \
c * subband1 ------------------> dilute -> high                   + out
c *                              by 2      pass \                /  =in
c * subband2 -> dilute -> high                   + dilute -> low
c *             by 2      pass \                /  by 2      pass
c *                             + dilute -> low
c * Start   .                  /  by 2      pass
c * here!   . -> dilute -> low
c *         .    by 2      pass
c *
c * In a real-time application, the filters introduce delays, so you need
c * to compensate them by adding additional delays to less-delayed higher
c * bands, to get the summation work as intended.
c *
c * For periodic signals or windowed operation, this problem doesn't exist -
c * a single subband transformation is a matrix multiplication, with wrapping
c * implemented in the matrix:
c *
c * Decomposition:
c *
c * |L0|   |C0  C1  C2  C3                | |I0|     L = lowpass output
c * |H0|   |C3 -C2  C1 -C0                | |I1|     H = highpass output
c * |L1|   |        C0  C1  C2  C3        | |I2|     I = input
c * |H1| = |        C3 -C2  C1 -C0        | |I3|     C = coefficients
c * |L2|   |                C0  C1  C2  C3| |I4|
c * |H2|   |                C3 -C2  C1 -C0| |I5|
c * |L3|   |C2  C3                  C0  C1| |I6|
c * |H3|   |C1 -C0                  C3 -C2| |I7|     Daubechies 4-coef:
c *
c *      1+sqrt(3)        3+sqrt(3)        3-sqrt(3)        1-sqrt(3)
c * C0 = ---------   C1 = ---------   C2 = ---------   C3 = ---------
c *      4 sqrt(2)        4 sqrt(2)        4 sqrt(2)        4 sqrt(2)
c *
c * Reconstruction:
c *
c * |I0|   |C0  C3                  C2  C1| |L0|
c * |I1|   |C1 -C2                  C3 -C0| |H0|
c * |I2|   |C2  C1  C0  C3                | |L1|
c * |I3| = |C3 -C0  C1 -C2                | |H1|
c * |I4|   |        C2  C1  C0  C3        | |L2|
c * |I5|   |        C3 -C0  C1 -C2        | |H2|
c * |I6|   |                C2  C1  C0  C3| |L3|
c * |I7|   |                C3 -C0  C1 -C2| |H3|
c *
c * This file contains the lowpass FIR filter coefficients. Highpass
c * coefficients you get by reversing tap order and multiplying by
c * sequence 1,-1, 1,-1, ... Because these are orthogonal wavelets, the
c * analysis and reconstruction coefficients are the same.
c *
c * A coefficient set convolved by its reverse is an ideal halfband lowpass
c * filter multiplied by a symmetric windowing function. This creates the
c * kind of symmetry in the frequency domain that enables aliasing-free
c * reconstruction. Daubechies wavelets are the minimum-phase, minimum
c * number of taps solutions for a number of vanishing moments (seven in
c * Daub7 etc), which determines their frequency selectivity.
c */
c ----------------------------------------------------------------------
      SUBROUTINE FiltCoef(h,g,M)
c     The subroutine returns h and g filter coefficients with
c     M vanishing moments
      INTEGER M,k,MaxM
      REAL*8 h(2*M),g(2*M)

      MaxM=38

      IF (M.GT.MaxM) THEN
      WRITE (*,*) " "
      WRITE (*,'(A,I2,A)') "FiltCoef : Warning : 
     &Wavelets with more than ",MaxM," vanishing moments"
      WRITE (*,'(A,I2,A)') "                     were not implemented.
     & Using M=",MaxM," instead!"
      WRITE (*,*) " "
      M=MaxM
      ENDIF

      IF (M.LT.1) THEN
      WRITE (*,*) " "
      WRITE (*,*) "FiltCoef : Error : 
     &Wavelets must have at least one vanishing moment."
      stop " "
      ENDIF


      IF (M.EQ.1) THEN
      h(1)=
     &7.071067811865475244008443621048490392848359376884740365883398d-01
      h(2)=
     &7.071067811865475244008443621048490392848359376884740365883398d-01
      ENDIF

        IF (M.EQ.2) THEN
      h(1)=
     &4.829629131445341433748715998644486838169524195042022752011715d-01
      h(2)=
     &8.365163037378079055752937809168732034593703883484392934953414d-01
      h(3)=
     &2.241438680420133810259727622404003554678835181842717613871683d-01
      h(4)=
     &-1.294095225512603811744494188120241641745344506599652569070016d-1
      ENDIF

        IF (M.EQ.3) THEN
      h(1)=
     &3.326705529500826159985115891390056300129233992450683597084705d-01
      h(2)=
     &8.068915093110925764944936040887134905192973949948236181650920d-01
      h(3)=
     &4.598775021184915700951519421476167208081101774314923066433867d-01
      h(4)=
     &-1.350110200102545886963899066993744805622198452237811919756862d-1
      h(5)=
     &-8.544127388202666169281916918177331153619763898808662976351748d-2
      h(6)=
     &3.522629188570953660274066471551002932775838791743161039893406d-02
        ENDIF

        IF (M.EQ.4) THEN
      h(1)=
     &2.303778133088965008632911830440708500016152482483092977910968d-01
      h(2)=
     &7.148465705529156470899219552739926037076084010993081758450110d-01
      h(3)=
     &6.308807679298589078817163383006152202032229226771951174057473d-01
      h(4)=
     &-2.798376941685985421141374718007538541198732022449175284003358d-2
      h(5)=
     &-1.870348117190930840795706727890814195845441743745800912057770d-1
      h(6)=
     &3.084138183556076362721936253495905017031482172003403341821219d-02
      h(7)=
     &3.288301166688519973540751354924438866454194113754971259727278d-02
      h(8)=
     &-1.059740178506903210488320852402722918109996490637641983484974d-2
        ENDIF


        IF (M.EQ.5) THEN
      h(1)=
     &1.601023979741929144807237480204207336505441246250578327725699d-01
      h(2)=
     &6.038292697971896705401193065250621075074221631016986987969283d-01
      h(3)=
     &7.243085284377729277280712441022186407687562182320073725767335d-01
      h(4)=
     &1.384281459013207315053971463390246973141057911739561022694652d-01
      h(5)=
     &-2.422948870663820318625713794746163619914908080626185983913726d-1
      h(6)=
     &-3.224486958463837464847975506213492831356498416379847225434268d-2
      h(7)=
     &7.757149384004571352313048938860181980623099452012527983210146d-02
      h(8)=
     &-6.241490212798274274190519112920192970763557165687607323417435d-3
      h(9)=
     &-1.258075199908199946850973993177579294920459162609785020169232d-2
      h(10)=
     &3.335725285473771277998183415817355747636524742305315099706428d-03
      ENDIF

        IF (M.EQ.6) THEN
      h(1)=
     &1.115407433501094636213239172409234390425395919844216759082360d-01
      h(2)=
     &4.946238903984530856772041768778555886377863828962743623531834d-01
      h(3)=
     &7.511339080210953506789344984397316855802547833382612009730420d-01
      h(4)=
     &3.152503517091976290859896548109263966495199235172945244404163d-01
      h(5)=
     &-2.262646939654398200763145006609034656705401539728969940143487d-1
      h(6)=
     &-1.297668675672619355622896058765854608452337492235814701599310d-1
      h(7)=
     &9.750160558732304910234355253812534233983074749525514279893193d-02
      h(8)=
     &2.752286553030572862554083950419321365738758783043454321494202d-02
      h(9)=
     &-3.158203931748602956507908069984866905747953237314842337511464d-2
      h(10)=
     &5.538422011614961392519183980465012206110262773864964295476524d-04
      h(11)=
     &4.777257510945510639635975246820707050230501216581434297593254d-03
      h(12)=
     &-1.077301085308479564852621609587200035235233609334419689818580d-3
      ENDIF

        IF (M.EQ.7) THEN
      h(1)=
     &7.785205408500917901996352195789374837918305292795568438702937d-02
      h(2)=
     &3.965393194819173065390003909368428563587151149333287401110499d-01
      h(3)=
     &7.291320908462351199169430703392820517179660611901363782697715d-01
      h(4)=
     &4.697822874051931224715911609744517386817913056787359532392529d-01
      h(5)=
     &-1.439060039285649754050683622130460017952735705499084834401753d-1
      h(6)=
     &-2.240361849938749826381404202332509644757830896773246552665095d-1
      h(7)=
     &7.130921926683026475087657050112904822711327451412314659575113d-02
      h(8)=
     &8.061260915108307191292248035938190585823820965629489058139218d-02
      h(9)=
     &-3.802993693501441357959206160185803585446196938467869898283122d-2
      h(10)=
     &-1.657454163066688065410767489170265479204504394820713705239272d-2
      h(11)=
     &1.255099855609984061298988603418777957289474046048710038411818d-02
      h(12)=
     &4.295779729213665211321291228197322228235350396942409742946366d-04
      h(13)=
     &-1.801640704047490915268262912739550962585651469641090625323864d-3
      h(14)=
     &3.537137999745202484462958363064254310959060059520040012524275d-04
      ENDIF

        IF (M.EQ.8) THEN
      h(1)=
     &5.441584224310400995500940520299935503599554294733050397729280d-02
      h(2)=
     &3.128715909142999706591623755057177219497319740370229185698712d-01
      h(3)=
     &6.756307362972898068078007670471831499869115906336364227766759d-01
      h(4)=
     &5.853546836542067127712655200450981944303266678053369055707175d-01
      h(5)=
     &-1.582910525634930566738054787646630415774471154502826559735335d-2
      h(6)=
     &-2.840155429615469265162031323741647324684350124871451793599204d-1
      h(7)=
     &4.724845739132827703605900098258949861948011288770074644084096d-04
      h(8)=
     &1.287474266204784588570292875097083843022601575556488795577000d-01
      h(9)=
     &-1.736930100180754616961614886809598311413086529488394316977315d-2
      h(10)=
     &-4.408825393079475150676372323896350189751839190110996472750391d-2
      h(11)=
     &1.398102791739828164872293057263345144239559532934347169146368d-02
      h(12)=
     &8.746094047405776716382743246475640180402147081140676742686747d-03
      h(13)=
     &-4.870352993451574310422181557109824016634978512157003764736208d-3
      h(14)=
     &-3.917403733769470462980803573237762675229350073890493724492694d-4
      h(15)=
     &6.754494064505693663695475738792991218489630013558432103617077d-04
      h(16)=
     &-1.174767841247695337306282316988909444086693950311503927620013d-4
      ENDIF

        IF (M.EQ.9) THEN
      h(1)=
     &3.807794736387834658869765887955118448771714496278417476647192d-02
      h(2)=
     &2.438346746125903537320415816492844155263611085609231361429088d-01
      h(3)=
     &6.048231236901111119030768674342361708959562711896117565333713d-01
      h(4)=
     &6.572880780513005380782126390451732140305858669245918854436034d-01
      h(5)=
     &1.331973858250075761909549458997955536921780768433661136154346d-01
      h(6)=
     &-2.932737832791749088064031952421987310438961628589906825725112d-1
      h(7)=
     &-9.684078322297646051350813353769660224825458104599099679471267d-2
      h(8)=
     &1.485407493381063801350727175060423024791258577280603060771649d-01
      h(9)=
     &3.072568147933337921231740072037882714105805024670744781503060d-02
      h(10)=
     &-6.763282906132997367564227482971901592578790871353739900748331d-2
      h(11)=
     &2.509471148314519575871897499885543315176271993709633321834164d-04
      h(12)=
     &2.236166212367909720537378270269095241855646688308853754721816d-02
      h(13)=
     &-4.723204757751397277925707848242465405729514912627938018758526d-3
      h(14)=
     &-4.281503682463429834496795002314531876481181811463288374860455d-3
      h(15)=
     &1.847646883056226476619129491125677051121081359600318160732515d-03
      h(16)=
     &2.303857635231959672052163928245421692940662052463711972260006d-04
      h(17)=
     &-2.519631889427101369749886842878606607282181543478028214134265d-4
      h(18)=
     &3.934732031627159948068988306589150707782477055517013507359938d-05
      ENDIF


        IF (M.EQ.10) THEN
      h(1)=
     &2.667005790055555358661744877130858277192498290851289932779975d-02
      h(2)=
     &1.881768000776914890208929736790939942702546758640393484348595d-01
      h(3)=
     &5.272011889317255864817448279595081924981402680840223445318549d-01
      h(4)=
     &6.884590394536035657418717825492358539771364042407339537279681d-01
      h(5)=
     &2.811723436605774607487269984455892876243888859026150413831543d-01
      h(6)=
     &-2.498464243273153794161018979207791000564669737132073715013121d-1
      h(7)=
     &-1.959462743773770435042992543190981318766776476382778474396781d-1
      h(8)=
     &1.273693403357932600826772332014009770786177480422245995563097d-01
      h(9)=
     &9.305736460357235116035228983545273226942917998946925868063974d-02
      h(10)=
     &-7.139414716639708714533609307605064767292611983702150917523756d-2
      h(11)=
     &-2.945753682187581285828323760141839199388200516064948779769654d-2
      h(12)=
     &3.321267405934100173976365318215912897978337413267096043323351d-02
      h(13)=
     &3.606553566956169655423291417133403299517350518618994762730612d-03
      h(14)=
     &-1.073317548333057504431811410651364448111548781143923213370333d-2
      h(15)=
     &1.395351747052901165789318447957707567660542855688552426721117d-03
      h(16)=
     &1.992405295185056117158742242640643211762555365514105280067936d-03
      h(17)=
     &-6.858566949597116265613709819265714196625043336786920516211903d-4
      h(18)=
     &-1.164668551292854509514809710258991891527461854347597362819235d-4
      h(19)=
     &9.358867032006959133405013034222854399688456215297276443521873d-05
      h(20)=
     &-1.326420289452124481243667531226683305749240960605829756400674d-5
      ENDIF


        IF (M.EQ.11) THEN
      h(1)=
     &1.869429776147108402543572939561975728967774455921958543286692d-02
      h(2)=
     &1.440670211506245127951915849361001143023718967556239604318852d-01
      h(3)=
     &4.498997643560453347688940373853603677806895378648933474599655d-01
      h(4)=
     &6.856867749162005111209386316963097935940204964567703495051589d-01
      h(5)=
     &4.119643689479074629259396485710667307430400410187845315697242d-01
      h(6)=
     &-1.622752450274903622405827269985511540744264324212130209649667d-1
      h(7)=
     &-2.742308468179469612021009452835266628648089521775178221905778d-1
      h(8)=
     &6.604358819668319190061457888126302656753142168940791541113457d-02
      h(9)=
     &1.498120124663784964066562617044193298588272420267484653796909d-01
      h(10)=
     &-4.647995511668418727161722589023744577223260966848260747450320d-2
      h(11)=
     &-6.643878569502520527899215536971203191819566896079739622858574d-2
      h(12)=
     &3.133509021904607603094798408303144536358105680880031964936445d-02
      h(13)=
     &2.084090436018106302294811255656491015157761832734715691126692d-02
      h(14)=
     &-1.536482090620159942619811609958822744014326495773000120205848d-2
      h(15)=
     &-3.340858873014445606090808617982406101930658359499190845656731d-3
      h(16)=
     &4.928417656059041123170739741708273690285547729915802418397458d-03
      h(17)=
     &-3.085928588151431651754590726278953307180216605078488581921562d-4
      h(18)=
     &-8.930232506662646133900824622648653989879519878620728793133358d-4
      h(19)=
     &2.491525235528234988712216872666801088221199302855425381971392d-04
      h(20)=
     &5.443907469936847167357856879576832191936678525600793978043688d-05
      h(21)=
     &-3.463498418698499554128085159974043214506488048233458035943601d-5
      h(22)=
     &4.494274277236510095415648282310130916410497987383753460571741d-06
      ENDIF


        IF (M.EQ.12) THEN
      h(1)=
     &1.311225795722951750674609088893328065665510641931325007748280d-02
      h(2)=
     &1.095662728211851546057045050248905426075680503066774046383657d-01
      h(3)=
     &3.773551352142126570928212604879206149010941706057526334705839d-01
      h(4)=
     &6.571987225793070893027611286641169834250203289988412141394281d-01
      h(5)=
     &5.158864784278156087560326480543032700677693087036090056127647d-01
      h(6)=
     &-4.476388565377462666762747311540166529284543631505924139071704d-2
      h(7)=
     &-3.161784537527855368648029353478031098508839032547364389574203d-1
      h(8)=
     &-2.377925725606972768399754609133225784553366558331741152482612d-2
      h(9)=
     &1.824786059275796798540436116189241710294771448096302698329011d-01
      h(10)=
     &5.359569674352150328276276729768332288862665184192705821636342d-03
      h(11)=
     &-9.643212009650708202650320534322484127430880143045220514346402d-2
      h(12)=
     &1.084913025582218438089010237748152188661630567603334659322512d-02
      h(13)=
     &4.154627749508444073927094681906574864513532221388374861287078d-02
      h(14)=
     &-1.221864906974828071998798266471567712982466093116558175344811d-2
      h(15)=
     &-1.284082519830068329466034471894728496206109832314097633275225d-2
      h(16)=
     &6.711499008795509177767027068215672450648112185856456740379455d-03
      h(17)=
     &2.248607240995237599950865211267234018343199786146177099262010d-03
      h(18)=
     &-2.179503618627760471598903379584171187840075291860571264980942d-3
      h(19)=
     &6.545128212509595566500430399327110729111770568897356630714552d-06
      h(20)=
     &3.886530628209314435897288837795981791917488573420177523436096d-04
      h(21)=
     &-8.850410920820432420821645961553726598738322151471932808015443d-5
      h(22)=
     &-2.424154575703078402978915320531719580423778362664282239377532d-5
      h(23)=
     &1.277695221937976658714046362616620887375960941439428756055353d-05
      h(24)=
     &-1.529071758068510902712239164522901223197615439660340672602696d-6
      ENDIF


        IF (M.EQ.13) THEN
      h(1)=
     &9.202133538962367972970163475644184667534171916416562386009703d-03
      h(2)=
     &8.286124387290277964432027131230466405208113332890135072514277d-02
      h(3)=
     &3.119963221604380633960784112214049693946683528967180317160390d-01
      h(4)=
     &6.110558511587876528211995136744180562073612676018239438526582d-01
      h(5)=
     &5.888895704312189080710395347395333927665986382812836042235573d-01
      h(6)=
     &8.698572617964723731023739838087494399231884076619701250882016d-02
      h(7)=
     &-3.149729077113886329981698255932282582876888450678789025950306d-1
      h(8)=
     &-1.245767307508152589413808336021260180792739295173634719572069d-1
      h(9)=
     &1.794760794293398432348450072339369013581966256244133393042881d-01
      h(10)=
     &7.294893365677716380902830610477661983325929026879873553627963d-02
      h(11)=
     &-1.058076181879343264509667304196464849478860754801236658232360d-1
      h(12)=
     &-2.648840647534369463963912248034785726419604844297697016264224d-2
      h(13)=
     &5.613947710028342886214501998387331119988378792543100244737056d-02
      h(14)=
     &2.379972254059078811465170958554208358094394612051934868475139d-03
      h(15)=
     &-2.383142071032364903206403067757739134252922717636226274077298d-2
      h(16)=
     &3.923941448797416243316370220815526558824746623451404043918407d-03
      h(17)=
     &7.255589401617566194518393300502698898973529679646683695269828d-03
      h(18)=
     &-2.761911234656862178014576266098445995350093330501818024966316d-3
      h(19)=
     &-1.315673911892298936613835370593643376060412592653652307238124d-3
      h(20)=
     &9.323261308672633862226517802548514100918088299801952307991569d-04
      h(21)=
     &4.925152512628946192140957387866596210103778299388823500840094d-05
      h(22)=
     &-1.651289885565054894616687709238000755898548214659776703347801d-4
      h(23)=
     &3.067853757932549346649483228575476236600428217237900563128230d-05
      h(24)=
     &1.044193057140813708170714991080596951670706436217328169641474d-05
      h(25)=
     &-4.700416479360868325650195165061771321650383582970958556568059d-6
      h(26)=
     &5.220035098454864691736424354843176976747052155243557001531901d-07
       ENDIF

        IF (M.EQ.14) THEN
      h(1)=
     &6.461153460087947818166397448622814272327159419201199218101404d-03
      h(2)=
     &6.236475884939889832798566758434877428305333693407667164602518d-02
      h(3)=
     &2.548502677926213536659077886778286686187042416367137443780084d-01
      h(4)=
     &5.543056179408938359926831449851154844078269830951634609683997d-01
      h(5)=
     &6.311878491048567795576617135358172348623952456570017289788809d-01
      h(6)=
     &2.186706877589065214917475918217517051765774321270432059030273d-01
      h(7)=
     &-2.716885522787480414142192476181171094604882465683330814311896d-1
      h(8)=
     &-2.180335299932760447555558812702311911975240669470604752747127d-1
      h(9)=
     &1.383952138648065910739939690021573713989900463229686119059119d-01
      h(10)=
     &1.399890165844607012492943162271163440328221555614326181333683d-01
      h(11)=
     &-8.674841156816968904560822066727795382979149539517503657492964d-2
      h(12)=
     &-7.154895550404613073584145115173807990958069673129538099990913d-2
      h(13)=
     &5.523712625921604411618834060533403397913833632511672157671107d-02
      h(14)=
     &2.698140830791291697399031403215193343375766595807274233284349d-02
      h(15)=
     &-3.018535154039063518714822623489137573781575406658652624883756d-2
      h(16)=
     &-5.615049530356959133218371367691498637457297203925810387698680d-3
      h(17)=
     &1.278949326633340896157330705784079299374903861572058313481534d-02
      h(18)=
     &-7.462189892683849371817160739181780971958187988813302900435487d-4
      h(19)=
     &-3.849638868022187445786349316095551774096818508285700493058915d-3
      h(20)=
     &1.061691085606761843032566749388411173033941582147830863893939d-03
      h(21)=
     &7.080211542355278586442977697617128983471863464181595371670094d-04
      h(22)=
     &-3.868319473129544821076663398057314427328902107842165379901468d-4
      h(23)=
     &-4.177724577037259735267979539839258928389726590132730131054323d-5
      h(24)=
     &6.875504252697509603873437021628031601890370687651875279882727d-05
      h(25)=
     &-1.033720918457077394661407342594814586269272509490744850691443d-5
      h(26)=
     &-4.389704901781394115254042561367169829323085360800825718151049d-6
      h(27)=
     &1.724994675367812769885712692741798523587894709867356576910717d-06
      h(28)=
     &-1.787139968311359076334192938470839343882990309976959446994022d-7
        ENDIF

        IF (M.EQ.15) THEN
      h(1)=
     &4.538537361578898881459394910211696346663671243788786997916513d-03
      h(2)=
     &4.674339489276627189170969334843575776579151700214943513113197d-02
      h(3)=
     &2.060238639869957315398915009476307219306138505641930902702047d-01
      h(4)=
     &4.926317717081396236067757074029946372617221565130932402160160d-01
      h(5)=
     &6.458131403574243581764209120106917996432608287494046181071489d-01
      h(6)=
     &3.390025354547315276912641143835773918756769491793554669336690d-01
      h(7)=
     &-1.932041396091454287063990534321471746304090039142863827937754d-1
      h(8)=
     &-2.888825965669656462484125009822332981311435630435342594971292d-1
      h(9)=
     &6.528295284877281692283107919869574882039174285596144125965101d-02
      h(10)=
     &1.901467140071229823484893116586020517959501258174336696878156d-01
      h(11)=
     &-3.966617655579094448384366751896200668381742820683736805449745d-2
      h(12)=
     &-1.111209360372316933656710324674058608858623762165914120505657d-1
      h(13)=
     &3.387714392350768620854817844433523770864744687411265369463195d-02
      h(14)=
     &5.478055058450761268913790312581879108609415997422768564244845d-02
      h(15)=
     &-2.576700732843996258594525754269826392203641634825340138396836d-2
      h(16)=
     &-2.081005016969308167788483424677000162054657951364899040996166d-2
      h(17)=
     &1.508391802783590236329274460170322736244892823305627716233968d-02
      h(18)=
     &5.101000360407543169708860185565314724801066527344222055526631d-03
      h(19)=
     &-6.487734560315744995181683149218690816955845639388826407928967d-3
      h(20)=
     &-2.417564907616242811667225326300179605229946995814535223329411d-4
      h(21)=
     &1.943323980382211541764912332541087441011424865579531401452302d-03
      h(22)=
     &-3.734823541376169920098094213645414611387630968030256625740226d-4
      h(23)=
     &-3.595652443624688121649620075909808858194202454084090305627480d-4
      h(24)=
     &1.558964899205997479471658241227108816255567059625495915228603d-04
      h(25)=
     &2.579269915531893680925862417616855912944042368767340709160119d-05
      h(26)=
     &-2.813329626604781364755324777078478665791443876293788904267255d-5
      h(27)=
     &3.362987181737579803124845210420177472134846655864078187186304d-06
      h(28)=
     &1.811270407940577083768510912285841160577085925337507850590290d-06
      h(29)=
     &-6.316882325881664421201597299517657654166137915121195510416641d-7
      h(30)=
     &6.133359913305752029056299460289788601989190450885396512173845d-08
        ENDIF

        IF (M.EQ.16) THEN
      h(1)=
     &3.189220925347738029769547564645958687067086750131428767875878d-03
      h(2)=
     &3.490771432367334641030147224023020009218241430503984146140054d-02
      h(3)=
     &1.650642834888531178991252730561134811584835002342723240213592d-01
      h(4)=
     &4.303127228460038137403925424357684620633970478036986773924646d-01
      h(5)=
     &6.373563320837888986319852412996030536498595940814198125967751d-01
      h(6)=
     &4.402902568863569000390869163571679288527803035135272578789884d-01
      h(7)=
     &-8.975108940248964285718718077442597430659247445582660149624718d-2
      h(8)=
     &-3.270633105279177046462905675689119641757228918228812428141723d-1
      h(9)=
     &-2.791820813302827668264519595026873204339971219174736041535479d-2
      h(10)=
     &2.111906939471042887209680163268837900928491426167679439251042d-01
      h(11)=
     &2.734026375271604136485245757201617965429027819507130220231500d-02
      h(12)=
     &-1.323883055638103904500474147756493375092287817706027978798549d-1
      h(13)=
     &-6.239722752474871765674503394120025865444656311678760990761458d-3
      h(14)=
     &7.592423604427631582148498743941422461530405946100943351940313d-02
      h(15)=
     &-7.588974368857737638494890864636995796586975144990925400097160d-3
      h(16)=
     &-3.688839769173014233352666320894554314718748429706730831064068d-2
      h(17)=
     &1.029765964095596941165000580076616900528856265803662208854147d-02
      h(18)=
     &1.399376885982873102950451873670329726409840291727868988490100d-02
      h(19)=
     &-6.990014563413916670284249536517288338057856199646469078115759d-3
      h(20)=
     &-3.644279621498389932169000540933629387055333973353108668841215d-3
      h(21)=
     &3.128023381206268831661202559854678767821471906193608117450360d-03
      h(22)=
     &4.078969808497128362417470323406095782431952972310546715071397d-04
      h(23)=
     &-9.410217493595675889266453953635875407754747216734480509250273d-4
      h(24)=
     &1.142415200387223926440228099555662945839684344936472652877091d-04
      h(25)=
     &1.747872452253381803801758637660746874986024728615399897971953d-04
      h(26)=
     &-6.103596621410935835162369150522212811957259981965919143961722d-5
      h(27)=
     &-1.394566898820889345199078311998401982325273569198675335408707d-5
      h(28)=
     &1.133660866127625858758848762886536997519471068203753661757843d-05
      h(29)=
     &-1.043571342311606501525454737262615404887478930635676471546032d-6
      h(30)=
     &-7.363656785451205512099695719725563646585445545841663327433569d-7
      h(21)=
     &2.308784086857545866405412732942006121306306735866655525372544d-07
      h(32)=
     &-2.109339630100743097000572623603489906836297584591605307745349d-8
       ENDIF



        IF (M.EQ.17) THEN
      h(1)=
     &2.241807001037312853535962677074436914062191880560370733250531d-03
      h(2)=
     &2.598539370360604338914864591720788315473944524878241294399948d-02
      h(3)=
     &1.312149033078244065775506231859069960144293609259978530067004d-01
      h(4)=
     &3.703507241526411504492548190721886449477078876896803823650425d-01
      h(5)=
     &6.109966156846228181886678867679372082737093893358726291371783d-01
      h(6)=
     &5.183157640569378393254538528085968046216817197718416402439904d-01
      h(7)=
     &2.731497040329363500431250719147586480350469818964563003672942d-02
      h(8)=
     &-3.283207483639617360909665340725061767581597698151558024679130d-1
      h(9)=
     &-1.265997522158827028744679110933825505053966260104086162103728d-1
      h(10)=
     &1.973105895650109927854047044781930142551422414135646917122284d-01
      h(11)=
     &1.011354891774702721509699856433434802196622545499664876109437d-01
      h(12)=
     &-1.268156917782863110948571128662331680384792185915017065732137d-1
      h(13)=
     &-5.709141963167692728911239478651382324161160869845347053990144d-2
      h(14)=
     &8.110598665416088507965885748555429201024364190954499194020678d-02
      h(15)=
     &2.231233617810379595339136059534813756232242114093689244020869d-02
      h(16)=
     &-4.692243838926973733300897059211400507138768125498030602878439d-2
      h(17)=
     &-3.270955535819293781655360222177494452069525958061609392809275d-3
      h(18)=
     &2.273367658394627031845616244788448969906713741338339498024864d-02
      h(19)=
     &-3.042989981354637068592482637907206078633395457225096588287881d-3
      h(20)=
     &-8.602921520322854831713706413243659917926736284271730611920986d-3
      h(21)=
     &2.967996691526094872806485060008038269959463846548378995044195d-03
      h(22)=
     &2.301205242153545624302059869038423604241976680189447476064764d-03
      h(23)=
     &-1.436845304802976126222890402980384903503674530729935809561434d-3
      h(24)=
     &-3.281325194098379713954444017520115075812402442728749700195651d-4
      h(25)=
     &4.394654277686436778385677527317841632289249319738892179465910d-04
      h(26)=
     &-2.561010956654845882729891210949920221664082061531909655178413d-5
      h(27)=
     &-8.204803202453391839095482576282189866136273049636764338689593d-5
      h(28)=
     &2.318681379874595084482068205706277572106695174091895338530734d-05
      h(29)=
     &6.990600985076751273204549700855378627762758585902057964027481d-06
      h(30)=
     &-4.505942477222988194102268206378312129713572600716499944918416d-6
      h(31)=
     &3.016549609994557415605207594879939763476168705217646897702706d-07
      h(32)=
     &2.957700933316856754979905258816151367870345628924317307354639d-07
      h(33)=
     &-8.423948446002680178787071296922877068410310942222799622593133d-8
      h(34)=
     &7.267492968561608110879767441409035034158581719789791088892046d-09
      ENDIF

        IF (M.EQ.18) THEN
      h(1)=
     &1.576310218440760431540744929939777747670753710991660363684429d-03
      h(2)=
     &1.928853172414637705921391715829052419954667025288497572236714d-02
      h(3)=
     &1.035884658224235962241910491937253596470696555220241672976224d-01
      h(4)=
     &3.146789413370316990571998255652579931786706190489374509491307d-01
      h(5)=
     &5.718268077666072234818589370900623419393673743130930561295324d-01
      h(6)=
     &5.718016548886513352891119994065965025668047882818525060759395d-01
      h(7)=
     &1.472231119699281415750977271081072312557864107355701387801677d-01
      h(8)=
     &-2.936540407365587442479030994981150723935710729035053239661752d-1
      h(9)=
     &-2.164809340051429711237678625668271471437937235669492408388692d-1
      h(10)=
     &1.495339755653777893509301738913667208804816691893765610261943d-01
      h(11)=
     &1.670813127632574045149318139950134745324205646353988083152250d-01
      h(12)=
     &-9.233188415084628060429372558659459731431848000144569612074508d-2
      h(13)=
     &-1.067522466598284855932200581614984861385266404624112083917702d-1
      h(14)=
     &6.488721621190544281947577955141911463129382116634147846137149d-02
      h(15)=
     &5.705124773853688412090768846499622260596226120431038524600676d-02
      h(16)=
     &-4.452614190298232471556143559744653492971477891439833592755034d-2
      h(17)=
     &-2.373321039586000103275209582665216110197519330713490233071565d-2
      h(18)=
     &2.667070592647059029987908631672020343207895999936072813363471d-02
      h(19)=
     &6.262167954305707485236093144497882501990325204745013190268052d-03
      h(20)=
     &-1.305148094661200177277636447600807169755191054507571666606133d-2
      h(21)=
     &1.186300338581174657301741592161819084544899417452317405185615d-04
      h(22)=
     &4.943343605466738130665529516802974834299638313366477765295203d-03
      h(23)=
     &-1.118732666992497072800658855238650182318060482584970145512687d-3
      h(24)=
     &-1.340596298336106629517567228251583609823044524685986640323942d-3
      h(25)=
     &6.284656829651457125619449885420838217551022796301582874349652d-04
      h(26)=
     &2.135815619103406884039052814341926025873200325996466522543440d-04
      h(27)=
     &-1.986485523117479485798245416362489554927797880264017876139605d-4
      h(28)=
     &-1.535917123534724675069770335876717193700472427021513236587288d-7
      h(29)=
     &3.741237880740038181092208138035393952304292615793985030731363d-05
      h(30)=
     &-8.520602537446695203919254911655523022437596956226376512305917d-6
      h(31)=
     &-3.332634478885821888782452033341036827311505907796498439829337d-6
      h(32)=
     &1.768712983627615455876328730755375176412501359114058815453100d-06
      h(33)=
     &-7.691632689885176146000152878539598405817397588156525116769908d-8
      h(34)=
     &-1.176098767028231698450982356561292561347579777695396953528141d-7
      h(35)=
     &3.068835863045174800935478294933975372450179787894574492930570d-08
      h(36)=
     &-2.507934454948598267195173183147126731806317144868275819941403d-9
       ENDIF

        IF (M.EQ.19) THEN
      h(1)=
     &1.108669763181710571099154195209715164245299677773435932135455d-03
      h(2)=
     &1.428109845076439737439889152950199234745663442163665957870715d-02
      h(3)=
     &8.127811326545955065296306784901624839844979971028620366497726d-02
      h(4)=
     &2.643884317408967846748100380289426873862377807211920718417385d-01
      h(5)=
     &5.244363774646549153360575975484064626044633641048072116393160d-01
      h(6)=
     &6.017045491275378948867077135921802620536565639585963293313931d-01
      h(7)=
     &2.608949526510388292872456675310528324172673101301907739925213d-01
      h(8)=
     &-2.280913942154826463746325776054637207093787237086425909534822d-1
      h(9)=
     &-2.858386317558262418545975695028984237217356095588335149922119d-1
      h(10)=
     &7.465226970810326636763433111878819005865866149731909656365399d-02
      h(11)=
     &2.123497433062784888090608567059824197077074200878839448416908d-01
      h(12)=
     &-3.351854190230287868169388418785731506977845075238966819814032d-2
      h(13)=
     &-1.427856950387365749779602731626112812998497706152428508627562d-1
      h(14)=
     &2.758435062562866875014743520162198655374474596963423080762818d-02
      h(15)=
     &8.690675555581223248847645428808443034785208002468192759640352d-02
      h(16)=
     &-2.650123625012304089901835843676387361075068017686747808171345d-2
      h(17)=
     &-4.567422627723090805645444214295796017938935732115630050880109d-2
      h(18)=
     &2.162376740958504713032984257172372354318097067858752542571020d-02
      h(19)=
     &1.937554988917612764637094354457999814496885095875825546406963d-02
      h(20)=
     &-1.398838867853514163250401235248662521916813867453095836808366d-2
      h(21)=
     &-5.866922281012174726584493436054373773814608340808758177372765d-3
      h(22)=
     &7.040747367105243153014511207400620109401689897665383078229398d-03
      h(23)=
     &7.689543592575483559749139148673955163477947086039406129546422d-04
      h(24)=
     &-2.687551800701582003957363855070398636534038920982478290170267d-3
      h(25)=
     &3.418086534585957765651657290463808135214214848819517257794031d-04
      h(26)=
     &7.358025205054352070260481905397281875183175792779904858189494d-04
      h(27)=
     &-2.606761356786280057318315130897522790383939362073563408613547d-4
      h(28)=
     &-1.246007917341587753449784408901653990317341413341980904757592d-4
      h(29)=
     &8.711270467219922965416862388191128268412933893282083517729443d-05
      h(30)=
     &5.105950487073886053049222809934231573687367992106282669389264d-06
      h(31)=
     &-1.664017629715494454620677719899198630333675608812018108739144d-5
      h(32)=
     &3.010964316296526339695334454725943632645798938162427168851382d-06
      h(33)=
     &1.531931476691193063931832381086636031203123032723477463624141d-06
      h(34)=
     &-6.862755657769142701883554613486732854452740752771392411758418d-7
      h(35)=
     &1.447088298797844542078219863291615420551673574071367834316167d-08
      h(36)=
     &4.636937775782604223430857728210948898871748291085962296649320d-08
      h(37)=
     &-1.116402067035825816390504769142472586464975799284473682246076d-8
      h(38)=
     &8.666848838997619350323013540782124627289742190273059319122840d-10
       ENDIF

        IF (M.EQ.20) THEN
      h(1)=
     &7.799536136668463215861994818889370970510722039232863880031127d-04
      h(2)=
     &1.054939462495039832454480973015641498231961468733236691299796d-02
      h(3)=
     &6.342378045908151497587346582668785136406523315729666353643372d-02
      h(4)=
     &2.199421135513970450080335972537209392121306761010882209298252d-01
      h(5)=
     &4.726961853109016963710241465101446230757804141171727845834637d-01
      h(6)=
     &6.104932389385938201631515660084201906858628924695448898824748d-01
      h(7)=
     &3.615022987393310629195602665268631744967084723079677894136358d-01
      h(8)=
     &-1.392120880114838725806970545155530518264944915437808314813582d-1
      h(9)=
     &-3.267868004340349674031122837905370666716645587480021744425550d-1
      h(10)=
     &-1.672708830907700757517174997304297054003744303620479394006890d-2
      h(11)=
     &2.282910508199163229728429126648223086437547237250290835639880d-01
      h(12)=
     &3.985024645777120219790581076522174181104027576954427684456660d-02
      h(13)=
     &-1.554587507072679559315307870562464374359996091752285157077477d-1
      h(14)=
     &-2.471682733861358401587992299169922262915151413349313513685587d-2
      h(15)=
     &1.022917191744425578861013681016866083888381385233081516583444d-01
      h(16)=
     &5.632246857307435506953246988215209861566800664402785938591145d-03
      h(17)=
     &-6.172289962468045973318658334083283558209278762007041823250642d-2
      h(18)=
     &5.874681811811826491300679742081997167209743446956901841959711d-03
      h(19)=
     &3.229429953076958175885440860617219117564558605035979601073235d-02
      h(20)=
     &-8.789324923901561348753650366700695916503030939283830968151332d-3
      h(21)=
     &-1.381052613715192007819606423860356590496904285724730356602106d-2
      h(22)=
     &6.721627302259456835336850521405425560520025237915708362002910d-03
      h(23)=
     &4.420542387045790963058229526673514088808999478115581153468068d-03
      h(24)=
     &-3.581494259609622777556169638358238375765194248623891034940330d-3
      h(25)=
     &-8.315621728225569192482585199373230956924484221135739973390038d-4
      h(26)=
     &1.392559619323136323905254999347967283760544147397530531142397d-03
      h(27)=
     &-5.349759843997695051759716377213680036185796059087353172073952d-5
      h(28)=
     &-3.851047486992176060650288501475716463266233035937022303649838d-4
      h(29)=
     &1.015328897367029050797488785306056522529979267572003990901472d-04
      h(30)=
     &6.774280828377729558011184406727978221295796652200819839464354d-05
      h(31)=
     &-3.710586183394712864227221271216408416958225264980612822617745d-5
      h(32)=
     &-4.376143862183996810373095822528607606900620592585762190542483d-6
      h(33)=
     &7.241248287673620102843105877497181565468725757387007139555885d-06
      h(34)=
     &-1.011994010018886150340475413756849103197395069431085005709201d-6
      h(35)=
     &-6.847079597000556894163334787575159759109091330092963990364192d-7
      h(36)=
     &2.633924226270001084129057791994367121555769686616747162262697d-07
      h(37)=
     &2.014322023550512694324757845944026047904414136633776958392681d-10
      h(38)=
     &-1.814843248299695973210605258227024081458531110762083371310917d-8
      h(39)=
     &4.056127055551832766099146230616888024627380574113178257963252d-09
      h(40)=
     &-2.99883648961931956640776707
     &8372705385732460052685621923178375d-10
       ENDIF

        IF (M.EQ.21) THEN
      h(1)=
     &5.488225098526837086776336675992521426750673054588245523834775d-04
      h(2)=
     &7.776639052354783754338787398088799862510779059555623704879234d-03
      h(3)=
     &4.924777153817727491399853378340056968104483161598320693657954d-02
      h(4)=
     &1.813596254403815156260378722764624190931951510708050516519181d-01
      h(5)=
     &4.196879449393627730946850609089266339973601543036294871772653d-01
      h(6)=
     &6.015060949350038975629880664020955953066542593896126705346122d-01
      h(7)=
     &4.445904519276003403643290994523601016151342743089878478478962d-01
      h(8)=
     &-3.572291961725529045922914178005307189036762547143966578066838d-2
      h(9)=
     &-3.356640895305295094832978867114363069987575282256098351499731d-1
      h(10)=
     &-1.123970715684509813515004981340306901641824212464197973490295d-1
      h(11)=
     &2.115645276808723923846781645238468659430862736248896128529373d-01
      h(12)=
     &1.152332984396871041993434411681730428103160016594558944687967d-01
      h(13)=
     &-1.399404249325472249247758764839776903226503657502071670245304d-1
      h(14)=
     &-8.177594298086382887387303634193790542522570670234556157566786d-2
      h(15)=
     &9.660039032372422070232189700372539681627783322249829842275517d-02
      h(16)=
     &4.572340574922879239251202944731235421034828710753381191345186d-02
      h(17)=
     &-6.497750489373232063332311106008616685748929419452249544690967d-2
      h(18)=
     &-1.865385920211851534093244412008141266131208093007217139232170d-2
      h(19)=
     &3.972683542785044175197464400756126818299918992482587866999707d-02
      h(20)=
     &3.357756390338110842532604766376200760791669954106679933144723d-03
      h(21)=
     &-2.089205367797907948785235479746212371728219866525211135343707d-2
      h(22)=
     &2.403470920805434762380632169785689545910525667396313550679652d-03
      h(23)=
     &8.988824381971911875349463398395464114417817949738911101372312d-03
      h(24)=
     &-2.891334348588901247375268718015882610844675931117463495551958d-3
      h(25)=
     &-2.958374038932831280750770228215510959830170264176955719827510d-3
      h(26)=
     &1.716607040630624138494506282569230126333308533535502799235333d-03
      h(27)=
     &6.394185005120302146432543767052865436099994387647359452249347d-04
      h(28)=
     &-6.906711170821016507268939228893784790518270744313525548714065d-4
      h(29)=
     &-3.196406277680437193708834220804640347636984901270948088339102d-5
      h(30)=
     &1.936646504165080615323696689856004910579777568504218782029027d-04
      h(31)=
     &-3.635520250086338309442855006186370752206331429871136596927137d-5
      h(32)=
     &-3.499665984987447953974079490046597240276268044409625722689849d-5
      h(33)=
     &1.535482509276049283124233498646050472096482329299719141107128d-05
      h(34)=
     &2.790330539814487046106169582691767916283793946025922387556917d-06
      h(35)=
     &-3.090017164545699197158555936852697325985864588418167982685400d-6
      h(36)=
     &3.166095442367030556603889009833954440058545355777781782000278d-07
      h(37)=
     &2.992136630464852794401294607536813682771292352506328096125857d-07
      h(38)=
     &-1.000400879030597332045460600516621971679363965166249211063755d-7
      h(39)=
     &-2.254014974673330131563184851456825991617915549643308754828159d-9
      h(40)=
     &7.058033541231121859020947976903904685464512825731230495144226d-09
      h(41)=
     &-1.471954197650365265189549600816698778213247061389470277337173d-9
      h(42)=
     &1.038805571023706553035373138760372703492942617518816122570050d-10
       ENDIF

        IF (M.EQ.22) THEN
      h(1)=
     &3.862632314910982158524358900615460368877852009576899680767316d-04
      h(2)=
     &5.721854631334539120809783403484493333555361591386208129183833d-03
      h(3)=
     &3.806993723641108494769873046391825574447727068953448390456335d-02
      h(4)=
     &1.483675408901114285014404448710249837385836373168215616427030d-01
      h(5)=
     &3.677286834460374788614690818452372827430535649696462720334897d-01
      h(6)=
     &5.784327310095244271421181831735444106385099957908657145590104d-01
      h(7)=
     &5.079010906221639018391523325390716836568713192498711562711282d-01
      h(8)=
     &7.372450118363015165570139016530653113725172412104955350368114d-02
      h(9)=
     &-3.127265804282961918033226222621788537078452535993545440716988d-1
      h(10)=
     &-2.005684061048870939324361244042200174132905844868237447130382d-1
      h(11)=
     &1.640931881067664818606223226286885712554385317412228836705888d-01
      h(12)=
     &1.799731879928913037252154295313083168387840791424988422757762d-01
      h(13)=
     &-9.711079840911470969274209179691733251456735137994201552926799d-2
      h(14)=
     &-1.317681376866834107513648518146838345477875022352088357523838d-1
      h(15)=
     &6.807631439273221556739202147004580559367442550641388181886023d-02
      h(16)=
     &8.455737636682607503362813659356786494357635805197410905877078d-02
      h(17)=
     &-5.136425429744413245727949984018884707909441768477091944584584d-2
      h(18)=
     &-4.653081182750671347875833607846979997825771277976548080904423d-2
      h(19)=
     &3.697084662069802057615318892988581825637896696876361343354380d-02
      h(20)=
     &2.058670762756536044060249710676656807281671451609632981487139d-02
      h(21)=
     &-2.348000134449318868560142854519364987363882333754753819791381d-2
      h(22)=
     &-6.213782849364658499069336123807608293122900450508440420104462d-3
      h(23)=
     &1.256472521834337406887017835495604463815382993214296088172221d-02
      h(24)=
     &3.001373985076435951229129255588255746904937042979316054485183d-04
      h(25)=
     &-5.455691986156717076595353163071679107868762395367234726592273d-3
      h(26)=
     &1.044260739186025323350755659184734060807432172611689413745029d-03
      h(27)=
     &1.827010495657279080112597436850157110235336772062961041154607d-03
      h(28)=
     &-7.706909881231196232880372722955519781655769913634565757339739d-4
      h(29)=
     &-4.237873998391800799531947768003976978197438302533528661825758d-4
      h(30)=
     &3.286094142136787341983758471405935405823323072829619248523697d-04
      h(31)=
     &4.345899904532003379046992625575076092823809665933575578710696d-05
      h(32)=
     &-9.405223634815760421845190098352673647881298980040512091599943d-5
      h(33)=
     &1.137434966212593172736144274866639210339820203135670505287250d-05
      h(34)=
     &1.737375695756189356163565074505405906859746605867772002320509d-05
      h(35)=
     &-6.166729316467578372152251668422979152169587307212708981768966d-6
      h(36)=
     &-1.565179131995160159307426993578204733378112742579926503832095d-6
      h(37)=
     &1.295182057318877573889711232345068147800395721925682566394936d-06
      h(38)=
     &-8.779879873361286276888117046153049053917243760475816789226764d-8
      h(39)=
     &-1.283336228751754417819693932114064887075096030264748079976736d-7
      h(40)=
     &3.761228749337362366156711648187743399164239397803629022612862d-08
      h(41)=
     &1.680171404922988885554331183691280245962290247654438114807112d-09
      h(42)=
     &-2.729623146632976083449327361739104754443221903317745768938846d-9
      h(43)=
     &5.335938821667489905169783227036804533253011117886586305435615d-10
      h(44)=
     &-3.602113484339554703794807810939301847299106
     &970237814334104274d-11
       ENDIF

        IF (M.EQ.23) THEN
      h(1)=
     &2.719041941282888414192673609703302357098336003920923958924757d-04
      h(2)=
     &4.202748893183833538390034372523511472345215563611003407984701d-03
      h(3)=
     &2.931000365788411514736204018929480427874317460676079959515131d-02
      h(4)=
     &1.205155317839719336306053895611899089004274336891709067958035d-01
      h(5)=
     &3.184508138528652363416527748460472152790575031409830417259640d-01
      h(6)=
     &5.449311478735204282674240672421984387504149924834544495466793d-01
      h(7)=
     &5.510185172419193913452724227212507720514144116478727269717859d-01
      h(8)=
     &1.813926253638400136259098302138614937264260737638175539416540d-01
      h(9)=
     &-2.613921480306441118856795735210118413900307577511142987337375d-1
      h(10)=
     &-2.714020986078430556604069575184718123763697177381058877113471d-1
      h(11)=
     &9.212540708241805260646030910734894258577648089100630012130261d-02
      h(12)=
     &2.235736582420402317149513960822561717689875252792817094811874d-01
      h(13)=
     &-3.303744709428937875006612792463031409461636228731285046551636d-2
      h(14)=
     &-1.640113215318759250156057837165276039181451149292112929401186d-1
      h(15)=
     &2.028307457564929974897286607551313323418860610791382310375731d-02
      h(16)=
     &1.122970436181072886950734465075645977754665593869789965874572d-01
      h(17)=
     &-2.112621235622724100704783293549467048999443844657058425212982d-2
      h(18)=
     &-7.020739157490110946204219011957565343899895499962369353294028d-2
      h(19)=
     &2.176585683449997560776882472168730165799461445156766923497545d-02
      h(20)=
     &3.849533252256919901057154320407596073180564628069920893870768d-02
      h(21)=
     &-1.852351365015615979794689960740674782817814176166333519597796d-2
      h(22)=
     &-1.753710100303584537915846117408613551147985251726558719415169d-2
      h(23)=
     &1.275194393152828646243157404474947115052750581861997731041018d-02
      h(24)=
     &6.031840650024162816289878206037841640814102314209075233751820d-03
      h(25)=
     &-7.075319273706152814194039481466556204493276773483821748740018d-3
      h(26)=
     &-1.134865473356251691289337120013286756337393784110786907825400d-3
      h(27)=
     &3.122876449818144997419144765125750522437659393621577492535411d-03
      h(28)=
     &-2.465014005163512031940473100375377210862560761576109755841161d-4
      h(29)=
     &-1.061231228886651321139357625683805642193648671030425010215075d-3
      h(30)=
     &3.194204927099011503676530359692366990929679170022583007683112d-04
      h(31)=
     &2.567624520078737205563856675376636092314813400664190770435450d-04
      h(32)=
     &-1.500218503490340967673163290447832236259277810659068637402668d-4
      h(33)=
     &-3.378894834120903434270962452674534330903724108906662510305045d-5
      h(34)=
     &4.426071203109246077621875303440935335701832843654692827539837d-05
      h(35)=
     &-2.635207889249186237209225933170897825432335273771458456888097d-6
      h(36)=
     &-8.347875567854625544366043748844183086765894974439245409223337d-6
      h(37)=
     &2.397569546840240057403739507525641239509517148980849889986407d-06
      h(38)=
     &8.147574834779447778085443041422881439860288287528356019216814d-07
      h(39)=
     &-5.339005405209421154584783682848780965053642859373536945701365d-7
      h(40)=
     &1.853091785633965019353699857864654181728710556702529908304185d-08
      h(41)=
     &5.417549179539278736503176166323685597634496102979977037271945d-08
      h(42)=
     &-1.399935495437998845130909687361847103274208993447892120341999d-8
      h(43)=
     &-9.472885901812050535221582074673490573
     &092096712822067564903012d-10
      h(44)=
     &1.050446453696543404071105111096438573423068913105255997908040d-09
      h(45)=
     &-1.932405111313417542192651899622541612
     &314066389643607507706887d-10
      h(46)=
     &1.250203302351040941433216718217504240541423430995137507404787d-11
       ENDIF

        IF (M.EQ.24) THEN
      h(1)=
     &1.914358009475513695026138336474115599435172088053846745168462d-04
      h(2)=
     &3.082081714905494436206199424544404720984720556128685270556458d-03
      h(3)=
     &2.248233994971641072358415157184825628226776692231940577581580d-02
      h(4)=
     &9.726223583362519663806545734008355914527504417674578571164300d-02
      h(5)=
     &2.729089160677263268706137134412557268751671263458895098625356d-01
      h(6)=
     &5.043710408399249919771876890402814109246866444441814540282099d-01
      h(7)=
     &5.749392210955419968460807901923407033144945935105622912839838d-01
      h(8)=
     &2.809855532337118833442626085115402941842959475929278883281409d-01
      h(9)=
     &-1.872714068851562376981887159775791469060265778441667840307934d-1
      h(10)=
     &-3.179430789993627375453948489797707550898087789160025182664299d-1
      h(11)=
     &4.776613684344728187950198323031360866349104994035553200788631d-03
      h(12)=
     &2.392373887803108551973268291945824822214858134512317715815616d-01
      h(13)=
     &4.252872964148383258147364472170645232684343235486951540533893d-02
      h(14)=
     &-1.711753513703468896897638515080572393949165942335556397917666d-1
      h(15)=
     &-3.877717357792001620177594726199572688446488033750771020190283d-2
      h(16)=
     &1.210163034692242362312637311149062286659377039046006801523826d-01
      h(17)=
     &2.098011370914481534980883827326017063121637262728447783605518d-02
      h(18)=
     &-8.216165420800166702291466006164189460916816748629968198028898d-2
      h(19)=
     &-4.578436241819221637997516339765068825260159169893967894877272d-3
      h(20)=
     &5.130162003998087915555334881398688958843078494595140394873884d-02
      h(21)=
     &-4.944709428125628299815920032649550811877887219282751174798211d-3
      h(22)=
     &-2.821310709490189098113895361900699228886900995412759197674058d-2
      h(23)=
     &7.661721881646585897329899904308764405384658404613669817843430d-03
      h(24)=
     &1.304997087108573583052494067883717533043101857128653233783396d-02
      h(25)=
     &-6.291435370018187780721843581169343900864298634085743861509767d-3
      h(26)=
     &-4.746568786323113800477796959513558401732252800905982385017245d-3
      h(27)=
     &3.736046178282523345179052160810332868725126356493155728625572d-03
      h(28)=
     &1.153764936839481504858282495202271984454410046682805375157566d-03
      h(29)=
     &-1.696456818974824394274534636412116243080312601322325642741589d-3
      h(30)=
     &-4.416184856141520063365958900079406737636243682138363561877750d-5
      h(31)=
     &5.861270593183109933716735450272894035425792347806515678695765d-04
      h(32)=
     &-1.181233237969554740613021227756568966806892308457221016257961d-4
      h(33)=
     &-1.460079817762616838924301818082729036314539476811023255670666d-4
      h(34)=
     &6.559388639305634085303738560455061974369354538271316071502698d-05
      h(35)=
     &2.183241460466558363365044032984257709791187640963509380549307d-05
      h(36)=
     &-2.022888292612697682860859987200455702614855595412267510558659d-5
      h(37)=
     &1.341157750809114719319937553186023660581084151828593222893663d-08
      h(38)=
     &3.901100338597702610409014129024223853127911530009766793352492d-06
      h(39)=
     &-8.980253143938407724149926669980791166378388013293887718404796d-7
      h(40)=
     &-4.032507756879971624098983247358983425236092110387724315244646d-7
      h(41)=
     &2.166339653278574639176393978510246335478946697396400359281412d-07
      h(42)=
     &-5.057645419792500308492508924343248979
     &317507866520688417567606d-10
      h(43)=
     &-2.255740388176086107368821674947175804005323153443170526520277d-8
      h(44)=
     &5.157776789671999638950774266313208715015419699643333784626363d-09
      h(45)=
     &4.748375824256231118094453549799175824526559994333227456737433d-10
      h(46)=
     &-4.024658644584379774251499574468195118
     &601698713554294941756559d-10
      h(47)=
     &6.991801157638230974132696433509625934021677793453732225542951d-11
      h(48)=
     &-4.3427825038037102472590375528867494579
     &51053124203814185811297d-12
       ENDIF

        IF (M.EQ.25) THEN
      h(1)=
     &1.348029793470188994578489247159356055370460656508881471268611d-04
      h(2)=
     &2.256959591854779520121391049628056149270016860666661928130747d-03
      h(3)=
     &1.718674125404015533817186914954848902241194002444696221013131d-02
      h(4)=
     &7.803586287213267559750659320481403668422052199257139168386084d-02
      h(5)=
     &2.316935078860218199900621518057089104946216881512075361624214d-01
      h(6)=
     &4.596834151460945937896973864539659944010260858049947396093277d-01
      h(7)=
     &5.816368967460577833534892038757085635755639698734580573323031d-01
      h(8)=
     &3.678850748029466984371319740855532278670733841012809062966976d-01
      h(9)=
     &-9.717464096463814276130048169040892607068486428294030952842447d-2
      h(10)=
     &-3.364730796417461309562110148848845218930261030262170601615289d-1
      h(11)=
     &-8.758761458765466140226687673880006154266689569025041229545538d-2
      h(12)=
     &2.245378197451017129525176510409543155930843160711989062118482d-01
      h(13)=
     &1.181552867199598604563067876819931882639429216001523151773895d-01
      h(14)=
     &-1.505602137505796309518094206831433270850173484773520730186277d-1
      h(15)=
     &-9.850861528996022153725952822686729410420350758543226219234795d-2
      h(16)=
     &1.066338050184779528831274540522414711301747903916268438037723d-01
      h(17)=
     &6.675216449401860666895983072443984697329752470942906490126865d-02
      h(18)=
     &-7.708411105657419356208567671699032054872853174701595359329826d-2
      h(19)=
     &-3.717396286112250887598137324046870459877639250821705817221557d-2
      h(20)=
     &5.361790939877949960629041419546536897037332284703545849594129d-02
      h(21)=
     &1.554260592910229163981295854603203625062268043511894295387375d-02
      h(22)=
     &-3.404232046065334099320628584033729153497903561399447916116575d-2
      h(23)=
     &-3.079836794847036661636693963570288706232460663070983852354326d-3
      h(24)=
     &1.892280447662762841086581178691039363674755753459524525886478d-02
      h(25)=
     &-1.989425782202736494289461896386235348901617760816745484282494d-3
      h(26)=
     &-8.860702618046368399013064252456556969199612331833605310278698d-3
      h(27)=
     &2.726936258738495739871469244610042793734119359765762028996059d-03
      h(28)=
     &3.322707773973191780118197357194829286271392998979276105842863d-03
      h(29)=
     &-1.842484290203331280837780430014195744813667655929909114672154d-3
      h(30)=
     &-8.999774237462950491085382524008429604309720852269895692000702d-4
      h(31)=
     &8.772581936748274843488806190175921376284150686011179612908221d-04
      h(32)=
     &1.153212440466300456460181455345639872216326644527860903202733d-04
      h(33)=
     &-3.098800990984697989530544245356271119416614147098459162436317d-4
      h(34)=
     &3.543714523276059005284289830559259809540337561365927850248007d-05
      h(35)=
     &7.904640003965528255137496303166001735463107762646364003487560d-05
      h(36)=
     &-2.733048119960041746353244004225286857636045649642652816856524d-5
      h(37)=
     &-1.277195293199783804144903848434605690990373526086311486716394d-5
      h(38)=
     &8.990661393062588905369930197413951232059323587543226269327396d-06
      h(39)=
     &5.232827708153076417963912065899772684403904504491727061662335d-07
      h(40)=
     &-1.779201332653634562565948556039009149458987774189389221295909d-6
      h(41)=
     &3.212037518862519094895005816661093988294166712919881121802831d-07
      h(42)=
     &1.922806790142371601278104244711267420759978799176017569693322d-07
      h(43)=
     &-8.656941732278507163388031517930974947984281611717187862530250d-8
      h(44)=
     &-2.611598556111770864259843089151782206922842627174274274741722d-9
      h(45)=
     &9.279224480081372372250073354726511359667401736947170444723772d-09
      h(46)=
     &-1.880415755062155537197782595740975189878162661203102565611681d-9
      h(47)=
     &-2.22847491022816889931479335206479595
     &7306403503495743572518755d-10
      h(48)=
     &1.535901570162657197021927739530721955859277615795931442682785d-10
      h(49)=
     &-2.5276251634656448110488642861697581
     &28142169484216932624854015d-11
      h(50)=
     &1.509692082823910867903367712096001664979004526477422347957324d-12
       ENDIF

        IF (M.EQ.26) THEN
      h(1)=
     &9.493795750710592117802731381148054398461637804818126397577999d-05
      h(2)=
     &1.650520233532988247022384885622071050555268137055829216839523d-03
      h(3)=
     &1.309755429255850082057770240106799154079932963479202407364818d-02
      h(4)=
     &6.227474402514960484193581705107415937690538641013309745983962d-02
      h(5)=
     &1.950394387167700994245891508369324694703820522489789125908612d-01
      h(6)=
     &4.132929622783563686116108686666547082846741228042232731476147d-01
      h(7)=
     &5.736690430342222603195557147853022060758392664086633396520345d-01
      h(8)=
     &4.391583117891662321931477565794105633815363384084590559889493d-01
      h(9)=
     &1.774076780986685727823533562031556893226571319881417676492595d-03
      h(10)=
     &-3.263845936917800216385340830055349953447745005769416287177497d-1
      h(11)=
     &-1.748399612893925042664835683606584215248582345438816346170042d-1
      h(12)=
     &1.812918323111226960705459766025430918716233584167982942044424d-01
      h(13)=
     &1.827554095896723746537533832033286839689931924709760567945595d-01
      h(14)=
     &-1.043239002859270439148009137202747658420968144330108510179290d-1
      h(15)=
     &-1.479771932752544935782314546369458188243947772922980064071205d-1
      h(16)=
     &6.982318611329236513756591683950208955110603212379412334701145d-02
      h(17)=
     &1.064824052498086303236593797715344405836015002929319291715777d-01
      h(18)=
     &-5.344856168148319149493577269390074213960237013099439431132086d-2
      h(19)=
     &-6.865475960403591525454725258715351280947435823354011140858001d-2
      h(20)=
     &4.223218579637203541206570902753288247790857760067894456114927d-02
      h(21)=
     &3.853571597111186425832144567362328142994885395255438867968781d-02
      h(22)=
     &-3.137811036306775484244644776337594435094096964336402798072360d-2
      h(23)=
     &-1.776090356835818354094298625884058170354129044259951019182732d-2
      h(24)=
     &2.073492017996382475887790073068984224515077665517103399898854d-02
      h(25)=
     &5.829580555318887971939315747596613038479561943085291072787359d-03
      h(26)=
     &-1.178549790619302893728624468402138072504226527540325463847390d-2
      h(27)=
     &-5.287383992626814439198630765217969804966319971038003993984480d-4
      h(28)=
     &5.601947239423804853206514239940474788977188460452053462770324d-03
      h(29)=
     &-9.390582504738289646165698675070641765810790863514339205205998d-4
      h(30)=
     &-2.145530281567620980305401403432221668847980295600748913748902d-3
      h(31)=
     &8.383488056543616046381924054554052104937784379435436426690560d-04
      h(32)=
     &6.161382204574344193703789012696411561214682388271673214197731d-04
      h(33)=
     &-4.319557074261807466712901913481943478521991611607433971794602d-4
      h(34)=
     &-1.060574748283803889966150803551837402553866816191659959347053d-4
      h(35)=
     &1.574795238607493590547765666590811258087715699737771458390360d-04
      h(36)=
     &-5.277795493037868976293566636015627609248847457646525246271036d-6
      h(37)=
     &-4.109673996391477816326502438997466532822639385119090230965252d-5
      h(38)=
     &1.074221540872195031273584409245060623104931330938273936484593d-05
      h(39)=
     &7.000078682964986734859102495210684809643657474253921074934684d-06
      h(40)=
     &-3.887400161856795187587790410706550576033603097954065074023128d-6
      h(41)=
     &-4.650463220640262639231145944536092973446596027469833860001618d-7
      h(42)=
     &7.939210633709952088373459255067360793370284788682979065122810d-07
      h(43)=
     &-1.079004237578671411922961583845716126060658213943840375162654d-7
      h(44)=
     &-8.904466370168590769052983362721567202750591914741016835071257d-8
      h(45)=
     &3.407795621290730008673832107214820587991557116806912418558069d-08
      h(46)=
     &2.169328259850323106986222296525930099935873861026310788086221d-09
      h(47)=
     &-3.776010478532324328184043667556576385639846460337894963138621d-9
      h(48)=
     &6.780047245828636668305808192607091517605349478677442468580825d-10
      h(49)=
     &1.002303191046526913509281844136258004034177309673269533418644d-10
      h(50)=
     &-5.8404081853411714684654924477998192629
     &05317576847426970757700d-11
      h(51)=
     &9.130510016371796243923232926650252570239054815939483900056681d-12
      h(52)=
     &-5.251871224244435037810503452564279
     &828539007071678724285717464d-13
       ENDIF

        IF (M.EQ.27) THEN
      h(1)=
     &6.687131385431931734918880680779563307675740731544063787599480d-05
      h(2)=
     &1.205531231673213234251999812212394463872002561229330125152073d-03
      h(3)=
     &9.952588780876619771874091297340545740163119816300838847749336d-03
      h(4)=
     &4.945259998290488004302995584228917712171023349013386944893643d-02
      h(5)=
     &1.629220275023933206396286389387812803673796872000118325233533d-01
      h(6)=
     &3.671102141253898226423388094379126394383458407087000700420400d-01
      h(7)=
     &5.538498609904800487605460395549044755068663194750017660900436d-01
      h(8)=
     &4.934061226779989979265447084358038959373468583404767251300717d-01
      h(9)=
     &1.028408550618229112710739475157388764479351682549490307668477d-01
      h(10)=
     &-2.897168033145948463175311101489473923261698802610323264603418d-1
      h(11)=
     &-2.482645819032605667810198368127693701263349361209208170092197d-1
      h(12)=
     &1.148230195177853576326445213787661879970642975306605349249036d-01
      h(13)=
     &2.272732884141708265275037216925482827043581894357907763081103d-01
      h(14)=
     &-3.878641863180231062443346843661817078060143110529946543683356d-2
      h(15)=
     &-1.780317409590085821070366277249759321269342801053489323888575d-1
      h(16)=
     &1.579939746024048431173907799261019471878724997312653292884660d-02
      h(17)=
     &1.311979717171553289711406975836688896451835867594492827800969d-01
      h(18)=
     &-1.406275155580876537026622167053147161846397735962817855782362d-2
      h(19)=
     &-9.102290652956591798241345515773322449830692586525337562864481d-2
      h(20)=
     &1.731101826549371089085675445961947677452358872325373949295769d-02
      h(21)=
     &5.796940573471798814748840657698008349462526768238833307489106d-02
      h(22)=
     &-1.851249356199807710545837861298826718763077900221574749342712d-2
      h(23)=
     &-3.273906663102087145481936428049519742538150452785563039743756d-2
      h(24)=
     &1.614696692239566682272152627542980896527822528487665111124260d-02
      h(25)=
     &1.566559564892457873003263983940819950829497022298967052103291d-02
      h(26)=
     &-1.157718645897628140054089958116866381056430680879332334217267d-2
      h(27)=
     &-5.862096345462925972966025215266179082657169806555503857975278d-3
      h(28)=
     &6.856635609684880675273184141746359000591385833807880272568038d-03
      h(29)=
     &1.342626877303679609082208800217479591902967766971379107017011d-03
      h(30)=
     &-3.332854469520006162763300141047111065412307706449049389557931d-3
      h(31)=
     &1.457529625931728587128588244152604734177322144376309490881599d-04
      h(32)=
     &1.301177450244135139135787970279897042994109161268159963884641d-03
      h(33)=
     &-3.418351226915427611946547437228006377896519777431057005796358d-4
      h(34)=
     &-3.879018574101327604369144470124819695479087900682219330965466d-4
      h(35)=
     &2.019719879690326857104208791272390315160018069955787875123234d-04
      h(36)=
     &7.660058387068576876674274961751262847965101108848090019821555d-05
      h(37)=
     &-7.711145517797584208411720507329584053382646435270054267102827d-5
      h(38)=
     &-3.517483614907445391752737841583832374184046409747387149129674d-6
      h(39)=
     &2.063442647736885318487206413360228908558806028468062177953960d-05
      h(40)=
     &-3.901164070638425528170558032557368703418425915665413541985623d-6
      h(41)=
     &-3.657500908187104997045760131046655906827644494899206692043298d-6
      h(42)=
     &1.634369624725637835424610743915128591988676092276368687669255d-06
      h(43)=
     &3.050880686251999094242671997731089918322345713516567387655763d-07
      h(44)=
     &-3.472468147394389269364673179891460601330730511237974736379548d-7
      h(45)=
     &3.286558968055159530983261866450459360074591641809187825408848d-08
      h(46)=
     &4.026255052866908637178682747490340533992340623231336911661711d-08
      h(47)=
     &-1.321332273990056558848617809101876846857728483295631388083263d-8
      h(48)=
     &-1.309465606856955151282041809232358209226373823424148862843577d-9
      h(49)=
     &1.521614984778521740775073159445241799352681846880808663329946d-09
      h(50)=
     &-2.415526928011130660506395791946234
     &018673860470542996426005750d-10
      h(51)=
     &-4.37498622429365439506994768201399635
     &1823060759948583134078918d-11
      h(52)=
     &2.213662088067662485181472969374945928903854605356443772873438d-11
      h(53)=
     &-3.29579012247658580706995397504309613
     &9541415768606924980926275d-12
      h(54)=
     &1.828188352882424933624530026056448539377272017834175009418822d-13
       ENDIF

        IF (M.EQ.28) THEN
      h(1)=
     &4.710807775014051101066545468288837625869263629358873937759173d-05
      h(2)=
     &8.794985159843870273564636742144073059158975665525081816488582d-04
      h(3)=
     &7.542650377646859177160195786201116927568410621050693986450538d-03
      h(4)=
     &3.909260811540534426092083794403768111329778710541126982205076d-02
      h(5)=
     &1.351379142536410450770749411679708279921694061092200363031937d-01
      h(6)=
     &3.225633612855224257318486139030596702170126503618082416187649d-01
      h(7)=
     &5.249982316303355562348293243640252929543774162151269406404636d-01
      h(8)=
     &5.305162934414858075256978195354516449402692654391295761050628d-01
      h(9)=
     &2.001761440459844380384404537971725815970574972480152145882083d-01
      h(10)=
     &-2.304989540475825257279397658067038304888129374484095837624889d-1
      h(11)=
     &-3.013278095326417816909366061441334075444383937588485826752087d-1
      h(12)=
     &3.285787916338710468450547883547348694255260871071954509422161d-02
      h(13)=
     &2.458081513737595535752949960866466132239832334168533456626848d-01
      h(14)=
     &3.690688531571127205290633425993077868843846977265847006108551d-02
      h(15)=
     &-1.828773307329849166920408764650763092868965221608724574218473d-1
      h(16)=
     &-4.683823374455167616514752420549419665215987106243491879971921d-2
      h(17)=
     &1.346275679102260877490923315484152662987698625205479167761416d-01
      h(18)=
     &3.447863127509970524678534595639646616244376966117385829345554d-02
      h(19)=
     &-9.768535580565244174963692133038973587005628990493154911133358d-2
      h(20)=
     &-1.734192283130589908795581592406238282930530566316914040035812d-2
      h(21)=
     &6.774789550190933956165341752699717255041141690153626336867769d-02
      h(22)=
     &3.448018955540951137600471926079622335842207388713342609755316d-03
      h(23)=
     &-4.333336861608628393863254980828284403766309203453808666888800d-2
      h(24)=
     &4.431732910062988320487418656322338284504389482966303454010563d-03    
      h(25)=
     &2.468806001015186586264188361362046240243934625858343309818244d-02
      h(26)=
     &-6.815549764552309639259447104811254179605050667281644254737890d-3
      h(27)=
     &-1.206359196821849005842466619530619474644989878503490321948471d-2
      h(28)=
     &5.838816627748944864497370576838809711476027837762897602935327d-03
      h(29)=
     &4.784863112454241718009916669120329848973107781600157214960003d-03
      h(30)=
     &-3.725461247074254799171427871442937099025589672466088044410521d-3
      h(31)=
     &-1.360373845639692436577650137133777929659265166644839235882291d-3
      h(32)=
     &1.875998668202795626152766912508562385106168761893900192731562d-03
      h(33)=
     &1.415672393140464257573780581396205840941849282748250523509874d-04
      h(34)=
     &-7.486749559114629991320679819683227355746847370960399216568306d-4
      h(35)=
     &1.154656063658921251969297916771881248142872975490882572741198d-04
      h(36)=
     &2.295790982233456202366621544054366855729175050420515776344878d-04
      h(37)=
     &-8.903901490044488099517361247378396756893227855233897357882978d-5
      h(38)=
     &-4.907713416190250858324783990436748073854807494400738311968278d-5
      h(39)=
     &3.641401211050802781223450761733180188911730291497201507086247d-05
      h(40)=
     &4.638664981394294654002871426476885751050837817671843706915388d-06
      h(41)=
     &-1.004326041333422601781848560432120920634648692782357855473103d-5
      h(42)=
     &1.247900317574834146052381692752796047052443265982232422642017d-06
      h(43)=
     &1.840363734517769191684379309039277810350620305330900536404818d-06
      h(44)=
     &-6.670215479954892588747450458085225880096882699397256774967304d-7
      h(45)=
     &-1.757461173209842779903676264971918635870906983281392939812547d-7
      h(46)=
     &1.490660013535362170989340065033061951960933954388633507264360d-07
      h(47)=
     &-8.262387315626556965966429243600984899650039704831080988658278d-9
      h(48)=
     &-1.784138690875710077191713941441263246560738410213624546116655d-8
      h(49)=
     &5.044047056383436444631252840057862002264087720676808580373667d-09
      h(50)=
     &6.944540328946226952976704718677697525410051405055662575530111d-10
      h(51)=
     &-6.07704124722901022476024530559630780383
     &0053533836849384680534d-10
      h(52)=
     &8.492220011056382105461206077240377024404404638947591299761197d-11
      h(53)=
     &1.867367263783390418963879146175452376940453585791428841004699d-11
      h(54)=
     &-8.365490471258800799349289794397908900767
     &054085216008197372193d-12
      h(55)=
     &1.188850533405901520842321749021089497203940688882364518455403d-12
      h(56)=
     &-6.36777235471485733563269209226725426636
     &8934590973693820942617d-14
      ENDIF

        IF (M.EQ.29) THEN
      h(1)=
     &3.318966279841524761813546359818075441349169975922439988843475d-05
      h(2)=
     &6.409516803044434540833706729120596322083061716935004987374676d-04
      h(3)=
     &5.702126517773375434760843998623507494914551464968126455168657d-03
      h(4)=
     &3.077358022140837676716707336516751814713312018344719150923618d-02
      h(5)=
     &1.113701169517405304762186166370327770191325772342190715118617d-01
      h(6)=
     &2.806534559709829376968881262770480606500920092398534229615289d-01
      h(7)=
     &4.897588047621993143592705932993573539235839610055331620240518d-01
      h(8)=
     &5.513744327583751951223746071670135992466984391233429663886536d-01
      h(9)=
     &2.891052383358291634605691113586264061513180158354460952469246d-01
      h(10)=
     &-1.540287344599000542466293779503370141731339982919280951230240d-1
      h(11)=
     &-3.300409489175880520295083779487012611959310539629627124613719d-1
      h(12)=
     &-5.570680007294085781514541931715795784309410235726214400350351d-2
      h(13)=
     &2.361052361530259415983110734054626770649468357328362426830433d-01
      h(14)=
     &1.124191748731883764769740670535880543076817816861518667898467d-01
      h(15)=
     &-1.608779885941877360771615465531852333085159940159968393590303d-1
      h(16)=
     &-1.078459499387214201077881957354707913786241153934264316589273d-1
      h(17)=
     &1.144722958938182579734135930060053286267822797640393386903440d-01
      h(18)=
     &8.322074716244975790297348835032537357891920536002627784941129d-02
      h(19)=
     &-8.512549261563550232832311331420804581881235448862834507281486d-2
      h(20)=
     &-5.502748952532572320924541450626650067707344725344841099873446d-2
      h(21)=
     &6.347916458421186633577789314698972361081611994794140119302163d-02
      h(22)=
     &3.053154327270413646637328212093941030592133225231728964047047d-02
      h(23)=
     &-4.518798127778834515979704475304405691390090327474972089790857d-2
      h(24)=
     &-1.291714255426679462966473962555410660387671182428076570686472d-2
      h(25)=
     &2.947043187174764111028122319949903667638786379520519899154373d-02
      h(26)=
     &2.648327307678167915542397563479749119673768286990136051577167d-03
      h(27)=
     &-1.704122457360668969234196743407615179099529206118693044741086d-2
      h(28)=
     &1.737880332720511164430027824345354801611373419264590068097416d-03
      h(29)=
     &8.469725493560752287772961661104710791306496373354237126998903d-03
      h(30)=
     &-2.550807127789472659145072247724735637183590942511858255354005d-3
      h(31)=
     &-3.473798989681100630649790255076233970957721666820195620598374d-3
      h(32)=
     &1.877120925723650133179338154344873477230567340668548016358682d-03
      h(33)=
     &1.087053942226062966738944397844498417945523630053411148182206d-03
      h(34)=
     &-1.000778327085680541055696707760062870925897014530348262794137d-3
      h(35)=
     &-2.000711363076779808296301110796026470163110202848894744316755d-4
      h(36)=
     &4.111283454742767033424740543004041500054889660665367490129376d-04
      h(37)=
     &-2.292018041214499897382298271438084577065170236103859181134525d-5
      h(38)=
     &-1.293044840080720609161466939678226852440475312744714379499074d-4
      h(39)=
     &3.645026068562774967665464216602750761690984830805534178557146d-05
      h(40)=
     &2.913344750169041218495787251929571015775436967652945386217480d-05
      h(41)=
     &-1.657328395306616289863396387854880512976861409870690029695161d-5
      h(42)=
     &-3.593644804025187638066915189731950450034629392522542962477168d-6
      h(43)=
     &4.750609246452552850197117564759363194953518317428400241629683d-06
      h(44)=
     &-3.029054592052818286474228294307141792053791695855058563299597d-7
      h(45)=
     &-8.975701750636280734511651941681818767895052287332471537510510d-7
      h(46)=
     &2.633898386997696553900967704111473475368019612368922599394214d-07
      h(47)=
     &9.387197411095863026484410601284876812292554863800653292318725d-08
      h(48)=
     &-6.286156922010786166768503252870590953166867739448102804392389d-8
      h(49)=
     &1.076591906619196137385201975028785139607670319821266803566785d-09
      h(50)=
     &7.768978854770062238895964639391324551611701293594055935346266d-09
      h(51)=
     &-1.893995386171984147774611076618946011337498790609031626697228d-9
      h(52)=
     &-3.426800863263089001811012278889864200550
     &342566386405676893537d-10
      h(53)=
     &2.407099453509342962399811991929330725186626582891090462239366d-10
      h(54)=
     &-2.94058925076453258288847397463827366424
     &4682541297835986306504d-11
      h(55)=
     &-7.8325097336278170323565565828194947948
     &84131433810848844709881d-12
      h(56)=
     &3.152762413370310423797539876893861621418382024668704492620948d-12
      h(57)=
     &-4.2856548700683441018981850733763
     &07686875386259541180967347399d-13
      h(58)=
     &2.219191311588302960934661700068023727737812918006011019184982d-14
       ENDIF

        IF (M.EQ.30) THEN
      h(1)=
     &2.338616172731421471474407279894891960011661146356580425400538d-05
      h(2)=
     &4.666379504285509336662000111055365140848987563882199035322085d-04
      h(3)=
     &4.300797165048069510045016757402827408493482974782286966500398d-03
      h(4)=
     &2.413083267158837895194919987958311943976725005113561262334092d-02
      h(5)=
     &9.123830406701570679321575555085899708564500191080751595642650d-02
      h(6)=
     &2.420206709402140994467599658342919512318194032687898436229538d-01
      h(7)=
     &4.504878218533178366981351802898336415314944375740699506554771d-01
      h(8)=
     &5.575722329128364304078082520999850413492571645754785374629734d-01
      h(9)=
     &3.662426833716279793144871151369089533016299234992584741629624d-01
      h(10)=
     &-6.618367077593731501909741041813726474911212544474895441395148d-2
      h(11)=
     &-3.329669750208556069196849320598850505877494561268613506392514d-1
      h(12)=
     &-1.419685133300829310219026267403758254954270602825020111483505d-1
      h(13)=
     &1.994621215806643032428990062111230223523226088131364328774921d-01
      h(14)=
     &1.778298732448367361280250921330425046260289700971176750362566d-01
      h(15)=
     &-1.145582194327077814891518778613672243404957549114393749173137d-1
      h(16)=
     &-1.572368179599938126878197378886501553251711910617673398124611d-1
      h(17)=
     &7.277865897036442699893544326605244235248713804556715604416632d-02
      h(18)=
     &1.227477460450093778691578797698150091624353365248212907325446d-01
      h(19)=
     &-5.380646545825707676022015051837304300338645984615639237930800d-2
      h(20)=
     &-8.765869003638366048026572053699028353846982304851342479893827d-2
      h(21)=
     &4.380166467141773250305407710250135373016604593736480428415303d-02
      h(22)=
     &5.671236574473569492590636983030617493807140224924978946302257d-02
      h(23)=
     &-3.567339749675960965780819743176056734137251336781389369397564d-2
      h(24)=
     &-3.226375891935220815954913483392725682165778426411705216010280d-2
      h(25)=
     &2.707861959529418272206848318420006522973840949600186710327776d-02
      h(26)=
     &1.528796076985739546052896626042375110302102640936712142026221d-02
      h(27)=
     &-1.839974386811734118728169880549148389603890445324127330811811d-2
      h(28)=
     &-5.296859666131086629169938675330494864053932988161015674773617d-3
      h(29)=
     &1.091563165830488927536881480211929049886878831313700460017968d-02
      h(30)=
     &6.196717564977244383592534999284255315694546230739551683085460d-04
      h(31)=
     &-5.530730148192003288871383856487027893918513053091795443517653d-3
      h(32)=
     &8.433845866620933982126003584365932145598126087481400294999080d-04
      h(33)=
     &2.324520094060099304385756339638431339131122661576649123053845d-03
      h(34)=
     &-8.609276968110423879660725173525347077801305237644122054954659d-4
      h(35)=
     &-7.678782504380918697963922441514742758516706160788123977340073d-4
      h(36)=
     &5.050948239033467796256544554086554367969638627715114003635557d-04
      h(37)=
     &1.724825842351709725545759714374272164367933578194910678479473d-04
      h(38)=
     &-2.161718301169633804271038862087964094429005266172702380483361d-4
      h(39)=
     &-8.548305467584070994787824796256108217987765582429940610377190d-6
      h(40)=
     &6.982008370808327851082027193100914402221658444151889697045071d-05
      h(41)=
     &-1.339716863293971629296314599448901465078920406443516550195793d-5
      h(42)=
     &-1.636152478725426488654528710478856195004608401773950511915162d-5
      h(43)=
     &7.252145535890469015723401169934327900622894130695550273452916d-06
      h(44)=
     &2.327549098493686509557358103785598216688723737824121617676858d-06
      h(45)=
     &-2.187267676996166416699555236143059249832615777542412142603694d-6
      h(46)=
     &1.099474338526203304286307383463498542376432972308342428764576d-08
      h(47)=
     &4.261662326011572446469849114416378817419458434583398455985144d-07
      h(48)=
     &-1.000414682354500898864979332965559934104686157639553850670490d-7
      h(49)=
     &-4.764379965139453357729154748688006975561934425368712852985388d-8
      h(50)=
     &2.605442754977625431940885841950955928085338672381046225838880d-08
      h(51)=
     &5.553397861397053982967618072672572206490972606026556946910028d-10
      h(52)=
     &-3.331105680467578245901976412732595596538702049437802824373020d-9
      h(53)=
     &6.984862691832182584221096665570313611280449991512869846064780d-10
      h(54)=
     &1.613622978270904360610418704685783656905979134344922647926295d-10
      h(55)=
     &-9.461387997276802120884525814092001871
     &993910062127702293573920d-11
      h(56)=
     &1.000105131393171192746337860330428369495110180346654025287492d-11
      h(57)=
     &3.239428638532286114355931428908079297696045600279108835760520d-12
      h(58)=
     &-1.18523759210158232825423149631058461
     &1948560976394420324137742d-12
      h(59)=
     &1.543997570847620046003616417646988780670333040868954794039905d-13
      h(60)=
     &-7.737942630954405708679963277418806436
     &871098329050829841696327d-15
       ENDIF

        IF (M.EQ.31) THEN
      h(1)=
     &1.648013386456140748122177817418358316441195236228590958603489d-05
      h(2)=
     &3.394122037769956699157160165352942212213928231154233571163033d-04
      h(3)=
     &3.236884068627721221829662672296912258338131668810067169630813d-03
      h(4)=
     &1.885369161298591269159568944275763468999829139547989648553486d-02
      h(5)=
     &7.433609301164788697908776495388047669378919816041031344650271d-02
      h(6)=
     &2.070128744852353286198055444111916450619762837756134323019573d-01
      h(7)=
     &4.091922000374278563928213235836188963704298775635493549519369d-01
      h(8)=
     &5.511398409142754983590484577074663132074992263886810324421617d-01
      h(9)=
     &4.294688082061372955430413148799008354573408538414331312236645d-01
      h(10)=
     &2.716921249736946422305354732634261873401679092095992827198308d-02
      h(11)=
     &-3.109551183195075186926560285811004715398678229333522634202008d-1
      h(12)=
     &-2.179784855235633521693544507220105631639547435903112747133934d-1
      h(13)=
     &1.401782887652732681656253206993073895422881511380152633441096d-01
      h(14)=
     &2.249667114737370933697297905066886078307490136415302624018330d-01
      h(15)=
     &-4.992634916046823977000579399730138693074543903234092797936484d-2
      h(16)=
     &-1.869623608957154494374577196258383009208655076187653847079167d-1
      h(17)=
     &1.543698842948893409652995335281236231845293548571166883219023d-02
      h(18)=
     &1.450895009319931981518942907854879059128872873116921504156674d-01
      h(19)=
     &-8.139832273469236863527708715566588550006680549152344840146851d-3
      h(20)=
     &-1.076127733234956326668605511648013952380301953590447106075614d-1
      h(21)=
     &1.094129745236496925725237900637802669504835743555466811796369d-02
      h(22)=
     &7.535361174328140695528289751109133941376701984419452638686226d-02
      h(23)=
     &-1.488002661810482202699555987503429289100801979910046913257306d-2
      h(24)=
     &-4.861907546485433003537603385831190109391263542044516048871113d-2
      h(25)=
     &1.615417156598591113619453864586701665635869166193865651960591d-02
      h(26)=
     &2.804761936675616906861927211659154977049392281479113764697785d-02
      h(27)=
     &-1.427627527776351943309800140756746087215016194775579070599004d-2
      h(28)=
     &-1.390055293926652880755898888934447671732373519028670201124816d-2
      h(29)=
     &1.051763948737184089128633441244991643331033825102031908858652d-02
      h(30)=
     &5.516163573310992566561289762241160214476622662764637181816550d-03
      h(31)=
     &-6.520852375874612553325469682628530079210293774541131381751695d-3
      h(32)=
     &-1.428264223218909891400516038687842292177211292295049238921068d-3
      h(33)=
     &3.393066776715931928419358796960612411097347419792355896915546d-03
      h(34)=
     &-6.397901106014600492881202314307290077992972755016494062875201d-5
      h(35)=
     &-1.459041741985160943114515221598080223845239255190055621901681d-3
      h(36)=
     &3.431398296904734438118401084929505912208229684629857530009147d-04
      h(37)=
     &4.998816175637222614896912406679513231966722440032799024979502d-04
      h(38)=
     &-2.396583469402949615285646688069476140260781708006174912535660d-4
      h(39)=
     &-1.243411617250228669409179807383399199879641177993453588807726d-4
      h(40)=
     &1.089584350416766882738651833752634206358441308880869184416670d-04
      h(41)=
     &1.501335727444532997071651937630983442758297688087711521441229d-05
      h(42)=
     &-3.631255157860086164261313773172162991107348698083164489165837d-5
      h(43)=
     &4.034520235184278839752741499546098778993926344831736074409765d-06
      h(44)=
     &8.795301342692987765440618030678349427367022581211855857458220d-06
      h(45)=
     &-3.035142365891509630069007852947057220760887215249503512783023d-6
      h(46)=
     &-1.369060230942940782050489751987123955074404782177163471279285d-6
      h(47)=
     &9.810015422044371573950976088058064384946146188110905321673802d-07
      h(48)=
     &5.327250656974915426977440959783080593776012130063170688309127d-08
      h(49)=
     &-1.975925129170206248152121156696590501303803187231928513867046d-7
      h(50)=
     &3.616826517331004805247567218405798591329788122337274956172315d-08
      h(51)=
     &2.328309713821409644308538888589329921141948539678106680777082d-08
      h(52)=
     &-1.061529602150252306500404266150823962402673780484965538270541d-8
      h(53)=
     &-6.47431168795986139870258153934195
     &4438747926255671605657095807d-10
      h(54)=
     &1.408568151025177427076547804944585301332087108125727813194374d-09
      h(55)=
     &-2.5240439541533533061836437029332183
     &08617979467184848456565837d-10
      h(56)=
     &-7.34893003248626390476691391965362437
     &9586487437915175106407348d-11
      h(57)=
     &3.692108808871129411604189196259677640440919369478263728899602d-11
      h(58)=
     &-3.3270089671259799299106362463371508
     &51642079794871116041187279d-12
      h(59)=
     &-1.324334917243963163878274345609465
     &717294426628053460151843705d-12
      h(60)=
     &4.445467096291932163298411852093011459626037560439178917611592d-13
      h(61)=
     &-5.559442050579014337641375730083534521
     &513818164827556763756543d-14
      h(62)=
     &2.699382879762665647295493928801387173921314576598505507855504d-15
        ENDIF

        IF (M.EQ.32) THEN
      h(1)=
     &1.161463302135014885567464100760659332951431420121048996305591d-05
      h(2)=
     &2.466566906380903352739104211274667134470169443886449124673996d-04
      h(3)=
     &2.431261919572266100780423071905958127811969678055971488060574d-03
      h(4)=
     &1.468104638141913563547809006402194831107662001343421893488086d-02
      h(5)=
     &6.025749912033537081745451975527967031851677384078997261920024d-02
      h(6)=
     &1.757507836394388988189299915753348505208376399651864661397588d-01
      h(7)=
     &3.675096285973496361995340339143234125206079560406868595968025d-01
      h(8)=
     &5.343179193409538322901117858552186425529774700290587495921679d-01
      h(9)=
     &4.778091637339484033555130814414794130354053753675509287934741d-01
      h(10)=
     &1.206305382656178269538098710665261299391507308342013788891222d-01
      h(11)=
     &-2.666981814766755535489784087869865024226542605534080371507405d-1
      h(12)=
     &-2.774215815584272153338153320303401666681294506143291967655666d-1
      h(13)=
     &6.471335480551623831000090095167664918448659157720155321560811d-02
      h(14)=
     &2.483106423568801736064852157222867588791898170114101300999760d-01
      h(15)=
     &2.466244483969740441701479334808723214802614938081258920635302d-02
      h(16)=
     &-1.921023447085468984341365278247990525863123891147783426068990d-1
      h(17)=
     &-4.899511718467173853355943225576377418394280156945986899417475d-2
      h(18)=
     &1.452320794752866460838830744051944832326998342053148426312341d-01
      h(19)=
     &4.440490819993974022640619534046603571086531544468421519143629d-02
      h(20)=
     &-1.094561131160893831027722774343269232755171130623890041619420d-1
      h(21)=
     &-2.962787250844770491204452379051215505049068645551070779367843d-2
      h(22)=
     &8.087414063848395744090831590426327690818854671836423275412813d-02
      h(23)=
     &1.410615151610660772869738802931740150275269382463799031013905d-02
      h(24)=
     &-5.692631406247843550478416271158537960555270097953330567652364d-2
      h(25)=
     &-2.380264464932573834443178362086503847328134994591954135879789d-3
      h(26)=
     &3.705145792354468010437633458013030898015496905609424004450953d-02
      h(27)=
     &-4.145907660827218781460700428862611061267328108653649653634276d-3
      h(28)=
     &-2.166282283639119347634778516947485598599029367518033869601702d-2
      h(29)=
     &6.167527310685675112579059689520105004744367282412921739811164d-03
      h(30)=
     &1.101740071540688116532806119564345712473051769079712407908648d-02
      h(31)=
     &-5.411568257275791208581502410752383050600045942275647685361370d-3
      h(32)=
     &-4.649216751184411528658094984504900172989190128905887602541396d-3
      h(33)=
     &3.627224640687864960122122984391704782343548385375321260251988d-03
      h(34)=
     &1.468955100468467772528811782840480639166582822577191079260543d-03
      h(35)=
     &-1.964740555821778254183647540656746450092725858126595984907304d-3
      h(36)=
     &-2.211678729579097916278097586914956834196749138610403102772710d-4
      h(37)=
     &8.673058518450555343925662389563539890596549655683386287799624d-04
      h(38)=
     &-1.024537310607396186949656796812972062290796122915930356634122d-4    
      h(39)=
     &-3.059654423826911750479261161552574500739091332121504634422577d-4
      h(40)=
     &1.053915461739828114700905192091104141076083602686374410146603d-04
      h(41)=
     &8.103678329134838389828091896334156224227821362491626044950428d-05
      h(42)=
     &-5.259809282684322782648914338377962890245975842272425408122506d-5
      h(43)=
     &-1.294045779405512723950480259110995722517019870286295908085366d-5
      h(44)=
     &1.824268401980691220603850117995712615809177092802967489081228d-05
      h(45)=
     &-6.361781532260254953363913076575914206506177493714496098327288d-7
      h(46)=
     &-4.558309576264423135123964145585288808181431652781253437738445d-6
      h(47)=
     &1.202889036321620990296134494079846952404216422923750605507047d-06
      h(48)=
     &7.560047625595947819392627283726711361273296630256477108501994d-07
      h(49)=
     &-4.285970693151457255418342315045357407199066350632593899896712d-7
      h(50)=
     &-5.003361868748230293692887222336390314786090450819216035110269d-8
      h(51)=
     &8.965966311957728376981484572655177545054433542721057470726361d-08
      h(52)=
     &-1.219924359483373093110396748985081720383992859961285213840740d-8
      h(53)=
     &-1.104383021722648979552131128575075255513372249283096583736746d-8
      h(54)=
     &4.250422311980592983740943309197245384991941251563471671065543d-09
      h(55)=
     &4.384387799940474369553236949848427579687147486892033587998023d-10
      h(56)=
     &-5.881091462634605628881794361152305
     &108432139465417759716875076d-10
      h(57)=
     &8.904723796221605490455387579189371137903330749397374037644960d-11
      h(58)=
     &3.263270741332907875981844980104948375955551273115386408552080d-11
      h(59)=
     &-1.43091876516920232018802221173975
     &0594608742928641485026836608d-11
      h(60)=
     &1.075610653501062115165734990153347111902874668945095034791947d-12
      h(61)=
     &5.361482229611801638107331379599434078296259332654994508124989d-13
      h(62)=
     &-1.6638004894334023698898181929622598
     &23988673359967722467427927d-13
      h(63)=
     &2.000715303810524954375796020597627467104635766752154321244151d-14
      h(64)=
     &-9.42101913953507842131465536229108
     &8223782497046057523323473331d-16
       ENDIF

        IF (M.EQ.33) THEN
      h(1)=
     &8.186358314175091939858945975190102731733968885547217619434602d-06
      h(2)=
     &1.791016153702791479424389068736094134247294413108336017758506d-04
      h(3)=
     &1.822709435164084208084617771787691709255513374281497713580568d-03
      h(4)=
     &1.139594337458160925830840619716397130445853638888472948832932d-02
      h(5)=
     &4.861466653171619508385707681587366397164931431125053574327899d-02
      h(6)=
     &1.481863131800528081784673514426737436792606299953305691300616d-01
      h(7)=
     &3.267181301177075783930752787756046348844272437670999719562429d-01
      h(8)=
     &5.093761725149396552227892926384090200953139820961482931291482d-01
      h(9)=
     &5.112547705832674655425831875568453973369927971748064975152374d-01
      h(10)=
     &2.095823507130554216526494469993023406452629154801126958766008d-01
      h(11)=
     &-2.042026223985421049629055102642279430174095014493415546881477d-1
      h(12)=
     &-3.159974107665602561905181464284910961862968513875028980451424d-1
      h(13)=
     &-1.927833943695275915600583425408664108893845271616240406358226d-2
      h(14)=
     &2.454206121192791114179964351253140999836791489738418857473689d-01
      h(15)=
     &9.985155868033815698139640215477639365289384281516885362929979d-02
      h(16)=
     &-1.714280990518593279308738113273443832545615219650436927029674d-1
      h(17)=
     &-1.108441331167107910806084983056783194189909198734302929909672d-1
      h(18)=
     &1.219678564037346149389134584371009777591763921148126952722200d-01
      h(19)=
     &9.478808805061595889263191779090571160237408179346345390888721d-02
      h(20)=
     &-9.114696835133148913093153757138373418923462847746880902676089d-2
      h(21)=
     &-7.030248505405615921453280814171665167171986608963193275084895d-2
      h(22)=
     &7.019114394099653254998935842432841393915841096633514680190145d-02
      h(23)=
     &4.573456189389667743139040427641638967843459421665709740086516d-02
      h(24)=
     &-5.347125133582228919431110824663168583260050383336359554980188d-2
      h(25)=
     &-2.524858297747649929258392207837724793937727346177294684700378d-2
      h(26)=
     &3.868706076024496481748675031852528047303323816250150793091832d-02
      h(27)=
     &1.070326582001954942654534968137727769698168853186071888736311d-02
      h(28)=
     &-2.572876175473297336123211392278301875687760837710204579628265d-2
      h(29)=
     &-2.167758617353607324783298657172830203896433848418061622436727d-3
      h(30)=
     &1.531695411585766548347442266431874060229304787191589430967538d-02
      h(31)=
     &-1.594288782414604768637856446111392724059836934455189837500244d-3
      h(32)=
     &-7.953540387057939240459305406538116220678495240302592677582773d-3
      h(33)=
     &2.389062408165908575935815973439728988151836094753689966108405d-03
      h(34)=
     &3.480800953405711999411461002429227385937942254778524257436278d-03
      h(35)=
     &-1.860718214455795912074482150710567824317228203897000129729967d-3
      h(36)=
     &-1.204309257604658876916644980097327372892008586047095719636829d-3
      h(37)=
     &1.074380696351291355073899234941719080473877020595209197706651d-03
      h(38)=
     &2.727305847336937211749282358350196461733595290569540045817329d-04
      h(39)=
     &-4.908329007590351474487792254066540683724948757382104652497458d-4
      h(40)=
     &4.393166251766185755059005296958129844094063524324718175254673d-06
      h(41)=
     &1.780431898251245351831728023200069586928513661382622116969992d-04
      h(42)=
     &-4.160438516273709306234368807933932360567787692918883118883736d-5
      h(43)=
     &-4.929564423417301834310231482621574127409950921583062559483686d-5
      h(44)=
     &2.423335398816890365621188379922041046073808819182024026589770d-05
      h(45)=
     &9.070805757828453800203677464921508178468256685438211818575040d-06
      h(46)=
     &-8.866121366757736169176034432364298134186929098274651022820760d-6
      h(47)=
     &-3.607516102879771631230351118595069330196155459105589342866625d-7
      h(48)=
     &2.288371276141527305481395545993763010565968667577768164201792d-06
      h(49)=
     &-4.426923407952870147984002129341809185622768353983550670755106d-7
      h(50)=
     &-3.985791291985944076942626511739220753169387460984290019185514d-7
      h(51)=
     &1.822443332571053437467128998002798233969112236553215291639303d-07
      h(52)=
     &3.377972703730854377516206663481869099376154259897212784144779d-08
      h(53)=
     &-3.987838198518880722819502850814936369197384392561970319349663d-8
      h(54)=
     &3.672863576838181340505563759379169099717712645283448779390320d-09
      h(55)=
     &5.111211857347453839549366593998758891130921028374576213256027d-09
      h(56)=
     &-1.671392677251932495173219614104411841891545601521784559793012d-9
      h(57)=
     &-2.4964021052461936480735192693701973
     &31176405371538404298745013d-10
      h(58)=
     &2.426833102305682309891302883361232297664099485514601790344279d-10
      h(59)=
     &-3.0495744539458634303612969314551
     &41500128170151643206937547928d-11
      h(60)=
     &-1.4202368598899367924370778449404127
     &49343225644487770840543290d-11
      h(61)=
     &5.509414720765524548752673631197714447818740985929081064907524d-12
      h(62)=
     &-3.3434812189532787659825327226899847
     &25170758193566174566492199d-13
      h(63)=
     &-2.1524883868333026185206035456859
     &94753329478275805993737095214d-13
      h(64)=
     &6.214740247174398315576214699577230693021307854673557214652751d-14
      h(65)=
     &-7.1965105453633224140336544707790705
     &92316600780697558361083151d-15
      h(66)=
     &3.289373678416306368625564108782095644036415401902518812978798d-16
       ENDIF

        IF (M.EQ.34) THEN
      h(1)=
     &5.770510632730285627466067796809329117324708919047900817738025d-06
      h(2)=
     &1.299476200679530037833484815390569400369432658207722720405084d-04
      h(3)=
     &1.364061390059049998200014449396877439591680435610837369411339d-03
      h(4)=
     &8.819889403884978803182764563095879335330977939541630862804757d-03
      h(5)=
     &3.904884135178594138905026219591569204043816577941517019631916d-02
      h(6)=
     &1.241524821113768081954449898210969172708199672428635378051285d-01
      h(7)=
     &2.877650592337145629334256618087718872558560120999651277991839d-01
      h(8)=
     &4.784787462793710621468610706120519466268010329031345843336104d-01
      h(9)=
     &5.305550996564631773133260223990794445605699030503652382795600d-01
      h(10)=
     &2.903663295072749510455945186199530115755664977934564128822650d-01
      h(11)=
     &-1.282468421744371672912377747048558427612774932943748628650824d-1
      h(12)=
     &-3.315253015083869417715548463087537345035828886426345397256876d-1
      h(13)=
     &-1.038919155156404718287260506925867970596448618647006698388596d-1
      h(14)=
     &2.169072201874275950610018667099322465619408030256534197819784d-01
      h(15)=
     &1.666017504122074437311574334509261366682993700573488534577890d-01
      h(16)=
     &-1.273373582238011562843862636988693890108793629966541695807247d-1
      h(17)=
     &-1.609249271778668063014799490429649196614628857267382976958607d-1
      h(18)=
     &7.799184693794810738265349531832015087096882277333968473726399d-02
      h(19)=
     &1.341259602711361284802399913977387999358280900708582462625539d-01
      h(20)=
     &-5.448296806413904636632671383140642554265865948686157271017286d-2
      h(21)=
     &-1.029475969928140852342073823689090498245496056845473569066667d-1
      h(22)=
     &4.357609464963129726428486610925800727137724136370669421246609d-02
      h(23)=
     &7.318523543679560555546221335452045680757998947493883124934567d-02
      h(24)=
     &-3.701283841786244960356402125554190040750079009127461655784927d-2
      h(25)=
     &-4.743855964527776247220681410983851377889756018716427358008296d-2
      h(26)=
     &3.073974657395934459931226513844134346305562928466993208164603d-02
      h(27)=
     &2.722835075635419610095839895805858855202745897718117731496534d-02
      h(28)=
     &-2.367173792282636485046786438094940427456079528043555566867110d-2
      h(29)=
     &-1.314398001665716086105827506126287041342680578404007359439612d-2
      h(30)=
     &1.640937419986519252112261495537409592363156309874473310057471d-02
      h(31)=
     &4.713649260999809905918876125437488856235874027077755004539205d-03
      h(32)=
     &-1.004550670836151917439146861146431000364858401181337134891421d-2
      h(33)=
     &-6.194748845153872839014356621835501857322345445234809347431098d-4
      h(34)=
     &5.334950768759936032170270195983921511565539100791906952901398d-03
      h(35)=
     &-7.692127975067836975989490900561029844887285335804349474993607d-4
      h(36)=
     &-2.399453943537055863933124827688081952701780599883067560501870d-3
      h(37)=
     &8.589959874363661955444898475746536583497522107459291718900058d-04
      h(38)=
     &8.751999064078688732610570055224339733760304773327228476255647d-04
      h(39)=
     &-5.527355762144197975516415296735124460550632283763688359649888d-4
      h(40)=
     &-2.326732140233531635428863212833942245597361085708567528230733d-4
      h(41)=
     &2.650772397558057819755811309071002543822145660933016957735937d-04
      h(42)=
     &2.660050018453441903046828468025589086403126180798464347801678d-05
      h(43)=
     &-9.914697770780134603580350758869378471802751837608461971022567d-5
      h(44)=
     &1.353117227249649581251887376414486225127346352042209141315562d-05
      h(45)=
     &2.844951419697807376503080001943765930601242225183893658540032d-05
      h(46)=
     &-1.057657494257950623848316304755218120233253479317574337409622d-5
      h(47)=
     &-5.710826510998303938275050074333400305512451419983646591762318d-6
      h(48)=
     &4.169871758547028398316761659984928804362023643629741358799744d-06
      h(49)=
     &4.979718101421307748081857636471761057429219265531618602960147d-07
      h(50)=
     &-1.116306534817008428597995070751765080383261658112656948526954d-6
      h(51)=
     &1.448195708333185127061180618150009526758658641231104901703561d-07
      h(52)=
     &2.025990666667859216690536885693725545344933235432307649205497d-07
      h(53)=
     &-7.526701740412589411177481797841044281662555785969415398369019d-8
      h(54)=
     &-1.990346501531736915866180448337614967570744211158241514589121d-8
      h(55)=
     &1.740423332936068076497051274445147160190783847854409836489662d-08
      h(56)=
     &-8.66574426136872221586474116624538588881
     &8567571145958531936939d-10
      h(57)=
     &-2.316501946995482751582294240136010067415084499025753117941001d-9
      h(58)=
     &6.446378210323402313101214894500231181606520211579581132442548d-10
      h(59)=
     &1.300410318609415248880403259300467720631189120978928377152233d-10
      h(60)=
     &-9.904774537632409015479530333979124540
     &183199174591377762845227d-11
      h(61)=
     &1.004208735461769864836516428998306778031143650101842361622330d-11
      h(62)=
     &6.080125354000167254059025929915591291115751734288584563131636d-12
      h(63)=
     &-2.1078791089153015462853703954437788
     &64676275235126044599683271d-12
      h(64)=
     &9.799451158211597727901178520526388692140586041163624252991805d-14
      h(65)=
     &8.579194051799733179793112298652600511486581216528683482143106d-14
      h(66)=
     &-2.31708370390640848107825708190308952
     &3234020423092175261925515d-14
      h(67)=
     &2.587338381935699555813538163144986688834142571207152879144731d-15
      h(68)=
     &-1.14894475448059012824481579431260
     &6245147888158018823490936280d-16
       ENDIF

        IF (M.EQ.35) THEN
      h(1)=
     &4.067934061148559026665247110206084571051201477121972612218005d-06
      h(2)=
     &9.421469475576740631603027533116630224451049736050903361458759d-05
      h(3)=
     &1.019122680375098109319314672751485080202557607467199213778085d-03
      h(4)=
     &6.807292884319132011971333979015625113494050642797397817625326d-03
      h(5)=
     &3.123628851149071453063391210769353068187088999495893257051179d-02
      h(6)=
     &1.034044558614783789938787754929279183985553322796063517049140d-01
      h(7)=
     &2.513073789944933128513251971488905042866779761014740192816902d-01
      h(8)=
     &4.435927392240354378183910489448494594782039032807956294826105d-01
      h(9)=
     &5.370084275091661028670690231716974547580034932361053607723887d-01
      h(10)=
     &3.603456405180473278744458573988718422538114217890792270621563d-01
      h(11)=
     &-4.388388187393404111343479394097224312100349011932028865098625d-2
      h(12)=
     &-3.238228649121161212147302807993176715625480327235512530593160d-1
      h(13)=
     &-1.817869767667278325788350264528191676841493369460849123538616d-1
      h(14)=
     &1.660413574907809195438433327470947940538097914525298064477785d-01
      h(15)=
     &2.172992893210892977675493456199559114036326358517672106972956d-01
      h(16)=
     &-6.526287131067753892154895911331108284007380738865652420304233d-2
      h(17)=
     &-1.919195892985939528760786800798636198516495957924798820500876d-1
      h(18)=
     &1.930954466601835091947734585938109944647435243484967057775110d-02
      h(19)=
     &1.552924803962371144206753760712566993987319378965231186477630d-01
      h(20)=
     &-4.752680834111350445288110998030979143710864689041902167119118d-3
      h(21)=
     &-1.205855226433935545076589480704957722635324456812322150437989d-1
      h(22)=
     &4.734229172641948763293980314992213293971770695480616789828384d-03
      h(23)=
     &8.991354757072954417865374195261962983644048998218233900481856d-02
      h(24)=
     &-9.318558949903924837875002823617504227246562152671894579504378d-3
      h(25)=
     &-6.335603744044346612098887534020545705731671718057964802006671d-2
      h(26)=
     &1.322854958503655524455929847605110719648746890497356808289302d-02
      h(27)=
     &4.125469306470509212749750814299126656151504805845417994651417d-02
      h(28)=
     &-1.436683978422007182104025173214012797788904894291716373493525d-2
      h(29)=
     &-2.416949780166026740294880681731084091264533168816746227537030d-2
      h(30)=
     &1.276645671565674419403918018742432714973656598227939824940035d-02
      h(31)=
     &1.228943600811871086161967625814297050611100200023898377949151d-02
      h(32)=
     &-9.577797899235709998147309703713518608283233882793489733491642d-3
      h(33)=
     &-5.085991649233429881797636583578921194675393807761154549733547d-3
      h(34)=
     &6.137754586740521089596801883631921221145712545042519987641234d-03
      h(35)=
     &1.428088794070762107355585870669842132609159040625895090070111d-03
      h(36)=
     &-3.357644380922383229567732565298665639037348585961127075507937d-3
      h(37)=
     &7.615969435172736546769649923895317451534703066016116257300160d-06
      h(38)=
     &1.549637469702362975561719246539787717204438637997824935787688d-03
      h(39)=
     &-3.346692164250854961608526121524596908041109918361306282201310d-4
      h(40)=
     &-5.864810318991817532175809224131456738367101035694188223408841d-4
      h(41)=
     &2.648328819961289039302810122699710966048565368047575218693134d-04
      h(42)=
     &1.700012283661249043584690194716767771204207742625746308522935d-04
      h(43)=
     &-1.365883072261161602559926714744746422567509177443594045709653d-4
      h(44)=
     &-2.976995962848509743944225866488519668585242655980656646544319d-5
      h(45)=
     &5.304143122913310222538317980686374696005605533475685587486683d-05
      h(46)=
     &-2.437001526827789860990429478540556752694389693432668831073769d-6
      h(47)=
     &-1.572442077270281693663288966405861215692805972737981986121447d-5
      h(48)=
     &4.308047861716731191350493437937513220737450410132878032163179d-06
      h(49)=
     &3.353345862871309889390877168046133657377105681618708355266688d-06
      h(50)=
     &-1.895929617693153288493891051875444439753318548105998166574535d-6
      h(51)=
     &-3.903931733287306166657519468494511920760767388397825775326745d-7
      h(52)=
     &5.302368616904760917074352633915743250769600635829229600812520d-07
      h(53)=
     &-3.700308378205124537986402644918879149894035910106489082512364d-8
      h(54)=
     &-9.990396944534900755781728477561240762191443422318249128866740d-8
      h(55)=
     &3.008188650719066928230268918661718274504955045022550217051301d-08
      h(56)=
     &1.084902733789934825266560240100449884702749303326571747323086d-08
      h(57)=
     &-7.458116552893037631192407611262788593505988638365840409367117d-9
      h(58)=
     &5.897951310384361575470355861162022501172491937837712969865619d-11
      h(59)=
     &1.030823345485433383811700481488557422005210168069163779730908d-09
      h(60)=
     &-2.433545573751672936168877250405
     &940817227367937230289801251648d-10
      h(61)=
     &-6.40793825650188901843060832323597
     &4406219193176918284664973727d-11
      h(62)=
     &4.000536627253744510742788201354093006471710416671002244302586d-11
      h(63)=
     &-3.1256393571085575405980982286781507
     &68528121565391376265627294d-12
      h(64)=
     &-2.567065476155081449204643852428401530
     &283519685638256074752850d-12
      h(65)=
     &8.015088533687900921948605418789324826115616416343391081288979d-13
      h(66)=
     &-2.597954328893848084315198205094389
     &145706680129208998638802995d-14
      h(67)=
     &-3.39772085679626743195678382565906959
     &6940335130100871912329556d-14
      h(68)=
     &8.624037434720089202680337663692777682810714650060805832406135d-15
      h(69)=
     &-9.29801252932418542092155566471986350
     &1848315099116725184370339d-16
      h(70)=
     &4.014628712333488654318569164614220308046021091178184654250982d-17
       ENDIF

        IF (M.EQ.36) THEN
      h(1)=
     &2.867925182755946334630479473029238615535511775894262711054705d-06
      h(2)=
     &6.826028678546358691748629102209605362240344266505035981791715d-05
      h(3)=
     &7.602151099668488285869792677106082100141275054892389379198545d-04
      h(4)=
     &5.240297377409884366201603524392995696042174937194435235003941d-03
      h(5)=
     &2.489056564482796484885927333959115579403023347044729739255255d-02
      h(6)=
     &8.565209259526409083864716995521111486437594750377856524772704d-02
      h(7)=
     &2.177569530979008149637945915719999746248969705650625533415876d-01
      h(8)=
     &4.064336977082553467407793990250384445903151630768558142125382d-01
      h(9)=
     &5.322668952607286914777444748641462027213554723153906901129337d-01
      h(10)=
     &4.178753356009697863620634559374236455222275302996931178265919d-01
      h(11)=
     &4.397519752934862993862182898358763783110745559238982179690132d-02
      h(12)=
     &-2.944210395891145711100715969898758940722458887377844633443675d-1
      h(13)=
     &-2.468070369781255270524798278622698446566520718230313889086016d-1
      h(14)=
     &9.811420416311477050518401371401568038943437322299913514049728d-02
      h(15)=
     &2.465372776089742110529709111809595434656418762898152706621356d-01
      h(16)=
     &7.278515095792229009687682299460382878643139026668958884429641d-03
      h(17)=
     &-1.993372056086496198603363400094784142714162256792182570541036d-1
      h(18)=
     &-4.586140074639271639145126228774831743002971373998329604574394d-2
      h(19)=
     &1.541062366276428841776316300420654875883842819413623395358262d-01
      h(20)=
     &5.027618007353842862036816972809884096761706036019748316890913d-02
      h(21)=
     &-1.188037543101356316801816931383547446073152951044444224449501d-1
      h(22)=
     &-3.988085357551317584091699967924044034100374257075864260934102d-2
      h(23)=
     &9.115678225801654406336059281306715151058903055370522031843771d-02
      h(24)=
     &2.503872144956848989919484296709846860569180993040383621980546d-02
      h(25)=
     &-6.820901663681751124880436344265538690580358108714540763125119d-2
      h(26)=
     &-1.131910031681742794381808082173695022123056280821611354577883d-2
      h(27)=
     &4.851308354780908538616267662315735632292989749013261207046367d-02
      h(28)=
     &1.424972661765391603147802607378542396323429657660009755652404d-03
      h(29)=
     &-3.198072067763969654470293513742344601172739688274251641873778d-2
      h(30)=
     &3.984040198717004857397179486790082321314291366656151213429068d-03
      h(31)=
     &1.906359478062535932877576164368198274858108513696832728889209d-02
      h(32)=
     &-5.657813245058818380424016973516714570499161434975761798379020d-3
      h(33)=
     &-9.990263473281372348001743806489172665465685056975652497503772d-3
      h(34)=
     &5.022989106665829004699819220796538830393945994687289792465541d-03
      h(35)=
     &4.413484835350575251918616780287775585471012556848037301025999d-03
      h(36)=
     &-3.484541445404883311209541395428535732697661971818727286003028d-3
      h(37)=
     &-1.503074066296643749549363655363411879858070202740814054964603d-3
      h(38)=
     &1.990793771851737270404293245701878186600899439513475823305914d-03
      h(39)=
     &2.776812795712026068152384207605140383490242756921936501940389d-04
      h(40)=
     &-9.463403823261101964604918059447913047725482130063492242779878d-4
      h(41)=
     &8.614565758992702032613879159402330909634737204578606399403107d-05
      h(42)=
     &3.693507284967510502620040341882236687749563414433432842567511d-04
      h(43)=
     &-1.155118895843527096848376999413102395191976350936666573818799d-4
      h(44)=
     &-1.131899468084665671727391922924411467938450743565106978099456d-4
      h(45)=
     &6.694741196930590257104231749283786251555566773398199990337698d-05
      h(46)=
     &2.375106683660860777161950832380341362257503761490580896617678d-05
      h(47)=
     &-2.731390824654337912922346414722045404779935825834384250023192d-5
      h(48)=
     &-1.183471059985615942783182762352360917304348034947412986608322d-6
      h(49)=
     &8.372218198160788432628056043217491552198857358432112275253310d-06
      h(50)=
     &-1.586145782434577495502614631566211839722879492827911790709498d-6
      h(51)=
     &-1.870811602859180713762972281154953528056257451900381097476968d-6
      h(52)=
     &8.311421279707778528163597405935375886855029592150424544500718d-07
      h(53)=
     &2.548423522556577831218519052844387478819866531902854523544709d-07
      h(54)=
     &-2.455377658434232699135878286794578515387138194247693201846263d-7
      h(55)=
     &2.753249073339512254085076456700241929492720457889076058451072d-09
      h(56)=
     &4.799043465450992009934526867650497683545716858606119786327559d-08
      h(57)=
     &-1.156093688817008406756913949175208452083765368825442482226093d-8
      h(58)=
     &-5.612784343327791397474114357094368557982413895802980814813369d-9
      h(59)=
     &3.138841695782424018351567952158415003571380699236147752239001d-09
      h(60)=
     &1.090815553713751810964713058800448676068475673611349566405716d-10
      h(61)=
     &-4.51254577856324963442520085608849019500
     &4077806062978067796020d-10
      h(62)=
     &8.962418203859611987065968320295929679774693465791367610044773d-11
      h(63)=
     &3.037429098112535221800013609576297196061786927734556635696416d-11
      h(64)=
     &-1.5997166892613571432003969224094485
     &15398648489795044468046420d-11
      h(65)=
     &8.876846287217374213524399682895564055949886050748321818411161d-13
      h(66)=
     &1.070969357114017002424433471621197579059927261727846375968378d-12
      h(67)=
     &-3.02928502697487726889613458976947
     &3854669758797446795757329862d-13
      h(68)=
     &5.542263182639804235231685861028995158694397223907295269180336d-15
      h(69)=
     &1.338071386299105896025578761458472955294763310766371178363783d-14
      h(70)=
     &-3.204628543401749860439316638848579
     &711789176444320134355253750d-15
      h(71)=
     &3.339971984818693213132578777712503670014459411167839211495237d-16
      h(72)=
     &-1.4032741753731906174898232091680
     &13922564353495443487431242610d-17
       ENDIF

        IF (M.EQ.37) THEN
      h(1)=
     &2.022060862498392121815038335333633351464174415618614893795880d-06
      h(2)=
     &4.942343750628132004714286117434454499485737947791397867195910d-05
      h(3)=
     &5.662418377066724013768394373249439163518654840493603575144737d-04
      h(4)=
     &4.024140368257286770702140124893772447952256842478891548092703d-03
      h(5)=
     &1.976228615387959153244055502205017461538589475705618414896893d-02
      h(6)=
     &7.058482597718160832030361890793007659963483925312132741868671d-02
      h(7)=
     &1.873263318620649448028843491747601576761901656888288838192023d-01
      h(8)=
     &3.684409724003061409445838616964941132670287724754729425204047d-01
      h(9)=
     &5.181670408556228873104519667534437205387109579265718071174178d-01
      h(10)=
     &4.622075536616057145505448401528172070050768534504278694229363d-01
      h(11)=
     &1.308789632330201726057701201017649601034381070893275586898075d-01
      h(12)=
     &-2.461804297610834132869018581145720710365433914584680691693717d-1
      h(13)=
     &-2.943759152626617722808219575932673733674290772235644691367427d-1
      h(14)=
     &1.967150045235938977077768648740052380288156507222647187301894d-02
      h(15)=
     &2.515232543602686933435224095078166291442923992611593827552710d-01
      h(16)=
     &8.180602838721862339029076982652411696000045533716726027662147d-02
      h(17)=
     &-1.819622917786080007408824256525225216444443143868752611284260d-1
      h(18)=
     &-1.084517138233017845554078812341876568514835176341639783558543d-1
      h(19)=
     &1.299296469598537527842528895259188653120602318620944502979726d-01
      h(20)=
     &1.017802968388141797470948228505865617480048287983176581607964d-01
      h(21)=
     &-9.660754061668439030915405045955772715988585374771282291315496d-2
      h(22)=
     &-8.233021190655740867404073660920379414988302492018783774702028d-2
      h(23)=
     &7.504761994836017933579005072594245435071674452882148228583865d-02
      h(24)=
     &5.956741087152995245435589042520108066877114768216272503684398d-02
      h(25)=
     &-5.925681563265897095153806724965924334077555174281436189512239d-2
      h(26)=
     &-3.825382947938424882011108885090442116802994193611884738133373d-2
      h(27)=
     &4.580794415126833246633256156110381805848138158784734496981778d-02
      h(28)=
     &2.097280059259754883313769469036393294461497749083921162354229d-02
      h(29)=
     &-3.352358406410096994358662875913243067234786296009238949920582d-2
      h(30)=
     &-8.833493890410232394064187990625563257107429109130726291528648d-3
      h(31)=
     &2.261865154459947356571431658958802912061105608212828675323452d-02
      h(32)=
     &1.690472383484423743663952859090705636512807161536954018400081d-03
      h(33)=
     &-1.376398196289478433857985486097070339786225136728067000591187d-2
      h(34)=
     &1.519305778833399218481261844599507408563295102235964076544334d-03
      h(35)=
     &7.387757452855583640107787619408806919082115520707105052944171d-03
      h(36)=
     &-2.248053187003824706127276829147166466869908326245810952521710d-3
      h(37)=
     &-3.394523276408398601988475786247462646314228994098320665709345d-3
      h(38)=
     &1.816871343801423525477184531347879515909226877688306010517914d-03
      h(39)=
     &1.263934258117477182626760951047019242187910977671449470318766d-03
      h(40)=
     &-1.111484865318630197259018233162929628309920117691177260742614d-3
      h(41)=
     &-3.280788470880198419407186455190899535706232295554613820907245d-4
      h(42)=
     &5.490532773373631230219769273898345809368332716288071475378651d-04
      h(43)=
     &1.534439023195503211083338679106161291342621676983096723309776d-05
      h(44)=
     &-2.208944032455493852493630802748509781675182699536797043565515d-4
      h(45)=
     &4.336726125945695214852398433524024058216834313839357806404424d-05
      h(46)=
     &7.055138782065465075838703109997365141906130284669094131032488d-05
      h(47)=
     &-3.098662927619930052417611453170793938796310141219293329658062d-5
      h(48)=
     &-1.639162496160583099236044020495877311072716199713679670940295d-5
      h(49)=
     &1.354327718416781810683349121150634031343717637827354228989989d-05
      h(50)=
     &1.849945003115590390789683032647334516600314304175482456338006d-06
      h(51)=
     &-4.309941556597092389020622638271988877959028012481278949268461d-6
      h(52)=
     &4.854731396996411681769911684430785681028852413859386141424939d-07
      h(53)=
     &1.002121399297177629772998172241869405763288457224082581829033d-06
      h(54)=
     &-3.494948603445727645895194867933547164628229076947330682199174d-7
      h(55)=
     &-1.509885388671583553484927666148474078148724554849968758642331d-7
      h(56)=
     &1.109031232216439389999036327867142640916239658806376290861690d-07
      h(57)=
     &5.350657515461434290618742656970344024396382191417247602674540d-09
      h(58)=
     &-2.252193836724805775389816424695618411834716065179297102428180d-8
      h(59)=
     &4.224485706362419268050011630338101126995607958955688879525896d-09
      h(60)=
     &2.793974465953982659829387370821677112004867350709951380622807d-09
      h(61)=
     &-1.297205001469435139867686007585972538983682739297235604327668d-9
      h(62)=
     &-1.0314111290969749656779506464981530
     &71722880698222864687038596d-10
      h(63)=
     &1.946164894082315021308714557636277980079559327508927751052218d-10
      h(64)=
     &-3.203398244123241367987902201268363088
     &933939831689591684670080d-11
      h(65)=
     &-1.3984157155376414879595516825574833486
     &61602836709278513081908d-11
      h(66)=
     &6.334955440973913249611879065201632922100533284261000819747915d-12
      h(67)=
     &-2.0963631942348005416147757427555557132795
     &49381264881030843258d-13
      h(68)=
     &-4.4216124098721053673335727348544
     &01373201808896976552663098518d-13
      h(69)=
     &1.138052830921439682522395208295427884729893377395129205716662d-13
      h(70)=
     &-4.518889607463726394454509623712773172513
     &778367070839294449849d-16
      h(71)=
     &-5.24302569188420583226035450374832533
     &4301994904062750850180233d-15
      h(72)=
     &1.189012387508252879928637969242590755033933791160383262132698d-15
      h(73)=
     &-1.199280335852879554967035114674445
     &327319437557227036460257649d-16
      h(74)=
     &4.906615064935203694857690087429901193139905690549533773201453d-18
       ENDIF

        IF (M.EQ.38) THEN
      h(1)=
     &1.425776641674131672055420247567865803211784397464191115245081d-06
      h(2)=
     &3.576251994264023012742569014888876217958307227940126418281357d-05
      h(3)=
     &4.211702664727116432247014444906469155300573201130549739553848d-04
      h(4)=
     &3.083088119253751774288740090262741910177322520624582862578292d-03
      h(5)=
     &1.563724934757215617277490102724080070486270026632620664785632d-02
      h(6)=
     &5.788994361285925649727664279317241952513246287766481213301801d-02
      h(7)=
     &1.600719935641106973482800861166599685169395465055048951307626d-01
      h(8)=
     &3.307757814110146511493637534404611754800768677041577030757306d-01
      h(9)=
     &4.965911753117180976599171147718708939352414838951726087564419d-01
      h(10)=
     &4.933560785171007975728485346997317064969513623594359091115804d-01
      h(11)=
     &2.130505713555785138286743353458562451255624665951160445122307d-01
      h(12)=
     &-1.828676677083358907975548507946239135218223185041410632924815d-1
      h(13)=
     &-3.216756378089978628483471725406916361929841940528189059002548d-1
      h(14)=
     &-6.226650604782432226643360160478765847565862101045597180310490d-2
      h(15)=
     &2.321259638353531085028708104285994998671615563662858079262996d-01
      h(16)=
     &1.499851196187170199586403453788927307298226028262603028635758d-01
      h(17)=
     &-1.417956859730596216710053144522330276392591055375830654519080d-1
      h(18)=
     &-1.599125651582443618288533214523534937804208844386102639177693d-1
      h(19)=
     &8.563812155615105741612217814369165313487129645536001850276987d-02
      h(20)=
     &1.414147340733826800884683119379170594092606174915755283496153d-01
      h(21)=
     &-5.658645863072738145681787657843320646815509410635114234947902d-2
      h(22)=
     &-1.147311707107443752394144019458942779715665489230169950201022d-1
      h(23)=
     &4.309589543304764288137871223616030624246568683595408792078602d-02
      h(24)=
     &8.720439826203975011910714164154456762073786124233088471855868d-02
      h(25)=
     &-3.660510340287429567372071039506772372567938710943432838908247d-2
      h(26)=
     &-6.176620870841315993604736705613246241897497782373337911398117d-2
      h(27)=
     &3.198987753153780630818381136366859026137035450576631134176875d-02
      h(28)=
     &4.005498110511594820952087086241114309038577379366732959648548d-02
      h(29)=
     &-2.689149388089451438550851767715967313417890393287236700072071d-2
      h(30)=
     &-2.311413402054931680856913553585621248925303865540203357180768d-2
      h(31)=
     &2.090464525565524340215982365351342094670261491526831672682244d-02
      h(32)=
     &1.129049727868596484270081487761544232851115891449843967151657d-02
      h(33)=
     &-1.470188206539868213708986402816605045648481224662435114088245d-2
      h(34)=
     &-4.131306656031089274123231103326745723188134548520938157995702d-3
      h(35)=
     &9.214785032197180512031534870181734003522861645903894504302286d-03
      h(36)=
     &5.625715748403532005741565594881148757066703437214522101740941d-04
      h(37)=
     &-5.071314509218348093935061417505663002006821323958752649640329d-3
      h(38)=
     &7.169821821064019257784165364894915621888541496773370435889585d-04
      h(39)=
     &2.400697781890973183892306914082592143984140550210130139535193d-03
      h(40)=
     &-8.448626665537775009068937851465856973251363010924003314643612d-4
      h(41)=
     &-9.424614077227377964015942271780098283910230639908018778588910d-4
      h(42)=
     &5.810759750532863662020321063678196633409555706981476723988312d-04
      h(43)=
     &2.817639250380670746018048967535608190123523180612961062603672d-04
      h(44)=
     &-3.031020460726611993600629020329784682496477106470427787747855d-4
      h(45)=
     &-4.555682696668420274688683005987764360677217149927938344795290d-5
      h(46)=
     &1.262043350166170705382346537131817701361522387904917335958705d-04
      h(47)=
     &-1.155409103833717192628479047983460953381959342642374175822863d-5
      h(48)=
     &-4.175141648540397797296325065775711309197411926289412468280801d-5
      h(49)=
     &1.334176149921350382547503457286060922218070031330137601427324d-05
      h(50)=
     &1.037359184045599795632258335010065103524959844966094870217687d-05
      h(51)=
     &-6.456730428469619160379910439617575420986972394137121953806236d-6
      h(52)=
     &-1.550844350118602575853380148525912999401292473185534395740371d-6
      h(53)=
     &2.149960269939665207789548199790770596890252405076394885606038d-06
      h(54)=
     &-8.487087586072593071869805266089426629606479876982221840833098d-8
      h(55)=
     &-5.187733738874144426008474683378542368066310000602823096009187d-7
      h(56)=
     &1.396377545508355481227961581059961184519872502493462010264633d-07
      h(57)=
     &8.400351046895965526933587176781279507953080669259318722910523d-08
      h(58)=
     &-4.884757937459286762082185411608763964041010392101914854918157d-8
      h(59)=
     &-5.424274800287298511126684174854414928447521710664476410973981d-9
      h(60)=
     &1.034704539274858480924046490952803937328239537222908159451039d-08
      h(61)=
     &-1.436329487795135706854539856979275911183628476521636251660849d-9
      h(62)=
     &-1.349197753983448821850381770889786301246741304307934955997111d-9
      h(63)=
     &5.261132557357598494535766638772624572100332209198979659077082d-10
      h(64)=
     &6.732336490189308685740626964182623159759767536724844030164551d-11
      h(65)=
     &-8.278256522538134727330692938158991
     &115335384611795874767521731d-11
      h(66)=
     &1.101692934599454551150832622160224231280195362919498540913658d-11
      h(67)=
     &6.291537317039508581580913620859140835852886308989584198166174d-12
      h(68)=
     &-2.484789237563642857043361214502760
     &723611468591833262675852242d-12
      h(69)=
     &2.626496504065252070488282876470525379851429538389481576454618d-14
      h(70)=
     &1.808661236274530582267084846343959377085922019067808145635263d-13
      h(71)=
     &-4.2498178195714630069666163715542065
     &72863122562744916796556474d-14
      h(72)=
     &-4.5633971621273731091016916430479237
     &47796563449194075621854491d-16
      h(73)=
     &2.045099676788988907802272564402310095398641092819367167252952d-15
      h(74)=
     &-4.40530704248346134244902713983830
     &1611006835285455050155842865d-16
      h(75)=
     &4.304596839558790016251867477122791508849697688058169053134463d-17
      h(76)=
     &-1.71615245108874418873240428173796
     &4277713026087224248235541071d-18
       ENDIF




            DO k=1,2*M
            g(k)=(-1)**(k-1)*h(2*M-k+1)
            ENDDO

      RETURN
      END      
      
      
