C -----------------------------------------------------------------------------
      SUBROUTINE NanofluidProp(inp,io,env,a,b,c)
C
C     $: Calculates nanofluid properties
C        w, diff clen : Pr * a
C        w, vzgon     : PrRa * b
C        T, diff clen : 1 * c
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(InputType) inp
      TYPE(IOtype)    :: io
      TYPE(penv) :: env

      REAL(8) muNFmuF
      REAL(8) rhoNFrhoF
      REAL(8) bNFbF
      REAL(8) kNFkF
      REAL(8) rhocpNFrhocpF
      REAL(8) a,b,c

c     mu_nf / mu_f  dinamicna viskoznost
      muNFmuF=1.0D0/(1.0D0-inp%nanovf)**2.5D0

c     rho_nf / rho_f  gostota
      rhoNFrhoF = 1.0D0-inp%nanovf+inp%nanovf*inp%nanoS_rho/inp%nanoF_rho

c     beta_nf / beta_f  temperaturni raztezek
      bNFbF = inp%nanoS_beta/inp%nanoF_beta/
     &        (1.0D0+(1.0D0-inp%nanovf)*inp%nanoF_rho/(inp%nanovf*inp%nanoS_rho))+
     &        1.0D0/(1.0D0+inp%nanovf*inp%nanoS_rho/((1.0D0-inp%nanovf)*inp%nanoF_rho))

c     k_nf / k_f  toplotna prevodnost
      kNFkF = (inp%nanoS_k+2.0D0*inp%nanoF_k-2.0D0*inp%nanovf*(inp%nanoF_k-inp%nanoS_k))/
     &        (inp%nanoS_k+2.0D0*inp%nanoF_k+      inp%nanovf*(inp%nanoF_k-inp%nanoS_k))

c     rho*cp_nf / rho*cp_f gostota krat specificna toplota
      rhocpNFrhocpF = 1.0D0-inp%nanovf+inp%nanovf*inp%nanoS_rho*inp%nanoS_cp/(inp%nanoF_rho*inp%nanoF_cp)

      a =  muNFmuF / rhoNFrhoF
      b =  bNFbF
      c =  kNFkF / rhocpNFrhocpF

      IF (env%myproc.EQ.1) THEN
        WRITE(io%l,'(A)') "Simulating flow of nanofluid"
        WRITE(io%l,*) "fluid  cp [J/kgK]  =",inp%nanoF_cp
        WRITE(io%l,*) "fluid rho [kg/m3]  =",inp%nanoF_rho
        WRITE(io%l,*) "fluid k   [W /mK]  =",inp%nanoF_k
        WRITE(io%l,*) "fluid beta [K^-1]  =",inp%nanoF_beta
        WRITE(io%l,*) "solid  cp [J/kgK]  =",inp%nanoS_cp
        WRITE(io%l,*) "solid rho [kg/m3]  =",inp%nanoS_rho
        WRITE(io%l,*) "solid k   [W /mK]  =",inp%nanoS_k
        WRITE(io%l,*) "solid beta [K^-1]  =",inp%nanoS_beta
        WRITE(io%l,*) "solid volume frac  =",inp%nanovf
        WRITE(io%l,'(A)') ""
        WRITE(io%l,*) "muNFmuF            =",muNFmuF
        WRITE(io%l,*) "rhoNFrhoF          =",rhoNFrhoF
        WRITE(io%l,*) "bNFbF              =",bNFbF
        WRITE(io%l,*) "kNFkF              =",kNFkF
        WRITE(io%l,*) "rhocpNFrhocpF      =",rhocpNFrhocpF
        WRITE(io%l,'(A)') ""
        WRITE(io%l,*) "nano parameter a   =",a
        WRITE(io%l,*) "nano parameter b   =",b
        WRITE(io%l,*) "nano parameter c   =",c
        WRITE(io%l,'(A)') ""
        WRITE(io%l,*) "Prandtl number     =",inp%prn
        WRITE(io%l,*) "fl. din. visc.[Pas]=",inp%prn*inp%nanoF_k/inp%nanoF_cp
        WRITE(io%l,*) "fl.kin. visc.[m2/s]=",inp%prn*inp%nanoF_k/inp%nanoF_cp/inp%nanoF_rho
        WRITE(io%l,'(A)') ""
      END IF

c      print *,inp%inano,inp%nanoF_cp,inp%nanoF_rho,inp%nanoF_k,inp%nanoF_beta,
c     &                                    inp%nanoS_cp,inp%nanoS_rho,inp%nanoS_k,inp%nanoS_beta
c      print *,muNFmuF,rhoNFrhoF,bNFbF,kNFkF,rhocpNFrhocpF
c      print *,a,b,c
      END

C -----------------------------------------------------------------------------
      SUBROUTINE NanofluidMaterialProperties(inp,mesh,gp,rho,cp,lambda,PartVolFrac)
C
C      Variable nanofluid properties
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      TYPE(InputType) inp
      TYPE(GaussType) :: gp

      REAL(8) rho(mesh%nnodes)
      REAL(8) cp(mesh%nnodes)
      REAL(8) lambda(mesh%nnodes)
      REAL(8) PartVolFrac(mesh%nnodes)

      REAL(8) muNFmuF
      REAL(8) rhoNFrhoF
      REAL(8) bNFbF
      REAL(8) kNFkF
      REAL(8) rhocpNFrhocpF

      REAL(8) arho,acp,alambda,apvf,vol

      INTEGER i

      arho=0.0D0
      acp=0.0D0
      alambda=0.0D0
      apvf=0.0D0

      DO i=1,mesh%nnodes

c       mu_nf / mu_f  dinamicna viskoznost
        muNFmuF=1.0D0/(1.0D0-PartVolFrac(i))**2.5D0

c       rho_nf / rho_f  gostota
        rhoNFrhoF = 1.0D0-PartVolFrac(i)+PartVolFrac(i)*inp%nanoS_rho/inp%nanoF_rho

c       beta_nf / beta_f  temperaturni raztezek
        bNFbF = inp%nanoS_beta/inp%nanoF_beta/
     &        (1.0D0+(1.0D0-PartVolFrac(i))*inp%nanoF_rho/(PartVolFrac(i)*inp%nanoS_rho))+
     &        1.0D0/(1.0D0+PartVolFrac(i)*inp%nanoS_rho/((1.0D0-PartVolFrac(i))*inp%nanoF_rho))

c       k_nf / k_f  toplotna prevodnost
        kNFkF = (inp%nanoS_k+2.0D0*inp%nanoF_k-2.0D0*PartVolFrac(i)*(inp%nanoF_k-inp%nanoS_k))/
     &        (inp%nanoS_k+2.0D0*inp%nanoF_k+      PartVolFrac(i)*(inp%nanoF_k-inp%nanoS_k))

c       rho*cp_nf / rho*cp_f gostota krat specificna toplota
        rhocpNFrhocpF = 1.0D0-PartVolFrac(i)+PartVolFrac(i)*inp%nanoS_rho*inp%nanoS_cp/(inp%nanoF_rho*inp%nanoF_cp)

        rho(i)   = rhoNFrhoF
        cp(i)    = rhocpNFrhocpF / rhoNFrhoF
        lambda(i)= kNFkF

c       calculate average values
        arho=arho+rho(i)
        acp=acp+cp(i)
        alambda=alambda+lambda(i)
        apvf=apvf+PartVolFrac(i)


      END DO

c      print *,arho/DBLE(mesh%nnodes)
c      print *,acp/DBLE(mesh%nnodes)
c      print *,alambda/DBLE(mesh%nnodes)
c      print *,apvf/DBLE(mesh%nnodes)

c      CALL DomainIntegral(mesh,gp,rho,arho)
c      CALL DomainIntegral(mesh,gp,cp,acp)
c      CALL DomainIntegral(mesh,gp,lambda,alambda)
c      CALL DomainIntegral(mesh,gp,PartVolFrac,apvf)

c      vol=0.0D0
c      do i=1,mesh%nicell
c        vol=vol+mesh%CellVolume(i)
c      end do

c      print *,arho/vol
c      print *,acp/vol
c      print *,alambda/vol
c      print *,apvf/vol

c      stop
      END
