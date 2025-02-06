C
C
C          subdomain (macro) kinematics for interal velocity
C
C
C     List of subroutines:
C
c       solv_kmMacro3sq
c       rhsKM3sq
c       solv_kmMacro3
c       rhsKM3
c       SetUpXknv3
c       sMat2crsSysRhsKM3
c       DetKMunk3
C
C
C -----------------------------------------------------------------------------
      SUBROUTINE solv_kmMacro3sq(eqn,inp,cpu,mesh,nit,sysm,rhsm,prec,kmp,
     &                         smatHtx,smatHty,smatHtz,
     &                         smatDx,smatDy,smatDz,velocity,vorticity)
C
C     Solve the kinematics equation for domain velocities
c     determined system of equations
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(meshtype) mesh
      TYPE(MATRIX) sysm,rhsm,prec
      TYPE(KMpointer) kmp
      TYPE(InputType) :: inp
      TYPE(CPUtype) cpu

      INTEGER ierr,i,nit

      REAL(8), ALLOCATABLE :: b(:),x(:),knv(:)

      REAL(8) smatHtx(mesh%nicnpoc,mesh%npoc) ! doktorat, enacba (4.7)
      REAL(8) smatHty(mesh%nicnpoc,mesh%npoc) ! to je za prvi integral
      REAL(8) smatHtz(mesh%nicnpoc,mesh%npoc) ! na desni

      REAL(8) smatDx(mesh%nicnsp,mesh%npoc) ! doktorat, enacba (4.7)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc) ! to je za drugi integral
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc) ! na desni

      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3)
      INTEGER eqn
      REAL(4) cptime,rcpu,xcpu

c     alociram spomin za vektor neznank in vektor na desni
      ALLOCATE (knv(kmp%nrhscol)) ! vector of known values
      ALLOCATE (x(kmp%nsmcol))    ! vector of unknown values
      ALLOCATE (b(sysm%neq))    ! r.h.s vector of the system of linear equations b=rhsm * knv

      rcpu=cptime(0.)

C     Set up vectors of known an unknown values
      CALL SetUpXknv3(eqn,mesh,kmp,velocity,knv,x)

C     Set up r.h.s. vector : b= rhsm * knv
      CALL CRSxV(rhsm,knv,kmp%nrhscol,b)
c
c     b = b + matrike * rhsv
c
      CALL rhsKM3sq(eqn,mesh,b,sysm,smatHtx,smatHty,smatHtz,
     &                  smatDx,smatDy,smatDz,velocity,vorticity)

      cpu%time((eqn-1)*6+3)=cpu%time((eqn-1)*6+3)+cptime(rcpu)
      rcpu=cptime(0.)

c
c      solve
c
      CALL SolvSLE(inp%sqrs_type,inp%sqrs_prec,2,inp%sqrs_maxit,5,inp%dlse(eqn),sysm%neq,sysm%nnz,
     &                   prec,sysm,b,x,nit,xcpu,ierr)

      cpu%time((eqn-1)*6+4)=cpu%time((eqn-1)*6+4)+cptime(rcpu)

      IF (ierr.NE.0) print *,"KM NAPAKA V SOVLERJU",ierr

      DO i=1,kmp%nsmcol
        IF (kmp%bunt(i).EQ.1) THEN  ! vx
          velocity(kmp%bunn(i),eqn)=(1.0D0-inp%urDv(eqn))*velocity(kmp%bunn(i),eqn)+inp%urDv(eqn)*x(i)
        ELSE
         print *,"error!"
        END IF
      END DO


      DEALLOCATE (x,b,knv)

C     Employ dv/dn=0 boundary condition
      DO i=1,mesh%nbnodes
        IF (mesh%kode(i,eqn).EQ.2.OR.mesh%kode(i,eqn).EQ.3) THEN  ! d(v_eqn)/dn=0
          velocity(mesh%gbn(i),eqn)=velocity(mesh%kobc(i),eqn)
        END IF
      END DO

C     Employ periodic boundary condition
      IF (mesh%iPeri.GT.0) THEN
        DO i=1,mesh%nPeri
          velocity(mesh%pwFT(i,2),eqn)=velocity(mesh%pwFT(i,1),eqn)
        END DO
      END IF


      END

C -----------------------------------------------------------------------------
      SUBROUTINE rhsKM3sq(eqn,mesh,b,sysm,smatHtx,smatHty,smatHtz,
     &                  smatDx,smatDy,smatDz,velocity,vorticity)
C
C      sets up right hand side for kinematics eqaution
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(meshtype) mesh
      TYPE(matrix) sysm

      REAL(8) smatHtx(mesh%nicnpoc,mesh%npoc) ! doktorat, enacba (4.7)
      REAL(8) smatHty(mesh%nicnpoc,mesh%npoc) ! to je za prvi integral
      REAL(8) smatHtz(mesh%nicnpoc,mesh%npoc) ! na desni

      REAL(8) smatDx(mesh%nicnsp,mesh%npoc) ! doktorat, enacba (4.7)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc) ! to je za drugi integral
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc) ! na desni

      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3)

      REAL(8) b(sysm%neq)

      INTEGER i,ii,eqn,ic,node,j
      INTEGER isqr,inode,isuma

      i=0
c
c     vx
c
      IF (eqn.EQ.1) THEN

        DO isqr=1,sysm%neq !nunk ! zanka po vrsticah kvadratne matrike
          inode=mesh%KMeql(isqr)  ! vozlisce v U  mrezi, kateremu pripada ta enacba
c         dodajamo prispevke enack, ki se sestejejo
          DO isuma=1,mesh%sqUlistNO(inode)
            i=mesh%sqUlist(inode,isuma)
            ii=mesh%sqUlistKM(inode,isuma)
            ic=mesh%sqUlistIC(inode,isuma)
            DO j=1,mesh%npoc
              node=mesh%idc(ic,j)
              b(isqr)=b(isqr)
     &                   -smatHtz(ii,j)*velocity(node,2)
     &                   +smatHty(ii,j)*velocity(node,3)
     &                   +smatDz(i,j) *vorticity(node,2)
     &                   -smatDy(i,j) *vorticity(node,3)
            END DO
          END DO
        END DO

c
c     vy
c
      ELSE IF (eqn.EQ.2) THEN

        DO isqr=1,sysm%neq !nunk ! zanka po vrsticah kvadratne matrike
          inode=mesh%KMeql(isqr)  ! vozlisce v U  mrezi, kateremu pripada ta enacba
c         dodajamo prispevke enack, ki se sestejejo
          DO isuma=1,mesh%sqUlistNO(inode)
            i=mesh%sqUlist(inode,isuma)
            ii=mesh%sqUlistKM(inode,isuma)
            ic=mesh%sqUlistIC(inode,isuma)
            DO j=1,mesh%npoc
              node=mesh%idc(ic,j)
              b(isqr)=b(isqr)
     &                   +smatHtz(ii,j)*velocity(node,1)
     &                   -smatHtx(ii,j)*velocity(node,3)
     &                   -smatDz(i,j) *vorticity(node,1)
     &                   +smatDx(i,j) *vorticity(node,3)
            END DO
          END DO
        END DO
c
c     vz
c
      ELSE IF (eqn.EQ.3) THEN

        DO isqr=1,sysm%neq !nunk ! zanka po vrsticah kvadratne matrike
          inode=mesh%KMeql(isqr)  ! vozlisce v U  mrezi, kateremu pripada ta enacba
c         dodajamo prispevke enack, ki se sestejejo
          DO isuma=1,mesh%sqUlistNO(inode)
            i=mesh%sqUlist(inode,isuma)
            ii=mesh%sqUlistKM(inode,isuma)
            ic=mesh%sqUlistIC(inode,isuma)
            DO j=1,mesh%npoc
              node=mesh%idc(ic,j)
              b(isqr)=b(isqr)
     &                   -smatHty(ii,j)*velocity(node,1)
     &                   +smatHtx(ii,j)*velocity(node,2)
     &                   +smatDy(i,j) *vorticity(node,1)
     &                   -smatDx(i,j) *vorticity(node,2)
            END DO
          END DO
        END DO

      END IF

      END





C -----------------------------------------------------------------------------
      SUBROUTINE solv_kmMacro3(eqn,inp,cpu,mesh,nit,sysm,rhsm,precv,kmp,
     &                         smatHtx,smatHty,smatHtz,
     &                         smatDx,smatDy,smatDz,velocity,vorticity)
C
C     Solve the kinematics equation for domain velocities
c     overdetermined system of equations
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(meshtype) mesh
      TYPE(MATRIX) sysm,rhsm
      TYPE(KMpointer) kmp      
      TYPE(InputType) :: inp   
      TYPE(CPUtype) cpu            
      
      INTEGER ierr,i,nit
      
      REAL(8), ALLOCATABLE :: b(:),x(:),knv(:)
      
      REAL(8) smatHtx(mesh%nicnpoc,mesh%npoc) ! doktorat, enacba (4.7)
      REAL(8) smatHty(mesh%nicnpoc,mesh%npoc) ! to je za prvi integral 
      REAL(8) smatHtz(mesh%nicnpoc,mesh%npoc) ! na desni
      
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc) ! doktorat, enacba (4.7)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc) ! to je za drugi integral 
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc) ! na desni        

      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3)      
      REAL(8) precv(kmp%nsmcol)
      INTEGER eqn 
      REAL(4) cptime,rcpu        

c     alociram spomin za vektor neznank in vektor na desni
      ALLOCATE (knv(kmp%nrhscol)) ! vector of known values
      ALLOCATE (x(kmp%nsmcol))    ! vector of unknown values
      ALLOCATE (b(sysm%neq))    ! r.h.s vector of the system of linear equations b=rhsm * knv

      rcpu=cptime(0.)

C     Set up vectors of known an unknown values      
      CALL SetUpXknv3(eqn,mesh,kmp,velocity,knv,x)

C     Set up r.h.s. vector : b= rhsm * knv
      CALL CRSxV(rhsm,knv,kmp%nrhscol,b)
c
c     b = b + matrike * rhsv
c
      CALL rhsKM3(eqn,mesh,b,smatHtx,smatHty,smatHtz,
     &                  smatDx,smatDy,smatDz,velocity,vorticity)
      cpu%time((eqn-1)*6+3)=cpu%time((eqn-1)*6+3)+cptime(rcpu)
      rcpu=cptime(0.)

c
c      solve
c    
      CALL slvlsqr2(sysm%neq,kmp%nsmcol,sysm%nnz,inp%lsqs_maxit,inp%dlse(eqn),nit,ierr,
     &              precv,sysm%i,sysm%j,sysm%v,b,x)
      cpu%time((eqn-1)*6+4)=cpu%time((eqn-1)*6+4)+cptime(rcpu)

      IF (ierr.NE.0) print *,"KM NAPAKA V SOVLERJU",ierr
         
      DO i=1,kmp%nsmcol
        IF (kmp%bunt(i).EQ.1) THEN  ! vx
          velocity(kmp%bunn(i),eqn)=(1.0D0-inp%urDv(eqn))*velocity(kmp%bunn(i),eqn)+inp%urDv(eqn)*x(i)
        ELSE 
         print *,"error!"
        END IF
      END DO

      
      DEALLOCATE (x,b,knv)

C     Employ dv/dn=0 boundary condition     
      DO i=1,mesh%nbnodes
        IF (mesh%kode(i,eqn).EQ.2.OR.mesh%kode(i,eqn).EQ.3) THEN  ! d(v_eqn)/dn=0
          velocity(mesh%gbn(i),eqn)=velocity(mesh%kobc(i),eqn)
        END IF
      END DO   
      
C     Employ periodic boundary condition
      IF (mesh%iPeri.GT.0) THEN
        DO i=1,mesh%nPeri
          velocity(mesh%pwFT(i,2),eqn)=velocity(mesh%pwFT(i,1),eqn)
        END DO
      END IF

      END

C -----------------------------------------------------------------------------
      SUBROUTINE rhsKM3(eqn,mesh,b,smatHtx,smatHty,smatHtz,
     &                  smatDx,smatDy,smatDz,velocity,vorticity)
C
C      sets up right hand side for kinematics eqaution
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(meshtype) mesh
      
      REAL(8) smatHtx(mesh%nicnpoc,mesh%npoc) ! doktorat, enacba (4.7)
      REAL(8) smatHty(mesh%nicnpoc,mesh%npoc) ! to je za prvi integral 
      REAL(8) smatHtz(mesh%nicnpoc,mesh%npoc) ! na desni
      
      REAL(8) smatDx(mesh%nicnsp,mesh%npoc) ! doktorat, enacba (4.7)
      REAL(8) smatDy(mesh%nicnsp,mesh%npoc) ! to je za drugi integral 
      REAL(8) smatDz(mesh%nicnsp,mesh%npoc) ! na desni        

      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) vorticity(mesh%nnodes,3) 
      
      REAL(8) b(mesh%nicell*mesh%npoc)
      
      INTEGER i,ii,eqn,ic,it,node,j
      
      i=0
c
c     vx
c      
      IF (eqn.EQ.1) THEN
        ii=0      
        DO ic=1,mesh%nicell
          DO it=1,mesh%npoc !mesh%nsp ! 51 izvornih tock v vsaki celici
            i=(ic-1)*mesh%nsp+it ! vrstica v matriki      
            ii=ii+1            
            DO j=1,mesh%npoc
              node=mesh%idc(ic,j)
              b(ii)=b(ii)-smatHtz(ii,j)*velocity(node,2)
     &                   +smatHty(ii,j)*velocity(node,3)
     &                   +smatDz(i,j) *vorticity(node,2)
     &                   -smatDy(i,j) *vorticity(node,3)
            END DO
          END DO
        END DO      
c
c     vy
c      
      ELSE IF (eqn.EQ.2) THEN
        ii=0
        DO ic=1,mesh%nicell
          DO it=1,mesh%npoc !mesh%nsp ! 51 izvornih tock v vsaki celici
            i=(ic-1)*mesh%nsp+it ! vrstica v matriki            
            ii=ii+1
            DO j=1,mesh%npoc
              node=mesh%idc(ic,j)
              b(ii)=b(ii)+smatHtz(ii,j)*velocity(node,1)
     &                   -smatHtx(ii,j)*velocity(node,3)
     &                   -smatDz(i,j) *vorticity(node,1)
     &                   +smatDx(i,j) *vorticity(node,3)
            END DO
          END DO
        END DO      
c
c     vz
c      
      ELSE IF (eqn.EQ.3) THEN
        ii=0
        DO ic=1,mesh%nicell
          DO it=1,mesh%npoc !mesh%nsp ! 51 izvornih tock v vsaki celici
            i=(ic-1)*mesh%nsp+it ! vrstica v matriki           
            ii=ii+1
            DO j=1,mesh%npoc
              node=mesh%idc(ic,j)
              b(ii)=b(ii)-smatHty(ii,j)*velocity(node,1)
     &                   +smatHtx(ii,j)*velocity(node,2)
     &                   +smatDy(i,j) *vorticity(node,1)
     &                   -smatDx(i,j) *vorticity(node,2)
            END DO
          END DO
        END DO      

      END IF

      END


C______________________________________________________________________C
      SUBROUTINE SetUpXknv3(eqn,mesh,kmp,velocity,knv,x)
C
C     Set up vectors of known an unknown values
C
      USE inc_types
      TYPE(meshtype) mesh
      TYPE(KMpointer) kmp      
      
      INTEGER i,eqn   
      
      REAL(8) velocity(mesh%nnodes,3)
      REAL(8) x(kmp%nsmcol)
      REAL(8) knv(kmp%nrhscol)      
            
      DO i=1,mesh%nnodes
      
        IF (kmp%u(i).LT.0) THEN ! vrednost hitrosti neznana
          x(abs(kmp%u(i)))=velocity(i,eqn)
        ELSE ! vrednost hitrosti  znana
          knv(kmp%u(i))=velocity(i,eqn)
        END IF
  
      END DO
      
      END



C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsKM3sq(mesh,smatH,sysm,rhsm,kmp)
C
C     $: Iz pravokotnih matrik h, ht, d naredi CRS sistemsko in rhs
c        sysm in rhsm sta enaki za vse tri enacbe
c        izvorne tocke samo v u nodeih
C
C -----------------------------------------------------------------------------
      USE inc_types

      TYPE(meshType) :: mesh
      INTEGER nunk,nb
      REAL(8) smatH(mesh%nicnsp,mesh%npoc)

      REAL(8), ALLOCATABLE :: FRsys(:),FRrhs(:)

      TYPE(matrix) :: sysm,rhsm
      TYPE(KMpointer) kmp

      INTEGER ic,row,crsrow,col,nuq,j,i,ks,kr
      INTEGER isuma,isqr,inode
c
c     dolzini vektorja neznank in vektorja na desni
c
      nunk=kmp%nsmcol
      nb=kmp%nrhscol
C
C     stevilo enacb
C
      sysm%neq=nunk !mesh%nicell*mesh%npoc !mesh%nicnsp
      rhsm%neq=nunk !mesh%nicell*mesh%npoc !mesh%nicnsp

      ALLOCATE (mesh%KMeql(nunk))
      ALLOCATE (FRsys(nunk),FRrhs(nb))

c     v EQL(i) je stevilka vozlisca, kateremu pripada i-ta enacba
c     zanka po vrstica v matriki
      DO i=1,nunk
        DO j=1,mesh%nnodes ! samo u u izvornih tockah, ker pri nimam q izvornih tock
          IF (kmp%u(j).EQ.-i) THEN
            mesh%KMeql(i)=j
          END IF
        END DO
      END DO

      sysm%nnz=0
      rhsm%nnz=0
      FRsys=0.0D0
      FRsys=0.0D0
      DO isqr=1,nunk ! zanka po vrsticah kvadratne matrike
        inode=ABS(mesh%KMeql(isqr))  ! vozlisce v U mrezi, kateremu pripada ta enacba
        FRsys=0
        FRrhs=0
c       dodajamo prispevke enacb, ki se sestejejo
        DO isuma=1,mesh%sqUlistNO(inode)
          row=mesh%sqUlist(inode,isuma)
          ic=mesh%sqUlistIC(inode,isuma)
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            nuq=kmp%u(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              FRsys(ABS(nuq))=FRsys(ABS(nuq))+smatH(row,col)
            ELSE ! znana vrednost - rhs
              FRrhs(nuq)=FRrhs(nuq)-smatH(row,col)
            END IF
          END DO
        END DO

c       prestejem nnz
        DO i=1,nunk
          IF (FRsys(i).NE.0.0D0) sysm%nnz=sysm%nnz+1
        END DO
        DO i=1,nb
          IF (FRrhs(i).NE.0.0D0) rhsm%nnz=rhsm%nnz+1
        END DO

      END DO
C
C     Alokacija pomnilnika
C
      ALLOCATE(sysm%v(sysm%nnz),sysm%i(sysm%neq+1),sysm%j(sysm%nnz),sysm%d(sysm%neq))
      ALLOCATE(rhsm%v(rhsm%nnz),rhsm%i(rhsm%neq+1),rhsm%j(rhsm%nnz),rhsm%d(rhsm%neq))
C
C     Polnim matriki
C
      ks=0
      kr=0
      crsrow=0   ! vrstica v sistemski in matrikah integralov
      DO isqr=1,nunk ! zanka po vrsticah kvadratne matrike
        inode=ABS(mesh%KMeql(isqr))  ! vozlisce v U mrezi, kateremu pripada ta enacba
        FRsys=0
        FRrhs=0
c       dodajamo prispevke enacb, ki se sestejejo
        DO isuma=1,mesh%sqUlistNO(inode)
          row=mesh%sqUlist(inode,isuma)
          ic=mesh%sqUlistIC(inode,isuma)
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H
            nuq=kmp%u(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              FRsys(ABS(nuq))=FRsys(ABS(nuq))+smatH(row,col)
            ELSE ! znana vrednost - rhs
              FRrhs(nuq)=FRrhs(nuq)-smatH(row,col)  ! minus zato, ker gre na drugo stran enacbe
            END IF
          END DO
        END DO
C
C       predelam v CRS
C
C       zacetek vrste
        crsrow=crsrow+1
        sysm%i(crsrow)=ks+1
        rhsm%i(crsrow)=kr+1
C       sys
        DO i=1,nunk
          IF (FRsys(i).NE.0.0D0) THEN
            ks=ks+1
            sysm%v(ks)=FRsys(i)
            sysm%j(ks)=i
            IF (i.EQ.isqr) sysm%d(isqr)=ks ! postavimo diagonalo
          END IF
        END DO
C       rhs
        DO i=1,nb
          IF (FRrhs(i).NE.0.0D0) THEN
            kr=kr+1
            rhsm%v(kr)=FRrhs(i)
            rhsm%j(kr)=i
          END IF
        END DO
      END DO
c     zadnji clen, zacetek neobstojece vrste
      sysm%i(crsrow+1)=ks+1
      rhsm%i(crsrow+1)=kr+1

      DEALLOCATE (FRsys,FRrhs)

      END




C -----------------------------------------------------------------------------
      SUBROUTINE sMat2crsSysRhsKM3(mesh,smatH,sysm,rhsm,kmp)      
C
C     $: Iz pravokotnih matrik h, ht, d naredi CRS sistemsko in rhs
c        sysm in rhsm sta enaki za vse tri enacbe
c        izvorne tocke samo v u nodeih
C
C -----------------------------------------------------------------------------
      USE inc_types 

      TYPE(meshType) :: mesh
      INTEGER nx,nb
      REAL(8) smatH(mesh%nicnsp,mesh%npoc)
     
      REAL(8), ALLOCATABLE :: FRsys(:),FRrhs(:)       
      
      TYPE(matrix) :: sysm,rhsm
      TYPE(KMpointer) kmp      
            
      INTEGER ic,it,row,crsrow,col,nuq,j,i,ks,kr
c
c     dolzini vektorja neznank in vektorja na desni
c
      nx=kmp%nsmcol
      nb=kmp%nrhscol
C      
C     stevilo enacb
C
      sysm%neq=mesh%nicell*mesh%npoc !mesh%nicnsp      
      rhsm%neq=mesh%nicell*mesh%npoc !mesh%nicnsp
C
C     stevilo nenicelnih clenov
C      
      sysm%nnz=0
      rhsm%nnz=0
      DO ic=1,mesh%nicell ! zanka po celicah
        DO it=1,mesh%npoc !mesh%nsp ! po izvornih tockah znotraj celice
          row=(ic-1)*mesh%nsp+it ! vrstica v sistemski matriki                  
          ! preskakuje, ker imam v SmatH tudi q izvorne tocke
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H = zaporedna stevilka v vektroju znank oz neznank
c           ** vx ali vy ali vz **
            nuq=kmp%u(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              IF (smatH(row,col).NE.0.0D0) sysm%nnz=sysm%nnz+1           
            ELSE ! znana vrednost - rhs
              IF (smatH(row,col).NE.0.0D0) rhsm%nnz=rhsm%nnz+1             
            END IF          
          END DO             
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
      ALLOCATE (FRsys(nx),FRrhs(nb))          
      ks=0
      kr=0
      crsrow=0   ! vrstica v sistemski in matrikah integralov = 3*vec
      DO ic=1,mesh%nicell ! zanka po celicah
        DO it=1,mesh%npoc !mesh%nsp ! po izvornih tockah znotraj celice
          row=(ic-1)*mesh%nsp+it ! vrstica v matrikah integralov
          FRsys=0.0D0
          FRrhs=0.0D0
C
C         obravnavano hkrati vse tri enacbe
C
          DO col=1,mesh%npoc ! zanka po stolpcih pravokotne matrike
            j=mesh%idc(ic,col) ! stolpec za H = zaporedna stevilka v vektroju znank oz neznank
c           ** vx **
            nuq=kmp%u(j) ! stevilka v vektorju neznank oziroma znank
            IF (nuq.LT.0) THEN ! neznaka - sys
              FRsys(ABS(nuq))= smatH(row,col)            
            ELSE ! znana vrednost - rhs
              FRrhs(nuq)=-smatH(row,col)  ! minus zato, ker gre na drugo stran enacbe             
            END IF
          END DO            
C               
C         predelam v CRS
C   
C         zacetek vrste
          crsrow=crsrow+1
          sysm%i(crsrow)=ks+1
          rhsm%i(crsrow)=kr+1          
C         sys
          DO i=1,nx
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
        END DO
      END DO     
c     zadnji clen, zacetek neobstojece vrste
      sysm%i(crsrow+1)=ks+1
      rhsm%i(crsrow+1)=kr+1          
      
      DEALLOCATE (FRsys,FRrhs)
      
      END


C______________________________________________________________________C
      SUBROUTINE DetKMunk3(mesh,kmp)
C
C     Determine kinematics unknowns, each equation separetely
C
      USE inc_types
      TYPE(meshtype) mesh
      TYPE(KMpointer) kmp
      
      INTEGER i,j

c     naredim pointerje za vx,vy,vz,wx,wy,wz kateri je znan, kateri neznan
      ALLOCATE (kmp%u(mesh%nnodes))
    
      kmp%nrhscol=0
      kmp%nsmcol=0

c     pointer na stevilko neznanke / znanke
c     hkrati stejem stevilo stoplcev v sistemski in na rhs
c     torej v kmp%vx(i) je zaporedna stevilka v vektorju 
c     (na desni ce +, na levi ce -) hitrosti vx v i-tem vozliscu

c
c     treba bi bilo upostevati robne pogoje, zaenkrat imam
c     hitrosti na robu znane
c     uporabim enega za vse tri enacbr

c     vx, vy, vz
      DO i=1,mesh%nnodes
        IF (mesh%lbn(i).NE.0) THEN ! tocka je na robu -> znana
          kmp%nrhscol=kmp%nrhscol+1
          kmp%u(i)=kmp%nrhscol
        ELSE ! tocka je v notranjosti -> neznanka             
          kmp%nsmcol=kmp%nsmcol+1
          kmp%u(i)=-kmp%nsmcol
        END IF
      END DO      
      
c     naredim se povraten indeks - iz neznanega vektorja v velocity, vorcitiy      
      ALLOCATE (kmp%bunt(kmp%nsmcol)) !  1=u
      ALLOCATE (kmp%bunn(kmp%nsmcol)) !  node number of the unkonwn      
      
      DO j=1,kmp%nsmcol
        DO i=1,mesh%nnodes    
          IF (kmp%u(i).EQ.-j) THEN
            kmp%bunt(j)=1
            kmp%bunn(j)=i
          END IF
        END DO
      END DO           
      
      END     

