C -----------------------------------------------------------------------------
      SUBROUTINE SetUpMagneticField(mesh,inp,mfs,hshg)
C
C     $: sets up magnetic field
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(InputType) :: inp
      TYPE(meshType) :: mesh

      REAL(8) mfs(mesh%nnodes,3) ! H
      REAL(8), ALLOCATABLE :: hsh(:) ! H skarano H
      REAL(8) hshg(mesh%nnodes,3) ! H skarano H

      INTEGER i
C     Magnetno polje tu je brez dimenzij
C     Deljeno je z H_0, ki je upostevan v magnetic Rayleigh number

      ALLOCATE (hsh(mesh%nnodes))

      mfs=0.0D0

c     Polje v smeri z
      IF (inp%iferro.EQ.1) THEN
        DO i=1,mesh%nnodes
          mfs(i,1)=1.0D0
        END DO
      ELSE IF (inp%iferro.EQ.2) THEN
        DO i=1,mesh%nnodes
          mfs(i,1)=mesh%x(i,1)
        END DO
      END IF

c     Izracunam H \cdot H
      DO i=1,mesh%nnodes
        hsh(i)=mfs(i,1)**2+mfs(i,2)**2+mfs(i,3)**2
      END DO

c     Izracunam gradient hsh
      CALL setGrad(mesh,hsh,hshg)

      DEALLOCATE (hsh)

      END

C -----------------------------------------------------------------------------
      SUBROUTINE CalKelvinBodyForce(mesh,inp,temp,mfs,hshg,kbf,tempGrad)
C
C     $: sets up magnetic field
C
C -----------------------------------------------------------------------------
      USE inc_types
      TYPE(InputType) :: inp
      TYPE(meshType) :: mesh

      REAL(8) mfs(mesh%nnodes,3) ! H
      REAL(8) kbf(mesh%nnodes,3) ! Kelvin Body Force
      REAL(8) hshg(mesh%nnodes,3) ! gradient H skarano H
      REAL(8) temp(mesh%nnodes) ! Temperature
      REAL(8) tempGrad(mesh%nnodes,3) ! Temperature

      REAL(8), ALLOCATABLE :: hsgt (:) ! H skalarno grad T

      REAL(8) enap
      INTEGER i,j


c     Izracunam H skalarno grad T
      ALLOCATE (hsgt(mesh%nnodes))
      DO i=1,mesh%nnodes
        hsgt(i)=mfs(i,1)*tempGrad(i,1)+mfs(i,2)*tempGrad(i,2)+mfs(i,3)*tempGrad(i,3)
      END DO

c     Izracunam \vec f'_K = kbf
      DO i=1,mesh%nnodes
        enap=1.0D0+temp(i)*inp%bdt
        DO j=1,3
          kbf(i,j)=0.5D0/(enap**2)*( (1.0D0+inp%hi0)/inp%bdt + temp(i))* hshg(i,j)+
     &             inp%hi0/(enap**3)*mfs(i,j)*hsgt(i)
        END DO
      END DO

      DEALLOCATE (hsgt)

      END
