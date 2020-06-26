
! tl_ad_tensor.f90
!
!> Tangent Linear (TL) and Adjoint (AD) model versions of MAOOAM.
!> Tensors definition module
!
!> @copyright                                                               
!> 2016-2020 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!                                                                           
!---------------------------------------------------------------------------!

MODULE tl_ad_tensor

  !-----------------------------------------------------!
  !                                                     !
  ! Preamble and variables declaration                  !
  !                                                     !
  !-----------------------------------------------------!

  USE aotensor_def
  USE tensor_def
  IMPLICIT NONE

  PRIVATE

  ! Epsilon to test equality with 0
  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16

  !> Tensor representation of the Tangent Linear tendencies.
  TYPE, PUBLIC :: TlTensor
    TYPE(Tensor) :: tensor                                          !< The TL tensor object
    INTEGER, DIMENSION(:), ALLOCATABLE :: count_elems               !< A list of the number of non-zero entries of the tensor component along the first index.
    PROCEDURE(tl_coeff), PRIVATE, POINTER :: tl_operation
    LOGICAL, PUBLIC :: initialized
  CONTAINS
    PROCEDURE :: init => init_tltensor
    PROCEDURE :: clean => delete_tensor
    PROCEDURE, PRIVATE :: compute_tensor => compute_tltensor
  END TYPE TlTensor

  !> Tensor representation of the Adjoint tendencies.
  TYPE, EXTENDS(TlTensor), PUBLIC :: AdTensor
    PROCEDURE(ad_coeff), PRIVATE, POINTER :: ad_operation
  CONTAINS
    PROCEDURE :: init => init_adtensor
    PROCEDURE, PRIVATE :: compute_tensor => compute_adtensor
  END TYPE AdTensor

  !-----------------------------------------------------!
  !                                                     !
  ! End of preamble                                     !
  !                                                     !
  !-----------------------------------------------------!

CONTAINS

  !> Subroutine to clean a TL tensor
  !> @param[in,out] tens The tensor to clean.
  SUBROUTINE delete_tensor(tens)
    CLASS(TlTensor), INTENT(INOUT) :: tens

    IF (allocated(tens%count_elems)) DEALLOCATE(tens%count_elems)

    CALL tens%tensor%clean

    tens%initialized = .FALSE.

  END SUBROUTINE delete_tensor


  !-----------------------------------------------------!
  !                                                     !
  ! Tangent linear model functions                      !
  !                                                     !
  !-----------------------------------------------------!
  
  !> Subroutine to initialise the TL tensor
  !> @param[in,out] tens The tensor to clean.
  !> @param[in] aot A Atmosphere-Ocean tensor to initialize the TL tensor with.
  SUBROUTINE init_tltensor(tens, aot)
    CLASS(TlTensor), INTENT(INOUT) :: tens
    CLASS(AtmOcTensor), INTENT(IN), TARGET :: aot

    INTEGER :: i
    INTEGER, POINTER :: ndim
    INTEGER :: AllocStat

    IF (.NOT.aot%initialized) THEN
      PRINT*, 'Provided AO tensor not yet initialized, impossible to initialize TL tensor!'
      RETURN
    END IF

    ndim => aot%ndim

    ALLOCATE(tens%count_elems(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) THEN
        PRINT*, "*** init_tltensor: Problem with allocation! ***"
        STOP "Exiting ..."
    END IF
    tens%count_elems=0

    CALL tens%tensor%init(ndim)

    tens%tl_operation => tl_add_count
    CALL tens%compute_tensor(aot)

    DO i=1,ndim
      ALLOCATE(tens%tensor%t(i)%elems(tens%count_elems(i)), STAT=AllocStat)
      IF (AllocStat /= 0) THEN
        PRINT*, "*** init_tltensor: Problem with allocation! ***"
        STOP "Exiting ..."
      END IF
    END DO

    tens%tl_operation => tl_coeff
    CALL tens%compute_tensor(aot)

    CALL tens%tensor%simplify

    tens%initialized = .TRUE.

  END SUBROUTINE init_tltensor

  ! Subroutine to compute the TL tensor
  SUBROUTINE compute_tltensor(tens, aot)
    CLASS(TlTensor), INTENT(INOUT) :: tens
    CLASS(AtmOcTensor), INTENT(IN), TARGET :: aot
    INTEGER :: i,j,k,n,nj
    INTEGER, POINTER :: ndim
    REAL(KIND=8) :: v
    ndim => aot%ndim
    DO i=1,ndim
      nj=aot%tensor%t(i)%nelems
      DO n=1,nj
        j=aot%tensor%t(i)%elems(n)%j
        k=aot%tensor%t(i)%elems(n)%k
        v=aot%tensor%t(i)%elems(n)%v
        CALL tens%tl_operation(i,j,k,v)
      ENDDO
    ENDDO
  END SUBROUTINE compute_tltensor

  ! Subroutine used to count the non-zero TL tensor entries
  SUBROUTINE tl_add_count(tens,i,j,k,v)
    CLASS(TlTensor), INTENT(INOUT) :: tens
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN)  :: v
    IF (abs(v) .ge. real_eps) THEN
      IF (j /= 0) tens%count_elems(i)=tens%count_elems(i)+1
      IF (k /= 0) tens%count_elems(i)=tens%count_elems(i)+1
    ENDIF
  END SUBROUTINE tl_add_count

  ! Subroutine used to compute the TL tensor entries
  SUBROUTINE tl_coeff(tens,i,j,k,v)
    CLASS(TlTensor), INTENT(INOUT) :: tens
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN) :: v
    INTEGER :: n
    IF (.NOT. allocated(tens%tensor%t)) THEN
      PRINT*, "Warning: Trying to compute a tltensor not previously allocated."
      PRINT*, "Aborting aotensor initialization."
      tens%initialized = .FALSE.
      RETURN
    END IF
    IF (.NOT. allocated(tens%tensor%t(i)%elems)) THEN
      PRINT*, "Warning: Trying to compute a tltensor not previously allocated."
      PRINT*, "Aborting tltensor initialization."
      tens%initialized = .FALSE.
      RETURN
    END IF
    IF (abs(v) .ge. real_eps) THEN
      IF (j /=0) THEN
        n=(tens%tensor%t(i)%nelems)+1
        tens%tensor%t(i)%elems(n)%j=j
        tens%tensor%t(i)%elems(n)%k=k
        tens%tensor%t(i)%elems(n)%v=v
        tens%tensor%t(i)%nelems=n
      END IF
      IF (k /=0) THEN
        n=(tens%tensor%t(i)%nelems)+1
        tens%tensor%t(i)%elems(n)%j=k
        tens%tensor%t(i)%elems(n)%k=j
        tens%tensor%t(i)%elems(n)%v=v
        tens%tensor%t(i)%nelems=n
      END IF
    END IF
  END SUBROUTINE tl_coeff

  !-----------------------------------------------------!
  !                                                     !
  ! Adjoint model functions                             !
  !                                                     !
  !-----------------------------------------------------!

  !> Subroutine to initialise the AD tensor
  !> @param[in,out] tens The tensor to clean.
  !> @param[in] aot A Atmosphere-Ocean tensor to initialize the AD tensor with.
  SUBROUTINE init_adtensor(tens, aot)
    CLASS(AdTensor), INTENT(INOUT) :: tens
    CLASS(AtmOcTensor), INTENT(IN), TARGET :: aot

    INTEGER :: i
    INTEGER, POINTER :: ndim
    INTEGER :: AllocStat

    IF (.NOT.aot%initialized) THEN
      PRINT*, 'Provided AO tensor not yet initialized, impossible to initialize AD tensor!'
      RETURN
    END IF

    ndim => aot%ndim

    ALLOCATE(tens%count_elems(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) THEN
        PRINT*, "*** init_adtensor: Problem with allocation! ***"
        STOP "Exiting ..."
    END IF
    tens%count_elems=0

    CALL tens%tensor%init(ndim)

    tens%ad_operation => ad_add_count
    CALL tens%compute_tensor(aot)

    DO i=1,ndim
      ALLOCATE(tens%tensor%t(i)%elems(tens%count_elems(i)), STAT=AllocStat)
      IF (AllocStat /= 0) THEN
        PRINT*, "*** init_adtensor: Problem with allocation! ***"
        STOP "Exiting ..."
      END IF
    END DO

    tens%ad_operation => ad_coeff
    CALL tens%compute_tensor(aot)

    CALL tens%tensor%simplify

    tens%initialized = .TRUE.

  END SUBROUTINE init_adtensor

  ! Subroutine to compute the AD tensor
  SUBROUTINE compute_adtensor(tens, aot)
    CLASS(AdTensor), INTENT(INOUT) :: tens
    CLASS(AtmOcTensor), INTENT(IN), TARGET :: aot
    INTEGER :: i,j,k,n,nj
    INTEGER, POINTER :: ndim
    REAL(KIND=8) :: v
    ndim => aot%ndim
    DO i=1,ndim
      nj=aot%tensor%t(i)%nelems
      DO n=1,nj
        j=aot%tensor%t(i)%elems(n)%j
        k=aot%tensor%t(i)%elems(n)%k
        v=aot%tensor%t(i)%elems(n)%v
        CALL tens%ad_operation(i,j,k,v)
      ENDDO
    ENDDO
  END SUBROUTINE compute_adtensor

  ! Subroutine used to count the number of non-zero AD tensor entries
  SUBROUTINE ad_add_count(tens,i,j,k,v)
    CLASS(AdTensor), INTENT(INOUT) :: tens
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN)  :: v
    IF ((abs(v) .ge. real_eps).AND.(i /= 0)) THEN
      IF (k /= 0) tens%count_elems(k)=tens%count_elems(k)+1
      IF (j /= 0) tens%count_elems(j)=tens%count_elems(j)+1
    ENDIF
  END SUBROUTINE ad_add_count

  ! Subroutine used to compute the AD tensor entries
  SUBROUTINE ad_coeff(tens,i,j,k,v)
    CLASS(AdTensor), INTENT(INOUT) :: tens
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN) :: v
    INTEGER :: n
    IF (.NOT. allocated(tens%tensor%t)) THEN
      PRINT*, "Warning: Trying to compute a tltensor not previously allocated."
      PRINT*, "Aborting aotensor initialization."
      tens%initialized = .FALSE.
      RETURN
    END IF

    IF ((abs(v) .ge. real_eps).AND.(i /=0)) THEN
      IF (k /=0) THEN
        IF (.NOT. allocated(tens%tensor%t(k)%elems)) THEN
          PRINT*, "Warning: Trying to compute a tltensor not previously allocated."
          PRINT*, "Aborting tltensor initialization."
          tens%initialized = .FALSE.
          RETURN
        END IF
        n=tens%tensor%t(k)%nelems+1
        tens%tensor%t(k)%elems(n)%j=i
        tens%tensor%t(k)%elems(n)%k=j
        tens%tensor%t(k)%elems(n)%v=v
        tens%tensor%t(k)%nelems=n
      END IF
      IF (j /=0) THEN
        IF (.NOT. allocated(tens%tensor%t(j)%elems)) THEN
          PRINT*, "Warning: Trying to compute a tltensor not previously allocated."
          PRINT*, "Aborting tltensor initialization."
          tens%initialized = .FALSE.
          RETURN
        END IF

        n=(tens%tensor%t(j)%nelems)+1
        tens%tensor%t(j)%elems(n)%j=i
        tens%tensor%t(j)%elems(n)%k=k
        tens%tensor%t(j)%elems(n)%v=v
        tens%tensor%t(j)%nelems=n
      END IF
    END IF
  END SUBROUTINE ad_coeff
END MODULE tl_ad_tensor
