
! aotensor_def.f90
!
!>  The equation tensor \f$\mathcal{T}_{i,j,k}\f$ for the coupled ocean-atmosphere model
!>  with temperature which allows for an extensible set of modes
!>  in the ocean and in the atmosphere.
!
!> @copyright                                                               
!> 2015-2020 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!                                                                           
!---------------------------------------------------------------------------!


MODULE aotensor_def

  !-----------------------------------------------------!
  !                                                     !
  ! Preamble and variables declaration                  !
  !                                                     !
  !-----------------------------------------------------!

  USE params
  USE inprod_analytic
  USE tensor_def
  IMPLICIT NONE

  PRIVATE

  ! Epsilon to test equality with 0
  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16

  !> Class to hold the tensor \f$\mathcal{T}_{i,j,k}\f$ representation of the tendencies.
  TYPE, PUBLIC :: AtmOcTensor
    TYPE(Tensor) :: tensor                                  !< The tensor object
    INTEGER, DIMENSION(:), ALLOCATABLE :: count_elems       !< A list of the number of non-zero entries of the tensor component along \f$i\f$.
    PROCEDURE(ao_coeff), PRIVATE, POINTER :: operation
    LOGICAL :: initialized
    INTEGER, POINTER :: noc, natm, ndim
  CONTAINS
    PROCEDURE :: init => init_aotensor
    PROCEDURE :: clean => delete_aotensor
    PROCEDURE, PRIVATE :: compute_tensor => compute_aotensor
    PROCEDURE, PRIVATE :: psi, theta, A, T
  END TYPE AtmOcTensor

  !-----------------------------------------------------!
  !                                                     !
  ! End of preamble                                     !
  !                                                     !
  !-----------------------------------------------------!

CONTAINS

  !-----------------------------------------------------!
  !                                                     !
  ! Function declarations                               !
  !                                                     !
  !-----------------------------------------------------!

  ! Translate the \f$\psi_{a,i}\f$ coefficients into effective coordinates
  FUNCTION psi(aot, i)
    CLASS(AtmOcTensor), INTENT(IN) :: aot
    INTEGER :: i,psi
    psi = i
  END FUNCTION psi

  ! Translate the \f$\theta_{a,i}\f$ coefficients into effective coordinates
  FUNCTION theta(aot,i)
    CLASS(AtmOcTensor), INTENT(IN) :: aot
    INTEGER :: i,theta
    theta = i + aot%natm
  END FUNCTION theta

  ! Translate the \f$\psi_{o,i}\f$ coefficients into effective coordinates
  FUNCTION A(aot,i)
    CLASS(AtmOcTensor), INTENT(IN) :: aot
    INTEGER :: i,A
    A = i + 2 * aot%natm
  END FUNCTION A

  ! Translate the \f$\delta T_{o,i}\f$ coefficients into effective coordinates
  FUNCTION T(aot,i)
    CLASS(AtmOcTensor), INTENT(IN) :: aot
    INTEGER :: i,T
    T = i + 2 * aot%natm + aot%noc
  END FUNCTION T

  ! Kronecker delta function
  FUNCTION kdelta(i,j)
    INTEGER :: i,j,kdelta
    kdelta=0
    IF (i == j) kdelta = 1
  END FUNCTION kdelta

  ! Subroutine to add element in the aotensor \f$\mathcal{T}_{i,j,k}\f$ structure.
  SUBROUTINE ao_coeff(aot,i,j,k,v)
    CLASS(AtmOcTensor), INTENT(INOUT) :: aot
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN) :: v
    INTEGER :: n
    IF (.NOT. allocated(aot%tensor%t)) THEN
      PRINT*, "Warning: Trying to compute an aotensor not previously allocated."
      PRINT*, "Aborting aotensor initialization."
      aot%initialized = .FALSE.
      RETURN
    END IF
    IF (.NOT. allocated(aot%tensor%t(i)%elems)) THEN
      PRINT*, "Warning: Trying to compute an aotensor not previously allocated."
      PRINT*, "Aborting aotensor initialization."
      aot%initialized = .FALSE.
      RETURN
    END IF
    IF (abs(v) .ge. real_eps) THEN
       n=(aot%tensor%t(i)%nelems)+1
       IF (j .LE. k) THEN
          aot%tensor%t(i)%elems(n)%j=j
          aot%tensor%t(i)%elems(n)%k=k
       ELSE
          aot%tensor%t(i)%elems(n)%j=k
          aot%tensor%t(i)%elems(n)%k=j
       END IF
       aot%tensor%t(i)%elems(n)%v=v
       aot%tensor%t(i)%nelems=n
    END IF
  END SUBROUTINE ao_coeff

  ! Subroutine to count the elements of the aotensor \f$\mathcal{T}_{i,j,k}\f$. Add +1 to count_elems(i) for each value that is added to the tensor i-th component.
  SUBROUTINE ao_add_count(aot,i,j,k,v)
    CLASS(AtmOcTensor), INTENT(INOUT) :: aot
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN)  :: v
    IF (abs(v) .ge. real_eps) aot%count_elems(i)=aot%count_elems(i)+1
  END SUBROUTINE ao_add_count

  ! Subroutine to compute the tensor aotensor
  SUBROUTINE compute_aotensor(aot, model_config, inprods)
    CLASS(AtmOcTensor), INTENT(INOUT) :: aot
    CLASS(ModelConfiguration), INTENT(IN), TARGET :: model_config
    CLASS(InnerProducts), INTENT(IN), TARGET :: inprods

    TYPE(PhysicsConfiguration), POINTER :: phys
    TYPE(InnerProducts), POINTER :: ips

    INTEGER :: i,j,k

    phys => model_config%physics
    ips => inprods

    CALL aot%operation(aot%theta(1),0,0,(phys%Cpa / (1 - ips%atmos%a(1,1) * phys%sig0)))
    DO i = 1, aot%natm
       DO j = 1, aot%natm
          CALL aot%operation(aot%psi(i),aot%psi(j),0,-(((ips%atmos%c(i,j) * phys%betp) / ips%atmos%a(i,i))) -&
               &(phys%kd * kdelta(i,j)) / 2 + ips%atmos%a(i,j)*phys%nuap)
          CALL aot%operation(aot%theta(i),aot%psi(j),0,(ips%atmos%a(i,j) * phys%kd * phys%sig0) &
               & / (-2 + 2 * ips%atmos%a(i,i) * phys%sig0))
          CALL aot%operation(aot%psi(i),aot%theta(j),0,(phys%kd * kdelta(i,j)) / 2)
          CALL aot%operation(aot%theta(i),aot%theta(j),0,(-((phys%sig0 * (2. * ips%atmos%c(i,j) * phys%betp +&
               & ips%atmos%a(i,j) * (phys%kd + 4. * phys%kdp)))) + 2. * (phys%LSBpa + phys%sc * phys%Lpa) &
               &* kdelta(i,j)) / (-2. + 2. * ips%atmos%a(i,i) * phys%sig0))
          DO k = 1, aot%natm
             CALL aot%operation(aot%psi(i),aot%psi(j),aot%psi(k),-((ips%atmos%b(i,j,k) / ips%atmos%a(i,i))))
             CALL aot%operation(aot%psi(i),aot%theta(j),aot%theta(k),-((ips%atmos%b(i,j,k) / ips%atmos%a(i,i))))
             CALL aot%operation(aot%theta(i),aot%psi(j),aot%theta(k),(ips%atmos%g(i,j,k) -&
                  & ips%atmos%b(i,j,k) * phys%sig0) / (-1 + ips%atmos%a(i,i) *&
                  & phys%sig0))
             CALL aot%operation(aot%theta(i),aot%theta(j),aot%psi(k),(ips%atmos%b(i,j,k) * phys%sig0) &
                  & / (1 - ips%atmos%a(i,i) * phys%sig0))
          END DO
       END DO
       DO j = 1, aot%noc
          CALL aot%operation(aot%psi(i),aot%A(j),0,phys%kd * ips%atmos%d(i,j) / (2 * ips%atmos%a(i,i)))
          CALL aot%operation(aot%theta(i),aot%A(j),0,phys%kd * (ips%atmos%d(i,j) * phys%sig0) &
                & / (2 - 2 * ips%atmos%a(i,i) * phys%sig0))
          CALL aot%operation(aot%theta(i),aot%T(j),0,ips%atmos%s(i,j) * (2 * phys%LSBpo + phys%Lpa) &
                & / (2 - 2 * ips%atmos%a(i,i) * phys%sig0))
       END DO
    END DO
    DO i = 1, aot%noc
       DO j = 1, aot%natm
          CALL aot%operation(aot%A(i),aot%psi(j),0,ips%ocean%K(i,j) * phys%dp / (ips%ocean%M(i,i) + phys%G))
          CALL aot%operation(aot%A(i),aot%theta(j),0,-(ips%ocean%K(i,j)) * phys%dp / (ips%ocean%M(i,i) + phys%G))
       END DO
       DO j = 1, aot%noc
          CALL aot%operation(aot%A(i),aot%A(j),0,-((ips%ocean%N(i,j) * phys%betp &
               & + ips%ocean%M(i,i) * (phys%rp + phys%dp) * kdelta(i,j)&
               & - ips%ocean%M(i,j)**2*phys%nuop)) / (ips%ocean%M(i,i) + phys%G))
          DO k = 1, aot%noc
             CALL aot%operation(aot%A(i),aot%A(j),aot%A(k),-(ips%ocean%C(i,j,k)) / (ips%ocean%M(i,i) + phys%G))
          END DO
       END DO
    END DO
    DO i = 1, aot%noc
       CALL aot%operation(aot%T(i),0,0,phys%Cpo * ips%ocean%W(i,1))
       DO j = 1, aot%natm
          CALL aot%operation(aot%T(i),aot%theta(j),0,ips%ocean%W(i,j) * (2 * phys%sc * phys%Lpo + phys%sBpa))
       END DO
       DO j = 1, aot%noc
          CALL aot%operation(aot%T(i),aot%T(j),0,-((phys%Lpo + phys%sBpo)) * kdelta(i,j))
          DO k = 1, aot%noc
             CALL aot%operation(aot%T(i),aot%A(j),aot%T(k),-(ips%ocean%O(i,j,k)))
          END DO
       END DO
    END DO
  END SUBROUTINE compute_aotensor

  !-----------------------------------------------------!
  !                                                     !
  ! Initialisation routine                              !
  !                                                     !
  !-----------------------------------------------------!

  !> Subroutine to initialise the AtmOcTensor tensor.
  !> @param[in,out] aot The AO tensor object to initialize.
  !> @param[in] model_configuration A model configuration object to initialize the model tensor with.
  !> @param[in] inprods A model inner products object to initialize the model with.
  SUBROUTINE init_aotensor(aot, model_configuration, inprods)
    CLASS(AtmOcTensor), INTENT(INOUT) :: aot
    CLASS(ModelConfiguration), INTENT(IN), TARGET :: model_configuration
    CLASS(InnerProducts), INTENT(IN), TARGET :: inprods

    INTEGER :: i
    INTEGER :: AllocStat 

    IF (.NOT.model_configuration%initialized) THEN
      PRINT*, "Warning: Model configuration not initialized."
      PRINT*, "Aborting aotensor initialization."
      RETURN
    END IF

    IF (.NOT.inprods%initialized) THEN
      PRINT*, "Warning: Inner products not initialized."
      PRINT*, "Aborting aotensor initialization."
      RETURN
    END IF

    aot%ndim => model_configuration%modes%ndim
    aot%natm => model_configuration%modes%natm
    aot%noc => model_configuration%modes%noc

    ALLOCATE(aot%count_elems(aot%ndim), STAT=AllocStat)
    IF (AllocStat /= 0) THEN
      PRINT*, "*** init_aotensor: Problem with allocation! ***"
      STOP "Exiting ..."
    END IF
    aot%count_elems=0

    CALL aot%tensor%init(aot%ndim)

    aot%operation => ao_add_count
    CALL aot%compute_tensor(model_configuration, inprods)

    DO i=1,aot%ndim
      ALLOCATE(aot%tensor%t(i)%elems(aot%count_elems(i)), STAT=AllocStat)
      IF (AllocStat /= 0) THEN
        PRINT*, "*** init_aotensor: Problem with allocation! ***"
      STOP "Exiting ..."
    END IF

    END DO

    aot%operation => ao_coeff
    CALL aot%compute_tensor(model_configuration, inprods)

    CALL aot%tensor%simplify

    aot%initialized = .TRUE.

  END SUBROUTINE init_aotensor

  !> Subroutine to clean a AtmOcTensor tensor.
  !> @param[in,out] aot The AtmOcTensor tensor object to initialize.
  SUBROUTINE delete_aotensor(aot)
    CLASS(AtmOcTensor), INTENT(INOUT) :: aot

    IF (allocated(aot%count_elems)) DEALLOCATE(aot%count_elems)

    CALL aot%tensor%clean
    NULLIFY(aot%ndim)
    NULLIFY(aot%natm)
    NULLIFY(aot%noc)

    aot%initialized = .FALSE.

  END SUBROUTINE delete_aotensor

END MODULE aotensor_def
      


