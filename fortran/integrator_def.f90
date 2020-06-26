
! integrator_def.f90
!
!>  Base class definition for the model's integrators.
!
!> @copyright
!> 2020 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.
!
!---------------------------------------------------------------------------!
!
!
!---------------------------------------------------------------------------

MODULE integrator_def
  USE model_def

  IMPLICIT NONE

  PRIVATE

  !> Base class to be subclassed to create a new integrator.
  TYPE, ABSTRACT, PUBLIC :: Integrator
    TYPE(Model), POINTER :: pmodel      !< A pointer to the model to integrate.
    REAL(KIND=8), POINTER :: dt         !< Time step of the integrator.
    INTEGER, POINTER :: ndim            !< Dimension of the phase space of the model to integrate.
  CONTAINS
    PROCEDURE(init_int), DEFERRED :: init
    PROCEDURE(step_int), DEFERRED :: step
    PROCEDURE(clean_int), DEFERRED :: clean
  END TYPE Integrator

  ABSTRACT INTERFACE
    !> Abstract interface for the procedures initializing the integrator objects.
    !> @param[in,out] integr Integrator object to initialize.
    !> @param[in] imodel Model object to initialize the integrator with.
    SUBROUTINE init_int(integr, imodel)
      IMPORT Integrator, Model
      CLASS(Integrator), INTENT(INOUT) :: integr
      CLASS(Model), INTENT(IN), TARGET :: imodel
    END SUBROUTINE init_int
  END INTERFACE

  ABSTRACT INTERFACE
    !> Abstract interface for the procedure to make the integrator compute a model's time step.
    !> @param[in,out] integr Integrator object to perform the step with.
    !> @param[in] y Initial point.
    !> @param[in] t Actual integration time
    !> @param[out] res Final point after the step.
    SUBROUTINE step_int(integr, y,t,res)
      IMPORT Integrator
      CLASS(Integrator), INTENT(INOUT) :: integr
      REAL(KIND=8), DIMENSION(0:integr%ndim), INTENT(IN) :: y
      REAL(KIND=8), INTENT(INOUT) :: t
      REAL(KIND=8), DIMENSION(0:integr%ndim), INTENT(OUT) :: res
    END SUBROUTINE step_int
  END INTERFACE

  ABSTRACT INTERFACE
    !> Abstract interface for the procedure to clean the integrator objects.
    !> @param[in,out] integr Integrator object to clean.
    SUBROUTINE clean_int(integr)
      IMPORT Integrator
      CLASS(Integrator), INTENT(INOUT) :: integr
    END SUBROUTINE clean_int
  END INTERFACE

END MODULE integrator_def
