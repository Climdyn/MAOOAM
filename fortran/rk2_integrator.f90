
! rk2_integrator.f90
!
!>  Module containing the second-order Runge-Kutta (RK2) integration routines.
!
!> @copyright                                                               
!> 2015-2020 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!
!>  @remark
!>  This module actually contains the Heun algorithm routines.
!                                                                           
!---------------------------------------------------------------------------

MODULE rk2_integrator
  USE model_def
  USE integrator_def
  IMPLICIT NONE

  PRIVATE

  !> Class for the Heun (RK2) integrator object.
  TYPE, EXTENDS(Integrator), PUBLIC :: RK2Integrator
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y1 !< Buffer to hold the intermediate position (Heun algorithm)
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_f0 !< Buffer to hold tendencies at the initial position
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_f1 !< Buffer to hold tendencies at the intermediate position
  CONTAINS
    PROCEDURE :: init
    PROCEDURE :: step
    PROCEDURE :: clean
  END TYPE RK2Integrator

CONTAINS
  
  !> Routine to initialise the integration buffers.
  !> @param[in,out] integr Integrator object to initialize.
  !> @param[in] imodel Model object to initialize the integrator with.
  SUBROUTINE init(integr, imodel)
    CLASS(RK2Integrator), INTENT(INOUT) :: integr
    CLASS(Model), INTENT(IN), TARGET :: imodel
    INTEGER :: AllocStat
    IF (.NOT.imodel%initialized) THEN
      PRINT*, 'Model not yet initialized, impossible to associate an integrator to an empty model!'
      RETURN
    END IF

    integr%pmodel => imodel
    integr%dt => imodel%model_configuration%integration%dt
    integr%ndim => imodel%model_configuration%modes%ndim

    ALLOCATE(integr%buf_y1(0:integr%ndim) ,STAT=AllocStat)
    IF (AllocStat /= 0) THEN
      PRINT*, "*** rk2integrator%init: Problem with allocation! ***"
      STOP "Exiting ..."
    END IF
    ALLOCATE(integr%buf_f0(0:integr%ndim) ,STAT=AllocStat)
    IF (AllocStat /= 0) THEN
      PRINT*, "*** rk2integrator%init: Problem with allocation! ***"
      STOP "Exiting ..."
    END IF
    ALLOCATE(integr%buf_f1(0:integr%ndim) ,STAT=AllocStat)
    IF (AllocStat /= 0) THEN
      PRINT*, "*** rk2integrator%init: Problem with allocation! ***"
      STOP "Exiting ..."
    END IF

  END SUBROUTINE init

  !> Routine to perform an integration step (Heun algorithm). The incremented time is returned.
  !> @param[in,out] integr Integrator object to perform the step with.
  !> @param[in] y Initial point.
  !> @param[in] t Actual integration time
  !> @param[out] res Final point after the step.
  SUBROUTINE step(integr, y,t,res)
    CLASS(RK2Integrator), INTENT(INOUT) :: integr
    REAL(KIND=8), DIMENSION(0:integr%ndim), INTENT(IN) :: y
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), DIMENSION(0:integr%ndim), INTENT(OUT) :: res
    
    CALL integr%pmodel%tendencies(t,y,integr%buf_f0)
    integr%buf_y1 = y+integr%dt*integr%buf_f0
    CALL integr%pmodel%tendencies(t+integr%dt,integr%buf_y1,integr%buf_f1)
    res=y+0.5*(integr%buf_f0+integr%buf_f1)*integr%dt
    t=t+integr%dt
  END SUBROUTINE step

  !> Routine to clean the integrator
  !> @param[in,out] integr Integrator object to clean.
  SUBROUTINE clean(integr)
      CLASS(RK2Integrator), INTENT(INOUT) :: integr

      IF (allocated(integr%buf_y1)) DEALLOCATE(integr%buf_y1)
      IF (allocated(integr%buf_f1)) DEALLOCATE(integr%buf_f1)
      IF (allocated(integr%buf_f0)) DEALLOCATE(integr%buf_f0)

  END SUBROUTINE clean

END MODULE rk2_integrator
