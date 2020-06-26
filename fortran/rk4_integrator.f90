
! rk4_integrator.f90
!
!>  Module containing the fourth-order Runge-Kutta (RK4) integration routines.
!
!> @copyright                                                               
!> 2015-2020 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!                                                                           
!---------------------------------------------------------------------------

MODULE rk4_integrator
  USE model_def
  USE integrator_def
  IMPLICIT NONE

  PRIVATE

  !> Class for the fourth-order Runge-Kutta (RK4) integrator object.
  TYPE, EXTENDS(Integrator), PUBLIC :: RK4Integrator
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y1 !< Buffer to hold the intermediate position
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_kA !< Buffer to hold tendencies at the initial position
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_kB !< Buffer to hold tendencies at the intermediate position
  CONTAINS
    PROCEDURE :: init
    PROCEDURE :: step
    PROCEDURE :: clean
  END TYPE RK4Integrator

CONTAINS

  !> Routine to initialise the integration buffers.
  !> @param[in,out] integr Integrator object to initialize.
  !> @param[in] imodel Model object to initialize the integrator with.
  SUBROUTINE init(integr, imodel)
    CLASS(RK4Integrator), INTENT(INOUT) :: integr
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
      PRINT*, "*** rk4integrator%init: Problem with allocation! ***"
      STOP "Exiting ..."
    END IF
    ALLOCATE(integr%buf_kA(0:integr%ndim) ,STAT=AllocStat)
    IF (AllocStat /= 0) THEN
      PRINT*, "*** rk4integrator%init: Problem with allocation! ***"
      STOP "Exiting ..."
    END IF
    ALLOCATE(integr%buf_kB(0:integr%ndim) ,STAT=AllocStat)
    IF (AllocStat /= 0) THEN
      PRINT*, "*** rk4integrator%init: Problem with allocation! ***"
      STOP "Exiting ..."
    END IF

  END SUBROUTINE init
  
  !> Routine to perform an integration step (RK4 algorithm). The incremented time is returned.
  !> @param[in,out] integr Integrator object to perform the step with.
  !> @param[in] y Initial point.
  !> @param[in] t Actual integration time
  !> @param[out] res Final point after the step.
  SUBROUTINE step(integr, y,t,res)
    CLASS(RK4Integrator), INTENT(INOUT) :: integr
    REAL(KIND=8), DIMENSION(0:integr%ndim), INTENT(IN) :: y
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), DIMENSION(0:integr%ndim), INTENT(OUT) :: res

    CALL integr%pmodel%tendencies(t,y,integr%buf_kA)
    integr%buf_y1 = y + 0.5*integr%dt*integr%buf_kA

    CALL integr%pmodel%tendencies(t+0.5*integr%dt,integr%buf_y1,integr%buf_kB)
    integr%buf_y1 = y + 0.5*integr%dt*integr%buf_kB
    integr%buf_kA = integr%buf_kA + 2*integr%buf_kB
    
    CALL integr%pmodel%tendencies(t+0.5*integr%dt,integr%buf_y1,integr%buf_kB)
    integr%buf_y1 = y + integr%dt*integr%buf_kB
    integr%buf_kA = integr%buf_kA + 2*integr%buf_kB
    
    CALL integr%pmodel%tendencies(t+integr%dt,integr%buf_y1,integr%buf_kB)
    integr%buf_kA = integr%buf_kA + integr%buf_kB
    
    t=t+integr%dt
    res=y+integr%buf_kA*integr%dt/6
  END SUBROUTINE step

  !> Routine to clean the integrator
  !> @param[in,out] integr Integrator object to clean.
  SUBROUTINE clean(integr)
      CLASS(RK4Integrator), INTENT(INOUT) :: integr

      IF (allocated(integr%buf_y1)) DEALLOCATE(integr%buf_y1)
      IF (allocated(integr%buf_kA)) DEALLOCATE(integr%buf_kA)
      IF (allocated(integr%buf_kB)) DEALLOCATE(integr%buf_kB)

  END SUBROUTINE clean

END MODULE rk4_integrator
