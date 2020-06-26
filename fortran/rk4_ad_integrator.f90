
! rk4_ad_integrator.f90
!
!> Adjoint (AD) model versions of MAOOAM.
!> Fourth-order Runge-Kutta (RK4) integrators module.
!
!> @copyright                                                               
!> 2020 Lesley De Cruz, Jonathan Demaeyer & Sebastian Schubert.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!                                                                           
!---------------------------------------------------------------------------



MODULE rk4_ad_integrator

  USE model_def
  USE rk4_integrator
  IMPLICIT NONE

  PRIVATE

  !> Class for the fourth-order Runge-Kutta (RK4) AD integrator object.
  TYPE, EXTENDS(RK4Integrator), PUBLIC :: RK4AdIntegrator
  CONTAINS
    PROCEDURE :: ad_step
  END TYPE RK4AdIntegrator

CONTAINS

  !-----------------------------------------------------!
  !                                                     !
  ! Adjoint model integrator step                       !
  !                                                     !
  !-----------------------------------------------------!

  !> Routine to perform an integration step (RK4 algorithm) of the adjoint model. The incremented time is returned.
  !> @param[in,out] integr Integrator object to perform the step with.
  !> @param[in] y Initial point.
  !> @param[in] ystar Evaluating the adjoint model at the point \f$\boldsymbol{y}^\ast\f$.
  !> @param[in] t Actual integration time
  !> @param[out] res Final point after the step.
  SUBROUTINE ad_step(integr,y,ystar,t,res)
    CLASS(RK4AdIntegrator), INTENT(INOUT) :: integr
    REAL(KIND=8), DIMENSION(0:integr%ndim), INTENT(IN) :: y,ystar
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), DIMENSION(0:integr%ndim), INTENT(OUT) :: res

    CALL integr%pmodel%ad_tendencies(t,ystar,y,integr%buf_kA)
    integr%buf_y1 = y+0.5*integr%dt*integr%buf_kA
    CALL integr%pmodel%ad_tendencies(t+0.5*integr%dt,ystar,integr%buf_y1,integr%buf_kB)
    integr%buf_y1 = y+0.5*integr%dt*integr%buf_kB
    integr%buf_kA = integr%buf_kA+2*integr%buf_kB
    CALL integr%pmodel%ad_tendencies(t+0.5*integr%dt,ystar,integr%buf_y1,integr%buf_kB)
    integr%buf_y1 = y+0.5*integr%dt*integr%buf_kB
    integr%buf_kA = integr%buf_kA+2*integr%buf_kB
    CALL integr%pmodel%ad_tendencies(t+integr%dt,ystar,integr%buf_y1,integr%buf_kB)
    integr%buf_kA = integr%buf_kA+integr%buf_kB
    res=y+integr%buf_kA*integr%dt/6
    t=t+integr%dt
  END SUBROUTINE ad_step

END MODULE rk4_ad_integrator
