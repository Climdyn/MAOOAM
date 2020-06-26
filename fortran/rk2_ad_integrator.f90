
! rk2_ad_integrator.f90
!
!> Adjoint (AD) model versions of MAOOAM.
!> Second-order Runge-Kutta (RK2) integrators module.
!
!> @copyright                                                               
!> 2020 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!>  @remark                                                                 
!>  This module actually contains the Heun algorithm routines.
!                                                                           
!---------------------------------------------------------------------------



MODULE rk2_ad_integrator

  USE model_def
  USE rk2_integrator
  IMPLICIT NONE

  PRIVATE

  !> Class for the Heun (RK2) AD integrator object.
  TYPE, EXTENDS(RK2Integrator), PUBLIC :: RK2AdIntegrator
  CONTAINS
    PROCEDURE :: ad_step
  END TYPE RK2AdIntegrator

CONTAINS

  !-----------------------------------------------------!
  !                                                     !
  ! Adjoint model integrator step                       !
  !                                                     !
  !-----------------------------------------------------!

  !> Routine to perform an integration step (Heun algorithm) of the adjoint model. The incremented time is returned.
  !> @param[in,out] integr Integrator object to perform the step with.
  !> @param[in] y Initial point.
  !> @param[in] ystar Evaluating the adjoint model at the point \f$\boldsymbol{y}^\ast\f$.
  !> @param[in] t Actual integration time
  !> @param[out] res Final point after the step.
  SUBROUTINE ad_step(integr,y,ystar,t,res)
    CLASS(RK2AdIntegrator), INTENT(INOUT) :: integr
    REAL(KIND=8), DIMENSION(0:integr%ndim), INTENT(IN) :: y,ystar
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), DIMENSION(0:integr%ndim), INTENT(OUT) :: res

    CALL integr%pmodel%ad_tendencies(t,ystar,y,integr%buf_f0)
    integr%buf_y1 = y+integr%dt*integr%buf_f0
    CALL integr%pmodel%ad_tendencies(t+integr%dt,ystar,integr%buf_y1,integr%buf_f1)
    res=y+0.5*(integr%buf_f0+integr%buf_f1)*integr%dt
    t=t+integr%dt
  END SUBROUTINE ad_step
  
END MODULE rk2_ad_integrator
