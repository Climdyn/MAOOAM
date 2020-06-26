
! rk2_tl_integrator.f90
!
!> Tangent Linear (TL) model versions of MAOOAM.
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



MODULE rk2_tl_integrator

  USE model_def
  USE rk2_integrator
  IMPLICIT NONE

  PRIVATE

  !> Class for the Heun (RK2) TL integrator object.
  TYPE, EXTENDS(RK2Integrator), PUBLIC :: RK2TlIntegrator
  CONTAINS
    PROCEDURE :: tl_step
  END TYPE RK2TlIntegrator

CONTAINS

  !-----------------------------------------------------!
  !                                                     !
  ! Tangent linear model integrator step                !
  !                                                     !
  !-----------------------------------------------------!

  !> Routine to perform an integration step (Heun algorithm) of the tangent linear model. The incremented time is returned.
  !> @param[in,out] integr Integrator object to perform the step with.
  !> @param[in] y Initial point.
  !> @param[in] ystar Evaluating the adjoint model at the point \f$\boldsymbol{y}^\ast\f$.
  !> @param[in] t Actual integration time
  !> @param[out] res Final point after the step.
  SUBROUTINE tl_step(integr,y,ystar,t,res)
    CLASS(RK2TlIntegrator), INTENT(INOUT) :: integr
    REAL(KIND=8), DIMENSION(0:integr%ndim), INTENT(IN) :: y,ystar
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), DIMENSION(0:integr%ndim), INTENT(OUT) :: res

    CALL integr%pmodel%tl_tendencies(t,ystar,y,integr%buf_f0)
    integr%buf_y1 = y+integr%dt*integr%buf_f0
    CALL integr%pmodel%tl_tendencies(t+integr%dt,ystar,integr%buf_y1,integr%buf_f1)
    res=y+0.5*(integr%buf_f0+integr%buf_f1)*integr%dt
    t=t+integr%dt
  END SUBROUTINE tl_step
  
END MODULE rk2_tl_integrator
