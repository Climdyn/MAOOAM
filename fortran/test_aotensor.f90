
! test_aotensor.f90
!
!> Small program to print the tendencies tensor.
!     
!> @copyright                                                               
!> 2015-2020 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!


PROGRAM test_aotensor

  USE params, only: ModelConfiguration
  USE inprod_analytic
  USE aotensor_def, only: AtmOcTensor
  USE util, only: str
  IMPLICIT NONE

  TYPE(ModelConfiguration) :: mc
  TYPE(InnerProducts) :: ips
  TYPE(AtmOcTensor) :: aot

  ! Epsilon to test equality with 0
  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16
  INTEGER :: i,j,k, n

  CALL mc%init
  CALL ips%init(mc)
  CALL aot%init(mc, ips)

  

  ! Program
  DO i=1,mc%modes%ndim
    DO n=1,aot%tensor%t(i)%nelems
      j=aot%tensor%t(i)%elems(n)%j
      k=aot%tensor%t(i)%elems(n)%k
      IF( abs(aot%tensor%t(i)%elems(n)%v) .GE. real_eps) THEN
        write(*,"(A,ES12.5)") "aotensor["//TRIM(str(i))//"]["//TRIM(str(j)) &
        &//"]["//TRIM(str(k))//"] = ",aot%tensor%t(i)%elems(n)%v
      END IF
    END DO
  END DO

END PROGRAM test_aotensor
