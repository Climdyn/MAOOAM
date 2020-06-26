
! test_inprod_analytic.f90
!
!> Small program to print the inner products.   
!     
!> @copyright                                                               
!> 2015-2020 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!
!>  @remark
!>  Print in the same order as test_inprod.lua
!
!---------------------------------------------------------------------------!


PROGRAM inprod_analytic_test

  USE params, only: ModelConfiguration
  USE inprod_analytic
  USE util, only: str

  IMPLICIT NONE
  INTEGER :: i,j,kk
  ! Epsilon to test equality with 0
  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16

  TYPE(ModelConfiguration) :: mc
  TYPE(InnerProducts) :: ips

  CALL mc%init
  CALL ips%init(mc)

  do i = 1, mc%modes%natm
    do j = 1, mc%modes%natm
      IF(abs(ips%atmos%a(i,j)) .GE. real_eps) write(*,"(A,ES12.5)") &
      & "a["//TRIM(str(i))//"]["//TRIM(str(j))// "] = ", ips%atmos%a(i,j)
      IF(abs(ips%atmos%c(i,j)) .GE. real_eps) write(*,"(A,ES12.5)") &
      & "c["//TRIM(str(i))//"]["//TRIM(str(j))// "] = ", ips%atmos%c(i,j)
      do kk = 1, mc%modes%natm
        IF( abs(ips%atmos%b(i,j,kk)) .GE. real_eps) write(*,"(A,ES12.5)") &
        & "b[" //TRIM(str(i))//"][" //TRIM(str(j))//"][" //TRIM(str(kk))//"] = ", ips%atmos%b(i,j,kk)
        IF( abs(ips%atmos%g(i,j,kk)) .GE. real_eps) write(*,"(A,ES12.5)") &
        & "g[" //TRIM(str(i))//"][" //TRIM(str(j))//"][" //TRIM(str(kk))//"] = ", ips%atmos%g(i,j,kk)
      end do
    end do
    do j = 1, mc%modes%noc
      IF(abs(ips%atmos%d(i,j)) .GE. real_eps) write(*,"(A,ES12.5)") &
      & "d[" //TRIM(str(i))//"][" //TRIM(str(j))//"] = ", ips%atmos%d(i,j)
      IF(abs(ips%atmos%s(i,j)) .GE. real_eps) write(*,"(A,ES12.5)") &
      & "s[" //TRIM(str(i))//"][" //TRIM(str(j))//"] = ", ips%atmos%s(i,j)
    end do
  end do
  do i = 1, mc%modes%noc
    do j = 1, mc%modes%noc
      IF(abs(ips%ocean%M(i,j)) .GE. real_eps) write(*,"(A,ES12.5)") &
      & "M[" //TRIM(str(i))//"][" //TRIM(str(j))//"] = ", ips%ocean%M(i,j)
      IF(abs(ips%ocean%N(i,j)) .GE. real_eps) write(*,"(A,ES12.5)") &
      & "N[" //TRIM(str(i))//"][" //TRIM(str(j))//"] = ", ips%ocean%N(i,j)
      do kk = 1, mc%modes%noc
        IF(abs(ips%ocean%O(i,j,kk)) .GE. real_eps) write(*,"(A,ES12.5)") &
        & "O[" //TRIM(str(i))//"][" //TRIM(str(j))//"][" //TRIM(str(kk))//"] = ", ips%ocean%O(i,j,kk)
        IF(abs(ips%ocean%C(i,j,kk)) .GE. real_eps) write(*,"(A,ES12.5)") &
        & "C[" //TRIM(str(i))//"][" //TRIM(str(j))//"][" //TRIM(str(kk))//"] = ", ips%ocean%C(i,j,kk)
      end do
    end do
    do j = 1, mc%modes%natm
      IF(abs(ips%ocean%K(i,j)) .GE. real_eps) write(*,"(A,ES12.5)") &
      & "K[" //TRIM(str(i))//"][" //TRIM(str(j))//"] = ", ips%ocean%K(i,j)
      IF(abs(ips%ocean%W(i,j)) .GE. real_eps) write(*,"(A,ES12.5)") &
      & "W[" //TRIM(str(i))//"][" //TRIM(str(j))//"] = ", ips%ocean%W(i,j)
    end do
  end do

END PROGRAM inprod_analytic_test

