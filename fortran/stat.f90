
! stat.f90
!
!>  Statistics accumulators
!
!> @copyright                                                               
!> 2015-2020 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!



MODULE stat
  IMPLICIT NONE

  PRIVATE
  
  !> Statistics accumulator objects class
  TYPE, PUBLIC :: StatAccumulator
    INTEGER :: i=0 !< Number of stats accumulated
  
    ! Vectors holding the stats
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m       !< Vector storing the inline mean
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: mprev   !< Previous mean vector
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: v       !< Vector storing the inline variance
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: mtmp
  CONTAINS
    PROCEDURE :: init => init_stat
    PROCEDURE :: accumulate => acc
    PROCEDURE :: mean
    PROCEDURE :: var
    PROCEDURE :: iter
    PROCEDURE :: reset
    PROCEDURE :: clean
  END TYPE StatAccumulator

  CONTAINS

    !> Initialize the accumulators
    !> @param[in,out] istat Statistical accumulator to initialize
    !> @param[in] ndim Dimension of the state space to accumulate statistics for.
    SUBROUTINE init_stat(istat, ndim)
      CLASS(StatAccumulator), INTENT(INOUT) :: istat
      INTEGER, INTENT(in) :: ndim
      INTEGER :: AllocStat
      
      ALLOCATE(istat%m(ndim), istat%mprev(ndim), STAT=AllocStat)
      IF (AllocStat /= 0) THEN
        PRINT*, "*** init_stat: Problem with allocation! ***"
        STOP "Exiting ..."
      END IF
      ALLOCATE(istat%v(ndim), istat%mtmp(ndim), STAT=AllocStat)
      IF (AllocStat /= 0) THEN
        PRINT*, "*** init_stat: Problem with allocation! ***"
        STOP "Exiting ..."
      END IF
      istat%m=0.D0
      istat%mprev=0.D0
      istat%v=0.D0
      istat%mtmp=0.D0
      
    END SUBROUTINE init_stat

    !> Accumulate one state
    !> @param[in,out] istat Statistical accumulator to initialize
    !> @param[in] x State to accumulate
    SUBROUTINE acc(istat, x)
      CLASS(StatAccumulator), INTENT(INOUT) :: istat
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x
      istat%i=istat%i+1
      istat%mprev=istat%m+(x-istat%m)/istat%i
      istat%mtmp=istat%mprev
      istat%mprev=istat%m
      istat%m=istat%mtmp
      istat%v=istat%v+(x-istat%mprev)*(x-istat%m)
    END SUBROUTINE acc

    !> Function returning the mean
    !> @param[in,out] istat Statistical accumulator to initialize
    !> @return The mean of the accumulated states
    FUNCTION mean(istat)
      CLASS(StatAccumulator), INTENT(IN) :: istat
      REAL(KIND=8), DIMENSION(size(istat%m)) :: mean
      mean=istat%m
    END FUNCTION mean

    !> Function returning the variance
    !> @param[in,out] istat Statistical accumulator to initialize
    !> @return The variance of the accumulated states
    FUNCTION var(istat)
      CLASS(StatAccumulator), INTENT(IN) :: istat
      REAL(KIND=8), DIMENSION(size(istat%m)) :: var
      var=istat%v/(istat%i-1)
    END FUNCTION var

    !> Function returning the number of data accumulated
    !> @param[in,out] istat Statistical accumulator to initialize
    !> @return The number of the accumulated states
    FUNCTION iter(istat)
      CLASS(StatAccumulator), INTENT(IN) :: istat
      INTEGER :: iter
      iter=istat%i
    END FUNCTION iter

    !> Routine resetting the accumulator
    !> @param[in,out] istat Statistical accumulator to initialize
    SUBROUTINE reset(istat)
      CLASS(StatAccumulator), INTENT(INOUT) :: istat
      istat%m=0.D0
      istat%mprev=0.D0
      istat%v=0.D0
      istat%i=0
    END SUBROUTINE reset

  !> Routine to clean the accumulator
  !> @param[in,out] istat Statistical accumulator to clean
  SUBROUTINE clean(istat)
      CLASS(StatAccumulator), INTENT(INOUT) :: istat

      IF (allocated(istat%m)) DEALLOCATE(istat%m)
      IF (allocated(istat%mprev)) DEALLOCATE(istat%mprev)
      IF (allocated(istat%v)) DEALLOCATE(istat%v)
      IF (allocated(istat%mtmp)) DEALLOCATE(istat%mtmp)

  END SUBROUTINE clean


END MODULE stat
