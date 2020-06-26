
! tensor_def.f90
!
!>  Tensor utility module. Contains class to represent sparse tensors.
!
!> @copyright                                                               
!> 2015-2020 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!


MODULE tensor_def
  USE util, only: str
  IMPLICIT NONE

  PRIVATE

  ! Epsilon to test equality with 0
  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16

  !> Coordinate list element type. Elementary elements of the sparse tensors.
  TYPE :: CoolistElem
    INTEGER :: j !< Index \f$j\f$ of the element
    INTEGER :: k !< Index \f$k\f$ of the element
    REAL(KIND=8) :: v !< Value of the element
  END TYPE CoolistElem

  !> Coordinate list. Type used to represent the sparse tensor.
  TYPE :: CooList
    TYPE(CoolistElem), DIMENSION(:), ALLOCATABLE :: elems !< Lists of elements tensor_def::coolist_elem
    INTEGER :: nelems = 0                                 !< Number of elements in the list.
  END TYPE CooList
  
  !> General class to represent a sparse tensor
  TYPE, PUBLIC :: Tensor
    TYPE(CooList), DIMENSION(:), ALLOCATABLE :: t         !< Sparse representation of the tensor as a tensor_def::coolist .
  CONTAINS
    PROCEDURE :: init
    PROCEDURE :: sparse_mul3
    PROCEDURE :: copy
    PROCEDURE, PASS(dst) :: from_mat
    PROCEDURE :: simplify
    PROCEDURE :: jsparse_mul
    PROCEDURE :: jsparse_mul_mat
    PROCEDURE :: add_elem
    PROCEDURE, PASS(dst) :: add_from_tensor
    PROCEDURE :: print_tensor
    PROCEDURE :: load_from_file => load_tensor_from_file
    PROCEDURE :: write_to_file => write_tensor_to_file
    PROCEDURE :: clean
    PROCEDURE :: allocated => test_alloc
    PROCEDURE :: empty
    PROCEDURE :: ndim => tensor_size
!    PROCEDURE :: add_check
  END TYPE Tensor

CONTAINS

  !> Function to test if the tensor is allocated.
  !> @param mtensor The tensor to test.
  !> @return A boolean indicating if the tensor is allocated.
  FUNCTION test_alloc(mtensor)
    CLASS(Tensor) :: mtensor
    LOGICAL :: test_alloc

    test_alloc = allocated(mtensor%t)

  END FUNCTION test_alloc

  !> Function to test if the tensor is empty.
  !> @param mtensor The tensor to test.
  !> @return A boolean indicating if the tensor is empty.
  FUNCTION empty(mtensor)
    CLASS(Tensor) :: mtensor
    LOGICAL :: empty
    INTEGER :: i

    empty = .TRUE.

    IF (.NOT. mtensor%allocated()) RETURN

    DO i=1,mtensor%ndim()
      IF (mtensor%t(i)%nelems/=0) empty = .FALSE.
    END DO

  END FUNCTION empty

  !> Routine to clean (deallocate) a tensor.
  !> @param[in,out] mtensor The tensor to clean.
  SUBROUTINE clean(mtensor)
    CLASS(Tensor), INTENT(INOUT) :: mtensor

    IF (mtensor%allocated()) DEALLOCATE(mtensor%t)

  END SUBROUTINE clean

  !> Routine to initialize a tensor.
  !> @param[in,out] mtensor The tensor to clean.
  !> @param[in] ndim The first dimension of the tensor.
  SUBROUTINE init(mtensor, ndim)
    CLASS(Tensor), INTENT(INOUT) :: mtensor
    INTEGER, INTENT(IN) :: ndim
    INTEGER :: AllocStat

    CALL mtensor%clean
    ALLOCATE(mtensor%t(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) THEN
        PRINT*, "*** tensor%init: Problem with allocation! ***"
        STOP "Exiting ..."
    END IF

  END SUBROUTINE init

  ! Function to return the size of the tensor.
  !> @param mtensor The tensor to return the size of.
  !> @result The size of the tensor
  FUNCTION tensor_size(mtensor) RESULT(ndim)
    CLASS(Tensor) :: mtensor
    INTEGER :: ndim
    ndim = size(mtensor%t)
  END FUNCTION tensor_size


  !> Routine to copy a tensor into another one.
  !> @param[in] src Source tensor.
  !> @param[out] dst Destination tensor.
  !> @warning The destination tensor will be reinitialized, erasing all previous content! Use with care...
  SUBROUTINE copy(src,dst)
    CLASS(Tensor), INTENT(IN) :: src
    CLASS(Tensor), INTENT(OUT) :: dst
    INTEGER :: i,j,AllocStat
    
    CALL dst%init(src%ndim())
    DO i=1,src%ndim()
      ALLOCATE(dst%t(i)%elems(src%t(i)%nelems), STAT=AllocStat)
      IF (AllocStat /= 0) THEN
        PRINT*, "*** tensor%copy: Problem with allocation! ***"
        STOP "Exiting ..."
      END IF
      DO j=1,src%t(i)%nelems
        dst%t(i)%elems(j)%j=src%t(i)%elems(j)%j
        dst%t(i)%elems(j)%k=src%t(i)%elems(j)%k
        dst%t(i)%elems(j)%v=src%t(i)%elems(j)%v
      ENDDO
      dst%t(i)%nelems=src%t(i)%nelems
    ENDDO
  END SUBROUTINE copy

  !> Routine to convert a matrix to a tensor, using only the fist two indices of the rank-3 tensor.
  !> @param[in] src Source matrix
  !> @param[out] dst Destination tensor.
  !> @warning The destination tensor will be reinitialized, erasing all previous content! Use with care...
  SUBROUTINE from_mat(src,dst)
    CLASS(Tensor), INTENT(INOUT) :: dst
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: src
    INTEGER :: i,j,n,AllocStat
    INTEGER :: ndim
    INTEGER, DIMENSION(2) :: sh

    sh = shape(src)
    ndim = sh(1)
    CALL dst%init(ndim)

    DO i=1,ndim
      n=0
      DO j=1,ndim
        IF (abs(src(i,j))>real_eps) n=n+1
      ENDDO
      ALLOCATE(dst%t(i)%elems(n), STAT=AllocStat)
      IF (AllocStat /= 0) THEN
        PRINT*, "*** tensor%from_mat: Problem with allocation! ***"
        STOP "Exiting ..."
      END IF
      n=0
      DO j=1,ndim
        IF (abs(src(i,j))>real_eps) THEN
          n=n+1
          dst%t(i)%elems(n)%j=j
          dst%t(i)%elems(n)%k=0
          dst%t(i)%elems(n)%v=src(i,j)
        ENDIF
      ENDDO
      dst%t(i)%nelems=n
    ENDDO
  END SUBROUTINE from_mat
  
  !> Sparse multiplication of a tensor with two vectors:  \f${\displaystyle \sum_{j,k=0}^{ndim}} \mathcal{T}_{i,j,k} \, a_j \,b_k\f$.
  !> @param[in] mtensor A sparse tensor of which index 2 and 3 will be contracted.
  !> @param[in] arr_j The vector \f$\boldsymbol{a}\f$ to be contracted with index 2 of the tensor.
  !> @param[in] arr_k The vector \f$\boldsymbol{b}\f$ to be contracted with index 3 of the tensor.
  !> @param[out] res Vector to store the result of the contraction.
  !> @remark Note that it is NOT safe to pass `arr_j` or `arr_k` as a result buffer,
  !> as this operation does multiple passes. However, passsing the same vector as `arr_j` and `arr_k` is safe.
  SUBROUTINE sparse_mul3(mtensor, arr_j, arr_k, res)
    CLASS(Tensor), INTENT(IN) :: mtensor
    REAL(KIND=8), DIMENSION(0:size(mtensor%t)), INTENT(IN)  :: arr_j, arr_k
    REAL(KIND=8), DIMENSION(0:size(mtensor%t)), INTENT(OUT) :: res
    INTEGER :: i,j,k,n
    res=0.D0
    DO i=1,mtensor%ndim()
      DO n=1,mtensor%t(i)%nelems
        j=mtensor%t(i)%elems(n)%j
        k=mtensor%t(i)%elems(n)%k
        res(i) = res(i) + mtensor%t(i)%elems(n)%v * arr_j(j)*arr_k(k)
      END DO
    END DO
  END SUBROUTINE sparse_mul3

  !> Routine to simplify a coolist (sparse tensor). For each index \f$i\f$, it upper triangularize the matrix
  !> \f[\mathcal{T}_{i,j,k} \qquad 0 \leq j,k \leq ndim.\f]
  !> @param[in,out] mtensor A sparse tensor which will be simplified.
  SUBROUTINE simplify(mtensor)
    CLASS(Tensor), INTENT(INOUT) :: mtensor
    INTEGER :: i,j,k
    INTEGER :: li,lii,liii,n
    DO i= 1,mtensor%ndim()
      n=mtensor%t(i)%nelems
      DO li=n,2,-1
        j=mtensor%t(i)%elems(li)%j
        k=mtensor%t(i)%elems(li)%k
        DO lii=li-1,1,-1
          IF (((j==mtensor%t(i)%elems(lii)%j).AND.(k==mtensor%t(i)&
            &%elems(lii)%k)).OR.((j==mtensor%t(i)%elems(lii)%k).AND.(k==mtensor%t(i)%elems(lii)%j))) THEN
            ! Found another entry with the same i,j,k: merge both into
            ! the one listed first (of those two).
            mtensor%t(i)%elems(lii)%v=mtensor%t(i)%elems(lii)%v+mtensor%t(i)%elems(li)%v
            IF (j>k) THEN
              mtensor%t(i)%elems(lii)%j=mtensor%t(i)%elems(li)%k
              mtensor%t(i)%elems(lii)%k=mtensor%t(i)%elems(li)%j
            ENDIF

            ! Shift the rest of the items one place down.
            DO liii=li+1,n
              mtensor%t(i)%elems(liii-1)%j=mtensor%t(i)%elems(liii)%j
              mtensor%t(i)%elems(liii-1)%k=mtensor%t(i)%elems(liii)%k
              mtensor%t(i)%elems(liii-1)%v=mtensor%t(i)%elems(liii)%v
            END DO
            mtensor%t(i)%nelems=mtensor%t(i)%nelems-1
            ! Here we should stop because the li no longer points to the
            ! original i,j,k element
            EXIT
          ENDIF
        ENDDO
      ENDDO
      n=mtensor%t(i)%nelems
      li=1
      DO WHILE (li<=mtensor%t(i)%nelems)
        ! Clear new "almost" zero entries and shift rest of the items one place down.
        ! Make sure not to skip any entries while shifting!
        DO WHILE (abs(mtensor%t(i)%elems(li)%v) < real_eps)
          DO liii=li+1,n
            mtensor%t(i)%elems(liii-1)%j=mtensor%t(i)%elems(liii)%j
            mtensor%t(i)%elems(liii-1)%k=mtensor%t(i)%elems(liii)%k
            mtensor%t(i)%elems(liii-1)%v=mtensor%t(i)%elems(liii)%v
          ENDDO
          mtensor%t(i)%nelems=mtensor%t(i)%nelems-1
          if (li > mtensor%t(i)%nelems) THEN
            EXIT
          ENDIF
        ENDDO
        li=li+1
      ENDDO

      n=mtensor%t(i)%nelems
      DO li=1,n
        ! Upper triangularize
        j=mtensor%t(i)%elems(li)%j
        k=mtensor%t(i)%elems(li)%k
        IF (j>k) THEN
          mtensor%t(i)%elems(li)%j=k
          mtensor%t(i)%elems(li)%k=j
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE simplify

  !> Sparse multiplication of two tensors to determine the Jacobian:
  !> \f[J_{i,j} = {\displaystyle \sum_{k=0}^{ndim}} \left( \mathcal{T}_{i,j,k} + \mathcal{T}_{i,k,j} \right) \, a_k.\f]
  !> It's implemented slightly differently: for every \f$\mathcal{T}_{i,j,k}\f$, we add to \f$J_{i,j}\f$ as follows:
  !> \f[J_{i,j} = J_{i,j} + \mathcal{T}_{i,j,k} \, a_k \\ J_{i,k} = J_{i,k} + \mathcal{T}_{i,j,k} \, a_j\f]
  !> This version return a sparse tensor.
  !> @param[in] mtensor A sparse tensor of which index 2 or 3 will be contracted.
  !> @param[in] arr_j The vector \f$\boldsymbol{a}\f$ to be contracted with.
  !> @param[out] jtensor A sparse tensor to store the result of the contraction
  !> @warning The output jtensor will be reinitialized, erasing all previous content! Use with care...
  SUBROUTINE jsparse_mul(mtensor, arr_j, jtensor)
    CLASS(Tensor), INTENT(IN) :: mtensor
    TYPE(Tensor), INTENT(INOUT):: jtensor
    REAL(KIND=8), DIMENSION(0:size(mtensor%t)), INTENT(IN)  :: arr_j
    REAL(KIND=8) :: v
    INTEGER :: i,j,k,n,nj,AllocStat
    CALL jtensor%init(mtensor%ndim())
    DO i=1,mtensor%ndim()
       nj=2*jtensor%t(i)%nelems
       ALLOCATE(jtensor%t(i)%elems(nj), STAT=AllocStat)
       IF (AllocStat /= 0) THEN
         PRINT*, "*** tensor%jsparse_mul: Problem with allocation! ***"
         STOP "Exiting ..."
       END IF
       nj=0
       DO n=1,mtensor%t(i)%nelems
          j=mtensor%t(i)%elems(n)%j
          k=mtensor%t(i)%elems(n)%k
          v=mtensor%t(i)%elems(n)%v
          IF (j /=0) THEN
             nj=nj+1
             jtensor%t(i)%elems(nj)%j=j
             jtensor%t(i)%elems(nj)%k=0
             jtensor%t(i)%elems(nj)%v=v*arr_j(k)
          END IF

          IF (k /=0) THEN
             nj=nj+1
             jtensor%t(i)%elems(nj)%j=k
             jtensor%t(i)%elems(nj)%k=0
             jtensor%t(i)%elems(nj)%v=v*arr_j(j)
          END IF
       END DO
       jtensor%t(i)%nelems=nj
    END DO
  END SUBROUTINE jsparse_mul

  !> Sparse multiplication of two tensors to determine the Jacobian:
  !> \f[J_{i,j} = {\displaystyle \sum_{k=0}^{ndim}} \left( \mathcal{T}_{i,j,k} + \mathcal{T}_{i,k,j} \right) \, a_k.\f]
  !> It's implemented slightly differently: for every \f$\mathcal{T}_{i,j,k}\f$, we add to \f$J_{i,j}\f$ as follows:
  !> \f[J_{i,j} = J_{i,j} + \mathcal{T}_{i,j,k} \, a_k \\ J_{i,k} = J_{i,k} + \mathcal{T}_{i,j,k} \, a_j\f]
  !> This version return a matrix.
  !> @param[in] mtensor A sparse tensor of which index 2 or 3 will be contracted.
  !> @param[in] arr_j The vector \f$\boldsymbol{a}\f$ to be contracted with.
  !> @param[out] jmatrix A matrix to store the result of the contraction.
  SUBROUTINE jsparse_mul_mat(mtensor, arr_j, jmatrix)
    CLASS(Tensor), INTENT(IN) :: mtensor
    REAL(KIND=8), DIMENSION(size(mtensor%t),size(mtensor%t)), INTENT(OUT):: jmatrix
    REAL(KIND=8), DIMENSION(0:size(mtensor%t)), INTENT(IN)  :: arr_j
    REAL(KIND=8) :: v
    INTEGER :: i,j,k,n
    jmatrix=0.D0
    DO i=1,mtensor%ndim()
       DO n=1,mtensor%t(i)%nelems
          j=mtensor%t(i)%elems(n)%j
          k=mtensor%t(i)%elems(n)%k
          v=mtensor%t(i)%elems(n)%v
          IF (j /=0) jmatrix(i,j)=jmatrix(i,j)+v*arr_j(k)
          IF (k /=0) jmatrix(i,k)=jmatrix(i,k)+v*arr_j(j)
       END DO
    END DO
  END SUBROUTINE jsparse_mul_mat

  !> Sparse multiplication of a 2d sparse tensor with a vector:  \f${\displaystyle \sum_{j=0}^{ndim}} \mathcal{T}_{i,j,k} \, a_j \f$.
  !> @param[in] mtensor A sparse tensor of which index 2 will be contracted.
  !> @param[in] arr_j The vector \f$\boldsymbol{a}\f$ to be contracted with.
  !> @param[out] res vector (buffer) to store the result of the contraction
  !> @remark Note that it is NOT safe to pass `arr_j` as a result buffer,
  !> as this operation does multiple passes.
  SUBROUTINE sparse_mul2(mtensor, arr_j, res)
    CLASS(Tensor), INTENT(IN) :: mtensor
    REAL(KIND=8), DIMENSION(0:size(mtensor%t)), INTENT(IN)  :: arr_j
    REAL(KIND=8), DIMENSION(0:size(mtensor%t)), INTENT(OUT) :: res
    INTEGER :: i,j,n
    res=0.D0
    DO i=1,mtensor%ndim()
      DO n=1,mtensor%t(i)%nelems
         j=mtensor%t(i)%elems(n)%j
         res(i) = res(i) + mtensor%t(i)%elems(n)%v * arr_j(j)
      END DO
    END DO
  END SUBROUTINE sparse_mul2

  !> Subroutine to add element to a coolist.
  !> @param[in,out] mtensor A tensor to add the element to.
  !> @param[in] i tensor \f$i\f$ index
  !> @param[in] j tensor \f$j\f$ index
  !> @param[in] k tensor \f$k\f$ index
  !> @param[in] v value to add
  SUBROUTINE add_elem(mtensor,i,j,k,v)
    CLASS(Tensor), INTENT(INOUT) :: mtensor
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN) :: v
    INTEGER :: n
    IF (abs(v) .ge. real_eps) THEN
       n=(mtensor%t(i)%nelems)+1
       mtensor%t(i)%elems(n)%j=j
       mtensor%t(i)%elems(n)%k=k
       mtensor%t(i)%elems(n)%v=v
       mtensor%t(i)%nelems=n
    END IF
  END SUBROUTINE add_elem

!  !> Subroutine to add element to a coolist and check for overflow.
!  !> Once the t buffer tensor is full, add it to the destination buffer.
!  !> @param t temporary buffer tensor for the destination tensor
!  !> @param i tensor \f$i\f$ index
!  !> @param j tensor \f$j\f$ index
!  !> @param k tensor \f$k\f$ index
!  !> @param v value to add
!  !> @param dst destination tensor
!  SUBROUTINE add_check(t,i,j,k,v,dst)
!    TYPE(CooList), DIMENSION(ndim), INTENT(INOUT) :: t
!    TYPE(CooList), DIMENSION(ndim), INTENT(INOUT) :: dst
!    INTEGER, INTENT(IN) :: i,j,k
!    REAL(KIND=8), INTENT(IN) :: v
!    INTEGER :: n
!    CALL add_elem(t,i,j,k,v)
!    IF (t(i)%nelems==size(t(i)%elems)) THEN
!       CALL add_to_tensor(t,dst)
!       DO n=1,ndim
!           t(n)%nelems=0
!       ENDDO
!    ENDIF
!  END SUBROUTINE add_check

  !> Routine to add the entries of a rank-3 tensor to another one.
  !> @param[in] src Tensor to add
  !> @param[in,out] dst Destination tensor
  SUBROUTINE add_from_tensor(src,dst)
    CLASS(Tensor), INTENT(IN) :: src
    CLASS(Tensor), INTENT(INOUT) :: dst
    TYPE(CoolistElem), DIMENSION(:), ALLOCATABLE :: celems
    INTEGER :: i,j,n,AllocStat

    DO i=1,dst%ndim()
       IF (src%t(i)%nelems/=0) THEN
          IF (dst%t(i)%nelems==0) THEN
             IF (allocated(dst%t(i)%elems)) THEN
                DEALLOCATE(dst%t(i)%elems, STAT=AllocStat)
                IF (AllocStat /= 0) THEN
                  PRINT*, "*** tensor%add_from_tensor: Problem with allocation! ***"
                  STOP "Exiting ..."
                END IF
             END IF
             ALLOCATE(dst%t(i)%elems(src%t(i)%nelems), STAT=AllocStat)
             IF (AllocStat /= 0) THEN
                  PRINT*, "*** tensor%add_from_tensor: Problem with allocation! ***"
                  STOP "Exiting ..."
             END IF
             n=0
          ELSE
             n=dst%t(i)%nelems
             ALLOCATE(celems(n), STAT=AllocStat)
             DO j=1,n
                celems(j)%j=dst%t(i)%elems(j)%j
                celems(j)%k=dst%t(i)%elems(j)%k
                celems(j)%v=dst%t(i)%elems(j)%v
             ENDDO
             IF (allocated(dst%t(i)%elems)) DEALLOCATE(dst%t(i)%elems, STAT=AllocStat)
             ALLOCATE(dst%t(i)%elems(src%t(i)%nelems+n), STAT=AllocStat)
             IF (AllocStat /= 0) THEN
                  PRINT*, "*** tensor%add_from_tensor: Problem with allocation! ***"
                  STOP "Exiting ..."
             END IF
             DO j=1,n
                dst%t(i)%elems(j)%j=celems(j)%j
                dst%t(i)%elems(j)%k=celems(j)%k
                dst%t(i)%elems(j)%v=celems(j)%v
             ENDDO
             IF (allocated(celems)) DEALLOCATE(celems, STAT=AllocStat)
          ENDIF
          DO j=1,src%t(i)%nelems
             dst%t(i)%elems(n+j)%j=src%t(i)%elems(j)%j
             dst%t(i)%elems(n+j)%k=src%t(i)%elems(j)%k
             dst%t(i)%elems(n+j)%v=src%t(i)%elems(j)%v
          ENDDO
          dst%t(i)%nelems=src%t(i)%nelems+n
       ENDIF
    ENDDO

  END SUBROUTINE add_from_tensor

  !> Routine to print a rank-3 tensor.
  !> @param[in] mtensor Tensor to print.
  !> @param[in] s String to put before tensor entries. Default to "t".
  SUBROUTINE print_tensor(mtensor,s)
    CLASS(Tensor), INTENT(IN) :: mtensor
    CHARACTER(LEN=*), INTENT(IN), TARGET, OPTIONAL :: s

    CHARACTER, TARGET :: sr = "t"
    CHARACTER, POINTER :: r
    INTEGER :: i,n,j,k
    IF (present(s)) THEN
       r => s
    ELSE
       r => sr
    END IF
    DO i=1,mtensor%ndim()
       DO n=1,mtensor%t(i)%nelems
          j=mtensor%t(i)%elems(n)%j
          k=mtensor%t(i)%elems(n)%k
          IF( abs(mtensor%t(i)%elems(n)%v) .GE. real_eps) THEN
             write(*,"(A,ES12.5)") r//"["//TRIM(str(i))//"]["//TRIM(str(j)) &
                  &//"]["//TRIM(str(k))//"] = ",mtensor%t(i)%elems(n)%v
          END IF
       END DO
    END DO
  END SUBROUTINE print_tensor

  !> Write a rank-3 tensor coolist to a file.
  !> @param mtensor The tensor to write
  !> @param[in] s Filename
  SUBROUTINE write_tensor_to_file(mtensor, s)
    CLASS(Tensor), INTENT(IN) :: mtensor
    CHARACTER (LEN=*), INTENT(IN) :: s
    INTEGER :: i,j,k,n
    OPEN(30,file=s)
    WRITE(30,*) mtensor%ndim()
    DO i=1,mtensor%ndim()
       WRITE(30,*) i,mtensor%t(i)%nelems
       DO n=1,mtensor%t(i)%nelems
          j=mtensor%t(i)%elems(n)%j
          k=mtensor%t(i)%elems(n)%k
          WRITE(30,*) i,j,k,mtensor%t(i)%elems(n)%v
       END DO
    END DO
    CLOSE(30)
  END SUBROUTINE write_tensor_to_file

  !> Load a rank-3 tensor coolist from a file definition
  !> @param[in,out] mtensor The tensor to load to.
  !> @param[in] s Filename of the tensor definition file.
  !> @remark The destination tensor have to be an empty tensor, i.e. with unallocated list of elements and nelems set to 0.
  SUBROUTINE load_tensor_from_file(mtensor,s)
    CLASS(Tensor), INTENT(INOUT) :: mtensor
    CHARACTER (LEN=*), INTENT(IN) :: s
    INTEGER :: i,ir,j,k,n,AllocStat, ndim
    REAL(KIND=8) :: v
    OPEN(30,file=s,status='old')
    READ(30, *) ndim
    ALLOCATE(mtensor%t(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) THEN
      PRINT*, "*** tensor%load_tensor_from_file: Problem with allocation! ***"
      STOP "Exiting ..."
    END IF
    DO i=1,ndim
       READ(30,*) ir,n
       IF (n /= 0) THEN
          ALLOCATE(mtensor%t(i)%elems(n), STAT=AllocStat)
          IF (AllocStat /= 0) THEN
            PRINT*, "*** tensor%load_tensor_from_file: Problem with allocation! ***"
            STOP "Exiting ..."
          END IF
          mtensor%t(i)%nelems=n
       ENDIF
       DO n=1,mtensor%t(i)%nelems
          READ(30,*) ir,j,k,v
          mtensor%t(i)%elems(n)%j=j
          mtensor%t(i)%elems(n)%k=k
          mtensor%t(i)%elems(n)%v=v
       ENDDO
    END DO
    CLOSE(30)
  END SUBROUTINE load_tensor_from_file


END MODULE tensor_def

