  module SparseSymmetricMatrix
  implicit none

  type :: SparseSymmetricMatrixType
    integer :: n
    integer, allocatable :: row(:)
    integer, allocatable :: col(:)
    real, allocatable :: val(:)
  end type SparseSymmetricMatrixType

contains

  subroutine StoreSymmetricMatrix(matrix, n, row, col, val)
    type(SparseSymmetricMatrixType), intent(out) :: matrix
    integer, intent(in) :: n
    integer, intent(in) :: row(:)
    integer, intent(in) :: col(:)
    real, intent(in) :: val(:)
    integer :: i, nnz

    ! Count non-zero elements in the upper triangular part
    nnz = 0
    do i = 1, size(row)
      if (row(i) <= col(i)) then
        nnz = nnz + 1
      end if
    end do

    matrix%n = n
    allocate(matrix%row(nnz))
    allocate(matrix%col(nnz))
    allocate(matrix%val(nnz))

    nnz = 0
    do i = 1, size(row)
      if (row(i) <= col(i)) then
        nnz = nnz + 1
        matrix%row(nnz) = row(i)
        matrix%col(nnz) = col(i)
        matrix%val(nnz) = val(i)
      end if
    end do
  end subroutine StoreSymmetricMatrix

  subroutine ReadSymmetricMatrix(matrix, n, row, col, val)
    type(SparseSymmetricMatrixType), intent(in) :: matrix
    integer, intent(out) :: n
    integer, allocatable, intent(out) :: row(:)
    integer, allocatable, intent(out) :: col(:)
    real, allocatable, intent(out) :: val(:)

    n = matrix%n
    allocate(row(size(matrix%row)))
    allocate(col(size(matrix%col)))
    allocate(val(size(matrix%val)))

    row = matrix%row
    col = matrix%col
    val = matrix%val
  end subroutine ReadSymmetricMatrix

  function GetElement(matrix, i, j) result(value)
    type(SparseSymmetricMatrixType), intent(in) :: matrix
    integer, intent(in) :: i, j
    real :: value
    integer :: k

    value = 0.0
    if (i > j) then
      ! Swap i and j to ensure we are always looking in the upper triangular part
      call swap(i, j)
    end if

    do k = 1, size(matrix%row)
      if (matrix%row(k) == i .and. matrix%col(k) == j) then
        value = matrix%val(k)
        return
      end if
    end do
  end function GetElement

  subroutine swap(a, b)
    integer, intent(inout) :: a, b
    integer :: temp

    temp = a
    a = b
    b = temp
  end subroutine swap

end module SparseSymmetricMatrix
