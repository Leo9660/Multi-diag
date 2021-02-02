! struct sptrsv
subroutine struct_sptrsv(n_size, nnz, nrhs, a, x)
    integer, intent(in) :: n_size
    integer, intent(in) :: nnz
    integer, intent(in) :: nrhs
    md_type, intent(in) :: a(nnz, n_size)
    md_type, intent(inout) :: x(nrhs, n_size)

    integer :: i, j, k
    integer :: lnz

    lnz = (nnz - 1) / 2
    !md_type :: res(nrhs, n_size)

    ! forward  LY = RHS
    do i = 1, n_size
        do j = 1, lnz
            k = i - lnz - 1 + j
            if (k >= 1) x(:, i) = x(:, i) - a(j, i) * x(:, k)
        end do
    end do
    ! backward UX = Y
    do i = n_size, 1, -1
        do j = 1, lnz
            k = i + j
            if (k <= n_size) x(:, i) = x(:, i) - a(lnz + 1 + j, i) * x(:, k)
        end do
        x(:, i) = x(:, i) / a(lnz + 1, i)
    end do
end subroutine
