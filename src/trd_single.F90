subroutine trd_single_solver_kernel(neqs, nsize, a, b, c, rhs, x)
    integer, intent(in) :: neqs !垂向三对角方程个数
    integer, intent(in) :: nsize !垂向三对角方程的行数
    md_type, intent(in) :: a(neqs, nsize), b(neqs, nsize), c(neqs, nsize) !代表对角元和两个非对角元
    md_type, intent(in) :: rhs(neqs, nsize) !求解右端向量
    md_type, intent(out) :: x(neqs, nsize) !求解结果

    md_type :: cp(neqs, nsize)!, dp(neqs, nsize)

    integer :: i, j

    !x(:, :) = rhs(:, :)

    cp(:, 1) = c(:, 1) / b(:, 1)
    x(:, 1) = rhs(:, 1) / b(:, 1)
    do i = 2, nsize
        cp(:, i) = c(:, i) / (b(:, i) - cp(:, i-1) * a(:, i))
        x(:, i) = (rhs(:, i) - x(:, i-1) * a(:, i)) / (b(:, i) - cp(:, i-1) * a(:, i))
    end do

    do i = nsize - 1, 1, -1
        x(:, i) = x(:, i) - cp(:, i) * x(:, i+1)
    end do

end subroutine