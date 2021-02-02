! circular multi-diagonal solver

subroutine md_cyclic_fact_kernel(neqs, ma)
    integer, intent(in) :: neqs !纬圈三对角方程个数
    type(md_matrix), intent(inout) :: ma(neqs) !预处理结构体，每个方程对应一个

    integer :: n, i, j, k
    md_type :: lij
    integer :: n_size, lnz, nnz

    integer :: buf_count, buf_size, pos_v, pos_w
    integer :: status(mpi_status_size), ierr
    md_type, allocatable, dimension(:) :: send_buf_v, recv_buf_v
    md_type, allocatable, dimension(:) :: send_buf_w, recv_buf_w

    ! LU factorization
    lnz = ma(1)%lnz
    nnz = ma(1)%nnz
    n_size = ma(1)%n_size

    do n = 2, neqs
        if (ma(n)%nnz .ne. nnz .or. ma(n)%n_size .ne. n_size) then
            print *, "Error: matrices of md_cyclic_fact must have the same size!"
            return
        end if
    end do

    buf_count = lnz * lnz * neqs
    buf_size = buf_count * md_size
    pos_v = 0; pos_w = 0

    allocate(send_buf_v(lnz*lnz*neqs))
    allocate(send_buf_w(lnz*lnz*neqs))
    allocate(recv_buf_v(lnz*lnz*neqs))
    allocate(recv_buf_w(lnz*lnz*neqs))

    do n = 1, neqs
        ! if (pid == 0) then
        !     print *, "matrix", ma(n)%matrix
        ! end if
        do j = 1, n_size - 1
            do i = j + 1, min(j + lnz, n_size)
                lij = ma(n)%matrix(lnz + 1 + j - i, i)/ ma(n)%matrix(lnz + 1, j)
                ma(n)%matrix(lnz + 1 + j - i, i) = lij
                do k = 1, nnz - 1 + j - i
                    ma(n)%matrix(lnz + 1 + j - i + k, i) = ma(n)%matrix(lnz + 1 + j - i + k, i) - &
                    lij * ma(n)%matrix(lnz + 1 + k, j)
                end do
            end do
        end do
        ! if (pid == 0) then
        !     print *, "a", ma(n)%matrix
        !     print *, "a(2, 2)", ma(n)%matrix(1, n_size) * ma(n)%matrix(3, n_size-1) + ma(n)%matrix(2, n_size)
        ! end if
        ma(n)%v(:, :) = 0
        ma(n)%w(:, :) = 0
        do i = 1, lnz
            do j = 1, i
                ma(n)%v(j, n_size-lnz+i) = ma(n)%matrix(nnz-i+j, n_size-lnz+i)
            end do
            do j = 1, lnz - i + 1
                ma(n)%w(j, i) = ma(n)%matrix(j, i)
            end do
        end do

        call struct_sptrsv(n_size, nnz, lnz, ma(n)%matrix, ma(n)%v)
        call struct_sptrsv(n_size, nnz, lnz, ma(n)%matrix, ma(n)%w)

        call mpi_pack(ma(n)%v(1, n_size-lnz+1), lnz*lnz, MD_MPI_TYPE, send_buf_v, buf_size, pos_v, comm, ierr)
        call mpi_pack(ma(n)%w(1, 1), lnz*lnz, MD_MPI_TYPE, send_buf_w, buf_size, pos_w, comm, ierr)
        
    end do

    ! MPI communication
    ! Send vb & wt
    call mpi_sendrecv(send_buf_v, buf_count, MD_MPI_TYPE, next_id, 0, &
    recv_buf_v, buf_count, MD_MPI_TYPE, prev_id, 0, comm, status, ierr)

    call mpi_sendrecv(send_buf_w, buf_count, MD_MPI_TYPE, prev_id, 0, &
    recv_buf_w, buf_count, MD_MPI_TYPE, next_id, 0, comm, status, ierr)

    pos_v = 0
    pos_w = 0
    do n = 1, neqs
        ma(n) = ma(n)
        call mpi_unpack(recv_buf_v, buf_size, pos_v, ma(n)%v_prev(1, 1), lnz*lnz, MD_MPI_TYPE, comm, ierr)
        call mpi_unpack(recv_buf_w, buf_size, pos_w, ma(n)%w_next(1, 1), lnz*lnz, MD_MPI_TYPE, comm, ierr)
    end do

    deallocate(send_buf_v)
    deallocate(send_buf_w)
    deallocate(recv_buf_v)
    deallocate(recv_buf_w)

    if (pid <= 3) then
        ! print *, pid, "wt", ma(1)%w_next
        ! print *, pid, "vb", ma(1)%v_prev
        ! print *, pid, "w", ma(1)%w
        ! print *, pid, "v", ma(1)%v
    end if

end subroutine

subroutine md_cyclic_solver_kernel(neqs, n_size, nrhs, ma, rhs, x)
    integer, intent(in) :: neqs !纬圈三对角方程个数
    integer, intent(in) :: n_size !当前进程对应三对角方程的行数
    integer, intent(in) :: nrhs !每个三对角方程的右端向量个数
    type(md_matrix), intent(in) :: ma(neqs) !预处理结构体，每个方程对应一个
    md_type, intent(in) :: rhs(nrhs, n_size, neqs) !求解右端向量
    md_type, intent(out) :: x(nrhs, n_size, neqs) !求解结果

    integer :: buf_count, buf_size, pos_b, pos_t
    integer :: status(mpi_status_size), ierr
    md_type, allocatable, dimension(:) :: send_buf_b, recv_buf_b
    md_type, allocatable, dimension(:) :: send_buf_t, recv_buf_t

    md_type, allocatable, dimension(:, :) :: x_next, x_prev

    !!! For trd solver       !!!
    !!! | 1 b1 | |x1| = |g1| !!!
    !!! | b2 1 | |x2| = |g2| !!!
    md_type :: b1, b2, g1(nrhs), g2(nrhs)

    integer :: n, i
    integer :: lnz, nnz

    x(:, :, :) = rhs(:, :, :)

    lnz = ma(1) % lnz
    nnz = ma(1) % nnz

    do n = 2, neqs
        if (ma(n)%nnz .ne. nnz .or. ma(n)%n_size .ne. n_size) then
            print *, "Error: matrices of md_cyclic_solver must have the same size!"
            return
        end if
    end do

    ! DY = B
    do n = 1, neqs
        call struct_sptrsv(n_size, nnz, nrhs, ma(n)%matrix, x(1, 1, n))
    end do

    allocate(send_buf_b(nrhs*lnz*neqs))
    allocate(recv_buf_b(nrhs*lnz*neqs))
    allocate(send_buf_t(nrhs*lnz*neqs))
    allocate(recv_buf_t(nrhs*lnz*neqs))
    allocate(x_next(nrhs, lnz))
    allocate(x_prev(nrhs, lnz))

    buf_count = nrhs * lnz * neqs
    buf_size = buf_count * md_size
    pos_b = 0; pos_t = 0

    ! MPI communication
    ! Send rhs_b & rhs_t
    do n = 1, neqs
        call mpi_pack(x(1, 1, n), nrhs*lnz, MD_MPI_TYPE, send_buf_t, buf_size, pos_t, comm, ierr)
        call mpi_pack(x(1, n_size-lnz+1, n), nrhs*lnz, MD_MPI_TYPE, send_buf_b, buf_size, pos_b, comm, ierr)
    end do

    call mpi_sendrecv(send_buf_b, buf_count, MD_MPI_TYPE, next_id, 0, &
    recv_buf_b, buf_count, MD_MPI_TYPE, prev_id, 0, comm, status, ierr)

    call mpi_sendrecv(send_buf_t, buf_count, MD_MPI_TYPE, prev_id, 0, &
    recv_buf_t, buf_count, MD_MPI_TYPE, next_id, 0, comm, status, ierr)

    if (lnz == 1) then
        do n = 1, neqs
            ! bottom
            b1 = ma(n)%v(1, n_size)
            b2 = ma(n)%w_next(1, 1)
            g1(:) = x(:, n_size, n)
            g2(:) = recv_buf_t((n-1)*nrhs+1:n*nrhs)
            x(:, n_size, n) = (g1(:) - b1 * g2(:)) / (1 - b1 * b2)
            x_next(:, 1) = (g2(:) - b2 * g1(:)) / (1 - b1 * b2)

            ! top
            b1 = ma(n)%v_prev(1, 1)
            b2 = ma(n)%w(1, 1)
            g1(:) = recv_buf_b((n-1)*nrhs+1:n*nrhs)
            g2(:) = x(:, 1, n)
            x_prev(:, 1) = (g1(:) - b1 * g2(:)) / (1 - b1 * b2)
            x(:, 1, n) = (g2(:) - b2 * g1(:)) / (1 - b1 * b2)

            ! Calculate X'
            do i = 2, n_size - 1
                x(:, i, n) = x(:, i, n) - &
                ma(n)%v(1, i) * x_next(:, 1) - ma(n)%w(1, i) * x_prev(:, 1)
            end do
        end do
    else
        print *, "Not implemented!"
    end if

    ! print *, pid, "b1", rhs(1, 1, 1)
    ! print *, pid, "test1", x_prev(1, 1)*test_a(1, 1, 1) + x(1, 1, 1)*test_a(2, 1, 1) + x(1, 2, 1)*test_a(3, 1, 1)
    ! print *, pid, "bn", rhs(1, n_size, 1)
    ! print *, pid, "test2", x(1, n_size-1, 1)*test_a(1, n_size, 1) + x(1, n_size, 1)*test_a(2, n_size, 1) + &
    ! x_next(1, 1)*test_a(3, n_size, 1)

    deallocate(send_buf_b)
    deallocate(recv_buf_b)
    deallocate(send_buf_t)
    deallocate(recv_buf_t)
    deallocate(x_next)
    deallocate(x_prev)

end subroutine

subroutine md_cyclic_multi_solver_kernel(nas, nbs, n_size, nrhs, pos, ma, rhs, x)
    integer, intent(in) :: nas !循环三对角矩阵A的个数
    integer, intent(in) :: nbs !右端向量b的个数
    integer, intent(in) :: n_size !当前进程对应三对角方程的行数
    integer, intent(in) :: nrhs !每个三对角方程的右端向量个数
    integer, intent(in) :: pos(nbs) !每个b对应的矩阵下标
    type(md_matrix), intent(in) :: ma(nas) !预处理结构体，每个方程对应一个
    md_type, intent(in) :: rhs(nrhs, n_size, nbs) !求解右端向量
    md_type, intent(out) :: x(nrhs, n_size, nbs) !求解结果

    integer :: buf_count, buf_size, pos_b, pos_t
    integer :: status(mpi_status_size), ierr
    md_type, allocatable, dimension(:) :: send_buf_b, recv_buf_b
    md_type, allocatable, dimension(:) :: send_buf_t, recv_buf_t

    md_type, allocatable, dimension(:, :) :: x_next, x_prev

    !!! For trd solver       !!!
    !!! | 1 b1 | |x1| = |g1| !!!
    !!! | b2 1 | |x2| = |g2| !!!
    md_type :: b1, b2, g1(nrhs), g2(nrhs)

    integer :: n, i
    integer :: lnz, nnz

    x(:, :, :) = rhs(:, :, :)

    lnz = ma(1) % lnz
    nnz = ma(1) % nnz

    do n = 2, nas
        if (ma(n)%nnz .ne. nnz .or. ma(n)%n_size .ne. n_size) then
            print *, "Error: matrices of md_cyclic_solver must have the same size!"
            return
        end if
    end do

    ! DY = B
    do n = 1, nbs
        call struct_sptrsv(n_size, nnz, nrhs, ma(pos(n))%matrix, x(1, 1, n))
    end do

    allocate(send_buf_b(nrhs*lnz*nbs))
    allocate(recv_buf_b(nrhs*lnz*nbs))
    allocate(send_buf_t(nrhs*lnz*nbs))
    allocate(recv_buf_t(nrhs*lnz*nbs))
    allocate(x_next(nrhs, lnz))
    allocate(x_prev(nrhs, lnz))

    buf_count = nrhs * lnz * nbs
    buf_size = buf_count * md_size
    pos_b = 0; pos_t = 0

    ! MPI communication
    ! Send rhs_b & rhs_t
    do n = 1, nbs
        call mpi_pack(x(1, 1, n), nrhs*lnz, MD_MPI_TYPE, send_buf_t, buf_size, pos_t, comm, ierr)
        call mpi_pack(x(1, n_size-lnz+1, n), nrhs*lnz, MD_MPI_TYPE, send_buf_b, buf_size, pos_b, comm, ierr)
    end do

    call mpi_sendrecv(send_buf_b, buf_count, MD_MPI_TYPE, next_id, 0, &
    recv_buf_b, buf_count, MD_MPI_TYPE, prev_id, 0, comm, status, ierr)

    call mpi_sendrecv(send_buf_t, buf_count, MD_MPI_TYPE, prev_id, 0, &
    recv_buf_t, buf_count, MD_MPI_TYPE, next_id, 0, comm, status, ierr)

    if (lnz == 1) then
        do n = 1, nbs
            ! bottom
            b1 = ma(pos(n))%v(1, n_size)
            b2 = ma(pos(n))%w_next(1, 1)
            g1(:) = x(:, n_size, n)
            g2(:) = recv_buf_t((n-1)*nrhs+1:n*nrhs)
            x(:, n_size, n) = (g1(:) - b1 * g2(:)) / (1 - b1 * b2)
            x_next(:, 1) = (g2(:) - b2 * g1(:)) / (1 - b1 * b2)

            ! top
            b1 = ma(pos(n))%v_prev(1, 1)
            b2 = ma(pos(n))%w(1, 1)
            g1(:) = recv_buf_b((n-1)*nrhs+1:n*nrhs)
            g2(:) = x(:, 1, n)
            x_prev(:, 1) = (g1(:) - b1 * g2(:)) / (1 - b1 * b2)
            x(:, 1, n) = (g2(:) - b2 * g1(:)) / (1 - b1 * b2)

            ! Calculate X'
            do i = 2, n_size - 1
                x(:, i, n) = x(:, i, n) - &
                ma(pos(n))%v(1, i) * x_next(:, 1) - ma(pos(n))%w(1, i) * x_prev(:, 1)
            end do
        end do
    else
        print *, "Not implemented!"
    end if

    deallocate(send_buf_b)
    deallocate(recv_buf_b)
    deallocate(send_buf_t)
    deallocate(recv_buf_t)
    deallocate(x_next)
    deallocate(x_prev)

end subroutine
