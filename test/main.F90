#define md_type real(8)
#define md_size 8

program test
    
    use multi_diag_mod
    use mpi

    implicit none

    integer :: n_size
    integer :: nnz
    integer :: nrhs
    integer :: neqs
    integer :: pid, pnum, ierr, next, prev
    integer :: comm

    type(md_matrix), allocatable, dimension(:) :: matrix
    md_type, allocatable, dimension(:, :, :) :: ma
    md_type, allocatable, dimension(:, :, :) :: b, x

    md_type, allocatable, dimension(:, :) :: a2, b2, c2, d2, x2
    md_type :: tmp

    integer :: i, j, n

    md_type :: test_b, test_x
    integer :: test_flag

    n_size = 60
    nnz = 3
    nrhs = 10
    neqs = 10

    comm = MPI_COMM_WORLD
    
    call mpi_init(ierr)
    call mpi_comm_rank(comm, pid, ierr)
    call mpi_comm_size(comm, pnum, ierr)
    call pos2id(pid, pnum, next, prev)

    call md_init(comm, pid, next, prev)

    !!!! 纬圈三对角测试

    allocate(matrix(neqs))
    allocate(ma(nnz, n_size, neqs))
    allocate(b(nrhs, n_size, neqs))
    allocate(x(nrhs, n_size, neqs))
    do i = 1, neqs
        call allocate_md_matrix(n_size, nnz, matrix(i))
        call init_data(nnz, n_size, nrhs, ma(1, 1, i), b(1, 1, i), pid * neqs + i)
        matrix(i)%matrix(:, :) = ma(:, :, i)
    end do
    call md_cyclic_fact(neqs, matrix)
    call md_cyclic_solver(neqs, n_size, nrhs, matrix, b, x)

    test_flag = 0

    do n = 1, neqs
        do i = 2, n_size - 1
            do j = 1, nrhs
                if (pid >= 0) then
                    test_x = x(j, i-1, n) * ma(1, i, n) + & 
                    x(j, i, n) * ma(2, i, n) + x(j, i+1, n) * ma(3, i, n)
                    test_b = b(j, i, n)
                    !print *, pid, "x backward", i, test_x - test_b
                    if (test_x - test_b >= 0.00000000001d0) then
                        test_flag = i
                    end if
                end if
            end do
        end do
    end do

    if (test_flag == 0) then
        print *, pid, "Parallel Correct!"
    else
        print *, pid, "Parallel error at", test_flag
    end if

    !!!!!
    
    !!!! 串行

    n_size = 64
    neqs = 1
    allocate(a2(neqs, n_size))
    allocate(b2(neqs, n_size))
    allocate(c2(neqs, n_size))
    allocate(d2(neqs, n_size))
    allocate(x2(neqs, n_size))

    do i = 1, n_size
        do j = 1, neqs
            call random_number(tmp)
            a2(j, i) = 0.05 + 0.25 * tmp
            call random_number(tmp)
            b2(j, i) = 1 + 0.2 * tmp
            call random_number(tmp)
            c2(j, i) = 0.05 + 0.25 * tmp
            call random_number(tmp)
            d2(j, i) = 100 * tmp
        end do
    end do

    call trd_solver(neqs, n_size, a2, b2, c2, d2, x2)

    test_flag = 0
    do i = 2, n_size - 1
        do j = 1, neqs
            test_x = x2(j, i-1) * a2(j, i) + x2(j, i) * b2(j, i) &
            + x2(j, i+1) * c2(j, i)
            test_b = d2(j, i)
            if (test_x - test_b >= 0.00000000001d0) then
                test_flag = i
            end if
        end do
    end do

    if (test_flag == 0) then
        print *, pid, "Sequential Correct!"
    else
        print *, pid, "Sequential error at", test_flag
    end if

    !!!!!


    do i = 1, neqs
        call deallocate_md_matrix(matrix(i))
    end do
    
contains
    subroutine init_data(nnz, p_size, nrhs, ma, rhs, seed_input)

        integer, intent(in) :: nnz, p_size, nrhs
        md_type, intent(out) :: ma(nnz, p_size)
        md_type, intent(out) :: rhs(nrhs, p_size)
        integer, intent(in) :: seed_input

        integer :: i, j
        integer :: lnz
        md_type :: tmp
        integer :: seed_put(1)

        lnz = (nnz - 1) / 2

        seed_put(1) = seed_input

        call random_seed(put = seed_put)
        call random_number(tmp)
        do i = 1, p_size
            do j = 1, lnz
                call random_number(tmp)
                ma(j, i) = 0.05 + 0.25 * tmp
            end do
            call random_number(tmp)
            ma(lnz + 1, i) = 1 + 0.2 * tmp
            do j = lnz + 2, nnz
                call random_number(tmp)
                ma(j, i) = 0.05 + 0.25 * tmp
            end do
        end do

        do i = 1, p_size
            do j = 1, nrhs
                call random_number(tmp)
                rhs(j, i) = 50 + 50 * tmp
            end do
        end do

        ! if (pid == 0) then
        !     print *, rhs(1, 1)
        ! end if

    end subroutine

    subroutine pos2id(pid, pnum, next, prev)

        integer, intent(in) :: pid, pnum
        integer, intent(out) :: next, prev

        next = mod(pid + 1 + pnum, pnum)
        prev = mod(pid - 1 + pnum, pnum)

    end subroutine

end program

#undef md_type
#undef md_size