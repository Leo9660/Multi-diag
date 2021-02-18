#define md_type real(8)
#define md_size 8
#if (md_size == 8)
#define MD_MPI_TYPE MPI_REAL8
#else
#define MD_MPI_TYPE MPI_REAL4
#endif

module multi_diag_mod

    use mpi

    implicit none

    private

    ! process id
    integer :: pid
    integer :: next_id
    integer :: prev_id

    ! whether cyclic solver is on a single process
    logical :: single_process

    ! communication group
    integer :: comm

    type md_matrix
        integer :: n_size !矩阵规模
        integer :: nnz !总对角元个数
        integer :: lnz !副对角元个数
        md_type, allocatable, dimension(:, :) :: matrix
        !矩阵数据(nnz, n_size)，1~lnz为左副对角元，第lnz+1个为主对角元，lnz+2~nnz为右副对角元
        md_type, allocatable, dimension(:, :) :: w, v
        !SPIKE分解用数据(lnz, n_size)，分别代表SPIKE分解A=DS后，S矩阵单位阵左/右的条状数据
        !单进程中，用v和w分别存储循环三对角LU分解后L矩阵右端的竖状数据和L矩阵下端的横状数据
        md_type, allocatable, dimension(:, :) :: w_next, v_prev
        !SPIKE分解用数据(lnz, lnz)，w_next为下个划分的wt，v_prev为上个划分的vb
    end type md_matrix

    !private pid, next_id, prev_id, single_process, comm
    public md_init, md_cyclic_fact, md_cyclic_solver, &
    md_cyclic_multi_solver, trd_solver, md_matrix, &
    allocate_md_matrix, deallocate_md_matrix

contains
    subroutine md_init(comm_in, rank, next_rank, prev_rank)

        integer, intent(in) :: comm_in !纬圈三对角通信子
        integer, intent(in) :: rank, next_rank, prev_rank!当前进程id，纬圈下一个/上一个进程id

        comm = comm_in

        pid = rank
        next_id = next_rank
        prev_id = prev_rank

        if (pid == next_id .and. pid == prev_rank) then
            single_process = .true.
        else if (pid .ne. next_id .and. pid .ne. prev_rank) then
            single_process = .false.
        else
            print *, "Error: Invalid init args!"
        end if

    end subroutine

    subroutine md_cyclic_fact(neqs, ma)
        integer, intent(in) :: neqs !纬圈三对角方程个数
        type(md_matrix), intent(inout) :: ma(neqs) !预处理结构体，每个方程对应一个

        call md_cyclic_fact_kernel(neqs, ma)

    end subroutine

    ! 以下两个接口均为循环三对角求解接口
    ! 根据初始化参数，多进程MPI并行求解循环三对角方程组Ax = b

    ! 批量求解neqs个循环三对角方程组
    subroutine md_cyclic_solver(neqs, n_size, nrhs, ma, rhs, x)
        integer, intent(in) :: neqs !循环三对角方程个数
        integer, intent(in) :: n_size !当前进程对应三对角方程的行数
        integer, intent(in) :: nrhs !每个三对角方程的右端向量个数
        type(md_matrix), intent(in) :: ma(neqs) !预处理结构体，每个方程对应一个
        md_type, intent(in) :: rhs(nrhs, n_size, neqs) !求解右端向量
        md_type, intent(out) :: x(nrhs, n_size, neqs) !求解结果

        call md_cyclic_solver_kernel(neqs, n_size, nrhs, ma, rhs, x)

    end subroutine

    ! 批量求解循环三对角方程组，但部分A矩阵可能相同
    ! 此时右端向量b的个数多于A矩阵的个数，用pos数组指定每个b对应的矩阵下标
    ! 比如，同时求解A1x1 = b1, A1x2 = b2, A2x3 = b3
    ! 则nas = 2，nbs = 3，pos数组为[1, 1, 2]
    subroutine md_cyclic_multi_solver(nas, nbs, n_size, nrhs, pos, ma, rhs, x)
        integer, intent(in) :: nas !循环三对角矩阵A的个数
        integer, intent(in) :: nbs !右端向量b的个数
        integer, intent(in) :: n_size !当前进程对应三对角方程的行数
        integer, intent(in) :: nrhs !每个三对角方程的右端向量个数
        integer, intent(in) :: pos(nbs) !每个b对应的矩阵下标
        type(md_matrix), intent(in) :: ma(nas) !预处理结构体，每个方程对应一个
        md_type, intent(in) :: rhs(nrhs, n_size, nbs) !求解右端向量
        md_type, intent(out) :: x(nrhs, n_size, nbs) !求解结果

        call md_cyclic_multi_solver_kernel(nas, nbs, n_size, nrhs, pos, ma, rhs, x)
    end subroutine

    subroutine trd_solver(neqs, nsize, a, b, c, rhs, x)
        integer, intent(in) :: neqs !垂向三对角方程个数
        integer, intent(in) :: nsize !垂向三对角方程的行数
        md_type, intent(in) :: a(neqs, nsize), b(neqs, nsize), c(neqs, nsize) !代表对角元和两个非对角元
        md_type, intent(in) :: rhs(neqs, nsize) !求解右端向量
        md_type, intent(out) :: x(neqs, nsize) !求解结果

        call trd_single_solver_kernel(neqs, nsize, a, b, c, rhs, x)
    end subroutine

#include "cyclic.F90"

#include "memory.F90"

#include "struct_sptrsv.F90"

#include "trd_single.F90"

end module

#undef md_type
#undef md_size
