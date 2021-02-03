# Multi-diag
Parallel multi-diagonal solver module in Fortran

## v1.0
md_mod模块提供多进程并行循环三对角方程组Ax = b的快速近似求解方法，注意此时矩阵大小n应该远大于对角线条数（目前支持为3），否则可能出现数值精度问题。md_type和md_size用于指定库的数据类型和数据类型种别值，注意应保持一致，比如md_type为rea(4)时，md_size为4。

## 数据结构

令矩阵规模为n，对角线条数为nnz，下三角部分对角线条数为lnz，注意lnz = (nnz - 1) / 2。在三对角矩阵中，nnz = 3，lnz = 1。目前只支持三对角情形。
md_mod模块定义了一个数据结构md_matrix，该数据结构用于存储三对角矩阵分解后的信息，包括每个子矩阵的LU分解和SPIKE划分的相关数据。用户需要修改的成员变量为matrix数组，该数组维度为(nnz, n)，即每次用连续的nnz个元素表示一行，其中第1 ~ lnz个元素为左副对角元，第lnz+1个元素为主对角元，第lnz+2 ~ nnz个元素为右副对角元。

## 内存申请

调用allocate_trd_prec(n_size, nnz, prec)对结构体进行内存申请，其中n_size为矩阵规模，nnz为对角线条数，prec为要申请的md_matrix数据结构。
调用deallocate_trd_prec(prec)释放对应的结构体内存。

## 模块初始化

调用md_init(comm_in, rank, next_rank, prev_rank)进行初始化，其中comm_in为并行通信子，rank为当前进程id，next_rank和prev_rank分别为下一个和上一个进程的id。后面提到的分解和求解接口均是多进程同时运行，每个进程对应循环三对角矩阵中的n_size个连续行，与其矩阵分块相邻的下/上进程对应next_rank和prev_rank。

## 矩阵分解

正式求解前，先调用md_lat_fact(neqs, ma)对循环三对角矩阵进行分解操作，其中neqs为三对角矩阵个数，ma是维度为(neqs)的md_matrix结构体数组。分解结果存在md_mod自定义结构体md_matrix中，需要用户预分配空间为并保存分解结果。每个三对角方程组对应一个md_matrix结构体。

## 求解接口

调用md_cyclic_solver(neqs, n_size, nrhs, ma, rhs, x)求解neqs个规模为n_size，右端向量个数为nrhs的并行循环三对角方程组，其中ma是维度为(neqs)的md_matrix结构体数组，rhs和x均是维度为(nrhs, n_size, neqs)的实数数组，对应右端向量输入和求解输出。
<code fortran>
subroutine md_cyclic_solver(neqs, n_size, nrhs, ma, rhs, x)
    integer, intent(in) :: neqs !循环三对角方程个数
    integer, intent(in) :: n_size !当前进程对应三对角方程的行数
    integer, intent(in) :: nrhs !每个三对角方程的右端向量个数
    type(md_matrix), intent(in) :: ma(neqs) !预处理结构体，每个方程对应一个
    md_type, intent(in) :: rhs(nrhs, n_size, neqs) !求解右端向量
    md_type, intent(out) :: x(nrhs, n_size, neqs) !求解结果
end subroutine
</code>

调用md_cyclic_multi_solver(nas, nbs, n_size, nrhs, pos, ma, rhs, x)求解矩阵A的个数少于右端向量b的个数的情况，其中nas是矩阵A的个数，nbs是右端向量b的个数，pos是维度为()
