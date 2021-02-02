# Multi-diag
Parallel multi-diagonal solver module in Fortran

v1.0
md_mod模块提供多进程并行循环三对角方程组Ax = b的快速近似求解方法，注意此时矩阵大小n应该远大于对角线条数（目前支持为3），否则可能出现数值精度问题。md_type和md_size用于指定库的数据类型和数据类型种别值，注意应保持一致，比如md_type为rea(4)时，md_size为4。

令矩阵规模为n，对角线条数为nnz，下三角部分对角线条数为lnz，注意lnz = (nnz - 1) / 2。在三对角矩阵中，nnz = 3，lnz = 1。目前只支持三对角情形。
md_mod模块定义了一个数据结构md_matrix，该数据结构用于存储三对角矩阵分解后的信息，包括每个子矩阵的LU分解和SPIKE划分的相关数据。用户需要修改的成员变量为matrix数组，该数组维度为(nnz, n)，即每次用连续的nnz个元素表示一行，其中第1 ~ lnz个元素为左副对角元，第lnz+1个元素为主对角元，第lnz+2 ~ nnz个元素为右副对角元。
## 申请
