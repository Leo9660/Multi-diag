! structure management
subroutine allocate_md_matrix(n_size, nnz, ma)
    integer, intent(in) :: n_size
    integer, intent(in) :: nnz
    type(md_matrix), intent(inout) :: ma

    ma%n_size = n_size
    ma%nnz = nnz
    ma%lnz = (nnz - 1) / 2
    allocate(ma%matrix(ma%nnz, n_size))
    allocate(ma%w(ma%lnz, n_size))
    allocate(ma%v(ma%lnz, n_size))
    allocate(ma%w_next(ma%lnz, ma%lnz))
    allocate(ma%v_prev(ma%lnz, ma%lnz))
    
end subroutine

subroutine deallocate_md_matrix(ma)
    type(md_matrix), intent(inout) :: ma
    deallocate(ma%matrix)
    deallocate(ma%w)
    deallocate(ma%v)
    deallocate(ma%w_next)
    deallocate(ma%v_prev)
end subroutine
