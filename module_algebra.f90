module module_algebra

contains
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------- inv_mat --------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine inv_mat(A,A_inv,nmat)
    ! computes matrix inverse
    implicit none
    integer :: nmat ! size of the matrix
    real(kind=8) :: A(nmat,nmat), A_inv(nmat,nmat)
    integer :: ipiv(nmat), info
    real(kind=8) :: work(nmat*nmat)
    integer :: lwork
    
    lwork = nmat*nmat
    A_inv = A

    ! Call LAPACK routine to invert the matrix
    call dgetrf(nmat, nmat, A_inv, nmat, ipiv, info)
    if (info /= 0) then
        print *, "Error: Matrix is singular."
        stop
    endif

    ! Compute the inverse of A_inv
    call dgetri(nmat, A_inv, nmat, ipiv, work, lwork, info)
    if (info /= 0) then
        print *, "Error: Unable to invert the matrix."
        stop
    endif
end subroutine inv_mat

end module module_algebra