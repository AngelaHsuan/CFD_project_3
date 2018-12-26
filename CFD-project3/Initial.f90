!---------------------------------------------------------
! Initial
subroutine Initial(N,del_x,del_y,conv,w,del_t,T,v_max,method,G)
    implicit NONE
    
    integer :: N,i,j,G
    real*8 :: Lx,Ly,del_x,del_y,conv,w,del_t,v_max
    real*8, dimension(N,N) :: T
    character :: method
    
    conv = 1D-6
    Lx = 8.d0
    Ly = 8.d0
    del_x = Lx/(N-1)
    del_y = Ly/(N-1)
    v_max = 0.385d0
    
    do j = 1,N,1
        T(:,j) = -tanh(0.5d0*(-4 + (j-1)*del_y ))
    end do
    !do i = 1,N,1
    !    do j = 1,N,1
    !        write(*,*)T(j,i),j,i
    !    end do
    !end do
    !pause
    !----------------------------------------------------
    w = 1.d0
    if (method == 'I' .AND. G.NE.0) then
        write(*,*)'Please enter the relaxation coefficient.'
        read(*,*)w
    end if
    !----------------------------------------------------
    write(*,*)'Please enter the time step.'
    read(*,*)del_t
    
end subroutine Initial