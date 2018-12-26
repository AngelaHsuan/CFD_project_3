!---------------------------------------------------------
! Boundary
subroutine Boundary(N,tt,del_t,del_x,del_y,T_left,T_right,T_bottom,T_top,v_max)
    implicit NONE
    
    integer :: N,i,j
    real*8 :: del_t,del_x,del_y,r,freq,v_max,tt
    real*8, dimension(N) :: T_left,T_right,T_bottom,T_top
    
    do i = 1,N,1
        do j = 1,N,1
            r = sqrt((-4 + (i-1)*del_x)**2 + (-4 + (j-1)*del_y)**2)
            freq = tanh(r)/((cosh(r))**2)/r/v_max
            
            if (j == 1) then
                !bottom
                T_bottom(i) = -tanh(-2.d0*(-4 + (j-1)*del_y )*cos(freq*tt) - 0.5d0*sin(freq*tt) )
            else if (j == N) then
                !top
                T_top(i) = -tanh(2.d0*(-4 + (j-1)*del_y )*cos(freq*tt) - 0.5d0*sin(freq*tt) )
            else if (i == 1) then
                !left
                T_left(j) = -tanh(0.5d0*(-4 + (j-1)*del_y )*cos(freq*tt) + 2.d0*sin(freq*tt) )
            else if (i == N) then
                !right
                T_right(j) = -tanh(0.5d0*(-4 + (j-1)*del_y )*cos(freq*tt) - 2.d0*sin(freq*tt) )
            else
            end if
            
            
        end do
        !if (tt == 120) then
        !    write(*,*)T_bottom(i),i
        !end if
    end do
    !if (tt == 120) then
    !    pause
    !end if
    
    !write(*,*)T_left
    !pause
    !write(*,*)T_right
    !pause
    !write(*,*)T_bottom
    !pause
    !write(*,*)T_top
    !pause
    
    
end subroutine Boundary