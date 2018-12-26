!------------------------------------------------
!velocity
subroutine velocity(N,u,v,del_x, del_y,v_max)
    implicit NONE
    
    integer :: i, j, N, Pad
    real*8 :: del_x, del_y,r,freq,v_max
    real*8, dimension(N,N) :: u, v
    
    do i = 1,N,1
        do j = 1, N, 1
            if (i == (N+1)/2 .AND. j == (N+1)/2) then
                u(i,j) = 0.d0
                v(i,j) = 0.d0
            else
                r = sqrt((-4 + (i-1)*del_x)**2 + (-4 + (j-1)*del_y)**2)
                freq = tanh(r)/((cosh(r))**2)/r/v_max
                u(i,j) = -freq*(-4.d0 + (j-1)*del_y)
                v(i,j) = freq*(-4.d0 + (i-1)*del_x)
            end if
        end do
    end do
    
    
    Pad = 100
    open(unit = Pad, file = "velocity.txt", status = "UNKNOWN")
    write(Pad, *)'variables="x","y","u","v"'
    write(Pad,*) 'zone i=', N,'j=',N,'DATAPACKING=POINT'
    
    do j = 1,N,1
        do i = 1,N,1
            write(Pad,'(4F20.12)') -4.d0+(i-1)*del_x, -4.d0+(j-1)*del_y, u(i,j), v(i,j)
        end do
    end do
    close (Pad,status = 'Keep')
    
end subroutine velocity