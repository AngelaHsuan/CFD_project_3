!------------------------------------------------
!animation
subroutine animation(del_x,del_y,N,T,printout,tt)
    implicit NONE
    
    integer :: i,j,N,Pad,printout
    real*8 :: del_x,del_y,tt
    real*8, dimension(N,N) :: T
    
    Pad = 100
    
    open (unit=Pad, file="animation.txt", status="UNKNOWN", ACCESS='APPEND')
    
    write(Pad,*) 'variables="x","y","T"'
    write(Pad,*) 'zone i=', N,'j=',N,'DATAPACKING=POINT'
    
    do j = 1,N,1
        do i = 1,N,1
            write(Pad,'(3F20.12)') -4.d0+(i-1)*del_x, -4.d0+(j-1)*del_y, T(i,j)
        end do
    end do
    close (Pad,status = 'Keep')
    
end subroutine animation