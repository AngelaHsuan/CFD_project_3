!------------------------------------------------
!output
subroutine output(del_x,del_y,N,T,printout,tt)
    implicit NONE
    
    integer :: i,j,N,Pad,printout
    real*8 :: del_x,del_y,tt
    real*8, dimension(N,N) :: T
    
    Pad = 100
    
    if (printout == 1) then
        open (unit=Pad, file="FOU.txt", status="UNKNOWN")
    else if (printout == 2) then
        open (unit=Pad, file="SOCD.txt", status="UNKNOWN")
    else if (printout == 3) then
        open (unit=Pad, file="QUICK.txt", status="UNKNOWN")
    else if (printout == 4) then
        open (unit=Pad, file="TOU.txt", status="UNKNOWN")
    else if (printout == 5) then
        open (unit=Pad, file="SOU.txt", status="UNKNOWN")
    else
        write(*,*)'there is error in printout'
        pause
    end if
    
    write(Pad,*) 'variables="x","y","T"'
    write(Pad,*) 'zone i=', N,'j=',N,'DATAPACKING=POINT'
    
    do j = 1,N,1
        do i = 1,N,1
            write(Pad,'(3F30.12)') -4.d0+(i-1)*del_x, -4.d0+(j-1)*del_y, T(i,j)
        end do
    end do
    close (Pad,status = 'Keep')
    
end subroutine output