! Written by Ting-Hsuan Hsu

! Start date 12 09 2017
! End date 12 20 2017

! Trying to use numerical method to solve the transient convection-diffusion problem
!------------------------------------------------------------------------------------------------------------------------------
program Main
    implicit NONE
    
    integer :: printout,N,G
    real*8 :: Pe,del_x, del_y,conv,w,beta,del_t,v_max,Time_Exe
    integer Time_Start,Time_End,Rate
    real*8, allocatable, dimension(:,:) :: T

    character :: again,method
    
    again = 'y'
    do while (again == 'y')
        
        call User(Pe,beta,N,G,printout,method)
        allocate(T(N,N))
        call Initial(N,del_x,del_y,conv,w,del_t,T,v_max,method,G)
        call system_clock(Time_Start,Rate)
        
        if (method == 'E') then
            call FTCS(del_t,del_x,del_y,G,T,N,v_max,Pe,beta,printout)
        else if (method == 'I') then
            call ADI(conv,w,del_t,del_x,del_y,G,T,N,v_max,Pe,beta,printout)
        else
            write(*,*)'There is mistake on method(Main)!!'
            pause
        end if
        
        call system_clock(Time_End,Rate)
        Time_Exe = (real(Time_End)-real(Time_Start))/real(Rate)
        write(*,*) '------------------------------------------------------------'
        Write(*,*)'The computational time(sec):', Time_Exe
        write(*,*) '------------------------------------------------------------'

        write(*,*)'Do you want to do it again?(y/n)'
        read(*,*)again
        deallocate(T)
    end do
    
end program Main