!---------------------------------------------------------
! User
subroutine User(Pe,beta,N,G,printout,method)
    implicit NONE
    
    integer :: Pe_select,N_select,method_select,N,G,printout
    real*8 :: Pe,beta
    character :: method
    
    
    write(*,*)'Choose the number you want for Pe.'
    write(*,*)'(1) 100'
    write(*,*)'(2) 500'
    write(*,*)'(3) 1000'
    write(*,*)'(4) 2000'
    read(*,*)Pe_select
    select case (Pe_select)
        case(1)
            Pe = 100.d0
        case(2)
            Pe = 500.d0
        case(3)
            Pe = 1000.d0
        case(4)
            Pe = 2000.d0
        case default
            write(*,*)'Please follow the instructions~'
            stop
    end select
        
!--------------------------------------------------------------
    write(*,*)'How many grid points do you want?'
    write(*,*)'(1) 81 by 81'
    write(*,*)'(2) else'
    read(*,*)N_select
    select case (N_select)
        case(1)
            N = 81
        case(2)
            write(*,*)'Enter the number of grid points you want.'
            read(*,*)N
        case default
            write(*,*)'Please follow the instructions~'
            stop
    end select
    
!---------------------------------------------------------------
    write(*,*)'Choose a method to solve the problem.'
    write(*,*)'(1) Explicit'
    write(*,*)'(2) Implicit'
    read(*,*)method_select
    select case (method_select)
        case(1)
            method = 'E'
            write(*,*)'Choose a method for explicit scheme.'
            write(*,*)'(1) FOU'
            write(*,*)'(2) SOCD'
            write(*,*)'(3) QUICK'
            write(*,*)'(4) TOU'
            write(*,*)'(5) SOU'
            read(*,*)method_select
            select case (method_select)
                case(1)
                    G = 0
                    beta = 0.d0
                    printout = 1
                case(2)
                    G = 1
                    beta = 0.d0
                    printout = 2
                case(3)
                    G = 1
                    beta = 1.d0/8.d0
                    printout = 3
                case(4)
                    G = 1
                    beta = 1.d0/6.d0
                    printout = 4
                case(5)
                    G = 1
                    beta = 1.d0/2.d0
                    printout = 5
                case default
                    write(*,*)'Please follow the instructions~'
                    stop
            end select
        case(2)
            method = 'I'
            write(*,*)'Choose a method for implicit scheme.'
            write(*,*)'(1) FOU'
            write(*,*)'(2) SOCD'
            write(*,*)'(3) QUICK'
            write(*,*)'(4) TOU'
            write(*,*)'(5) SOU'
            read(*,*)method_select
            select case (method_select)
                case(1)
                    G = 0
                    beta = 0.d0
                    printout = 1
                case(2)
                    G = 1
                    beta = 0.d0
                    printout = 2
                case(3)
                    G = 1
                    beta = 1.d0/8.d0
                    printout = 3
                case(4)
                    G = 1
                    beta = 1.d0/6.d0
                    printout = 4
                case(5)
                    G = 1
                    beta = 1.d0/2.d0
                    printout = 5
                case default
                    write(*,*)'Please follow the instructions~'
                    stop
            end select
        case default
        write(*,*)'Please follow the instructions~'
        stop
    end select
    
        
end subroutine User