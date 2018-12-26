!---------------------------------------------------------
! FTCS
subroutine FTCS(del_t,del_x,del_y,G,T,N,v_max,Pe,beta,printout)
    implicit NONE
    
    integer :: i,j,k,nn,G,N,printout,Pad
    real*8 :: del_t,del_x,del_y,v_max,Pe,beta,tt
    real*8 :: C_w_FOU,C_px_FOU,C_py_FOU,C_e_FOU,C_s_FOU,C_n_FOU,C_px_gen,C_py_gen,C_w_gen,C_e_gen,C_s_gen,C_n_gen,C_ww_HO,C_ee_HO,C_nn_HO,C_ss_HO
    real*8, dimension(N) :: T_left,T_right,T_bottom,T_top
    real*8, dimension(N,N) :: T,T_after,u,v
    
    call velocity (N,u,v,del_x, del_y,v_max)
    
    T_after = T
    
    !open (unit=Pad, file="T_bottom.txt", status="UNKNOWN")
    !write(Pad,*) 'variables="t","Temp"'
    !write(Pad,*) 'zone i=', 4.d0/del_t + 1,'DATAPACKING=POINT'
    !close (Pad,status = 'Keep')
    
    do nn = 1,4.d0/del_t + 1,1
        tt = (nn-1)*del_t
        !write(*,'(1F6.4)')tt
        call Boundary(N,tt,del_t,del_x,del_y,T_left,T_right,T_bottom,T_top,v_max)
        
        !Pad = 100
        !open (unit=Pad, file="T_bottom.txt", status="UNKNOWN", ACCESS='APPEND')
        !write(Pad,'(3F25.12)') (nn-1)*del_t, T_bottom(N/2)
        !close (Pad,status = 'Keep')
        
        do j = 1,N,1
            do i = 1,N,1
                if (i == 1) then
                    !left
                    T(i,j) = T_left(j)
                else if(i == N) then
                    !right
                    T(i,j) = T_right(j)
                else if(j == 1) then
                    !bottom
                    T(i,j) = T_bottom(i)
                else if(j == N) then
                    !top
                    T(i,j) = T_top(i)
                else
                    call Coeff_E(i,j,Pe,beta,u,v,N,C_w_FOU,C_px_FOU,C_py_FOU,C_e_FOU,C_s_FOU,C_n_FOU,C_px_gen,C_py_gen,C_w_gen,C_e_gen,C_s_gen,C_n_gen,C_ww_HO,C_ee_HO,C_nn_HO,C_ss_HO,del_x,del_y)
                    if (i == 2 .OR. i == N-1 .OR. j == 2 .OR. j == N-1) then
                        T_after(i,j) = T(i,j)&
                            & + del_t*(-C_Px_FOU*T(i,j) - C_Py_FOU*T(i,j) - C_E_FOU*T(i+1,j) - C_W_FOU*T(i-1,j) - C_N_FOU*T(i,j+1) - C_S_FOU*T(i,j-1))
                    else
                        T_after(i,j) = T(i,j)&
                            & + del_t*(-C_Px_FOU*T(i,j) - C_Py_FOU*T(i,j) - C_E_FOU*T(i+1,j) - C_W_FOU*T(i-1,j) - C_N_FOU*T(i,j+1) - C_S_FOU*T(i,j-1)&
                            & - G*((C_Px_gen + C_Py_gen)*T(i,j) + C_E_gen*T(i+1,j) + C_W_gen*T(i-1,j) + C_N_gen*T(i,j+1) + C_S_gen*T(i,j-1)&
                            & + C_EE_HO*T(i+2,j) + C_WW_HO*T(i-2,j) + C_NN_HO*T(i,j+2) + C_SS_HO*T(i,j-2)) )
                    end if
                    
                end if
            end do
        end do
        T = T_after
        !if ( tt == 0 .OR. tt == 1 .OR. tt == 2 .OR. tt == 3 .OR. tt == 4) then
        !    call output(del_x,del_y,N,T,printout,tt)
        !    pause
        !end if
        !call animation(del_x,del_y,N,T,printout,tt)
        !pause
    end do
    call output(del_x,del_y,N,T,printout,tt)
end subroutine FTCS