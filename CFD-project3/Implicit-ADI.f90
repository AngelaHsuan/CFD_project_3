!---------------------------------------------------------
! ADI
subroutine ADI(conv,w,del_t,del_x,del_y,G,T,N,v_max,Pe,beta,printout)
    implicit NONE
    
    integer :: i,j,k,nn,G,N,printout,N_iter
    real*8 :: del_t,del_x,del_y,v_max,Pe,beta,tt,conv,w
    real*8, dimension(N) :: C_w_FOU,C_px_FOU,C_py_FOU,C_e_FOU,C_s_FOU,C_n_FOU,C_px_gen,C_py_gen,C_w_gen,C_e_gen,C_s_gen,C_n_gen,C_ww_HO,C_ee_HO,C_nn_HO,C_ss_HO
    real*8, dimension(N) :: T_left,T_right,T_bottom,T_top,T_after_j,T_after_i,RHS
    real*8, dimension(N,N) :: T,T_after,T_after_new,T_err,u,v
    character :: circle
    
    T_after = T
    
    call velocity (N,u,v,del_x, del_y,v_max)
    
    do nn = 1,4.d0/del_t + 1,1
        tt = (nn-1)*del_t
        !write(*,'(1F6.4)')tt
        
        call Boundary(N,tt,del_t,del_x,del_y,T_left,T_right,T_bottom,T_top,v_max)
        
        T_err = 1.d0
        N_iter = 0
        !Row---------------------------------------------------------
        circle = 'r'
        do while (maxval(T_err) > conv)
            do j = 1,N,1
                if (j == 1) then
                    !Bottom
                    do i = 1,N,1
                        C_w_FOU(i) = 0
                        C_px_FOU(i) = 1
                        C_e_FOU(i) = 0
                        RHS(i) = T_bottom(i)
                        !T_after(i,j) = T_bottom(i)
                    end do
                else if (j == N) then
                    !Top
                    do i = 1,N,1
                        C_w_FOU(i) = 0
                        C_px_FOU(i) = 1
                        C_e_FOU(i) = 0
                        RHS(i) = T_top(i)
                        !T_after(i,j) = T_top(i)
                    end do
                else
                    do i = 1,N,1
                        if (i == 1) then
                            !left
                            C_w_FOU(i) = 0
                            C_px_FOU(i) = 1
                            C_e_FOU(i) = 0
                            RHS(i) = T_left(j)
                            !T_after(i,j) = T_left(j)
                        else if (i == N) then
                            !right
                            C_w_FOU(i) = 0
                            C_px_FOU(i) = 1
                            C_e_FOU(i) = 0
                            RHS(i) = T_right(j)
                            !T_after(i,j) = T_right(j)
                        else
                            !center
                            call Coeff_I(circle,i,j,Pe,beta,u,v,N,C_w_FOU,C_px_FOU,C_py_FOU,C_e_FOU,C_s_FOU,C_n_FOU,C_px_gen,C_py_gen,C_w_gen,C_e_gen,C_s_gen,C_n_gen,C_ww_HO,C_ee_HO,C_nn_HO,C_ss_HO,del_x,del_y)
                            C_px_FOU(i) = C_px_FOU(i)+(2.d0/del_t)
                            if (i == 2 .OR. i == N-1 .OR. j == 2 .OR. j == N-1) then
                                RHS(i) = (2.d0/del_t - C_py_FOU(i))*T(i,j) - C_n_FOU(i)*T(i,j+1) - C_s_FOU(i)*T(i,j-1)
                            else
                                RHS(i) = (2.d0/del_t - C_py_FOU(i))*T(i,j) - C_n_FOU(i)*T(i,j+1) - C_s_FOU(i)*T(i,j-1)&
                                & - G*( C_px_gen(i)*T_after(i,j) + C_py_gen(i)*T(i,j) + C_e_gen(i)*T_after(i+1,j) + C_w_gen(i)*T_after(i-1,j) + C_n_gen(i)*T(i,j+1) + C_s_gen(i)*T(i,j-1)&
                                & + C_ee_HO(i)*T_after(i+2,j) + C_ww_HO(i)*T_after(i-2,j) + C_nn_HO(i)*T(i,j+2) + C_ss_HO(i)*T(i,j-2) )
                            end if
                        end if
                    end do
                end if
                !if (j == 2) then
                    !write(*,*)C_e_FOU,'yo',j
                    !pause
                !end if
                call TDMA(N,C_w_FOU,C_px_FOU,C_e_FOU,RHS,T_after_j)
                
                !if (j == 2) then
                !    write(*,*)C_e_FOU,'yo',j
                !    pause
                !end if
                
                do i = 1,N,1
                    T_after_new(i,j) = T_after_j(i)
                    T_err(i,j) = abs(T_after_new(i,j)-T_after(i,j))
                    !T_after(i,j) = T_after_new(i,j)
                    T_after(i,j) = w*T_after_new(i,j) + (1-w)*T_after(i,j)
                end do 
            end do
            !write(*,'(A22,I7,A11,1E12.5,A13,1E12.5)') 'number of iteration :', N_iter,'error=', maxval(T_err)
            !N_iter = N_iter + 1
        end do
        !pause
        T = T_after
        !call output(del_x,del_y,N,T,printout,tt)
        !write(*,*)'here'
        !pause
        
        T_err = 1.d0
        
        circle = 'c'
        ! Column-----------------------------------------------------
        do while (maxval(T_err) > conv)
            do i = 1,N,1
                if (i == 1) then
                    !left
                    do j = 1,N,1
                        C_s_FOU(j) = 0
                        C_py_FOU(j) = 1
                        C_n_FOU(j) = 0
                        RHS(j) = T_left(j)
                        !T_after(i,j) = T_bottom(i)
                    end do
                else if (i == N) then
                    !Right
                    do j = 1,N,1
                        C_s_FOU(j) = 0
                        C_py_FOU(j) = 1
                        C_n_FOU(j) = 0
                        RHS(j) = T_right(j)
                        !T_after(i,j) = T_top(i)
                    end do
                else
                    do j = 1,N,1
                        if (j == 1) then
                            !bottom
                            C_s_FOU(j) = 0
                            C_py_FOU(j) = 1
                            C_n_FOU(j) = 0
                            RHS(j) = T_bottom(i)
                            !T_after(i,j) = T_bottom(j)
                        else if (i == N) then
                            !top
                            C_s_FOU(j) = 0
                            C_py_FOU(j) = 1
                            C_n_FOU(j) = 0
                            RHS(j) = T_top(i)
                            !T_after(i,j) = T_top(j)
                        else
                            !center
                            call Coeff_I(circle,i,j,Pe,beta,u,v,N,C_w_FOU,C_px_FOU,C_py_FOU,C_e_FOU,C_s_FOU,C_n_FOU,C_px_gen,C_py_gen,C_w_gen,C_e_gen,C_s_gen,C_n_gen,C_ww_HO,C_ee_HO,C_nn_HO,C_ss_HO,del_x,del_y)
                            C_py_FOU(j) = C_py_FOU(j)+(2.d0/del_t)
                            if (i == 2 .OR. i == N-1 .OR. j == 2 .OR. j == N-1) then
                                RHS(j) = (2.d0/del_t - C_px_FOU(j))*T(i,j) - C_e_FOU(j)*T(i+1,j) - C_w_FOU(j)*T(i-1,j)
                            else
                                RHS(j) = (2.d0/del_t - C_px_FOU(j))*T(i,j) - C_e_FOU(j)*T(i+1,j) - C_w_FOU(j)*T(i-1,j)&
                                & - G*( C_px_gen(j)*T(i,j) + C_py_gen(j)*T_after(i,j) + C_e_gen(j)*T(i+1,j) + C_w_gen(j)*T(i-1,j) + C_n_gen(j)*T_after(i,j+1) + C_s_gen(j)*T_after(i,j-1)&
                                & + C_ee_HO(j)*T(i+2,j) + C_ww_HO(j)*T(i-2,j) + C_nn_HO(j)*T_after(i,j+2) + C_ss_HO(j)*T_after(i,j-2) )
                            end if
                        end if
                    end do
                end if
                call TDMA(N,C_s_FOU,C_py_FOU,C_n_FOU,RHS,T_after_i)
                do j = 1,N,1
                    T_after_new(i,j) = T_after_i(j)
                    T_err(i,j) = abs(T_after_new(i,j)-T_after(i,j))
                    T_after(i,j) = w*T_after_new(i,j) + (1-w)*T_after(i,j)
                end do
            end do
            !write(*,'(A22,I7,A11,1E12.5,A13,1E12.5)') 'number of iteration :', N_iter,'error=', maxval(T_err)
            !N_iter = N_iter + 1
            
        end do
        
        T = T_after
        !if ( tt == 1 .OR. tt == 2 .OR. tt == 3 .OR. tt == 4) then
            !call output(del_x,del_y,N,T,printout,tt)
        !    pause
        !end if
        !call animation(del_x,del_y,N,T,printout,tt)
        !pause
    end do
    call output(del_x,del_y,N,T,printout,tt)
    
end subroutine ADI