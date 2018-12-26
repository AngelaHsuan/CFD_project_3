!---------------------------------------------------------
! Coefficient-implicit
subroutine Coeff_I(circle,i,j,Pe,beta,u,v,N,C_w_FOU,C_px_FOU,C_py_FOU,C_e_FOU,C_s_FOU,C_n_FOU,C_px_gen,C_py_gen,C_w_gen,C_e_gen,C_s_gen,C_n_gen,C_ww_HO,C_ee_HO,C_nn_HO,C_ss_HO,del_x,del_y)
    implicit NONE
    
    integer :: N,i,j,pad
    real*8 :: del_x,del_y,r,freq,v_max,Pe,beta,u_EI_p, u_EI_n, u_WI_p, u_WI_n, v_NI_p, v_NI_n, v_SI_p, v_SI_n,u_p,v_p
    real*8, dimension(N) :: C_w_FOU,C_px_FOU,C_py_FOU,C_e_FOU,C_s_FOU,C_n_FOU,C_px_gen,C_py_gen,C_w_gen,C_e_gen,C_s_gen,C_n_gen
    real*8, dimension(N) :: C_ww_HO,C_ee_HO,C_nn_HO,C_ss_HO,C_w_HO,C_e_HO,C_n_HO,C_s_HO,C_px_HO,C_py_HO
    real*8, dimension(N,N) :: u, v
    character :: circle
    
    call uv(N,i,j,u_EI_p, u_EI_n, u_WI_p, u_WI_n, v_NI_p, v_NI_n, v_SI_p, v_SI_n,u_p,v_p,u,v)
    
    if (circle == 'r') then
        C_w_HO(i) = 0
        C_px_HO(i) = 0
        C_py_HO(i) = 0
        C_e_HO(i) = 0
        C_s_HO(i) = 0
        C_n_HO(i) = 0
        C_ww_HO(i) = 0
        C_ee_HO(i) = 0
        C_ss_HO(i) = 0
        C_nn_HO(i) = 0
    
        C_w_FOU(i) = -1.d0/del_x*u_WI_p - 1.d0/Pe/del_x**2
        C_px_FOU(i) = 1.d0/del_x*(u_EI_p - u_WI_n) + 1.d0/Pe*2.d0/del_x**2
        C_py_FOU(i) = 1.d0/del_y*(v_NI_p - v_SI_n) + 1.d0/Pe*2.d0/del_y**2
        C_e_FOU(i) = 1.d0/del_x*u_EI_n - 1.d0/Pe/del_x**2
        C_s_FOU(i) = -1.d0/del_y*v_SI_p - 1.d0/Pe/del_y**2 
        C_n_FOU(i) = 1.d0/del_y*v_NI_n - 1.d0/Pe/del_y**2
    
        C_w_HO(i) = 1.d0/del_x*(-beta*u_EI_p - u_WI_p*(0.5d0 + 2*beta) - u_WI_n*(0.5d0 - beta)) - 1.d0/Pe/del_x**2
        C_px_HO(i) = 1.d0/del_x*(u_EI_p*(0.5d0 + 2*beta) + u_EI_n*(0.5d0 - beta) - u_WI_p*(0.5d0 - beta) - u_WI_n*(0.5d0 + 2*beta))&
                & + 1.d0/Pe*2.d0/del_x**2
        C_py_HO(i) = 1.d0/del_y*(v_NI_p*(0.5d0 + 2*beta) + v_NI_n*(0.5d0 - beta) - v_SI_p*(0.5d0 - beta) - v_SI_n*(0.5d0 + 2*beta))&
                & + 1.d0/Pe*2.d0/del_y**2
        C_e_HO(i) = 1.d0/del_x*(beta*u_WI_n + u_EI_p*(0.5d0 - beta) + u_EI_n*(0.5d0 + 2*beta)) - 1.d0/Pe/del_x**2
        C_s_HO(i) = 1.d0/del_y*(-beta*v_NI_p - v_SI_p*(0.5d0 + 2*beta) - v_SI_n*(0.5d0 - beta)) - 1.d0/Pe/del_y**2
        C_n_HO(i) = 1.d0/del_y*(beta*v_SI_n + v_NI_p*(0.5d0 - beta) + v_NI_n*(0.5d0 + 2*beta)) - 1.d0/Pe/del_y**2
        
        C_ww_HO(i) = 1.d0/del_x*(beta*u_WI_p)
        C_ee_HO(i) = -1.d0/del_x*(beta*u_EI_n)
        C_ss_HO(i) = 1.d0/del_y*(beta*v_SI_p)
        C_nn_HO(i) = -1.d0/del_y*(beta*v_NI_n)
    
        C_px_gen(i) = C_px_HO(i) - C_px_FOU(i)
        C_py_gen(i) = C_py_HO(i) - C_py_FOU(i)
        C_e_gen(i) = C_e_HO(i) - C_e_FOU(i)
        C_w_gen(i) = C_w_HO(i) - C_w_FOU(i)
        C_n_gen(i) = C_n_HO(i) - C_n_FOU(i)
        C_s_gen(i) = C_s_HO(i) - C_s_FOU(i)
        !if (j == 2 ) then
        !    write(*,*)C_e_FOU,u_EI_n,'yA'
        !    pause
        !end if
        !write(*,*)C_px_gen(i),C_py_gen(i),C_w_gen(i),C_e_gen(i),C_s_gen(i),C_n_gen(i),C_ww_HO(i),C_ee_HO(i),C_nn_HO(i),C_ss_HO(i)
        !pause
    else if (circle == 'c') then
        C_w_HO(j) = 0
        C_px_HO(j) = 0
        C_py_HO(j) = 0
        C_e_HO(j) = 0
        C_s_HO(j) = 0
        C_n_HO(j) = 0
        C_ww_HO(j) = 0
        C_ee_HO(j) = 0
        C_ss_HO(j) = 0
        C_nn_HO(j) = 0
    
        C_w_FOU(j) = -1.d0/del_x*u_WI_p - 1.d0/Pe/del_x**2
        C_px_FOU(j) = 1.d0/del_x*(u_EI_p - u_WI_n) + 1.d0/Pe*2.d0/del_x**2
        C_py_FOU(j) = 1.d0/del_y*(v_NI_p - v_SI_n) + 1.d0/Pe*2.d0/del_y**2
        C_e_FOU(j) = 1.d0/del_x*u_EI_n - 1.d0/Pe/del_x**2
        C_s_FOU(j) = -1.d0/del_y*v_SI_p - 1.d0/Pe/del_y**2
        C_n_FOU(j) = 1.d0/del_y*v_NI_n - 1.d0/Pe/del_y**2
    
        C_w_HO(j) = 1.d0/del_x*(-beta*u_EI_p - u_WI_p*(0.5d0 + 2*beta) - u_WI_n*(0.5d0 - beta)) - 1.d0/Pe/del_x**2
        C_px_HO(j) = 1.d0/del_x*(u_EI_p*(0.5d0 + 2*beta) + u_EI_n*(0.5d0 - beta) - u_WI_p*(0.5d0 - beta) - u_WI_n*(0.5d0 + 2*beta))&
                & + 1.d0/Pe*2.d0/del_x**2
        C_py_HO(j) = 1.d0/del_y*(v_NI_p*(0.5d0 + 2*beta) + v_NI_n*(0.5d0 - beta) - v_SI_p*(0.5d0 - beta) - v_SI_n*(0.5d0 + 2*beta))&
                & + 1.d0/Pe*2.d0/del_y**2
        C_e_HO(j) = 1.d0/del_x*(beta*u_WI_n + u_EI_p*(0.5d0 - beta) + u_EI_n*(0.5d0 + 2*beta)) - 1.d0/Pe/del_x**2
        C_s_HO(j) = 1.d0/del_y*(-beta*v_NI_p - v_SI_p*(0.5d0 + 2*beta) - v_SI_n*(0.5d0 - beta)) - 1.d0/Pe/del_y**2
        C_n_HO(j) = 1.d0/del_y*(beta*v_SI_n + v_NI_p*(0.5d0 - beta) + v_NI_n*(0.5d0 + 2*beta)) - 1.d0/Pe/del_y**2
    
        C_ww_HO(j) = 1.d0/del_x*(beta*u_WI_p)
        C_ee_HO(j) = -1.d0/del_x*(beta*u_EI_n)
        C_ss_HO(j) = 1.d0/del_y*(beta*v_SI_p)
        C_nn_HO(j) = -1.d0/del_y*(beta*v_NI_n)
    
        C_px_gen(j) = C_px_HO(j) - C_px_FOU(j)
        C_py_gen(j) = C_py_HO(j) - C_py_FOU(j)
        C_e_gen(j) = C_e_HO(j) - C_e_FOU(j)
        C_w_gen(j) = C_w_HO(j) - C_w_FOU(j)
        C_n_gen(j) = C_n_HO(j) - C_n_FOU(j)
        C_s_gen(j) = C_s_HO(j) - C_s_FOU(j)
        !write(*,*)C_px_gen(j),C_py_gen(j),C_w_gen(j),C_e_gen(j),C_s_gen(j),C_n_gen(j),C_ww_HO(j),C_ee_HO(j),C_nn_HO(j),C_ss_HO(j)
        !pause
    else
        write(*,*)'There is mistake in coefficient-Im'
        pause
    end if
    
end subroutine Coeff_I