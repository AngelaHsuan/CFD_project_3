!---------------------------------------------------------
! Coefficient-explicit
subroutine Coeff_E(i,j,Pe,beta,u,v,N,C_w_FOU,C_px_FOU,C_py_FOU,C_e_FOU,C_s_FOU,C_n_FOU,C_px_gen,C_py_gen,C_w_gen,C_e_gen,C_s_gen,C_n_gen,C_ww_HO,C_ee_HO,C_nn_HO,C_ss_HO,del_x,del_y)
    implicit NONE
    
    integer :: N,i,j
    real*8 :: del_x,del_y,r,freq,v_max,Pe,beta,u_EI_p, u_EI_n, u_WI_p, u_WI_n, v_NI_p, v_NI_n, v_SI_p, v_SI_n,u_p,v_p
    real*8 :: C_w_FOU,C_px_FOU,C_py_FOU,C_e_FOU,C_s_FOU,C_n_FOU,C_px_gen,C_py_gen,C_w_gen,C_e_gen,C_s_gen,C_n_gen
    real*8 :: C_ww_HO,C_ee_HO,C_nn_HO,C_ss_HO,C_w_HO,C_e_HO,C_n_HO,C_s_HO,C_px_HO,C_py_HO
    real*8, dimension(N,N) :: u, v
    
    call uv(N,i,j,u_EI_p, u_EI_n, u_WI_p, u_WI_n, v_NI_p, v_NI_n, v_SI_p, v_SI_n,u_p,v_p,u,v)
    
    C_w_FOU = -1.d0/del_x*u_WI_p - 1.d0/Pe/del_x**2
    C_px_FOU = 1.d0/del_x*(u_EI_p - u_WI_n) + 1.d0/Pe*2.d0/del_x**2
    C_py_FOU = 1.d0/del_y*(v_NI_p - v_SI_n) + 1.d0/Pe*2.d0/del_y**2
    C_e_FOU = 1.d0/del_x*u_EI_n - 1.d0/Pe/del_x**2
    C_s_FOU = -1.d0/del_y*v_SI_p - 1.d0/Pe/del_y**2 
    C_n_FOU = 1.d0/del_y*v_NI_n - 1.d0/Pe/del_y**2
    
    C_w_HO = 1.d0/del_x*(-beta*u_EI_p - u_WI_p*(0.5d0 + 2*beta) - u_WI_n*(0.5d0 - beta)) - 1.d0/Pe/del_x**2
    C_px_HO = 1.d0/del_x*(u_EI_p*(0.5d0 + 2*beta) + u_EI_n*(0.5d0 - beta) - u_WI_p*(0.5d0 - beta) - u_WI_n*(0.5d0 + 2*beta))&
            & + 1.d0/Pe*2.d0/del_x**2
    C_py_HO = 1.d0/del_y*(v_NI_p*(0.5d0 + 2*beta) + v_NI_n*(0.5d0 - beta) - v_SI_p*(0.5d0 - beta) - v_SI_n*(0.5d0 + 2*beta))&
            & + 1.d0/Pe*2.d0/del_y**2
    C_e_HO = 1.d0/del_x*(beta*u_WI_n + u_EI_p*(0.5d0 - beta) + u_EI_n*(0.5d0 + 2*beta)) - 1.d0/Pe/del_x**2
    C_s_HO = 1.d0/del_y*(-beta*v_NI_p - v_SI_p*(0.5d0 + 2*beta) - v_SI_n*(0.5d0 - beta)) - 1.d0/Pe/del_y**2
    C_n_HO = 1.d0/del_y*(beta*v_SI_n + v_NI_p*(0.5d0 - beta) + v_NI_n*(0.5d0 + 2*beta)) - 1.d0/Pe/del_y**2
    
    C_ww_HO = 1.d0/del_x*(beta*u_WI_p)
    C_ee_HO = -1.d0/del_x*(beta*u_EI_n)
    C_ss_HO = 1.d0/del_y*(beta*v_SI_p)
    C_nn_HO = -1.d0/del_y*(beta*v_NI_n)
    
    C_px_gen = C_px_HO - C_px_FOU
    C_py_gen = C_py_HO - C_py_FOU
    C_e_gen = C_e_HO - C_e_FOU
    C_w_gen = C_w_HO - C_w_FOU
    C_n_gen = C_n_HO - C_n_FOU
    C_s_gen = C_s_HO - C_s_FOU
    !write(*,*)C_px_gen,C_py_gen,C_w_gen,C_e_gen,C_s_gen,C_n_gen,C_ww_HO,C_ee_HO,C_nn_HO,C_ss_HO
    !pause
    
end subroutine Coeff_E