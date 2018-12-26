!-------------------------------------------------------------------
subroutine uv(N,i,j,u_EI_p, u_EI_n, u_WI_p, u_WI_n, v_NI_p, v_NI_n, v_SI_p, v_SI_n,u_p,v_p,u,v)
    implicit NONE
    
    integer :: i, j, N
    real*8 :: u_EI_p,u_EI_n,u_WI_p,u_WI_n,v_NI_p,v_NI_n,v_SI_p,v_SI_n,u_p,v_p
    real*8, dimension(N,N) :: u,v
    
    if (i.ne.N) then
        u_EI_p = (  0.5d0*( u(i+1,j)+u(i,j) ) + abs( 0.5d0*( u(i+1,j)+u(i,j) ) )  )/2
        u_EI_n = (  0.5d0*( u(i+1,j)+u(i,j) ) - abs( 0.5d0*( u(i+1,j)+u(i,j) ) )  )/2
    end if
    if (i.ne.1) then
        u_WI_p = (  0.5d0*( u(i-1,j)+u(i,j) ) + abs( 0.5d0*( u(i-1,j)+u(i,j) ) )  )/2
        u_WI_n = (  0.5d0*( u(i-1,j)+u(i,j) ) - abs( 0.5d0*( u(i-1,j)+u(i,j) ) )  )/2
    end if
    if (j.ne.N) then
        v_NI_p = (  0.5d0*( v(i,j+1)+v(i,j) ) + abs( 0.5d0*( v(i,j+1)+v(i,j) ) )  )/2
        v_NI_n = (  0.5d0*( v(i,j+1)+v(i,j) ) - abs( 0.5d0*( v(i,j+1)+v(i,j) ) )  )/2
    end if
    if (j.ne.1) then
        v_SI_p = (  0.5d0*( v(i,j-1)+v(i,j) ) + abs( 0.5d0*( v(i,j-1)+v(i,j) ) )  )/2
        v_SI_n = (  0.5d0*( v(i,j-1)+v(i,j) ) - abs( 0.5d0*( v(i,j-1)+v(i,j) ) )  )/2
    end if
    u_p = u(i,j)
    v_p = v(i,j)
    
end subroutine uv