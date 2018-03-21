subroutine prepare_coefficients_of_field_equations()
  use precision,only: p_
  use constants,only:one_half,two,one,elementary_charge,mu0,twopi
  use normalizing,only:bn,ln,tn_i,vn_e
  use perturbation_field_matrix,only:coeff1_ex, coeff1_ey, coeff1_ex_yy, coeff1_ex_xy, coeff1_ey_xx, coeff1_ey_xy !output, coefficients of the 1nd perpendicular equation
  use perturbation_field_matrix,only:coeff2_ex, coeff2_ey, coeff2_ex_yy, coeff2_ex_xy, coeff2_ey_xx, coeff2_ey_xy !output, coefficients of the 2nd perpendicular equation
  use perturbation_field_matrix,only:coeff3_epar_xx, coeff3_epar_yy,coeff3_epar_xy,&
     & coeff3_epar_x,coeff3_epar_y,coeff3_epar !output, coefficients of the parallel field equation
  use perturbation_field_matrix,only:coeff1_ey_implicit, coeff2_ex_implicit !output
  use domain_decomposition,only: theta_interval,theta_start !the 1D domain decomposion along theta coordinates, theta_start is the starting poloidal location of the cell treated by No. myid proc
  use magnetic_coordinates,only: mpoloidal,nflux,theta_1d_array,radcor_1d_array
  use array_in_mc,only: b_mc_matrix,grad_psi_matrix, grad_alpha_matrix,grad_psi_dot_grad_alpha_matrix !as input
  use electrons_module,only:charge_e,mass_e,ne0
  use flux_tube_model,only: radcor_fixed
  use interpolate_module
  implicit none

  real(p_):: grad_psi_val,grad_alpha_val, grad_psi_dot_grad_alpha_val,bval
  real(p_):: theta0,constant,term,dzerox,dzeroy
  real(p_):: gs_psi_prime !function name

  constant=tn_i*bn/(mu0*elementary_charge*ne0*ln**2) !a factor of dtao is included in the field solver
!  theta0=theta_start+theta_interval/two !the value of theta at the center of the poloidal cell
  theta0=theta_start
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_psi_matrix,&
       & theta0,radcor_fixed,grad_psi_val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_alpha_matrix,&
       & theta0,radcor_fixed,grad_alpha_val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_psi_dot_grad_alpha_matrix,&
       & theta0,radcor_fixed,grad_psi_dot_grad_alpha_val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,b_mc_matrix,theta0,radcor_fixed,bval)
!  write(*,*) 'constant=',constant, 'grad_psi_val=',grad_psi_val, 'grad_alpha_val=',grad_alpha_val
!  write(*,*) 'grad_psi_dot_grad_alpha_val=',grad_psi_dot_grad_alpha_val

  term=grad_psi_dot_grad_alpha_val**2-grad_psi_val**2*grad_alpha_val**2
  coeff1_ey_implicit=-twopi*gs_psi_prime(radcor_fixed)*term
  coeff2_ex_implicit=-twopi*gs_psi_prime(radcor_fixed)*(-term)
!  coeff1_ey_implicit=0. !for testing
!  coeff2_ex_implicit=0. !for testing

  !for the 1st perpendicular equation, at integer time-step
  coeff1_ex=grad_psi_val**2
  coeff1_ey=grad_psi_dot_grad_alpha_val
!  coeff1_ex=0. !for testing
!  coeff1_ey=0. !for testing

  coeff1_ex_yy=constant*bval**2/gs_psi_prime(radcor_fixed)*grad_psi_dot_grad_alpha_val
  coeff1_ex_xy=constant*bval**2/gs_psi_prime(radcor_fixed)*grad_psi_val**2
  coeff1_ey_xx=-coeff1_ex_xy
  coeff1_ey_xy=-coeff1_ex_yy

!!$  coeff1_ex_yy=0._p_ !for testing
!!$  coeff1_ex_xy=0._p_ !for testing
!!$  coeff1_ey_xx=0._p_ !for testing
!!$  coeff1_ey_xy=0._p_ !for testing

  !for the 2nd perpendicular equation, at integer time-step
  coeff2_ex=grad_psi_dot_grad_alpha_val
  coeff2_ey=grad_alpha_val**2
!  coeff2_ex=0._p_ !for testing
!  coeff2_ey=0._p_ !for testing

  coeff2_ex_yy=constant*bval**2/gs_psi_prime(radcor_fixed)*grad_alpha_val**2
  coeff2_ex_xy=constant*bval**2/gs_psi_prime(radcor_fixed)*grad_psi_dot_grad_alpha_val
  coeff2_ey_xx=-coeff2_ex_xy
  coeff2_ey_xy=-coeff2_ex_yy

!!$  coeff2_ex_yy=0._p_ !for testing
!!$  coeff2_ex_xy=0._p_ !for testing
!!$  coeff2_ey_xx=0._p_ !for testing
!!$  coeff2_ey_xy=0._p_ !for testing


  !for the 3rd equation (parallel direction field equation)
  coeff3_epar_xx=-grad_psi_val**2
  coeff3_epar_xy=-two*grad_psi_dot_grad_alpha_val
  coeff3_epar_yy=-grad_alpha_val**2

!!$  coeff3_epar_xx=0._p_ !for testing
!!$  coeff3_epar_xy=0._p_ !for testing
!!$  coeff3_epar_yy=0._p_ !for testing


  coeff3_epar=charge_e**2*mu0*ln**2*ne0/mass_e
  term=-tn_i*charge_e*mu0*vn_e**2*ne0/bn*gs_psi_prime(radcor_fixed)/bval
  dzerox=0._p_ !d0_dot_grad_psi, set to zero for testing
  dzeroy=0._p_ !d0_dot_grad_alpha, set to zero for testing
  coeff3_epar_x=term*(grad_psi_dot_grad_alpha_val*dzerox-grad_psi_val**2*dzeroy)
  coeff3_epar_y=term*(grad_alpha_val**2*dzerox-grad_psi_dot_grad_alpha_val*dzeroy)

!  write(*,*)  ' coeff1_ex, coeff1_ey, coeff2_ex,coeff2_ey, coeff1_ey_implicit,  coeff2_ex_implicit,=',&
 !      &  coeff1_ex, coeff1_ey, coeff2_ex,coeff2_ey, coeff1_ey_implicit,  coeff2_ex_implicit
end subroutine prepare_coefficients_of_field_equations
