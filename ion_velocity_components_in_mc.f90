subroutine ion_velocity_components_in_mc(r,phi,z,radcor,theta,vr,vphi,vz,&
     & v,vpar,vx,vy,grad_psi,grad_alpha,grad_psi_dot_grad_alpha,bval,tmp) !given velocity components in cylindrical coordinates, obtain the components in magnetic coordinates
  use precision,only: p_
  use constants,only: one
  use magnetic_coordinates,only: mpoloidal,nflux,radcor_1d_array,theta_1d_array,r_mag_surf,z_mag_surf
  use array_in_mc,only:grad_psi_matrix,grad_alpha_r_matrix,grad_alpha_z_matrix,grad_alpha_matrix,grad_psi_dot_grad_alpha_matrix
  use array_in_mc,only:grad_psi_r_matrix,grad_psi_z_matrix
  use flux_tube_model,only: radcor_fixed 
  use interpolate_module
  implicit none
  real(p_),intent(in):: r,phi,z,vr,vphi,vz !in cylindrical coordinates
  real(p_),intent(in):: radcor,theta
  real(p_),intent(out):: v,vpar,vx,vy !vx is defined by vx=v_dot_grad_x, vy is defined by vy=v_dot_grad_y,note that grad_x and grad_y are not perpendicular to each other
  real(p_),intent(out):: grad_psi,grad_alpha, grad_psi_dot_grad_alpha,bval,tmp

  real(p_)::br_si,bphi_si,bz_si !function names
  real(p_):: br,bphi,bz
  real(p_):: psi_r_func,psi_z_func
  real(p_):: psi_r,psi_z
  real(p_):: ddelta_dr_func,ddelta_dz_func !function names
  real(p_):: delta_r,delta_z
  real(p_):: grad_alpha_r_val,grad_alpha_z_val
  real(p_):: grad_psi_r,grad_psi_z
  real(p_):: vper_sq,rval,zval
  real(p_):: radcor0,theta0

  v=sqrt(vr*vr+vz*vz+vphi*vphi)

  radcor0=radcor_fixed !the radial location is fixed to a specified location (without fixing the radial location (especially for calculating vpar), there are numerical instabilities, this is the last revision I made before the code generates ion cyclontron waves, 2017.Oct.13), this must be changed when the code is extended to the global case
  !radcor0=radcor !global case
  theta0=theta

  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
       & r_mag_surf,theta0,radcor0,rval)  !in magnetic coordinates

  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
       & z_mag_surf,theta0,radcor0,zval)  !in magnetic coordinates

  br=br_si(rval,zval)
  bphi=bphi_si(rval,zval)
  bz=bz_si(rval,zval)
  bval=sqrt(br*br+bphi*bphi+bz*bz)

  vpar=(br*vr+bz*vz+bphi*vphi)/bval

  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
       & grad_psi_r_matrix,theta0,radcor0,grad_psi_r)  !in magnetic coordinates

  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
       & grad_psi_z_matrix,theta0,radcor0,grad_psi_z)  !in magnetic coordinates

  vx=(vr*grad_psi_r+vz*grad_psi_z) !/grad_psi

  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
       & grad_alpha_r_matrix,theta0,radcor0,grad_alpha_r_val)  !in magnetic coordinates
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
       & grad_alpha_z_matrix,theta0,radcor0,grad_alpha_z_val)  !in magnetic coordinates

  vy=vphi/r+vr*grad_alpha_r_val+vz*grad_alpha_z_val

  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
       & grad_psi_matrix,theta0,radcor0,grad_psi)  !in magnetic coordinates

  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
       & grad_alpha_matrix,theta0,radcor0,grad_alpha)  !in magnetic coordinates

  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
       & grad_psi_dot_grad_alpha_matrix,theta0,radcor0,grad_psi_dot_grad_alpha)  !in magnetic coordinates

  !for testing
  !vper_sq=(vx/grad_psi)**2+(vy/grad_alpha)**2 +2*(vx/grad_psi)*(vy/grad_alpha)*grad_psi_dot_grad_alpha/(grad_psi*grad_alpha) !wrong formula
  !tmp=sqrt(vpar**2+vper_sq)
end subroutine ion_velocity_components_in_mc

subroutine ion_velocity_components_in_mc_at_half_time_step() !wrapper subroutine
  use precision,only: p_
  use ions_module,only: nmarker_i,r_i_mid,phi_i_mid,z_i_mid,radcor_i_mid,theta_i_mid
  use ions_module,only: vr_i_mid,vphi_i_mid,vz_i_mid,active_i !in cylindrical coordinates, as input
  use ions_module,only: v_i,vpar_i,vx_i,vy_i, grad_psi_i,grad_alpha_i,grad_psi_dot_grad_alpha_i,bval_i !as output
  implicit none
  real(p_):: tmp(nmarker_i)
  integer:: k

  do k=1,nmarker_i
     if(active_i(k).eqv..true.)  call ion_velocity_components_in_mc(r_i_mid(k),phi_i_mid(k),z_i_mid(k),&
          & radcor_i_mid(k),theta_i_mid(k),vr_i_mid(k),vphi_i_mid(k),vz_i_mid(k),&
          & v_i(k),vpar_i(k),vx_i(k),vy_i(k),grad_psi_i(k),grad_alpha_i(k),grad_psi_dot_grad_alpha_i(k),bval_i(k),tmp(k))
  enddo

end subroutine


subroutine ion_velocity_components_in_mc_at_integer_time_step() !wrapper subroutine
  use precision,only: p_
  use ions_module,only: nmarker_i,r_i,phi_i,z_i,radcor_i,theta_i
  use ions_module,only:vr_i_integer,vphi_i_integer,vz_i_integer,active_i !in cylindrical coordinates, as input
  use ions_module,only: v_i,vpar_i,vx_i,vy_i, grad_psi_i,grad_alpha_i,grad_psi_dot_grad_alpha_i,bval_i !as output
  implicit none
  real(p_):: tmp(nmarker_i)
  integer:: k

  do k=1,nmarker_i
     if(active_i(k).eqv..true.) call ion_velocity_components_in_mc(r_i(k),phi_i(k),z_i(k),&
          & radcor_i(k),theta_i(k), vr_i_integer(k),vphi_i_integer(k),vz_i_integer(k),&
         & v_i(k),vpar_i(k),vx_i(k),vy_i(k), grad_psi_i(k),grad_alpha_i(k),grad_psi_dot_grad_alpha_i(k),bval_i(k),tmp(k)) 
  enddo
!for testing
!k=nmarker_i/4
!write(*,'(8(1pe14.5))') vx_i(k)/grad_psi_i(k),vy_i(k)/grad_alpha_i(k),vpar_i(k),v_i(k),tmp(k),r_i(k),phi_i(k),z_i(k)
end subroutine


