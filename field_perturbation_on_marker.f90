subroutine field_perturbation_on_marker(radcor,theta,alpha,active,&
     & epar_val,ex_val,ey_val,mf_par_val,mf_x_val,mf_y_val)
  use precision,only:p_
  use constants,only: one,zero,twopi
  use perturbation_field_matrix,only: epar_left,ex_left,ey_left,mf_par_left,mf_x_left,mf_y_left
  use perturbation_field_matrix,only: epar_right,ex_right,ey_right,mf_par_right,mf_x_right,mf_y_right !field at nearby (larger) poloidal location
  use magnetic_coordinates,only: mtoroidal,nflux2,tor_1d_array,radcor_1d_array2,nsegment
  use domain_decomposition,only: theta_interval,theta_start
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  real(p_),intent(in)::radcor,theta,alpha
  logical,intent(in):: active
  real(p_),intent(out)::epar_val,ex_val,ey_val,mf_par_val,mf_x_val,mf_y_val
  real(p_):: coeff1,coeff2,tmp1,tmp2

  if(active.eqv..false.) then !force on particles outside the computational region is set to zero
     epar_val=0._p_
     ex_val=0._p_
     ey_val=0._p_
     mf_par_val=0._p_
     mf_x_val=0._p_
     mf_y_val=0._p_
  else
     coeff1=(theta-theta_start)/theta_interval
     coeff2=one-coeff1

     !  if(alpha>twopi .or. alpha<zero) write(*,*) 'field_on_marker,alpha=',alpha
     if(alpha>twopi/nsegment .or. alpha<zero) write(*,*) 'field_on_marker,alpha is out of range, alpha=',alpha
     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,epar_left,alpha,radcor,tmp1)  !uniform xarray and zarray are assumed
     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,epar_right,alpha,radcor,tmp2)  !uniform xarray and zarray are assumed
     epar_val=tmp1*coeff2+tmp2*coeff1

     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,ex_left,alpha,radcor,tmp1)  
     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,ex_right,alpha,radcor,tmp2)  
     ex_val=tmp1*coeff2+tmp2*coeff1

     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,ey_left,alpha,radcor,tmp1)  
     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,ey_right,alpha,radcor,tmp2)  
     ey_val=tmp1*coeff2+tmp2*coeff1

     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,mf_par_left,alpha,radcor,tmp1)
     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,mf_par_right,alpha,radcor,tmp2)
     mf_par_val=tmp1*coeff2+tmp2*coeff1

     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,mf_x_left,alpha,radcor,tmp1)  
     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,mf_x_right,alpha,radcor,tmp2)  
     mf_x_val=tmp1*coeff2+tmp2*coeff1

     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,mf_y_left,alpha,radcor,tmp1)
     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,mf_y_right,alpha,radcor,tmp2)
     mf_y_val=tmp1*coeff2+tmp2*coeff1
  endif
end subroutine field_perturbation_on_marker


subroutine field_perturbation_on_marker2(radcor,theta,alpha,&
     & ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val) !for adiabatic electron model
  use precision,only:p_
  use constants,only: one,zero,twopi
  use perturbation_field_matrix,only:ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left
  use perturbation_field_matrix,only: ef_cyl_r_right,ef_cyl_z_right,ef_cyl_phi_right !field at nearby (larger) poloidal location
  use magnetic_coordinates,only: mtoroidal,nflux2,tor_1d_array,radcor_1d_array2
  use domain_decomposition,only: theta_interval,theta_start !as input
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  real(p_),intent(in)::radcor,theta,alpha
  real(p_),intent(out)::ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val
  real(p_):: coeff1,coeff2,tmp1,tmp2

     coeff1=(theta-theta_start)/theta_interval
     coeff2=one-coeff1

     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_r_left,alpha,radcor,tmp1)  !uniform xarray and zarray are assumed
     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_r_right,alpha,radcor,tmp2)  !uniform xarray and zarray are assumed
     ef_cyl_r_val=tmp1*coeff2+tmp2*coeff1

     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_z_left,alpha,radcor,tmp1)  
     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_z_right,alpha,radcor,tmp2)  
     ef_cyl_z_val=tmp1*coeff2+tmp2*coeff1

     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_phi_left,alpha,radcor,tmp1)
     call linear_2d_interpolation(mtoroidal+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_phi_right,alpha,radcor,tmp2)  
     ef_cyl_phi_val=tmp1*coeff2+tmp2*coeff1
end subroutine field_perturbation_on_marker2
