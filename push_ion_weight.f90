subroutine push_ion_weight_for_adiabatic_electron_model(dtao,radcor_i,theta_i,alpha_i,&
     & r_i,z_i,phi_i,vr_i,vz_i,vphi_i,active_i,w_i,w_i_star,nmarker_i)
  use precision,only:p_
  use constants,only: two,twopi,one_half,kev
  use normalizing,only: vn_i,ln
  use ions_module,only:ps_vol_i,mass_i,kappa_ni,kappa_ti,ti0,v_i,ni0 ,normalizing_factor
  use array_in_mc,only: grad_psi_r_matrix,grad_psi_z_matrix
  use magnetic_coordinates,only:mpoloidal,nflux,theta_1d_array,radcor_1d_array
  use interpolate_module
  implicit none
  real(p_),intent(in):: dtao
  integer,intent(in):: nmarker_i
  real(p_),intent(in):: radcor_i(nmarker_i),theta_i(nmarker_i),alpha_i(nmarker_i),w_i(nmarker_i)
  real(p_),intent(in):: r_i(nmarker_i),z_i(nmarker_i),phi_i(nmarker_i)
  real(p_),intent(in):: vr_i(nmarker_i),vz_i(nmarker_i),vphi_i(nmarker_i)
  real(p_),intent(out):: w_i_star(nmarker_i)
  logical,intent(in):: active_i(nmarker_i)
  real(p_):: eq_particle_number,wiprime
  real(p_):: ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val
  real(p_):: tmp,tmp2,fi0
  real(p_):: minor_r_prime,gs_psi_prime
  real(p_):: br_val,bz_val,bphi_val,br_si,bz_si,bphi_si, grad_psi_r,grad_psi_z,bval
  integer:: k,ierr

  do k=1,nmarker_i
     v_i(k)=sqrt(vr_i(k)**2+vz_i(k)**2+vphi_i(k)**2)
  enddo
  tmp=ti0*kev/(mass_i*vn_i**2)

  eq_particle_number=ni0*sqrt((mass_i/(twopi*ti0*kev))**3)*ln**3*vn_i**3/normalizing_factor !explicitly cancel the exp(-v^2/vt^2) dependence, valid only for the case that both physical equilibrium distribution and marker distribution are Maxwellian.
  do k=1,nmarker_i
     if( active_i(k).eqv..false.) then
        w_i_star(k)=w_i(k) !weight keep unchanged, this is actually set to be zero in the routine dealing with the radial particle boundary condition
     else
        !eq_particle_number=ps_vol_i(k)*fi0(v_i(k)) !general case
        tmp2=kappa_ni+(v_i(k)**2/(two*tmp)-1.5_p_)*kappa_ti

        bz_val=bz_SI(r_i(k),z_i(k))
        br_val=br_SI(r_i(k),z_i(k))
        bphi_val=bphi_SI(r_i(k),z_i(k))
        bval=sqrt(bz_val*bz_val+br_val*br_val+bphi_val*bphi_val)

        call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
             & grad_psi_r_matrix,theta_i(k),radcor_i(k),grad_psi_r)  !in magnetic coordinates

        call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
             & grad_psi_z_matrix, theta_i(k),radcor_i(k),grad_psi_z)  !in magnetic coordinates

        call field_perturbation_on_marker2(radcor_i(k),theta_i(k),alpha_i(k),&
             & ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val)

        wiprime=eq_particle_number*twopi/tmp*(ef_cyl_r_val*vr_i(k)+ef_cyl_z_val*vz_i(k)+ef_cyl_phi_val*vphi_i(k)) &
             &+eq_particle_number*tmp2*minor_r_prime(radcor_i(k))/bval**2*&
             (grad_psi_z*ef_cyl_r_val*bphi_val-grad_psi_r*ef_cyl_z_val*bphi_val&
             & -grad_psi_z*ef_cyl_phi_val*br_val+grad_psi_r*ef_cyl_phi_val*bz_val)

        w_i_star(k)=w_i(k)+wiprime*dtao
     endif
  enddo
end subroutine push_ion_weight_for_adiabatic_electron_model


subroutine push_ion_weight0(dtao,radcor_i,theta_i,alpha_i,w_i,w_i_star,nmarker_i)
  use precision,only:p_
  use constants,only: two,twopi,one_half,kev
  use normalizing,only: vn_i
  use ions_module,only: vpar_i,vx_i,vy_i,v_i,active_i
  use ions_module,only:grad_psi_i,grad_alpha_i,grad_psi_dot_grad_alpha_i,bval_i
  use ions_module,only:ps_vol_i,mass_i,kappa_ni,kappa_ti,ti0
  !  use magnetic_coordinates,only:mpoloidal,nflux,theta_1d_array,radcor_1d_array
  !  use array_in_magnetic_coordinates_for_guiding_center_pusher,only: b_mc_matrix
  implicit none
  real(p_),intent(in):: dtao
  integer,intent(in):: nmarker_i
  real(p_),intent(in):: radcor_i(nmarker_i),theta_i(nmarker_i),alpha_i(nmarker_i),w_i(nmarker_i)
  real(p_),intent(out):: w_i_star(nmarker_i)

  real(p_):: eq_particle_number,wiprime
  real(p_):: ef_par_val,ef_x_val,ef_y_val,mf_par_val,mf_x_val,mf_y_val
  real(p_):: tmp,tmp2,fi0
  real(p_):: minor_r_prime,gs_psi_prime
  integer:: k

  tmp=ti0*kev/(mass_i*vn_i**2)

  !call ion_velocity_components_in_mc_at_half_time_step() 

  do k=1,nmarker_i
     if( active_i(k).eqv..false.) cycle
     eq_particle_number=ps_vol_i(k)*fi0(v_i(k))
     tmp2=kappa_ni+(v_i(k)**2/(two*tmp)-1.5_p_)*kappa_ti

     call field_perturbation_on_marker(radcor_i(k),theta_i(k),alpha_i(k),active_i(k),&
          & ef_par_val,ef_x_val,ef_y_val,mf_par_val,mf_x_val,mf_y_val)
     wiprime=eq_particle_number*twopi/tmp*(ef_par_val*vpar_i(k)+ef_x_val*vx_i(k)+ef_y_val*vy_i(k)) &
          &+eq_particle_number*tmp2*gs_psi_prime(radcor_i(k))*minor_r_prime(radcor_i(k))/bval_i(k)**2 &
          & *(grad_psi_i(k)**2*grad_alpha_i(k)**2-grad_psi_dot_grad_alpha_i(k)**2)*ef_y_val&
          &-eq_particle_number*tmp2*minor_r_prime(radcor_i(k))/bval_i(k)*(vx_i(k)*mf_par_val&
          & -vpar_i(k)*(mf_x_val*grad_psi_i(k)**2+mf_y_val*grad_psi_dot_grad_alpha_i(k)))
     !wiprime=eq_particle_number*twopi/tmp*(ef_par_val*vpar_i(k)+0._p_+0._p_) 
     w_i_star(k)=w_i(k)+wiprime*dtao
  enddo
end subroutine push_ion_weight0


subroutine push_ion_weight_without_eper(dtao,radcor_i,theta_i,alpha_i,w_i,w_i_star,nmarker_i)
  use precision,only:p_
  use constants,only: two,twopi,one_half,kev
  use normalizing,only: vn_i
  use ions_module,only: vpar_i,vx_i,vy_i,v_i,active_i
  use ions_module,only:grad_psi_i,grad_alpha_i,grad_psi_dot_grad_alpha_i,bval_i
  use ions_module,only:ps_vol_i,mass_i,kappa_ni,kappa_ti,ti0
  use flux_tube_model,only:radcor_fixed
!  use magnetic_coordinates,only:mpoloidal,nflux,theta_1d_array,radcor_1d_array
!  use array_in_magnetic_coordinates_for_guiding_center_pusher,only: b_mc_matrix
  implicit none
  real(p_),intent(in):: dtao
  integer,intent(in):: nmarker_i
  real(p_),intent(in):: radcor_i(nmarker_i),theta_i(nmarker_i),alpha_i(nmarker_i),w_i(nmarker_i)
  real(p_),intent(out):: w_i_star(nmarker_i)

  real(p_):: eq_particle_number,wiprime
  real(p_):: ef_par_val,ef_x_val,ef_y_val,mf_par_val,mf_x_val,mf_y_val
  real(p_):: tmp,tmp2,fi0
  real(p_):: minor_r_prime,gs_psi_prime
  integer:: k

  tmp=ti0*kev/(mass_i*vn_i**2)

  !call ion_velocity_components_in_mc_at_half_time_step() 

  do k=1,nmarker_i
    if( active_i(k).eqv..false.) cycle
     eq_particle_number=ps_vol_i(k)*fi0(v_i(k))
     tmp2=kappa_ni+(v_i(k)**2/(two*tmp)-1.5_p_)*kappa_ti

     call field_perturbation_on_marker(radcor_i(k),theta_i(k),alpha_i(k),active_i(k),&
          & ef_par_val,ef_x_val,ef_y_val,mf_par_val,mf_x_val,mf_y_val)

!     wiprime=eq_particle_number*twopi/tmp*(ef_par_val*vpar_i(k)+0._p_+0._p_) &
     wiprime=eq_particle_number*twopi/tmp*(0._p_+0._p_+0._p_) &
          &+eq_particle_number*tmp2*gs_psi_prime(radcor_fixed)*minor_r_prime(radcor_fixed)/bval_i(k)**2 &
          & *(grad_psi_i(k)**2*grad_alpha_i(k)**2-grad_psi_dot_grad_alpha_i(k)**2)*ef_y_val&
          &-eq_particle_number*tmp2*minor_r_prime(radcor_fixed)/bval_i(k)*(vx_i(k)*mf_par_val&
          & -vpar_i(k)*(mf_x_val*grad_psi_i(k)**2+mf_y_val*grad_psi_dot_grad_alpha_i(k)))
     w_i_star(k)=w_i(k)+wiprime*dtao
  enddo
end subroutine push_ion_weight_without_eper

subroutine update_ion_weight(dtao,w_i_star)      !using values of (r,v,E,B) at t_{n+1/2} to evaluate the rhs of the weight evolution equation
  use precision,only:p_
  use constants,only: two,twopi,one_half,kev
  use normalizing,only: vn_i
  use ions_module,only: nmarker_i, mass_i,kappa_ni,kappa_ti,ti0,ps_vol_i,active_i
!  use ions_module,only:radcor_i,theta_i,alpha_i
  use ions_module,only:radcor_i,theta_i,alpha_i
  use ions_module,only: vpar_i,vx_i,vy_i,v_i,grad_psi_i,grad_alpha_i,grad_psi_dot_grad_alpha_i,bval_i
!  use ions_module,only: w_i_star
  use ions_module,only: w_i !as input and output
!  use magnetic_coordinates,only:mpoloidal,nflux,theta_1d_array,radcor_1d_array
!  use array_in_magnetic_coordinates_for_guiding_center_pusher,only: b_mc_matrix
  implicit none
  real(p_),intent(in):: dtao
  real(p_),intent(in):: w_i_star(nmarker_i)
  real(p_):: eq_particle_number,wiprime
  real(p_):: ef_par_val,ef_x_val,ef_y_val,mf_par_val,mf_x_val,mf_y_val
  real(p_):: tmp,tmp2,fi0
  real(p_):: minor_r_prime,gs_psi_prime
  integer:: k

  tmp=ti0*kev/(mass_i*vn_i**2)

  !call ion_velocity_components_in_mc_at_half_time_step() 

  do k=1,nmarker_i
    if( active_i(k).eqv..false.) cycle
     eq_particle_number=ps_vol_i(k)*fi0(v_i(k))
     tmp2=kappa_ni+(v_i(k)**2/(two*tmp)-1.5_p_)*kappa_ti

     call field_perturbation_on_marker(radcor_i(k),theta_i(k),alpha_i(k),active_i(k),&
          & ef_par_val,ef_x_val,ef_y_val,mf_par_val,mf_x_val,mf_y_val)
!!$     wiprime=eq_particle_number*twopi/tmp*(ef_par_val*vpar_i(k)+ef_x_val*vx_i(k)+ef_y_val*vy_i(k)) &
!!$          &+eq_particle_number*tmp2*gs_psi_prime(radcor_i(k))*minor_r_prime(radcor_i(k))/bval_i(k)**2 &
!!$          & *(grad_psi_i(k)**2*grad_alpha_i(k)**2-grad_psi_dot_grad_alpha_i(k)**2)*ef_y_val&
!!$          &-eq_particle_number*tmp2*minor_r_prime(radcor_i(k))/bval_i(k)*(vx_i(k)*mf_par_val&
!!$          & -vpar_i(k)*(mf_x_val*grad_psi_i(k)**2+mf_y_val*grad_psi_dot_grad_alpha_i(k)))
! wiprime=eq_particle_number*twopi/tmp*(ef_x_val*vx_i(k)+ef_y_val*vy_i(k))
 wiprime=eq_particle_number*twopi/tmp*(ef_par_val*vpar_i(k)+ef_x_val*vx_i(k)+ef_y_val*vy_i(k))
!     wiprime=0._p_ !for testing
!     w_i(k)=w_i(k)+wiprime*dtao
w_i(k)=w_i_star(k) +wiprime*dtao
  enddo
end subroutine update_ion_weight



function fi0(v) result (z) ! v in unit of vn, fi0 in unit 1/(Ln**3*vn**3)
  use precision,only: p_
  use constants,only: two,twopi,kev
  use normalizing,only: vn_i,Ln
  use ions_module,only: mass_i,ti0,ni0
  implicit none
  real(p_):: v,z
  real(p_):: v_si
  v_si=v*vn_i
  z=ni0*sqrt((mass_i/(twopi*ti0*kev))**3)*exp(-mass_i*v_si**2/(two*ti0*kev))
  z=z*(vn_i**3*Ln**3)
end function




