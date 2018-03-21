subroutine push_electron_orbit_weight_half_step(dtao) !Euler step (to estimate the values at the half-time-step)
  use precision,only:p_
  use constants,only: one,two,twopi,one_half,kev
  use normalizing,only: vn_e,vn_i
  use array_in_magnetic_coordinates_for_guiding_center_pusher,only:w1,w2,w3,w4,w5,w6,w7,w8,w9,w10 !as input
  use array_in_mc,only: b_mc_matrix,grad_psi_matrix, grad_alpha_matrix,&
       & grad_psi_dot_grad_alpha_matrix !as input
  use magnetic_coordinates,only: mpoloidal,nflux,theta_1d_array,radcor_1d_array,radcor_low1,radcor_upp1
  use electrons_module,only: nmarker_e,mass_e,te0,kappa_ne,kappa_te !as input
  use electrons_module,only: radcor_e,theta_e,alpha_e,vpar_e,mu_e,ps_vol_e,w_e !input
  use electrons_module,only: touch_bdry_e,active_e !both input and output
  use electrons_module,only: radcor_e_mid,theta_e_mid,alpha_e_mid,vpar_e_mid,w_e_mid !output
  use pputil !sorting markers among different processors
  use interpolate_module
  implicit none
  real(p_),intent(in)::  dtao
  real(p_):: radial_drift,theta_drift,alpha_drift,mirror_force,factor1,factor2,half_step
  real(p_):: w1val,w2val,w3val,w4val,w5val,w6val,w7val,w8val,w9val,w10val
  real(p_)::grad_psi_val,grad_alpha_val,grad_psi_dot_grad_alpha_val
  real(p_)::ef_par_val,ef_x_val,ef_y_val,mf_par_val,mf_x_val,mf_y_val
  real(p_):: b_val,v_val,eq_particle_number,tmp,tmp2,weprime
  real(p_):: fe0,gs_psi_prime !function names
  integer:: i,np_old,np_new,ierr

  half_step=dtao*one_half
  do i=1,nmarker_e
     if( touch_bdry_e(i).eqv..true.) cycle
     !the equilibrium field is 2D
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w1,theta_e(i),radcor_e(i),w1val) 
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w2,theta_e(i),radcor_e(i),w2val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w3,theta_e(i),radcor_e(i),w3val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w4,theta_e(i),radcor_e(i),w4val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w5,theta_e(i),radcor_e(i),w5val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w6,theta_e(i),radcor_e(i),w6val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w7,theta_e(i),radcor_e(i),w7val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w8,theta_e(i),radcor_e(i),w8val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w9,theta_e(i),radcor_e(i),w9val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w10,theta_e(i),radcor_e(i),w10val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_psi_matrix,&
          & theta_e(i),radcor_e(i),grad_psi_val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_alpha_matrix,&
          & theta_e(i),radcor_e(i),grad_alpha_val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_psi_dot_grad_alpha_matrix,&
          & theta_e(i),radcor_e(i),grad_psi_dot_grad_alpha_val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,b_mc_matrix,theta_e(i),radcor_e(i),b_val)

     call field_perturbation_on_marker(radcor_e(i),theta_e(i),alpha_e(i),active_e(i),&
          & ef_par_val,ef_x_val,ef_y_val,mf_par_val,mf_x_val,mf_y_val)

     factor1=one+vpar_e(i)/twopi*w1val
     radial_drift=vpar_e(i)**2/twopi*w3val/factor1 +mu_e(i)/(twopi*factor1)*w6val
     theta_drift=vpar_e(i)*(w2val+vpar_e(i)/twopi*w4val)/factor1 + mu_e(i)/(twopi*factor1)*w7val

     alpha_drift=vpar_e(i)**2/(twopi*factor1)*w5val+mu_e(i)/(twopi*factor1)*w8val
     mirror_force=-mu_e(i)/factor1*(w9val+vpar_e(i)/twopi*w10val)

     v_val=sqrt(vpar_e(i)**2+mu_e(i)*two*b_val)
     eq_particle_number=ps_vol_e(i)*fe0(v_val)
!write(*,*) 'eq_particle_number=',eq_particle_number
     tmp=te0*kev/mass_e/vn_e**2
     tmp2=kappa_ne+(v_val**2/(two*tmp)-1.5_p_)*kappa_te
     factor2=gs_psi_prime(radcor_e(i))*(grad_psi_val**2*grad_alpha_val**2-grad_psi_dot_grad_alpha_val**2)&
          & /(grad_psi_val*b_val**2)
 

    weprime=-(ef_y_val*factor2*vn_i/vn_e)*tmp2*eq_particle_number &
          & -(mf_x_val*grad_psi_val*vpar_e(i)/b_val &
          & +mf_y_val*grad_psi_dot_grad_alpha_val*vpar_e(i)/(grad_psi_val*b_val))&
          & *tmp2*eq_particle_number &
         & -vn_i/vn_e*twopi/tmp*(vpar_e(i)*ef_par_val+ef_x_val*radial_drift+ef_y_val*alpha_drift)*eq_particle_number !&
!         & -vn_i/vn_e*twopi/tmp*(ef_x_val*radial_drift+ef_y_val*alpha_drift)*eq_particle_number

!weprime=0._p_ !for testing

     radcor_e_mid(i)=radcor_e(i)+radial_drift*half_step
     theta_e_mid(i) =theta_e(i) + theta_drift*half_step
     alpha_e_mid(i) =alpha_e(i)  + alpha_drift*half_step
     vpar_e_mid(i)  =vpar_e(i)+   mirror_force*half_step
     w_e_mid(i)     =w_e(i)+           weprime*half_step

!     theta_e_mid(i)=theta_e_mid(i)-int(theta_e_mid(i)/twopi)*twopi !shift into the range [0:twopi]
!     if(theta_e_mid(i).lt.0) theta_e_mid(i)=theta_e_mid(i)+twopi !shift into the range [0:twopi]
!     call shift_to_zero_twopi_range(theta_e_mid(i))
     if(theta_e_mid(i)>twopi .or. theta_e_mid(i)<0._p_) &
          & call shift_particle_theta_then_alpha(radcor_e_mid(i),theta_e_mid(i),alpha_e_mid(i))

!     alpha_e_mid(i)=alpha_e_mid(i)-int(alpha_e_mid(i)/twopi)*twopi !shift  into the range [0:twopi]
!     if(alpha_e_mid(i).lt.0) alpha_e_mid(i)=alpha_e_mid(i)+twopi !shift into the range [0:twopi]
!    call shift_to_zero_twopi_range(alpha_e_mid(i))
        call shift_to_specified_toroidal_range(alpha_e_mid(i))
!     if(radcor_e_mid(i).ge.radcor_upp1 .or. radcor_e_mid(i).le.radcor_low1) touch_bdry_e(i)=.true.
      call check_whether_marker_in_boundary(radcor_e_mid(i),touch_bdry_e(i),active_e(i))
!if(abs(w_e_mid(i))>0.0000001) write(*,*) 'w_e(i) is nonzero', w_e_mid(i),i, mf_x_val,mf_y_val
  enddo


  !sorting particles according to thier theta and then assign them to the corresponding processors, using the subroutines provided in pputil_yj.f90
  np_old=nmarker_e
  call init_pmove(theta_e_mid(:),np_old,twopi,ierr)
  call pmove(theta_e_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(radcor_e_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(alpha_e_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vpar_e_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(w_e_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit 

  call pmove(theta_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(radcor_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(alpha_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vpar_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(mu_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(w_e(:),  np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit 
  call pmove(ps_vol_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove2(touch_bdry_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove2(active_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call end_pmove(ierr)
  nmarker_e=np_new

end subroutine push_electron_orbit_weight_half_step


subroutine push_electron_orbit_weight_full_step(dtao) !push electrons from t{n} to t{n+1} (using the values of field perturbations at the half-time-step)
  use precision,only:p_
  use constants,only: one,two,twopi,one_half,kev
  use normalizing,only: vn_e,vn_i
  use array_in_magnetic_coordinates_for_guiding_center_pusher,only:w1,w2,w3,w4,w5,w6,w7,w8,w9,w10 
  use array_in_mc,only: b_mc_matrix,grad_psi_matrix, grad_alpha_matrix,&
       & grad_psi_dot_grad_alpha_matrix 
  use magnetic_coordinates,only: mpoloidal,nflux,theta_1d_array,radcor_1d_array,radcor_low1,radcor_upp1
  use electrons_module,only: nmarker_e,mu_e,mass_e,te0,kappa_ne,kappa_te,ps_vol_e 
  use electrons_module,only: radcor_e_mid,theta_e_mid,alpha_e_mid,vpar_e_mid 
  use electrons_module,only: radcor_e,theta_e,alpha_e,vpar_e,w_e,touch_bdry_e,active_e !both input and output
  use pputil !sorting markers among different processors
  use interpolate_module
  implicit none
  real(p_),intent(in)::  dtao
  real(p_):: radial_drift,theta_drift,alpha_drift,mirror_force,factor1,factor2
  real(p_):: w1val,w2val,w3val,w4val,w5val,w6val,w7val,w8val,w9val,w10val
  real(p_)::grad_psi_val,grad_alpha_val,grad_psi_dot_grad_alpha_val
  real(p_)::ef_par_val,ef_x_val,ef_y_val,mf_par_val,mf_x_val,mf_y_val
  real(p_):: b_val,v_val,eq_particle_number,tmp,tmp2,weprime
  real(p_):: fe0,gs_psi_prime !function names
  integer:: i,np_old,np_new,ierr

  do i=1,nmarker_e
     if( touch_bdry_e(i).eqv..true.) cycle
     !the equilibrium field is 2D
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w1,theta_e_mid(i),radcor_e_mid(i),w1val) 
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w2,theta_e_mid(i),radcor_e_mid(i),w2val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w3,theta_e_mid(i),radcor_e_mid(i),w3val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w4,theta_e_mid(i),radcor_e_mid(i),w4val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w5,theta_e_mid(i),radcor_e_mid(i),w5val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w6,theta_e_mid(i),radcor_e_mid(i),w6val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w7,theta_e_mid(i),radcor_e_mid(i),w7val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w8,theta_e_mid(i),radcor_e_mid(i),w8val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w9,theta_e_mid(i),radcor_e_mid(i),w9val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w10,theta_e_mid(i),radcor_e_mid(i),w10val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_psi_matrix,&
          & theta_e_mid(i),radcor_e_mid(i),grad_psi_val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_alpha_matrix,&
          & theta_e_mid(i),radcor_e_mid(i),grad_alpha_val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_psi_dot_grad_alpha_matrix,&
          & theta_e_mid(i),radcor_e_mid(i),grad_psi_dot_grad_alpha_val)
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,b_mc_matrix,theta_e_mid(i),radcor_e_mid(i),b_val)

     call field_perturbation_on_marker(radcor_e_mid(i),theta_e_mid(i),alpha_e_mid(i),active_e(i),&
          & ef_par_val,ef_x_val,ef_y_val,mf_par_val,mf_x_val,mf_y_val)

     factor1=one+vpar_e_mid(i)/twopi*w1val
     radial_drift=vpar_e_mid(i)**2/twopi*w3val/factor1 +mu_e(i)/(twopi*factor1)*w6val
     theta_drift=vpar_e_mid(i)*(w2val+vpar_e_mid(i)/twopi*w4val)/factor1 + mu_e(i)/(twopi*factor1)*w7val

     alpha_drift=vpar_e_mid(i)**2/(twopi*factor1)*w5val+mu_e(i)/(twopi*factor1)*w8val
     mirror_force=-mu_e(i)/factor1*(w9val+vpar_e_mid(i)/twopi*w10val)

     v_val=sqrt(vpar_e_mid(i)**2+mu_e(i)*two*b_val)
     eq_particle_number=ps_vol_e(i)*fe0(v_val)
     tmp=te0*kev/(mass_e*vn_e**2)
     tmp2=kappa_ne+(v_val**2/(two*tmp)-1.5_p_)*kappa_te
   factor2=gs_psi_prime(radcor_e(i))*(grad_psi_val**2*grad_alpha_val**2-grad_psi_dot_grad_alpha_val**2)&
          & /(grad_psi_val*b_val**2)


    weprime=-(ef_y_val*factor2*vn_i/vn_e)*tmp2*eq_particle_number&
          & -(mf_x_val*grad_psi_val*vpar_e_mid(i)/b_val&
          & +mf_y_val*grad_psi_dot_grad_alpha_val*vpar_e_mid(i)/(grad_psi_val*b_val))&
          & *tmp2*eq_particle_number &
          & -vn_i/vn_e*twopi/tmp*(vpar_e_mid(i)*ef_par_val+ef_x_val*radial_drift+ef_y_val*alpha_drift)*eq_particle_number !&
!          & -vn_i/vn_e*twopi/tmp*(ef_x_val*radial_drift+ef_y_val*alpha_drift)*eq_particle_number
!weprime=0._p_ !for testing
     radcor_e(i)=radcor_e(i)+radial_drift*dtao
     theta_e(i) =theta_e(i) + theta_drift*dtao
     alpha_e(i)=alpha_e(i)  + alpha_drift*dtao
     vpar_e(i) =vpar_e(i)+   mirror_force*dtao
     w_e(i)=w_e(i)+               weprime*dtao

!!$     theta_e(i)=theta_e(i)-int(theta_e(i)/twopi)*twopi !shift into the range [0:twopi]
!!$     if(theta_e(i).lt.0) theta_e(i)=theta_e(i)+twopi !shift into the range [0:twopi]
!    call shift_to_zero_twopi_range(theta_e(i))
 if(theta_e(i)>twopi .or. theta_e(i)<0._p_) &
& call shift_particle_theta_then_alpha(radcor_e(i),theta_e(i),alpha_e(i))
!!$     alpha_e(i)=alpha_e(i)-int(alpha_e(i)/twopi)*twopi !shift  into the range [0:twopi]
!!$     if(alpha_e(i).lt.0) alpha_e(i)=alpha_e(i)+twopi !shift into the range [0:twopi]
!    call shift_to_zero_twopi_range(alpha_e(i))
        call shift_to_specified_toroidal_range(alpha_e(i))
     !if(radcor_e(i).ge.radcor_upp1 .or. radcor_e(i).le.radcor_low1) touch_bdry_e(i)=.true.
      call check_whether_marker_in_boundary(radcor_e(i),touch_bdry_e(i),active_e(i))
  enddo

  !sorting particles according to thier theta and then assign them to the corresponding processors, using the subroutines provided in pputil_yj.f90
  np_old=nmarker_e
  call init_pmove(theta_e(:),np_old,twopi,ierr)
  call pmove(theta_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(radcor_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(alpha_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(mu_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vpar_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(w_e(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit 
  call pmove(ps_vol_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove2(touch_bdry_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove2(active_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call end_pmove(ierr)
  nmarker_e=np_new
end subroutine push_electron_orbit_weight_full_step


function fe0(v) result (z) ! v in unit of vn, fe0 in unit 1/(Ln**3*vn**3)
  use precision,only: p_
  use constants,only: two,twopi,kev
  use normalizing,only: vn_e,Ln
  use electrons_module,only: mass_e,te0,ne0
  implicit none
  real(p_):: v,z
  real(p_):: v_si
  v_si=v*vn_e
  z=ne0*sqrt((mass_e/(twopi*te0*kev))**3)*exp(-mass_e*v_si**2/(two*te0*kev))
  z=z*(vn_e**3*Ln**3)
end function fe0


subroutine shift_particle_theta_then_alpha(radcor,theta,alpha) !shifting theta and make sure that phi is not changed by accordingly shifting alpha, i.e., keep the particle in the same sptial location
  use precision,only:p_
  use constants,only: twopi
  implicit none
  real(p_),intent(in):: radcor
  real(p_),intent(inout):: theta,alpha
  real(p_):: q_func_pfn !function name
  integer:: ishift
!!$  a=a-int(a/twopi)*twopi !shift into the range [0:twopi]
!!$  if(a.lt.0) a=a+twopi !shift into the range [0:twopi]

  ishift=floor(theta/twopi)
  theta=theta-ishift*twopi
  alpha=alpha+ishift*twopi*q_func_pfn(radcor) !this shift is needed to make the cylindrical toroidal angle phi be the same as the value before the theta-shifting, i.e., make the particle lie at the same sptaial location when doing the theta-shifting
end subroutine shift_particle_theta_then_alpha
