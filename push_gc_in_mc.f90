module push_gc_in_mc_module
  private
  public push_gc_in_mc
contains
subroutine push_gc_in_mc(charge_sign,dtao,radcor,theta,alpha,vpar,mu,touch_bdry) !push a single guiding-center orbit in magnetic coordinates (field-line following coordinates) using 2nd Runge-Kutta
  use precision,only:p_
  use constants,only: one,twopi,one_half
  use array_in_magnetic_coordinates_for_guiding_center_pusher,only:w1,w2,w3,w4,w5,w6,w7,w8,w9,w10
  use magnetic_coordinates,only: mpoloidal,nflux,theta_1d_array,radcor_1d_array,radcor_low1,radcor_upp1
  use interpolate_module
  implicit none
  real(p_),intent(in)::charge_sign,dtao,mu
  real(p_),intent(inout)::  radcor,theta,alpha,vpar
  logical,intent(out):: touch_bdry
  real(p_):: radial_drift,theta_drift,alpha_drift,mirror_force,factor1
  real(p_):: w1_func,w2_func,w3_func,w4_func,w5_func,w6_func,w7_func,w8_func,w9_func,w10_func
  real(p_):: w1val,w2val,w3val,w4val,w5val,w6val,w7val,w8val,w9val,w10val
  real(p_):: radcor_mid,theta_mid,alpha_mid,vpar_mid
  integer:: i

  if(radcor.ge.radcor_upp1 .or. radcor.le.radcor_low1)  then !check whether particle location is within the specified region
     touch_bdry=.true.
     return
  endif

  call shift_to_minus_pi_positive_pi_range(theta)
!!$if(theta.gt.twopi .and. theta.lt.0._p_) then
!!$write(*,*) '*********theta=', theta
!!$stop
!!$endif

  !write(*,*) theta,radcor
  !2nd Runge-Kutta, first step
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w1,theta,radcor,w1val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w2,theta,radcor,w2val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w3,theta,radcor,w3val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w4,theta,radcor,w4val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w5,theta,radcor,w5val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w6,theta,radcor,w6val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w7,theta,radcor,w7val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w8,theta,radcor,w8val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w9,theta,radcor,w9val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w10,theta,radcor,w10val)
!!$w1val=w1_func(theta,radcor)
!!$w2val=w2_func(theta,radcor)
!!$w3val=w3_func(theta,radcor)
!!$w4val=w4_func(theta,radcor)
!!$w5val=w5_func(theta,radcor)
!!$w6val=w6_func(theta,radcor)
!!$w7val=w7_func(theta,radcor)
!!$w8val=w8_func(theta,radcor)
!!$w9val=w9_func(theta,radcor)
!!$w10val=w10_func(theta,radcor)

  factor1=one+charge_sign*vpar/twopi*w1val
  radial_drift=charge_sign*vpar**2/twopi*w3val/factor1 +charge_sign*mu/(twopi*factor1)*w6val
  theta_drift=vpar*(w2val+charge_sign*vpar/twopi*w4val)/factor1 + charge_sign*mu/(twopi*factor1)*w7val
  alpha_drift=charge_sign*vpar**2/(twopi*factor1)*w5val+charge_sign*mu/(twopi*factor1)*w8val
  mirror_force=-mu/factor1*(w9val+charge_sign*vpar/twopi*w10val)

  radcor_mid=radcor+radial_drift*dtao*one_half
  theta_mid =theta + theta_drift*dtao*one_half
  alpha_mid=alpha  + alpha_drift*dtao*one_half
  vpar_mid =vpar+   mirror_force*dtao*one_half

  !write(*,*) 'drift=',radial_drift,theta_drift,alpha_drift,mirror_force
  call shift_to_minus_pi_positive_pi_range(theta_mid)

  if(radcor_mid.ge.radcor_upp1 .or. radcor_mid.le.radcor_low1)  then
     touch_bdry=.true.
     return
  endif

!!$  if(theta_mid>twopi .or.theta_mid<0) then
!!$
!!$     write(*,*) 'theta=',theta,'theta_drift*dtao',theta_drift*dtao*one_half,'theta_mid=,radcor_mid=',&
!!$          & theta_mid,radcor_mid, 'warning****************'
!!$  endif


  !2nd Runge-Kutta, second step
  !write(*,*) 'theta_mid,radcor_mid=',theta_mid,radcor_mid
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w1,theta_mid,radcor_mid,w1val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w2,theta_mid,radcor_mid,w2val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w3,theta_mid,radcor_mid,w3val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w4,theta_mid,radcor_mid,w4val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w5,theta_mid,radcor_mid,w5val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w6,theta_mid,radcor_mid,w6val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w7,theta_mid,radcor_mid,w7val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w8,theta_mid,radcor_mid,w8val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w9,theta_mid,radcor_mid,w9val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w10,theta_mid,radcor_mid,w10val)

!!$w1val=w1_func(theta_mid,radcor_mid)
!!$w2val=w2_func(theta_mid,radcor_mid)
!!$w3val=w3_func(theta_mid,radcor_mid)
!!$w4val=w4_func(theta_mid,radcor_mid)
!!$w5val=w5_func(theta_mid,radcor_mid)
!!$w6val=w6_func(theta_mid,radcor_mid)
!!$w7val=w7_func(theta_mid,radcor_mid)
!!$w8val=w8_func(theta_mid,radcor_mid)
!!$w9val=w9_func(theta_mid,radcor_mid)
!!$w10val=w10_func(theta_mid,radcor_mid)

  factor1=one+charge_sign*vpar_mid/twopi*w1val
  radial_drift=charge_sign*vpar_mid**2/twopi*w3val/factor1 +charge_sign*mu/(twopi*factor1)*w6val
  theta_drift=vpar_mid*(w2val+charge_sign*vpar_mid/twopi*w4val)/factor1 +charge_sign*mu/(twopi*factor1)*w7val
  alpha_drift=charge_sign*vpar_mid**2/(twopi*factor1)*w5val+charge_sign*mu/(twopi*factor1)*w8val
  mirror_force=-mu/factor1*(w9val+charge_sign*vpar_mid/twopi*w10val)

  !write(*,*) 'mirror_foce=',mirror_force, w9val
  radcor=radcor+radial_drift*dtao
  theta =theta + theta_drift*dtao
  alpha=alpha  + alpha_drift*dtao
  vpar =vpar+mirror_force*dtao

end subroutine push_gc_in_mc
end module push_gc_in_mc_module
