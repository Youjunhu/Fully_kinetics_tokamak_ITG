subroutine forward_half_step_for_boris(dtao,r0,z0,phi0,vr0,vz0,vphi0,r1,z1,phi1,vr1,vz1,vphi1) 
  !input: initial condition of the orbit: r0,z0,phi0,vr0,vz0,vphi0, at the same time t_{0}
  !Output: the instanteous value of (r,z,phi) after half_dtao, i.e., t_{1/2}, and the projection of velocity at t_{0} onto the basis vectors at t_{1/2}
  use precision,only:p_
  use constants,only:zero,one,two,one_half,twopi
  implicit none
  real(p_),intent(in):: dtao
  real(p_),intent(in):: r0,z0,phi0,vr0,vz0,vphi0  !instantaneous value of orbit at t0
  real(p_),intent(out):: r1,z1,phi1 !location at t0+half_dtao
  real(p_),intent(out):: vr1,vz1,vphi1   !projection of the old velocity (t=t0) on the new basis vector (determined by the new (r,z,phi))
  real(p_):: step,dt,dr,dz,dphi,dvr,dvz,dvphi
  real(p_):: r_fo_dot,z_fo_dot,phi_fo_dot,vr_fo_dot,vz_fo_dot,vphi_fo_dot !equations of motion
  real(p_):: kr1,kz1,kphi1,kvr1,kvz1,kvphi1 !Runge-Kutta steps
!  real(p_):: kr2,kz2,kphi2,kvr2,kvz2,kvphi2 !Runge-Kutta steps
  integer,parameter::m=100 ! if dtao is comparable or larger than the gyro-period, then the first backward half step needs to be finished with multiple steps
  integer:: k
  real(p_):: r,z,phi,vr,vz,vphi !working variables

  r=r0
  z=z0
  phi=phi0
  vr=vr0
  vz=vz0
  vphi=vphi0

  step=0.5_p_*dtao
  dt=step/m
  do k=1,m
     !2nd order Rung-Kuta method
     kr1=     one_half*dt*r_fo_dot(r,z,phi,vr,vz,vphi)
     kz1=     one_half*dt*z_fo_dot(r,z,phi,vr,vz,vphi)
     kphi1=   one_half*dt*phi_fo_dot(r,z,phi,vr,vz,vphi)
     kvr1=    one_half*dt*vr_fo_dot(r,z,phi,vr,vz,vphi)
     kvz1=    one_half*dt*vz_fo_dot(r,z,phi,vr,vz,vphi)
     kvphi1=  one_half*dt*vphi_fo_dot(r,z,phi,vr,vz,vphi)

     dr=    dt*r_fo_dot    (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dz=    dt*z_fo_dot    (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dphi=  dt*phi_fo_dot  (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dvr=   dt*vr_fo_dot   (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dvz=   dt*vz_fo_dot   (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dvphi= dt*vphi_fo_dot (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)

     !update
     r=r+      dr
     z=z+      dz
     phi=phi+  dphi
     vr=vr+dvr
     vz=vz+dvz
     vphi=vphi+dvphi
     !write(*,*) 'dvphi=',dvphi
  enddo

r1=r
z1=z
phi1=phi

dphi=phi1-phi0

vz1=vz0
vr1=  vphi0*sin(dphi)+vr0*cos(dphi) !projection of the old velocity on the new basis vector 
vphi1=vphi0*cos(dphi)-vr0*sin(dphi) !projection of the old velocity on the new basis vector 

!!$phi1=phi1-int(phi1/twopi)*twopi !shift into the range [0:twopi]
!!$if(phi1.lt.0) phi1=phi1+twopi !shift into the range [0:twopi]
!   call shift_to_zero_twopi_range(phi1)
   call shift_to_specified_toroidal_range(phi1)
end subroutine forward_half_step_for_boris



function r_fo_dot(r,z,phi,vr,vz,vphi)
  use precision,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: r_fo_dot,r,z,phi,vr,vz,vphi

  r_fo_dot=vr
end function r_fo_dot

function phi_fo_dot(r,z,phi,vr,vz,vphi)
  use precision,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: phi_fo_dot,r,z,phi,vr,vz,vphi

  phi_fo_dot=vphi/r
end function phi_fo_dot

function z_fo_dot(r,z,phi,vr,vz,vphi)
  use precision,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: z_fo_dot,r,z,phi,vr,vz,vphi

  z_fo_dot=vz
end function z_fo_dot

function vr_fo_dot(r,z,phi,vr,vz,vphi)
  use precision,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: vr_fo_dot,r,z,phi,vr,vz,vphi
  real(p_):: bphi,bz !function names
!  vr_fo_dot=twopi*(-vz*bphi(r,z,phi)+vphi*bz(r,z,phi)) !wrong, term vphi**2/r**2 is missing
!  vr_fo_dot=twopi*(-vz*bphi(r,z,phi)+vphi*bz(r,z,phi))+ vphi**2/r**2 !2017-Aug.12, a bug found, where r**2 should be replaced by r. kinetic energy is not well conserved by the buggy code (compared with results given by the Cartesian version of orbit integrator), which forced me to examine the codes to find possible bugs, and I finally found this bug. After correcting this code, the conservation of kinetic energy given by this code is as good as that of the Cartesian version.
  vr_fo_dot=twopi*(-vz*bphi(r,z,phi)+vphi*bz(r,z,phi))+ vphi**2/r !correct
end function vr_fo_dot



function vphi_fo_dot(r,z,phi,vr,vz,vphi) 
  use precision,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: vphi_fo_dot,r,z,phi,vr,vz,vphi
  real(p_):: br,bz  !function names

  vphi_fo_dot=twopi*(-vr*bz(r,z,phi)+vz*br(r,z,phi))-vphi*vr/r !the last term is inertial force

end function vphi_fo_dot

function vz_fo_dot(r,z,phi,vr,vz,vphi)
  use precision,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: vz_fo_dot,r,z,phi,vr,vz,vphi
  real(p_):: bphi,br !function names
  vz_fo_dot=twopi*(vr*bphi(r,z,phi)-vphi*br(r,z,phi))
end function vz_fo_dot
