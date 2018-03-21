subroutine test_full_orbit(dtao0)
  use precision,only:p_
  use constants,only:zero,kev,pi,twopi
  use normalizing,only: vn=>vn_i,tn=>tn_i
  use ions_module,only:mass=>mass_i,charge=>charge_i
  implicit none
  real(p_),intent(in):: dtao0
integer,parameter::ln=1.0_p_
  real(p_)::r,z,phi,vr,vz,vphi,t !,r0,z0,phi0
  integer:: kk
  integer,parameter:: maxstep=29000
  real(p_):: kin_eng,pphi
  real(p_):: b,psi_func
  real(p_):: br_SI, bz_SI,bphi_SI,b_SI !function names
  real(p_):: bval,brval,bzval,bphival
  real(p_):: v,energy,sai,pitch_angle
  real(p_):: rg,zg,phig,rg1,zg1,phig1
  real(p_):: b_dot_v,omega_local,dtao
  integer,parameter:: n_tor_period=1
  logical,parameter:: check_boundary_loss=.true.
  character(100), parameter::orbit_file="fo_go.txt"
  real(p_):: r0,z0,phi0,vr0,vphi0,vz0
  r= 2.1_p_
  z= 0._p_
  phi=0._p_
  vr=1.0d6
  vz=1.0d6
  vphi=5d5
!  dtao=1._p_

dtao=1.0*dtao0

  !--guiding-center orbit, to roughly verify the full orbit--
  bval=b_SI(r,z)
  brval=br_SI(r,z)
  bzval=bz_SI(r,z)
  bphival=bphi_SI(r,z)

!!$  call particle_to_guiding_center_location(r,phi,z,vr,vphi,vz,brval,bphival,bzval,rg,phig,zg)
!!$
!!$  v=sqrt(vr**2+vz**2+vphi**2)
!!$  energy=0.5_p_*mass*v*v/kev
!!$  write(*,*) 'kinetic energy (kev)=', energy
!!$
!!$  b_dot_v=brval*vr+bphival*vphi+bzval*vz
!!$  sai=b_dot_v/(bval*v)
!!$  pitch_angle=acos(sai)/pi*180._p_
!!$
!!$  call orbit(mass,charge,energy,pitch_angle,phig,rg,zg,dtao,n_tor_period,check_boundary_loss,orbit_file)

  !----
  !then calculate full orbit

  r=r/ln
  z=z/ln
  vr=vr/vn
  vz=vz/vn
  vphi=vphi/vn
!r=1.6182173678310701;z=-0.17532304528354634;phi=0.51458392767352690
!vr= 2.2153551525963856E-002; vz= -3.2974526832309299E-003; vphi=-3.1659760504767286E-002

!   omega_local=b_SI(r*ln,z*ln)*charge/mass
!   dtao=twopi/omega_local/tn/8_p_ !the time-step is chosen as in terms of the local gyro-period
t=0._p_
write(*,*) 'local gyro-period=',twopi/omega_local
 
     write(*,*) 'Using Boris algorithm to push full orbit'
     open(163,file='full_orbit_boris.txt')
     !call backward_half_step_cartesian(dtao,r,z,phi,vr,vz,vphi) !push only velocity, to set initial condition for the first step of boris algorithm
     r0=r;z0=z;phi0=phi; vr0=vr; vz0=vz; vphi0=vphi
     call forward_half_step_for_boris(dtao,r0,z0,phi0,vr0,vz0,vphi0,r,z,phi,vr,vz,vphi) !for testing

!     do kk=1,maxstep
     do kk=1,200
        call push_full_orbit_cylindrical_boris(dtao,r,phi,z,vr,vphi,vz)

        t=t+dtao
        kin_eng=0.5_p_*mass*(vr**2+vz**2+vphi**2)*vn**2/kev
        pphi=mass*r*ln*vphi*vn+charge*psi_func(r*ln,z*ln)
        !write(163,*) t, t*tn,r*ln,z*ln,phi,vr,vz,vphi,kin_eng, (vr**2+vz**2)/b(r,z),pphi   
        write(163,*) t+dtao/2, (t+dtao/2)*tn,r*ln,z*ln,phi,vr,vz,vphi,kin_eng, (vr**2+vz**2)/b(r,z),pphi !for the case of using forward initialization

     enddo
     write(*,*) 'kk=',kk
     close(163)

end subroutine test_full_orbit



!!$subroutine backward_half_step_cartesian(dtao,r,z,phi,vr,vz,vphi) !push only velocity, to set initial condition for the first step of boris algorithm, !using multi-steps, instead of one step, considering dtao used in Boris may be comparable to the gyro-period
!!$!actually working in Cartesian coordinates (i.e., constant basis vectors)
!!$  use precision,only:p_
!!$  use constants,only:zero,one,two,one_half,three,six,twopi,kev
!!$  use normalizing,only: vn=>vn_i
!!$  use ions_module,only:mass=>mass_i
!!$  implicit none
!!$  real(p_),intent(in):: dtao
!!$  real(p_),intent(inout):: r,z,phi,vr,vz,vphi  !instantaneous value of orbit
!!$  real(p_):: dvx,dvy,dvz
!!$  real(p_):: x_fo_dot,y_fo_dot,z_cartesian_fo_dot
!!$  real(p_):: vx_fo_dot,vy_fo_dot,vz_cartesian_fo_dot
!!$  real(p_):: kx1,ky1,kz1,kvx1,kvy1,kvz1 !Runge-Kutta steps
!!$!  real(p_):: kx2,ky2,kz2,kvx2,kvy2,kvz2 !Runge-Kutta steps
!!$  real(p_)::vx,vy,x,y,dx,dy,dz,z0,dt
!!$  real(p_):: step
!!$  integer,parameter::m=100 !if dtao is comparable or larger than the gyro-period, then the first backward half step needs to be finished with multiple rk steps
!!$  integer::k
!!$
!!$write(*,*)  "energy calculated in Cartesian ,before evolution=", 0.5_p_*mass*(vr**2+vz**2+vphi**2)*vn**2/kev
!!$  x=r
!!$  y=0._p_
!!$  z0=z
!!$  vx=vr
!!$  vy=vphi
!!$
!!$step=-0.5_p_*dtao
!!$dt=step/m
!!$do k=1,m
!!$  kx1=one_half*dt*x_fo_dot(x,y,z,vx,vy,vz)
!!$  ky1=one_half*dt*y_fo_dot(x,y,z,vx,vy,vz)
!!$  kz1=one_half*dt*z_cartesian_fo_dot(x,y,z,vx,vy,vz)
!!$  kvx1=   one_half*dt*vx_fo_dot(x,y,z,vx,vy,vz)
!!$  kvy1=   one_half*dt*vy_fo_dot(x,y,z,vx,vy,vz)
!!$  kvz1=   one_half*dt*vz_cartesian_fo_dot(x,y,z,vx,vy,vz)
!!$
!!$  dx=dt*x_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dy=dt*y_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dz=dt*z_cartesian_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dvx=   dt*vx_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dvy=   dt*vy_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dvz=   dt*vz_cartesian_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$
!!$  !update
!!$  x=x+dx
!!$  y=y+dy
!!$  z=z+dz
!!$  vx=vx+dvx
!!$  vy=vy+dvy
!!$  vz=vz+dvz
!!$enddo
!!$
!!$z=z0 !resume to the original value
!!$vr=vx
!!$vphi=vy
!!$
!!$write(*,*)  "energy calculated in Cartesian after evolution=", 0.5_p_*mass*(vr**2+vz**2+vphi**2)*vn**2/kev
!!$!write(*,*) 'dvphi=',dvphi
!!$end subroutine backward_half_step_cartesian
!!$
!!$
!!$function x_fo_dot(x,y,z,vx,vy,vz)
!!$  use precision,only:p_
!!$  use constants,only:twopi
!!$  implicit none
!!$  real(p_):: x_fo_dot,x,y,z,vx,vy,vz
!!$
!!$  x_fo_dot=vx
!!$end function 
!!$
!!$
!!$function y_fo_dot(x,y,z,vx,vy,vz)
!!$  use precision,only:p_
!!$  use constants,only:twopi
!!$  implicit none
!!$  real(p_):: y_fo_dot,x,y,z,vx,vy,vz
!!$
!!$  y_fo_dot=vy
!!$end function 
!!$
!!$
!!$function z_cartesian_fo_dot(x,y,z,vx,vy,vz)
!!$  use precision,only:p_
!!$  use constants,only:twopi
!!$  implicit none
!!$  real(p_):: z_cartesian_fo_dot,x,y,z,vx,vy,vz
!!$
!!$  z_cartesian_fo_dot=vz
!!$end function
!!$
!!$ 
!!$function vx_fo_dot(x,y,z,vx,vy,vz) 
!!$  use precision,only:p_
!!$  use constants,only:twopi
!!$  implicit none
!!$  real(p_):: vx_fo_dot,x,y,z,vx,vy,vz
!!$  real(p_):: bphi,bz,br !function names
!!$  real(p_):: brval,bphival,bzval,by
!!$  real(p_)::r,phi
!!$  r=sqrt(x*x+y*y)
!!$  phi=acos(x/r)
!!$  if(y<0) phi=-phi  !phi is assumed in the range [-pi,pi]
!!$  brval=br(r,z,phi)
!!$  bphival=bphi(r,z,phi)
!!$  bzval=bz(r,z,phi)
!!$  by=brval*sin(phi)+bphival*cos(phi)
!!$  vx_fo_dot=twopi*(-vz*by+vy*bzval) 
!!$end function
!!$
!!$
!!$function vy_fo_dot(x,y,z,vx,vy,vz) !without inertial force
!!$  use precision,only:p_
!!$  use constants,only:twopi
!!$  implicit none
!!$  real(p_):: vy_fo_dot,x,y,z,vx,vy,vz
!!$  real(p_):: br,bz,bphi  !function names
!!$  real(p_):: brval,bphival,bzval,bx
!!$  real(p_)::r,phi
!!$  r=sqrt(x*x+y*y)
!!$  phi=acos(x/r)
!!$  if(y<0) phi=-phi  !phi is assumed in the range [-pi,pi]
!!$  brval=br(r,z,phi)
!!$  bphival=bphi(r,z,phi)
!!$  bzval=bz(r,z,phi)
!!$  
!!$  bx=brval*cos(phi)-bphival*sin(phi)
!!$
!!$  vy_fo_dot=twopi*(-vx*bzval+vz*bx)
!!$
!!$end function
!!$
!!$
!!$function vz_cartesian_fo_dot(x,y,z,vx,vy,vz)
!!$  use precision,only:p_
!!$  use constants,only:twopi
!!$  implicit none
!!$  real(p_):: vz_cartesian_fo_dot,x,y,z,vx,vy,vz
!!$  real(p_):: bphi,br ,bz !function names
!!$  real(p_):: brval,bphival,bzval,bx,by
!!$  real(p_)::r,phi
!!$  r=sqrt(x*x+y*y)
!!$  phi=acos(x/r)
!!$  if(y<0) phi=-phi  !phi is assumed in the range [-pi,pi]
!!$  brval=br(r,z,phi)
!!$  bphival=bphi(r,z,phi)
!!$  bzval=bz(r,z,phi)
!!$
!!$  bx=brval*cos(phi)-bphival*sin(phi)
!!$  by=brval*sin(phi)+bphival*cos(phi)
!!$ vz_cartesian_fo_dot=twopi*(vx*by-vy*bx)
!!$
!!$end function
!!$
!!$
!!$
!!$
!!$
!!$subroutine forward_half_step_cartesian(dtao,r0,z0,phi0,vr0,vz0,vphi0,r1,z1,phi1,vr1,vz1,vphi1) 
!!$  use precision,only:p_
!!$  use constants,only:zero,one,two,one_half,three,six,twopi,kev
!!$  use normalizing,only: vn=>vn_i
!!$  use ions_module,only:mass=>mass_i
!!$  implicit none
!!$  real(p_),intent(in):: dtao
!!$  real(p_),intent(in):: r0,z0,phi0,vr0,vz0,vphi0
!!$  real(p_),intent(out):: r1,z1,phi1,vr1,vz1,vphi1
!!$  real(p_):: x_fo_dot,y_fo_dot,z_cartesian_fo_dot
!!$  real(p_):: vx_fo_dot,vy_fo_dot,vz_cartesian_fo_dot
!!$  real(p_):: kx1,ky1,kz1,kvx1,kvy1,kvz1 !Runge-Kutta steps
!!$!  real(p_):: kx2,ky2,kz2,kvx2,kvy2,kvz2 !Runge-Kutta steps
!!$  real(p_):: x,y,z,vx,vy,vz,dx,dy,dz,dvx,dvy,dvz !working variables
!!$  real(p_):: step,dt,alpha
!!$  integer,parameter::m=100 !if dtao is comparable or larger than the gyro-period, then the first backward half step needs to be finished with multiple rk steps
!!$  integer::k
!!$
!!$write(*,*)  "energy calculated in Cartesian ,before evolution=", 0.5_p_*mass*(vr0**2+vz0**2+vphi0**2)*vn**2/kev
!!$  x=r0
!!$  y=0._p_
!!$  z=z0
!!$  vx=vr0
!!$  vy=vphi0
!!$  vz=vz0
!!$
!!$step=0.5_p_*dtao
!!$dt=step/m
!!$do k=1,m
!!$  kx1=one_half*dt*x_fo_dot(x,y,z,vx,vy,vz)
!!$  ky1=one_half*dt*y_fo_dot(x,y,z,vx,vy,vz)
!!$  kz1=one_half*dt*z_cartesian_fo_dot(x,y,z,vx,vy,vz)
!!$  kvx1=   one_half*dt*vx_fo_dot(x,y,z,vx,vy,vz)
!!$  kvy1=   one_half*dt*vy_fo_dot(x,y,z,vx,vy,vz)
!!$  kvz1=   one_half*dt*vz_cartesian_fo_dot(x,y,z,vx,vy,vz)
!!$
!!$  dx=dt*x_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dy=dt*y_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dz=dt*z_cartesian_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dvx=   dt*vx_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dvy=   dt*vy_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dvz=   dt*vz_cartesian_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$
!!$  !update
!!$  x=x+dx
!!$  y=y+dy
!!$  z=z+dz
!!$  vx=vx+dvx
!!$  vy=vy+dvy
!!$  vz=vz+dvz
!!$enddo
!!$
!!$ r1=sqrt(x*x+y*y)
!!$
!!$  alpha=asin(y/r1)
!!$
!!$  phi1=phi0+alpha
!!$z1=z
!!$
!!$
!!$  vr1=cos(alpha)*vr0+sin(alpha)*vphi0
!!$  vphi1=-sin(alpha)*vr0+cos(alpha)*vphi0
!!$vz1=vz0
!!$write(*,*)  "energy calculated in Cartesian after evolution=", 0.5_p_*mass*(vr1**2+vz1**2+vphi1**2)*vn**2/kev
!!$!write(*,*) 'dvphi=',dvphi
!!$end subroutine forward_half_step_cartesian
