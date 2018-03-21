subroutine push_gc_in_cylindrical(dtao,mu,r,z,phi,vpar)
  !input: initial condition of the orbit: mu,r,z,phi,vpar
  !Output: the instanteous value of the orbit after dtao: mu,r,z,phi,vpar
  use precision,only:p_
  use constants,only:zero,one,two,one_half,three,six,twopi
  implicit none
  real(p_),intent(in):: dtao,mu
  real(p_),intent(inout):: r,z,phi,vpar  !instantaneous value of orbit
  real(p_):: dr,dz,dphi,dvpar
  real(p_):: r_dot,z_dot,phi_dot,vpar_dot !equations of motion
  real(p_):: kr1,kz1,kvpar1,kphi1,kr2,kz2,kvpar2,kphi2,kr3,kz3,kvpar3,kphi3,kr4,kz4,kvpar4,kphi4 !Runge-Kutta steps

  !4nd order Rung-Kutta method
  kr1=dtao*r_dot(r,z,phi,vpar,mu)
  kz1=dtao*z_dot(r,z,phi,vpar,mu)
  kvpar1=dtao*vpar_dot(r,z,phi,vpar,mu)
  kphi1=dtao*phi_dot(r,z,phi,vpar,mu)

  kr2=dtao*r_dot      (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu)
  kz2=dtao*z_dot      (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu)
  kvpar2=dtao*vpar_dot(r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu)
  kphi2=dtao*phi_dot  (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu)

  kr3=dtao*r_dot      (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu)
  kz3=dtao*z_dot      (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu)
  kvpar3=dtao*vpar_dot(r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu)
  kphi3=dtao*phi_dot  (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu)

  kr4=dtao*r_dot      (r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu)
  kz4=dtao*z_dot      (r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu)
  kvpar4=dtao*vpar_dot(r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu)
  kphi4=dtao*phi_dot  (r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu)

  dr=kr1/six+kr2/three+kr3/three+kr4/six
  dz=kz1/six+kz2/three+kz3/three+kz4/six
  dphi=kphi1/six+kphi2/three+kphi3/three+kphi4/six
  dvpar=kvpar1/six+kvpar2/three+kvpar3/three+kvpar4/six
  !update
  r=r+      dr
  z=z+      dz
  vpar=vpar+dvpar
  phi=phi+  dphi
  !write(*,*) 'dvpar=',dvpar
end subroutine push_gc_in_cylindrical


function r_dot(r,z,phi,vpar,mu) result(funcval)
  use precision,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: funcval,r,z,phi,vpar,mu
  real(p_):: bstar_r,bstar_parallel,bstar_parallelval
  real(p_):: b,bphi,b_z,bz,b_phi

  bstar_parallelval=bstar_parallel(r,z,phi,vpar)
  funcval=bstar_r(r,z,phi,vpar)/bstar_parallelval*vpar +&
       & mu/(twopi*b(r,z,phi)*bstar_parallelval)*&
       & (bphi(r,z,phi)*b_z(r,z,phi)-bz(r,z,phi)*b_phi(r,z,phi)/r)

end function r_dot


function z_dot(r,z,phi,vpar,mu) result(funcval)
  use precision,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: funcval,r,z,phi,vpar,mu
  real(p_):: bstar_z,bstar_parallel,bstar_parallelval
  real(p_)::b,bphi,b_r,br,b_phi

  bstar_parallelval=bstar_parallel(r,z,phi,vpar)
  funcval= bstar_z(r,z,phi,vpar)/bstar_parallelval*vpar + &
       & mu/(twopi*b(r,z,phi)*bstar_parallelval)*(br(r,z,phi)*b_phi(r,z,phi)/r-bphi(r,z,phi)*b_r(r,z,phi))
end function z_dot

function phi_dot(r,z,phi,vpar,mu) result(funcval)
  use precision,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: funcval,r,z,phi,vpar,mu
  real(p_):: bstar_phi,bstar_parallel,bstar_parallelval
  real(p_):: b,bz,b_r,br,b_z

  bstar_parallelval=bstar_parallel(r,z,phi,vpar)
  funcval=bstar_phi(r,z,phi,vpar)/bstar_parallelval*vpar+&
       & mu/(twopi*b(r,z,phi)*bstar_parallelval)*(bz(r,z,phi)*b_r(r,z,phi)-br(r,z,phi)*b_z(r,z,phi))
  funcval=funcval/r
end function phi_dot

function vpar_dot(r,z,phi,vpar,mu) result(funcval)
  use precision,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: funcval,r,z,phi,vpar,mu
  real(p_):: bstar_r,bstar_z,bstar_phi,bstar_parallel,bstar_parallelval
  real(p_):: b_r,b_z,b_phi

  bstar_parallelval=bstar_parallel(r,z,phi,vpar)
  funcval=-mu*(bstar_r(r,z,phi,vpar)/bstar_parallelval*b_r(r,z,phi)+bstar_z(r,z,phi,vpar)/bstar_parallelval*b_z(r,z,phi) &
       & +bstar_phi(r,z,phi,vpar)/bstar_parallelval*b_phi(r,z,phi)/r )
end function


function bstar_r(r,z,phi,vpar) result(funcval)
  use precision,only:p_
  use constants,only:twopi,one
  implicit none
  real(p_):: funcval,r,z,phi,vpar
  real(p_):: br,unitbphi_z,unitbz_phi

  funcval=br(r,z,phi)+vpar/twopi*(unitbz_phi(r,z,phi)/r-unitbphi_z(r,z,phi))

end function bstar_r


function bstar_z(r,z,phi,vpar) result(funcval)
  use precision,only:p_
  use constants,only:twopi,one
  implicit none
  real(p_):: funcval,r,z,phi,vpar
  real(p_):: b,bphi,bz,unitbphi_r,unitbr_phi
  real(p_):: unitbphi

  unitbphi=bphi(r,z,phi)/b(r,z,phi)

  funcval=bz(r,z,phi)+vpar/twopi*(unitbphi_r(r,z,phi)+unitbphi/r-unitbr_phi(r,z,phi)/r)

end function bstar_z

function bstar_phi(r,z,phi,vpar) result(funcval)
  use precision,only:p_
  use constants,only:twopi,one
  implicit none
  real(p_)::funcval,r,z,phi,vpar
  real(p_):: bphi,unitbr_z,unitbz_r

  funcval=bphi(r,z,phi)+vpar/twopi*(unitbr_z(r,z,phi)-unitbz_r(r,z,phi))

end function bstar_phi


function bstar_parallel(r,z,phi,vpar) result(funcval)
  use precision,only:p_
  use constants,only:twopi,one
  implicit none
  real(p_):: funcval,r,z,phi,vpar
  real(p_)::unitb_dot_curl_unitb
  real(p_)::b,br,bz,bphi
  real(p_)::unitbr_z,unitbz_r,unitbphi_r,unitbphi_z
  real(p_):: unitbr_phi,unitbz_phi
  real(p_):: unitbr,unitbphi,unitbz,bval

  bval=b(r,z,phi)
  unitbr=br(r,z,phi)/bval
  unitbz=bz(r,z,phi)/bval
  unitbphi=bphi(r,z,phi)/bval

  unitb_dot_curl_unitb=unitbr*(unitbz_phi(r,z,phi)/r-unitbphi_z(r,z,phi))&
       & +unitbphi*(unitbr_z(r,z,phi)-unitbz_r(r,z,phi))&
       & +unitbz*(unitbphi_r(r,z,phi)+unitbphi/r-unitbr_phi(r,z,phi)/r)
  funcval=bval+vpar/twopi*unitb_dot_curl_unitb
end function bstar_parallel






