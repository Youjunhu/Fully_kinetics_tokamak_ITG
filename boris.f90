subroutine push_full_orbit_cylindrical_boris(dtao,r,phi,z,vr,vphi,vz)
  !input: initial condition of the orbit: r,phi,z,vr,vphi,vz
  !Output: the instanteous value of the orbit after dtao: r,phi,z,vr,vphi,vz
!in essence, the algorithm uses Cartesian coordinates, and then takes into account the rotation of the the basis vectors due to the change of particle's location
  use precision,only:p_
  use constants,only:zero,one,two,twopi
  implicit none
  real(p_),intent(in):: dtao
  real(p_),intent(inout):: r,phi,z,vr,vphi,vz  !instantaneous value of orbit
  real(p_):: br,bz,bphi !function names
  real(p_):: er,ez,ephi !function names
  real(p_):: tx,ty,tz,t,factor,sx,sy,sz
  real(p_):: vx_minus,vy_minus,vz_minus
  real(p_):: vx_prime,vy_prime,vz_prime
  real(p_):: vx_plus,vy_plus,vz_plus
  real(p_):: cx,cy,cz
  real(p_):: x,y,vx,vy
  real(p_):: alpha


  vx_minus=vr  +  er  (r,z,phi)*dtao/two*(twopi)
  vy_minus=vphi+  ephi(r,z,phi)*dtao/two*(twopi)
  vz_minus=vz  +  ez  (r,z,phi)*dtao/two*(twopi)

  tx=br(r,z,phi)*dtao/two*(twopi)
  ty=bphi(r,z,phi)*dtao/two*(twopi)
  tz=bz(r,z,phi)*dtao/two*(twopi)

  call cross_product_in_cartesian(vx_minus,vy_minus,vz_minus,tx,ty,tz,cx,cy,cz)

  vx_prime=  vx_minus+ cx
  vy_prime=  vy_minus+ cy
  vz_prime=  vz_minus+ cz

  t=sqrt(tx*tx+tz*tz+ty*ty)
  factor=two/(one+t*t)
  sx=tx*factor
  sy=ty*factor
  sz=tz*factor

  call cross_product_in_cartesian(vx_prime,vy_prime,vz_prime,sx,sy,sz,cx,cy,cz)

  vx_plus =  vx_minus+ cx
  vy_plus =  vy_minus+ cy
  vz_plus =  vz_minus+ cz

  vx=vx_plus+  er  (r,z,phi)*dtao/two*(twopi)
  vy=vy_plus+  ephi(r,z,phi)*dtao/two*(twopi)
  vz=vz_plus+  ez  (r,z,phi)*dtao/two*(twopi)


  x=r+vx*dtao
  y=vy*dtao
  z=z+vz*dtao


  r=sqrt(x*x+y*y)

  alpha=asin(y/r)

  phi=phi+alpha

  vr=cos(alpha)*vx+sin(alpha)*vy
  vphi=-sin(alpha)*vx+cos(alpha)*vy
end subroutine push_full_orbit_cylindrical_boris


subroutine push_full_orbit_cylindrical_boris_with_additional_input_output(dtao,r,phi,z,vr,vphi,vz,phi_mid,vr_mid,vphi_mid,vz_mid)
!derived from "push_full_orbit_cylindrical_boris", the modification is that: additional input phi_mid, which is value of phi at t_{n+1/2}, and additional output: vr_mid,vphi_mid,vz_mid, which are projection of velocity at t_{n+1/2} onto the basis vector at t_{n+1/2}, instead of t_{n+1}. Here, for general cases, n can also be an half-integer
  use precision,only:p_
  use constants,only:zero,one,two,twopi
  implicit none
  real(p_),intent(in):: dtao
  real(p_),intent(inout):: r,phi,z,vr,vphi,vz  !instantaneous value of (r,phi,z)  at t_{n} (as input) and t_{n+1} (as output), (vr,vphi,vz) is the projection of velocity at t_{n-1/2} at basis vector at t_{n} (as input); and the projection of velocity at t_{n+1/2} on the basis vector at t_{n+1} (as output}
  real(p_),intent(in):: phi_mid !value of phi at t_{n+1/2}
  real(p_),intent(out):: vr_mid,vphi_mid,vz_mid !projection velocity at t_{n+1/2} on the basis vector at t_{n+1/2}
  real(p_):: br,bz,bphi !function names
  real(p_):: er,ez,ephi !function names
  real(p_):: tx,ty,tz,t,factor,sx,sy,sz
  real(p_):: vx_minus,vy_minus,vz_minus
  real(p_):: vx_prime,vy_prime,vz_prime
  real(p_):: vx_plus,vy_plus,vz_plus
  real(p_):: cx,cy,cz
  real(p_):: x,y,vx,vy
  real(p_):: dphi,phi0

phi0=phi !record the old value of phi
  vx_minus=vr  +  er  (r,z,phi)*dtao/two*(twopi)
  vy_minus=vphi+  ephi(r,z,phi)*dtao/two*(twopi)
  vz_minus=vz  +  ez  (r,z,phi)*dtao/two*(twopi)

  tx=br(r,z,phi)*dtao/two*(twopi)
  ty=bphi(r,z,phi)*dtao/two*(twopi)
  tz=bz(r,z,phi)*dtao/two*(twopi)

  call cross_product_in_cartesian(vx_minus,vy_minus,vz_minus,tx,ty,tz,cx,cy,cz)

  vx_prime=  vx_minus+ cx
  vy_prime=  vy_minus+ cy
  vz_prime=  vz_minus+ cz

  t=sqrt(tx*tx+tz*tz+ty*ty)
  factor=two/(one+t*t)
  sx=tx*factor
  sy=ty*factor
  sz=tz*factor

  call cross_product_in_cartesian(vx_prime,vy_prime,vz_prime,sx,sy,sz,cx,cy,cz)

  vx_plus =  vx_minus+ cx
  vy_plus =  vy_minus+ cy
  vz_plus =  vz_minus+ cz

  vx=vx_plus+  er  (r,z,phi)*dtao/two*(twopi)
  vy=vy_plus+  ephi(r,z,phi)*dtao/two*(twopi)
  vz=vz_plus+  ez  (r,z,phi)*dtao/two*(twopi)


  x=r+vx*dtao
  y=vy*dtao
  z=z+vz*dtao


  r=sqrt(x*x+y*y)

  dphi=asin(y/r)

  phi=phi+dphi

  vr=cos(dphi)*vx+sin(dphi)*vy
  vphi=-sin(dphi)*vx+cos(dphi)*vy

dphi=phi_mid-phi0
  vr_mid=cos(dphi)*vx+sin(dphi)*vy !the projection of v_{n+1/2} on basis vector at t_{n+1/2}
  vphi_mid=-sin(dphi)*vx+cos(dphi)*vy
  vz_mid=vz

end subroutine push_full_orbit_cylindrical_boris_with_additional_input_output


function er(r,z,phi) !R component of electric field
  use precision,only:p_
  implicit none
  real(p_):: er,r,z,phi

  er=0._p_

end function er


function ez(r,z,phi) !R component of electric field
  use precision,only:p_
  implicit none
  real(p_):: ez,r,z,phi

  ez=0._p_

end function ez


function ephi(r,z,phi) !R component of electric field
  use precision,only:p_
  implicit none
  real(p_):: ephi,r,z,phi

  ephi=0._p_

end function ephi



subroutine cross_product_in_cartesian(ax,ay,az,bx,by,bz,cx,cy,cz)
  use precision,only:p_
  implicit none

  real(p_),intent(in):: ax,ay,az,bx,by,bz
  real(p_),intent(out)::cx,cy,cz

  cx=ay*bz-az*by
  cy=az*bx-ax*bz
  cz=ax*by-ay*bx
end subroutine cross_product_in_cartesian


subroutine normalize_full_orbit_variables(nmarker,r,phi,z,vr,vphi,vz)
 use precision,only:p_
  use normalizing,only: Ln,vn_i

  implicit none
  integer,intent(in):: nmarker
  real(p_),intent(inout):: r(nmarker),phi(nmarker),z(nmarker)
  real(p_),intent(inout):: vr(nmarker),vphi(nmarker),vz(nmarker)

 !normalization
  r=r/Ln !convert to unit Ln
  z=z/Ln !convert to unit Ln
!  phi=phi
  vr=vr/vn_i !normalized by vn_i
  vphi=vphi/vn_i
  vz=vz/vn_i
end subroutine normalize_full_orbit_variables


