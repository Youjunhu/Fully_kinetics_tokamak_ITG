subroutine magnetic_coordinates_derivatives_in_cylindrical() !this subroutine calculates dradcor_dR, dradcor_dZ, dtheta_dR,dtheta_dZ, ddelta_dR,ddelta_dZ on uniform grids in cylindrical coordinates (R,Z). These arrays are to be used as interpolating table in the corresponding interpolation functions
  use precision,only:p_
  use constants,only:twopi,two,pi
  use normalizing,only: Ln
  use magnetic_coordinates,only:r_mag_surf,z_mag_surf,nflux,mpoloidal
  use mapping_module,only: nx_mapping,nz_mapping,r_cyl,z_cyl,radcor,theta_a,tor_shift_a,tor_shift_b,i0,j0 !input
  use mapping_module,only: dradial_dr,dradial_dz,dtheta_dr,dtheta_dz,ddelta_dr_a,ddelta_dz_a,ddelta_dr_b,ddelta_dz_b !as output, which will be used to construct interpolating functions
  use domain_decomposition,only: myid
  use interpolate_module
  implicit none
  real(p_):: r(nx_mapping),z(nz_mapping),ddelta_dz_j0_plus_half,ddelta_dz_j0_minus_half
  integer:: i,j,ierr,inout1(nx_mapping,nz_mapping),inout2(nx_mapping,nz_mapping)
  logical:: within_region(nx_mapping,nz_mapping)

  r=r_cyl/Ln
  z=z_cyl/Ln

  call calculate_2d_partial_derivatives(nx_mapping,nz_mapping,r,z,radcor,dradial_dr,dradial_dz) !after this, functions dradcor_dr_func and dradcor_dz_func are ready to be used

!!$  call calculate_2d_partial_derivatives(nx_mapping,nz_mapping,r,z,theta,dtheta_dr,dtheta_dz) !the following special treatment at theta cut is needed, 
!!$  do i=i0+1,nx_mapping !for low-field-side, and j0 corresponds to midplane, !this is important, to get correct values of dtheta_dz at these special locations
!!$     dtheta_dz(i,j0)= (theta(i,j0+1)+twopi-theta(i,j0-1))/(z(j0+1)-z(j0-1)) !correct
!!$     dtheta_dz(i,j0-1)= (twopi-theta(i,j0-2))/(z(j0)-z(j0-2))
!!$  enddo
!!$  do i=1,i0-1 !j0 corresponds to midplane
!!$     dtheta_dz(i,j0)= (theta(i,j0+1)-twopi-theta(i,j0-1))/(z(j0+1)-z(j0-1)) !correct
!!$     dtheta_dz(i,j0-1)= (-pi-theta(i,j0-2))/(z(j0)-z(j0-2))
!!$  enddo
  !after this, functions dtheta_dr_func and dtheta_dz_func are ready to be used

  call calculate_2d_partial_derivatives(nx_mapping,nz_mapping,r,z,theta_a,dtheta_dr,dtheta_dz) !the following special treatment at theta cut is needed, !theta cut was changed from the low-field-side to high-field side on 2017-Dec. 5
  do i=1,i0-1 !for high-field-side (theta cut is here), j0 corresponds to midplane
     dtheta_dz(i,j0)= (theta_a(i,j0+1)-(theta_a(i,j0-1)+twopi))/(z(j0+1)-z(j0-1)) !correct
     !dtheta_dz(i,j0-1)= (-pi-theta(i,j0-2))/(z(j0)-z(j0-2))
     dtheta_dz(i,j0+1)= (theta_a(i,j0+2)-pi)/(z(j0+2)-z(j0))
  enddo
  !after this, functions dtheta_dr_func and dtheta_dz_func are ready to be used

!!$  ddelta_dr=0._p_
!!$  ddelta_dz=0._p_
!!$  ddelta_dz_b=0._p_

  call calculate_2d_partial_derivatives(nx_mapping,nz_mapping,r,z,tor_shift_a,ddelta_dr_a,ddelta_dz_a)
  call calculate_2d_partial_derivatives(nx_mapping,nz_mapping,r,z,tor_shift_b,ddelta_dr_b,ddelta_dz_b) !the following special treatment at the cut (low-field midplane) is needed before the function ddleta_dr_func, ddelta_dz_func can be used
!!$  do i=i0+1,nx_mapping !for low-field-side, and j0 corresponds to midplane !this is important, to get correct values of ddelta_dz at these special locations
!!$     !ddelta_dz(i,j0)= (tor_shift(i,j0+1)+tor_shift_b(i,j0)-tor_shift(i,j0-1))/(z(j0+1)-z(j0-1)) !wrong, this center-difference is equivalent to averaging the values of derivative at two nearby half-grid of j0, however the derivative itself is not continous accros j0, thus the average does not make sence. This is different from calculating dtheta_dz, which is continuous acros the j0 (midplane) and it makes sence to average the values of derivative at two nearby half-grid of j0 to approximate the value at j0. 
!!$     ddelta_dz_j0_plus_half= (tor_shift(i,j0+1)-tor_shift(i,j0))/(z(j0+1)-z(j0))
!!$     ddelta_dz(i,j0)=two*ddelta_dz_j0_plus_half-ddelta_dz(i,j0+1) !linear extrapolate of ddelta_dz along Z direction (vertical direction)
!!$     ddelta_dz(i,j0-1)= (tor_shift_b(i,j0)-tor_shift_b(i,j0-2))/(z(j0)-z(j0-2))
!!$  enddo
!!$
!!$  ddelta_dz_b=ddelta_dz !ddelta_dz_b is identical to ddelta_dz, except on the midplane
!!$  do i=i0+1,nx_mapping !for low-field-side, and j0 corresponds to midplane !this is important, to get correct values of ddelta_dz at these special locations
!!$     !ddelta_dz_b(i,j0)= (tor_shift(i,j0+1)+tor_shift_b(i,j0)-tor_shift_b(i,j0-1))/(z(j0+1)-z(j0-1)) !wrong, the reason is mentioned above
!!$     ddelta_dz_j0_minus_half= (tor_shift_b(i,j0)-tor_shift_b(i,j0-1))/(z(j0)-z(j0-1))
!!$     ddelta_dz_b(i,j0)=two*ddelta_dz_j0_minus_half-ddelta_dz(i,j0-1) !linear extrapolate of ddelta_dz along Z direction (vertical direction)
!!$  enddo !after this, function ddelta_dz_func are ready to be used
!!$
!!$  ddelta_dr_b=ddelta_dr !ddelta_dr_b is identical to ddelta_dr, except on the midplane, the ddelta_dr has a jump on low-field-side-midplane.
!!$  do i=i0+1,nx_mapping !for low-field-side, and j0 corresponds to midplane
!!$     ddelta_dr_b(i,j0)=two*ddelta_dr(i,j0-1)-ddelta_dr(i,j0-2) !linear extrapolate of ddelta_dr along the Z direction (vertical direction)
!!$  enddo !after this, function ddelta_dr_func are ready to be used

!the above "calculate_2d_partial_derivatives" does not take into account the cut (high-field-side midplane), so we need manually provide the correct values of the derivatives near the cut
  do i=1,i0-1 !for high-field-side midplane (theta cut), and j0 corresponds to midplane
     ddelta_dz_a(i,j0)=two*ddelta_dz_a(i,j0-1)-ddelta_dz_a(i,j0-2) !linear extrapolate of ddelta_dz along Z direction (vertical direction), to get the left (below the midplane) derivative of delta over z
     ddelta_dz_a(i,j0+1)=(tor_shift_b(i,j0+2)-tor_shift_b(i,j0))/(z(j0+2)-z(j0))
  enddo

  do i=1,i0-1 !for high-field-side midplane (theta cut), and j0 corresponds to midplane
     ddelta_dz_b(i,j0)=two*ddelta_dz_b(i,j0+1)-ddelta_dz_b(i,j0+2) !linear extrapolate of ddelta_dz along Z direction (vertical direction), to get the right (above the midplane) derivative of delta over z
     ddelta_dz_b(i,j0-1)=(tor_shift_a(i,j0)-tor_shift_a(i,j0-2))/(z(j0)-z(j0-2))
  enddo

  do i=1,i0-1 !for high-field-side midplane (theta cut), j0 corresponds to midplane
     ddelta_dr_a(i,j0)=two*ddelta_dr_a(i,j0-1)-ddelta_dz_a(i,j0-2) !linear extrapolate of ddelta_dr along Z direction (vertical direction), to get the left (below the midplane) derivative of delta over r
  enddo

  do i=1,i0-1 !for high-field-side midplane (theta cut), j0 corresponds to midplane
     ddelta_dr_b(i,j0)=two*ddelta_dr_b(i,j0+1)-ddelta_dz_b(i,j0+2) !linear extrapolate of ddelta_dr along Z direction (vertical direction), to get the right (above the midplane) derivative of delta over r
  enddo
!Note that ddelta_dr_a(i,j0) is value of the left (i.e., below the midplane) derivative of delta at the midplane while ddelta_dr_b(i,j0) is value of the right (i.e., above the midplane) derivative of delta at the midplane, these two limits turn out to be differnt, which implies that the derivative of delta at the theta cut does not exist.

!!$  if(myid.eq.0) write(*,*) maxval(ddelta_dz),maxval(ddelta_dz_b),maxloc(ddelta_dz),maxloc(ddelta_dz_b)
  if(myid.eq.0) then !for test
     within_region=.true.
     do i=1,nx_mapping
        do j=1,nz_mapping !check whether a point is within the specifed region or not.
           call PNPOLY(r_cyl(i),z_cyl(j),r_mag_surf(:,nflux),z_mag_surf(:,nflux),mpoloidal,INOUT1(i,j))
           call PNPOLY(r_cyl(i),z_cyl(j),r_mag_surf(:,1),z_mag_surf(:,1),mpoloidal,INOUT2(i,j))
           if((inout1(i,j).eq. -1) .or. (inout2(i,j).eq.1)) within_region(i,j)=.false.
        enddo
     enddo
     open(124,file='r_z_ddelta_dr_dz.txt')
     do i=1,nx_mapping
        do j=1,nz_mapping 
           if(within_region(i,j).eqv..true.) then
              write(124,*) r_cyl(i),z_cyl(j),ddelta_dr_a(i,j),ddelta_dz_a(i,j),ddelta_dr_b(i,j),ddelta_dz_b(i,j)
           else
              write(124,*) r_cyl(i),z_cyl(j), 'NaN',' NaN', ' NAN'
           endif
        enddo
        write(124,*) 
     enddo
     close(124)
  endif
end subroutine magnetic_coordinates_derivatives_in_cylindrical

!tested, dradial_dr_fun and dradial_dr_func2 agree with each other.Also dradial_dz_fun and dradial_dz_func2 agree with each other. 
function dradial_dr_func2(r,z) result(funcval) 
  use precision,only:p_
  use normalizing,only: Ln
  use mapping_module,only:nx_mapping,nz_mapping,r_cyl,z_cyl,dradial_dr
  use interpolate_module
  implicit none
  real(p_):: r,z,funcval
  call linear_2d_interpolation(nx_mapping,nz_mapping,r_cyl,z_cyl,dradial_dr,r*Ln,z*Ln,funcval)  
end function

function dradial_dz_func2(r,z) result(funcval)
  use precision,only:p_
  use normalizing,only: Ln
  use mapping_module,only:nx_mapping,nz_mapping,r_cyl,z_cyl,dradial_dz
  use interpolate_module
  implicit none
  real(p_):: r,z,funcval
  call linear_2d_interpolation(nx_mapping,nz_mapping,r_cyl,z_cyl,dradial_dz,r*Ln,z*Ln,funcval)  

end function

function dradial_dr_func(r,z) result(funcval)
  use precision,only:p_
  use normalizing,only: ln
  use radial_module,only:psi_lcfs,psi_axis
  implicit none
  real(p_):: r,z,funcval
  real(p_):: psi_r_func
  real(p_):: psi_prime

  psi_prime=psi_lcfs-psi_axis
  funcval=psi_r_func(r*Ln,z*Ln)/psi_prime

end function


function dradial_dz_func(r,z) result(funcval)
  use precision,only:p_
  use normalizing,only: ln
  use radial_module,only:psi_lcfs,psi_axis
  implicit none
  real(p_):: r,z,funcval
  real(p_):: psi_z_func
  real(p_):: psi_prime
  psi_prime=psi_lcfs-psi_axis
  funcval=psi_z_func(r*Ln,z*Ln)/psi_prime
end function

function dtheta_dr_func(r,z) result(funcval)
  use precision,only:p_
  use normalizing,only: Ln
  use mapping_module,only:nx_mapping,nz_mapping,r_cyl,z_cyl,dtheta_dr
  use interpolate_module
  implicit none
  real(p_):: r,z,funcval
!  if(r<minval(r_cyl)) write(*,*)  'r, minval(r_cyl)=',r, minval(r_cyl)
!  if(z<minval(z_cyl)) write(*,*)  'z, minval(z_cyl)=',z, minval(z_cyl)
  call linear_2d_interpolation(nx_mapping,nz_mapping,r_cyl,z_cyl,dtheta_dr,r*Ln,z*Ln,funcval)  
end function dtheta_dr_func

function dtheta_dz_func(r,z) result(funcval)
  use precision,only:p_
  use normalizing,only: Ln
  use mapping_module,only:nx_mapping,nz_mapping,r_cyl,z_cyl,dtheta_dz
  use interpolate_module
  implicit none
  real(p_):: r,z,funcval
  call linear_2d_interpolation(nx_mapping,nz_mapping,r_cyl,z_cyl,dtheta_dz,r*Ln,z*Ln,funcval)  
end function dtheta_dz_func

function ddelta_dr_func(r,z) result(funcval)
  use precision,only:p_
  use normalizing,only: Ln
  use mapping_module,only:nx_mapping,nz_mapping,r_cyl,z_cyl,ddelta_dr_a,ddelta_dr_b
  use interpolate_module
  implicit none
  real(p_):: r,z,funcval
!  call linear_2d_interpolation(nx_mapping,nz_mapping,r_cyl,z_cyl,ddelta_dr,r*Ln,z*Ln,funcval)  
 call linear_2d_interpolation_special(nx_mapping,nz_mapping,r_cyl,z_cyl,ddelta_dr_a,ddelta_dr_b,r*Ln,z*Ln,funcval)  
end function ddelta_dr_func

function ddelta_dz_func(r,z) result(funcval)
  use precision,only:p_
  use normalizing,only: Ln
  use mapping_module,only:nx_mapping,nz_mapping,r_cyl,z_cyl,ddelta_dz_a,ddelta_dz_b
  use interpolate_module
  implicit none
  real(p_):: r,z,funcval
  call linear_2d_interpolation_special(nx_mapping,nz_mapping,r_cyl,z_cyl,ddelta_dz_a,ddelta_dz_b,r*Ln,z*Ln,funcval)  
end function ddelta_dz_func


subroutine calculate_2d_partial_derivatives(nx,nz,xarray,zarray,psi,psi_x,psi_z)
  use precision,only:p_
  implicit none
  integer,intent(in)::nx,nz
  real(p_),intent(in):: psi(nx,nz)
  real(p_),intent(in):: xarray(nx),zarray(nz)
  real(p_),intent(out):: psi_x(nx,nz),psi_z(nx,nz)

  integer:: i,j,i1,i2,j1,j2

  !first-order partial derivatives
  do i=1,nx
     do j=1,nz
        i2=i+1
        i1=i-1
        j2=j+1
        j1=j-1
        if(i.eq.1) i1=i
        if(j.eq.1) j1=j
        if(i.eq.nx) i2=i
        if(j.eq.nz) j2=j
        psi_x(i,j)=(psi(i2,j)-psi(i1,j))/(xarray(i2)-xarray(i1))
        psi_z(i,j)=(psi(i,j2)-psi(i,j1))/(zarray(j2)-zarray(j1))
     enddo
  enddo

end subroutine calculate_2d_partial_derivatives

function ddelta_dz_lsf_midplane_twopi(r) result(funcval)
  use precision,only:p_
  use mapping_module,only:nx_mapping,r_cyl,ddelta_dz_b
  use mapping_module,only: i0, j0 !index of the point at magnetic axis
  use interpolate_module
  implicit none
  real(p_):: r,funcval
call linear_1d_interpolation(nx_mapping,r_cyl,ddelta_dz_b(:,j0),r,funcval)
end function ddelta_dz_lsf_midplane_twopi

function ddelta_dr_lsf_midplane_twopi(r) result(funcval)
  use precision,only:p_
  use mapping_module,only:nx_mapping,r_cyl,ddelta_dr_b
  use mapping_module,only: i0, j0 !index of the point at magnetic axis
  use interpolate_module
  implicit none
  real(p_):: r,funcval
call linear_1d_interpolation(nx_mapping,r_cyl,ddelta_dr_b(:,j0),r,funcval)
end function ddelta_dr_lsf_midplane_twopi
