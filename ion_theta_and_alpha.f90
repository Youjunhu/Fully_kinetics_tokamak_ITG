subroutine ion_theta_and_alpha(nmarker_i,r_i,phi_i,z_i, radcor_i,active_i,touch_bdry_i,theta_i,alpha_i)
  use precision,only:p_
  use magnetic_coordinates,only:mpoloidal,nflux,theta_1d_array,radcor_1d_array,tor_shift_mc
  use domain_decomposition,only:myid
  use interpolate_module

  implicit none
  integer,intent(in):: nmarker_i
  real(p_),intent(in):: r_i(nmarker_i),z_i(nmarker_i),phi_i(nmarker_i)
  real(p_),intent(in):: radcor_i(nmarker_i)
  logical,intent(in):: active_i(nmarker_i),touch_bdry_i(nmarker_i)
  real(p_),intent(out):: theta_i(nmarker_i),alpha_i(nmarker_i)
  real(p_):: tor_shift
  integer:: k

  do k=1,nmarker_i
     !if(active_i(k).eqv..false.) cycle !wrong, inactive markers' theta must be computed so that they can be sorted by the sorting subroutine
     if(touch_bdry_i(k).eqv..true.) cycle
     call interpolate_from_cylindrical_to_magnetic_coordinates1(r_i(k),z_i(k),theta_i(k))
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,tor_shift_mc,&
          & theta_i(k),radcor_i(k),tor_shift) !interpolating in magnetic coordinates to get tor_shift
     alpha_i(k)=phi_i(k)-tor_shift !generalized toroidal angle
     call shift_to_specified_toroidal_range(alpha_i(k))
  enddo

end subroutine ion_theta_and_alpha


subroutine check_whether_marker_in_boundary(radcor,touch_bdry,active)
  use precision,only:p_
  use magnetic_coordinates,only: radcor_low1,radcor_upp1,radcor_low2,radcor_upp2 !as input
  implicit none
  real(p_),intent(in)::  radcor
  logical,intent(out):: touch_bdry,active

  if(radcor.gt.radcor_upp1 .or. radcor.lt.radcor_low1)  then
     touch_bdry=.true.
     active=.false.
  else if(radcor.gt.radcor_upp2 .or. radcor.lt.radcor_low2)  then
     touch_bdry=.false.
     active=.false.
  else
     touch_bdry=.false.
     active=.true.
  endif

!!$  if(radcor.le.0.) then 
!!$     write(*,*) 'in checking subroutine,warning****************', 'radcor',radcor
!!$     write(*,*) 'radcor_low2,radcor_up2=',radcor_low2,radcor_upp2
!!$     write(*,*) 'touch_bdry=',touch_bdry
!!$  endif
end subroutine check_whether_marker_in_boundary
