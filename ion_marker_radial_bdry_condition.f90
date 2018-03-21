subroutine ion_marker_radial_bdry_condition_original(nmarker_i,r_i,z_i,radcor_i,active_i,w_i) !get the radial coordinates of markers, and place a marker at a random location within the buffer region if the marker is outside of the buffer region
  use precision,only:p_
  use constants,only:pi,twopi
  use poloidal_flux_2d,only:xarray,zarray,nx,nz
  use magnetic_coordinates,only: radcor_low1,radcor_upp1,radcor_low2,radcor_upp2
  implicit none
  integer,intent(in):: nmarker_i
  real(p_),intent(inout):: r_i(nmarker_i),z_i(nmarker_i)
  real(p_),intent(out):: radcor_i(nmarker_i),w_i(nmarker_i)
  logical,intent(out):: active_i(nmarker_i)
  integer:: k,flag,next_seed
  real(p_):: theta,rannum1,rannum2,radcor_as_func_of_pfn,pfn_func
  logical:: outside_box

  active_i=.false.
  do k=1,nmarker_i
     flag=0
     outside_box=r_i(k).ge.xarray(nx) .or.r_i(k).le.xarray(1) .or. z_i(k).ge.zarray(nz) .or.  z_i(k).le.zarray(1)
     if(outside_box) then
        flag=1      
     else
        radcor_i(k)=radcor_as_func_of_pfn(pfn_func(r_i(k),z_i(k))) !get radial coordinate
        if(radcor_i(k).gt.radcor_upp1) flag=1     
        if(radcor_i(k).lt.radcor_low1) flag=-1 
        if(radcor_i(k).gt.radcor_low2 .and. radcor_i(k).lt.radcor_upp2) active_i(k)=.true.
     endif

     if(flag.ne.0) then !marker is outside the buffer region, re-introducing is needed
        w_i(k)=0._p_ !set weight to zero
        call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
        call sub_random_yj(0,next_seed,rannum2) !use last random number as iseed
        if(flag.eq.1)  radcor_i(k)=radcor_upp2+rannum1*(radcor_upp1-radcor_upp2) !place the marker at a random location in the outer buffer region
        if(flag.eq.-1) radcor_i(k)=radcor_low1+rannum1*(radcor_low2-radcor_low1) !place the marker at a random location in the inner buffer region
        theta=-pi+rannum2*twopi
        call magnetic_coordinates_to_cylindrical_coordinates(theta,radcor_i(k),r_i(k),z_i(k)) !given (theta,radcor), return (R,Z)
        !write(*,*) 'marker is outside of buffer region, place it at a random location in the buffer region',r_i(k),z_i(k)
!write(*,*) '00000',r_i(k),z_i(k)
     endif
  enddo
end subroutine ion_marker_radial_bdry_condition_original


subroutine ion_marker_radial_bdry_condition(nmarker_i,r_i,z_i,radcor_i,active_i,touch_bdry_i,w_i) !get the radial coordinates of markers, and place a marker at a random location within the buffer region if the marker is outside of the buffer region
  use precision,only:p_
  use constants,only:pi,twopi
  use poloidal_flux_2d,only:xarray,zarray,nx,nz
  use magnetic_coordinates,only:radcor_low2,radcor_upp2,radcor_low0,radcor_upp0
  implicit none
  integer,intent(in):: nmarker_i
  real(p_),intent(inout):: r_i(nmarker_i),z_i(nmarker_i)
  real(p_),intent(out):: radcor_i(nmarker_i),w_i(nmarker_i)
  logical,intent(out):: active_i(nmarker_i),touch_bdry_i(nmarker_i)
  integer:: k,flag,next_seed
  real(p_):: theta,rannum1,rannum2,radcor_as_func_of_pfn,pfn_func
  logical:: outside_box

  do k=1,nmarker_i
     outside_box=r_i(k).ge.xarray(nx) .or.r_i(k).le.xarray(1) .or. z_i(k).ge.zarray(nz) .or.  z_i(k).le.zarray(1) 
     if(outside_box.eqv..true.) then
        touch_bdry_i(k)=.true. !marker is lost forever, will never be re-introduced
        active_i(k)=.false.
     else
        radcor_i(k)=radcor_as_func_of_pfn(pfn_func(r_i(k),z_i(k))) !get radial coordinate
        if(radcor_i(k).gt.radcor_upp0 .or. radcor_i(k).lt.radcor_low0) then 
           touch_bdry_i(k)=.true. !marker is lost forever, will never be re-introduced
           active_i(k)=.false.
        else if(radcor_i(k).lt.radcor_low2 .or. radcor_i(k).gt.radcor_upp2) then
           touch_bdry_i(k)=.false.
           active_i(k)=.false.
           z_i(k)=-z_i(k) !relocate the marker by up-down reversing
           radcor_i(k)=radcor_as_func_of_pfn(pfn_func(r_i(k),z_i(k))) !re-calculate the radial coordinate, which is different from the original value if the flux-surface is not up-down symmetric
        else
           touch_bdry_i(k)=.false.
           active_i(k)=.true.
        endif
     endif
  enddo

  do k=1,nmarker_i
     if(active_i(k).eqv..false.) w_i(k)=0._p_
  enddo

end subroutine ion_marker_radial_bdry_condition
