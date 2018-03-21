subroutine construct_poloidal_coordinate(r_old,z_old,np_lcfs,r_new,z_new,mpoloidal,theta_uniform) !to construct poloidal angle coordinate
  use precision,only: p_
  use constants,only: two,pi,twopi
  use control_parameters,only: poloidal_angle_type
  implicit none
  integer,intent(in):: np_lcfs,mpoloidal
  real(p_),intent(in):: r_old(np_lcfs),z_old(np_lcfs)
  real(p_),intent(in):: theta_uniform(mpoloidal)
  real(p_),intent(out):: r_new(mpoloidal),z_new(mpoloidal)
  real(p_):: theta_old(np_lcfs),dl(np_lcfs-1)
  real(p_):: rmid,zmid,psi_gradient_func
  real(p_):: y2(np_lcfs),y_tmp
  integer:: i

  call arc_length(r_old,z_old,np_lcfs,dl)
  !calculate poloidal angle 
  theta_old(1)=0._p_
  if(poloidal_angle_type .eq. 'equal-arc') then
     do i=2,np_lcfs
        theta_old(i)=theta_old(i-1)+dl(i-1) !equal-arc-length poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'equal-volume') then
     do i=2,np_lcfs
        rmid=0.5_p_*(r_old(i-1)+r_old(i))
        zmid=0.5_p_*(z_old(i-1)+z_old(i))
        theta_old(i)=theta_old(i-1)+dl(i-1)*rmid/psi_gradient_func(rmid,zmid) !equal-volume poloidal angle
     enddo
  else
     stop 'please choose poloidal angle type between equal-arc and equal-volume'
  endif

  !theta_old=theta_old*twopi/theta_old(np_lcfs) !normalized to the range [0:twopi]
   theta_old=theta_old*twopi/theta_old(np_lcfs)-pi !normalized to the range [-pi:pi]

  !The theta array obtained above is usually not uniform. Next, interpolate R and Z to uniform theta grids on every magnetic surfac.
  !interpolate R  to uniform theta grid points
  call spline(theta_old,r_old,np_lcfs,2.d30,2.d30,y2) !prepare the second order derivative needed in the cubic spline interpolation
  do i=2,mpoloidal-1
     call splint(theta_old,r_old,y2,np_lcfs,theta_uniform(i),y_tmp) !to get R corresponding uniform theta grids
     r_new(i)=y_tmp
  enddo
  !interpolate Z  to uniform theta grid points
  call spline(theta_old,z_old,np_lcfs,2.d30,2.d30,y2) !prepare the second order derivative needed in the cubic spline interpolation
  do i=2,mpoloidal-1
     call splint(theta_old,z_old,y2,np_lcfs,theta_uniform(i),y_tmp) !to get Z corresponding uniform theta grids
     z_new(i)=y_tmp
  enddo

  r_new(1)=r_old(1) !ending points are not included in the above interpolation since we know the anwser (as given here)
  z_new(1)=z_old(1)
  r_new(mpoloidal)=r_old(np_lcfs)
  z_new(mpoloidal)=z_old(np_lcfs)

end subroutine construct_poloidal_coordinate
