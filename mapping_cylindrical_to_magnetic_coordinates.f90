subroutine mapping_cylindrical_to_magnetic_coordinates()
  !used to use r_mag_surf, z_mag_surf in meter, later r_mag_surf z_mag_surf was normalized by Ln before entering this subroutine, since presently I use Ln=1meter, therefore this subroutine remains correct
  !  use mpi
  use precision,only:p_
  !  use magnetic_coordinates,only: mpoloidal,nflux, r_mag_surf, z_mag_surf
  !  use radial_module,only: r_axis,z_axis
  use constants,only:pi
  use boundary,only: np_lcfs !,x_lcfs,z_lcfs
  use mapping_module,only:nx_mapping,nz_mapping
  use mapping_module,only:r_cyl,z_cyl,dr,dz,radcor,i0,j0 !as output
  use mapping_module,only: theta_a,theta_b,tor_shift_a,tor_shift_b !as output
  use domain_decomposition,only:myid
  implicit none
  logical:: within_region(nx_mapping,nz_mapping)
  integer:: i,j,inout1(nx_mapping,nz_mapping),inout2(nx_mapping,nz_mapping)
  real(p_):: theta(nx_mapping,nz_mapping),tor_shift(nx_mapping,nz_mapping)
  real(p_):: psi_func,q_func
  !  integer:: ierr
  real(p_):: r_inner_surf(np_lcfs),z_inner_surf(np_lcfs),r_outer_surf(np_lcfs),z_outer_surf(np_lcfs)

  call choose_boundary_magnetic_surfaces_for_the_mapping(r_inner_surf,z_inner_surf,r_outer_surf,z_outer_surf)
  call create_cylindrical_grids(r_outer_surf,z_outer_surf,np_lcfs,nx_mapping,nz_mapping,r_cyl,z_cyl,dr,dz,i0,j0)

  do i=1,nx_mapping
     do j=1,nz_mapping
!!$        call PNPOLY(r_cyl(i),z_cyl(j),r_mag_surf(:,nflux),z_mag_surf(:,nflux),mpoloidal,INOUT1(i,j))
!!$        call PNPOLY(r_cyl(i),z_cyl(j),r_mag_surf(:,2),z_mag_surf(:,2),mpoloidal,INOUT2(i,j))
        call PNPOLY(r_cyl(i),z_cyl(j),r_outer_surf,z_outer_surf,np_lcfs,INOUT1(i,j))
        call PNPOLY(r_cyl(i),z_cyl(j),r_inner_surf,z_inner_surf,np_lcfs,INOUT2(i,j))
     enddo
  enddo

  within_region=.true.
  do i=1,nx_mapping !check whether a point is within the specifed region or not.
     do j=1,nz_mapping
        if((inout1(i,j).eq. -1) .or. (inout2(i,j).eq.1)) within_region(i,j)=.false.
     enddo
  enddo
  !  if(myid.eq.0) write(*,*) 'q_edge=',q_func(psi_func(r_mag_surf(1,nflux),z_mag_surf(1,nflux))) !value of q at boundary
!!$  if(myid.eq.0) then
!!$     open(11,file='mapping_in.txt')
!!$     open(12,file='mapping_out.txt')
!!$     do i=1,nx_mapping
!!$        do j=1,nz_mapping
!!$           if(within_region(i,j) .eqv. .false.) then 
!!$              write(12,*) r_cyl(i),z_cyl(j) !out
!!$           else
!!$              write(11,*) r_cyl(i),z_cyl(j) !in
!!$           endif
!!$        enddo
!!$     enddo
!!$     close(11)
!!$     close(12)
!!$  endif

  !  call arrange_lcfs(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis)

  !     open(123,file='t_map.txt')

  radcor=0._p_ !initialized
  theta=0._p_
  tor_shift=0._p_

  do i=1,nx_mapping
     do j=1,nz_mapping
        if(i<i0 .and. j.eq.j0) then !theta cut, special treatment is needed here
           !done at the end of this subroutine
        else
           if(within_region(i,j).eqv..true.)  call mapping(r_cyl(i),z_cyl(j),radcor(i,j),theta(i,j),tor_shift(i,j))  !calculate (radcor, theta, tor_shift) of the point (r_cyl(i),z_cyl(j)).
        endif
     enddo
  enddo
  !close(123)

  theta_a=theta
  theta_b=theta
  tor_shift_a=tor_shift
  tor_shift_b=tor_shift

  do i=1,nx_mapping
     do j=1,nz_mapping
        if(i<i0 .and. j.eq.j0) then !theta cut, special treatment is needed here
           theta_a(i,j)=-pi
           theta_b(i,j)=pi
           call mapping2(r_cyl(i),z_cyl(j),radcor(i,j),tor_shift_a(i,j),tor_shift_b(i,j))
        endif
     enddo
  enddo

!!$
!!$  tor_shift_b=tor_shift
!!$  do i=i0,nx_mapping
!!$     j=j0 !re-calculate tor_shift_angle at low-field-side midplane (2pi*q, instead of zero) and store them in tor_shift_b array
!!$     if(within_region(i,j).eqv..true.)   call mapping2(r_cyl(i),z_cyl(j),tor_shift_b(i,j))
!!$  enddo


!!$  tor_shift_b=tor_shift
!!$  do i=1,nx_mapping
!!$     do j=1,nz_mapping
!!$        if(i<i0 .and. j.eq.j0) then !the value of theta and tor_shift at the high-field-side midplane
!!$           if(within_region(i,j).eqv..true.)   call mapping2(r_cyl(i),z_cyl(j),tor_shift_b(i,j))
!!$           if(within_region(i,j).eqv..true.)   call mapping2(r_cyl(i),z_cyl(j),tor_shift_b(i,j))
!!$        endif
!!$     enddo
!!$  enddo

  if(myid.eq.0) then
     open(123,file='r_z_theta_tor_shift.txt')
     do i=1,nx_mapping
        do j=1,nz_mapping 
           !do j=j0,j0
           if(within_region(i,j).eqv..true.) then
              write(123,*) r_cyl(i),z_cyl(j),theta_a(i,j),theta_b(i,j),tor_shift_a(i,j),tor_shift_b(i,j)
              !write(123,*) r_cyl(i),z_cyl(j),psi_func(r_cyl(i),z_cyl(j))
           else
              write(123,*) r_cyl(i),z_cyl(j), 'NaN',' NaN',' NaN',' NaN'
           endif
        enddo
        write(123,*) 
     enddo
     close(123)
  endif
!  write(*,*) 'maximum of tor_shift=',maxval(tor_shift_a),'minimum of tor_shift=',minval(tor_shift_a)
!  write(*,*) 'maximum of tor_shift=',maxval(tor_shift_b),'minimum of tor_shift=',minval(tor_shift_b)
end subroutine mapping_cylindrical_to_magnetic_coordinates


subroutine mapping(r,z,radcor,theta,tor_shift)
  !given (R,Z), this subroutine finds the magnetic surface that passes through the point and calculates its poloidal angle and toroidal shift
  use precision,only:p_
  use constants,only:zero,one,two,twopi,pi
  use boundary, only: x_lcfs,z_lcfs,np_lcfs
  use radial_module,only:r_axis,z_axis,psi_axis,psi_lcfs
  use control_parameters,only: poloidal_angle_type
  use domain_decomposition,only:myid
  implicit none
  real(p_),intent(in):: r,z
  real(p_),intent(out):: radcor,theta,tor_shift

  real(p_)::psival,psi_func,radcor_as_func_of_pfn
  real(p_):: x_contour(np_lcfs+1),z_contour(np_lcfs+1)
  real(p_)::dl(np_lcfs), sum

  real(p_),parameter:: xacc=1.0d-6 !tolerance used in bi-section root-finder
  real(p_):: x1,x2,z1,z2

  real(p_):: slope(np_lcfs),slope2(np_lcfs)
  real(p_):: rtbis !function name of the root finder using the bisection method
  real(p_):: zfunc,xfunc !equation of the straight line (in poloidal plane) that passing throught the magnetic axis point and one point on LCFS
  real(p_):: one_dim_psi_func,one_dim_psi_func2 !one dimension function [psi(x,z(x)) and psi(x(z),z)]on the straight line mentioned in the above.
  external:: one_dim_psi_func,one_dim_psi_func2 !this two function will be passed to a root-finding subroutine
  integer:: i,end_i !,ierr
  real(p_):: value1,value2,value3, rmid,zmid,normalization,psi_gradient_func
  real(p_),parameter:: large_number=1d30

  psival=psi_func(r,z)
  radcor=radcor_as_func_of_pfn((psival-psi_axis)/(psi_lcfs-psi_axis))

!!$  do i=1,np_lcfs
!!$     if(x_lcfs(i)-r_axis .eq. 0._p_) then !since I use compiler option which catches all erroneous arithmetic operation, I need to avoid dividing by zero
!!$        slope(i)= large_number
!!$     else
!!$        slope(i)= (z_lcfs(i)-z_axis)/(x_lcfs(i)-r_axis) !the slope for function Z=Z(X)
!!$     endif
!!$     if(z_lcfs(i)-z_axis .eq. 0._p_) then
!!$        slope2(i)=large_number
!!$     else
!!$        slope2(i)=(x_lcfs(i)-r_axis)/(z_lcfs(i)-z_axis) !the slope for function X=X(Z)
!!$        !write(*,*) i,slope(i),slope2(i)
!!$     endif
!!$  enddo
!!$
!!$  do i=1,np_lcfs
!!$     if(abs(slope(i)).le.1.0_p_) then !use Z=Z(X) function, the reason that I switch between using function X=X(Z) and Z=Z(X) is to aviod large slope.
!!$        x1=r_axis
!!$        x2=x_lcfs(i) !+0.01 !shift left a little to gurrantee that the range is enough for a root to lie in
!!$        x_contour(i)=rtbis(one_dim_psi_func,x1,x2,xacc,r_axis,z_axis,slope(i),psival)
!!$        z_contour(i)=zfunc(r_axis,z_axis,slope(i),x_contour(i))
!!$     else !switch to using X=X(Z) function
!!$        z1=z_axis
!!$        z2=z_lcfs(i)
!!$        z_contour(i)=rtbis(one_dim_psi_func2,z1,z2,xacc,r_axis,z_axis,slope2(i),psival)
!!$        x_contour(i)=xfunc(r_axis,z_axis,slope2(i),z_contour(i)) 
!!$     endif
!!$  enddo

  call contour(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis,psival,x_contour,z_contour)

  call arc_length(x_contour,z_contour,np_lcfs,dl)
!!$     sum=0.
!!$     do i=1,np_lcfs-1
!!$        sum=sum+dl(i)
!!$     enddo
!!$     circumference=sum

  normalization=0._p_
  if(poloidal_angle_type .eq. 'equal-arc') then
     do i=2,np_lcfs !finish a full poloidal circle integration to get the normalization factor
        normalization=normalization+dl(i-1) !equal-arc-length poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'equal-volume') then
     do i=2,np_lcfs !finish a full poloidal circle integration to get the normalization factor
        rmid=0.5_p_*(x_contour(i-1)+x_contour(i))
        zmid=0.5_p_*(z_contour(i-1)+z_contour(i))
        normalization=normalization+dl(i-1)*rmid/psi_gradient_func(rmid,zmid) !equal-volume poloidal angle
     enddo
  else
     stop 'please choose poloidal angle type between equal-arc and equal-volume'
  endif

  call locate_poloidal_index(r,z,x_lcfs,z_lcfs,np_lcfs,end_i) !poloidal index of point (r,z) is between end_i and end_i+1
!!$ if(myid.eq.0) then
!!$     do i=1,end_i
!!$        write(123,*) x_contour(i),z_contour(i)
!!$     enddo
!!$     write(123,*) 
!!$     write(123,*) 
!!$  endif

  x_contour(end_i+1)=r
  z_contour(end_i+1)=z
  dl(end_i)=sqrt((x_contour(end_i+1)-x_contour(end_i))**2+(z_contour(end_i+1)-z_contour(end_i))**2)

  !calculate poloidal angle 
  theta=0._p_

  if(poloidal_angle_type .eq. 'equal-arc') then
     do i=2,end_i+1
        theta=theta+dl(i-1) !equal-arc-length poloidal angle
     enddo

  elseif(poloidal_angle_type .eq. 'equal-volume') then
     do i=2,end_i+1
        rmid=0.5_p_*(x_contour(i-1)+x_contour(i))
        zmid=0.5_p_*(z_contour(i-1)+z_contour(i))
        theta=theta+dl(i-1)*rmid/psi_gradient_func(rmid,zmid) !equal-volume poloidal angle
     enddo
  else
     stop 'please choose poloidal angle type between equal-arc and equal-volume'
  endif

  !theta=theta*twopi/normalization !normalized to the range [0:twopi]
   theta=theta*twopi/normalization-pi !normalized to the range [-pi:pi]

  call cal_toroidal_shift(psival,x_contour,z_contour,np_lcfs,end_i,tor_shift) !calculate toroidal shift which is needed in the definition of the generalized toroidal angle

end subroutine mapping

subroutine mapping2(r,z,radcor,tor_shift_a,tor_shift_b)
  !given (R,Z), this subroutine find the magnetic surface that passes throught the point and calculates its toroidal angle shift
  !the same as mapping(), the difference is that this only takes care of the total_tor_shift at the theta cut, see the comments where this subroutine is called
  use precision,only:p_
  use constants,only:zero,one,two,twopi,pi
  use boundary, only: x_lcfs,z_lcfs,np_lcfs
  use radial_module,only:r_axis,z_axis,psi_axis,psi_lcfs

  implicit none
  real(p_),intent(in):: r,z
  real(p_),intent(out):: radcor,tor_shift_a,tor_shift_b
  real(p_)::psival,psi_func,radcor_as_func_of_pfn
  real(p_):: x_contour(np_lcfs),z_contour(np_lcfs)
  real(p_)::dl(np_lcfs-1)
  integer:: i,end_i

  psival=psi_func(r,z)
  radcor=radcor_as_func_of_pfn((psival-psi_axis)/(psi_lcfs-psi_axis))
  call contour(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis,psival,x_contour,z_contour)

  !  call cal_toroidal_shift_total(psival,x_contour,z_contour,np_lcfs,dl,tor_shift_a) 
  call cal_toroidal_shift_at_theta_cut(psival,x_contour,z_contour,np_lcfs,tor_shift_a,tor_shift_b) !calculate toroidal shift which is needed in the definition of the generalized toroidal angle
end subroutine mapping2


subroutine locate_poloidal_index(r,z,x_lcfs,z_lcfs,np_lcfs,end_i) !poloidal index of point (r,z) is between end_i and end_i+1
  use precision,only:p_
  use constants,only:twopi,pi
  use radial_module,only:r_axis,z_axis
  implicit none
  real(p_),intent(in):: r,z
  integer,intent(in):: np_lcfs
  real(p_),intent(in):: x_lcfs(np_lcfs),z_lcfs(np_lcfs)
  integer,intent(out):: end_i !!poloidal index of point (r,z) is between end_i and end_i+1
  real(p_):: xn,zn,xn0,zn0,angle0,angle1,angle2
  integer:: i

  xn0=r-r_axis
  zn0=z-z_axis
  angle0=acos(xn0/sqrt(xn0*xn0+zn0*zn0))
  !if(zn0<0.0) angle0=twopi-angle0 !theta in [0:twopi]
  if(zn0.le.0.0) angle0=-angle0 !theta in [-pi:pi]


  do i=1,np_lcfs-1
     xn=x_lcfs(i)-r_axis
     zn=z_lcfs(i)-z_axis
     angle1=acos(xn/sqrt(xn**2+zn**2)) !theta in [0:twopi]
     !if(zn<0.0) angle1=twopi-angle1 !theta in [-pi:pi]
     if(zn<0.0) angle1=-angle1
     if(i.eq.1) angle1=-pi

     xn=x_lcfs(i+1)-r_axis
     zn=z_lcfs(i+1)-z_axis
     angle2=acos(xn/sqrt(xn**2+zn**2))
     !if(zn<0.0) angle2=twopi-angle2 !theta in [0:twopi]
     if(zn<0.0) angle2=-angle2 !theta in [-pi:pi]
     !if(i+1.eq.np_lcfs) angle2=twopi !special case, should be twopi, instead of zero. missing this generates a wrong ending-point
      if(i+1.eq.np_lcfs) angle2=pi

     if((angle0-angle1)*(angle0-angle2).le.0._p_) exit

  enddo

  end_i=i

end subroutine locate_poloidal_index

subroutine choose_boundary_magnetic_surfaces_for_the_mapping(r_inner_surf,z_inner_surf,r_outer_surf,z_outer_surf)
  use precision,only:p_
  use constants,only: two
  use magnetic_coordinates,only: pfn_inner,pfn_bdry !boundary for the region in which magnetic_coordinates are constructed
  use radial_module,only: r_axis,z_axis,psi_axis,psi_lcfs
  use boundary,only: x_lcfs,z_lcfs,np_lcfs

  implicit none
  real(p_),intent(out):: r_inner_surf(np_lcfs),z_inner_surf(np_lcfs),r_outer_surf(np_lcfs),z_outer_surf(np_lcfs)
  real(p_):: psi_val,psi_bdry,psi_inner
!  integer:: ierr

  psi_inner=psi_axis+pfn_inner*(psi_lcfs-psi_axis) 
  psi_val=(psi_axis+psi_inner)/two !the inner boundary for doing the mapping
  call contour(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis,psi_val,r_inner_surf,z_inner_surf)

  psi_bdry=psi_axis+pfn_bdry*(psi_lcfs-psi_axis)
  psi_val=(psi_lcfs+psi_bdry)/two !the outer boundary for doing the mapping, which is chosen to be between lcfs and the outer boundary of the region in which magnetic_coordinates are availabe
  call contour(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis,psi_val,r_outer_surf,z_outer_surf)

end subroutine choose_boundary_magnetic_surfaces_for_the_mapping

subroutine create_cylindrical_grids(r_outer_surf,z_outer_surf,np_lcfs,nx,nz,r,z,dr,dz,i0,j0)
!create rectangular box (with cylindrical grids in it) on poloidal plane with the boundary flux surface within the box and (r_axis,z_axis) is exactly on a grid point
  use precision,only:p_
!  use magnetic_coordinates,only:r_mag_surf, z_mag_surf,mpoloidal,nflux
  use radial_module,only: r_axis,z_axis
  implicit none
  integer,intent(in):: np_lcfs,nx,nz
  real(p_),intent(out):: r(nx),z(nz),dr,dz
  integer,intent(out):: i0,j0 !index of the point at magnetic axis
  real(p_):: r_min, r_max,z_min,z_max,r_width,z_width
!  real(p_):: r_mag_surf1(mpoloidal), z_mag_surf1(mpoloidal)
  real(p_):: r_outer_surf(np_lcfs),z_outer_surf(np_lcfs)
  integer:: i,j,nxp,nzp
  real(p_):: rp(nx-1),zp(nz-1)


nxp=nx-1 !using a reduced number, so that I can append the array with an additional element
nzp=nz-1 !using a reduced number, so that I can append the array with an additional element

!!$  do i=1,mpoloidal !select the boundary magnetic surface
!!$     r_mag_surf1(i)=r_mag_surf(i,nflux)
!!$     z_mag_surf1(i)=z_mag_surf(i,nflux)
!!$  enddo

  r_min=minval(r_outer_surf)
  r_max=maxval(r_outer_surf)
  r_width=r_max-r_min

  z_min=minval(z_outer_surf)
  z_max=maxval(z_outer_surf)
  z_width=z_max-z_min

  dr=r_width/(nxp-1)
  dz=z_width/(nzp-1)

!i0=nx/2 !i index at magnetic axis
i0=floor((r_axis-r_min)/dr)+1 !i index near the magnetic axis
j0=floor((z_axis-z_min)/dz)+1 !j index near the magnetic axis
  rp(i0)=r_axis !set (i0,j0) to be on the magnetic axis, since there is a floor operation when getting (i0,j0), this shifts the box to the low-field-side (lfs) and upward
  zp(j0)=z_axis 

 !I do the above steps because I want the magnetic axis to lie exactly on a grid so that I can get grids that exactly represent the midplane.
  do i=i0-1,1,-1 !the grids at the high-field-side (hfs) of the magnetic axis
     rp(i)=rp(i+1)-dr
  enddo
  do i=i0+1,nxp,1 !the girds at the low-field-side (lfs) of the magnetic axis
     rp(i)=rp(i-1)+dr
  enddo

  do j=j0-1,1,-1 !the grids below the midplane
     zp(j)=zp(j+1)-dz
  enddo
  do j=j0+1,nzp,1 !the grids above the midplane
     zp(j)=zp(j-1)+dz
  enddo

!add an additional element. This is needed because setting the magnetic axis to be on a grid usually shifts the box to the lfs, which will make some points on the hfs not included in the box
r(1)=rp(1)-dr
do i=1,nxp
r(i+1)=rp(i)
enddo

z(1)=zp(1)-dz !similar reason as mentioned above
do j=1,nzp
z(j+1)=zp(j)
enddo

i0=i0+1 !the i index of the magnetic axis at the new array r
j0=j0+1 !the j index of the magnetic axis at the new array z

end subroutine create_cylindrical_grids
