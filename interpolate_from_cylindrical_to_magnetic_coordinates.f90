subroutine interpolate_from_cylindrical_to_magnetic_coordinates(r0,z0,theta0,tor_shift0) !calculate the magnetic coordinats (theta0,tor_shift0) of (R0,Z0) by using interpolation, radcor0 is not calculated in this subroutine because we can use my reliable way to calculate radcor0 from (r0,z0), i.e., radcor0=radcor_as_func_of_pfn(pfn_func(r0,z0))
!to use this iterpolation, first to make sure that (r0,z0) is within the specfied region
use constants,only:p_,two,twopi
use mapping_module,only: r_cyl,z_cyl,dr,dz,radcor,theta_a,theta_b,tor_shift_a,i0,j0,tor_shift_b
use interpolate_module,only: linear_2d_interpolation_kernel
implicit none
real(p_),intent(in):: r0,z0
real(p_),intent(out):: theta0,tor_shift0
real(p_):: radcor_tmp(2,2),theta_tmp(2,2),tor_shift_tmp(2,2),q_func2,qval
integer::i,j,ii,jj
real(p_):: radcor_as_func_of_pfn,pfn_func
!locate
i=(r0-r_cyl(1))/dr+1
j=(z0-z_cyl(1))/dz+1

!  2D interpolations to get radcor0, theta0, and tor_shift0
!!$  do ii=1,2
!!$     do jj=1,2
!!$        radcor_tmp(ii,jj)=radcor(i+ii-1,j+jj-1)
!!$     enddo
!!$  enddo
!!$
!!$  call linear_2d_interpolation_kernel(r_cyl(i),z_cyl(j),radcor_tmp,r0,z0,radcor0)

!  radcor0=radcor_as_func_of_pfn(pfn_func(r0,z0))


  do ii=1,2
     do jj=1,2
        theta_tmp(ii,jj)=theta_a(i+ii-1,j+jj-1)
     enddo
  enddo

!!$ if(i>i0 .and. j+1.eq.j0) then !handle the boundary at the low-field-side midplane
!!$theta_tmp(1,2)=twopi 
!!$theta_tmp(2,2)=twopi 
!!$endif
 if(i<i0 .and. j.eq.j0) then !handle the boundary at the high-field-side (theta cut)
     theta_tmp(1,1)=theta_b(i,j0)
     theta_tmp(2,1)=theta_b(i+1,j0)
  endif
  call linear_2d_interpolation_kernel(r_cyl(i),z_cyl(j),theta_tmp,r0,z0,theta0)

  do ii=1,2
     do jj=1,2
        tor_shift_tmp(ii,jj)=tor_shift_a(i+ii-1,j+jj-1)
     enddo
  enddo

!!$if(i>i0 .and. j+1.eq.j0) then !handle the boundary at the low-field-side midplne
!!$!qval=q_func2(radcor0)
!!$!qval=2.3074808101594884*1.0004
!!$!tor_shift_tmp(1,2)=two*tor_shift(i,j)-tor_shift(i,j-1) 
!!$!tor_shift_tmp(1,2)=qval*twopi
!!$!tor_shift_tmp(2,2)=qval*twopi
!!$tor_shift_tmp(1,2)=tor_shift_b(i,j0)
!!$tor_shift_tmp(2,2)=tor_shift_b(i+1,j0)
!!$endif

 if(i<i0 .and. j.eq.j0) then !handle the boundary at the high-field-side midplane (theta cut)
     tor_shift_tmp(1,1)=tor_shift_b(i,j0)
     tor_shift_tmp(2,1)=tor_shift_b(i+1,j0)
  endif

  call linear_2d_interpolation_kernel(r_cyl(i),z_cyl(j),tor_shift_tmp,r0,z0,tor_shift0)

!alpha0=phi0-tor_shift0

end subroutine interpolate_from_cylindrical_to_magnetic_coordinates


subroutine interpolate_from_cylindrical_to_magnetic_coordinates1(r0,z0,theta0) !calculate the poloidal angle theta0 of (R0,Z0) by using interpolation, radcor0 is not calculated in this subroutine because we can use my reliable way to calculate radcor0 from (r0,z0), i.e., radcor0=radcor_as_func_of_pfn(pfn_func(r0,z0))
!to use this iterpolation, first to make sure that (r0,z0) is within the specfied region
use constants,only:p_,two,twopi
use mapping_module,only: r_cyl,z_cyl,dr,dz,radcor,theta_a,theta_b,i0,j0
  use interpolate_module

implicit none
real(p_),intent(in):: r0,z0
real(p_),intent(out):: theta0
real(p_):: radcor_tmp(2,2),theta_tmp(2,2),tor_shift_tmp(2,2),q_func2,qval
integer::i,j,ii,jj
real(p_):: radcor_as_func_of_pfn,pfn_func
!locate
i=(r0-r_cyl(1))/dr+1
j=(z0-z_cyl(1))/dz+1

!  2D interpolations to get radcor0, theta0, and tor_shift0
!!$  do ii=1,2
!!$     do jj=1,2
!!$        radcor_tmp(ii,jj)=radcor(i+ii-1,j+jj-1)
!!$     enddo
!!$  enddo
!!$
!!$  call linear_2d_interpolation_kernel(r_cyl(i),z_cyl(j),radcor_tmp,r0,z0,radcor0)

!  radcor0=radcor_as_func_of_pfn(pfn_func(r0,z0))

  do ii=1,2
     do jj=1,2
        theta_tmp(ii,jj)=theta_a(i+ii-1,j+jj-1)
     enddo
  enddo

 if(i<i0 .and. j.eq.j0) then !handle the boundary at the high-field-side (theta cut)
     theta_tmp(1,1)=theta_b(i,j0)
     theta_tmp(2,1)=theta_b(i+1,j0)
  endif
  call linear_2d_interpolation_kernel(r_cyl(i),z_cyl(j),theta_tmp,r0,z0,theta0)

!!$  do ii=1,2
!!$     do jj=1,2
!!$        tor_shift_tmp(ii,jj)=tor_shift_a(i+ii-1,j+jj-1)
!!$     enddo
!!$  enddo
!!$
!!$ if(i<i0 .and. j.eq.j0) then !handle the boundary at the high-field-side midplane (theta cut)
!!$     tor_shift_tmp(1,1)=tor_shift_b(i,j0)
!!$     tor_shift_tmp(2,1)=tor_shift_b(i+1,j0)
!!$  endif
!!$
!!$  call linear_2d_interpolation_kernel(r_cyl(i),z_cyl(j),tor_shift_tmp,r0,z0,tor_shift0)

end subroutine interpolate_from_cylindrical_to_magnetic_coordinates1
