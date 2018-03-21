subroutine cal_toroidal_shift_at_theta_cut(psival,x_contour,z_contour,np_lcfs,tor_shift_a,tor_shift_b)
  use precision,only:p_
  use constants,only:one_half
  ! use poloidal_flux_2d,only: nx,nz,xarray,zarray
  use magnetic_coordinates,only:sign_of_jacobian,sign_of_gs_psi_prime
  use radial_module,only:psi_axis,psi_lcfs
  use domain_decomposition,only:myid
  implicit none
  real(p_),intent(in):: psival
  integer,intent(in):: np_lcfs
  real(p_),intent(in):: x_contour(np_lcfs),z_contour(np_lcfs)
  real(p_),intent(out):: tor_shift_a,tor_shift_b
  real(p_):: g_value,g_func !toroidal field function
  real(p_):: psi_z_func,psi_r_func,q_func_pfn
  real(p_):: x_mid,z_mid,grad_psi,dl
  integer:: i,i_theta_zero

  g_value=g_func(psival)
  i_theta_zero=(np_lcfs+1)/2

  tor_shift_a=0._p_
  do i=i_theta_zero-1,1,-1 ! for location near the cut below the midplane
     x_mid=(x_contour(i)+x_contour(i+1))*one_half
     z_mid=(z_contour(i)+z_contour(i+1))*one_half
     grad_psi=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
     dl=-sqrt((x_contour(i)-x_contour(i+1))**2+(z_contour(i)-z_contour(i+1))**2)
     tor_shift_a=tor_shift_a+g_value/(x_mid*grad_psi)*dl
  enddo

  tor_shift_b=0._p_
  do i=i_theta_zero+1,np_lcfs ! for location near the cut above the midplane
     x_mid=(x_contour(i)+x_contour(i-1))*one_half
     z_mid=(z_contour(i)+z_contour(i-1))*one_half
     grad_psi=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
     dl=sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
     tor_shift_b=tor_shift_b+g_value/(x_mid*grad_psi)*dl
  enddo

  tor_shift_a=tor_shift_a*(-sign_of_jacobian)/sign_of_gs_psi_prime !include the correct sign
  tor_shift_b=tor_shift_b*(-sign_of_jacobian)/sign_of_gs_psi_prime !include the correct sign
end subroutine cal_toroidal_shift_at_theta_cut

subroutine cal_toroidal_shift(psival,x_contour,z_contour,np_lcfs,end_i,tor_shift)
  use precision,only:p_
  use constants,only:one_half
  ! use poloidal_flux_2d,only: nx,nz,xarray,zarray
  use magnetic_coordinates,only:sign_of_jacobian,sign_of_gs_psi_prime
  use radial_module,only:psi_axis,psi_lcfs
  use domain_decomposition,only:myid
  implicit none
  real(p_),intent(in):: psival
  integer,intent(in):: np_lcfs,end_i
  !real(p_),intent(in):: x_contour(end_i+1),z_contour(end_i+1),dl(end_i)
  real(p_),intent(in):: x_contour(np_lcfs),z_contour(np_lcfs) !,dl(np_lcfs-1)
  real(p_),intent(out):: tor_shift
  real(p_):: g_value,g_func !toroidal field function
  real(p_):: psi_z_func,psi_r_func,q_func_pfn
  !real(p_):: grad_psi(end_i+1)
  real(p_):: x_mid,z_mid,grad_psi,dl
  real(p_):: pfn,mr,costh,sinth,r0=1.32_p_
  integer:: i,i_theta_zero

  g_value=g_func(psival)

!!$  do i=1,end_i
!!$     x_mid=(x_contour(i)+x_contour(i+1))*one_half
!!$     z_mid=(z_contour(i)+z_contour(i+1))*one_half
!!$     grad_psi(i)=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
!!$  enddo

  i_theta_zero=(np_lcfs+1)/2
  tor_shift=0._p_
  if(end_i.ge.i_theta_zero) then
     do i=i_theta_zero+1,end_i+1
        x_mid=(x_contour(i)+x_contour(i-1))*one_half
        z_mid=(z_contour(i)+z_contour(i-1))*one_half
        grad_psi=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
        dl=sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
        tor_shift=tor_shift+g_value/(x_mid*grad_psi)*dl
     enddo

  else if (end_i.lt.i_theta_zero) then 
     do i=i_theta_zero-1,end_i+1,-1
        x_mid=(x_contour(i)+x_contour(i+1))*one_half
        z_mid=(z_contour(i)+z_contour(i+1))*one_half
        grad_psi=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
        dl=-sqrt((x_contour(i)-x_contour(i+1))**2+(z_contour(i)-z_contour(i+1))**2)
        tor_shift=tor_shift+g_value/(x_mid*grad_psi)*dl
     enddo
  endif
  !--for testing----analytical formula, for concentric circular configuration
!!$  mr=sqrt((x_contour(end_i+1)-r0)**2+(z_contour(end_i+1)-0._p_)**2)
!!$  costh=(x_contour(end_i+1)-r0)/mr
!!$  sinth=z_contour(end_i+1)/mr
!!$  pfn=(psival-psi_axis)/(psi_lcfs-psi_axis)
!!$  tor_shift=2*q_func_pfn(pfn)*atan((r0-mr)/sqrt(r0**2-mr**2)*sinth/(costh+1))
  !----for tesing----------------------

!!$  do i=2,end_i+1
!!$     x_mid=(x_contour(i-1)+x_contour(i))*one_half
!!$     tor_shift=tor_shift+g_value/(x_mid*grad_psi(i-1))*dl(i-1)
!!$  enddo
  tor_shift=tor_shift*(-sign_of_jacobian)/sign_of_gs_psi_prime !include the correct sign
end subroutine cal_toroidal_shift


subroutine cal_toroidal_shift2(g_value,x_contour,z_contour,m,tor_shift) !cf. cal_toroidal_shift, the differences are (1) g_value instead of psival is provided as input; (2) compute a series of tor_shift on a magnetic surface. Dec-5-2017: change the range of theta to [-pi:pi]
  use precision,only:p_
  use constants,only:one_half
  use magnetic_coordinates,only:sign_of_jacobian,sign_of_gs_psi_prime,i_theta_zero
  implicit none
  real(p_),intent(in):: g_value
  integer,intent(in):: m
  real(p_),intent(in):: x_contour(m),z_contour(m)
  real(p_),intent(out):: tor_shift(m)
  real(p_):: psi_r_func,psi_z_func,x_mid,z_mid
  real(p_):: grad_psi,dl
  integer:: i


  tor_shift(i_theta_zero)=0._p_

  do i=i_theta_zero+1,m !for theta in [0:pi]
     x_mid=(x_contour(i-1)+x_contour(i))*one_half
     z_mid=(z_contour(i-1)+z_contour(i))*one_half
     grad_psi=sqrt((psi_r_func(x_mid,z_mid))**2+(psi_z_func(x_mid,z_mid))**2) !can also use the value of grad_psi calculated in magnetic coordinates
     dl=sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
     tor_shift(i)=tor_shift(i-1)+g_value/(x_mid*grad_psi)*dl
  enddo

  do i=i_theta_zero-1,1,-1 !for theta in [0:-pi]
     x_mid=(x_contour(i+1)+x_contour(i))*one_half
     z_mid=(z_contour(i+1)+z_contour(i))*one_half
     grad_psi=sqrt((psi_r_func(x_mid,z_mid))**2+(psi_z_func(x_mid,z_mid))**2) !can also use the value of grad_psi calculated in magnetic coordinates
     dl=-sqrt((x_contour(i)-x_contour(i+1))**2+(z_contour(i)-z_contour(i+1))**2)
     tor_shift(i)=tor_shift(i+1)+g_value/(x_mid*grad_psi)*dl
  enddo

!another way of calculating tor_shift
!!$ tor_shift(1)=0._p_
!!$  do i=1,m-1 
!!$     x_mid=(x_contour(i+1)+x_contour(i))*one_half
!!$     z_mid=(z_contour(i+1)+z_contour(i))*one_half
!!$     grad_psi=sqrt((psi_r_func(x_mid,z_mid))**2+(psi_z_func(x_mid,z_mid))**2)
!!$     dl=sqrt((x_contour(i+1)-x_contour(i))**2+(z_contour(i+1)-z_contour(i))**2)
!!$     tor_shift(i+1)=tor_shift(i)+g_value/(x_mid*grad_psi)*dl
!!$  enddo
!!$  tor_shift=tor_shift-tor_shift(i_theta_zero)

!!$  tor_shift(1)=0._p_
!!$  do i=2,m
!!$     x_mid=(x_contour(i-1)+x_contour(i))*one_half
!!$     z_mid=(z_contour(i-1)+z_contour(i))*one_half
!!$     grad_psi=sqrt((psi_r_func(x_mid,z_mid))**2+(psi_z_func(x_mid,z_mid))**2)
!!$     dl=sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
!!$     tor_shift(i)=tor_shift(i-1)+g_value/(x_mid*grad_psi)*dl
!!$  enddo

  tor_shift=tor_shift*(-sign_of_jacobian)/sign_of_gs_psi_prime !include the correct sign
end subroutine cal_toroidal_shift2


