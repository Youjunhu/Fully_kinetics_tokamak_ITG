subroutine calculate_metric(mpoloidal,nflux)
  use precision,only:p_
  use magnetic_coordinates,only: r_mag_surf,z_mag_surf,dtheta,dradcor
  use magnetic_coordinates,only: tor_shift_mc,jacobian,abs_jacobian_max,abs_jacobian_min !as output
 use magnetic_coordinates,only: pfn_inner,pfn_bdry,gs_psi_array
  use magnetic_coordinates,only: theta_1d_array,radcor_1d_array
  use radial_module,only: psi_lcfs,psi_axis
  use magnetic_coordinates,only:sign_of_jacobian,sign_of_gs_psi_prime
  use array_in_mc,only: grad_psi_matrix, grad_alpha_matrix,grad_theta_matrix,&
       & grad_psi_dot_grad_alpha_matrix,grad_psi_dot_grad_theta_matrix,grad_alpha_dot_grad_theta_matrix
  use normalizing,only: Ln
  use domain_decomposition,only: myid
  implicit none
  integer,intent(in):: mpoloidal,nflux
  real(p_):: g_value,g_func,q_func_pfn !function name
  integer:: i,j,file_unit
  character(8):: filename

  !  allocate(jacobian(mpoloidal,nflux))
  call calculate_metric1(mpoloidal,nflux,r_mag_surf,z_mag_surf,dtheta,dradcor) !computing gradients of psi and theta in magnetic coordinate

  call plasma_volume(nflux,mpoloidal,jacobian,dradcor,dtheta,myid) !calculate volume element
  if(myid.eq.0) write(*,*) 'max_jacobian=',maxval(jacobian),'min_jacobian=',minval(jacobian)
  abs_jacobian_max=maxval(abs(jacobian))
  abs_jacobian_min=minval(abs(jacobian))
  !write(*,*) 'max_abs(jacobian)=',maxval(abs(jacobian)),'min_abs(jacobian)=',minval(abs(jacobian))

  sign_of_jacobian=sign(1._p_,jacobian(mpoloidal/3,nflux/3))
  sign_of_gs_psi_prime=sign(1._p_,psi_lcfs-psi_axis) !radial coordinator is assumed to be pfn
  call cal_q_with_correct_sign() ! after this call, function "q_func_pfn" is ready to be used
  if(myid.eq.0) write(*,*) 'q_inner=',q_func_pfn(pfn_inner),'q_bdry=',q_func_pfn(pfn_bdry)

  allocate(tor_shift_mc(mpoloidal,nflux))
  do j=1,nflux
     g_value=g_func(gs_psi_array(j))
     call cal_toroidal_shift2(g_value,r_mag_surf(:,j),z_mag_surf(:,j),mpoloidal,tor_shift_mc(:,j))
  enddo
if(myid.eq.0) write(*,*) 'maximum of tor_shift_mc=',maxval(tor_shift_mc),'minimum of tor_shift_mc=',minval(tor_shift_mc)
  if(myid.eq.0) then
     open(113,file='theta_tor_shift2d.txt')
     do j=1,nflux
        do i=1,mpoloidal
           write(113,*) r_mag_surf(i,j),z_mag_surf(i,j),theta_1d_array(i),tor_shift_mc(i,j)
        enddo
        write(113,*)
        write(113,*)
     enddo
  endif

!!$ if(myid.eq.0) call draw_grids_on_theta_isosurface(mpoloidal,nflux,tor_shift_mc,r_mag_surf,z_mag_surf)
!!$ if(myid.eq.0) call draw_alpha_isosurface(mpoloidal,nflux,tor_shift_mc,r_mag_surf,z_mag_surf)
!!$ if(myid.eq.0) call draw_alpha_contours_on_a_magnetic_surface(mpoloidal,nflux,tor_shift_mc,r_mag_surf,z_mag_surf)

  call calculate_metric2(mpoloidal,nflux,r_mag_surf,z_mag_surf,tor_shift_mc,dtheta,dradcor) !computing gradient of generalized toroidal angle alpha


  if(myid.eq.0) then
     write(filename,'(a4,i4.4)') 'compa',myid
     file_unit=myid+311
     open(file_unit,file=filename)
     do i=1,mpoloidal
        do j=1,nflux
           write(file_unit,'(2i8.4,9(1pe14.5))')  i,j,grad_alpha_matrix(i,j),grad_psi_dot_grad_alpha_matrix(i,j),&
                & grad_psi_matrix(i,j),jacobian(i,j),grad_alpha_dot_grad_theta_matrix(i,j),&
                & tor_shift_mc(i,j),grad_psi_dot_grad_theta_matrix(i,j),theta_1d_array(i),radcor_1d_array(j)
        enddo
        write(file_unit,*)
     enddo
     close(file_unit)
  endif

  !some tests
!!$  if(myid.eq.0) then
!!$     do j=1,nflux
!!$        write(*,*) j,radcor_1d_array(j),minor_r_radcor(radcor_1d_array(j)),r_mag_surf0(1,j)-r_axis,&
!!$             & minor_r_prime(radcor_1d_array(j))
!!$     enddo
!!$  endif

end subroutine calculate_metric

subroutine calculate_metric1(m,n,r,z,dtheta,dpsi)
  use precision,only:p_
  use constants,only:zero,one,two
  use magnetic_coordinates,only: jacobian !as output
  use array_in_mc,only: rpsi,zpsi,rth,zth,grad_psi_matrix,grad_theta_matrix,grad_psi_dot_grad_theta_matrix !as output
  use array_in_mc,only:grad_psi_r_matrix,grad_psi_z_matrix,grad_theta_r_matrix,grad_theta_z_matrix !as output
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: r(m,n),z(m,n)
  real(p_),intent(in):: dtheta, dpsi
  integer:: i,j,u
  real(p_):: minor_r_radcor,minor_r_val
  real(p_)::minor_r_prime, minor_r_prime_val

  allocate(jacobian(m,n))
  allocate(grad_psi_matrix(m,n))
  allocate(grad_theta_matrix(m,n))
  allocate(grad_psi_dot_grad_theta_matrix(m,n))
  allocate(grad_psi_r_matrix(m,n))
  allocate(grad_psi_z_matrix(m,n))

  allocate(rpsi(m,n))
  allocate(zpsi(m,n))
  allocate(rth(m,n))
  allocate(zth(m,n))

  call partial_derivative_in_mc(m,n,r,z,dtheta,dpsi,rpsi,rth,zpsi,zth,jacobian) 

  grad_psi_matrix=abs(r/jacobian)*sqrt(zth**2+rth**2)
  grad_theta_matrix=abs(r/jacobian)*sqrt(zpsi**2+rpsi**2)
  grad_psi_dot_grad_theta_matrix=-(r/jacobian)**2*(zth*zpsi+rth*rpsi)

  grad_psi_r_matrix=-r/jacobian*zth
  grad_psi_z_matrix=r/jacobian*rth
  grad_theta_r_matrix=r/jacobian*zpsi
  grad_theta_z_matrix=-r/jacobian*rpsi
  !call plot_psi_r_z_theta_r_z_mc()
contains
  subroutine plot_psi_r_z_theta_r_z_mc()
    use magnetic_coordinates,only: theta_1d_array,radcor_1d_array
    use domain_decomposition,only: myid
    integer:: u
    real(p_):: minor_r_radcor,minor_r_val
    real(p_)::minor_r_prime, minor_r_prime_val
    if(myid.eq.0) then
       open(newunit=u, file='psi_r_z_theta_r_z_mc')
       do i=1,m
          do j=1,n
             !write(u,*) r(i,j),z(i,j),grad_theta_matrix(i,j),jacobian(i,j)
             minor_r_prime_val=minor_r_prime(radcor_1d_array(j))
             write(u,*) r(i,j),z(i,j),grad_psi_r_matrix(i,j)*minor_r_prime_val, grad_psi_z_matrix(i,j)*minor_r_prime_val,&
                  & grad_theta_r_matrix(i,j),grad_theta_z_matrix(i,j)
          enddo
          write(u,*)
       enddo
       close(u)

       open(newunit=u, file='psi_r_z_theta_r_z_mc_analytic') !for concentric-circular magnetic field (equal-arc-length theta and geometric minor radius as the flux surface label)
       do i=1,m
          do j=1,n
             minor_r_val=minor_r_radcor(radcor_1d_array(j))
             write(u,*) r(i,j),z(i,j),cos(theta_1d_array(i)),sin(theta_1d_array(i)),&
                  &-sin(theta_1d_array(i))/minor_r_val,cos(theta_1d_array(i))/minor_r_val
          enddo
          write(u,*)
       enddo
       close(u)
    endif
  end subroutine plot_psi_r_z_theta_r_z_mc
end subroutine calculate_metric1


subroutine calculate_metric2(m,n,r,z,tor_shift_mc,dtheta,dpsi) !compute the gradient of the generalized toroidal angle alpha
  use precision,only:p_
  use constants,only:zero,one,two
  use magnetic_coordinates,only: jacobian !as input
  use array_in_mc,only: rpsi,zpsi,rth,zth !as input
  use array_in_mc,only: grad_alpha_r_matrix,grad_alpha_z_matrix,grad_alpha_phi_matrix !as output 
  use array_in_mc,only: grad_alpha_matrix,grad_psi_dot_grad_alpha_matrix,grad_alpha_dot_grad_theta_matrix !as output
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: r(m,n),z(m,n),tor_shift_mc(m,n)
  real(p_),intent(in):: dtheta, dpsi
  real(p_):: tor_shift_psi(m,n), tor_shift_th(m,n) !temporary array
  integer:: i,j

  allocate(grad_alpha_r_matrix(m,n))
  allocate(grad_alpha_z_matrix(m,n))
  allocate(grad_alpha_phi_matrix(m,n))
  allocate(grad_alpha_matrix(m,n))
  allocate(grad_psi_dot_grad_alpha_matrix(m,n))
  allocate(grad_alpha_dot_grad_theta_matrix(m,n))

  call partial_derivative_of_tor_shift_in_mc(m,n,tor_shift_mc,dtheta,dpsi,tor_shift_psi,tor_shift_th) 
  grad_alpha_r_matrix=tor_shift_psi*r/jacobian*zth-tor_shift_th*r/jacobian*zpsi
  grad_alpha_z_matrix=tor_shift_th*r/jacobian*rpsi-tor_shift_psi*r/jacobian*rth
  grad_alpha_phi_matrix=one/r
  grad_alpha_matrix=sqrt(grad_alpha_phi_matrix**2+grad_alpha_r_matrix**2+grad_alpha_z_matrix**2)
  grad_psi_dot_grad_alpha_matrix=-r/jacobian*zth*grad_alpha_r_matrix+r/jacobian*rth*grad_alpha_z_matrix
  grad_alpha_dot_grad_theta_matrix=grad_alpha_r_matrix*r/jacobian*zpsi-grad_alpha_z_matrix*r/jacobian*rpsi
  !call plot_alpha_r_z_mc()
contains
  subroutine plot_alpha_r_z_mc()
    use magnetic_coordinates,only: radcor_1d_array,theta_1d_array
    use domain_decomposition,only: myid
    integer:: u
    real(p_):: theta,minor_r,major_r0,major_r,q,local_q,qprime,factor
    real(p_):: ddelta_minor_r,alpha_r,alpha_z
    real(p_):: q_func_pfn, minor_r_radcor

    if(myid.eq.0) then
       open(newunit=u,file='alpha_r_z_mc')
       do i=1,m
          do j=1,n
             write(u,*) r(i,j),z(i,j),grad_alpha_r_matrix(i,j),grad_alpha_z_matrix(i,j) !,grad_alpha_matrix(i,j)
          enddo
          write(u,*)
       enddo
       close(u)
       !write(*,*) 'maxval(grad_alpha_matrix)=',maxval(grad_alpha_matrix),'minval(grad_alpha_matrix)=',minval(grad_alpha_matrix)
       open(newunit=u,file='alpha_r_z_mc_analytic') !for concentric-circular magnetic field (equal-arc-length theta and geometric minor radius as the flux surface label)
       do i=1,m
          theta=theta_1d_array(i)
          do j=1,n
             major_r0=1.32_p_!for DIII-D cyclone base case parameter
             minor_r=minor_r_radcor(radcor_1d_array(j))
             major_r=major_r0+minor_r*cos(theta)
             q=q_func_pfn(radcor_1d_array(j))
             local_q=q*sqrt(major_r0**2-minor_r**2)/major_r
             qprime=0.78*1.4/0.24_p_ !for DIII-D cyclone base case parameter, a linear q profile is assumed, the same as what Ben did
             factor=(major_r0-minor_r)/sqrt(major_r0**2-minor_r**2)*tan(theta/two)
             ddelta_minor_r=two*qprime*atan(factor)+two*q/(1+factor**2)*tan(theta/2._p_)*&
                  & (-major_r0)/((major_r0+minor_r)*sqrt(major_r0**2-minor_r**2))
             alpha_r=-ddelta_minor_r*cos(theta)+local_q*sin(theta)/minor_r
             alpha_z=-local_q*cos(theta)/minor_r-ddelta_minor_r*sin(theta)
             write(u,*) r(i,j),z(i,j), alpha_r,alpha_z
          enddo
          write(u,*)
       enddo
       close(u)

    endif
  end subroutine plot_alpha_r_z_mc
end subroutine calculate_metric2


subroutine partial_derivative_in_mc(m,n,r,z,dtheta,dpsi,rpsi,rth,zpsi,zth,jacob)
  !calculate the partial derivative of R and Z with respect to the magnetic cooordinates (psi,theta)
  !jacob is also calculated in this subroutine
  use precision,only:p_
  use constants,only:zero,one,two,twopi,one_half
  implicit none

  integer,intent(in):: m,n
  real(p_),intent(in):: r(m,n),z(m,n)
  real(p_),intent(in):: dtheta,dpsi
  real(p_),intent(out):: rpsi(m,n),rth(m,n),zpsi(m,n),zth(m,n),jacob(m,n)
  real(p_):: tmp0
  integer:: i,j

  do i=1,m  
     do j=2,n-1 !use center difference scheme for inner points
        rpsi(i,j)=(r(i,j+1)-r(i,j-1))/(two*dpsi)
        zpsi(i,j)=(z(i,j+1)-z(i,j-1))/(two*dpsi)
     enddo

     !use linear interpolation to get the value  j=n
     tmp0=(r(i,n)-r(i,n-1))/dpsi
     rpsi(i,n)=two*tmp0-rpsi(i,n-1)
     tmp0=(z(i,n)-z(i,n-1))/dpsi
     zpsi(i,n)=two*tmp0-zpsi(i,n-1)

     !use linear interpolation to get the value j=1
     tmp0=(r(i,2)-r(i,1))/dpsi
     rpsi(i,1)=two*tmp0-rpsi(i,2)

     tmp0=(z(i,2)-z(i,1))/dpsi
     zpsi(i,1)=two*tmp0-zpsi(i,2)

  enddo

  do j=1,n
     do i=2,m-1 !use center difference scheme for inner points
        rth(i,j)= (r(i+1,j)-r(i-1,j))/(two*dtheta)
        zth(i,j)=(z(i+1,j)-z(i-1,j))/(two*dtheta)
     enddo

     !use peroidic property of r and z to calculate the partial derivative for boundary points at theta=0 and 2pi
     rth(1,j)=(r(2,j)-r(m-1,j))/(two*dtheta)
     zth(1,j)=(z(2,j)-z(m-1,j))/(two*dtheta)
     rth(m,j)=rth(1,j)
     zth(m,j)=zth(1,j)
  enddo

  !calculate the Jacobian:
  do i=1,m
!          do j=2,n !the jacobian at the magnetic axis is zero
     do j=1,n
        jacob(i,j)=r(i,j)*(rth(i,j)*zpsi(i,j)-rpsi(i,j)*zth(i,j))   !Jacobain of coordinate system (psi,theta,fai)
     enddo
     !jacob(i,n)=two*jacob(i,n-1)-jacob(i,n-2)
     !use linear interpolation to get the value of jacobian at the magnetic surface near the magnetic axis
  enddo

  !write(*,*) jacob(1:m,5)
end subroutine partial_derivative_in_mc


subroutine partial_derivative_of_tor_shift_in_mc(m,n,tor_shift,dtheta,dpsi,tor_shift_psi,tor_shift_th)
  !calculate the partial derivative with respect to the magnetic cooordinates (psi,theta)
  use precision,only:p_
  use constants,only:zero,one,two,twopi,one_half
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: tor_shift(m,n)
  real(p_),intent(in):: dtheta,dpsi
  real(p_),intent(out):: tor_shift_psi(m,n),tor_shift_th(m,n)
  real(p_):: tmp0,twopi_q
  integer:: i,j

  do i=1,m  
     do j=2,n-1 !use center difference scheme for inner points
        tor_shift_psi(i,j)=(tor_shift(i,j+1)-tor_shift(i,j-1))/(two*dpsi)
     enddo
     !use linear interpolation to get the value at j=n
     tmp0=(tor_shift(i,n)-tor_shift(i,n-1))/dpsi
     tor_shift_psi(i,n)=two*tmp0-tor_shift_psi(i,n-1)
     !use linear interpolation to get the value j=1
     tmp0=(tor_shift(i,2)-tor_shift(i,1))/dpsi
     tor_shift_psi(i,1)=two*tmp0-tor_shift_psi(i,2)
  enddo
  do j=1,n
     do i=2,m-1 !use center difference scheme for inner points
        tor_shift_th(i,j)= (tor_shift(i+1,j)-tor_shift(i-1,j))/(two*dtheta)
     enddo
     !for boundary points at theta cut
!!$     twopi_q=tor_shift(m,j)
!!$     tor_shift_th(1,j)=(tor_shift(2,j)+twopi_q-tor_shift(m-1,j))/(two*dtheta)
!!$     tor_shift_th(m,j)=tor_shift_th(1,j)

   !use linear interpolation to get the value at i=1 and i=m
     tor_shift_th(1,j)=two*tor_shift_th(2,j)-tor_shift_th(3,j)
     tor_shift_th(m,j)=two*tor_shift_th(m-1,j)-tor_shift_th(m-2,j)
  enddo
end subroutine partial_derivative_of_tor_shift_in_mc



subroutine cal_q_with_correct_sign()
  use radial_module,only:npsi,qpsi
  use magnetic_coordinates,only:sign_of_jacobian,sign_of_gs_psi_prime
  use radial_module,only: sign_bphi
  use radial_module,only: q_with_sign !as output
  implicit none

  allocate(q_with_sign(npsi))

  q_with_sign=abs(qpsi)*(-sign_bphi)*sign_of_jacobian*sign_of_gs_psi_prime !include the correct sign

end subroutine cal_q_with_correct_sign




subroutine plasma_volume(nflux,mpoloidal,jacobian,dradcor,dtheta,myid)
  !to calculate the volume within the LCFS
  use precision,only: p_
  use constants,only: two,twopi
  use magnetic_coordinates,only: dv,vol0 !as output
  implicit none
  integer,intent(in)::nflux,mpoloidal,myid
  real(p_),intent(in):: jacobian(mpoloidal,nflux),dradcor,dtheta
  real(p_):: vol
  integer:: i,j

  allocate(dv(mpoloidal,nflux))

  do i=1,mpoloidal
     do j=1,nflux
        dv(i,j)=abs(jacobian(i,j))*dradcor*dtheta*twopi
     enddo
  enddo

  vol=0._p_
  do i=1,mpoloidal-1
     do j=1,nflux-1
        vol=vol+dv(i,j)
     enddo
  enddo
  vol0=vol


  if(myid.eq.0) write(*,*) 'Plasma volume within boundary flux surface is (m^3)',vol0

end subroutine plasma_volume




function jacobian_func(theta,radcor) result (z) 
  use precision,only: p_
  use magnetic_coordinates,only: mpoloidal,nflux,theta_1d_array,radcor_1d_array,jacobian !as input
  use interpolate_module,only: linear_2d_interpolation

  implicit none
  real(p_)::radcor,theta,z
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,jacobian,theta,radcor,z)  !uniform 1darray is assumed
end function jacobian_func

function jacobian_func_normalized(theta,radcor) result (z) !!used in generating non-uniformly distributed random numbers that satisfied a probability density function proportional to the |jacobian|
  use precision,only: p_
  use magnetic_coordinates,only: mpoloidal,nflux,theta_1d_array,radcor_1d_array,jacobian !as input
  use magnetic_coordinates,only: nflux1,j_low1
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  real(p_)::radcor,theta,z
  real(p_):: jacobian0(mpoloidal,nflux)
  real(p_):: jac1(mpoloidal,nflux1),abs_jacobian_min,abs_jacobian_max
  integer:: i,j,jshift

  jacobian0=abs(jacobian)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,jacobian0,theta,radcor,z)  !uniform 1darray is assumed

  do i=1,mpoloidal
     do j=1,nflux1
        jshift=j-1+j_low1
        jac1(i,j)=jacobian0(i,jshift)
     enddo
  enddo
  abs_jacobian_max=maxval(jac1)
!  abs_jacobian_min=minval(jac1)

  !  z=(z-abs_jacobian_min)/(abs_jacobian_max-abs_jacobian_min) !turns out to be wrong! a subtle bug which I spent much time in finding, the shift is wrong. This bug is due to my misunderstanding about the rejection method, not due to programming mistakes.
    z=z/abs_jacobian_max
end function jacobian_func_normalized
