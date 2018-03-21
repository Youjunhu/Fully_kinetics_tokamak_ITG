subroutine adiabatic_electron_field_solver()
  use mpi
  use constants,only: one,two,one_half
  use precision,only:p_
  use magnetic_coordinates,only: mtoroidal,nflux2,dradcor,dtor,j_low2,theta_1d_array,dtheta,r_mag_surf
  use magnetic_coordinates,only: radcor_1d_array2,tor_1d_array
  use perturbation_field_matrix,only: den_left,potential
  use domain_decomposition,only: theta_start
  use array_in_mc,only: grad_psi_r_matrix,grad_psi_z_matrix,grad_theta_r_matrix,grad_theta_z_matrix
  use array_in_mc,only: grad_alpha_r_matrix,grad_alpha_z_matrix,grad_alpha_phi_matrix
  use perturbation_field_matrix,only: ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left !as output
  use perturbation_field_matrix,only: filter_toroidal,filter_radial
  use electrons_module,only: charge_e,te0,ne0
  use constants,only:kev
  use normalizing,only:bn,ln,vn_i
  use domain_decomposition,only:myid
  use delta_ne_module,only: delta_ne, delta_ne_theta,delta_ne_psi,delta_ne_alpha !for testing
  use transform_module,only:oned_fourier_transform1, oned_backward_fourier_transform1,&
       &oned_sine_transform2,oned_inverse_sine_transform2
  use filter_module,only: toroidal_filter,radial_sine_filter
  use derivatives_in_field_line_following_coordinates,only: radial_derivative,toroidal_derivative,theta_derivative
  implicit none

  real(p_):: potential_psi(mtoroidal,nflux2),potential_alpha(mtoroidal,nflux2),potential_theta(mtoroidal,nflux2)
  real(p_):: ef_cyl_r(mtoroidal,nflux2),ef_cyl_z(mtoroidal,nflux2),ef_cyl_phi(mtoroidal,nflux2)
  !complex(p_):: ef_cyl_r_dft(mtoroidal,nflux2),ef_cyl_z_dft(mtoroidal,nflux2),ef_cyl_phi_dft(mtoroidal,nflux2)
  complex(p_):: potential_dft(mtoroidal,nflux2),out(mtoroidal,nflux2)
  real(p_):: potential_dst(mtoroidal,nflux2) !discrete sine transform
  integer:: i,j,i_left,i_right,j_left,j_right,ierr
  integer:: jeq,ipoloidal
  real(p_):: normal
  real(p_):: tmp1(mtoroidal,nflux2),tmp2(mtoroidal,nflux2),tmp3(mtoroidal,nflux2)

  !----testing
!!$  do i=1,mtoroidal
!!$     do j=1,nflux2
!!$        den_left(i,j)=delta_ne(radcor_1d_array2(j),theta_start,tor_1d_array(i))/ne0
!!$        tmp1(i,j)=delta_ne_theta(radcor_1d_array2(j),theta_start,tor_1d_array(i))/ne0
!!$        tmp2(i,j)=delta_ne_psi(radcor_1d_array2(j),theta_start,tor_1d_array(i))/ne0
!!$        tmp3(i,j)=delta_ne_alpha(radcor_1d_array2(j),theta_start,tor_1d_array(i))/ne0
!!$     enddo
!!$  enddo !---testing end

  potential=den_left !assuming adiabatic electrons and quasineutrality, potential normalized by Te/e

  !smoothing is moved outside of this subroutine
  if(filter_toroidal.eqv..true.) then !filter over the toroidal mode number, keeping the perturbation with desired toroidal mode number
     call oned_fourier_transform1(potential,potential_dft,mtoroidal,nflux2) !calculating 1d DFT of s(:,:) along the first dimension
     call toroidal_filter(potential_dft,out,mtoroidal,nflux2)
     potential_dft=out
     call oned_backward_fourier_transform1(potential_dft,potential,mtoroidal,nflux2)
  endif

  if(filter_radial.eqv..true.) then !filter over the radial mode number, keeping only low-radial-harmonics of the perturbation
!!$     call oned_fourier_transform2(potential,potential_dft,mtoroidal,nflux2) !calculating 1d DFT of s(:,:) along the second dimension
!!$     call radial_fourier_filter(potential_dft,out,mtoroidal,nflux2)
!!$     potential_dft=out
!!$     call oned_backward_fourier_transform2(potential_dft,potential,mtoroidal,nflux2)
     call oned_sine_transform2(potential,potential_dst,mtoroidal,nflux2) !calculating 1d DST of s(:,:) along the second dimension
     call radial_sine_filter(potential_dst,mtoroidal,nflux2)
     call oned_inverse_sine_transform2(potential_dst,potential,mtoroidal,nflux2)
  endif

  call radial_derivative  (potential,potential_psi,dradcor)
  call toroidal_derivative(potential,potential_alpha,dtor)
  call theta_derivative(potential,potential_theta,mtoroidal,nflux2) !partial derivative with respect to theta with psi and alpha fixed, i.e., along the magnetic field line
  !write(*,*) potential_theta(10,40),tmp1(10,40),'myid=',myid
  !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ipoloidal=1+nint((theta_start-theta_1d_array(1))/dtheta)
  do i=1,mtoroidal
     do j=1,nflux2
        jeq=j-1+j_low2
        ef_cyl_r(i,j)=-potential_psi(i,j)*grad_psi_r_matrix(ipoloidal,jeq)&
             & -potential_theta(i,j)*grad_theta_r_matrix(ipoloidal,jeq)-potential_alpha(i,j)*grad_alpha_r_matrix(ipoloidal,jeq)

        ef_cyl_z(i,j)=-potential_psi(i,j)*grad_psi_z_matrix(ipoloidal,jeq) &
             & -potential_theta(i,j)*grad_theta_z_matrix(ipoloidal,jeq)-potential_alpha(i,j)*grad_alpha_z_matrix(ipoloidal,jeq)
        !ef_cyl_phi(i,j)=-potential_alpha(i,j)/r_mag_surf(ipoloidal,jeq)*grad_alpha_phi_matrix(ipoloidal,jeq) !wrong, an additional 1/R factor is wrongly included, a bug found on 2018-Jan.-4 evening
        ef_cyl_phi(i,j)=-potential_alpha(i,j)/r_mag_surf(ipoloidal,jeq) !corrected
     enddo
  enddo

  normal=te0*kev/(ln*bn*vn_i*charge_e)
  !write(*,*) 'normal=',normal
  ef_cyl_r=ef_cyl_r*normal !normalized by the unit used in the code (Bn*vn_i)
  ef_cyl_z=ef_cyl_z*normal
  ef_cyl_phi=ef_cyl_phi*normal

  do j=1,nflux2
     do i=1,mtoroidal
        ef_cyl_r_left(i,j)=ef_cyl_r(i,j) !store the data in the proper arrays
        ef_cyl_z_left(i,j)=ef_cyl_z(i,j)
        ef_cyl_phi_left(i,j)=ef_cyl_phi(i,j)
     enddo
     ef_cyl_r_left(mtoroidal+1,j)=ef_cyl_r_left(1,j)  !peroidic toroidal boundary condition
     ef_cyl_z_left(mtoroidal+1,j)=ef_cyl_z_left(1,j) 
     ef_cyl_phi_left(mtoroidal+1,j)=ef_cyl_phi_left(1,j)
  enddo

!!$  do i=1,mtoroidal+1 !fixed zero radial boundary condition, this is not necessary because the sine transform used later will enforce that the field is set to zero at j=1-1 and j=nflux2+1
!!$     ef_cyl_r_left(i,1)=0._p_
!!$     ef_cyl_z_left(i,1)=0._p_
!!$     ef_cyl_phi_left(i,1)=0._p_
!!$
!!$     ef_cyl_r_left(i,nflux2)=  0._p_
!!$     ef_cyl_z_left(i,nflux2)=  0._p_
!!$     ef_cyl_phi_left(i,nflux2)= 0._p_
!!$  enddo
  call communicate_field_value_between_neighbour_cells2() !for adiabatic electrons model
end subroutine adiabatic_electron_field_solver


subroutine average_electric_field()
use perturbation_field_matrix
implicit none
  ef_cyl_r_left=0.5_p_*(ef_cyl_r_left+ef_cyl_r_left_old)
  ef_cyl_z_left=0.5_p_*(ef_cyl_z_left+ef_cyl_z_left_old)
  ef_cyl_phi_left=0.5_p_*(ef_cyl_phi_left+ef_cyl_phi_left_old)

  ef_cyl_r_right=0.5_p_*(ef_cyl_r_right+ef_cyl_r_right_old)
  ef_cyl_z_right=0.5_p_*(ef_cyl_z_right+ef_cyl_z_right_old)
  ef_cyl_phi_right=0.5_p_*(ef_cyl_phi_right+ef_cyl_phi_right_old)

end subroutine average_electric_field

subroutine store_old_electric_field()
use perturbation_field_matrix
implicit none
  ef_cyl_r_left_old=ef_cyl_r_left
  ef_cyl_z_left_old=ef_cyl_z_left
  ef_cyl_phi_left_old=ef_cyl_phi_left

  ef_cyl_r_right_old=ef_cyl_r_right
  ef_cyl_z_right_old= ef_cyl_z_right
  ef_cyl_phi_right_old=ef_cyl_phi_right

end subroutine store_old_electric_field
