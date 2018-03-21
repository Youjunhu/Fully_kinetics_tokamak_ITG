subroutine cal_magnetic_coordinates(mpoloidal,nflux) !to obtain R(radcor,theta), Z(radcor,theta)
  use precision,only:p_
  use constants,only: twopi,pi
  use boundary,only:x_lcfs,z_lcfs,np_lcfs
  use radial_module,only: r_axis,z_axis,psi_lcfs,psi_axis
  use magnetic_coordinates,only: r_mag_surf,z_mag_surf,pfn,gs_psi_array !as output
  use magnetic_coordinates,only: pfn_inner,pfn_bdry,theta_1d_array,radcor_1d_array,dtheta,dradcor !as output
  use magnetic_coordinates,only:minor_r_array,minor_r_prime_array
  use magnetic_coordinates,only:i_theta_zero
  use normalizing,only: Ln
  use domain_decomposition,only: myid
  implicit none
  integer,intent(in):: mpoloidal,nflux
  real(p_),dimension(:,:),allocatable::r_mag_surf0, z_mag_surf0

  real(p_):: theta_uniform(mpoloidal)
  real(p_):: radcor_as_func_of_pfn,minor_r_radcor,minor_r_prime !function name
  integer:: i,j,ierr

  allocate(gs_psi_array(nflux))
  allocate(pfn(nflux))  
  allocate(radcor_1d_array(nflux))
  call create_gs_psi_array(nflux,pfn_inner,pfn_bdry,gs_psi_array,pfn,radcor_1d_array,myid) !select some flux surfaces (labeld by gs_psi_array)
  allocate(minor_r_array(nflux))
  allocate(minor_r_prime_array(nflux))
  allocate(r_mag_surf0(np_lcfs,nflux))
  allocate(z_mag_surf0(np_lcfs,nflux))
  allocate(r_mag_surf(mpoloidal,nflux))
  allocate(z_mag_surf(mpoloidal,nflux))

  !$omp parallel do  
  do j=1,nflux
     call contour(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis,gs_psi_array(j),r_mag_surf0(:,j),z_mag_surf0(:,j))
  enddo
  !$omp end parallel do

  if(myid.eq.0) then !diagnostic
     open(113,file='mag_surf_shape0.txt')
     do j=1,nflux
        do i=1,np_lcfs
           write(113,*) r_mag_surf0(i,j),z_mag_surf0(i,j)
        enddo
        write(113,*)
        write(113,*)
     enddo
     close(113)
  endif

  do j=1,nflux
     minor_r_array(j)=abs(r_mag_surf0(1,j)-r_axis) !minor radius on the high-field-side midplane, this array is to be used in the interpolation appearing in the function minor_r_radcor
  enddo
  call one_dimensional_derivative(nflux,pfn,minor_r_array,minor_r_prime_array) !compute minor_r_prime_array, this array is to be used as interpolating table in the function minor_r_prime

  do i=1,mpoloidal !First, constructe a uniform theta grids
     !theta_uniform(i)=0._p_+twopi/(mpoloidal-1)*(i-1) ! theta ranging from 0 to twopi
     theta_uniform(i)=-pi+twopi/(mpoloidal-1)*(i-1)  ! theta ranging from -pi to +pi
  enddo
  allocate(theta_1d_array(mpoloidal))
  theta_1d_array=theta_uniform
  i_theta_zero=(mpoloidal+1)/2 !poloidal index corresponding to theta=0, which is the low-field-side midplane for up-down symmetric flux surfaces
  !$omp parallel do  
  do j=1,nflux
     call construct_poloidal_coordinate(r_mag_surf0(:,j),z_mag_surf0(:,j),np_lcfs, &
          & r_mag_surf(:,j),z_mag_surf(:,j),mpoloidal,theta_uniform)
  enddo
  !$omp end parallel do

  if (myid.eq.0) then !diagnostics
     open(113,file='theta_line.txt')
     do i=1,mpoloidal
        do j=1,nflux
           write(113,*) r_mag_surf(i,j),z_mag_surf(i,j)
        enddo
        write(113,*)
        write(113,*)
     enddo
     close(113)

     open(113,file='mag_surf_shape.txt')
     do j=1,nflux
        do i=1,mpoloidal
           write(113,*) r_mag_surf(i,j),z_mag_surf(i,j)
        enddo
        write(113,*)
        write(113,*)
     enddo
     close(113)
  endif

  !transform to the Ln unit
  r_mag_surf=r_mag_surf/Ln
  z_mag_surf=z_mag_surf/Ln

  dtheta=twopi/(mpoloidal-1) !grid interval, uniform theta grid is assumed
  dradcor=radcor_1d_array(2)-radcor_1d_array(1) !radial grid interval, uniform grid is assumed

end subroutine cal_magnetic_coordinates

subroutine create_gs_psi_array(nflux,pfn_inner,pfn_bdry,gs_psi_array,pfn,radcor_1d_array,myid)
  !select some poloidal_flux_values that are used in finding magnetic surfaces
  use precision,only:p_
  use radial_module,only: psi_axis,psi_lcfs !as input
  implicit none
  integer,intent(in)::nflux,myid
  real(p_),intent(in):: pfn_inner,pfn_bdry
  real(p_),intent(out):: gs_psi_array(nflux),pfn(nflux),radcor_1d_array(nflux)

  real(p_):: radcor_as_func_of_pfn,dpfn,q_func2
  integer:: j

  dpfn=(pfn_bdry-pfn_inner)/(nflux-1)
  do j=1,nflux
     pfn(j)=pfn_inner+dpfn*(j-1)
  enddo

  do j=1,nflux
     gs_psi_array(j)=psi_axis+pfn(j)*(psi_lcfs-psi_axis)
  enddo

  do j=1,nflux
     radcor_1d_array(j)=radcor_as_func_of_pfn(pfn(j))
  enddo
end subroutine create_gs_psi_array


function radcor_as_func_of_pfn(pfn) result (z)
  use precision,only: p_
  implicit none
  real(p_):: pfn,z

  z=pfn

end function radcor_as_func_of_pfn

function gs_psi_prime(radcor) result (z) !radcor is assumed to be pfn (normalized poloidal flux)
  use precision,only: p_
  use normalizing,only: bn,Ln
  use radial_module,only:psi_axis,psi_lcfs
  implicit none
  real(p_):: radcor,z

  z=(psi_lcfs-psi_axis)/(bn*Ln**2)

end function gs_psi_prime

subroutine magnetic_coordinates_to_cylindrical_coordinates(theta,radcor,r,z) !given (theta,radcor), return (R,Z)
  use precision,only:p_
  use magnetic_coordinates,only: mpoloidal,nflux,r_mag_surf,z_mag_surf,theta_1d_array,radcor_1d_array
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  real(p_),intent(in):: theta,radcor
  real(p_),intent(out):: r,z
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,r_mag_surf,theta,radcor,R)  !uniform 1darray is assumed
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,z_mag_surf,theta,radcor,Z)  !uniform 1darray is assumed
end subroutine magnetic_coordinates_to_cylindrical_coordinates

function minor_r_radcor(radcor) result (z)
  use precision,only: p_
  use constants,only: two,twopi
  use magnetic_coordinates,only: nflux,radcor_1d_array,minor_r_array
  use interpolate_module,only: linear_1d_interpolation
  implicit none
  real(p_):: radcor,z

  call linear_1d_interpolation(nflux,radcor_1d_array,minor_r_array,radcor,z)  
end function minor_r_radcor

function minor_r_prime(radcor) result (z) !derivative of minor_r with respect to the radial coordinate
 use precision,only: p_
  use constants,only: two
  use magnetic_coordinates,only: nflux,radcor_1d_array,minor_r_prime_array
  use interpolate_module,only: linear_1d_interpolation
  implicit none
!  real(p_),parameter:: a=0.45 !minor radius of the LFCS, in unit of Ln, which is 1 meter
real(p_):: radcor,z

!minor_r=sqrt(radcor)*a
!z=sqrt(radcor)*a !needs to be revised
!z=a/(two*sqrt(radcor)) !needs to be revised
  call linear_1d_interpolation(nflux,radcor_1d_array,minor_r_prime_array,radcor,z)  
end function minor_r_prime
