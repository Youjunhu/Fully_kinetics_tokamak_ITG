module deposit_ions_module
private
public deposit_ions
contains
subroutine deposit_ions(nmarker_i,active_i,radcor_i,theta_i,alpha_i,vr_i,vphi_i,vz_i,w_i) !given (r,v,w) of markers, do the deposition to get perpendicular currents and number density on spatial grids
!periodic toroidal boundary condition and the poloidal boundary condtion (connection condition along the field line) are taken into account
  use mpi
  use precision,only:p_
  use constants,only: one,two
  use ions_module,only: ni0
  use magnetic_coordinates,only:nflux2,mtoroidal,radcor_1d_array2,theta_1d_array,j_low2,tor_1d_array,dtor,dradcor,dtheta,jacobian
  use magnetic_coordinates,only:vol1
  use perturbation_field_matrix,only: jr_left,jphi_left,jz_left,den_left !as output
  use domain_decomposition,only: GRID_COMM,TUBE_COMM,theta_interval,theta_start,numprocs,myid,GCLR,ntube,GCLR_cut,my_left,my_right
  use connection_condition
  implicit none
  integer,intent(in)::nmarker_i
  logical,intent(in):: active_i(nmarker_i)
  real(p_),intent(in)::radcor_i(nmarker_i),theta_i(nmarker_i),alpha_i(nmarker_i)
  real(p_),intent(in)::vr_i(nmarker_i),vphi_i(nmarker_i),vz_i(nmarker_i),w_i(nmarker_i)
  real(p_):: jr_right(mtoroidal,nflux2),jphi_right(mtoroidal,nflux2),jz_right(mtoroidal,nflux2),den_right(mtoroidal,nflux2)
  real(p_):: jr_left_tmp(mtoroidal,nflux2),jphi_left_tmp(mtoroidal,nflux2),jz_left_tmp(mtoroidal,nflux2),&
       & den_left_tmp(mtoroidal,nflux2) !temporary arrays

  real(p_):: my_jr_left(mtoroidal,nflux2),my_jr_right(mtoroidal,nflux2)
  real(p_):: my_jz_left(mtoroidal,nflux2),my_jz_right(mtoroidal,nflux2)
  real(p_):: my_jphi_left(mtoroidal,nflux2),my_jphi_right(mtoroidal,nflux2)
  real(p_):: my_den_left(mtoroidal,nflux2),my_den_right(mtoroidal,nflux2)

  real(p_):: coeff_theta_1,coeff_theta_2,coeff_alpha_1,coeff_alpha_2,coeff_radcor_1,coeff_radcor_2
  real(p_):: kernel,dv,mysum,sum
  integer:: k,i,j, jshift,ipoloidal,i_plus_one,j_plus_one
  integer:: status(MPI_STATUS_SIZE),ierr

  my_jr_left=0._p_ !before the deposition, set the array to zero
  my_jr_right=0._p_
  my_jphi_left=0._p_
  my_jphi_right=0._p_
  my_jz_left=0._p_
  my_jz_right=0._p_
  my_den_left=0._p_
  my_den_right=0._p_

  do k=1,nmarker_i !for each marker, deposit it to the corresponding grids
     if(active_i(k).eqv..false.) cycle ! markers outside the computational region do not contribute density or current to any grids
     !determine the interpolating coeefficeint according to the location of a marker
     coeff_theta_1=(theta_i(k)-theta_start)/theta_interval
     coeff_theta_2=one-coeff_theta_1

     !call location(mtoroidal,tor_1d_array,alpha_i(k),i)
     i=floor((alpha_i(k)-tor_1d_array(1))/dtor+1) !uniform xarray is assumed, otherwise we need to call location() subroutine to locate xval
     coeff_alpha_1=(alpha_i(k)-tor_1d_array(i))/dtor
     coeff_alpha_2=one-coeff_alpha_1
     !if(myid.eq.2) write(*,*) 'alpha, i=',i, 'k=',k
     !call location(nflux2, radcor_1d_array2, radcor_i(k),j)
     j=floor((radcor_i(k)-radcor_1d_array2(1))/dradcor+1)
     coeff_radcor_1= (radcor_i(k)-radcor_1d_array2(j))/dradcor
     coeff_radcor_2=one-coeff_radcor_1

     i_plus_one=i+1
     if(i.eq.mtoroidal) i_plus_one=1 !periodic condition in the toroidal direction
     j_plus_one=j+1
     if(j.eq.nflux2) j_plus_one=j !all the active markers are actually within the boundary flux surface labeled by nflux2, if j.eq.nflux2, then the marker is exactly on the boundary flux surface, this code line is needed to avoid exeeding array bounds in case that the marker is exactly on the boundary flux surface or slightly outside of it due to numerical truncation errors
     !if(j.eq.nflux2) j_plus_one=1 !periodic condition. is disabled later in this subroutine, now using fixed zero boundary condition along the radial direction. In the past, I used the periodic bounary condition because I want to use Fourier transform along the radial direction. However this periodic radial boundary condition is not reasonable. The reason is as follows: Here the radial direction is d/dpsi in (psi,theta,alpha) coordinates. This direction is a combinition of the toroidal and the usual radial direction and thus is an artificial direction introduced when using (psi,theta,alpha) coordinates, and does not has any physical reason to justify that a peroidic condition should be satisfied along this artificial direction.

     !deposit markers to the grids, using the above coefficients
     kernel=w_i(k)*vr_i(k) !neglecting the difference between basis vectors at particle location and those at the grid where the particle are deposited
     !kernel=w_i(k) !for testing
     my_jr_left(i,j)= my_jr_left(i,j)+ kernel*coeff_theta_2*coeff_alpha_2*coeff_radcor_2
     my_jr_right(i,j)=my_jr_right(i,j)+kernel*coeff_theta_1*coeff_alpha_2*coeff_radcor_2

     my_jr_left(i_plus_one,j)=my_jr_left(i_plus_one,j)+kernel*coeff_theta_2*coeff_alpha_1*coeff_radcor_2
     my_jr_right(i_plus_one,j)=my_jr_right(i_plus_one,j)+kernel*coeff_theta_1*coeff_alpha_1*coeff_radcor_2

     my_jr_left(i,j_plus_one)=my_jr_left(i,j_plus_one)+kernel*coeff_theta_2*coeff_alpha_2*coeff_radcor_1
     my_jr_right(i,j_plus_one)=my_jr_right(i,j_plus_one)+kernel*coeff_theta_1*coeff_alpha_2*coeff_radcor_1

     my_jr_left(i_plus_one,j_plus_one)=my_jr_left(i_plus_one,j_plus_one)+kernel*coeff_theta_2*coeff_alpha_1*coeff_radcor_1
     my_jr_right(i_plus_one,j_plus_one)=my_jr_right(i_plus_one,j_plus_one)+kernel*coeff_theta_1*coeff_alpha_1*coeff_radcor_1

     kernel=w_i(k)*vphi_i(k) !for another kernel
     !kernel=w_i(k) !for testing
     my_jphi_left(i,j)=my_jphi_left(i,j) +kernel*coeff_theta_2*coeff_alpha_2*coeff_radcor_2
     my_jphi_right(i,j)=my_jphi_right(i,j)+kernel*coeff_theta_1*coeff_alpha_2*coeff_radcor_2

     my_jphi_left(i_plus_one,j)= my_jphi_left(i_plus_one,j)+kernel*coeff_theta_2*coeff_alpha_1*coeff_radcor_2
     my_jphi_right(i_plus_one,j)=my_jphi_right(i_plus_one,j)+kernel*coeff_theta_1*coeff_alpha_1*coeff_radcor_2

     my_jphi_left(i,j_plus_one)=my_jphi_left(i,j_plus_one)+kernel*coeff_theta_2*coeff_alpha_2*coeff_radcor_1
     my_jphi_right(i,j_plus_one)=my_jphi_right(i,j_plus_one)+kernel*coeff_theta_1*coeff_alpha_2*coeff_radcor_1

     my_jphi_left(i_plus_one,j_plus_one)=my_jphi_left(i_plus_one,j_plus_one)+kernel*coeff_theta_2*coeff_alpha_1*coeff_radcor_1
     my_jphi_right(i_plus_one,j_plus_one)=my_jphi_right(i_plus_one,j_plus_one)+kernel*coeff_theta_1*coeff_alpha_1*coeff_radcor_1

     kernel=w_i(k)*vz_i(k) !for another kernel
     my_jz_left(i,j)=my_jz_left(i,j) +kernel*coeff_theta_2*coeff_alpha_2*coeff_radcor_2
     my_jz_right(i,j)=my_jz_right(i,j)+kernel*coeff_theta_1*coeff_alpha_2*coeff_radcor_2

     my_jz_left(i_plus_one,j)= my_jz_left(i_plus_one,j)+kernel*coeff_theta_2*coeff_alpha_1*coeff_radcor_2
     my_jz_right(i_plus_one,j)=my_jz_right(i_plus_one,j)+kernel*coeff_theta_1*coeff_alpha_1*coeff_radcor_2

     my_jz_left(i,j_plus_one)=my_jz_left(i,j_plus_one)+kernel*coeff_theta_2*coeff_alpha_2*coeff_radcor_1
     my_jz_right(i,j_plus_one)=my_jz_right(i,j_plus_one)+kernel*coeff_theta_1*coeff_alpha_2*coeff_radcor_1

     my_jz_left(i_plus_one,j_plus_one)=my_jz_left(i_plus_one,j_plus_one)+kernel*coeff_theta_2*coeff_alpha_1*coeff_radcor_1
     my_jz_right(i_plus_one,j_plus_one)=my_jz_right(i_plus_one,j_plus_one)+kernel*coeff_theta_1*coeff_alpha_1*coeff_radcor_1

     kernel=w_i(k) !density, only used in isotropic fluid electrons model or adiabatic electrons model to provide delta_ne, which is assumed to be equal to delta_ni
     my_den_left(i,j)=my_den_left(i,j) +kernel*coeff_theta_2*coeff_alpha_2*coeff_radcor_2
     my_den_right(i,j)=my_den_right(i,j)+kernel*coeff_theta_1*coeff_alpha_2*coeff_radcor_2

     my_den_left(i_plus_one,j)= my_den_left(i_plus_one,j)+kernel*coeff_theta_2*coeff_alpha_1*coeff_radcor_2
     my_den_right(i_plus_one,j)=my_den_right(i_plus_one,j)+kernel*coeff_theta_1*coeff_alpha_1*coeff_radcor_2

     my_den_left(i,j_plus_one)=my_den_left(i,j_plus_one)+kernel*coeff_theta_2*coeff_alpha_2*coeff_radcor_1
     my_den_right(i,j_plus_one)=my_den_right(i,j_plus_one)+kernel*coeff_theta_1*coeff_alpha_2*coeff_radcor_1

     my_den_left(i_plus_one,j_plus_one)=my_den_left(i_plus_one,j_plus_one)+kernel*coeff_theta_2*coeff_alpha_1*coeff_radcor_1
     my_den_right(i_plus_one,j_plus_one)=my_den_right(i_plus_one,j_plus_one)+kernel*coeff_theta_1*coeff_alpha_1*coeff_radcor_1

  enddo

  !summing over all those procs that are in the same cell (ntube procs in each cell)
  call MPI_ALLREDUCE(my_jr_left,   jr_left,   mtoroidal*nflux2,MPI_REAL8,MPI_SUM,GRID_COMM,ierr)
  call MPI_ALLREDUCE(my_jr_right,  jr_right,  mtoroidal*nflux2,MPI_REAL8,MPI_SUM,GRID_COMM,ierr)
  call MPI_ALLREDUCE(my_jz_left,   jz_left,   mtoroidal*nflux2,MPI_REAL8,MPI_SUM,GRID_COMM,ierr)
  call MPI_ALLREDUCE(my_jz_right,  jz_right,  mtoroidal*nflux2,MPI_REAL8,MPI_SUM,GRID_COMM,ierr)
  call MPI_ALLREDUCE(my_jphi_left, jphi_left, mtoroidal*nflux2,MPI_REAL8,MPI_SUM,GRID_COMM,ierr)
  call MPI_ALLREDUCE(my_jphi_right,jphi_right,mtoroidal*nflux2,MPI_REAL8,MPI_SUM,GRID_COMM,ierr)
  call MPI_ALLREDUCE(my_den_left,   den_left,   mtoroidal*nflux2,MPI_REAL8,MPI_SUM,GRID_COMM,ierr)
  call MPI_ALLREDUCE(my_den_right,  den_right,  mtoroidal*nflux2,MPI_REAL8,MPI_SUM,GRID_COMM,ierr)

  !communication between neighbour cells
  call MPI_Sendrecv(jr_right,mtoroidal*nflux2, MPI_real8, my_right, 2,&
       &            jr_left_tmp, mtoroidal*nflux2,  MPI_real8, my_left, 2,Tube_comm,status,ierr)
  call MPI_Sendrecv(jz_right,mtoroidal*nflux2, MPI_real8, my_right, 3,&
       &            jz_left_tmp,mtoroidal*nflux2, MPI_real8, my_left, 3,Tube_comm,status,ierr)
  call MPI_Sendrecv(jphi_right,mtoroidal*nflux2, MPI_real8, my_right, 4,&
       &            jphi_left_tmp,mtoroidal*nflux2, MPI_real8, my_left, 4,Tube_COMM,status,ierr)
  call MPI_Sendrecv(den_right,mtoroidal*nflux2, MPI_real8, my_right, 3,&
       &            den_left_tmp,mtoroidal*nflux2, MPI_real8, my_left, 3,Tube_comm,status,ierr)
  !  if(GCLR.eq.0) then !special treatment at theta cut (theta=0), current density received from the left-neigbour of No. 0 cell needs to be interpolated to the grids on the left-plane of the No. 0 cell
  !  if(GCLR.eq.GCLR_cut) then !special treatment at the theta cut, wrong! a bug found
  if(GCLR.eq.GCLR_cut+1) then !special treatment at the theta cut
     call connection_condition_at_theta_cut2(jr_left_tmp) 
     call connection_condition_at_theta_cut2(jz_left_tmp) 
     call connection_condition_at_theta_cut2(jphi_left_tmp)
     call connection_condition_at_theta_cut2(den_left_tmp) 
  endif

  jr_left=jr_left+jr_left_tmp !add the contribution from the neighbour cell
  jz_left=jz_left+jz_left_tmp !add the contribution from the neighbour cell
  jphi_left=jphi_left+jphi_left_tmp !add the contribution from the neighbour cell
  den_left=den_left+den_left_tmp !add the contribution from the neighbour cell

  do j=1,nflux2  !divided by the space volume of a cell, to give the current denisty, note that this 'cell' is the cell defined by PIC (i.e., grid is the center of the cell). (while grids are the boundaries of the "mpi cell")
     jshift=j_low2+(j-1)
     ipoloidal=1+nint((theta_start-theta_1d_array(1))/dtheta) !index of begining theta angle of this subdomain in the original magnetic coordinate grids
     dv=abs(jacobian(ipoloidal,jshift))*dradcor*theta_interval*dtor !volume of the cell (the center of the cell is the grid)
     do i=1,mtoroidal
        jr_left(i,j)=jr_left(i,j)/dv
        jphi_left(i,j)=jphi_left(i,j)/dv
        jz_left(i,j)=jz_left(i,j)/dv
        den_left(i,j)=den_left(i,j)/dv
     enddo
  enddo

  jr_left=jr_left/ni0 !current density normalized by ni0*vni*qi
  jz_left=jz_left/ni0 !current density normalized by ni0*vni*qi
  jphi_left=jphi_left/ni0 !current density normalized by ni0*vni*qi
  den_left=den_left/ni0 !number density normalized by ni0

!!$do i=1,mtoroidal !set the radial boundary condition: fixed zero boundary condition, not necessary, so removed these codes
!!$  den_left(i,1)=0._p_ !fixed (zero) boundary condition
!!$  !den_left(i,nflux2)=(den_left(i,nflux2-1)+den_left(i,1))/two !periodic boundary condition, den(j=1)=den(j=nflux2+1), not reasonable to use peroidic boundary along the radial direction. I used this in the past because I used the Fourier transform in the radial direction. Later I switched to the sine transform and thus the peroidic boundary condition became unnecessary
!!$  den_left(i,nflux2)=0._p_ !fixed (zero) boundary condition
!!$  jr_left(i,1)=0._p_
!!$  jr_left(i,nflux2)=0._p_
!!$  jz_left(i,1)=0._p_
!!$  jz_left(i,nflux2)=0._p_
!!$  jphi_left(i,1)=0._p_
!!$  jphi_left(i,nflux2)=0._p_
!!$enddo

  !call MPI_Barrier( MPI_Comm_world, ierr)
  !some tests
!!$call exact_value2(myid, numprocs)
!!$jr_left_tmp=jr_left*ni0
!!$call record_data(myid, numprocs,jr_left_tmp)
!  write(*,*) 'den_left=',den_left(5,5), 'myid=',myid
!  write(*,*) 'den_left=',den_left(10,8), 'myid=',myid
!  write(*,*) 'averaged den_left=',sum(den_left(:,:))/(mtoroidal*nflux2), 'myid=',myid
!  write(*,*) 'den_left, max,min=',maxval(den_left(:,:)),minval(den_left(:,:)), 'myid=',myid
  !write(*,*) 'nmarker_i=', nmarker_i, 'myid=',myid
!!$  mysum=0.
!!$  do k=1,nmarker_i
!!$     mysum=mysum+w_i(k)
!!$  enddo
!!$  call MPI_ALLREDUCE(mysum, sum, 1,MPI_REAL8,MPI_SUM,mpi_COMM_world,ierr)
!!$  sum=sum/vol1
!!$  write(*,*) 'myid=',myid,'sum=',sum
  !    call check_domain_particles(theta_i,nmarker_i)
end subroutine deposit_ions
end module deposit_ions_module
