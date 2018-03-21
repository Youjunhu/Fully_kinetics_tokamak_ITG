subroutine deposit_electrons(nmarker_e,active_e,radcor_e,theta_e,alpha_e,mu_e,vpar_e,w_e)!with markers' (v,x,w) given, do the deposition to get pressure on grids
  use mpi
  use precision,only:p_
  use constants,only: one,kev
  use normalizing,only: vn_e
  use electrons_module,only: te0,mass_e
  use perturbation_field_matrix,only: den_left
  use array_in_mc,only: b_mc_matrix
  use magnetic_coordinates,only:radcor_1d_array,theta_1d_array,mpoloidal,nflux
  use magnetic_coordinates,only:nflux2,mtoroidal,radcor_1d_array2,j_low2,tor_1d_array,dtor,dradcor,dtheta,jacobian
  use perturbation_field_matrix,only: pper_e_left,ppar_e_left !as output
  use domain_decomposition,only: numprocs,myid,GRID_COMM,TUBE_COMM,theta_interval,theta_start,GCLR,ntube,GCLR_cut,my_left,my_right
use electrons_module,only:ne0
use connection_condition
  use interpolate_module
  implicit none
  integer,intent(in)::nmarker_e
  logical,intent(in):: active_e(nmarker_e)
  real(p_),intent(in):: radcor_e(nmarker_e),theta_e(nmarker_e),alpha_e(nmarker_e),w_e(nmarker_e),mu_e(nmarker_e),vpar_e(nmarker_e)
  real(p_):: my_pper_e_left(mtoroidal,nflux2),my_pper_e_right(mtoroidal,nflux2)
  real(p_):: my_ppar_e_left(mtoroidal,nflux2),my_ppar_e_right(mtoroidal,nflux2)
  real(p_):: pper_e_left_tmp(mtoroidal,nflux2),pper_e_right(mtoroidal,nflux2)
  real(p_):: ppar_e_left_tmp(mtoroidal,nflux2),ppar_e_right(mtoroidal,nflux2)
  real(p_):: coeff_theta_1,coeff_theta_2,coeff_alpha_1,coeff_alpha_2,coeff_radcor_1,coeff_radcor_2
  real(p_):: kernel,dv,b_val
  integer:: i,j,k,jshift,ipoloidal,i_plus_one,j_plus_one
  integer:: status(MPI_STATUS_SIZE),ierr

  my_pper_e_left =0._p_ !before the deposition, set the array to zero
  my_pper_e_right=0._p_ !before the deposition, set the array to zero
  my_ppar_e_left =0._p_ !before the deposition, set the array to zero
  my_ppar_e_right=0._p_ !before the deposition, set the array to zero

do k=1,nmarker_e
     if(active_e(k).eqv..false.) cycle ! markers outside the computational region do not contribute anything to any grids
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,b_mc_matrix,theta_e(k),radcor_e(k),b_val)
     kernel=w_e(k)*mu_e(k)*b_val
!     kernel=w_e(k) !for testing
!if(abs(w_e(k))>0.0000001) write(*,*) 'w_e(k) is nonzero', w_e(k),k

!     if(myid.eq.0) write(*,*) 'elctron bval,mu,w*mu*bval=',b_val,mu_e(k),kernel
     !determine the interpolating coeefficeint according to the location of a marker
     coeff_theta_1=(theta_e(k)-theta_start)/theta_interval
     coeff_theta_2=one-coeff_theta_1

!     call location(mtoroidal,tor_1d_array,alpha_e(k),i)
     i=floor((alpha_e(k)-tor_1d_array(1))/dtor+1)!uniform xarray is assumed, otherwise we need to call location() subroutine to locate xval
     coeff_alpha_1=(alpha_e(k)-tor_1d_array(i))/dtor
     coeff_alpha_2=one-coeff_alpha_1

!     call location(nflux2, radcor_1d_array2, radcor_e(k),j)
     j=floor((radcor_e(k)-radcor_1d_array2(1))/dradcor+1) !uniform xarray is assumed, otherwise we need to call location() subroutine to locate xval
     coeff_radcor_1= (radcor_e(k)-radcor_1d_array2(j))/dradcor
     coeff_radcor_2=one-coeff_radcor_1

     !deposit markers to the No. (i,j) grid, using the above coefficients
     my_pper_e_left(i,j)=  my_pper_e_left(i,j)+kernel*coeff_theta_2*coeff_alpha_2*coeff_radcor_2
     my_pper_e_right(i,j)=my_pper_e_right(i,j)+kernel*coeff_theta_1*coeff_alpha_2*coeff_radcor_2

     i_plus_one=i+1
     if(i.eq.mtoroidal) i_plus_one=1 !periodic condition
     j_plus_one=j+1
     if(j.eq.nflux2) j_plus_one=1 !periodic condition
     my_pper_e_left(i_plus_one,j)=my_pper_e_left(i_plus_one,j)  +kernel*coeff_theta_2*coeff_alpha_1*coeff_radcor_2
     my_pper_e_right(i_plus_one,j)=my_pper_e_right(i_plus_one,j)+kernel*coeff_theta_1*coeff_alpha_1*coeff_radcor_2

     my_pper_e_left(i,j_plus_one)= my_pper_e_left(i,j_plus_one) +kernel*coeff_theta_2*coeff_alpha_2*coeff_radcor_1
     my_pper_e_right(i,j_plus_one)=my_pper_e_right(i,j_plus_one)+kernel*coeff_theta_1*coeff_alpha_2*coeff_radcor_1

     my_pper_e_left(i_plus_one,j_plus_one)=my_pper_e_left(i_plus_one,j_plus_one)+  kernel*coeff_theta_2*coeff_alpha_1*coeff_radcor_1
     my_pper_e_right(i_plus_one,j_plus_one)=my_pper_e_right(i_plus_one,j_plus_one)+kernel*coeff_theta_1*coeff_alpha_1*coeff_radcor_1

     kernel=w_e(k)*vpar_e(k)**2 !for computing perturbed electron parallel pressure

     my_ppar_e_left(i,j)=  my_ppar_e_left(i,j)+kernel*coeff_theta_2*coeff_alpha_2*coeff_radcor_2
     my_ppar_e_right(i,j)=my_ppar_e_right(i,j)+kernel*coeff_theta_1*coeff_alpha_2*coeff_radcor_2
     my_ppar_e_left(i_plus_one,j)=my_ppar_e_left(i_plus_one,j)  +kernel*coeff_theta_2*coeff_alpha_1*coeff_radcor_2
     my_ppar_e_right(i_plus_one,j)=my_ppar_e_right(i_plus_one,j)+kernel*coeff_theta_1*coeff_alpha_1*coeff_radcor_2
     my_ppar_e_left(i,j_plus_one)= my_ppar_e_left(i,j_plus_one) +kernel*coeff_theta_2*coeff_alpha_2*coeff_radcor_1
     my_ppar_e_right(i,j_plus_one)=my_ppar_e_right(i,j_plus_one)+kernel*coeff_theta_1*coeff_alpha_2*coeff_radcor_1
     my_ppar_e_left(i_plus_one,j_plus_one)=my_ppar_e_left(i_plus_one,j_plus_one)+  kernel*coeff_theta_2*coeff_alpha_1*coeff_radcor_1
     my_ppar_e_right(i_plus_one,j_plus_one)=my_ppar_e_right(i_plus_one,j_plus_one)+kernel*coeff_theta_1*coeff_alpha_1*coeff_radcor_1
enddo

 do j=1,nflux2  !divided by the spatial volume of a cell, to give the pressure, note that this 'cell' is the cell defined by PIC (i.e., grid is the center of the cell)
     jshift=j_low2+(j-1)
     ipoloidal=1+nint((theta_start-theta_1d_array(1))/dtheta) !index of begining theta angle of this subdomain in the original magnetic coordinate grids
     dv=abs(jacobian(ipoloidal,jshift))*dradcor*theta_interval*dtor !volume of the cell (the center of the cell is the grid)
     do i=1,mtoroidal
        my_pper_e_left(i,j)= my_pper_e_left (i,j)/dv
        my_ppar_e_left(i,j)= my_ppar_e_left (i,j)/dv
     enddo
  enddo

  do j=1,nflux2
     jshift=j_low2+(j-1)
     ipoloidal=1+nint((theta_start+theta_interval-theta_1d_array(1))/dtheta) !index of ending theta angle of this subdomain in the original magnetic coordinate grids
     dv=abs(jacobian(ipoloidal,jshift))*dradcor*theta_interval*dtor !volume of the cell (the center of the cell is the grid)
     do i=1,mtoroidal
        my_pper_e_right(i,j)=my_pper_e_right(i,j)/dv
        my_ppar_e_right(i,j)=my_ppar_e_right(i,j)/dv
     enddo
  enddo

my_pper_e_left=my_pper_e_left/ne0 !normalization
my_pper_e_right=my_pper_e_right/ne0 !normalization
my_ppar_e_left=my_ppar_e_left/ne0 !normalization
my_ppar_e_right=my_ppar_e_right/ne0 !normalization

!summing over all those procs that are in the same cell (ntube procs in each cell)
     call MPI_ALLREDUCE(my_pper_e_left, pper_e_left, mtoroidal*nflux2,MPI_REAL8,MPI_SUM,GRID_COMM,ierr)
     call MPI_ALLREDUCE(my_pper_e_right,pper_e_right,mtoroidal*nflux2,MPI_REAL8,MPI_SUM,GRID_COMM,ierr)
     call MPI_ALLREDUCE(my_ppar_e_left, ppar_e_left, mtoroidal*nflux2,MPI_REAL8,MPI_SUM,GRID_COMM,ierr)
     call MPI_ALLREDUCE(my_ppar_e_right,ppar_e_right,mtoroidal*nflux2,MPI_REAL8,MPI_SUM,GRID_COMM,ierr)

 !communication between neighbour cells, single proc in a cell
!!$  my_right=myid+1
!!$  if(myid.eq.numprocs-1) my_right=0
!!$  my_left=myid-1
!!$  if(myid.eq.0) my_left=numprocs-1

!!$  my_right=GCLR+1 !these values are now pre-computed and stored in a module, no need to re-compute this in each subroutine
!!$  if(GCLR.eq.numprocs/ntube-1) my_right=0
!!$  my_left=GCLR-1
!!$  if(GCLR.eq.0) my_left=numprocs/ntube-1

  call MPI_Sendrecv(pper_e_right,   mtoroidal*nflux2, MPI_real8, my_right, 1,&
       &            pper_e_left_tmp,mtoroidal*nflux2, MPI_real8, my_left,  1,Tube_COMM,status,ierr)
  call MPI_Sendrecv(ppar_e_right,   mtoroidal*nflux2, MPI_real8, my_right, 2,&
       &            ppar_e_left_tmp,mtoroidal*nflux2, MPI_real8, my_left,  2,Tube_COMM,status,ierr)

! if(GCLR.eq.0) then !special treatment at the theta cut
 if(GCLR.eq.GCLR_cut+1) then !special treatment at the theta cut
     call connection_condition_at_theta_cut2(pper_e_left_tmp)
     call connection_condition_at_theta_cut2(ppar_e_left_tmp)
  endif

  pper_e_left=pper_e_left+pper_e_left_tmp !add the contribution from the neighbour cell
  ppar_e_left=ppar_e_left+ppar_e_left_tmp !add the contribution from the neighbour cell


!some tests
!call exact_value2(myid, numprocs)
!call record_data(myid, numprocs,pper_e_left)

end subroutine deposit_electrons





