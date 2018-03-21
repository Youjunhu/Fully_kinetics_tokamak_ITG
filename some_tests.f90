subroutine some_tests()
use precision,only:p_
use magnetic_coordinates,only:mpoloidal,nflux,theta_1d_array,radcor_1d_array
implicit none
integer :: i,j

open(11,file='theta_1d_array')
do i=1,mpoloidal
write(11,*) i, theta_1d_array(i)

enddo
close(11)

open(11,file='radcor_1d_array')
do j=1,nflux
write(11,*) j, radcor_1d_array(j)
enddo
close(11)

end subroutine some_tests


subroutine some_test2(myid,numprocs)
  use mpi
  use precision,only:p_
  use electrons_module,only:nmarker_e,  touch_bdry_e, active_e,ntouch_bdry_e, total_ntouch_bdry_e
  use electrons_module,only: radcor_e,theta_e,alpha_e,vpar_e,mu_e,ps_vol_e,ne0
  use array_in_mc,only: b_mc_matrix
  use magnetic_coordinates,only: mpoloidal,nflux,theta_1d_array,radcor_1d_array,vol2
  use interpolate_module
  implicit none
  integer,intent(in):: myid,numprocs
  real(p_):: sum,my_sum,b_val
  real(p_):: fe0 !function name
  integer :: k,ierr

  my_sum=0.
  do k=1,nmarker_e
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,b_mc_matrix,theta_e(k),radcor_e(k),b_val)
     my_sum=my_sum+ps_vol_e(k)*fe0(sqrt(vpar_e(k)**2+2*mu_e(k)*b_val))
  enddo

  !my_sum=0.1
  call  MPI_Reduce(my_sum, sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
  if(myid.eq.0) then
     write(*,*) 'electron sum=',sum, 'myid=',myid, 'exact total_number=', ne0*vol2
  endif
end subroutine some_test2


subroutine some_test3(myid,numprocs)
  use mpi
  use precision,only:p_
  use ions_module,only: nmarker_i,w_i,ps_vol_i,ni0,vr_i,vz_i,vphi_i
 use magnetic_coordinates,only: vol2
  implicit none
  integer,intent(in):: myid,numprocs
  real(p_):: sum,my_sum,v
  real(p_):: fi0 !function name
  integer :: k,ierr

  my_sum=0.
  do k=1,nmarker_i
     v=sqrt(vr_i(k)**2+vz_i(k)**2+vphi_i(k)**2)
     my_sum=my_sum+ps_vol_i(k)*fi0(v)
  enddo

  !my_sum=0.1
  call  MPI_Reduce(my_sum, sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
  if(myid.eq.0) then
     write(*,*) 'ion, sum=',sum, 'exact total_number=', ni0*vol2
  endif
end subroutine some_test3


subroutine some_test4(myid,numprocs)
  use mpi
  use precision,only:p_
  use ions_module,only: nmarker_i,w_i,ps_vol_i,vr_i,vz_i,vphi_i,radcor_i
! use magnetic_coordinates,only: vol2
  implicit none
  integer,intent(in):: myid,numprocs
  real(p_):: sum,my_sum,v
  real(p_):: fi_arbitrary !function name
  integer :: k,ierr

  my_sum=0.
  do k=1,nmarker_i
     v=sqrt(vr_i(k)**2+vz_i(k)**2+vphi_i(k)**2)
     my_sum=my_sum+ps_vol_i(k)*fi_arbitrary(radcor_i(k),v)
  enddo

  !my_sum=0.1
  call  MPI_Reduce(my_sum, sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD,ierr)

  if(myid.eq.0) then
     write(*,*) 'ion, sum=',sum
     call exact_value() !to compare
  endif
end subroutine some_test4

subroutine some_test5(myid,numprocs)
  use mpi
  use precision,only:p_
  use normalizing,only: vn_i,vn_e
  use electrons_module,only: nmarker_e,w_e,ps_vol_e,mu_e,vpar_e,radcor_e,theta_e
  use magnetic_coordinates,only:mpoloidal,nflux,radcor_1d_array,theta_1d_array
  use array_in_mc,only: b_mc_matrix
  use interpolate_module
! use magnetic_coordinates,only: vol2
  implicit none
  integer,intent(in):: myid,numprocs
  real(p_):: sum,my_sum,v
  real(p_):: fe_arbitrary !function name
  integer :: k,ierr
  real(p_):: b_val

  my_sum=0.
  do k=1,nmarker_e
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,b_mc_matrix,theta_e(k),radcor_e(k),b_val)
     v=sqrt(mu_e(k)*b_val*2+vpar_e(k)**2)
     my_sum=my_sum+ps_vol_e(k)*fe_arbitrary(radcor_e(k),v)
  enddo

  !my_sum=0.1
  call  MPI_Reduce(my_sum, sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD,ierr)

  if(myid.eq.0) then
     write(*,*) 'electron, sum=',sum
     call exact_value() !to compare
  endif
end subroutine some_test5


subroutine exact_value()
  use precision,only:p_
  use magnetic_coordinates,only: mpoloidal,nflux,j_low2,j_upp2,radcor_1d_array,dv
  implicit none
  real(p_):: ni_arbitrary !function name
  real(p_):: sum
integer:: i,j
  sum=0._p_
  do i=1,mpoloidal
     do j=j_low2,j_upp2
        sum=sum+ni_arbitrary(radcor_1d_array(j))*dv(i,j)
     enddo
  enddo
  write(*,*)  'exact total_number=', sum
end subroutine exact_value

function fi_arbitrary(radcor,v) result (z) ! for tesing Monte-Carlo integration,v in unit of vn, fi0 in unit 1/(Ln**3*vn**3)
  use precision,only: p_
  use constants,only: two,twopi,kev
  use normalizing,only: vn_i,Ln
  use ions_module,only: mass_i
  implicit none
  real(p_):: radcor,v,z
  real(p_):: v_si
  real(p_):: ti0
  real(p_):: ni_arbitrary !function name
  ti0=2.0*kev
!  ni0=1.0d19 !1/m^3
  v_si=v*vn_i
  z=ni_arbitrary(radcor)*sqrt((mass_i/(twopi*ti0))**3)*exp(-mass_i*v_si**2/(two*ti0))
!  z=ni0*sqrt((mass_i/(twopi*ti0))**3)*exp(-mass_i*v_si**2/(two*ti0))
  z=z*(vn_i**3*Ln**3)
end function

function ni_arbitrary(radcor) result (z) ! for tesing Monte-Carlo integration,v in unit of vn, fi0 in unit 1/(Ln**3*vn**3)
  use precision,only: p_
  implicit none
  real(p_):: radcor,z
  real(p_):: ni0
  ni0=1.0d19 !1/m^3
!z=ni0
!  z=ni0*(1-sqrt(radcor))
 !z=ni0*(1-radcor)**2
z=ni0*(1-radcor)**2

end function


function fe_arbitrary(radcor,v) result (z) ! for tesing Monte-Carlo integration,v in unit of vn_e, z in unit 1/(Ln**3*vn**3)
  use precision,only: p_
  use constants,only: two,twopi,kev
  use normalizing,only: vn_e,Ln
  use electrons_module,only: mass_e
  implicit none
  real(p_):: radcor,v,z
  real(p_):: v_si
  real(p_):: te0
  real(p_):: ni_arbitrary !function name
  te0=2.0*kev
!  ni0=1.0d19 !1/m^3
  v_si=v*vn_e
  z=ni_arbitrary(radcor)*sqrt((mass_e/(twopi*te0))**3)*exp(-mass_e*v_si**2/(two*te0))
  z=z*(vn_e**3*Ln**3)
end function


function fe2_arbitrary(radcor,theta,alpha,v) result (z) ! for tesing Monte-Carlo integration,v in unit of vn_e, z in unit 1/(Ln**3*vn**3)
  use precision,only: p_
  use constants,only: two,twopi,kev
  use normalizing,only: vn_e,Ln
  use electrons_module,only: mass_e
  implicit none
  real(p_):: radcor,theta,alpha,v,z
  real(p_):: v_si
  real(p_):: ni_arbitrary !function name
  real(p_):: te0
  te0=2.0*kev
!  ni0=1.0d19 !1/m^3
  v_si=v*vn_e
  z=ni_arbitrary(radcor)*sqrt((mass_e/(twopi*te0))**3)*exp(-mass_e*v_si**2/(two*te0))*cos(theta)*sin(alpha)
  z=z*(vn_e**3*Ln**3)
end function




subroutine record_data(myid, numprocs,a)
  use precision,only:p_
  use magnetic_coordinates,only: mpoloidal,nflux2,mtoroidal,j_low2,j_upp2,radcor_1d_array,dv
  use magnetic_coordinates,only: tor_1d_array
  use domain_decomposition,only: theta_interval,theta_start !as input 
  implicit none
  real(p_):: ni_arbitrary !function name
  ! real(p_):: sum
  integer,intent(in):: myid,numprocs
  real(p_),intent(in):: a(mtoroidal,nflux2)
  character(7):: filename
  integer:: itor,j,file_unit

  write(filename,'(a3,i4.4)') 'num',myid
  file_unit=myid+310
  open(file_unit,file=filename)

     do itor=1,mtoroidal
  do j=1,nflux2
        write(file_unit,'(2i8.4,1(1pe14.5))')  itor,j,a(itor,j)
     enddo
     write(file_unit,'(4(1pe14.5))')
  enddo
  close(file_unit)

  !  write(*,*)  'exact total_number=', sum
end subroutine record_data


subroutine exact_value2(myid, numprocs) !for testing, compare the value given by this subroutine with that given by the deposition subroutine
  use precision,only:p_
  use constants,only: two,twopi,fourpi,kev
  use normalizing,only:vn_e
  use magnetic_coordinates,only: mpoloidal,nflux2,mtoroidal,j_low2,j_upp2,radcor_1d_array
  use magnetic_coordinates,only: tor_1d_array
  use domain_decomposition,only: theta_interval,theta_start !as input 
  use ions_module,only: vmin_i,vmax_i
  !  use electrons_module,only:mass_e
  !use perturbation_field_matrix,only: ntor=>toroidal_mode_number_included
  implicit none
  !  real(p_):: ni_arbitrary !function name
  real(p_):: delta_ne !function name
  real(p_):: initial_delta_f_i !function name
   real(p_):: sum
  integer,intent(in):: myid,numprocs
  integer:: itor,j,jshift,file_unit,jv
  real(p_):: a(mtoroidal,nflux2)
  character(7):: filename
  real(p_):: te0,v,dvelocity,dvol_v
  integer,parameter:: nv=100

  !te0=2.0*kev/(mass_e*vn_e**2)
  !  do i=1,mpoloidal
  do j=1,nflux2
     jshift=j_low2+(j-1)
     do itor=1,mtoroidal
!        a(itor,j)=delta_ne(radcor_1d_array(jshift),theta_start,tor_1d_array(itor)) !*cos(theta_start)*sin(ntor*tor_1d_array(itor))
        sum=0.
        dvelocity=(vmax_i-vmin_i)/(nv-1)
        do jv=1,nv !integration in velocity space
           v=vmin_i+dvelocity*(jv-1)
           dvol_v=fourpi*v*v*dvelocity
           sum=sum+initial_delta_f_i(radcor_1d_array(jshift),theta_start,tor_1d_array(itor),v)*dvol_v
        enddo
        a(itor,j)=sum
     enddo
  enddo

  write(filename,'(a3,i4.4)') 'ana',myid
  file_unit=myid+200
  open(file_unit,file=filename)

  do itor=1,mtoroidal
     do j=1,nflux2
        write(file_unit,'(2i8.4,1(1pe14.5))')  itor,j,a(itor,j)
     enddo
     write(file_unit,'(4(1pe14.5))')
  enddo

  close(file_unit)
  !  write(*,*)  'exact total_number=', sum
end subroutine exact_value2


