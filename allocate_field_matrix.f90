subroutine allocate_field_matrix()
  use precision,only:p_
  use magnetic_coordinates,only: mtoroidal,nflux2,radcor_1d_array,j_low2,j_upp2 !as input
  use perturbation_field_matrix,only: epar_left,ex_left,ey_left, mf_par_left,mf_x_left,mf_y_left
  use perturbation_field_matrix,only: epar_right,ex_right,ey_right, mf_par_right,mf_x_right,mf_y_right
  use perturbation_field_matrix,only: mf_par_old,mf_x_old,mf_y_old
  use perturbation_field_matrix,only: ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left
  use perturbation_field_matrix,only: ef_cyl_r_right,ef_cyl_z_right,ef_cyl_phi_right
  use perturbation_field_matrix,only: ef_cyl_r_left_old,ef_cyl_z_left_old,ef_cyl_phi_left_old
  use perturbation_field_matrix,only: ef_cyl_r_right_old,ef_cyl_z_right_old,ef_cyl_phi_right_old

!  use perturbation_field_matrix,only: jper_x_i_left,jper_y_i_left
  use perturbation_field_matrix,only: jr_left,jphi_left,jz_left,den_left,potential
  use perturbation_field_matrix,only: pper_e_left,ppar_e_left
  use perturbation_field_matrix,only:source_e1,source_e2,source_e3, source_faraday1,source_faraday2,source_faraday3
  use perturbation_field_matrix,only:source_i1,source_i2,source_i3
  use perturbation_field_matrix,only: source1,source2,source3
  use flux_tube_model,only: j_fixed,radcor_fixed
  use domain_decomposition,only: numprocs,myid
  implicit none
  integer:: i,j

  if(myid.eq.0) write(*,*) 'mtoroidal,nflux2=',mtoroidal,nflux2
  allocate(epar_left(mtoroidal+1,nflux2)) !field at left-boundary (smaller theta) of the cell for which the present processor is responsible
  allocate(ex_left(mtoroidal+1,nflux2))
  allocate(ey_left(mtoroidal+1,nflux2))
  allocate(mf_par_left(mtoroidal+1,nflux2))
  allocate(mf_x_left(mtoroidal+1,nflux2))
  allocate(mf_y_left(mtoroidal+1,nflux2))

  allocate(epar_right(mtoroidal+1,nflux2)) !field at right-boundary of the cell for which the present processor is responsible
  allocate(ex_right(mtoroidal+1,nflux2))
  allocate(ey_right(mtoroidal+1,nflux2))
  allocate(mf_par_right(mtoroidal+1,nflux2))
  allocate(mf_x_right(mtoroidal+1,nflux2))
  allocate(mf_y_right(mtoroidal+1,nflux2))

  allocate(mf_par_old(mtoroidal+1,nflux2))
  allocate(mf_x_old(mtoroidal+1,nflux2))
  allocate(mf_y_old(mtoroidal+1,nflux2))

allocate(ef_cyl_r_left(mtoroidal+1,nflux2))
allocate(ef_cyl_z_left(mtoroidal+1,nflux2))
allocate(ef_cyl_phi_left(mtoroidal+1,nflux2))
allocate(ef_cyl_r_right(mtoroidal+1,nflux2))
allocate(ef_cyl_z_right(mtoroidal+1,nflux2))
allocate(ef_cyl_phi_right(mtoroidal+1,nflux2))
ef_cyl_r_left=0._p_
ef_cyl_z_left=0._p_
ef_cyl_phi_left=0._p_
ef_cyl_r_right=0._p_
ef_cyl_z_right=0._p_
ef_cyl_phi_right=0._p_

allocate(ef_cyl_r_left_old(mtoroidal+1,nflux2))
allocate(ef_cyl_z_left_old(mtoroidal+1,nflux2))
allocate(ef_cyl_phi_left_old(mtoroidal+1,nflux2))
allocate(ef_cyl_r_right_old(mtoroidal+1,nflux2))
allocate(ef_cyl_z_right_old(mtoroidal+1,nflux2))
allocate(ef_cyl_phi_right_old(mtoroidal+1,nflux2))



!field equations are solved at the left-boundary of the cell, thus only source at the left-boundary need to be stored
!  allocate(jper_x_i_left(mtoroidal,nflux2))
!  allocate(jper_y_i_left(mtoroidal,nflux2))
  allocate(pper_e_left(mtoroidal,nflux2)) 
  allocate(ppar_e_left(mtoroidal,nflux2))

  allocate(jr_left(mtoroidal,nflux2))
  allocate(jphi_left(mtoroidal,nflux2))
  allocate(jz_left(mtoroidal,nflux2))
  allocate(den_left(mtoroidal,nflux2))
  allocate(potential(mtoroidal,nflux2))

  allocate(source_i1(mtoroidal,nflux2))
  allocate(source_i2(mtoroidal,nflux2))
  allocate(source_i3(mtoroidal,nflux2))

  allocate(source_e1(mtoroidal,nflux2))
  allocate(source_e2(mtoroidal,nflux2))
  allocate(source_e3(mtoroidal,nflux2))

  allocate(source_faraday1(mtoroidal,nflux2))
  allocate(source_faraday2(mtoroidal,nflux2))
  allocate(source_faraday3(mtoroidal,nflux2))

  allocate(source1(mtoroidal,nflux2)) !soure term at left-boundary (smaller theta) of the cell for which the present processor is responsible
  allocate(source2(mtoroidal,nflux2))
  allocate(source3(mtoroidal,nflux2))

  j_fixed=(j_low2+j_upp2)/2
  radcor_fixed=radcor_1d_array(j_fixed) !the radcor of the center of computational region, used in flux tube model

call set_field_initial_value(numprocs,myid)
end subroutine allocate_field_matrix


subroutine set_field_initial_value(numprocs,myid) !for testing whislter waves
  use precision,only:p_
  use constants,only: twopi
  use magnetic_coordinates,only: mtoroidal,nflux2,tor_1d_array,dtor,radcor_1d_array,radcor_upp2,radcor_low2
 use magnetic_coordinates,only: radcor_1d_array2
  use perturbation_field_matrix,only:   nincluded=>toroidal_mode_number_included
  use perturbation_field_matrix,only: epar_left,ex_left,ey_left, mf_par_left,mf_x_left,mf_y_left
  use perturbation_field_matrix,only: ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left
  use perturbation_field_matrix,only: ef_cyl_r_right,ef_cyl_z_right,ef_cyl_phi_right
!  use perturbation_field_matrix,only: epar_right,ex_right,ey_right, mf_par_right,mf_x_right,mf_y_right
implicit none
integer,intent(in)::numprocs, myid
  real(p_),parameter:: eps=1.0d-3
  real(p_)::random_yj,rannum1,tmp
  integer:: iseed,next_seed
  real(p_):: kr
integer:: i,j

  iseed=-(2777+myid*3) !set the iseed in different procs, when using this, it is a parallel generator, but the random numbers in different procs may be related if the iseed chosen for differnt procs is not good enough
  call sub_random_yj(iseed,next_seed,tmp) !just to trigger the use of the iseed, the generated random number is not used, 

  kr=10*twopi/(radcor_upp2-radcor_low2) !radial wave number
  !initial value of the perturbed field is set to be a small value, for testing
  ex_left=0._p_
  ey_left=0._p_
  epar_left=0._p_
  do i=1,mtoroidal
     do j=1,nflux2
!        call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
 !       mf_par_left(i,j)=eps*rannum1
        mf_par_left(i,j)=eps*cos(kr*radcor_1d_array2(j))*cos(real(myid-1)/numprocs*twopi)&
             & *sin(nincluded*tor_1d_array(i))
!        mf_par_left(i,j)=eps*cos(kr*radcor)*cos(theta)*
     enddo
  enddo
  mf_x_left=0._p_
  mf_y_left=0._p_
 mf_par_left=0._p_
!!$  ex_right=0._p_
!!$  ey_right=0._p_
!!$  epar_right=0._p_
!!$  do i=1,mtoroidal
!!$     do j=1,nflux2
!!$        mf_par_right(i,j)=eps*random_yj(0._p_)
!!$     enddo
!!$  enddo
!!$  mf_x_right=0._p_
!!$  mf_y_right=0._p_

ef_cyl_r_left=0._p_
ef_cyl_z_left=0._p_
ef_cyl_phi_left=0._p_
ef_cyl_r_right=0._p_
ef_cyl_z_right=0._p_
ef_cyl_phi_right=0._p_

end subroutine set_field_initial_value

