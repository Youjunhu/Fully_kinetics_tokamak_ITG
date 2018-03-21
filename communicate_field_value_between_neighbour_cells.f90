subroutine communicate_field_value_between_neighbour_cells()
  use mpi
  use precision,only:p_
  use perturbation_field_matrix,only: ex_left,ey_left,epar_left,mf_x_left,mf_y_left,mf_par_left !already known before entering this subroutine
  use perturbation_field_matrix,only: ex_right,ey_right,epar_right,mf_x_right,mf_y_right,mf_par_right !as output
  use magnetic_coordinates,only: m=>mtoroidal,n=>nflux2
  use domain_decomposition,only:myid,GCLR,Tube_comm,ntube,GCLR_cut,my_left,my_right
  use connection_condition,only:connection_condition_at_theta_cut,components_transformation_at_theta_cut
  implicit none

  integer:: status(MPI_STATUS_SIZE),ierr
  !communication between neighbour cells: Every proc is response for one cell which has two boundary grids, only the field on left-boundary-grid is computed by the present proc, the field on the right-boundary is received from the neighbour proc. Note that the definition of cell here is different from the the definition of cell in PIC: the grids are the centers of the cells defined in PIC while the grids are the boundaries of the cells mentioned here.

  call MPI_Sendrecv(ex_left,  (m+1)*n,  MPI_real8, my_left,  4,&
       &            ex_right, (m+1)*n,  MPI_real8, my_right, 4,Tube_COMM,status,ierr)
  call MPI_Sendrecv(ey_left,  (m+1)*n,  MPI_real8, my_left,  5,&
       &            ey_right, (m+1)*n,  MPI_real8, my_right, 5,Tube_COMM,status,ierr)
  call MPI_Sendrecv(epar_left, (m+1)*n,  MPI_real8, my_left,  6,&
       &            epar_right,(m+1)*n,  MPI_real8, my_right, 6,Tube_COMM,status,ierr)
  call MPI_Sendrecv(mf_x_left,  (m+1)*n,  MPI_real8, my_left,  7,&
       &            mf_x_right, (m+1)*n,  MPI_real8, my_right, 7,Tube_COMM,status,ierr)
  call MPI_Sendrecv(mf_y_left,  (m+1)*n,  MPI_real8, my_left,  8,&
       &            mf_y_right, (m+1)*n,  MPI_real8, my_right, 8,Tube_COMM,status,ierr)
  call MPI_Sendrecv(mf_par_left,  (m+1)*n,  MPI_real8, my_left,  9,&
       &            mf_par_right, (m+1)*n,  MPI_real8, my_right, 9,Tube_COMM,status,ierr)

  if(GCLR.eq.GCLR_cut) then !special treatment at theta cut
     call connection_condition_at_theta_cut(ex_right) 
     call connection_condition_at_theta_cut(ey_right) 
     call connection_condition_at_theta_cut(epar_right) 
     call connection_condition_at_theta_cut(mf_x_right) 
     call connection_condition_at_theta_cut(mf_y_right) 
     call connection_condition_at_theta_cut(mf_par_right)
     !the basis vectors grad_alpha, in term of which the pertubed electric field and magnetic field are decomposed, is discontinuous across the theta-cut. Therefore the components at one side need to be transformed to the components on the other side
     call  components_transformation_at_theta_cut(ex_right,ey_right)
     call  components_transformation_at_theta_cut(mf_x_right,mf_y_right)
  endif

  !   if(myid.eq.0) write(*,*) 'field=, after',ex_left(m/3,n/3), ey_left(m/3,n/3)
end subroutine communicate_field_value_between_neighbour_cells

subroutine communicate_field_value_between_neighbour_cells2() !for adiabatic electrons model
  use mpi
  use precision,only:p_
  use perturbation_field_matrix,only: ef_cyl_r_left, ef_cyl_z_left, ef_cyl_phi_left !already known before entering this subroutine
  use perturbation_field_matrix,only: ef_cyl_r_right, ef_cyl_z_right, ef_cyl_phi_right !as output
  use magnetic_coordinates,only: m=>mtoroidal,n=>nflux2
  use domain_decomposition,only:myid,GCLR,Tube_comm,ntube,GCLR_cut,my_left,my_right
  use connection_condition,only:  connection_condition_at_theta_cut
  implicit none

  integer:: status(MPI_STATUS_SIZE),ierr
  !communication between neighbour cells: Every proc is response for one cell which has two boundary grids, only the field on left-boundary-grid is computed by the present proc, the field on the right-boundary is received from the neighbour proc. Note that the definition of cell here is different from the the definition of cell in PIC: the grids are the centers of the cells defined in PIC while the grids are the boundaries of the cells mentioned here.

  call MPI_Sendrecv(ef_cyl_r_left,  (m+1)*n,  MPI_real8, my_left,  4,&
       &            ef_cyl_r_right, (m+1)*n,  MPI_real8, my_right, 4,Tube_COMM,status,ierr)
  call MPI_Sendrecv(ef_cyl_z_left,  (m+1)*n,  MPI_real8, my_left,  5,&
       &            ef_cyl_z_right, (m+1)*n,  MPI_real8, my_right, 5,Tube_COMM,status,ierr)
  call MPI_Sendrecv(ef_cyl_phi_left, (m+1)*n,  MPI_real8, my_left,  6,&
       &            ef_cyl_phi_right,(m+1)*n,  MPI_real8, my_right, 6,Tube_COMM,status,ierr)

  if(GCLR.eq.GCLR_cut) then !special treatment at theta cut
     call connection_condition_at_theta_cut(ef_cyl_r_right) 
     call connection_condition_at_theta_cut(ef_cyl_z_right) 
     call connection_condition_at_theta_cut(ef_cyl_phi_right) 
  endif
  !   if(myid.eq.0) write(*,*) 'field=, after',ex_left(m/3,n/3), ey_left(m/3,n/3)
end subroutine communicate_field_value_between_neighbour_cells2





