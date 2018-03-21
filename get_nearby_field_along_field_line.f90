subroutine get_nearby_field_along_field_line(a,a_left,a_right,m,n) !thoughly tested 2018-Jan.7
  use precision,only:p_
  use domain_decomposition,only: GCLR,ntube,TUBE_COMM,my_left,my_right,GCLR_cut
  use mpi
  use connection_condition,only:connection_condition_at_theta_cut,connection_condition_at_theta_cut3
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: a(m,n)
  real(p_),intent(out):: a_left(m,n),a_right(m,n)
  integer:: status(MPI_STATUS_SIZE),ierr

  call MPI_Sendrecv(a, m*n, MPI_real8, my_right, 1,&
       &            a_left,m*n, MPI_real8, my_left, 1,Tube_COMM,status,ierr)

  call MPI_Sendrecv(a, m*n, MPI_real8, my_left, 2,&
       &            a_right,m*n, MPI_real8, my_right, 2,Tube_COMM,status,ierr)

  if(GCLR.eq.GCLR_cut) then
     call connection_condition_at_theta_cut(a_right) !a_right array received by GCLR_cut from its right neighbour can not be directly used by GCLR_cut processor because the original a_right array is defined on a set grids which are different from the grids desired by GCLR processor. Interpolation is needed. !Interpolating the numerical table defined on the theta=-pi plane to get values on the grids defined on the theta=+pi plane (theta=+pi and -pi correspond to the same place, but the grids on the planes are different)
  endif

  if(GCLR.eq.GCLR_cut+1) then
     call connection_condition_at_theta_cut3(a_left) !a_left array received by GCLR_cut+1 from its left neighbour can not be directly used by GCLR_cut+1 processor because the original a_left array is defined on a set grids which are different from the grids desired by the GCLR_cut+ processor. Interpolation is needed. !Interpolating the numerical table defined on the theta=+pi-theta_interval plane to get values on the grids defined on the theta=-pi-theta_interval plane (theta=+pi-theta_interval and -pi-theta_interval correspond to the same place, but the grids on the planes are different)
  endif

end subroutine get_nearby_field_along_field_line
