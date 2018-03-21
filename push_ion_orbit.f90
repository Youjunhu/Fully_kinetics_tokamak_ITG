subroutine push_ion_orbit_first_step(dtao_i) !wrapper of "push_full_orbit_cylindrical_boris_with_additional_output" subroutine
  use precision,only:p_
  use constants,only: twopi
  use ions_module,only: nmarker_i,touch_bdry_i
  use ions_module,only: r_i,phi_i,z_i,vr_i,vphi_i,vz_i !as input and output
  use ions_module,only: phi_i_mid
  use ions_module,only: vr_i_mid,vphi_i_mid,vz_i_mid !as output
  implicit none
  real(p_),intent(in):: dtao_i
  integer:: k

  !$omp parallel do
  do k=1,nmarker_i
     if(touch_bdry_i(k).eqv..true.) cycle
     !call push_full_orbit_cylindrical_boris(dtao_i,r_i(k),phi_i(k),z_i(k),vr_i(k),vphi_i(k),vz_i(k))
     call push_full_orbit_cylindrical_boris_with_additional_input_output(dtao_i,r_i(k),phi_i(k),z_i(k),vr_i(k),vphi_i(k),vz_i(k),&
          & phi_i_mid(k),vr_i_mid(k),vphi_i_mid(k),vz_i_mid(k)) !obtain velocity at t_{n+1/2}, outputting the projections of this velocity onto both the basis vectors at t_{n+1/2} and those at t_{n+1}
!!$     !{r,z,phi}_i_mid is known before entering this subroutine ({r,z,phi}_i_mid is ethier the initial codintion of the second-pusher or the output of it.)
  enddo
  !$omp end parallel do
end subroutine push_ion_orbit_first_step


subroutine push_ion_orbit_second_step(dtao_i) 
  !calculate the velocity at t_{n+1}, giving the projections of this velocity onto (1) the basis vectors at t_{n+1} (to prepare for the deposition process) and (2) the basis vectors at t_{n+3/2} (to prepare input for the next secod-pusher)
  !output the location at t_{n+3/2} (to prepare input for the next secod-pusher)
  !(the spatial location at t_{n+1} is already computed by the first boris-pusher)
  use precision,only:p_
  use constants,only: twopi
  use ions_module,only: nmarker_i,touch_bdry_i_mid
  use ions_module,only: r_i_mid,phi_i_mid,z_i_mid !as input (t_{n+1/2}) and output (t_{n+3/2})
  use ions_module,only: vr_i_integer_mid,vphi_i_integer_mid,vz_i_integer_mid !input (projection of v at t_{n} onto basis vector at t_{n+1/2}) and output (projection of v at t_{n+1} onto basis vector at t_{n+3/2})
  use ions_module,only: vr_i_integer,vphi_i_integer,vz_i_integer !as output, projection of velocity at t_{n+1} onto the basis vectors at t_{n+1}
  use ions_module,only: phi_i !as input, location at t_{n+1}
  implicit none
  real(p_),intent(in):: dtao_i
  integer:: k

  do k=1,nmarker_i
     if(touch_bdry_i_mid(k).eqv..true.) cycle
     call push_full_orbit_cylindrical_boris_with_additional_input_output(dtao_i,r_i_mid(k),phi_i_mid(k),z_i_mid(k),&
          & vr_i_integer_mid(k),vphi_i_integer_mid(k),vz_i_integer_mid(k),&
          & phi_i(k),vr_i_integer(k),vphi_i_integer(k),vz_i_integer(k))
  enddo
end subroutine push_ion_orbit_second_step
