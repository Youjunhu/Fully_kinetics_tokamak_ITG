program lorentz_ions !Particle in cell simulation of Ion Temperature Gradient (ITG) mode in tokamak with fully kinetic ions and drift-kinetic electrons, project (part of GEM project) started in May of 2017
  use mpi
  use control_parameters,only:kstart,kend,dtao_omega_i_axis,niter,iplot_mode_structure
  use precision,only:p_
  use constants,only: twopi,pi,zero,one_half,kev,two,mu0
  use normalizing,only:ln,bn,tn_i
  use ions_module,only: mass_i,charge_i,dtao_i,nmarker_i,total_nmarker_i,vt_i,ti0,ni0,ti0,kappa_ti
  use ions_module,only: r_i,z_i,phi_i, radcor_i,theta_i,alpha_i,tor_shift_i,ps_vol_i, touch_bdry_i,active_i
  use ions_module,only: touch_bdry_i_mid,active_i_mid
  use ions_module,only: radcor_i_mid,theta_i_mid,alpha_i_mid
  use ions_module,only: r_i_old,z_i_old,phi_i_old
  use ions_module,only: vr_i,vz_i,vphi_i,vpar_i,vx_i,vy_i, ntouch_bdry_i, total_ntouch_bdry_i
  use ions_module,only: vr_i_integer,vz_i_integer,vphi_i_integer
  use ions_module,only: vr_i_mid,vz_i_mid,vphi_i_mid
  use ions_module,only: r_i_mid,z_i_mid,phi_i_mid,vr_i_integer_mid,vz_i_integer_mid,vphi_i_integer_mid
  use ions_module,only: w_i,w_i_mid,w_i_star
  use electrons_module,only:fluid_electron,mass_e, nmarker_e,total_nmarker_e
  use electrons_module,only:  touch_bdry_e, active_e,ntouch_bdry_e, total_ntouch_bdry_e
  use electrons_module,only: radcor_e,theta_e,alpha_e,tor_shift_e,rg_e,zg_e,phig_e,vpar_e,mu_e,ps_vol_e,w_e
  use electrons_module,only:radcor_e_mid,theta_e_mid,alpha_e_mid,vpar_e_mid,w_e_mid
  use magnetic_coordinates,only: mpoloidal,nflux,radcor_low1,radcor_upp1,radcor_low2,radcor_upp2,toroidal_range,nflux2,mtoroidal
  use magnetic_coordinates,only: nflux2,mtoroidal,tor_1d_array,radcor_1d_array,theta_1d_array,tor_shift_mc
  use array_in_mc,only: grad_psi_matrix
  use radial_module,only: baxis,r_axis
  use flux_tube_model,only:radcor_fixed,j_fixed
  use perturbation_field_matrix,only:mf_par_left,mf_x_left,mf_y_left,toroidal_mode_number_included,ey_left,den_left
  use perturbation_field_matrix,only:mf_par_old,mf_x_old,mf_y_old,ef_cyl_phi_left,potential
  use domain_decomposition,only: myid,numprocs,tube_comm,grid_comm,ntube,gclr,tclr,GCLR_cut,&
       & theta_interval,theta_start,my_right,my_left
  use pputil
  use mode_structure
  use deposit_ions_module,only: deposit_ions
  use transform_module,only: mode_evolution_analysis3
  use smoothing_module,only: smoothing_along_field_line_for_adiabatic_electron_model
  implicit none
  integer::   ierr
  real(p_):: dtao_e,t_omegai
  integer:: k,kt,i_new,iter

  real(p_):: b_si,omega_local,omega_i_axis,rho_i,k_binormal1,k_binormal2
  real(p_):: minor_r_radcor,gs_psi_prime,q_func2,pfn_func,radcor_as_func_of_pfn
  real(p_):: minor_r_min,minor_r_max,minor_r_width
  real(p_):: tarray(2) !store the cpu clock time
  INTEGER :: count1,count2, count3, count_rate, count_max
  real(p_):: z_old,theta_old,real_shift, real_shift_e, tor_shift_val

  integer:: np_old,np_new
  integer:: file_unit_i,file_unit_e,file_unit
integer :: u_evolution
  character(6):: filename_i,filename_e,filename

  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  call cpu_time(tarray(1))    !cpu_time is a f95 intrinsic subroutine

  call read_parameters()
  call ppinit(myid,numprocs,ntube,tube_comm,grid_comm)
  GCLR=INT(myid/ntube)
  TCLR=MOD(myid,ntube)

  theta_interval=twopi/(numprocs/ntube)  
  my_right=GCLR+1
  if(GCLR.eq.numprocs/ntube-1) my_right=0
  my_left=GCLR-1
  if(GCLR.eq.0) my_left=numprocs/ntube-1
  GCLR_cut=numprocs/ntube/2-1 !the value of GCLR at the theta cut (this case is at the high-field-side), GCLR=0 corresponds to low-field-side midplane
  theta_start=GCLR*theta_interval
  if(GCLR.gt.GCLR_cut) theta_start=-pi+(GCLR-GCLR_cut-1)*theta_interval

  !if(TCLR.eq.0) write(*,*) GCLR, theta_start
  call construct_numerical_tokamak_equilibrium()

  call cal_magnetic_coordinates(mpoloidal,nflux) !this to prepare numerical tables R(radcor,theta) and Z(radcor,theta)
  call calculate_metric(mpoloidal,nflux)
  call set_computational_radial_region()
  call mapping_cylindrical_to_magnetic_coordinates()  !this to prepare numerical tables radcor(R,Z), theta(R,Z), and tor_shift(R,Z). !Boris pusher works in cylindrical coordinates, we need these pre-mapping data to interpolate the Boris orbit into magnetic coordinates. Later I decide that tor_shift(R,Z) is not going to be used in the rest of the program. The value of tor_shift and the various gradients of psi, theta, and alpha will be computed by interpolating the numerical table in magnetic coordinates.
  !  if(myid.eq.0) call field_lines_analyse()
  call magnetic_coordinates_derivatives_in_cylindrical() !Later I decide that the data from this subroutine will not be used in the rest of the program. The reason is given above.
  call cal_array_in_magnetic_coordinates_for_guiding_center_pusher(mpoloidal,nflux) !numerical table (e.g., grad_alpha) in magnetic coordinates
  !  call cal_arc_length_along_field_line() !these arc lengths are needed when calculating the parallel derivative of perturbations, update: not used now

  call set_toroidal_range()
  call plasma_volume_of_computational_region()
  call plasma_volume_to_buffer()

  call load_ions() !spatial location in magnetic coordinates (psi,theta,phi) and then transform to cylindrical coordinates
  call ion_marker_radial_bdry_condition(nmarker_i,r_i,z_i,radcor_i,active_i,touch_bdry_i,w_i) !get the radial coordinates of markers, and place a marker at a random location within the buffer region if the marker is outside of the buffer region
  call ion_theta_and_alpha(nmarker_i,r_i,phi_i,z_i, radcor_i,active_i,touch_bdry_i,theta_i,alpha_i)

  call load_electrons() !in magnetic coordinates
  if(myid.eq.0) write(*,*) 'total number of electrons=',total_nmarker_e
  if(myid.eq.0) write(*,*) 'total number of ions=',total_nmarker_i

  call initial_weight_e()

  call create_toroidal_grids()
  call allocate_field_matrix()

!!$ call field_solver_at_full_step() 

  !omega_local=b_SI(r_i(1)*ln,z_i(1)*ln)*charge_i/mass_i
  omega_local=baxis*charge_i/mass_i
  !  call normalize_full_orbit_variables(nmarker_i,r_i,phi_i,z_i,vr_i,vphi_i,vz_i)
  !  dtao_i=twopi/omega_local/tn_i/8._p_/40._p_ !the time-step is chosen in terms of the local gyro-period
  !  dtao_i=twopi/omega_local/tn_i/64._p_ !the time-step is chosen in terms of the local gyro-period
  omega_i_axis=baxis*charge_i/mass_i
  rho_i=sqrt(ti0*kev/mass_i)/omega_i_axis
  if(myid.eq.0) write(*,*) 'omega_i_axis (kHz)=', omega_i_axis/1000._p_
  if(myid.eq.0) write(*,*) 'kappa_ti*rho_i=', kappa_ti*rho_i

  dtao_i=dtao_omega_i_axis/omega_i_axis/tn_i !time step in unit of tn_i
  minor_r_min=minor_r_radcor(radcor_low2)
  minor_r_max=minor_r_radcor(radcor_upp2)
  minor_r_width=minor_r_max-minor_r_min
  if(myid.eq.0) write(*,*) 'minor_r_min,minor_r_max,minor_r_width=',minor_r_min,minor_r_max,minor_r_width
  if(myid.eq.0) write(*,*) 'first radial sine harmonic: kr*rho_i=',pi/minor_r_width*rho_i
  !rho_i=vt_i/sqrt(2.0_p_)/omega_i_axis
  if(myid.eq.0) write(*,*) 'R0/rho_i=',r_axis/rho_i
  if(myid.eq.0) write(*,*) 'minor_r_width/rho_i=',minor_r_width/rho_i

  k_binormal1=toroidal_mode_number_included*q_func2(radcor_fixed)/minor_r_radcor(radcor_fixed)*rho_i
  k_binormal2=toroidal_mode_number_included*baxis/(gs_psi_prime(radcor_fixed)*&
       & (grad_psi_matrix(1,j_fixed)+grad_psi_matrix(mpoloidal/2,j_fixed))/two)*rho_i
  if(myid.eq.0) write(*,*) 'k_binorm*rhoi1=',k_binormal1,'k_binorm*rhoi2=',k_binormal2
  if(myid.eq.0) write(*,*) 'number of radial harmonics that should be included (pi*shear0*k_binormal/kr1)=',&
       & pi*0.78*k_binormal1/(pi/minor_r_width*rho_i)
  if(myid.eq.0) write(*,*) 'omega_star1/twopi (kHz)=', k_binormal1*kappa_ti*rho_i*omega_i_axis/twopi/1000._p_
  if(myid.eq.0) write(*,*) 'omega_star2/twopi (kHz)=', k_binormal2*kappa_ti*rho_i*omega_i_axis/twopi/1000._p_
  if(myid.eq.0) write(*,*) 'Ion beta=',ti0*kev*ni0/(baxis**2/(two*mu0)), 'ti0 (keV)=',ti0,'ni0 (m^-3)=',ni0

  vr_i_integer=vr_i !vr_i_integer is the projection of velocity at t_{n} to the basis vectors at t_{n}
  vz_i_integer=vz_i 
  vphi_i_integer=vphi_i 
  do k=1,nmarker_i !vr_i is initially the the projection of velocity at t_{0} to basis vectors at t_{0}, after the backward pushing, vr_i is the projection of velocity at t_{-1/2} to basis vector at t_{0}
     call backward_half_step_for_boris(dtao_i,r_i(k),z_i(k),phi_i(k),vr_i(k),vz_i(k),vphi_i(k)) !push only velocity, to set initial condition for the first step of boris algorithm, using mutiple steps (2nd RK)
     call forward_half_step_for_boris(dtao_i,r_i(k),z_i(k),phi_i(k),vr_i_integer(k),vz_i_integer(k),vphi_i_integer(k),&
          & r_i_mid(k),z_i_mid(k),phi_i_mid(k),vr_i_integer_mid(k),vz_i_integer_mid(k),vphi_i_integer_mid(k)) !push location half-step and then the velocity is projected onto the new local vector basis
  enddo

  !  open(21,file='full_orbit.txt')
  !  open(22,file='gc_electron.txt')
!!$  write(filename_i,'(a2,i4.4)') 'ib',myid
!!$  write(filename_e,'(a2,i4.4)') 'eb',myid
!!$  file_unit_i=myid+256
!!$  file_unit_e=myid+256*2
!!$  open(file_unit_i,file=filename_i)
  !!  open(file_unit_e,file=filename_e)

  dtao_e=dtao_i*(mass_i/mass_e) !time step in unit of electron gyro-period, tn_e=twopi/omegan_e
  if(myid.eq.0) write(*,*) 'dtao_i=',dtao_i,'dtao_i= (s)',dtao_i*tn_i,'dtao_e=',dtao_e,&
       &  'total_evolution_time (ms)=',dtao_i*tn_i*(kend-kstart+1)*1000.

  call initial_weight_i(w_i)
!!$  call ion_velocity_components_in_mc_at_integer_time_step()
  call prepare_coefficients_of_field_equations() ! time-step length dtao_i is involved in the coefficients, which are all independent of time
!!$  write(filename,'(a1,i4.4)') 'd',myid
!!$  file_unit=myid+118
!!$  open(file_unit,file=filename)

  !include implicit.in

  !--initial field--
  call sort_ions_according_to_poloidal_location(theta_i)
  call deposit_ions(nmarker_i,active_i,radcor_i,theta_i,alpha_i,vr_i_integer,vphi_i_integer,vz_i_integer,w_i)
  call adiabatic_electron_field_solver()
  kt=0
  if(TCLR.eq.0) then
!     call mode_structure_on_theta_plane(kt,GCLR,den_left,'den')
     call mode_structure_on_theta_plane(kt,GCLR,potential,'potential')
!     call mode_structure_on_poloidal_plane(kt,den_left,'den')
     call mode_structure_on_poloidal_plane(kt,potential,'potential')
!!$     call mode_structure_on_xz_plane(kt,den_left,'den')
     call mode_structure_on_xz_plane(kt,potential,'potential')
!!$     call mode_structure_on_yz_plane(kt,den_left,'den')
     call mode_structure_on_yz_plane(kt,potential,'potential')
  endif
  !--end initial field

if(myid.eq.2)open(newunit=u_evolution, file="evolution.txt")

  do kt=kstart,kend 
     !call some_test3(myid,numprocs)
     t_omegai=dtao_i*(kt-1)*tn_i*omega_i_axis
     if(myid.eq.1) call report(t_omegai)
     !if(myid.eq.2) call mode_evolution_analysis2(t_omegai)
     if(myid.eq.2) call mode_evolution_analysis3(t_omegai,potential,mtoroidal,nflux2,u_evolution)

     mf_x_old=mf_x_left
     mf_y_old=mf_y_left
     mf_par_old=mf_par_left

     r_i_old=r_i
     z_i_old=z_i
     phi_i_old=phi_i
     call store_old_electric_field()

     do iter=1,niter !if only two iterations are made, then this corresponds to a predictor-corrector scheme called Heun's method
        if(iter.eq.1) then !In linear case or nonlinear case with explicit trajectory pusher, the pusher needs to invoked only once.
           call push_ion_orbit_first_step(dtao_i) !get particle velocity at t_{n+1/2}, the location at t_{n+1} is also obtained
           r_i_mid=(r_i+r_i_old)/two
           z_i_mid=(z_i+z_i_old)/two
           phi_i_mid=(phi_i+phi_i_old)/two
           call ion_marker_radial_bdry_condition(nmarker_i,r_i_mid,z_i_mid,&
                & radcor_i_mid,active_i_mid,touch_bdry_i_mid,w_i_mid) !get the radial coordinates of markers, and place a marker at a random location within the buffer region if the marker is outside of the buffer region
           call ion_theta_and_alpha(nmarker_i,r_i_mid,phi_i_mid,z_i_mid, &
                & radcor_i_mid,active_i_mid,touch_bdry_i_mid,theta_i_mid,alpha_i_mid)
        endif

        call sort_ions_according_to_poloidal_location(theta_i_mid)
        call average_electric_field()
        call push_ion_weight_for_adiabatic_electron_model(dtao_i,radcor_i_mid,theta_i_mid,alpha_i_mid,&
             & r_i_mid,z_i_mid,phi_i_mid,vr_i_mid,vz_i_mid,vphi_i_mid,active_i_mid,w_i,w_i_star,nmarker_i)

!!$  write(filename,'(a1,i4.4)') 'i',myid
!!$  file_unit=myid+100
!!$  open(file_unit,file=filename)
!!$  do k=1,nmarker_i
!!$     write(file_unit,'(4(1pe14.5),i6,i3)')  radcor_i_mid(k) ,theta_i_mid(k),r_i_mid(k),z_i_mid(k),k,myid
!!$  enddo
!!$  close(file_unit)


!!$     call ion_velocity_components_in_mc_at_half_time_step()
        !     call sort_ions_according_to_poloidal_location(theta_i_mid)
        !     call deposit_ions(radcor_i_mid,theta_i_mid,alpha_i_mid,vr_i_mid,vphi_i_mid,vz_i_mid,w_i_star)
        !     call adiabatic_electron_field_solver()
!!$     if(fluid_electron.eqv..true.) then
!!$        call fluid_electron_pressure()
!!$     else
!!$        call push_electron_orbit_weight_half_step(dtao_e) !using field at t_{n} 
!!$        call deposit_electrons(nmarker_e,active_e,radcor_e_mid,theta_e_mid,alpha_e_mid,mu_e,vpar_e_mid,w_e_mid)!given the (v,x,w) of markers, do the deposition to get currents on grids
!!$     endif
!!$     call prepare_electron_source_terms()
!!$
!!$     call field_solver_electromagnetic_case(dtao_i*0.5_p_)
!!$
!!$     call evolve_parallel_magnetic_field(dtao_i*one_half)
!!$     call evolve_perpendicular_magnetic_field(dtao_i*one_half)
!!$     call smoothing_along_field_line() 
!!$     call communicate_field_value_between_neighbour_cells()
        !----------------------------second step--------------------------------------------------
!!$     call push_ion_weight_without_eper(dtao_i,radcor_i_mid,theta_i_mid,alpha_i_mid,w_i,w_i_star,nmarker_i) 
        !     call push_ion_weight_for_adiabatic_electron_model(dtao_i,radcor_i_mid,theta_i_mid,alpha_i_mid,&
        !          & r_i_mid,z_i_mid,phi_i_mid,vr_i_mid,vz_i_mid,vphi_i_mid,active_i,w_i,w_i_star,nmarker_i)
        !     w_i=w_i_star !update the weight

        !m if(iter.eq.1) call push_ion_orbit_second_step(dtao_i) !get velocity at t_{n+1} (and also the location at t_{n+3/2}, the location at t_{n+1} has already been calculated in "push_ion_orbit_first_step"
!!$     call ion_velocity_components_in_mc_at_integer_time_step()
        !yj     call sort_ions_according_to_poloidal_location(theta_i)
        !    call deposit_ions(radcor_i,theta_i,alpha_i,vr_i,vphi_i,vz_i,w_i_star)
        !    call adiabatic_electron_field_solver()
        if(iter.eq.1) call ion_marker_radial_bdry_condition(nmarker_i,r_i,z_i,radcor_i,active_i,touch_bdry_i,w_i) !get the radial coordinates of markers, and place a marker at a random location within the buffer region if the marker is outside of the buffer region
        if(iter.eq.1) call ion_theta_and_alpha(nmarker_i,r_i,phi_i,z_i, radcor_i,active_i,touch_bdry_i,theta_i,alpha_i)

        call sort_ions_according_to_poloidal_location(theta_i)
        call deposit_ions(nmarker_i,active_i,radcor_i,theta_i,alpha_i,vr_i_integer,vphi_i_integer,vz_i_integer,w_i_star)
        call adiabatic_electron_field_solver()

     enddo !iteration for the implicit weight pusher
     w_i=w_i_star !update the weight
     call smoothing_along_field_line_for_adiabatic_electron_model()

     !  write(file_unit,'(4(1pe14.4))') dtao_i*kt*tn_i*omega_i_axis,den_left(5,5),den_left(5,55),den_left(10,30)

!!$     if(fluid_electron.eqv..true.) then
!!$        call fluid_electron_pressure()
!!$     else
!!$        call push_electron_orbit_weight_full_step(dtao_e) !using values of (r,v,w) at t{n+1/2} to compute the rhs of evolution equation
!!$        call deposit_electrons(nmarker_e,active_e,radcor_e,theta_e,alpha_e,mu_e,vpar_e,w_e)
!!$     endif
!!$     call prepare_electron_source_terms()
!!$
!!$     call field_solver_electromagnetic_case(dtao_i)
!!$
!!$     mf_x_left=mf_x_old
!!$     mf_y_left=mf_y_old
!!$     mf_par_left=mf_par_old
!!$
!!$     call evolve_parallel_magnetic_field(dtao_i)
!!$     call evolve_perpendicular_magnetic_field(dtao_i)
!!$     call smoothing_along_field_line()
!!$     call communicate_field_value_between_neighbour_cells()

!!$     call update_ion_weight(dtao_i,w_i_star)
!!$
!!$     !call some_test2(myid,numprocs)
!!$     !call some_test4(myid,numprocs)
!!$     !call some_test5(myid,numprocs)
     if(TCLR.eq.0 .and. mod((kt-1),iplot_mode_structure).eq.0) then
!        call mode_structure_on_theta_plane(kt,GCLR,den_left,'den')
        call mode_structure_on_theta_plane(kt,GCLR,potential,'potential')
 !       call mode_structure_on_poloidal_plane(kt,den_left,'den')
        call mode_structure_on_poloidal_plane(kt,potential,'potential')
!!$        call mode_structure_on_xz_plane(kt,den_left,'den')
       call mode_structure_on_xz_plane(kt,potential,'potential')
!!$        call mode_structure_on_yz_plane(kt,den_left,'den')
        call mode_structure_on_yz_plane(kt,potential,'potential')
     endif
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  enddo !------------------------------------!main time advancing loop
if(myid.eq.2) close(u_evolution)
  !  close(file_unit)
  !close(file_unit_i)
  !  close(file_unit_e)
  !  close(21)
  !  close(22)

  call save_data_for_restarting()


!!$  write(*,*) 'myid=',myid,'ntouch_bdry_i,e=',ntouch_bdry_i,ntouch_bdry_e
!!$  call MPI_Allreduce( ntouch_bdry_i,  total_ntouch_bdry_i, 1, MPI_integer, MPI_sum, MPI_Comm_World, ierr)
!!$  call MPI_Allreduce( ntouch_bdry_e,  total_ntouch_bdry_e, 1, MPI_integer, MPI_sum, MPI_Comm_World, ierr)
!!$  if(myid.eq.0) write(*,*) 'total_ntouch_bdry_i,e in all processes=',total_ntouch_bdry_i,total_ntouch_bdry_e, ', fraction=',&
!!$       & real(total_ntouch_bdry_i)/total_nmarker_i,real(total_ntouch_bdry_e)/total_nmarker_e

  !for testing
  !     call field_solver() 

!!$  write(filename,'(a1,i4.4)') 'e',myid
!!$  file_unit=myid+100
!!$  open(file_unit,file=filename)
!!$  do k=1,nmarker_e
!!$     write(file_unit,*)  radcor_e(k) ,theta_e(k),rg_e(k),zg_e(k)
!!$  enddo
!!$  close(file_unit)


!!$  write(filename,'(a1,i4.4)') 'k',myid
!!$  file_unit=myid+100
!!$  open(file_unit,file=filename)
!!$  do k=1,nmarker_i
!!$    if(touch_bdry_i(k).eqv..true.)  write(file_unit,*)  radcor_i(k) ,theta_i(k),r_i(k),z_i(k)
!!$!     write(file_unit,*)  vr_i(k),vphi_i(k),vz_i(k)
!!$  enddo
!!$  close(file_unit)

  call cpu_time(tarray(2))
  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  if(myid.eq.0) write(*,*) 'CPU time used (seconds)', tarray(2)-tarray(1), &
       & 'Wall time used (seconds) ',  (count2-count1)/count_rate

  call MPI_FINALIZE(ierr)
end program lorentz_ions


