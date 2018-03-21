  subroutine read_parameters()
!    use mpi
    use control_parameters,only:kstart,kend,dtao_omega_i_axis,niter,ion_spatial_loading_scheme,ion_velocity_loading_scheme,&
         &  electron_spatial_loading_scheme,electron_velocity_loading_scheme,poloidal_angle_type,iplot_mode_structure
    use constants,only: twopi
    use normalizing,only:ln,bn,omegan_i,tn_i,vn_i,omegan_e,tn_e,vn_e
    use ions_module,only: mass_i,charge_i,ti0,ni0,kappa_ni,kappa_ti,total_nmarker_i ! total number of ion markers (including all the particles in all theprocessors).
    use electrons_module,only: mass_e,charge_e,te0,ne0,kappa_ne,kappa_te,total_nmarker_e,fluid_electron ! total number of ion markers (including all the particles in all theprocessors).
    use magnetic_coordinates,only:nflux,mpoloidal,mtoroidal,pfn_inner, pfn_bdry,nsegment
    use perturbation_field_matrix,only:filter_toroidal,filter_radial,toroidal_mode_number_included,radial_harmonics_included
    use domain_decomposition,only: ntube
    implicit none
!    integer,intent(in):: myid,numprocs
    integer:: ierr
    namelist/normalizing_nmlt/ln,bn
    namelist/control_nmlt/kstart,kend,dtao_omega_i_axis,niter,ion_spatial_loading_scheme,ion_velocity_loading_scheme,&
         & electron_spatial_loading_scheme,electron_velocity_loading_scheme,iplot_mode_structure,&
         & filter_toroidal,filter_radial,toroidal_mode_number_included,radial_harmonics_included,&
         & poloidal_angle_type,nsegment,nflux,mpoloidal,mtoroidal,pfn_inner, pfn_bdry,ntube,fluid_electron
    namelist/ions_nmlt/mass_i,charge_i,ni0,ti0,kappa_ni,kappa_ti,total_nmarker_i
    namelist/electrons_nmlt/mass_e,charge_e,te0,ne0,kappa_ne,kappa_te,total_nmarker_e

!    integer:: file_free
 !   integer:: status(MPI_STATUS_SIZE)
    ! ---a file can be safely accessed by only one process at a time, so the following codes gaurantee that the inputfile is read in turn by each processor
!!$    if ( myid .eq. 0 ) then ! master gets permission to read parameter file first, others wait in line
!!$       file_free = 1
!!$    else 
!!$       call MPI_Recv(file_free, 1, MPI_INT, myid-1, 1, MPI_COMM_WORLD, status,ierr)
!!$    endif

!    if (file_free .eq. 1) then ! this process read the input file
       open(31,file='input.nmlt')
       read(31,normalizing_nmlt)
       read(31,ions_nmlt)
       read(31,electrons_nmlt)
       read(31,control_nmlt)
       close(31)
 !   endif

  !  if (myid .ne. numprocs-1) call MPI_Send(file_free, 1, MPI_INT, myid+1, 1, MPI_COMM_WORLD,ierr)  !give reading file permission to the next process 

!!$    if(myid .eq.0)  then
!!$       !write(*,normalizing_nmlt)
!!$       !write(*,ions_nmlt)
!!$       !write(*,electrons_nmlt)
!!$    endif

    omegan_i=bn*charge_i/mass_i !cyclotron angular frequency in Hz
    tn_i=twopi/omegan_i !time unit used in this program
    vn_i=Ln/tn_i !the value of the normalizing velocity in SI unit m/s

    omegan_e=bn*charge_e/mass_e !cyclotron angular frequency in Hz
    tn_e=twopi/omegan_e
    vn_e=ln/tn_e
  end subroutine read_parameters




