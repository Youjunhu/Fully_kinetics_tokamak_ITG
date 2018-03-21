subroutine load_electrons ()
  use precision,only:p_
  use constants,only:twopi,pi,two,kev,fourpi
  use normalizing,only: Ln,vn_e
  use control_parameters,only: electron_spatial_loading_scheme  !1=>uniform in (psi,theta,alpha) coordinates; 2=>uniform in real space
  use control_parameters,only: electron_velocity_loading_scheme !1=>uniform in (v,theta_v,phi_v) coordinates; 2=>Isotropic Gaussian in (v,theta_v,phi_v)
  use magnetic_coordinates,only: radcor_low2,radcor_upp2,vol2 !as input
  use magnetic_coordinates,only: mpoloidal,nflux,jacobian,radcor_1d_array,theta_1d_array,toroidal_range
  use electrons_module,only: total_nmarker_e,mass_e,te0 !as input
  use electrons_module,only: nmarker_e,radcor_e,theta_e,alpha_e,touch_bdry_e,active_e ! as output
  use electrons_module,only: tor_shift_e !only allocate these array, their values are not assigned in the present subroutine
  use electrons_module,only: rg_e,zg_e,phig_e !as output, loading using magnetic coordinates, the Boris pusher works in cylindrical coordinates, therefore, we need to transform from mag. cor. to cylin. cor.
  use electrons_module,only: vpar_e,mu_e
  use electrons_module,only:ps_vol_e,w_e
  use electrons_module,only:radcor_e_mid,theta_e_mid,alpha_e_mid,vpar_e_mid,w_e_mid !vlaues at half-time-step
  use magnetic_coordinates,only: mpoloidal,nflux,theta_1d_array,radcor_1d_array,tor_shift_mc
  use domain_decomposition,only: numprocs,myid
  use interpolate_module
  use pputil !containing subroutines that sort particles into differnt processors
  implicit none
  integer:: fixed_large_size
  integer:: iseed,next_seed
  real(p_),allocatable:: v_e(:),vper_e(:) !used only when doing initial loading

  integer,parameter:: max_try=10000 !used in rejection method

  real(p_):: radcor_val,theta_val,rannum1,rannum2,tmp
  !  real(p_):: random_yj
  integer:: i,ierr,j,file_unit
!  integer:: status(MPI_STATUS_SIZE)
  character(5):: filename
  real(p_):: pos1,pos2
  real(p_):: jacobian_func_normalized !func name
  real(p_):: jacobian_val
  integer:: np_old,np_new
  real(p_):: vt,vmin,vmax,v_val,maxwellian_func_electron
  real(p_),allocatable:: theta_v(:),phi_v(:)

  real(p_)::maxwellian_min,maxwellian_max
  real(p_):: br_mc_func,bz_mc_func,bphi_mc_func !br_mc_func is function name, normalized br as a function of magnetic coordinates (theta,radcor)
  real(p_):: bval,bxval,byval,bzval
  real(p_) vx,vy,vz
  real(p_):: normalizing_factor

  nmarker_e=total_nmarker_e/numprocs !nmarker_e initially store the number of markers initially loaded per processor (i.e.total_nmarker_e/numprocs), latter actual number of markers per proc will be assigned to nmarker_e, the value of which will be differnt for differnt processors and at differnt time
  fixed_large_size=(total_nmarker_e/numprocs)*2 !the number of particle per proc after re-arranging the particles between the processors may exceed the number of original loaded particles per proc (i.e., total_nmarker_e/numprocs), increasing the array length by a factor of 2 is needed to make sure that the array is big enough to contain all the particles that belong to the domain for which the processor is responsible.
  allocate(radcor_e(fixed_large_size)) 
  allocate(theta_e(fixed_large_size))
  allocate(alpha_e(fixed_large_size))
  allocate(tor_shift_e(fixed_large_size))
  allocate(v_e(fixed_large_size))
  allocate(vpar_e(fixed_large_size))
  allocate(vper_e(fixed_large_size))
  allocate(mu_e(fixed_large_size))
  allocate(w_e(fixed_large_size)) !only allocate the array, its value is set in other subroutines
  allocate(ps_vol_e(fixed_large_size))

  allocate(radcor_e_mid(fixed_large_size)) !only allocate the array, the values are set in subroutine "push_electron_orbit_weight_half_step"
  allocate(theta_e_mid(fixed_large_size))!only allocate the array
  allocate(alpha_e_mid(fixed_large_size))  !only allocate the array
  allocate(vpar_e_mid(fixed_large_size))
  allocate(w_e_mid(fixed_large_size)) !only allocate the array, its value is set in other subroutines

  !  allocate(pitch_angle_e(fixed_large_size))
  !  allocate(gyro_angle_e(fixed_large_size))

  allocate(rg_e(fixed_large_size))
  allocate(zg_e(fixed_large_size))
  allocate(phig_e(fixed_large_size))

  allocate(theta_v(nmarker_e)) !local array
  allocate(phi_v(nmarker_e)) !local array

  allocate(touch_bdry_e(fixed_large_size))
  allocate(active_e(fixed_large_size))
  touch_bdry_e=.false. ! initially, all markers are considered as not touching the boundary
  active_e=.true. !initially, all markers are within the computational boundary

  !  radcor_min=minval(radcor_1d_array)
  !  radcor_max=maxval(radcor_1d_array)


  ! ---random generator, when use MPI_send to generate iseed for other processes, it is actual a sequence generator,instead of parallel generator
!!$  if ( myid .eq. 0 ) then ! master generates random numbers first, others wait in line
!!$     iseed = 0
!!$  else 
!!$     call MPI_Recv(iseed, 1, MPI_INT, myid-1, 1, MPI_COMM_WORLD, status,ierr) !other processes wait to receive the iseed
!!$  endif

  iseed=-(2777+myid*3) !set the iseed in different procs, when using this, it is a parallel generator, but the random numbers in different procs may be related if the iseed chosen for differnt procs is not good enough
  !  write(*,*) 'myid=',myid, 'iseed=',iseed

  ! now generate the random numbers
  call sub_random_yj(iseed,next_seed,tmp) !just to trigger the use of the iseed, the generated random number is not used, 

  if(electron_spatial_loading_scheme.eq.2) then
     do i=1,nmarker_e
        do j=1,max_try !rejection method to generate nonuniform random numbers
           call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
           call sub_random_yj(0,next_seed,rannum2) !use last random number as iseed
           radcor_val=radcor_low2+(radcor_upp2-radcor_low2)*rannum1 !scale the random number to the range [radcor_low2: radcor_upp2]
           !theta_val=rannum2*twopi
           theta_val=(rannum2-0.5_p_)*twopi
           pos1=jacobian_func_normalized(theta_val,radcor_val)
           !       write(*,*) 'jacobian_func_normalized(theta_val,radcor_val)=',pos1
           call sub_random_yj(0,next_seed,pos2) !use last random number as iseed   
           if(pos1<pos2) then
              cycle
           else
              radcor_e(i)=radcor_val
              theta_e(i)=theta_val
              exit
           endif
        enddo !uniform space distiribution, which non-uniform in (radcor,theta) coordinates
        !     if(j.eq.max_try+1) stop "***stop**, rejection method is not successful in generating distribution"
        !     write(*,*) 'j=',j
     enddo
  elseif(electron_spatial_loading_scheme.eq.1) then    !use uniform loading in (radcor,theta), the testing indicate this loading can give more accurate Monte-Color integration than uniform loading in real space
     do i=1,nmarker_e     
        call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
        call sub_random_yj(0,next_seed,rannum2) !use last random number as iseed
        radcor_val=radcor_low2+(radcor_upp2-radcor_low2)*rannum1 !scale the random number to the range [radcor_low2: radcor_upp2]
        !theta_val=rannum2*twopi
        theta_val=(rannum2-0.5_p_)*twopi
        radcor_e(i)=radcor_val
        theta_e(i)=theta_val
     enddo
  else
     stop 'please specify a loading scheme for the spatial distribution of electron markers'
  endif


  do i=1,nmarker_e
     call magnetic_coordinates_to_cylindrical_coordinates(theta_e(i),radcor_e(i),rg_e(i),zg_e(i)) !to get the corresponding (R,Z) coordinates
  enddo

  rg_e=rg_e/Ln !changed from SI unit to Ln unit
  zg_e=zg_e/Ln !changed from SI unit to Ln unit


  do i=1,nmarker_e !setting toroidal coordinate of particles
     call sub_random_yj(0,next_seed,rannum1) !use last random number as iseed
     !phig_e(i)=twopi*rannum1
phig_e(i)=toroidal_range*rannum1
  enddo

  do i=1,nmarker_e
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,tor_shift_mc,theta_e(i),radcor_e(i),tor_shift_e(i)) 
     alpha_e(i)=phig_e(i)-tor_shift_e(i)
!!$     alpha_e(i)=alpha_e(i)-int(alpha_e(i)/twopi)*twopi !shift into the range [0:twopi]
!!$     if(alpha_e(i).lt.0) alpha_e(i)=alpha_e(i)+twopi !shift into the range [0:twopi]
!    call shift_to_zero_twopi_range(alpha_e(i))
         call shift_to_specified_toroidal_range(alpha_e(i))
  enddo

  !setting velocity
  vt=sqrt(two*te0*kev/mass_e)
  vmin=0._p_
  vmax=3._p_*vt
  maxwellian_min=maxwellian_func_electron(vmax)
  maxwellian_max=maxwellian_func_electron(vmin)
  vmin=vmin/vn_e
  vmax=vmax/vn_e

  if(electron_velocity_loading_scheme.eq.1) then
     do i=1,nmarker_e
        call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
        v_val=vmin+rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
        v_e(i)=v_val
     enddo
  elseif(electron_velocity_loading_scheme.eq.2) then
     do i=1,nmarker_e
        do j=1,max_try !rejection method to generate nonuniform random numbers
           call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
           v_val=vmin+rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
           pos1=maxwellian_func_electron(v_val*vn_e)
           call sub_random_yj(0,next_seed,pos2) !0 means using last random number as iseed   
           pos2=maxwellian_min+pos2*(maxwellian_max-maxwellian_min) !scaled to [maxwellian_min,maxwellian_max]
           if(pos1<pos2) then
              cycle
           else
              v_e(i)=v_val
              exit
           endif
        enddo
     enddo
  else
     stop 'please specify a loading scheme for the velocity distribution of electron markers'
  endif


  !  v_e=v_e/vn_e !normalized by vn_e, vn_e=ln/tn_e, tn_e=twopi/omegan_e, omegan_e=bn*charge_e/mass_e

  do i=1,nmarker_e !setting direction of velocity, assuming a global spherical coordinates, !can not assume a local poloidal coordinates with the magnetic field defining the axis, because the Jacobian will be spatially dependendent in doing this, which is not desired
     call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed
     theta_v(i)=pi*rannum1
     call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed
     phi_v(i)=twopi*rannum1
  enddo


  do i=1,nmarker_e
     vz=v_e(i)*cos(theta_v(i)) !velocity components in a constant Cartesian coordinate system
     vx=v_e(i)*sin(theta_v(i))*cos(phi_v(i)) !the same as above
     vy=v_e(i)*sin(theta_v(i))*sin(phi_v(i)) !the same as above
     !gyro-radius is assumed to be nearly zero, and thus the particle location is approximated by the guiding-center location
     bxval=br_mc_func(theta_e(i),radcor_e(i))*cos(phig_e(i)) +bphi_mc_func(theta_e(i),radcor_e(i))*(-sin(phig_e(i))) !b components in a constant Cartesian coor. system
     byval=br_mc_func(theta_e(i),radcor_e(i))*sin(phig_e(i)) +bphi_mc_func(theta_e(i),radcor_e(i))*cos(phig_e(i)) !b components in a constant Cartesian coor. system
     bzval=bz_mc_func(theta_e(i),radcor_e(i)) !b components in a constant Cartesian coor. system
     bval=sqrt(bxval**2+byval**2+bzval**2) 
     vpar_e(i)=(vx*bxval+vy*byval+vz*bzval)/bval !scalar product
     vper_e(i)=sqrt((vy*bzval-byval*vz)**2+(vz*bxval-vx*bzval)**2+(vx*byval-bxval*vy)**2)/bval !cross product
     mu_e(i)=vper_e(i)**2/(two*bval) !normalized magnetic moment
  enddo

  if(electron_spatial_loading_scheme.eq.1 .and. electron_velocity_loading_scheme.eq.1) then
     normalizing_factor=total_nmarker_e/(twopi*toroidal_range*(radcor_upp2-radcor_low2)*twopi*pi*(vmax-vmin))
     do i=1,nmarker_e
        call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,jacobian,theta_e(i),radcor_e(i),jacobian_val)
        ps_vol_e(i)=v_e(i)**2*sin(theta_v(i))*abs(jacobian_val)/(normalizing_factor)
     enddo

  elseif(electron_spatial_loading_scheme.eq.1 .and. electron_velocity_loading_scheme.eq.2) then
     normalizing_factor=total_nmarker_e/(twopi*toroidal_range*(radcor_upp2-radcor_low2)*twopi*pi*vt/vn_e*sqrt(pi)/two)
     do i=1,nmarker_e
        call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,jacobian,theta_e(i),radcor_e(i),jacobian_val)
        ps_vol_e(i)=v_e(i)**2*sin(theta_v(i))*abs(jacobian_val)/(normalizing_factor*maxwellian_func_electron(v_e(i)*vn_e))
     enddo
  elseif(electron_spatial_loading_scheme.eq.2 .and. electron_velocity_loading_scheme.eq.1) then
   normalizing_factor=total_nmarker_e/(vol2*twopi*pi*(vmax-vmin)) 
     do i=1,nmarker_e
        ps_vol_e(i)=v_e(i)**2*sin(theta_v(i))/(normalizing_factor)
     enddo
  elseif(electron_spatial_loading_scheme.eq.2 .and. electron_velocity_loading_scheme.eq.2) then
     !  normalizing_factor=total_nmarker_e/(vol2*fourpi*vt/vn_e*sqrt(pi)/two) !wrong
     normalizing_factor=total_nmarker_e/(vol2*twopi*pi*vt/vn_e*sqrt(pi)/two) !for loading particles unifrom in real-space
     do i=1,nmarker_e
        ps_vol_e(i)=v_e(i)**2*sin(theta_v(i))/(normalizing_factor*maxwellian_func_electron(v_e(i)*vn_e))
     enddo
  endif

!!$   iseed=next_seed
!!$  if (myid .ne. numprocs-1) then
!!$     call MPI_Send(iseed, 1, MPI_INT, myid+1, 1, MPI_COMM_WORLD,ierr)  !send the iseed to next process
!!$  endif

  !assign the loaded particles to the corresponding processors, using the subroutines provided in pputil_yj.f90
  np_old=nmarker_e
  call init_pmove(theta_e(:),np_old,twopi,ierr)
  call pmove(theta_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(radcor_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(alpha_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(rg_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(zg_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(phig_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vpar_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(mu_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(ps_vol_e(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call end_pmove(ierr)
  nmarker_e=np_new

!!$  write(filename,'(a1,i4.4)') 'e',myid
!!$  file_unit=myid+100
!!$  open(file_unit,file=filename)
!!$  do i=1,nmarker_e
!!$     write(file_unit,*)  radcor_e(i) ,theta_e(i),rg_e(i),zg_e(i),i
!!$  enddo
!!$  close(file_unit)

!!$  !     call some_test2(myid,numprocs)
end subroutine load_electrons


function maxwellian_func_electron(v) result(z) !v in SI unit
use precision,only:p_
  use constants,only:two,kev
  use electrons_module,only: mass_e,te0 !as input
implicit none
real(p_):: v,z

z=exp(-mass_e*v*v/(two*te0*kev))
end function maxwellian_func_electron
