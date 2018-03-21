subroutine load_ions()
  use precision,only:p_
  use constants,only:one,twopi,pi,two,kev,fourpi
  use control_parameters,only: ion_spatial_loading_scheme  !1=>uniform in (psi,theta,alpha) coordinates; 2=>uniform in real space
  use control_parameters,only: ion_velocity_loading_scheme !1=>uniform in (v,theta_v,phi_v) coordinates; 2=>Isotropic Gaussian in (v,theta_v,phi_v)
  use normalizing,only: vn_i
  use magnetic_coordinates,only: radcor_low=>radcor_low2,radcor_upp=>radcor_upp2,vol=>vol2 !as input
  use magnetic_coordinates,only: mpoloidal,nflux,radcor_1d_array,theta_1d_array,jacobian,toroidal_range
  use ions_module,only: total_nmarker_i,fixed_large_size,mass_i,ti0 !as input
  use ions_module,only: nmarker_i,radcor_i,theta_i,active_i ! as output
  use ions_module,only:active_i_mid,touch_bdry_i_mid,touch_bdry_i
  use ions_module,only: alpha_i,tor_shift_i !only allocate these array, their values are not assigned in the present subroutine
  use ions_module,only: r_i,z_i,phi_i,v_i,vr_i,vz_i,vphi_i !as output, loading using magnetic coordinates, the Boris pusher works in cylindrical coordinates, therefore, we need to transform from mag. cor. to cylin. cor.
  use ions_module,only: r_i_old,z_i_old,phi_i_old
  use ions_module,only: ps_vol_i,normalizing_factor !as output
  use ions_module,only: w_i,w_i_mid,w_i_star !only allocate the array, the value is not set in this subroutine
  use ions_module,only: vpar_i,vx_i,vy_i, grad_psi_i,grad_alpha_i,grad_psi_dot_grad_alpha_i,bval_i !only allocate the array, the value is not set in this subroutine
  use ions_module,only: r_i_mid,z_i_mid,phi_i_mid,radcor_i_mid,theta_i_mid,alpha_i_mid !only allocate the array, the value is not set in this subroutine
  use ions_module,only: vr_i_mid,vz_i_mid,vphi_i_mid !only allocate the array, the value is not set in this subroutine
  use ions_module,only: vr_i_integer,vz_i_integer,vphi_i_integer !only allocate the array, the value is not set in this subroutine
  use ions_module,only: vr_i_integer_mid,vz_i_integer_mid,vphi_i_integer_mid
  use ions_module,only: vt_i,vmin_i,vmax_i !as output, vt_i in SI unit
  use pputil !containing subroutines that sort particles into differnt processors
  use domain_decomposition,only: numprocs,myid
  use interpolate_module
  implicit none

  integer:: iseed,next_seed
  integer,parameter:: max_try=10000

  real(p_):: radcor_val,theta_val,rannum1,rannum2,rannum3,tmp
  !  real(p_):: random_yj
  integer:: i,ierr,j,file_unit
  !  integer:: status(MPI_STATUS_SIZE)
  character(5):: filename
  real(p_):: pos1,pos2,jacobian_val
  real(p_):: jacobian_func_normalized !func name
  integer:: np_old,np_new
  real(p_):: vt,vmin,vmax,v_val,maxwellian_func_ion
  !real(p_),allocatable:: theta_v(:),phi_v(:)
  real(p_),allocatable::vx(:),vy(:),vz(:),tmp_array(:)
  real(p_)::maxwellian_min,maxwellian_max

  nmarker_i=total_nmarker_i/numprocs !nmarker_i initially store the number of markers initially loaded per processor (i.e.total_nmarker_i/numprocs), latter actual number of markers per proc will be assigned to nmarker_i, the value of which will be differnt for differnt processors and at differnt time
  fixed_large_size=(total_nmarker_i/numprocs)*2 !the number of particle per proc after re-arranging the particles between the processors may exceed the number of original loaded particles per proc (i.e., total_nmarker_i/numprocs), increasing the array length by a factor of 2 is needed to make sure that the array is big enough to contain all the particles that belong to the domain for which the processor is responsible.
  allocate(radcor_i(fixed_large_size)) 
  allocate(theta_i(fixed_large_size))
  allocate(alpha_i(fixed_large_size))  !only allocate the array,the values are set in other places
  allocate(tor_shift_i(fixed_large_size))

  allocate(radcor_i_mid(fixed_large_size)) !only allocate the array, the values are set in other places
  allocate(theta_i_mid(fixed_large_size))!only allocate the array
  allocate(alpha_i_mid(fixed_large_size))  !only allocate the array

  allocate(v_i(fixed_large_size))
  allocate(vr_i(fixed_large_size))
  allocate(vz_i(fixed_large_size))
  allocate(vphi_i(fixed_large_size))

  !  allocate(pitch_angle_i(fixed_large_size))
  !  allocate(gyro_angle_i(fixed_large_size))

  allocate(r_i(fixed_large_size))
  allocate(z_i(fixed_large_size))
  allocate(phi_i(fixed_large_size))
  allocate(r_i_mid(fixed_large_size)) !only allocate the array, the values are set in other places
  allocate(z_i_mid(fixed_large_size)) !only allocate the array, the values are set in other places
  allocate(phi_i_mid(fixed_large_size))!only allocate the array, the values are set in other places

  allocate(w_i(fixed_large_size)) !only allocate the array, the initialization is done in another subroutine
  allocate(w_i_mid(fixed_large_size)) !only allocate the array, the initialization is done in another subroutine
  allocate(w_i_star(fixed_large_size)) !only allocate the array, the initialization is done in another subroutine

  allocate(ps_vol_i(fixed_large_size))
  allocate(active_i(fixed_large_size)) !whether particles are within computational boundary
  allocate(active_i_mid(fixed_large_size)) !whether particles are within computational boundary
  allocate(touch_bdry_i(fixed_large_size)) !whether particles are within computational boundary
  allocate(touch_bdry_i_mid(fixed_large_size)) !whether particles are within computational boundary

  !  allocate(theta_v(nmarker_i)) !local array
  !  allocate(phi_v(nmarker_i)) !local array
  allocate(vx(nmarker_i)) !local array
  allocate(vy(nmarker_i)) !local array
  allocate(vz(nmarker_i)) !local array
  allocate(tmp_array(3*nmarker_i))

  !  radcor_min=minval(radcor_1d_array)
  !  radcor_max=maxval(radcor_1d_array)


  ! ---random generator, when use MPI_send to generate iseed for other processes, it is actual a sequence generator,instead of parallel generator
!!$  if ( myid .eq. 0 ) then ! master generates random numbers first, others wait in line
!!$     iseed = 0
!!$  else 
!!$     call MPI_Recv(iseed, 1, MPI_INT, myid-1, 1, MPI_COMM_WORLD, status,ierr) !other processes wait to receive the iseed
!!$  endif

  iseed=-(1777+myid*3) !set the iseed in different procs, when using this, it is a parallel generator, but the random numbers in different procs may be related if the iseed chosen for differnt procs is not good enough
  !  write(*,*) 'myid=',myid, 'iseed=',iseed

  ! now generate the random numbers
  call sub_random_yj(iseed,next_seed,tmp) !just to trigger the use of the iseed, the generated random number is not used, 

  if(ion_spatial_loading_scheme.eq.1) then
     do i=1,nmarker_i
        call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
        call sub_random_yj(0,next_seed,rannum2) !use last random number as iseed
        radcor_val=radcor_low+(radcor_upp-radcor_low)*rannum1 !scale the random number to the range [radcor_low: radcor_upp]
        !theta_val=rannum2*twopi
        theta_val=(rannum2-0.5_p_)*twopi
        radcor_i(i)=radcor_val
        theta_i(i)=theta_val
     enddo
  elseif (ion_spatial_loading_scheme.eq.2) then

     do i=1,nmarker_i     
        do j=1,max_try !rejection method to generate nonuniform random numbers
           call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
           call sub_random_yj(0,next_seed,rannum2) !use last random number as iseed
           call sub_random_yj(0,next_seed,rannum3) !use last random number as iseed
           radcor_val=radcor_low+(radcor_upp-radcor_low)*rannum1 !scale the random number to the range [radcor_low: radcor_upp]
           theta_val=-pi+rannum2*twopi
           pos1=jacobian_func_normalized(theta_val,radcor_val)
           !       write(*,*) 'jacobian_func_normalized(theta_val,radcor_val)=',pos1
           call sub_random_yj(0,next_seed,pos2) !use last random number as iseed   
           if(pos1<pos2) then
              cycle
           else
              radcor_i(i)=radcor_val
              theta_i(i)=theta_val
              phi_i(i)=toroidal_range*rannum3
              exit
           endif
        enddo
        !     if(j.eq.max_try+1) stop "***stop**, rejection method is not successful in generating distribution"
        !     write(*,*) 'j=',j
     enddo
  else
     stop 'please specify a loading scheme for the spatial distribution ion markers'
  endif

  do i=1,nmarker_i
     call magnetic_coordinates_to_cylindrical_coordinates(theta_i(i),radcor_i(i),r_i(i),z_i(i)) !to get the corresponding (R,Z) coordinates
  enddo

!!$  do i=1,nmarker_i !setting toroidal coordinate of particles
!!$     call sub_random_yj(0,next_seed,rannum1) !use last random number as iseed
!!$     phi_i(i)=toroidal_range*rannum1
!!$  enddo

  !setting velocity
  vt=sqrt(two*ti0*kev/mass_i)
  vmin=-3._p_*vt/vn_i !normalized by vn_i
  vmax=3._p_*vt/vn_i !normalized by vn_i
  maxwellian_max=maxwellian_func_ion(0._p_)

!if(myid.eq.0) write(*,*)  'vt=',vt, 'vmin=',vmin,'vmax=',vmax, 'vmin*vn_i=',vmin*vn_i, 'vt/vn_i=',vt/vn_i

  if(ion_velocity_loading_scheme.eq.1) then !using uniform loading in v, instead of Gaussian
     do i=1,3*nmarker_i        
        call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
        v_val=vmin+rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
        tmp_array(i)=v_val
     enddo
  elseif (ion_velocity_loading_scheme.eq.2) then
     do i=1,3*nmarker_i        
        do j=1,max_try !rejection method to generate nonuniform random numbers
           call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
           v_val=vmin+rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
           pos1=maxwellian_func_ion(v_val*vn_i)
           call sub_random_yj(0,next_seed,pos2) !0 means using last random number as iseed   
           pos2=pos2*maxwellian_max !scaled to [0,maxwellian_max]
           if(pos1<pos2) then
              cycle
           else
              tmp_array(i)=v_val
              exit
           endif
        enddo
!        if(myid.eq.0) write(*,*) 'j=',j
     enddo
  else
     stop 'please specify a loading scheme for the velocity distribution ion markers'
  endif

  do i=1,nmarker_i
     vx(i)=tmp_array(i)
     vy(i)=tmp_array(i+nmarker_i)
     vz(i)=tmp_array(i+2*nmarker_i)
  v_i(i)=sqrt(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
  enddo



!if(myid.eq.0) call calculate_possibility_density(vz,nmarker_i,100,vmin,vmax)

  do i=1,nmarker_i
     vz_i(i)=vz(i)
     vr_i(i)=vx(i)*cos(phi_i(i))+vy(i)*sin(phi_i(i))
     vphi_i(i)=-vx(i)*sin(phi_i(i))+vy(i)*cos(phi_i(i))
  enddo

!  v_i=sqrt(vr_i*vr_i+vz_i*vz_i+vphi_i*vphi_i)
  !if(myid.eq.3) call calculate_possibility_density(v_i,nmarker_i,100,vmin,vmax)

  !  v_i=v_i/vn_i !normalized by vn_i
!!$  do i=1,nmarker_i !setting direction of velocity
!!$     call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed
!!$     theta_v(i)=pi*rannum1
!!$     call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed
!!$     phi_v(i)=twopi*rannum1
!!$  enddo
  !transform to components in cylindrical coordinates
!!$  do i=1,nmarker_i !velocity components in a constant Cartesian coordinate system
!!$     vz_i(i)=v_i(i)*cos(theta_v(i))
!!$     vx(i)=v_i(i)*sin(theta_v(i))*cos(phi_v(i))
!!$     vy(i)=v_i(i)*sin(theta_v(i))*sin(phi_v(i))
!!$  enddo

!!$  do i=1,nmarker_i !projected onto the basis vectors of cylindrical coordinates
!!$     vr_i(i)=vx(i)*cos(phi_i(i))+vy(i)*sin(phi_i(i))
!!$     vphi_i(i)=vy(i)*cos(phi_i(i))-vx(i)*sin(phi_i(i))
!!$  enddo

  if(ion_spatial_loading_scheme.eq.1 .and. ion_velocity_loading_scheme.eq.1) then
     !normalizing_factor=total_nmarker_i/(twopi*toroidal_range*(radcor_upp-radcor_low)*twopi*pi*(vmax-vmin))
     normalizing_factor=total_nmarker_i/(twopi*toroidal_range*(radcor_upp-radcor_low)*(vmax-vmin)**3)
     do i=1,nmarker_i
        call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,jacobian,theta_i(i),radcor_i(i),jacobian_val)
        !ps_vol_i(i)=v_i(i)**2*sin(theta_v(i))*abs(jacobian_val)/(normalizing_factor)
        ps_vol_i(i)=abs(jacobian_val)/(normalizing_factor)
     enddo
  elseif(ion_spatial_loading_scheme.eq.1 .and. ion_velocity_loading_scheme.eq.2) then
     !normalizing_factor=total_nmarker_i/((radcor_upp-radcor_low)*twopi*toroidal_range*twopi*pi*vt/vn_i*sqrt(pi)/two)
     normalizing_factor=total_nmarker_i/((radcor_upp-radcor_low)*twopi*toroidal_range*(sqrt(pi)*vt/vn_i)**3)
     do i=1,nmarker_i
        call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,jacobian,theta_i(i),radcor_i(i),jacobian_val)
        !ps_vol_i(i)=v_i(i)**2*sin(theta_v(i))*abs(jacobian_val)/(normalizing_factor*maxwellian_func_ion(v_i(i)*vn_i))
        ps_vol_i(i)=abs(jacobian_val)/(normalizing_factor*maxwellian_func_ion(v_i(i)*vn_i))
     enddo

  elseif(ion_spatial_loading_scheme.eq.2 .and. ion_velocity_loading_scheme.eq.1) then
     !normalizing_factor=total_nmarker_i/(vol*twopi*pi*(vmax-vmin))
     normalizing_factor=total_nmarker_i/(vol*(vmax-vmin)**3)
     do i=1,nmarker_i
        !ps_vol_i(i)=v_i(i)**2*sin(theta_v(i))/(normalizing_factor)
        ps_vol_i(i)=one/(normalizing_factor)
     enddo

  elseif(ion_spatial_loading_scheme.eq.2 .and. ion_velocity_loading_scheme.eq.2) then
     !  normalizing_factor=total_nmarker_i/(vol*fourpi*vt/vn_i*sqrt(pi)/two) !wrong
     !normalizing_factor=total_nmarker_i/(vol*twopi*pi*vt/vn_i*sqrt(pi)/two) !wrong again
     normalizing_factor=total_nmarker_i/(vol*(sqrt(pi)*vt/vn_i)**3) !corrected
     do i=1,nmarker_i
        !ps_vol_i(i)=v_i(i)**2*sin(theta_v(i))/(normalizing_factor*maxwellian_func_ion(v_i(i)*vn_i)) !wrong
        ps_vol_i(i)=one/(normalizing_factor*maxwellian_func_ion(v_i(i)*vn_i))
     enddo
  endif


!!$   iseed=next_seed
!!$  if (myid .ne. numprocs-1) then
!!$     call MPI_Send(iseed, 1, MPI_INT, myid+1, 1, MPI_COMM_WORLD,ierr)  !send the iseed to next process
!!$  endif

  !assign the loaded particles to the corresponding processors, using the subroutines provided in pputil_yj.f90
  np_old=nmarker_i
  call init_pmove(theta_i(:),np_old,twopi,ierr)
  call pmove(theta_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(radcor_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(r_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(phi_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(v_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vr_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vz_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vphi_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(ps_vol_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call end_pmove(ierr)
  nmarker_i=np_new

  !  call some_test3(myid,numprocs)

!!$  write(filename,'(a1,i4.4)') 'i',myid
!!$  file_unit=myid+100
!!$  open(file_unit,file=filename)
!!$  do i=1,nmarker_i
!!$     write(file_unit,'(4(1pe14.5),i6,i3)')  radcor_i(i) ,theta_i(i),r_i(i),z_i(i),i,myid
!!$     !     write(file_unit,*)  vr_i(i),vphi_i(i),vz_i(i)
!!$  enddo
!!$  close(file_unit)

  active_i=.true. ! initially, all markers are active, i.e., within the computational region

  allocate(vr_i_integer(fixed_large_size)) !only allocate, value is not set in this subroutine
  allocate(vz_i_integer(fixed_large_size)) !only allocate, value is not set in this subroutine
  allocate(vphi_i_integer(fixed_large_size)) !only allocate, value is not set in this subroutine

  allocate(vr_i_integer_mid(fixed_large_size)) !only allocate, value is not set in this subroutine
  allocate(vz_i_integer_mid(fixed_large_size)) !only allocate, value is not set in this subroutine
  allocate(vphi_i_integer_mid(fixed_large_size)) !only allocate, value is not set in this subroutine

  allocate(vr_i_mid(fixed_large_size)) !only allocate, value is not set in this subroutine
  allocate(vz_i_mid(fixed_large_size)) !only allocate, value is not set in this subroutine
  allocate(vphi_i_mid(fixed_large_size)) !only allocate, value is not set in this subroutine

  allocate(vpar_i(fixed_large_size)) !velocity components in magnetic coordinates,  only allocate, value is not set in this subroutine
  allocate(vx_i(fixed_large_size)) !only allocate, value is not set in this subroutine
  allocate(vy_i(fixed_large_size)) !only allocate, value is not set in this subroutine
  allocate(grad_psi_i(fixed_large_size)) !only allocate, value is not set in this subroutine
  allocate(grad_alpha_i(fixed_large_size)) !only allocate, value is not set in this subroutine
  allocate(grad_psi_dot_grad_alpha_i(fixed_large_size)) !only allocate, value is not set in this subroutine
  allocate(bval_i(fixed_large_size)) !only allocate, value is not set in this subroutine

  allocate(r_i_old(fixed_large_size)) !only allocate, value is not set in this subroutine
  allocate(z_i_old(fixed_large_size)) !only allocate, value is not set in this subroutine
  allocate(phi_i_old(fixed_large_size)) !only allocate, value is not set in this subroutine


  vmin_i=vmin
  vmax_i=vmax
  vt_i=vt
end subroutine load_ions


function maxwellian_func_ion(v) result(z)
use precision,only:p_
  use constants,only:two,kev
  use ions_module,only: mass_i,ti0 !as input
implicit none
real(p_):: v,z

z=exp(-mass_i*v*v/(two*ti0*kev)) !the normalizing factor (mi/(twopi*Ti*kev))^(3/2) is not included


end function maxwellian_func_ion



function random_yj(seed) result (z) !return a random number uniform distribute in [0:1]
  !linear congruental method to genereate random number
  !This corresponds to Park and Miller's choice implemented using Schrange's algorithm to aviod integer overflow
  !refer to R. Fitzpatrick's book "Computational physics, an introduction course" for the details
  use precision,only:p_
  implicit none
  integer:: seed
  real(p_):: z 
!!$  integer, parameter:: a=106,c=1283,m=6075
!!$    z=mod(a*seed+c,m)
  integer, parameter:: a=16807,m=2147483647 !m=2^31-1,c=0, this choice is called Park and Miller method
  integer,parameter:: q=127773 !q=m/a
  integer,parameter:: r=2836 !r=mod(m,a)
  real(p_),parameter:: RANDMAX=2147483646._p_
  integer,save:: next=1
  !  write(*,*) m
  if (seed .ne.0) next=seed
  next=a*mod(next,q)-r*(next/q)
  if(next<0) next=next+m
  z=next/RANDMAX
end function random_yj


subroutine sub_random_yj(seed,next_seed,randnum)  !return a random number uniform distribute in [0:1] 
!modified form random_yj(), also return the seed for next generator which may be needed by another generator in another proc
  !linear congruental method to genereate random number
  !This corresponds to Park and Miller's choice implemented using Schrange's algorithm to aviod integer overflow
  !refer to R. Fitzpatrick's book "Computational physics, an introduction course" for the details
  use precision,only:p_
  implicit none
  integer,intent(in):: seed
  real(p_),intent(out):: randnum 
  integer,intent(out):: next_seed
  integer, parameter:: a=16807,m=2147483647 !m=2^31-1,c=0, this choice is called Park and Miller method
  integer,parameter:: q=127773 !q=m/a
  integer,parameter:: r=2836 !r=mod(m,a)
  real(p_),parameter:: RANDMAX=2147483646._p_
  integer,save:: next=1
  !  write(*,*) m
  if (seed .ne.0) next=seed
  next=a*mod(next,q)-r*(next/q)
  if(next<0) next=next+m
  randnum=next/RANDMAX
  next_seed=next
end subroutine sub_random_yj


subroutine calculate_possibility_density(v,total_number,sample_number,starting_value,ending_value)
  use  precision,only:p_
  use  constants,only:one,two,four,five,twopi,eight
  implicit none

  integer,intent(in):: total_number,sample_number
  real(p_),intent(in):: v(total_number)
  real(p_),intent(in):: starting_value,ending_value
  real(p_):: xcenter(sample_number-1)
  real(p_):: possibility_density(sample_number-1)
  real(p_):: bin_npt(sample_number-1)
  real(p_):: interval
  real(p_):: x(sample_number)
  integer:: i,j

  interval=(ending_value-starting_value)/(sample_number-1)
  do i=1,sample_number
     x(i)=starting_value+interval*(i-1)
  enddo

  bin_npt=0
  do i=1,sample_number-1
     do j=1,total_number
        if (v(j).ge.x(i) .and. v(j).le.x(i+1)) bin_npt(i)=bin_npt(i)+1
     enddo
  enddo

  do i=1,sample_number-1
     xcenter(i)=(x(i)+x(i+1))/two
     possibility_density(i)=bin_npt(i)/(total_number)/interval
     write(*,*) xcenter(i), possibility_density(i)
  enddo

end subroutine calculate_possibility_density



