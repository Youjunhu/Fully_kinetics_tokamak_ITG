subroutine construct_numerical_tokamak_equilibrium()
  use mpi
  use precision,only:p_
  use constants,only:zero,one,two,three,four,five,twopi 
  use normalizing,only: Ln,bn
  use poloidal_flux_2d,only:xarray,zarray,nx,nz,psi,psi_gradient,psi_x,psi_z,psi_xx,psi_zz,psi_xz,psi_zx, &
       & y2a_psi,y2a_gradient,y2a_psi_x,y2a_psi_z,y2a_psi_xx,y2a_psi_zz,y2a_psi_xz,y2a_psi_zx !y2a_psi_* is a tempory array used in 2d cubic spline interpolation of psi_*, psi_x is the partial derivative with respect to x, and similar meaning for psi_x, psi_z etc.
  use radial_module,only:psi_axis,psi_lcfs,npsi,psi_1d,fpsi,ffprime,fprime,qpsi,pfn_npsi,tfn_npsi,baxis,r_axis,z_axis !!location of magnetic axis
  use boundary,only: x_lcfs,z_lcfs,np_lcfs
  use domain_decomposition,only:myid,numprocs
!  use magnetic_coordinates,only:mpoloidal
  !  use radial_module,only: y2_fpsi,y2_fprime

  implicit none
  character(100):: gfile_name
  namelist /magnetic_configuration/gfile_name
  integer,parameter:: nwmax=501 !assumed maximum length of the arrays, the actual length will be determined by nx, nz, or, np_lcfs
  real(p_):: tmp_psi(nwmax,nwmax),tmp_fpsi(nwmax),tmp_qpsi(nwmax),tmp_press(nwmax),tmp_pprime(nwmax),tmp_ffprime(nwmax)
  real(p_):: rleft,zmid,xdim,zdim !specification of the rectangular compuational within which the value of psi is known
  !  real(p_):: x_lcfs0(nwmax),z_lcfs0(nwmax) !the x and z coordinates of the points on the LCFS, nwmax is the assumed maximum length of the two arrays, the actual length will be determined by np_lcfs
  !  integer:: np_lcfs ! the actual length of arrays x_lcfs and z_lcfs

  real(p_):: psi_func !the interpolating poloidal flux function of two variables x and z (constructed by 2d cubic spline interpolation)
  real(p_):: psi_gradient_func,psi_x_func,psi_z_func,psi_xx_func,psi_zz_func,psi_xz_func,psi_zx_func !names of interpolating functions
  integer:: i,j,file_free,ierr
  integer:: status(MPI_STATUS_SIZE)

  ! ---a file can be safely accessed by only one process at a time, so the following codes gaurantee that the inputfile is read in turn by each processor
!!$  if ( myid .eq. 0 ) then ! master gets permission to read parameter file first, others wait in line
!!$     file_free = 1
!!$  else 
!!$     call MPI_Recv(file_free, 1, MPI_INT, myid-1, 1, MPI_COMM_WORLD, status,ierr)
!!$  endif

!  if (file_free .eq. 1) then ! this process read the input file
     open(31,file='input.nmlt')
     read(31,magnetic_configuration)
     close(31)
     call read_gfile(gfile_name,nwmax,tmp_psi,nx,nz,rleft,zmid,xdim,zdim,psi_axis,psi_lcfs, &
       & r_axis,z_axis,baxis,tmp_fpsi,tmp_qpsi,tmp_press,tmp_pprime,tmp_ffprime,myid)
     if(myid .eq.0)  write(*,magnetic_configuration)
!  endif

 ! if (myid .ne. numprocs-1) call MPI_Send(file_free, 1, MPI_INT, myid+1, 1, MPI_COMM_WORLD,ierr)  !give reading file permission to the next process 
!------------

!  write(*,*) '(R,Z) grids in G-file are ', 'nx=',nx,'nz=',nz
  !now the value of nx and nz is known, we allocate the arrays using the actual lenght of the corresponding array:
  !Note that nx in g-file is also used to define the number of radial grid points.
  allocate(psi(nx,nz))
  allocate(xarray(nx))
  allocate(zarray(nz))

  !now we use the actual length of the array, instead of the assumed maximum length
  do i=1,nx
     do j=1,nz
        psi(i,j)=tmp_psi(i,j) !use the actual length of the psi array, instead of the assumed maximum length
     enddo
  enddo

  do i=1,nx !construct the X array
     xarray(i)=rleft+xdim/(nx-1)*(i-1)
  enddo

  do j=1,nz !construct the Z array
     zarray(j)=(zmid-zdim/two)+zdim/(nz-1)*(j-1)
  enddo

  allocate(y2a_psi(nx,nz)) !this is an intermedial array needed in cubic spline interpolation.
  call splie2(xarray,zarray,psi,nx,nz,y2a_psi) !after this call, the function psi_func (defined later in this file) is ready to be used.

  allocate(psi_x(nx,nz))
  allocate(psi_z(nx,nz))
  allocate(psi_xx(nx,nz))
  allocate(psi_zz(nx,nz))
  allocate(psi_xz(nx,nz))
  allocate(psi_zx(nx,nz))
  allocate(psi_gradient(nx,nz))
  call calculate_poloidal_flux_partial_derivatives(nx,nz,xarray,zarray,psi,psi_x,psi_z,psi_xx,psi_zz,psi_xz,psi_zx,psi_gradient)
!!$  allocate(y2a_gradient(nx,nz)) ! an array in spline interpolation to store 2nd derivatives
!!$  allocate(y2a_psi_x(nx,nz))
!!$  allocate(y2a_psi_z(nx,nz))
!!$  allocate(y2a_psi_xx(nx,nz))
!!$  allocate(y2a_psi_zz(nx,nz))
!!$  allocate(y2a_psi_xz(nx,nz))
!!$  allocate(y2a_psi_zx(nx,nz))
!!$  call splie2(xarray,zarray,psi_gradient,nx,nz,y2a_gradient) !after this call, the function psi_gradient_func is ready to be used.
!!$  call splie2(xarray,zarray,psi_x,nx,nz,y2a_psi_x) !after this call, the function psi_x_func is ready to be used.
!!$  call splie2(xarray,zarray,psi_z,nx,nz,y2a_psi_z) !after this call, the function psi_z_func is ready to be used.
!!$  call splie2(xarray,zarray,psi_xx,nx,nz,y2a_psi_xx) !after this call, the function psi_xx_func is ready to be used.
!!$  call splie2(xarray,zarray,psi_zz,nx,nz,y2a_psi_zz) !after this call, the function psi_zz_func is ready to be used.
!!$  call splie2(xarray,zarray,psi_xz,nx,nz,y2a_psi_xz) !after this call, the function psi_xz_func is ready to be used.
!!$  call splie2(xarray,zarray,psi_zx,nx,nz,y2a_psi_zx) !after this call, the function psi_zx_func is ready to be used.

  !  write(*,*) 'minimum value of psi determined directly from the discrete psi data', minval(psi)
  ! write(*,*) 'value of psi at magnetic axis calculated from interpolating function',psi_func(r_axis,z_axis)
  !write(*,*) 'value of psi at magnetic axis specified in g-file, psi_axis=',psi_axis

  npsi=nx !nx in g-file is also used to define the number of radial grid points.
  allocate(qpsi(npsi))
  allocate(psi_1d(npsi))
  allocate(fpsi(npsi))
  allocate(ffprime(npsi))
  allocate(fprime(npsi))
  do i=1,npsi !uniform psi array, which is the assumed radial coordinator by G-file for magnetic surface function.
     psi_1d(i)=psi_axis+(psi_lcfs-psi_axis)/(npsi-1)*(i-1) 
  enddo

  do i=1,npsi !use the actual length of the radial array, instead of the assumed maximum length
     fpsi(i)=tmp_fpsi(i) 
     ffprime(i)=tmp_ffprime(i)
     qpsi(i)=tmp_qpsi(i)
  enddo
  fprime=ffprime/fpsi

  allocate(tfn_npsi(npsi))
  call calculate_tfn(npsi,psi_1d,qpsi,tfn_npsi)
  allocate(pfn_npsi(npsi))
  pfn_npsi=(psi_1d-psi_1d(1))/(psi_1d(npsi)-psi_1d(1))


!!$  allocate(y2_fpsi(npsi))
!!$  allocate(y2_fprime(npsi))
!!$  call spline(psi_1d,fpsi,npsi,2.d30,2.d30,y2_fpsi) !prepare the second order derivative needed in the cubic spline interpolation
!!$  call spline(psi_1d,fprime,npsi,2.d30,2.d30,y2_fprime) !prepare the second order derivative needed in the cubic spline interpolation
  call draw_rect_region(nx,nz,xarray,zarray)  !draw the compuational box used in G-file
  call draw_3d_tokamak()

  call arrange_lcfs(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis,myid)

end subroutine construct_numerical_tokamak_equilibrium


function psi_func(x,z) result(funcval)!SI units, poloidal magnetic flux function
  use precision,only:p_
  use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi !as input
  use interpolate_module,only: linear_2d_interpolation

  implicit none
  real(p_):: x,z,funcval
  call linear_2d_interpolation(nx,nz,xarray,zarray,psi,x,z,funcval)  
end function

function pfn_func(x,z) result(funcval)!SI units, poloidal magnetic flux function
  use precision,only:p_
  use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi !as input
  use radial_module,only:psi_axis,psi_lcfs
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  real(p_):: x,z,funcval
  call linear_2d_interpolation(nx,nz,xarray,zarray,psi,x,z,funcval)  
funcval=(funcval-psi_axis)/(psi_lcfs-psi_axis)
end function


function psi_func_spline(xval,zval) result(funcval)!SI units, poloidal magnetic flux function
  use precision,only:p_
  use poloidal_flux_2d,only: xarray,zarray,psi,y2a_psi,nx,nz
  implicit none
  real(p_)::xval,zval,funcval
  call splin2(xarray,zarray,psi,y2a_psi,nx,nz,xval,zval,funcval) !use spline interpolation
end function


function psi_r_func(x,z)
  use precision,only:p_
  use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_x
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  real(p_):: x,z,psi_r_func
  call linear_2d_interpolation(nx,nz,xarray,zarray,psi_x,x,z,psi_r_func)  
end function psi_r_func


function psi_z_func(x,z)
  use precision,only:p_
  use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_z !as input
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  real(p_):: x,z,psi_z_func
  call linear_2d_interpolation(nx,nz,xarray,zarray,psi_z,x,z,psi_z_func)  
end function psi_z_func


function psi_rr_func(x,z)
  use precision,only:p_
  use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_xx !as input
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  real(p_):: x,z,psi_rr_func
  call linear_2d_interpolation(nx,nz,xarray,zarray,psi_xx,x,z,psi_rr_func)  
end function psi_rr_func

function psi_zz_func(x,z)
  use precision,only:p_
  use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_zz !as input
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  real(p_):: x,z,psi_zz_func
  call linear_2d_interpolation(nx,nz,xarray,zarray,psi_zz,x,z,psi_zz_func)  
end function psi_zz_func


function psi_rz_func(x,z)
  use precision,only:p_
  use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_xz !as input
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  real(p_):: x,z,psi_rz_func
  call linear_2d_interpolation(nx,nz,xarray,zarray,psi_xz,x,z,psi_rz_func)  
end function psi_rz_func


function psi_zr_func(x,z)
  use precision,only:p_
  use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_zx !as input
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  real(p_):: x,z,psi_zr_func
  call linear_2d_interpolation(nx,nz,xarray,zarray,psi_zx,x,z,psi_zr_func)  
end function psi_zr_func


function q_func(psival) result(z) !safety factor
  !not used in computing the orbits, only used to do analytical estimation of some quantities, such as bounce frequency
  use precision,only:p_
  use radial_module,only:npsi,psi_1d,qpsi
  use interpolate_module,only: linear_1d_interpolation
  implicit none
  real(p_):: z,psival
  call linear_1d_interpolation(npsi,psi_1d,qpsi,psival,z)  
end function q_func


function q_func2(pfn) result(z) !safety factor, use pfn instead of psival as the independent variable
  !not used in computing the orbits, only used to do analytical estimation of some quantities, such as bounce frequency
  use precision,only:p_
  use radial_module,only:npsi,pfn_npsi,qpsi
  use interpolate_module,only: linear_1d_interpolation
  implicit none
  real(p_):: z,pfn
  call linear_1d_interpolation(npsi,pfn_npsi,qpsi,pfn,z)  
end function q_func2

function q_func_pfn(pfn) result(z) !safety factor, with correct sign
  use precision,only:p_
  use radial_module,only:npsi,pfn_npsi,q_with_sign
  use interpolate_module,only: linear_1d_interpolation
  implicit none
  real(p_):: z,pfn
  call linear_1d_interpolation(npsi,pfn_npsi,q_with_sign,pfn,z)  
end function q_func_pfn


function g_func(psival) result(z) !all quantities are in S.I units.
  use precision,only:p_
  use radial_module,only:npsi,psi_1d,fpsi
  use interpolate_module,only: linear_1d_interpolation
  implicit none
  real(p_):: z,psival
  call linear_1d_interpolation(npsi,psi_1d,fpsi,psival,z)  
end function g_func


!!$function g_func(psival) result(z)
!!$  use precision,only:p_
!!$  use radial_module,only:npsi,psi_1d,fpsi,y2_fpsi
!!$  implicit none
!!$  real(p_):: z,psival
!!$  call splint(psi_1d,fpsi,y2_fpsi,npsi,psival,z) 
!!$end function g_func


function gprime(psival) result(z)
  use precision,only:p_
  use radial_module,only:npsi,psi_1d,fprime
  use interpolate_module,only: linear_1d_interpolation
  implicit none
  real(p_):: z,psival
  call linear_1d_interpolation(npsi,psi_1d,fprime,psival,z)  
end function gprime

!!$function gprime(psival) result(z)
!!$  use precision,only:p_
!!$  use radial_module,only:npsi,psi_1d,fprime,y2_fprime
!!$  implicit none
!!$  real(p_):: z,psival
!!$  call splint(psi_1d,fprime,y2_fprime,npsi,psival,z) 
!!$end function gprime


function g_r(r,z)
  use precision,only:p_
  implicit none
  real(p_):: r,z,g_r
  real(p_):: psi_func,psi_r_func,gprime
  real(p_):: psival
  psival=psi_func(r,z)
  g_r=gprime(psival)*psi_r_func(r,z)
end function g_r

function g_z(r,z)
  use precision,only:p_
  implicit none
  real(p_):: r,z,g_z
  real(p_):: psi_func,psi_z_func,gprime
  real(p_):: psival
  psival=psi_func(r,z)
  g_z=gprime(psival)*psi_z_func(r,z)
end function g_z



!!$function psi_func(xval,zval)
!!$  use precision,only:p_
!!$  use poloidal_flux_2d,only: xarray,zarray,psi,y2a_psi,nx,nz
!!$  implicit none
!!$  real(p_):: psi_func
!!$  real(p_)::xval,zval,psival
!!$  call splin2(xarray,zarray,psi,y2a_psi,nx,nz,xval,zval,psival)
!!$
!!$  psi_func=psival
!!$end function psi_func
!!$


function psi_gradient_func(xval,zval) result (z)
  use precision,only:p_
  use poloidal_flux_2d,only: xarray,zarray,psi_gradient,y2a_gradient,nx,nz
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  real(p_):: z,xval,zval
!  call splin2(xarray,zarray,psi_gradient,y2a_gradient,nx,nz,xval,zval,z)
  call linear_2d_interpolation(nx,nz,xarray,zarray,psi_gradient,xval,zval,z)
end function psi_gradient_func


!!$function psi_r_func(xval,zval) result (z)
!!$  use precision,only:p_
!!$  use poloidal_flux_2d,only: xarray,zarray,psi_x,y2a_psi_x,nx,nz
!!$  implicit none
!!$  real(p_):: z,xval,zval
!!$  call splin2(xarray,zarray,psi_x,y2a_psi_x,nx,nz,xval,zval,z)
!!$end function psi_r_func

!!$function psi_z_func(xval,zval) result (z)
!!$  use precision,only:p_
!!$  use poloidal_flux_2d,only: xarray,zarray,psi_z,y2a_psi_z,nx,nz
!!$  implicit none
!!$  real(p_):: z,xval,zval
!!$  call splin2(xarray,zarray,psi_z,y2a_psi_z,nx,nz,xval,zval,z)
!!$end function psi_z_func
!!$
!!$function psi_rr_func(xval,zval) result (z)
!!$  use precision,only:p_
!!$  use poloidal_flux_2d,only: xarray,zarray,psi_xx,y2a_psi_xx,nx,nz
!!$  implicit none
!!$  real(p_):: z,xval,zval
!!$  call splin2(xarray,zarray,psi_xx,y2a_psi_xx,nx,nz,xval,zval,z)
!!$end function psi_rr_func
!!$
!!$function psi_zz_func(xval,zval) result (z)
!!$  use precision,only:p_
!!$  use poloidal_flux_2d,only: xarray,zarray,psi_zz,y2a_psi_zz,nx,nz
!!$  implicit none
!!$  real(p_):: z,xval,zval
!!$  call splin2(xarray,zarray,psi_zz,y2a_psi_zz,nx,nz,xval,zval,z)
!!$end function psi_zz_func
!!$
!!$
!!$function psi_rz_func(xval,zval) result (z)
!!$  use precision,only:p_
!!$  use poloidal_flux_2d,only: xarray,zarray,psi_xz,y2a_psi_xz,nx,nz
!!$  implicit none
!!$  real(p_):: z,xval,zval
!!$  call splin2(xarray,zarray,psi_xz,y2a_psi_xz,nx,nz,xval,zval,z)
!!$end function psi_rz_func
!!$
!!$function psi_zr_func(xval,zval) result (z)
!!$  use precision,only:p_
!!$  use poloidal_flux_2d,only: xarray,zarray,psi_zx,y2a_psi_zx,nx,nz
!!$  implicit none
!!$  real(p_):: z,xval,zval
!!$  call splin2(xarray,zarray,psi_zx,y2a_psi_zx,nx,nz,xval,zval,z)
!!$end function psi_zr_func


subroutine calculate_poloidal_flux_partial_derivatives(nx,nz,xarray,zarray,psi,psi_x,psi_z,&
     & psi_xx,psi_zz,psi_xz,psi_zx,psi_gradient)
  !subroutine calculate_poloidal_flux_gradient(nx,nz,xarray,zarray,psi,psi_gradient)
  use precision,only:p_
  use constants,only:one,two,twopi
  implicit none
  integer,intent(in)::nx,nz
  real(p_),intent(in):: psi(nx,nz)
  real(p_),intent(in):: xarray(nx),zarray(nz)
  real(p_),intent(out):: psi_x(nx,nz),psi_z(nx,nz),psi_xx(nx,nz),psi_zz(nx,nz),psi_xz(nx,nz),psi_zx(nx,nz)
  real(p_),intent(out):: psi_gradient(nx,nz)

  integer:: i,j,i1,i2,j1,j2

  !first-order partial derivatives
  do i=1,nx
     do j=1,nz
        i2=i+1
        i1=i-1
        j2=j+1
        j1=j-1
        if(i.eq.1) i1=i
        if(j.eq.1) j1=j
        if(i.eq.nx) i2=i
        if(j.eq.nz) j2=j
        psi_x(i,j)=(psi(i2,j)-psi(i1,j))/(xarray(i2)-xarray(i1))
        psi_z(i,j)=(psi(i,j2)-psi(i,j1))/(zarray(j2)-zarray(j1))
     enddo
  enddo

  !second-order partial derivatives
  do i=1,nx
     do j=1,nz
        i2=i+1
        i1=i-1
        j2=j+1
        j1=j-1
        if(i.eq.1) i1=i
        if(j.eq.1) j1=j
        if(i.eq.nx) i2=i
        if(j.eq.nz) j2=j
        psi_xx(i,j)=(psi_x(i2,j)-psi_x(i1,j))/(xarray(i2)-xarray(i1))
        psi_zz(i,j)=(psi_z(i,j2)-psi_z(i,j1))/(zarray(j2)-zarray(j1))
        psi_xz(i,j)=(psi_x(i,j2)-psi_x(i,j1))/(zarray(j2)-zarray(j1))
        psi_zx(i,j)=(psi_z(i2,j)-psi_z(i1,j))/(xarray(i2)-xarray(i1))
     enddo
  enddo


  psi_gradient=sqrt(psi_x**2+psi_z**2)


end subroutine calculate_poloidal_flux_partial_derivatives


subroutine calculate_tfn(npsi,psi_1d,qpsi,tfn_npsi)
  !calculate the toroidal magnetic flux
  use precision,only:p_
  use constants,only: two,twopi
  implicit none
  integer,intent(in):: npsi
  real(p_),intent(in):: psi_1d(npsi),qpsi(npsi)
  real(p_),intent(out):: tfn_npsi(npsi)
  real(p_):: dpsi,tf_npsi(npsi),total_tf
  integer:: j

  dpsi=(psi_1d(npsi)-psi_1d(1))/(npsi-1)
  tf_npsi(1)=0._p_
  do j=2,npsi 
     tf_npsi(j)=tf_npsi(j-1)+(qpsi(j)+qpsi(j-1))/two*twopi*dpsi  !using the formula dtf=q*dpf=q*d(pf_gs)*twopi
  enddo

  open(213,file='pf_tf.txt')
  do j=1,npsi
     write(213,*) psi_1d(j),tf_npsi(j)
  enddo
  close(213)

  total_tf=tf_npsi(npsi)
  tfn_npsi= tf_npsi/total_tf !normalized toroidal magnetic flux
end subroutine calculate_tfn


subroutine draw_rect_region(nx,nz,r_1d,z_1d)
  use precision,only:p_
  implicit none
  integer,intent(in):: nx,nz
  real(p_),intent(in):: r_1d(nx),z_1d(nz)
  integer:: i,j

  open(12,file='rectangular')
  do j=1,nz
     write(12,*) r_1d(1),z_1d(j)
  enddo
  do i=1,nx
     write(12,*) r_1d(i),z_1d(nz)
  enddo
  do j=1,nz
     write(12,*) r_1d(nx),z_1d(nz-j+1)
  enddo

  do i=1,nx
     write(12,*) r_1d(nx-i+1),z_1d(1)
  enddo
  close(12)
end subroutine draw_rect_region


subroutine draw_3d_tokamak()
  use precision,only:p_
  use constants,only: twopi
  use boundary,only: nlim,rlim,zlim,np_lcfs,x_lcfs,z_lcfs !as input
  implicit none
  !  integer,intent(in):: np_lcfs
  ! real(p_),intent(in):: x_lcfs(np_lcfs),z_lcfs(np_lcfs)
  integer:: i,k
  integer,parameter:: nphi=100
  real(p_):: phi
  open(33,file='3d_lcfs.txt')

  do k=1,nphi
     phi=0.+(k-1)*twopi/(nphi-1)
     do i=1,np_lcfs
        write(33,*)  phi, z_lcfs(i), x_lcfs(i),1.0
     enddo
     write(33,*)
  enddo

  open(33,file='3d_limiter.txt')

  do k=1,nphi
     phi=0.+(k-1)*twopi/2./(nphi-1)
     do i=1,nlim
        write(33,*)  phi, zlim(i), rlim(i),1.0
     enddo
     write(33,*)
  enddo
  close(33)
end subroutine draw_3d_tokamak


