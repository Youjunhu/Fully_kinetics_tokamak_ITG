module mode_structure
  implicit none
  private
  public mode_structure_on_theta_plane,mode_structure_on_poloidal_plane,mode_structure_on_xz_plane,&
       & mode_structure_on_yz_plane

contains
  subroutine mode_structure_on_theta_plane(kt,GCLR,a,partial_file_name)
    use precision,only:p_
    use magnetic_coordinates,only:radcor_1d_array2,tor_1d_array
    use domain_decomposition,only:myid
    integer,intent(in)::kt,GCLR
    real(p_),intent(in):: a(:,:)
    character(len=*),intent(in)::partial_file_name
    character(100):: full_file_name
    integer:: i,j,m,n,u

    full_file_name='ms/tor_rad_polxxx_txxxxxx'//partial_file_name
    write(full_file_name(15:17),'(i3.3)') GCLR
    write(full_file_name(20:25),'(i6.6)') kt
    open(newunit=u,file=full_file_name)

    m=size(a,1)
    n=size(a,2)
    do j=1,n
       do i=1,m
          write(u,*) radcor_1d_array2(j),tor_1d_array(i),a(i,j)
       enddo
       write(u,*)
    enddo
    close(u)
  end subroutine mode_structure_on_theta_plane

  subroutine mode_structure_on_poloidal_plane(kt,a,partial_file_name)
    use precision,only:p_
    use constants,only:pi
    use magnetic_coordinates,only: tor_1d_array,r_mag_surf, z_mag_surf,tor_shift_mc,mpoloidal,nsegment,&
         & dtheta,theta_1d_array
    use domain_decomposition,only:numprocs,myid,ntube,tube_comm,GCLR,theta_start,theta_interval,GCLR_cut
    use constants,only:twopi
  use interpolate_module
    use mpi
    implicit none
    integer,intent(in)::kt
    real(p_),intent(in):: a(:,:)
    character(len=*),intent(in)::partial_file_name
    integer:: itor,ipol,ipol_eq,j,ierr,m,n
    real(p_):: my_a_poloidal_plane(size(a,2)),a_poloidal_plane(size(a,2),0:numprocs/ntube-1)
    real(p_):: phi_array(size(a,1)),a_phi(size(a,1))
    real(p_)::phi_val,my_theta
    character(100)::full_file_name
    integer:: u !file unit number

    m=size(a,1)
    n=size(a,2)

  !  allocate(my_a_poloidal_plane(n))
  !  allocate(a_poloidal_plane(n,0:numprocs/ntube-1))
  !  allocate(phi_array(m))
  !  allocate(a_phi(m))

    phi_val=0._p_+0.5_p_*twopi/nsegment !choose a fixed cylindrical toroidal angle, this is the poloidal plane on which the mode structure is computed
    !ipol_eq=1+GCLR*(mpoloidal-1)/(numprocs/ntube) !the poloidal sequence (in the equilibrium grids) of proc GCLR, wrong!
    ipol_eq=1+nint((theta_start-theta_1d_array(1))/dtheta) !index of begining theta angle of this subdomain in the original magnetic coordinate grids
    do j=1,n
       do itor=1,m
          phi_array(itor)=tor_1d_array(itor)+tor_shift_mc(ipol_eq,j)
          call shift_to_specified_toroidal_range(phi_array(itor))
          a_phi(itor)=a(itor,j)
       enddo
       call bubble_sort(phi_array,a_phi,m) !make phi_array change monotoically, which is required (by the interpolating subroutine) for phi_array to be used in the interpolating subroutine, the dependent variable a_phi should be also be changed accordingly if the sequence of phi_array array is changed
       call linear_1d_interpolation_general_case(m,phi_array,a_phi,phi_val,my_a_poloidal_plane(j)) !interpolated to the same cylindrical toroidal angle
    enddo
    call MPI_gather(my_a_poloidal_plane, n, MPI_real8,&
                       a_poloidal_plane, n, MPI_real8, 0, tube_COMM, ierr)
    if(myid.eq.0) then
       full_file_name='ms/poloidal_plane_txxxxxx'//partial_file_name
       write(full_file_name(20:25),'(i6.6)') kt
       open(newunit=u,file=full_file_name)
       do j=1,n
          do ipol=0,numprocs/ntube-1
             !ipol_eq=1+ipol*(mpoloidal-1)/(numprocs/ntube) !the poloidal sequence (in the equilibrium grids), wrong!
             my_theta=ipol*theta_interval
             if(ipol.gt.GCLR_cut) my_theta=-pi+(ipol-GCLR_cut-1)*theta_interval
             ipol_eq=1+nint((my_theta-theta_1d_array(1))/dtheta) !index of begining theta angle of this subdomain in the original magnetic coordinate grids
             write(u,*) r_mag_surf(ipol_eq,j), z_mag_surf(ipol_eq,j),a_poloidal_plane(j,ipol)
          enddo
          write(u,*) 
       enddo
       close(u)
    endif
  end subroutine mode_structure_on_poloidal_plane


SUBROUTINE Bubble_Sort(a,b,n) !I learned this algorithm when I was a undergraduate student. The following code is adapted from the code at http://rosettacode.org/wiki/Sorting_algorithms/Bubble_sort#Fortran
 use precision,only:p_
  implicit none
  integer,intent(in)::n
!  REAL(p_), INTENT(in out), DIMENSION(:) :: a,b
  REAL(p_), INTENT(in out), DIMENSION(n) :: a,b
  REAL(p_) :: temp
  INTEGER :: i, j
  LOGICAL :: swapped
 
!  DO j = SIZE(a)-1, 1, -1
  DO j = n-1, 1, -1
    swapped = .FALSE.
    DO i = 1, j
      IF (a(i) > a(i+1)) THEN
        temp = a(i) !swap array a
        a(i) = a(i+1)
        a(i+1) = temp
        swapped = .TRUE.
        temp = b(i) !also swap array b
        b(i) = b(i+1)
        b(i+1) = temp
      END IF
    END DO
    IF (.NOT. swapped) EXIT
  END DO
END SUBROUTINE Bubble_Sort



subroutine mode_structure_on_xz_plane(kt,a,partial_file_name)
  use precision,only:p_
  use constants,only:pi
  use magnetic_coordinates,only: radcor_1d_array2,theta_1d_array,dtheta,mtoroidal
  use domain_decomposition,only:numprocs,myid,ntube,tube_comm,GCLR,GCLR_cut,theta_interval
  use constants,only:twopi
  use mpi
  implicit none
  integer,intent(in)::kt
  real(p_),intent(in):: a(:,:)
  character(len=*),intent(in)::partial_file_name
  integer:: itor,ipol,ipol_eq,j,ierr,m,n
  real(p_):: my_a_xz_plane(size(a,2)),a_xz_plane(size(a,2),0:numprocs/ntube-1)
  real(p_)::my_theta
  character(100)::full_file_name
  integer:: u !file unit number

   m=size(a,1)
   n=size(a,2)
!  allocate(my_a_xz_plane(n))
!  allocate(a_xz_plane(n,0:numprocs/ntube-1))

  itor=mtoroidal/2 !choose a alpha (i.e., y) grid
  do j=1,n
     my_a_xz_plane(j)=a(itor,j)
  enddo
  call MPI_gather(my_a_xz_plane, n, MPI_real8,&
                     a_xz_plane, n, MPI_real8, 0, tube_COMM, ierr)
  if(myid.eq.0) then
     full_file_name='ms/xz_txxxxxx'//partial_file_name
     write(full_file_name(8:13),'(i6.6)') kt
     open(newunit=u,file=full_file_name)
     do ipol=0,numprocs/ntube-1 !poloidal direction
        my_theta=ipol*theta_interval
        if(ipol.gt.GCLR_cut) my_theta=-pi+(ipol-GCLR_cut-1)*theta_interval
        ipol_eq=1+nint((my_theta-theta_1d_array(1))/dtheta) !index of begining theta angle of this subdomain in the original magnetic coordinate grids
        do j=1,n !radial direction
           write(u,*) radcor_1d_array2(j), theta_1d_array(ipol_eq),a_xz_plane(j,ipol)
        enddo
        write(u,*) 
     enddo
     close(u)
  endif
end subroutine mode_structure_on_xz_plane


subroutine mode_structure_on_yz_plane(kt,a,partial_file_name)
  use precision,only:p_
  use constants,only:pi
  use magnetic_coordinates,only: mtoroidal,tor_1d_array,theta_1d_array,dtheta,nflux2
  use domain_decomposition,only:numprocs,myid,ntube,tube_comm,GCLR,GCLR_cut,theta_interval
  use constants,only:twopi
  use mpi
  implicit none
  integer,intent(in)::kt
  real(p_),intent(in):: a(:,:)
  character(len=*),intent(in)::partial_file_name
  integer:: itor,ipol,ipol_eq,jrad,ierr,m,n
  real(p_):: my_a_yz_plane(size(a,1)),a_yz_plane(size(a,1),0:numprocs/ntube-1)
  real(p_)::my_theta
  character(100)::full_file_name
  integer:: u !file unit number

  m=size(a,1)
   n=size(a,2)
 
!  allocate(my_a_yz_plane(m))
!  allocate(a_yz_plane(m,0:numprocs/ntube-1))

  jrad=nflux2/2 !choose a radial (i.e., x) grid
  do itor=1,mtoroidal
     my_a_yz_plane(itor)=a(itor,jrad)
  enddo
  call MPI_gather(my_a_yz_plane, m, MPI_real8,&
                     a_yz_plane, m, MPI_real8, 0, tube_COMM, ierr)
  if(myid.eq.0) then
     full_file_name='ms/yz_txxxxxx'//partial_file_name
     write(full_file_name(8:13),'(i6.6)') kt
     open(newunit=u,file=full_file_name)
     do ipol=0,numprocs/ntube-1 !poloidal direction
        my_theta=ipol*theta_interval
        if(ipol.gt.GCLR_cut) my_theta=-pi+(ipol-GCLR_cut-1)*theta_interval
        ipol_eq=1+nint((my_theta-theta_1d_array(1))/dtheta) !index of begining theta angle of this subdomain in the original magnetic coordinate grids
        do itor=1,mtoroidal !toroidal direction
           write(u,*) tor_1d_array(itor), theta_1d_array(ipol_eq), a_yz_plane(itor,ipol)
        enddo
        write(u,*) 
     enddo
     close(u)
  endif
end subroutine mode_structure_on_yz_plane
end module mode_structure






