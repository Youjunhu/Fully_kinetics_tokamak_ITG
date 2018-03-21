subroutine set_toroidal_range
  use domain_decomposition,only:myid
  use constants,only: twopi
  use magnetic_coordinates,only: nsegment
  use magnetic_coordinates,only: toroidal_range !output
  implicit none

  toroidal_range=twopi/nsegment
if(myid.eq.0) write(*,*) 'toroidal range=',toroidal_range
end subroutine set_toroidal_range


subroutine create_toroidal_grids()
  use precision,only:p_
  use constants,only: twopi
  use magnetic_coordinates,only: mtoroidal,toroidal_range
  use magnetic_coordinates,only: tor_1d_array,dtor !as output
  implicit none
  integer:: i

  allocate(tor_1d_array(mtoroidal+1)) 
  !dtor=twopi/mtoroidal
  dtor=toroidal_range/mtoroidal
  do i=1,mtoroidal+1 ! No. 1 point is phi=0, No. mtoroidal+1 is phi=toroidal_range
     tor_1d_array(i)=0._p_+dtor*(i-1)
  enddo
end subroutine create_toroidal_grids
