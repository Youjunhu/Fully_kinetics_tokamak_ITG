subroutine set_computational_radial_region()
  use precision,only:p_
  use magnetic_coordinates,only: dradcor,radcor_1d_array !as input
  use magnetic_coordinates,only: radcor_low0,radcor_upp0, radcor_low1,radcor_upp1,radcor_low2,radcor_upp2,j_low2,j_upp2,nflux2 !as output
  use magnetic_coordinates,only: j_low1,j_upp1,nflux1, radcor_1d_array2 !as output
  use domain_decomposition,only:myid
  implicit none

  integer::j,j_shift

  radcor_low0=minval(radcor_1d_array)
  radcor_upp0=maxval(radcor_1d_array)

  radcor_low1=radcor_low0+(radcor_upp0-radcor_low0)*0.05!decrease radcor_low1 to allow for a buffer region for particle motion
  radcor_upp1=radcor_upp0-(radcor_upp0-radcor_low0)*0.05 !increase radcor_upp1 to allow for a buffer region for particle motion
  j_low1=1+ceiling((radcor_low1-radcor_1d_array(1))/dradcor) !to locate a nearby grid
  j_upp1=1+  floor((radcor_upp1-radcor_1d_array(1))/dradcor)  !to locate a nearby grid
  nflux1=j_upp1-j_low1+1
  radcor_low1=radcor_1d_array(j_low1)
  radcor_upp1=radcor_1d_array(j_upp1)

  radcor_low2=radcor_low1+(radcor_upp1-radcor_low1)*0.05!further shrink the radial range, to be used as the computational boundary, to allow for a buffer region for the interpolation when marker go out of this region
  radcor_upp2=radcor_upp1-(radcor_upp1-radcor_low1)*0.05 !furthre shrink the radial range

  j_low2=1+ceiling((radcor_low2-radcor_1d_array(1))/dradcor) !to locate a nearby grid
  j_upp2=1+floor  ((radcor_upp2-radcor_1d_array(1))/dradcor)  !to locate a nearby grid
  !if(mod(j_upp2-j_low2+1,2).ne.0) j_low2=j_low2+1 !make j_upp2-j_low2+1  an even number
  nflux2=j_upp2-j_low2+1 !the total number of radial grids

  radcor_low2=radcor_1d_array(j_low2)
  radcor_upp2=radcor_1d_array(j_upp2)

  if(myid.eq.0) write(*,*) 'j_low1=',j_low1,'j_upp1=',j_upp1,'j_low2=',j_low2,'j_upp2=',j_upp2
  if(myid.eq.0) write(*,*) 'radcor_low1=',radcor_low1,'radcor_upp1=',radcor_upp1
  if(myid.eq.0) write(*,*) 'radcor_low2=',radcor_low2,'radcor_upp2=',radcor_upp2

  allocate(radcor_1d_array2(nflux2))
  do j=1,nflux2
     j_shift=j_low2+j-1
     radcor_1d_array2(j)=radcor_1d_array(j_shift)
  enddo
end subroutine set_computational_radial_region

subroutine plasma_volume_of_computational_region()
  !to calculate the spatial volume of computational region
  use precision,only: p_
  use constants,only: two,twopi
  use domain_decomposition,only: myid
  use magnetic_coordinates,only: dv,j_low2,j_upp2,mpoloidal,toroidal_range
  use magnetic_coordinates,only: vol2 !as output
  implicit none
  integer:: i,j

  vol2=0._p_
  do i=1,mpoloidal-1
     do j=j_low2,j_upp2
        vol2=vol2+dv(i,j)
     enddo
  enddo

  if(myid.eq.0) write(*,*) 'volume of full torus of the computational domain=', vol2

  vol2=vol2/twopi*toroidal_range
  !if(myid.eq.0) write(*,*) 'volume of the toroidal wedge=', vol2
end subroutine plasma_volume_of_computational_region

subroutine plasma_volume_to_buffer()
  !to calculate the spatial volume of domain including buffer region
  use precision,only: p_
  use constants,only: two,twopi
  use domain_decomposition,only: myid
  use magnetic_coordinates,only: dv,j_low1,j_upp1,mpoloidal,toroidal_range
  use magnetic_coordinates,only: vol1 !as output
  implicit none
  integer:: i,j

  vol1=0._p_
  do i=1,mpoloidal-1
     do j=j_low1,j_upp1
        vol1=vol1+dv(i,j)
     enddo
  enddo

  if(myid.eq.0) write(*,*) 'volume of full torus of domain including buffer=', vol1

  vol1=vol1/twopi*toroidal_range
  !if(myid.eq.0) write(*,*) 'volume of the toroidal wedge including buffer=', vol1
end subroutine plasma_volume_to_buffer
