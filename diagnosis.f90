subroutine draw_grids_on_theta_isosurface(mpoloidal,nflux,tor_shift_mc,r_mag_surf,z_mag_surf) !on theta=constant surface, the subroutine name is wrong, not necessarily top view
  use precision,only: p_
  use constants,only: two,twopi
  !  use magnetic_coordinates,only:mtoroidal
  implicit none
  integer,intent(in):: mpoloidal,nflux
  real(p_),intent(in)::tor_shift_mc(mpoloidal,nflux),r_mag_surf(mpoloidal,nflux),z_mag_surf(mpoloidal,nflux)
  real(p_):: phi,alpha,dalpha
  integer:: i,j,itor,mtoroidal

  i=10
  i=20
  i=60
  i=1
  !  i=mpoloidal
  mtoroidal=20
  dalpha=twopi/mtoroidal

  open(113,file='grids_on_theta_isosurface.txt')
  do itor=1,mtoroidal+1
     !  do itor=1,2
     alpha=dalpha*(itor-1)
     !     do j=1,nflux,10
     !    do j=1,nflux,1
     do j=1,1
        phi=alpha+tor_shift_mc(i,j) !phi is changing due to the radial dependene of tor_shift
        write(113,*) phi,r_mag_surf(i,j),z_mag_surf(i,j)
     enddo
     write(113,*)
     write(113,*)
  enddo
  close(113)

  open(113,file='grids_on_theta_isosurface2.txt')
  do j=1,nflux,10
     do itor=1,mtoroidal+1
        alpha=dalpha*(itor-1)
        phi=alpha+tor_shift_mc(i,j)
        write(113,*) phi,r_mag_surf(i,j),z_mag_surf(i,j)
     enddo
     write(113,*)
     write(113,*)
  enddo
  close(113)
end subroutine draw_grids_on_theta_isosurface


subroutine draw_alpha_isosurface(mpoloidal,nflux,tor_shift_mc,r_mag_surf,z_mag_surf) !on a nonzero radial range
  use precision,only: p_
  use constants,only: two,twopi
  !  use magnetic_coordinates,only:mtoroidal
  implicit none
  integer,intent(in):: mpoloidal,nflux
  real(p_),intent(in)::tor_shift_mc(mpoloidal,nflux),r_mag_surf(mpoloidal,nflux),z_mag_surf(mpoloidal,nflux)
  real(p_):: phi,alpha0
  integer:: i,j

  alpha0=twopi/8._p_
  open(113,file='alpha_isosurface.txt')
  do j=1,101,5
     do i=1,mpoloidal
        phi=alpha0+tor_shift_mc(i,j) 
        !phi=alpha0
        write(113,*) phi,r_mag_surf(i,j),z_mag_surf(i,j)
     enddo
     write(113,*)
     write(113,*)
  enddo
  close(113)

end subroutine draw_alpha_isosurface


subroutine draw_alpha_contours_on_a_magnetic_surface(mpoloidal,nflux,tor_shift_mc,r_mag_surf,z_mag_surf) !field lines on a magnetic surface
  use precision,only: p_
  use constants,only: two,twopi
  !  use magnetic_coordinates,only:mtoroidal
  implicit none
  integer,intent(in):: mpoloidal,nflux
  real(p_),intent(in)::tor_shift_mc(mpoloidal,nflux),r_mag_surf(mpoloidal,nflux),z_mag_surf(mpoloidal,nflux)
  real(p_):: phi,alpha0,tor_range
  integer:: i,j,ialpha,ishift
  integer::nalpha=10

  open(113,file='alpha_contours_on_magnetic_surface.txt')
  !  do j=1,30,2
  j=15
  tor_range=twopi/4
  do ialpha=1,nalpha
     alpha0=0._p_+tor_range/(nalpha-1)*(ialpha-1)
     do i=1,mpoloidal
        phi=alpha0+tor_shift_mc(i,j)
        if(phi<0. .or. phi>tor_range) then !shift into the range [0:tor_range]
           ishift=floor(phi/tor_range)
           phi=phi-ishift*tor_range 
           alpha0=alpha0-ishift*tor_range
           write(113,*)
           write(113,*)
        endif
        !phi=alpha0
        write(113,*) phi,r_mag_surf(i,j),z_mag_surf(i,j)
     enddo
     write(113,*)
     write(113,*)
  enddo
  close(113)

end subroutine draw_alpha_contours_on_a_magnetic_surface
