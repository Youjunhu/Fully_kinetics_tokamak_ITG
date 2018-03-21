subroutine draw_magnetic_surface(r0,z0,filename) !draw the magnetic surface which passes through (r0,z0)
  use precision,only:p_
  use constants,only:zero,one,two,twopi
  use boundary,only:np_lcfs,x_lcfs,z_lcfs
  use radial_module,only: r_axis,z_axis
  implicit none

  real(p_),intent(in):: r0,z0
  character(*),intent(in)::  filename
  real(p_):: psival
  real(p_):: x_contour(np_lcfs),z_contour(np_lcfs)
  real(p_):: psi_func !function name
  integer:: i
  
  psival=psi_func(r0,z0)

  call contour(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis,psival,x_contour,z_contour)

  open(381,file=filename)
  do i=1,np_lcfs
     write(381,*) x_contour(i),z_contour(i)
  enddo
  close(381)
end subroutine draw_magnetic_surface





