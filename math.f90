subroutine shift_to_zero_twopi_range(a) !shift a into the range [0:twopi]
  use precision,only:p_
  use constants,only: twopi
  implicit none
  real(p_),intent(inout):: a
  integer:: ishift
!!$  a=a-int(a/twopi)*twopi !shift into the range [0:twopi]
!!$  if(a.lt.0) a=a+twopi !shift into the range [0:twopi]

   ishift=floor(a/twopi)
   a=a-ishift*twopi

end subroutine shift_to_zero_twopi_range


subroutine shift_to_specified_toroidal_range(a) !shift "a" into the range [0:toroidal_range]
  use precision,only:p_
!  use constants,only: twopi
  use magnetic_coordinates,only: toroidal_range
  implicit none
  real(p_),intent(inout):: a
  integer:: ishift
!!$  a=a-int(a/twopi)*twopi !shift into the range [0:twopi]
!!$  if(a.lt.0) a=a+twopi !shift into the range [0:twopi]
!!$  a=a-int(a/toroidal_range)*toroidal_range !shift into the range [0:toroidal_range]
!!$  if(a.lt.0) a=a+toroidal_range !shift into the range [0:toroidal_range]

 ishift=floor(a/toroidal_range)
 a=a-ishift*toroidal_range

end subroutine shift_to_specified_toroidal_range


subroutine shift_to_minus_pi_positive_pi_range(a) !shift "a" into the range [-pi:pi]
  use precision,only:p_
  use constants,only: twopi,pi
  implicit none
  real(p_),intent(inout):: a
  integer:: ishift

 ishift=floor(a/twopi)
 a=a-ishift*twopi
if(a>pi) a=a-twopi
end subroutine shift_to_minus_pi_positive_pi_range

subroutine one_dimensional_derivative(n,x,y,dydx)
  use precision,only:p_
  use constants,only:zero,one,two,twopi,one_half
  implicit none

  integer,intent(in):: n
  real(p_),intent(in):: x(n),y(n)
  real(p_),intent(out):: dydx(n)
  real(p_):: tmp0,dx
  integer:: j

  dx=x(2)-x(1) !uniform interval is assumed

  do j=2,n-1 !use center difference scheme for inner points
     dydx(j)=(y(j+1)-y(j-1))/(two*dx)
  enddo

  !use linear interpolation to get the value  j=n
  tmp0=(y(n)-y(n-1))/dx
  dydx(n)=two*tmp0-dydx(n-1)
  !use linear interpolation to get the value j=1
  tmp0=(y(2)-y(1))/dx
  dydx(1)=two*tmp0-dydx(2)
end subroutine one_dimensional_derivative

