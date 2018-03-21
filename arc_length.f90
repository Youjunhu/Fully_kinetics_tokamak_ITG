subroutine arc_length(x_contour,z_contour,np_lcfs,dl)
  !calculate the poloidal arc length between neighbour points on every contour line.
  use precision,only:p_
  implicit none
  integer,intent(in):: np_lcfs
  real(p_),intent(in):: x_contour(np_lcfs),z_contour(np_lcfs) 
  real(p_),intent(out):: dl(np_lcfs-1)
  integer:: i
  integer:: i_plus_one,i_plus_two,i_minus_one
  real(p_):: x(4),z(4),tmp
 
     do i=1,np_lcfs-1
        i_plus_one=i+1  !i_plus_one indicates the right point
        i_minus_one=i-1 !i_minus_one indicates the left point
        i_plus_two=i+2
        if (i .eq. np_lcfs) i_plus_one=2 !deal with boundary points
        if (i .eq. 1)       i_minus_one=np_lcfs-1 !deal with boundary points
        if (i .eq. np_lcfs) i_plus_two=3 !deal with boundary points
        if (i .eq. np_lcfs-1) i_plus_two=2 !deal with boundary points
        x(1)=x_contour(i_minus_one)
        x(2)=x_contour(i)
        x(3)=x_contour(i_plus_one)
        x(4)=x_contour(i_plus_two)
        z(1)=z_contour(i_minus_one)
        z(2)=z_contour(i)
        z(3)=z_contour(i_plus_one)
        z(4)=z_contour(i_plus_two)
        
        call arc_between_two_points(x,z,tmp)
        dl(i)=tmp

        !dl(i)=sqrt((x_contour(i)-x_contour(i+1))**2+(z_contour(i)-z_contour(i+1))**2) !use simple formula to calculate the arc length
     enddo

end subroutine arc_length



subroutine   arc_between_two_points(x,z,dl)
!calculate the arc length between point (x(2),z(2)) and point (x(3),z(3))
  use precision,only:p_
  use constants,only:one,two
  implicit none
  real(p_),intent(in):: x(4),z(4)
  real(p_),intent(out):: dl
  real(p_):: ds,a_length,b_length
  real(p_):: dot_a_and_ds,dot_b_and_ds, cos_tha,cos_thb,m1,m2

  !ds is the length of straight-line segment passing through (x(2),z(2)) and (x(3),z(3))
  ds=sqrt((x(3)-x(2))**2+(z(3)-z(2))**2) 

  a_length=sqrt((x(3)-x(1))**2+(z(3)-z(1))**2)
  b_length=sqrt((x(4)-x(2))**2+(z(4)-z(2))**2)
  dot_a_and_ds=(x(3)-x(1))*(x(3)-x(2)) &
       +(z(3)-z(1))*(z(3)-z(2))
  dot_b_and_ds=(x(4)-x(2))*(x(3)-x(2)) &
       +(z(4)-z(2))*(z(3)-z(2))
  cos_tha=dot_a_and_ds/(a_length*ds)
  cos_thb=dot_b_and_ds/(b_length*ds)

  m1=sqrt(one-cos_tha**2)/cos_tha 
  m2=sqrt(one-cos_thb**2)/cos_thb
  !the value of m1 and m2 should be positive for most cases
  dl=ds*(one+(two*m1**2+two*m2**2+m1*m2)/30._p_) !calculate arc length using Eq. (5.38) in  S. Jardin's book. Here I assume that the dot product of the slope is negative, so the -m1*m2 term is replaced with +abs(slope1)*abs(slope2), need checking the correctness for general case
  !dl=ds*(1._p_+0.) !use linear function to approximate the arc length
end subroutine arc_between_two_points
