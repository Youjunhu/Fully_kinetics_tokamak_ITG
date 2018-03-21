module derivatives_in_field_line_following_coordinates
  use precision,only:p_
  use constants,only: two
  implicit none
  private
  public radial_derivative,toroidal_derivative, theta_derivative

contains

  subroutine radial_derivative(field,field_x,dradcor) !calculate the derivative of field with respect to psi in field-line-following-coordinates (psi,theta,alpha)
    real(p_),intent(in)::field(:,:)
    real(p_),intent(out)::field_x(:,:)
    real(p_),intent(in):: dradcor
    integer:: m,n,i,j

    m=size(field,1)
    n=size(field,2)
    do i=1,m
       do j=2,n-1
          field_x(i,j)= (field(i,j+1)-field(i,j-1))/(two*dradcor)
       enddo
       field_x(i,1)=(field(i,2)-0._p_)/(two*dradcor) !zero boundary condition for the field is used
       field_x(i,n)=(0._p_-field(i,n-1))/(two*dradcor)!zero boundary condition for the field is used
    enddo

  end subroutine radial_derivative

  subroutine toroidal_derivative(field,field_y,dtor) !calculate the derivative of pper_e with respect to y:
    real(p_),intent(in)::field(:,:)
    real(p_),intent(out)::field_y(:,:)
    real(p_),intent(in):: dtor
    integer:: m,n,i,j,i_left,i_right

    m=size(field,1)
    n=size(field,2)
    do i=1,m   
       do j=1,n
          i_left=i-1
          if(i.eq.1) i_left=m !periodic boundary condition
          i_right=i+1
          if(i.eq.m) i_right=1 !periodic boundary condition
          field_y(i,j)=(field(i_right,j)-field(i_left,j))/(two*dtor)
       enddo
    enddo

  end subroutine toroidal_derivative

  subroutine theta_derivative(a,a_theta,m,n) !partial derivative with respect to theta with psi and alpha fixed, i.e., along the magnetic field line
    use constants,only:two
    use domain_decomposition,only: theta_interval,myid
    implicit none
    integer,intent(in):: m,n
    real(p_),intent(in):: a(m,n)
    real(p_),intent(out):: a_theta(m,n)
    real(p_):: a_left(m,n),a_right(m,n)

    call get_nearby_field_along_field_line(a,a_left,a_right,m,n)
    a_theta=(a_right-a_left)/(two*theta_interval) !centered difference

    !write(*,*) a_left(10,10),a_right(10,10),a_theta(10,10),'myid=',myid
  end subroutine theta_derivative
end module derivatives_in_field_line_following_coordinates

