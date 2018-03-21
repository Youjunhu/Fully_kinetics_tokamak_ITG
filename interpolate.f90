module interpolate_module
  use precision,only:p_
  use constants,only:one
  implicit none
  private
  public linear_2d_interpolation_kernel,linear_2d_interpolation, linear_1d_interpolation, linear_1d_interpolation_tmp,&
       & linear_2d_interpolation_special, linear_1d_interpolation_general_case
contains
  subroutine linear_2d_interpolation(nx,nz,xarray,zarray,psi,x,z,psival)  !uniform xarray and zarray are assumed
    integer,intent(in):: nx,nz
    real(p_),intent(in):: xarray(nx),zarray(nz),psi(nx,nz)
    real(p_),intent(in):: x,z
    real(p_),intent(out)::psival

    real(p_):: dx,dz
    integer:: i,j,ii,jj
    real(p_):: psi_tmp(2,2)

    dx=xarray(2)-xarray(1)
    i=floor(one+(x-xarray(1))/dx) !uniform xarray is assumed, otherwise we need to call location() subroutine to locate xval
    dz=zarray(2)-zarray(1)
    j=floor(one+(z-zarray(1))/dz)

    do ii=1,2   !define a 2x2 array near (x,z)
       do jj=1,2
          psi_tmp(ii,jj)=psi(i+ii-1,j+jj-1)
       enddo
    enddo

    call linear_2d_interpolation_kernel(xarray(i),zarray(j),psi_tmp,x,z,psival)
  end subroutine linear_2d_interpolation


  subroutine linear_2d_interpolation_kernel(x1a,x2a,ya,x1,x2,y)
    real(p_)::x1a(2),x2a(2),ya(2,2),x1,x2,y
    real(p_):: ytmp(2),slope
    integer:: j

    do j=1,2
       slope=(ya(2,j)-ya(1,j))/(x1a(2)-x1a(1))
       ytmp(j)=ya(1,j)+slope*(x1-x1a(1))
    enddo
    slope=(ytmp(2)-ytmp(1))/(x2a(2)-x2a(1))
    y=ytmp(1)+slope*(x2-x2a(1))

  end subroutine linear_2d_interpolation_kernel


  subroutine linear_1d_interpolation(n,x,y,xval,yval)
    use precision,only:p_
    use constants,only:one
    implicit none
    integer,intent(in):: n
    real(p_),intent(in):: x(n),y(n)
    real(p_),intent(in):: xval
    real(p_),intent(out):: yval

    real(p_):: dx,slope
    integer:: i

    dx=x(2)-x(1)
    i=floor(one+(xval-x(1))/dx) !uniform xarray is assumed, otherwise we need to call location() subroutine to locate xval
    !call location(n,x,xval,i)
    if(i.ge.n) i=n-1
    slope=(y(i+1)-y(i))/(x(i+1)-x(i))
    yval=y(i)+slope*(xval-x(i))

  end subroutine linear_1d_interpolation


  subroutine linear_1d_interpolation_general_case(n,x,y,xval,yval)
    use precision,only:p_
    use constants,only:one
    implicit none
    integer,intent(in):: n
    real(p_),intent(in):: x(n),y(n)
    real(p_),intent(in):: xval
    real(p_),intent(out):: yval

    real(p_):: dx,slope
    integer:: i

    !dx=x(2)-x(1)
    !  i=floor(one+(xval-x(1))/dx) !uniform xarray is assumed, otherwise we need to call location() subroutine to locate xval
    call location(n,x,xval,i)
    if(i.ge.n) write(*,*) 'warning*****, x is is exceed the range'
    if(i.ge.n) i=n-1
    slope=(y(i+1)-y(i))/(x(i+1)-x(i))
    yval=y(i)+slope*(xval-x(i))

  end subroutine linear_1d_interpolation_general_case

  subroutine linear_1d_interpolation_tmp(n,x,y,xval,yval)  
    use precision,only:p_
    use constants,only:one
    implicit none
    integer,intent(in):: n
    real(p_),intent(in):: x(n),y(n)
    real(p_),intent(in):: xval
    real(p_),intent(out):: yval

    real(p_):: slope
    integer:: i

    !dx=x(2)-x(1)
    !i=floor(one+(xval-x(1))/dx) !this for uniform x, otherwise we need to call location() subroutine to locate xval
    call location(n,x,xval,i)
    if(i.ge.n) i=n-1

    slope=(y(i+1)-y(i))/(x(i+1)-x(i))
    yval=y(i)+slope*(xval-x(i))

  end subroutine linear_1d_interpolation_tmp


  subroutine location(n,x,xval,k) !use bisection method to locate xval in an array
    !return k (xval is located between x(k) and x(k+1)
    use precision,only:p_
    implicit none
    integer,intent(in):: n
    real(p_),intent(in):: x(n),xval
    integer,intent(out)::k
    integer:: kl,ku,km

    if(xval.gt.x(n) .or.xval.lt.x(1)) write(*,*) "***warning****, x provided is not in the range"
    kl=1
    ku=n
30  if(ku-kl .gt. 1) then  !use bisection method to search location of theta
       km=(ku+kl)/2
       if((x(n).ge.x(1)).eqv.(xval.ge.x(km))) then
          kl=km
       else
          ku=km
       endif
       goto 30
    endif
    k=kl
  end subroutine location


  subroutine linear_2d_interpolation_special(nx,nz,xarray,zarray,psi_a,psi_b,x,z,psival)  !uniform xarray and zarray are assumed
    use precision,only:p_
    use constants,only:one
    use mapping_module,only: i0, j0 !index of the point at magnetic axis
    implicit none
    integer,intent(in):: nx,nz
    real(p_),intent(in):: xarray(nx),zarray(nz),psi_a(nx,nz),psi_b(nx,nz)
    real(p_),intent(in):: x,z
    real(p_),intent(out)::psival
    real(p_):: psi(nx,nz)

    real(p_):: dx,dz
    integer:: i,j,ii,jj
    real(p_):: psi_tmp(2,2)

    dx=xarray(2)-xarray(1)
    i=floor(one+(x-xarray(1))/dx) !uniform xarray is assumed, otherwise we need to call location() subroutine to locate xval
    dz=zarray(2)-zarray(1)
    j=floor(one+(z-zarray(1))/dz)

    psi=psi_a
    !if(i>i0 .and. j+1.eq.j0) psi=psi_b !special treatment for ponts near midplane at the low-field-side 
    if(i<i0 .and. j.eq.j0) psi=psi_b !special treatment for points near the theta cut (now at the high-field-side )

    !define a 2x2 array near (x,z)
    do ii=1,2
       do jj=1,2
          psi_tmp(ii,jj)=psi(i+ii-1,j+jj-1)
       enddo
    enddo

    call linear_2d_interpolation_kernel(xarray(i),zarray(j),psi_tmp,x,z,psival)
  end subroutine linear_2d_interpolation_special
end module interpolate_module

