subroutine contour(x_lcfs,z_lcfs,np_lcfs,x_axis,z_axis,psival,x_contour,z_contour)
  !given a value of the poloidal flux, psival, this subroutine find the magnetic surface corresponding to this poloidal flux
  use precision,only:p_
  use constants,only:zero,one,two,twopi
  implicit none

  integer,intent(in):: np_lcfs
  real(p_),intent(in):: x_lcfs(np_lcfs),z_lcfs(np_lcfs),x_axis,z_axis,psival
  real(p_),intent(out):: x_contour(np_lcfs),z_contour(np_lcfs)

  real(p_),parameter:: xacc=1.0d-6 !tolerance used in bi-section root-finder
  real(p_):: x1,x2,z1,z2

  real(p_):: slope(np_lcfs),slope2(np_lcfs)
  real(p_):: rtbis !function name of the root finder using the bisection method
  real(p_):: zfunc,xfunc !equation of the straight line (in poloidal plane) that passing throught the magnetic axis point and one point on LCFS
  real(p_):: one_dim_psi_func,one_dim_psi_func2 !one dimension function [psi(x,z(x)) and psi(x(z),z)]on the straight line mentioned in the above.
  external:: one_dim_psi_func,one_dim_psi_func2 !this two function will be passed to a root-finding subroutine
  integer:: i,ierr
  real(p_),parameter:: large_number=1d30


  !do i=1,np_lcfs-1
  do i=1,np_lcfs
     if(x_lcfs(i)-x_axis .eq. 0._p_) then !since I use compiler option which catches all erroneous arithmetic operation, I need to avoid dividing by zero
     slope(i)= large_number
     else
     slope(i)= (z_lcfs(i)-z_axis)/(x_lcfs(i)-x_axis) !the slope for function Z=Z(X)
     endif
     if(z_lcfs(i)-z_axis .eq. 0._p_) then
     slope2(i)=large_number
     else
     slope2(i)=(x_lcfs(i)-x_axis)/(z_lcfs(i)-z_axis) !the slope for function X=X(Z)
     endif
     !write(*,*) i,slope(i),slope2(i)
  enddo

!!$  write(*,*) maxval(slope),minval(slope)
!!$  write(*,*) maxloc(slope),minloc(slope)
!!$  write(*,*) x_axis,x_lcfs( maxloc(slope)),x_lcfs(minloc(slope))
  do i=1,np_lcfs-1  !exclude i=np_lcfs because it is identical to i=1
     if(abs(slope(i)).le.1.0_p_) then !use Z=Z(X) function, the reason that I switch between using function X=X(Z) and Z=Z(X) is to aviod large slope.
        x1=x_axis
        x2=x_lcfs(i) !+0.01 !shift left a little to gurrantee that the range is enough for a root to lie in
        x_contour(i)=rtbis(one_dim_psi_func,x1,x2,xacc,x_axis,z_axis,slope(i),psival)
        z_contour(i)=zfunc(x_axis,z_axis,slope(i),x_contour(i))
     else !switch to using X=X(Z) function
        z1=z_axis
        z2=z_lcfs(i)
        z_contour(i)=rtbis(one_dim_psi_func2,z1,z2,xacc,x_axis,z_axis,slope2(i),psival)
        x_contour(i)=xfunc(x_axis,z_axis,slope2(i),z_contour(i)) 
     endif
  enddo

  x_contour(np_lcfs)=x_contour(1) !i=1 and i=np_lcfs are identical
  z_contour(np_lcfs)=z_contour(1) !i=1 and i=np_lcfs are identical


end subroutine contour


function one_dim_psi_func(x_axis,z_axis,slope,psival,x) 
  !poloidal flux as a function of x on a straight line with slope "slope" in poloidal plane
  use precision,only:p_
  implicit none
  real(p_):: one_dim_psi_func,x,x_axis,z_axis,slope,psival
  real(p_):: zfunc,psi_func
  one_dim_psi_func=psi_func(x,zfunc(x_axis,z_axis,slope,x))-psival
end function one_dim_psi_func


function zfunc(x_axis,z_axis,slope,x) !straight line Z=Z(x) with slope "slope" in poloidal plane starting from the location of magnetic axis
  use precision,only:p_
  implicit none
  real(p_):: zfunc,x,x_axis,z_axis,slope
  zfunc=z_axis+slope*(x-x_axis)
end function zfunc


function one_dim_psi_func2(x_axis,z_axis,slope,psival,z) result(fun_val)
  !poloidal flux as a function of z on a straight line with slope "slope" in poloidal plane
  use precision,only:p_
  implicit none
  real(p_):: fun_val,z
  real(p_):: x_axis,z_axis,slope,psival
  real(p_):: xfunc,psi_func
  fun_val=psi_func(xfunc(x_axis,z_axis,slope,z),z)-psival
end function one_dim_psi_func2


function xfunc(x_axis,z_axis,slope,z) !straight line X=X(Z) with slope "slope" in poloidal plane starting from the location of magnetic axis
  use precision,only:p_
  implicit none
  real(p_):: xfunc,z
  real(p_)::x_axis,z_axis,slope
  xfunc=x_axis+slope*(z-z_axis)
end function xfunc


FUNCTION rtbis(func,x1,x2,xacc,xmaxis,zmaxis,slope,psival)
  !find a root of func by using the bisection method
  use precision,only: p_
  implicit none
  INTEGER JMAX
  REAL(p_) rtbis,x1,x2,xacc,func,xmaxis,zmaxis,slope,psival
  EXTERNAL func
  PARAMETER (JMAX=40)
  INTEGER j
  REAL(p_) dx,f,fmid,xmid
  fmid=func(xmaxis,zmaxis,slope,psival,x2)
  f=   func(xmaxis,zmaxis,slope,psival,x1)
  !      write(*,*) 'f1=', f, 'f2=',fmid
  if(f*fmid.ge.0.) stop 'root must be bracketed in rtbis'

  if(f.lt.0.)then
     rtbis=x1
     dx=x2-x1
  else
     rtbis=x2
     dx=x1-x2
  endif
  do  j=1,JMAX
     dx=dx*.5
     xmid=rtbis+dx
     fmid=func(xmaxis,zmaxis,slope,psival,xmid)
     if(fmid.le.0.)rtbis=xmid
     if(abs(dx).lt.xacc .or. fmid.eq.0.) return
  enddo
  stop 'too many bisections in rtbis'
end FUNCTION rtbis

