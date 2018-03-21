subroutine field_lines_analyse()
  use precision,only:p_
  use constants,only: two
  use radial_module,only: z_axis
  implicit none

  integer:: n_tor_loop,max_npt_along_field_line,krad !used in field line tracing module
  namelist /field_line_tracing_nl/  n_tor_loop,max_npt_along_field_line,krad
  real(p_),allocatable:: r_start(:),z_start(:),phi_start(:) !starting point of field lines
  real(p_),allocatable:: r_poincare(:,:),z_poincare(:,:),phi_poincare(:,:) !pointcare points
  integer:: j,k
  integer,allocatable::nloop_actual(:)

  open(31,file='input.nmlt')
  read(31,field_line_tracing_nl)
  close(31)
  write(*,field_line_tracing_nl)
  allocate(r_start(krad))
  allocate(z_start(krad))
  allocate(phi_start(krad))
  allocate(r_poincare(n_tor_loop+1,krad))
  allocate(z_poincare(n_tor_loop+1,krad))
  allocate(phi_poincare(n_tor_loop+1,krad))
  allocate(nloop_actual(krad))

!  !$omp parallel do
  do k=1,krad
     r_start(k)=2.15_p_+(k-1)*0.15d0/(20-1)
     z_start(k)=z_axis
     phi_start(k)=0._p_
     call field_line_tracing(r_start(k),z_start(k),phi_start(k), max_npt_along_field_line,n_tor_loop,&
          & r_Poincare(:,k),z_Poincare(:,k),phi_Poincare(:,k),nloop_actual(k))
  enddo
!  !$omp end parallel do

  open(26,file='poincare.txt')
  do k=1,krad
     do j=1,nloop_actual(k)
        write(26,*) r_Poincare(j,k),z_Poincare(j,k),phi_Poincare(j,k)
     enddo
     write(26,*)
     write(26,*)
  enddo
  close(26)



!!$  do k=1,krad
!!$     call draw_magnetic_surface(r_start(k),z_start(k),'ref_field_line.txt') !draw the magnetic surface which passes through (r0,z0)
!!$  enddo


end subroutine field_lines_analyse



subroutine field_line_tracing(r0,z0,phi0,npt,n_tor_loop,r_poincare,z_poincare,phi_poincare,nloop_actual)
  !given coordinates (R,Z,phi), this subroutine follows the field lines passing through this point until it has finish n_tor_loop toroidal loop or exceeds the specifed maximum number of points along the field-line, npt. This subroutine also calculates the safety factor of the field-line found.
  use precision,only:p_
  use constants,only: two,twopi,one_half
  use boundary,only: np_lcfs,x_lcfs,z_lcfs,nlim,rlim,zlim !use to check whether field line touch the boundary
  implicit none

  real(p_),intent(in):: r0,z0,phi0
  integer,intent(in)::npt,n_tor_loop
  real(p_):: r(npt),z(npt),phi(npt)
  real(p_),intent(out):: r_poincare(n_tor_loop+1),z_poincare(n_tor_loop+1),phi_poincare(n_tor_loop+1)
  integer,intent(out):: nloop_actual
  real(p_),parameter:: step=1d-3 !meter, trial of dr or dz step
  real(p_):: brval,bzval,bphival,bpolval,dr,dz,dphi
  real(p_):: br_SI,bz_SI,bphi_SI
  real(p_):: r_mid,z_mid,dl_pol,qval
  logical:: loss
  integer:: j,k,jj


  k=1 !Poincare points
  r_poincare(k)=r0 !Poincare points
  z_poincare(k)=z0
  phi_poincare(k)=phi0


  r(1)=r0
  z(1)=z0
  phi(1)=phi0

  loss=.false.
  do j=1,npt-1
     !2nd Runge-Kutta
     brval=    br_SI(r(j),z(j),phi(j))
     bzval=    bz_SI(r(j),z(j),phi(j))
     bphival=bphi_SI(r(j),z(j),phi(j))
     bpolval=sqrt(brval**2+bzval**2)

     if(abs(bzval).lt.abs(brval)) then
        dr=step*one_half
        if(brval.lt.0._p_) dr=-step*one_half
        dz=bzval/brval*dr

     else
        dz=step*one_half
        if(bzval.lt.0._p_) dz=-step*one_half
        dr=brval/bzval*dz
     endif

     dl_pol=sqrt(dr**2+dz**2)
     dphi=bphival/bpolval*dl_pol/r(j)

     !first step:
     r_mid=r(j)+dr
     z_mid=z(j)+dz
     brval=    br_SI(r_mid,z_mid,phi(j)+dphi)
     bzval=    bz_SI(r_mid,z_mid,phi(j)+dphi)
     bphival=bphi_SI(r_mid,z_mid,phi(j)+dphi)
     bpolval=sqrt(brval**2+bzval**2)


     if(abs(bzval).lt.abs(brval)) then
        dr=step
        if(brval.lt.0._p_) dr=-step
        dz=bzval/brval*dr
     else
        dz=step
        if(bzval.lt.0._p_) dz=-step
        dr=brval/bzval*dz
     endif

     dl_pol=sqrt(dr**2+dz**2)
     dphi=bphival/(bpolval*r_mid)*dl_pol

     r(j+1)=r(j)+dr
     z(j+1)=z(j)+dz
     phi(j+1)=phi(j)+dphi

     call  check_whether_field_line_touch_boundary(r(j+1),z(j+1),phi(j+1),x_lcfs,z_lcfs,np_lcfs,loss)
     if (loss.eqv. .true.) exit

     if(abs(floor(abs(phi(j)-phi0)/twopi)-floor(abs(phi(j+1)-phi0)/twopi)).eq.1) then ! finish one toroidal turn
        !   write(*,*) 'j=',j,'k=',k, phi(j),phi(j+1)
        k=k+1
        r_poincare(k)=(r(j)+r(j+1))/two
        z_poincare(k)=(z(j)+z(j+1))/two
        phi_poincare(k)=((phi(j)+phi(j+1))/two)/twopi
     endif

     if(abs(phi0-phi(j+1))/twopi.ge.n_tor_loop) exit

  enddo

 if(j.eq.npt) then
     open(76,file='bad_line.txt')
     do jj=1,j-1
     write(76,*) r(jj),z(jj),phi(jj)
     enddo
     close(76)
     call safety_factor_a_field_line(r,z,phi,j,qval)
     call draw_magnetic_surface(r0,z0,'ref_field_line.txt') !draw the magnetic surface which passes through (r0,z0)
     stop 'max number of tracing steps of field line is exceeded before achiving the specified number of toroidal loop'
endif




  nloop_actual=k
  write(*,*) 'nloop_actual=',nloop_actual, 'actual step along field line=',j



  call safety_factor_a_field_line(r,z,phi,j+1,qval)

!!$  open(76,file='field_line.txt')
!!$  do jj=1,j+1
!!$     write(76,*) r(jj),z(jj),phi(jj)
!!$  enddo
!!$  close(76)


call check_field_line_in_field_aligned_coordinates(r,z,phi,j,qval)

end subroutine field_line_tracing

subroutine check_field_line_in_field_aligned_coordinates(r,z,phi,npt,qval)
  use precision,only:p_
  use constants,only: two,twopi
  use radial_module,only:z_axis
  use mapping_module,only: nx_mapping ,j0,r_cyl,tor_shift_b
  implicit none
  integer,intent(in):: npt
  real(p_),intent(in):: r(npt),z(npt),phi(npt),qval
  real(p_):: radcor(npt),theta(npt),alpha(npt),tor_shift(npt)
  real(p_)::x1,x2,y1,y2,z1,z2,dx,dy,dz,dl(npt)
  real(p_):: q_func,psi_func,q_func2,pfn_func,radcor_as_func_of_pfn
  real(p_):: sum=0._p_,real_shift,total_shift,tmp_array(nx_mapping)
  integer:: j,kk



  do j=1,npt
     radcor(j)=radcor_as_func_of_pfn(pfn_func(r(j),z(j))) !get radial coordinate
     call interpolate_from_cylindrical_to_magnetic_coordinates(r(j),z(j),theta(j),tor_shift(j))
  enddo

  dl(1)=0._p_
  do j=2,npt
     x1=r(j-1)*cos(phi(j-1))
     x2=r(j)*cos(phi(j))
     y1=r(j-1)*sin(phi(j-1))
     y2=r(j)*sin(phi(j))
     z1=z(j-1)
     z2=z(j)
     dx=x2-x1
     dy=y2-y1
     dz=z2-z1
     dl(j)=dl(j-1)+sqrt(dx*dx+dy*dy+dz*dz)
  enddo

  open(11,file='field_line_in_field_aligned_co.txt')
  do j=2,npt

!!$     if(abs(theta(j)-theta(j-1)) .ge. twopi*0.9) then !indecate finishing one poloidal loop
!!$        !total_shift=tor_shift(j)*(z_axis-z(j-1))/(z(j)-z_axis)+tor_shift(j-1)
!!$        !total_shift=q_func2(radcor(j))*twopi
!!$        !total_shift=qval*twopi*1.0004
!!$        !total_shift=2.35227*twopi
!!$        do kk=1,nx_mapping
!!$           tmp_array(kk)= tor_shift_b(kk,j0)
!!$        enddo
!!$        call linear_1d_interpolation(nx_mapping,r_cyl,tmp_array,r(j),total_shift)
!!$        sum=sum+total_shift
!!$     endif
!!$     alpha(j)=phi(j)-(tor_shift(j)+sum)
   call accumulate_tor_shift(theta(j-1),theta(j),r(j),tor_shift(j),real_shift)
    alpha(j)=phi(j)-real_shift 
     write(11,*) dl(j),radcor(j),theta(j),alpha(j), tor_shift(j), phi(j),r(j),z(j)
  enddo
  close(11)
end subroutine check_field_line_in_field_aligned_coordinates


subroutine safety_factor_a_field_line(r,z,phi,npt,qval)
  !given a field line, calculate its safety factor
  use precision,only:p_
  use constants,only: two,twopi
  implicit none
  integer,intent(in):: npt
  real(p_),intent(in):: r(npt),z(npt),phi(npt)
  real(p_),intent(out):: qval
  real(p_):: phi_old
  real(p_):: q_func,psi_func
  integer:: j,npass_midplane

  npass_midplane=0
  do j=1,npt-1
     if(z(j)*z(j+1).lt.0) then !indicates one midplane-crossing
        npass_midplane=npass_midplane+1
        if(npass_midplane.eq.1) phi_old=(phi(j)+phi(j+1))/two
     endif
     if(npass_midplane.eq.3) then !indicates that the line has finished one poloidal period
        qval=abs(phi_old-phi(j))/twopi
        write(*,*) 'safety factor of field line passing (r,z)',r(1),z(1),'is', qval,&
             & 'q value specified in gfile =', q_func(psi_func(r(1),z(1)))
        exit
     endif
  enddo

  write(*,*) 'toroidal loops the field line travels=',(phi(npt)-phi(1))/twopi

end subroutine safety_factor_a_field_line


subroutine check_whether_field_line_touch_boundary(r,z,phi,rlim,zlim,nlim,loss)
  use precision,only:p_
!  use boundary,only: nlim,rlim,zlim 
  implicit none
  real(p_),intent(in):: r,z,phi
  integer,intent(in):: nlim
  real(p_),intent(out):: rlim(nlim),zlim(nlim)
  logical,intent(out):: loss
  integer:: inout

  call PNPOLY(r,z,rlim,zlim,nlim,INOUT) !find out wheter the point (r,z) is within the limiter
  !        if (inout.eq.1) then !within the LCFS
  if (inout.eq.-1 .or.inout.eq.0) then !the particle is out of the limiter
     write(*,*) '==>This field line touches the limiter at (R,Z,phi)=', r,z,phi
     loss=.true.
     !stop
  else
     loss=.false.
  endif

end subroutine


  subroutine accumulate_tor_shift(theta_old,theta_new,r,tor_shift,real_shift) !,kt)
    use precision,only:p_
    use constants,only: twopi
    use mapping_module,only: nx_mapping,tor_shift_b,j0,r_cyl
  use interpolate_module,only: linear_1d_interpolation
    implicit none
    real(p_),intent(in):: theta_old,theta_new,r,tor_shift
    real(p_),intent(out)::real_shift
    real(p_),save::sum=0._p_
    real(p_):: tmp_array(nx_mapping),twopi_q
!integer,intent(in):: kt
    integer::kk
    if(abs(theta_old-theta_new) .ge. twopi*0.9) then !indicate finishing one poloidal loop
       do kk=1,nx_mapping
          tmp_array(kk)= tor_shift_b(kk,j0)
       enddo
       call linear_1d_interpolation(nx_mapping,r_cyl,tmp_array,r,twopi_q) !the result is twopi*q, I use this instead of directly using twopi*q because the latter may cause some cancellation problem
       sum=sum+twopi_q*sign(1._p_,theta_old-theta_new)
!write(*,*) 'sum=',sum, 'twopi_q=',twopi_q,'tor_shift=',tor_shift,'sum+tor_shift/ real_shift',sum+tor_shift, 'kt=',kt
    endif

    real_shift=sum+tor_shift


    !real_shift=sum-(twopi_q-tor_shift)
  end subroutine accumulate_tor_shift
