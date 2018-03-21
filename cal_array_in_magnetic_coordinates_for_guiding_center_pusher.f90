subroutine cal_array_in_magnetic_coordinates_for_guiding_center_pusher(mpoloidal,nflux)
  use precision,only:p_
  use constants,only: twopi,one
  use radial_module,only:psi_lcfs,psi_axis
  use magnetic_coordinates,only: r_mag_surf,z_mag_surf,jacobian !as intput
  use array_in_magnetic_coordinates_for_guiding_center_pusher,only: w1,w2,w3,w4,w5,w6,w7,w8,w9,w10 !as output
  use array_in_mc,only: b_mc_matrix, br_mc_matrix,bz_mc_matrix,bphi_mc_matrix,bp_mc_matrix !as output
  use domain_decomposition,only: myid
  !  use array_in_mc,only: grad_psi_matrix, grad_alpha_matrix,grad_psi_dot_grad_alpha_matrix !as output
  implicit none
  integer,intent(in):: mpoloidal,nflux
  real(p_):: br,bz,bphi,b !function names
  real(p_):: b_r,b_z,unitbr_z,unitbz_r,unitbphi_r,unitbphi_z !function names
  real(p_):: brval,bzval,bphival,bval
  real(p_):: b_rval,b_zval
  real(p_):: dradial_dr_func,dradial_dz_func !function names
  real(p_):: dradial_dr_func2,dradial_dz_func2 !function names
  real(p_):: dtheta_dr_func,dtheta_dz_func !function names
  real(p_):: ddelta_dr_func,ddelta_dz_func,ddelta_dr_lsf_midplane_twopi,ddelta_dz_lsf_midplane_twopi !function names
  real(p_):: dradial_dr_val,dradial_dz_val
  real(p_):: dtheta_dr_val,dtheta_dz_val
  real(p_):: ddelta_dr_val,ddelta_dz_val
  real(p_):: unitbr,unitbz,unitbphi
  real(p_):: curl_unitb_rcomp,curl_unitb_zcomp,curl_unitb_phicomp
  real(p_):: unitb_dot_curl_unitb
  real(p_):: r,z
  integer:: i,j,file_unit
  character(8)::filename


  allocate(w1(mpoloidal,nflux))
  allocate(w2(mpoloidal,nflux))
  allocate(w3(mpoloidal,nflux))
  allocate(w4(mpoloidal,nflux))
  allocate(w5(mpoloidal,nflux))
  allocate(w6(mpoloidal,nflux))
  allocate(w7(mpoloidal,nflux))
  allocate(w8(mpoloidal,nflux))
  allocate(w9(mpoloidal,nflux))
  allocate(w10(mpoloidal,nflux))
  allocate(b_mc_matrix(mpoloidal,nflux))
  allocate(br_mc_matrix(mpoloidal,nflux))
  allocate(bz_mc_matrix(mpoloidal,nflux))
  allocate(bphi_mc_matrix(mpoloidal,nflux))
  allocate(bp_mc_matrix(mpoloidal,nflux))

   if(myid.eq.0) open(81,file='grad_theta_alpha') !for testing
  do i=1,mpoloidal
     !     do j=2,nflux-2
     do j=1,nflux
        r=r_mag_surf(i,j)
        z=z_mag_surf(i,j)
        !        write(*,*) 'i,j=',i,j,'r,z=',r*ln,z*ln
        brval=br(r,z)
        bzval=bz(r,z)
        bphival=bphi(r,z)
        bval=sqrt(brval**2+bzval**2+bphival**2)
        b_zval=b_z(r,z)
        b_rval=b_r(r,z)
        unitbr=brval/bval
        unitbz=bzval/bval
        unitbphi=bphival/bval
        curl_unitb_rcomp=-unitbphi_z(r,z)
        curl_unitb_phicomp=unitbr_z(r,z)-unitbz_r(r,z)
        curl_unitb_zcomp=unitbphi_r(r,z)+unitbphi/r

        unitb_dot_curl_unitb=unitbr*curl_unitb_rcomp +unitbphi*curl_unitb_phicomp&
             & +unitbz*curl_unitb_zcomp

        dradial_dr_val=dradial_dr_func(r,z) !can be replaced by data already known in magnetic coordinates
        dradial_dz_val=dradial_dz_func(r,z) !can be replaced by data already known in magnetic coordinates
        dtheta_dr_val=dtheta_dr_func(r,z) !can be replaced by data already known in magnetic coordinates
        dtheta_dz_val=dtheta_dz_func(r,z) !can be replaced by data already known in magnetic coordinates
        ddelta_dr_val=ddelta_dr_func(r,z) !can be replaced by data already known in magnetic coordinates, do this later
        ddelta_dz_val=ddelta_dz_func(r,z) !can be replaced by data already known in magnetic coordinates

!        if(i.eq.mpoloidal)    ddelta_dr_val=ddelta_dr_lsf_midplane_twopi(r) !at theta=twopi cut
!        if(i.eq.mpoloidal)    ddelta_dz_val=ddelta_dz_lsf_midplane_twopi(r) !at theta=twopi cut

!!$!write(*,*) 'i,j=',dradial_dr_val, dradial_dr_func2(r,z),dradial_dz_val, dradial_dz_func2(r,z)
        w1(i,j)=unitb_dot_curl_unitb/bval
        w2(i,j)=unitbr*dtheta_dr_val+unitbz*dtheta_dz_val
        !      w2(i,j)=-(psi_lcfs-psi_axis)/(bval*jacobian(i,j))
        w3(i,j)=curl_unitb_rcomp/bval*dradial_dr_val+curl_unitb_zcomp/bval*dradial_dz_val
        w4(i,j)=curl_unitb_rcomp/bval*dtheta_dr_val+curl_unitb_zcomp/bval*dtheta_dz_val
        w5(i,j)=curl_unitb_phicomp/(bval*r)-curl_unitb_rcomp/bval*ddelta_dr_val-curl_unitb_zcomp/bval*ddelta_dz_val
        w6(i,j)=(bphival*b_zval)*dradial_dr_val+(-bphival*b_rval)*dradial_dz_val
        w6(i,j)=w6(i,j)/bval**2
        w7(i,j)=(bphival*b_zval)*dtheta_dr_val+(-bphival*b_rval)*dtheta_dz_val
        w7(i,j)=w7(i,j)/bval**2
        w8(i,j)=(bzval*b_rval-brval*b_zval)/r-(bphival*b_zval)*ddelta_dr_val-(-bphival*b_rval)*ddelta_dz_val
        w8(i,j)=w8(i,j)/bval**2

        w9(i,j)=unitbr*b_rval+unitbz*b_zval
        w10(i,j)=(b_rval*curl_unitb_rcomp+b_zval*curl_unitb_zcomp)/bval

!!$        grad_psi_matrix(i,j)=sqrt(dradial_dr_val**2+dradial_dz_val**2) !now calculated in another subroutine
!!$        grad_alpha_matrix(i,j)=sqrt(one/r**2+ddelta_dr_val**2+ddelta_dz_val**2)
!!$        grad_psi_dot_grad_alpha_matrix(i,j)=dradial_dr_val*ddelta_dr_val+dradial_dz_val*ddelta_dz_val

        b_mc_matrix(i,j)=bval
        br_mc_matrix(i,j)=brval
        bz_mc_matrix(i,j)=bzval
        bphi_mc_matrix(i,j)=bphival
        bp_mc_matrix(i,j)=sqrt(brval**2+bzval**2)

        if(myid.eq.0) write(81,*) r_mag_surf(i,j), z_mag_surf(i,j), dtheta_dr_val,dtheta_dz_val,&
             & sqrt(dtheta_dz_val**2+dtheta_dr_val**2), ddelta_dr_val,  ddelta_dz_val,&
             sqrt(ddelta_dr_val**2+ddelta_dz_val**2), sqrt(one/r**2+ddelta_dr_val**2+ddelta_dz_val**2)

!!$        if (isnan(w6(i,j))) then
!!$           write(*,*) 'i=,j=',i,j, w6
!!$           stop
!!$        endif
     enddo
     if(myid.eq.0) write(81,*)
  enddo
  if(myid.eq.0) close(81)
!!$if(myid.eq.0) then
!!$  write(filename,'(a4,i4.4)') 'metr',myid
!!$  file_unit=myid+311
!!$  open(file_unit,file=filename)
!!$  do i=1,mpoloidal
!!$     do j=1,nflux
!!$        write(file_unit,'(2i8.4,3(1pe14.5))')  i,j,grad_alpha_matrix(i,j),grad_psi_dot_grad_alpha_matrix(i,j), grad_psi_matrix(i,j)
!!$     enddo
!!$     write(file_unit,*)
!!$  enddo
!!$  close(file_unit)
!!$endif

  ! write(*,*) 'min max w1-4=', maxval(w1),minval(w1),maxval(w2),minval(w2),maxval(w3),minval(w3),maxval(w4),minval(w4)
  !write(*,*) 'min max w5-8=', maxval(w5),minval(w5),maxval(w6),minval(w6),maxval(w7),minval(w7),maxval(w8),minval(w8)
  !write(*,*) 'min max w9-10=', maxval(w9),minval(w9),maxval(w10),minval(w10)
end subroutine cal_array_in_magnetic_coordinates_for_guiding_center_pusher

function b_mc_func(theta,radcor) result(funcval)
  use precision,only:p_
  use array_in_mc,only: b_mc_matrix
  use magnetic_coordinates,only:mpoloidal,nflux,theta_1d_array,radcor_1d_array
  use interpolate_module,only: linear_2d_interpolation

  implicit none
  real(p_)::theta,radcor,funcval
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,b_mc_matrix,theta,radcor,funcval)
end function

function br_mc_func(theta,radcor) result(funcval)
  use precision,only:p_
  use array_in_mc,only: br_mc_matrix
  use magnetic_coordinates,only:mpoloidal,nflux,theta_1d_array,radcor_1d_array
  use interpolate_module,only: linear_2d_interpolation

  implicit none
  real(p_)::theta,radcor,funcval
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,br_mc_matrix,theta,radcor,funcval)
end function

function bz_mc_func(theta,radcor) result(funcval)
  use precision,only:p_
  use array_in_mc,only: bz_mc_matrix
  use magnetic_coordinates,only:mpoloidal,nflux,theta_1d_array,radcor_1d_array
  use interpolate_module,only: linear_2d_interpolation

  implicit none
  real(p_)::theta,radcor,funcval
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,bz_mc_matrix,theta,radcor,funcval)
end function


function bphi_mc_func(theta,radcor) result(funcval)
  use precision,only:p_
  use array_in_mc,only: bphi_mc_matrix
  use magnetic_coordinates,only:mpoloidal,nflux,theta_1d_array,radcor_1d_array
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  real(p_)::theta,radcor,funcval
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,bphi_mc_matrix,theta,radcor,funcval)
end function
