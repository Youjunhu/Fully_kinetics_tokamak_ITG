module connection_condition
  use interpolate_module,only: linear_1d_interpolation
  implicit none
  private
  public connection_condition_at_theta_cut,connection_condition_at_theta_cut2,connection_condition_at_theta_cut3,&
       & components_transformation_at_theta_cut

contains
  subroutine connection_condition_at_theta_cut(u) !interpolating data on theta=-pi plane to get values on theta=+pi plane
    !the two adjacent cells separated by theta=-pi or pi share the interface. However the grids on the interface of one cell is different from another,
    !specificly, the toroidal location of point (itor,jrad) on theta=+pi is different from the corresponding point (itor,jrad) on theta=-pi. Therefore an interpolation over the toroidal angle is needed.
    use precision,only:p_
    use magnetic_coordinates,only:tor_1d_array,tor_shift_mc,mpoloidal,j_low2,mtoroidal
    real(p_),intent(inout):: u(:,:)
    real(p_):: phi_old(mtoroidal+1),u_old(mtoroidal+1)
    real(p_):: phi_new
    integer:: i,j,jeq,m,n

    m=size(u,1)
    n=size(u,2)
    !write(*,*) 'm,n=',m,n
    do j=1,n
       jeq=j-1+j_low2
       do i=1,mtoroidal+1
          phi_old(i)=tor_1d_array(i)+tor_shift_mc(1,jeq) !at theta=-pi
          phi_old(i)=phi_old(i)-tor_shift_mc(1,jeq) !substracting by tor_shift_mc(1,jeq), which is a constant for the i=1,m+1 loop, is to re-define phi. If we do not re-define phi_old and use shift_to_specified_toroidal_rangle subroutine to shift the phi_old to the specified range, then phi_old array is changing in the way like 0.8pi->end_value-->starting_value--0.8pi, which is non-monotonical, which is not compatible with the interpolating function which requires the numerical table to be monotonical. The phi of the points that are to be interpolated should also be redefined in this way, i.e., subtracted by tor_shift_mc(1,jeq) , yj2017-Oct.18
          !call shift_to_specified_toroidal_range(phi_old(i)) !this shift is not necessary because phi_old array is already in the specified range after we re-define it.
          if(i.eq.m+1) then
             u_old(i)=u(1,j) !using toroidal periodic condition, u(m+1,j)=u(1,j).
          else
             u_old(i)=u(i,j)
          endif
       enddo

       do i=1,m
          phi_new=tor_1d_array(i)+tor_shift_mc(mpoloidal,jeq) !at theta=+pi
          phi_new=phi_new-tor_shift_mc(1,jeq) !re-define the toroidal angle in the same way that the above numerical table is defined, so that the following interpolating can work correctly
          call shift_to_specified_toroidal_range(phi_new) !this call is necessary
          call linear_1d_interpolation(mtoroidal+1,phi_old,u_old,phi_new,u(i,j))  
       enddo
    enddo
  end subroutine connection_condition_at_theta_cut


  subroutine connection_condition_at_theta_cut2(u) !interpolating data on theta=+pi plane to get values on theta=-pi plane
    use precision,only:p_
    use constants,only: twopi
    use magnetic_coordinates,only:tor_1d_array,tor_shift_mc,mpoloidal,j_low2,mtoroidal
    implicit none
    real(p_),intent(inout):: u(:,:)
    real(p_):: phi_old(mtoroidal+1),u_old(mtoroidal+1)
    real(p_):: phi_new
    integer:: i,j,jeq,m,n

    m=size(u,1)
    n=size(u,2)

    do j=1,n
       jeq=j-1+j_low2
       do i=1,mtoroidal+1
          phi_old(i)=tor_1d_array(i)+tor_shift_mc(mpoloidal,jeq) !at theta cut: theta=+pi
          phi_old(i)=phi_old(i)-tor_shift_mc(mpoloidal,jeq) !re-define phi_old, the reason is explained at the place labeled by yj2017-Oct.18 in this file
          if(i.eq.m+1) then
             u_old(i)=u(1,j) !using toroidal periodic condition, u(m+1,j)=u(1,j).
          else
             u_old(i)=u(i,j)
          endif
       enddo

       do i=1,m
          phi_new=tor_1d_array(i)+tor_shift_mc(1,jeq) !at theta cut: theta=-pi.
          phi_new=phi_new-tor_shift_mc(mpoloidal,jeq) !The subtraction by tor_shift_mc(mpoloidal,jeq) is to re-define phi,the reason is explained in "connection_condition_at_theta_cut"
          call shift_to_specified_toroidal_range(phi_new)
          call linear_1d_interpolation(mtoroidal+1,phi_old,u_old,phi_new,u(i,j))  
       enddo
    enddo
  end subroutine connection_condition_at_theta_cut2

  subroutine connection_condition_at_theta_cut3(u) !interpolating data on theta=+pi-theta_interval plane to get data on grids defined on the plane of theta=-pi-theta_interval. The grids on theta=-pi-theta_interval are determined by following the magnetic field line starting from the grids defined on theta=-pi plane. The intersecting points of these field lines with theta=-pi-theta_interval plane define the grids on this plane
    !see the comments in subroutine "connection_condition_at_theta_cut"
    use precision,only:p_
    use constants,only: twopi
    use magnetic_coordinates,only:tor_1d_array,tor_shift_mc,mpoloidal,dtheta,j_low2,mtoroidal
    use domain_decomposition,only: theta_interval
    implicit none
    real(p_),intent(inout):: u(:,:)
    real(p_):: phi_old(mtoroidal+1),u_old(mtoroidal+1)
    real(p_):: phi_at_midplane,phi_new
    integer:: i,j,multi_eq_cells,jeq,m,n

    m=size(u,1)
    n=size(u,2)

    multi_eq_cells=NINT(theta_interval/dtheta)
    do j=1,n
       jeq=j-1+j_low2
       do i=1,mtoroidal+1
          phi_old(i)=tor_1d_array(i)+tor_shift_mc(mpoloidal-multi_eq_cells,jeq)
          phi_old(i)=phi_old(i)-tor_shift_mc(mpoloidal-multi_eq_cells,jeq) !re-define phi_old by shifting phi_old by a contant, the reason is explained at the place labeled by mark2017-11-7-1 in this file
          if(i.eq.m+1) then
             u_old(i)=u(1,j) !using toroidal periodic condition, u(m+1,j)=u(1,j).
          else
             u_old(i)=u(i,j)
          endif
       enddo

       do i=1,m
          phi_at_midplane=tor_1d_array(i)+tor_shift_mc(1,jeq) !phi angle of point (i,j) on theta=-pi plane
          phi_new=phi_at_midplane+(tor_shift_mc(mpoloidal-multi_eq_cells,jeq)-tor_shift_mc(mpoloidal,jeq)) !the phi angle of field-lines when it gos backward to theta=-pi-theata_interval
          phi_new=phi_new-tor_shift_mc(mpoloidal-multi_eq_cells,jeq) !shift by a constant, so it is consistent with the definition of phi_old
          call shift_to_specified_toroidal_range(phi_new)
          call linear_1d_interpolation(mtoroidal+1,phi_old,u_old,phi_new,u(i,j))  
       enddo
    enddo
  end subroutine connection_condition_at_theta_cut3

  subroutine components_transformation_at_theta_cut(ex,ey) !due to the basis vector grad_alpha is discontinuous across the theta-cut.
    use precision,only:p_
    use magnetic_coordinates,only: mpoloidal
    use flux_tube_model,only:radcor_fixed,j_fixed
    use array_in_mc,only: grad_psi_matrix, grad_alpha_matrix,&
         & grad_psi_dot_grad_alpha_matrix,grad_psi_dot_grad_theta_matrix
    use array_in_mc,only:grad_psi_r_matrix,grad_psi_z_matrix, &
         & grad_alpha_r_matrix,grad_alpha_z_matrix,grad_alpha_phi_matrix
    implicit none
    real(p_),intent(inout):: ex(:,:),ey(:,:)
    !  real(p_),allocatable:: source1(:,:),source2(:,:),det_x(:,:),det_y(:,:)
    real(p_)::source1(size(ex,1),size(ex,2)),source2(size(ex,1),size(ex,2))
    real(p_)::det_x(size(ex,1),size(ex,2)),det_y(size(ex,1),size(ex,2))

    real(p_):: det_c,a11,a12,a21,a22
    !  real(p_):: theta0
    real(p_):: grad_psi_val,grad_alpha_val,grad_psi_dot_grad_alpha_val
    real(p_)::grad_alpha_dot_grad_alpha_positive,grad_psi_dot_grad_alpha_positive !_positive indicates value at theta=+pi, i.e., i=mpoloidal
    integer:: i,j 

!!$  allocate(source1(m,n))
!!$  allocate(source2(m,n))
!!$  allocate(det_x(m,n))
!!$  allocate(det_y(m,n))

    !  theta0=0._p_
!!$  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_psi_matrix,&
!!$       & theta0,radcor_fixed,grad_psi_val)
!!$  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_alpha_matrix,&
!!$       & theta0,radcor_fixed,grad_alpha_val)
!!$  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_psi_dot_grad_alpha_matrix,&
!!$       & theta0,radcor_fixed,grad_psi_dot_grad_alpha_val)

    i=1
    j=j_fixed

    grad_psi_val=grad_psi_matrix(i,j)
    grad_alpha_val=grad_alpha_matrix(i,j)
    grad_psi_dot_grad_alpha_val=grad_psi_dot_grad_alpha_matrix(i,j)

    source1= ex*grad_psi_val**2+ey*grad_psi_dot_grad_alpha_val
    source2= ex*grad_psi_dot_grad_alpha_val+ey*grad_alpha_val**2

    grad_psi_dot_grad_alpha_positive=grad_psi_r_matrix(i,j)*grad_alpha_r_matrix(mpoloidal,j)+&
         & grad_psi_z_matrix(i,j)*grad_alpha_z_matrix(mpoloidal,j)

    grad_alpha_dot_grad_alpha_positive=grad_alpha_r_matrix(i,j)*grad_alpha_r_matrix(mpoloidal,j)+&
         & grad_alpha_z_matrix(i,j)*grad_alpha_z_matrix(mpoloidal,j)+&
         & grad_alpha_phi_matrix(i,j)*grad_alpha_phi_matrix(mpoloidal,j)

    a11=grad_psi_val**2
    a12=grad_psi_dot_grad_alpha_positive
    a21=grad_psi_dot_grad_alpha_val
    a22=grad_alpha_dot_grad_alpha_positive

    det_c=a11*a22-a12*a21
    det_x=source1*a22-a12*source2
    det_y=a11*source2-source1*a21
    ex=det_x/det_c !Cramer's Rule to solve linear system of equations
    ey=det_y/det_c !Cramer's Rule to solve linear system of equations
  end subroutine components_transformation_at_theta_cut
end module connection_condition
