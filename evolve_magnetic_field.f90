subroutine evolve_parallel_magnetic_field(dtao_i)
  use precision,only:p_
  use constants,only:two
  use perturbation_field_matrix,only:ex=>ex_left,ey=>ey_left
  use perturbation_field_matrix,only:mf_par=>mf_par_left !input and output
  use magnetic_coordinates,only: mpoloidal,nflux,radcor_1d_array,theta_1d_array,mtoroidal,nflux2,dradcor,dtor,dtheta
  use flux_tube_model,only: radcor_fixed
  use array_in_mc,only: b_mc_matrix
  use domain_decomposition,only:theta_start,myid,numprocs
  use interpolate_module,only: linear_2d_interpolation
  use derivatives_in_field_line_following_coordinates,only:radial_derivative,toroidal_derivative
  implicit none
  real(p_),intent(in):: dtao_i
  integer:: i,j,j_left,j_right,i_left,i_right !,ipoloidal
  real(p_):: ex0(mtoroidal,nflux2),ey0(mtoroidal,nflux2)
  real(p_):: ex_y(mtoroidal,nflux2),ey_x(mtoroidal,nflux2)
  real(p_):: bval,theta0,gs_psi_prime

  do j=1,nflux2
     do i=1,mtoroidal   
        ex0(i,j)=ex(i,j) !select part of the 2d array
        ey0(i,j)=ey(i,j) !select part of the 2d array
     enddo
  enddo

  call radial_derivative  (ey0,ey_x,dradcor)
  call toroidal_derivative(ex0,ex_y,dtor)

  theta0=theta_start
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,b_mc_matrix,theta0,radcor_fixed,bval)

  do i=1,mtoroidal   
     do j=1,nflux2
        mf_par(i,j)=mf_par(i,j)-dtao_i*(ey_x(i,j)-ex_y(i,j))*bval/gs_psi_prime(radcor_fixed)
     enddo
  enddo
  mf_par=0._p_ !for testing
  !if(myid.eq.0) write(*,*) ey_x(mtoroidal/3,nflux2/3),ex_y(mtoroidal/3,nflux2/3)
  !if(myid.eq.0) write(*,*) ex(mtoroidal/3,nflux2/3),ey(mtoroidal/3,nflux2/3)
end subroutine evolve_parallel_magnetic_field


subroutine evolve_perpendicular_magnetic_field(dtao_i) 
  use precision,only:p_
  use constants,only:two
  use perturbation_field_matrix,only:epar=>epar_left
  use perturbation_field_matrix,only:mf_x=>mf_x_left,mf_y=>mf_y_left !as input and output
  use magnetic_coordinates,only: radcor_1d_array,theta_1d_array,mpoloidal,nflux,mtoroidal,nflux2,dradcor,dtor
  use flux_tube_model,only: radcor_fixed
  use array_in_mc,only: b_mc_matrix,grad_psi_matrix, grad_alpha_matrix,&
       & grad_psi_dot_grad_alpha_matrix
  use domain_decomposition,only:theta0=>theta_start,myid,numprocs
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  real(p_),intent(in):: dtao_i
  real(p_):: epar_x(mtoroidal,nflux2),epar_y(mtoroidal,nflux2)
  real(p_):: rhs1(mtoroidal,nflux2),rhs2(mtoroidal,nflux2)
  real(p_):: gs_psi_prime
  real(p_):: grad_psi_val,grad_alpha_val, grad_psi_dot_grad_alpha_val,bval
  integer:: i,j
  integer:: j_left,j_right,i_left,i_right
  integer,parameter:: nvar=2 !mf_x,mf_y
  real(p_)::  a(nvar,nvar),b(nvar)
  INTEGER:: ipiv(nvar),info

  do i=1,mtoroidal   !calculate the derivative of epar with respect to x:
     do j=1,nflux2
        j_left=j-1
        if(j_left.eq.0) j_left=nflux2 !periodic boundary condition
        j_right=j+1
        if(j_right.eq.nflux2+1) j_right=1 !periodic boundary condition
        epar_x(i,j)= (epar(i,j_right)-epar(i,j_left))/(two*dradcor)
     enddo
  enddo

  do i=1,mtoroidal   !calculate the derivative of epar with respect to y:
     do j=1,nflux2
        i_left=i-1
        if(i.eq.1) i_left=mtoroidal !periodic boundary condition
        i_right=i+1
        if(i.eq.mtoroidal) i_right=1 !periodic boundary condition
        epar_y(i,j)=(epar(i_right,j)-epar(i_left,j))/(two*dtor)
     enddo
  enddo

 call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_psi_matrix,&
       & theta0,radcor_fixed,grad_psi_val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_alpha_matrix,&
       & theta0,radcor_fixed,grad_alpha_val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_psi_dot_grad_alpha_matrix,&
       & theta0,radcor_fixed,grad_psi_dot_grad_alpha_val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,b_mc_matrix,theta0,radcor_fixed,bval)
!

  do i=1,mtoroidal   
     do j=1,nflux2
        rhs1(i,j)=mf_x(i,j)*grad_psi_val**2+mf_y(i,j)*grad_psi_dot_grad_alpha_val&
             & -dtao_i*epar_y(i,j)*bval/gs_psi_prime(radcor_fixed)
        rhs2(i,j)=mf_x(i,j)*grad_psi_dot_grad_alpha_val+mf_y(i,j)*grad_alpha_val**2 &
             &+dtao_i*epar_x(i,j)*bval/gs_psi_prime(radcor_fixed)
     enddo
  enddo

 do i=1,mtoroidal   
     do j=1,nflux2
        b(1)=rhs1(i,j) !right-hand-side
        b(2)=rhs2(i,j)

        a(1,1)=grad_psi_val**2
        a(1,2)=grad_psi_dot_grad_alpha_val
        a(2,1)=grad_psi_dot_grad_alpha_val
        a(2,2)=grad_alpha_val**2

        call dgesv(nvar,1,a,nvar,ipiv,b,nvar,info) !lapack routine to slove linear system of equations, ref: http://www.netlib.org/lapack/explore-3.1.1-html/zgesv.f.html
        mf_x(i,j)=b(1)
        mf_y(i,j)=b(2)
     enddo
  enddo
end subroutine evolve_perpendicular_magnetic_field
