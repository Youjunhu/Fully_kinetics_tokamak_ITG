subroutine prepare_ion_source_terms()
  use precision,only:p_
  use constants,only: one,two
  use magnetic_coordinates,only:nflux,mpoloidal,r_mag_surf,dtheta,j_low2,radcor_1d_array,theta_1d_array
  use magnetic_coordinates,only:nflux2,mtoroidal !,radcor_1d_array2
  use array_in_mc,only:grad_psi_matrix,grad_alpha_matrix,grad_psi_dot_grad_alpha_matrix
  use array_in_mc,only:grad_psi_r_matrix,grad_psi_z_matrix,grad_alpha_r_matrix,grad_alpha_z_matrix
  use perturbation_field_matrix,only: jr=>jr_left,jphi=>jphi_left,jz=>jz_left
  use perturbation_field_matrix,only: source_i1,source_i2,source_i3 !as output
  use domain_decomposition,only:numprocs,myid,theta_interval,theta_start
  use flux_tube_model,only:radcor_fixed !the radcor of the center of computational region
  use interpolate_module
  implicit none
  real(p_):: jx(mtoroidal,nflux2),jy(mtoroidal,nflux2)
  real(p_):: theta,radcor,rval
  real(p_):: grad_psi_r_val,grad_psi_z_val,grad_alpha_r_val,grad_alpha_z_val
  real(p_):: grad_psi_val,grad_alpha_val,grad_psi_dot_grad_alpha_val
  real(p_):: gs_psi_prime
  integer:: i,j,ipoloidal

  theta=theta_start
  radcor=radcor_fixed
  do j=1,nflux2
     !radcor=radcor_1d_array2(j) !for global simulation
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
          & grad_psi_r_matrix,theta,radcor,grad_psi_r_val)  !in magnetic coordinates
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
          & grad_psi_z_matrix,theta,radcor,grad_psi_z_val)  !in magnetic coordinates
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
          & grad_alpha_r_matrix,theta,radcor,grad_alpha_r_val)  !in magnetic coordinates
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
          & grad_alpha_z_matrix,theta,radcor,grad_alpha_z_val)  !in magnetic coordinates
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
          & grad_psi_matrix,theta,radcor,grad_psi_val)  !in magnetic coordinates
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
          & grad_alpha_matrix,theta,radcor,grad_alpha_val)  !in magnetic coordinates
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
          & grad_psi_dot_grad_alpha_matrix,theta,radcor,grad_psi_dot_grad_alpha_val)  !in magnetic coordinates
     call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,&
          & r_mag_surf,theta,radcor,rval)  !in magnetic coordinates

     do i=1,mtoroidal
        jx(i,j)=(jr(i,j)*grad_psi_r_val+jz(i,j)*grad_psi_z_val) !jx is defined by jx=J_dot_grad_x      
        jy(i,j)=jphi(i,j)/rval+jr(i,j)*grad_alpha_r_val+jz(i,j)*grad_alpha_z_val !jy is defined by jy=J_dot_grad_y
        source_i1(i,j)=gs_psi_prime(radcor)*(jx(i,j)*grad_psi_dot_grad_alpha_val-jy(i,j)*grad_psi_val**2)
        source_i2(i,j)=gs_psi_prime(radcor)*(jx(i,j)*grad_alpha_val**2-jy(i,j)*grad_psi_dot_grad_alpha_val)
        source_i3(i,j)=0._p_ !ions do not contribute source term to the parallel field equation
     enddo
  enddo
end subroutine prepare_ion_source_terms


subroutine prepare_electron_source_terms()
  use precision,only:p_
  use normalizing,only: vn_e,vn_i,ln,bn
  use electrons_module,only:charge_e,ne0
  use constants,only: one,two,twopi,mu0
  use perturbation_field_matrix,only: pper_e=>pper_e_left, ppar_e=>ppar_e_left
  use perturbation_field_matrix,only: source_e1,source_e2,source_e3 !as output
  use magnetic_coordinates,only: nflux,mpoloidal,mtoroidal,radcor_1d_array,theta_1d_array,nflux2,dtheta,dradcor,dtor,jacobian
  use array_in_magnetic_coordinates_for_guiding_center_pusher,only: w9
  use array_in_mc,only: grad_psi_matrix, grad_alpha_matrix,grad_psi_dot_grad_alpha_matrix,b_mc_matrix
  use domain_decomposition,only: theta_start,myid,numprocs,theta_interval
  use flux_tube_model,only: radcor_fixed
  use derivatives_in_field_line_following_coordinates,only: radial_derivative,toroidal_derivative,theta_derivative
  use interpolate_module
  !  use flux_tube_model,only: j_fixed
  implicit none

  real(p_):: pper_e_x(mtoroidal,nflux2),pper_e_y(mtoroidal,nflux2) !derivatives of electron perpendicular pressure with respect to x and y, dpper_e/dx, dpper_e/dy
  real(p_):: ppar_e_theta(mtoroidal,nflux2),ppar_e_par(mtoroidal,nflux2)
  integer:: i,j,i_left,i_right,j_left,j_right
  integer::ipoloidal
  real(p_):: coeff,factor,bval,theta0
  real(p_):: grad_psi_val,grad_alpha_val, grad_psi_dot_grad_alpha_val,w9val,jacobian_val
  real(p_):: gs_psi_prime !function name
!!$do i=1,mtoroidal
!!$do j=1,nflux2
!!$if(abs(pper_e(i,j))>0.0000001) write(*,*) 'pper_e is nonzero', pper_e(i,j)
!!$enddo
!!$enddo

call radial_derivative  (pper_e,pper_e_x,dradcor)
call toroidal_derivative(pper_e,pper_e_y,dtor)


!!$  do i=1,mtoroidal   !calculate the derivative of pper_e with respect to y:
!!$     do j=1,nflux2
!!$        i_left=i-1
!!$        if(i.eq.1) i_left=mtoroidal !periodic boundary condition
!!$        i_right=i+1
!!$        if(i.eq.mtoroidal) i_right=1 !periodic boundary condition
!!$        pper_e_y(i,j)=(pper_e(i_right,j)-pper_e(i_left,j))/(two*dtor)
!!$     enddo
!!$  enddo

  call theta_derivative(ppar_e,ppar_e_theta,mtoroidal,nflux2)
  factor=-vn_e/vn_i*charge_e*mu0*ln*vn_e*ne0/bn

  theta0=theta_start+theta_interval*0.5_p_
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_psi_matrix,&
       & theta0,radcor_fixed,grad_psi_val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_alpha_matrix,&
       & theta0,radcor_fixed,grad_alpha_val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,grad_psi_dot_grad_alpha_matrix,&
       & theta0,radcor_fixed,grad_psi_dot_grad_alpha_val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,b_mc_matrix,theta0,radcor_fixed,bval)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,w9,theta0,radcor_fixed,w9val)
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,jacobian,theta0,radcor_fixed,jacobian_val) 

  do j=1,nflux2    !for flux tube model, fixing the radial location to j_fixed to evaluate the equilibrium quantities
     do i=1,mtoroidal
        source_e1(i,j)=pper_e_x(i,j)*grad_psi_val**2&
             & +pper_e_y(i,j)*grad_psi_dot_grad_alpha_val
        source_e2(i,j)=pper_e_x(i,j)*grad_psi_dot_grad_alpha_val&
             &+ pper_e_y(i,j)*grad_alpha_val**2

        ppar_e_par(i,j)=-gs_psi_prime(radcor_fixed)/(bval*jacobian_val)*ppar_e_theta(i,j)
        source_e3(i,j)=factor*(w9val*pper_e(i,j)/bval+ppar_e_par(i,j)-ppar_e(i,j)/bval*w9val) !parallel term, partially implemented
        !        source_e3(i,j)=factor*(w9val*pper_e(i,j)/bval-ppar_e(i,j)/bval*w9val) !parallel term, partially implemented
        !source_e3(i,j)=factor*(w9val*pper_e(i,j)/bval) !parallel term, partially implemented
     enddo
  enddo

  coeff=-vn_e/vn_i/twopi
  source_e1=source_e1*coeff
  source_e2=source_e2*coeff
  source_e3=source_e3
end subroutine prepare_electron_source_terms

subroutine prepare_source_term_due_to_Faraday_law()
  !source due to Faraday's law in the perpendicular field equations. in the high-n approximation, only delta_B_parallel contribute to this source
  use precision,only:p_
  use constants,only: one,two,mu0
  use electrons_module,only: charge_e,ne0
  use normalizing,only:Ln,bn,vn_i
  use magnetic_coordinates,only: mtoroidal,nflux2,theta_1d_array,dtor,dradcor,dtheta
  use perturbation_field_matrix,only: mf_par=>mf_par_left
  use perturbation_field_matrix,only:source_faraday1,source_faraday2,source_faraday3 !as output
  use array_in_mc,only: grad_psi_matrix, grad_alpha_matrix,grad_psi_dot_grad_alpha_matrix,b_mc_matrix
  use domain_decomposition,only: myid,numprocs,theta_start 
  use flux_tube_model,only: j_fixed
  use derivatives_in_field_line_following_coordinates,only: radial_derivative,toroidal_derivative
  implicit none
  integer:: i,j,ipoloidal,j_left,j_right,i_left,i_right
  real(p_):: mf_par0(mtoroidal,nflux2),mf_par_x(mtoroidal,nflux2),mf_par_y(mtoroidal,nflux2)
  real(p_):: bval

  !  call partial_derivative_in_mc0(mtoroidal,nflux2,mf_par,dtor,dradcor,mf_par_x,mf_par_y) !wrong, the shape of the matrix is not consistent between actual argument and formal argument, a bug which it took almost one day for me to find

  do j=1,nflux2
     do i=1,mtoroidal   
        mf_par0(i,j)=mf_par(i,j) !select part of the 2d array mf_par(mtoroidal+1,nflux2)
     enddo
  enddo

  call radial_derivative  (mf_par0,mf_par_x,dradcor)
  call toroidal_derivative(mf_par0,mf_par_y,dtor)

  ipoloidal=1+nint((theta_start-theta_1d_array(1))/dtheta) !index of begining theta angle of this subdomain in the original magnetic coordinate grids
  bval=b_mc_matrix(ipoloidal,j_fixed)
  do j=1,nflux2    !for flux tube model, fixing the radial location to j_fixed to evaluate the equilibrium quantities
     do i=1,mtoroidal
        source_faraday1(i,j)=mf_par_x(i,j)*bval*grad_psi_matrix(ipoloidal,j_fixed)**2 &
             & +mf_par_y(i,j)*bval*grad_psi_dot_grad_alpha_matrix(ipoloidal,j_fixed)&
             & +mf_par(i,j)*0._p_ !presently the last term is neglected
        source_faraday2(i,j)= mf_par_x(i,j)*bval*grad_psi_dot_grad_alpha_matrix(ipoloidal,j_fixed)&
             & +mf_par_y(i,j)*bval*grad_alpha_matrix(ipoloidal,j_fixed)**2&
             & +mf_par(i,j)*0._p_ !presently the last term is neglected
        source_faraday3(i,j)=0._p_ !parallel term
     enddo
  enddo
  source_faraday1= source_faraday1*(-bn/(mu0*vn_i*charge_e*ne0*ln))
  source_faraday2= source_faraday2*(-bn/(mu0*vn_i*charge_e*ne0*ln))
  source_faraday3= source_faraday3*(-bn/(mu0*vn_i*charge_e*ne0*ln))
end subroutine prepare_source_term_due_to_Faraday_law

subroutine prepare_total_source_terms()
  use precision,only:p_
  use constants,only: one,two
  !use perturbation_field_matrix,only: jper_x_i=>jper_x_i_left,jper_y_i=>jper_y_i_left
  use perturbation_field_matrix,only:source_e1,source_e2,source_e3 !as input
  use perturbation_field_matrix,only: source_i1,source_i2,source_i3 !as input
  use perturbation_field_matrix,only:source_faraday1,source_faraday2,source_faraday3
  use perturbation_field_matrix,only: source1,source2,source3 !as output
  use magnetic_coordinates,only: mtoroidal,nflux2
  use domain_decomposition,only: myid,numprocs
  implicit none

  integer:: i,j

  do i=1,mtoroidal
     do j=1,nflux2
!!$        source1(i,j)=jper_x_i(i,j) 
!!$        source2(i,j)=jper_y_i(i,j) 
!!$        source3(i,j)=0._p_         

!!$        source1(i,j)= source_faraday1(i,j)
!!$        source2(i,j)= source_faraday2(i,j)
!!$        source3(i,j)= source_faraday3(i,j)

!!$        source1(i,j)= source_faraday1(i,j)+source_e1(i,j)
!!$        source2(i,j)= source_faraday2(i,j) +source_e2(i,j)
!!$        source3(i,j)= source_faraday3(i,j) +source_e3(i,j)

!!$        source1(i,j)= source_faraday1(i,j)+source_e1(i,j)+jper_x_i(i,j)
!!$        source2(i,j)= source_faraday2(i,j)+source_e2(i,j)+jper_y_i(i,j)
!!$        source3(i,j)= source_faraday3(i,j)+source_e3(i,j)+0._p_

!!$        source1(i,j)= source_e1(i,j) +source_i1(i,j)
!!$        source2(i,j)= source_e2(i,j) +source_i2(i,j)
!!$        source3(i,j)= source_e3(i,j) +source_i3(i,j)

        source1(i,j)= source_i1(i,j) +source_e1(i,j)+source_faraday1(i,j)
        source2(i,j)= source_i2(i,j) +source_e2(i,j)+source_faraday2(i,j)
        source3(i,j)= source_i3(i,j) +source_e3(i,j)+source_faraday3(i,j)

!!$        source1(i,j)=source_e1(i,j)
!!$        source2(i,j)=source_e2(i,j)
!!$        source3(i,j)=source_e3(i,j)

!!$        source1(i,j)=jper_x_i(i,j) 
!!$        source2(i,j)=jper_y_i(i,j) 
!!$        source3(i,j)=0._p_         
     enddo
  enddo
end subroutine prepare_total_source_terms


