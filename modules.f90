module  precision
implicit none
!integer,parameter:: p_=kind(1.0)
integer,parameter:: p_=kind(1.0d0)
end module precision

module constants
  use precision,only:p_
  implicit none
  real(p_),parameter:: coulomb_log=15._p_ !assumed to be a constant
  real(p_),parameter:: kev=1.6022d-16   !unit J
  real(p_),parameter:: elementary_charge=1.6022d-19   !unit C
  real(p_),parameter:: electron_mass=9.1094d-31 !in unit of kg
  real(p_),parameter:: epsilon0=8.8542d-12 !Permittivity of free space 
  real(p_),parameter:: atom_mass_unit=1.6726d-27  !kg
  real(p_),parameter:: pi=3.1415926_p_
  real(p_),parameter:: twopi=pi*2.0_p_
  real(p_),parameter:: fourpi=pi*4.0_p_
  real(p_),parameter:: half_pi=pi*0.5_p_
  real(p_),parameter:: mu0=fourpi*1.0d-7 !permeability in SI unit
  real(p_),parameter:: zero=0.0_p_
  real(p_),parameter:: one=1.0_p_
  real(p_),parameter:: two=2.0_p_
  real(p_),parameter:: three=3.0_p_
  real(p_),parameter:: four=4.0_p_
  real(p_),parameter:: five=5.0_p_
  real(p_),parameter:: six=6.0_p_
  real(p_),parameter:: seven=7.0_p_
  real(p_),parameter:: eight=8.0_p_
  real(p_),parameter:: nine=9.0_p_
  real(p_),parameter:: ten=10.0_p_
  real(p_),parameter:: eleven=11.0_p_
  real(p_),parameter:: twelve=12.0_p_
  real(p_),parameter:: thirteen=13.0_p_
  real(p_),parameter:: fourteen=14.0_p_
  real(p_),parameter:: fifteen=15.0_p_
  real(p_),parameter:: one_half=0.5_p_
  real(p_),parameter:: one_third=one/three
  real(p_),parameter:: one_fifth=0.2_p_
  real(p_),parameter:: three_halfs=1.5_p_
end module constants

module normalizing
  use precision,only:p_
  implicit none
  real(p_):: Ln !a characteristic length (in unit of meter) chosen by users
  real(p_):: bn !a characteristic magnetic field strength (in unit of Tesla) chosen by users
  real(p_):: omegan_i,tn_i,vn_i !omegan_i is the ions gyro-angular-frequency in magnetic field bn
  real(p_):: omegan_e,tn_e,vn_e !omegan_e is the electron gyro-angular-frequency in magnetic field bn
end module normalizing


module  poloidal_flux_2d
  !poloidal flux and its partial derivatives on (R,Z) plane, where (R,phi,Z) is the cylindrical coordinate system
  use precision,only:p_
  implicit none
  integer:: nx,nz !nx,nz are respectively the numbers of grids in R and Z directions (specified in G-file)
  real(p_),dimension(:),allocatable ::xarray,zarray ! R and Z array
  real(p_),dimension(:,:),allocatable ::psi,y2a_psi !psi, 2D array on (R,Z) plane, is the poloidal flux function appearing in Grad-Shafranov equation. 
  !psi is related to poloidal magnetic field by Bp=\nabal{psi}\times\nabla{fai},where fai is the usual toroidal angle;  y2a_psi is the 2nd order erivative array used in the 2D spline interpolation of psi
  real(p_),dimension(:,:),allocatable ::psi_x,psi_z,psi_xx,psi_zz,psi_xz,psi_zx !psi_x is the partial derivative with respect to x, similar meaning for others
  real(p_),dimension(:,:),allocatable ::psi_gradient,y2a_gradient !strength of the gradient of the poloidal flux, y2a_gradient is the second derivative array used in 2D spline interpolation
  real(p_),dimension(:,:),allocatable ::y2a_psi_x,y2a_psi_z,y2a_psi_xx,y2a_psi_zz,y2a_psi_xz,y2a_psi_zx ! y2a_* is the second derivative array used in 2D spline interpolation
end module poloidal_flux_2d


module radial_module
  use precision,only:p_
  implicit none
  integer:: npsi
  real(p_),dimension(:),allocatable:: psi_1d,fpsi,ffprime,fprime,qpsi
  real(p_),dimension(:),allocatable:: pfn_npsi,tfn_npsi
  real(p_):: r_axis,z_axis,baxis
  real(p_):: psi_axis,psi_lcfs
  real(p_):: sign_bphi !sign of the toroidal component of the magnetic field
  real(p_),dimension(:),allocatable:: q_with_sign
end module radial_module


module magnetic_coordinates
  use precision,only:p_
  implicit none
  integer:: mpoloidal,nflux
  real(p_),dimension(:,:),allocatable:: r_mag_surf, z_mag_surf !SI units
  real(p_),dimension(:,:),allocatable:: tor_shift_mc
  real(p_),dimension(:,:),allocatable:: jacobian,dv !SI units
  real(p_)::sign_of_jacobian,sign_of_gs_psi_prime
  real(p_),dimension(:),  allocatable:: pfn,gs_psi_array,radcor_1d_array,minor_r_array,minor_r_prime_array !gs_psi_array is Grad-Shafranov poloidal flux in SI units
  real(p_),dimension(:),  allocatable:: theta_1d_array
  real(p_):: pfn_inner, pfn_bdry
  real(p_):: dtheta,dradcor
  real(p_)::radcor_low0,radcor_upp0
  real(p_)::radcor_low1,radcor_upp1
  real(p_)::radcor_low2,radcor_upp2
  integer:: j_low2,j_upp2,nflux2
  integer:: j_low1,j_upp1,nflux1
  integer:: i_theta_zero !poloidal index of 2D array corresponding to the low-field-side midplane
  real(p_),dimension(:),  allocatable:: radcor_1d_array2 !shranked radial array
  integer::  mtoroidal !grid number along the toroidal direction
  integer:: nsegment !computational region is 1/nsegment of the full torus
  real(p_):: toroidal_range !toroidal_range=twopi/nsegment
  real(p_),dimension(:),  allocatable:: tor_1d_array !toroidal grids for perturbations
  real(p_):: dtor
  real(p_):: abs_jacobian_min,abs_jacobian_max, vol0,vol1,vol2
end  module magnetic_coordinates

module array_in_mc
  use precision,only:p_
  implicit none
  real(p_),dimension(:,:),allocatable::rth,zth, rpsi,zpsi
  real(p_),dimension(:,:),allocatable:: bp_mc_matrix, b_mc_matrix !normalized by Bn
  real(p_),dimension(:,:),allocatable:: br_mc_matrix,bz_mc_matrix,bphi_mc_matrix !normalized by Bn
  real(p_),dimension(:,:),allocatable:: grad_psi_matrix, grad_alpha_matrix,grad_theta_matrix,&
       & grad_psi_dot_grad_alpha_matrix,grad_psi_dot_grad_theta_matrix,grad_alpha_dot_grad_theta_matrix
  real(p_),dimension(:,:),allocatable::grad_psi_r_matrix,grad_psi_z_matrix,grad_theta_r_matrix,grad_theta_z_matrix
  real(p_),dimension(:,:),allocatable::grad_alpha_r_matrix,grad_alpha_z_matrix,grad_alpha_phi_matrix
  real(p_),dimension(:,:),allocatable:: arc_length_along_field_line
end module array_in_mc

module domain_decomposition
  use precision,only:p_
  implicit none
  integer::numprocs,myid
  INTEGER :: GRID_COMM,TUBE_COMM
  integer:: GCLR,TCLR,ntube,my_left,my_right
  integer:: GCLR_cut !the value of GCLR at theta cut
  real(p_):: theta_interval,theta_start !the 1D domain decomposion along theta coordinates, theta_start is the starting poloidal location of the cell treated by No. myid proc
end module domain_decomposition

module control_parameters
  use precision,only:p_
 implicit none
  integer:: kstart,kend,niter
  real(p_):: dtao_omega_i_axis
  character(100):: poloidal_angle_type
integer:: ion_spatial_loading_scheme  !1=>uniform in (psi,theta,alpha) coordinates; 2=>uniform in real space
integer:: ion_velocity_loading_scheme !1=>uniform in (v,theta_v,phi_v) coordinates; 2=>Isotropic Gaussian in (v,theta_v,phi_v)
integer:: electron_spatial_loading_scheme !1=>uniform in (psi,theta,alpha) coordinates; 2=>uniform in real space
integer:: electron_velocity_loading_scheme !1=>uniform in (v,theta_v,phi_v) coordinates; 2=>Isotropic Gaussian in (v,theta_v,phi_v)
integer:: iplot_mode_structure
end module control_parameters

module perturbation_field_matrix
  use precision,only:p_
 implicit none
real(p_),dimension(:,:),allocatable:: epar_left,ex_left,ey_left !perturbed electric field at the left-grid of the cell for which the present proc is resposible
real(p_),dimension(:,:),allocatable:: mf_par_left,mf_x_left,mf_y_left !perturbed magnetic field
real(p_),dimension(:,:),allocatable:: mf_par_old,mf_x_old,mf_y_old

real(p_),dimension(:,:),allocatable:: epar_right,ex_right,ey_right !perturbed electric field at the right-grid of the cell for which the present proc is resposible
real(p_),dimension(:,:),allocatable:: mf_par_right,mf_x_right,mf_y_right !perturbed magnetic field

real(p_),dimension(:,:),allocatable:: ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left
real(p_),dimension(:,:),allocatable:: ef_cyl_r_right,ef_cyl_z_right,ef_cyl_phi_right
real(p_),dimension(:,:),allocatable:: ef_cyl_r_left_old,ef_cyl_z_left_old,ef_cyl_phi_left_old
real(p_),dimension(:,:),allocatable:: ef_cyl_r_right_old,ef_cyl_z_right_old,ef_cyl_phi_right_old

!real(p_),dimension(:,:),allocatable:: jper_x_i_left,  jper_y_i_left !perturbed perpendicular ion current at the left-grid of the cell for which the present proc is resposible 
!real(p_),dimension(:,:),allocatable:: jper_x_i_right, jper_y_i_right !perturbed perpendicular ion current at the right-grid of the cell for which the present proc is resposible 
!jper_i_x is defined as grad_x dot (B0 cross v),jper_i_y is defined as grad_y dot (B0 cross v)

real(p_),dimension(:,:),allocatable:: jr_left,jphi_left,jz_left
real(p_),dimension(:,:),allocatable:: den_left,potential !used in the adiabatic electron model
real(p_),dimension(:,:),allocatable:: pper_e_left !electron perpendicular pressure
real(p_),dimension(:,:),allocatable:: ppar_e_left !electron parallel pressure
real(p_),dimension(:,:),allocatable:: source_e1,source_e2,source_e3 !electron terms on the right-hand-side of the two perpendicular field equations and the parallel field equation
real(p_),dimension(:,:),allocatable:: source_i1,source_i2,source_i3
real(p_),dimension(:,:),allocatable:: source1, source2,source3
real(p_),dimension(:,:),allocatable:: source_faraday1, source_faraday2,source_faraday3

real(p_)::  coeff1_ex, coeff1_ey, coeff1_ex_yy, coeff1_ex_xy, coeff1_ey_xx, coeff1_ey_xy !coefficients of the 1nd perpendicular equation
real(p_)::  coeff2_ex, coeff2_ey, coeff2_ex_yy, coeff2_ex_xy, coeff2_ey_xx, coeff2_ey_xy !coefficients of the 2nd perpendicular equation
real(p_):: coeff3_epar_xx, coeff3_epar_yy,coeff3_epar_xy,&
     & coeff3_epar_x,coeff3_epar_y,coeff3_epar !coefficients of the parallel field equation
real(p_):: coeff1_ey_implicit, coeff2_ex_implicit
logical:: filter_toroidal,filter_radial
integer:: toroidal_mode_number_included,radial_harmonics_included
end module perturbation_field_matrix

module flux_tube_model
use precision,only:p_
implicit none
integer:: j_fixed
real(p_):: radcor_fixed !the radcor of the center of computational region
end module flux_tube_model

module array_in_magnetic_coordinates_for_guiding_center_pusher
  use precision,only:p_
  implicit none
  real(p_),dimension(:,:),allocatable:: w1,w2,w3,w4,w5,w6,w7,w8,w9,w10 !for guiding-center pusher
end module array_in_magnetic_coordinates_for_guiding_center_pusher


module boundary
  use precision,only:p_
  integer:: nlim,np_lcfs
  real(p_),dimension(:),allocatable ::rlim,zlim,x_lcfs,z_lcfs
end module boundary


module mapping_module !from cylindrical coordinates to magnetic coordinates
  use precision,only:p_
  implicit none
  integer,parameter:: nx_mapping=100,nz_mapping=100
  real(p_):: r_cyl(nx_mapping),z_cyl(nz_mapping)
  real(p_):: radcor(nx_mapping,nz_mapping)
  real(p_):: theta_a(nx_mapping,nz_mapping),theta_b(nx_mapping,nz_mapping)
  real(p_):: tor_shift_a(nx_mapping,nz_mapping),tor_shift_b(nx_mapping,nz_mapping)
  real(p_):: dr,dz
  integer:: i0,j0 !index of the point at magnetic axis
  real(p_):: dtheta_dr(nx_mapping,nz_mapping),dtheta_dz(nx_mapping,nz_mapping)
  real(p_):: ddelta_dr_a(nx_mapping,nz_mapping),ddelta_dz_a(nx_mapping,nz_mapping)
  real(p_):: ddelta_dr_b(nx_mapping,nz_mapping),ddelta_dz_b(nx_mapping,nz_mapping)
  real(p_):: dradial_dr(nx_mapping,nz_mapping),dradial_dz(nx_mapping,nz_mapping)
end module mapping_module

module ions_module
  use precision,only:p_
  implicit none
  real(p_)::mass_i,charge_i,dtao_i

  integer:: total_nmarker_i  ! total number of ion markers (including the particles in all the processors).
  integer:: nmarker_i !particle number in a single processor, its value will be differnt for differnt processors and at differnt time
  integer:: fixed_large_size
  real(p_):: ni0,ti0,kappa_ni,kappa_ti !ti0 is a typical value of ion temperature in the compuational region, in flux tube model, ti0 is the temperature at the reference magnetic surface
  real(p_):: vt_i,vmin_i,vmax_i
  real(p_),allocatable:: ps_vol_i(:) !determined by the initial loading, is constant for each marker over the time-evoultion
  real(p_)::normalizing_factor

  real(p_),allocatable:: r_i(:),z_i(:),phi_i(:) !Cylindrical coordinates at time t_{n},  unit Ln, rad
  real(p_),allocatable:: r_i_old(:),z_i_old(:),phi_i_old(:) !Cylindrical coordinates at time t_{n},  unit Ln, rad, temporary working arrays
  real(p_),allocatable:: r_i_mid(:),z_i_mid(:),phi_i_mid(:) !Cylindrical coordinates at time t_{n+1/2}, unit Ln, rad, 
  real(p_),allocatable:: radcor_i(:),theta_i(:), alpha_i(:),tor_shift_i(:) !magnetic coordinates at integer-time-step, alpha is the generalized toroidal angle
  real(p_),allocatable:: radcor_i_mid(:),theta_i_mid(:),alpha_i_mid(:) !magnetic coordinates at half time-step

  real(p_),allocatable:: vr_i(:),vz_i(:),vphi_i(:) ! projection of velocity at t_{n-1/2} to the basis cylindrical vector at integer-time-step t_{n}, unit vn_i=ln/tn_i, 
  real(p_),allocatable:: vr_i_integer_mid(:),vz_i_integer_mid(:),vphi_i_integer_mid(:) !projection of velocity at t_{n} to the basis cylindrical vector at t_{n+1/2}, unit vn_i=ln/tn_i, 

  real(p_),allocatable:: vr_i_mid(:),vz_i_mid(:),vphi_i_mid(:) !projection of velocity at t_{n+1/2} to the local basis vectors at t_{n+1/2}
  real(p_),allocatable:: vr_i_integer(:),vz_i_integer(:),vphi_i_integer(:) !projection of velocity at t_{n} to the local basis vectors at t_{n}

  real(p_),allocatable:: w_i(:)  !weight of ion markers
  real(p_),allocatable:: w_i_mid(:)  !weight of ion markers at t_{n+1/2}, presently, not actually used in the compuation 
 real(p_),allocatable:: w_i_star(:)  !weight of ion markers

  real(p_),allocatable:: v_i(:),vpar_i(:),vx_i(:),vy_i(:) !velocity components in magnetic coordinates,vx is defined by vx=v_dot_grad_x, vy is defined by vy=v_dot_grad_y,note that grad_x and grad_y are not perpendicular to each other, can be at half-time step or integer-time-step, depending on the context
  real(p_),allocatable:: grad_psi_i(:),grad_alpha_i(:), grad_psi_dot_grad_alpha_i(:),bval_i(:) ! can be at half-time step or integer-time-step, depending on the context

  logical,allocatable:: touch_bdry_i(:),active_i(:) !indicates whether the orbit of a marker touches the boundary
  logical,allocatable:: touch_bdry_i_mid(:),active_i_mid(:) !indicates whether the orbit of a marker touches the boundary
  integer:: ntouch_bdry_i=0, total_ntouch_bdry_i !numbe of markers that touch the boundary in each process and all processes, respectively
end module ions_module


module electrons_module
  use precision,only:p_
  implicit none
  logical::fluid_electron
  real(p_)::mass_e,charge_e 
  integer:: total_nmarker_e  ! total number of markers (including all the particles in all theprocessors).
  integer:: nmarker_e !particle number in a single processor, its value will be differnt for differnt processors and at differnt time
  real(p_):: ne0,te0, kappa_ne,kappa_te !te0 is a typical value of temperature in the compuational region, in flux tube model, it is the temperature at the reference magnetic surface

  real(p_),allocatable:: radcor_e(:),theta_e(:),alpha_e(:) !in magnetic coordinates, alpha_e is the generalized toroidal angle
  real(p_),allocatable::vpar_e(:),mu_e(:),w_e(:)

  real(p_),allocatable::radcor_e_mid(:),theta_e_mid(:),alpha_e_mid(:) !vlaues at half-time-step
  real(p_),allocatable::vpar_e_mid(:), w_e_mid(:)

  real(p_),allocatable:: rg_e(:),zg_e(:),phig_e(:) !in cylindrical coordinates, unit Ln, rad, for testing
  real(p_),allocatable:: tor_shift_e(:) !used when transformation from cylindrical coordinates to magnetic coordinates

  real(p_),allocatable:: ps_vol_e(:) !weight of electron markers
  logical,allocatable:: touch_bdry_e(:),active_e(:) !indicates whether the orbit of a marker touches the boundary
  integer:: ntouch_bdry_e=0, total_ntouch_bdry_e !numbe of markers that touch the boundary in each process and all processes, respectively

end module electrons_module
