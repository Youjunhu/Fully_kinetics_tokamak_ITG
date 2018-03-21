module delta_ne_module
  use precision,only: p_
  use constants,only: twopi,two,pi
  use electrons_module,only:ne0 !1/m^3
  use magnetic_coordinates,only: radcor_low2,radcor_upp2
  use perturbation_field_matrix,only: ntor=>toroidal_mode_number_included
  implicit none
  private
  public delta_ne,delta_ne_theta,delta_ne_psi,delta_ne_alpha
  real(p_):: kradcor0,delta_ne0
  real(p_):: radcor_center,radcor_width
  real(p_)::mprime
  real(p_):: q,q_func_pfn,qmax,qmin
  integer:: integer_val
  external:: q_func_pfn


contains
  function delta_ne(radcor,theta,alpha) result (z) ! 1/m^3
    real(p_)::radcor,theta,alpha,z
    delta_ne0=ne0*0.001_p_ !amplitude of the initial perturbation is chosen 0.05 of the equilibrium density
    !kradcor0=twopi/(radcor_upp2-radcor_low2) !radial wave number
    kradcor0=pi/(radcor_upp2-radcor_low2) !the fundament radial wave number in sine expansion
    radcor_center=(radcor_upp2+radcor_low2)/two
    radcor_width=(radcor_upp2-radcor_low2)*0.2_p_
    q=q_func_pfn(radcor)
    !mprime=mod(ntor*q,1.0_p_) !this value is chosen in order to satisfy perioidic condition on (theta,alpha) plane when using the following spatial distribution
    qmax=q_func_pfn(radcor_upp2)
    qmin=q_func_pfn(radcor_low2)
    integer_val=nint(ntor*(qmax+qmin)/2)
    mprime=ntor*q-integer_val

    !  write(*,*) 'mprime=',mprime
    !  if(abs(mprime)>0.5) mprime=mprime-sign(1._p_,mprime) !choose the smallest abs(mprime) that is possible, can this help solve the radial resolution problem? to be verified, seems unlikely.
    !  write(*,*) 'mprime=',mprime, 'after'

    !z=delta_ne0 !for testing
    !z=delta_ne0*cos(kradcor0*radcor+mprime*theta+ntor*alpha) !incorrect
    z=delta_ne0*sin(kradcor0*(radcor-radcor_low2))*cos(mprime*theta+ntor*alpha)
    !z=delta_ne0*exp(-(radcor-radcor_center)**2/radcor_width**2)*cos(mprime*theta+ntor*alpha) !set exponential radial dependence so that the perturbation near the radial boundary is small
    !z=delta_ne0*cos(ntor*alpha) !for testing
  end function delta_ne

  function delta_ne_theta_orgin(radcor,theta,alpha) result (z) ! temporary, for testing
    real(p_)::radcor,theta,alpha,z
    delta_ne0=ne0*0.01_p_ !amplitude of the initial perturbation is chosen 0.05 of the equilibrium density
    kradcor0=twopi/(radcor_upp2-radcor_low2) !radial wave number
    radcor_center=(radcor_upp2+radcor_low2)/two
    radcor_width=(radcor_upp2-radcor_low2)*0.1_p_
    q=q_func_pfn(radcor)
    !mprime=mod(ntor*q,1.0_p_) !this value is chosen in order to satisfy perioidic condition on (theta,alpha) plane when using the following spatial distribution
    qmax=q_func_pfn(radcor_upp2)
    qmin=q_func_pfn(radcor_low2)
    integer_val=nint(ntor*(qmax+qmin)/2)
    mprime=ntor*q-integer_val

    z=-delta_ne0*sin(kradcor0*(radcor-radcor_low2))*sin(mprime*theta+ntor*alpha)*mprime
    !z=-delta_ne0*exp(-(radcor-radcor_center)**2/radcor_width**2)*sin(mprime*theta+ntor*alpha)*mprime
  end function


 function delta_ne_theta(radcor,theta,alpha) result (z) ! temporary, for testing
    real(p_)::radcor,theta,alpha,z
    real(p_):: del_theta
    del_theta=twopi/8._p_
    z=(delta_ne(radcor,theta+del_theta,alpha)-delta_ne(radcor,theta-del_theta,alpha))/(2*del_theta) !using numerical difference
  end function



  function delta_ne_psi(radcor,theta,alpha) result (z) ! temporary, for testing
    real(p_)::radcor,theta,alpha,z
    real(p_):: del_psi

    !z=kradcor0*delta_ne0*cos(kradcor0*(radcor-radcor_low2))*cos(mprime*theta+ntor*alpha) !wrong, there is a radial dependence in mprime
    del_psi=0.01*(radcor_upp2-radcor_low2)
    z=(delta_ne(radcor+del_psi,theta,alpha)-delta_ne(radcor-del_psi,theta,alpha))/(2*del_psi) !using numerical difference
  end function delta_ne_psi

  function delta_ne_alpha(radcor,theta,alpha) result (z) ! temporary, for testing
    real(p_)::radcor,theta,alpha,z
    delta_ne0=ne0*0.01_p_ !amplitude of the initial perturbation is chosen 0.05 of the equilibrium density
    kradcor0=twopi/(radcor_upp2-radcor_low2) !radial wave number
    radcor_center=(radcor_upp2+radcor_low2)/two
    radcor_width=(radcor_upp2-radcor_low2)*0.1_p_
    q=q_func_pfn(radcor)
    !mprime=mod(ntor*q,1.0_p_) !this value is chosen in order to satisfy perioidic condition on (theta,alpha) plane when using the following spatial distribution
    qmax=q_func_pfn(radcor_upp2)
    qmin=q_func_pfn(radcor_low2)
    integer_val=nint(ntor*(qmax+qmin)/2)
    mprime=ntor*q-integer_val

    z=-delta_ne0*sin(kradcor0*(radcor-radcor_low2))*sin(mprime*theta+ntor*alpha)*ntor
    !z=-delta_ne0*exp(-(radcor-radcor_center)**2/radcor_width**2)*sin(mprime*theta+ntor*alpha)*ntor
  end function delta_ne_alpha

  function delta_ne_alpha_tmp(radcor,theta,alpha) result (z) ! temporary, for testing
    use magnetic_coordinates,only:toroidal_range
   real(p_)::radcor,theta,alpha,z
   real(p_):: del_alpha
   del_alpha=toroidal_range*0.01_p_
    !z=-delta_ne0*sin(kradcor0*(radcor-radcor_low2))*sin(mprime*theta+ntor*alpha)*ntor
    !z=-delta_ne0*exp(-(radcor-radcor_center)**2/radcor_width**2)*sin(mprime*theta+ntor*alpha)*ntor
 z=(delta_ne(radcor,theta,alpha+del_alpha)-delta_ne(radcor,theta,alpha-del_alpha))/(2*del_alpha) !using numerical difference
  end function
end module delta_ne_module


subroutine initial_weight_i(w_i)
  use precision,only:p_
  use ions_module,only: nmarker_i, ps_vol_i 
  use ions_module,only: radcor_i,theta_i,alpha_i,v_i,ni0
  use domain_decomposition,only:myid
  implicit none
  real(p_),intent(out):: w_i(nmarker_i)
  integer:: k
  real(p_):: rannum,tmp
  integer:: iseed,next_seed
  ! real(p_),parameter:: eps=1.0d20
  real(p_):: initial_delta_f_i !function name
!!$  do k=1,nmarker_i
!!$     w_i(k)=ps_vol_i(k)*eps !initial weight of markers
!!$  enddo  

  do k=1,nmarker_i
     w_i(k)=ps_vol_i(k)*initial_delta_f_i(radcor_i(k),theta_i(k),alpha_i(k),v_i(k)) !initial weight of markers
  enddo

!!$  iseed=-(3777+myid*3) !set the iseed in different procs, when using this, it is a parallel generator, but the random numbers in different procs may be related if the iseed chosen for differnt procs is not good enough
!!$  ! now generate the random numbers
!!$  call sub_random_yj(iseed,next_seed,tmp) !just to trigger the use of the iseed, the generated random number is not used, 
!!$  do k=1,nmarker_i
!!$     call sub_random_yj(0,next_seed,rannum) !0 means using last random number as iseed 
!!$     w_i(k)=ni0*1.0d-10*(rannum-0.5_p_)*2._p_ !random, for testing
!!$  enddo


!!$if(myid.eq.1) then
!!$   do k=1,nmarker_i
!!$      write(*,*) alpha_i(k),initial_delta_f_i(radcor_i(k),theta_i(k),alpha_i(k),v_i(k)),w_i(k)
!!$   enddo
!!$endif

end subroutine initial_weight_i


subroutine initial_weight_e()
  use precision,only:p_
  use constants,only: two
  use electrons_module,only: nmarker_e,ps_vol_e !input
  use electrons_module,only: mu_e,vpar_e,radcor_e,theta_e,alpha_e !!input
  use magnetic_coordinates,only:mpoloidal,nflux,radcor_1d_array,theta_1d_array !input
  use array_in_mc,only: b_mc_matrix !input
  use electrons_module,only: w_e !output
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  integer:: k
!real(p_),parameter:: eps=1.0d20
real(p_):: v,b_val
!real(p_):: fe2_arbitrary
real(p_):: initial_delta_f_e

!!$  do k=1,nmarker_e
!!$     w_e(k)=ps_vol_e(k)*eps !initial weight of markers
!!$  enddo  

  do k=1,nmarker_e
  call linear_2d_interpolation(mpoloidal,nflux,theta_1d_array,radcor_1d_array,b_mc_matrix,theta_e(k),radcor_e(k),b_val)
     v=sqrt(two*mu_e(k)*b_val*2+vpar_e(k)**2) !in unit of vn_e
!     w_e(k)=ps_vol_e(k)*(mu_e(k)*b_val*2+vpar_e(k)**2)/3.0*fe2_arbitrary(radcor_e(k),theta_e(k),alpha_e(k),v) !initial weight of markers, for testing
     w_e(k)=ps_vol_e(k)*initial_delta_f_e(radcor_e(k),theta_e(k),alpha_e(k),v) !initial weight of markers, for testing
  enddo  

end subroutine initial_weight_e

function initial_delta_f_e(radcor,theta,alpha,v) result (z) ! for tesing Monte-Carlo integration,v in unit of vn_e, z in unit 1/(Ln**3*vn**3)
  use precision,only: p_
  use constants,only: two,twopi,kev
  use normalizing,only: vn_e,Ln
  use electrons_module,only: mass_e,te0
  use delta_ne_module,only: delta_ne
  implicit none
  real(p_):: radcor,theta,alpha,v,z
  real(p_):: v_si

!  ni0=1.0d19 !1/m^3
  v_si=v*vn_e
  z=delta_ne(radcor,theta,alpha)*sqrt((mass_e/(twopi*te0*kev))**3)*exp(-mass_e*v_si**2/(two*te0*kev))
  z=z*(vn_e**3*Ln**3)
end function


function initial_delta_f_i(radcor,theta,alpha,v) result (z) ! for tesing Monte-Carlo integration, v in unit of vn, z in unit 1/(Ln**3*vn**3)
  use precision,only: p_
  use constants,only: two,twopi,kev
  use normalizing,only: vn_i,Ln
  use ions_module,only: mass_i,ti0
  use delta_ne_module,only: delta_ne !function name
  implicit none
  real(p_):: radcor,theta,alpha,v,z
  real(p_):: v_si

  v_si=v*vn_i
  z=delta_ne(radcor,theta,alpha)*sqrt((mass_i/(twopi*ti0*kev))**3)*exp(-mass_i*v_si**2/(two*ti0*kev))
  z=z*(vn_i**3*Ln**3)
end function







