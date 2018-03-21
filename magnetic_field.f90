!all independent and dependent variables are in SI units

function br_SI(r,z) !R component of magnetic field
  use precision,only:p_
  implicit none
  real(p_):: br_SI,r,z
  real(p_):: psi_z_func
  br_SI=-psi_z_func(r,z)/r
end function 

function bz_SI(r,z) !Z component of magnetic field
  use precision,only:p_
  implicit none
  real(p_):: bz_SI,r,z
  real(p_):: psi_r_func
  bz_SI=psi_r_func(r,z)/r
end function 

function bphi_SI(r,z) !phi component of magnetic field
  use precision,only:p_
  implicit none
  real(p_):: bphi_SI,r,z
  real(p_)::g_func,psi_func
  bphi_SI=g_func(psi_func(r,z))/r
end function 


function b_SI(r,z) !strength of magnetic field
  use precision,only:p_
  implicit none
  real(p_):: b_SI,r,z
  real(p_):: br_SI,bz_SI,bphi_SI
  b_SI=sqrt(br_SI(r,z)**2+bz_SI(r,z)**2+bphi_SI(r,z)**2)
end function 

function b_r_SI(r,z) !partial derivative of  magnetic field
  use precision,only:p_
  use constants,only:one
  implicit none
  real(p_):: b_r_SI,r,z
  real(p_):: b_SI,psi_r_func,psi_rr_func,psi_z_func,psi_rz_func,g_func,psi_func, g_r
  b_r_SI=-b_SI(r,z)/r+one/(b_SI(r,z)*r*r)*(psi_r_func(r,z)*psi_rr_func(r,z) &
       & +psi_z_func(r,z)*psi_rz_func(r,z)+g_func(psi_func(r,z))*g_r(r,z))
end function 

function b_z_SI(r,z) !partial derivative of magnetic field
  use precision,only:p_
  use constants,only:one
  implicit none
  real(p_):: b_z_SI,r,z
  real(p_):: b_SI,psi_r_func,psi_rz_func,psi_z_func,psi_zz_func,g_func,psi_func, g_z
  b_z_SI=one/(b_SI(r,z)*r*r)*(psi_r_func(r,z)*psi_rz_func(r,z)+&
       & psi_z_func(r,z)*psi_zz_func(r,z)+g_func(psi_func(r,z))*g_z(r,z))
end function 



function br_r_SI(r,z) !partial derivative of component of magnetic field
  use precision,only:p_
  use constants,only:one
  implicit none
  real(p_):: br_r_SI,r,z
  real(p_):: psi_z_func,psi_rz_func
  br_r_SI=psi_z_func(r,z)/r**2-psi_rz_func(r,z)/r
end function 


function br_z_SI(r,z) !partial derivative of component of magnetic field
  use precision,only:p_
  use constants,only:one
  implicit none
  real(p_):: br_z_SI,r,z
  real(p_):: psi_zz_func
  br_z_SI=-psi_zz_func(r,z)/r
end function 

function bz_r_SI(r,z) !partial derivative of component of magnetic field
  use precision,only:p_
  use constants,only:one
  implicit none
  real(p_):: bz_r_SI,r,z
  real(p_):: psi_r_func,psi_rr_func
  bz_r_SI=-psi_r_func(r,z)/r**2+psi_rr_func(r,z)/r
end function 


function bz_z_SI(r,z) !partial derivative of component of magnetic field
  use precision,only:p_
  use constants,only:one
  implicit none
  real(p_):: bz_z_SI,r,z
  real(p_):: psi_rz_func
  bz_z_SI=psi_rz_func(r,z)/r
end function 


function bphi_r_SI(r,z) !partial derivative of component of magnetic field
  use precision,only:p_
  use constants,only:one
  implicit none
  real(p_):: bphi_r_SI,r,z
  real(p_):: g_func,gprime,psi_func,psi_r_func
  bphi_r_SI=gprime(psi_func(r,z))*psi_r_func(r,z)/r-g_func(psi_func(r,z))/r**2
end function 

function bphi_z_SI(r,z) !partial derivative of component of magnetic field
  use precision,only:p_
  use constants,only:one,zero
  implicit none
  real(p_):: bphi_z_SI,r,z
  real(p_):: gprime,psi_func,psi_z_func
  bphi_z_SI=gprime(psi_func(r,z))*psi_z_func(r,z)/r
end function 


!!$function unitbr_SI(r,z) !component of the unit vector of magnetic field
!!$  use precision,only:p_
!!$  implicit none
!!$  real(p_):: unitbr_SI,r,z
!!$  real(p_):: b_SI,br_SI
!!$  unitbr_SI=br_SI(r,z)/b_SI(r,z)
!!$end function 
!!$
!!$function unitbz_SI(r,z) !component of the unit vector of magnetic field
!!$  use precision,only:p_
!!$  implicit none
!!$  real(p_):: unitbz_SI,r,z
!!$  real(p_):: b_SI,bz_SI
!!$  unitbz_SI=bz_SI(r,z)/b_SI(r,z)
!!$end function 
!!$
!!$function unitbphi_SI(r,z) !component of the unit vector of magnetic field
!!$  use precision,only:p_
!!$  implicit none
!!$  real(p_):: unitbphi_SI,r,z
!!$  real(p_):: b_SI,bphi_SI
!!$  unitbphi_SI=bphi_SI(r,z)/b_SI(r,z)
!!$end function 


!!$function unitbr_r_SI(r,z) !partial derivative of the component of magnetic unit vector
!!$  use precision,only:p_
!!$  implicit none
!!$  real(p_):: unitbr_r_SI,r,z
!!$  real(p_):: b_SI,br_r_SI,b_r_SI,br_SI
!!$  unitbr_r_SI=(br_r_SI(r,z)*b_SI(r,z)-b_r_SI(r,z)*br_SI(r,z))/b_SI(r,z)**2
!!$end function 
!!$
!!$function unitbr_z_SI(r,z) !partial derivative of the component of magnetic unit vector
!!$  use precision,only:p_
!!$  implicit none
!!$  real(p_):: unitbr_z_SI,r,z
!!$  real(p_):: b_SI,br_z_SI,b_z_SI,br_SI
!!$  unitbr_z_SI=(br_z_SI(r,z)*b_SI(r,z)-b_z_SI(r,z)*br_SI(r,z))/b_SI(r,z)**2
!!$end function 
!!$
!!$function unitbz_r_SI(r,z) !partial derivative of the component of magnetic unit vector
!!$  use precision,only:p_
!!$  implicit none
!!$  real(p_):: unitbz_r_SI,r,z
!!$  real(p_):: b_SI,bz_r_SI,b_r_SI,bz_SI
!!$  unitbz_r_SI=(bz_r_SI(r,z)*b_SI(r,z)-b_r_SI(r,z)*bz_SI(r,z))/b_SI(r,z)**2
!!$end function 
!!$
!!$function unitbz_z_SI(r,z) !partial derivative of the component of magnetic unit vector
!!$  use precision,only:p_
!!$  implicit none
!!$  real(p_):: unitbz_z_SI,r,z
!!$  real(p_):: b_SI,bz_z_SI,b_z_SI,bz_SI
!!$  unitbz_z_SI=(bz_z_SI(r,z)*b_SI(r,z)-b_z_SI(r,z)*bz_SI(r,z))/b_SI(r,z)**2
!!$end function 
!!$
!!$function unitbphi_r_SI(r,z) !partial derivative of the component of the unit vector of the magnetic field
!!$  use precision,only:p_
!!$  implicit none
!!$  real(p_):: unitbphi_r_SI,r,z
!!$  real(p_):: b_SI,bphi_SI,b_r_SI,bphi_r_SI
!!$  unitbphi_r_SI=(bphi_r_SI(r,z)*b_SI(r,z)-b_r_SI(r,z)*bphi_SI(r,z))/b_SI(r,z)**2
!!$end function 
!!$
!!$
!!$function unitbphi_z_SI(r,z) !partial derivative of the component of magnetic unit vector
!!$  use precision,only:p_
!!$  implicit none
!!$  real(p_):: unitbphi_z_SI,r,z
!!$  real(p_):: b_SI,bphi_z_SI,b_z_SI,bphi_SI
!!$  unitbphi_z_SI=(bphi_z_SI(r,z)*b_SI(r,z)-b_z_SI(r,z)*bphi_SI(r,z))/b_SI(r,z)**2
!!$end function 
!!$
!!$
