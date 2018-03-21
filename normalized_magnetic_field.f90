!wrapper functions about the magnetic field
!what the the wrapper do is to change the unit of equilibrium magnetic field function B(R,Z), in which B is given in unit of bn, and R and Z coordinates are given in unit of Ln

function br(r,z) !R component of magnetic field
  use precision,only:p_
  use normalizing,only: Ln,bn
  implicit none
  real(p_):: br,r,z
  real(p_):: br_SI

     br=br_SI(r*Ln,z*Ln)/bn

end function br


function bz(r,z) !Z component of magnetic field
  use precision,only:p_
  use normalizing,only: Ln,bn
  implicit none
  real(p_):: bz,r,z
  real(p_):: bz_SI

  bz=bz_SI(r*Ln,z*Ln)/bn

end function bz

function bphi(r,z) !phi component of magnetic field
  use precision,only:p_
  use normalizing,only: Ln,bn
  implicit none
  real(p_):: bphi,r,z
  real(p_):: bphi_SI
     bphi=bphi_SI(r*Ln,z*Ln)/bn

end function bphi

function b(r,z) !strength of magnetic field
  use precision,only:p_
  use normalizing,only: Ln,bn

  implicit none
  real(p_):: b,r,z
  real(p_):: b_SI

  b=b_SI(r*Ln,z*Ln)/bn

end function b

function b_r(r,z) !partial derivative of magtitude of magnetic field
  use precision,only:p_
  use normalizing,only: Ln,bn
  implicit none
  real(p_):: b_r,r,z
  real(p_):: b_r_SI

  b_r=b_r_SI(r*Ln,z*Ln)/bn*Ln

end function b_r

function b_z(r,z) !partial derivative of magtitude of magnetic field
  use precision,only:p_
  use normalizing,only: Ln,bn
  implicit none
  real(p_):: b_z,r,z
  real(p_):: b_z_SI

  b_z=b_z_SI(r*Ln,z*Ln)/bn*Ln

end function b_z

function b_phi(r,z) result(funcval)
  use precision,only:p_
  use normalizing,only: Ln,bn
  implicit none
  real(p_):: funcval,r,z

  funcval=0._p_

end function


function br_r(r,z) !partial derivative of component of magnetic field
  use precision,only:p_
  use normalizing,only: Ln,bn
  implicit none
  real(p_):: br_r,r,z
  real(p_):: br_r_SI

  br_r=br_r_SI(r*Ln,z*Ln)/bn*Ln

end function br_r


function br_z(r,z) !partial derivative of component of magnetic field
  use precision,only:p_
  use normalizing,only: Ln,bn

  implicit none
  real(p_):: br_z,r,z
  real(p_):: br_z_SI

  br_z=br_z_SI(r*Ln,z*Ln)/bn*Ln
end function br_z


function br_phi(r,z) result(funcval)
  use precision,only:p_
  use normalizing,only: Ln,bn
  implicit none
  real(p_):: funcval,r,z

  funcval=0._p_

end function 


function bz_r(r,z) !partial derivative of component of magnetic field
  use precision,only:p_
  use normalizing,only: Ln,bn
  implicit none
  real(p_):: bz_r,r,z
  real(p_):: bz_r_SI

  bz_r=bz_r_SI(r*Ln,z*Ln)/bn*Ln

end function bz_r


function bz_z(r,z) !partial derivative of component of magnetic field
  use precision,only:p_
  use normalizing,only: Ln,bn
  implicit none
  real(p_):: bz_z,r,z
  real(p_):: bz_z_SI

  bz_z=bz_z_SI(r*Ln,z*Ln)/bn*Ln

end function bz_z

function bz_phi(r,z) result(funcval)
  use precision,only:p_
  use normalizing,only: Ln,bn
  implicit none
  real(p_):: funcval,r,z

  funcval=0._p_

end function 


function bphi_r(r,z) !partial derivative of component of magnetic field
  use precision,only:p_
  use normalizing,only: Ln,bn

  implicit none
  real(p_):: bphi_r,r,z
  real(p_):: bphi_r_SI
  bphi_r= bphi_r_SI(r*Ln,z*Ln)/bn*Ln

end function bphi_r

function bphi_z(r,z) !partial derivative of component of magnetic field
  use precision,only:p_
  use normalizing,only: Ln,bn

  implicit none
  real(p_):: bphi_z,r,z
  real(p_):: bphi_z_SI

  bphi_z=bphi_z_SI(r*Ln,z*Ln)/bn*Ln

end function bphi_z


function unitbr_r(r,z)  result(funcval)!partial derivative of the component of magnetic unit vector
  use precision,only:p_
  use normalizing,only:Ln
  implicit none
  real(p_):: funcval,r,z
  !real(p_):: unitbr_r_SI
  !unitbr_r= unitbr_r_SI(r*Ln,z*Ln)*Ln
  real(p_):: b,bval,br_r,b_r,br
  bval=b(r,z)
  funcval=(br_r(r,z)*bval-b_r(r,z)*br(r,z))/bval**2
end function 

function unitbr_z(r,z) result(funcval)!partial derivative of the component of magnetic unit vector
  use precision,only:p_
  use normalizing,only:Ln
  implicit none
  real(p_):: funcval,r,z
 ! real(p_):: unitbr_z_SI
  !unitbr_z= unitbr_z_SI(r*Ln,z*Ln)*Ln
  real(p_):: b,bval,br_z,b_z,br

  bval=b(r,z)
  funcval=(br_z(r,z)*bval-b_z(r,z)*br(r,z))/bval**2
end function


function unitbr_phi(r,z) result(funcval)!partial derivative of the component of magnetic unit vector
  use precision,only:p_
  use normalizing,only:Ln
  implicit none
  real(p_):: funcval,r,z
  real(p_):: b,bval,br_phi,b_phi,br

  bval=b(r,z)
  funcval=(br_phi(r,z)*bval-b_phi(r,z)*br(r,z))/bval**2
end function 


function unitbz_r(r,z) result(funcval) !partial derivative of the component of magnetic unit vector
  use precision,only:p_
  use normalizing,only:Ln
  implicit none
  real(p_):: funcval,r,z
  !real(p_):: unitbz_r_SI
  !unitbz_r=unitbz_r_SI(r*Ln,z*Ln)*Ln
  real(p_)::b,bval,bz_r,b_r,bz
  bval=b(r,z)
  funcval=(bz_r(r,z)*bval-b_r(r,z)*bz(r,z))/bval**2
end function 

function unitbz_z(r,z) result (funcval)!partial derivative of the component of magnetic unit vector
  use precision,only:p_
!  use normalizing,only:Ln
  implicit none
  real(p_):: funcval,r,z
 ! real(p_):: unitbz_z_SI
  !unitbz_z=unitbz_z_SI(r*Ln,z*Ln)*Ln
real(p_):: b,bval,bz_z,b_z,bz

  bval=b(r,z)
  funcval=(bz_z(r,z)*bval-b_z(r,z)*bz(r,z))/bval**2
end function 

function unitbz_phi(r,z) result (funcval)!partial derivative of the component of magnetic unit vector
  use precision,only:p_
!  use normalizing,only:Ln
  implicit none
  real(p_):: funcval,r,z
 ! real(p_):: unitbz_z_SI
  !unitbz_z=unitbz_z_SI(r*Ln,z*Ln)*Ln
real(p_):: b,bval,bz_phi,b_phi,bz

  bval=b(r,z)
  funcval=(bz_phi(r,z)*bval-b_phi(r,z)*bz(r,z))/bval**2
end function 


function unitbphi_r(r,z) result (funcval)!partial derivative of the component of the unit vector of the magnetic field
  use precision,only:p_
!  use normalizing,only:Ln
  implicit none
  real(p_):: funcval,r,z
  !real(p_):: unitbphi_r_SI
!  unitbphi_r=unitbphi_r_SI(r*Ln,z*Ln)*Ln
real(p_):: b,bval,bphi_r,b_r,bphi

  bval=b(r,z)
  funcval=(bphi_r(r,z)*bval-b_r(r,z)*bphi(r,z))/bval**2
end function 


function unitbphi_z(r,z) result (funcval)!partial derivative of the component of magnetic unit vector
  use precision,only:p_
!  use normalizing,only:Ln
  implicit none
  real(p_):: funcval,r,z
  !real(p_):: unitbphi_z_SI
!  unitbphi_z=unitbphi_z_SI(r*Ln,z*Ln)*Ln

real(p_):: b,bval,bphi_z,b_z,bphi

  bval=b(r,z)
  funcval=(bphi_z(r,z)*bval-b_z(r,z)*bphi(r,z))/bval**2
end function 






