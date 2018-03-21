subroutine fluid_electron_pressure()
  use mpi
  use precision,only:p_
  use constants,only: one,kev
  use normalizing,only: vn_e
  use electrons_module,only: te0,mass_e
  use perturbation_field_matrix,only: den_left
  use perturbation_field_matrix,only: pper_e_left,ppar_e_left !as output
  implicit none

!for isothermal fluid electron model:
  pper_e_left=den_left*(te0*kev)/(mass_e*vn_e**2) !delta_ne=delta_ni is assumed, here den_left is from deposite_ions subroutine
  ppar_e_left=pper_e_left

end subroutine fluid_electron_pressure
