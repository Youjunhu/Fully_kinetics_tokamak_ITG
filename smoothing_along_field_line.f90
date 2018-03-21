module smoothing_module
private
public smoothing_along_field_line, smoothing_along_field_line_for_adiabatic_electron_model
contains
subroutine smoothing_along_field_line_core(a,m,n) !smoothing along theta with psi and alpha fixed, i.e., along the magnetic field line
  use precision,only:p_
  use constants,only:one,two,six
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(inout):: a(m,n)
  real(p_):: a_left(m,n),a_right(m,n)
  real(p_),parameter:: weight1=0.5_p_,weight2=-one/six !got to know these values from Dr. Yang Chen (GEM code)

  call get_nearby_field_along_field_line(a,a_left,a_right,m,n) !get value of field on the two grids that are to the left/rightof the present grid

  a=(weight1*a_right+a+weight1*a_left)/(one+two*weight1) !smoothing using weight1, get to know this smoothing scheme from Dr. Yang Chen

  call get_nearby_field_along_field_line(a,a_left,a_right,m,n) !get value of field on the two grids that are to the left/right of the present grid

  a=(weight2*a_right+a+weight2*a_left)/(one+two*weight2) !smoothing using weight2

end subroutine smoothing_along_field_line_core


subroutine smoothing_along_field_line()
  use precision,only:p_
  use perturbation_field_matrix,only: epar_left,ex_left,ey_left, mf_par_left,mf_x_left,mf_y_left
  use magnetic_coordinates,only: m=>mtoroidal,n=>nflux2 
  implicit none
  call smoothing_along_field_line_core(epar_left,m+1,n) 
  call smoothing_along_field_line_core(ex_left,m+1,n) 
  call smoothing_along_field_line_core(ey_left,m+1,n) 
  call smoothing_along_field_line_core(mf_par_left,m+1,n) 
  call smoothing_along_field_line_core(mf_x_left,m+1,n) 
  call smoothing_along_field_line_core(mf_y_left,m+1,n) 

end subroutine smoothing_along_field_line


subroutine smoothing_along_field_line_for_adiabatic_electron_model()
  use precision,only:p_
  use perturbation_field_matrix,only: ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left
  use magnetic_coordinates,only: m=>mtoroidal,n=>nflux2 
  implicit none
  call smoothing_along_field_line_core(ef_cyl_r_left,m+1,n) 
  call smoothing_along_field_line_core(ef_cyl_z_left,m+1,n) 
  call smoothing_along_field_line_core(ef_cyl_phi_left,m+1,n)
  call communicate_field_value_between_neighbour_cells2() !for adiabatic electrons model
end subroutine smoothing_along_field_line_for_adiabatic_electron_model

end module smoothing_module
