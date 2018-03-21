# YouJunHu's makefile
program_name=Lorentz_ions_model
#main_src=main_DK_electrons.f90
main_src=main_adiabatic_electrons.f90
#main_src=test1.f90
ifeq  (${NERSC_HOST},edison)
  COMPILER=	ftn
   #OPTION = -FR -r8 -heap-arrays -O2 -g -traceback -check bounds
OPTION = -FR -r8 -heap-arrays -O2 -g -traceback 
else ifeq (${NERSC_HOST},cori)
  COMPILER=	ftn
  #OPTION = -FR -r8 -heap-arrays -O2 -g -traceback -check bounds -check all
OPTION = -FR -r8 -heap-arrays -O2 -g -traceback -qopenmp
else
  COMPILER=	mpif90
  #OPTION= -fbounds-check -fopenmp -fimplicit-none
  #OPTION= -Wall -Wextra -Wconversion -pedantic  -ffpe-trap=zero,overflow,underflow -fbounds-check 
  #OPTION=-Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace
  OPTION= -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -fimplicit-none -fbounds-check 
  lapack_location:=/usr/lib/lapack/liblapack.so.3
  blas_location:=/usr/lib/libblas/libblas.so.3
endif

fftw_include=/usr/include

COMPILE=	$(COMPILER) $(OPTION) -c 

f90sources= modules.f90 interpolate.f90 mode_structure.f90 pputil_yj.f90  construct_numerical_tokamak_equilibrium.f90 read_gfile.f90 read_parameters.f90 magnetic_field.f90 normalized_magnetic_field.f90  contour.f90 \
 draw_magnetic_surface.f90 calculate_metric.f90  cal_magnetic_coordinates.f90 construct_poloidal_coordinate.f90 mapping_cylindrical_to_magnetic_coordinates.f90  \
connection_condition_at_theta_cut.f90 \
cal_toroidal_shift.f90 field_lines_analyse.f90 arrange_lcfs.f90 arc_length.f90 interpolate_from_cylindrical_to_magnetic_coordinates.f90 boris.f90 backward_half_step_for_boris.f90 forward_half_step_for_boris.f90 test_full_orbit.f90 push_ion_orbit.f90 sort_ions.f90 deposit_ions.f90  deposit_electrons.f90 push_electron_orbit_weight.f90 push_gc_in_mc.f90 push_electron_cylindrical.f90 fluid_electron_pressure.f90 load_ions.f90  load_electrons.f90  initial_weight.f90 derivatives_in_field_line_following_coordinates.f90 filter.f90 transform.f90 field_solving.f90 cal_array_in_magnetic_coordinates_for_guiding_center_pusher.f90 \
magnetic_coordinates_derivatives_in_cylindrical.f90 \
push_ion_weight.f90 ion_velocity_components_in_mc.f90 field_perturbation_on_marker.f90 prepare_coefficients_of_field_equations.f90 prepare_source_terms.f90 some_tests.f90 set_computational_radial_region.f90 spline.f90 spline2.f90 allocate_field_matrix.f90 evolve_magnetic_field.f90 communicate_field_value_between_neighbour_cells.f90  \
get_nearby_field_along_field_line.f90    smoothing_along_field_line.f90 adiabatic_electron_field_solver.f90 \
create_toroidal_grids.f90  math.f90  ion_theta_and_alpha.f90 ion_marker_radial_bdry_condition.f90 diagnosis.f90 \
report.f90 save_data_for_restarting.f90   $(main_src)
f77sources=  pnpoly.for  

f90objs= $(f90sources:.f90=.o)
f77objs= $(f77sources:.for=.o)
main_obj=$(main_src:.f90=.o)
$(program_name): $(f90objs) $(f77objs)
	$(COMPILER)  $(OPTION) $(f90objs) $(f77objs)  $(lapack_location) $(blas_location) -lfftw3 -lm -o $(program_name) 

modules.o: modules.f90
	$(COMPILE) $< -o $@
pputil_yj.o: pputil_yj.f90
	$(COMPILE) $< -o $@
connection_condition_at_theta_cut.o: connection_condition_at_theta_cut.f90 modules.f90
	 $(COMPILE) $< -o $@
mode_structure.o:mode_structure.f90 modules.f90
	 $(COMPILE) $< -o $@
$(main_obj): ${main_src} modules.f90
	$(COMPILE) $< -o $@
construct_numerical_tokamak_equilibrium.o: construct_numerical_tokamak_equilibrium.f90 modules.f90
	$(COMPILE) $< -o $@
read_gfile.o: read_gfile.f90 modules.f90
	 $(COMPILE) $< -o $@
magnetic_field.o: magnetic_field.f90 modules.f90
	 $(COMPILE) $< -o $@
fluid_electron_pressure.o: fluid_electron_pressure.f90 modules.f90
	 $(COMPILE) $< -o $@
construct_poloidal_coordinate.o: construct_poloidal_coordinate.f90 modules.f90
	 $(COMPILE) $< -o $@
cal_toroidal_shift.o: cal_toroidal_shift.f90 modules.f90
	 $(COMPILE) $< -o $@
read_parameters.o:read_parameters.f90 modules.f90
	 $(COMPILE) $< -o $@
push_ion_orbit.o: push_ion_orbit.f90 modules.f90
	 $(COMPILE) $< -o $@
boris.o: boris.f90 modules.f90
	 $(COMPILE) $< -o $@
backward_half_step_for_boris.o: backward_half_step_for_boris.f90 modules.f90
	 $(COMPILE) $< -o $@
forward_half_step_for_boris.o: forward_half_step_for_boris.f90 modules.f90
	 $(COMPILE) $< -o $@
test_full_orbit.o: test_full_orbit.f90 modules.f90
	 $(COMPILE) $< -o $@
ion_velocity_components_in_mc.o: ion_velocity_components_in_mc.f90 modules.f90
	$(COMPILE) $< -o $@
push_ion_weight.o: push_ion_weight.f90 modules.f90
	 $(COMPILE) $< -o $@
ion_theta_and_alpha.o:ion_theta_and_alpha.f90 modules.f90
	 $(COMPILE) $< -o $@
ion_marker_radial_bdry_condition.o:ion_marker_radial_bdry_condition.f90 modules.f90
	 $(COMPILE) $< -o $@
sort_ions.o: sort_ions.f90 modules.f90
	 $(COMPILE) $< -o $@
push_electron_orbit_weight.o: push_electron_orbit_weight.f90 modules.f90
	 $(COMPILE) $< -o $@
field_perturbation_on_marker.o: field_perturbation_on_marker.f90 modules.f90
	 $(COMPILE) $< -o $@
deposit_ions.o: deposit_ions.f90 modules.f90
	 $(COMPILE) $< -o $@
deposit_electrons.o: deposit_electrons.f90 modules.f90
	 $(COMPILE) $< -o $@
push_gc_in_mc.o: push_gc_in_mc.f90 modules.f90
	 $(COMPILE) $< -o $@
push_electron_cylindrical.o: push_electron_cylindrical.f90 modules.f90
	 $(COMPILE) $< -o $@
cal_array_in_magnetic_coordinates_for_guiding_center_pusher.o: cal_array_in_magnetic_coordinates_for_guiding_center_pusher.f90 modules.f90
	 $(COMPILE) $< -o $@
magnetic_coordinates_derivatives_in_cylindrical.o: magnetic_coordinates_derivatives_in_cylindrical.f90 modules.f90
	 $(COMPILE) $< -o $@
set_computational_radial_region.o: set_computational_radial_region.f90 modules.f90
	 $(COMPILE) $< -o $@
load_ions.o:load_ions.f90 modules.f90
	 $(COMPILE) $< -o $@
load_electrons.o:load_electrons.f90 modules.f90
	$(COMPILE) $< -o $@
initial_weight.o: initial_weight.f90
	 $(COMPILE) $< -o $@
field_solving.o: field_solving.f90 modules.f90
	 $(COMPILE)  -I$(fftw_include) $< -o $@
transform.o: transform.f90 modules.f90
	$(COMPILE)  -I$(fftw_include) $< -o $@
filter.o:filter.f90 modules.f90
	 $(COMPILE)   $< -o $@

communicate_field_value_between_neighbour_cells.o: communicate_field_value_between_neighbour_cells.f90 modules.f90
	 $(COMPILE) $< -o $@
adiabatic_electron_field_solver.o:adiabatic_electron_field_solver.f90 modules.f90
	 $(COMPILE) $< -o $@
create_toroidal_grids.o: create_toroidal_grids.f90 modules.f90
	 $(COMPILE) $< -o $@
report.o: report.f90 modules.f90
	 $(COMPILE)  -I$(fftw_include) $< -o $@
some_tests.o: some_tests.f90 modules.f90
	 $(COMPILE) $< -o $@
math.o: math.f90 modules.f90
	 $(COMPILE) $< -o $@
diagnosis.o: diagnosis.f90 modules.f90
	$(COMPILE) $< -o $@
prepare_coefficients_of_field_equations.o: prepare_coefficients_of_field_equations.f90 modules.f90
	 $(COMPILE)   $< -o $@
prepare_source_terms.o: prepare_source_terms.f90 modules.f90
	 $(COMPILE)  $< -o $@
normalized_magnetic_field.o: normalized_magnetic_field.f90 modules.f90
	 $(COMPILE) $< -o $@
field_lines_analyse.o: field_lines_analyse.f90 modules.f90
	  $(COMPILE) $< -o $@
contour.o: contour.f90 modules.f90
	  $(COMPILE) $< -o $@
cal_magnetic_coordinates.o: cal_magnetic_coordinates.f90 modules.f90
	  $(COMPILE) $< -o $@
mapping_cylindrical_to_magnetic_coordinates.o: mapping_cylindrical_to_magnetic_coordinates.f90 modules.f90
	 $(COMPILE) $< -o $@
interpolate_from_cylindrical_to_magnetic_coordinates.o: interpolate_from_cylindrical_to_magnetic_coordinates.f90 modules.f90
	 $(COMPILE) $< -o $@
interpolate.o: interpolate.f90 modules.f90
	 $(COMPILE) $< -o $@
arrange_lcfs.o: arrange_lcfs.f90 modules.f90
	 $(COMPILE) $< -o $@
arc_length.o: arc_length.f90 modules.f90
	  $(COMPILE) $< -o $@
draw_magnetic_surface.o: draw_magnetic_surface.f90 modules.f90
	  $(COMPILE) $< -o $@
pnpoly.o: pnpoly.for modules.f90
	 $(COMPILE) $< -o $@
spline.o: spline.f90 modules.f90
	 $(COMPILE) $< -o $@
spline2.o: spline2.f90 modules.f90
	  $(COMPILE) $< -o $@
calculate_metric.o: calculate_metric.f90 modules.f90
	   $(COMPILE) $< -o $@
allocate_field_matrix.o:allocate_field_matrix.f90 modules.f90
	 $(COMPILE) $< -o $@
evolve_magnetic_field.o: evolve_magnetic_field.f90 modules.f90
	 $(COMPILE) $< -o $@
derivatives_in_field_line_following_coordinates.o:derivatives_in_field_line_following_coordinates.f90 modules.f90
	 $(COMPILE) $< -o $@
get_nearby_field_along_field_line.o:get_nearby_field_along_field_line.f90 modules.f90
	 $(COMPILE) $< -o $@
smoothing_along_field_line.o:smoothing_along_field_line.f90 modules.f90
	 $(COMPILE) $< -o $@
save_data_for_restarting.o: save_data_for_restarting.f90 modules.f90
	 $(COMPILE) $< -o $@
.PHONY : clean run tarfile
clean :
	 rm -f program_name  $(f90objs) $(f77objs) *.mod  $(program_name)
tarfile:
	mkdir $(program_name)_version`date +"%Y-%m-%d"` && cp -r $(f90sources) $(f77sources) makefile input.nmlt $(program_name)_version`date +"%Y-%m-%d"` && tar -cf $(program_name)_version`date +"%Y-%m-%d"`.tar $(program_name)_version`date +"%Y-%m-%d"` && rm -r $(program_name)_version`date +"%Y-%m-%d"` 
merge_to_one_file:
	cat $(f90sources) >one_file.f90
	cat $(f77sources) >one_file.for
	$(COMPILER)  $(OPTION) one_file.f90 one_file.for  $(lapack_location) $(blas_location) -I$(fftw_include) -lfftw3 -lm -o ./a.out
#sync_to_cluster:
#	 make clean && rsync -avz --delete /home/yj/project_new/lorentz_ions_drift_e_model/ yj@202.127.204.22:/scratch/yj/lorentz_ions_drift_e_model/
sync_to_cori:
	 make clean && rsync -avz --delete /home/yj/project_new/lorentz_ions_drift_e_model/ yjhu@cori.nersc.gov:~/lorentz_ions_drift_e_model/
run: $(program_name)
	date +%H:%M:%S
	mpiexec -n 8 ./$(program_name)
	date +%H:%M:%S
