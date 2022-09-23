#!/bin/tcsh
setenv OMP_NUM_THREADS 2

echo -n "Please enter model (i.e., 'wd01_abs'): "
set model=$<
echo "$model"

if ( -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr1.in") then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr1.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr1.in
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr1.in
endif
if ( -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr2.in") then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr2.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr2.in
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr2.in
endif
if ( -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr3.in") then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr3.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr3.in
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr3.in
endif
if ( -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr4.in") then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr4.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr4.in
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr4.in
endif
if ( -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr5.in") then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr5.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr5.in
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr5.in
endif
if ( -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr6.in") then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr6.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr6.in
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr6.in
endif
if ( -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_HII.in") then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_HII.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_HII.in
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr_HII.in
endif
if ( -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_HII4.in") then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_HII4.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_HII4.in
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr_HII4.in
endif
if ( -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_tot.in") then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_tot.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_tot.in
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_tot.in
endif
