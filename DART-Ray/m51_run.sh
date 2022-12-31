#!/bin/tcsh
setenv OMP_NUM_THREADS 2

set model="wd01_abs"
echo "$model"

if ( -f "./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr1.in") then
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr1.in
endif
if ( -f "./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr2.in") then
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr2.in
endif
if ( -f "./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr3.in") then
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr3.in
endif
if ( -f "./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr4.in") then
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr4.in
endif
if ( -f "./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr5.in") then
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr5.in
endif
if ( -f "./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr6.in") then
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr6.in
endif
if ( -f "./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr_HII.in") then
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr_HII.in
endif
if ( -f "./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr_HII4.in") then
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr_HII4.in
endif
if ( -f "./2DTO3D_GRIDS/$model/input_2dto3d_m51_tot.in") then
	mpirun -n 2 dartray_galaxy ./2DTO3D_GRIDS/$model/input_2dto3d_m51_tot.in
endif
