#!/bin/bash
export OMP_NUM_THREADS=8

read -p "Choose your model (i.e., 'wd01_abs'): " model
echo "model: $model"

if test -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr1.in"; then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr1.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr1.in
	qsub ./SCRIPTS/script_dartray_galaxy_openmpi8.sh ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr1.in
fi

if test -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr2.in"; then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr2.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr2.in
	qsub ./SCRIPTS/script_dartray_galaxy_openmpi8.sh ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr2.in
fi

if test -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr3.in"; then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr3.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr3.in
	qsub ./SCRIPTS/script_dartray_galaxy_openmpi8.sh ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr3.in
fi

if test -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr4.in"; then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr4.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr4.in
	qsub ./SCRIPTS/script_dartray_galaxy_openmpi8.sh ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr4.in
fi

if test -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr5.in"; then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr5.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr5.in
	qsub ./SCRIPTS/script_dartray_galaxy_openmpi8.sh ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr5.in
fi

if test -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr6.in"; then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr6.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr6.in
	qsub ./SCRIPTS/script_dartray_galaxy_openmpi8.sh ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr6.in
fi

if test -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_HII.in"; then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_HII.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_HII.in
	qsub ./SCRIPTS/script_dartray_galaxy_openmpi8.sh ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr_HII.in
fi

if test -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_HII4.in"; then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_HII4.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_HII4.in
	qsub ./SCRIPTS/script_dartray_galaxy_openmpi8.sh ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr_HII4.in
fi

if test -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_HII6.in"; then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_HII6.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_HII6.in
	qsub ./SCRIPTS/script_dartray_galaxy_openmpi8.sh ./2DTO3D_GRIDS/$model/input_2dto3d_m51_irr_HII6.in
fi

if test -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_clump1.in"; then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/input_grid_2dto3d_m51_irr_clump1.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/input_grid_2dto3d_m51_irr_clump1.in
	qsub ./SCRIPTS/script_dartray_galaxy_openmpi8.sh ./2DTO3D_GRIDS/input_2dto3d_m51_irr_clump1.in
fi

if test -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_clump3.in"; then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/input_grid_2dto3d_m51_irr_clump3.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/input_grid_2dto3d_m51_irr_clump3.in
	qsub ./SCRIPTS/script_dartray_galaxy_openmpi8.sh ./2DTO3D_GRIDS/input_2dto3d_m51_irr_clump3.in
fi

if test -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_irr_clump5.in"; then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/input_grid_2dto3d_m51_irr_clump5.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/input_grid_2dto3d_m51_irr_clump5.in
	qsub ./SCRIPTS/script_dartray_galaxy_openmpi8.sh ./2DTO3D_GRIDS/input_2dto3d_m51_irr_clump5.in
fi

if test -f "./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_tot.in"; then
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_tot.in
	./create_adap_grid_2dto3d ./2DTO3D_GRIDS/$model/input_grid_2dto3d_m51_tot.in
	qsub ./SCRIPTS/script_dartray_galaxy_openmpi8.sh ./2DTO3D_GRIDS/$model/input_2dto3d_m51_tot.in
fi

watch qstat -u sedmodel
