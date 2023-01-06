# DART-Ray

DART-Ray is a 3D ray-tracing dust radiative transfer code. It can produce the surface brightness maps for a precalculated grid of dust emissivity.

For reference, the following nomenclature is used to denote the different components:
- 'irr1' is the main thick dust disk.
- 'irr2' is the main thin dust disk.
- 'irr3' is the inner thick dust disk.
- 'irr4' is the inner thin dust disk.
- 'irr5' is the outer thick dust disk.
- 'irr6' is the outer thin dust disk.
- 'irr_HII' is the main HII component.
- 'irr_HII4' is the inner HII component.
- 'tot' is the combination of the above components.


0. Run `create_directories.sh` to create directories for later steps.

1. In /2DTO3D_GRIDS/wd01_abs/ there are files named `input_2dto3d_m51_*.in` for each dust component. In these files find the string beginning "label_model_2d= ", this variable should correspond to the model you are running and should be identical to the strings produced in /emission_NUrad/outdata_intlum/. For example, for the file in /emission_NUrad/outdata_intlum/ named: "grid_irr_wd01_q06_t635_s800_no210_bd1_hd8000_zd160_hd1_4700_zd1_90_hs4000_zs190_hs1_4700_zs1_90_reff500_ell175_abs.xdr", label_model_2d for the "tot" component should be "label_model_2d= 'wd01_q06_t635_s800_no210_bd1_hd8000_zd160_hd1_4700_zd1_90_hs4000_zs190_hs1_4700_zs1_90_reff500_ell175_abs_tot',".
The quickest way to change these strings is using the sed command, e.g. `> sed -i 's/old_string/new_string/g' input_grid_2dto3d_m51_*.in`.

1. Run `m51_writegrids.pro`. This reads the files from /emission_NUrad/outdata_intlum and writes files readable by DART-Ray.

2. Run `m51_gridmaster.sh`. This program creates the main and lambda grids containing the spatial coordinates and size of each grid element, along with the corresponding volume emissivity and extinction factor. The former is defined in 2DTO3D_GRIDS/wd01_abs/input_grid_2dto3d_m51_*.in and the latter is determined by the outputs of /emission_NUrad/.

3. Run `m51_run.sh`. Now the grids have been made, `m51_run.sh` executes the radiative transfer calculations according to the contents of 2DTO3D_GRIDS/wd01_abs/input_2dto3d_m51_*.in. This code can be run in parallel and the number of cores used is specified by `setenv OMP_NUM_THREADS N` where N is the number of CPUs. Note that in the command `mpirun -n X...` X must be an integer multiple of N.

4. Run `m51_fits.pro`. The previous step has now produced the data cubes containing the surface brightness maps as viewed at the specified angle, this program now writes these outputs as .fits files (in units of [MJy/sr]).

5. Run `m51_phot.pro`. Now the surface brightness maps have been made, the photometry can be done to produce the azimuthally averaged surface brightness profiles.
