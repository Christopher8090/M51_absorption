# DART-Ray

DART-Ray is a 3D ray-tracing dust radiative transfer code. It can produce the surface brightness maps for a precalculated grid of dust emissivity.

0. Run `create_directories.sh` to create directories for later steps.

1. Run `m51_writegrids.pro`. This reads the files from /emission_NUrad/outdata_intlum and writes files readable by DART-Ray.

2. Run `m51_gridmaster.sh`. This program creates the main and lambda grids containing the spatial coordinates and size of each grid element, along with the corresponding volume emissivity and extinction factor. The former is defined in 2DTO3D_GRIDS/wd01_abs/input_grid_2dto3d_m51_*.in and the latter is determined by the outputs of /emission_NUrad/.

3. Run `m51_run.sh`. Now the grids have been made, `m51_run.sh` executes the radiative transfer calculations according to the contents of 2DTO3D_GRIDS/wd01_abs/input_2dto3d_m51_*.in. This code can be run in parallel and the number of cores used is specified by `setenv OMP_NUM_THREADS N` where N is the number of CPUs. Note that in the command `mpirun -n X...` X must be an integer multiple of N.

4. Run `m51_fits.pro`. The previous step has now produced the data cubes containing the surface brightness maps as viewed at the specified angle, this program now writes these outputs as .fits files (in units of [MJy/sr]).

5. Run `m51_phot.pro`. Now the surface brightness maps have been made, the photometry can be done to produce the azimuthally averaged surface brightness profiles.
