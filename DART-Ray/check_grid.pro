pro check_grid
compile_opt idl2

gridname = './GALAXY_GRIDS/M33/grid_model_M33_all_l0.220um.h5'

grid = h5f_open(gridname)
dens_data = h5d_open(grid, 'dens')
dens_stars_data = h5d_open(grid, 'dens_stars')

dens = h5d_read(dens_data)
dens_stars = h5d_read(dens_stars_data)

h5d_close, dens_data
h5d_close, dens_stars_data
h5f_close, grid

help, dens
print, 'dens = ', total(dens)
print, 'dens_stars = ', total(dens_stars)
stop

end

