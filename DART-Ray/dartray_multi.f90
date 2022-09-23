PROGRAM DARTRAY_MULTI
  use smooth_grid_routines
  use user_routines_multi
  use rt_routines_multi
  use io_routines
  use dartray_hub
  IMPLICIT NONE
  
  !! Initialize MPI 
  call initialize_mpi
  
  !! INPUT VARIABLE INITIALIZE  
  call input_initialize

!!! Read input file 
  call read_input_file

!!! Read wavelength grid 
  call read_lambda_list
  
!!! CHECK INPUT
  call check_input 
  
!!! GRID INITIALIZATION 
  call grid_initialize
  
!!! write info file
  call write_file_info
  
!!! RT calculation
  call dartray_switch

!!! Terminate mpi 
  call terminate_mpi
  
CONTAINS
  
END PROGRAM DARTRAY_MULTI
