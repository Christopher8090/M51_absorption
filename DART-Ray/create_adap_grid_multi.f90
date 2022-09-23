!> Creates the main 3D grid for an axysimmetric galaxy model.
PROGRAM create_adap_grid_multi
  use iso_fortran_env
  use smooth_grid_routines
  use user_routines_multi
  use io_routines
  use sed_routines
  IMPLICIT NONE
  integer :: il,i0(1:1),i1(1:1)
  logical :: file_exists
  integer :: ierr 

  !! Initialize MPI 
  call initialize_mpi

  !! INPUT VARIABLE INITIALIZE  
  call input_initialize

  ! read input file 
  call read_input_multi

  ! read lambda grid 
  call read_lambda_list
  call set_lambda_arr_si

  ! check lambda_list 
  call check_lambda_list_multi

  ! set opacity values 
  call prepare_dust_model
  
  print *, 'Are the output filenames OK? (Press RETURN to continue)'
  print *, trim(adjustl(dir_grid))//grid_file, trim(adjustl(dir_grid))//grid_info_file
  print *, trim(adjustl(dir_grid))//'grid_'//trim(adjustl(label_model_lambda_grid))//'_grid_type_wavelength.h5'
   
  read(*,*)
    
  ! check if main grid already exists. 
  
  inquire(file=trim(adjustl(dir_grid))//grid_file, exist=file_exists)

 ! IF not create it !!!

  select case (file_exists)

  case (.FALSE.)
     print *, 'Main grid does not exist. Create it !'
     ! --------------------------------------------------------------------! 
     ! calculate scaling factors of the disks  (this can be changed depending on specific model) or set up input for dust emission grid 
     
     lambda_in = lambda_ref_SI
     call set_hs_disk  
     call set_multi_input

      call print_profiles_stars
      !stop
     
     ! --------------------------------------------------------------------! 
     ! create grid arrays     
     call create_grid_arrays()
     
     ! create extra arrays for galaxy model 
     call create_grid_arrays_multi()
     
     
     call set_base
     
     ! --------------------------------------------------------------------!
     ! start subdivision loop 
     
     call subdivision_loop()
     
  ! ---------------------------------------------------------------------!
     
     ! fix scaling disks and bulge 
     call fix_dens_stars_arrays
     
     ! --------------------------------------------------------------------!
     ! print 3D grid on H5DF file  
     call print_3d_grid_file()
     
     ! print extra dens stars arrays 
     call print_3d_grid_file_multi()
     
     ! print info file 
     call print_info_file 
     
     print *, 'Main grid creation completed!'
     stop
case(.TRUE.)
   print *, 'Main grid does already exist. Read it !'
   call read_main_grid
   
   call create_grid_arrays_multi()
        
end select


! ---------------------------------------------------------------------!
! loop on lambdas to make monochromatic emissivity/ density tables

! rename label_model

label_model= trim(adjustl(label_model_lambda_grid))//'_'//trim(adjustl(grid_type))

! make lambda grids
i0=minloc(abs(lambda_arr-lambda_min)/lambda_min)-1
i1=minloc(abs(lambda_arr-lambda_max)/lambda_max)-1

do il=i0(1), i1(1)
    
   lambda_in=lambda_arr_SI(il)  
   lambda=lambda_in
   print *, 'lambda= ', lambda_arr(il)
   write(label_wave,'(F9.3)') lambda_arr(il) ! this converts number to string 

   ! set filename lambda grid
   grid_file_lambda='grid_'//trim(adjustl(label_model))//'_l'//trim(adjustl(label_wave))//'um.h5'
   print *, trim(adjustl(grid_file_lambda))
     
   ! set old disk hs_disk
   call set_hs_disk  !!! check you are using the right wavelenghts

   ! set star and dust coefficients 
   call set_multi_input

   !! skip wavelengths if no emission there 
   if ((eta_disk0 == 0).and.(eta_tdisk0 == 0).and.&
		(eta_disk03 == 0).and.(eta_tdisk04 == 0).and.&
		(eta_disk05 == 0).and.(eta_tdisk06 == 0).and.(eta_bulge0 == 0)) cycle 
     
   !! print profiles 
   !call print_profiles_stars
   !stop
  
   ! create grids
   call make_lambda_grid    

   ! check on total luminosities 
   call fix_dens_stars_arrays
     
   ! assign values to parent cells 
   call assign_dens_to_parent
     
   ! print them
   call print_lambda_grid

   !print *, 'grid for lambda = ', lambda_arr(il), ' printed' 

end do

CONTAINS

!> It sets the scale length for the old stellar disk. 
  subroutine set_hs_disk
    integer :: i 
    
    do i = 0, lnum

       if (abs(lambda_in-lambda_arr_SI(i))/lambda_arr_SI(i) < 1E-4) exit 

    end do

    if (i > lnum -1) then 
       print *, 'ERROR(set_hs_disk): lambda_in not found in lambda_arr_SI array!'
       stop
    endif

    hs_disk = hs_disk_arr(i)
	 hs_disk3 = hs_disk3_arr(i)
	 hs_disk5 = hs_disk5_arr(i)
    
 end subroutine set_hs_disk

!> It calculates lambda grids at the selected wavelength.  
   subroutine make_lambda_grid
    
    real(kind=real64) :: x,y,z,cellsize,av_rho_dust,av_rho_disk, av_rho_tdisk, av_rho_bulge
    integer :: cc,i
    real(kind=real64) :: av_rho_dust_disk, av_rho_dust_tdisk,av_rho_dustem
	 real(kind=real64) :: av_rho_disk3, av_rho_tdisk4, av_rho_dust_disk3, av_rho_dust_tdisk4
	 real(kind=real64) :: av_rho_disk5, av_rho_tdisk6, av_rho_dust_disk5, av_rho_dust_tdisk6,av_rho_tdisk8

    grid_creation_lambda = .TRUE. 

    !print *, 'calculating grid for lambda =', lambda 
    ! initialize arrays 
    dens = 0            
    dens_disk = 0
    dens_tdisk = 0
    dens_disk3 = 0
    dens_tdisk4 = 0
    dens_disk5 = 0
    dens_tdisk6 = 0
    dens_tdisk8 = 0
    dens_bulge = 0

    call OMP_SET_NESTED(.true.)
    call omp_set_num_threads(nproc)
    
    !$OMP PARALLEL DEFAULT(NONE), PRIVATE(i,x,y,z,cellsize,av_rho_dust,av_rho_disk, av_rho_tdisk,av_rho_disk3, av_rho_tdisk4,av_rho_disk5, av_rho_tdisk6,av_rho_tdisk8, av_rho_bulge,av_rho_dust_disk,av_rho_dust_tdisk,av_rho_dust_disk3,av_rho_dust_tdisk4,av_rho_dust_disk5,av_rho_dust_tdisk6,av_rho_dustem), &
!$OMP SHARED(tot_ncell,csize,dens,dens_disk,dens_tdisk,dens_disk3,dens_tdisk4,dens_disk5,dens_tdisk6,dens_tdisk8,dens_bulge,ccoord,cchild,max_rad,max_z,dens_stars)

    !$OMP DO SCHEDULE(DYNAMIC,50)

    do i=0, tot_ncell-1
             
       if (cchild(i) /= -1 ) cycle 
       !print *, 'i=',i,'of ',tot_ncell-1
       x=ccoord(1,i)
       y=ccoord(2,i)
       z=ccoord(3,i) 
       cellsize=csize(i)

       
       !call calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_disk,av_rho_tdisk, av_rho_bulge,av_rho_dust_disk, av_rho_dust_tdisk)
   	 call calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_disk, av_rho_tdisk, av_rho_disk3, av_rho_tdisk4, av_rho_disk5, av_rho_tdisk6,av_rho_tdisk8, &
			av_rho_bulge, av_rho_dust_disk,av_rho_dust_tdisk, av_rho_dust_disk3,av_rho_dust_tdisk4, av_rho_dust_disk5,av_rho_dust_tdisk6)

       
!!$       if ((sqrt(x**2+y**2+z**2) < max_rad).and.(abs(z) < max_z)) then 
!!$       print *, av_rho_dust, av_rho_disk, av_rho_tdisk, av_rho_bulge
!!$       read(*,*)
!!$       endif 

     
       dens(i)=av_rho_dust                
       !dens_stars(tot_ncell)=av_rho_stars
       dens_disk(i) = av_rho_disk
       dens_tdisk(i) = av_rho_tdisk
!---
       dens_disk3(i) = av_rho_disk3
       dens_tdisk4(i) = av_rho_tdisk4
       dens_disk5(i) = av_rho_disk5
       dens_tdisk6(i) = av_rho_tdisk6
       dens_tdisk8(i) = av_rho_tdisk8
!---
       dens_bulge(i) = av_rho_bulge
       
    enddo

    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

end subroutine make_lambda_grid


!> This subroutine obtains the average dust extinction coefficient and the average stellar luminosity of the current cell from the user-defined routines.
!  subroutine calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_disk, av_rho_tdisk, av_rho_bulge, av_rho_dust_disk,av_rho_dust_tdisk)
!  
!    real(Kind=real64) :: x,y,z,cellsize,av_rho_dust,av_rho_disk, av_rho_tdisk, av_rho_bulge,rad,sigtruncate,truncate
!    real(kind=real64) :: av_rho_dust_disk,av_rho_dust_tdisk,sigtruncate1,truncate1
!  character(LEN=lcar_type) :: disk_comp1, disk_comp2
  
  !Thick disk truncation
!  sigtruncate=sha*rtrun
!  rad=sqrt(x*x+y*y)      
!  truncate=0.5*(1.0-derf((rad-rtrun)/sigtruncate))   

  !thin disk truncation
!  sigtruncate1=sha1*rtrun
!  truncate1=0.5*(1.0-derf((rad-rtrun)/sigtruncate1))

!  disk_comp1='dust_disk'
!  disk_comp2='dust_tdisk'

!  av_rho_dust_disk=av_disk(x,y,z,cellsize,disk_comp1,thick_disk_type_ID)*truncate  
!  av_rho_dust_tdisk=av_disk(x,y,z,cellsize,disk_comp2,thin_disk_type_ID)*truncate1
!  av_rho_dust=av_rho_dust_disk+av_rho_dust_tdisk
  
!  disk_comp1='disk'
!  disk_comp2='tdisk'

!  av_rho_disk=av_disk(x,y,z,cellsize,disk_comp1,old_disk_type_ID)*truncate
!  av_rho_tdisk=av_disk(x,y,z,cellsize,disk_comp2,young_disk_type_ID)*truncate1
!  av_rho_bulge= av_star_bulge(x,y,z,cellsize)*truncate

!end subroutine calc_dens

!> This subroutine obtains the average dust extinction coefficient and the average stellar luminosity of the current cell from the user-defined routines.
  subroutine calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_disk, av_rho_tdisk, av_rho_disk3, av_rho_tdisk4, av_rho_disk5, av_rho_tdisk6, av_rho_tdisk8, &
		av_rho_bulge, av_rho_dust_disk,av_rho_dust_tdisk, av_rho_dust_disk3,av_rho_dust_tdisk4, av_rho_dust_disk5,av_rho_dust_tdisk6)
  
    real(Kind=real64) :: x,y,z,cellsize,av_rho_dust,av_rho_disk, av_rho_tdisk, av_rho_bulge,rad,sigtruncate,truncate,sigtruncate_dust,truncate_dust
    real(kind=real64) :: av_rho_dust_disk,av_rho_dust_tdisk,sigtruncate1,truncate1,sigtruncate_tdust1,truncate_tdust1
!---
    real(Kind=real64) ::av_rho_disk3, av_rho_dust_disk3, sigtruncate3,truncate3, sigtruncate_dust3,truncate_dust3
    real(kind=real64) :: av_rho_tdisk4, av_rho_dust_tdisk4,sigtruncate4,truncate4, sigtruncate_tdust4,truncate_tdust4

    real(Kind=real64) ::av_rho_disk5, av_rho_dust_disk5, sigtruncate5,truncate5, sigtruncate_dust5,truncate_dust5
    real(kind=real64) :: av_rho_tdisk6, av_rho_dust_tdisk6,sigtruncate6,truncate6, sigtruncate_tdust6,truncate_tdust6
    real(Kind=real64) :: av_rho_tdisk8,sigtruncate8,truncate8
!---
    character(LEN=lcar_type) :: disk_comp1, disk_comp2,disk_comp3, disk_comp4,disk_comp5, disk_comp6, disk_comp8
  
  !Old disk truncation
  sigtruncate=sha*rtrun_disk
  rad=sqrt(x*x+y*y)      
  truncate=0.5*(1.0-derf((rad-rtrun_disk)/sigtruncate))   

  !Young disk truncation
  sigtruncate1=sha1*rtrun_tdisk
  truncate1=0.5*(1.0-derf((rad-rtrun_tdisk)/sigtruncate1))

  !Thick disk truncation
  sigtruncate_dust=sha*rtrun_dust_disk
  truncate_dust=0.5*(1.0-derf((rad-rtrun_dust_disk)/sigtruncate_dust))   

  !Thin disk truncation
  sigtruncate_tdust1=sha1*rtrun_dust_tdisk
  truncate_tdust1=0.5*(1.0-derf((rad-rtrun_dust_tdisk)/sigtruncate_tdust1))

!---
  !Old disk3 truncation
  sigtruncate3=sha*rtrun_disk3   
  truncate3=0.5*(1.0-derf((rad-rtrun_disk3)/sigtruncate3))   

  !Young disk4 truncation
  sigtruncate4=sha1*rtrun_tdisk4
  truncate4=0.5*(1.0-derf((rad-rtrun_tdisk4)/sigtruncate4))

  !Thick disk3 truncation
  sigtruncate_dust3=sha*rtrun_dust_disk3
  truncate_dust3=0.5*(1.0-derf((rad-rtrun_dust_disk3)/sigtruncate_dust3))   

  !Thin disk4 truncation
  sigtruncate_tdust4=sha1*rtrun_dust_tdisk4
  truncate_tdust4=0.5*(1.0-derf((rad-rtrun_dust_tdisk4)/sigtruncate_tdust4))
!--
  !Old disk5 truncation
  sigtruncate5=sha*rtrun_disk5   
  truncate5=0.5*(1.0-derf((rad-rtrun_disk5)/sigtruncate5))   

  !Young disk6 truncation
  sigtruncate6=sha1*rtrun_tdisk6
  truncate6=0.5*(1.0-derf((rad-rtrun_tdisk6)/sigtruncate6))

  !Thick disk5 truncation
  sigtruncate_dust5=sha*rtrun_dust_disk5
  truncate_dust5=0.5*(1.0-derf((rad-rtrun_dust_disk5)/sigtruncate_dust5))   

  !Thin disk6 truncation
  sigtruncate_tdust6=sha1*rtrun_dust_tdisk6
  truncate_tdust6=0.5*(1.0-derf((rad-rtrun_dust_tdisk6)/sigtruncate_tdust6))
!---
  !Young disk8 truncation
  sigtruncate8=sha1*rtrun_tdisk8
  truncate8=0.5*(1.0-derf((rad-rtrun_tdisk8)/sigtruncate8))
!---


  disk_comp1='dust_disk'
  disk_comp2='dust_tdisk'
  disk_comp3='dust_disk3'
  disk_comp4='dust_tdisk4'
  disk_comp5='dust_disk5'
  disk_comp6='dust_tdisk6'

  av_rho_dust_disk=av_disk(x,y,z,cellsize,disk_comp1,thick_disk_type_ID)*truncate_dust  
  av_rho_dust_tdisk=av_disk(x,y,z,cellsize,disk_comp2,thin_disk_type_ID)*truncate_tdust1

!---
  av_rho_dust_disk3=av_disk(x,y,z,cellsize,disk_comp3,thick_disk3_type_ID)*truncate_dust3  
  av_rho_dust_tdisk4=av_disk(x,y,z,cellsize,disk_comp4,thin_disk4_type_ID)*truncate_tdust4

  av_rho_dust_disk5=av_disk(x,y,z,cellsize,disk_comp5,thick_disk5_type_ID)*truncate_dust5  
  av_rho_dust_tdisk6=av_disk(x,y,z,cellsize,disk_comp6,thin_disk6_type_ID)*truncate_tdust6
!---
  av_rho_dust=av_rho_dust_disk+av_rho_dust_tdisk+av_rho_dust_disk3+av_rho_dust_tdisk4+av_rho_dust_disk5+av_rho_dust_tdisk6
  
  disk_comp1='disk'
  disk_comp2='tdisk'
  disk_comp3='disk3'
  disk_comp4='tdisk4'
  disk_comp5='disk5'
  disk_comp6='tdisk6'
  disk_comp8='tdisk8'


  av_rho_disk=av_disk(x,y,z,cellsize,disk_comp1,old_disk_type_ID)*truncate
  av_rho_tdisk=av_disk(x,y,z,cellsize,disk_comp2,young_disk_type_ID)*truncate1

!---
  av_rho_disk3=av_disk(x,y,z,cellsize,disk_comp3,old_disk3_type_ID)*truncate3
  av_rho_tdisk4=av_disk(x,y,z,cellsize,disk_comp4,young_disk4_type_ID)*truncate4

  av_rho_disk5=av_disk(x,y,z,cellsize,disk_comp5,old_disk5_type_ID)*truncate5
  av_rho_tdisk6=av_disk(x,y,z,cellsize,disk_comp6,young_disk6_type_ID)*truncate6

  av_rho_tdisk8=av_disk(x,y,z,cellsize,disk_comp8,young_disk8_type_ID)*truncate8
!---
  av_rho_bulge= av_star_bulge(x,y,z,cellsize)*truncate

end subroutine calc_dens







!!$subroutine calc_dens_dustem(x,y,z,cellsize,av_rho_dustem)
!!$  ! this subroutine obtains the average dust average luminosity of the current cell from the user-defined routines.
!!$    real(Kind=real64) :: x,y,z,cellsize,av_rho_dustem,rad,sigtruncate,truncate
!!$  
!!$  av_rho_dustem=av_dustem(x,y,z,cellsize)
!!$
!!$end subroutine calc_dens_dustem

!> Returns true if the cell with central coordinates x,y,z and size cellsize has to be subdivided. 
logical function subdivision(x,y,z,cellsize) 
real(kind=real64) :: x,y,z,cellsize,tau,lum_stars_disk
integer :: ib=1 

  tau=dens(ncurr)*cellsize   ! tau=K_ext*rho_N*length  here rho_N is in H atoms per volume

  lum_stars_disk=(dens_disk(ncurr)+dens_tdisk(ncurr)+dens_disk3(ncurr)+dens_tdisk4(ncurr)+dens_disk5(ncurr)+dens_tdisk6(ncurr)+dens_tdisk8(ncurr))*cellsize**3  
  
! subdivision for RT calculations

select case(subdivision_criteria)
case('standard') 
   subdivision=(((nlevel < max_lvl).and.((tau > max_dtau).or.& 
(lum_stars_disk > max_dlum*lnu_tot) & 
.or.(nlevel < min_lvl).or.((abs(z) < z_subd_lim).and.(sqrt(x**2+y**2) < R_subd_lim)) )).or.(lvl(ncurr) < nlevel))  ! do not remove term (lvl(ncurr) < nlevel). Needed for smooth cell refinement algorithm. 
case('milky_way')
   subdivision= (((nlevel < max_lvl).and.((tau > max_dtau).or.& 
(lum_stars_disk > max_dlum*lnu_tot) & 
.or.(nlevel < min_lvl).or.((abs(z) < z_subd_lim).and.(sqrt(x**2+y**2) < R_subd_lim)).or.((omega_cell >= omega_max).and.((abs(z)-cellsize/2.) < z_subd_lim2)) ) ).or.(lvl(ncurr) < nlevel))
! note: the term (abs(z)-cellsize/2.) < z_subd_lim2) is correct. No need to use external abs!!!! otherwise big cells with centre on the plane are not subdivided  
case('fixed') 
   subdivision= ((cellsize/2 > res_min).or.(lvl(ncurr) < nlevel))
case default
   print *, 'STOP(subdivision): subdivision_criteria not recognized!'
   print *, 'subdivision_criteria =',subdivision_criteria
   stop
end select
  
end function subdivision

  subroutine print_info_file()
! this subroutine calculates maximum and minimum tau/lum cell and print info file 
integer :: k,i
real(kind=real32) :: max_dtau_eff,max_dlum_eff,mean_dtau,mean_dlum  

! calculate max and mean 
k=0
max_dtau_eff=-999.
max_dlum_eff=-999.
mean_dtau=0
mean_dlum=0

do i=0, tot_ncell-1 
if (cchild(i) == -1) then  
k=k+1
   if (dens(i)*csize(i) > max_dtau_eff) max_dtau_eff=dens(i)*csize(i)
mean_dtau=mean_dtau+dens(i)*csize(i)
  if (dens_stars(i)*csize(i)**3 > max_dlum_eff) max_dlum_eff=dens_stars(i)*csize(i)**3
mean_dlum=mean_dlum+dens_stars(i)*csize(i)**3

endif 

end do 

mean_dtau=mean_dtau/k
mean_dlum=mean_dlum/k

! print grid info on a file
open(42, file=trim(adjustl(dir_grid))//trim(adjustl(grid_info_file)), status='unknown') 

   write(42,*) 'GRID TYPE '
   write(42,*) 'grid_type= ', grid_type
   write(42,*)  
   write(42,*)  'OLD STELLAR DISK'
   write(42,*)  'old_disk_type= ', old_disk_type
   write(42,*)  'old = ', old
   write(42,*)  'hs_disk_b= ', hs_disk_b 
   write(42,*)  'zs_disk= ',zs_disk
   write(42,*)  'zs_disk_r1= ', zs_disk_r1
   write(42,*)  'zs_disk_rsun= ', zs_disk_rsun   
   write(42,*)  'hsin= ',hsin 
   write(42,*)  'hs_disk= ', hs_disk
   write(42,*)	 'intrunc_disk= ', intrunc_disk
	write(42,*)  'rtrun_disk= ', rtrun_disk
   write(42,*)
   write(42,*)   '! YOUNG STELLAR DISK'
   write(42,*)   'young_disk_type= ', young_disk_type
   write(42,*)   'sfr= ', sfr
   write(42,*)   'hs_tdisk= ', hs_tdisk     
   write(42,*)   'zs_tdisk= ', zs_tdisk
   write(42,*)   'zs_tdisk_r1= ', zs_tdisk_r1
   write(42,*)   'zs_tdisk_rsun= ', zs_tdisk_rsun
   write(42,*)   'hs1in= ', hs1in
   write(42,*)	 'intrunc_tdisk= ', intrunc_tdisk
	write(42,*)  'rtrun_tdisk= ', rtrun_tdisk
   write(42,*) 

!---
   write(42,*)  'OLD STELLAR DISK3'
   write(42,*)  'old_disk3_type= ', old_disk3_type
   write(42,*)  'old3 = ', old3
   write(42,*)  'hs_disk3_b= ', hs_disk3_b 
   write(42,*)  'zs_disk3= ',zs_disk3
   write(42,*)  'zs_disk3_r1= ', zs_disk3_r1
   write(42,*)  'zs_disk3_rsun= ', zs_disk3_rsun   
   write(42,*)  'hs3in= ',hs3in 
   write(42,*)  'hs_disk3= ', hs_disk3
   write(42,*)	 'intrunc_disk3= ', intrunc_disk3
	write(42,*)  'rtrun_disk3= ', rtrun_disk3
   write(42,*)
   write(42,*)   '! YOUNG STELLAR DISK4'
   write(42,*)   'young_disk4_type= ', young_disk4_type
   write(42,*)   'sfr4= ', sfr4
   write(42,*)   'hs_tdisk4= ', hs_tdisk4     
   write(42,*)   'zs_tdisk4= ', zs_tdisk4
   write(42,*)   'zs_tdisk4_r1= ', zs_tdisk4_r1
   write(42,*)   'zs_tdisk4_rsun= ', zs_tdisk4_rsun
   write(42,*)   'hs4in= ', hs4in
   write(42,*)	 'intrunc_tdisk4= ', intrunc_tdisk4
	write(42,*)  'rtrun_tdisk4= ', rtrun_tdisk4
   write(42,*) 

   write(42,*)  'OLD STELLAR DISK5'
   write(42,*)  'old_disk5_type= ', old_disk5_type
   write(42,*)  'old5 = ', old5
   write(42,*)  'hs_disk5_b= ', hs_disk5_b 
   write(42,*)  'zs_disk5= ',zs_disk5
   write(42,*)  'zs_disk5_r1= ', zs_disk5_r1
   write(42,*)  'zs_disk5_rsun= ', zs_disk5_rsun   
   write(42,*)  'hs5in= ',hs5in 
   write(42,*)  'hs_disk5= ', hs_disk5
   write(42,*)	 'intrunc_disk5= ', intrunc_disk5
	write(42,*)  'rtrun_disk5= ', rtrun_disk5
   write(42,*)
   write(42,*)   '! YOUNG STELLAR DISK6'
   write(42,*)   'young_disk6_type= ', young_disk6_type
   write(42,*)   'sfr6= ', sfr6
   write(42,*)   'hs_tdisk6= ', hs_tdisk6     
   write(42,*)   'zs_tdisk6= ', zs_tdisk6
   write(42,*)   'zs_tdisk6_r1= ', zs_tdisk6_r1
   write(42,*)   'zs_tdisk6_rsun= ', zs_tdisk6_rsun
   write(42,*)   'hs6in= ', hs6in
   write(42,*)	 'intrunc_tdisk6= ', intrunc_tdisk6
	write(42,*)  'rtrun_tdisk6= ', rtrun_tdisk6
   write(42,*) 
   write(42,*)   '! YOUNG STELLAR DISK8'
   write(42,*)   'young_disk8_type= ', young_disk8_type
   write(42,*)   'sfr8= ', sfr8
   write(42,*)   'hs_tdisk8= ', hs_tdisk8     
   write(42,*)   'zs_tdisk8= ', zs_tdisk8
   write(42,*)   'zs_tdisk8_r1= ', zs_tdisk8_r1
   write(42,*)   'zs_tdisk8_rsun= ', zs_tdisk8_rsun
   write(42,*)   'hs8in= ', hs8in
   write(42,*)	 'intrunc_tdisk8= ', intrunc_tdisk8
	write(42,*)  'rtrun_tdisk8= ', rtrun_tdisk8
   write(42,*) 
!---
   
   write(42,*)   '! BULGE'
   write(42,*)   'reff= ', reff
   write(42,*)   'acap_bulge= ',acap_bulge 
   write(42,*)   'ellipt= ', ellipt
   write(42,*)   'mtrunc= ', mtrunc
   write(42,*)   'bd_ratio= ', bd_ratio
   write(42,*)   'nsersic= ', nsersic
   write(42,*)
   write(42,*)   '! THICK DUST DISK '
   write(42,*)   'thick_disk_type= ', thick_disk_type
   write(42,*)   'tau1= ', tau1
   write(42,*)   'hd_disk= ', hd_disk
   write(42,*)   'zd_disk= ', zd_disk
   write(42,*)   'zd_disk_r1= ', zd_disk_r1
   write(42,*)   'zd_disk_rsun= ', zd_disk_rsun
   write(42,*)	 'intrunc_dust_disk= ', intrunc_dust_disk
	write(42,*)  'rtrun_dust_disk= ', rtrun_dust_disk
   
   write(42,*)   'hdin= ', hdin
   write(42,*)
   write(42,*)   '! THIN DUST DISK'
   write(42,*)   'thin_disk_type= ', thin_disk_type
   write(42,*)   'tau2= ', tau2    
   write(42,*)   'hd_tdisk= ', hd_tdisk
   write(42,*)   'zd_tdisk= ', zd_tdisk
   write(42,*)   'zd_tdisk_r1= ', zd_tdisk_r1
   write(42,*)   'zd_tdisk_rsun= ', zd_tdisk_rsun
   write(42,*)   'hd1in= ', hd1in
   write(42,*)	 'intrunc_dust_tdisk= ', intrunc_dust_tdisk
	write(42,*)  'rtrun_dust_tdisk= ', rtrun_dust_tdisk

!---
   write(42,*)
   write(42,*)   '! THICK DUST DISK3 '
   write(42,*)   'thick_disk3_type= ', thick_disk3_type
   write(42,*)   'tau3= ', tau3
   write(42,*)   'hd_disk3= ', hd_disk3
   write(42,*)   'zd_disk3= ', zd_disk3
   write(42,*)   'zd_disk3_r1= ', zd_disk3_r1
   write(42,*)   'zd_disk3_rsun= ', zd_disk3_rsun
   write(42,*)   'hd3in= ', hd3in
   write(42,*)	 'intrunc_dust_disk3= ', intrunc_dust_disk3
	write(42,*)  'rtrun_dust_disk3= ', rtrun_dust_disk3
   write(42,*)
   write(42,*)   '! THIN DUST DISK4'
   write(42,*)   'thin_disk4_type= ', thin_disk4_type
   write(42,*)   'tau4= ', tau4    
   write(42,*)   'hd_tdisk4= ', hd_tdisk4
   write(42,*)   'zd_tdisk4= ', zd_tdisk4
   write(42,*)   'zd_tdisk4_r1= ', zd_tdisk4_r1
   write(42,*)   'zd_tdisk4_rsun= ', zd_tdisk4_rsun
   write(42,*)   'hd4in= ', hd4in
   write(42,*)	 'intrunc_dust_tdisk4= ', intrunc_dust_tdisk4
	write(42,*)  'rtrun_dust_tdisk4= ', rtrun_dust_tdisk4

   write(42,*)
   write(42,*)   '! THICK DUST DISK5 '
   write(42,*)   'thick_disk5_type= ', thick_disk5_type
   write(42,*)   'tau5= ', tau5
   write(42,*)   'hd_disk5= ', hd_disk5
   write(42,*)   'zd_disk5= ', zd_disk5
   write(42,*)   'zd_disk5_r1= ', zd_disk5_r1
   write(42,*)   'zd_disk5_rsun= ', zd_disk5_rsun
   write(42,*)   'hd5in= ', hd5in
   write(42,*)	 'intrunc_dust_disk5= ', intrunc_dust_disk5
	write(42,*)  'rtrun_dust_disk5= ', rtrun_dust_disk5
   write(42,*)
   write(42,*)   '! THIN DUST DISK6'
   write(42,*)   'thin_disk6_type= ', thin_disk6_type
   write(42,*)   'tau6= ', tau6    
   write(42,*)   'hd_tdisk6= ', hd_tdisk6
   write(42,*)   'zd_tdisk6= ', zd_tdisk6
   write(42,*)   'zd_tdisk6_r1= ', zd_tdisk6_r1
   write(42,*)   'zd_tdisk6_rsun= ', zd_tdisk6_rsun
   write(42,*)   'hd6in= ', hd6in
   write(42,*)	 'intrunc_dust_tdisk6= ', intrunc_dust_tdisk6
	write(42,*)  'rtrun_dust_tdisk6= ', rtrun_dust_tdisk6

!---
  
   write(42,*)
   write(42,*)   '! Other geometrical parameters' 
   write(42,*)   'rsun= ', rsun
   write(42,*)   'max_z= ', max_z
   write(42,*)   'max_rad= ', max_rad
   write(42,*)   'sha= ', sha
   write(42,*)   'omega_max= ', omega_max
   write(42,*)
   write(42,*)   '! Wavelength '
   write(42,*)   'lambda_in= ', lambda_in
   write(42,*)
   write(42,*)  '! 3D grid parameters '  
   write(42,*)  'modelsize = ', modelsize
   write(42,*) ' base =', base  
   write(42,*) ' max_ncell', max_ncell 
   write(42,*) ' MAX DTAU PARAM (input)=', max_dtau
   write(42,*) ' MAX DLUM PARAM (input)=', max_dlum
   write(42,*) ' MIN NLVL (input)=', min_lvl
   write(42,*) ' MAX NLVL (input)=', max_lvl
   write(42,*) ' MAX DTAU (output)=', max_dtau_eff
   write(42,*) ' MEAN DTAU (output)=', mean_dtau
   write(42,*) ' MAX DLUM (output)=', max_dlum_eff
   write(42,*) ' MEAN DLUM (output)=', mean_dlum
   write(42,*)

close(42)

end subroutine print_info_file


!!!!-----------------------------------!!!!!
!!! DO NOT MODIFY FOLLOWING ROUTINES !!!! 
!!!------------------------------------!!!!   

  subroutine subdivision_loop()
      ! This subroutine performs the subdivision loop until all cells satisfy the input subdivision criteria 
    real(kind=real64) :: x,y,z, cellsize, av_rho_dust,av_rho_disk,av_rho_tdisk, av_rho_bulge
    real(kind=real64) :: lum_stars, tau
    real(kind=real64) :: sun_pos(3),cell_pos(3),rel_vec(3),dist
    integer :: i,outcube(3)
    integer(kind=int32) :: a,b
    real(kind=real64) :: av_rho_dust_disk,av_rho_dust_tdisk,av_rho_dustem
!---
    real(Kind=real64) ::av_rho_disk3, av_rho_dust_disk3
    real(kind=real64) :: av_rho_tdisk4, av_rho_dust_tdisk4

    real(Kind=real64) ::av_rho_disk5, av_rho_dust_disk5
    real(kind=real64) :: av_rho_tdisk6, av_rho_dust_tdisk6
    real(kind=real64) :: av_rho_tdisk8
!---
    
    sun_pos=(/rsun,0._real64,0._real64/)  ! observer position

    ! FIRST CELL initialization
      x=ccoord(1,0)
      y=ccoord(2,0)
      z=ccoord(3,0)

      call calc_cellsize(cellsize,nlevel)

      
!      call calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_disk,av_rho_tdisk, av_rho_bulge,av_rho_dust_disk,av_rho_dust_tdisk)
        call calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_disk, av_rho_tdisk, av_rho_disk3, av_rho_tdisk4, av_rho_disk5, av_rho_tdisk6, av_rho_tdisk8, &
        	av_rho_bulge, av_rho_dust_disk,av_rho_dust_tdisk, av_rho_dust_disk3,av_rho_dust_tdisk4, av_rho_dust_disk5,av_rho_dust_tdisk6)

         dens(0) = av_rho_dust
         !dens_stars(0) = av_rho_stars
         dens_disk(0)= av_rho_disk
         dens_tdisk(0) = av_rho_tdisk
         dens_disk3(0)= av_rho_disk3
         dens_tdisk4(0) = av_rho_tdisk4
         dens_disk5(0)= av_rho_disk5
         dens_tdisk6(0) = av_rho_tdisk6
         dens_tdisk8(0) = av_rho_tdisk8
         dens_bulge(0) = av_rho_bulge
         dens_dust_disk(0)= av_rho_dust_disk
         dens_dust_tdisk(0)=av_rho_dust_tdisk
         dens_dust_disk3(0)= av_rho_dust_disk3
         dens_dust_tdisk4(0)=av_rho_dust_tdisk4
         dens_dust_disk5(0)= av_rho_dust_disk5
         dens_dust_tdisk6(0)=av_rho_dust_tdisk6
         
         
      do 

         subdiv:  do ncurr=nstart(nlevel),nstart(nlevel+1)-1    
! check no parent cell in the list (apart from first step) 
   if ((ncurr /= 0).and.(cchild(ncurr) /= (-1))) then 
      print *, ncurr, cchild(ncurr)
      STOP 'there is something weird going on: subdivision loop ' 
      endif    ! this is needed to check problems with neighbour cell subdivision 

      x=ccoord(1,ncurr)
      y=ccoord(2,ncurr)
      z=ccoord(3,ncurr)

      call calc_cellsize(cellsize,nlevel)

        ! these lines are for creating a good grid from sun position 
      cell_pos=(/x,y,z/)
      rel_vec=sun_pos-cell_pos

      outcube=0
      where(abs(rel_vec) > cellsize/2.)  
      outcube=1
      endwhere 

      if (sum(outcube) > 0) then
        dist=(sum((sun_pos-cell_pos)**2))**0.5
        omega_cell=cellsize/dist ! linear angular size !!!  

         else 
            omega_cell=omega_max/2. 

         endif



      !call calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_stars)
    
      ! dens(ncurr)=av_rho_dust
     
       !dens_stars(ncurr)=av_rho_stars

         
!!$      av_rho_dust=dens(ncurr)
!!$      !av_rho_stars=dens_stars(ncurr)
!!$      av_rho_disk=dens_disk(ncurr)
!!$      av_rho_tdisk=dens_tdisk(ncurr)
!------------------------30 MARCH 20
!!$      av_rho_disk3=dens_disk3(ncurr)
!!$      av_rho_tdisk4=dens_tdisk4(ncurr)
!!$      av_rho_disk5=dens_disk5(ncurr)
!!$      av_rho_tdisk6=dens_tdisk6(ncurr)
!------------------------------------



!!$      av_rho_bulge=dens_bulge(ncurr)
      
      
 if(subdivision(x,y,z,cellsize)) then 

     ! Taking care of neighbour cells if necessary 

    if (lvl(ncurr) > 1) then

     call subdivide_neighbour_cells()

    endif
!!$
!!$! Subdivide current cell  
     
   if (lvl(ncurr) == nlevel) then  
 call subdivide_cell(ncurr, nlevel) 
   endif
!!$
  endif 

 end do subdiv

! update current level
      nlevel=nlevel+1
      nstart(nlevel+1)=tot_ncell+1
 
      if (nlevel < max_lvl) then 
      print *,'level up to',nlevel
      print *,(nstart(i),i=0,nlevel)
      
        cycle
      else          
         exit
      endif 

   end do

! ---------------------------------------------
! Check if there are subdivided neighbours in the list and check if they have neighbours which need to be subdivided
! -------------------------------------------------

   a=nstart(nlevel)
   b=nstart(nlevel+1)

print *, 'A   B'
print *, a,b
!read(*,*)

do

   do ncurr=a,b-1  
      if (lvl(ncurr) < nlevel) then 
         if (lvl(ncurr) > 1) then

            call subdivide_neighbour_cells
               
         endif
      endif
   end do

   a=b
   b=tot_ncell+1
   print *, 'A   B'
   print *, a,b
   ! read(*,*)

   if ((b-a) > 1) then      
      cycle
   else          
      exit
   endif

end do

tot_ncell=tot_ncell+1  ! this is important for io_routines

print *, 'grid subdivision completed'
   

 end subroutine subdivision_loop



subroutine subdivide_neighbour_cells
! this routine subdivide neighbour cells to cell ncurr if necessary (that is if the neighbours of ncurr do not have same nesting level of ncurr) 
  integer :: isel, inc, out, flag_jump
  integer(kind=int32) :: cc,nlevel_cc 
  
 dir_loop: do isel=1,3
  
 inc_loop: do inc=+1,-1,-2    
    
   call find_neighbours(isel,inc,cc,out)

   if (out /= 1) then

      call check_level_jump(cc,flag_jump)
      if (flag_jump == 1) then 
          
         nlevel_cc=lvl(cc)
       
         call subdivide_cell(cc,nlevel_cc) 
 
      endif

   endif

   end do inc_loop
   end do dir_loop


  end subroutine subdivide_neighbour_cells

 subroutine subdivide_cell(cc,nlevel_cc) 
!This subroutine splits a cell in many child cells and creates appropriate links 
IMPLICIT NONE 

integer, pointer  :: ix,iy,iz
integer, target :: ii(3)
integer :: i,j,ib
integer(kind=int32) :: cc,nlevel_cc 
real(Kind=real64) :: cellsize,x,y,z,av_rho_dust,av_rho_disk, av_rho_tdisk, av_rho_bulge
real(kind=real64) :: av_rho_dust_disk,av_rho_dust_tdisk,av_rho_dustem
!---
    real(Kind=real64) ::av_rho_disk3, av_rho_dust_disk3
    real(kind=real64) :: av_rho_tdisk4, av_rho_dust_tdisk4

    real(Kind=real64) ::av_rho_disk5, av_rho_dust_disk5
    real(kind=real64) :: av_rho_tdisk6, av_rho_dust_tdisk6
    real(kind=real64) :: av_rho_tdisk8
!---

ix => ii(1)
iy => ii(2)
iz => ii(3)

   cchild(cc)=tot_ncell+1

   ib=1
   if (nlevel_cc > 0 ) ib=2

if (tot_ncell+base(ib)**3 > max_ncell) then
   print *, 'too many cells! Raise max_ncell and try again' 
   stop
endif

   
   do iz=0,base(ib)-1
   do iy=0,base(ib)-1
   do ix=0,base(ib)-1

     tot_ncell=tot_ncell+1
     ncell(tot_ncell)=tot_ncell
    
     !print *, cindex(cc)

     if (nlevel_cc >0) then
   cindex(tot_ncell)=ior(cindex(cc),((iz*base(ib)+iy)*base(ib)+ix+1)*basediv(1)*basediv(2)**(nlevel_cc-1))
   else 
     cindex(tot_ncell)=ior(cindex(cc),((iz*base(ib)+iy)*base(ib)+ix+1))
     endif

   if (cindex(tot_ncell) < 0) then 
      print *, 'CINDEX problem'
      print *, '< 0'
      print *, cindex(cc)
      stop
      endif 

      call calc_cellsize(cellsize, nlevel_cc+1)
    !cellsize=modelsize*dble(baseinv(1)*baseinv(2)**(nlevel_cc))
     
    csize(tot_ncell)=cellsize
    lvl(tot_ncell)=(nlevel_cc+1)
    

    do j=1,3
       if (mod(base(ib),2) == 1) then 
    ccoord(j,tot_ncell)=ccoord(j,cc)+cellsize*dble(ii(j)-base(ib)/2)
        
       elseif (mod(base(ib),2) == 0) then
    ccoord(j,tot_ncell)=ccoord(j,cc)+cellsize*dble((ii(j)-base(ib)/2+0.5))
         
       else
          STOP 'something wrong here: subdivide routine'
       endif
          
 end do

      x=ccoord(1,tot_ncell)
      y=ccoord(2,tot_ncell)
      z=ccoord(3,tot_ncell)

      
         
         !call calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_disk, av_rho_tdisk, av_rho_bulge, av_rho_dust_disk, av_rho_dust_tdisk)
        call calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_disk, av_rho_tdisk, av_rho_disk3, av_rho_tdisk4, av_rho_disk5, av_rho_tdisk6,av_rho_tdisk8, &
		av_rho_bulge, av_rho_dust_disk,av_rho_dust_tdisk, av_rho_dust_disk3,av_rho_dust_tdisk4, av_rho_dust_disk5,av_rho_dust_tdisk6)
       
         dens(tot_ncell)=av_rho_dust    
         !dens_stars(tot_ncell)=av_rho_stars
         dens_disk(tot_ncell) = av_rho_disk
         dens_tdisk(tot_ncell) = av_rho_tdisk
         dens_disk3(tot_ncell) = av_rho_disk3
         dens_tdisk4(tot_ncell) = av_rho_tdisk4
         dens_disk5(tot_ncell) = av_rho_disk5
         dens_tdisk6(tot_ncell) = av_rho_tdisk6
         dens_tdisk8(tot_ncell) = av_rho_tdisk8
         dens_bulge(tot_ncell) = av_rho_bulge
         dens_dust_disk(tot_ncell)= av_rho_dust_disk
         dens_dust_tdisk(tot_ncell) = av_rho_dust_tdisk
         dens_dust_disk3(tot_ncell)= av_rho_dust_disk3
         dens_dust_tdisk4(tot_ncell) = av_rho_dust_tdisk4
         dens_dust_disk5(tot_ncell)= av_rho_dust_disk5
         dens_dust_tdisk6(tot_ncell) = av_rho_dust_tdisk6
         
     cchild(tot_ncell)=-1

print *, tot_ncell, (ccoord(j,tot_ncell),j=1,3) !, dens(tot_ncell),dens_stars(tot_ncell)


end do
end do
end do

end subroutine subdivide_cell  
  



END PROGRAM create_adap_grid_multi
