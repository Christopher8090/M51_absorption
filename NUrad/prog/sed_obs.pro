PRO sed_obs
compile_opt idl2

; Calculates the spatially integrated flux densities for M51a up to the radius 'limit_rad'
root = '/nfs/d58/vlk/sedmodel/cinman/m51a/NUrad_M51a/'
save_dir = root+'saves/obs/'

pc_m = 3.0857D16        ; conversion between [pc] and [m].
dist_pc = 8.58D6        ; distance to the galaxy in [pc].
dist_m = dist_pc * pc_m ; distance to the galaxy in [m].
limit_rad = 24.  ; radius in [kpc] out to which the luminosity is integrated
print, "Integrate out to "+strtrim(limit_rad,1)+" kpc."

wavelength = STRARR(20)
wavelength[0]  = '0.1542'       ; GALEX.FUV     ; 
wavelength[1]  = '0.2274'       ; GALEX.NUV     ; 
wavelength[2]  = '0.3562'       ; SDSS.u        ; 
wavelength[3]  = '0.4719'       ; SDSS.g        ; 
wavelength[4]  = '0.6185'       ; SDSS.r        ; 
wavelength[5]  = '0.7500'       ; SDSS.i        ; 
wavelength[6]  = '0.8961'       ; SDSS.z        ; 
wavelength[7]  = '1.200'        ; 2MASS.J       ; 
wavelength[8]  = '1.600'        ; 2MASS.H       ; 
wavelength[9]  = '2.200'        ; 2MASS.K       ; 
wavelength[10] = '3.507'        ; IRAC.I1       ; MJy/sr
wavelength[11] = '4.437'        ; IRAC.I2       ; MJy/sr
wavelength[12] = '5.739'        ; IRAC.I3       ; MJy/sr
wavelength[13] = '7.927'        ; IRAC.I4       ; MJy/sr
wavelength[14] = '24.000'       ; MIPS.24       ; MJy/sr
wavelength[15] = '70.000'       ; PACS.70       ; Jy/pixel
wavelength[16] = '160.000'      ; PACS.160      ; Jy/pixel
wavelength[17] = '250.000'      ; SPIRE.250     ; Jy/beam
wavelength[18] = '350.000'      ; SPIRE.350     ; Jy/beam
wavelength[19] = '500.000'      ; SPIRE.500     ; Jy/beam

llambda = FLOAT(wavelength)	; converts from string to float, units: [um]
int = DBLARR(N_ELEMENTS(wavelength))	; units: [MJy]
int_lam = DBLARR(N_ELEMENTS(wavelength))
int_error = int
print, "wavelength    flux [Jy]     error [Jy]"
FOR i = 0, 19 DO BEGIN
	wave = DOUBLE(wavelength[i]) * 1D-6     ; converts the wavelength from [um] to [m]
	conv_factor = (299792458 / wave^2)  ; conversion used for CpS to Jy [/m /s] or [Hz /m]
	RESTORE, save_dir+'datafile_'+wavelength[i]+'um.save'
	int[i] = MAX(int_rad[WHERE(x_rad LE limit_rad)], /NaN)
	int_lam[i] = int[i]	; MJy
	int_lam[i] = int_lam[i] * 1E6	; Jy
	int_lam[i] = int_lam[i] * 1E-26	; W /m^2 /Hz
	int_lam[i] = int_lam[i] * conv_factor	; W /m^2 /m
	int_lam[i] = int_lam[i] * wave	; W /m^2	; integrated flux in lambda * f_lambda

	int_error[i] = max(cal_error_tot[where(x_rad le limit_rad)])
	PRINT, wavelength[i], int[i]*1e6, int_error[i]*1e6
;continue
ENDFOR
wavelength = llambda
int_jy = int*1E6	; integrated flux in [Jy]
int_jy[16] = int_jy[16] -2.418	; [CII] correction calculated from De Looze 2011
int_error = int_error*1e6	; converts from [MJy] to [Jy]
int_watts = (int_jy * (4*!DPI*dist_m^2)) *1E-26	;[W/Hz]
wave=[912.,1350.,1500.,1650.,2000.,2200.,2500.,2800.,3650.,4430.,5640.,8090.,12590.,22000.,36000.,45000.,58000.]*1e-4
lum_lunit_nu = 2.241d+36
lum_lunit_diff_nu = 1.437d+36 ;in W
;----------------------------------------- INTERPOLATE TO MODEL WAVELENGTHS -------------------------------
int1 = INTERPOL(int_jy, llambda, wave, /SPLINE)
obs_lum_interpolated = [int1, int_jy[13:19]]
wavelength_interpolated = [wave, llambda[13:19]]
;-----------------------------------------------------------------------------------------------------------
save_name = save_dir+'obs_sed.save'
save, wavelength, int_jy, int_watts, int_lam, int_error, filename=save_name
;SAVE, wavelength, int_jy, int_watts, int_lam, wavelength_interpolated, obs_lum_interpolated, int_error, FILENAME=save_name
PRINT, 'WRITTEN: '+save_name
stop

END
