PRO obs_phot
compile_opt idl2
; The purpose of this code is to carry out photometry on the observation data of M51a of the wavelengths listed.
; The input should be the raw maps from NED with a complete mask including masking NGC 5195.
; This code first calculates map size and resolution from the header, and then converts the map units into [MJy/sr].
; Then creates an ellipse centred on NGC 5194, and then calculates the average background flux.
; Then carries out azimuthally averaged radial photometry on the map with background removed.
; Then corrects for the foreground extinction, and plots the radial profile and the curve of growth. CJI 03/12/2020

root = '/nfs/d58/vlk/sedmodel/'
root = root+'cinman/m51a/NUrad_M51a/'
mapdir = root+'maps/obs/'
maskdir = mapdir+'masks/'
plotdir = mapdir+'plotting/'
kerneldir = mapdir+'kernels/'
figdir = root+'figures/obs/'
savedir = root+'saves/obs/'

m_per_pc = 3.0857D+16	; conversion between [pc] and [m].
distance_pc = 8.58D+6	; distance to object in [pc]
c = 299792458d+0
distance_m = distance_pc * m_per_pc	; distance to the galaxy in [m].
flux_to_lum = 4d * !DPI * distance_m^2d * 1d-26	; conversion factor for [Jy] -> [W/hz]
ra = DOUBLE(202.46967)	;right ascension        ; from Turner & Ho (1994)
dec = DOUBLE(47.194989)	;declination    ; from Turner & Ho (1994)
bg_aperture = 40d+0	; ellipticity of the background aperture in degrees
aperture = 20.3d+0	; ellipticity of the photometry aperture in degrees
skip = 0	; 1= skip background calculation, 0= include background calculation
ratio = (COS(aperture * !DPI/180))^(-1)	; ratio between major and minor axis, incl = 20.3 from Hu, 2013
ratio_background = (COS(bg_aperture * !DPI/180))^(-1)	; inclination of 40 degrees used for background calculation
pos_ang = 12.0d+0	; position angle of 170 from Tully (1974) or 158 from Brown, 2014 or 12 from Hu , 2013
print_text = 1	; =1 to print steps to terminal

dim_wave = 20
wavelength = STRARR(dim_wave)
wavelength[0]  = '0.1542'	; GALEX.FUV	; counts/s	uv15
wavelength[1]  = '0.2274'	; GALEX.NUV	; counts/s	uv22
wavelength[2]  = '0.3562'	; SDSS.u	; nMgy (nanoMaggies)	uv36 	
wavelength[3]  = '0.4719'	; SDSS.g	; nMgy (nanoMaggies)	b
wavelength[4]  = '0.6185'	; SDSS.r	; nMgy (nanoMaggies)	v
wavelength[5]  = '0.7500'	; SDSS.i	; nMgy (nanoMaggies)	i
wavelength[6]  = '0.8961'	; SDSS.z	; nMgy (nanoMaggies)	-
wavelength[7]  = '1.200'	; 2MASS.J	; data-number units [DN]	j
wavelength[8]  = '1.600'	; 2MASS.H	; data-number units [DN]	-
wavelength[9]  = '2.200'	; 2MASS.K	; data-number units [DN]	k
wavelength[10] = '3.507'	; IRAC.I1	; MJy/sr	ir36
wavelength[11] = '4.437'	; IRAC.I2	; MJy/sr	ir45
wavelength[12] = '5.739'	; IRAC.I3	; MJy/sr	ir58
wavelength[13] = '7.927'	; IRAC.I4	; MJy/sr
wavelength[14] = '24.000'	; MIPS.24	; MJy/sr
wavelength[15] = '70.000'	; PACS.70	; Jy/pixel
wavelength[16] = '160.000'	; PACS.160	; Jy/pixel
wavelength[17] = '250.000'	; SPIRE.250	; Jy/beam
wavelength[18] = '350.000'	; SPIRE.350	; Jy/beam
wavelength[19] = '500.000'	; SPIRE.500	; Jy/beam

;Initialise arrays for total flux/luminosity as a function of wavelength
tot_flux_jy = dblarr(dim_wave)
tot_lum_whz = dblarr(dim_wave)

;Loop through wavelengths
FOR i = 0, dim_wave-1 DO BEGIN
map = mapdir+wavelength[i]+'um.fits'
mask = maskdir+'mask_'+wavelength[i]+'um.fits'
;--------------------------------------- MAP SIZE AND SCALES ----------------------------------------------
data = READFITS(map, hd)	; reads in the map for photometry to be done on
mask_data = READFITS(mask, hdm)
IF print_text EQ 1 THEN PRINT, 'READ - Map: '+map+', Mask: '+mask
ADXY, hd, ra, dec, xc, yc	; defines the x and y centre from the coordinates: ra and dec.
xc = FLOAT(xc[0])	; converts xc to float to avoid errors.
yc = FLOAT(yc[0])	; converts yc to flost to avoid errors.
dim = SIZE(data,/DIMENSIONS)	; 2d array with the size of the x and y axes.
nx = dim[0]	; number of pixels in the x axis.
ny = dim[1]	; number of pixels in the y axis.
X = FINDGEN(nx) # REPLICATE(1.0, ny)	;map of X pixel coordinates.
Y = REPLICATE(1.0, nx) # FINDGEN(ny)	;map of Y pixel coordinates.

pixsize_deg = ABS(fxpar(hd, 'CDELT1', count=count))	; x-axis resolution [degrees/pixel].
IF count EQ 0 THEN pixsize_deg = ABS(fxpar(hd, 'CD1_1', count=count))
pixsize_arcsec = pixsize_deg*3600.	; pixel size [arcsec]
pixsize_pc = distance_pc * TAN(!DPI/180 * pixsize_deg)     ; [pc/pixel].
pixsize_sr = (pixsize_pc/distance_pc)^2.	; converts pixsize from [pc/pixel] to [sr/pix] same as [deg^2/pix].

mapsize = pixsize_pc * nx	; size of the map in [pc].
IF print_text EQ 1 THEN PRINT, 'Pixel size of map: '+STRTRIM(pixsize_pc,1)+' [pc/pixel]'
IF print_text EQ 1 THEN PRINT, ''
;-----------------------------------------------------------------------------------------------------------
;---------------------------------------- CONVOLVE MAP -----------------------------------------------------
convolve_affix = ''
if i ge 2 and i le 13 then begin
convolve_affix = '_convolved'
if i le 6 then kernel=kerneldir+'PSF_Original_Gauss_03.0.fits' else kernel=kerneldir+'PSF_Original_Gauss_01.5.fits'
print, kernel
kernel = readfits(kernel, hdk)
pixsize_kernel = abs(fxpar(hdk,'CD1_1')*3600.)	; pixel size of kernel [arcsec]
if abs(pixsize_arcsec-pixsize_kernel)/pixsize_kernel gt 0.05 then begin
	size_kernel = (size(kernel))[1]
	size_new = round(double(size_kernel)*pixsize_kernel/pixsize_arcsec)
	if (size_new mod 2) eq 0 then size_new = size_new+1
endif
kernel = congrid(kernel+0.0, size_new, size_new, cubic=-0.5, /center)
kernel = kernel/total(kernel)	; normalises kernel
data[*,*] = convol(data[*,*], kernel)
print, 'Convolved '+map
endif
;----------------------------------------------------------------------------------------------------------
;-------------------------------------------- UNIT CONVERSION ----------------------------------------------
wave = DOUBLE(wavelength[i]) * 1D-6	; converts the wavelength from [um] to [m].
conv_factor = (299792458 / wave^2) * 1E-23 * 1E-10	; conversion between frequency and wavelength dependence[/s /A]/[Hz/A]
FUV_corr = 1.40E-15     ; from https://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html.
NUV_corr = 2.06E-16     ; from https://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html.
JMAGZP = 20.86689949    ; instrumental zero point for 2MASS data.
HMAGZP = 20.65309906    ; instrumental zero point for 2MASS data.
KMAGZP = 20.04360008    ; instrumental zero point for 2MASS data.
F0_J = 1594     ; from the 2MASS explanatory https://old.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html.
F0_H = 1024	; from the 2MASS explanatory https://old.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html.
F0_K = 666.7	; from the 2MASS explanatory https://old.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html.

IF i EQ 0 THEN data = data * ((FUV_corr / conv_factor) / (1D6 * pixsize_sr))	; converts from [counts/sec] to [MJy/sr].
IF i EQ 1 THEN data = data * ((NUV_corr / conv_factor) / (1D6 * pixsize_sr))	; converts from [counts/sec] to [MJy/sr].
IF i GE 2 AND i LE 6 THEN data = data * ((3.631 * 1D-6) / (1D6 * pixsize_sr))	; converts from [nMgy] to [MJy/sr].
IF i EQ 7 THEN data = data * F0_J * 10^(-JMAGZP / 2.5) / (1D6 * pixsize_sr)	; converts from [DN] to [MJy/sr].
IF i EQ 8 THEN data = data * F0_H * 10^(-HMAGZP / 2.5) / (1D6 * pixsize_sr)	; converts from [DN] to [MJy/sr].
IF i EQ 9 THEN data = data * F0_K * 10^(-KMAGZP / 2.5) / (1D6 * pixsize_sr)	; converts from [DN] to [MJy/sr].
IF i EQ 15 OR i EQ 16 THEN data = data / (1D6 * pixsize_sr)	; converts from [Jy/pix] to [MJy/sr].
IF i EQ 17 THEN data = data * (((pixsize_deg*3600)^2)/ 469.3542) / (1D6 * pixsize_sr)	; converts from [Jy/beam] to [MJy/sr] from SPIRE drg.
IF i EQ 18 THEN data = data * (((pixsize_deg*3600)^2)/ 831.275) / (1D6 * pixsize_sr)	; converts from [Jy/beam] to [MJy/sr] from SPIRE drg.
IF i EQ 19 THEN data = data * (((pixsize_deg*3600)^2)/ 1804.3058) / (1D6 * pixsize_sr)	; converts from [Jy/beam] to [MJy/sr] from SPIRE drg.
data_mjy = data
;writefits, mapdir+wavelength[i]+'um'+convolve_affix+'_MJy.fits', data_mjy, hd
;-----------------------------------------------------------------------------------------------------------
;----------------------------------------- C[II] CORRECTION ------------------------------------------------
;if wavelength[i] = '160.000' then begin
;dist_ellipse, el_map, dim, xc, yc, ratio, pos_ang
;map_extent = where(el_map lt 21000/pixsize_pc, count)
;data(*) = data(*) - ((24.18/(pixsize_sr*1e6))/count)
;endif
;-----------------------------------------------------------------------------------------------------------
;--------------------------------------- APPLY MASK TO THE MAP ---------------------------------------------
data[WHERE(mask_data EQ 1)] = !VALUES.F_NaN	; where regions are being masked, set the map data to NaN.
IF print_text EQ 1 THEN PRINT, 'Mask succesfully applied to the map.'
;------------------------------------- BACKGROUND CALCULATION ---------------------------------------------
av_bg = 0.	; initialise values
av_error_bg = 0.	; initialise values
bg_error = 0.	; initialise values
n_pix_bg = 0.
IF skip EQ 1 THEN GOTO, skip_background
IF i GE 7 AND i NE 10 AND i NE 11 AND i NE 12 AND i NE 15 THEN GOTO, skip_background
dist_ellipse, el_map, dim, xc, yc, ratio_background, pos_ang	; creates mask array for elliptical aperture photometry
back_in = 25000 / pixsize_pc	; inner radius of the ellipse for the background in units [pc]/[pixel/px]->[pix]
back_out = 30000 / pixsize_pc	; outer radius of the ellipse for the background in units [pc]/[pixel/px]->[pix]
IF print_text EQ 1 THEN PRINT, 'Begin background subraction...'
back_step = (back_out-back_in) / 6	; the distance  between the two limits dividing it into 6

iq = WHERE(el_map GT back_in AND el_map LT back_out AND FINITE(data) EQ 1, count_tot)
iq1 = WHERE(el_map GT back_in AND el_map LT back_in+back_step AND FINITE(data) EQ 1, count1)
iq2 = WHERE(el_map GT back_in+back_step AND el_map LT back_in+(2*back_step) AND FINITE(data) EQ 1, count2)
iq3 = WHERE(el_map GT back_in+(2*back_step) AND el_map LT back_in+(3*back_step) AND FINITE(data) EQ 1, count3)
iq4 = WHERE(el_map GT back_in+(3*back_step) AND el_map LT back_in+(4*back_step) AND FINITE(data) EQ 1, count4)
iq5 = WHERE(el_map GT back_in+(4*back_step) AND el_map LT back_in+(5*back_step) AND FINITE(data) EQ 1, count5)
iq6 = WHERE(el_map GT back_in+(5*back_step) AND el_map LT back_in+(6*back_step) AND FINITE(data) EQ 1, count6)
n_pix_bg = [count1, count2, count3, count4, count5, count6]

av_bg_total = TOTAL(data[iq],/NaN)/count_tot
av_bg1 = TOTAL(data[iq1],/NaN)/count1	; average background within the annulus, eqn (A1)
av_bg2 = TOTAL(data[iq2],/NaN)/count2
av_bg3 = TOTAL(data[iq3],/NaN)/count3
av_bg4 = TOTAL(data[iq4],/NaN)/count4
av_bg5 = TOTAL(data[iq5],/NaN)/count5
av_bg6 = TOTAL(data[iq6],/NaN)/count6
av_bg = (av_bg1 + av_bg2 + av_bg3 + av_bg4 + av_bg5 + av_bg6)/6	; eqn (A4) in Thirlwall2020
av_bg_stddev = STDDEV([av_bg1, av_bg2, av_bg3, av_bg4, av_bg5, av_bg6])	; eqn (A5) (sample standard deviation)
bg_error = av_bg_stddev / SQRT(6)	; eqn (A6) [MJy/sr]

bg_stddev1 = STDDEV(data[iq1],/NaN)	; standard deviation within the annulus, eqn (A2)
bg_stddev2 = STDDEV(data[iq2],/NaN)
bg_stddev3 = STDDEV(data[iq3],/NaN)
bg_stddev4 = STDDEV(data[iq4],/NaN)
bg_stddev5 = STDDEV(data[iq5],/NaN)
bg_stddev6 = STDDEV(data[iq6],/NaN)

error_bg1 = bg_stddev1 / SQRT(count1)	; error in standard deviation for he annulus, eqn (A3)
error_bg2 = bg_stddev2 / SQRT(count2)
error_bg3 = bg_stddev3 / SQRT(count3)
error_bg4 = bg_stddev4 / SQRT(count4)
error_bg5 = bg_stddev5 / SQRT(count5)
error_bg6 = bg_stddev6 / SQRT(count6)
av_error_bg = (error_bg1 + error_bg2 + error_bg3 + error_bg4 + error_bg5 + error_bg6)/6

print_background_map = 0
IF print_background_map EQ 1 THEN BEGIN
	background_aperture = (data * !VALUES.F_NaN)
	;background_aperture = data * 0
	width_pc = 100	; half the width of the outline of the background aperture
	IF i GE 17 THEN width_pc = 300
	width = width_pc / pixsize_pc
	background_aperture[WHERE(el_map GE back_in-width AND el_map LE back_in+width)] = 1.
	background_aperture[WHERE(el_map GE back_out-width AND el_map LE back_out+width)] = 1.
	back_name = plotdir+'background_aperture_'+wavelength[i]+'um.fits'
	writefits, back_name, background_aperture, hd
	PRINT, 'WRITTEN: '+back_name
ENDIF
skip_background:
IF print_text EQ 1 THEN PRINT, 'Average background: '+STRTRIM(av_bg,1)+' [MJy/sr]'
IF print_text EQ 1 THEN PRINT, 'Average background error: '+STRTRIM(av_error_bg * 1E3,1)+' [kJy/sr]'
IF print_text EQ 1 THEN PRINT, ''
;--------------------------------------------------------------------------------------------------------------
;--------------------------------------------- ANNULI FOR ANALYSIS --------------------------------------------
dist_ellipse, el_map, dim, xc, yc, ratio, pos_ang       ; creates mask array for elliptical aperture photometry
print_aperture_map = 0
IF print_aperture_map EQ 1 THEN BEGIN
	rad_lim_pc = [0, 5, 10, 15] *1E3      ; radii in [pc] from the centre for which a ring is to be marked for analysis purposes.
	width_pc = 100	; half of the width of the annulus from the line above.
	IF i GE 17 THEN width_pc = 300
	rad_arr = rad_lim_pc / pixsize_pc    ; converts the annulus distance to [pixels].
	width = width_pc / pixsize_pc	; converts width_pc from [pc] to [pixels].
	rad_lim_map = el_map * !VALUES.F_NaN
	ii=0
	FOR ii = 0, N_ELEMENTS(rad_arr)-1 DO BEGIN
		mark = WHERE(el_map GE rad_arr[ii]-width AND el_map LE rad_arr[ii]+width)
		IF mark[0] EQ -1 THEN CONTINUE
		rad_lim_map[mark] = 1.
	ENDFOR
	map_mark_name = plotdir+'map_marks_'+wavelength[i]+'um.fits'
	writefits, map_mark_name, rad_lim_map, hd
	PRINT, 'WRITTEN: '+map_mark_name
ENDIF
;-------------------------------------------------------------------------------------------------------------
;----------------------------------------------- SURFACE PHOTOMETRY ------------------------------------------
stepsize_kpc = DBLARR(250)
stepsize_kpc[*] = 0.1
IF i GE 15 THEN BEGIN	; changes step size for PACS & SPIRE data due to low resolution
	stepsize_kpc = DBLARR(100)
	stepsize_kpc[*] = 0.3
ENDIF
stepsize_pc = stepsize_kpc * 1E3	; units [pc]
stepsize_pix = stepsize_pc / pixsize_pc	; units [pix]
;PRINT, 'stepsize_kpc total = '+STRTRIM(TOTAL(stepsize_kpc),1)+' kpc or '+STRTRIM(TOTAL(stepsize_pix),1)+' pix'

nsteps = N_ELEMENTS(stepsize_kpc)	; the total number of annuli used for the photometry
av_int_rad = DBLARR(nsteps)	; the average flux within an annulus
av_int_rad1 = DBLARR(nsteps)
av_int_rad2 = DBLARR(nsteps)
av_int_rad3 = DBLARR(nsteps)
av_int_rad4 = DBLARR(nsteps)
int_rad = DBLARR(nsteps)	; the total flux within an annulus
n_pix = DBLARR(nsteps)	; array used to contain the number of pixels used for each annulus in the photometry
x_rad_pix = DBLARR(nsteps)	; radial position in the galaxy in [pixels]
x_rad = DBLARR(nsteps)	; radial position in the galaxy in [kpc]
x_rad_pix[0] = 0
x_rad[0] = 0

sb_bg_error = DBLARR(nsteps)
av_int_rad_mean = DBLARR(nsteps)
flux_stddev = DBLARR(nsteps)
FOR j = 1, nsteps-1 DO BEGIN
	stepdiff = stepsize_pix[j] - stepsize_pix[j-1]	; defines the midpoint of the annuli
        x_rad_pix[j] = x_rad_pix[j-1] + stepsize_pix[j]    ;defines array of pixels in units of stepsize
        x_rad[j] = x_rad_pix[j-1] + stepdiff ; later to be converted from [pixels] to [kpc] (element j is at the midpoint of the annulus)

        iq = WHERE((el_map GE x_rad_pix[j-1]) AND (el_map LT x_rad_pix[j]) AND FINITE(data), count)      ;pixels that lie within the annulus
	IF MIN(iq) EQ -1 THEN CONTINUE
	av_int_rad[j] = ((TOTAL(data[iq],/NaN)-(count*av_bg)) / count)	; average flux within annulus [MJy/sr/pix]
	;int_rad[j] = int_rad[j-1] + (av_int_rad[j] * count)	; summed averaged fluxes [MJy/sr]
	int_rad[j] = int_rad[j-1] + total(data[iq], /NaN)

iq1 = WHERE((el_map GE x_rad_pix[j-1]) AND (el_map LT x_rad_pix[j]) AND Y GT yc AND x LT xc AND FINITE(data), count1)
IF iq1[0] NE -1 THEN av_int_rad1[j] = ((TOTAL(data[iq1],/NaN)-(count1*av_bg)) / count1)      ; eqn. A10

iq2 = WHERE((el_map GE x_rad_pix[j-1]) AND (el_map LT x_rad_pix[j]) AND Y GT yc AND x GT xc AND FINITE(data), count2)
IF iq2[0] NE -1 THEN av_int_rad2[j] = ((TOTAL(data[iq2],/NaN)-(count2*av_bg)) / count2)

iq3 = WHERE((el_map GE x_rad_pix[j-1]) AND (el_map LT x_rad_pix[j]) AND Y LT yc AND x LT xc AND FINITE(data), count3)
IF iq3[0] NE -1 THEN av_int_rad3[j] = ((TOTAL(data[iq3],/NaN)-(count3*av_bg)) / count3)

iq4 = WHERE((el_map GE x_rad_pix[j-1]) AND (el_map LT x_rad_pix[j]) AND Y LT yc AND x GT xc AND FINITE(data), count4)
IF iq4[0] NE -1 THEN av_int_rad4[j] = ((TOTAL(data[iq4],/NaN)-(count4*av_bg)) / count4)

	n_pix[j] = count	; the number of pixels within annulus j [pix].
	sb_bg_error[j] = bg_error * SQRT(TOTAL(n_pix_bg) / n_pix[j])	; eqn. A9, error in SB due to background.

	av_int_rad_mean[j] = (av_int_rad1[j] + av_int_rad2[j] + av_int_rad3[j] + av_int_rad4[j]) / 4	; eqn. A11
	flux_stddev[j] = STDDEV([av_int_rad1[j], av_int_rad2[j], av_int_rad3[j], av_int_rad4[j]],/NaN)	; eqn. A12
ENDFOR
av_int_rad[0] = !VALUES.F_NAN
x_rad = x_rad * pixsize_pc / 1E3        ; converts x_rad from [pixels] to [kpc]
;---------------------------------------------------------------------------------------------------------------
;------------------------------------ FOREGROUND EXTINCTION CORRECTION -----------------------------------------
ebv = 0.0308	; value for E(B-V) from NASA/IPAC irsa
fm_unred, wavelength[i]*10000., 1., Ebv, flux_unred       ; IDL function to de-redden light. wavelength must be in angstroms
a_corr = 2.5 * alog10(flux_unred)       ; converts the corrected flux into magnitudes
av_int_rad = av_int_rad * flux_unred	; dereddens average flux.
int_rad = int_rad * flux_unred * pixsize_sr	; dereddens and multiplies by [sr] to give units =>[MJy]
av_bg = av_bg * flux_unred	; dereddens the calculated background value.
data_mjy = (data_mjy * flux_unred) - av_bg
IF print_text EQ 1 THEN PRINT, 'All flux values de-reddened.'
;---------------------------------------------------------------------------------------------------------------
;------------------------------------------ ERROR CALCULATIONS--------------------------------------------------
IF i LE 1 THEN cal_percent = 0.1	; GALEX
IF i GE 2 AND i LE 6 THEN cal_percent = 0.02	; SDSS
IF i GE 7 AND i LE 13 THEN cal_percent = 0.03	; 2MASS & IRAC
IF i EQ 14 THEN cal_percent = 0.01	; MIPS
IF i EQ 15 THEN cal_percent = 0.1	; PACS70
IF i EQ 16 THEN cal_percent = 0.2	; PACS160
IF i ge 17 THEN cal_percent = 0.15	; SPIRE250, 350, 500

cal_error = av_int_rad[*] * cal_percent
cal_error_tot = int_rad[*] * cal_percent

bg_error_tot = bg_error * TOTAL(n_pix_bg, /NaN) * pixsize_sr      ; eqn. A7, error in total flux density due to background [MJy/sr]->[MJy]
n_pix[where(x_rad ge 17.5)] = 0.        ; ignores annuli beyond a critical radius from being considered
poisson = SQRT(total(n_pix))    ; Poisson noise
poisson = 0.

flux_density_error = SQRT(cal_error_tot^2 +bg_error_tot^2 + poisson^2)  ; eqn. A8, total error in flux densities
print, 'cal_error_tot: ', cal_error_tot[(where(x_rad ge 17.5))[0]]
print, 'bg_error_tot: ', bg_error_tot
print, 'poisson: ', poisson

config_error = flux_stddev / SQRT(4)    ; eqn. A13
total_error = SQRT(cal_error^2. + sb_bg_error[*]^2. + config_error[*]^2.)       ; eqn. A14
total_error[WHERE(x_rad GE 17.5)] = 0.
int_rad_stddev = STDDEV(int_rad,/NaN)
PRINT, 'Calibration error: '+STRTRIM(MEAN(cal_error, /NaN),1)
PRINT, 'Background error: '+STRTRIM(MEAN(sb_bg_error, /NaN),1)
PRINT, 'Configuration error: '+STRTRIM(MEAN(config_error, /NaN),1)
;for kk = 0, n_elements(total_error)-1 do begin
;	print, x_rad(kk), cal_error(kk), sb_bg_error(kk), config_error(kk)
;endfor
;---------------------------------------------------------------------------------------------------------------
;---------------------------------------------------------------------------------------------------------------
limit = WHERE(x_rad GE 21.)	; array of positions beyond limit raidus (limit(0) = furtherest extend of photometry)
tot_lum = int_rad[limit[0]]*1d+6	; total flux of M51a at limit radius
tot_flux_jy[i] = int_rad[limi[0]]*1d+6	; total flux of M51a at limit radius
wave = DOUBLE(wavelength[i]) * 1D-6     ; converts the wavelength from [um] to [m]
conv_factor = (299792458 / wave^2)  ; conversion used for CpS to Jy [/m /s] or [Hz /m]
tot_lum_watts = tot_lum * 1E6 * 1E-26 * conv_factor * wave	; converts from [MJy] => [W/m^2]
IF print_text EQ 1 THEN PRINT, 'Total luminosity: '+STRTRIM(tot_lum,1)+' [Jy], or '+STRTRIM(tot_lum_watts,1)+' [W/m^2]'
;---------------------------------------------------------------------------------------------------------------
IF wavelength[i] EQ '0.1542'  THEN wave_save = 'FUV'	; GALEX.FUV	; counts/s 
IF wavelength[i] EQ '0.2274'  THEN wave_save = 'NUV'	; GALEX.NUV	; counts/s
IF wavelength[i] EQ '0.3562'  THEN wave_save = 'u'
IF wavelength[i] EQ '0.4719'  THEN wave_save = 'g'
IF wavelength[i] EQ '0.6185'  THEN wave_save = 'r'
IF wavelength[i] EQ '0.7500'  THEN wave_save = 'i'
IF wavelength[i] EQ '0.8961'  THEN wave_save = 'z'
IF wavelength[i] EQ '1.200'   THEN wave_save = 'J'
IF wavelength[i] EQ '1.600'   THEN wave_save = 'H'
IF wavelength[i] EQ '2.200'   THEN wave_save = 'K'
IF wavelength[i] EQ '3.507'   THEN wave_save = '3.5'
IF wavelength[i] EQ '4.437'   THEN wave_save = '4.4'
IF wavelength[i] EQ '5.739'   THEN wave_save = '5.8'
IF wavelength[i] EQ '7.927'   THEN wave_save = '8'
IF wavelength[i] EQ '24.000'  THEN wave_save = '24'
IF wavelength[i] EQ '70.000'  THEN wave_save = '70'
IF wavelength[i] EQ '160.000' THEN wave_save = '160'
IF wavelength[i] EQ '250.000' THEN wave_save = '250'
IF wavelength[i] EQ '350.000' THEN wave_save = '350'
IF wavelength[i] EQ '500.000' THEN wave_save = '500'

wavelength[2]  = '0.3562'       ; SDSS.u        ; nMgy (nanoMaggies)    uv36    
wavelength[3]  = '0.4719'       ; SDSS.g        ; nMgy (nanoMaggies)    b
wavelength[4]  = '0.6185'       ; SDSS.r        ; nMgy (nanoMaggies)    v
wavelength[5]  = '0.7500'       ; SDSS.i        ; nMgy (nanoMaggies)    i
wavelength[6]  = '0.8961'       ; SDSS.z        ; nMgy (nanoMaggies)    
wavelength[7]  = '1.200'        ; 2MASS.J       ; data-number units [DN]        j
wavelength[8]  = '1.600'        ; 2MASS.H       ; data-number units [DN]        
wavelength[9]  = '2.200'        ; 2MASS.K       ; data-number units [DN]        k
wavelength[10] = '3.507'        ; IRAC.I1       ; MJy/sr        ir36
wavelength[11] = '4.437'        ; IRAC.I2       ; MJy/sr        ir45
wavelength[12] = '5.739'        ; IRAC.I3       ; MJy/sr        ir58
wavelength[13] = '7.927'        ; IRAC.I4       ; MJy/sr        
wavelength[14] = '24.000'       ; MIPS.24       ; MJy/sr        
wavelength[15] = '70.000'       ; PACS.70       ; Jy/pixel      
wavelength[16] = '160.000'      ; PACS.160      ; Jy/pixel      
wavelength[17] = '250.000'      ; SPIRE.250     ; Jy/beam       
wavelength[18] = '350.000'      ; SPIRE.350     ; Jy/beam       
wavelength[19] = '500.000'      ; SPIRE.500     ; Jy/beam 

;IF i GE 2 AND i LE 13 THEN BEGIN
;        WRITEFITS, mapdir+wavelength[i]+'um_convolved_MJy.fits', data_mjy, hd
;ENDIF ELSE BEGIN
;        WRITEFITS, mapdir+wavelength[i]+'um_MJy.fits', data_mjy, hd
;ENDELSE
save_name = savedir+'datafile_'+wavelength[i]+'um.save'
save_name = './TEMP/obs_tot_profile_'+wave_save+'.save'
;IF skip EQ 1 THEN save_name = savedir+'datafile_no_background_'+wavelength[i]+'um.save'

SAVE, x_rad, av_int_rad, int_rad, int_rad_stddev, total_error, bg_error, cal_error_tot, $
	pixsize_deg, distance_pc, ra, dec, ratio, pos_ang, ebv, flux_unred, $
	DESCRIPTION='Radius [kpc]	Average intensity [MJy/sr]', FILENAME=save_name
PRINT, 'SAVED FILE: '+save_name
PRINT, '-----------------------------------------------------------------------------'
ENDFOR
tot_lum_whz = tot_lum_jy * flux_to_lum	; convert total flux to total luminosity [W/Hz]
;convert wavlength to frequency
nu = c / double(wavelength)
;difference between frequency n - n+1
dnu = nu[0:dim_wave-2] - nu[1:dim_wave-1]
;average lum of these two points
av_lum = (tot_lum_whz[1:dim_wave-1] + tot_lum_whz[0:dim_wave-2]) * 2d
int_lum = av_lum *dnu	; THIS CONSIDERS UNSED BANDS IN THE INTEGRATION
PRINT, '------------------------------- DONE ----------------------------------------'
stop
END
