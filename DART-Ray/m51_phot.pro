PRO m51_phot
compile_opt idl2
; The purpose of this code is to carry produce azimuthally averaged SB profiles of the dust emission of M51a.
; Maps produced from 2Dto3D projection in DART-Ray.
; This code reads in the model maps for every component produced from DART-Ray and, where necessary, convolves all components
; with the appropriate PSF. Then carries out the azimuthally averaged photometry and combines the profiles with the stellar
; emission for the NIR wavelengths. A chi^2 value is calculated then all relevant data is saved.

root = '../'
nurad_dir = root+'NUrad/'
mapdir = root+'DART-Ray/MAPS/M51a/'	; directory containing surface brightness profiles.
griddir = root+'DART-Ray/RUNS/GALAXY_M51a/'
savedir = root+'saves/model/'	; directory for the .save files containing the SB profiles.
kernel_dir = nurad_dir+'maps/obs/kernels/'
figdir = root+'figures/'
print_text = 0	; =1 for text to be printed to terminal.

read_scaling, model, qyear, scaabs, tau, nsersic, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5, f_bd

distance = 8.58E6    ; from gal_param.in 12/02/21.
inclination = 20.3*!DPI/180.	; Same as the value Theta, in dir_out_galaxy.dat.
ratio = (COS(inclination))^(-1)	; ratio of major and minor axis.
pos_ang = 90.   ;position angle of major axis (not the same as observation, specific for DART-Ray output).

map = griddir+'grid_m2dto3d_m51a_tot_maps_part2_int.h5'	; read in grid for finding model size [pc]
data_map = H5F_OPEN(map)
map_size = H5D_READ(H5D_OPEN(data_map, 'size_map'))     ; size of the model map in [pc]
H5F_CLOSE, data_map

wavelength=STRARR(11)
wavelength[0]  = '3.6'	;stellar + dust ;IRAC_3.4
wavelength[1]  = '4.5'	;stellar + dust ;IRAC_4.5
wavelength[2]  = '5.8'	;stellar + dust
wavelength[3]  = '8.0'	;dust
wavelength[4]  = '24'	;dust
wavelength[5]  = '70'	;dust
wavelength[6]  = '100'	;dust
wavelength[7]  = '160'	;dust
wavelength[8]  = '250'	;dust
wavelength[9]  = '350'	;dust
wavelength[10] = '500'	;dust

;read, convolve_model, prompt="Convolve? [1]=yes, [0]=no: "
convolve_model = 0
FOR i = 0, 10 DO BEGIN
IF wavelength[i] EQ '100' THEN CONTINUE	; skips 100um maps as there is no data at this wavelength.
morph = ['tot','irr1','irr2','irr3','irr4','irr5','irr6','irr_HII','irr_HII4','irr_HII6','irr_HII7']
map = STRARR(N_ELEMENTS(morph))	; defines an array to fit all morphological components
morph_check = LONARR(N_ELEMENTS(morph))	; checks if maps for a given morphology is present
FOR k = 0, N_ELEMENTS(morph)-1 DO BEGIN
	map[k] = mapdir+model+'_map_'+morph[k]+'_'+wavelength[i]+'um_'+scaabs+'.fits'
	morph_check[k] = file_test(map[k])	; checks if maps for a given morphology is present
ENDFOR

irr1_data = READFITS(map[WHERE(morph EQ 'irr1')],hd)
irr2_data = READFITS(map[WHERE(morph EQ 'irr2')],hd)
irr3_data = READFITS(map[WHERE(morph EQ 'irr3')],hd)
irr4_data = READFITS(map[WHERE(morph EQ 'irr4')],hd)
irr5_data = READFITS(map[WHERE(morph EQ 'irr5')],hd)
irr6_data = READFITS(map[WHERE(morph EQ 'irr6')],hd)
IF morph_check[7] EQ 1 THEN HII_data  = READFITS(map[WHERE(morph EQ 'irr_HII')],hd) ELSE HII_data = 0
IF morph_check[8] EQ 1 THEN HII4_data = READFITS(map[WHERE(morph EQ 'irr_HII4')],hd) ELSE HII4_data = 0
IF morph_check[9] EQ 1 THEN HII6_data = READFITS(map[WHERE(morph EQ 'irr_HII6')],hd) ELSE HII6_data = 0
IF morph_check[10] EQ 1 THEN HII7_data = READFITS(map[WHERE(morph EQ 'irr_HII7')],hd) ELSE HII7_data = 0
tot_data  = READFITS(map[WHERE(morph EQ 'tot')],hd)
PRINT, 'READ: '+model+' '+wavelength[i]+'um maps'
PRINT, 'Components found: ', morph[WHERE(morph_check EQ 1)]
PRINT, ''
;----------------------------------------------------------------------------------------------------------------
;------------------------------------ DEFINE MAPS SIZES AND RESOLUTIONS -----------------------------------------
dim = SIZE(tot_data,/DIMENSIONS)    ; 2d array with the size of the x and y axes.
nx = dim[0]     ; number of pixels in the x axis.
ny = dim[1]     ; number of pixels in the y axis.
xc = (nx /2) - 1	;determines the x centre of the map.
yc = (ny /2) - 1	;determines the y centre of the map.
X = findgen(nx) # replicate(1.0, ny)  ;map of X pixel coordinates.
Y = replicate(1.0, nx) # findgen(ny)  ;map of Y pixel coordinates.
pixsize_pc = FLOAT(map_size[0] / nx)	; pixel size of the map in [pc/pixel].
pixsize_arcsec = ATAN(pixsize_pc/distance) * 3600. * 180. / !DPI	; pixel size in [arcsec/pixel]
pix_area_sr = (pixsize_pc/distance)^2.; solid_angle=Area/radius^2
IF print_text EQ 1 THEN PRINT, 'Pixel size of the map: '+STRTRIM(pixsize_pc,1)+' pc, '+STRTRIM(pixsize_arcsec,1)+' arcsec'
IF print_text EQ 1 THEN PRINT, 'Map size: '+STRTRIM(map_size,1)+' pc'
IF print_text EQ 1 THEN PRINT, ''
;-------------------------------- CONVOLVE MODEL MAP -------------------------------------
;convolve_model = 1	; set to 1 to convolve the model, set to 0 to skip this step.
conv_affix = ''
IF convolve_model NE 1 THEN GOTO, skip_convol
kernel = ''
IF wavelength[i] EQ '24' THEN kernel = kernel_dir+'PSF_resized_MIPS_24.fits'
IF wavelength[i] EQ '70' THEN kernel = kernel_dir+'PSF_resized_PACS_70.fits'
IF wavelength[i] EQ '160' THEN kernel = kernel_dir+'PSF_resized_PACS_160.fits'
IF wavelength[i] EQ '250' THEN kernel = kernel_dir+'PSF_resized_SPIRE_250.fits'
IF wavelength[i] EQ '350' THEN kernel = kernel_dir+'PSF_resized_SPIRE_350.fits'
IF wavelength[i] EQ '500' THEN kernel = kernel_dir+'PSF_resized_SPIRE_500.fits'
IF kernel EQ '' THEN GOTO, skip_convol
IF print_text EQ 1 THEN PRINT, 'Convolving '+wavelength[i]+'um maps with '+kernel

conv_affix = '_convolved'	;  adds this affix to some save files.
IF i LT 80 THEN kernel_image = READFITS(kernel, hdk)	; reads in the kernel and its header.
;IF i GE 8 THEN kernel_image = READFITS(kernel, hdk, /EXTEN)	; reads in the kernel and its header.
pixsize_kernel = ABS(fxpar(hdk,'CDELT1',count=count)*3600.)	; pixel size of the kernel in [''].
IF count eq 0 THEN pixsize_kernel = ABS(fxpar(hdk,'CD1_1',count=count)*3600.)	; finds alternate pixel size.

IF (ABS(pixsize_arcsec-pixsize_kernel)/pixsize_kernel) GT 0.05 THEN BEGIN	; checks if difference between scales is 5%.
        size_kernel = (SIZE(kernel_image))[1]
        size_new = ROUND(FLOAT(size_kernel)*pixsize_kernel/pixsize_arcsec)
        PRINT, 'Pixel size of kernel - '+STRTRIM(pixsize_kernel,1)+' model - '+STRTRIM(pixsize_arcsec,1)
        PRINT, 'Resizing kernel from '+STRTRIM(size_kernel,1)+' to '+STRTRIM(size_new,1)
        IF (size_new MOD 2) EQ 0 THEN size_new = size_new+1	; ensures the size of the kernel is odd.
ENDIF ELSE BEGIN
	PRINT, 'Kernel and image are on the same pixel scale of '+STRTRIM(pixsize_kernel,1)+' [arcsec/pix]!'
	size_new = (SIZE(kernel_image))[1]
	PRINT, 'Kernel size is '+STRTRIM(size_new,1)+' x '+STRTRIM(size_new,1)+' [pix]'
ENDELSE

kernel_image = CONGRID(kernel_image+0.0, size_new, size_new, cubic=-0.5,/center)	; re-scales the kernel to new pixel size.
kernel_image = kernel_image/TOTAL(kernel_image)	; normalises the kernel to converve flux.
nx_pad = nx+200	; add pixels to the map to provide a buffer for the convolution.
ny_pad = ny+200
xc_pad = nx_pad /2
yc_pad = ny_pad /2
tot_temp  = DBLARR(nx_pad, ny_pad)	; creates a temporary padded array for the map
irr1_temp = DBLARR(nx_pad, ny_pad)
irr2_temp = DBLARR(nx_pad, ny_pad)
irr3_temp = DBLARR(nx_pad, ny_pad)
irr4_temp = DBLARR(nx_pad, ny_pad)
irr5_temp = DBLARR(nx_pad, ny_pad)
irr6_temp = DBLARR(nx_pad, ny_pad)
IF morph_check[7] EQ 1 THEN HII_temp  = DBLARR(nx_pad, ny_pad)
IF morph_check[8] EQ 1 THEN HII4_temp = DBLARR(nx_pad, ny_pad)
IF morph_check[9] EQ 1 THEN HII6_temp = DBLARR(nx_pad, ny_pad)
IF morph_check[10] EQ 1 THEN HII7_temp = DBLARR(nx_pad, ny_pad)
tot_zero  = WHERE(tot_data EQ 0)		; finds where the map is zero to restore proper truncation after convolution.
irr1_zero = WHERE(irr1_data EQ 0)
irr2_zero = WHERE(irr2_data EQ 0)
irr3_zero = WHERE(irr3_data EQ 0)
irr4_zero = WHERE(irr4_data EQ 0)
irr5_zero = WHERE(irr5_data EQ 0)
irr6_zero = WHERE(irr6_data EQ 0)
IF morph_check[7] EQ 1 THEN HII_zero  = WHERE(HII_data EQ 0)
IF morph_check[8] EQ 1 THEN HII4_zero = WHERE(HII4_data EQ 0)
IF morph_check[9] EQ 1 THEN HII6_zero = WHERE(HII6_data EQ 0)
IF morph_check[10] EQ 1 THEN HII7_zero = WHERE(HII7_data EQ 0)
FOR kk = 0, nx-1 DO BEGIN
	FOR jj = 0, ny-1 DO BEGIN
		tot_temp[kk+xc_pad-xc, jj+yc_pad-yc] = tot_data[kk,jj]
		irr1_temp[kk+xc_pad-xc, jj+yc_pad-yc] = irr1_data[kk,jj]
		irr2_temp[kk+xc_pad-xc, jj+yc_pad-yc] = irr2_data[kk,jj]
		irr3_temp[kk+xc_pad-xc, jj+yc_pad-yc] = irr3_data[kk,jj]
		irr4_temp[kk+xc_pad-xc, jj+yc_pad-yc] = irr4_data[kk,jj]
		irr5_temp[kk+xc_pad-xc, jj+yc_pad-yc] = irr5_data[kk,jj]
		irr6_temp[kk+xc_pad-xc, jj+yc_pad-yc] = irr6_data[kk,jj]
		IF morph_check[7] EQ 1 THEN HII_temp[kk+xc_pad-xc, jj+yc_pad-yc] = HII_data[kk,jj]
		IF morph_check[8] EQ 1 THEN HII4_temp[kk+xc_pad-xc, jj+yc_pad-yc] = HII4_data[kk,jj]
		IF morph_check[9] EQ 1 THEN HII6_temp[kk+xc_pad-xc, jj+yc_pad-yc] = HII6_data[kk,jj]
		IF morph_check[10] EQ 1 THEN HII7_temp[kk+xc_pad-xc, jj+yc_pad-yc] = HII7_data[kk,jj]
	ENDFOR
ENDFOR
tot_temp[*,*]  = CONVOL(tot_temp[*,*], kernel_image)
irr1_temp[*,*]  = CONVOL(irr1_temp[*,*], kernel_image)
irr2_temp[*,*]  = CONVOL(irr2_temp[*,*], kernel_image)
irr3_temp[*,*]  = CONVOL(irr3_temp[*,*], kernel_image)
irr4_temp[*,*]  = CONVOL(irr4_temp[*,*], kernel_image)
irr5_temp[*,*]  = CONVOL(irr5_temp[*,*], kernel_image)
irr6_temp[*,*]  = CONVOL(irr6_temp[*,*], kernel_image)
IF morph_check[7] EQ 1 THEN HII_temp[*,*]  = CONVOL(HII_temp[*,*], kernel_image)
IF morph_check[8] EQ 1 THEN HII4_temp[*,*]  = CONVOL(HII4_temp[*,*], kernel_image)
FOR kk = 0, nx-1 DO BEGIN
        FOR jj = 0, ny-1 DO BEGIN
                tot_data[kk,jj] = tot_temp[kk+xc_pad-xc, jj+yc_pad-yc]
		irr1_data[kk,jj] = irr1_temp[kk+xc_pad-xc, jj+yc_pad-yc]
		irr2_data[kk,jj] = irr2_temp[kk+xc_pad-xc, jj+yc_pad-yc]
		irr3_data[kk,jj] = irr3_temp[kk+xc_pad-xc, jj+yc_pad-yc]
		irr4_data[kk,jj] = irr4_temp[kk+xc_pad-xc, jj+yc_pad-yc]
		irr5_data[kk,jj] = irr5_temp[kk+xc_pad-xc, jj+yc_pad-yc]
		irr6_data[kk,jj] = irr6_temp[kk+xc_pad-xc, jj+yc_pad-yc]
		IF morph_check[7] EQ 1 THEN HII_data[kk,jj] = HII_temp[kk+xc_pad-xc, jj+yc_pad-yc]
		IF morph_check[8] EQ 1 THEN HII4_data[kk,jj] = HII4_temp[kk+xc_pad-xc, jj+yc_pad-yc]
		IF morph_check[9] EQ 1 THEN HII6_data[kk,jj] = HII6_temp[kk+xc_pad-xc, jj+yc_pad-yc]
		IF morph_check[10] EQ 1 THEN HII7_data[kk,jj] = HII7_temp[kk+xc_pad-xc, jj+yc_pad-yc]
        ENDFOR
ENDFOR
tot_data[tot_zero] = 0		; restores zeroes to regions that have been truncated.
irr1_data[irr1_zero] = 0
irr2_data[irr2_zero] = 0
irr3_data[irr3_zero] = 0
irr4_data[irr4_zero] = 0
irr5_data[irr5_zero] = 0
irr6_data[irr6_zero] = 0
IF morph_check[7] EQ 1 THEN HII_data[HII_zero] = 0
IF morph_check[8] EQ 1 THEN HII4_data[HII4_zero] = 0
IF morph_check[9] EQ 1 THEN HII6_data[HII6_zero] = 0
IF morph_check[10] EQ 1 THEN HII7_data[HII7_zero] = 0

IF print_text EQ 1 THEN PRINT, 'All components successfully convolved!' 
write_conv_maps = 0
IF write_conv_maps EQ 1 THEN BEGIN
	sxaddpar, hd_out, 'NAXIS', 2
	sxaddpar, hd_out, 'NAXIS1', nx ; changes the value 'NAXIS1' in the header to the correct value.
	sxaddpar, hd_out, 'NAXIS2', ny ; changes the value 'NAXIS2' in the header to the correct value.
	sxaddpar, hd_out, 'CDELT1', pixsize_arcsec / 3600.
	sxaddpar, hd_out, 'CDELT2', pixsize_arcsec / 3600.
	writefits, mapdir+model+'_map_tot_'+wavelength[i]+'um_conv.fits', tot_data, hd_out
	;writefits, mapdir+model+'_map_irr1_'+wavelength[i]+'um_conv.fits', irr1_data
	;writefits, mapdir+model+'_map_irr2_'+wavelength[i]+'um_conv.fits', irr2_data
	;writefits, mapdir+model+'_map_irr3_'+wavelength[i]+'um_conv.fits', irr3_data
	;writefits, mapdir+model+'_map_irr4_'+wavelength[i]+'um_conv.fits', irr4_data
	;writefits, mapdir+model+'_map_irr5_'+wavelength[i]+'um_conv.fits', irr5_data
	;writefits, mapdir+model+'_map_irr6_'+wavelength[i]+'um_conv.fits', irr6_data
	;IF morph_check[7] EQ 1 THEN writefits, mapdir+model+'_map_HII_'+wavelength[i]+'um_conv.fits', HII_data
	;IF morph_check[8] EQ 1 THEN writefits, mapdir+model+'_map_HII4_'+wavelength[i]+'um_conv.fits', HII4_data
	;IF morph_check[9] EQ 1 THEN writefits, mapdir+model+'_map_HII6_'+wavelength[i]+'um_conv.fits', HII6_data
	;IF morph_check[10] EQ 1 THEN writefits, mapdir+model+'_map_HII7_'+wavelength[i]+'um_conv.fits', HII7_data
	PRINT, 'Written convolved maps to '+mapdir
ENDIF
IF print_text EQ 1 THEN PRINT, ''
skip_convol:
;-----------------------------------------------------------------------------------------
dist_ellipse, el_map, dim, xc, yc, ratio, pos_ang

print_aperture_map = 0
IF print_aperture_map EQ 1 THEN BEGIN
        rad_lim_pc = [5., 10., 15., 20., 25.] *1E3      ; radii in [pc] from the centre for which a ring is to be marked for analysis purposes.
        width_pc = 75   ; half of the width of the annulus from the line above.
        rad_arr = rad_lim_pc / pixsize_pc    ; converts the annulus distance to [pixels].
        width = width_pc / pixsize_pc   ; converts width_pc from [pc] to [pixels].
        rad_lim_map = el_map * !VALUES.F_NaN
        ii=0
        FOR ii = 0, N_ELEMENTS(rad_arr)-1 DO BEGIN
                mark = WHERE(el_map GE rad_arr[ii]-width AND el_map LE rad_arr[ii]+width)
                rad_lim_map[mark] = 1.
        ENDFOR
        map_mark_name = 'map_marks.fits'
        writefits, map_mark_name, rad_lim_map, hd
        PRINT, 'WRITTEN: '+map_mark_name
ENDIF

nsteps = 250	; number of steps used in the photometry
stepsize_kpc = DBLARR(nsteps)      ; array for size of the annuli steps
stepsize_kpc[*] = 0.1   ; assigns the size of every step in [kpc]
stepsize_pc = stepsize_kpc * 1E3	; stepsize in [pc]
stepsize_pix = stepsize_pc / pixsize_pc	; stepsize in [pixels]

IF print_text EQ 1 THEN PRINT, 'Step size: '+STRTRIM(MEAN(stepsize_pc),1)+' pc'
IF print_text EQ 1 THEN PRINT, 'Total stepsize: '+STRTRIM(TOTAL(stepsize_kpc),1)+' kpc'
IF print_text EQ 1 THEN PRINT, ''

av_int_rad_irr1 = DBLARR(nsteps)
int_rad_irr1 = DBLARR(nsteps)
av_int_rad_irr2 = DBLARR(nsteps)
int_rad_irr2 = DBLARR(nsteps)
av_int_rad_irr3 = DBLARR(nsteps)
int_rad_irr3 = DBLARR(nsteps)
av_int_rad_irr4 = DBLARR(nsteps)
int_rad_irr4 = DBLARR(nsteps)
av_int_rad_irr5 = DBLARR(nsteps)
int_rad_irr5 = DBLARR(nsteps)
av_int_rad_irr6 = DBLARR(nsteps)
int_rad_irr6 = DBLARR(nsteps)
av_int_rad_tot = DBLARR(nsteps)
int_rad_tot = DBLARR(nsteps)
av_int_rad_HII = DBLARR(nsteps)
int_rad_HII = DBLARR(nsteps)
av_int_rad_HII4 = DBLARR(nsteps)
int_rad_HII4 = DBLARR(nsteps)
av_int_rad_HII6 = DBLARR(nsteps)
int_rad_HII6 = DBLARR(nsteps)
av_int_rad_HII7 = DBLARR(nsteps)
int_rad_HII7 = DBLARR(nsteps)

x_rad_pix = DBLARR(nsteps)	; [pix]
x_rad_model = DBLARR(nsteps)	; currently in [pix] later to be converted to [kpc]
area_model = dblarr(nsteps)	; array to contain the surface area corresponding to radius
x_rad_pix[0] = 0
x_rad_model[0] = 0

FOR j = 1, nsteps-1 DO BEGIN
	stepdiff = stepsize_pix[j]-stepsize_pix[j-1]
	x_rad_pix[j] = x_rad_pix[j-1] + stepsize_pix[j]	;defines array of pixels in units of stepsize
	x_rad_model[j] = x_rad_pix[j-1] + stepdiff	; later to be converted from [pixels] to [kpc] (element j is at the midpoint of the annulus)
	
	iq = WHERE((el_map GE x_rad_pix[j-1]) AND (el_map LT x_rad_pix[j]), count)	;pixels that lie within the annulus
	IF iq[0] LT 0 THEN CONTINUE
	av_int_rad_irr1[j] = (TOTAL(irr1_data[iq])) / count
	int_rad_irr1[j] = int_rad_irr1[j-1] + (av_int_rad_irr1[j] * count)

        av_int_rad_irr2[j] = (TOTAL(irr2_data[iq])) / count
        int_rad_irr2[j] = int_rad_irr2[j-1] + (av_int_rad_irr2[j] * count)

        av_int_rad_irr3[j] = (TOTAL(irr3_data[iq])) / count
        int_rad_irr3[j] = int_rad_irr3[j-1] + (av_int_rad_irr3[j] * count)

        av_int_rad_irr4[j] = (TOTAL(irr4_data[iq])) / count
        int_rad_irr4[j] = int_rad_irr4[j-1] + (av_int_rad_irr4[j] * count)

        av_int_rad_irr5[j] = (TOTAL(irr5_data[iq])) / count
        int_rad_irr5[j] = int_rad_irr5[j-1] + (av_int_rad_irr5[j] * count)

        av_int_rad_irr6[j] = (TOTAL(irr6_data[iq])) / count
        int_rad_irr6[j] = int_rad_irr6[j-1] + (av_int_rad_irr6[j] * count)

        av_int_rad_tot[j] = (TOTAL(tot_data[iq])) / count
        int_rad_tot[j] = int_rad_tot[j-1] + (av_int_rad_tot[j] * count)

	av_int_rad_HII[j] = (TOTAL(HII_data[iq])) / count
	int_rad_HII[j] = int_rad_HII[j-1] + (av_int_rad_HII[j] * count)
	
	av_int_rad_HII4[j] = (TOTAL(HII4_data[iq])) / count
	int_rad_HII4[j] = int_rad_HII4[j-1] + (av_int_rad_HII4[j] * count)
	
	av_int_rad_HII6[j] = (TOTAL(HII6_data[iq])) / count
	int_rad_HII6[j] = int_rad_HII6[j-1] + (av_int_rad_HII6[j] * count)

	av_int_rad_HII7[j] = (TOTAL(HII7_data[iq])) / count
	int_rad_HII7[j] = int_rad_HII7[j-1] + (av_int_rad_HII7[j] * count)
ENDFOR

x_rad_model = x_rad_model * pixsize_pc * 1E-3	;converts form [pixels] to [kpc]
area_model[*] = (!DPI * (x_rad_model[*]^2)) / ratio     ; [kpc^2]
int_rad_tot = int_rad_tot * pix_area_sr	; convert total lum to [MJy]
av_int_rad_irr1[0] = !VALUES.F_NaN
av_int_rad_irr2[0] = !VALUES.F_NaN
av_int_rad_irr3[0] = !VALUES.F_NaN
av_int_rad_irr4[0] = !VALUES.F_NaN
av_int_rad_irr5[0] = !VALUES.F_NaN
av_int_rad_irr6[0] = !VALUES.F_NaN
av_int_rad_tot[0]  = !VALUES.F_NaN
av_int_rad_HII[0]  = !VALUES.F_NaN
av_int_rad_HII4[0] = !VALUES.F_NaN
av_int_rad_HII6[0] = !VALUES.F_NaN
av_int_rad_HII7[0] = !VALUES.F_NaN
;------------------------------------------ COMBINE NIR PROFILES -----------------------------------------------------
combine = 0	; MUST BE KEPT AT 0!! If the necessary file is found, this is automatically set to 1!
IF wavelength[i] EQ '3.6' THEN BEGIN
	model_data = root+"saves/model/"+model+"_datafile_ir36_"+scaabs+".save"
	combine = 1
ENDIF
IF wavelength[i] EQ '4.5' THEN BEGIN
	model_data = root+"saves/model/"+model+"_datafile_ir45_"+scaabs+".save"
	combine = 1
ENDIF
IF wavelength[i] EQ '5.8' THEN BEGIN
	model_data = root+"saves/model/"+model+"_datafile_ir58_"+scaabs+".save"
	combine = 1
ENDIF
IF combine NE 1 THEN GOTO, skip_combine
RESTORE, model_data
av_int_rad_tot  = av_int_rad_tot + av_int_rad_model
int_rad_tot = int_rad_tot + int_rad_model
IF print_text EQ 1 THEN PRINT, 'Combined dust and stellar components for '+wavelength[i]+'um maps.'
skip_combine:
;----------------------------------------------------------------------------------------------------------------
;----------------------------------------------------------------------------------------------------------------
limit = WHERE(x_rad_model GE 24.)
tot_lum_model = int_rad_tot[limit[0]]
wave = DOUBLE(wavelength[i]) * 1D-6     ; converts the wavelength from [um] to [m]
IF print_text EQ 1 THEN PRINT, 'tot_lum = '+STRTRIM(tot_lum_model,1)+' [MJy]'
;---------------------------------------------------------------------------------------------------------------
IF wavelength[i] EQ '3.6' THEN wave_obs = '3.507'
IF wavelength[i] EQ '4.5' THEN wave_obs = '4.437'
IF wavelength[i] EQ '5.8' THEN wave_obs = '5.739'
IF wavelength[i] EQ '8.0' THEN wave_obs = '7.927'
IF wavelength[i] EQ '24'  THEN wave_obs = '24.000'
IF wavelength[i] EQ '70'  THEN wave_obs = '70.000'
IF wavelength[i] EQ '160' THEN wave_obs = '160.000'
IF wavelength[i] EQ '250' THEN wave_obs = '250.000'
IF wavelength[i] EQ '350' THEN wave_obs = '350.000'
IF wavelength[i] EQ '500' THEN wave_obs = '500.000'
obs_data = root+'saves/obs/datafile_'+wave_obs+'um.save'
RESTORE, obs_data	;profiles from observation data
;---------------------------------- CHI-SQUARED TEST -------------------------------------
nsteps_obs = n_elements(x_rad)
chi_sqr = DBLARR(nsteps_obs)
chi_count = 0	; initialises number of items used in chi^2 calculation.
FOR k = 0, nsteps_obs-2 DO BEGIN
        xpos_obs = x_rad[k]
	if x_rad[k] GE 7. then continue
	IF k GE N_ELEMENTS(av_int_rad) THEN CONTINUE
        xpos_model = WHERE(MIN(ABS(x_rad_model-xpos_obs)) EQ ABS(x_rad_model-xpos_obs))
        chi_sqr[k] = (av_int_rad[k] - av_int_rad_tot[xpos_model[0]])^2. / (total_error[k])^2.
	chi_count = chi_count +1
ENDFOR
chi_sqr_tot = TOTAL(chi_sqr, /NaN)
chi_sqr_r = chi_sqr_tot / chi_count
PRINT, 'Reduced X^2 = '+STRTRIM(chi_sqr_r,1)
;-------------------------------- RESIDUALS -------------------------------------
if obs_data ne '' then begin
resi = dblarr(n_elements(x_rad))
for k = 0, n_elements(x_rad)-1 do begin
        xposi = where(abs(x_rad_model-x_rad[k]) eq min(abs(x_rad_model-x_rad[k])))
        resi[k] = (av_int_rad[k] - av_int_rad_tot[xposi]) / av_int_rad[k] *100
endfor
endif
;-----------------------------------------------------------------------------------------
;--------------------------- PLOTTING ----------------------------------------------------
if obs_data ne '' then begin    ; only plots if obs data is present
xmin = 0
xmax = 17.5
ymax = max(av_int_rad_tot, /NaN) *1.5
ymin = ymax * 1e-3
plotname = figdir+model+'_sb_profile_'+wavelength[i]+'_'+scaabs+'.ps'
aspect = cgpswindow()
thisdevice = !D.NAME
set_plot, 'ps', /copy
device, filename=plotname, xsize=aspect.xsize, ysize=aspect.ysize, xoffset=aspect.xoffset, yoffset=aspect.yoffset,$
        color=1, encapsulated=encapsulated, inches=aspect.inches, $
        bits=8, set_character_size=[180,200], set_font='HELVETICA', /landscape

        !p.thick=6
        !x.thick=3
        !y.thick=3
        !p.charthick=2
        !p.charsize=1.5

cgplot, x_rad, av_int_rad, linestyle=0, color='black', ytitle='I, [MJy /sr]', xrange=[xmin,xmax], yrange=[ymin,ymax],$
        /ylog, position=[0.17,0.26,0.95,0.95], xtickformat='(A1)', title=model+' '+wavelength[i]+'$\mu$m, '+scaabs
if wavelength[i] eq 3.6 or wavelength[i] eq 4.5 or wavelength[i] eq 5.8 then begin
	cgplot, x_rad_model, av_int_rad_b, linestyle=1, color='purple', /overplot
	cgplot, x_rad_model, av_int_rad_tdi, linestyle=2, color='green', /overplot
	cgplot, x_rad_model, av_int_rad_td,  linestyle=2, color='orange', /overplot
	cgplot, x_rad_model, av_int_rad_tdo, linestyle=2, color='blue', /overplot
	cgplot, x_rad_model, av_int_rad_di,  linestyle=1, color='green', /overplot
	cgplot, x_rad_model, av_int_rad_d,   linestyle=1, color='orange', /overplot
	cgplot, x_rad_model, av_int_rad_do,  linestyle=1, color='blue', /overplot
endif
cgplot, x_rad_model, av_int_rad_irr1, linestyle=1, color='grey', /overplot
cgplot, x_rad_model, av_int_rad_irr2, linestyle=2, color='grey', /overplot
cgplot, x_rad_model, av_int_rad_irr3, linestyle=1, color='grey', /overplot
cgplot, x_rad_model, av_int_rad_irr4, linestyle=2, color='grey', /overplot
cgplot, x_rad_model, av_int_rad_irr5, linestyle=1, color='grey', /overplot
cgplot, x_rad_model, av_int_rad_irr6, linestyle=2, color='grey', /overplot
cgplot, x_rad_model, av_int_rad_HII,  linestyle=1, color='cyan', /overplot
cgplot, x_rad_model, av_int_rad_HII4, linestyle=1, color='cyan', /overplot
cgplot, x_rad_model, av_int_rad_tot,  linestyle=0, color='red', /overplot

if wavelength[i] eq 3.6 or wavelength[i] eq 4.5 or wavelength[i] eq 5.8 then begin
	cgLegend, colors=['black','red','green','orange','blue','green','orange','blue','cyan'], $
        linestyle=[0,0,2,2,2,1,1,1,1], $
        title=['obs','model','inner thin','main thin','outer thin','inner thick','main thick','outer thick','HII'],$
        length=0.03, location=[0.75,0.9], vspace=4
endif else begin
	cgLegend, colors=['black','red','grey','grey','cyan'],$
	linestyle=[0,0,1,2,1], $
	title=['obs','model','Thick dust','Thin dust','HII'],$
	length=0.03, location=[0.75,0.9], vspace=4
endelse

cgplot, [xmin,xmax], [20,20], color="blue", linestyle=2, position=[0.17,0.1,0.95,0.25], $
        ytitle="R [%]", xrange=[xmin,xmax], xtitle='Radius, [kpc]', /NoErase, yrange=[-50,50], yticks=1
cgplot, [xmin,xmax], [-20,-20], color="blue", linestyle=2, /overplot
cgplot, [xmin,xmax], [0,0], color="black", /overplot
cgplot, x_rad, resi, color='red', /overplot

device, /close_file
cgfixps, plotname
print, 'SAVED: '+plotname
endif
;-----------------------------------------------------------------------------------------
;-----------------------------------------------------------------------------------------
PRINT, 'Max intensity profile value - obs: ', STRTRIM(MAX((int_rad*1e6),/NaN),1), ' [Jy]'
print, '                            model: ', STRTRIM(MAX(int_rad_tot*1e6,/NaN),1), ' [Jy]'

save_name = savedir+model+'_datafile_'+wavelength[i]+'um_'+scaabs+'.save'
IF wavelength[i] EQ '3.6' OR wavelength[i] EQ '4.5' OR wavelength[i] EQ '5.8' THEN BEGIN
	SAVE, x_rad_model, av_int_rad_irr1, av_int_rad_irr2, av_int_rad_irr3, av_int_rad_irr4, av_int_rad_irr5, $
	av_int_rad_irr6, av_int_rad_tot, int_rad_tot, $
	av_int_rad_HII, av_int_rad_HII4, av_int_rad_HII6, av_int_rad_HII7, $
	av_int_rad_b, av_int_rad_d, av_int_rad_td, av_int_rad_di, av_int_rad_tdi, av_int_rad_do, av_int_rad_tdo, $
	chi_sqr_r, resi, chi_sqr, chi_count, FILENAME=save_name
ENDIF ELSE BEGIN
	SAVE, x_rad_model, av_int_rad_irr1, av_int_rad_irr2, av_int_rad_irr3, av_int_rad_irr4, av_int_rad_irr5, $
	av_int_rad_irr6, av_int_rad_tot, int_rad_tot, $
	av_int_rad_HII, av_int_rad_HII4, av_int_rad_HII6, av_int_rad_HII7, $
	chi_sqr_r, resi, chi_sqr, chi_count, FILENAME=save_name
ENDELSE
PRINT, 'SAVED: '+save_name
PRINT, '-----------------------------------------------------------------------------'
ENDFOR
PRINT, '------------------------------- DONE ----------------------------------------'
stop
END
