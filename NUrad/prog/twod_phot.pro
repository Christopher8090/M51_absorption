PRO twod_phot,model,qyear,scaabs,tau,sfr,sfr4,sfr6,sfr7,old,old3,old5,bd,f_uv,f_uv4,f_uv6,f_uv7,f_BVIK,f_BVIK3,f_BVIK5
compile_opt idl2
; The purpose of this routine is to produce the surface brightness profiles of the UV/optical output of NUrad for
; M51a. 
; The required inputs include:
;	- geometry.in
;	- gal_param.in
;	- scaling.in
;	- ff_scaling.in
;
; Looping through the different wavelengths, the maps for each morphological component are resized such that they
; spatially coincide according to 'geometry.in'. The maps are convolved if necessary, then scaled according to 
; 'scaling.in' and 'ff_scaling.in'. Then, the radial SB profiles are calculated, and the corresponding residuals
; and chi^2 is calculated. Lastly, a plot is produced and a '.save' file is saved.
;	
start
read_scaling, model, qyear, scaabs, tau, nsersic, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5, f_bd

print_text = 0	; to print information about each step, set to 1.
root = '../../'
nurad_dir = root+'NUrad/'
map_dir = nurad_dir+'out/'
save_dir = root+'saves/model/'
figdir = root+'figures/'
emiss_dir = root+'/emission_NUrad/indata/'
kernel_dir = root+'maps/obs/kernels/'

ss=''
fname = nurad_dir+'indata/geometry.in'
openr, unit, fname, /get_lun
readf, unit, ss
readf, unit, tau1
readf, unit, ss
readf, unit, tau2
readf, unit, ss
readf, unit, hd
readf, unit, ss
readf, unit, zd
readf, unit, ss
readf, unit, hdin
readf, unit, ss
readf, unit, zdin
readf, unit, ss
readf, unit, hdsolar
readf, unit, ss
readf, unit, zdsolar
readf, unit, ss
readf, unit, hd1
readf, unit, ss
readf, unit, zd1
readf, unit, ss
readf, unit, hd1in
readf, unit, ss
readf, unit, zd1in
readf, unit, ss
readf, unit, hd1solar
readf, unit, ss
readf, unit, zd1solar
readf, unit, ss
readf, unit, h_bdisc
readf, unit, ss
readf, unit, h_vdisc
readf, unit, ss
readf, unit, h_idisc
readf, unit, ss
readf, unit, h_jdisc
readf, unit, ss
readf, unit, h_kdisc
readf, unit, ss
readf, unit, h_ir34disc
readf, unit, ss
readf, unit, h_ir45disc
readf, unit, ss
readf, unit, h_ir58disc
readf, unit, ss
readf, unit, zs
readf, unit, ss
readf, unit, hsin
readf, unit, ss
readf, unit, zsin
readf, unit, ss
readf, unit, hssolar
readf, unit, ss
readf, unit, zssolar
readf, unit, ss
readf, unit, hs1
readf, unit, ss
readf, unit, zs1
readf, unit, ss
readf, unit, hs1in
readf, unit, ss
readf, unit, zs1in
readf, unit, ss
readf, unit, hs1solar
readf, unit, ss
readf, unit, zs1solar
readf, unit, ss
readf, unit, rtrun
readf, unit, ss
readf, unit, sharp
readf, unit, ss
readf, unit, rtrund
readf, unit, ss
readf, unit, sharpd
readf, unit, ss
readf, unit, rtrun1
readf, unit, ss
readf, unit, sharp1
readf, unit, ss
readf, unit, rtrund1
readf, unit, ss
readf, unit, sharpd1
readf, unit, ss
readf, unit, reff
readf, unit, ss
readf, unit, ellipt
readf, unit, ss
readf, unit, nsersic
readf, unit, ss
readf, unit, xis0
readf, unit, ss
readf, unit, xis1
readf, unit, ss
readf, unit, xid0
readf, unit, ss
readf, unit, xid1
readf, unit, ss
readf, unit, idisc1
readf, unit, ss
readf, unit, idisc2
FOR skip = 0,190 DO BEGIN
	readf, unit, ss
ENDFOR
readf, unit, ss
readf, unit, hstin
readf, unit, ss
readf, unit, hdtin
readf, unit, ss
readf, unit, hs1tin
readf, unit, ss
readf, unit, hd1tin
readf, unit, ss
readf, unit, hs3tin
readf, unit, ss
readf, unit, hd3tin
readf, unit, ss
readf, unit, hs4tin
readf, unit, ss
readf, unit, hd4tin
readf, unit, ss
readf, unit, hs5tin
readf, unit, ss
readf, unit, hd5tin
readf, unit, ss
readf, unit, hs6tin
readf, unit, ss
readf, unit, hd6tin

free_lun, unit
hs = h_bdisc
shd = strcompress(string(fix(hd*1000.)),/remove_all)
szd = strcompress(string(fix(zd*1000.)),/remove_all)
shd1 = strcompress(string(fix(hd1*1000.)),/remove_all)
szd1 = strcompress(string(fix(zd1*1000.)),/remove_all)
shs = strcompress(string(fix(hs*1000.)),/remove_all)
szs = strcompress(string(fix(zs*1000.)),/remove_all)
shs1 = strcompress(string(fix(hs1*1000.)),/remove_all)
szs1 = strcompress(string(fix(zs1*1000.)),/remove_all)
sreff = strcompress(string(fix(reff*1000.)),/remove_all)
sellipt = strcompress(string(fix(ellipt*100.)),/remove_all)
strsfr = strcompress(string(fix(sfr*100.)),/remove_all)
strold = strcompress(string(fix(old*100.)),/remove_all)
strbd = strcompress(string(fix(bd*100.)),/remove_all)
stau = strcompress(string(fix(tau*10.)),/remove_all)
snsersic = strcompress(string(fix(nsersic)),/remove_all)

fname = nurad_dir+'indata/gal_param.in'
openr, unit, fname, /get_lun
readf, unit, ss
readf, unit, ss
readf, unit, ss
readf, unit, ss
readf, unit, ss
readf, unit, ss
readf, unit, ss
readf, unit, ss
readf, unit, ss
readf, unit, dist_in
readf, unit, ss
readf, unit, inclination
readf, unit, ss
readf, unit, ss
readf, unit, pix_size
readf, unit, ss
readf, unit, nx_b
readf, unit, ss
readf, unit, ny_b
readf, unit, ss
readf, unit, nx_n
readf, unit, ss
readf, unit, ny_n
readf, unit, ss
readf, unit, nx_i
readf, unit, ss
readf, unit, ny_i
readf, unit, ss
readf, unit, nx_m
readf, unit, ss
readf, unit, ny_m
readf, unit, ss
readf, unit, nx_o
readf, unit, ss
readf, unit, ny_o
free_lun, unit

sinclination=strcompress(string(fix(inclination)),/remove_all)
wave_option=['uv09','uv13','uv15','uv16','uv20','uv22','uv25','uv28','uv36','b','v','i','j','k','ir36','ir45','ir58']
lambda=[912.,1350.,1500.,1650.,2000.,2200.,2500.,2800.,3650.,4430.,5640.,8090.,12590.,22000.,36000.,45000.,58000.]*1.d ;angstrom
distance_pc = dist_in * 1d+6  ;converts from Mpc to pc, "dist[Mpc]" in gal_param
pc_m = 3.0857D16        ; conversion between [pc] and [m]
dist_m = distance_pc * pc_m ; distance to the galaxy in [m]
flux_to_lum = (4d*!DPI*dist_m^2d) *1D-26  ; factor converting from flux to luminosity [Jy]->[W/Hz]
inc = STRCOMPRESS(STRING(FIX(ROUND(inclination))), /REMOVE_ALL)	; inclination from gal_param, truncated for later files
pixsize_model = pix_size
read, intrin, prompt="Do photometry on intrinsic maps? [1]=yes [0]=no: "

FOR i = 0, N_ELEMENTS(lambda)-1 DO BEGIN
wavelength = wave_option[i]
obs_data=''
IF lambda[i] EQ 912. THEN obs_data = ''
IF lambda[i] EQ 1350. THEN obs_data = ''
IF lambda[i] EQ 1500. THEN obs_data = '0.1542'	; i = 2
IF lambda[i] EQ 1650. THEN obs_data = ''
IF lambda[i] EQ 2000. THEN obs_data = ''
IF lambda[i] EQ 2200. THEN obs_data = '0.2274'	; i = 5
IF lambda[i] EQ 2500. THEN obs_data = ''
IF lambda[i] EQ 2800. THEN obs_data = ''
IF lambda[i] EQ 3650. THEN obs_data = '0.3562'	; i = 8
IF lambda[i] EQ 4430. THEN obs_data = '0.4719'	; i = 9
IF lambda[i] EQ 5640. THEN obs_data = '0.6185'	; i = 10
IF lambda[i] EQ 8090. THEN obs_data = '0.7500'	; i = 11
IF lambda[i] EQ 12590. THEN obs_data = '1.200'	; i = 12
IF lambda[i] EQ 22000. THEN obs_data = '2.200'	; i = 13
IF lambda[i] EQ 36000. THEN obs_data = '3.507'	; i = 14
IF lambda[i] EQ 45000. THEN obs_data = '4.437'	; i = 15
IF lambda[i] EQ 58000. THEN obs_data = '5.739'	; i = 16
if obs_data ne '0.4719' then continue
savewave = ['uv09','uv13','GALEX_FUV','uv16','uv20','GALEX_NUV','uv25','uv28','SDSS_u','SDSS_g','SDSS_r','SDSS_i','2MASS_J','2MASS_K','IRAC_3.6','IRAC_4.5','IRAC_5.8']

IF obs_data NE '' THEN BEGIN
	RESTORE, root+'saves/obs/datafile_'+obs_data+'um.save'	; restores obs data file
	PRINT, ''
	PRINT,'Obs_data: '+obs_data+' Model_data: '+wavelength
ENDIF
;-------------------------------------- SINGLE MODEL -----------------------------------------------------
if intrin eq 1 then stau='0'
read_scaling, model, qyear, scaabs, tau, nsersic, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5, f_bd
;define filenames with model data
bulge_param = wavelength+'_'+model+'_q'+qyear+'_i'+inc+'_t'+stau+$
                '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
                '_reff'+sreff+'_ell'+sellipt+'_n'+snsersic+'_'+scaabs+'.fits'

model_param = wavelength+'_'+model+'_q'+qyear+'_i'+inc+'_t'+stau+$
                '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
                '_hs1_'+shs1+'_zs1_'+szs1+'_'+scaabs+'.fits'

old_stellar = wavelength+'_'+model+'_q'+qyear+'_i'+inc+'_t'+stau+$
                '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
                '_hs'+shs+'_zs'+szs+'_'+scaabs+'.fits'

b_model_map   = map_dir+'map_mb_'+bulge_param   ; bulge map
d_model_map   = map_dir+'map_md_'+old_stellar   ; main old stellar map
td_model_map  = map_dir+'map_mtd_'+model_param  ; main young stellar map
di_model_map  = map_dir+'map_mdi_'+old_stellar  ; inner old stellar map
tdi_model_map = map_dir+'map_mtdi_'+model_param ; inner young stellar map
do_model_map  = map_dir+'map_mdo_'+old_stellar  ; outer old stellar map
tdo_model_map = map_dir+'map_mtdo_'+model_param ; outer young stellar map
ntd_model_map = map_dir+'map_mntd_'+model_param ; nuclear thin stellar map
print, model+' maps read.'
;--------------------------------------------------------------------------------------------------------
map_nx = [nx_n, nx_i, nx_m, nx_o, nx_b]	; array of map sizes for each component in [pixels]
len = MAX(map_nx)	; the number of different components in map_nx
result_b = FILE_TEST(b_model_map)	; check if a bulge map is present
IF (result_b EQ 1) THEN BEGIN
	b_data = readfits(b_model_map,hd)	; reads the model map for this component
	b_data_rot = REVERSE(b_data)	; creates a mirror of the map half
	b_tot = DBLARR(nx_b,nx_b)	; defines a square array to fit both halves of the model map
	b_tot[0:(nx_b/2)-1,*] = b_data_rot	; sets the left side of the map to be the mirrored values
	b_tot[nx_b/2:nx_b-1,*] = b_data		; sets the right side of the map to be the un-mirrored values
	b_data = b_tot[*,*]
ENDIF

result_d = FILE_TEST(d_model_map)
IF (result_d EQ 1) THEN BEGIN
	d_data = readfits(d_model_map,hd)
	d_data_rot = REVERSE(d_data)
	d_tot = DBLARR(nx_m,nx_m)
	d_tot[0:(nx_m/2)-1,*] = d_data_rot
	d_tot[nx_m/2:nx_m-1,*] = d_data
	d_data = d_tot[*,*]
ENDIF

result_td = FILE_TEST(td_model_map)
IF (result_td EQ 1) THEN BEGIN
	td_data = readfits(td_model_map,hd)
	td_data_rot = REVERSE(td_data)
	td_tot = DBLARR(nx_m,nx_m)
	td_tot[0:(nx_m/2)-1,*] = td_data_rot
	td_tot[nx_m/2:nx_m-1,*] = td_data
	td_data = td_tot[*,*]
ENDIF

result_di = FILE_TEST(di_model_map)
IF (result_di EQ 1) THEN BEGIN
	di_data = readfits(di_model_map,hd)
	di_data_rot = REVERSE(di_data)
	di_tot = DBLARR(nx_i,nx_i)
	di_tot[0:(nx_i/2)-1,*] = di_data_rot
	di_tot[nx_i/2:nx_i-1,*] = di_data
	di_data = di_tot[*,*]
ENDIF

result_tdi = FILE_TEST(tdi_model_map)
IF (result_tdi EQ 1) THEN BEGIN
	tdi_data = readfits(tdi_model_map,hd)
	tdi_data_rot = REVERSE(tdi_data)
	tdi_tot = DBLARR(nx_i,nx_i)
	tdi_tot[0:(nx_i/2)-1,*] = tdi_data_rot
	tdi_tot[nx_i/2:nx_i-1,*] = tdi_data
	tdi_data = tdi_tot[*,*]
ENDIF

result_do = FILE_TEST(do_model_map)
IF (result_do EQ 1) THEN BEGIN
	do_data = readfits(do_model_map,hd)
	do_data_rot = REVERSE(do_data)
	do_tot = DBLARR(nx_o,nx_o)
	do_tot[0:(nx_o/2)-1,*] = do_data_rot
	do_tot[nx_o/2:nx_o-1,*] = do_data
	do_data = do_tot[*,*]
ENDIF

result_tdo = FILE_TEST(tdo_model_map)
IF (result_tdo EQ 1) THEN BEGIN
	tdo_data = readfits(tdo_model_map,hd)
	tdo_data_rot = REVERSE(tdo_data)
	tdo_tot = DBLARR(nx_o,nx_o)
	tdo_tot[0:(nx_o/2)-1,*] = tdo_data_rot
	tdo_tot[nx_o/2:nx_o-1,*] = tdo_data
	tdo_data = tdo_tot[*,*]
ENDIF

result_ntd = FILE_TEST(ntd_model_map)
IF (result_ntd EQ 1) THEN BEGIN
	ntd_data = readfits(ntd_model_map,hd)
	ntd_data_rot = REVERSE(ntd_data)
	ntd_tot = DBLARR(nx_n,nx_n)
	ntd_tot[0:(nx_n/2)-1,*] = ntd_data_rot
	ntd_tot[nx_n/2:nx_n-1,*] = ntd_data
	ntd_data = ntd_tot[*,*]
ENDIF
IF print_text EQ 1 THEN BEGIN
	PRINT, 'Main thick,   d   ', d_model_map
	PRINT, 'Main thin ,   td  ', td_model_map
	PRINT, 'Inner thick,  di  ', di_model_map
	PRINT, 'Inner thin,   tdi ', tdi_model_map
	PRINT, 'Outer thick,  do  ', do_model_map
	PRINT, 'Outer thin,   tdo ', tdo_model_map
	PRINT, 'Bulge,        b   ', b_model_map
	PRINT, 'Nuclear thin, ntd ', ntd_model_map
ENDIF
;-------------------- RESIZE COMPONENT MAPS TO BE ON THE SAME SCALE ----------------------
map_dim = MAX(map_nx)	; finds the maximum maps size in [pixels].
centre = (map_dim/2) -1	; fins the centre of the map.
b_temp = DBLARR(map_dim, map_dim)	; creates a null array of maximum map size for bulge.
ntd_temp = DBLARR(map_dim, map_dim)	; creates a null array of maximum map size for nuclear disk.
tdi_temp = DBLARR(map_dim, map_dim)	; creates a null array of maximum map size for thin inner disk.
td_temp = DBLARR(map_dim, map_dim)	; creates a null array of maximum map size for thin main disk.
di_temp = DBLARR(map_dim, map_dim)	; creates a null array of maximum map size for thick inner disk.
d_temp = DBLARR(map_dim, map_dim)	; creates a null array of maximum map size for thick main disk.
do_temp = dblarr(map_dim, map_dim)

IF nx_b LT map_dim AND result_b EQ 1 THEN BEGIN
	IF print_text EQ 1 THEN PRINT, 'Resizing bulge map to '+STRTRIM(CEIL(map_dim),1)+'x'+STRTRIM(CEIL(map_dim),1)
        FOR ii = 0, nx_b-1 DO BEGIN
                FOR jj = 0, nx_b-1 DO BEGIN
                        b_temp[centre-(nx_b/2)+ii,centre-(nx_b/2)+jj] = b_data[ii,jj]
                ENDFOR
        ENDFOR
ENDIF
IF nx_n LT map_dim AND result_ntd EQ 1 THEN BEGIN
        FOR ii = 0, nx_n-1 DO BEGIN
                FOR jj = 0, nx_n-1 DO BEGIN
                        ntd_temp[centre-(nx_n/2)+ii,centre-(nx_n/2)+jj] = ntd_data[ii,jj]
                ENDFOR
        ENDFOR
ENDIF
IF nx_i LT map_dim AND result_tdi EQ 1 THEN BEGIN
        FOR ii = 0, nx_i-1 DO BEGIN
                FOR jj = 0, nx_i-1 DO BEGIN
                        tdi_temp[centre-(nx_i/2)+ii,centre-(nx_i/2)+jj] = tdi_data[ii,jj]
			IF result_di EQ 1 THEN di_temp[centre-(nx_i/2)+ii,centre-(nx_i/2)+jj] = di_data[ii,jj]
                ENDFOR
        ENDFOR
ENDIF
IF nx_m LT map_dim AND result_td EQ 1 THEN BEGIN
        FOR ii = 0, nx_m-1 DO BEGIN
                FOR jj = 0, nx_m-1 DO BEGIN
                        td_temp[centre-(nx_m/2)+ii,centre-(nx_m/2)+jj] = td_data[ii,jj]
			IF result_d EQ 1 THEN d_temp[centre-(nx_m/2)+ii,centre-(nx_m/2)+jj] = d_data[ii,jj]
                ENDFOR
        ENDFOR
ENDIF
b_data = b_temp
ntd_data = ntd_temp
tdi_data = tdi_temp
td_data = td_temp
di_data = di_temp
d_data = d_temp
if file_test(do_model_map) eq 0 then do_data = do_temp

;cal_map_tot = d_data + td_data + do_data + tdo_data + b_data + tdi_data + di_data + ntd_data
;unit_str = wave_option[i]+'_'+model+'_q'+qyear+'_i'+inc+'_t'+stau+'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+'_hs1_'+shs1+'_zs1_'+szs1+'_'+scaabs+'.fits'
;writefits, './temp/'+unit_str, cal_map_tot

tdi_zero = WHERE(tdi_data LT 1E-11)
td_zero = WHERE(td_data LT 1E-11)
tdo_zero = WHERE(tdo_data LT 1E-11)
IF result_di EQ 1 THEN di_zero = WHERE(di_data LT 1E-11)
IF result_d EQ 1 THEN d_zero = WHERE(d_data LT 1E-11)
IF result_do EQ 1 THEN do_zero = WHERE(do_data LT 1E-11)
ntd_zero = WHERE(ntd_data LT 1E-11)
b_zero = WHERE(b_data LT 1E-11)
;-----------------------------------------------------------------------------------------
print_test_fits = 0
IF print_test_fits EQ 1 THEN BEGIN
	writefits, 'b_test.fits', b_data
	writefits, 'td_test.fits', td_data
	writefits, 'ntd_test.fits', ntd_data
	writefits, 'tdi_test.fits', tdi_data
	writefits, 'tdo_test.fits', tdo_data
	PRINT, 'Printed test maps.'
ENDIF
;-----------------------------------------------------------------------------------------
;-------------------------------- CONVOLVE MODEL MAP -------------------------------------
convolve_model = 0      ; set to 1 to convolve the model, set to 0 to skip this step.
conv_affix = ''
IF convolve_model EQ 1 THEN BEGIN
kernel = ''
IF obs_data EQ '0.1542' THEN kernel = kernel_dir+'PSF_Original_GALEX_FUV.fits'
IF obs_data EQ '0.2274' THEN kernel = kernel_dir+'PSF_Original_GALEX_NUV.fits'
IF kernel EQ '' THEN GOTO, skip_convol
IF print_text EQ 1 THEN PRINT, ''
IF print_text EQ 1 THEN PRINT, 'Convolving '+wavelength+' maps with '+kernel

conv_affix = '';'_convolved'       ;  adds this affix to some save files.
IF i LT 8 THEN kernel_image = READFITS(kernel, hdk)     ; reads in the kernel and its header.
IF i GE 8 THEN kernel_image = READFITS(kernel, hdk, /EXTEN)     ; reads in the kernel and its header.
kernel_image = CONGRID(kernel_image, 100, 100)
pixsize_kernel = ABS(fxpar(hdk,'CDELT1',count=count)*3600.)     ; pixel size of the kernel in [''].
IF count eq 0 THEN pixsize_kernel = ABS(fxpar(hdk,'CD1_1',count=count)*3600.)   ; finds alternate pixel size.

IF (ABS(pixsize_model-pixsize_kernel)/pixsize_kernel) GT 0.05 THEN BEGIN        ; checks if difference between scales is 5%.
        size_kernel = (SIZE(kernel_image))[1]
        size_new = ROUND(FLOAT(size_kernel)*pixsize_kernel/pixsize_model)
        PRINT, 'pixel size of kernel - '+STRTRIM(pixsize_kernel,1)+' model - '+STRTRIM(pixsize_model,1)
        PRINT, 'Resizing kernel from '+STRTRIM(size_kernel,1)+' to '+STRTRIM(size_new,1)
        IF (size_new MOD 2) EQ 0 THEN size_new = size_new+1     ; ensures the size of the kernel is odd.
ENDIF

kernel_image = CONGRID(kernel_image+0.0, size_new, size_new, cubic=-0.5,/center)        ; re-scales the kernel to new pixel size.
kernel_image = kernel_image/TOTAL(kernel_image) ; normalises the kernel to converve flux.
IF result_b   EQ 1 THEN b_data[*,*] = CONVOL(b_data[*,*], kernel_image)
IF result_ntd EQ 1 THEN ntd_data[*,*] = CONVOL(ntd_data[*,*], kernel_image)
IF result_tdi EQ 1 THEN tdi_data[*,*] = CONVOL(tdi_data[*,*], kernel_image)
IF result_td  EQ 1 THEN td_data[*,*] = CONVOL(td_data[*,*], kernel_image)
IF result_tdo EQ 1 THEN tdo_data[*,*] = CONVOL(tdo_data[*,*], kernel_image)
IF result_di  EQ 1 THEN di_data[*,*] = CONVOL(di_data[*,*], kernel_image)
IF result_d   EQ 1 THEN d_data[*,*] = CONVOL(d_data[*,*], kernel_image)
IF result_do  EQ 1 THEN do_data[*,*] = CONVOL(do_data[*,*], kernel_image)
IF print_text EQ 1 THEN PRINT, 'All components successfully convolved!'
IF print_text EQ 1 THEN PRINT, ''
ENDIF
skip_convol:
;-----------------------------------------------------------------------------------------
;-------------------------------- SCALE TEMPLATE SED -------------------------------------
;apply old, sfr, ffactor, and bulge/disc ratio taken from radiation_fields_template.pro
;define a wavelength dependence f factor for the young uv disc
filter_opt=['uv36','b','v','i','j','k','ir36','ir45','ir58']
dim_opt = N_ELEMENTS(filter_opt)
filter_uv = ['uv09','uv13','uv15','uv16','uv20','uv22','uv25','uv28']
dim_uv = N_ELEMENTS(filter_uv)
filter = ['uv09','uv13','uv15','uv16','uv20','uv22','uv25','uv28','uv36','b','v','i','j','k','ir36','ir45','ir58']
dim = N_ELEMENTS(filter)

;define a wavelength dependence f factor for the young uv disk
f=[0.573,0.516,0.473,0.455,0.372,0.324,0.261,0.206,0.108,0.068,0.038,0.016,0.009,0.006,0.001,0.001,0.001] ;Fcal=0.35
f_diff = 1.- f	;f diffuse

;standard template
corr_old=1.
corr_old3=1.
corr_old5=1.
lum_lunit_nu = 2.241d+36
lum_lunit_diff_nu = 1.437d+36 ;in W

ff = dblarr(dim)
ff_td = sfr * f_diff * (lum_lunit_nu/lum_lunit_diff_nu) * f_uv  ; scale factor for main young disk. eqn (26) PT11
ff_td4 = sfr4 * f_diff * (lum_lunit_nu/lum_lunit_diff_nu) * f_uv4       ; scale factor for inner young disk
ff_td6 = sfr6 * f_diff * (lum_lunit_nu/lum_lunit_diff_nu) * f_uv6       ; scale factor for young outer disk
ff_td7 = sfr7 * f_diff * (lum_lunit_nu/lum_lunit_diff_nu) * f_uv7       ; scale factor for nuclear thin disk
ff_d = ff       ; scale factor for main old disk
ff_d3 = ff      ; scale factor for inner old disk
ff_d5 = ff      ; scale factor for outer old disk
ff_b = ff       ; scale factor for bulge
ff_d[(dim-dim_opt):(dim-1)]=(old/corr_old) * f_BVIK     ; scale factor for main old disk. eqn (25) PT11
ff_d3[(dim-dim_opt):(dim-1)]=(old3/corr_old) * f_BVIK3  ; scale factor for inner old disk
ff_d5[(dim-dim_opt):(dim-1)]=(old5/corr_old) * f_BVIK5  ; scale factor for outer old disk
ff_b[(dim-dim_opt):(dim-1)]= old * f_bd * bd    ; scale factor for bulge

IF (result_b EQ 1) THEN BEGIN
	b_data = b_data[*,*] * ff_b[i]	; apply correction factor to the bulge
	b_data[b_zero] = 0
ENDIF ELSE BEGIN
	b_data = 0
ENDELSE
IF (result_d EQ 1) THEN BEGIN
	d_data = d_data[*,*] * ff_d[i]	; apply correction factor to the main old disk
	d_data[d_zero] = 0
ENDIF ELSE BEGIN
	d_data = 0
ENDELSE
IF (result_di EQ 1) THEN BEGIN
	di_data = di_data[*,*] * ff_d3[i]	; apply correction factor to the inner old disk
	di_data[di_zero] = 0
ENDIF ELSE BEGIN
	di_data = 0
ENDELSE

IF (result_do EQ 1) THEN BEGIN
	do_data = do_data[*,*] * ff_d5[i]	; apply correction factor to the outer old disk
	do_data[do_zero] = 0
ENDIF ELSE BEGIN
	do_data = 0
ENDELSE

IF (result_ntd EQ 1) THEN BEGIN
	ntd_data = ntd_data[*,*] * ff_td7[i]	; apply correction factor to the nuclear thin disk
	ntd_data[ntd_zero] = 0
ENDIF ELSE BEGIN
	ntd_data = 0
ENDELSE

td_data = td_data[*,*] * ff_td[i]	; apply correction factor to the main thin disk
tdi_data = tdi_data[*,*] * ff_td4[i]	; apply correction factor to the inner thin disk
tdo_data = tdo_data[*,*] * ff_td6[i]	; apply correction factor to the outer thin disk
tdi_data[tdi_zero] = 0
td_data[td_zero] = 0
tdo_data[tdo_zero] = 0
;-----------------------------------------------------------------------------------------

cal_map_tot = d_data + td_data + do_data + tdo_data + b_data + tdi_data + di_data + ntd_data
main_map_tot = d_data + td_data + do_data + tdo_data
inner_map_tot = b_data + tdi_data + di_data + ntd_data
;cal_name = 'calmap_'+model+'_t'+stau+'_'+wave_option[i]+'.fits'
;writefits, cal_name, cal_map_tot

unit_str = wave_option[i]+'_'+model+'_q'+qyear+'_i'+inc+'_t'+stau+'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+'_hs1_'+shs1+'_zs1_'+szs1+'_'+scaabs+'.fits'

dim = SIZE(tdo_data,/dimensions)
nx = dim[0]
ny = dim[1]

;PRINT,'Dimensions: nx:', STRTRIM(STRING(nx),1), ' ny:', STRTRIM(STRING(ny),1)
X = FINDGEN(nx) # REPLICATE(1.0, ny)  ;map of X pixel coordinates
Y = REPLICATE(1.0, nx) # FINDGEN(ny)  ;map of Y pixel coordinates

xc=(nx/2.)-1.
yc=(ny/2.)-1.

;define ellipse parameters
ratio = (cos(inclination*!DPI/180.))^(-1.)   ;ratio between major and minor axis
pos_ang = 90. ;position angle 90deg as this is ensures that the annuli match the shape of the model. Not the same as observation due to difference in orientation.     
dist_ellipse, el_map, dim, xc, yc, ratio, pos_ang;, /double  ; map of elliptical distances (major axis)

pixsize_deg = pixsize_model / 3600	; pixel size converted from [''] to [degrees]
pixsize_pc = distance_pc * TAN(!DPI/180 * pixsize_deg)     ; [pc/pixel]
pixsize_sr = (pixsize_pc/distance_pc)^2.   ; converts pixsize from [pc/pixel] to [sr/pix] same as [deg^2/pix]
mapsize = pixsize_pc * nx
;PRINT, 'pixsize_pc = '+STRTRIM(pixsize_pc,1)+' pc'
;PRINT, 'mapsize = '+STRTRIM(mapsize,1)+' pc'

stepsize_kpc = DBLARR(250)	; array for size of the annuli steps
nsteps = N_ELEMENTS(stepsize_kpc)	; the total number of steps used
stepsize_kpc[*] = 0.1	; assigns the size of every step in [kpc]
stepsize_kpc = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $	; 2kpc
		0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $	; 4kpc
		0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $	; 6
		0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $	; 8
		0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $	; 10
		0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $	; 12
		0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $	; 14
		0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $	; 16
		0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $	; 18
		0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $	; 20
		0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $	; 22
		0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $	; 24
		0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]						; 25

stepsize_kpc = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.15,0.15,0.2,0.3,0.3, $      ; 2kpc
                0.3,0.3,0.2,0.15,0.15,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $      ; 4kpc
                0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $      ; 6
                0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $      ; 8
                0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $      ; 10
                0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $      ; 12
                0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $      ; 14
                0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $      ; 16
                0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $      ; 18
                0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $      ; 20
                0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $      ; 22
                0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, $      ; 24
                0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
;print, total(stepsize_kpc)
;stop
stepsize_kpc = DBLARR(250)      ; array for size of the annuli steps
stepsize_kpc[*] = 0.1   ; assigns the size of every step in [kpc]
nsteps = N_ELEMENTS(stepsize_kpc)       ; the total number of steps used
stepsize_pc = stepsize_kpc * 1E3	; stepsize in [pc]
stepsize_pix = stepsize_pc / pixsize_pc	; stepsize in [pixels]
;PRINT, 'Total stepsize = '+STRTRIM(TOTAL(stepsize_kpc),1)+' kpc'

int_rad_model=DBLARR(nsteps) 
av_int_rad_model=DBLARR(nsteps)

int_rad_b=DBLARR(nsteps)
av_int_rad_b=DBLARR(nsteps)

int_rad_d=DBLARR(nsteps)
av_int_rad_d=DBLARR(nsteps)

int_rad_td=DBLARR(nsteps)
av_int_rad_td=DBLARR(nsteps)

int_rad_di=DBLARR(nsteps)
av_int_rad_di=DBLARR(nsteps)

int_rad_tdi=DBLARR(nsteps)
av_int_rad_tdi=DBLARR(nsteps)

int_rad_do=DBLARR(nsteps)
av_int_rad_do=DBLARR(nsteps)

int_rad_tdo=DBLARR(nsteps)
av_int_rad_tdo=DBLARR(nsteps)

int_rad_ntd=DBLARR(nsteps)
av_int_rad_ntd=DBLARR(nsteps)

x_rad_pix=DBLARR(nsteps)
x_rad_model=DBLARR(nsteps)	; array containing the radius
area_model = dblarr(nsteps)	; array to contain the surface area corresponding to radius
n_pixels = DBLARR(nsteps)	; array containing the number of pixels used to average in each annulus.

x_rad_pix[0] = 0	; [pix]
x_rad_model[0] = 0	; [pc]
FOR j = 1, nsteps-1 DO BEGIN
	step_diff = stepsize_pix[j] - stepsize_pix[j-1]	; defines the midpoint of the annuli
	x_rad_pix[j] = x_rad_pix[j-1] + stepsize_pix[j]	; defines array of pixels in units of stepsize
	x_rad_model[j] = x_rad_pix[j-1] + step_diff	; later to be converted from [pixels] to [kpc] (element j is at the midpoint of the annulus)
	
	iq = WHERE(el_map GE x_rad_pix[j-1] AND el_map LT x_rad_pix[j], count)	; pixels that lie within the annulus
	IF iq[0] EQ -1 THEN CONTINUE
	n_pixels[j] = count

	av_int_rad_b[j] = (TOTAL(b_data[iq]))/count  ; average flux within each annulus
        int_rad_b[j] = int_rad_b[j-1] + av_int_rad_b[j]*count  ; curve of growth (summed flux

	av_int_rad_td[j] = (TOTAL(td_data[iq]))/count  ; average flux within each annulus
	int_rad_td[j] = int_rad_td[j-1] + av_int_rad_td[j]*count  ; curve of growth (summed flux)
		
	av_int_rad_d[j] = (TOTAL(d_data[iq]))/count  ; average flux within each annulus
	int_rad_d[j] = int_rad_d[j-1] + av_int_rad_d[j]*count  ; curve of growth (summed flux
	
	av_int_rad_di[j] = (TOTAL(di_data[iq]))/count  ; average flux within each annulus
	int_rad_di[j] = int_rad_di[j-1] + av_int_rad_di[j]*count  ; curve of growth (summed flux
	
	av_int_rad_do[j] = (TOTAL(do_data[iq]))/count  ; average flux within each annulus
	int_rad_do[j] = int_rad_do[j-1] + av_int_rad_do[j]*count  ; curve of growth (summed flux
	
	av_int_rad_tdi[j] = (TOTAL(tdi_data[iq]))/count  ; average flux within each annulus
	int_rad_tdi[j] = int_rad_tdi[j-1] + av_int_rad_tdi[j]*count  ; curve of growth (summed flux)
	
	av_int_rad_tdo[j] = (TOTAL(tdo_data[iq]))/count  ; average flux within each annulus
	int_rad_tdo[j] = int_rad_tdo[j-1] + av_int_rad_tdo[j]*count  ; curve of growth (summed flux)
	
	av_int_rad_ntd[j] = (TOTAL(ntd_data[iq]))/count  ; average flux within each annulus
	int_rad_ntd[j] = int_rad_ntd[j-1] + av_int_rad_ntd[j]*count  ; curve of growth (summed flux)

	av_int_rad_model[j] = (total(cal_map_tot[iq]))/count
	int_rad_model[j] = int_rad_model[j-1] + total(cal_map_tot[iq])
ENDFOR
av_lum_model = av_int_rad_model * flux_to_lum	; average luminosity as a function of radius [W/Hz/kpc]
int_rad_model = int_rad_model * 1d-6	; converts int_rad to MJy
int_rad_b = int_rad_b * 1d-6
int_rad_di = int_rad_di * 1d-6
int_rad_d = int_rad_d * 1d-6
int_rad_do = int_rad_do * 1d-6
int_rad_tdi = int_rad_tdi * 1d-6
int_rad_td = int_rad_td * 1d-6
int_rad_tdo = int_rad_tdo * 1d-6

x_rad_model = x_rad_model * pixsize_pc * 1d-3	; converts form [pixels] to [kpc]
area_model[*] = (!DPI * (x_rad_model[*]^2)) / ratio	; [kpc^2]
;area of a pixel in sr
pix_area_sr_model = (pixsize_pc/distance_pc)^2.; solid_angle=Area/radius^2 UNITS: [sr]
;print, pix_area_sr_model

;remove 0 value from r=0
av_int_rad_model[0] = !VALUES.F_NaN
av_int_rad_b[0] = !VALUES.F_NaN
av_int_rad_d[0] = !VALUES.F_NaN
av_int_rad_td[0] = !VALUES.F_NaN 
av_int_rad_di[0] = !VALUES.F_NaN
av_int_rad_tdi[0] = !VALUES.F_NaN 
av_int_rad_do[0] = !VALUES.F_NaN
av_int_rad_tdo[0] = !VALUES.F_NaN 
av_int_rad_ntd[0] = !VALUES.F_NaN 

;convert model from Jy to MJy/sr
av_int_rad_model[*] = av_int_rad_model[*]/pix_area_sr_model * 1d-6
av_int_rad_b[*] = av_int_rad_b[*]/pix_area_sr_model * 1d-6
av_int_rad_d[*] = av_int_rad_d[*]/pix_area_sr_model * 1d-6
av_int_rad_td[*] = av_int_rad_td[*]/pix_area_sr_model * 1d-6
av_int_rad_di[*] = av_int_rad_di[*]/pix_area_sr_model * 1d-6
av_int_rad_tdi[*] = av_int_rad_tdi[*]/pix_area_sr_model * 1d-6
av_int_rad_do[*] = av_int_rad_do[*]/pix_area_sr_model * 1d-6
av_int_rad_tdo[*] = av_int_rad_tdo[*]/pix_area_sr_model * 1d-6
av_int_rad_ntd[*] = av_int_rad_ntd[*]/pix_area_sr_model * 1d-6
;print, 'Max obs: '+strtrim(max(int_rad,/NaN),1)+' [MJy]'
print, 'Max model: '+strtrim(max(int_rad_model,/NaN),1)+' [MJy]'

cal_map_tot = ((b_data)+(d_data)+(td_data)+(di_data)+(tdi_data)+(do_data)+(tdo_data)+(ntd_data))/pix_area_sr_model * 1d-6

sxaddpar, hd_out, 'NAXIS', 2
sxaddpar, hd_out, 'NAXIS1', nx ; changes the value 'NAXIS1' in the header to the correct value.
sxaddpar, hd_out, 'NAXIS2', ny ; changes the value 'NAXIS2' in the header to the correct value.
sxaddpar, hd_out, 'CDELT1', pixsize_deg
sxaddpar, hd_out, 'CDELT2', pixsize_deg
;sxaddpar, hd_out, 'END', 'END'
cal_map_name = map_dir+model+'_calibrated_map_'+wavelength+conv_affix+'.fits'
;writefits, cal_map_name, cal_map_tot;, hd_out
;print, 'Write: '+cal_map_name
;-----------------------------------------------------------------------------------------
if savewave[i] eq 'GALEX_NUV' or savewave[i] eq 'SDSS_g' or savewave[i] eq '2MASS_Ks' then begin
map_lim = 30d+3
map_lim_pix = map_lim / pixsize_pc
cal_map_tot_temp = dblarr(map_lim_pix, map_lim_pix)
for ii = 0, map_lim_pix-1 do begin
for jj = 0, map_lim_pix-1 do begin
        xpos = xc - (map_lim_pix/2) +ii
        ypos = yc - (map_lim_pix/2) +jj
        cal_map_tot_temp[ii,jj] = cal_map_tot[xpos,ypos]
endfor
endfor
cal_map_tot_temp = rot(cal_map_tot_temp, 78., /interp)
writefits, '../../maps/model_'+savewave[i]+'.fits', cal_map_tot_temp
print, 'written: ../../maps/model_'+savewave[i]+'.fits'
help, cal_map_tot_temp
endif
stop
;-----------------------------------------------------------------------------------------
;---------------------------------- CHI-SQUARED TEST -------------------------------------
IF obs_data EQ '' THEN BEGIN
	chi_sqr_r = 0
	chi_sqr_tot = 0
	chi_count = 0
	int_rad = 0
ENDIF ELSE BEGIN
nsteps_obs = n_elements(x_rad)
chi_sqr = DBLARR(nsteps_obs)
chi_count = 0
FOR k = 0, nsteps_obs-2 DO BEGIN
	xpos_obs = x_rad[k]
	IF xpos_obs GE 7 THEN CONTINUE
	xpos_model = WHERE(MIN(ABS(x_rad_model-xpos_obs)) EQ ABS(x_rad_model-xpos_obs))
	chi_sqr[k] = ((av_int_rad[k] - av_int_rad_model[k]) / total_error[k])^2.
	chi_count = chi_count +1
ENDFOR
chi_sqr_tot = TOTAL(chi_sqr, /NaN)
chi_sqr_r = chi_sqr_tot / chi_count
PRINT, 'Reduced X^2 = '+STRTRIM(chi_sqr_r,1)
ENDELSE
;-------------------------------- RESIDUALS -------------------------------------
if obs_data ne '' then begin
resi = dblarr(n_elements(x_rad_model))
for k = 0, n_elements(x_rad_model)-1 do begin
	xposi = where(abs(x_rad_model[k]-x_rad) eq min(abs(x_rad_model[k]-x_rad)))
	resi[k] = (av_int_rad[xposi] - av_int_rad_model[k]) / av_int_rad[xposi] *100
endfor
endif else resi = 0
;--------------------------- PLOTTING ----------------------------------------------------
if obs_data ne '' then begin	; only plots if obs data is present
xmin = 0
xmax = 17.5
ymax = max(av_int_rad_model, /NaN) *1.5
ymin = ymax * 1e-3
plotname = figdir+model+'_sb_profile_'+wavelength+'_'+scaabs+'.ps'
IF stau EQ 0 THEN plotname = figdir+model+'_sb_profile_t0_'+wavelength+'_'+scaabs+'.ps'
;plotname = figdir+'sb_profile_'+model+'_q'+qyear+'_t'+stau+'_sfr'+strsfr+'_old'+strold+'_bd'+strbd+'_'+savewave[i]+'_'+scaabs+'.ps'
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
	/ylog, position=[0.17,0.26,0.95,0.95], xtickformat='(A1)', title=model+' '+wavelength+' '+scaabs+' '+strtrim(chi_sqr_r,1)
cgplot, x_rad_model, av_int_rad_b, linestyle=1, color='purple', /overplot
cgplot, x_rad_model, av_int_rad_tdi, linestyle=2, color='green', /overplot
cgplot, x_rad_model, av_int_rad_td,  linestyle=2, color='orange', /overplot
cgplot, x_rad_model, av_int_rad_tdo, linestyle=2, color='blue', /overplot
cgplot, x_rad_model, av_int_rad_di,  linestyle=1, color='green', /overplot
cgplot, x_rad_model, av_int_rad_d,   linestyle=1, color='orange', /overplot
cgplot, x_rad_model, av_int_rad_do,  linestyle=1, color='blue', /overplot
cgplot, x_rad_model, av_int_rad_model,linestyle=0, color='red', /overplot

cgLegend, colors=['black','red','purple','green','orange','blue','green','orange','blue'], $
	linestyle=[0,0,1,2,2,2,1,1,1], $
	title=['obs','model','bulge','inner thin','main thin','outer thin','inner thick','main thick','outer thick'], $
	length=0.03, location=[0.75,0.9], vspace=4

cgplot, [xmin,xmax], [20,20], color="blue", linestyle=2, position=[0.17,0.1,0.95,0.25], $
	ytitle="R [%]", xrange=[xmin,xmax], xtitle='Radius, [kpc]', /NoErase, yrange=[-50,50], yticks=1
cgplot, [xmin,xmax], [-20,-20], color="blue", linestyle=2, /overplot
cgplot, [xmin,xmax], [0,0], color="black", /overplot
cgplot, x_rad_model, resi, color='red', /overplot

device, /close_file
cgfixps, plotname
print, 'SAVED: '+plotname
endif
;-----------------------------------------------------------------------------------------
PRINT, 'Max intensity profile value - obs: ', STRTRIM(MAX(int_rad,/NaN),1), ' [MJy]  model: ', STRTRIM(MAX(int_rad_model,/NaN),1), ' [MJy]'

save_name = save_dir+'model_'+savewave[i]+'_'+model+'_q'+qyear+'_i'+inc+'_t'+stau+'_hd'+shd+'_hd1_'+shd1+'_reff'+sreff+'_ell'+sellipt+scaabs+'_sfr'+strsfr+'_old'+strold+'_bd'+strbd+'.save'
save_name = save_dir+model+'_datafile_'+wavelength+'_'+scaabs+'.save'
IF stau EQ 0 THEN save_name = save_dir+model+'_t'+stau+'_datafile_'+wavelength+'_'+scaabs+'.save'
SAVE, x_rad_model, av_int_rad_model, int_rad_model, area_model, $
	av_int_rad_b, av_int_rad_di, av_int_rad_d, av_int_rad_do, $
	av_int_rad_ntd, av_int_rad_tdi, av_int_rad_td, av_int_rad_tdo, $
	int_rad_b, int_rad_di, int_rad_d, int_rad_do, $
	int_rad_ntd, int_rad_tdi, int_rad_td, int_rad_tdo, $
	resi, chi_sqr_r, chi_sqr, chi_count, $
	n_pixels, pixsize_sr, av_lum_model, FILENAME=save_name
PRINT, 'SAVED: '+save_name
PRINT, ''
;---------------------------------------------------------------------------------------------
ENDFOR
stop
END
