PRO m51_fits
compile_opt idl2
;Converts the maps output form dartray_multi in the .h5 format to .fits

read_scaling, model, qyear, scaabs, tau, nsersic, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5, f_bd

root = '/nfs/d58/vlk/sedmodel/cinman/DART-Ray/'
nudir = root+'../m51a/NUrad_M51a/'
numap_dir = nudir+'out/'
griddir = root+'RUNS/GALAXY_M51a/'	; directory of .h5 files to be read
outdir = root+'MAPS/M51a/'
print_to_term = 1
print_check = 0
distance = 8.58E6

ss=''
name = nudir+'indata/geometry.in'
openr,unit,name,/get_lun
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

wavelength=STRARR(11)
wavelength[0]  = '3.6'   ;stellar + dust ;IRAC_3.4
wavelength[1]  = '4.5'   ;stellar + dust ;IRAC_4.5
wavelength[2]  = '5.8'   ;stellar + dust
wavelength[3]  = '8.0'   ;dust
wavelength[4]  = '24'  ;dust
wavelength[5]  = '70'  ;dust
wavelength[6]  = '100' ;dust
wavelength[7]  = '160' ;dust
wavelength[8]  = '250' ;dust
wavelength[9]  = '350' ;dust
wavelength[10] = '500' ;dust

grid = '2DTO3D_GRIDS/M51a_GRIDS/grid_m2dto3d_m51a_tot_main.h5'
grid_data = H5F_OPEN(grid)
csize = H5D_READ(H5D_OPEN(grid_data, 'csize'))
H5F_CLOSE, grid_data

;morph = ['tot','irr1','irr2','irr3','irr4','irr5','irr6']	; different morphologies
morph = ['tot','irr1','irr2','irr3','irr4','irr5','irr6','irr_HII','irr_HII4','irr_HII6','irr_HII7']
FOR i = 0, N_ELEMENTS(morph)-1 DO BEGIN
	map = griddir+'grid_m2dto3d_m51a_'+morph[i]+'_maps_part2_int.h5'
	IF FILE_TEST(map) EQ 0 THEN BEGIN
		PRINT, map+' not found! Skipping this component.'
		GOTO, skip
	ENDIF ELSE PRINT, 'READ: '+map
	data_map = H5F_OPEN(map)

	data = H5D_READ(H5D_OPEN(data_map, 'map_arr_out'))
	lambda = H5D_READ(H5D_OPEN(data_map, 'lambda_arr_maps'))
	map_size = H5D_READ(H5D_OPEN(data_map, 'size_map'))
	H5F_CLOSE, data_map
	
	nx = N_ELEMENTS(data[*,0,0])
	ny = N_ELEMENTS(data[0,*,0])
	xc = nx/2
	yc = ny/2
	pixsize_pc = map_size / nx
	pixsize_deg = FLOAT(ATAN(pixsize_pc/distance) *(180./!DPI))
	pixsize_arcsec = pixsize_deg * 3600.
	map = DBLARR(nx,ny)
	IF print_to_term EQ 1 AND print_check EQ 0 THEN BEGIN
		PRINT, 'Size of map read: '+STRTRIM(map_size/1E3,1)+' x '+STRTRIM(map_size/1E3,1)+' [kpc]'
		PRINT, 'Pixel size of map: '+STRTRIM(pixsize_pc,1)+' [pc/pixel], '+STRTRIM(pixsize_arcsec,1)+' [arcsec/pixel]'
		PRINT, 'Resolution of model grid: '+STRTRIM(MIN(csize),1)+' [pc]'
		PRINT, ''
		print_check = 1
	ENDIF

	FOR k = 0, N_ELEMENTS(wavelength)-1 DO BEGIN
		map[*,*] = data[*,*,k]
		outputname = outdir+model+'_map_'+morph[i]+'_'+wavelength[k]+'um_'+scaabs+'.fits'      ;directory and name of output
;outputname = '../../M51_lmc_abs/maps_dust/map_'+model+'_q'+qyear+'_t'+stau+'_s'+strsfr+'_no'+strold+'_bd'+strbd+'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+'_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+'_reff'+sreff+'_ell'+sellipt+'_'+scaabs+'_l'+wavelength[k]+'_'+morph[i]+'.fits'
;if morph[i] eq 'tot' then outputname = '../../M51_lmc_abs/maps_dust_total/map_'+model+'_q'+qyear+'_t'+stau+'_s'+strsfr+'_no'+strold+'_bd'+strbd+'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+'_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+'_reff'+sreff+'_ell'+sellipt+'_'+scaabs+'_l'+wavelength[k]+'_'+morph[i]+'.fits'
		FXADDPAR, header, 'SIMPLE', 'T'
		FXADDPAR, header, 'BITPIX', -32
		FXADDPAR, header, 'NAXIS', SIZE(map, /N_DIMENSIONS) 
		FXADDPAR, header, 'NAXIS1', nx
		FXADDPAR, header, 'NAXIS2', ny
		FXADDPAR, header, 'CDELT1', pixsize_deg[0]
		FXADDPAR, header, 'CDELT2', pixsize_deg[0]
	        writefits, outputname, map, header      ;writes the fits file
	        PRINT, 'WRITTEN: ' + outputname
	ENDFOR
	skip:
	PRINT, ''
ENDFOR
END
