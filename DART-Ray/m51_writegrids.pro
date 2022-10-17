PRO m51_writegrids, tau, sfr, old, bd, ffactor 
compile_opt idl2

; This program takes the IDL emissivity files for the different dust
; components (calculated using the 2D RT code) and extracts the values
; only for a set of input wavelengths. Depending on the keyword
; type_grid, one can make tables of  

; INPUT
; type_grid [string]: dust emission component. It can be
; 'irr1': first dust disc
; 'irr2': second dust disc
; 'irr_HII': HII regions dust emission
; 'tot' : all dust emission components
; 
; lambda_arr [float]: array containing the wavelengths [um] for the
; output files.
;
; label_model [string]: name of the model in the output file
; file_irr1/irr2/hii [string]: names of the dust emissivities IDL
; files (within directory ./INPUT_GRIDS ). 

; OUTPUT 
; monochromatic 2D grid files within directory './INPUT_GRIDS/.  

root = '../'	; root directory.
dir_nurad = root+'NUrad/'	; NUrad directory.
dir_nurad_out = dir_nurad+'out/'	; directory containing output data from NUrad.
dir_nurad_indata = dir_nurad+'indata/'	; directory containing geometry input paramter files.
dir_emiss_intlum = root+'emission_NUrad/outdata_intlum/'     ; directory for emission data.
dir_out = './2DTO3D_GRIDS/M51a_GRIDS/'	; directory for grids to be outputted.

read_scaling, model, qyear, scaabs, tau, nsersic, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5, f_bd

ss = ''
filename_param = 'gal_param.in'
name_param = dir_nurad_indata+filename_param
openr,unit,name_param,/get_lun
print, 'read ', name_param
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
readf, unit, ss
readf, unit, ss
readf, unit, mopt_b
readf, unit, ss
readf, unit, mopt_n
readf, unit, ss
readf, unit, mopt_i
readf, unit, ss
readf, unit, mopt_m
readf, unit, ss
readf, unit, mopt_o
readf, unit, ss
readf, unit, ss
readf, unit, mstep1_b
readf, unit, ss
readf, unit, mlength1_b
readf, unit, ss
readf, unit, mstep2_b
readf, unit, ss
readf, unit, mlength2_b
readf, unit, ss
readf, unit, mstep3_b
readf, unit, ss
readf, unit, ss
readf, unit, mstep1_n
readf, unit, ss
readf, unit, mlength1_n
readf, unit, ss
readf, unit, mstep2_n
readf, unit, ss
readf, unit, mlength2_n
readf, unit, ss
readf, unit, mstep3_n
readf, unit, ss
readf, unit, ss
readf, unit, mstep1_i
readf, unit, ss
readf, unit, mlength1_i
readf, unit, ss
readf, unit, mstep2_i
readf, unit, ss
readf, unit, mlength2_i
readf, unit, ss
readf, unit, mstep3_i
readf, unit, ss
readf, unit, ss
readf, unit, mstep1_m
readf, unit, ss
readf, unit, mlength1_m
readf, unit, ss
readf, unit, mstep2_m
readf, unit, ss
readf, unit, mlength2_m
readf, unit, ss
readf, unit, mstep3_m
readf, unit, ss
readf, unit, ss
readf, unit, mstep1_o
readf, unit, ss
readf, unit, mlength1_o
readf, unit, ss
readf, unit, mstep2_o
readf, unit, ss
readf, unit, mlength2_o
readf, unit, ss
readf, unit, mstep3_o
readf, unit, ss
readf, unit, ss
readf, unit, elm_b
readf, unit, ss
readf, unit, elm_n
readf, unit, ss
readf, unit, elm_i
readf, unit, ss
readf, unit, elm_m
readf, unit, ss
readf, unit, elm_o
free_lun, unit

name = dir_nurad_indata+'/geometry.in'
openr,unit,name,/get_lun
print, 'read ', name
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
readf, unit, h_bdisk
readf, unit, ss
readf, unit, h_vdisk
readf, unit, ss
readf, unit, h_idisk
readf, unit, ss
readf, unit, h_jdisk
readf, unit, ss
readf, unit, h_kdisk
readf, unit, ss
readf, unit, h_ir36disk
readf, unit, ss
readf, unit, h_ir45disk
readf, unit, ss
readf, unit, h_ir58disk
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
readf, unit,ss
readf, unit,rtrund
readf, unit,ss
readf, unit,sharpd1
readf, unit, ss
readf, unit, rtrun1
readf, unit, ss
readf, unit, sharp1
readf, unit,ss
readf, unit,rtrund1
readf, unit,ss
readf, unit,sharpd1
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
readf, unit, idisk1
readf, unit, ss
readf, unit, idisk2
;read inner component
readf, unit, ss
readf, unit, ss
readf, unit, tau3
readf, unit, ss
readf, unit, tau4
readf, unit, ss
readf, unit, hd3
readf, unit, ss
readf, unit, zd3
readf, unit, ss
readf, unit, hd3in
readf, unit, ss
readf, unit, zd3in
readf, unit, ss
readf, unit, hd3solar
readf, unit, ss
readf, unit, zd3solar
readf, unit, ss
readf, unit, hd4
readf, unit, ss
readf, unit, zd4
readf, unit, ss
readf, unit, hd4in
readf, unit, ss
readf, unit, zd4in
readf, unit, ss
readf, unit, hd4solar
readf, unit, ss
readf, unit, zd4solar
readf, unit, ss
readf, unit, h_bdisk3
readf, unit, ss
readf, unit, h_vdisk3
readf, unit, ss
readf, unit, h_idisk3
readf, unit, ss
readf, unit, h_jdisk3
readf, unit, ss
readf, unit, h_kdisk3
readf, unit, ss
readf, unit, h_ir36disk3
readf, unit, ss
readf, unit, h_ir45disk3
readf, unit, ss
readf, unit, h_ir58disk3
readf, unit, ss
readf, unit, zs3
readf, unit, ss
readf, unit, hs3in
readf, unit, ss
readf, unit, zs3in
readf, unit, ss
readf, unit, hs3solar
readf, unit, ss
readf, unit, zs3solar
readf, unit, ss
readf, unit, hs4
readf, unit, ss
readf, unit, zs4
readf, unit, ss
readf, unit, hs4in
readf, unit, ss
readf, unit, zs4in
readf, unit, ss
readf, unit, hs4solar
readf, unit, ss
readf, unit, zs4solar
readf, unit, ss
readf, unit, rtrun3
readf, unit, ss
readf, unit, sharp3
readf, unit,ss
readf, unit,rtrund3
readf, unit,ss
readf, unit,sharpd3
readf, unit, ss
readf, unit, rtrun4
readf, unit, ss
readf, unit, sharp4
readf, unit,ss
readf, unit,rtrund4
readf, unit,ss
readf, unit,sharpd4
readf, unit, ss
readf, unit, xis3
readf, unit, ss
readf, unit, xis4
readf, unit, ss
readf, unit, xid3
readf, unit, ss
readf, unit, xid4
readf, unit, ss
readf, unit, idisk3
readf, unit, ss
readf, unit, idisk4
;read outer component
readf, unit, ss
readf, unit, ss
readf, unit, tau5
readf, unit, ss
readf, unit, tau6
readf, unit, ss
readf, unit, hd5
readf, unit, ss
readf, unit, zd5
readf, unit, ss
readf, unit, hd5in
readf, unit, ss
readf, unit, zd5in
readf, unit, ss
readf, unit, hd5solar
readf, unit, ss
readf, unit, zd5solar
readf, unit, ss
readf, unit, hd6
readf, unit, ss
readf, unit, zd6
readf, unit, ss
readf, unit, hd6in
readf, unit, ss
readf, unit, zd6in
readf, unit, ss
readf, unit, hd6solar
readf, unit, ss
readf, unit, zd6solar
readf, unit, ss
readf, unit, h_bdisk5
readf, unit, ss
readf, unit, h_vdisk5
readf, unit, ss
readf, unit, h_idisk5
readf, unit, ss
readf, unit, h_jdisk5
readf, unit, ss
readf, unit, h_kdisk5
readf, unit, ss
readf, unit, h_ir36disk5
readf, unit, ss
readf, unit, h_ir45disk5
readf, unit, ss
readf, unit, h_ir58disk5
readf, unit, ss
readf, unit, zs5
readf, unit, ss
readf, unit, hs5in
readf, unit, ss
readf, unit, zs5in
readf, unit, ss
readf, unit, hs5solar
readf, unit, ss
readf, unit, zs5solar
readf, unit, ss
readf, unit, hs6
readf, unit, ss
readf, unit, zs6
readf, unit, ss
readf, unit, hs6in
readf, unit, ss
readf, unit, zs6in
readf, unit, ss
readf, unit, hs6solar
readf, unit, ss
readf, unit, zs6solar
readf, unit, ss
readf, unit, rtrun5
readf, unit, ss
readf, unit, sharp5
readf, unit,ss
readf, unit,rtrund5
readf, unit,ss
readf, unit,sharpd5
readf, unit, ss
readf, unit, rtrun6
readf, unit, ss
readf, unit, sharp6
readf, unit,ss
readf, unit,rtrund6
readf, unit,ss
readf, unit,sharpd6
readf, unit, ss
readf, unit, xis5
readf, unit, ss
readf, unit, xis6
readf, unit, ss
readf, unit, xid5
readf, unit, ss
readf, unit, xid6
readf, unit, ss
readf, unit, idisk5
readf, unit, ss
readf, unit, idisk6
readf, unit, ss
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
;read nuclear component
readf, unit, ss
readf, unit, ss
readf, unit, hs7
readf, unit, ss
readf, unit, zs7
readf, unit, ss
readf, unit, hs7in
readf, unit, ss
readf, unit, zs7in
readf, unit, ss
readf, unit, hs7solar
readf, unit, ss
readf, unit, zs7solar
readf, unit, ss
readf, unit, rtrun7
readf, unit, ss
readf, unit, sharp7
readf, unit, ss
readf, unit, xis7
readf, unit, ss
readf, unit, idisk7
readf, unit, ss
readf, unit, hs7tin

free_lun, unit
hs  =  h_bdisk
shd = strcompress(string(round(hd*1000.)),/REMOVE_ALL)
szd = strcompress(string(round(zd*1000.)),/REMOVE_ALL)
shd1 = strcompress(string(round(hd1*1000.)),/REMOVE_ALL)
szd1 = strcompress(string(round(zd1*1000.)),/REMOVE_ALL)
shs = strcompress(string(round(hs*1000.)),/REMOVE_ALL)
szs = strcompress(string(round(zs*1000.)),/REMOVE_ALL)
shs1 = strcompress(string(round(hs1*1000.)),/REMOVE_ALL)
szs1 = strcompress(string(round(zs1*1000.)),/REMOVE_ALL)
sreff = strcompress(string(round(reff*1000.)),/REMOVE_ALL)
sellipt = strcompress(string(round(ellipt*100.)),/REMOVE_ALL)
strsfr = strcompress(string(round(sfr*100.)),/REMOVE_ALL)
strold = strcompress(string(round(old*100.)),/REMOVE_ALL)
strbd = strcompress(string(round(bd*100.)),/REMOVE_ALL)
stau = strcompress(string(round(tau*10.)),/REMOVE_ALL)
sffactor = strcompress(string(round(ffactor*100.)),/REMOVE_ALL)
sinclination = strcompress(string(round(inclination)),/REMOVE_ALL)

snsersic = strcompress(string(round(nsersic)),/REMOVE_ALL)
;------------------------------------ FILE NAMES IN EMISSION/OUT_LUM --------------------------------------------
type_grid_arr = ['tot']
label_model = model+'_q'+qyear+ $
   '_t'+stau+'_s'+strsfr+'_no'+strold+'_bd'+strbd+'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
   '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+'_reff'+sreff+'_ell'+sellipt+'_'+scaabs
print,'label_model = ',label_model
file_irr1 = 'grid_irr1_'+label_model+'.xdr'
if file_test(dir_emiss_intlum+file_irr1) eq 1 then begin 
	type_grid_arr = [type_grid_arr, 'irr1']
	print,file_irr1
endif
file_irr2 = 'grid_irr2_'+label_model+'.xdr'
if file_test(dir_emiss_intlum+file_irr2) then begin
	type_grid_arr = [type_grid_arr, 'irr2']
	print,file_irr2
endif
file_irr3 = 'grid_irr3_'+label_model+'.xdr'
if file_test(dir_emiss_intlum+file_irr3) then begin
	type_grid_arr = [type_grid_arr, 'irr3']
	print,file_irr3
endif
file_irr4 = 'grid_irr4_'+label_model+'.xdr'
if file_test(dir_emiss_intlum+file_irr4) then begin
	type_grid_arr = [type_grid_arr, 'irr4']
	print,file_irr4
endif
file_irr5 = 'grid_irr5_'+label_model+'.xdr'
if file_test(dir_emiss_intlum+file_irr5) then begin
	type_grid_arr = [type_grid_arr, 'irr5']
	print,file_irr5
endif
file_irr6 = 'grid_irr6_'+label_model+'.xdr'
if file_test(dir_emiss_intlum+file_irr6) then begin
	type_grid_arr = [type_grid_arr, 'irr6']
	print,file_irr6
endif
file_hii = 'grid_irr_HII_'+model+'_q'+qyear+ '_t'+stau+'_s'+strsfr+'_no'+strold+'_bd'+strbd+$
	   '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+'_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
	   '_reff'+sreff+'_ell'+sellipt+'_f'+sffactor+'.xdr'
if file_test(dir_emiss_intlum+file_hii) then begin
	type_grid_arr = [type_grid_arr, 'irr_HII', 'irr_HII4', 'irr_HII6', 'irr_HII7']
	print,file_hii
endif
print, type_grid_arr
;------------------------------------------------------------------------------------------------------------
for k = 0, n_elements(type_grid_arr)-1 do begin
type_grid = type_grid_arr[k]

iq = where(type_grid EQ type_grid_arr, count)
if (count ne 1) then print, 'type_grid not recognized'

lambda_arr = [3.6, 4.5, 5.8, 8., 24., 70., 100., 160., 250., 350., 500.]	; wavelengths for grids to be produced at.

if (idisk1 EQ 1) then begin
	restore, dir_emiss_intlum+file_irr1
	lum_diff1 = lumdouble1
endif ELSE begin
        lum_diff1 = 0
endELSE

if (idisk2 EQ 1) then begin
        restore, dir_emiss_intlum+file_irr2
        lum_diff2 = lumdouble2
endif ELSE begin
        lum_diff2 = 0
endELSE

if (idisk3 EQ 1) then begin
        restore, dir_emiss_intlum+file_irr3
        lum_diff3 = lumdouble3
endif ELSE begin
        lum_diff3 = 0
endELSE

if (idisk4 EQ 1) then begin
        restore, dir_emiss_intlum+file_irr4
        lum_diff4 = lumdouble4
endif ELSE begin
        lum_diff4 = 0
endELSE

if (idisk5 EQ 1) then begin
        restore, dir_emiss_intlum+file_irr5
        lum_diff5 = lumdouble5
endif ELSE begin
        lum_diff5 = 0
endELSE

if (idisk6 EQ 1) then begin
        restore, dir_emiss_intlum+file_irr6
        lum_diff6 = lumdouble6
endif ELSE begin
        lum_diff6 = 0
endELSE

r_arr = rrr	; array of radial positions.
z_arr = zzz	; array of height positions.

restore, dir_emiss_intlum+file_hii	; restores the file for all HII (clumpy) components.
lum_hii = eta
lum_hii4 = eta4
lum_hii6 = eta6
lum_hii7 = eta7
if type_grid eq 'irr_HII4' and max(lum_hii4) eq 0. then continue
if type_grid eq 'irr_HII6' and max(lum_hii6) eq 0. then continue
if type_grid eq 'irr_HII7' and max(lum_hii7) eq 0. then continue
lambda = lambda *1d-4

print, type_grid, count
if (type_grid ne 'tot') then begin	; zeroes every other component if the total isn't being calculated.
	if (type_grid ne 'irr_HII') then lum_hii = 0
	if (type_grid ne 'irr_HII4') then lum_hii4 = 0
	if (type_grid ne 'irr_HII6') then lum_hii6 = 0
	if (type_grid ne 'irr_HII7') then lum_hii7 = 0
	if (type_grid ne 'irr1') then lum_diff1 = 0
	if (type_grid ne 'irr2') then lum_diff2 = 0
	if (type_grid ne 'irr3') then lum_diff3 = 0
	if (type_grid ne 'irr4') then lum_diff4 = 0
	if (type_grid ne 'irr5') then lum_diff5 = 0
	if (type_grid ne 'irr6') then lum_diff6 = 0
endif

nr = n_elements(r_arr[*,0])	; number of elements in the radial sampling.
nz = n_elements(r_arr[0,*])	; number of elements in the height sampling.
n_lambda = n_elements(lambda)	; number of wavelengths emission is calculated for.

lum_tot = lum_diff1 + lum_diff2 + lum_diff3 + lum_diff4 + lum_diff5 + lum_diff6 +$
		lum_hii + lum_hii4 + lum_hii6 + lum_hii7
lum_tot_ext = dblarr(nr,nz)

for il = 0, n_elements(lambda_arr)-1 do begin	; loop through each wavelength in lambda_arr
	lambda_grid = lambda_arr[il]   ; [microns]
	label = strcompress(string(lambda_grid, format = ("(F9.3)")), /remove_all)
	file_2d_grid = 'grid_'+label_model+'_'+type_grid+'_l'+label+'um.dat'	; name of the output grid.
	
	;lambda_grid = lambda_grid*10d^4	; converts wavelength input from [um] to [A]
	lum_tot_ext[*,*] = 0.	; intialises the array.
	
	; extract wavelength of interest and interpolate fluxes
	iq0 = value_locate(lambda, lambda_grid)	; finds where lambda is equal to lambda_grid
	iq1 = iq0+1
	if (iq0 LT 0 OR iq1 GT n_lambda) then begin
		print, 'wavelength outside value range'
		stop
	endif
	
	for i = 0, nr-1 do begin	; loops through radial positions
		for j = 0, nz-1 do begin	; loops through height positions
			arr = lum_tot[i,j,*]	; array of luminosities at a given R, z position for all wavelengths.
			const = interpol(arr, lambda, lambda_grid)	; interpolates the value of arr at the specific wavelength lambda_grid.
			lum_tot_ext[i,j] = const	; assigns this interpolated value to the respective R, z position.
		endfor
	endfor
	openw, unit10, dir_out+file_2d_grid, /GET_LUN
	printf, unit10, ' R [pc]      Z [pc]     j_nu   [W/Hz/pc^3]    krho  [pc^-1]'
	for i = 0,nr -1 do begin
		for j = 0,nz -1 do begin
			printf, unit10, r_arr[i,j], z_arr[i,j], lum_tot_ext[i,j], 0.
		endfor
	endfor
	free_lun, unit10
endfor
endfor
return
end
