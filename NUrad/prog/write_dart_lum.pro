pro write_dart_lum
compile_opt idl2
; The purpose of this routine is to produce the apparent sed (the dust attenuated stellar sed) from the nurad output for each morphological component.

;------------------------------------------------ DEFINE DIRECTORIES ----------------------------------------
print_text = 0  ; to print information about each step, set to 1.
root = '/nfs/d58/vlk/sedmodel/cinman/'
emiss_indata = root+'m51a/emission_NUrad_M51a/outdata_intlum/'
nurad_indata = root+'m51a/NUrad_M51a/indata/'
save_dir = root+'m51a/NUrad_M51a/saves/model/'
obsdir = root+'m51a/NUrad_M51a/saves/obs/'
figdir = root+'m51a/NUrad_M51a/figures/model/'
;------------------------------------------------------------------------------------------------------------
read_scaling, model, qyear, scaabs, tau, nsersic, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5, f_bd
ss=''
name = nurad_indata+'geometry.in'
openr, unit, name, /get_lun
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
readf, unit, rtruncate
readf, unit, ss
readf, unit, sharp
readf, unit, ss
readf, unit, rtruncated
readf, unit, ss
readf, unit, sharpd
readf, unit, ss
readf, unit, rtruncate1
readf, unit, ss
readf, unit, sharp1
readf, unit, ss
readf, unit, rtruncated1
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
readf, unit, rtruncate3
readf, unit, ss
readf, unit, sharp3
readf, unit, ss
readf, unit, rtruncated3
readf, unit, ss
readf, unit, sharpd3
readf, unit, ss
readf, unit, rtruncate4
readf, unit, ss
readf, unit, sharp4
readf, unit, ss
readf, unit, rtruncated4
readf, unit, ss
readf, unit, sharpd4
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
readf, unit, rtruncate5
readf, unit, ss
readf, unit, sharp5
readf, unit, ss
readf, unit, rtruncated5
readf, unit, ss
readf, unit, sharpd5
readf, unit, ss
readf, unit, rtruncate6
readf, unit, ss
readf, unit, sharp6
readf, unit, ss
readf, unit, rtruncated6
readf, unit, ss
readf, unit, sharpd6
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
free_lun, unit

hs = h_bdisc
shd = strcompress(string(round(hd*1000.)),/remove_all)
szd = strcompress(string(round(zd*1000.)),/remove_all)
shd1 = strcompress(string(round(hd1*1000.)),/remove_all)
szd1 = strcompress(string(round(zd1*1000.)),/remove_all)
shs = strcompress(string(round(hs*1000.)),/remove_all)
szs = strcompress(string(round(zs*1000.)),/remove_all)
shs1 = strcompress(string(round(hs1*1000.)),/remove_all)
szs1 = strcompress(string(round(zs1*1000.)),/remove_all)
sreff = strcompress(string(round(reff*1000.)),/remove_all)
sellipt = strcompress(string(round(ellipt*100.)),/remove_all)
strsfr = strcompress(string(round(sfr*100.)),/remove_all)
strold = strcompress(string(round(old*100.)),/remove_all)
strbd = strcompress(string(round(bd*100.)),/remove_all)
stau = strcompress(string(round(tau*10.)),/remove_all)
snsersic = strcompress(string(round(nsersic)),/remove_all)
sffactor = strcompress(string(round(ffactor*100.)),/remove_all)
;-------------------------------------------------------------------------------------------------------------
inclination = 20.3d	; inclination of galaxy [degrees]
ratio = inclination*!DPI/180.d	;ratio between major and minor axis
pc_m = 3.0857D16	; conversion between [pc] and [m]
dist_pc = 8.58D6	; distance to the galaxy in [pc]
c = 299792458d	; speed of light [m /s]
dist_m = dist_pc * pc_m	; distance to the galaxy in [m]
flux_to_lum = (4*!DPI*dist_m^2) *1E-26	; factor converting from flux to luminosity [Jy]->[W/Hz]
limit_rad = 20.d	; radius in [kpc] out to which the luminosity is integrated
print, 'Consider galaxy where radius < '+strtrim(limit_rad,1)+' [kpc]'
;-------------------------------------------------------------------------------------------------------------
wave_option=['uv09','uv13','uv15','uv16','uv20','uv22','uv25','uv28','uv36','b','v','i','j','k','ir36','ir45','ir58']
wavelength_nurad=[912.,1350.,1500.,1650.,2000.,2200.,2500.,2800.,3650.,4430.,5640.,8090.,12590.,22000.,36000.,45000.,58000.]*1e-4	; [um]

PT11_young = [0.344,0.905,0.844,0.863,0.908,0.926,0.843,0.910,1.842,2.271,3.837,5.734,0.931,0.728,0.141,0.141,0.141]*1d21
PT11_old = [0,0,0,0,0,0,0,0,0,4.771, 9.382, 19.54, 72.20, 64.97, 12.58, 12.58, 12.58]*1d21
PT11_young = PT11_young / flux_to_lum	;converts from [W/Hz] to [Jy]
PT11_old = PT11_old / flux_to_lum	;converts from [W/Hz] to [Jy]

filter_opt=['uv36','b','v','i','j','k','ir36','ir45','ir58']
dim_opt = n_elements(filter_opt)
filter_uv = ['uv09','uv13','uv15','uv16','uv20','uv22','uv25','uv28']
dim_uv = n_elements(filter_uv)
dim = n_elements(wave_option)
opt_nan = dblarr(dim-dim_opt) * 0.
f_BVIK = [opt_nan,f_BVIK]
f_BVIK3 = [opt_nan,f_BVIK3]
f_BVIK5 = [opt_nan,f_BVIK5]
f_bd = [opt_nan,f_bd]
F_cal = 0.35	; fiducial calibration value set to 0.35 from PT11
f_lambda = [1.637,1.474,1.351,1.300,1.063,0.926,0.746,0.589,0.309,0.194,0.109,0.0457,0.0257,0.0171,0.00286,0.00286,0.00286]
	; f_lambda from Table A.1. Tuffs et al (2004)
f=[0.573,0.516,0.473,0.455,0.372,0.324,0.261,0.206,0.108,0.068,0.038,0.016,0.009,0.006,0.001,0.001,0.001]	; table A.1. Tuffs et al. (2004) multiplied by F_cal=0.35
f_diff = 1.- f	; Wavelength dependence of the fraction of photons escaping the clumpy component (Table E.4. PT11)[dimensionless].

lum_lunit_nu = 2.241d+36	; "Luminosity of the non-ionising UV photons" [W].
lum_lunit_diff_nu = 1.437d+36	; Luminosity of the non-ionising UV photons of the diffuse component [W].
model='wd01'
combined=0
;--------------------------------------------------------------------------------------------------------------
;------------------------------------ CALCULATE APPARENT SED --------------------------------------------------
morph_comp = ['tot','tdi','td','tdo','b','di','d','do']	; morphological components
bulge_i = where(morph_comp eq 'b')
inner_i = where(morph_comp eq 'tdi' or morph_comp eq 'di')
main_i = where(morph_comp eq 'td' or morph_comp eq 'd')
outer_i = where(morph_comp eq 'tdo' or morph_comp eq 'do')
young_i = where(morph_comp eq 'tdi' or morph_comp eq 'td' or morph_comp eq 'tdo')
old_i = where(morph_comp eq 'di' or morph_comp eq 'd' or morph_comp eq 'do')
tot_i = where(morph_comp eq 'tot')
tdi_i = where(morph_comp eq 'tdi')
td_i = where(morph_comp eq 'td')
tdo_i = where(morph_comp eq 'tdo')
b_i = where(morph_comp eq 'b')
di_i = where(morph_comp eq 'di')
d_i = where(morph_comp eq 'd')
do_i = where(morph_comp eq 'do')

dim_morph = n_elements(morph_comp)	; number of morphological components
app_flux = dblarr(dim_morph, dim)
if combined eq 1 then begin
for i = 0, dim-1 do begin	; apparent components.
	fname = save_dir+"wd01_datafile_"+wave_option[i]+"_"+scaabs+".save"     ; save file from NUrad photometry.
        restore, fname
        app_flux[td_i,i] = max(int_rad_td[where(x_rad_model LE limit_rad)]) * 1e6
        app_flux[tdo_i,i] = max(int_rad_tdo[where(x_rad_model LE limit_rad)]) * 1e6
        app_flux[d_i,i] = max(int_rad_d[where(x_rad_model LE limit_rad)]) * 1e6
        app_flux[do_i,i] = max(int_rad_do[where(x_rad_model LE limit_rad)]) * 1e6
        fname = save_dir+"lmc1_datafile_"+wave_option[i]+"_"+scaabs+".save"   ; save file from NUrad photometry.
        restore, fname
        app_flux[tdi_i,i] = max(int_rad_tdi[where(x_rad_model LE limit_rad)]) * 1e6
        app_flux[b_i,i] = max(int_rad_b[where(x_rad_model LE limit_rad)]) * 1e6
        app_flux[di_i,i] = max(int_rad_di[where(x_rad_model LE limit_rad)]) * 1e6
        app_flux[tot_i,i] = total(app_flux[1:7,i])      ; total of all components per wavelength
endfor
endif else begin
for i = 0, dim-1 do begin       ; apparent components looping through wavelength.
	fname = save_dir+model+"_datafile_"+wave_option[i]+"_"+scaabs+".save"   ; save file from NUrad photometry.
        restore, fname
        app_flux[tot_i,i] = max(int_rad_model[where(x_rad_model LE limit_rad)]) * 1e6
        app_flux[tdi_i,i] = max(int_rad_tdi[where(x_rad_model LE limit_rad)]) * 1e6
        app_flux[td_i,i] = max(int_rad_td[where(x_rad_model LE limit_rad)]) * 1e6
        app_flux[tdo_i,i] = max(int_rad_tdo[where(x_rad_model LE limit_rad)]) * 1e6
        app_flux[b_i,i] = max(int_rad_b[where(x_rad_model LE limit_rad)]) * 1e6
        app_flux[di_i,i] = max(int_rad_di[where(x_rad_model LE limit_rad)]) * 1e6
        app_flux[d_i,i] = max(int_rad_d[where(x_rad_model LE limit_rad)]) * 1e6
        app_flux[do_i,i] = max(int_rad_do[where(x_rad_model LE limit_rad)]) * 1e6
endfor
endelse
app_flux_tot = app_flux[tot_i,*]
app_flux_inner = app_flux[tdi_i,*] + app_flux[di_i,*]
app_flux_main = app_flux[td_i,*] + app_flux[d_i,*]
app_flux_outer = app_flux[tdo_i,*] + app_flux[do_i,*]
app_flux_young = app_flux[tdi_i,*] + app_flux[td_i,*] + app_flux[tdo_i,*]
app_flux_old = app_flux[di_i,*] + app_flux[d_i,*] + app_flux[do_i,*]
app_flux_bulge = app_flux[b_i,*]
;---------------------------------------------------------------------------------------------------------------
;--------------------------------------- CALCULATE INTRINSIC SED -----------------------------------------------
dim_morph = n_elements(morph_comp)      ; number of morphological components
int_lum = dblarr(dim_morph, dim)
if combined eq 1 then begin
for i = 0, dim-1 do begin       ; apparent components.
	fname = save_dir+"wd01_t0_datafile_"+wave_option[i]+"_sca.save"   ; save file from NUrad photometry.
	restore, fname
	int_lum[td_i,i] = max(int_rad_td[where(x_rad_model LE limit_rad)]) * 1e6* flux_to_lum
	int_lum[tdo_i,i] = max(int_rad_tdo[where(x_rad_model LE limit_rad)]) * 1e6* flux_to_lum
	int_lum[d_i,i] = max(int_rad_d[where(x_rad_model LE limit_rad)]) * 1e6* flux_to_lum
	int_lum[do_i,i] = max(int_rad_do[where(x_rad_model LE limit_rad)]) * 1e6* flux_to_lum
	fname = save_dir+"lmc1_t0_datafile_"+wave_option[i]+"_sca.save"   ; save file from NUrad photometry.
	restore, fname
	int_lum[tdi_i,i] = max(int_rad_tdi[where(x_rad_model LE limit_rad)]) * 1e6* flux_to_lum
	int_lum[b_i,i] = max(int_rad_b[where(x_rad_model LE limit_rad)]) * 1e6* flux_to_lum
	int_lum[di_i,i] = max(int_rad_di[where(x_rad_model LE limit_rad)]) * 1e6* flux_to_lum
	int_lum[tot_i,i] = total(int_lum[1:7,i])      ; total of all components per wavelength
endfor
endif else begin
for i = 0, dim-1 do begin
	fname = save_dir+model+"_t0_datafile_"+wave_option[i]+"_"+scaabs+".save" ; save file from NUrad photometry.
	restore, fname
	int_lum[tot_i,i] = max(int_rad_model[where(x_rad_model LE limit_rad)]) * 1e6* flux_to_lum
	int_lum[tdi_i,i] = max(int_rad_tdi[where(x_rad_model LE limit_rad)]) * 1e6* flux_to_lum
	int_lum[td_i,i] = max(int_rad_td[where(x_rad_model LE limit_rad)]) * 1e6* flux_to_lum
	int_lum[tdo_i,i] = max(int_rad_tdo[where(x_rad_model LE limit_rad)]) * 1e6* flux_to_lum
	int_lum[b_i,i] = max(int_rad_b[where(x_rad_model LE limit_rad)]) * 1e6* flux_to_lum
	int_lum[di_i,i] = max(int_rad_di[where(x_rad_model LE limit_rad)]) * 1e6* flux_to_lum
	int_lum[d_i,i] = max(int_rad_d[where(x_rad_model LE limit_rad)]) * 1e6* flux_to_lum
	int_lum[do_i,i] = max(int_rad_do[where(x_rad_model LE limit_rad)]) * 1e6 * flux_to_lum
endfor
endelse
dartdir = root+'DART-Ray/GALAXY_GRIDS/M51/'
lumtype = ['i_young', 'm_young', 'o_young', 'b_old', 'i_old', 'm_old', 'o_old']
scale = [sfr4, sfr, sfr6, bd*old, old3, old, old5]

for kk = 0, n_elements(lumtype)-1 do begin
	filename = dartdir+'M51'+lumtype[kk]+'_star_sed.dat'
	openw, unit, filename, /get_lun
	printf, unit, 'WAVELENGTH [um]      LUMINOSITY [W/Hz]'
	print, lumtype[kk], scale[kk]
	for nn =0, n_elements(wavelength_nurad)-1 do begin
		if int_lum[kk+1,nn] eq 0. then continue
		;printf, unit, wavelength_nurad[nn], (int_lum[kk+1,nn] *flux_to_lum)
		;printf, unit, wavelength_nurad[nn], int_lum[kk+1,nn] / scale[kk]
		print, wavelength_nurad[nn], int_lum[kk+1,nn], int_lum[kk+1,nn] / scale[kk]
	endfor
	free_lun, unit
	print, 'WRITTEN: '+filename
endfor
stop
end
