pro sed_model, model, qyear, tau, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, $
	ffactor, scaabs, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7,$
	 f_BVIK, f_BVIK3, f_BVIK5, f_bd
compile_opt idl2
; The purpose of this routine is to produce the apparent sed (the dust attenuated stellar sed) from the nurad output for each morphological component.

read_scaling, model, qyear, scaabs, tau, nsersic, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5, f_bd
;-------------------------------------- DEFINE DIRECTORIES ------------------------------
print_text = 0  ; to print information about each step, set to 1.
root = '../../'
nurad_indata = root+'NUrad/indata/'
save_dir = root+'saves/model/'
obsdir = root+'saves/obs/'
figdir = root+'figures/'
results_dir = root+'results/'
;-----------------------------------------------------------------------------------------
;filename = nurad_indata+'ff_scaling.in'
;READCOL, filename, wave, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5, f_bd, SKIPLINE=1, FORMAT="(A,D,D,D,D,D,D,D,D)"

;f_BVIK = f_BVIK[WHERE(f_BVIK NE 0)]
;f_BVIK3 = f_BVIK3[WHERE(f_BVIK5 NE 0)]
;f_BVIK5 = f_BVIK5[WHERE(f_BVIK5 NE 0)]
;f_bd = f_bd[WHERE(f_bd NE 0)]

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
limit_rad = 24.d	; radius in [kpc] out to which the luminosity is integrated
print, 'Consider galaxy where radius < '+strtrim(limit_rad,1)+' [kpc]'
;-------------------------------------------------------------------------------------------------------------
wave_option=['uv09','uv13','uv15','uv16','uv20','uv22','uv25','uv28','uv36','b','v','i','j','k','ir36','ir45','ir58']
wavelength_nurad=[912.,1350.,1500.,1650.,2000.,2200.,2500.,2800.,3650.,4430.,5640.,8090.,12590.,22000.,36000.,45000.,58000.]*1e-4	; [um]
wavelengths_m = wavelength_nurad * 1e-6 ; [m]
sf_wavelengths = [912.,1350.,1500.,1650.,2000.,2200.,2500.,2800.,3650.,4430.]	; wavelengths for SFR calculation

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
int_flux = dblarr(dim_morph, dim)
int_lum = dblarr(dim_morph, dim)
openw, t0_lun, results_dir+'lum_intrinsic_'+model+'.dat', /get_lun
printf, t0_lun, 'Band	tdi[W/Hz]  td		  tdo	       di           d            do           b'
for i = 0, dim-1 do begin
	fname = save_dir+model+"_t0_datafile_"+wave_option[i]+"_"+scaabs+".save" ; save file from NUrad photometry.
	restore, fname
	int_flux[tot_i,i] = max(int_rad_model[where(x_rad_model LE limit_rad)]) * 1e6
	int_flux[tdi_i,i] = max(int_rad_tdi[where(x_rad_model LE limit_rad)]) * 1e6
	int_flux[td_i,i] = max(int_rad_td[where(x_rad_model LE limit_rad)]) * 1e6
	int_flux[tdo_i,i] = max(int_rad_tdo[where(x_rad_model LE limit_rad)]) * 1e6
	int_flux[b_i,i] = max(int_rad_b[where(x_rad_model LE limit_rad)]) * 1e6
	int_flux[di_i,i] = max(int_rad_di[where(x_rad_model LE limit_rad)]) * 1e6
	int_flux[d_i,i] = max(int_rad_d[where(x_rad_model LE limit_rad)]) * 1e6
	int_flux[do_i,i] = max(int_rad_do[where(x_rad_model LE limit_rad)]) * 1e6
	int_lum[*,i] = int_flux[*,i] * flux_to_lum
	printf, t0_lun, format='(A4,7(e13.3))',wave_option[i], (int_lum[tdi_i,i])[0], (int_lum[td_i,i])[0], (int_lum[tdo_i,i])[0], $
			(int_lum[di_i,i])[0], (int_lum[d_i,i])[0], (int_lum[do_i,i])[0], (int_lum[b_i,i])[0]
endfor
free_lun, t0_lun
int_flux_tot = int_flux[tot_i,*]
int_flux_inner = int_flux[tdi_i,*] + int_flux[di_i,*]
int_flux_main = int_flux[td_i,*] + int_flux[d_i,*]
int_flux_outer = int_flux[tdo_i,*] + int_flux[do_i,*]
int_flux_young = int_flux[tdi_i,*] + int_flux[td_i,*] + int_flux[tdo_i,*]	; intrinsic diffuse only lum
int_flux_old = int_flux[di_i,*] + int_flux[d_i,*] + int_flux[do_i,*]
int_flux_bulge = int_flux[b_i,*]
int_flux_tot = int_flux_young + int_flux_old + int_flux_bulge

sfr_err = [-0.0966, 0.102]	; +10.2% error wd01
sfr_err[0] = -0.0966	; -9.66% error wd01
;sfr_err[1] = 0.0467	; +4.67% error lmc1
;sfr_err[0] = -0.0657	; -6.57% error lmc1
;-------------------------------------- CALCULATE STELLAR MASSES ----------------------------------------------
mass_light_Meidt = 0.6	; mass-to-light ratio found in Meidt et al. (2014) +/-0.06 dex [M_sol/solar_lum]
solar_lum = 3.828e26	; luminosity of the sun [W]
Meidt_nu = c/3.6e-6	; frequency corresponding to wavelength 3.6um [Hz]
stellar_mass_conv = (mass_light_Meidt * Meidt_nu * flux_to_lum)/solar_lum	; conversion to solar mass following Meidt et al. (2014) [M_sol/solar_lum]
tot_stellar_mass = int_flux_tot[where(wave_option eq 'ir36')] * stellar_mass_conv
inner_stellar_mass = int_flux_inner[where(wave_option eq 'ir36')] * stellar_mass_conv
main_stellar_mass = int_flux_main[where(wave_option eq 'ir36')] * stellar_mass_conv
outer_stellar_mass = int_flux_outer[where(wave_option eq 'ir36')] * stellar_mass_conv

;print, '-----------------------------------------------------------------------------------------'
;print, 'Total stellar mass: '+strtrim(tot_stellar_mass,1)+' [M_sol]'
;print, '+'+strtrim(tot_stellar_mass*sfr_err[1],1)+', '+strtrim(tot_stellar_mass*sfr_err[0],1)
;print, 'Inner stellar mass: '+strtrim(inner_stellar_mass,1)+' [M_sol]'
;print, '+'+strtrim(inner_stellar_mass*sfr_err[1],1)+', '+strtrim(inner_stellar_mass*sfr_err[0],1)
;print, 'Main stellar mass: '+strtrim(main_stellar_mass,1)+' [M_sol]'
;print, '+'+strtrim(main_stellar_mass*sfr_err[1],1)+', '+strtrim(main_stellar_mass*sfr_err[0],1)
;print, 'Outer stellar mass: '+strtrim(outer_stellar_mass,1)+' [M_sol]'
;print, '+'+strtrim(outer_stellar_mass*sfr_err[1],1)+', '+strtrim(outer_stellar_mass*sfr_err[0],1)
;print, ''
;--------------------------------------------------------------------------------------------------------------
;------------------------------------------ Attenuation Curve -------------------------------------------------
morph_comp = ['tot','tdi','td','tdo','b','di','d','do']
app_lum = app_flux * flux_to_lum        ; convert from flux [Jy] to luminosity [W/Hz]
int_lum = int_flux * flux_to_lum        ; convert from flux [Jy] to luminosity [W/Hz]

atten_global = int_lum[tot_i,*] / app_lum[tot_i,*]
atten_bulge = (int_lum[b_i,*] / app_lum[b_i,*])
atten_inner = (int_lum[tdi_i,*] + int_lum[di_i,*]) / (app_lum[tdi_i,*] + app_lum[di_i,*])
atten_main = (int_lum[td_i,*] + int_lum[d_i,*]) / (app_lum[td_i,*] + app_lum[d_i,*])
atten_outer = (int_lum[tdo_i,*] + int_lum[do_i,*]) / (app_lum[tdo_i,*] + app_lum[do_i,*])
nwave = 89
atten_lim = where(wave_option eq 'uv36')
dwave = (wavelength_nurad[atten_lim] - min(wavelength_nurad))/nwave	; interpolated wavelengths between limits
wave_interp = wavelength_nurad[0] + dwave[0]*findgen(nwave)   ; sets the interpolated wavelengths

;atten_global_interp = spline(wavelength_nurad, atten_global, wave_interp, 50)
;atten_bulge_interp = spline(wavelength_nurad, atten_bulge, wave_interp, 50)
;atten_inner_interp = spline(wavelength_nurad, atten_inner, wave_interp, 50)
;atten_main_interp = spline(wavelength_nurad, atten_main, wave_interp, 50)
;atten_outer_interp = spline(wavelength_nurad, atten_outer, wave_interp, 50)

;wave_interp = [wave_interp[*], wavelength_nurad[atten_lim+1:*]]
;atten_global_interp = [atten_global_interp[*], atten_global[atten_lim+1:*]]
;atten_bulge_interp = [atten_bulge_interp[*], atten_inner[atten_lim+1:*]]
;atten_inner_interp = [atten_inner_interp[*], atten_inner[atten_lim+1:*]]
;atten_main_interp = [atten_main_interp[*], atten_main[atten_lim+1:*]]
;atten_outer_interp = [atten_outer_interp[*], atten_outer[atten_lim+1:*]]

atten_outname = save_dir+model+'_atten_curve.save'
wave = wavelength_nurad
save, wave, atten_global, atten_bulge, atten_inner, atten_main, atten_outer, $
      wave_interp, atten_global_interp, atten_bulge_interp, atten_inner_interp, atten_main_interp, atten_outer_interp, $
	filename=atten_outname
print, 'SAVED: '+atten_outname
;--------------------------------------------------------------------------------------------------------------
;------------------------------------- CALCULATE SFR AND DENSITIES --------------------------------------------

min_wave = 912 * 1d-10
max_wave = 4430 * 1d-10
sfr = dblarr(dim_morph)
for i = 0, dim_morph-1 do sfr[i] = calc_sfr(wavelengths_m, min_wave, max_wave, int_lum[i,*])
print, sfr

items = n_elements(int_lum[*,0])	; number of disk components
int_integ_lum = dblarr(items,2)	; intrinsic integrated luminosity for each disk component [W]
app_integ_lum = dblarr(items,2)	; apparent integrated luminosity for each disk component [W]
int_integ_lum_young = dblarr(items,2)
int_integ_lum_old = dblarr(items,2)
int_integ_lum_bulge = dblarr(items,2)
app_integ_lum_young = dblarr(items,2)
app_integ_lum_old = dblarr(items,2)
app_integ_lum_bulge = dblarr(items,2)
wave = wavelength_nurad * 1e-6	; all nurad wavelengths in [m]
nu = c / wave	; all nurad wavelengths converted to frequency [Hz]
dnu = abs(nu[1:(dim-1)] - nu[0:(dim-2)])	; frequency element for integration [Hz]

;print, 'Integrate over: '+strtrim(wave(0)*1e6,1)+'um to '+strtrim(wave(dim-1)*1e6,1)+'um.'
for i = 0, items-1 do begin	; loops through all disks and integrates over all wavelengths
	int_val = int_lum[i,*]
	int_lum_mean = 0.5 * (int_val[1:dim-1] + int_val[0:dim-2])
	int_integ_lum[i,0] = total(int_lum_mean * dnu)

	app_val = app_lum[i,*]
	app_lum_mean = 0.5 * (app_val[1:dim-1] + app_val[0:dim-2])
        app_integ_lum[i,0] = total(app_lum_mean * dnu)
	
	if morph_comp[i] eq 'tdi' or morph_comp[i] eq 'td' or morph_comp[i] eq 'tdo' then begin
                int_integ_lum_young[i,0] = int_integ_lum[i,0]
                app_integ_lum_young[i,0] = app_integ_lum[i,0]
        endif
        if morph_comp[i] eq 'd' or morph_comp[i] eq 'di' or morph_comp[i] eq 'do' then begin
                int_integ_lum_old[i,0] = int_integ_lum[i,0]
                app_integ_lum_old[i,0] = app_integ_lum[i,0]
        endif
        if morph_comp[i] eq 'b' then begin
                int_integ_lum_bulge[i,0] = int_integ_lum[i,0]
		app_integ_lum_bulge[i,0] = app_integ_lum[i,0]
        endif
endfor
;print, total(int_integ_lum_bulge[*,0]
wave = wavelength_nurad[0:dim_uv] * 1e-6  ; all nurad wavelengths in [m]
nu = c / wave   ; all nurad wavelengths converted to frequency [Hz]
dnu = abs(nu[1:(dim_uv)] - nu[0:(dim_uv-1)])        ; frequency element for integration [Hz]
;print, 'Integrate over: '+strtrim(wave(0)*1e6,1)+'um to '+strtrim(wave(dim_uv)*1e6,1)+'um.'
for i = 0, items-1 do begin     ; loops through all disks and integrates over ONLY UV wavelengths
        int_val = int_lum[i,0:dim_uv]
        int_lum_mean = 0.5 * (int_val[1:dim_uv] + int_val[0:dim_uv-1])
        int_integ_lum[i,1] = total(int_lum_mean * dnu)

        app_val = app_lum[i,0:dim_uv]
        app_lum_mean = 0.5 * (app_val[1:dim_uv] + app_val[0:dim_uv-1])
        app_integ_lum[i,1] = total(app_lum_mean * dnu)
	
	if morph_comp[i] eq 'tdi' or morph_comp[i] eq 'td' or morph_comp[i] eq 'tdo' then begin
		int_integ_lum_young[i,1] = int_integ_lum[i,1]
		app_integ_lum_young[i,1] = app_integ_lum[i,1]
	endif
	if morph_comp[i] eq 'd' or morph_comp[i] eq 'di' or morph_comp[i] eq 'do' then begin
		int_integ_lum_old[i,1] = int_integ_lum[i,1]
                app_integ_lum_old[i,1] = app_integ_lum[i,1]
	endif
	if morph_comp[i] eq 'b' then begin
		int_integ_lum_bulge[i,1] = int_integ_lum[i,1]
		app_integ_lum_bulge[i,1] = app_integ_lum[i,1]
	endif
endfor
tot_intrinsic_lum = int_integ_lum[0,0]
young_intrinsic_lum = total(int_integ_lum_young[*,0])
old_intrinsic_lum = total(int_integ_lum_old[*,0])
bulge_intrinsic_lum = total(int_integ_lum_bulge[*,0])
print, 'd lum [W]: '+strtrim(int_integ_lum[d_i,0],1)
print, 'di lum [W]: '+strtrim(int_integ_lum[di_i,0],1)
print, 'do lum [W]: '+strtrim(int_integ_lum[do_i,0],1)
print, 'td lum [W]: '+strtrim(int_integ_lum[td_i,0],1)
print, 'tdi lum [W]: '+strtrim(int_integ_lum[tdi_i,0],1)
print, 'tdo lum [W]: '+strtrim(int_integ_lum[tdo_i,0],1)
print, 'b lum [W]: '+strtrim(int_integ_lum[b_i,0],1)


;print, '% of light due to young stellar: '+strtrim((young_intrinsic_lum/tot_intrinsic_lum)*100,1)+' %'
;print, '% of light due to old stellar: '+strtrim((old_intrinsic_lum/tot_intrinsic_lum)*100,1)+' %'
;print, '% of light due to bulge: '+strtrim((bulge_intrinsic_lum/tot_intrinsic_lum)*100,1)+' %'

stellar_emission_tot = total(int_integ_lum[0,0])
stellar_inner_emission = total(int_integ_lum[inner_i,0]) / stellar_emission_tot
stellar_bulge_emission = total(int_integ_lum[bulge_i,0]) / stellar_emission_tot
stellar_main_emission = total(int_integ_lum[main_i,0]) / stellar_emission_tot
stellar_outer_emission = total(int_integ_lum[outer_i,0]) / stellar_emission_tot
;print, 'bulge stellar contribution: '+strtrim(stellar_bulge_emission*100,1)+' %'
;print, 'inner stellar contribution: '+strtrim(stellar_inner_emission*100,1)+' %'
;print, 'main stellar contribution: '+strtrim(stellar_main_emission*100,1)+' %'
;print, 'outer stellar contribution: '+strtrim(stellar_outer_emission*100,1)+' %'
;print, ''
bd_ratio = int_integ_lum_bulge[b_i,0]/int_integ_lum_old[old_i,0]
print, 'b/d = '+strtrim(bd_ratio,1)

young_light_abs = total(app_integ_lum_young[*,0])/total(int_integ_lum_young[*,0])
old_light_abs = total(app_integ_lum_old[*,0])/total(int_integ_lum_old[*,0])
bulge_light_abs = total(app_integ_lum_bulge[*,0])/total(int_integ_lum_bulge[*,0])
young_light_abs = total(app_integ_lum_young[*,1])/total(int_integ_lum_young[*,1])
old_light_abs = total(app_integ_lum_old[*,1]+app_integ_lum_bulge[*,1])/$
                total(int_integ_lum_old[*,1]+int_integ_lum_bulge[*,1])
;print, 'young/old absorbed: '+strtrim(young_light_abs/old_light_abs*100,1)+' %'
;print, 'Stellar light absorbed: '+strtrim(100-(app_integ_lum[0,0]/int_integ_lum[0,0]*100),1)+' %'
;print, ''

sfr_global = total(int_integ_lum[0,1]) / lum_lunit_nu	; eqn (17) in PT11.
sfr_inner = total(int_integ_lum[inner_i,1]) / lum_lunit_nu	; eqn (17) in PT11.
sfr_main = total(int_integ_lum[main_i,1]) / lum_lunit_nu	; eqn (17) in PT11.
sfr_outer = total(int_integ_lum[outer_i,1]) / lum_lunit_nu	; eqn (17) in PT11.
print, 'Total global UV SFR = '+strtrim(sfr_global,1)+' [M_sol /yr]'
;print, '+'+strtrim(sfr_global*sfr_err[1],1)+', '+strtrim(sfr_global*sfr_err[0],1)
print, 'Total inner UV SFR = '+strtrim(sfr_inner,1)+' [M_sol /yr]'
;print, '+'+strtrim(sfr_inner*sfr_err[1],1)+', '+strtrim(sfr_inner*sfr_err[0],1)
print, 'Total main UV SFR = '+strtrim(sfr_main,1)+' [M_sol /yr]'
;print, '+'+strtrim(sfr_main*sfr_err[1],1)+', '+strtrim(sfr_main*sfr_err[0],1)
print, 'Total outer UV SFR = '+strtrim(sfr_outer,1)+' [M_sol /yr]'
;print, '+'+strtrim(sfr_outer*sfr_err[1],1)+', '+strtrim(sfr_outer*sfr_err[0],1)
print, ''

area_global = max(area_model[where(x_rad_model le limit_rad)])	; surface area of map at limit radius [kpc^2]
area_inner = area_model[(where(x_rad_model ge rtruncate3))[0]]-area_model[(where(x_rad_model le hs3tin))[0]]
area_main = area_model[(where(x_rad_model ge rtruncate1))[0]]-area_model[(where(x_rad_model le hs1tin))[0]]
area_outer = area_model[(where(x_rad_model ge limit_rad))[0]]-area_model[(where(x_rad_model le hs5tin))[0]]
surface_dens_sfr_tot = sfr_global/area_global	; [M_sol /yr /kpc^2]
surface_dens_sfr_tdi = sfr_inner/area_inner	; [M_sol /yr /kpc^2]
surface_dens_sfr_td  = sfr_main/area_main	; [M_sol /yr /kpc^2]
surface_dens_sfr_tdo = sfr_outer/area_outer	; [M_sol /yr /kpc^2]
print, 'Surface density in SFR, global: '+strtrim(surface_dens_sfr_tot,1)+' [M_sol /yr /kpc^2]'
;print, '+'+strtrim(surface_dens_sfr_tot*sfr_err[1],1)+', '+strtrim(surface_dens_sfr_tot*sfr_err[0],1)
print, 'Surface density in SFR, inner: '+strtrim(surface_dens_sfr_tdi,1)+' [M_sol /yr /kpc^2]'
;print, '+'+strtrim(surface_dens_sfr_tdi*sfr_err[1],1)+', '+strtrim(surface_dens_sfr_tdi*sfr_err[0],1)
print, 'Surface density in SFR, main: '+strtrim(surface_dens_sfr_td,1)+' [M_sol /yr /kpc^2]'
;print, '+'+strtrim(surface_dens_sfr_td*sfr_err[1],1)+', '+strtrim(surface_dens_sfr_td*sfr_err[0],1)
print, 'Surface density in SFR, outer: '+strtrim(surface_dens_sfr_tdo,1)+' [M_sol /yr /kpc^2]'
;print, '+'+strtrim(surface_dens_sfr_tdo*sfr_err[1],1)+', '+strtrim(surface_dens_sfr_tdo*sfr_err[0],1)
print, '-----------------------------------------------------------------------------------------'
surface_dens_mass_tot = tot_stellar_mass/area_global
surface_dens_mass_tdi = inner_stellar_mass/area_inner
surface_dens_mass_td = main_stellar_mass/area_main
surface_dens_mass_tdo = outer_stellar_mass/area_outer
print, 'Surface density in stellar mass, global: '+strtrim(tot_stellar_mass,1)+' [M_sol/kpc^2]'
print, 'Surface density in stellar mass, inner: '+strtrim(inner_stellar_mass,1)+' [M_sol/kpc^2]'
print, 'Surface density in stellar mass, main: '+strtrim(main_stellar_mass,1)+' [M_sol/kpc^2]'
print, 'Surface density in stellar mass, outer: '+strtrim(outer_stellar_mass,1)+' [M_sol/kpc^2]'
print, ''
print, '-----------------------------------------------------------------------------------------'
print, ''
sfr_outname = save_dir+model+'_sfr_data.save'
sfr_upper_err = sfr_err[0]
sfr_lower_err = sfr_err[1]
save,	sfr_global, sfr_inner, sfr_outer, sfr_lower_err, sfr_upper_err, $
	surface_dens_sfr_tot, surface_dens_sfr_tdi, surface_dens_sfr_td, surface_dens_sfr_tdo, $
	surface_dens_mass_tot, surface_dens_mass_tdi, surface_dens_mass_td, surface_dens_mass_tdo, $
	filename=sfr_outname
;--------------------------------------------------------------------------------------------------------------



;--------------------------------------------------- ANALYSIS -------------------------------------------------
sed_filename = obsdir+'obs_sed.save'	; corresponding observed SED
restore, sed_filename
resi = dblarr(n_elements(wavelength))
for i = 0, n_elements(wavelength)-1 do begin	; calculate residuals at available wavelengths
        wave_posi = where(min(abs(wavelength_nurad - wavelength[i])) eq abs(wavelength_nurad - wavelength[i]))
        resi[i] = ((int_jy[i] - app_flux[wave_posi]) /int_jy[i]) * 100
	print, wavelength[i], int_jy[i], app_flux[wave_posi], resi[i]
        if abs(resi[i]) gt 50. then resi[i]=!values.f_NaN
endfor
max_resi = where(max(abs(resi)) eq abs(resi))
print, 'Maximum residual: '+strtrim(resi[max_resi],1)+' %, occurs at '+strtrim(wavelength[max_resi],1)+' um'
print, 'Mean residual: '+strtrim(mean(resi),1)+' %'
print, ''

;--------------------------------------------------------------------------------------------------------------
print, resi
;------------------------------------------------ PLOT --------------------------------------------------------
plot_check = 1
if plot_check eq 1 then begin
xmin = 0.09
xmax = 6
ymax = max(app_flux) *4.0
ymin = ymax*2e-4
plotname = figdir+'sed_'+model+'_'+scaabs+'.ps'
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

cgplot, wavelength_nurad, app_flux_tot, linestyle=0, color='red', ytitle='Flux [Jy]', xrange=[xmin,xmax], yrange=[ymin,ymax],$
        /xlog, /ylog, position=[0.17,0.26,0.95,0.95], xtickformat='(A1)',$
	title=model+'_'+scaabs+' SED'
;cgplot, wavelength_dust, dust_flux_diff[0,*], color='red', linestyle=5, /overplot
cgplot, wavelength_nurad, int_flux_tot, color='purple', linestyle=0, /overplot
;cgplot, wavelength_nurad, app_flux_young, color='blue', linestyle=5, /overplot
;cgplot, wavelength_nurad, app_flux_old, color='green', linestyle=5, /overplot
;cgplot, wavelength_nurad, app_flux, color='black', linestyle=0, /overplot
cgscatter2d, wavelength, int_jy, fit=0, /overplot	; scatter plot of corresponding SED

cglegend, colors=['red','purple','black'], linestyle=[0,0,2], title=['Apparent flux','Intrinsic flux','Data'], $
        length = 0.05, location=[0.2,0.9], vspace=4

cgplot, [xmin,xmax], [20,20], color="black", linestyle=2, position=[0.17,0.1,0.95,0.25], $
        ytitle="R [%]", xrange=[xmin,xmax], xtitle='Wavelength, $\lambda$ [$\mu$m]', /NoErase, yrange=[-50,50], yticks=1, /xlog
cgplot, [xmin,xmax], [-20,-20], color="black", linestyle=2, /overplot
cgplot, [xmin,xmax], [0,0], color="black", /overplot
cgscatter2d, wavelength, resi, fit=0, /overplot

device, /close_file
cgfixps, plotname
print, 'SAVED: '+plotname
endif
;--------------------------------------------------------------------------------------------------------------
outname = save_dir+model+'_sed_'+scaabs+'.save'
;SAVE, tot_lum, dust_lum, wavelength_nurad, wavelength_nurad, wavelength_dust, int_lum_combined, app_lum_combined, $
;	app_lum_tot, app_lum_young, app_lum_old, app_b, app_di, app_d, app_do, app_tdi, app_td, app_tdo, app_ntd, $
;	int_lum_tot, int_lum_young_diff, int_lum_old, int_b, int_di, int_d, int_do, int_tdi, int_td, int_tdo, int_ntd, $
;	tot_sfr_tdisk, tot_sfr_tdi, tot_sfr_td, tot_sfr_tdo, int_lum_young_local, int_lum_young, $
;	llambda, tot_diff, tot_local, dust_inner, dust_main, dust_outer, tot_int_lum, FILENAME=outname
save, morph_comp, wavelength_nurad, app_flux_combined, int_flux_combined, dust_flux, app_flux, tot_int_flux, $
	wavelength_dust, dust_flux_diff, dust_flux_local, wavelength_nurad, app_flux, int_flux, resi, $
	surface_dens_sfr_tot, surface_dens_sfr_tdi, surface_dens_sfr_td, surface_dens_sfr_tdo, $
	surface_dens_mass_tot, surface_dens_mass_tdi, surface_dens_mass_td, surface_dens_mass_tdo, $
	filename=outname
print,''
print, 'SAVED: '+outname
end
