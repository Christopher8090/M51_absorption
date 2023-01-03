PRO diag_emissivity

; The purpose of this code is to produce diagnostic plots of the dust emission.
; Produces radial profiles at a height z=0 for the wavelengths specified in the `wavelength` array.

;------------------------------------------------ DEFINE DIRECTORIES -------------------------------------------
root = '../../'
em_dir = root+'emission_NUrad/outdata_intlum/'
out_dir = root+'figures/'
savedir = root+'saves/model/'
;---------------------------------------------------------------------------------------------------------------
;---------------------------------------------- DEFINE INPUT PARAMETERS ----------------------------------------
read_scaling, model, qyear, scaabs, tau, nsersic, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5

ss=''
name = root+'NUrad/indata/geometry.in'
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
;sinclination=strcompress(string(round(inclination)),/remove_all)

wavelength = ['3.6','4.5','5.8','8.0','24','70','100','160','250','350','500','850']
item=['irr1','irr2','irr3','irr4','irr5','irr6']
n_items = n_elements(item)
;----------------------------------------------------------------------------------------------------------------
;-------------------------------------------- GET DUST EMISSION DATA---------------------------------------------
FOR j = 0, N_ELEMENTS(wavelength)-1 DO BEGIN	; loop through wavelengths of interest
if wavelength[j] ne '500' then continue
FOR i = 0, n_elements(item)-1 DO BEGIN	; loop through morphological components
RESTORE, em_dir+'grid_'+item[i]+'_'+model+'_q'+qyear+'_t'+stau+'_s'+strsfr+'_no'+strold+'_bd'+strbd+'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+'_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+'_reff'+sreff+'_ell'+sellipt+'_'+scaabs+'.xdr'
endfor
z_height = 0.d		;defines height of z being looked at
radius = WHERE(zzz EQ z_height)	;finds the all radii at the given z_height
lam = lambda/1E4	;converts lambda values to microns
wave_arr = MIN(ABS(lam-wavelength[j]),posi)	;finds the closest wavelength to the input wavelength

emiss_irr1 = lumdouble1[radius,z_height,posi]	; main thick disk
emiss_irr2 = lumdouble2[radius,z_height,posi]	; main thin disk
emiss_irr3 = lumdouble3[radius,z_height,posi]	; inner thick disk
emiss_irr4 = lumdouble4[radius,z_height,posi]	; inner thin disk
emiss_irr5 = lumdouble5[radius,z_height,posi]	; outer thick disk
emiss_irr6 = lumdouble6[radius,z_height,posi]	; outer thin disk
emiss_tot = emiss_irr1+emiss_irr2+emiss_irr3+emiss_irr4+emiss_irr5+emiss_irr6	; all disks
;---------------------------------------------------------------------------------------------------------------
;-------------------------------------------------- SAVE & PLOT ------------------------------------------------
ymax = MAX(emiss_tot) * 1.5	;determines 150% maximum value of emiss_tot for plotting range
ymin = ymax * 1e-7
r_range = rrr[radius]*1E-3	;defines the plotting range of the radius and converts from pc to kpc

save_name = savedir+model+'_diag_'+wavelength[j]+'um.save'
SAVE, r_range, emiss_tot, emiss_irr1, emiss_irr2, emiss_irr3, emiss_irr4, emiss_irr5, emiss_irr6, FILENAME=save_name
PRINT, 'Saved: '+save_name
filename = out_dir+'radial_emissivity_'+STRTRIM(wavelength[j],1)+'um.ps'
header = STRTRIM(wavelength[j],1)+'$\mu$m at z='+STRTRIM(fix(z_height),1)+'pc'
y_name = 'Luminosity [Arb. units]'

aspect = cgPSWINDOW()
thisdevice = !D.NAME
SET_PLOT,'ps', /Copy
DEVICE,FILENAME = filename, XSIZE=aspect.xsize, YSIZE=aspect.ysize, XOFFSET=aspect.xoffset, YOFFSET=aspect.yoffset, COLOR=1, $
ENCAPSULATED=encapsulated, Inches=aspect.inches,bits=8,set_character_size=[180,200],set_font='HELVETICA',/LANDSCAPE
	!p.thick=6
	!x.thick=3
	!y.thick=3
        !p.charthick=3
       	!p.charsize=2

cgPLOT, r_range, emiss_irr1, LINESTYLE=4, color='orange', XTITLE='Radius [kpc]', YTITLE=y_name, TITLE=header, XRANGE=[0,25], YRANGE=[ymin,ymax], /YLOG
cgPLOT, r_range, emiss_irr2, linestyle=5, color='orange', /overplot
cgPLOT, r_range, emiss_irr3, linestyle=4, color='green', /overplot
cgPLOT, r_range, emiss_irr4, linestyle=5, color='green', /overplot
cgPLOT, r_range, emiss_irr5, linestyle=4, color='blue', /overplot
cgPLOT, r_range, emiss_irr6, linestyle=5, color='blue', /overplot
cgPLOT, r_range, emiss_tot, linestyle=0, color='black', /overplot
	
cgLegend, Colors=['black','orange','orange','green','green','blue','blue'], Location= [0.75,0.85],$
	Linestyle=[0,4,5,4,5,4,5], PSyms=[0,0,0,0,0,0,0], Titles=['Total','irr1','irr2','irr3','irr4','irr5','irr6'], VSpace=2.75, /Background

;cgPLOT, r_range, emiss_tot*0, linestyle=1, /overplot

DEVICE,/CLOSE_FILE
cgFIXPS, filename
PRINT, 'Plotted: '+filename
stop
;---------------------------------------------------------------------------------------------------------------
ENDFOR
END
