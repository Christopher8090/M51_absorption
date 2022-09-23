;this program scales the radiation fields corresponding to a certain
;SFR, old and B/D ratio 
;it uses a different template for the old
;JJT 6.10.17 - 13.10.17 - updated to deal
;with extra components
;JJT 19MAR18 decouple dust and stellar truncation radii


pro radiation_fields_template,swdisk3,swdisk4,swdisk5,swdisk6,swdisk7

;input parameters:
; model - the dust model used (e.g. 'wd01' or 'pop00' refers to the model of
;         Weingartner and Draine 2001 or Popescu et al. 2000)
; qyear - year of last updated the extinction emissivities (e.g. '06' or'01')
; tau   - the central face-on opacity in the B band
; sfr   - the equivalent sfr corresponding to the non-ionising UV luminosity
;         that escapes in the diffuse component (in units of 
;         1Msolar/yr=2.495x1.e+36 W, as given by the standard population 
;         synthesis)
; sfr4/6 - sfr for extra thin disk components 
; old- a scaling factor for the optical luminosity of the old stellar disk
;         where one unit of optical luminosity corresponds to 10 times the
;         luminosity of the non-ionisng UV photons for a total (real)
;         SFR=1Msolar/yr=2.495x1.e+36W
; old3/5 - old for extra disk components    
; bd        the bulge-to-disk ratio
; swdisk - do you want to take a particular disk into account 'yes'/'no'
start_time = systime(/seconds)
;common dirdef,rootdir
close,/all
common scaling
;get directories from ./mydirectories.in BTC 2019-11-18
;ss=''
;urad_dir=''
;emission_dir=''
;OPENR, dirunit, 'mydirectories.in', /GET_LUN
;	READF, dirunit, ss
;	READF, dirunit, urad_dir
;	READF, dirunit, ss
;	READF, dirunit, emission_dir
;FREE_LUN, dirunit

dir = urad_dir+'indata/'
dir1 = urad_dir+'unit/'
dir2 = emission_dir+'indata/'

rootdir = '../../'
common param, shd, szd, shd1, szd1, shs, szs, shs1, szs1, sreff, sellipt

filter_opt=['uv36','b','v','i','j','k','ir36','ir45','ir58'];,'nir']
dim_opt = n_elements(filter_opt)
filter_uv = ['uv09','uv13','uv15','uv16','uv20','uv22','uv25','uv28'];,'uv36']
dim_uv = n_elements(filter_uv)
filter = ['uv09','uv13','uv15','uv16','uv20','uv22','uv25','uv28','uv36','b','v','i','j','k','ir36','ir45','ir58'];,'nir']
dim = n_elements(filter)

;lambda=[912.,1350.,1500.,1650.,2000.,2200.,2500.,2800.,3650.,4430.,5640.,8090.,12590.,22000.,50000.] * 1.d ;in AA
lambda=[912.,1350.,1500.,1650.,2000.,2200.,2500.,2800.,3650.,4430.,5640.,8090.,12590.,22000.,36000.,45000.,58000.] * 1.d ;in AA

;define a wavelength dependence f factor for the young uv disk
;f=[0.573,0.516,0.473,0.455,0.372,0.324,0.261,0.206,0.108,0.068,0.038,0.016,0.009,0.006,0.001] ;Fcal=0.35
f=[0.573,0.516,0.473,0.455,0.372,0.324,0.261,0.206,0.108,0.068,0.038,0.016,0.009,0.006,0.001,0.001,0.001] ;Fcal=0.35
f_diff = 1.- f

;define a correction in the B,V,I bands to account for a different
;intrinsic SED 
;final old SED 20NOV18
;   f_BVIK = [0.3,1.,0.7,0.6,0.09,0.1,0.24,0.3,0.4]   ;old=0.2
;   f_BVIK3 = [0.1,1.,0.67,0.5,0.1,0.23,0.72,0.5,0.4]   ;old3=0.0025
;   f_BVIK5 = [0.04,1.,0.9,1.5,0.4,0.4,2.491,1.62,1.08]   ;old5=0.005

;      ['uv09','uv13','uv15','uv16','uv20','uv22','uv25','uv28','uv36','b','v','i','j','k','nir']
;final young SED 20NOV18
;   f_uv = [1.,0.9,0.9,1.,1.,1.,1.,1.,0.25,0.2,0.1,0.05,0.11,0.1,0.46,0.32,0.28] ;sfr=0.58
;   f_uv4 = [1.,0.73,0.73,1.,1.,1.,1.,1.,0.3,0.2,0.12,0.1,0.3,0.01,0.012,0.0125,0.01] ;sfr4=0.025
;   f_uv6 = [1.,0.8,0.8,1.,1.,1.,1.,1.,0.35,0.14,0.06,0.03,0.1,0.01,0.043,0.022,0.022] ;sfr6=0.025
;   f_uv7 = [1.,0.58,0.58,1.,1.,1.,1.,1.,0.6,0.9,0.9,0.6,5.,4.5,6.,4.,4.] ;sfr7=0.002


;test
;   f_uv = [0.001,0.9,0.9,1.,1.,1.,1.,1.,0.25,0.2,0.1,0.05,0.11,0.1,0.46,0.32,0.28] ;sfr=0.58
;   f_uv4 = [0.001,0.73,0.73,1.,1.,1.,1.,1.,0.3,0.2,0.12,0.1,0.3,0.01,0.012,0.0125,0.01] ;sfr4=0.025
;   f_uv6 = [0.001,0.8,0.8,1.,1.,1.,1.,1.,0.35,0.14,0.06,0.03,0.1,0.01,0.043,0.022,0.022] ;sfr6=0.025
;   f_uv7 = [0.001,0.58,0.58,1.,1.,1.,1.,1.,0.6,0.9,0.9,0.6,5.,4.5,6.,4.,4.] ;sfr7=0.002

;read_scaling, model, qyear, scaabs, tau, nsersic, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5, f_bd
;define a correction factor to derive the value of old from the
;standard template
;corr_old = 1.099
corr_old=1.
lum_lunit_nu = 2.241d+36 ;in W
lum_lunit_diff_nu = 1.437d+36 ;in W

ss = ''
idisk1 = 0L
idisk2 = 0L
idisk3 = 0L
idisk4 = 0L
idisk5 = 0L
idisk6 = 0L
idisk7 = 0L

filename = 'geometry.in'
name = dir+filename
openr,unit,name,/get_lun
;print, 'read ', name
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
readf, unit, idisk1
readf, unit, ss
readf, unit, idisk2

free_lun, unit

hs = h_bdisk
shd=strcompress(string(fix(hd*1000.)),/remove_all)
szd=strcompress(string(fix(zd*1000.)),/remove_all)
shd1=strcompress(string(fix(hd1*1000.)),/remove_all)
szd1=strcompress(string(fix(zd1*1000.)),/remove_all)
shs=strcompress(string(fix(hs*1000.)),/remove_all)
szs=strcompress(string(fix(zs*1000.)),/remove_all)
shs1=strcompress(string(fix(hs1*1000.)),/remove_all)
szs1=strcompress(string(fix(zs1*1000.)),/remove_all)
sreff=strcompress(string(fix(reff*1000.)),/remove_all)
sellipt=strcompress(string(fix(ellipt*100.)),/remove_all)

combine = 0
;Allows total unit radiation fields to be written for the combined dust model
if combine eq 1 then begin
	;Over-write read scaling factors to allow for combining lmc+wd01
	sfr = 5.2
	sfr4 = 3.5
	sfr6 =0.7 
	sfr7 = 0.
	old = 1.7
	lmc_old = 1.5
	old3 =0.8  
	old5 = 0.55
	bd = 0.012
	f_BVIK = [0.18,1.,0.75,0.37,0.17,0.19,0.45,0.28,0.72]
	f_BVIK3 = [0.15,1.,0.78,0.36,0.15,0.17,0.46,0.31,0.45]
	f_BVIK5 = [0.22,1.,0.69,0.39,0.15,0.15,0.39,0.36,1.1]
	f_bd = [0.25,1.,0.78,0.38,0.17,0.16,0.42,0.38,0.28]
	
	f_uv = [0.45,0.5,0.6 ,0.66,0.78,1.,0.88,0.72,0.26,0.2 ,0.1 ,0.05,0.2 ,0.2 ,0.5 ,0.5 ,0.5]
	f_uv4 = [0.2,0.3,0.58,0.6 ,0.78,1.,0.78,0.6 ,0.16,0.13,0.06,0.03,0.15,0.15,0.15,0.15,0.45]
	f_uv6 = [1. ,1. ,1.  ,1.  ,1.  ,1.,1.  ,1.  ,0.45,0.3 ,0.25,0.15,0.3 ,0.8 ,1.  ,1.5 ,1.5]
	f_uv7 = [1  ,1  ,1   ,1   ,1   ,1 ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1] ;sfr7=0.002
endif
;enable/disable the below for calculations of energy absorbed by different populations
;sfr = 0.
;sfr4 = 0.
;sfr6 = 0.
;old = 0.
;old3 = 0.
;old5 = 0.
;bd = 0.

ff = dblarr(dim)
ff_td = sfr * f_diff * (lum_lunit_nu/lum_lunit_diff_nu) *f_uv ;added by JJT 10/8/17
ff_d = ff
ff_td4 = sfr4 * f_diff * (lum_lunit_nu/lum_lunit_diff_nu) *f_uv4 ;added by JJT 6/10/17
ff_d3 = ff ;added by JJT 6/10/17
ff_td6 = sfr6 * f_diff * (lum_lunit_nu/lum_lunit_diff_nu) *f_uv6 ;added by JJT 6/10/17
ff_td7 = sfr7 * f_diff * (lum_lunit_nu/lum_lunit_diff_nu) *f_uv7
ff_d5 = ff ;added by JJT 6/10/17
ff_b = ff
ff_d[(dim-dim_opt):(dim-1)]=(old/corr_old) * f_BVIK
ff_d3[(dim-dim_opt):(dim-1)]=(old3/corr_old) * f_BVIK3 ;added by JJT 6/10/17
ff_d5[(dim-dim_opt):(dim-1)]=(old5/corr_old) * f_BVIK5 ;added by JJT 6/10/17
ff_b[(dim-dim_opt):(dim-1)]=(old/corr_old) * f_bd * bd
if combine eq 1 then ff_b[(dim-dim_opt):(dim-1)]=(lmc_old/corr_old) * f_bd * bd

;define idisk and ibulge in the uv
idisk='no'
ibulge='no'
idisk3='no' ;added by JJT 6/10/17
idisk5='no' ;added by JJT 6/10/17
for ii = 0, dim_uv-1 do begin
readwrite_rf_template,filter[ii],ibulge,idisk,idisk3,idisk5,$
ff_b[ii],ff_d[ii],ff_td[ii],ff_d3[ii],ff_td4[ii],ff_d5[ii],ff_td6[ii],ff_td7[ii],$
dir1,dir2,swdisk3,swdisk4,swdisk5,swdisk6,swdisk7
endfor

;define idisk and ibulge in the optical
if old eq 0. then idisk='no' else idisk='yes'
if old3 eq 0. then idisk3='no' else idisk3='yes' ;added by JJT 6/10/17
if old5 eq 0. then idisk5='no' else idisk5='yes' ;added by JJT 6/10/17
if bd eq 0. then ibulge='no' else ibulge='yes'
for ii = dim_uv, dim-1 do begin
readwrite_rf_template,filter[ii],ibulge,idisk,idisk3,idisk5,$
ff_b[ii],ff_d[ii],ff_td[ii],ff_d3[ii],ff_td4[ii],ff_d5[ii],ff_td6[ii],ff_td7[ii],$
dir1,dir2,swdisk3,swdisk4,swdisk5,swdisk6,swdisk7
    
endfor
print, 'Done: radiation_fields_unit_template.pro ('+strtrim(ceil(systime(/seconds)-start_time),1)+' s)'
end
