PRO read_scaling_geom, model, scaabs, qyear, tau, nsersic, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5, f_bd, hs, shd, szd, shd1, szd1, shs, szs, shs1, szs1, sreff, sellipt, strsfr, strold, strbd, stau, snsersic, sinclination, dist_in, pix_size, nx_b, nx_n, nx_i, nx_o, nx_m

; The purpose of this code is to read in the scaling factors defined by "scaling.in" and "ff_scaling.in".
; "scaling.in" contains the following information:
;	- model: string of the desired dust model.
;	- qyear: 
;	- scaabs: either "sca" for scattering, or "abs" for absorption.
;	- tau: the value in the string output from NUrad divided by 10, i.e.: for *_t254_*, tau=25.4.
;	- nersic: the sersic index defined in "geometry.in".
;	- sfr,sfr4.sfr6,sfr7: the star formation rate for uv22 for each respective morphological component.
;	- old,old3,old5: the amplitude scaling for the b-band for each respective morphological component.
;	- bd: the amplitude scaling for the bulge in the b band.
;	- ffactor,ffactor4,ffacotr6,ffactor7: scaling factor for the dust emission in the HII regions.
;
; "ff_scaling.in" contains the following information:
;	- f_uv,f_uv4,f_uv6,f_uv7: wavelength dependant scaling factors for young stars (scaled relative to sfr)
;	- f_BVIK,f_BVIK3,f_BVIK5: wavelength dependant scaling factors for old stars (scaled relative to old)
;	- f_bd: wavelength dependant scaling factors for the bulge (scaled relative to bd)
;
; To use this program, include the following line in any program you wish to call these values;
;
;	- read_scaling, model, qyear, scaabs, tau, nsersic, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5
;
; Please ensure that this program is included in the directory you wish to use it in.
;
root = '/nfs/d58/vlk/sedmodel/cinman/'	; root directory.
indata = root+'m51a/NUrad_M51a/indata/'+model+'_'+scaabs+'/'	; directory where both "scaling.in" and "ff_scaling.in" are contained.
print, "Read from: "+indata

filename = indata+'scaling.in'
READCOL, filename, x, y, SKIPLINE=0, FORMAT="(A,D)"     ;read the data file, skip the first line (which is the number of lines in the file) and assign the data types to each entry
;model = x(0)
qyear = x(1)
;scaabs = x(2)
tau = (y(WHERE(x EQ 'tau')))(0)
nsersic = (y(WHERE(x EQ 'nsersic')))(0)
sfr = (y(WHERE(x EQ 'sfr')))(0)
sfr4 = (y(WHERE(x EQ 'sfr4')))(0)
sfr6 = (y(WHERE(x EQ 'sfr6')))(0)
sfr7 = (y(WHERE(x EQ 'sfr7')))(0)
old = (y(WHERE(x EQ 'old')))(0)
old3 = (y(WHERE(x EQ 'old3')))(0)
old5 = (y(WHERE(x EQ 'old5')))(0)
bd = (y(WHERE(x EQ 'bd')))(0)
ffactor = (y(WHERE(x EQ 'ffactor')))(0)
ffactor4 = (y(WHERE(x EQ 'ffactor4')))(0)
ffactor6 = (y(WHERE(x EQ 'ffactor6')))(0)
ffactor7 = (y(WHERE(x EQ 'ffactor7')))(0)
;PRINT, 'model: '+STRTRIM(model,1)
;PRINT, 'qyear: '+STRTRIM(qyear,1)
;PRINT, 'scaabs: '+STRTRIM(scaabs,1)
;PRINT, 'tau: '+STRTRIM(tau,1)
;PRINT, 'nsersic: '+STRTRIM(nsersic,1)
;PRINT, 'sfr: '+STRTRIM(sfr,1)
;PRINT, 'sfr4: '+STRTRIM(sfr4,1)
;PRINT, 'sfr6: '+STRTRIM(sfr6,1)
;PRINT, 'sfr7: '+STRTRIM(sfr7,1)
;PRINT, 'old: '+STRTRIM(old,1)
;PRINT, 'old3: '+STRTRIM(old3,1)
;PRINT, 'old5: '+STRTRIM(old5,1)
;PRINT, 'bd: '+STRTRIM(bd,1)
;PRINT, 'ffactor: '+STRTRIM(ffactor,1)
;PRINT, 'ffactor4: '+STRTRIM(ffactor4,1)
;PRINT, 'ffactor6: '+STRTRIM(ffactor6,1)
;PRINT, 'ffactor7: '+STRTRIM(ffactor7,1)
;PRINT, ''

filename = indata+'ff_scaling.in'
READCOL, filename, wave, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5, f_bd, SKIPLINE=1, FORMAT="(A,D,D,D,D,D,D,D,D)"

f_BVIK = f_BVIK(WHERE(f_BVIK NE 0))
f_BVIK3 = f_BVIK3(WHERE(f_BVIK5 NE 0))
f_BVIK5 = f_BVIK5(WHERE(f_BVIK5 NE 0))
f_bd = f_bd(WHERE(f_bd NE 0))

ss=''
name = indata+'geometry.in'
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

filename_param = indata+'gal_param.in'
name_param = filename_param
openr,unit,name_param,/get_lun
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
sinclination=strcompress(string(fix(inclination)),/remove_all)

END
