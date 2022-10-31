PRO read_scaling, model, qyear, scaabs, tau, nsersic, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5, f_bd

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
root = '../../'	; root directory.
indata = root+'NUrad/indata/'	; directory where both "scaling.in" and "ff_scaling.in" are contained.

filename = indata+'scaling.in'
READCOL, filename, x, y, SKIPLINE=0, FORMAT="(A,D)"     ;read the data file, skip the first line (which is the number of lines in the file) and assign the data types to each entry
model = x[0]
qyear = x[1]
scaabs = x[2]
tau = (y[WHERE(x EQ 'tau')])[0]
nsersic = (y[WHERE(x EQ 'nsersic')])[0]
sfr = (y[WHERE(x EQ 'sfr')])[0]
sfr4 = (y[WHERE(x EQ 'sfr4')])[0]
sfr6 = (y[WHERE(x EQ 'sfr6')])[0]
sfr7 = (y[WHERE(x EQ 'sfr7')])[0]
old = (y[WHERE(x EQ 'old')])[0]
old3 = (y[WHERE(x EQ 'old3')])[0]
old5 = (y[WHERE(x EQ 'old5')])[0]
bd = (y[WHERE(x EQ 'bd')])[0]
ffactor = (y[WHERE(x EQ 'ffactor')])[0]
ffactor4 = (y[WHERE(x EQ 'ffactor4')])[0]
ffactor6 = (y[WHERE(x EQ 'ffactor6')])[0]
ffactor7 = (y[WHERE(x EQ 'ffactor7')])[0]

PRINT, 'model: '+STRTRIM(model,1)
;PRINT, 'qyear: '+STRTRIM(qyear,1)
PRINT, 'scaabs: '+STRTRIM(scaabs,1)
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
; 17 wavelengths in total, uv36 is at index 8
f_BVIK = f_BVIK[8:16]
f_BVIK3 = f_BVIK3[8:16]
f_BVIK5 = f_BVIK5[8:16]
f_bd = f_bd[8:16]
END
