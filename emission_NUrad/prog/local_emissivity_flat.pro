;this program produces a file of infrared emissivity as a function of
;R,z corresponding to the dust emission from the HII regions

pro local_emissivity_flat,model,qyear,tau,sfr,sfr4,sfr6,sfr7,old,bd,scaabs

start_time = systime(/seconds)
common scaling
close,/all

dir=emission_dir
dir2=urad_dir+'indata/'

suv=strcompress(string(round(sfr*100)),/remove_all)
sbd=strcompress(string(round(bd*100)),/remove_all)
sold=strcompress(string(round(old*100)),/remove_all)
stau=strcompress(string(round(tau*10)),/remove_all)
sffactor = strcompress(string(round(ffactor*100)),/remove_all)

ss = ''

idisk1 = 0L
idisk2 = 0L
idisk3 = 0L
idisk4 = 0L
idisk5 = 0L
idisk6 = 0L
idisk7 = 0L

filename = dir2+'geometry.in'
OPENR, unit, filename, /GET_LUN
READF, unit, ss
READF, unit, tau1
READF, unit, ss
READF, unit, tau2
READF, unit, ss
READF, unit, hd
READF, unit, ss
READF, unit, zd
READF, unit, ss
READF, unit, hdin
READF, unit, ss
READF, unit, zdin
READF, unit, ss
READF, unit, hdsolar
READF, unit, ss
READF, unit, zdsolar
READF, unit, ss
READF, unit, hd1
READF, unit, ss
READF, unit, zd1
READF, unit, ss
READF, unit, hd1in
READF, unit, ss
READF, unit, zd1in
READF, unit, ss
READF, unit, hd1solar
READF, unit, ss
READF, unit, zd1solar
READF, unit, ss
READF, unit, h_bdisk
READF, unit, ss
READF, unit, h_vdisk
READF, unit, ss
READF, unit, h_idisk
READF, unit, ss
READF, unit, h_jdisk
READF, unit, ss
READF, unit, h_kdisk
READF, unit, ss
READF, unit, h_ir36disk
READF, unit, ss
READF, unit, h_ir45disk
READF, unit, ss
READF, unit, h_ir58disk
READF, unit, ss
READF, unit, zs
READF, unit, ss
READF, unit, hsin
READF, unit, ss
READF, unit, zsin
READF, unit, ss
READF, unit, hssolar
READF, unit, ss
READF, unit, zssolar
READF, unit, ss
READF, unit, hs1
READF, unit, ss
READF, unit, zs1
READF, unit, ss
READF, unit, hs1in
READF, unit, ss
READF, unit, zs1in
READF, unit, ss
READF, unit, hs1solar
READF, unit, ss
READF, unit, zs1solar
READF, unit, ss
READF, unit, rtrun
READF, unit, ss
READF, unit, sharp
READF, unit, ss
READF, unit, rtrund
READF, unit, ss
READF, unit, sharpd
READF, unit, ss
READF, unit, rtrun1
READF, unit, ss
READF, unit, sharp1
READF, unit, ss
READF, unit, rtrund1
READF, unit, ss
READF, unit, sharpd1
READF, unit, ss
READF, unit, reff
READF, unit, ss
READF, unit, ellipt
READF, unit, ss
READF, unit, nsersic
READF, unit, ss
READF, unit, xis0
READF, unit, ss
READF, unit, xis1
READF, unit, ss
READF, unit, xid0
READF, unit, ss
READF, unit, xid1
READF, unit, ss
READF, unit, idisk1
READF, unit, ss
READF, unit, idisk2
;read inner component
READF, unit, ss
READF, unit, ss
READF, unit, tau3
READF, unit, ss
READF, unit, tau4
READF, unit, ss
READF, unit, hd3
READF, unit, ss
READF, unit, zd3
READF, unit, ss
READF, unit, hd3in
READF, unit, ss
READF, unit, zd3in
READF, unit, ss
READF, unit, hd3solar
READF, unit, ss
READF, unit, zd3solar
READF, unit, ss
READF, unit, hd4
READF, unit, ss
READF, unit, zd4
READF, unit, ss
READF, unit, hd4in
READF, unit, ss
READF, unit, zd4in
READF, unit, ss
READF, unit, hd4solar
READF, unit, ss
READF, unit, zd4solar
READF, unit, ss
READF, unit, h_bdisk3
READF, unit, ss
READF, unit, h_vdisk3
READF, unit, ss
READF, unit, h_idisk3
READF, unit, ss
READF, unit, h_jdisk3
READF, unit, ss
READF, unit, h_kdisk3
READF, unit, ss
READF, unit, h_ir36disk3
READF, unit, ss
READF, unit, h_ir45disk3
READF, unit, ss
READF, unit, h_ir58disk3
READF, unit, ss
READF, unit, zs3
READF, unit, ss
READF, unit, hs3in
READF, unit, ss
READF, unit, zs3in
READF, unit, ss
READF, unit, hs3solar
READF, unit, ss
READF, unit, zs3solar
READF, unit, ss
READF, unit, hs4
READF, unit, ss
READF, unit, zs4
READF, unit, ss
READF, unit, hs4in
READF, unit, ss
READF, unit, zs4in
READF, unit, ss
READF, unit, hs4solar
READF, unit, ss
READF, unit, zs4solar
READF, unit, ss
READF, unit, rtrun3
READF, unit, ss
READF, unit, sharp3
READF, unit, ss
READF, unit, rtrund3
READF, unit, ss
READF, unit, sharpd3
READF, unit, ss
READF, unit, rtrun4
READF, unit, ss
READF, unit, sharp4
READF, unit, ss
READF, unit, rtrund4
READF, unit, ss
READF, unit, sharpd4
READF, unit, ss
READF, unit, xis3
READF, unit, ss
READF, unit, xis4
READF, unit, ss
READF, unit, xid3
READF, unit, ss
READF, unit, xid4
READF, unit, ss
READF, unit, idisk3
READF, unit, ss
READF, unit, idisk4
;read outer component
READF, unit, ss
READF, unit, ss
READF, unit, tau5
READF, unit, ss
READF, unit, tau6
READF, unit, ss
READF, unit, hd5
READF, unit, ss
READF, unit, zd5
READF, unit, ss
READF, unit, hd5in
READF, unit, ss
READF, unit, zd5in
READF, unit, ss
READF, unit, hd5solar
READF, unit, ss
READF, unit, zd5solar
READF, unit, ss
READF, unit, hd6
READF, unit, ss
READF, unit, zd6
READF, unit, ss
READF, unit, hd6in
READF, unit, ss
READF, unit, zd6in
READF, unit, ss
READF, unit, hd6solar
READF, unit, ss
READF, unit, zd6solar
READF, unit, ss
READF, unit, h_bdisk5
READF, unit, ss
READF, unit, h_vdisk5
READF, unit, ss
READF, unit, h_idisk5
READF, unit, ss
READF, unit, h_jdisk5
READF, unit, ss
READF, unit, h_kdisk5
READF, unit, ss
READF, unit, h_ir36disk5
READF, unit, ss
READF, unit, h_ir45disk5
READF, unit, ss
READF, unit, h_ir58disk5
READF, unit, ss
READF, unit, zs5
READF, unit, ss
READF, unit, hs5in
READF, unit, ss
READF, unit, zs5in
READF, unit, ss
READF, unit, hs5solar
READF, unit, ss
READF, unit, zs5solar
READF, unit, ss
READF, unit, hs6
READF, unit, ss
READF, unit, zs6
READF, unit, ss
READF, unit, hs6in
READF, unit, ss
READF, unit, zs6in
READF, unit, ss
READF, unit, hs6solar
READF, unit, ss
READF, unit, zs6solar
READF, unit, ss
READF, unit, rtrun5
READF, unit, ss
READF, unit, sharp5
READF, unit, ss
READF, unit, rtrund5
READF, unit, ss
READF, unit, sharpd5
READF, unit, ss
READF, unit, rtrun6
READF, unit, ss
READF, unit, sharp6
READF, unit,ss
READF, unit, rtrund6
READF, unit, ss
READF, unit,sharpd6
READF, unit, ss
READF, unit, xis5
READF, unit, ss
READF, unit, xis6
READF, unit, ss
READF, unit, xid5
READF, unit, ss
READF, unit, xid6
READF, unit, ss
READF, unit, idisk5
READF, unit, ss
READF, unit, idisk6
;read inner truncation radius
READF, unit, ss
READF, unit, ss
READF, unit, hstin
READF, unit, ss
READF, unit, hdtin
READF, unit, ss
READF, unit, hs1tin
READF, unit, ss
READF, unit, hd1tin
READF, unit, ss
READF, unit, hs3tin
READF, unit, ss
READF, unit, hd3tin
READF, unit, ss
READF, unit, hs4tin
READF, unit, ss
READF, unit, hd4tin
READF, unit, ss
READF, unit, hs5tin
READF, unit, ss
READF, unit, hd5tin
READF, unit, ss
READF, unit, hs6tin
READF, unit, ss
READF, unit, hd6tin
;read nuclear component
READF, unit, ss
READF, unit, ss
READF, unit, hs7
READF, unit, ss
READF, unit, zs7
READF, unit, ss
READF, unit, hs7in
READF, unit, ss
READF, unit, zs7in
READF, unit, ss
READF, unit, hs7solar
READF, unit, ss
READF, unit, zs7solar
READF, unit, ss
READF, unit, rtrun7
READF, unit, ss
READF, unit, sharp7
READF, unit, ss
READF, unit, xis7
READF, unit, ss
READF, unit, idisk7
READF, unit, ss
READF, unit, hs7tin

free_lun, unit
hs = h_bdisk
shd=strcompress(string(round(hd*1000.)),/remove_all)
szd=strcompress(string(round(zd*1000.)),/remove_all)
shd1=strcompress(string(round(hd1*1000.)),/remove_all)
szd1=strcompress(string(round(zd1*1000.)),/remove_all)
shs=strcompress(string(round(hs*1000.)),/remove_all)
szs=strcompress(string(round(zs*1000.)),/remove_all)
shs1=strcompress(string(round(hs1*1000.)),/remove_all)
szs1=strcompress(string(round(zs1*1000.)),/remove_all)
sreff=strcompress(string(round(reff*1000.)),/remove_all)
sellipt=strcompress(string(round(ellipt*100.)),/remove_all)

hii_filename = dir2+'hii_geometry.in'	; unique geometry for hii regions. Added CJI 30/04/21
IF FILE_TEST(hii_filename) THEN BEGIN
print, 'Using the geometry defined by hii_geometry.in.'
OPENR, unit, hii_filename, /GET_LUN
READF, unit, ss
READF, unit, hs1
READF, unit, ss
READF, unit, zs1
READF, unit, ss
READF, unit, hs1in
READF, unit, ss
READF, unit, zs1in
READF, unit, ss
READF, unit, hs1solar
READF, unit, ss
READF, unit, zs1solar
READF, unit, ss
READF, unit, rtrun1
READF, unit, ss
READF, unit, sharp1
READF, unit, ss
READF, unit, xis1
READF, unit, ss
READF, unit, hs1tin
READF, unit, ss
READF, unit, ss
READF, unit, hs4
READF, unit, ss
READF, unit, zs4
READF, unit, ss
READF, unit, hs4in
READF, unit, ss
READF, unit, zs4in
READF, unit, ss
READF, unit, hs4solar
READF, unit, ss
READF, unit, zs4solar
READF, unit, ss
READF, unit, rtrun4
READF, unit, ss
READF, unit, sharp4
READF, unit, ss
READF, unit, xis4
READF, unit, ss
READF, unit, hs4tin
READF, unit, ss
READF, unit, ss
READF, unit, hs6
READF, unit, ss
READF, unit, zs6
READF, unit, ss
READF, unit, hs6in
READF, unit, ss
READF, unit, zs6in
READF, unit, ss
READF, unit, hs6solar
READF, unit, ss
READF, unit, zs6solar
READF, unit, ss
READF, unit, rtrun6
READF, unit, ss
READF, unit, sharp6
READF, unit, ss
READF, unit, xis6
READF, unit, ss
READF, unit, hs6tin
READF, unit, ss
READF, unit, ss
READF, unit, hs7
READF, unit, ss
READF, unit, zs7
READF, unit, ss
READF, unit, hs7in
READF, unit, ss
READF, unit, zs7in
READF, unit, ss
READF, unit, hs7solar
READF, unit, ss
READF, unit, zs7solar
READF, unit, ss
READF, unit, rtrun7
READF, unit, ss
READF, unit, sharp7
READF, unit, ss
READF, unit, xis7
READF, unit, ss
READF, unit, hs7tin
ENDIF

ztrun1 = 84000. ;in pc

hs1 = hs1 * 1000.d
zs1 = zs1 * 1000.d
hs1in = hs1in * 1000.d
rtrun1 = rtrun1 * 1000.d
zs1in = zs1in * 1000.d
hs1solar = hs1solar * 1000.d
zs1solar = zs1solar * 1000.d

hs4 = hs4 * 1000.d
zs4 = zs4 * 1000.d
hs4in = hs4in * 1000.d
rtrun4 = rtrun4 * 1000.d
zs4in = zs4in * 1000.d
hs4solar = hs4solar * 1000.d
zs4solar = zs4solar * 1000.d

hs6 = hs6 * 1000.d
zs6 = zs6 * 1000.d
hs6in = hs6in * 1000.d
rtrun6 = rtrun6 * 1000.d
zs6in = zs6in * 1000.d
hs6solar = hs6solar * 1000.d
zs6solar = zs6solar * 1000.d

hs7 = hs7 * 1000.d
zs7 = zs7 * 1000.d
hs7in = hs7in * 1000.d
rtrun7 = rtrun7 * 1000.d
zs7in = zs7in * 1000.d
hs7solar = hs7solar * 1000.d
zs7solar = zs7solar * 1000.d

hs1tin=hs1tin * 1000.d
hs4tin=hs4tin * 1000.d
hs6tin=hs6tin * 1000.d
hs7tin=hs7tin * 1000.d

RESTORE, filename= dir+'outdata_intlum/localised_luminosity_'+model+'_q'+qyear+$
       '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_f'+sffactor+'.xdr'


;lambda in Angstrom
;lum_template_scaled in W/Hz
dim_lambda = N_ELEMENTS(lambda)

termz1 = 1.-exp(-ztrun1/zs1)
termr1 = exp(-hs1in/hs1)-exp(-rtrun1/hs1)+(hs1in/hs1)*exp(-hs1in/hs1)-(rtrun1/hs1)*exp(-rtrun1/hs1)

termz4 = 1.-exp(-ztrun1/zs4)
termr4 = exp(-hs4in/hs4)-exp(-rtrun4/hs4)+(hs4in/hs4)*exp(-hs4in/hs4)-(rtrun4/hs4)*exp(-rtrun4/hs4)

termz6 = 1.-exp(-ztrun1/zs6)
termr6 = exp(-hs6in/hs6)-exp(-rtrun6/hs6)+(hs6in/hs6)*exp(-hs6in/hs6)-(rtrun6/hs6)*exp(-rtrun6/hs6)

termz7 = 1.-exp(-ztrun1/zs7)
termr7 = exp(-hs7in/hs7)-exp(-rtrun7/hs7)+(hs7in/hs7)*exp(-hs7in/hs7)-(rtrun7/hs7)*exp(-rtrun7/hs7)

;termr1in = (4./3.) * !Dpi *(1.+xis1/2.)*hs1in^2*zs1*termz*exp(-hs1in/hs1)
termr1in = (4./3.)*!Dpi*((1.+xis1/2.)*hs1in^2-(1-xis1) * hs1tin^3/hs1in - (3./2.) * xis1 * hs1tin^2)*zs1*termz1*exp(-hs1in/hs1)
termr4in = (4./3.)*!Dpi*((1.+xis4/2.)*hs4in^2-(1-xis4) * hs4tin^3/hs4in - (3./2.) * xis4 * hs4tin^2)*zs4*termz4*exp(-hs4in/hs4)
termr6in = (4./3.)*!Dpi*((1.+xis6/2.)*hs6in^2-(1-xis6) * hs6tin^3/hs6in - (3./2.) * xis6 * hs6tin^2)*zs6*termz6*exp(-hs6in/hs6)
termr7in = (4./3.)*!Dpi*((1.+xis7/2.)*hs7in^2-(1-xis7) * hs7tin^3/hs7in - (3./2.) * xis7 * hs7tin^2)*zs7*termz7*exp(-hs7in/hs7)

if hs1in eq 0 then begin
        termr1in(*)=0.
endif
if hs4in eq 0 then begin
        termr4in(*)=0.
endif
if hs6in eq 0 then begin
        termr6in(*)=0.
endif
if hs7in eq 0 then begin
        termr7in(*)=0.
endif

;calculate the emissivity in the centre in W/Hz/pc^3
rho  = lum_template_scaled /(4. * !Dpi * zs1 * hs1^2*termz1*termr1+termr1in)
rho4 = lum_template_scaled4/(4. * !Dpi * zs4 * hs4^2*termz4*termr4+termr4in)
rho6 = lum_template_scaled6/(4. * !Dpi * zs6 * hs6^2*termz6*termr6+termr6in)
rho7 = lum_template_scaled7/(4. * !Dpi * zs7 * hs7^2*termz7*termr7+termr7in)
;help, rho

;read the spatial grid
dim_r = 0L
dim_z = 0L
filename1 = 'urad.in'
name1 = dir2+filename1
OPENR, unit1, name1, /GET_LUN
READF, unit1, pass
READF, unit1, pass
READF, unit1, dim_r
r = DBLARR(dim_r)
for ii = 0, dim_r-1 do begin
	READF, unit1, passr
	r(ii) = passr
endfor
READF, unit1, dim_z
z = DBLARR(dim_z)
for jj = 0, dim_z-1 do begin
	READF, unit1, passz
	z(jj) = passz
endfor 
;transform the spatial position from kpc to pc
r = r * 1000.
z = z * 1000.


kr = 0L
kz = 0L
k = 0L

eta = DBLARR(dim_r,dim_z,dim_lambda)
eta4 = DBLARR(dim_r,dim_z,dim_lambda)
eta6 = DBLARR(dim_r,dim_z,dim_lambda)
eta7 = DBLARR(dim_r,dim_z,dim_lambda)

;      flared disk 
;Atentie, zs11 trebuie declarat vector - vezi linia 283
;Note that zs11 must be declared vector - see line 283
;zs11 = zs1
zs11=DBLARR(dim_r)
zs41=DBLARR(dim_r)
zs61=DBLARR(dim_r)
zs71=DBLARR(dim_r)

if (zs1in gt zs1) then flares11 = alog10((zs1solar - zs1)/(zs1in - zs1)) $
else flares11 = 0.

if (zs4in gt zs4) then flares41 = alog10((zs4solar - zs4)/(zs4in - zs4)) $
else flares41 = 0.

if (zs6in gt zs6) then flares61 = alog10((zs6solar - zs6)/(zs6in - zs6)) $
else flares61 = 0.

if (zs7in gt zs7) then flares71 = alog10((zs7solar - zs7)/(zs7in - zs7)) $
else flares71 = 0.

if (hs1in gt 0.) then begin         
	flares12 = alog10(hs1solar/hs1in)
	flares11 = flares11/flares12
endif else flares11 = 0.

if (hs4in gt 0.) then begin         
	flares42 = alog10(hs4solar/hs4in)
	flares41 = flares41/flares42
endif else flares41 = 0.

if (hs6in gt 0.) then begin         
	flares62 = alog10(hs6solar/hs6in)
	flares61 = flares61/flares62
endif else flares61 = 0.

if (hs7in gt 0.) then begin         
    flares72 = alog10(hs7solar/hs7in)
    flares71 = flares71/flares72
endif else flares11 = 0.


;zs11 = zs1 + (zs1in - zs1) * (r/hs1in)^flares11

;Atentie aici! Trateaza separat cazul hs1in=0
;Be careful here! Treat the case hs1in = 0 separately
for kr = 0L,dim_r-1 do begin
   if (hs1in gt 0.) then zs11[kr] = zs1 + (zs1in - zs1) * (r[kr]/hs1in)^flares11 $
   else zs11[kr]=zs1

   if (hs4in gt 0.) then zs41[kr] = zs4 + (zs4in - zs4) * (r[kr]/hs4in)^flares41 $
   else zs41[kr]=zs4

   if (hs6in gt 0.) then zs61[kr] = zs6 + (zs6in - zs6) * (r[kr]/hs6in)^flares61 $
   else zs61[kr]=zs6

   if (hs7in gt 0.) then zs71[kr] = zs7 + (zs7in - zs7) * (r[kr]/hs7in)^flares71 $
   else zs71[kr]=zs7
endfor

for kr = 0L,dim_r-1 do begin
 for kz = 0L,dim_z-1 do begin
  if idisk2 eq 1 then begin ;double exponential distribution
    if r[kr] ge hs1in then begin
        eta(kr,kz,*) = rho(*) * (zs1/zs11[kr]) * exp(-r[kr]/hs1-z(kz)/zs11[kr])
    endif
    if r[kr] lt hs1in then begin
         eta(kr,kz,*) = rho(*) * (zs1/zs11[kr]) * $
                        ((r[kr]/hs1in)*(1. - xis1) + xis1) * $
                        exp(-hs1in/hs1-z(kz)/zs11[kr])
    endif
  endif ;end double exponential distribution

  if idisk4 eq 1 then begin ;double exponential distribution
    if r[kr] ge hs4in then begin
        eta4(kr,kz,*) = rho4(*) * (zs4/zs41[kr]) * exp(-r[kr]/hs4-z(kz)/zs41[kr])
    endif
    if r[kr] lt hs4in then begin
         eta4(kr,kz,*) = rho4(*) * (zs4/zs41[kr]) * $
                        ((r[kr]/hs4in)*(1. - xis4) + xis4) * $
                        exp(-hs4in/hs4-z(kz)/zs41[kr])
    endif
  endif ;end double exponential distribution

  if idisk6 eq 1 then begin ;double exponential distribution
    if r[kr] ge hs6in then begin
        eta6(kr,kz,*) = rho6(*) * (zs6/zs61[kr]) * exp(-r[kr]/hs6-z(kz)/zs61[kr])
    endif
    if r[kr] lt hs6in then begin
         eta6(kr,kz,*) = rho6(*) * (zs6/zs61[kr]) * $
                        ((r[kr]/hs6in)*(1. - xis6) + xis6) * $
                        exp(-hs6in/hs6-z(kz)/zs61[kr])
    endif
  endif ;end double exponential distribution

  if idisk7 eq 1 then begin ;double exponential distribution
    if r[kr] ge hs7in then begin
        eta7(kr,kz,*) = rho7(*) * (zs7/zs71[kr]) * exp(-r[kr]/hs7-z(kz)/zs71[kr])
    endif
    if r[kr] lt hs7in then begin
         eta7(kr,kz,*) = rho7(*) * (zs7/zs71[kr]) * $
                        ((r[kr]/hs7in)*(1. - xis7) + xis7) * $
                        exp(-hs7in/hs7-z(kz)/zs71[kr])
    endif
  endif ;end double exponential distribution

  sech11 = 1.d0/cosh(z(kz)/zs11[kr])
  sech12=sech11*sech11
  sech41 = 1.d0/cosh(z(kz)/zs41[kr])
  sech42=sech41*sech41
  sech61 = 1.d0/cosh(z(kz)/zs61[kr])
  sech62=sech61*sech61
  sech71 = 1.d0/cosh(z(kz)/zs71[kr])
  sech72=sech71*sech71

  if idisk2 eq 3 then begin ;sech2 vertical distribution
    if r[kr] ge hs1in then begin
        eta(kr,kz,*) = rho(*) * (zs1/zs11[kr]) * exp(-r[kr]/hs1) * sech12
    endif
    if r[kr] lt hs1in then begin
         eta(kr,kz,*) = rho(*) * (zs1/zs11[kr]) * $
                        ((r[kr]/hs1in)*(1. - xis1) + xis1) * $
                        exp(-hs1in/hs1) * sech12
    endif
  endif ;end sech2 vertical distribution

  if idisk4 eq 3 then begin ;sech2 vertical distribution
    if r[kr] ge hs4in then begin
        eta4(kr,kz,*) = rho4(*) * (zs4/zs41[kr]) * exp(-r[kr]/hs4) * sech42
    endif
    if r[kr] lt hs4in then begin
         eta4(kr,kz,*) = rho4(*) * (zs4/zs41[kr]) * $
                        ((r[kr]/hs4in)*(1. - xis4) + xis4) * $
                        exp(-hs4in/hs4) * sech42
    endif
  endif ;end sech2 vertical distribution

  if idisk6 eq 3 then begin ;sech2 vertical distribution
    if r[kr] ge hs6in then begin
        eta6(kr,kz,*) = rho6(*) * (zs6/zs61[kr]) * exp(-r[kr]/hs6) * sech62
    endif
    if r[kr] lt hs6in then begin
         eta6(kr,kz,*) = rho6(*) * (zs6/zs61[kr]) * $
                        ((r[kr]/hs6in)*(1. - xis6) + xis6) * $
                        exp(-hs6in/hs6) * sech62
    endif
  endif ;end sech2 vertical distribution

  if idisk7 eq 3 then begin ;sech2 vertical distribution
    if r[kr] ge hs7in then begin
        eta7(kr,kz,*) = rho7(*) * (zs7/zs71[kr]) * exp(-r[kr]/hs7) * sech72
    endif
    if r[kr] lt hs7in then begin
         eta7(kr,kz,*) = rho7(*) * (zs7/zs71[kr]) * $
                        ((r[kr]/hs7in)*(1. - xis7) + xis7) * $
                        exp(-hs7in/hs7) * sech72
    endif
  endif ;end sech2 vertical distribution

 endfor 
endfor

;;;inserted by JJT 13/04/17;;;
for kr = 0L,dim_r-1 do begin
 for kz = 0L,dim_z-1 do begin
;truncating HII emissivity
      rplane=sqrt(r[kr]*r[kr]+z(kz)*z(kz))
      rplane=double(rplane)
      sigtruncate1=sharp1*rtrun1
      truncate1=0.5*(1.0-erf((rplane-rtrun1)/sigtruncate1))

      sigtruncate4=sharp4*rtrun4
      truncate4=0.5*(1.0-erf((rplane-rtrun4)/sigtruncate4))

      sigtruncate6=sharp6*rtrun6
      truncate6=0.5*(1.0-erf((rplane-rtrun6)/sigtruncate6))

      sigtruncate7=sharp7*rtrun7
      truncate7=0.5*(1.0-erf((rplane-rtrun7)/sigtruncate7))

      eta(kr,kz,*)=eta(kr,kz,*)*truncate1
      eta4(kr,kz,*)=eta4(kr,kz,*)*truncate4
      eta6(kr,kz,*)=eta6(kr,kz,*)*truncate6
      eta7(kr,kz,*)=eta7(kr,kz,*)*truncate7

; added inner truncation radius
        if(r[kr] lt hs1tin) then begin
         eta(kr,kz,*)=0
        endif

        if(r[kr] lt hs4tin) then begin
         eta4(kr,kz,*)=0
        endif


        if(r[kr] lt hs6tin) then begin
         eta6(kr,kz,*)=0
        endif

        if(r[kr] lt hs7tin) then begin
         eta7(kr,kz,*)=0
        endif

 endfor
endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
IF ffactor EQ 0. THEN eta[*,*,*] = 0
IF ffactor4 EQ 0. THEN eta4[*,*,*] = 0
IF ffactor6 EQ 0. THEN eta6[*,*,*] = 0
IF ffactor7 EQ 0. THEN eta7[*,*,*] = 0
print,'max eta: '+strtrim(max(eta),1)
print,'max eta4: '+strtrim(max(eta4),1)
print,'max eta6: '+strtrim(max(eta6),1)
print,'max eta7: '+strtrim(max(eta7),1)

eta_plot = eta + eta4; + eta6; + eta7
set_plot,'ps'
name=dir+'figures/local_emissivity_r.ps'
PRINT, 'plot ', name
device,filename=name
    !p.thick=3
    !x.thick=4
    !y.thick=4
    !p.charthick=3
    !p.charsize=1.5

plot_io,r/1.d+3,eta_plot(*,0,233),$
        xtitle='R [kpc]', ytitle='L[W Hz!U-1!Npc!U-3]',linestyle=0, yrange=[max(eta_plot)*1e-3,max(eta_plot)], xrange=[0,7]
device,/close
;PRINT, eta(*,0,300)

set_plot,'ps'
name=dir+'figures/local_emissivity_z.ps'
PRINT, 'plot ', name
device,filename=name
    !p.thick=3
    !x.thick=4
    !y.thick=4
    !p.charthick=3
    !p.charsize=1.5

plot_io,z/1.d+3,eta(0,*,233),$
        xtitle='z [kpc]', ytitle='L[W Hz!U-1!Npc!U-3]',linestyle=0, yrange=[max(eta_plot)*1e-3,max(eta_plot)]
device,/close

out_name = dir+$
'outdata_intlum/grid_irr_HII_'+model+'_q'+qyear+$
'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_f'+sffactor+'.xdr'
save,lambda,r,z,eta,eta4,eta6,eta7,file=out_name
PRINT, 'Written: '+out_name

mark1:
print, 'DONE: local_emissivity_flat.pro ('+strtrim(ceil(systime(/seconds)-start_time),1)+' s)'
end
