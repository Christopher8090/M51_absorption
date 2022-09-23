;this program read the files with the probability distribution in temperature
;for each grain compositions and calculates the infrared emissivity in
;W/m^2/Sr/AA 
pro luminosity_ec_template,check 

; input:
; model - dust model, e.g. 'wd01'
; qyear - e.g.'06'
; tau - total B band face-on central opacity
; sfr - star formation rate in Msolar/yr
; old - scaling factor for the optical emission
; bd  - bulge-to-disk ratio
; check - a logical switch if check='yes' then the ISRF is considered;
;                          if check='no' then the radiation fields in the 
;                                        galaxy are considered


;subroutines
;kappagraincomp_H.pro
;planck_l.pro

start_time = systime(/seconds)
common scaling
;common dirdef,rootdir
close,/all

;ss=''
;urad_dir=''
;emission_dir=''
;OPENR, dirunit, 'mydirectories.in', /GET_LUN
;	READF, dirunit, ss
;	READF, dirunit, urad_dir
;	READF, dirunit, ss
;	READF, dirunit, emission_dir
;FREE_LUN, dirunit

dir = emission_dir+'outdata/'
dir2 = emission_dir+'outdata_temp/'
dir3 = emission_dir+'/outdata_lum/'
dir4 = urad_dir+'indata/'

filter = 4430.     ;B band wavelength in AA
dim=800		   ;dimension of the wavelength vector
zeichen=''
Zahl=0
T_min=0.D
T_max=0.D
pass1=0.D
pass2=0.D
pass3=0.D
inorm = 0.D0
iafinal = dblarr(dim)
kappa = dblarr(dim)
newdim_positionsr = 0L ;number of positions in r
newdim_positionsz = 0L ;number of positions in z
newdim_positions = 0L ;total number of positions rxz in the galaxy


;produce the wavelength vector between 10^3 and 5x10^7 Angstroem, with equal 
;bins in a logarithmic scale
a1=Alog10(1.0D3)
a2=ALOG10(5.D7)
a3=(a2-a1)/(dim-1.)
lambda=10.D^(findgen(dim)*a3+a1)

;define some strings for file names
suv=strcompress(string(round(sfr*100)),/remove_all)
sbd=strcompress(string(round(bd*100)),/remove_all)
sold=strcompress(string(round(old*100)),/remove_all)
stau=strcompress(string(round(tau*10)),/remove_all)
;stau='86096'

;read file with geometry
ss = ''

idisk1 = 0L
idisk2 = 0L

filename = 'geometry.in'
name = dir4+filename
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

;read the file with the dust model parameters and grain size distribution
if model eq 'wd01' then namer1 = dir+'grain_sizeswd01_q'+qyear+'.dat'
if model eq 'lmc1' then namer1 = dir+'grain_sizeslmc1_q'+qyear+'.dat'
if model eq 'wd01_c60' then namer1 = dir+'grain_sizeswd01_q'+qyear+'.dat'
if model eq 'wd01_c50' then namer1 = dir+'grain_sizeswd01_q'+qyear+'.dat'
if model eq 'wd01_c40' then namer1 = dir+'grain_sizeswd01_q'+qyear+'.dat'
if model eq 'wd01_c30' then namer1 = dir+'grain_sizeswd01_q'+qyear+'.dat'
openr,unit11,namer1,/get_lun
ss = ' '
dim_size = 0L
dim_comp=0L
dim_Draine=0L
ig_amin=0L
ig_amax=0L
readf,unit11,ss
readf,unit11,dim_comp
comp_array = strarr(dim_comp)
readf,unit11,ss
for ii = 0, dim_comp-1 do begin
    readf, unit11, ss
    comp_array[ii] = ss
endfor

comp = ''
qcomp_year = ''
for j=0L,dim_comp-1 do begin	; begin loop in composition
readf,unit11,ss                 ;
readf,unit11,comp               ;grain composition
readf,unit11,ss                 ;
readf,unit11,qcomp_year         ;model for Q (year when updated)
readf,unit11,ss                 ;
readf,unit11,dim_Draine         ;number of grain sizes tabulated by Draine
readf,unit11,ss                 ;
readf,unit11,amin               ;minimum grain size in micron
readf,unit11,ss                 ;
readf,unit11,amax               ;maximum grain size in micron
readf,unit11,ss                 ;
readf,unit11,ig_amin            ;index of amin in Draine sizes
readf,unit11,ss                 ;
readf,unit11,ig_amax            ;index of amax in Draine sizes
readf,unit11,ss                 ;
readf,unit11,s_grain            ;grain density in g/cm^3
readf,unit11,ss                 ;
readf,unit11,dim_size           ;number of grain sizes
readf,unit11,ss                 ;
readf,unit11,weight             ;
readf,unit11,ss                 ;
readf,unit11,ss                 ;
readf,unit11,ss                 ;
readf,unit11,ss                 ;
a = dblarr(dim_size)           ;grain sizes
na = dblarr(dim_size)          ;grain size distribution

for ii = 0L,dim_size-1 do begin
	readf,unit11,format='(2(e13.6))',pass1,pass2
	a[ii] = pass1 ;cm
	na[ii] = pass2 ;in cm^-1 H-1
endfor

kappagraincomp_H,comp,qcomp_year,dim_Draine,amin,amax,ig_amin,ig_amax,$
	dim_size,a,na,lambda,dim,kappa_ab,kappa_sc,kappa_ph,g,Q_abs
kappa = kappa_ab+kappa_sc

da = a[1:dim_size-1]-a[0:dim_size-2] ;in cm

;derive the extinction cross section in the B band
minlambda = min(lambda)
maxlambda = max(lambda)
if filter lt maxlambda and filter gt minlambda then begin 
	kappa_b = interpol(kappa,lambda,filter)
endif else begin
	print, 'filter outside wavelength range; program stops'
	goto, mark3
endelse


;define a wavelength array to cover only the infrared part

iwave = where(lambda ge 1.d+4 and lambda le 5.d+7, nwave)
if nwave gt 0 then begin
	nlambda = dblarr(nwave)
	nlambda = lambda[iwave] 
endif else begin
	print,'input wavelength does not contain the IR range; program stops'
	goto, mark3
endelse 
dim_wave = nwave
Q_abs_new = dblarr(dim_size,dim_wave)
Q_abs_new[*,*] = Q_abs[*,iwave]
 

 if check eq 'no' then begin
  ;read the file with the temperature distributions
  ;name of the file to read the PT curves
  Name='PT_'+comp+'_'+model+'_'+qyear+$
            '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+'.dat'
  ;name of the file to write the infrared
  ;emissivity (integrated over size distribution) 
  Namew='lum_'+comp+'_'+model+'_q'+qyear+$
              '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+'.dat'
  ;name of the file to write the brightness of each grain size
  ;Namew1='lum_grain_'+comp+'_'+model+'_q'+qyear+$
  ;       '_t'+stau+'_s'+suv+'_o'+sold+'_bd'+sbd+'.dat'
  ;name of the file to write the energy emitted by each grain
  Namew2='energy_emitted_grain_'+comp+'_'+model+'_q'+qyear+$
               '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+'.dat'
  ;name of the file to write the energy emitted and absorbed integrated over
  ;grain size
  Namew3='energy_'+comp+'_'+model+'_q'+qyear+$
                '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+'.dat'
  ;name of the file to read the energy absorbed by each grain
  ;needed to check energy conservation
  Name1='energy_absorbed_grain_'+comp+'_'+model+'_'+qyear+$
             '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+'.dat' 
 endif

 if check eq 'yes' then begin
  ;read the file with the temperature distributions
  ;name of the file to read with the PT curves
  Name='PT_'+comp+'_'+model+'_'+qyear+'_isrf.dat'
  ;name of the file to write with the infrared
  ;emissivity (integrated over size distribution) 
  Namew='lum_'+comp+'_'+model+'_q'+qyear+'_isrf.dat'
  ;name of the file to write with the brightness of each grain size
  ;Namew1='lum_grain_'+comp+'_'+model+'_q'+qyear+'_isrf.dat'
  ;name of the file to write with the energy emitted by each grain
  Namew2='energy_emitted_grain_'+comp+'_'+model+'_q'+qyear+'_isrf.dat'
  ;name of the file to write the energy emitted integrated over grain size
  Namew3='energy_'+comp+'_'+model+'_q'+qyear+'_isrf.dat'
  ;name of the file to read the energy absorbed by each grain
  ;needed to check energy conservation
  Name1='energy_absorbed_grain_'+comp+'_'+model+'_'+qyear+'_isrf.dat' 
 endif

  openr,unit1,dir2+Name,/get_lun
  openw, unit2, dir3+Namew,/get_lun
  ;openw,unit3,dir3+Namew1,/get_lun
  openw,unit4,dir3+Namew2,/get_lun
  openw,unit6,dir3+Namew3,/get_lun
  openr,unit5,dir2+Name1,/get_lun
  print,'write ',dir3+Namew
  ;print,'write ',dir3+Namew1
  print,'write ',dir3+Namew2
  print,'write ',dir3+Namew3

  ss=''

  ;read the header of the input file
  ;the first line contain the input parameters
  readf, unit1, ss
  readf, unit1, newtau
  readf, unit1, ss
  readf, unit1, newsfr
  readf, unit1, ss
  readf, unit1, newold
  readf, unit1, ss
  readf, unit1, newbd
  ;check that the input parameters read from the file are the same as those
  ;inputed to the program (since we have a multiple definition of parameters)
  if double(newtau)-tau ge 0.001 then begin
     print, 'unexpected tau in PT file; program stops'
     goto, mark3
  endif
  if double(newsfr)-sfr ge 0.001 then begin
     print, 'unexpected sfr in PT file; program stops'
     goto, mark3
  endif
  if double(newbd)-bd ge 0.001 then begin
     print, 'unexpected bd in PT file; program stops'
     goto, mark3
  endif
  if double(newold)-old ge 0.001 then begin
     print, 'unexpected old in PT file; program stops'
     goto, mark3
  endif

  readf, unit1, ss
  readf, unit1, newtau1
  if j eq 0L then begin 
     tau1 = newtau1 
  endif else begin
     if double(newtau1)-tau1 ge 0.001 then begin
        print, 'unexpected tau1; program stops'
        goto, mark3
     endif
  endelse 
  readf, unit1, ss
  readf, unit1, newhd ;in pc
  if j eq 0L then begin 
     hd = newhd 
  endif else begin
     if double(newhd)-hd ge 0.001 then begin
        print, 'unexpected hd; program stops'
        goto, mark3
     endif
  endelse 
  readf, unit1, ss
  readf, unit1, newzd ;in pc
  if j eq 0L then begin 
     zd = newzd 
  endif else begin
     if double(newzd)-zd ge 0.001 then begin
        print, 'unexpected zd; program stops'
        goto, mark3
     endif
  endelse 
  readf, unit1, ss
  readf, unit1, hdin
  readf, unit1, ss
  readf, unit1, zdin
  readf, unit1, ss
  readf, unit1, hdsolar
  readf, unit1, ss
  readf, unit1, zdsolar
  readf, unit1, ss
  readf, unit1, newtau2
  if j eq 0L then begin 
     tau2 = newtau2 
  endif else begin
     if double(newtau2)-tau2 ge 0.001 then begin
        print, 'unexpected tau2; program stops'
        goto, mark3
     endif
  endelse 
  readf, unit1, ss
  readf, unit1, newhd1 ;in pc
  if j eq 0L then begin 
     hd1 = newhd1 
  endif else begin
     if double(newhd1)-hd1 ge 0.001 then begin
        print, 'unexpected hd1; program stops'
        goto, mark3
     endif
  endelse 
  readf, unit1, ss
  readf, unit1, newzd1 ;in pc
  if j eq 0L then begin 
     zd1 = newzd1 
  endif else begin
     if double(newzd1)-zd1 ge 0.001 then begin
        print, 'unexpected zd1; program stops'
        goto, mark3
     endif
  endelse 
  readf, unit1, ss
  readf, unit1, hd1in
  readf, unit1, ss
  readf, unit1, zd1in
  readf, unit1, ss
  readf, unit1, hd1solar
  readf, unit1, ss
  readf, unit1, zd1solar
  readf, unit1, ss
  readf, unit1, hs
  readf, unit1, ss
  readf, unit1, zs
  readf, unit1, ss
  readf, unit1, hsin
  readf, unit1, ss
  readf, unit1, zsin
  readf, unit1, ss
  readf, unit1, hssolar
  readf, unit1, ss
  readf, unit1, zssolar
  readf, unit1, ss
  readf, unit1, hs1
  readf, unit1, ss
  readf, unit1, zs1
  readf, unit1, ss
  readf, unit1, hs1in
  readf, unit1, ss
  readf, unit1, zs1in
  readf, unit1, ss
  readf, unit1, hs1solar
  readf, unit1, ss
  readf, unit1, zs1solar
  readf, unit1, ss
  readf, unit1, rtrun
  readf, unit1, ss
  readf, unit1, rtrun1
  readf, unit1, ss
  readf, unit1, xis0
  readf, unit1, ss
  readf, unit1, xis1
  readf, unit1, ss
  readf, unit1, xid0
  readf, unit1, ss
  readf, unit1, xid1
  readf, unit1, ss
  readf, unit1, idisk1
  readf, unit1, ss
  readf, unit1, idisk2
  readf, unit1, ss
  readf, unit1, newdim_positionsr
  if j eq 0L then begin 
     dim_positionsr = newdim_positionsr
  endif else begin
     if newdim_positionsr ne dim_positionsr then begin
        print, 'unexpected dim_positionsr; program stops'
        goto, mark3
     endif
  endelse 
  readf, unit1, ss
  readf, unit1, newdim_positionsz
  if j eq 0L then begin 
     dim_positionsz = newdim_positionsz 
  endif else begin
     if newdim_positionsz ne dim_positionsz then begin
        print, 'unexpected dim_positionsz; program stops'
        goto, mark3
     endif
  endelse 
  readf, unit1, ss
  readf,unit1, ss
  readf,unit1,newdim_positions
  if j eq 0L then begin 
     dim_positions = newdim_positions 
  endif else begin
     if newdim_positions ne dim_positions then begin
        print, 'unexpected dim_positions; program stops'
        goto, mark3
     endif
  endelse 
  readf,unit1,ss
  readf,unit1,ss
  nsize=0L
  readf,unit1,nsize
  if nsize ne dim_size then begin
     print,'unexpected number of grain sizes; program stops'
     goto, mark3
  endif 
  readf,unit1,ss

  printf, unit2, 'central face-on tau B'
  printf, unit2, tau
  printf, unit2, 'sfr'
  printf, unit2, sfr
  printf, unit2, 'old'
  printf, unit2, old
  printf, unit2, 'bd'
  printf, unit2, bd
  printf, unit2, 'tau1'
  printf, unit2, tau1
  printf, unit2, 'hd [pc]'
  printf, unit2, hd
  printf, unit2, 'zd [pc]'
  printf, unit2, zd
  printf, unit2, 'hdin [pc]'
  printf, unit2, hdin
  printf, unit2, 'zdin [pc]'
  printf, unit2, zdin
  printf, unit2, 'hdsolar [pc]'
  printf, unit2, hdsolar
  printf, unit2, 'zdsolar [pc]'
  printf, unit2, zdsolar
  printf, unit2, 'tau2'
  printf, unit2, tau2
  printf, unit2, 'hd1 [pc]'
  printf, unit2, hd1
  printf, unit2, 'zd1 [pc]'
  printf, unit2, zd1
  printf, unit2, 'hd1in [pc]'
  printf, unit2, hd1in
  printf, unit2, 'zd1in [pc]'
  printf, unit2, zd1in
  printf, unit2, 'hd1solar [pc]'
  printf, unit2, hd1solar
  printf, unit2, 'zd1solar [pc]'
  printf, unit2, zd1solar
  printf, unit2, 'hs [pc]'
  printf, unit2, hs
  printf, unit2, 'zs [pc]'
  printf, unit2, zs
  printf, unit2, 'hsin [pc]'
  printf, unit2, hsin
  printf, unit2, 'zsin [pc]'
  printf, unit2, zsin
  printf, unit2, 'hssolar [pc]'
  printf, unit2, hssolar
  printf, unit2, 'zssolar [pc]'
  printf, unit2, zssolar
  printf, unit2, 'hs1 [pc]'
  printf, unit2, hs1
  printf, unit2, 'zs1 [pc]'
  printf, unit2, zs1
  printf, unit2, 'hs1in [pc]'
  printf, unit2, hs1in
  printf, unit2, 'zs1in [pc]'
  printf, unit2, zs1in
  printf, unit2, 'hs1solar [pc]'
  printf, unit2, hs1solar
  printf, unit2, 'zs1solar [pc]'
  printf, unit2, zs1solar
  printf, unit2, 'rtrun [pc]'
  printf, unit2, rtrun
  printf, unit2, 'rtrun1 [pc]'
  printf, unit2, rtrun1
  printf, unit2, 'xis0'
  printf, unit2, xis0
  printf, unit2, 'xis1'
  printf, unit2, xis1
  printf, unit2, 'xid0'
  printf, unit2, xid0
  printf, unit2, 'xid1'
  printf, unit2, xid1
  printf, unit2, 'idisk1'
  printf, unit2, idisk1
  printf, unit2, 'idisk2'
  printf, unit2, idisk2
  printf, unit2, 'number of r positions'
  printf, unit2, dim_positionsr
  printf, unit2, 'number of z positions'
  printf, unit2, dim_positionsz
  printf, unit2, ss
  printf, unit2,  'number of positions in the galaxy:'
  printf, unit2, dim_positions
  printf,unit2
  printf, unit2, 'composition'
  printf, unit2, comp
  printf,unit2
  printf, unit2, 'weight'
  printf, unit2, weight
  printf,unit2
  printf, unit2, 'Number of wavelengths:'
  printf, unit2, dim_wave
  printf,unit2
  printf, unit2, 'C_ext(B) [cm^2/H]'
  printf, unit2, kappa_b 
  ;printf, unit3, 'central face-on tau B'
  ;printf, unit3, tau
  ;printf, unit3, 'sfr'
  ;printf, unit3, sfr
  ;printf, unit3, 'old'
  ;printf, unit3, old
  ;printf, unit3, 'bd'
  ;printf, unit3, bd
  ;printf,unit3
  ;printf,unit3,'composition'
  ;printf,unit3,comp
  ;printf,unit3
  ;printf,unit3,'weight'
  ;printf,unit3,weight
  ;printf,unit3
  ;printf,unit3,'Number of positions in the galaxy:'
  ;printf,unit3,dim_positions
  ;printf,unit3
  ;printf,unit3,'Number of wavelengths:'
  ;printf,unit3,dim
  ;printf,unit3
  ;printf,unit3,'Number of grain sizes'
  ;printf,unit3,dim_size
  ;printf,unit3

  printf, unit4, 'central face-on tau B'
  printf, unit4, tau
  printf, unit4, 'sfr'
  printf, unit4, sfr
  printf, unit4, 'old'
  printf, unit4, old
  printf, unit4, 'bd'
  printf, unit4, bd
  printf,unit4
  printf,unit4,'composition'
  printf,unit4,comp
  printf,unit4
  printf,unit4,'weight'
  printf,unit4,weight
  printf,unit4
  printf,unit4,'Number of positions in the galaxy:'
  printf,unit4,dim_positions
  printf,unit4
  printf,unit4,'Number of wavelengths:'
  printf,unit4,dim_wave
  printf,unit4
  printf,unit4,'Number of grain sizes'
  printf,unit4,dim_size
  printf,unit4

  printf, unit6, 'central face-on tau B'
  printf, unit6, tau
  printf, unit6, 'sfr'
  printf, unit6, sfr
  printf, unit6, 'old'
  printf, unit6, old
  printf, unit6, 'bd'
  printf, unit6, bd
  printf,unit6
  printf,unit6,'composition'
  printf,unit6,comp
  printf,unit6
  printf,unit6,'weight'
  printf,unit6,weight
  printf,unit6
  printf,unit6,'Number of positions in the galaxy:'
  printf,unit6,dim_positions
  printf,unit6
  printf,unit6,'Number of wavelengths:'
  printf,unit6,dim_wave
  printf,unit6
  printf,unit6,'Number of grain sizes'
  printf,unit6,dim_size
  printf,unit6
  printf,unit6, '      r            z         energy_abs      energy_emitted'
  printf,unit6, '      pc           pc        W/H             W/H'

  readf, unit5, ss
  readf, unit5, newtau
  readf, unit5, ss
  readf, unit5, newsfr
  readf, unit5, ss
  readf, unit5, newold
  readf, unit5, ss
  readf, unit5, newbd
  ;check that the input parameters read from the file are the same as those
  ;inputed to the program (since we have a multiple definition of parameters)
  if double(newtau)-tau ge 0.001 then begin
     print, 'unexpected tau in energy absorbed file; program stops'
     goto, mark3
  endif
  if double(newsfr)-sfr ge 0.001 then begin
     print, 'unexpected sfr in energy absorbed file; program stops'
     goto, mark3
  endif
  if double(newbd)-bd ge 0.001 then begin
     print, 'unexpected bd in energy absorbed file; program stops'
     goto, mark3
  endif
  if double(newold)-old ge 0.001 then begin
     print, 'unexpected old in energy absorbed file; program stops'
     goto, mark3
  endif
  readf,unit5, ss
  readf,unit5,dim_pos
  if dim_pos ne dim_positions then begin
     print, 'unexpected number of positions in the galaxy; program stops'
     goto, mark3
  endif
  readf,unit5, ss
  readf,unit5, ss
  readf,unit5,dim_s
  if dim_s ne dim_size then begin
     print, 'unexpected number of grain sizes; program stops'
     goto, mark3
  endif
  readf,unit5, ss

  ig = lonarr(dim_size)
  aa = dblarr(dim_size)
  aaa = dblarr(dim_size)
  ta=dblarr(dim_size,dim_wave)
  energy_abs = dblarr(dim_size)
  ia_grain=dblarr(dim_size,dim_wave)

  for k=0L,dim_positions-1 do begin  ;begin loop in position
      readf,unit1,ss
      readf,unit1, rr, zz

      readf,unit5,ss
      readf,unit5, rrr, zzz
      if abs(rrr-rr) or abs(zzz-zz) gt 0.1 then begin
         print, 'unexpected position in the galaxy; program stops'
         goto, mark3
      endif
      readf,unit5,ss
      readf,unit5,ss
      readf,unit5,ss 

      for ii=0L,dim_size-1 do begin  ;loop for different grain sizes
        jj = dim_size - 1 - ii
        readf,unit1,ss
        readf,unit1,pass
        ig[jj] = pass
        readf,unit1,ss
        readf,unit1,ss
        readf,unit1,ss
        readf,unit1,pass
        aa[jj] = pass ;in micron
        if abs(aa[jj]*1.d-4 - a[jj]) gt 0.1 then begin
           print, 'unexpected grain size in the PT file; program stops'
           goto, mark3
        endif
        readf,unit1,ss
        readf,unit1,T_min,T_max
        readf,unit1,ss
        readf,unit1,Zahl
        readf,unit1,ss
        readf,unit1,ss
        T_m=dblarr(Zahl-1)   
        prob=dblarr(Zahl-1)
        bb=dblarr(dim_wave,Zahl-1)
        for kk=0L,Zahl-2 do begin
            readf,unit1,pass1,pass2,pass3
            T_m(kk)=pass1
            prob(kk)=pass2
            bb(*,kk) = planck_l(nlambda,T_m(kk))
        endfor

        ;make the integration over temperature
        for l=0L,dim_wave-1 do begin
         f=dblarr(Zahl-1)
         f[*]=bb[l,*]*prob[*]
         integ = total(f)
         ta[ig[jj],l] = integ*Q_abs_new[ig[jj],l] ;in W/m^2/A/sr
         ia_grain[ig[jj],l] = !Dpi * ta[ig[jj],l] * a[ig[jj]]^2  
        endfor

        ;read the energy absorbed
         readf,unit5,pass1,pass2
         aaa[jj] = pass1  ;in micron
         if abs(aaa[jj]-aa[jj]) gt 0.1 then begin
            print, 'unexpected grain size; program stops'
            goto, mark3
         endif         
         energy_abs[jj] = pass2
      endfor ;end loop in grain size ii orjj 

      ia = dblarr(dim_size)
      ;make the integration over the size distribution
      for l=0L,dim_wave-1 do begin
        ia(*) = !Dpi * na * ta[*,l] * a^2 
        iia = 0.5*(ia[1:dim_size-1]+ia[0:dim_size-2])
        y = iia*da
        iafinal[l] = TOTAL(y)      
      endfor
      
      ;the 1.D-4 comes from the transformation of a^2 from cm to m
      ;the 1.d+7 comes from the transformation from W in erg/sec
      iafinal = 1.D-4 * iafinal * 1.d+7 ; in W/AA/sr/H
      ia_grain = 1.D-4 * ia_grain ;in W/AA/sr


      ;check the energy emitted by a single grain

      dl = nlambda[1:dim_wave-1]-nlambda[0:dim_wave-2]
      energy_emitted = dblarr(dim_size)
      for ii = 0L, dim_size-1 do begin
          energ = 0.5*(ia_grain[ii,1:dim_wave-1]+ia_grain[ii,0:dim_wave-2])
          energy_emitted[ii] = 4. * !Dpi * total(energ*dl) ; in W or J/s
          ;print,'a, energy emitted [J/s]: ',$ 
          ;      a[ii],energy_emitted[ii]
      endfor


      ;integrate the energy emitted over the grain size distribution
       iee = na * energy_emitted
       ieem = 0.5*(iee[1:dim_size-1]+iee[0:dim_size-2])
       int_energy_emitted = total(ieem*da) ;in W/H

     ;integrate the energy absorbed over the grain size distribution
       iea = na * energy_abs
       ieam = 0.5*(iea[1:dim_size-1]+iea[0:dim_size-2])
       int_energy_abs = total(ieam*da) ;in W/H

       
      ;check energy conservation
      ;energy_emitted_m = 0.5*(energy_emitted[0:dim_size-2]+$
      ;                   energy_emitted[1:dim_size-1])    
      ;energy_abs_m = 0.5*(energy_abs[0:dim_size-2]+$
      ;                   energy_abs[1:dim_size-1])
      ;print, 'energy_emitted_m', energy_emitted_m
      ;print, 'energy_abs_m', energy_abs_m
      ;total_energy_emitted = total(energy_emitted_m * da) 
      ;total_energy_abs = total(energy_abs_m * da) * 1.d
      ;print, 'total energy emitted', int_energy_emitted
      ;print, 'total energy absorbed', int_energy_abs
      dif = abs(int_energy_abs - int_energy_emitted)/int_energy_abs
      if dif gt 0.025 then begin
		print, ' '
		print, 'comp, rr, zz:  ', comp, rr, zz
 		print, '**** energy conservation worse than 2.5 percent **** ', dif 
      endif

      printf, unit2
      printf, unit2, '; R z in pc'
      printf, unit2,  rr, zz
      printf, unit2
      printf, unit2, '       lambda(A)   I(erg s^-1 sr^-1 H^-1 A^-1)' 
      for i=0,dim_wave-1 do begin
      	printf, unit2, nlambda[i],iafinal[i]
      endfor
       
      ;printf,unit3
      ;printf,unit3,'; R z in pc'
      ;printf,unit3, rr, zz
      ;for ii = 0, dim_size-1 do begin
      ;    printf,unit3
      ;    printf,unit3,'grain size in cm'
      ;    printf,unit3,a[ii]
      ;    printf,unit3
      ;    printf,unit3,'       lambda(A)   I(W m^-2 sr^-1 A^-1)'  
      ;    for i=0,dim-1 do begin
      ;        printf,unit3,nlambda[i],ta[ii,i]
      ;    endfor
      ;endfor

      printf,unit4
      printf,unit4,'; R z in pc'
      printf,unit4, rr, zz
      printf,unit4
      printf,unit4,';energy conservation accuracy'
      printf,unit4, dif
      printf,unit4,'grain size      energy emitted'
      printf,unit4,'cm              W'
      for ii = 0, dim_size-1 do begin
		printf,unit4,a[ii],energy_emitted[ii]
      endfor
      printf,unit6
      printf,unit6, rr, zz, int_energy_abs, int_energy_emitted
  endfor ;loop in position


free_lun, unit1
free_lun, unit2
;free_lun, unit3
free_lun, unit4
free_lun, unit5
free_lun, unit6
endfor ;end loop in composition
mark3:
free_lun, unit11
print, 'DONE: luminosity_ec_template.pro ('+strtrim(ceil(systime(/seconds)-start_time),1)+' s)'
end
