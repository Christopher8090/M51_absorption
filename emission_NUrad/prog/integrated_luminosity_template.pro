;this program reads the infrared emissivity in erg/s/sr/H/A at each 
;R,z position and resamples this on a regular grid using cubic spline
;interpolation. Then it integrates the emission and derives the total 
;luminosity/Hz of the galaxy. 
;11/10/17 updated by JJT to include extra components 
;          geometrical parameters for extra components are multiplied by
;          10d^3. as input parameters for original components are read in
;          twice. Once in [kpc], and then overwritten in [pc].
;27/11/17 JJT updated to deal with inner truncation radii
;10/07/18 BTC updated to be able to disable extra components.
;29/11/18 JJT updated to allow calculation of heating due to different stellar components 

pro integrated_luminosity_template,model,qyear,tau,sfr,sfr4,sfr6,sfr7,old,old3,old5,bd,scaabs,dr_req,dz_req,swdisk3,swdisk4,swdisk5,swdisk6,swheatr
;integrated_luminosity_template,'wd01','06',3.8,2.28,0.2,0.09,'sca',100.,25.,'yes','yes','yes','yes' ;call for MW
;integrated_luminosity_template,'wd01','06',25.4,7.8,5.8,1.,0.,2.5,1.45,0.7,0.01,'abs',100.,25.,'yes','yes','yes','yes','no'



;input:
; model - dust model, e.g. 'wd01'
; qyear - e.g.'06'
; tau - total B band face-on central opacity
; sfr/sfr4/sfr6 - star formation rate in Msolar/yr
; old/old3/old5 - scaling factor for the optical emission
; bd  - bulge-to-disk ratio
; dr_req,dz_req - the required resampling; the actual resampling may be a bit
;                 finer than this
;                 the 100., 25. pc was tested to give the best compromise 
;                 between accuracy and speed
;sw*  -  switches for various components or heating source contribution calculation values of 'yes' or 'no' only
;		-heatr tells program if you want to calculate dust emission due to old stellar population only.

;subroutines:
;regrid.pro

start_time = systime(/seconds)
common scaling
close,/all

COMMON gridinfo, lambda, rrr, zzz, dim

cl = 2.99792458d8                ;speed of light in m/sec
newdim_positionsr = 0L  ;number of r positions
newdim_positionsz = 0L  ;number of z position
newdim_positions = 0L   ;total number of rxz positions in the galaxy 
newdim = 0L            ;number of wavelengths
radius = 0.         ;r position in the galaxy in units of pc
height = 0.         ;z position in the galaxy in units of pc
ss=''
dim_comp = 0L       ;number of grain compositions
dim_size = 0L       ;number of grain sizes

;define some strings for file names
suv=strcompress(string(round(sfr*100)),/remove_all)
sbd=strcompress(string(round(bd*100)),/remove_all)
sold=strcompress(string(round(old*100)),/remove_all)
stau=strcompress(string(round(tau*10)),/remove_all)

dir = emission_dir+'outdata/'
dir3 = emission_dir+'outdata_lum/'
dir4 = emission_dir+'outdata_intlum/'
dir2= urad_dir+'indata/'

;set up for new identifier for grids
ratoldident=''
if (swheatr eq 'yes') then ratoldident='_old'

;read file with geometry
ss = ''

idisk1 = 0L
idisk2 = 0L
; JJT 16/10/17
if (swdisk3 eq 'yes') then idisk3 = 0L
if (swdisk4 eq 'yes') then idisk4 = 0L
if (swdisk5 eq 'yes') then idisk5 = 0L
if (swdisk6 eq 'yes') then idisk6 = 0L

filename = 'geometry.in'

name = dir2+filename
openr,unit,name,/get_lun
print, 'read ', name

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
h_ir36disk=h_ir36disk*10d^3.
readf, unit, ss
readf, unit, h_ir45disk
h_ir45disk=h_ir45disk*10d^3.
readf, unit, ss
readf, unit, h_ir58disk
h_ir58disk=h_ir58disk*10d^3
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
;read inner component
readf, unit, ss
readf, unit, ss
readf, unit, tau3
readf, unit, ss
readf, unit, tau4
readf, unit, ss
readf, unit, hd3
hd3=hd3*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zd3
zd3=zd3*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hd3in
hd3in=hd3in*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zd3in
zd3in=zd3in*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hd3solar
hd3solar=hd3solar*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zd3solar
zd3solar=zd3solar*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hd4
hd4=hd4*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zd4
zd4=zd4*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hd4in
hd4in=hd4in*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zd4in
zd4in=zd4in*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hd4solar
hd4solar=hd4solar*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zd4solar
zd4solar=zd4solar*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, h_bdisk3
h_bdisk3=h_bdisk3*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, h_vdisk3
h_vdisk3=h_vdisk3*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, h_idisk3
h_idisk3=h_idisk3*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, h_jdisk3
h_jdisk3=h_jdisk3*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, h_kdisk3
h_kdisk3=h_kdisk3*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, h_ir36disk3
h_ir36disk3=h_ir36disk3*10d^3.
readf, unit, ss
readf, unit, h_ir45disk3
h_ir45disk3=h_ir45disk3*10d^3.
readf, unit, ss
readf, unit, h_ir58disk3
h_ir58disk3=h_ir58disk3*10d^3.
readf, unit, ss
readf, unit, zs3
zs3=zs3*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hs3in
hs3in=hs3in*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zs3in
zs3in=zs3in*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hs3solar
hs3solar=hs3solar*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zs3solar
zs3solar=zs3solar*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hs4
hs4=hs4*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zs4
zs4=zs4*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hs4in
hs4in=hs4in*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zs4in
zs4in=zs4in*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hs4solar
hs4solar=hs4solar*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zs4solar
zs4solar=zs4solar*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, rtrun3
readf, unit, ss
readf, unit, sharp3
readf, unit, ss
readf, unit, rtrund3
readf, unit, ss
readf, unit, sharpd3
readf, unit, ss
readf, unit, rtrun4
readf, unit, ss
readf, unit, sharp4
readf, unit, ss
readf, unit, rtrund4
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
;read outer component
readf, unit, ss
readf, unit, ss
readf, unit, tau5
readf, unit, ss
readf, unit, tau6
readf, unit, ss
readf, unit, hd5
hd5=hd5*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zd5
zd5=zd5*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hd5in
hd5in=hd5in*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zd5in
zd5in=zd5in*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hd5solar
hd5solar=hd5solar*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zd5solar
zd5solar=zd5solar*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hd6
hd6=hd6*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zd6
zd6=zd6*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hd6in
hd6in=hd6in*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zd6in
zd6in=zd6in*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hd6solar
hd6solar=hd6solar*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zd6solar
zd6solar=zd6solar*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, h_bdisk5
h_bdisk5=h_bdisk5*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, h_vdisk5
h_vdisk5=h_vdisk5*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, h_idisk5
h_idisk5=h_idisk5*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, h_jdisk5
h_jdisk5=h_jdisk5*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, h_kdisk5
h_kdisk5=h_kdisk5*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, h_ir36disk5
h_ir36disk5=h_ir36disk5*10d^3.
readf, unit, ss
readf, unit, h_ir45disk5
h_ir45disk5=h_ir45disk5*10d^3.
readf, unit, ss
readf, unit, h_ir58disk5
h_ir58disk5=h_ir58disk5*10d^3.
readf, unit, ss
readf, unit, zs5
zs5=zs5*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hs5in
hs5in=hs5in*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zs5in
zs5in=zs5in*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hs5solar
hs5solar=hs5solar*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zs5solar
zs5solar=zs5solar*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hs6
hs6=hs6*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zs6
zs6=zs6*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hs6in
hs6in=hs6in*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zs6in
zs6in=zs6in*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, hs6solar
hs6solar=hs6solar*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, zs6solar
zs6solar=zs6solar*10d^3.;[pc] conversion
readf, unit, ss
readf, unit, rtrun5
readf, unit, ss
readf, unit, sharp5
readf, unit, ss
readf, unit, rtrund5
readf, unit, ss
readf, unit, sharpd5
readf, unit, ss
readf, unit, rtrun6
readf, unit, ss
readf, unit, sharp6
readf, unit, ss
readf, unit, rtrund6
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
;read inner truncation radius
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

;must multiply rtrun3/4/5/6 by 10d^3 as rtrun and rtrun1
;are read in again from a different file in [pc] rather than [kpc]
rtrun3=rtrun3*10d^3.
rtrun4=rtrun4*10d^3.
rtrun5=rtrun5*10d^3.
rtrun6=rtrun6*10d^3.
rtrund=rtrund*10d^3.
rtrund1=rtrund1*10d^3.
rtrund3=rtrund3*10d^3.
rtrund4=rtrund4*10d^3.
rtrund5=rtrund5*10d^3.
rtrund6=rtrund6*10d^3.
hstin=hstin*10d^3.
hdtin=hdtin*10d^3.
hs1tin=hs1tin*10d^3.
hd1tin=hd1tin*10d^3.
hs3tin=hs3tin*10d^3.
hd3tin=hd3tin*10d^3.
hs4tin=hs4tin*10d^3.
hd4tin=hd4tin*10d^3.
hs5tin=hs5tin*10d^3.
hd5tin=hd5tin*10d^3.
hs6tin=hs6tin*10d^3.
hd6tin=hd6tin*10d^3.

;read the information about grain composition
if model eq 'wd01' then namer1=rootdir+dir+'grain_sizeswd01_q'+qyear+'.dat'
if model eq 'lmc1' then namer1=rootdir+dir+'grain_sizeslmc1_q'+qyear+'.dat'
if model eq 'wd01_c60' then namer1=rootdir+dir+'grain_sizeswd01_q'+qyear+'.dat'
if model eq 'wd01_c50' then namer1=rootdir+dir+'grain_sizeswd01_q'+qyear+'.dat'
if model eq 'wd01_c40' then namer1=rootdir+dir+'grain_sizeswd01_q'+qyear+'.dat'
if model eq 'wd01_c30' then namer1=rootdir+dir+'grain_sizeswd01_q'+qyear+'.dat'
openr, unit11, namer1,/get_lun
readf, unit11, ss
readf, unit11, dim_comp
comp_array = strarr(dim_comp)
readf, unit11, ss
for ii = 0L, dim_comp-1 do begin
    readf, unit11, ss
    comp_array(ii) = ss
endfor
free_lun,unit11

Cext_b = 0.D
kk = 0L
k = 0L
lu=0.D
ll=0.D
compgrain = ' '

for jj = 0L,dim_comp-1 do begin  ;loop in composition
comp = comp_array[jj]

Name1='lum_'+comp+'_'+model+'_q'+qyear
Name1 = Name1 + '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
            '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
            '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
            '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+'.dat'
Name2='energy_'+comp+'_'+model+'_q'+qyear
Name2 = Name2 + '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
            '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
            '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
            '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+'.dat'
openr, unit1, dir3+Name1,/get_lun
openr,unit5,dir3+Name2,/get_lun
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
if newtau ne tau then begin
	print, 'unexpected tau in lum file; program ;stops'
	goto, mark1
endif
if newsfr ne sfr then begin
	print, 'unexpected sfr in lum file; program ;stops'
	goto, mark1
endif
if newbd ne bd then begin
	print, 'unexpected bd in lum file; program ;stops'
	goto, mark1
endif
if newold ne old then begin
	print, 'unexpected old in lum file; program ;stops'
	goto, mark1
endif
readf, unit1, ss
readf, unit1, newtau1
if jj eq 0L then begin 
	tau1 = newtau1 
endif else begin
	if newtau1 ne tau1 then begin
		print, 'unexpected tau1; program ;stops'
		goto, mark1
	endif
endelse 
readf, unit1, ss
readf, unit1, newhd ;in pc
if jj eq 0L then begin 
	hd = newhd 
endif else begin
	if newhd ne hd then begin
		print, 'unexpected hd1; program ;stops'
		goto, mark1
	endif
endelse 
readf, unit1, ss
readf, unit1, newzd ;in pc
if jj eq 0L then begin 
	zd = newzd 
endif else begin
	if newzd ne zd then begin
		print, 'unexpected zd1; program ;stops'
		goto, mark1
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
if jj eq 0L then begin 
	tau2 = newtau2 
endif else begin
	if newtau2 ne tau2 then begin
		print, 'unexpected tau2; program ;stops'
		goto, mark1
	endif
endelse 
readf, unit1, ss
readf, unit1, newhd1 ;in pc
if jj eq 0L then begin 
	hd1 = newhd1 
endif else begin
	if newhd1 ne hd1 then begin
		print, 'unexpected hd2; program ;stops'
		goto, mark1
	endif
endelse 
readf, unit1, ss
readf, unit1, newzd1 ;in pc
if jj eq 0L then begin 
	zd1 = newzd1 
endif else begin
	if newzd1 ne zd1 then begin
		print, 'unexpected zd2; program ;stops'
		goto, mark1
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
if jj eq 0L then begin 
	dim_positionsr = newdim_positionsr
endif else begin
	if newdim_positionsr ne dim_positionsr then begin
		print, 'unexpected dim_positionsr; program ;stops'
		goto, mark1
	endif
endelse 
readf, unit1, ss
readf, unit1, newdim_positionsz
if jj eq 0L then begin 
	dim_positionsz = newdim_positionsz 
endif else begin
	if newdim_positionsz ne dim_positionsz then begin
		print, 'unexpected dim_positionsz; program ;stops'
		goto, mark1
	endif
endelse 
readf, unit1, ss
readf, unit1,  ss
readf, unit1, newdim_positions
if jj eq 0L then begin 
	dim_positions = newdim_positions 
endif else begin
	if newdim_positions ne dim_positions then begin
		print, 'unexpected dim_positions; program ;stops'
		goto, mark1
	endif
endelse 
readf, unit1, ss
readf, unit1, ss
readf, unit1, compgrain
if compgrain ne comp then begin
	print, 'unexpected grain composition; program ;stops'
	goto, mark1
endif 
readf, unit1, ss
readf, unit1, ss
readf, unit1, weight
readf, unit1, ss
readf, unit1, ss
readf, unit1, newdim
if jj eq 0L then begin 
	dim = newdim 
endif else begin
	if newdim ne dim then begin
		print, 'unexpected wavelength dimension; program ;stops'
		goto, mark1
	endif
endelse  
readf, unit1, ss
readf, unit1, ss
readf, unit1, cb
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
if newtau ne tau then begin
	print, 'unexpected tau in lum file; program ;stops'
	goto, mark1
endif
if newsfr ne sfr then begin
	print, 'unexpected sfr in lum file; program ;stops'
	goto, mark1
endif
if newbd ne bd then begin
	print, 'unexpected bd in lum file; program ;stops'
	goto, mark1
endif
if newold ne old then begin
	print, 'unexpected old in lum file; program ;stops'
	goto, mark1
endif
readf,unit5,ss
readf,unit5,ss
readf,unit5,compgrain
if compgrain ne comp then begin
	print, 'unexpected grain composition; program ;stops'
	goto, mark1
endif 
readf,unit5,ss
readf,unit5,ss
readf,unit5,weight
readf,unit5,ss
readf,unit5,ss
readf,unit5,newdim_positions
if newdim_positions ne dim_positions then begin
	print, 'unexpected number of positions; program ;stops'
	goto, mark1
endif 
readf,unit5,ss
readf,unit5,ss
readf,unit5,newdim
if newdim ne dim then begin
	print, 'unexpected number of wavelengths; program ;stops'
	goto, mark1
endif 
readf,unit5,ss
readf,unit5,ss
readf,unit5,dim_size
readf,unit5,ss
readf,unit5,ss
readf,unit5,ss
 
;input the extinction coefficient in the B band in cm^2/H
Cext_B = Cext_B + cb

;the number density of H atoms in the centre of the galaxy multiplied to the
;extinction cross section in the B band [pc^-1]

nHc1Cext = tau1/(2.*zd)
nHc2Cext = tau2/(2.*zd1)
nHc3Cext = tau3/(2.*zd3)
nHc4Cext = tau4/(2.*zd4)
nHc5Cext = tau5/(2.*zd5)
nHc6Cext = tau6/(2.*zd6)

lum_comp = dblarr(dim_positions,dim)
lum1_comp = dblarr(dim_positions,dim)
lum2_comp = dblarr(dim_positions,dim)
dust1 = dblarr(dim_positions)
dust2 = dblarr(dim_positions)
dust = dblarr(dim_positions)
lambda = dblarr(dim) 

if (swdisk3 eq 'yes') then lum3_comp = dblarr(dim_positions,dim)
if (swdisk4 eq 'yes') then lum4_comp = dblarr(dim_positions,dim)
if (swdisk3 eq 'yes') then dust3 = dblarr(dim_positions)
if (swdisk4 eq 'yes') then dust4 = dblarr(dim_positions)
if (swdisk5 eq 'yes') then lum5_comp = dblarr(dim_positions,dim)
if (swdisk6 eq 'yes') then lum6_comp = dblarr(dim_positions,dim)
if (swdisk5 eq 'yes') then dust5 = dblarr(dim_positions)
if (swdisk6 eq 'yes') then dust6 = dblarr(dim_positions)

if jj eq 0 then begin 
	lum = lum_comp
	lum1 = lum1_comp
	lum2 = lum2_comp
	if (swdisk3 eq 'yes') then lum3 = lum3_comp
	if (swdisk4 eq 'yes') then lum4 = lum4_comp
	if (swdisk5 eq 'yes') then lum5 = lum5_comp
	if (swdisk6 eq 'yes') then lum6 = lum6_comp
endif
rr = dblarr(dim_positions)
zz = dblarr(dim_positions)
zd01 = dblarr(dim_positions)
zd11 = dblarr(dim_positions)
if (swdisk3 eq 'yes') then zd31 = dblarr(dim_positions)
if (swdisk4 eq 'yes') then zd41 = dblarr(dim_positions)
if (swdisk5 eq 'yes') then zd51 = dblarr(dim_positions)
if (swdisk6 eq 'yes') then zd61 = dblarr(dim_positions)

rrr = dblarr(dim_positionsr,dim_positionsz)
zzz = dblarr(dim_positionsr,dim_positionsz)
lumdouble = dblarr(dim_positionsr,dim_positionsz,dim)
lumdouble1 = dblarr(dim_positionsr,dim_positionsz,dim)
lumdouble2 = dblarr(dim_positionsr,dim_positionsz,dim)
energy_abs1 = dblarr(dim_positions)
energy_abs2 = dblarr(dim_positions)

lumdouble3 = dblarr(dim_positionsr,dim_positionsz,dim)
lumdouble4 = dblarr(dim_positionsr,dim_positionsz,dim)
energy_abs3 = dblarr(dim_positions)
energy_abs4 = dblarr(dim_positions)
lumdouble5 = dblarr(dim_positionsr,dim_positionsz,dim)
lumdouble6 = dblarr(dim_positionsr,dim_positionsz,dim)
energy_abs5 = dblarr(dim_positions)
energy_abs6 = dblarr(dim_positions)

energy_abs_comp = dblarr(dim_positions)
if jj eq 0 then energy_abs = energy_abs_comp
energy_abs_double = dblarr(dim_positionsr,dim_positionsz)

if (zdin gt zd) then flared01 = alog10((zdsolar - zd)/(zdin - zd)) else flared01 = 0.
if (hdin gt 0.) then begin
	flared02 = alog10(hdsolar/hdin)
	flared00 = flared01/flared02
endif else flared00 = 0.

if (swdisk3 eq 'yes') then begin
	if (zd3in gt zd3) then flared31 = alog10((zd3solar - zd3)/(zd3in - zd3)) else flared31 = 0.
	if (hd3in gt 0.) then begin
		flared32 = alog10(hd3solar/hd3in)
		flared30 = flared31/flared32
	endif else flared30 = 0.
endif

if (swdisk5 eq 'yes') then begin
	if (zd5in gt zd5) then flared51 = alog10((zd5solar - zd5)/(zd5in - zd5)) else flared51 = 0.
	if (hd5in gt 0.) then begin
		flared52 = alog10(hd5solar/hd5in)
		flared50 = flared51/flared52
	endif else flared50 = 0.
endif

if (zd1in gt zd1) then flared11 = alog10((zd1solar - zd1)/(zd1in - zd1)) else flared11 = 0.
if (hd1in gt 0.) then begin
	flared12 = alog10(hd1solar/hd1in)
	flared11 = flared11/flared12
endif else flared11 = 0.

if (swdisk4 eq 'yes') then begin
	if (zd4in gt zd4) then flared41 = alog10((zd4solar - zd4)/(zd4in - zd4)) else flared41 = 0.
	if (hd4in gt 0.) then begin
		flared42 = alog10(hd4solar/hd4in)
		flared41 = flared41/flared42
	endif else flared41 = 0.
endif

if (swdisk6 eq 'yes') then begin
	if (zd6in gt zd6) then flared61 = alog10((zd6solar - zd6)/(zd6in - zd6)) else flared61 = 0.
	if (hd6in gt 0.) then begin
		flared62 = alog10(hd6solar/hd6in)
		flared61 = flared61/flared62
	endif else flared61 = 0. 
endif

for k=0L,dim_positions-1 do begin  ;loop for different locations in the galaxy
	readf, unit1, ss
	readf, unit1, ss
	readf, unit1,  radius, height
	;read radial and vertical position in pc
	rr[k] = radius
	zz[k] = height
	readf, unit1, ss
	readf, unit1, ss
	readf, unit5, ss
	;read energy absorbed and emitted in W/H
	readf, unit5, pass1,pass2,pass3,pass4
	if pass1 ne radius then begin
		print, 'unexpected radial position in energy file; program ;stops'
		goto, mark1
	endif
	if pass2 ne height then begin
		print, 'unexpected vertical position in energy file; program ;stops'
		goto, mark1
	endif
  
	;Atentie aici! Tratam separat cazurile hdin=0 si hd1in=0 pentru care 
	;avem zdin=zd si zd1in=zd1
	if (zdin gt zd) then zd01[k] = zd + (zdin - zd) * (rr[k]/hdin)^flared00 else zd01[k] = zd
	if (zd1in gt zd1) then zd11[k] = zd1 + (zd1in - zd1) * (rr[k]/hd1in)^flared11 else zd11[k] = zd1
	if (swdisk3 eq 'yes') then begin
		if (zd3in gt zd3) then zd31[k] = zd3 + (zd3in - zd3) * (rr[k]/hd3in)^flared30 else zd31[k] = zd3
	endif
	if (swdisk4 eq 'yes') then begin
		if (zd4in gt zd4) then zd41[k] = zd4 + (zd4in - zd4) * (rr[k]/hd4in)^flared41 else zd41[k] = zd4
	endif
	if (swdisk5 eq 'yes') then begin
		if (zd5in gt zd5) then zd51[k] = zd5 + (zd5in - zd5) * (rr[k]/hd5in)^flared50 else zd51[k] = zd5
	endif
	if (swdisk6 eq 'yes') then begin
		if (zd6in gt zd6) then zd61[k] = zd6 + (zd6in - zd6) * (rr[k]/hd6in)^flared61 else zd61[k] = zd6
	endif
     
if idisk1 eq 1 then begin ;double exponential distribution first disk
      if rr[k] ge hdin then begin
           dust1[k]=nHc1Cext * (zd/zd01[k])*$
                    exp(-rr[k]/hd) * exp(-zz[k]/zd01[k])
      endif
 
      if rr[k] lt hdin then begin
          dust1[k]=nHc1Cext * (zd/zd01[k])*$
          ((rr[k]/hdin)*(1.-xid0)+xid0)*$
          exp(-hdin/hd) * exp(-zz[k]/zd01[k])
       endif
     endif ;end double exponential distribution first disk

if (swdisk3 eq 'yes') then begin
     if idisk3 eq 1 then begin ;double exponential distribution third disk
      if rr[k] ge hd3in then begin
           dust3[k]=nHc3Cext * (zd3/zd31[k])*$
                    exp(-rr[k]/hd3) * exp(-zz[k]/zd31[k])
      endif
 
      if rr[k] lt hd3in then begin
          dust3[k]=nHc3Cext * (zd3/zd31[k])*$
          ((rr[k]/hd3in)*(1.-xid3)+xid3)*$
          exp(-hd3in/hd3) * exp(-zz[k]/zd31[k])
       endif
     endif ;end double exponential distribution third disk
endif

if (swdisk5 eq 'yes') then begin
     if idisk5 eq 1 then begin ;double exponential distribution fifth disk
     if rr[k] ge hd5in then begin
           dust5[k]=nHc5Cext * (zd5/zd51[k])*$
                    exp(-rr[k]/hd5) * exp(-zz[k]/zd51[k])
      endif
 
      if rr[k] lt hd5in then begin
          dust5[k]=nHc5Cext * (zd5/zd51[k])*$
          ((rr[k]/hd5in)*(1.-xid5)+xid5)*$
          exp(-hd5in/hd5) * exp(-zz[k]/zd51[k])
       endif
     endif ;end double exponential distribution fifth disk
endif
      if idisk2 eq 1 then begin ;double exponential distribution second disk 
       if rr[k] ge hd1in then begin
          dust2[k]=nHc2Cext * (zd1/zd11[k]) * $
                    exp(-rr[k]/hd1) * exp(-zz[k]/zd11[k])
       endif
       if rr[k] lt hd1in then begin
          dust2[k]=nHc2Cext * (zd1/zd11[k])*$
          ((rr[k]/hd1in)*(1.-xid1)+xid1)*$
          exp(-hd1in/hd1) * exp(-zz[k]/zd11[k])
       endif
      endif ;end double exponential distribution second disk

if (swdisk4 eq 'yes') then begin
      if idisk4 eq 1 then begin ;double exponential distribution forth disk 
       if rr[k] ge hd4in then begin
          dust4[k]=nHc4Cext * (zd4/zd41[k]) * $
                    exp(-rr[k]/hd4) * exp(-zz[k]/zd41[k])
       endif
       if rr[k] lt hd4in then begin
          dust4[k]=nHc4Cext * (zd4/zd41[k])*$
          ((rr[k]/hd4in)*(1.-xid4)+xid4)*$
          exp(-hd4in/hd4) * exp(-zz[k]/zd41[k])
       endif
      endif ;end double exponential distribution forth disk
endif

if (swdisk6 eq 'yes') then begin
      if idisk6 eq 1 then begin ;double exponential distribution sixth disk 
       if rr[k] ge hd6in then begin
          dust6[k]=nHc6Cext * (zd6/zd61[k]) * $
                    exp(-rr[k]/hd6) * exp(-zz[k]/zd61[k])
       endif
       if rr[k] lt hd6in then begin
          dust6[k]=nHc6Cext * (zd6/zd61[k])*$
          ((rr[k]/hd6in)*(1.-xid6)+xid6)*$
          exp(-hd6in/hd6) * exp(-zz[k]/zd61[k])
       endif
      endif ;end double exponential distribution sixth disk
endif

      if idisk1 eq 3 then begin ;sech2 vertical distribution first disk
        sech01 = 1.d0/cosh(zz[k]/zd01[k])
        sech02 = sech01 * sech01
        if rr[k] ge hdin then begin
          dust1[k]=nHc1Cext *(zd/zd01[k])*$
                    exp(-rr[k]/hd) * sech02
        endif 
        if rr[k] lt hdin then begin
          dust1[k]=nHc1Cext * (zd/zd01[k])*$
          ((rr[k]/hdin)*(1.-xid0)+xid0)*$
          exp(-hdin/hd) * sech02
        endif
       endif ;end sech2 vertical  distribution first disk

if (swdisk3 eq 'yes') then begin
      if idisk3 eq 3 then begin ;sech2 vertical distribution third disk
        sech31 = 1.d0/cosh(zz[k]/zd31[k])
        sech32 = sech31 * sech31
        if rr[k] ge hd3in then begin
          dust3[k]=nHc3Cext *(zd3/zd31[k])*$
                    exp(-rr[k]/hd3) * sech32
        endif 
        if rr[k] lt hd3in then begin
          dust3[k]=nHc3Cext * (zd3/zd31[k])*$
          ((rr[k]/hd3in)*(1.-xid3)+xid3)*$
          exp(-hd3in/hd3) * sech32
        endif
       endif ;end sech2 vertical  distribution third disk
endif

if (swdisk5 eq 'yes') then begin
      if idisk5 eq 3 then begin ;sech2 vertical distribution fifth disk
        sech51 = 1.d0/cosh(zz[k]/zd51[k])
        sech52 = sech51 * sech51
        if rr[k] ge hd5in then begin
          dust5[k]=nHc5Cext *(zd5/zd51[k])*$
                    exp(-rr[k]/hd5) * sech52
        endif 
        if rr[k] lt hd5in then begin
          dust5[k]=nHc5Cext * (zd5/zd51[k])*$
          ((rr[k]/hd5in)*(1.-xid5)+xid5)*$
          exp(-hd5in/hd5) * sech52
        endif
       endif ;end sech2 vertical  distribution fifth disk
endif

       if idisk2 eq 3 then begin ;sech2 vertical  distribution second disk 
         sech11 = 1.d0/cosh(zz[k]/zd11[k])
         sech12=sech11*sech11
        if rr[k] ge hd1in then begin
          dust2[k]=nHc2Cext * (zd1/zd11[k]) * $
                    exp(-rr[k]/hd1) * sech12
        endif
        if rr[k] lt hd1in then begin
          dust2[k]=nHc2Cext * (zd1/zd11[k])*$
          ((rr[k]/hd1in)*(1.-xid1)+xid1)*$
          exp(-hd1in/hd1) * sech12
        endif
       endif ;end sech2 vertical distribution second disk

if (swdisk4 eq 'yes') then begin
       if idisk4 eq 3 then begin ;sech2 vertical  distribution forth disk 
         sech41 = 1.d0/cosh(zz[k]/zd41[k])
         sech42=sech41*sech41
        if rr[k] ge hd4in then begin
          dust4[k]=nHc4Cext * (zd4/zd41[k]) * $
                    exp(-rr[k]/hd4) * sech42
        endif
        if rr[k] lt hd4in then begin
          dust4[k]=nHc4Cext * (zd4/zd41[k])*$
          ((rr[k]/hd4in)*(1.-xid4)+xid4)*$
          exp(-hd4in/hd4) * sech42
        endif
       endif ;end sech2 vertical distribution forth disk
endif

if (swdisk6 eq 'yes') then begin
       if idisk6 eq 3 then begin ;sech2 vertical  distribution sixth disk 
         sech61 = 1.d0/cosh(zz[k]/zd61[k])
         sech62=sech61*sech61
        if rr[k] ge hd6in then begin
          dust6[k]=nHc6Cext * (zd6/zd61[k]) * $
                    exp(-rr[k]/hd6) * sech62
        endif
        if rr[k] lt hd6in then begin
          dust6[k]=nHc6Cext * (zd6/zd61[k])*$
          ((rr[k]/hd6in)*(1.-xid6)+xid6)*$
          exp(-hd6in/hd6) * sech62
        endif
       endif ;end sech2 vertical distribution sixth disk
endif
;truncating dust disks
      rplane=sqrt(rr[k]*rr[k]+zz[k]*zz[k])
      rplane=double(rplane)
;     first disk
      rtrund=double(rtrund)
      sigtruncate=sharpd*rtrund
      truncate=0.5*(1.0-erf((rplane-rtrund)/sigtruncate))
;     second disk
      rtrund1=double(rtrund1)
      sigtruncate1=sharpd1*rtrund1
      truncate1=0.5*(1.0-erf((rplane-rtrund1)/sigtruncate1))
;     third disk
if (swdisk3 eq 'yes') then begin
      rtrund3=double(rtrund3)
      sigtruncate3=sharpd3*rtrund3
      truncate3=0.5*(1.0-erf((rplane-rtrund3)/sigtruncate3))
endif
;     forth disk
if (swdisk4 eq 'yes') then begin
      rtrund4=double(rtrund4)
      sigtruncate4=sharpd4*rtrund4
      truncate4=0.5*(1.0-erf((rplane-rtrund4)/sigtruncate4))
endif
;     fifth disk
if (swdisk5 eq 'yes') then begin
      rtrund5=double(rtrund5)
      sigtruncate5=sharpd5*rtrund5
      truncate5=0.5*(1.0-erf((rplane-rtrund5)/sigtruncate5))
endif
;     sixth disk
if (swdisk6 eq 'yes') then begin
      rtrund6=double(rtrund6)
      sigtruncate6=sharpd6*rtrund6
      truncate6=0.5*(1.0-erf((rplane-rtrund6)/sigtruncate6))
endif

dust1[k]=dust1[k]*truncate
dust2[k]=dust2[k]*truncate1
if (swdisk3 eq 'yes') then dust3[k]=dust3[k]*truncate3
if (swdisk4 eq 'yes') then dust4[k]=dust4[k]*truncate4
if (swdisk5 eq 'yes') then dust5[k]=dust5[k]*truncate5
if (swdisk6 eq 'yes') then dust6[k]=dust6[k]*truncate6

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	if rr[k] lt hdtin then begin
		dust1[k]=0.
	endif
        if rr[k] lt hd1tin then begin
                dust2[k]=0.
        endif
  if (swdisk3 eq 'yes') then begin
        if rr[k] lt hd3tin then begin
                dust3[k]=0.
        endif
  endif
  if (swdisk4 eq 'yes') then begin      
	if rr[k] lt hd4tin then begin
                dust4[k]=0.
        endif
  endif
  if (swdisk5 eq 'yes') then begin
        if rr[k] lt hd5tin then begin
                dust5[k]=0.
        endif
  endif
  if (swdisk6 eq 'yes') then begin
        if rr[k] lt hd6tin then begin
                dust6[k]=0.
        endif
  endif

;;; JJT edit - 4 September 17 ;;;
;Set negative dust values to 0
;First dust disk
	if dust1[k] lt 0. then begin
	   dust1[k] = 0.
	endif

;Second dust disk
	if dust2[k] lt 0. then begin
	   dust2[k] = 0.
	endif

;JJT 16/10/17
;Third dust disk
  if (swdisk3 eq 'yes') then begin
	if dust3[k] lt 0. then begin
	   dust3[k] = 0.
	endif
  endif
;JJT 16/10/17
;Forth dust disk
  if (swdisk4 eq 'yes') then begin
	if dust4[k] lt 0. then begin
	   dust4[k] = 0.
	endif
  endif
;JJT 16/10/17
;Fifth dust disk
  if (swdisk5 eq 'yes') then begin
	if dust5[k] lt 0. then begin
	   dust5[k] = 0.
	endif
  endif
;JJT 16/10/17
;Sixth dust disk
  if (swdisk6 eq 'yes') then begin
	if dust6[k] lt 0. then begin
	   dust6[k] = 0.
	endif
  endif

dust[k] = dust1[k] + dust2[k]
if (swdisk3 eq 'yes') then dust[k] = dust[k] + dust3[k]
if (swdisk4 eq 'yes') then dust[k] = dust[k] + dust4[k]
if (swdisk5 eq 'yes') then dust[k] = dust[k] + dust5[k]
if (swdisk6 eq 'yes') then dust[k] = dust[k] + dust6[k]

     for i=0L,dim-1 do begin
         ;read wavelength in AA and luminosity in erg/s/sr/H/A
        readf, unit1, ll,lu
         ;if k eq 263 and i eq 628 then print, ll, lu
        lambda[i] = ll
         ;the luminosity for the first and second disk in units of 
         ;W/A/H/pc (not yet devided to the extinction cross section)
         ;1.d-7 comes from the transformation from erg/sec into W
        lum1_comp[k,i]=4.*!Dpi*1.d-7 * dust1[k]*lu
        lum2_comp[k,i]=4.*!Dpi*1.d-7 * dust2[k]*lu

       if (swdisk3 eq 'yes') then lum3_comp[k,i]=4.*!Dpi*1.d-7 * dust3[k]*lu
       if (swdisk4 eq 'yes') then lum4_comp[k,i]=4.*!Dpi*1.d-7 * dust4[k]*lu
       if (swdisk5 eq 'yes') then lum5_comp[k,i]=4.*!Dpi*1.d-7 * dust5[k]*lu
       if (swdisk6 eq 'yes') then lum6_comp[k,i]=4.*!Dpi*1.d-7 * dust6[k]*lu

;although dust*[k] will make lum*_comp[k,i]=0 if rr[k] lt hd*tin
;however this is not necessarily the same as hs*tin.
;set lum*_comp=0 if rr[k] lt hs*tin
	if rr[k] lt hstin then begin
		lum1_comp[k,i]=0
	endif
        if rr[k] lt hs1tin then begin
                lum2_comp[k,i]=0
        endif
   if (swdisk3 eq 'yes') then begin
        if rr[k] lt hs3tin then begin
                lum3_comp[k,i]=0
        endif
   endif
   if (swdisk4 eq 'yes') then begin
        if rr[k] lt hs4tin then begin
                lum4_comp[k,i]=0
        endif
   endif
   if (swdisk5 eq 'yes') then begin
        if rr[k] lt hs5tin then begin
                lum5_comp[k,i]=0
        endif
   endif
   if (swdisk6 eq 'yes') then begin
        if rr[k] lt hs6tin then begin
                lum6_comp[k,i]=0
        endif
   endif

	lum_comp[k,i] = lum1_comp[k,i] + lum2_comp[k,i]
	if (swdisk3 eq 'yes') then lum_comp[k,i] = lum_comp[k,i] + lum3_comp[k,i]
	if (swdisk4 eq 'yes') then lum_comp[k,i] = lum_comp[k,i] + lum4_comp[k,i]
	if (swdisk5 eq 'yes') then lum_comp[k,i] = lum_comp[k,i] + lum5_comp[k,i]
	if (swdisk6 eq 'yes') then lum_comp[k,i] = lum_comp[k,i] + lum6_comp[k,i]


     endfor ;end loop in wavelength

       energy_abs1[k] = pass3 * dust1[k]  ;in W/H/pc
       energy_abs2[k] = pass3 * dust2[k]  ;in W/H/pc

;JJT 16/10/17
       if (swdisk3 eq 'yes') then energy_abs3[k] = pass3 * dust3[k]  ;in W/H/pc
       if (swdisk4 eq 'yes') then energy_abs4[k] = pass3 * dust4[k]  ;in W/H/pc
       if (swdisk5 eq 'yes') then energy_abs5[k] = pass3 * dust5[k]  ;in W/H/pc
       if (swdisk6 eq 'yes') then energy_abs6[k] = pass3 * dust6[k]  ;in W/H/pc

;    energy_abs_comp[k] = energy_abs1[k] + energy_abs2[k]  ;in W/H/pc
energy_abs_comp[k] = energy_abs1[k] + energy_abs2[k] + $
	energy_abs3[k] + energy_abs4[k] + $
	energy_abs5[k] + energy_abs6[k]  ;in W/H/pc
endfor ;end loop in position

;the total luminosity (summation over composition) 
;in units of W/A/H/pc 
lum[*,*] = lum[*,*] + lum_comp[*,*]
lum1[*,*] = lum1[*,*] + lum1_comp[*,*]
lum2[*,*] = lum2[*,*] + lum2_comp[*,*]
if (swdisk3 eq 'yes') then lum3[*,*] = lum3[*,*] + lum3_comp[*,*]
if (swdisk4 eq 'yes') then lum4[*,*] = lum4[*,*] + lum4_comp[*,*]
if (swdisk5 eq 'yes') then lum5[*,*] = lum5[*,*] + lum5_comp[*,*]
if (swdisk6 eq 'yes') then lum6[*,*] = lum6[*,*] + lum6_comp[*,*]

energy_abs[*] = energy_abs[*] + energy_abs_comp[*] ;in W/H/pc
free_lun, unit5
free_lun, unit1

endfor ;end loop in composition

;transform from cm^2/H in pc^2/H
;print, 'Cext_B', Cext_B
Cext_B = Cext_B/(3.086 * 1.d+18)^2
;print,'Cext_B', Cext_B

;luminosity/volume in W/A/pc^3
lum = lum/Cext_B
lum1 = lum1/Cext_B
lum2 = lum2/Cext_b

if (swdisk3 eq 'yes') then lum3 = lum3/Cext_B
if (swdisk4 eq 'yes') then lum4 = lum4/Cext_b
if (swdisk5 eq 'yes') then lum5 = lum5/Cext_B
if (swdisk6 eq 'yes') then lum6 = lum6/Cext_b
energy_abs = energy_abs/Cext_B  ;in W/pc^3

;luminosity/volume in W/Hz/pc^3
for k = 0L,dim_positions-1 do begin
	lum[k,*] = lum[k,*] * 1.d-10 * lambda[*] * lambda[*]/cl
	lum1[k,*] = lum1[k,*] * 1.d-10 * lambda[*] * lambda[*]/cl
	lum2[k,*] = lum2[k,*] * 1.d-10 * lambda[*] * lambda[*]/cl
	if (swdisk3 eq 'yes') then lum3[k,*] = lum3[k,*] * 1.d-10 * lambda[*] * lambda[*]/cl
	if (swdisk4 eq 'yes') then lum4[k,*] = lum4[k,*] * 1.d-10 * lambda[*] * lambda[*]/cl
	if (swdisk5 eq 'yes') then lum5[k,*] = lum5[k,*] * 1.d-10 * lambda[*] * lambda[*]/cl
	if (swdisk6 eq 'yes') then lum6[k,*] = lum6[k,*] * 1.d-10 * lambda[*] * lambda[*]/cl
endfor
;transform the input two dimensional array of luminosity into a 3 dimensional
;array 
k=0L

for kr = 0L,dim_positionsr-1 do begin
	for kz = 0L,dim_positionsz-1 do begin    
		lumdouble[kr,kz,*] = lum[k,*]
		lumdouble1[kr,kz,*] = lum1[k,*]
		lumdouble2[kr,kz,*] = lum2[k,*]
		if (swdisk3 eq 'yes') then lumdouble3[kr,kz,*] = lum3[k,*]
		if (swdisk4 eq 'yes') then lumdouble4[kr,kz,*] = lum4[k,*]
		if (swdisk5 eq 'yes') then lumdouble5[kr,kz,*] = lum5[k,*]
		if (swdisk6 eq 'yes') then lumdouble6[kr,kz,*] = lum6[k,*]
		
		energy_abs_double[kr,kz] = energy_abs[k]
		rrr[kr,kz] = rr[k]
		zzz[kr,kz] = zz[k]
		k = k + 1
	endfor
endfor

;Remove negative values from luminosity
;NB the where command will cycle through all
;parts of the array and change values acordingly

;lumdouble
;locate where lumdouble < 0
	iq=where(lumdouble lt 0)
;set values at iq to 0
	if iq ne -1 then begin
	 lumdouble[iq]=0
	endif
;lumdouble1
;locate where lumdouble1 < 0
        iq=where(lumdouble1 lt 0)
;set values at iq to 0
	if iq ne -1 then begin
         lumdouble1[iq]=0
	endif
;lumdouble2
;locate where lumdouble2 < 0
        iq=where(lumdouble2 lt 0)
;set values at iq to 0
	if iq ne -1 then begin
         lumdouble2[iq]=0
	endif

;lumdouble3
;locate where lumdouble3 < 0
if (swdisk3 eq 'yes') then begin
        iq=where(lumdouble3 lt 0)
;set values at iq to 0
	if iq ne -1 then begin
         lumdouble3[iq]=0
	endif
endif
;lumdouble4
;locate where lumdouble4 < 0
if (swdisk4 eq 'yes') then begin
        iq=where(lumdouble4 lt 0)
;set values at iq to 0
	if iq ne -1 then begin
         lumdouble4[iq]=0
	endif
endif

;lumdouble5
;locate where lumdouble5 < 0
if (swdisk5 eq 'yes') then begin
        iq=where(lumdouble5 lt 0)
;set values at iq to 0
	if iq ne -1 then begin
         lumdouble5[iq]=0
	endif
endif
;lumdouble6
;locate where lumdouble6 < 0
if (swdisk6 eq 'yes') then begin
        iq=where(lumdouble6 lt 0)
;set values at iq to 0
	if iq ne -1 then begin
         lumdouble6[iq]=0
	endif
endif

;integrate luminosity over all input positions and calculate energy emitted

pass = dblarr(dim)
pass_energy = dblarr(dim)

for kz = 0L,dim_positionsz-2 do begin
	for kr = 0L,dim_positionsr-2 do begin
		pass[*] = !Dpi * (rrr[kr+1,kz]^2-rrr[kr,kz]^2) * (zzz[kr,kz+1]-zzz[kr,kz]) * $
			(lumdouble[kr,kz,*]+lumdouble[kr+1,kz+1,*])/2.
		pass_energy[*] = pass_energy[*] + pass[*]  ;in W/Hz
	endfor
endfor

pass_abs = 0.d
pass_energy_abs = 0.d
;integrate energy abs over all input positions and calculate energy absprbed
for kz = 0L,dim_positionsz-2 do begin
	for kr = 0L,dim_positionsr-2 do begin
		pass_abs = !Dpi * (rrr[kr+1,kz]^2-rrr[kr,kz]^2) * $
			(zzz[kr,kz+1]-zzz[kr,kz]) * $
			(energy_abs_double[kr,kz]+energy_abs_double[kr+1,kz+1])/2.
		pass_energy_abs = pass_energy_abs + pass_abs  ;in W
	endfor
endfor
print, 'pass_energy_abs [W]: ',total(pass_energy_abs *2d)

;the energy is derived only for half of the galaxy. One needs to multiply with
;2 to get the energy emitted by the whole galaxy
pass_energy = pass_energy * 2.

pass_energy_abs = pass_energy_abs * 2. * 1.d+7
print, 'energy absorbed [erg/s]', pass_energy_abs
print, 'energy_absorbed [W]', pass_energy_abs * 1d-7

namew ='total_luminosity_raw'+'_'+model+'_q'+qyear+$
               '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.dat'
 openw,unit3,dir4+Namew,/get_lun
 print, 'write ', dir4+Namew
 printf, unit3, 'central face-on tau B'
 printf, unit3, tau
 printf, unit3, 'sfr'
 printf, unit3, sfr
 printf, unit3, 'old'
 printf, unit3, old
 printf, unit3, 'bd'
 printf, unit3, bd
 printf,unit3
 printf,unit3,'Number of wavelengths:'
 printf,unit3,dim
 printf,unit3
 printf,unit3,'       lambda(A)   L(W Hz^-1)'
 for i = 0, dim-1 do begin
    printf, unit3, lambda[i], pass_energy[i]
 endfor
;transform from W/Hz in W/A
wave_test = lambda
pc_m = 3.0857D16        ; conversion between [pc] and [m]
dist_pc = 8.58D6        ; distance to the galaxy in [pc]
dist_m = dist_pc * pc_m ; distance to the galaxy in [m]
flux_to_lum = (4*!DPI*dist_m^2) *1E-26  ; factor converting from flux to luminosity [Jy]->[W/Hz]
pass_jy = pass_energy / flux_to_lum
save, wave_test, pass_jy, filename="irr.save"
pass_energy = pass_energy * cl * 1.e+10/lambda^2

;integrate over wavelength
dlambda = lambda[1:(dim-1)] - lambda[0:(dim-2)]
pass_energym = 0.5 * (pass_energy[1:(dim-1)] + pass_energy[0:(dim-2)])
energy = total(pass_energym * dlambda) ;in W

energy = energy * 1.d+7  ;in erg/sec
print, 'energy raw [erg/sec]', energy

Namedust='dust'+'_'+model+'_q'+qyear+'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.xdr'
;JJT 16/10/17
if (swdisk3 eq 'no' and swdisk4 eq 'no' and swdisk5 eq 'no' and swdisk6 eq 'no') then $ 
    save, dust1, dust2, rr, zz ,file=dir4+Namedust

if (swdisk3 eq 'yes' and swdisk4 eq 'yes' and swdisk5 eq 'yes' and swdisk6 eq 'yes') then $
    save, dust1, dust2, dust3, dust4, dust5, dust6, rr, zz ,file=dir4+Namedust

Nameidl='grid_irr'+'_'+model+'_q'+qyear+'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.xdr'
save,lumdouble, lambda, rrr, zzz, dim,file=dir4+Nameidl

Nameidl1='grid_irr1'+'_'+model+'_q'+qyear+'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.xdr'
save,lumdouble1, lambda, rrr, zzz, dim,file=dir4+Nameidl1


Nameidl2='grid_irr2'+'_'+model+'_q'+qyear+'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.xdr'
save,lumdouble2, lambda, rrr, zzz, dim,file=dir4+Nameidl2

;JJT 16/10/17
if (swdisk3 eq 'yes') then begin
  Nameidl3='grid_irr3'+'_'+model+'_q'+qyear+'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.xdr'
  save,lumdouble3, lambda, rrr, zzz, dim,file=dir4+Nameidl3
endif

;JJT 16/10/17
if (swdisk4 eq 'yes') then begin
  Nameidl4='grid_irr4'+'_'+model+'_q'+qyear+'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.xdr'
  save,lumdouble4, lambda, rrr, zzz, dim,file=dir4+Nameidl4
endif

;JJT 16/10/17
if (swdisk5 eq 'yes') then begin
  Nameidl5='grid_irr5'+'_'+model+'_q'+qyear+'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.xdr'
  save,lumdouble5, lambda, rrr, zzz, dim,file=dir4+Nameidl5
endif

;JJT 16/10/17
if (swdisk6 eq 'yes') then begin
  Nameidl6='grid_irr6'+'_'+model+'_q'+qyear+'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.xdr'
  save,lumdouble6, lambda, rrr, zzz, dim,file=dir4+Nameidl6
endif

regrid,lumdouble,dr_req,dz_req,dr_out,dz_out,lumdouble_reg,rrr_reg,zzz_reg,nr,nz
Nameidl='grid'+'_'+model+'_q'+qyear+'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.xdr'
;save,lumdouble_reg, lambda, rrr_reg, zzz_reg, dim,file=dir4+Nameidl

regrid,lumdouble1,dr_req,dz_req,dr_out,dz_out,lumdouble1_reg,rrr_reg,zzz_reg,nr,nz
Nameidl1='grid'+'_'+model+'_q'+qyear+'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.xdr'
;save,lumdouble1_reg, lambda, rrr_reg, zzz_reg, dim,file=dir4+Nameidl1

regrid,lumdouble2,dr_req,dz_req,dr_out,dz_out,lumdouble2_reg,rrr_reg,zzz_reg,nr,nz
Nameidl2='grid'+'_'+model+'_q'+qyear+'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.xdr'
;save,lumdouble2_reg, lambda, rrr_reg, zzz_reg, dim,file=dir4+Nameidl2

regrid,lumdouble3,dr_req,dz_req,dr_out,dz_out,lumdouble3_reg,rrr_reg,zzz_reg,nr,nz
Nameidl3='grid'+'_'+model+'_q'+qyear+'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.xdr'
;save,lumdouble3_reg, lambda, rrr_reg, zzz_reg, dim,file=dir4+Nameidl3

regrid,lumdouble4,dr_req,dz_req,dr_out,dz_out,lumdouble4_reg,rrr_reg,zzz_reg,nr,nz
Nameidl4='grid'+'_'+model+'_q'+qyear+'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.xdr'
save,lumdouble4_reg, lambda, rrr_reg, zzz_reg, dim,file=dir4+Nameidl4

regrid,lumdouble5,dr_req,dz_req,dr_out,dz_out,lumdouble5_reg,rrr_reg,zzz_reg,nr,nz
Nameidl5='grid'+'_'+model+'_q'+qyear+'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.xdr'
;save,lumdouble5_reg, lambda, rrr_reg, zzz_reg, dim,file=dir4+Nameidl5

regrid,lumdouble6,dr_req,dz_req,dr_out,dz_out,lumdouble6_reg,rrr_reg,zzz_reg,nr,nz
Nameidl6='grid'+'_'+model+'_q'+qyear+'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.xdr'
;save,lumdouble6_reg, lambda, rrr_reg, zzz_reg, dim,file=dir4+Nameidl6

save,lumdouble_reg,lumdouble1_reg,lumdouble2_reg,lumdouble3_reg,lumdouble4_reg,lumdouble5_reg,lumdouble6_reg,lambda,rrr_reg,zzz_reg,file=dir4+Nameidl


;print,'dimensions of resampled arrays', nr,nz

;integrate the resampled luminosity over all resampled positions 
;and recalculate energy emitted

pass = dblarr(dim)
pass_energy = dblarr(dim)

pass1 = dblarr(dim)
pass_energy1 = dblarr(dim)

pass2 = dblarr(dim)
pass_energy2 = dblarr(dim)

pass3 = dblarr(dim)
pass_energy3 = dblarr(dim)

pass4 = dblarr(dim)
pass_energy4 = dblarr(dim)

pass5 = dblarr(dim)
pass_energy5 = dblarr(dim)

pass6 = dblarr(dim)
pass_energy6 = dblarr(dim)

for kz = 0L,nz-2 do begin
	for kr = 0L,nr-2 do begin
		pass[*] = !pi * (rrr_reg[kr+1,kz]^2-rrr_reg[kr,kz]^2) * $
		(zzz_reg[kr,kz+1]-zzz_reg[kr,kz]) * $
		(lumdouble_reg[kr,kz,*]+lumdouble_reg[kr+1,kz+1,*])/2.
		pass1[*] = !pi * (rrr_reg[kr+1,kz]^2-rrr_reg[kr,kz]^2) * $
		(zzz_reg[kr,kz+1]-zzz_reg[kr,kz]) * $
		(lumdouble1_reg[kr,kz,*]+lumdouble1_reg[kr+1,kz+1,*])/2.
		pass2[*] = !pi * (rrr_reg[kr+1,kz]^2-rrr_reg[kr,kz]^2) * $
		(zzz_reg[kr,kz+1]-zzz_reg[kr,kz]) * $
		(lumdouble2_reg[kr,kz,*]+lumdouble2_reg[kr+1,kz+1,*])/2.
		pass3[*] = !pi * (rrr_reg[kr+1,kz]^2-rrr_reg[kr,kz]^2) * $
		(zzz_reg[kr,kz+1]-zzz_reg[kr,kz]) * $
		(lumdouble3_reg[kr,kz,*]+lumdouble3_reg[kr+1,kz+1,*])/2.
		pass4[*] = !pi * (rrr_reg[kr+1,kz]^2-rrr_reg[kr,kz]^2) * $
		(zzz_reg[kr,kz+1]-zzz_reg[kr,kz]) * $
		(lumdouble4_reg[kr,kz,*]+lumdouble4_reg[kr+1,kz+1,*])/2.
		pass5[*] = !pi * (rrr_reg[kr+1,kz]^2-rrr_reg[kr,kz]^2) * $
		(zzz_reg[kr,kz+1]-zzz_reg[kr,kz]) * $
		(lumdouble5_reg[kr,kz,*]+lumdouble5_reg[kr+1,kz+1,*])/2.
		pass6[*] = !pi * (rrr_reg[kr+1,kz]^2-rrr_reg[kr,kz]^2) * $
		(zzz_reg[kr,kz+1]-zzz_reg[kr,kz]) * $
		(lumdouble6_reg[kr,kz,*]+lumdouble6_reg[kr+1,kz+1,*])/2.
		
		;IF rrr_reg[kr+1,kz] GT 7000. THEN BEGIN	; ZEROES LUMINOSITY AFTER 7000pc!!!
		;	pass(*) = 0.
		;	pass1(*) = 0.
		;	pass2(*) = 0.
		;	pass3(*) = 0.
		;	pass4(*) = 0.
		;	pass5(*) = 0.
		;	pass6(*) = 0.
		;ENDIF
		
		pass_energy[*] = pass_energy[*] + pass[*]  ;in W/Hz
		pass_energy1[*] = pass_energy1[*] + pass1[*]  ;in W/Hz
		pass_energy2[*] = pass_energy2[*] + pass2[*]  ;in W/Hz
		pass_energy3[*] = pass_energy3[*] + pass3[*]  ;in W/Hz
		pass_energy4[*] = pass_energy4[*] + pass4[*]  ;in W/Hz
		pass_energy5[*] = pass_energy5[*] + pass5[*]  ;in W/Hz
		pass_energy6[*] = pass_energy6[*] + pass6[*]  ;in W/Hz
	endfor
endfor

;the energy is derived only for half of the galaxy. One needs to multiply with
;2 to get the energy emitted by the whole galaxy
pass_energy = pass_energy * 2.
pass_energy1 = pass_energy1 * 2.
pass_energy2 = pass_energy2 * 2.
pass_energy3 = pass_energy3 * 2.
pass_energy4 = pass_energy4 * 2.
pass_energy5 = pass_energy5 * 2.
pass_energy6 = pass_energy6 * 2.

save,pass_energy1,pass_energy2,pass_energy3,pass_energy4,pass_energy5,pass_energy6,pass_energy,file=dir4+'reg_pass_energy'+ratoldident+'.xdr'

namew ='total_luminosity'+'_'+model+'_q'+qyear+$
               '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.dat'
 openw,unit4,dir4+Namew,/get_lun
 print, 'write ', dir4+Namew
 printf, unit4, 'central face-on tau B'
 printf, unit4, tau
 printf, unit4, 'sfr'
 printf, unit4, sfr
 printf, unit4, 'old'
 printf, unit4, old
 printf, unit4, 'bd'
 printf, unit4, bd
 printf, unit4, 'tau1'
 printf, unit4, tau1
 printf, unit4, 'hd [pc]'
 printf, unit4, hd
 printf, unit4, 'zd [pc]'
 printf, unit4, zd
 printf, unit4, 'hdin [pc]'
 printf, unit4, hdin
 printf, unit4, 'zdin [pc]'
 printf, unit4, zdin
 printf, unit4, 'hdsolar [pc]'
 printf, unit4, hdsolar
 printf, unit4, 'zdsolar [pc]'
 printf, unit4, zdsolar
 printf, unit4, 'tau2'
 printf, unit4, tau2
 printf, unit4, 'hd1 [pc]'
 printf, unit4, hd1
 printf, unit4, 'zd1 [pc]'
 printf, unit4, zd1
 printf, unit4, 'hd1in [pc]'
 printf, unit4, hd1in
 printf, unit4, 'zd1in [pc]'
 printf, unit4, zd1in
 printf, unit4, 'hd1solar [pc]'
 printf, unit4, hd1solar
 printf, unit4, 'zd1solar [pc]'
 printf, unit4, zd1solar
 printf, unit4, 'hs [pc]'
 printf, unit4, hs
 printf, unit4, 'zs [pc]'
 printf, unit4, zs
 printf, unit4, 'hsin [pc]'
 printf, unit4, hsin
 printf, unit4, 'zsin [pc]'
 printf, unit4, zsin
 printf, unit4, 'hssolar [pc]'
 printf, unit4, hssolar
 printf, unit4, 'zssolar [pc]'
 printf, unit4, zssolar
 printf, unit4, 'hs1 [pc]'
 printf, unit4, hs1
 printf, unit4, 'zs1 [pc]'
 printf, unit4, zs1
 printf, unit4, 'hs1in [pc]'
 printf, unit4, hs1in
 printf, unit4, 'zs1in [pc]'
 printf, unit4, zs1in
 printf, unit4, 'hs1solar [pc]'
 printf, unit4, hs1solar
 printf, unit4, 'zs1solar [pc]'
 printf, unit4, zs1solar
 printf, unit4, 'rtrun [pc]'
 printf, unit4, rtrun
 printf, unit4, 'rtrun1 [pc]'
 printf, unit4, rtrun1
 printf, unit4, 'xis0'
 printf, unit4, xis0
 printf, unit4, 'xis1'
 printf, unit4, xis1
 printf, unit4, 'xid0'
 printf, unit4, xid0
 printf, unit4, 'xid1'
 printf, unit4, xid1 
 printf, unit4, 'idisk1'
 printf, unit4, idisk1
 printf, unit4, 'idisk2'
 printf, unit4, idisk2
 printf,unit4
 printf,unit4,'Number of wavelengths:'
 printf,unit4,dim
 printf,unit4
 printf,unit4,'       lambda(A)   L(W Hz^-1)'
 for i = 0, dim-1 do begin
     printf, unit4, lambda[i], pass_energy[i]
 endfor

nu = cl / (lambda*1d-10)
dnu = nu[0:dim-2] - nu[1:dim-1]
pass_e = 0.5d * (pass_energy[1:(dim-1)] + pass_energy[0:(dim-2)])
enery = total(pass_energy * dnu) ;in W
print, 'energy [W] wavelength to freq: ',enery

;transform from W/Hz in W/A
pass_energy = pass_energy * cl * 1.e+10/lambda^2

;integrate over wavelength
dlambda = lambda[1:(dim-1)] - lambda[0:(dim-2)]
pass_energym = 0.5 * (pass_energy[1:(dim-1)] + pass_energy[0:(dim-2)])
energy = total(pass_energym * dlambda) ;in W

print, 'energy [W]:',energy
energy = energy * 1.d+7  ;in erg/sec
print, 'energy [erg/sec]', energy

free_lun, unit3
free_lun, unit4

;namew ='total_luminosity_irr1_'+model+'_q'+qyear+$
;               '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
;             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
;             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
;             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.dat'
; openw,unit5,dir4+Namew,/get_lun
; printf,unit5,'lambda(A)	L(W Hz^-1)'
; for i = 0, dim-1 do begin
;     printf, unit5, lambda[i], pass_energy1[i]
; endfor
;
;namew ='total_luminosity_irr2_'+model+'_q'+qyear+$
;               '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
;             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
;             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
;             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.dat'
; openw,unit6,dir4+Namew,/get_lun
; printf,unit6,'lambda(A)	L(W Hz^-1)'
; for i = 0, dim-1 do begin
;     printf, unit6, lambda[i], pass_energy2[i]
; endfor
;
;namew ='total_luminosity_irr3_'+model+'_q'+qyear+$
;               '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
;             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
;             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
;             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.dat'
; openw,unit7,dir4+Namew,/get_lun
; printf,unit7,'lambda(A)	L(W Hz^-1)'
; for i = 0, dim-1 do begin
;     printf, unit7, lambda[i], pass_energy3[i]
; endfor
;
;namew ='total_luminosity_irr4_'+model+'_q'+qyear+$
;               '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
;             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
;             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
;             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.dat'
; openw,unit8,dir4+Namew,/get_lun
; printf,unit8,'lambda(A)	L(W Hz^-1)'
; for i = 0, dim-1 do begin
;     printf, unit8, lambda[i], pass_energy4[i]
; endfor
;
;namew ='total_luminosity_irr5_'+model+'_q'+qyear+$
;               '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
;             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
;             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
;             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.dat'
; openw,unit9,dir4+Namew,/get_lun
; printf,unit9,'lambda(A)	L(W Hz^-1)'
; for i = 0, dim-1 do begin
;     printf, unit9, lambda[i], pass_energy5[i]
; endfor
;
;namew ='total_luminosity_irr6_'+model+'_q'+qyear+$
;               '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
;             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
;             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
;             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.dat'
; openw,unit10,dir4+Namew,/get_lun
; printf,unit10,'lambda(A)	L(W Hz^-1)'
; for i = 0, dim-1 do begin
;     printf, unit10, lambda[i], pass_energy6[i]
; endfor

mark1:
print, 'DONE: integrated_luminosity_template.pro ('+strtrim(ceil(systime(/seconds)-start_time),1)+' s)'
end
