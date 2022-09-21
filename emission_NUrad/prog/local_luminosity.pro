;this program produces a HII template scaled to the sfr and ffactor required
pro local_luminosity,model,qyear,tau,sfr,sfr4,sfr6,sfr7,old,bd,scaabs

start_time = systime(/seconds)
common scaling
close,/all
c = 2.99792458d+8 ;speed of light in m/s
L_unit_young = 4.235e+36 ;in W

suv=strcompress(string(round(sfr*100)),/remove_all)
sbd=strcompress(string(round(bd*100)),/remove_all)
sold=strcompress(string(round(old*100)),/remove_all)
stau=strcompress(string(round(tau*10)),/remove_all)
sffactor = strcompress(string(round(ffactor*100)),/remove_all)

dir=emission_dir
dir2=urad_dir+'indata/'

;read file with geometry
ss = ''

idisk1 = 0L
idisk2 = 0L

filename = 'geometry.in'
name = dir2+filename
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

;read the file with diffuse luminosity (W Hz^-1)
namer ='total_luminosity'+'_'+model+'_q'+qyear+$
            '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+'.dat'
openr, unit1, dir+'outdata_intlum/'+namer,/get_lun
;print, 'read ', dir+'outdata_intlum/'+namer
ss = ''
readf, unit1, ss
readf, unit1, newtau
readf, unit1, ss
readf, unit1, newsfr
readf, unit1, ss
readf, unit1, newold
readf, unit1, ss
readf, unit1, newbd
;check that the input parameters read from the file are the same as those
if newtau ne tau then begin
	print, 'unexpected tau in total lum file; program stops'
	goto, mark1
endif
if newsfr ne sfr then begin
	print, 'unexpected sfr in total lum file; program stops'
	goto, mark1
endif
if newbd ne bd then begin
	print, 'unexpected bd in total lum file; program stops'
	goto, mark1
endif
if newold ne old then begin
	print, 'unexpected old in total lum file; program stops'
	goto, mark1
endif
readf, unit1, ss
readf, unit1, tau1
readf, unit1, ss
readf, unit1, hd
readf, unit1, ss
readf, unit1, zd
readf, unit1, ss
readf, unit1, hdin
readf, unit1, ss
readf, unit1, zdin
readf, unit1, ss
readf, unit1, hdsolar
readf, unit1, ss
readf, unit1, zdsolar
readf, unit1, ss
readf, unit1, tau2
readf, unit1, ss
readf, unit1, hd1
readf, unit1, ss
readf, unit1, zd1
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
readf, unit1, ss
readf, unit1, dim
readf, unit1, ss
readf, unit1, ss
lambda = dblarr(dim)
lum_integ = dblarr(dim)
for i = 0, dim-1 do begin
	readf, unit1, ll, lu
	lambda[i] = ll
	lum_integ[i] = lu
endfor

lambdaout = lambda/1.d+4
dopita1, lambdaout, sedout, 1.0e7, 1.0, 0.1, 10000., 1

lum_template=sedout

;transform luminosity from W/Hz in W/A
lum_template_lambda = lum_template * c * 1.e+10/lambda^2

dlambda = lambda[1:(dim-1)] - lambda[0:(dim-2)]	; interval between lambda(n) and lambda(n+1)

;calculate the luminosity of the localised component
lum_templatem =  0.5 * (lum_template_lambda[1:(dim-1)] + lum_template_lambda[0:(dim-2)])	; average lum in intervals
int_lum_template = total(lum_templatem * dlambda)	; integrate luminosity over all lambda
print,'int_lum_template [W]', int_lum_template,'  (integrated luminosity of the localised component).'

;scaled to unit luminosity
lum_template_unit = lum_template * L_unit_young/int_lum_template
lum_template_unit_lambda = lum_template_lambda * L_unit_young/int_lum_template

;check scaling
lum_templatem =  0.5 * (lum_template_unit_lambda[1:(dim-1)] + lum_template_unit_lambda[0:(dim-2)])
int_lum_template = total(lum_templatem * dlambda)
print,'int_lum_template [W]', int_lum_template

;scaled to the sfr
sfr_total = sfr/(1-ffactor)
sfr_local = ffactor * sfr_total
sfr_local = ffactor	; Added CJI 07/05/21

sfr_total4 = sfr4/(1-ffactor4)
sfr_local4 = ffactor4 * sfr_total4
sfr_local4 = ffactor4	; Added CJI 07/05/21

sfr_total6 = sfr6/(1-ffactor6)
sfr_local6 = ffactor6 * sfr_total6
;sfr_local6 = ffactor6	; Added CJI 11/05/21

sfr_total7 = sfr7/(1-ffactor7)
sfr_local7 = ffactor7 * sfr_total7

print,'sfr_local', sfr_local
print,'sfr_local4', sfr_local4
print,'sfr_local6', sfr_local6
print,'sfr_local7', sfr_local7

;luminosity of the HII regions in W/Hz
lum_template_scaled = lum_template_unit * sfr_local 
lum_template_scaled4 = lum_template_unit * sfr_local4 
lum_template_scaled6 = lum_template_unit * sfr_local6 
lum_template_scaled7 = lum_template_unit * sfr_local7

name_local='localised_luminosity_'+model+'_q'+qyear+$
       '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
             '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
             '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
             '_reff'+sreff+'_ell'+sellipt+'_f'+sffactor+'.xdr'
print,'save ', dir+'outdata_intlum/'+name_local
save,lambda,lum_template_scaled,lum_template_scaled4,lum_template_scaled6,lum_template_scaled7,file=dir+'outdata_intlum/'+name_local

mark1:
free_lun, unit1
print, 'DONE: local_luminosity.pro ('+strtrim(ceil(systime(/seconds)-start_time),1)+' s)'
end
