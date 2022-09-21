;this program plots the integrated FIR SED of a galaxy 

pro plot_total_luminosity_template,model,qyear,tau,sfr,old,bd,scaabs,swheatr
;plot_total_luminosity_template, 'wd01','06',25.4,7.9,2.6,0.01,'abs','no'

start_time = systime(/seconds)
common scaling
close,/all

;define where the input file geometry.in is
geometry_dir=urad_dir+'indata/geometry.in'

;set up for new identifier for grids
ratoldident=''
if (swheatr eq 'yes') then ratoldident='_old'

dist = 8.58 ;distance in Mpc
c = 2.99792458d+8 ;speed of light in m/s 
suv=strcompress(string(round(sfr*100)),/remove_all)
sbd=strcompress(string(round(bd*100)),/remove_all)
sold=strcompress(string(round(old*100)),/remove_all)
stau=strcompress(string(round(tau*10)),/remove_all)
sffactor = strcompress(string(round(ffactor*100)),/remove_all)

ss=''
openr, gunit,geometry_dir,/get_lun	;open geometry_ext.in
readf, gunit, ss	;read geometry parameters
readf, gunit, tau1
readf, gunit, ss
readf, gunit, tau2
readf, gunit, ss
readf, gunit, hd
readf, gunit, ss
readf, gunit, zd
readf, gunit, ss
readf, gunit, hdin
readf, gunit, ss
readf, gunit, zdin
readf, gunit, ss
readf, gunit, hdsolar
readf, gunit, ss
readf, gunit, zdsolar
readf, gunit, ss
readf, gunit, hd1
readf, gunit, ss
readf, gunit, zd1
readf, gunit, ss
readf, gunit, hd1in
readf, gunit, ss
readf, gunit, zd1in
readf, gunit, ss
readf, gunit, hd1solar
readf, gunit, ss
readf, gunit, zd1solar
readf, gunit, ss
readf, gunit, h_bdisk
readf, gunit, ss
readf, gunit, h_vdisk
readf, gunit, ss
readf, gunit, h_idisk
readf, gunit, ss
readf, gunit, h_jdisk
readf, gunit, ss
readf, gunit, h_kdisk
readf, gunit, ss
readf, gunit, h_ir36disk
readf, gunit, ss
readf, gunit, h_ir45disk
readf, gunit, ss
readf, gunit, h_ir58disk
readf, gunit, ss
readf, gunit, zs
readf, gunit, ss
readf, gunit, hsin
readf, gunit, ss
readf, gunit, zsin
readf, gunit, ss
readf, gunit, hssolar
readf, gunit, ss
readf, gunit, zssolar
readf, gunit, ss
readf, gunit, hs1
readf, gunit, ss
readf, gunit, zs1
readf, gunit, ss
readf, gunit, hs1in
readf, gunit, ss
readf, gunit, zs1in
readf, gunit, ss
readf, gunit, hs1solar
readf, gunit, ss
readf, gunit, zs1solar
readf, gunit, ss
readf, gunit, rtruncate
readf, gunit, ss
readf, gunit, sharp
readf, gunit, ss
readf, gunit, rtruncated
readf, gunit, ss
readf, gunit, sharpd
readf, gunit, ss
readf, gunit, rtruncate1
readf, gunit, ss
readf, gunit, sharp
readf, gunit, ss
readf, gunit, rtruncated1
readf, gunit, ss
readf, gunit, sharpd1
readf, gunit, ss
readf, gunit, reff
readf, gunit, ss
readf, gunit, ellipt
readf, gunit, ss
readf, gunit, nsersic
readf, gunit, ss
readf, gunit, xis0
readf, gunit, ss
readf, gunit, xis1
readf, gunit, ss
readf, gunit, xid0
readf, gunit, ss
readf, gunit, xid1
readf, gunit, ss
readf, gunit, idisk1
readf, gunit, ss
readf, gunit, idisk2
free_lun, gunit

hs =strcompress(string(round(h_bdisk*1000.)),/remove_all)	    ;convert geometry parameters to strings BTC
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

dir = emission_dir
;read the file with integrated luminosity (W Hz^-1)
namer = 'total_luminosity'+'_'+model+'_q'+qyear+$
       '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
	'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+$
	'_zd1_'+szd1+'_hs'+hs+'_zs'+$
	szs+'_hs1_'+shs1+'_zs1_'+szs1+$
	'_reff'+sreff+'_ell'+sellipt+'_'+scaabs+ratoldident+'.dat'

openr, unit1, dir+'outdata_intlum/'+namer,/get_lun
print, 'read ', dir+'outdata_intlum/'+namer
ss=''
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

for i = 1, 70 do readf, unit1, ss	;skip to dim in input file BTC
readf, unit1, dim
readf, unit1, ss
readf, unit1, ss
lambda = dblarr(dim)
lum_integ = dblarr(dim)

restore, dir+'outdata_intlum/reg_pass_energy'+ratoldident+'.xdr'

for i = 0, dim-1 do begin
	readf, unit1, ll, lu
	lambda[i] = ll
	lum_integ[i] = lu
endfor
dist = dist * 3.086d+22                   ;in meter
flux = lum_integ/(4. * !pi * dist * dist) ;in W Hz^-1 m^-2
flux1 = pass_energy1/(4. * !pi * dist * dist) ;in W Hz^-1 m^-2
flux2 = pass_energy2/(4. * !pi * dist * dist) ;in W Hz^-1 m^-2
flux3 = pass_energy3/(4. * !pi * dist * dist) ;in W Hz^-1 m^-2
flux4 = pass_energy4/(4. * !pi * dist * dist) ;in W Hz^-1 m^-2
flux5 = pass_energy5/(4. * !pi * dist * dist) ;in W Hz^-1 m^-2
flux6 = pass_energy6/(4. * !pi * dist * dist) ;in W Hz^-1 m^-2
flux = flux * 1.d+26                      ;in Jy
flux1 = flux1 * 1.d+26                      ;in Jy
flux2 = flux2 * 1.d+26                      ;in Jy
flux3 = flux3 * 1.d+26                      ;in Jy
flux4 = flux4 * 1.d+26                      ;in Jy
flux5 = flux5 * 1.d+26                      ;in Jy
flux6 = flux6 * 1.d+26                      ;in Jy

xx1 = dblarr(1000)
funct1=dblarr(1000)

lum_diff = lum_integ 

name_diff='diffuse_luminosity_'+model+'_q'+qyear+$
       '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+'_f'+sffactor+ratoldident+'.xdr'
print, 'save ', dir+'outdata_intlum/'+name_diff
save, lambda,lum_diff,flux,flux1,flux2,flux3,flux4,flux5,flux6, file=dir+'outdata_intlum/'+name_diff

name = dir+'figures/sed'+'_'+model+'_q'+qyear+$
       '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+'_f'+sffactor+ratoldident+'.ps'

set_plot, 'ps'
print, name
device, filename=name
maxlum = max(lum_diff)
minlum = min(lum_diff)
!p.charsize=1.5
!p.charsize=1.5
plot_oo,lambda/1.d+4,lum_diff,xrange=[2.d,3.d+3],$
        yrange=[minlum,maxlum + 0.2 * maxlum],$
        xstyle=1,ystyle=1,xtitle='!4k!3 [!4l!3m]',ytitle='L!B!4m!3!N [W/Hz]',$
        linestyle=0
;oplot,lambda/1.d+4,total_flux_old
;oplot,lambda/1.d+4,lum_diff, linestyle=1
;oplot,lambda/1.d+4,0.36 * int_flux, linestyle=2
;oploterror,lambda_meas/1.d+4,flux_meas,e_lambda,e_flux,psym=4
;oploterror,lambda_meas/1.d+4,flux_meas1,e_lambda,e_flux,psym=2
;oplot,lambda/1.d+4,lum_template, linestyle=2
device,/close


;transform luminosity from W/Hz in W/A
lum_diff = lum_diff * c * 1.e+10/lambda^2
;lum_template = lum_template * c * 1.e+10/lambda^2
;total_lum = total_lum * c * 1.e+10/lambda^2
         

;plot results in units of energy
set_plot,'ps'
name=dir+'figures/sed1'+'_'+model+'_q'+qyear+$
       '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+'_f'+sffactor+ratoldident+'.ps'
print, name
device,filename=name
;emissivity = [11.,15.,10.,7.,5.,15.,30.,20.] * 1.d+33
;lambdaplot = [5.8,8.,8.,12.,25.,70.,100.,160.]

plot_oo,lambda/1.d+4,lum_diff*lambda,$;*0.77D+43,$
        xrange=[2.d,3.d+3],yrange=[3.d+33,3.d+36],$
        ;xrange=[0.1,1.d+3],yrange=[1.d+33,1.d+36],$
        xstyle=1,ystyle=1,$;
        xtitle='!4k!3 [!4l!3m]',ytitle='!4k!3 L!B!4k!3!N [W]',$
        linestyle=0
;        title='10*Tau= '+stau+'100*SFR= '+suv
;oplot,lambdaplot,emissivity,psym=4
;oplot,lambda/1.d+4,flux*lambda, linestyle=1
;oploterror,lambda_meas/1.d+4,flux_meas*lambda_meas,e_lambda,e_flux,psym=4
;oplot,lambda/1.d+4,ffactor * sedout * lambda, linestyle=2

mark1:
device,/close


;calculate the luminosity of the diffuse component
dlambda = lambda[1:(dim-1)] - lambda[0:(dim-2)]

lum_diffm = 0.5 * (lum_diff[1:(dim-1)] + lum_diff[0:(dim-2)])
int_lum_diff = total(lum_diffm * dlambda)  ;in W
print,'int_lum_diff [W]', int_lum_diff


;1 Msolar per year produces 0.267 * 1.d+43 erg/sec ionising radiation
;shortwards of 912 A, and 2.495 * 1.d+43 erg/sec of non-ionising radiation from
;912 A to 4430 A. Only 0.3 of the ionising radiation goes into heating the
;dust, the rest goes into the recombination lines
;sfr_loc = lum_loc/(2.495 * 1.d+43 * 1.d-7 + 0.3 * 0.267 * 1.d+43 * 1.d-7)

free_lun, unit1

urad_save_dir = urad_dir+'saves/model/'
i=0
wavelengths=['uv15']
restore, urad_save_dir+model+'_datafile_'+wavelengths[i]+'_'+scaabs+'.save'


print, 'DONE: plot_total_luminosity_template.pro ('+strtrim(ceil(systime(/seconds)-start_time),1)+' s)'
end
