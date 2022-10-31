;this routine calculates the extinction coefficients per unit column density 
;it also outputs the extinction cross section in the B band in cm^2/H

pro kappagrain_H, model, qyear, kappatf, filters, kappatf_ab, kappatf_sc,$
		lambda, kappat_ab 
;kappagrain_H,'wd01','06'
;kappagrain_H,'pop00','01'
;kappagrain_H,'test','bb'

;input parameters
;model - name of input dust model, e.g. 'wd01'
;qyear - e.g. '06'

;output parameters
;kappatf -    the extinction coefficients in cm^2/H given for selected filters
;filters -    the wavelength in AA of selected filters
;kappatf_ab - the absorption coefficients in cm^2/H given for selected filters
;kappatf_sc - the scatterign coefficients given for selected filters
;lambda     - the wavelength array in AA ranging from UV to submm where the 
;             absorption coefficients are calculated
;kappat_ab  - the absorption coefficients in cm^2/H calculated for the lambda 
;             array


;subroutines:
;fitzpatrick.pro
;kappagraincomp_H.pro

common dirdef,rootdir

rootdir = '../../'
mydir = rootdir+'cinman/clump/'

;Initialise dummy variables
ss = ''
dim_comp = 0L
namer1 = mydir+'/emission/outdata/grain_sizes'+model+'_q'+qyear+'.dat'
openr, unit1, namer1, /get_lun
readf, unit1, ss
readf, unit1, dim_comp	; Number of different grain composition
comp_array = strarr(dim_comp)
readf, unit1, ss
for ii = 0, dim_comp-1 do begin
	readf, unit1, ss
	comp_array(ii) = ss	; populate comp_array with the different grain compositions
endfor

dim = 1000	;dimension of wavelength

;specific wavelengths to output extinction coefficients
filters = [912.,1350.,1500.,1650.,2000.,2200.,2500.,2800.,3150, 3650.,4430.,$
	5640.,8090.,12590.,22000.,36000,45000,58000];,50000.] ;in AA
filters_nurad = [912.,1500.,2800.,3150.,4430.,5640.,8090.,12590.,22000.,2200.,$
	1350.,1650.,2000.,2500.,3650.,36000.,45000.,58000.]

dimfilt = n_elements(filters)

index_reorder = INDGEN(dimfilt)
; reorder the indexing to be printed with in correct order
for jj = 0 , dimfilt-1 do index_reorder(jj) = WHERE(filters EQ filters_nurad(jj))

;Initialise arrays for extinction for each composition for the fine wavelength sampling
kappa_ab_array = dblarr(dim_comp,dim)	; absoprtion coefficient
kappa_sc_array = dblarr(dim_comp,dim)	; scattering coefficient
kappa_ph_array = dblarr(dim_comp,dim)	; phase coefficient
g_array = dblarr(dim_comp,dim)	; anisotropy factor
kappa_pah = dblarr(dim)	;

;Initialise arrays for extinction for each composition for the wavelengths sampled by NUrad
kappaf_ab_array = dblarr(dim_comp,dimfilt)
kappaf_sc_array = dblarr(dim_comp,dimfilt)
kappaf_ph_array = dblarr(dim_comp,dimfilt)
gf_array = dblarr(dim_comp,dimfilt)

weight_array = fltarr(dim_comp)           ;weight for each comp

;Total extinction coefficients summed across grain composition
kappat_ab = dblarr(dim)
kappat_sc = dblarr(dim)
kappat_ph = dblarr(dim)

kappatf_ab = dblarr(dimfilt) ;as above but for wavelengths sampled by NUrad
kappatf_sc = dblarr(dimfilt)
kappatf_ph = dblarr(dimfilt)
gf = dblarr(dimfilt)

kappat = dblarr(dim)       ;total extinction
kappatf = dblarr(dimfilt)

;produce the wavelength vector between 1.12D2 and 10^7 Angstroem, with equal 
;bins in a logarithmic scale
a1 = Alog10(1.120D2)
a2 = ALOG10(1.d+7)
a3 = (a2-a1) / (dim-1)
lambda = 10.D^(findgen(dim)*a3+a1)

;Produce the Fitzpatrick extinction curve according to wavelengths in lambda
fitzpatrick,lambda/1.d+4,Kfitz
Kfitz = Kfitz/(5.8d+21 * 1.086d)

;initialise name of grain composition variable
comp = ''
qcomp_year = ''
;Loop through grain compositions to calculate extinction coefficients
for j = 0, dim_comp-1 do begin	;begin loop in grain composition
	readf,unit1,ss
	readf,unit1,comp	; name of grain composition
	readf,unit1,ss
	readf,unit1,qcomp_year	; qyear (year when model updated)
	readf,unit1,ss
	readf,unit1,dim_Draine	; number of grain sizes tabualted by Draine
	readf,unit1,ss
	readf,unit1,amin	; min grain size in microns
	readf,unit1,ss
	readf,unit1,amax	; max grain size in microns
	readf,unit1,ss
	readf,unit1,ig_amin	; index of amin in Draine sizes
	readf,unit1,ss
	readf,unit1,ig_amax	; index of amax in Draine sizes
	readf,unit1,ss
	readf,unit1,s_grain	; grain density in g/cm^3
	readf,unit1,ss
	readf,unit1,dim_size	; number of grain sizes
	readf,unit1,ss
	readf,unit1,weight	; weighting factor
	readf,unit1,ss
	readf,unit1,ss
	readf,unit1,ss
	readf,unit1,ss
	;initialise arrays holding grain size a and grain size distribution na
	a = dblarr(dim_size)	; grain size in [cm]
	na = dblarr(dim_size)	; number density of grains in [cm^-1 H^-1]
	
	for ii = 0,dim_size-1 do begin
		readf, unit1, format='(2(e13.6))', pass1, pass2
		a(ii) = pass1	;in cm
		na(ii) = pass2	;in (cm * H)^-1
	endfor
	
	kappagraincomp_H, comp, qcomp_year, dim_Draine, amin, amax, ig_amin, ig_amax,$
		dim_size, a, na, lambda, dim, kappa_ab, kappa_sc, kappa_ph, g, Q_abs

	kappa_ab_array(j,*) = kappa_ab(*)
	kappa_sc_array(j,*) = kappa_sc(*)
	kappa_ph_array(j,*) = kappa_ph(*)
	g_array(j,*) = g(*)
	
	inb = where(filters eq 4430., nb)
	
	pass_ab = interpol(kappa_ab,lambda,filters)
	kappaf_ab_array(j,*) = pass_ab(*)
	pass_sc = interpol(kappa_sc,lambda,filters)
	kappaf_sc_array(j,*) = pass_sc(*)
	pass_ph = interpol(kappa_ph,lambda,filters)
	kappaf_ph_array(j,*) = pass_ph(*)
	pass_g = interpol(g,lambda,filters)
	gf_array(j,*) = pass_g(*)
	weight_array(j) = weight
endfor	; end loop in grain composition
free_lun, unit1

;sumation over composition

for j = 0, dim_comp-1 do begin 
	kappat_ab(*) = kappat_ab(*) + weight_array(j) * kappa_ab_array(j,*) 
	kappat_sc(*) = kappat_sc(*) + weight_array(j) * kappa_sc_array(j,*)
	kappat_ph(*) = kappat_ph(*) + weight_array(j) * kappa_ph_array(j,*)
	;same but at specific wavelengths
	kappatf_ab(*) = kappatf_ab(*) + weight_array(j) * kappaf_ab_array(j,*)
	kappatf_sc(*) = kappatf_sc(*) + weight_array(j) * kappaf_sc_array(j,*)
	kappatf_ph(*) = kappatf_ph(*) + weight_array(j) * kappaf_ph_array(j,*)
endfor
 
;total extinction
kappat = kappat_ab + kappat_sc
kappatf = kappatf_ab + kappatf_sc

;write dustmodel file to the nurad/indata/ directory.
Name = mydir+'/nurad/indata/dustmodel_'+model+'_q'+qyear+'_clump.in'
openw,unit21,Name,/get_lun
printf,unit21,format='(I2.0)', dimfilt
for jj = 0,dimfilt-1 do begin
	kk = index_reorder(jj)
	printf,unit21,format='(f6.0,3(d7.3))', filters(kk), kappatf(kk)/kappatf(10),$
		kappatf_sc(kk)/kappatf(kk), kappatf_ph(kk)/kappatf_sc(kk)
endfor
free_lun,unit21
print, Name

savename = mydir+'/nurad/saves/clump_kappa_component.save'
wave = lambda
save, wave, kappat, kappa_ab_array, kappa_sc_array, kappa_pah, Kfitz, filename=savename

!p.thick=3
!x.thick=4
!y.thick=4
!p.charthick=3
!p.charsize=1.5
set_plot,'ps'
Name = mydir+'/emission/figures/kappa_component_grain'+model+'_q'+qyear+'.ps'
print,Name
device,filename=Name

plot_oo,lambda/1.d+4,1.d+21 * kappat,$
        linestyle=0, xrange=[912.e-4,3.d],xstyle=1,ystyle=1,yrange=[0.03,2.4],$
        xtitle='!4k!3 [!4l!3m]',$
        ytitle='C!3!Iext!N [10!E-21!N N!IH!N!E-1!N cm!E2!N]'
oplot,lambda/1.d+4,Kfitz*1.d+21,linestyle=3
for j = 0,dim_comp-1 do begin
if j lt 2 then begin
	oplot,lambda/1.d+4,(kappa_ab_array(j,*)+ kappa_sc_array(j,*))* 1.d+21,$
		linestyle=j+1
endif
plots,[0.8d,1.1d],[2.d,2.d],linestyle=3,/data
xyouts, 1.2d,1.85d, 'Fitzpatrick'
plots,[0.8d,1.1d],[1.4d,1.4d],linestyle=0,/data
xyouts,1.2,1.4,'C!3!Iext!N'
plots,[0.8d,1.1d],[1.0d,1.0d],linestyle=1,/data
xyouts,1.2,1.0,'C!3!Iext!N Si'
plots,[0.8d,1.1d],[0.7d,0.7d],linestyle=2,/data
xyouts,1.2,0.7,'C!3!Iext!N Gra'
if j gt 1 then begin
	kappa_pah = kappa_pah + kappa_ab_array(j,*)+ kappa_sc_array(j,*)
endif
if j eq 3 then begin
	oplot, lambda/1.d+4, kappa_pah * 1.d+21,linestyle=j+1
	plots,[0.8d,1.1d],[0.5d,0.5d],linestyle=j+1,/data
	xyouts,1.2,0.5,'C!3!Iext!N PAH'
endif
endfor

device,/close

print,'  wave      abs     scat    ext    ext/B albedo   g'
print,'  AA     cm^2/H   cm^2/H   cm^2/H             '
for jj = 0,dimfilt-1 do begin
    print,format='(f6.0,3(e9.2),1x,3(f6.3,1x))',filters(jj),$
          kappatf_ab(jj),kappatf_sc(jj),kappatf(jj),kappatf(jj)/kappatf(10),$
          kappatf_sc(jj)/kappatf(jj),kappatf_ph(jj)/kappatf_sc(jj)
endfor

Name = mydir+'/emission/outdata/kappagrain'+model+'_q'+qyear+'.dat'
openw,unit20,Name,/get_lun
printf,unit20,'  wave      abs     scat    ext   ext/B albedo    g'
printf,unit20,'  AA     cm^2/H   cm^2/H   cm^2/H             '
for jj = 0,dimfilt-1 do begin
    printf,unit20,format='(f6.0,3(e9.2),1x,3(f6.3,1x))',filters(jj),$
          kappatf_ab(jj),kappatf_sc(jj),kappatf(jj),kappatf(jj)/kappatf(10),$
          kappatf_sc(jj)/kappatf(jj),kappatf_ph(jj)/kappatf_sc(jj)
endfor
free_lun, unit20

set_plot,'ps'
Name = mydir+'/emission/figures/kappagrain'+model+'_q'+qyear+'.ps'
print,mydir+'/emission/figures/kappagrain'+model+'_q'+qyear+'.ps'
device,filename=Name

plot_oo, lambda/1.d+4, kappat*1.d+21, xtitle='!4k!3 [!4l!3m]',$
         ytitle='C!3!Iext!N [10!E-21!N N!IH!N!E-1!N cm!E2!N]',$
         linestyle=0, xrange=[912.e-4,3.d],xstyle=1,ystyle=1
oplot,  lambda/1.d+4, kappat_ab*1.d+21,linestyle=1
oplot, lambda/1.d+4, kappat_sc*1.d+21,linestyle=2       
oplot, lambda/1.d+4,Kfitz*1.d+21,linestyle=3
plots,[0.6d,0.8d],[2.d,2.d],linestyle=3,/data
plots,[0.6d,0.8d],[1.4d,1.4d],linestyle=0,/data
plots,[0.6d,0.8d],[1.0d,1.0d],linestyle=1,/data
plots,[0.6d,0.8d],[0.7d,0.7d],linestyle=2,/data
xyouts, 0.9d,1.85d, 'Fitzpatrick'
xyouts,0.9,1.4,'C!3!Iext!N'
xyouts,0.9,1.0,'C!3!Iabs!N'
xyouts,0.9,0.7,'C!3!Iscat!N'
device,/close
set_plot,'ps'
Name = mydir+'/emission/figures/kappagrain'+model+'_q'+qyear+'_1_lambda.ps'
print, mydir+'/emission/figures/kappagrain'+model+'_q'+qyear+'_1_lambda.ps'
device,filename=Name

;output curves
save,lambda,kappat,kappat_ab,kappat_sc,Kfitz,filename='../outdata_intlum/attenuation_curves.xdr'

;plot_oo, 1./(lambda/1.d+4), kappat*1.d+21, $
plot, 1./(lambda/1.d+4), kappat*1.d+21, $
         xtitle='!4k!3!E-1!N [!4l!3m!E-1!N]',$
         ytitle='C!3!Iext!N [10!E-21!N N!IH!N!E-1!N cm!E2!N]',$
         yrange=[0.03,2.5],linestyle=0, xrange=[0.3d,10.d],xstyle=1,ystyle=1
oplot,  1./(lambda/1.d+4), kappat_ab*1.d+21,linestyle=1
oplot, 1./(lambda/1.d+4), kappat_sc*1.d+21,linestyle=2       
oplot, 1./(lambda/1.d+4),Kfitz*1.d+21,linestyle=3
;oplot,1./(lambda_in),0.93 * cext*1.d+21,linestyle=1
plots,[0.6d,0.8d],[2.d,2.d],linestyle=3,/data
plots,[0.6d,0.8d],[1.4d,1.4d],linestyle=0,/data
plots,[0.6d,0.8d],[1.0d,1.0d],linestyle=1,/data
plots,[0.6d,0.8d],[0.7d,0.7d],linestyle=2,/data
xyouts, 0.9d,1.85d, 'Fitzpatrick'
xyouts,0.9,1.4,'C!3!Iext!N'
xyouts,0.9,1.0,'C!3!Iabs!N'
xyouts,0.9,0.7,'C!3!Iscat!N'
device,/close
set_plot,'x'
!p.thick=1
!x.thick=1
!y.thick=1
!p.charthick=1
!p.charsize=1.

end
