;this routines outputs the absorption, scattering and phase-function coss
;section for a specific grain composition
;output: kappa_ab,kappa_sc,kappa_ph,g,Q_abs
pro kappagraincomp_H, comp, qcomp_year, dim_Draine, amin, amax, ig_amin,$
			ig_amax, dim_size, a, na, lambda, dim,$
			kappa_ab, kappa_sc, kappa_ph, g, Q_abs

;input parameters:
;comp       - composition e.g. 'Si'
;qcomp_year  
;dim_Draine - the total numers of grains in the Draine's tables
;amin       - the minimum grain size in micron of the input dust model
;amax       - the maximum grain size in micron of the input dust model
;ig_amin    - the index of the minimum grain size of the input dust model in 
;             the original sizes tabulated by Draine
;ig_amax    - the index of the maximum grain size of the input dust model in 
;             the original sizes tabulated by Draine
;dim_size   - the total number of grains in the input dust model
;a          - the grain sizes in [cm] of the input dust model
;na         - the grain size distribution in (cm * H)^-1 of the input dust 
;             model
;lambda     - the input wavelength array
;dim        - the dimenstion of the input wavelength array lambda

;subroutines:
;q_grain.pro

kappa_ab = dblarr(dim)	;array with the absorption coefficient 
kappa_sc = dblarr(dim)
kappa_ph = dblarr(dim)

QQ = dblarr(dim_Draine,dim)
QQ_sc = dblarr(dim_Draine,dim)
QQ_ph = dblarr(dim_Draine,dim)
QQ_grain_const = dblarr(dim_Draine)
QQ_sc_const = dblarr(dim_Draine)
QQ_ph_const = dblarr(dim_Draine)


;calls the subroutine that reads the absorption, scaterring and phase 
;function coefficients as a function of grain size and wavelength from Draine
q_grain, comp, qcomp_year, lambda_Q, QQ_grain, QQ_grain_sc, QQ_grain_ph

ncheck = n_elements(QQ_grain[*,0])
if ncheck ne dim_Draine then begin
	print,'unexpected number of grain sizes; program stops'
	goto, mark1
endif
lambda_Q = lambda_Q * 1.E+4     ; in Angstroem

;calculate the absorption coefficient for our lambda
for i = 0, dim_Draine-1 do begin
	; Q for lambda>1000 microns
	QQ_grain_const[i] = lambda_Q[0]^2 * QQ_grain[i,0]
	;check where lambda < 1e7 \AA or 1000 microns
	nbr = WHERE(lambda lt 1.E7, n)
	;if any wavelengths < 1e7 interpolate QQ_grain interpolate to given lambda
	if n gt 0 then QQ[i,nbr] = interpol(QQ_grain[i,*],lambda_Q,lambda[nbr])
	;because lambda has the values monotonically increasing
	if n lt dim then QQ[i,n:(dim-1)] = QQ_grain_const[i]/lambda[n:(dim-1)]^2
endfor

;calculate the scattering coefficient for our lambda
for i = 0, dim_Draine-1 do begin
	; Q for lambda>1000 microns
	QQ_sc_const[i] = lambda_Q[0]^2 * QQ_grain_sc[i,0]
	nbr = WHERE(lambda lt 1.E7, n)
	if n gt 0 then QQ_sc[i,nbr] = interpol(QQ_grain_sc[i,*],lambda_Q,lambda[nbr])
	;because lambda has the values monotonically increasing
	if n lt dim then QQ_sc[i,n:(dim-1)] = QQ_sc_const[i]/lambda[n:(dim-1)]^2
endfor

;calculate the scattering phase function for our lambda
for i = 0, dim_Draine-1 do begin
	; Q for lambda>1000 microns
	QQ_ph_const[i] = lambda_Q[0]^2 * QQ_grain_ph[i,0]
	nbr = WHERE (lambda lt 1.E7, n)
	if n gt 0 then QQ_ph[i,nbr] = interpol(QQ_grain_ph[i,*],lambda_Q,lambda[nbr])
	;because lambda has the values monotonically increasing
	if n lt dim then QQ_ph[i,n:(dim-1)] = QQ_ph_const[i]/lambda[n:(dim-1)]^2
endfor


Q_abs = dblarr(dim_size,dim)
Q_sc = dblarr(dim_size,dim)
Q_ph = dblarr(dim_size,dim)
Q_abs = QQ[ig_amin:ig_amax,*]
Q_sc = QQ_sc[ig_amin:ig_amax,*]
Q_ph = QQ_ph[ig_amin:ig_amax,*]

;element in grain size [cm]
da = a[1:dim_size-1] - a[0:dim_size-2]
nna = 0.5 * (na[1:dim_size-1] + na[0:dim_size-2])
;print, comp 
;print, 'total: ',total(nna * da), ' H^-1'
;print, '40%:   ',total(nna * da) *0.4, ' H^-1'
;print, '60%:   ',total(nna * da) *0.6, ' H^-1'
;print, ''
;make the integration of the absorption over the size distribution
;loop through wavelength
for l = 0, dim-1 do begin
	ia = na * Q_abs[*,l] * a^2
	iia = 0.5*(ia[1:dim_size-1]+ia[0:dim_size-2])
	integ = total(iia*da)
	kappa_ab[l] = !pi * integ ;in (nH)^-1 * cm^2
endfor

;make the integration of the scaterring over the size distribution
for l=0,dim-1 do begin
	ia_sc = na * Q_sc[*,l] * a^2
	iia = 0.5*(ia_sc[1:dim_size-1]+ia_sc[0:dim_size-2])
	integ = total(iia*da)
	kappa_sc[l] = !pi * integ  ;in (nH^-1) * cm^2
endfor

;make the integration of the scaterring phase function over the size 
;distribution
for l=0,dim-1 do begin
	ia_ph = na * Q_ph[*,l] * Q_sc[*,l] * a^2
	iia = 0.5*(ia_ph[1:dim_size-1]+ia_ph[0:dim_size-2])
	integ = total(iia*da)
	kappa_ph[l] = !pi * integ ;in nH^-1 * cm^2 
endfor

g = kappa_ph/kappa_sc

mark1:
end
