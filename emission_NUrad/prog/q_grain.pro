;the program gives you the emissivity coefficient as a function of size and
;wavelength (from Draine)

pro q_grain,comp,qcomp_year,lambda,Q_abs,Q_sc,Q_ph

;common dirdef,rootdir
rootdir='../../'
Staub_dir=rootdir+'/Draine/'
name = Staub_dir+comp+qcomp_year
openr,unit,name,/get_lun

w=''
for i=0,3 do begin
	readf,unit,w
endfor
;Number of grain sizes
dim_Draine = long(strmid(w,0,4))
readf,unit,w
;Number of wavelengths sampled
dim_wave = long(strmid(w,0,4))

;Initialise arrays
lambda=fltarr(dim_wave)
Q_abs=fltarr(dim_Draine,dim_wave)
Q_sc=fltarr(dim_Draine,dim_wave)
Q_ph=fltarr(dim_Draine,dim_wave)

;Begin loop through grain size 
for i = 0, dim_Draine-1 do begin
	;skip header
	for k = 0, 2 do begin
		readf,unit,w
	endfor
	;Initialise variables
	b=0.
	c=0.
	d=0.
	e=0.
	;Loop through wavelengths sampled
	for j = 0, dim_wave-1 do begin
		readf, unit, b, c, d, e
		lambda[j] = b
		Q_abs[i,j] = c
		Q_sc[i,j] = d
		Q_ph[i,j] = e
	endfor
endfor

free_lun, unit
end
