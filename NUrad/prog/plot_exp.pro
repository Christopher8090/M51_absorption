pro plot_exp, hs
compile_opt idl2

; This program plots the function form of the exponential disk geometry. Currently all parameters are hard-coded
; EXCEPT for 'hs' which is part of the calling sequence of this program. The user may decide comment-out the
; hard-coded geometry parameters below in favour of adding them to the calling sequence. Likewise, the user
; may remove 'hs' from the calling sequence in favour of including it below.

; calling options:
;	- hs, scale length (gradient of the pure exponential region)

; restore SB profile to fit an exponential to
data = 'datafile_0.2274um.save'
restore, data

; geometry parameters
A = 1.		; this is an arbitrary value used to adjust the plots to overlap
;hs = 1.		; scale length (gradient of the pure exponential region)
rtruncate = 5.	; outer truncation radius
hstin = 1.	; inner truncation radius
hsin = 2.	; inner radius
xis = 0.	; defines gradient in the region hstin <= r < rtruncate
sharp = 0.01	; smoothing factor to define curve at r = rtruncate

; intialse arrays to hold the exponential profile
r = dindgen(200)/10.	; radial position from 0-20kpc in steps of 0.1kpc
y = dblarr(n_elements(r))

; define each value of y for each r given the above parameters
for i = 0, n_elements(r)-1 do begin ; begin loop in r
	; function to smooth the transition of the profile after rtruncate
	truncate = 0.5 * (1 - erf( (r[i]-rtruncate) / (sharp*rtruncate) ))

	; for the region hstin <= r < hsin use the following formula
	if r[i] ge hstin and r[i] lt hsin then begin
		y[i] = ( (r[i]/hsin)*(1.-xis) + xis ) * exp(-hsin/hs)
	endif

	; for the region r >= hsin use the pure exponential
	if r[i] ge hsin then begin
		y[i] = exp(-r[i]/hs)
	endif
	
	; apply smoothing factor calculated for position i
	y[i] = y[i] * truncate

endfor ; end loop in r

; define maximum y value to be plotted
y_max = max(av_int_rad,/NaN) * 1.2

plot, x_rad, av_int_rad, /ylog, max_value=y_max, min_value=y_max*1e-3
oplot, r, y

stop




end
