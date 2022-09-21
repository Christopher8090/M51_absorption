function integrate_lum, wavelengths, lambda1, lambda2, lums
compile_opt idl2

;This function calculates the SFR approximating the integration of 
;luminosity as a function of wavelength by taking a Riemann Sum.
;Inputs:
;	- Wavelengths, array of wavelengths in [m] corresponding to some luminosity (must have the same size)
;	- lambda1, the wavelength at which integration begins
;	- lambda2, the wavelength at which integration ends
;	- lums, luminosities in [W/Hz] corresponding to the "wavelengths" array

error=-1
if n_elements(lums) ne n_elements(wavelengths) then error=0
if lambda1 lt min(wavelengths) then error=1
if lambda2 gt max(wavelengths) then error=2
if error eq 0 then begin
	print, "Error - wavelengths and Lums are different sizes."
	return, -999
	goto, abort
endif
if error eq 1 then begin
	print, "Error - lower wavelength is too small."
	goto, abort
endif
if error eq 2 then begin
	print,"Error - upper wavelength is too large."
	return, -999
	goto, abort
endif

c = 299792458d

;Redefine wavelength array to be bounded by the upper and lower limits input
wave = wavelengths[where(abs(wavelengths-lambda1) eq min(abs(wavelengths-lambda1))):$
		   where(abs(wavelengths-lambda2) eq min(abs(wavelengths-lambda2)))]
;Redefine luminosity array to match the upper and lower bounds
lums = lums[where(abs(wavelengths-lambda1) eq min(abs(wavelengths-lambda1))):$
	    where(abs(wavelengths-lambda2) eq min(abs(wavelengths-lambda2)))]
dim_wave = n_elements(wave)
;Convert wavelengths to frequencies [Hz]
freq = c / wave

;Difference between frequencies i and i+1
d_freq = freq[0:dim_wave-2] - freq[1:dim_wave-2]
;Mean luminosity between elements i and i+1
av_lum = (lums[0:dim_wave-2] + lums[1:dim_wave-1]) / 2d
;Take the Riemann sum to approximate integration
integrated_lum = total(av_lum * d_freq)

return, integrated_lum
abort:
end
