;+
;Get the maps produced by NUrad and the data produced by the emission codes in master_clump and plot the integrated SED
;-
PRO plot_sed
;Standard compiler options
COMPILE_OPT idl2
;Get clump parameters and filenames
COMMON clump
;COMMON dirdef, rootdir
;Get directories from mydirectories.in
ss=''
urad_dir=''
emission_dir=''
rootdir=''
OPENR, dirunit, 'mydirectories.in', /GET_LUN
        READF, dirunit, ss
        READF, dirunit, urad_dir
        READF, dirunit, ss
        READF, dirunit, emission_dir
	readf, dirunit, ss
	readf, dirunit, rootdir
FREE_LUN, dirunit

;Define the parsec to meters conversion factor
m_per_pc=3.086D+16

;Convert the distance to the galaxy from Mpc to meters
distance_m = distance * 1D+6 * m_per_pc

;Conversion factor for flux [Jy] to luminosity [W/Hz]
flux_to_lum = 4d*!DPI*distance_m^2d * 1d-26     ; units [W/Hz]/[Jy]

;Get the IR luminosities from total_luminosity_*.dat:
OPENR, totlum, rootdir+emission_dir+'/outdata_intlum/'+totlum_filename, /GET_LUN
	;Skip to the number of wavelengths
	WHILE ss NE 'Number of wavelengths' DO READF, totlum, ss
	;Read the number of wavelengths
	READF, totlum, n_wavelengths_ir
	;Skip to the data
	READF, totlum, ss
	READF, totlum, ss
	
	;Ensure the number of wavelengths is an integer
	n_wavelengths_ir = LONG(n_wavelengths_ir)
	
	;Prepare arrays to store IR wavelengths and total luminosities
	ir_wavelengths = DBLARR(n_wavelengths_ir)
	ir_totlums = DBLARR(n_wavelengths_ir)
	;Read the IR data
	FOR j=0, n_wavelengths_ir-1 DO BEGIN
		READF, totlum, ir_wavelength, ir_totlum
		ir_wavelengths[j] = ir_wavelength
		ir_totlums[j] = ir_totlum
	ENDFOR
;Close the IR data file
FREE_LUN, totlum

;Concatonate UV-opt and IR data arrays
wavelengths=ir_wavelengths
totlums=ir_totlums
;energy in W/Hz * Hz
totenergy=totlums * 2.998D+18/wavelengths

;Save the results
SAVE, wavelengths, totlums, FILE='../outdata_intlum/'+sed_filename+'.xdr', DESCRIPTION="Wavelengths in Angstroms, total luminosities in W/Hz"

tot_flux = totlums / flux_to_lum

save, wavelengths, totlums, totenergy, tot_flux, filename = '../../nurad/saves/dust_sed.save'

;Set the filename for the plot
plot_name='../figures/'+sed_filename

;Make a title for the plot
;Convert parameters to strings
src=STRTRIM(STRING(FIX(rc*1E+3)),2)+'pc'
stau=STRTRIM(STRING(FIX(tau_clump)),2)
strunc=STRTRIM(STRING(FIX(rc * 1E+3 * clump_trunc)),2)+'pc'
sring=STRTRIM(STRING(FIX(rc * 1E+3 * clump_trunc * ring_factor)),2)+'pc'
plot_title="tau: "+stau+", rc: "+src+", truncation: "+strunc+", shell radius: "+sring

;Plot the result
@plot_begin.inc 
	cgPlot, wavelengths/1D+4, totenergy,$
		/XLOG, /YLOG, $
		TITLE=plot_title,$
		XTITLE="Wavelength, $\mu$m", $
		YTITLE="Luminosity, W/Hz x Hz"
@plot_end.inc

END

