PRO twod_all
compile_opt idl2

; call 'read_scaling.pro' routine to read in the following variables
read_scaling, model, qyear, scaabs, tau, nsersic, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5

; run 'start.pro' to allow for the Coyote graphics
start

; list of morphological components to include
morph_op=['td','tdi','tdo','d','di','do','b']
print, '------------------------------------------------maps.pro-------------------------------------------------'
; run 'maps.pro' for each morphological component
for i_op=0, 6 do begin
	; takes the .dat map file and print 2D map based on viewing angle in 'gal_param.in'
	maps, model, qyear, morph_op[i_op], tau, scaabs
endfor
print, 'Done: maps.pro'
print, ''
print, '----------------------------------------radiation_fields_unit.pro----------------------------------------'
for i_op=0, 6 do begin
	; Calibrates the radaiation field energy density files in /out/ and prints them to /unit/
	radiation_fields_unit, model, qyear, morph_op[i_op], tau, scaabs, nsersic
endfor
print,'Done: radiation_fields_unit.pro'
print, ''

; produces the azimuthally averaged surface brightness profiles from 2D maps
twod_phot,model,qyear,scaabs,tau,sfr,sfr4,sfr6,sfr7,old,old3,old5,bd,f_uv,f_uv4,f_uv6,f_uv7,f_BVIK,f_BVIK3,f_BVIK5

end
