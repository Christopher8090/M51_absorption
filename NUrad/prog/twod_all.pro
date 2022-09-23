PRO twod_all;,tau,sfr,sfr4,sfr6,sfr7,old,old3,old5,bd
compile_opt idl2

read_scaling, model, qyear, scaabs, tau, nsersic, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5

start
morph_op=['td','tdi','tdo','d','di','do','b','ntd']
print, '------------------------------------------------maps.pro-------------------------------------------------'
for i_op=0, 6 do begin
	maps, model, qyear, morph_op[i_op], tau, scaabs
endfor
print, 'Done: maps.pro'
print, ''
print, '----------------------------------------radiation_fields_unit.pro----------------------------------------'
for i_op=0, 6 do begin                                  ;every component
	radiation_fields_unit, model, qyear, morph_op[i_op], tau, scaabs, nsersic
endfor
print,'Done: radiation_fields_unit.pro'
print, ''

;twod_phot,model,qyear,scaabs,tau,sfr,sfr4,sfr6,sfr7,old,old3,old5,bd,f_uv,f_uv4,f_uv6,f_uv7,f_BVIK,f_BVIK3,f_BVIK5

end
