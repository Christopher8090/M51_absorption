PRO twod_all;,tau,sfr,sfr4,sfr6,sfr7,old,old3,old5,bd

root = '../'

read_scaling, model, qyear, scaabs, tau, nsersic, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5

start
morph_op=['td','tdi','tdo','d','di','do','b','ntd']
PRINT, '------------------------------------------------maps.pro-------------------------------------------------'
FOR i_op=0, 6 DO BEGIN
	maps, model, qyear, morph_op(i_op), tau, scaabs
ENDFOR
print, 'Done: maps.pro'
PRINT, ''

morph_op=['td','d','b','tdi','tdo','di','do','ntd']     ;every component
PRINT, '----------------------------------------radiation_fields_unit.pro----------------------------------------'
FOR i_op=0, 6 DO BEGIN                                  ;every component
   radiation_fields_unit, model, qyear, morph_op(i_op), tau, scaabs, nsersic
ENDFOR
PRINT,'Done: radiation_fields_unit.pro'
PRINT, ''
twod_phot,model,qyear,scaabs,tau,sfr,sfr4,sfr6,sfr7,old,old3,old5,bd,f_uv,f_uv4,f_uv6,f_uv7,f_BVIK,f_BVIK3,f_BVIK5

END
