;this routine give the spectrum of the energy density field as given by angelos

; big bang radiation can be included...

; HISTORY, rewritten by CRISTINA POPESCU
;         15.5.2002 updated by Joerg Fischera, make routine shorter

function u_l_ph2,wave,lu,ul,inclbigbang=inclbigbang

; velocity of light
c_e = 2.99792458D    ;E+8  m/s

fact = 3.4d-57   ;factor that transform erg/pc^3 in J/m^3

;transform the energy density from erg/pc^3/A into J/m^3/A
energy = ul * fact

;produce the spectrum of the energy density by interpolating
lambda = lu 

sp = wave*0.D0
nr =  WHERE(wave GE MIN(lambda) AND wave LE MAX(LAMBDA),NN)
IF NN GT 0 THEN BEGIN
  ;make the interpolation for the input wave and set to zero any negative values
  sp[nr] = interpol(energy,lambda,wave[nr])>0.D0
  sp = sp*1.d+10 ;in units of 10^-10 J/m^3
ENDIF

IF KEYWORD_SET(inclbigbang) THEN BEGIN
   sp = sp + 100.*4.*!PI/c_e*planck_l(wave,2.7)
ENDIF

return, sp
end
