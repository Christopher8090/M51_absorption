FUNCTION u_grain,T,comp,a,alpha=alpha

; 8.1.2002 Joerg Fischera

;this function calculates the integral over C(T)dT, namely the internal energy
;per unit volume
; the calculations are done for silicates and graphites

  el_e = 1.60217733D0 ; 1.d-19 C
  ; Si is default

  C_a=[1.4E3,2.2E4,4.8E5,3.41E7]      ;the heat capacity in erg/cm^3/K 
  C_p=[2.,1.3,0.68,0.]                ;the power coefficient of the power law
  C_b=[0.E0,50.E0,150.E0,500.E0,1.E6] ; Temp intervalls


IF comp EQ 'Gra' THEN BEGIN

  C_a=[3.84E2,2.32E3,5.61E3,7.74E5,4.14E7]
  C_p=[2.,1.56,1.37,0.57,0.]
  C_b=[0.E0,60.E0,100.E0,470.E0,1070.E0,1.E6]

ENDIF 

IF comp EQ 'SiC' THEN BEGIN

  ; choose alpha SiC
  IF KEYWORD_SET(alpha) THEN BEGIN
    C_a = [1.82E0, 4.07E-1, 3.02E2, 1.25E4, 5.63E5, 1.04E7]
    C_p = [3.03, 3.50, 2.03, 1.31, 0.641, 0.185]
    C_b = [0.E0, 25.E0, 90.E0, 175.E0, 300.E0, 600.E0, 1.E6]
  ENDIF ELSE BEGIN
    ; choose beta SiC
    C_a = [2.18E3, 5.57E-1, 3.84E2, 1.03E4, 5.55E5, 1.10E7]
    C_p = [0.670, 3.43, 1.98, 1.34, 0.644, 0.177]
    C_b = [0.E0, 20.E0, 90.E0, 175.E0, 300.E0, 600.E0, 1.E6]
  ENDELSE

ENDIF

IF comp eq 'PAHion' or comp eq 'PAHneu' then begin
   U=u_pah(T,a)
ENDIF ELSE BEGIN

;IF comp EQ 'Iron' THEN BEGIN
;  ; call extra routine
;  U = u_grain_iron(T);

;ENDIF ELSE BEGIN

  NI = N_ELEMENTS(T)
  CC = DBLARR(NI)

  FOR j=0,NI-1 DO BEGIN
    b = C_b
    nr = WHERE( C_b LT T(j), N)
    b(N) = T(j)

    FOR i = 0,N-1 DO BEGIN
      p = C_p(i)+1
      CC(j) = CC(j) + C_a(i)/p*(b(i+1)^p-b(i)^p)
    ENDFOR
  ENDFOR
  U = CC/10.  
  U = (1./el_e)*U*!Dpi*4.*a^3.*10./3.
ENDELSE

RETURN, U ;(in J/Coulomb and some factor)

end
