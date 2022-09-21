function c_grain,T,comp

; 15.1.1999 Joerg Fischera, function written for arrays
; 10.5.1999 Joerg Fischera, iron included


if comp eq 'Si' then begin 
  ;(ergs/cm^3/K)
  C_a=[1.4E3,2.2E4,4.8E5,3.41E7]
  C_p=[2.,1.3,0.68,0.]
  C_b=[0.,50,150,500,100000]
  CC=0.
endif

if comp eq 'Gra' then begin
  ; (ergs/cm^3/K)
  C_a=[3.84E2,2.32E3,5.61E3,7.74E5,4.14E7]
  C_p=[2.,1.56,1.37,0.57,0.]
  C_b=[0.,60,100,470,1070,100000]
endif

if comp eq 'PAHion' then begin
    ; assume the same heat capacity as for graphite
  ; (ergs/cm^3/K)
  C_a=[3.84E2,2.32E3,5.61E3,7.74E5,4.14E7]
  C_p=[2.,1.56,1.37,0.57,0.]
  C_b=[0.,60,100,470,1070,100000]
endif

if comp eq 'PAHneu' then begin
    ; assume the same heat capacity as for graphite
  ; (ergs/cm^3/K)
  C_a=[3.84E2,2.32E3,5.61E3,7.74E5,4.14E7]
  C_p=[2.,1.56,1.37,0.57,0.]
  C_b=[0.,60,100,470,1070,100000]
endif

IF comp EQ 'SiC' THEN BEGIN
  ; ergs/cm^3/K
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

IF comp EQ 'Iron' THEN BEGIN
  ; call extra routine
  CC = c_grain_iron(T)

ENDIF ELSE BEGIN

  N = N_ELEMENTS(T)
  CC = FLTARR(N)

  FOR i = 0,N-1 DO BEGIN

    n=0
    repeat n=n+1 until (T(i) lt C_b(n))
    CC(i) = C_a(n-1)*T(i)^C_p(n-1)

  ENDFOR

  CC = CC / 10. 
ENDELSE

  return, CC  ;(J/m^3/K)

end
