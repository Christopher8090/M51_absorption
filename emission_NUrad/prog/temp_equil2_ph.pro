PRO temp_equil2_ph, a, COMP, lu, ul, T_equil, lambda_Q, Q_a, T_Q_pl, Q_pl_a, $
	R1, Tmax=Tmax, energy, HEATING=HEATING, R0=R0,R2=R2,inclbigbang=inclbigbang

; *****************************************************************************
; ****  THIS PROGRAM CALCULATES THE EQULILIBRIUM TEMPERATURE FOR PHOTON  ******
; ** HEATING (CONTINUUM) OF A GRAIN WITH RELATIVE MOTION IN RESPECT TO A GAS **
; *****************************************************************************
; INPUT:
;	A	: GRAINSIZE IN MICRONS
;	COMP	: DEFINES THE GRAIN: 'Gra', 'Si', 'SiC', 'Iron'
;	LAMBDA_Q: WAVELENGTH FOR ABSORBTION COEFFICIENT
;	Q_A	: ABSORBTION COEFFICIENT OF GRAIN SIZE A
;	T_Q_PL	: TEMPERATURE OF PLANCK AVERAGED EMISSION COEFFICIENT
;	Q_PL_A	: PLANCK AVERAGED EMISSION COEFFICIENT OF GRAIN SIZE A
;	lu      : WAVELENGTH OF PHOTON DENSITY (ANGSTROEM)
;	ul	: PHOTON DENSITY
; KEYWORDS:
;	inclbigbang: RADIATION OF "BIG-BANG" IS INCLUDED
; 	Tmax	: MAXIMUM TEMPERATURE USED FOR INTRAPOLATION
; OUTPUT:
;	T_EQUIL	: EQUILIBRIUM TEMPERATURE
;	R1	; HEATING RATE [eV/s]
; ******************************************************************************
; HISTORY: WRITTEN BY JOERG FISCHERA (14.6.2001)
;	   WRITTEN FOR PHOTON HEATING (CONTINUUM) ONLY, BY JOERG FISCHERA
;		   15.5.2002
; ******************************************************************************


; ---------------------------------------------------------
; CONSTANT VALUES
; ---------------------------------------------------------

COMMON temp_const, m_el_e, el_e, stb_e, c_e, h_e, k_e


; -------------------------------------------------------------
; heating processes
; -------------------------------------------------------------

; energygains

R0_isrf=0.D
R1_isrf=0.D
R2_isrf=0.D
energyloss =0.D
T_best =0.D


  ; **** PHOTONHEATING (INTERSTELLAR BACKGROUND RADIATION) *****


    ; ------ ABSORBTION COEFFIENT Q -------
   
    ;calls the subroutine that reads the absorption coefficient
    ;Q_si_a,a,lambda_Q,QQ_si
    ;lambda_Q=lambda_Q*1.E+4     ; in Angstroem
    ;plot_oo,lambda_Q,QQ_si*lambda_Q/a
    QQ_a=Q_a

    ;  Q for lambda>1000 microns
    QQ_const=lambda_Q(0)^2*QQ_a(0)

    ; ------ CALCULATE THE ENERGY GAIN PER SECOND ------------
    ; METHODE: SIMPSONS RULE

    ; DEFINE THE WAVELENGTH ARRAY
    l_min = MIN(lu)
    a1 = Alog10(l_min)
    l_max = MAX(lu)
    IF KEYWORD_SET(inclbigbang) THEN l_max = 2.9D7/2.7*100.D0
    a2 = Alog10(l_max)
    ANZ = 2000
    ;the widths of the intervals in wavelength, in a logarithmic scale
    a3=(a2-a1)/(double(anz)-1.)
    ;produce a vector with wavelength for integration
    la=10.D^(findgen(anz)*a3+a1)
    la_m = (LA(0:ANZ-2)+LA(0:ANZ-1)) /2. 

    ; energies in eV
    En = h_e * c_e / el_e / la *1000.
    En_m = h_e * c_e / el_e / la_m*1000.

    ;call the subroutine that calculates the energy density of the radiation
    ;field for the corresponding wavelength vector 
    kurve = u_l_ph2(la,lu,ul,inclbigbang=inclbigbang)

    QQ = dblarr(anz)
    nbr = WHERE( la LT 1.E7, n)
    if n GT 0 THEN QQ[nbr] = interpol(QQ_a,lambda_Q,la[nbr])
    IF n LT anz THEN QQ[n:anz-1] = QQ_const/la[n:anz-1]^2
    f = QQ * 1.E5 * kurve * c_e / el_e / En

    ;check energy absorbed 
     deltala = dblarr(anz)
     deltala = la[1:anz-1]-la[0:anz-2]
     ener = QQ * kurve
     energ = 0.5*(ener[1:anz-1]+ener[0:anz-2])
     ;energy = pi * a^2[micron] * 1.d-12 * c_e * 1.d+8 * u * QQ * wave;
     ;in 1.E-10 J/m^3 * m^2 * m/s = 1.E-10 * J/s
     energy = total(energ * deltala) * !Dpi * a^2 * c_e * 1.d-14
     ;print, energy

    ;plot_oo,la,f
    kurve = u_l_ph2(la_m,lu,ul,inclbigbang=inclbigbang)

    QQ = dblarr(anz-1)
    nbr = WHERE( LA_m LT 1.E7, n)
    if n GT 0 THEN QQ(nbr) = interpol(QQ_a,lambda_Q,la_m[nbr])
    IF n LT anz-1 THEN QQ[n:anz-2] = QQ_const/la_m[n:anz-2]^2
    f_m = QQ * 1.E5 * kurve * c_e / el_e / En_m
    dla = (la[1:anz-1] - la[0:anz-2])/6.

    
    ; derive R1
    f1 = f*En
    f1_m = f_m*En_m
    ff = f1[0:ANZ-2]+f1[1:ANZ-1]+4.*f1_m
      
    ;pass the value of the inverse of the collisional time for the j bin of
    ;energy to a vector 
    R1_isrf = TOTAL(ff*dla) * !Dpi * a^2 
    ;print,'heating isrf',heating_R1.isrf
    IF N_ELEMENTS(R2) GT 0 THEN BEGIN
      ; derive R2
      f1 = f*En^2
      f1_m = f_m*En_m^2
      ff = f1[0:ANZ-2]+f1[1:ANZ-1]+4.*f1_m
      R2_isrf = TOTAL(ff*dla) * !Dpi * a^2 
    ENDIF

    IF N_ELEMENTS(R0) GT 0 THEN BEGIN
      ; derive R2
      f1 = f*En^0
      f1_m = f_m*En_m^0
      ff = f1[0:ANZ-2]+f1[1:ANZ-1]+4.*f1_m
      R0_isrf = TOTAL(ff*dla) * !Dpi * a^2 
    ENDIF

energygain = R1_isrf

if energygain GT 0 THEN BEGIN

; ---------------------------------------------------------
; fit the right temperature for cooling
; ---------------------------------------------------------


; planck averaged emissivities

;call the routine that provides the averaged emissivity
;planck_grain_a,a,'Si',TT_Q,Qp_a
Qp_a = Q_pl_a
TT_Q=T_q_pl
Tmin=min(TT_Q)
;Qp_si_a contain the averaged emissivity and TT_Q contains the respective
;dust temperature 

; define an array for interpolation that is better than the data
T_ARR = fltarr(1000)
T_log = findgen(1000)*(ALOG10(1200)-ALOG10(Tmin))/999.+ALOG10(Tmin)

; T =< min(TT_Q) K 
const = Qp_a(0) / min(TT_Q)^2.0

; T = min(TT_Q) .. 1500 K
; make a spline
Nii = WHERE(Qp_a GT 0, NNN)
Qp_a_log = spline(ALOG10(TT_Q(0:NNN-1)),ALOG10(Qp_a(0:NNN-1)),T_log)

;plot,T_log,Qp_a_log
;oplot,ALOG10(TT_Q(0:NNN-1)),ALOG10(Qp_a(0:NNN-1)),psym=5
; now the fitting 
; interpolation method
T_min=1.

IF KEYWORD_SET(Tmax) THEN T_max = 100 > Tmax < 1500 ELSE T_max=TT_Q(NNN-1)
T_max = T_max*1.
diff = 10000.
slope = 0.
T_best = T_min+(T_max-T_min)/2.


  Qp = const * T_min^2.0
  y0 = 4.D-1 * !Dpi * a^2 * stb_e * T_min^4 *Qp / el_e
  xx = ALOG10(T_max)
  yy = interpol(Qp_a_log,T_log,xx)
  Qp = 10^yy(0)
  y1 = 4.D-1 * !Dpi * a^2 * stb_e * T_max^4 *Qp / el_e

  f0 = Alog(y0)-ALOG(energygain)
  f1 = ALOG(y1)-ALOG(energygain)

REPEAT BEGIN

  VERGL1=diff
  slope = slope+1.

  IF f0 LT f1 THEN BEGIN
    TT = T_min - (T_max-T_min) * f0 / (f1 - f0)

    IF TT le min(TT_Q) THEN BEGIN
      Qp = const * TT^2.0
    ENDIF ELSE BEGIN
      xx = ALOG10(T_max)
      yy = interpol(Qp_a_log,T_log,xx)
      Qp = 10.^yy(0)
    ENDELSE

    y2 = 4.D-1 * !Dpi * a^2 * stb_e * TT^4 *Qp / el_e
    f2 = ALOG(y2) - ALOG(energygain)

    IF f2 LT 0 THEN BEGIN
      f0=f2
      T_min=TT
    ENDIF ELSE BEGIN
      f1=f2
      T_max=TT
    ENDELSE

    diff = abs(y2-energygain)/energygain
    IF diff LT vergl1 THEN BEGIN
      T_best = TT
      energyloss = y2
    ENDIF

  ENDIF ELSE SLOPE=500

ENDREP UNTIL slope EQ 500 OR diff LT 1.E-6

;print,'grain',a,'gain:',energygain,' lost (FIT):',energyloss
;print,'   temp',T_best,'needed slopes',slope
;print
;IF slope eq 500 then print,'not converged, diff:',diff/y2*100.,'%' 

ENDIF ELSE print,'grain',a,' zero heating'
;print,energygain,y2,'  diff:',diff/y2*100.,'%'

;T_equil=0.
T_equil = T_best
R1 = energygain
R2 = R2_isrf
R0 = R0_isrf
;PRINT,'R0,R1,R2 ',R0,R1,R2
;HEATING=HEATING_R1

return
END
