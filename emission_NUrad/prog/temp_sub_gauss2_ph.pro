PRO temp_sub_gauss2_ph, a, COMP, lu, ul, T_equil, $
	 lambda_Q, Q_a, T_Q_pl, Q_pl_a, $
	R2=R2, sig, E_equil,INCLBIGBANG=INCLBIGBANG

; *****************************************************************
; THIS PROGRAM CALCULATES THE TEMPERATURE DISTRIBUTION FOR PHOTON 
; HEATING FOR THE LARGE GRAIN LIMIT (VOIT, ApJ, 379, 122)
; *****************************************************************
; INPUT   : A   = GRAINSIZE IN MICRONS
;           COMP= DEFINES THE GRAIN: 'Gra', 'Si'
;	    lu  = WAVELENGTH OF PHOTON DENSITIES (VECTOR)
;	    ul  = PHOTON DENSITIES (VECTOR)
;           T_equil = equilibrium temperature
;	    R2       IF SET THEN INPUT VALUE IS USED, OTHERWISE R2
;	             WILL BE CALCULATED
; KEYWORD : INCLBIGBANG = IF SET, THEN RADIATION OF "BIG-BANG" 
;			  INCLUDED 
; OUTPUT  : TEMPERATURE DISTRIBUTION (GAUSS)
; ******************************************************************
; HISTORY: 15.5.2002 written by Joerg Fischera, based on the routine
;		     temp_sub_gauss2.pro
; ******************************************************************


; ---------------------------------------------------------
; CONSTANT VALUES
; ---------------------------------------------------------

COMMON temp_const, m_el_e, el_e, stb_e, c_e, h_e, k_e


; -------------------------------------------------------------
; calculation of R2 in VOIT, second energy moment
; -------------------------------------------------------------


IF KEYWORD_SET(R2) NE 1 THEN BEGIN

  ; energygains
  ; calculation in eV^2/s

  ; **** PHOTONHEATING (CONTINUUM RADIATION) *****


    ; ------ ABSORBTION COEFFIENT Q -------
   
    QQ_a=Q_a

    ;  Q for lambda>1000 microns
    QQ_const=lambda_Q[0]^2*QQ_a[0]

    ; ------ CALCULATE THE ENERGY GAIN PER SECOND ------------
    ; METHODE: SIMPSONS RULE

    ; DEFINE THE WAVELENGTH ARRAY
    l_min = 912.
    a1 = Alog10(l_min)
    l_max = max(lu)
    IF KEYWORD_SET(inclbigbang) THEN l_max = 2.9D7/2.7*100.
    a2 = Alog10(l_max)
    ANZ = 2000
    ;the widths of the intervals in wavelength, in a logarithmic scale
    a3=(a2-a1)/(double(anz)-1.)
    ;produce a vector with wavelength for integration
    la=10.D^(findgen(anz)*a3+a1)

    ;call the subroutine that calculate the energy density of the radiation
    ;field for the corresponding wavelength vector 
    kurve = u_l_ph2(la, lu, ul, inclbigbang=inclbigbang)

    QQ = dblarr(anz)
    nbr = WHERE( la LT 1.E7, n)
    if n GT 0 THEN QQ(nbr) = interpol(QQ_a,lambda_Q,la(nbr))
    IF n LT anz THEN QQ(n:anz-1) = QQ_const/la(n:anz-1)^2
    f = QQ * kurve / la

    LA1 = (LA(0:ANZ-2)+LA(0:ANZ-1)) /2. 
    kurve = u_l_ph2(la1, lu, ul, inclbigbang=inclbigbang)

    QQ = dblarr(anz-1)
    nbr = WHERE( LA1 LT 1.E7, n)
    if n GT 0 THEN QQ(nbr) = interpol(QQ_a,lambda_Q,LA1(nbr))
    IF n LT anz THEN QQ(n:anz-2) = QQ_const/LA1(n:anz-2)^2
    f1 = QQ * kurve / la1
 
    ff = f(0:ANZ-2)+F(1:ANZ-1)+4.*F1
    dla = (la(1:anz-1) - la(0:anz-2))/6.
      
    ;pass the value of the inverse of the collisional time for the j bin of
    ;energy to a vector 
    R2 = TOTAL(ff*dla) *!Dpi*a^2*c_e^2*h_e*1.D8 / el_e^2

  ; ************** END CONTINUUM RAIATION ***************

  R2 = R2

ENDIF
 

; -----------------------------------------------------------
; calculate sigma of gauss distribution
; -----------------------------------------------------------

  ; calculate the heat capacity of the equilibrium temperature
  ct = c_grain( t_equil, comp)
  
  ct=ct[0]
  ; calculate the planck averaged emissivity

  Qp_a = Q_pl_a
  TT_Q=T_q_pl
  ; T =< 10 K 
  const = Qp_a[0] / 10.^2.0

  IF T_equil le 10 THEN BEGIN
      Qp = const * T_equil^2.0    
  ENDIF ELSE BEGIN
      Qp = interpol(Qp_a,TT_Q,T_equil)
      Qp = Qp(0)
  ENDELSE

  ; make derivation of Qp
  IF T_equil LT 50 THEN dQp = 2.0 * const * T_equil^1. ELSE BEGIN
    ; easy extrapolation to the neighbours
    T_log = findgen(1000)*(ALOG10(1400)-ALOG10(10))/999.+ALOG10(10)
    ; T = 10 .. 500 K
    ; make a spline
    Nii = WHERE(Qp_a GT 0,NNN)
    Qp_a_log = spline(ALOG10(TT_Q(0:NNN-1)),ALOG10(Qp_a(0:NNN-1)),T_log)
    n = 0L
    REPEAT n=n+1 UNTIL T_log(n) GT ALOG10(T_equil)
    ;print,n,Qp_a(n)
    ;print,pot
    dQp = (Qp/T_equil)*(Qp_a_log(n)-Qp_a_log(n-1)) / (T_log(n)-T_log(n-1))
  ENDELSE

  dE = 4.*stb_e * T_equil^3 *( 4.*Qp + T_equil*dQp)

  sig = 0.5 * R2 * a * ct *1.D2 * 4./3. / dE
  sig = SQRT(sig)

  ; calculate the equilibrium energy
  ;Em = 4./3.*!Dpi*a^3 * u_grain(t_equil,comp) * 1.D1 / el_e
  Em = u_grain(t_equil,comp,a)  

  ;print,'sig: ',sig
  ;print,'Em : ',Em

E_equil = Em[0]


END
