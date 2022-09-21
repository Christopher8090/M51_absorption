pro temp_vert_plus_ph, ii, a, comp, lu, ul, T_min, T_max, lambda_Q, QQ_grain, $
	TT_Q, Qp_grain_a, T_P, U_P, PROP, dT, dE, energy, $
	VOIT=VOIT,inclbigbang=inclbigbang
; The program calculates the temperature distribution of a dust grain with
; radius  a
; UPDATED PROGRAM: temp_vert_plus.pro
; INPUT:    NUMBER OF INTERVALLS (ii-1)
;           GRAINSIZE a IN MICRON
;           COMPOSITION comp, STRING ('Si', 'Gra', 'SiC')
;           TEMPERATURE INTERVALL [T_min, T_max]
;           lambda_Q EMISSIVITY WAVELENGTH IN ANGSTROEM
;           QQ_grain EMISSIVITY
;           TT_Q TEMPERATURE ARRAY FOR PLANCK AVERAGED EMISSIVITY
;           Qp_grain_a PLANCK AVERAGED EMISSIVIY
;	    lu	       wavelength of photon density
;	    ul	       photon density
; KEYWORDS: CHOOSE THE HEATING PROCESS
;	    VOIT: CHOOSE VOIT APPROXIMATION (DEFAULT)
;	    INCLBIGBANG: INCLUDE RADIATION (2.7 K) FROM "BIG-BANG"
; OUTPUT:   DISTRIBUTION PROP
;           TEMPERATURE T_P CORRESPONDING TO P
; HISTORY:  8.10.1998 Joerg Fischera, calculation for planck averaged 
;			emissivities improved
;	    26.6.2001 Joerg Fischera, make routine faster
;           15.5.2002 Joerg Fischera, make routine running for photon 
;		      heating only (continuum)


; ------------------------------------------------------
; constants
; ------------------------------------------------------

COMMON temp_const, m_el_e, el_e, stb_e, c_e, h_e, k_e
AMU    = 1.66055D4    ;E-31 kg


; ------------------------------------------------------------------------
; ------------------------------- Model ----------------------------------
; ------------------------------------------------------------------------


; --------------------------------
; the temperature and energy scala
; --------------------------------

;call u_grain subroutine to calculate U_min and U_max
u_gr_min = u_grain(T_min,comp,a)
u_gr_max = u_grain(T_max,comp,a)
U_min = u_gr_min[0]
U_max = u_gr_max[0]


; set table for temperatures and thermal energies
a1=alog10(T_min)
a2=alog10(T_max)
;calculate the width of a bin in a logarithmic scale
a3=(1./DOUBLE(II-1))*(a2-a1)
Tarr = 10^(dindgen(II)*a3+a1)
UTarr = u_grain(Tarr,comp,a)

; set energy bin
U = 10.D0^(DINDGEN(ii)/(ii-1.)*(ALOG10(U_max)-ALOG10(U_min))+ALOG10(U_min))
U_m = 0.5*(U[0:ii-2]+U[1:ii-1])
dU = U[1:II-1]-U[0:II-2]

; set temperature bin
TT = INTERPOL(Tarr,UTarr,U)
TT_m = INTERPOL(Tarr,UTarr,U_m)


; ---------------------
; emission coefficients 
; ---------------------

;calculates the internal energy (eV) that corresponds to TT_Q
;U_T = t_nach_u_grain(TT_Q,comp)
;U_T = (1./el_e)*U_T*!Dpi*4.*a^3.*10./3.

;determine the averaged emissivities corresponding to our bins in temperature,
;by interpolating the values from Qp_grain_a 
Qp_grain_a_U = dblarr(ii)

; make interpolation only for points with T<Tmax(com)
  nr = WHERE(TT_m LT 10,N)
  Qp_grain_a_U = interpol(ALOG10(Qp_grain_a),ALOG10(TT_Q),ALOG10(TT_m))
  Qp_grain_a_U = 10.^Qp_grain_a_U 

; for T less than 10K there are no more tabulated values for the average
; emissivity, because they can be approximated by T^2.0 law, independent
; on grain size
; scale the average emissivities, by comparing with the values at 10 K. 

vergl1=interpol(Qp_grain_a,TT_Q,10.)
IF N GT 0 THEN BEGIN
  Qp_grain_a_U[nr] = vergl1[0] * (TT_m[nr]/10.)^2
ENDIF

; ------------------------------------------
; absorption coefficient at long wavelength
; (needed for heating by photons)
; ------------------------------------------

;  Q for lambda>1000 microns
QQ_grain_const=lambda_Q[0]^2*QQ_grain[0]
; QQ_grain=QQ_grain_const/lambda_Q^2


; --------------------------------------------------------------------------
;  calculation of the transition matrix BB
; --------------------------------------------------------------------------

; --------------------------------------------------------------------------
;***************************************************************************
;***************************** DEFINITIONS *********************************
;***************************************************************************
;* epsilon	 ARRAY OF MINIMUM ENERGIES THE INTEGRAL EQUATION
;*		 IS SOLVED USING A FIRST ORDER INTEGRATION (see Voit '91)
;* Rn (n=0,1,2)	 MOMENTS OF THE HEATING RATE (see Voit '91)
;* Ren (n=0,1,2) MOMENTS OF THE HEATING RATE INTEGRATED TILL epsilon
;* BB		 TRANSITION MATRIX (first order)
;* Uarr		 
;* Rdarr	 


AA = DBLARR(ii-2,ii-1)
BB = DBLARR(ii-2,ii-1)
R0 = 0.D0
R1 = 0.D0
R2 = 0.D0
Re0 = DBLARR(ii-2)
Re1 = DBLARR(ii-2)
Re2 = DBLARR(ii-2)


epsilon=U_m[1:ii-2]-U[1:ii-2]

Uarrmin = MIN([0.01,MIN(epsilon)/20.])
Uarrmax = MAX(epsilon)
anz = (FIX(ALOG10(Uarrmax)-ALOG10(Uarrmin))*50+2)<600
Uarr1 = 10.D^(DINDGEN(anz)/DOUBLE(anz-1.)*(ALOG10(Uarrmax)-ALOG10(Uarrmin))+ALOG10(Uarrmin))
Uarr2 = U_m[0:ii-2]-U[0]
;Uarr2 = U_m[1:ii-2]-U_m[0]
nr=WHERE(Uarr2 GT Uarrmax,NN)
Uarr = [Uarr1,Uarr2[nr]]
anz=anz+NN
Rdarr = DBLARR(anz)
; minimum energy step taken into account, default (eV):
Edep_min = 1.
; more general: minimum deposited heating energy

; ----------------------------------------------------------------

;----------- the contribution of continuum radiation -----------

;if KEYWORD_SET(CONT) then begin

    ; --------- create a table for heating by ISRF ----------
    ; the lowest energy used for integration (if bigbang is included)
    l_max = MAX(lu)
    IF KEYWORD_SET(inclbigbang) EQ 1 THEN BEGIN 
      ll_max = 2.9D7/2.7*100.  ; wavelength in angstroem
      U_dep_min = (U_m[1]-U[1])/20.
      l_max = 1.D3*h_e*c_e/(U_dep_min*el_e)
      l_max = MAX([l_max,ll_max])
    ENDIF
    ; the maximum energy
    ;U_dep_max = 13.6
    ;l_min = 1.D3*h_e*c_e/(13.65*el_e)
    l_min = MIN(lu)

    anz = 500
    la = 10.D^(dindgen(anz)/(anz-1.)*(ALOG10(l_max)-ALOG10(l_min))+ALOG10(l_min))
    la_m = 0.5*(la[0:anz-2]+la[1:anz-1])
    y0 = u_l_ph2(la,lu,ul,inclbigbang=inclbigbang)
    y0_m = u_l_ph2(la_m,lu,ul,inclbigbang=inclbigbang)
    ;y0 = u_l_isrf(la,chi=1)
    ;y0_m = u_l_isrf(la_m,chi=1)

    ; calculate Q_abs
    Q = DBLARR(anz)
    nbr = WHERE( la LT 1.E7, n)
    Q[nbr] = interpol(QQ_grain,lambda_Q,la[nbr])
    IF n LT ANZ THEN Q[n:anz-1] = QQ_grain_const/la[n:anz-1]^2
    Q_m = INTERPOL(Q,la,la_m)

    y = y0 * Q * la
    y_m = y0_m * Q_m * la_m
    yy = 4.*y_m + y[0:anz-2]+y[1:anz-1]
    dla = (la[1:anz-1]-la[0:anz-2])/6.

        ;check energy absorbed
    deltala = dblarr(anz)
    deltala = la[1:anz-1]-la[0:anz-2]
    ener = Q * y0
    energ = 0.5*(ener[1:anz-1]+ener[0:anz-2])
    energy = total(energ * deltala)
    ;print,'a, energy absorbed [E-10 *J/m^3]: ', a, energy

    G = DBLARR(anz-1)
    FOR i=0,anz-2 DO G[i] = TOTAL(yy[0:i]*dla[0:i])
    G = G * 100. / h_e
    Umin = 1.D3*h_e*c_e/(la_m*el_e)
    ;PLOT_OO,Umin,G,yr=[0.0000001,2]

    BB1 = DBLARR(ii-2,ii-1)
    AA1 = DBLARR(ii-2,ii-1)

    ;FOR j=0,ii-3 DO BEGIN
    ;  Umin1 = U_m[j+1:ii-2]-U_m[j]
    ;  GG = INTERPOL(G,Umin,Umin1)>0.
    ;  BB1[j:ii-3,j] = GG
    ;ENDFOR

    ; USE THE TABLE FOR ALL OTHER TRANSITIONS
    FOR j=1,ii-2 DO BEGIN
      Umin1 = U_m[j]-U[0:j]
      GG = INTERPOL(G,Umin,Umin1) >0.
      AA1[j-1,0:j-1] = (GG[1:j]-GG[0:j-1])/dU[0:j-1]

      ;Umin1 = U_m[j]-U_m[0:j-1]
      ;GG = INTERPOL(G,Umin,Umin1) >0.
      ;BB1[j-1,0:j-1] = GG
      BB1[j-1,0:j-1] = 0.5*(GG[0:j-1]*Umin1[0:j-1]+GG[1:j]*Umin1[1:j])/(U_m[j]-U_m[0:j-1])
      ;BB1[j-1,0:j-1] = 0.5*(GG[0:j-1]+GG[1:j])
    ENDFOR

    BB1 = -BB1* !DPI * a^2
    AA1 =  AA1* !DPI * a^2
    BB  = BB + BB1
    AA  = AA + AA1

    Rdarr = Rdarr + (INTERPOL(G,Umin,Uarr)>0.)*!DPI*a^2
    
  ; moments of the heating rates
      la_E=(h_e*c_e/(el_e*la))*1.E+3      ; Energy in eV
      la_m_E=(h_e*c_e/(el_e*la_m))*1.E+3      ; Energy in eV  
  
      ; first moment (1/s)
      R0 = R0+G[anz-2]*!DPI*a^2

   ;
   ;   Re0 = Re0+(G[anz-2]-INTERPOL(G,Umin,epsilon)>0.)*!DPI*a^2  
   ;
      ; secound moment (eV/s)
      y = y0 * Q * la * la_E * 100. / h_e
      y_m = y0_m * Q_m * la_m * la_m_E * 100. / h_e
      yy = 4.*y_m + y[0:anz-2]+y[1:anz-1]
      FOR i=0,anz-2 DO G[i] = TOTAL(yy[0:i]*dla[0:i])
      R1 = R1+G[anz-2]*!DPI*a^2
   ;   Re1 = Re1+(G[anz-2]-INTERPOL(G,Umin,epsilon)>0.)*!DPI*a^2 
   ;
      ; third moment (eV^2/s)
      y = y0 * Q * la * la_E^2 * 100. / h_e
      y_m = y0_m * Q_m * la_m * la_m_E^2 * 100. / h_e
      yy = 4.*y_m + y[0:anz-2]+y[1:anz-1]
      FOR i=0,anz-2 DO G[i] = TOTAL(yy[0:i]*dla[0:i])
      R2 = R2+G[anz-2]*!DPI*a^2
   ;   Re2 = Re2+(G[anz-2]-INTERPOL(G,Umin,epsilon)>0.)*!DPI*a^2

  ;PRINT,'R0,R1,R2 ',R0,R1,R2
  Edep_min = min(UARR)

;ENDIF


; -----------------------------------------------------
; 2. Energy gain from emission
; -----------------------------------------------------

   TRANS2 = stb_e*Qp_grain_a_U[1:ii-2]*TT_m[1:ii-2]^4
   TRANS2 = 4.D-1*!Dpi*a^2*TRANS2
   TRANS2 = TRANS2/(el_e*dU[1:ii-2])
   nr = (LINDGEN(ii-2)+1)*(ii-1)-1
   BB[nr] = TRANS2

; --------------------------------------------
; Calculation of the distribution
; --------------------------------------------

IF KEYWORD_SET(VOIT) NE 1 THEN BEGIN

  temp_calc_matrix,BB,P

ENDIF ELSE BEGIN

  ; moments of the heating rates 

  ; the rate per deposited energy Rd
  anz = N_ELEMENTS(Rdarr)
  RRd0 = [R0,Rdarr]
  Umin0= [0.,Uarr]

  dER  = [Uarr[0],Uarr[1:anz-1]-Uarr[0:anz-2]]
  Rd0  = (RRd0[0:anz-1] - RRd0[1:anz])/dER
  ERd  = 0.5*dER+Umin0[0:anz-1]

  ; derive the moments

  yy0 = DBLARR(anz)
  yy1 = DBLARR(anz)
  yy2 = DBLARR(anz)

  y = Rd0*dER
  FOR i = 0,anz-1 DO yy0[i] = TOTAL(y[0:i])
  Re0 = INTERPOL(yy0,Uarr,epsilon)>0.
  y = Rd0*ERd*dER
  FOR i = 0,anz-1 DO yy1[i] = TOTAL(y[0:i])
  Re1 = INTERPOL(yy1,Uarr,epsilon)>0.
  y = Rd0*ERd^2*dER
  FOR i = 0,anz-1 DO yy2[i] = TOTAL(y[0:i])
  Re2 = INTERPOL(yy2,Uarr,epsilon)>0.

  ;PLOT_OI,Uarr,yy0  ;,yr=[1.D-0,20]
  ;PLOT_OI,ERd,Rd0
  ;plot,Uarr,yy0/max(yy0),xr=[1,20]
 ; PLOT_OO,Umin0,RRd0/MAX(RRD0);,xr=[0.00001,100]
 ; OPLOT,epsilon,Re0/R0,psym=4
 ; OPLOT,epsilon,Re1/R1,psym=3
 ; OPLOT,epsilon,Re2/R2,psym=2
  ;PLOT_OO,U,dU/1.D7
;dfg
  ; cooling 
  Edot = stb_e*Qp_grain_a_U[0:ii-2]*TT_m[0:ii-2]^4
  Edot = 4.D-1*!Dpi*a^2*Edot/el_e
  Edotm= 0.5*(Edot[0:ii-3]+Edot[1:ii-2])

  ; calculate the constants of the analytic solution
  a1 = DBLARR(ii-2)
  b1 = DBLARR(ii-2)
  c1 = DBLARR(ii-2)
  l1 = DBLARR(ii-2)
  l2 = DBLARR(ii-2)

  nr = WHERE(Re2/R2 GT 1.D-6, NN)
  IF NN GT 0 THEN BEGIN
     a1[nr] =  (Edot[nr+1]-Re1[nr])/(Re2[nr]/2.)
     b1[nr] = -(R0-Re0[nr])/(Re2[nr]/2.)
     c1[nr] =  2.D0/Re2[nr]
     l1[nr] = -0.5D0*a1[nr]+SQRT(0.25*a1[nr]^2-b1[nr])
     l2[nr] = -0.5D0*a1[nr]-SQRT(0.25*a1[nr]^2-b1[nr])
  ENDIF  
  nr = WHERE(Re2/R2 LE 1.D-6, NN)
  IF NN GT 0 THEN BEGIN
     a1[nr] = Edot[nr+1]/R0
     l1[nr] = R0/Edot[nr+1]
  ENDIF

  ; DEFINE PROBABILITY VECTOR
  N =  N_ELEMENTS(BB(0,*))
  P = DBLARR(N)

  ; startpoint
  P[0] = 1.D0
  FOR j = 1,N-1 DO BEGIN

      ; ------ FIRST APPROXIMATION (Guhatakurta & Draine) -------
      P[j] = - double(TOTAL(BB[j-1,0:j-1]*P[0:j-1],2)) / BB[j-1,j] 

      ; -------- analytical approximation (Voit) ---------

      ; ---- source function Se(E) ------
      ; heating transition
      ; Utrans = U_m[j]-U_m[0:j-1]
      Aji = AA[j-1,0:j-1]
      ;OPLOT,Utrans,Aji,psym=5
      SeE = TOTAL(P[0:j-1]*Aji)

      ; for Re2/R2 ~ 0 solution of DEQ of 1. order
      ; for Re2/R2 > 0 solution of DEQ of 2. order

      ; set the integrand values
      ; choose the minimum stepsize smaller than the minimum 
      ; deposited energy
        ; energy intervall
        from = U[j]
        to   = U[j]+epsilon[j-1]

        ; choose a first minimum energy difference (dlog10)
        delta0 = ABS(ALOG10(1.+Edep_min/4./from))

        ; total number of intervalls
        NN = ROUND((ALOG10(to)-ALOG10(from))/delta0)
        NN = MAX([NN,4])<250

	; now the corresponding minimum energy difference (dlog10)
        delta = (ALOG10(to)-ALOG10(from))/DOUBLE(NN-1.)

        ; the energies of the energy intervall
        E  = 10.^(DINDGEN(NN)*delta+ALOG10(U[j]))
        ; the average eneries
        Em = 0.5 * (E[0:NN-2] + E[1:NN-1])

	; the energy bin sizes
        dE = E[1:NN-1] - E[0:NN-2]

        ; minimum energies for transitions to higher energies
        diff = E[NN-1] - Em

      IF Re2[j-1]/R2 LE 1.D-6 THEN BEGIN
        ; THE APPROXIMATE DISTRIBUTION FUNCTION
        f0 = p[j-1]/dU[j-1]
        fm = (f0-SeE/R0) * EXP((l1[j-1]*(Em-U_m[j-1]))<40.)+SeE/R0

        ;PRINT,j-1,'*',f0,l1[j-1],f0,SeE/R0
        ;sdf =0
      ENDIF ELSE BEGIN

        ; secound border condition: first derivertive at the lower border
        ff = 2./Re2[j-1]*(p[j]*BB[j-1,j]+(epsilon[j-1])*SeE+(Re1[j-1]-Edot[j])*p[j-1]/dU[j-1])
        ;if abs(a - 0.000630928) lt 1.e-5 then PRINT, j, ff
        ; COEFFICIENTS A,B,C OF THE SOLUTION
        C2 = -c1[j-1]/(b1[j-1]+1.D-12)*SeE
        f0 = p[j-1]/dU[j-1]
        D2 = l2[j-1]-l1[j-1]
        A2 = (l2[j-1]*(f0-C2)-ff)/D2
        B2 = (l1[j-1]*(f0-C2)-ff)/D2

        ; THE APPROXIMATE DISTRIBUTION FUNCTION
        fm = (A2*EXP((l1[j-1]*(Em-U_m[j-1]))<40.)-B2*EXP((l2[j-1]*(Em-U_m[j-1])))<40.)+C2

	fm = ( f0 * (l2[j-1]*EXP(l1[j-1]*(Em-U_m[j-1])<40)-l1[j-1]*EXP(l2[j-1]*(Em-U_m[j-1])<40) )/D2 - $
	     ff * ( EXP(l1[j-1]*(Em-U_m[j-1])<40) - EXP(l2[j-1]*(Em-U_m[j-1]))<40 ) / D2 )  + $ 
	     C2 * (1.-(l2[j-1]*EXP(l1[j-1]*(Em-U_m[j-1])<40)-l1[j-1]*EXP(l2[j-1]*(Em-U_m[j-1])<40))/D2 )
        ;PRINT,j-1,f0,l2[j-1],ff,SeE

        sdf=1
     
      ENDELSE


      yy = INTERPOL(RRd0,Umin0,diff)
      ;oplot,diff,yy,psym=4
      ;oplot,Em,fm*dU[j-1]
      ;PLOTS,U_m[j-1],p[j-1],psym=4
      ;PRINT,j-1,NN,fm[0],diff[0]

      yy = yy*fm*dE
      yy = double(TOTAL(yy))/BB[j-1,j]


      ; -------- THE FINAL SOLUTION --------

      p[j] = p[j]+yy

      ; renormalize everything to the highest probability
      P = P/(MAX(P)>1.D-12)

  ENDFOR

ENDELSE

P = P / double(TOTAL(P))

; average temperature
TT_av=double(TOTAL(TT_m*P))


; ------------------------------------------------------------
; ------------------------ the result ------------------------
; ------------------------------------------------------------

T_P = TT_m
U_P = U_m
PROP = p
dT = TT[1:ii-1] - TT[0:ii-2]
dE = U[1:ii-1] - U[0:ii-2]

; ------------------------------------------------------------
; ------------------- end of program -------------------------
; ------------------------------------------------------------

return
end








































