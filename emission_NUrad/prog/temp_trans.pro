pro temp_trans,comp,model,qyear,qcomp_year,dim_size,aa,lu,ul,lambda_q,qq_abs,$
               tt_q,qp_grain,$
               stau,suv,sold,sbd,dir,inclbigbang=inclbigbang,energy_abs,$
               check


; *********************************************************
; *** program to calculate temperature distribution of grains ***
; ---------------------------------------------------------
; INPUT: comp: grain composition
;        model - dust model, e.g. 'wd01'
;        qyear - e.g.'06'
;        qyear_comp
;        dim_size - number of grain sizes in the input dust model
;        aa - grain size in micron
;	 lu  : wavelength of photon density (in ANGSTROEM)
;	 ul  : photon density in erg/pc^3/AA
;        lambda_q - wavelength of q array in AA
;        qq_abs - absorption emissivities
;        sbd - string identifier for the bulge-to-disk ratio to use for the
;              file name 
;        suv - string identifier for the fraction of non-ionizing UV photons to
;              use for the file name
;        sold - string identifier for the scaling of the old stellar populations 
;                   stellar populations to use for the file name
;        stau - string identifier for opacity to use for the file name        
;        dir  - save directory for the files
;	 inclbigbang : if set, then radiation of "big-bang" is included

; ---------------------------------------------------------
; HISTORY
; 15.5.2002 written by Joerg Fischera
;	    based on temp_vert.pro, written by Cristina Popescu
; 11.05.2006 updated by C. Popescu
; *********************************************************
; METHOD:
; PROGRAM DERIVES FOR ALL GRAIN SIZES THE TEMPERATURE
; DISTRIBUTION, STARTING WITH LARGEST GRAIN SIZES.
; IF THE GAUSSIAN DISTRIBUTION IS A GOOD APROXIMATION TO THE P_E DISTRIBUTION
; AND IF THE SIGMA>MINVAR (E.G. 0.5 K), THEN THE DISTRIBUTIONS ARE CALCULATED
; USING THE GAUSSIAN APPROXIMATIONS (SEE VOIT '91).
; IF THE GAUSSIAN DISTRIBUTION IS NOT A GOOD APROXIMATION TO THE P_E DISTRIBUTION
; THEN THE DISTRIBUTIONS ARE CALCULATED BY A NUMERICAL INTEGRATION METHOD
; INVOLVING ALSO ANALYTICAL SOLUTIONS

; A SAVE CONDITION CAN BE CHOSEN, SO THAT (IN CASE THAT
; THE GAUSSIAN DISTRIBUTION IS APPROPRIATE) ONLY
; TEMPERATURE DISTRIBUTIONS WITH A VARIANCE LARGER THAN
; A GIVEN VALUE ARE WRITTEN ON DISC. THE DISTRIBUTIONS,
; NUMERICALLY DERIVED, ARE SAVED ANYWAY.
; *********************************************************
; USED ROUTINES:
; temp_equil2_ph.pro
; temp_vert_plus_ph.pro
; temp_sub_gauss2_ph.pro
; temp_sub_gauss_distr.pro
; q_grain.pro
; planck_grain.pro
; u_l_ph2.pro

;get the root directory
common dirdef,rootdir

; DIRECTORY OF THE FILES
dirfile=rootdir+dir

; EXTRA SAVE CONDITION
MINVAR=0.5

ii=61L                    ; number of bins+1 for the temperature distribution

COMMON temp_const, m_el_e, el_e, stb_e, c_e, h_e, k_e

c_e=2.99792458D    ;E+8  m/s
el_e=1.602189D     ;E-19 Coulomb
m_el_e=9.10953     ;E-31 kg
h_e=6.626176D      ;E-34 J/s
k_e = 1.38066D     ;E-23 J/K
stb_e=5.6703D      ;E-8  W/m^2/K^4



i = dim_size
energy_abs = dblarr(dim_size)
perc=0.1   ; used for search of the minimum/ maximum temperature
	   ; (dT1=perc*Tmin, dT2=perc*Tmax)
Tmax=500.  ; maximum temperature used for interpolation, to derive
	   ; the equilibrium temperature
gaussian='yes' ; in the beginning of the calculation
;automatic = 1  ; find automatically where gaussian approximation not 
	       ; appropriate
;SAVECOND=0     ; savecondition, in the beginning 0, automatically
	       ; changed to 1, if VARIANCE of pT curve > MINVAR,
	       ; DISTRIBUTION HAS TO BE GAUSSIAN LIKE, OR IF GAUSSIAN
	       ; not appropriate
VOIT = 1       ; for numerical integration use analytical approximations
R0=0.D0
R2=0.D0


; ------- NAME OF THE FILE --------
Name=''
Name = 'PT_'+comp+'_'+model+'_'+qyear
if check eq 'no' then begin
   Name = Name+'_t'+stau+'_s'+suv+'_o'+sold+'_bd'+sbd+'.temp'
endif
if check eq 'yes' then begin
   Name = Name+'_isrf.temp'
endif



; -------- WRITE NEW FILE ---------
openw,unit,dirfile+Name,/get_lun
;print,'write ',dirfile+Name


REPEAT BEGIN
  ;if comp eq 'PAHneu' then print,'i',' gaussian',i, gaussian
  ; counting counter clockwise
  ; start with largest grain
  i = i - 1
  a = aa[i]

  ; the emissivities
  q_a = qq_abs[i,*]
  qp_grain_a = qp_grain[i,*]
  ;if i eq 25 then print,a,qp_grain[i,*]

   


    ; ----------s----------------------------------------------------
    ;               CALCULATE EQUILIBRIUM TEMPERATURE
    ; --------------------------------------------------------------

    temp_equil2_ph, a, comp, lu, ul, t_equil, lambda_q, q_a, tt_q, qp_grain_a, $
	R1, R0=R0, R2=R2, Tmax=Tmax, energy, inclbigbang=inclbigbang
    ;print,'equilibrium temperature', a, t_equil  ;,R0,R2
    ;if i eq 25 then print,a, t_equil,Tmax
     energy_abs(i) = energy
  IF gaussian EQ 'yes' THEN BEGIN

    ; --------------------------------------------------------------
    ;                GAUSSIAN DISTRIBUTION AFTER VOIT
    ; --------------------------------------------------------------

      ; ----------- calculate dispersion for gaussian -------------
      temp_sub_gauss2_ph, a, COMP, lu,ul, T_equil,$ 
        lambda_q, q_a, tt_q, qp_grain_a,$
	sig, E_equil, R2=R2,inclbigbang=inclbigbang

        vol = a^3*4./3.*!Dpi*10./el_e
        sigT = sig / c_grain(T_equil,comp) /vol
        sigT = sigT[0]
        ;print,'dispersion in temperature ',sigT
	; CONDITION FOR WRITING ON DISC
;	IF sigT GT MINVAR THEN SAVECOND=1

;      IF ((R1/R0 / E_equil) LT 0.01 AND 2*sigT/T_equil LT .05) 
;      and KEYWORD_SET(automatic) THEN BEGIN
     IF ((R1/R0 / E_equil) LT 0.01 AND 2*sigT/T_equil LT .05) THEN BEGIN
     
        temp_sub_gauss_distr, ii, a, comp, T_equil, E_equil, sig, TT_m, $
		U_m,prob, dT, dE
        T_min = TT_m(0)-0.5*(TT_m(1)-TT_m(0))
        T_max = TT_m(ii-2)+0.5*(TT_m(1)-TT_m(0))
        ; calcname='voit'

      ENDIF ELSE BEGIN

        gaussian = 'NO'
;	SAVECOND = 1
        ;print,a,t_equil
        sig = sig / c_grain(T_equil,comp) /vol
        sig=sig(0)
        T_min = (t_equil - sig*5) > t_equil/2
        T_max = (t_equil + sig*5) ;< t_equil*3./2.        
        dT2 = T_max * perc
        dT1 = T_min * perc
       
      ENDELSE

  ENDIF


  IF gaussian EQ 'NO' THEN BEGIN

      ; -----------------------------------------------------------
      ;             CALCULATE STOCHASTIC HEATING
      ; -----------------------------------------------------------

 
        jump2:
        ;print,'a, sig, T_min,T_max',a, sig, T_min,T_max   ;,pl,pr
        ;call_procedure, "temp_vert_plus_ph", $
	;  ii, aa[i], comp, lu, ul, T_min, T_max,$ 
  	;  lambda_q, q_a, tt_q, qp_grain_a, TT_m, U_m, prob, dT, dE, $
	;  energy, VOIT=VOIT, inclbigbang=inclbigbang
        call_procedure, "temp_vert_plus_ph_neu", $
	  ii, aa[i], comp, lu, ul, T_min, T_max,$ 
  	  lambda_q, q_a, tt_q, qp_grain_a, TT_m, U_m, prob, dT, dE, $
          energy,VOIT=VOIT, inclbigbang=inclbigbang,Tback=2.7
          ;if comp eq 'PAHion' and i eq 0 then print, 'TT_m', a[i], TT_m, U_m, prob, dT, dE
          ;print, 'a, prob', a, prob[0], prob[ii-2]
          energy_abs(i) = energy
          ;print,energy_abs(i) in J/s
          dT1 = perc*T_min
          dT2 = perc*T_max

          ;print,'left,right',prob[0],prob[ii-2]

          IF PROB[0] GT 1.D-12 THEN BEGIN
            T_min = T_min-dT1
	    goto, jump2
          ENDIF

          if (prob(0) lt 1.E-60 ) AND (T_min+dT1 LT T_EQUIL) then begin
            T_min = (T_min+dT1)
            goto, jump2
          endif

          IF PROB[ii-2] GT 1.D-12 THEN BEGIN
            T_max = T_max+dT2
	    GOTO, jump2
          ENDIF


  ENDIF

  ;PLOT_OO,TT_m,prob

  ; ------------------------------------------------------------
  ; ----------------------write data----------------------------
  ; ------------------------------------------------------------
;first write the data starting with the biggest grains, which exhibit
; equilibrium temperature. The format is kept the same as those for stochastic
; heating. 
 if (gaussian EQ 'yes') then begin
    printf,unit
    printf,unit,i
    printf,unit,'; Temperature of grains ('+comp+')'
    printf,unit
    printf,unit,'; grain size in micron'
    printf,unit, aa[i]
    printf,unit,'; T-Interval (Kelvin)'
    printf,unit,t_equil,t_equil
    printf,unit,'; number of points : '
    ;two and one are introduced to have the same format as for stochastic
    ;heating 
    two = 2
    printf,unit,two
    printf,unit
    printf,unit,'    <T>     dp(T)'
    one = 1.0
    for k=0,two-2 do begin
      printf,unit,t_equil,one,one
    endfor    
  endif

;continue with the stochastic heating
  if (gaussian EQ 'NO') then begin
    printf,unit
    printf,unit,i
    printf,unit,'; Temperature of grains ('+comp+')'
    printf,unit
    printf,unit,'; grain size in micron'
    printf,unit, aa[i]
    printf,unit,'; T-Interval (Kelvin)'
    printf,unit,T_min,T_max
    printf,unit,'; number of points : '
    printf,unit,ii
    printf,unit
    printf,unit,'    <T>     dp(T)'
    for k=0,ii-2 do begin
      printf,unit,TT_m[k],prob[k],dT[k]
    endfor    

  endif

ENDREP UNTIL i EQ 0
;help,energy
free_lun,unit
return

end


