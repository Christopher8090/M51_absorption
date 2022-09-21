;this function calculates the energy density of the radiation field in the
;solar neighborhood, following Mathis,Metzger and Panagia (1983)
; THE ROUTINE GIVES THE ENERGY DENSITY PER WAVELENGTH BIN, IF NUMBER
; OF WAVELENGTH LARGER 1.
; -------------------------------------------------------------------------
; HISTORY: 3.6.03 WRITTEN BY JOERG FISCHERA BASED ON u_l_ph.pro


function u_l_isrf,lambda,chi=chi                 ;UV_flux         ;,UV_wave
  ; Result in E-10 J/m^3/A

  ; CONSTANTS
  DELTA = 0.001   ; MINIMAL BIN SIZE FOR INTEGRATION
  c_e=2.99792458D    ;E+8  m/s speed of light
  el_e = 1.60217733D-19  ; electron charge


  IF N_ELEMENTS(chi) EQ 0 THEN BEGIN
    chi=1.D0 ;gives the strength of the ISRF in units of the solar neighborhood
           ;chi=1 corresponds to the solar neighborhood
  ENDIF


  ;the UV component taken from MMP (1982)
  UV_flux=[7.6,6.9,8.26,9.17,10.4,12.4,15.1,20.4,20.4,17.4,14.8,12.4,10.8]
  UV_wave=[0.25,0.246,0.23,0.216,0.20,0.18,0.16,0.134,0.110,0.105,0.100,0.095,$
	0.0912]
  UV_flux=UV_flux*1.D-8       
  UV_wave=UV_wave*1.D4

  ;calculate the 3 diluted blackbody radiation fields plus the 2.9 K component
  ;give the dilution factor vector
  Wi=[1.D+2,1.D-12,1.D-11,4.D-11]
  chiarr = [1.,chi+FLTARR(3)]
  ;give the black body temperature vector
  Ti=[2.7,7500.,4000.,3000.]
  NLAMBDA=N_elements(lambda)
  Ul=DBLARR(NLAMBDA)

  DLAMBDA = ABS(LAMBDA[1:NLAMBDA-1] - LAMBDA[0:NLAMBDA-2])
  loglambda = ALOG10(LAMBDA)
  dloglambda = loglambda[1:NLAMBDA-1] - loglambda[0:NLAMBDA-2]

  FOR i=0,NLAMBDA-2 DO BEGIN

    Ni = LONG(ABS(dloglambda[i]/DELTA))>2

    LAi = 10.^(FINDGEN(Ni)/(Ni-1.D0)*DLOGLAMBDA[i]+LOGLAMBDA[i])
    dLAi = LAi[1:Ni-1]-LAi[0:Ni-2]
    ULi = DBLARR(Ni)
    nr = WHERE(LAi GT 912 AND LAi LT 2500.D0, N)
    IF N GT 0 THEN BEGIN
       ULi[nr] = Wi[0]*planck_l(LAi[nr],Ti[0])
       ULi[nr] = ULi[nr]+chi*interpol(UV_flux,UV_wave,LAi[nr])/(4.*!pi)
    ENDIF
    nr = WHERE(LAi GE 2500.D0, N)
    IF N GT 0 THEN BEGIN
       FOR j=0,3 DO ULi[nr]=ULi[nr]+chiarr[j]*Wi[j]*planck_l(LAi[nr],Ti[j])
    ENDIF
    
    ULi = 0.5*(ULi[0:Ni-2]+ULi[1:Ni-1])
    Ul[i] = TOTAL(ULi*dLAi)        

  ENDFOR

  ;print,TOTAL(Ul)*4.*!DPI/c_e/el_e*1.D-16
  ; transform to E-10 J/m^3/A
  Ul = Ul / DLAMBDA
  Ul = [Ul,Ul[NLAMBDA-2]]

  ;the total energy density of the radiation field
  Ul=Ul*4.*!pi/c_e
  ;nr = WHERE(lambda LE 4000.,N)
  ;IF N GT 0 THEN Ul(nr)=0.
  return,Ul

end

