function u_pah,T,a

;number of carbon atoms
nc = 468.*(a/0.001)^3.


  ; h-c abundance in molecule
  IF nc LE 25 THEN hc=0.5 ELSE BEGIN
    IF nc LE 100 THEN hc=0.5/SQRT(nc*1./25.) ELSE hc=0.25
  ENDELSE

  ; frequencies [cm^-1]
  wj = [886.,1161.,3030.]
  ; temperatures
  ; out-of the plane
  theta_op = 863.
  ; in plane mode
  theta_ip = 2504.

  h_e = 6.6260755D0   ; 1E-34 Js , Planck constant
  c_e = 2.9979D0      ; 1.d8 m , 1d10 cm
  k_e = 1.380658D0    ; 1E-23 J/K, Boltzmann constant
  el_e = 1.60217733D0 ; 1.d-19 C

  const = h_e*c_e*1.d-24

  temp1 = T*0.
  FOR j = 0,2 DO BEGIN

     x = const * wj[j]/k_e*1.d23
     temp1 = temp1 + wj[j]*const / (exp((x/T)<40.)-1.)
  ENDFOR
  temp1 = temp1/el_e*1.d19

  ;temp2 = (nc-2.)*(k_e*theta_op*fn(T/theta_op,2))*1.D-4/el_e
  ;temp3 = (nc-2.)*2.*k_e*theta_ip*fn(T/theta_ip,2)*1.D-4/el_e

  ; for fn function
  NN = N_ELEMENTS(T)
  x1 =  T/theta_op
  x2 =  T/theta_ip
  fn1 = x1*0.
  fn2 = x2*0.

  nr1 = WHERE(x1 LT 0.1, n1)
  nr2 = WHERE(x1 GT 5.0, n2)
  IF n1 GT 0 THEN fn1[0:n1-1] = exp(0.183003)*x1[0:n1-1]^2.99973
  IF n2 GT 0 THEN $
      fn1[NN-1-n2+1:NN-1] = exp(-1.48293)*x1[NN-1-n2+1:NN-1]^1.0184
  IF NN-n1-n2 GT 0 THEN BEGIN
      xx = x1[n1:NN-1-n2]
      y = -1.73830+1.37605*ALOG(xx)-0.219485*ALOG(xx)^2+0.0537687*ALOG(xx)^3
      fn1[n1:NN-1-n2] = EXP(y)
  ENDIF

  nr1 = WHERE(x2 LT 0.1, n1)
  nr2 = WHERE(x2 GT 5.0, n2)
  IF n1 GT 0 THEN fn2[0:n1-1] = exp(0.183003)*x2[0:n1-1]^2.99973
  IF n2 GT 0 THEN $
      fn2[NN-1-n2+1:NN-1] = exp(-1.48293)*x2[NN-1-n2+1:NN-1]^1.0184
  IF NN-n1-n2 GT 0 THEN BEGIN
      xx = x2[n1:NN-1-n2]
      y = -1.73830+1.37605*ALOG(xx)-0.219485*ALOG(xx)^2+0.0537687*ALOG(xx)^3
      fn2[n1:NN-1-n2] = EXP(y)
  ENDIF

  temp2 = (nc-2.)*k_e*theta_op*fn1*1.D-4/el_e
  temp3 = (nc-2.)*2.*k_e*theta_ip*fn2*1.D-4/el_e

  ;print,NC*hc*temp1,temp2,temp3
  ;the factor 4 is needed to correct the wrong formula from Draine & Li 2001
  Epah = nc*hc*temp1+4.0* (temp2+temp3)
  ;PLOT_OO,T,nc*hc*temp1
  ;OPLOT,T,temp2,color=255
  ;OPLOT,T,temp3,color=256L^2*255

  return, Epah ;
end
