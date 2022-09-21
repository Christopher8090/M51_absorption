;This function calculates the Planck function for a given wavelength and
;temperature 

; 8.1.2002 Joerg Fischera

; INPUT Array lambda in Angstroem
;       temperature T in Kelvin

function planck_l,lambda,T
  ;lambda in E-10 m
  c_e=2.99792458D0    ;E+8  m/s
  h_e=6.626176D0      ;E-34 J/s
  k_e=1.38066D0       ;E-23 J/K
  ;        s=N_elements(lambda)
  ;        l_exp=double(fix(Alog10(lambda)))
  ;lambda=lambda*1.D0
  ll=lambda/1.E10    ; in meter        ;*10.D^(-l_exp)
  TT = T *1.D0
  xx=h_e*c_e/(k_e*ll*TT)*1.D-3
  ;xx = 1.43883/ll/TT

;  l_exp=-34.+16.+50.-5.*l_exp-10.
;  ;print,l_exp
;  ;print,xx
;  yy=xx*0.D0
;  for i=0,s-1 do begin
;    if (xx(i) lt 50) then begin
;      ;print,i,xx(i),l_exp(i)
;      yy(i)=1./(exp(xx(i))-1.)
;      ;yy(i)=yy(i)*1.
;      ;print,i,xx(i),yy(i)
;    endif else yy(i)=0.
;    if (l_exp(i) lt 30) then begin
;      yy(i)=yy(i)*10.D^(l_exp(i))
;    endif else yy(i)=0
;  endfor
;  ;print,ll


  nr = WHERE((xx LT 87) ,number)
  yy = xx*0.
  IF number GT 0 THEN BEGIN 
    yy(nr) =  1. / ( (100.*ll(nr))^5 * (exp(xx(nr))-1.) )
    yy(nr) = yy(nr)*2.*h_e*c_e^2 
    ;yy(nr) = yy(nr)*3.74185E-5
    yy(nr) = yy(nr)*1.E-18
  ENDIF
  ;print,lambda
  ;print,yy
  ;print,'hui'
  return,yy  ;J/m^2/s/A
end

