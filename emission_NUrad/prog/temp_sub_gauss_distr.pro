pro temp_sub_gauss_distr, ii, a, comp, T_equil, E_equil, sig, Tdistr, E,distr,$
	dT, dE

; 8.1.2002 Joerg Fischera

; calculates the distribution after VOIT in the range of 8 sigma

el_e=1.602189D     ;E-19 Coulomb

E = findgen(ii-1)*(sig*16.)/(ii-2.) + E_equil - 8.*sig
dE = E(1)-E(0)
E1 = findgen(ii)*dE + E(0) - 0.5*dE

distr = 1.D0/(SQRT(2.*!Dpi)*sig)*EXP(-0.5*(E-E_equil)^2/sig^2)

; now calculate the distribution in temperature
T = findgen(1000)*(0.6*T_equil)/999. + T_equil*(1.-.20)
ut = u_grain(T,comp,a) 
tE = interpol(T,ut,E)
tE1 = interpol(T,ut,E1)

;FOR i = 0,ii-2 DO distr(i) = distr(i) * c_grain(tE(i),comp)

distr = distr/double(TOTAL(distr))
Tdistr = tE
;wset,1
;plot,tE,distr
;wset,0

; intervalls
dT = fltarr(ii-2)
dE = fltarr(ii-2)
dT = TE1(1:ii-1)-TE1(0:ii-2)
dE = E1(1:ii-1) - E1(0:ii-2)

END
