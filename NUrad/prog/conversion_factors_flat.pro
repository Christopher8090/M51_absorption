;this routine calculates the conversion factors needed to scale
;radiation fields calculated for a given surface brightness in
;magnitudes per arcsquare to a total magnitude corresponding to SFR=1
;and old=1
;the second stellar disk has a flat radial distribution out to a inner radius r1in

pro conversion_factors_flat,factord_old, factord1_young,factord3_old, factord4_young,factord5_old, factord6_young, factord7_young

common param, dist, factor, zs, hs, zs1, hs1, hsin, hs1in, ztrun, rtrun, ztrun1, rtrun1, xis0, xis1, ll, ll1, S0, S0_opt, ms, $
                            zs3, hs3, zs4, hs4, hs3in, hs4in, rtrun3, rtrun4, xis3, xis4, $
                            zs5, hs5, zs6, hs6, hs5in, hs6in, rtrun5, rtrun6, xis5, xis6, $
                            hstin, hdtin, hs1tin, hd1tin, hs3tin, hd3tin, hs4tin, hd4tin, $
                            hs5tin, hd5tin, hs6tin, hd6tin, $
                            rtrund, rtrund1,rtrund3, rtrund4, rtrund5, rtrund6, $
                            zs7, hs7, hs7in, rtrun7, xis7, hs7tin

dist = dist * 1.d+6 * 3.086 * 1.d+16 ;in metres
zs = zs * factor ;in metres
hs = hs * factor ;in metres
zs1 = zs1 * factor ;in metres
hs1 = hs1 * factor ;in metres
hsin = hsin * factor ;in metres
hs1in = hs1in * factor ;in metres
ztrun = ztrun * factor ;in metres
rtrun = rtrun * factor ;in metres
ztrun1 = ztrun1 * factor ;in metres
rtrun1 = rtrun1 * factor ;in metres
;inner component
zs3 = zs3 * factor ;in metres
hs3 = hs3 * factor ;in metres
zs4 = zs4 * factor ;in metres
hs4 = hs4 * factor ;in metres
hs3in = hs3in * factor ;in metres
hs4in = hs4in * factor ;in metres
rtrun3 = rtrun3 * factor ;in metres
rtrun4 = rtrun4 * factor ;in metres
;outer component
zs5 = zs5 * factor ;in metres
hs5 = hs5 * factor ;in metres
zs6 = zs6 * factor ;in metres
hs6 = hs6 * factor ;in metres
hs5in = hs5in * factor ;in metres
hs6in = hs6in * factor ;in metres
rtrun5 = rtrun5 * factor ;in metres
rtrun6 = rtrun6 * factor ;in metres

;dust truncation
rtrund = rtrund * factor ;in metres
rtrund1 = rtrund1 * factor ;in metres
rtrund3 = rtrund3 * factor ;in metres
rtrund4 = rtrund4 * factor ;in metres
rtrund5 = rtrund5 * factor ;in metres
rtrund6 = rtrund6 * factor ;in metres

;inner trunc in metres
hstin = hstin * factor
hdtin = hdtin * factor
hs1tin = hs1tin * factor
hd1tin = hd1tin * factor
hs3tin = hs3tin * factor
hd3tin = hd3tin * factor
hs4tin = hs4tin * factor
hd4tin = hd4tin * factor
hs5tin = hs5tin * factor
hd5tin = hd5tin * factor
hs6tin = hs6tin * factor
hd6tin = hd6tin * factor

;nuclear component
zs7 = zs7 * factor ;in metres
hs7 = hs7 * factor ;in metres
hs7in = hs7in * factor ;in metres
rtrun7 = rtrun7 * factor ;in metres
hs7tin = hs7tin * factor

termz = 1.-exp(-ztrun/zs)
termr = exp(-hsin/hs)-exp(-rtrun/hs)+(hsin/hs)*exp(-hsin/hs)-(rtrun/hs)*exp(-rtrun/hs)
termz1 = 1.-exp(-ztrun1/zs1)
termr1 = exp(-hs1in/hs1)-exp(-rtrun1/hs1)+(hs1in/hs1)*exp(-hs1in/hs1)-(rtrun1/hs1)*exp(-rtrun1/hs1)

termz3 = 1.-exp(-ztrun/zs3)
termr3 = exp(-hs3in/hs3)-exp(-rtrun3/hs3)+(hs3in/hs3)*exp(-hs3in/hs3)-(rtrun3/hs3)*exp(-rtrun3/hs3)
termz4 = 1.-exp(-ztrun1/zs4)
termr4 = exp(-hs4in/hs4)-exp(-rtrun4/hs4)+(hs4in/hs4)*exp(-hs4in/hs4)-(rtrun4/hs4)*exp(-rtrun4/hs4)
termz5 = 1.-exp(-ztrun/zs5)
termr5 = exp(-hs5in/hs5)-exp(-rtrun5/hs5)+(hs5in/hs5)*exp(-hs5in/hs5)-(rtrun5/hs5)*exp(-rtrun5/hs5)
termz6 = 1.-exp(-ztrun1/zs6)
termr6 = exp(-hs6in/hs6)-exp(-rtrun6/hs6)+(hs6in/hs6)*exp(-hs6in/hs6)-(rtrun6/hs6)*exp(-rtrun6/hs6)
termz7 = 1.-exp(-ztrun1/zs7)
termr7 = exp(-hs7in/hs7)-exp(-rtrun7/hs7)+(hs7in/hs7)*exp(-hs7in/hs7)-(rtrun7/hs7)*exp(-rtrun7/hs7)

termrin = (4./3.)*!Dpi*((1.+xis0/2.)*hsin^2-(1-xis0) * hstin^3/hsin - (3./2.) * xis0 * hstin^2)*zs*termz*exp(-hsin/hs)
termr1in = (4./3.)*!Dpi*((1.+xis1/2.)*hs1in^2-(1-xis1) * hs1tin^3/hs1in - (3./2.) * xis1 * hs1tin^2)*zs1*termz1*exp(-hs1in/hs1)
termr3in = (4./3.)*!Dpi*((1.+xis3/2.)*hs3in^2-(1-xis3) * hs3tin^3/hs3in - (3./2.) * xis3 * hs3tin^2)*zs3*termz3*exp(-hs3in/hs3)
termr4in = (4./3.)*!Dpi*((1.+xis4/2.)*hs4in^2-(1-xis4) * hs4tin^3/hs4in - (3./2.) * xis4 * hs4tin^2)*zs4*termz4*exp(-hs4in/hs4)
termr5in = (4./3.)*!Dpi*((1.+xis5/2.)*hs5in^2-(1-xis5) * hs5tin^3/hs5in - (3./2.) * xis5 * hs5tin^2)*zs5*termz5*exp(-hs5in/hs5)
termr6in = (4./3.)*!Dpi*((1.+xis6/2.)*hs6in^2-(1-xis6) * hs6tin^3/hs6in - (3./2.) * xis6 * hs6tin^2)*zs6*termz6*exp(-hs6in/hs6)
termr7in = (4./3.)*!Dpi*((1.+xis7/2.)*hs7in^2-(1-xis7) * hs7tin^3/hs7in - (3./2.) * xis7 * hs7tin^2)*zs7*termz7*exp(-hs7in/hs7)

if hsin eq 0 then begin
	termrin(*)=0.
endif
if hs1in eq 0 then begin
        termr1in(*)=0.
endif
if hs3in eq 0 then begin
        termr3in(*)=0.
endif
if hs4in eq 0 then begin
        termr4in(*)=0.
endif
if hs5in eq 0 then begin
        termr5in(*)=0.
endif
if hs6in eq 0 then begin
        termr6in(*)=0.
endif
if hs7in eq 0 then begin
        termr7in(*)=0.
endif

;calculate the luminosity density per unit volume in the centre of the galaxy 
;in W/m^3/Hz

ll0 = ll/(4. * !Dpi * zs * hs^2*termz*termr + termrin)
ll01 = ll1/(4. * !Dpi * zs1 * hs1^2*termz1*termr1 + termr1in)
ll03 = ll/(4. * !Dpi * zs3 * hs3^2*termz3*termr3 + termr3in)
ll04 = ll1/(4. * !Dpi * zs4 * hs4^2*termz4*termr4 + termr4in)
ll05 = ll/(4. * !Dpi * zs5 * hs5^2*termz5*termr5 + termr5in)
ll06 = ll1/(4. * !Dpi * zs6 * hs6^2*termz6*termr6 + termr6in)
ll07 = ll1/(4. * !Dpi * zs7 * hs7^2*termz7*termr7 + termr7in)

;calculate the face-on luminosity density per unit projected area in the centre
;of the galaxy in W/m^2/Hz

lls0 = ll0 * zs * 2. * termz
lls01 = ll01 * zs1 * 2. * termz1
lls03 = ll03 * zs3 * 2. * termz3
lls04 = ll04 * zs4 * 2. * termz4
lls05 = ll05 * zs5 * 2. * termz5
lls06 = ll06 * zs6 * 2. * termz6
lls07 = ll07 * zs7 * 2. * termz7

;calculte the face-on luminosity density per unit projected area (corresponding
;to 1 arcsec^2) in W/Hz/arcsec^2
;the radius in metres corresponding to 1 arcsec
rr = 1. * dist/206265.	; 206265 is one [radian] in arcsec.

;the area in metres^2 corresponding to 1 arcsec^2
aa = rr^2
lls0 = lls0 * aa
lls01 = lls01 * aa
lls03 = lls03 * aa
lls04 = lls04 * aa
lls05 = lls05 * aa
lls06 = lls06 * aa
lls07 = lls07 * aa

;calculate the face-on brightness in W/Hz/m^2/arcsec^2
lls0 = lls0/(4. * !Dpi * dist^2)
lls01 = lls01/(4. * !Dpi * dist^2)
lls03 = lls03/(4. * !Dpi * dist^2)
lls04 = lls04/(4. * !Dpi * dist^2)
lls05 = lls05/(4. * !Dpi * dist^2)
lls06 = lls06/(4. * !Dpi * dist^2)
lls07 = lls07/(4. * !Dpi * dist^2)

;transform in Jy/arcsec^2
lls0 = lls0 * 1.d+26
lls01 = lls01 * 1.d+26
lls03 = lls03 * 1.d+26
lls04 = lls04 * 1.d+26
lls05 = lls05 * 1.d+26
lls06 = lls06 * 1.d+26
lls07 = lls07 * 1.d+26

;transform in magnitudes
m = -2.5 * alog10(lls0/S0_opt)
m1 = -2.5 * alog10(lls01/S0)	; magnitude relative to the zero point.
m3 = -2.5 * alog10(lls03/S0_opt)
m4 = -2.5 * alog10(lls04/S0)
m5 = -2.5 * alog10(lls05/S0_opt)
m6 = -2.5 * alog10(lls06/S0)
m7 = -2.5 * alog10(lls07/S0)

factord_old = 10.^(-0.4 * (m - ms))
factord1_young = 10.^(-0.4 * (m1 - ms))	; ratio of flux of m1/ms.
factord3_old = 10.^(-0.4 * (m3 - ms))
factord4_young = 10.^(-0.4 * (m4 - ms))
factord5_old = 10.^(-0.4 * (m5 - ms))
factord6_young = 10.^(-0.4 * (m6 - ms))
factord7_young = 10.^(-0.4 * (m7 - ms))
end
