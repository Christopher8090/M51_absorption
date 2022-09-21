pro dopita1,lambdaout,sedout,p,tdisp,lam1,lam2,ifudge
;
; calculates template SEDOUT of a star formation region
; based on model HII region emission calculated by Dopita et al 2006 (paper I).
; as tabulated in table 6 of that paper.
;
; Eg dopita1,lambdaout,sedout,1.0e7,1.0,0.1,10000.,1
;
; input: 
; LAMBDAOUT (array specifying wavelengths in micron, starting at shortest,
;            at which the SED is to be sampled). The available line-free range
;            in micron accommodated in table 6 is 0.0831792 to 813112. 
;            Since SED is monochromatic the user should ensure that
;            lambdaout is sufficiently finely sampled (lambda/deltalambda of
;            at least 100). 
; P          ISM pressure (P/k in cm-3.K), in range 10^4 to 10^8
;            Note that in Mike's formulation external pressure alone
;            determines the dynamics of the HII region, its radius and thus
;            the dust illumination (for a given central star cluster).
;            In reality it is the external density of
;            the medium into which an HII region expands,
;            rather than external pressure, which determines these quantities.
;            Mike's calculations were all done with an external density of
;            1cm^-3 and these yield too high radii and thus too low dust 
;            temperatures compared to the real situation of an HII region
;            expanding in a molecular environment with n = 10^4.
;            To mimic this real situation P of 10^7 or even 10^8 must
;            be used in Mike's model, despite the fact that the
;            actual pressure in the typical galaxies concerned is much lower.
;
; TDISP      exponential timescale for cloud dispersion (Myr). This
;            should be in the range 1 to 32. Setting it to a large
;            value means that the local dust emission is powered by
;            somewhat lower mass stars produced in the starburst, and
;            that the intensity of the UV light incident on the grains
;            in the HII region is less (corresponding to larger HII region
;            radii).  
;
; LAM1 LAM2  wavelength range in micron of a plot of calculated SED
; IFUDGE     if IFUDEGE = 1 then multiply Mike's curve with a function
;            to ensure it fits the data on galactic HII region complexes
;
; output:
; SEDOUT     SED in units of flux density (flux per frquency interval)
;            scaled such that the integral of the oversampled SED
;            from which SED is generated is unity over the 
;            0.08 micron to 81.3cm wavelength range. Note that SED
;            is a monochromatic function
;            
;
prepareplot
;
inname='dopita1.dat'
outname='sfr_sed.dat'
openr,inunit,inname,/get_lun
;print,'DOPITA1: reading ',inname
ndata=3162L ; number of input wavelength sampling points
tdispin=[1.,2.,4.,8.,16.,32.] ; timescales for cloud dispersion (Myr)
ntdisp=n_elements(tdispin)
pin=10000.0d0*[1.,100.,1000.] ; pressure (P/k in cm-3.K)
np=n_elements(pin)
nlambda=ndata/np ; number of input wavelength points - should be 1054
freq=fltarr(nlambda) ; to contain input frequency in Hz
lambda=freq ; to contain input wavelength in micron
s=fltarr(nlambda,ntdisp,np) ; flux density in mW/m2/Hz (exp. dispersion model)
fs=s ; to contain freq * flux density in mW/m2 (exp. dispersion model)
slin=fltarr(nlambda,np) ; flux density in mW/m2/Hz (linear dispersion model)
fslin=slin ; to contain freq * flux density in mW/m2 (linear dispersion model)
;
;
;
if(p lt min(pin) or p gt max(pin))then begin
 print,'DOPITA1: aborting - P out of range'
 print,' P, min(pin), max(pin) = ',p,min(pin), max(pin)
 goto, c999
endif ; p
if(tdisp lt min(tdispin) or tdisp gt max(tdispin))then begin
 print,'DOPITA1: aborting - TDISP out of range'
 print,' TDISP, min(tdispin), max(tdispin) = ',tdisp,min(tdispin), max(tdispin)
 goto, c999
endif ; tdisp
;
w=''
for i=0,30 do begin
 readf,inunit,w
endfor ; i
d0=0.0d0
d1=d0
d2=d0
d3=d0
d4=d0
d5=d0
dlin=d0
e0=d0
e1=d0
e2=d0
e3=d0
e4=d0
e5=d0
elin=d0
for ip=0,np-1 do begin ; loop in pressure
 for ilambda=0,nlambda-1 do begin
  readf,inunit,pdum,f,l,d0,d1,d2,d3,d4,d5,dlin,e0,e1,e2,e3,e4,e5,elin
  if(pdum ne pin[ip])then begin
   print,'Unexpected pressure found - aborting'
   print,'ip,pin[ip],ilambda,l = ',ip,pin[ip],ilambda,l
   goto, c999
  endif ; pdum ne p(ip)
  freq[ilambda]=f
  lambda[ilambda]=l
;  s(ilambda,0,ip)=d0 ; these get overflows
;  s(ilambda,1,ip)=d1
;  s(ilambda,2,ip)=d2
;  s(ilambda,3,ip)=d3
;  s(ilambda,4,ip)=d4
;  s(ilambda,5,ip)=d5
;  slin(ilambda,ip)=dlin
  s(ilambda,0,ip)=e0/f
  s(ilambda,1,ip)=e1/f
  s(ilambda,2,ip)=e2/f
  s(ilambda,3,ip)=e3/f
  s(ilambda,4,ip)=e4/f
  s(ilambda,5,ip)=e5/f
  slin(ilambda,ip)=elin/f
  fs(ilambda,0,ip)=e0
  fs(ilambda,1,ip)=e1
  fs(ilambda,2,ip)=e2
  fs(ilambda,3,ip)=e3
  fs(ilambda,4,ip)=e4
  fs(ilambda,5,ip)=e5
  fslin(ilambda,ip)=elin
 endfor ; ilambda
endfor ; ip
;
; mask known lines included in mappings (as listed in table 2 of paper 3)
;
mask=lonarr(nlambda)
lamline=[0.12157,0.19091,0.19112,0.23252,0.27979,0.3726,0.37287,0.37979]
lamline=[lamline,0.38354,0.38687,0.38886,0.38891,0.39674,0.39701,0.41047]
lamline=[lamline,0.43405,0.44715,0.48613,0.49588,0.50068,0.58756,0.63002]
lamline=[lamline,0.6548,0.65628,0.65833,0.66782,0.67163,0.67307,0.71357]
lamline=[lamline,0.77510,0.90693,1.00494,1.0830,1.0833,1.09381,1.28181]
lamline=[lamline,1.8751,2.6252,4.05,6.98,8.99,10.52,12.81,15.55,18.68,26.2] 
lamline=[lamline,21.0,33.6366,34.7941,51.7972,57.3845,63.184,88.3017,145.525]
lamline=[lamline,157.741,205.178,610.0]
; as well as the false PAHquery
lamline=[lamline,22.0]
; and other lines seen in Dopita spectra (some seem to be lines
; from above list with inaccurate wavelengths)
lamline=[lamline,0.9089,0.9568,1.0859]
lamline=[lamline,33.9,35.15,52.4,58.0,89.5,124.,160.,208.7,377.,611.]
lamline=[lamline,11200.,12200.,14600.,79700.,104200.,115000.,148800.]
lamline=[lamline,233000.,254500.,364000.,476000.,569000.,680000.,742000.]
;
nlamline=n_elements(lamline)
for ilamline=0,nlamline-1 do begin
 test=abs(lambda-lamline(ilamline)) - min(abs(lambda-lamline(ilamline)))
 ii=where(test lt 0.01*lamline(ilamline),nii) ; for 1 percent resolution
 if(nii gt 0)then mask[ii]=1L
endfor ; ilamline
;
ii=where(mask eq 0,nii)
lambda=lambda[ii]
ss=s
sslin=slin
s=double(s)
slin=double(slin)
nlambda=nii
s=fltarr(nlambda,ntdisp,np)
slin=fltarr(nlambda,np)
for ip=0,np-1 do begin
 for itdisp=0,ntdisp-1 do begin
  s(*,itdisp,ip)=ss(ii,itdisp,ip)
  slin(*,ip)=sslin(ii,ip)
 endfor ; itdisp
endfor ; ip
;
; subtract off synchrotron emission calculated
; as asynch*lambda^(alpha) with alpha and asynch 
; determined from the spectrum produced by dopita from 1cm to 1m
; (which we take to be free-free free, since there is no curvature).
;
alpha=1.093
asynch=fltarr(ntdisp,np)
synch=fltarr(nlambda,ntdisp,np)
asynch(0,0)=1.90
asynch(1,0)=1.795
asynch(2,0)=1.695
asynch(3,0)=1.623
asynch(4,0)=1.602
asynch(5,0)=1.555
asynch(0,1)=2.12
asynch(1,1)=2.107
asynch(2,1)=2.095
asynch(3,1)=2.089
asynch(4,1)=2.120
asynch(5,1)=2.098
asynch(0,2)=1.897
asynch(1,2)=1.979
asynch(2,2)=2.062
asynch(3,2)=2.126
asynch(4,2)=2.168
asynch(5,2)=2.193
asynch=asynch*(10.0^(19.0))
for itdisp=0,ntdisp-1 do begin
 for ip=0,np-1 do begin
  synch(*,itdisp,ip)=asynch(itdisp,ip)*((lambda/1000000.)^alpha)
 endfor ; ip
endfor ; itdisp
iiii=where(lambda lt 100.0,niiii)
synch(iiii,*,*)=0.0d0 ; don't subtract synchrotron shorter than 100 micron
                  ; since we expect a spectral break before this and
                  ; anticipate Mike would not have included short wavelength
                  ; synchrotron
s=s-synch
;
; In the residual there is a step at 1000 micron
; which seems incompatible with a pure dust emission. We deal with
; this by forcing the emission longwards of 500 micron to decrease
; as a modified (m=2) planck function with some (low) temperature
; tplanck, set to 10, 12 and 15 K for the three pressures so that
; the SED is continuous with that for lambda < 500 micron. This was found
; to give a better looking extrapolation than that provided by
; dust emitting with emission peop lambda^(-4) corresponding to warm
; dust from the HII region. In reality Mike's model cannot handle the
; cold dust emission from the optically thick cloud interiors since 
; photons are only traced to Av 20 in cloud skin. In any case it is
; too difficult (in practice impossible) to calculate this interior 
; emisison     a priori since
; (i) the ratio of skin emission to interior emission will depend on
; the sizes of the clouds that are present in the star formation complex and
; [ii] the grains in the cloud interior will have ice mantles of
; unknown composition so we don't know the optical properties of these grains. 
; Thus the only way of dealing with the long wavelength emission
; of the SF region template is through the use of empirical data.
; Thus, the solution provided here should be later adjusted so that 
; it reproduces the observed submm/FIR colours of the integrated
; emission of SF complexes in the MW and in neaby galaxies. In this
; way we might hope to retain the diagnosic power of the link between
; the colours of the localised warm dust emission and ISM pressure provided
; by Mike's model while still properly accounting for the part of
; the submm emission which is provided by the localised sources.
; 
;lambda_rj=500.0
lambda_rj=300.0 ; better fits submm data
ii=where(lambda ge lambda_rj,nii)
lambdalong=lambda[ii]
factorlong=fltarr(np,nii)
dummy=fltarr(nii)
ii0=ii[nii-1] ; shortest wavelength in lambdalong greater or equal to lambda_rj
for ip = 0,np-1 do begin
	; factorlong(*,ip)=(lambdalong/lambda[ii0])^(-4.0d0)
	if(ip eq 0)then tplanck=10.0
	if(ip eq 1)then tplanck=12.0
	if(ip eq 2)then tplanck=15.0
	tplanck=tplanck*2.0 ; to better fit submm data
	dummy[*]=(planck(10000.*lambdalong,tplanck)/planck(10000.*lambda[ii0],tplanck))
	factorlong[ip,*]=dummy[*]
endfor ; ip=0,np-1
;help,factorlong
;print,factorlong
;print,lambdalong
slong=fltarr(nii,ntdisp,np)
for itdisp=0,ntdisp-1 do begin
 for ip=0,np-1 do begin
  slong[*,itdisp,ip]=s[ii0,itdisp,ip]*factorlong[ip,*]
  s[ii,itdisp,ip]=slong[*,itdisp,ip]
 endfor ; ip
endfor ; itdisp
;
; at this point the s array contains dust and PAH emission and stellar
; emission, with no spectral lines. We can integrate to find the
; total bolometric luminosities
;
lumtot=fltarr(ntdisp,np)
for itdisp=0,ntdisp-1 do begin
 for ip=0,np-1 do begin
  lumtot[itdisp,ip]=total(s[*,itdisp,ip]/lambda[*])
;  print,'itdisp,ip, lumtot = ',itdisp,ip,lumtot(itdisp,ip)
 endfor ; ip
endfor ; itdisp
;


;
; finally we subtract the continuum emission from the optical/NIR/MIR
; side. We do this by assuming (i) that 90 percent of 
; all the emission in Mike's SED at 2.5
; micron is stellar emission, [ii] that the stellar emission 
; varies as a black body with 4000K.
; The dust/PAH emission at 3.0 micron and longer wavelengths is taken to be 
; Mike's calculated emission minus this stellar continuum.
; The dust/PAH emission shorter than 3.0 micron is taken to 
; be the dust/PAH continuum at 3.0 micron (ie the remaining continuum
; after subtraction of the stellar continuum) decreasing towards
; shorter wavelengths as a Back Body with temperature 500K
; 
; 2.5 micron stellar continuum
;save,pin,tdispin,lambda,s,slin,file='dopita1.xdr'
;
lamref=2.5
ii=where(abs(lambda-lamref) lt 0.05,nii)
ilamref=ii(0)
astar=fltarr(ntdisp,np)
tstar=astar
starlight=fltarr(nlambda,ntdisp,np)
for itdisp=0,ntdisp-1 do begin
 for ip=0,np-1 do begin
  astar[itdisp,ip]=total(s[ii,itdisp,ip])/float(nii)
  tstar[itdisp,ip]=4000.0
  starlight[*,itdisp,ip]=$
astar[itdisp,ip]*planck(10000.0*lambda,tstar[itdisp,ip])
  starlight[*,itdisp,ip]=starlight[*,itdisp,ip]*lambda*lambda/(lamref*lamref)
  starlight[*,itdisp,ip]=starlight[*,itdisp,ip]/planck(10000.0*lambda[ilamref],tstar[itdisp,ip])
  s[*,itdisp,ip]=s[*,itdisp,ip]-starlight[*,itdisp,ip]
 endfor ; ip
endfor ; itdisp
;
; 3.0 micron PAH continuum
;
lamref=3.0
ii=where(abs(lambda-lamref) lt 0.05,nii)
ilamref=ii[0]
apah=fltarr(ntdisp,np)
tpah=apah
iii=where(lambda le lamref,niii)
pahcont=fltarr(niii,ntdisp,np)
lampah=lambda[iii]
for itdisp=0,ntdisp-1 do begin
 for ip=0,np-1 do begin
  apah[itdisp,ip]=total(s[ii,itdisp,ip])/float(nii)
  tpah[itdisp,ip]=500.0
  pahcont[*,itdisp,ip]=$
apah[itdisp,ip]*planck(10000.0*lampah,tpah[itdisp,ip])
  pahcont[*,itdisp,ip]=pahcont[*,itdisp,ip]*lampah*lampah/(lamref*lamref)
  pahcont[*,itdisp,ip]=pahcont[*,itdisp,ip]/planck(10000.0*lambda[ilamref],tpah[itdisp,ip])
  s[iii,itdisp,ip]=pahcont[*,itdisp,ip]
 endfor ; ip
endfor ; itdisp
;
; patch region from 3.4 to 4.5 micron
;
lamref1=3.4
lamref2=4.5
ii1=where(abs(lambda-lamref1) lt 0.05,nii1)
ilamref1=ii1(0)
ii2=where(abs(lambda-lamref2) lt 0.05,nii2)
ilamref2=ii2(0)
cont1=fltarr(ntdisp,np)
cont2=fltarr(ntdisp,np)
for itdisp=0L,ntdisp-1 do begin
 for ip=0L,np-1 do begin
  cont1(itdisp,ip)=median(s(ii1,itdisp,ip))
  cont2(itdisp,ip)=median(s(ii2,itdisp,ip))
  i1=max(ii1)
  i2=min(ii2) ; longest wavelength
  i1=long(double(i1))
  i2=long(double(i2))
  di=i1-i2+1L
  imid=i1-di/2
;  help,i1,i2,imid
  x1=lambda[i1]
  x2=lambda[i2]
  xm=lambda[imid]
  s1=cont1[itdisp,ip]
  s2=cont2[itdisp,ip]
  sm=(s1+s2)/3.
  dum1=(x1*x1-xm*xm)*(x1-x2)
  dum2=(x1*x1-x2*x2)*(x1-xm)
  dum3=(s1-sm)*(x1-x2)
  dum4=(s1-s2)*(x1-xm)
  a=(dum3-dum4)/(dum1-dum2)
  b=((s1-s2)-a*(x1*x1-x2*x2))/(x1-x2)
  c=s1-(a*x1*x1)-(b*x1)
  for i=i2,i1 do begin
   x=lambda[i]
;   help,s,i,itdisp,ip,x,a,b,c
   s[i,itdisp,ip]=a*x*x + b*x + c
endfor ; i
 endfor ; ip
endfor ; itdisp
;
;save,pin,tdispin,starlight,synch,lambda,s,slin,file='dopita1.xdr'
;
;
; at this point we have a pure dust/PAH spectrum in s. We can now derive the
; time-averaged average F factor corresponding to each combination
; of itdisp and ip
;
fav=fltarr(ntdisp,np)
lumdusttot=fltarr(ntdisp,np)
for itdisp=0,ntdisp-1 do begin
 for ip=0,np-1 do begin
  lumdusttot[itdisp,ip]=total(s[*,itdisp,ip]/lambda[*])
  fav[itdisp,ip]=lumdusttot[itdisp,ip]/lumtot[itdisp,ip]
;  print,'itdisp,ip,lumtot,lumdusttot,fav = ',itdisp,ip,lumtot(itdisp,ip),$
;lumdusttot(itdisp,ip),fav(itdisp,ip)
; fav found to be in the range 0.69 (for 1 MYr and 10^7) 
; to 0.978 (for 32 Myr and 10^4). Thus in this model most of the radiation
; is absorbed locally, which is differerent from the real situation.
; Mike must have done something wrong here as the escape fraction should
; not depend on P. Perhaps, in addition, he has put in some attenuation
; into the UV-optical values
; 
 endfor ; ip
endfor ; itdisp
;
; create SED interpolated to input values of P and TDISP 
;
; first do interpolation in TDISP to produce SED1
;
sed1=fltarr(nlambda,np)
dum=fltarr(ntdisp)
work=fltarr(ntdisp,np)
tdispsampin=findgen(ntdisp)
tdispout=[tdisp]
tdispsampout=interpol(tdispsampin,tdispin,tdispout)
for ip=0L,np-1 do begin
 for ilambda=0L,nlambda-1 do begin
  work(*,ip)=s(ilambda,*,ip)
  dum(*)=work(*,ip)
  sed1(ilambda,ip)=interpolate(dum,tdispsampout,cubic=-0.5)
 endfor ; ilambda
endfor ; ip
;
;save,pin,tdispin,starlight,synch,lambda,s,sed1,file='dopita1.xdr'
;
; now interpolate in log(p)
;
sed=fltarr(nlambda)
dum=fltarr(np)
work=fltarr(np)
psampin=findgen(np)
pout=[p]
poutlog=alog10(p)
pinlog=alog10(pin)
;pout=10000.*(1.0*findgen(1000L)) ; for test
psampout=interpol(psampin,pinlog,poutlog)
 for ilambda=0L,nlambda-1 do begin
  work(*)=alog10(sed1(ilambda,*))
  result=interpolate(work,psampout,cubic=-0.5)
  sed(ilambda)=10.0^(result)
;  if(result lt -18.0)then sed(ilambda)=0.0
 endfor ; ilambda
;
; remove -NaN
seddum=sed
ii=where(seddum gt 0,nii)
sed(*)=0.0d0
sed[ii]=seddum[ii]
;
;window,0
;plot_oo,lambda,s(*,0,0),$
;xst=1,yst=1,xran=[lam1,lam2],yran=max(s)*[1.0e-06,1.1],/nodata
;itdisp=3
;for ip=0,np-1 do begin
;oplot,lambda,s(*,itdisp,ip)
;endfor
;oplot,lambda,sed,linest=1
;
;save,pin,tdispin,starlight,synch,lambda,s,sed,file='dopita1.xdr'
;
; plot result
;window,0
;plot_oo,lambda,s(*,0,0),$
;xst=1,yst=1,xran=[lam1,lam2],yran=max(s)*[1.0e-06,1.1],/nodata
;itdisp=0
;for ip=0,np-1 do begin
;oplot,lambda,s(*,itdisp,ip)
;endfor
;oplot,lambda,sed,linest=1
;
; normalise SED
;
lambdamin=min(lambda)
lambdamax=max(lambda) ; range of available lambda
;print,'DOPITA1: range of available wavelengths (micron) = ',lambdamin,lambdamax
if(float(min(lambdaout)) lt lambdamin or float(max(lambdaout)) gt lambdamax)then begin
 print,'DOPITA1: aborting - requested wavelength range unavailable'
 print,' min(lambdaout), max(lambdaout) = ',min(lambdaout), max(lambdaout)
 goto,c999
endif ; min(lambdaout) lt min(lambda) or max(lambdaout) gt max(lambda)
;
nlambdagrid=100000L ; create a highly oversampled interpolated array
lambdagrid=alog10(lambdamin)+$
(alog10(lambdamax)-alog10(lambdamin))*$
findgen(nlambdagrid)/float(nlambdagrid-1L)
lambdagrid=10.0d0^(lambdagrid)
lambdaident=findgen(nlambda)
lambdasamp=interpol(lambdaident,lambda,lambdagrid)
sedgrid=interpolate(sed,lambdasamp,cubic=-0.5)
; oversampled lambdagrid, sedgrid is now set up. This can be used
; for future applications which require oversampling
; such as convolution with filterpass functions
; but is only used here to obtain the normalisation.
; 
dlambdagrid=lambdagrid(1:nlambdagrid-1)-lambdagrid(0:nlambdagrid-2)
sgrid=(sedgrid(1:nlambdagrid-1)+sedgrid(0:nlambdagrid-2))/2.0
lgrid=(lambdagrid(1:nlambdagrid-1)+lambdagrid(0:nlambdagrid-2))/2.0
dnu=dlambdagrid/(lgrid*lgrid)
flux=total(sgrid*dnu)
sedgrid=sedgrid/flux
sed=sed/flux
;
; now interpolate the (now normalised) SED (sampled by LAMBDA) 
; onto the user given wavelength sampling LAMBDAOUT
;
nlambdagrid=n_elements(lambdaout)
lambdagrid=lambdaout
lambdaident=findgen(nlambda)
lambdasamp=interpol(lambdaident,lambda,lambdagrid)
sedout=interpolate(sed,lambdasamp,cubic=-0.5)
;
;save,lambda,s,sed,lambdaout,sedout,file='dopita1.xdr'
;
;window,0
plot_oo,lambdaout,sedout,$
xst=1,yst=1,xran=[0.8,3000.],yran=max(sedout)*[3.0e-05,5.0],$
xtit='Wavelength (micron)',ytit='Predicted flux density (Jy)'
;
; overplot data normalised to 100 micron point
;
lambdapivot=100.0
i100=where(abs(lambdaout-lambdapivot) eq min(abs(lambdaout-lambdapivot)),ni100)
sedoutref=sedout(i100(0))
;
hiidata,source,catsource,lambdasource,fluxsource,lumsource
;help,source,catsource,lambdasource,fluxsource,lumsource
nsource=n_elements(source) ; number of sources
nlambdasource=n_elements(lambdasource) ; number of sampled wavelengths
shii=fltarr(nsource,nlambdasource) ; fluxes at each wavelength
nhii=lonarr(nlambdasource) ; number of measurements at each wavelength
i100source=where(lambdasource eq lambdapivot,ni100source)
for isource=0,nsource-1 do begin
 ii=where(fluxsource(isource,*) gt 0.001,nii)
; print,source(isource),alog10(fluxsource(isource,ii))
 factor=sedoutref/fluxsource(isource,i100source) ; fluxratio data point 100um
; print,'sedoutref,factor= ',sedoutref,factor
 for i=0,nii-1 do begin ; lambdas for which source seen
  idum=ii(i)
;  help,idum,nhii
  nhii(idum)=nhii(idum)+1L ; number of samples so far at wavelength sample idum
  iii=nhii(idum)
;  help,iii
  shii(iii-1,idum)=factor(0)*fluxsource(isource,idum)
  xplot=lambdasource(idum)
  yplot=shii(iii-1,idum)
;  help,xplot,yplot
  oplot,[xplot],[yplot],psym=1
endfor ; i
endfor ; isource
;save,nhii,shii,file='hii.xdr'
;
; now multiply Mike's curve with the fudge function such that the
; model curve fits the data
;
if(ifudge eq 1)then begin
 nmin=1L
 ifit=where(nhii ge nmin,nfit)
 wfit=fltarr(nfit)
 sfit=wfit
 lamfit=sfit
 model=sfit
 ii=-1L
 for i=0,nlambdasource-1 do begin ; wavelength
  nlam=nhii(i)
  if(nlam ge nmin)then begin
   ii=ii+1L
   wfit[ii]=nlam
   sfit[ii]=median(shii(0:nlam-1,i))
   lamfit[ii]=lambdasource(i)
   iref=where(abs(lambdaout-lambdasource(i)) eq $
        min(abs(lambdaout-lambdasource(i))),niref)
   model[ii]=sedout(iref(0))
  endif ; nlam
 endfor ; i
;save,sfit,wfit,lamfit,model,file='sedmeas.xdr'

llamfit=alog10(lamfit)
lsfit=alog10(sfit)
lsedout=alog10(sedout)
llambdapivot=alog10(lambdapivot)
llambdaout=alog10(lambdaout)
wfit(*)=1.0
; polynomial log
;pivotfit,llamfit,lsfit,wfit,lsedout,llambdapivot,2,llambdaout,result
;result=10.0^(result)
; polynomial linear
;pivotfit,lamfit,sfit,wfit,sedout,lambdapivot,3,lambdaout,result
;gammaupmax=1.2 ; for log 1.2 - must be greater than unity
;gammaupmax=20.0
;lgammaupmax=alog10(gammaupmax)
; exp log
;exppivotfit1,llamfit,sfit,wfit,sedout,lgammaupmax,llambdapivot,llambdaout,result
; exp linear
;exppivotfit1,lamfit,sfit,wfit,sedout,gammaupmax,lambdapivot,lambdaout,result
;
; sedout=result
;
; ad hoc solution for wavelengths less than 60 micron
;
; find the values of MODEL M at the reference wavelengths at which the
; data is sampled
;
factor=sfit(0)/model(0)
findalpha1,lambdaout,12.0,38.0,factor,alpha
sedout=sedout*alpha
oplot,lambdaout,sedout,linest=1
endif ; ifudge
;
returnplot
;
free_lun, inunit
close,/all
;
c999:
;
return
end
;
pro hiidata,source,cat,lambda,flux,lum
;
; outputs:
; SOURCE: (string array of dimension NSOURCE giving name of source)
; CAT: (string array of dimension NSOURCE giving category of source
;       GHII UCHII)
; LAMBDA: (1D array of dimension NLAMBDA containing wavelengths in micron)
; FLUX (2D array of dimension NSOURCE, NLAMBDA cotaining flux densities in Jy
;       such that wavelengths for which no flux is available are set to REF)
; LUM (1D array of dimension NSOURCE:
;      IR luminosity in units of solar luminosities) 
;
ref=alog10(10.0^(-6.))
lambda=[8.,12.,21.,25.,60.,70.,100.,160.,240.,350.,850.,1200.]
nlambda=n_elements(lambda)
nsource=1000L
source=strarr(nsource)
cat=strarr(nsource)
flux=fltarr(nsource,nlambda)
lum=fltarr(nsource)
;
isource=0L
source(isource)='0.361'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,1.82,2.56,2.63,3.34,ref,3.09,ref,ref,ref,ref,ref]
lum(isource)=10.0^5.49
;
isource=1L
source(isource)='0.394'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.14,2.93,2.87,3.62,ref,3.34,ref,ref,ref,ref,ref]
lum(isource)=10.0^5.79
;
isource=2L
source(isource)='0.489'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,1.31,2.11,2.23,3.43,ref,3.78,ref,ref,ref,ref,ref]
lum(isource)=10.0^5.60
;
isource=3L
source(isource)='0.518'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.82,3.55,3.58,4.63,ref,4.72,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.76
;
isource=4L
source(isource)='0.572'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.17,2.82,2.96,3.66,ref,3.80,ref,ref,ref,ref,ref]
lum(isource)=10.0^5.84
;
isource=5L
source(isource)='1.149'
cat(isource)='ghii'
; Hill et al. 2005 SEST 1.2mm map
flux(isource,*)=10.0^[ref,2.55,3.01,3.10,3.94,ref,4.30,ref,ref,ref,ref,1.02531]
lum(isource)=10.0^6.16
;
isource=6L
source(isource)='2.303'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,1.60,2.28,2.61,3.44,ref,3.53,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.06
;
isource=7L
source(isource)='3.270'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,1.63,2.31,2.78,3.78,ref,3.95,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.39
;
isource=8L
source(isource)='4.412'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,1.65,2.28,2.75,3.67,ref,3.96,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.34
;
isource=9L
source(isource)='5.973'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,3.05,3.69,4.11,4.62,ref,4.55,ref,ref,ref,ref,ref]
lum(isource)=10.0^5.84
;
isource=10L
source(isource)='8.137'
cat(isource)='ghii'
; thompson et al. 2006 scuba
flux(isource,*)=10.0^[ref,2.26,2.80,2.93,3.79,ref,3.95,ref,ref,2.502,1.55,ref]
lum(isource)=10.0^6.40
;
isource=11L
source(isource)='10.159'
cat(isource)='ghii'
; hill et al 2005 SEST 1.2mm map
flux(isource,*)=10.0^[ref,3.23,3.75,3.77,4.52,ref,4.62,ref,ref,ref,ref,1.382]
lum(isource)=10.0^6.22
;
isource=12L
source(isource)='10.315'
cat(isource)='ghii'
; thompson et al. 2006 scuba pointing
; Hill et al 2004 SEST map 1.2mm
flux(isource,*)=10.0^[ref,2.75,3.37,3.35,4.06,ref,4.17,ref,ref,2.443,1.509,1.356]
lum(isource)=10.0^6.82
;
isource=13L
source(isource)='15.032' ; M17
cat(isource)='ghii'
; Hill et al 2005 SEST map 1.2mm 80Jy
; Gordon & Jewell NRAO 12-m 1.3mm 110 Jy => 152 Jy at 1.2mm
; take average of 116 Jy
flux(isource,*)=10.0^[ref,4.35,4.95,4.92,5.33,ref,5.35,ref,ref,ref,ref,2.065]
lum(isource)=10.0^6.61
;
isource=14L
source(isource)='20.733'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.28,2.92,3.25,4.15,ref,4.30,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.61
;
isource=15L
source(isource)='25.382'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,3.00,3.63,3.64,4.32,ref,4.34,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.84
;
isource=16L
source(isource)='29.944'
cat(isource)='ghii'
; Hill et al 2005 SEST map 1.2mm == 23.3 Jy
; Mooney et al 1995 1.3mm MRT == 20.0Jy which is the same as 27.5 at 1.2 mm
; we take the average = 25.4Jy
flux(isource,*)=10.0^[ref,2.79,3.39,3.33,4.18,ref,4.28,ref,ref,ref,ref,1.405]
lum(isource)=10.0^6.15
;
isource=17L
source(isource)='30.776'
cat(isource)='ghii'
; W43
; Hill et al 2005 SEST map 1.2mm 50.9 Jy
; Mooney et al. 1995 SEST and MRT 66.1 Jy at 1300 == 91.0 at 1.2mm
; we take the average of 80.0Jy
flux(isource,*)=10.0^[ref,3.19,3.71,3.82,4.61,ref,4.72,ref,ref,ref,ref,1.9031]
lum(isource)=10.0^6.56
;
isource=18L
source(isource)='32.797'
cat(isource)='ghii'
; thompson et al. 2006 scuba. But, assuming a nu^4 law only gets 
;                             67 percent of SEST flux
; Hill et al 2005 SEST map 1.2mm
flux(isource,*)=10.0^[ref,1.66,2.28,2.42,3.54,ref,3.72,ref,ref,2.497,1.377,0.964]
lum(isource)=10.0^6.09
;
isource=19L
source(isource)='43.169'
; W49A
cat(isource)='ghii'
; Sievers et al 1991 SEST 1270 = 138 => 1200 = 152
; Ward-Thompson & Robson 1990 185 at 1100 => 131
; Buckley & Ward-Thompson 1996 155 at 1100 316 at 800 2883 at 450
;                            equivalent to 109, 248 and 2883 at 1200 850 450 
; Gordon 1987 70 Jy at 1295 => 95 at 1200
; average for 1200 is 129 (Sievers & Ward-Thompson x 2)
; Ward-Thompson & Robson 1990 403 at 800 => 316 at 850
; average of W-T * 2 at 850 is 282
; Ward-Thompson & Robson 1990 5500 at 350
; Gordon 1987 4700 at 400 => 8000 at 350 and 2934 at 450
; average gordon and  W-T at 450 is 
flux(isource,*)=10.0^[ref,2.90,3.47,3.59,4.55,ref,4.74,ref,ref,3.464,2.45,2.11]
lum(isource)=10.0^7.04
;
isource=20L
source(isource)='48.596'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.22,2.63,2.85,3.93,ref,4.11,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.23
;
isource=21L
source(isource)='48.930'
cat(isource)='ghii'
; Mooney et al. 1995 SEST and MRT map of all W51 sources in 21um field
;                 10.1 +5.8 Jy 1300 == 21.9Jy at 1.2mm
flux(isource,*)=10.0^[ref,2.91,3.54,3.66,4.23,ref,4.18,ref,ref,ref,ref,1.34]
lum(isource)=10.0^6.10
;
isource=22L
source(isource)='49.384'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.72,3.27,3.34,4.23,ref,4.43,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.08
;
isource=23L
source(isource)='49.486'
; W51A
cat(isource)='ghii'
; Hill et al 2005 SEST map 1.2mm 126Jy - we ignore this
; Sievers et al 1991 SEST 1270 = 300 => 1200 = 376
; Gordon 1987 8840 at 400 => 15080 at 350
flux(isource,*)=10.0^[ref,3.29,3.94,4.16,4.84,ref,4.97,ref,ref,4.178,ref,2.575]
lum(isource)=10.0^6.68
;
isource=24L
source(isource)='70.300'
cat(isource)='ghii'
; thompson et al. 2006 scuba
flux(isource,*)=10.0^[ref,2.78,3.31,3.35,4.21,ref,4.27,ref,ref,2.408,1.575,ref]
lum(isource)=10.0^6.43
;
isource=25L
source(isource)='79.293'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.16,2.73,3.05,3.90,ref,3.97,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.05
;
isource=26L
source(isource)='133.720'
; W3
cat(isource)='ghii'
; Moore et al. 2007 SCUBA scanmap 850um = 339Jy
; Gordon 1987 1.3mm 46Jy == 63.4Jy at 1200um NRAO 12-m
flux(isource,*)=10.0^[ref,3.25,3.98,4.06,4.70,ref,4.84,ref,ref,ref,2.530,1.802]
lum(isource)=10.0^6.35
;
isource=27L
source(isource)='274.013'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.81,3.42,3.52,4.30,ref,4.35,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.27
;
isource=28L
source(isource)='282.023'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.73,3.30,3.43,4.17,ref,4.18,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.06
;
isource=29L
source(isource)='284.301'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,3.75,4.37,4.43,5.09,ref,5.08,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.82
;
isource=30L
source(isource)='287.379'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,3.38,4.12,4.49,5.00,ref,5.03,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.14
;
isource=31L
source(isource)='289.066'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,1.97,2.42,2.83,3.82,ref,4.04,ref,ref,ref,ref,ref]
lum(isource)=10.0^5.93
;
isource=32L
source(isource)='291.284'
cat(isource)='ghii'
; Hill et al. 2005 SEST 1.2mm map
flux(isource,*)=10.0^[ref,3.28,4.05,4.03,4.64,ref,4.80,ref,ref,ref,ref,1.79934]
lum(isource)=10.0^6.07
;
isource=33L
source(isource)='291.610'
cat(isource)='ghii'
; Hill et al. 2005 SEST 1.2mm map
flux(isource,*)=10.0^[ref,3.88,4.55,4.59,5.06,ref,5.07,ref,ref,ref,ref,1.44091]
lum(isource)=10.0^6.07
;
isource=34L
source(isource)='298.227'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,3.00,3.69,3.74,4.30,ref,4.28,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.75
;
isource=35L
source(isource)='298.862'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.73,3.39,3.63,4.42,ref,4.46,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.78
;
isource=36L
source(isource)='305.359'
cat(isource)='ghii'
; Hill et al. 2005 1.2mm map
flux(isource,*)=10.0^[ref,2.64,3.17,3.52,4.35,ref,4.40,ref,ref,ref,ref,1.512]
lum(isource)=10.0^5.75
;
isource=37L
source(isource)='319.158'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,1.85,2.39,2.71,3.67,ref,3.86,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.12
;
isource=38L
source(isource)='319.392'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.26,2.93,3.00,3.87,ref,3.96,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.34
;
isource=39L
source(isource)='320.327'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.27,2.72,2.83,3.75,ref,3.95,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.31
;
isource=40L
source(isource)='327.304'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,3.05,3.65,3.82,4.60,ref,4.79,ref,ref,ref,ref,ref]
lum(isource)=10.0^5.92
;
isource=41L
source(isource)='327.993'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.12,2.48,2.64,3.59,ref,3.71,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.04
;
isource=42L
source(isource)='330.868'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.12,2.62,3.08,4.12,ref,4.30,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.49
;
isource=43L
source(isource)='331.324'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.15,2.71,2.90,3.67,ref,3.81,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.09
;
isource=44L
source(isource)='331.354'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.12,2.66,2.86,3.69,ref,3.84,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.11
;
isource=45L
source(isource)='331.529'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.92,3.44,3.52,4.48,ref,4.61,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.89
;
isource=46L
source(isource)='333.122'
cat(isource)='ghii'
; 1.2mm Bains et al 2006 SEST clump fit
flux(isource,*)=10.0^[ref,2.78,3.38,3.64,4.52,ref,4.66,ref,ref,ref,ref,2.162]
lum(isource)=10.0^5.93
;
isource=47L
source(isource)='333.293'
cat(isource)='ghii'
; 1.2mm Bains et al 2006 SEST clump fit
flux(isource,*)=10.0^[ref,2.86,3.50,3.60,4.39,ref,4.54,ref,ref,ref,ref,1.887]
lum(isource)=10.0^5.84
;
isource=48L
source(isource)='333.610'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,3.65,4.10,4.12,4.55,ref,4.60,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.03
;
isource=49L
source(isource)='338.398'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,1.93,2.65,3.45,4.23,ref,4.27,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.73
;
isource=50L
source(isource)='338.400'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,1.28,1.87,2.13,3.30,ref,3.48,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.00
;
isource=51L
source(isource)='345.555'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,1.96,2.54,2.77,3.74,ref,3.90,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.43
;
isource=52L
source(isource)='345.645'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,1.56,2.17,2.49,3.64,ref,3.75,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.30
;
isource=53L
source(isource)='347.611'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.16,2.91,3.38,4.17,ref,4.28,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.11
;
isource=54L
source(isource)='351.467'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.07,2.43,2.69,3.65,ref,3.81,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.25
;
isource=55L
source(isource)='359.429'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.16,2.93,3.30,4.53,ref,4.66,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.63
;
isource=56L
source(isource)='30dor'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,2.50,3.21,3.32,3.98,ref,3.93,ref,ref,ref,ref,ref]
lum(isource)=10.0^7.74
;
isource=57L
source(isource)='ngc346'
cat(isource)='ghii'
flux(isource,*)=10.0^[ref,0.35,1.23,1.51,2.43,ref,2.48,ref,ref,ref,ref,ref]
lum(isource)=10.0^6.29
;
nsource=isource+1L
source=source(0:nsource-1)
cat=cat(0:nsource-1)
lum=lum(0:nsource-1)
fluxout=fltarr(nsource,nlambda)
for ilambda=0,nlambda-1 do fluxout(*,ilambda)=flux(0:nsource-1,ilambda)
flux=fluxout
;
return
end
;
pro prepareplot
;
set_plot,'ps'
;
;print,'setting line and character widths for ps file ...'
!p.charsize=1.3 ; annotation size for axes
!p.charthick=2.5 ; thickness for annotation 
!p.thick=1.0  ; thickness for plotted points
!x.thick=2.5  ; thickness of x axis
!y.thick=2.5  ; thickness of y axis
;!p.symsize=10.0   ; plotting symbol size
;
;xsize=18.0
;ysize=10.0
;
device,filename='dopita1.ps';,xs=xsize,ys=ysize 
device,/close
;
return
end
;
pro returnplot
;

set_plot,'x'
;print,'written plot to file dopita1.ps'
;print,'re-setting line and character widths for screen plot...'
!p.charsize=1.0 ; annotation size for axes
!p.charthick=1.0 ; thickness for annotation 
!p.thick=1.0  ; thickness for plotted points
!x.thick=1.0  ; thickness of x axis
!y.thick=1.0  ; thickness of y axis
;
return
end
pro pivotfit,lambda,s,weight,model,lambdapivot,npoly,lambdafit,sfit
;
; Finds an NPOLY order polynomial fit SFIT sampled at LAMBDAFIT,
; to input data S sampled at LAMBDA and with weight WEIGHT.
; such that (i)   SFIT is an input function MODEL multiplied by
;                 an internally used polynomial function ALPHA
;           [ii]  ALPHA takes unity value at a pivot wavelength
;                 LAMBDAPIVOT, which should also be one of the sampling
;                 points in LAMBDA
;
; For a cubic fit set NPOLY to 4. The constraint about the pivot
; means that there are then only NPOLY-1 polynomial coefficients
; to be determined, since the offset in the fitting function ALPHA
; has to be unity.
;
w=weight
;w(*)=1.0
x=lambda-lambdapivot
if(min(abs(x))/lambdapivot gt 0.01)then begin
 print,'PIVOTFIT: aborting: pivot wavelength not sampled in input array'
 print,'lambdapivot = ', lambdapivot
 goto, c999
endif ; abs(min(test))/lambdapivot
;
ndata=n_elements(s)
if(ndata le npoly)then begin
 print,'PIVOTFIT: aborting: requested polynomial of too high order'
 print,'npoly = ', npoly
 goto, c999
endif ; 
;
; find the values of MODEL M at the reference wavelengths at which the
; data is sampled
;
m=fltarr(ndata)
for idata=0,ndata-1 do begin
 test=abs(lambdafit-lambda(idata))
 ii=where(test eq min(test),nii)
 m(idata)=model[ii]
endfor ; idata
;
a=fltarr(npoly,npoly) ; matrix elements
b=fltarr(npoly)
for j=1,npoly-1L do begin
 xj=float(j)
 b(j)=total(w*m*(s-m)*(x^xj))
 for k=1,npoly-1L do begin
  xk=float(k)
  a(j,k)=total(w*m*m*(x^xk)*(x^xj))
 endfor ; j
endfor ; k
;print,'new'
;
; solve equations for alpha
;
asolve=a(1:npoly-1,1:npoly-1)
bsolve=b(1:npoly-1)
asolveinvert=invert(asolve)
c=fltarr(npoly)
c(0)=1.0d0
alpha=lambdafit
alpha(*)=0.0d0
for k=1,npoly-1L do begin
 for j=1,npoly-1L do begin
  c(k)=c(k)+asolveinvert(j-1,k-1)*bsolve(j-1)
 endfor ; j
endfor ; k
for k=0,npoly-1L do begin
 xk=float(k)
 alpha=alpha+c(k)*((lambdafit-lambdapivot)^xk)
endfor ; k
;
sfit=model*alpha
;
;save,lambda,s,w,m,model,lambdapivot,npoly,lambdafit,sfit,alpha,file='pivotfit.xdr'
;
c999:
;
return
end
;
pro testpivotfit
;
;
;
;
lambda=[12.0000,21.0000,25.0000,60.0000,100.000,350.000,850.000,1200.00]
nlambda=n_elements(lambda)
;s=planck(10000.0*lambda,25.0)
;s=2.0*(lambda/100.)+1.0*(lambda/100.)^2.
lambdafit=12.0+float(findgen(1200))
modeltarget=5.0*planck(10000.0*lambdafit,25.0)
s=fltarr(nlambda)
for i=0,nlambda-1 do begin
 test=abs(lambdafit-lambda(i))
 ii=where(test eq min(test),nii)
 ii=ii(0)
 s(i)=modeltarget[ii]
endfor ; i
model=modeltarget*(lambdafit/100.)
;model=3.0*(lambdafit/100.)
npoly=4L ; if use pivotfit
gammaupmax=100.0 ; if use exppivotfit1
lambdapivot=100.0
;test=abs(lambdafit-lambdapivot)
;ii=where(test eq min(test),nii)
;ii=ii(0)
;model=model-model[ii]+s(4) ; ensure data and model match at 100 micron
weight=s
weight(*)=1.0
;
llambda=alog10(lambda)
ls=alog10(s)
lmodel=alog10(model)
llambdapivot=alog10(lambdapivot)
llambdafit=alog10(lambdafit)
; polynomial fit
;pivotfit,lambda,s,weight,model,lambdapivot,npoly,lambdafit,sfit
;pivotfit,llambda,ls,weight,lmodel,llambdapivot,npoly,llambdafit,lsfit 
; exppivotfit
lgammaupmax=alog10(gammaupmax)
exppivotfit1,lambda,s,weight,model,gammaupmax,lambdapivot,lambdafit,sfit
exppivotfit1,llambda,ls,weight,lmodel,lgammaupmax,llambdapivot,llambdafit,lsfit
; for log
sfit=10.0d0^(lsfit)
;
window,0
plot_oo,lambdafit,model
oplot,lambdafit,sfit,linest=1
oplot,lambda,s,psym=1
;
return
end
;
pro exppivotfit1,lambda,s,weight,model,gammaupmax,lambdapivot,lambdafit,sfit
;
; Finds a fit SFIT sampled at LAMBDAFIT,
; to input data S sampled at LAMBDA and with weight WEIGHT.
; This fit is specific to the particular dataset used to
; constrain the Dopita model. 
; SFIT is an input function MODEL multiplied by
; an internally used function ALPHA where
; alpha = 1.0 - betaup * [1 - exp(-abs(deltalambda / gammaup)) ]
; for positive deltalambda and
; alpha = 1.0 + betadown * [1 - exp(+abs(deltalambda / gammadown)) ]
; for negative deltalambda where deltalambda = lambda - lambdapivot and
; LAMBDAPIVOT is a pivot wavelength which should also be one of
; the sampling points in LAMBDA. The 4 internal parameters are fixed
; as follows:
; - 1.0 - betaup (the asymptotic value at high positive deltalambda)
;   is fixed to be the value s/model at the highest measured deltalambda
; - 1.0 + betadown (the asymptotic value at high positive deltalambda)
;   is fixed to be the value s/model at the most negative measured deltalambda
; - gammaup and gammadown are found by direct search by minimising
;   sum [ weight (s-model*alpha) ], where the constraint guaranteeing
;   continuity of the gradient at deltalambda=0, 
;   betaup/gammaup = betadown/gammadown, is enforced, and the search range
;   for gammaup is limited to the range 0 to the input parameter GAMMAUPMAX
;   
w=weight
x=lambda-lambdapivot
if(min(abs(x))/lambdapivot gt 0.01)then begin
 print,'EXPPIVOTFIT1: aborting: pivot wavelength not sampled in input array'
 print,'lambdapivot = ', lambdapivot
 goto, c999
endif ; abs(min(test))/lambdapivot
;
; find the values of MODEL M at the reference wavelengths at which the
; data is sampled
;
ndata=n_elements(lambda)
m=fltarr(ndata)
for idata=0,ndata-1 do begin
 test=abs(lambdafit-lambda(idata))
 ii=where(test eq min(test),nii)
 m(idata)=model[ii]
endfor ; idata
;
;save,m,s,lambda,file='dump.xdr'
imin=where(lambda eq min(lambda))
;betadown=1.0-(m(imin)/s(imin))
betadown=(s(imin)/m(imin))-1.0
if(betadown lt 0.0)then begin
 print,'EXPPIVOTFIT1: aborting - betadown lt 0.0; betadown = ',betadown
goto,c999
endif ; betadown
imax=where(lambda eq max(lambda))
;betaup=(m(imax)/s(imax))-1.0
betaup1=1.0-(s(imax)/m(imax))
betaup2=1.0-(s(imax-1)/m(imax-1)) ; include nehboring points
betaup=0.5*(betaup1+betaup2)
if(betaup lt 0.0 or betaup gt 1.0)then begin
print,'EXPPIVOTFIT1: aborting - betaup lt 0.0; betaup = ',betaup
goto, c999
endif ; betaup
;
print,'EXPPIVOTFIT1: betadown,betaup = ',betadown,betaup
;
; find gammaup
;
tol=gammaupmax/100.0 ; tolerance and minimum search value for gammaup in micron
ng=5L ; must be odd
gstep=(gammaupmax-tol)/float(ng-1L)
gtrial=tol+(gammaupmax-tol)*findgen(ng)/float(ng-1L)
findalpha,lambda,lambdapivot,betaup,betadown,tol,alphatrial
chi=w*(s-m*alphatrial)
chi=total(chi*chi)
chimin=chi
gmin=tol
alpha=alphatrial
itrial=-1L
itrialmin=itrial+1
c200:
itrial=itrial+1
g=gtrial(itrial)
findalpha,lambda,lambdapivot,betaup,betadown,g,alphatrial
chi=w*(s-m*alphatrial)
chi=total(chi*chi)
print,'itrial,g,chi = ',itrial,g,chi
if(chi lt chimin)then begin
 chimin=chi
 itrialmin=itrial
 gmin=g
 alpha=alphatrial
 print,'*** chimin,gmin,itrialmin = ',chimin,gmin,itrialmin
endif ; chi lt chimin
if(itrial ge ng-1L)then begin
 gstep=gstep/1.5
 if(gstep lt tol)then goto,c300
 g1=gmin-gstep*float(ng/2)
 g2=gmin+gstep*float(ng/2)
 if(g2 gt gammaupmax)then begin
  g2=gammaupmax
  g1=g2-gstep*float(ng-1L)
 endif ; g2 gt gammaupmax
 gtrial=g1+(g2-g1)*findgen(ng)/float(ng-1L)
 itrial=0L
endif ; itrial ge ng
goto, c200
c300:
gammaup=gmin
print,'EXPPIVOTFIT1: gammaup = ',gammaup
;
findalpha,lambdafit,lambdapivot,betaup,betadown,gammaup,alpha
sfit=model*alpha
;
;save,lambda,s,w,m,model,lambdapivot,gammaupmax,lambdafit,sfit,alpha,$
betaup,betadown,file='exppivotfit1.xdr'
;
c999:
;
return
end
;
pro findalpha,lambda,lambdapivot,betaup,betadown,gammaup,alpha
;
; calculates a function ALPHA for each value of an input array LAMBDA
; where alpha = 1.0 - betaup * [1 - exp(-abs(deltalambda / gammaup)) ]
; for positive deltalambda 
; and alpha = 1.0 + betadown * [1 - exp(-abs(deltalambda / gammadown)) ]
; for negative deltalambda 
; where deltalambda = lambda - lambdapivot and
; LAMBDAPIVOT is a pivot wavelength at which alpha should be unity.
; The 4 internal parameters are fixed as follows:
; - +1.0 - betaup (the asymptotic value at high positive deltalambda)
;   is fixed to be the value s/model at the highest measured deltalambda
; - +1.0 + betadown (the asymptotic value at high positive deltalambda)
;   is fixed to be the value s/model at the most negative measured deltalambda
; - gammaup and gammadown are symptotes for high positive and low negative
;   values of deltalambda 
;
alpha=lambda
alpha[*]=0.0d0
deltalambda=lambda-lambdapivot
gammadown=gammaup*(betadown/betaup)
ipos=where(deltalambda ge 0.0,nipos)
if(nipos gt 0)then begin
 for i=0,nipos-1L do begin
  ii=ipos[i]
  test=-abs(deltalambda[ii]/gammaup)
  alpha[ii]=1.0d0-betaup*(1.0-exp(test))
 endfor ; ipos
; test=-abs(deltalambda(ipos)/gammaup)
; ii=where(test lt -87.0,nii)
; if(nii gt 0)then test[ii]=-87.0
; alpha(ipos)=1.0d0-betaup*(1.0-exp(test))
endif ; nipos gt 0
ineg=where(deltalambda lt 0.0,nineg)
if(nineg gt 0)then begin
 for i=0,nineg-1L do begin
  ii=ineg[i]
  test=-abs(deltalambda[ii]/gammadown)
  alpha[ii]=1.0d0+betadown*(1.0-exp(test))
 endfor ; ineg
; test=-abs(deltalambda(ineg)/gammadown)
; ii=where(test lt -87.0,nii)
; if(nii gt 0)then test[ii]=-87.0
; alpha(ineg)=1.0d0+betadown*(1.0-exp(test))
endif ; nineg gt 0
;
return
end
pro findalpha1,lambda,lambda1,lambda2,beta,alpha
;
; calculates a function ALPHA for each value of an input array LAMBDA where 
; alpha = beta+1 for lambda <= lambda1
;       = 1.0  for lambda >= lambda2
;       = 1.0+(beta-1.0)*sin[0.5*!pi * -(lambda-lambda2)/(lambda2-lambda1) ]
;         for  lambda1 < lambda < lambda2
alpha=lambda
alpha[*]=0.0d0
ii1=where(lambda le lambda1,nii1)
ii2=where(lambda gt lambda1 and lambda lt lambda2,nii2)
ii3=where(lambda ge lambda2,nii3)
if(nii1 gt 0)then alpha[ii1]=beta+1.0
if(nii3 gt 0)then alpha[ii3]=1.0d0
if(nii2 gt 0)then begin
dum1=(lambda2-lambda[ii2])/lambda2
dum2=(lambda2-lambda1)/lambda2
dum1=(20.0^dum1)-1.0
dum2=(20.0^dum2)-1.0
dum=dum1/dum2
;dum=(lambda2-lambda(ii2))/(lambda2-lambda1)
alpha[ii2]=$
1.0+0.5*beta*(1.0+sin(!pi*(dum-0.5)))
; alpha(ii2)=$
;1.0+0.5*beta*(1.0+sin(!pi*(((lambda2-lambda(ii2))/(lambda2-lambda1))-0.5)))
endif ; nii2
;
;print,'new'
return
end
