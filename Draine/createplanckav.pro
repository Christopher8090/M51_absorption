pro createplanckav,comp,qcomp_year
;
; creates file containing planck averaged emissivities for
; composition COMP and model year QCOMP_YEAR
;
; eg createplanckav,'Gra','01'
; eg createplanckav,'Gra','06'
; eg createplanckav,'Si','01'
;
common dirdef,rootdir
;
tmaxemission=1500.0d0 ; maximum temperature for calculation of emission
;
Staub_dir=rootdir+'/emission1/Draine/'
inname=Staub_dir+comp+qcomp_year
;outname=Staub_dir+'planck_'+comp+qcomp_year
outname=Staub_dir+'planck_'+'test_'+comp+qcomp_year ; for debugging
;
nsizeout=0
if(comp+qcomp_year eq 'Gra01')then begin
 nsizeout=25
 sizeout=[0.001,0.002,0.003,0.005,0.006,0.007]
 sizeout=[sizeout,0.01,0.015,0.02,0.03,0.05,0.07]
 sizeout=[sizeout,0.1,0.15,0.2,0.3,0.5,0.7]
 sizeout=[sizeout,1.,1.5,2.,3.,5.,7.,10.]
 ntout=38L
 tmin=10.
 tmax=5.012E+04
 tout=alog10(tmin)+(alog10(tmax)-alog10(tmin))*$
   findgen(ntout)/(double(ntout-1L))
 tout=10.0d0^(tout)
endif
;
if(comp+qcomp_year eq 'Gra06')then begin
 nsizeout=25
 sizeout=[0.001,0.002,0.003,0.005,0.006,0.007]
 sizeout=[sizeout,0.01,0.015,0.02,0.03,0.05,0.07]
 sizeout=[sizeout,0.1,0.15,0.2,0.3,0.5,0.7]
 sizeout=[sizeout,1.,1.5,2.,3.,5.,7.,10.]
 ntout=38L
 tmin=10.
 tmax=5.012E+04
 tout=alog10(tmin)+(alog10(tmax)-alog10(tmin))*$
   findgen(ntout)/(double(ntout-1L))
 tout=10.0d0^(tout)
 endif
;
if(comp+qcomp_year eq 'Si01')then begin
 nsizeout=25
 sizeout=[0.001,0.002,0.003,0.005,0.006,0.007]
 sizeout=[sizeout,0.01,0.015,0.02,0.03,0.05,0.07]
 sizeout=[sizeout,0.1,0.15,0.2,0.3,0.5,0.7]
 sizeout=[sizeout,1.,1.5,2.,3.,5.,7.,10.]
 ntout=38L
 tmin=10.
 tmax=5.012E+04
 tout=alog10(tmin)+(alog10(tmax)-alog10(tmin))*$
   findgen(ntout)/(double(ntout-1L))
 tout=10.0d0^(tout)
endif
;
if(nsizeout eq 0)then begin
 print,'CREATEPLANCKAV: ',comp+qcomp_year,' not supported'
 goto, c999
endif ; nsizeout eq 0
openr,inunit,inname,/get_lun
print,'CREATEPLANCKAV: reading ',inname
q_grain_s,comp,qcomp_year,lambda,Q_abs,Q_sc,Q_ph,sizes
;help,sizes
;save,sizes,file='sizes.xdr'
nlambda=n_elements(lambda)
nsize=n_elements(q_abs(*,0))
print,'CREATEPLANCKAV: input nlambda, nsize  = ',nlambda,nsize
print,'CREATEPLANCKAV: input min, max lambda = ',min(lambda),max(lambda)
;
; find emissivities for input sizes, for each temperature TOUT
;
qe=fltarr(nsize,ntout) ; Planck-averaged Q_em for T_d=T_rad
qa=qe ; Planck-averaged Q_abs for T_d=25.0K
      ; at present these quanties are the same (except at 
      ; temperatures greater than tmaxemission) since we
      ; do not have the temperature dependent emissivities
q=fltarr(nlambda)
for isize=0,nsize-1L do begin ; input sizes
 for itout=0,ntout-1L do begin ; output temperatures
  trad=tout(itout)
  q(*)=q_abs(isize,*) ; for T_d=T_rad=tout(itout)
  planckav,q,lambda,trad,value
  qe(isize,itout)=value
  qa(isize,itout)=value
  if(trad gt tmaxemission)then qe(isize,itout)=0.0
;  q(*)=q_abs(isize,*) ; for T_d=25K
;  planckav,q,lambda,trad,tmaxemission,value
;  qa(isize,itout)=value
 endfor ; itemp
endfor ; isize
;
; interpolate emissivities for input sizes to output sizes and write output
;
openw,outunit,outname,/get_lun
printf,outunit,comp+qcomp_year
printf,outunit,'Planck-averaged Q_em, Q_abs calculated by CREATEPLANCKAV'
printf,outunit,' '
printf,outunit,' '
printf,outunit,'Emissivity <Qem> is computed only for T<',tmaxemission
printf,outunit,'Planck-averaged Q_em for T_d=T_rad'
printf,outunit,'Planck-averaged Q_abs and Q_pr for T_d= 25.00K'
printf,outunit,ntout,' temperature bins',format='(E9.3,a17)'
printf,outunit,nsizeout,' grain sizes:',format='(E9.3,a13)'
for isizeout=0,nsizeout-1L do begin
 printf,outunit,sizeout(isizeout),format='(E9.3)'
endfor ; isize
;
qeout=fltarr(nsizeout,ntout)
qaout=qeout
qedum=fltarr(nsize)
qeo=fltarr(nsizeout)
qadum=fltarr(nsize)
qao=fltarr(nsizeout)
dummy=0.0
nnn=n_elements(sizes)
sizesamp=findgen(nnn)
sizesampout=interpol(sizesamp,sizes,sizeout) ; sampling points of sizeout
                                             ; in units of input size sampling
 for itout=0,ntout-1L do begin ; output temperatures
  qedum(*)=qe(*,itout)
;  help,qeo,qedum,sizes,sizeout
;  qeo(*)=interpol(qedum,sizes,sizeout)
  qeo(*)=interpolate(qedum,sizesampout,cubic=-0.5)
  qeout(*,itout)=qeo(*)
  if(tout(itout) gt tmaxemission)then qeout(*,itout)=0.0
  qadum(*)=qa(*,itout)
;  qao(*)=interpol(qadum,sizes,sizeout)
  qao(*)=interpolate(qedum,sizesampout,cubic=-0.5)
  qaout(*,itout)=qao(*)
 endfor ; itout
;
for isizeout=0,nsizeout-1L do begin
 printf,outunit,' '
 printf,outunit,comp+qcomp_year
 printf,outunit,sizeout(isizeout),'= radius (micron)',format='(E9.3,a17)'
 printf,outunit,'  T_rad    <Qem>/a   <Qabs>/a  dummy  (micron-1)'
 for itout=0,ntout-1L do begin
  printf,outunit,tout(itout),qeout(isizeout,itout)/sizeout(isizeout),$
    qaout(isizeout,itout)/sizeout(isizeout),dummy,format='(E9.3,3E10.3)'
 endfor ; itout
endfor ; isizeout
print,'CREATEPLANCKAV: output written to: ',outname
;
c999:
;
free_lun,inunit
free_lun,outunit
;
return
end
;
;------------------------------------------------------------------------------
pro q_grain_s,comp,qcomp_year,lambda,Q_abs,Q_sc,Q_ph,sizes
;
; same as q_grain but also outputs the grain sizes (radii) SIZES in micron
;
common dirdef,rootdir
;print,rootdir
;
Staub_dir=rootdir+'/emission1/Draine/'
name = Staub_dir+comp+qcomp_year
;print,name
openr,unit,name,/get_lun

w=''
for i=0,3 do begin
  readf,unit,w
  ;print,w
endfor
dim_Draine = long(strmid(w,0,4))
readf,unit,w
dim_wave = long(strmid(w,0,4))
sizes=fltarr(dim_Draine)

lambda=fltarr(dim_wave)
Q_abs=fltarr(dim_Draine,dim_wave)
Q_sc=fltarr(dim_Draine,dim_wave)
Q_ph=fltarr(dim_Draine,dim_wave)

for i=0,dim_Draine-1 do begin
  for k=0,1 do begin
    readf,unit,w
    ;print,w
endfor
  size = strmid(w,0,9)
  sizes(i)=size
  readf,unit,w
  b=0.
  c=0.
  d=0.
  e=0.
  for j=0,dim_wave-1 do begin
    readf,unit,b,c,d,e
    lambda(j) = b
    Q_abs(i,j) = c
    Q_sc(i,j) = d
    Q_ph(i,j) = e
  endfor
endfor

free_lun, unit
end
;
pro planckav,q,lambda,t,x
save,q,lambda,t,file='planckav.xdr'
;
; Finds Planck-averaged emissivities 
; x = integral[q(lambda)*planck(t,lambda)] / integral[planck(t,lambda)]
;
; input parameters
; LAMBDA in micron
; T in kelvin
; Q dimensionless
;
;integral can be done over lambda using IDL planck function
;which is in units of power/area/wavelength
;
; 
;iswitch=1 ; to use input course sampling in lambda for integration
          ; with linear interpolation over wavelength bin
iswitch=0 ; to use finer sampling in lambda for integration
          ; with cubic interpolation between sample points
          ; and with an extrapolation to 10 times maximum input wavelength
;
nlambda=n_elements(lambda) ; dimensions of input wavelength array
;
; create alternative set of lambdas for integration
;
if(iswitch eq 0)then begin
 lmin=min(lambda)
; lmin=0.0912 ;
 lmaxin=max(lambda)
 lmax=10.0*lmaxin
 lambdaextend=[lmax,lambda]
 nl=10001L
 l=alog10(lmin)+findgen(nl)*(alog10(lmax)-alog10(lmin))/float(nl-1)
 l=10.0d0^(l) ; in micron
 lsamp_in=findgen(nlambda+1) ;sampling points of extended input wavelen array
 dlam=10000.*(l(1:nl-1)-l(0:nl-2)) ; in A
 lamav=0.5*(l(1:nl-1)+l(0:nl-2)) ; in micron
 lsamp=interpol(lsamp_in,lambdaextend,lamav) ; sampling points of lamav
                                   ; in units of extended input lambda sampling
 lamav=lamav*10000. ; in A
 qdum=q(0)*((lmaxin/lmax)^2.0d0)
 qextend=[qdum,q]
 qav=interpolate(qextend,lsamp,cubic=-0.5)
endif ; iswitch eq 0
;
; ... or use input lambda for integration
;
if(iswitch ne 0)then begin
 dlam=10000.*(lambda(1:nlambda-1)-lambda(0:nlambda-2)) ; in A
 lamav=0.5*10000.*(lambda(1:nlambda-1)+lambda(0:nlambda-2))
 qav=0.5*(q(1:nlambda-1)+q(0:nlambda-2))
endif ; iswitch ne 0
;
x=total(dlam*qav*planck(lamav,t))
x=x/total(dlam*planck(lamav,t))
;
return
end
