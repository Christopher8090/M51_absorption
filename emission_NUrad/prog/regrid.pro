pro regrid,lumdouble,dr_req,dz_req,dr_out,dz_out,$
lumdouble_reg,rrr_reg,zzz_reg,nr,nz
;
; converts irregularly spatially sampled array LUMDOUBLE
; into a regularly spatially sampled arrray LUMDOUBLE_REG
;
COMMON gridinfo, lambda, rrr, zzz, dim

;input:
; DR_REQ,DZ_REQ requested maximum sampling interval of
;               of emissivity in r,z
; LUMDOUBLE (3D array in r,z,lambda): emissivity 
; 
; RRR (2D array in r,z): value of r position
;   it is assumed that values are monotonically increasing
;   along rows and each row is identical
; ZZZ (2D array in r,z): value of z position
;   it is assumed that values are monotonically increasing
;   along columns and each column is identical
; LAMBDA (1D array) wavelengths
; DIM number of wavelengths - dimension of LUMDOUBLE(0,0,*) and LAMBDA

;output:
; LUMDOUBLE_REG (3D array in r,z,lambda): emissivity 
; RRR_REG (2D array in r,z): value of r position
; ZZZ_REG (2D array in r,z): value of z position
; NR,NZ: dimensions of RRR_REG, ZZZ_REG, LUMDOUBLE_REG(*,*,0)
; DR_OUT,DZ_OUT: sampling interval in r,z of LUMDOUBLE_REG
;     These will be the largest possible numbers smaller
;     than DR_REQ,DZ_REQ for which the range in r,z of RRR_REG,ZZZ_REG
;     exactly matches the range in r,z of LUMDOUBLE
;
;restore,/verb,'gridinfo.xdr'
;
nrrr=n_elements(rrr(*,0)) ; dimensions of input arrays
nzzz=n_elements(zzz(0,*))
min_rrr=min(rrr)
max_rrr=max(rrr)
min_zzz=min(zzz)
max_zzz=max(zzz)
if(min_rrr ne 0. or min_zzz ne 0.)then begin
 print,'RGRID: aborting - input minimum values of rrr zzz not zero'
 goto,c999
endif

; find dimensions and spatial sampling of output arrays
nr=1L+floor(max_rrr/dr_req)
if(float(nr)*dr_req ne max_rrr)then nr=nr+1L
dr_out=max_rrr/float(nr-1L)
nz=1L+floor(max_zzz/dz_req)
if(float(nz)*dz_req ne max_zzz)then nz=nz+1L
dz_out=max_zzz/float(nz-1L)

lumdouble_reg=dblarr(nr,nz,dim)


; peform spline interpolation to each column,
; followed by spline interpolation to each created row

work1=dblarr(nzzz)
work2=dblarr(nrrr,nz)
z_in=work1
z_in(*)=zzz(0,*)
z_out=dz_out*findgen(nz)

work3=dblarr(nrrr)
r_in=work3
r_in(*)=rrr(*,0)
r_out=dr_out*findgen(nr)

for ilambda=0,dim-1L do begin
; print,'ilambda=',ilambda
 for irrr=0,nrrr-1L do begin ; input array counter
  work1(*)=lumdouble(irrr,*,ilambda)
  result=interpol(work1,z_in,z_out,/spline) ; result has dimension nz
  work2(irrr,*)=result(*)
 endfor ; irrr
 ;if ilambda eq 700 then save,work2,file='regrid_work.xdr'
 for iz=0,nz-1L do begin
  work3(*)=work2(*,iz)
  result=interpol(work3,r_in,r_out,/spline) ; result has dimension nr
  lumdouble_reg(*,iz,ilambda)=result(*)
 endfor ; iz
endfor ; ilambda


;rrr_reg=dblarr(nr,nz)
;zzz_reg=rrr_reg
;for iz=0,nz-1 do begin
; rrr_reg(*,iz)=r_out(*)
;endfor ; iz
;for ir=0,nr-1 do begin
; zzz_reg(ir,*)=z_out(*)
;endfor ; ir
;---------------------------------------------------

lumdouble_reg_swap=dblarr(nr,nz,dim)


; peform spline interpolation to each column,
; followed by spline interpolation to each created row

work1=dblarr(nrrr)
work2=dblarr(nr,nzzz)

work3=dblarr(nzzz)

for ilambda=0,dim-1L do begin
; print,'ilambda=',ilambda
 for izzz=0,nzzz-1L do begin ; input array counter
  work1(*)=lumdouble(*,izzz,ilambda)
  result=interpol(work1,r_in,r_out,/spline) ; result has dimension nr
  work2(*,izzz)=result(*)
 endfor ; izzz
 if ilambda eq 700 then save,work2,file='regrid_work.xdr'
 for ir=0,nr-1L do begin
  work3(*)=work2(ir,*)
  result=interpol(work3,z_in,z_out,/spline) ; result has dimension nz
  lumdouble_reg_swap(ir,*,ilambda)=result(*)
 endfor ; ir
endfor ; ilambda
;------------------------------------------------------------------

lumdouble_reg = (lumdouble_reg + lumdouble_reg_swap)/2. 

rrr_reg=dblarr(nr,nz)
zzz_reg=rrr_reg
for iz=0,nz-1 do begin
 rrr_reg(*,iz)=r_out(*)
endfor ; iz
for ir=0,nr-1 do begin
 zzz_reg(ir,*)=z_out(*)
endfor ; ir
;
c999:
;
;
;save,lumdouble_reg_swap,rrr_reg,zzz_reg,nr,nz,dr_out,dz_out,file='regrid.xdr'
;
return
end
