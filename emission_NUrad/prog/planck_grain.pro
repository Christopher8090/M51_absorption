pro planck_grain,comp,qcomp_year,aa,T,Q_pl_grain

common scaling

;Planck averaged emissivities
; liefert planckgemitteltes Emissionsvermoegen aller Koernchen, von denen
; die Emissionsvermoegen Q(lambda,a) berechnet wurden (Draine '93)

Staub_dir=rootdir+'/Draine/'
name ='planck_'+comp+qcomp_year

openr,unit,Staub_dir+name,/get_lun

dim = n_elements(aa)
;print, 'dim', dim

w=''
for i=0,6 do begin
	readf,unit,w
	;print,w
endfor
readf,unit,w
;print, w
dim_temp = long(strmid(w,0,2))
;print,dim_temp
readf,unit,w
dim_size = long(strmid(w,0,2))
;print,dim_size
a_pl = dblarr(dim_size)
for i = 0,dim_size-1 do begin
	readf,unit,dum
	a_pl(i) = dum
	;print,a_pl(i)
endfor

a1=0.
b1=0.
c1=0.
d1=0.
n1=0
T=fltarr(dim_temp)
Q_grain=fltarr(dim_size,dim_temp)

for n=0,dim_size-1 do begin
	for k=0,3 do begin
		readf,unit,w
		;print,w
	endfor
	for m=0,dim_temp-1 do begin
		readf,unit,a1,b1,c1,d1
		;print, a1,b1,c1,d1
		T(m)=a1
		if b1 eq 0 then Q_grain(n,m)=c1 else Q_grain(n,m)=b1
	endfor
endfor

free_lun, unit


; interpolation to the grain sizes from Draine
Q_pl_grain = DBLARR(dim,dim_temp)

FOR i=0,dim_temp-1 DO BEGIN
	Q_pl_grain(*,i) = interpol(Q_grain(*,i)*a_pl,a_pl,aa)
ENDFOR
end
