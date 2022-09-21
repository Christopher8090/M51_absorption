;the program read the tables for Si and produce a file for black bodies

pro readwritefile2,comp
;readwritefile2,'Si01'

common dirdef,rootdir

Staub_dir=rootdir+'/emission1/Draine/'
name = Staub_dir+comp
name1 = Staub_dir+'Sibb'
print,name
openr,unit,name,/get_lun
openw,unit1,name1,/get_lun

w=''
for i=0,3 do begin
  readf,unit,w
  printf,unit1,w
  ;print,w
endfor
dim_Draine = long(strmid(w,0,4))
readf,unit,w
printf,unit1,w
;print,w
dim_wave = long(strmid(w,0,4))
;print,dim_wave

lambda=fltarr(dim_wave)
Q_abs=fltarr(dim_Draine,dim_wave)
Q_sc=fltarr(dim_Draine,dim_wave)
Q_ph=fltarr(dim_Draine,dim_wave)

for i=0,dim_Draine-1 do begin
  for k=0,2 do begin
    readf,unit,w
    printf,unit1,w
    ;print,w
  endfor
  b=0.
  c=0.
  d=0.
  e=0.
  for j=0,dim_wave-1 do begin
    readf,unit,b,c,d,e
    lambda(j) = b
    Q_abs(i,j) = 1.d
    Q_sc(i,j) = 0.d
    Q_ph(i,j) = 0.d
    printf,unit1,lambda(j),Q_abs(i,j),Q_sc(i,j),Q_ph(i,j)
  endfor
endfor
free_lun, unit
free_lun, unit1
end
