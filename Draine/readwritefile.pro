;the program read the tables for PAHs and put them in the same format as for Si
;and Gra
;one of the changes is that it removes the column with Q_ext
;after outp
uting the new files, go with the editor and delete the name of the
;column that was deleted (e.g. string-replace 'Q_ext' with ' ')

pro readwritefile,comp
;readwritefile,'PAHion01'

common dirdef,rootdir

Staub_dir=rootdir+'/emission1/Draine/'
name = Staub_dir+comp+'_old'
name1 = Staub_dir+comp
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
  a=0
  b=0.
  c=0.
  d=0.
  e=0.
  for j=0,dim_wave-1 do begin
    readf,unit,b,a,c,d,e
    printf,unit1,b,c,d,e
    lambda(j) = b
    Q_abs(i,j) = c
    Q_sc(i,j) = d
    Q_ph(i,j) = e
  endfor
endfor

free_lun, unit
free_lun, unit1
end
