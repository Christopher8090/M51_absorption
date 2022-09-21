pro plot_emissivities,comp
;plot_emissivities,'Gra06'

common dirdef,rootdir

dir = rootdir+'/emission1/Draine/'
name1 = dir+comp
print,'read ', name1

openr,unit,name1,/get_lun

ss = strmid(comp,0,1)
print,ss
w=''
for i=0,3 do begin
  readf,unit,w
  ;printf,unit1,w
  ;print,w
endfor
dim_Draine = long(strmid(w,0,4))
;print, dim_Draine
readf,unit,w
;printf,unit1,w
dim_wave = long(strmid(w,0,4))
;print,dim_wave

lambda=fltarr(dim_wave)
Q_abs=fltarr(dim_Draine,dim_wave)
Q_sc=fltarr(dim_Draine,dim_wave)
Q_ph=fltarr(dim_Draine,dim_wave)

for i=0,dim_Draine-1 do begin
  for k=0,2 do begin
    readf,unit,w
    ;printf,unit1,w
    ;print,w
  endfor
  a=0
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

;print, Q_abs(0,*)
dirfig = '/emission1/figures/'
set_plot,'ps'
name2 = 'emissivities_'+comp+'.ps'
name2 = rootdir+dirfig+name2
print, 'plot ', name2
device,filename=name2
!p.thick=3
!x.thick=4
!y.thick=4
!p.charthick=3
!p.charsize=1.5
if ss eq 'S' then begin
   xmin = 0.003
   xmax = 1.
   ymin = 0.01
   ymax = 2.
endif
if ss eq 'G' then begin
   xmin = 0.003
   xmax = 1.
   ymin = 0.03
   ymax = 3.
endif
if ss eq 'P' then begin
   xmin = 0.005
   xmax = 1.
   ymin = 0.007
   ymax = 3.
endif
if ss eq 'S' or ss eq 'G' then begin
  plot_oo, lambda,Q_abs(10,*),$
         xrange=[xmin,xmax],yrange=[ymin,ymax],$
         xstyle=1,ystyle=1,linestyle=0,$
         xtitle='!4k!3 [!4l!3m]',ytitle='Q!Iabs!N' ;0.003 micron

   oplot,lambda,Q_abs(20,*) ;0.01 micron
   oplot,lambda,Q_abs(30,*) ;0.03 micron
   oplot,lambda,Q_abs(40,*) ; 0.1 micron
   oplot,lambda,Q_abs(50,*) ;0.3 micron
   oplot,lambda,Q_abs(60,*) ;1 micron
endif
if ss eq 'P' then begin
  plot_oo, lambda,Q_abs(0,*),$
         xrange=[xmin,xmax],yrange=[ymin,ymax],$
         xstyle=1,ystyle=1,linestyle=0,$
         xtitle='!4k!3 [!4l!3m]',ytitle='Q!Iabs!N' ;0.00035 micron

   oplot,lambda,Q_abs(9,*) ;0.001 micron
   oplot,lambda,Q_abs(23,*) ;0.005 micron
   ;oplot,lambda,Q_abs(24,*) ; 0.006 micron
   ;oplot,lambda,Q_abs(26,*) ;0.007 micron
   oplot,lambda,Q_abs(29,*) ;0.01 micron
endif
device,/close
set_plot,'ps'
name3 = 'emissivities_ir_'+comp+'.ps'
name3 = rootdir+dirfig+name3
print, 'plot ', name3
device,filename=name3
!p.thick=3
!x.thick=4
!y.thick=4
!p.charthick=3
!p.charsize=1.5
if ss eq 'S' then begin
   xmin = 1.
   xmax = 1000.
   ymin = 0.1
   ymax = 20.
endif
if ss eq 'G' then begin
   xmin = 1.
   xmax = 1000.
   ymin = 0.2
   ymax = 10.
endif
if ss eq 'P' then begin
   xmin = 1.
   xmax = 1000.
   ymin = 0.02
   ymax = 4.
endif
if ss eq 'S' or ss eq 'G' then begin
  plot_oo, lambda,lambda * Q_abs(20,*)/0.01,$
         xrange=[xmin,xmax],yrange=[ymin,ymax],$
         xstyle=1,ystyle=1,linestyle=0,$
         xtitle='!4k!3 [!4l!3m]',ytitle='!4k!3Q!Iabs!N/a'

   ;oplot,lambda,lambda * Q_abs(20,*)/0.01
   ;oplot,lambda,lambda * Q_abs(30,*)/0.03
   oplot,lambda,lambda * Q_abs(40,*)/0.1
   oplot,lambda,lambda * Q_abs(50,*)/0.3
   oplot,lambda,lambda * Q_abs(60,*)/1.0
endif
if ss eq 'P' then begin
  plot_oo, lambda,lambda * Q_abs(23,*)/0.00501,$
         xrange=[xmin,xmax],yrange=[ymin,ymax],$
         xstyle=1,ystyle=1,linestyle=0,$
         xtitle='!4k!3 [!4l!3m]',ytitle='!4k!3Q!Iabs!N/a'
   ;oplot,lambda,lambda * Q_abs(9,*)/0.00100
   ;oplot,lambda,lambda * Q_abs(0,*)/0.00035
   oplot,lambda,lambda * Q_abs(24,*)/0.00562
   oplot,lambda,lambda * Q_abs(26,*)/0.00708
   oplot,lambda,lambda * Q_abs(29,*)/0.01
endif
device,/close
free_lun, unit
end

