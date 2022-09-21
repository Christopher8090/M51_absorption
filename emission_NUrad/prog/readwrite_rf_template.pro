pro readwrite_rf_template,filter,ibulge,idisk,idisk3,idisk5,model,qyear,tau,sfr,sfr4,sfr6,sfr7,newold,old3,old5,$
bd,nsersic,scaabs,factor_bulge,factor_disk,factor_disk1,factor_disk3,factor_disk4,factor_disk5,factor_disk6,$
factor_disk7,dir1,dir2,swdisk3,swdisk4,swdisk5,swdisk6,swdisk7

close,/all
common param, shd, szd, shd1, szd1, shs, szs, shs1, szs1, sreff, sellipt
common scaling

if idisk eq 'no' and ibulge eq 'yes' then begin
	print,'wrong combination of morphologies for the given wavelength'
	print, 'idisk,ibulge,filter', idisk,ibulge,filter
	print, 'aborting:'
	goto, mark11
endif

;define some strings for file names
stau=strcompress(string(round(tau*10)),/remove_all)
suv=strcompress(string(round(sfr*100)),/remove_all)
sbd=strcompress(string(round(bd*100)),/remove_all)
sold=strcompress(string(round(newold*100)),/remove_all)
snsersic=strcompress(string(round(nsersic)),/remove_all)

thin_param = '_'+model+'_q'+qyear+'_t'+stau+'_hd'+shd+'_zd'+szd+$
	     '_hd1_'+shd1+'_zd1_'+szd1+'_hs1_'+shs1+'_zs1_'+szs1+'_'+scaabs+'.dat'
thick_param = '_'+model+'_q'+qyear+'_t'+stau+'_hd'+shd+'_zd'+szd+$
	      '_hd1_'+shd1+'_zd1_'+szd1+'_hs'+shs+'_zs'+szs+'_'+scaabs+'.dat'
bulge_param = '_'+model+'_q'+qyear+'_t'+stau+'_hd'+shd+'_zd'+szd+$
	      '_hd1_'+shd1+'_zd1_'+szd1+'_reff'+sreff+'_ell'+sellipt+'_n'+snsersic+'_'+scaabs+'.dat'

combine = 0.
if combine eq 1 then begin
	wd_td = '_wd01_q06_t254_hd6000_zd160_hd1_4300_zd1_90_hs1_4300_zs1_90_sca.dat'
	wd_d = '_wd01_q06_t254_hd6000_zd160_hd1_4300_zd1_90_hs3800_zs190_sca.dat'
	wd_b = '_lmc1_q06_t177_hd6300_zd160_hd1_4300_zd1_90_reff350_ell145_n4_sca.dat'
	lmc_td = '_lmc1_q06_t177_hd6300_zd160_hd1_4300_zd1_90_hs1_5000_zs1_90_sca.dat'
	lmc_d = '_lmc1_q06_t177_hd6300_zd160_hd1_4300_zd1_90_hs4400_zs190_sca.dat'
	lmc_b = '_wd01_q06_t254_hd6000_zd160_hd1_4300_zd1_90_reff350_ell145_n4_sca.dat'
endif
components = []
;open file with radiation fields for the second disk
filenamed1 = 'unit_mtd_'+filter+thin_param
if combine eq 1 then filenamed1 = 'unit_mtd_'+filter+wd_td
named1 = dir1+filenamed1
openr,unit10,named1,/get_lun
components = [components, 'td'] 
;inner tdisk
if swdisk4 eq 'yes' then begin
	filenamed4 = 'unit_mtdi_'+filter+thin_param
	if combine eq 1 then filenamed4 = 'unit_mtdi_'+filter+lmc_td
	named4 = dir1+filenamed4
	openr, unit14, named4, /get_lun
	components = [components, 'tdi']
endif
;outer tdisk
if swdisk6 eq 'yes' then begin
	filenamed6 = 'unit_mtdo_'+filter+thin_param
	if combine eq 1 then filenamed6 = 'unit_mtdo_'+filter+wd_td
	named6 = dir1+filenamed6
	openr,unit16,named6,/get_lun
	components = [components, 'tdo']
endif
;nuclear tdisk
if swdisk7 eq 'yes' then begin
	filenamed7 = 'unit_mntd_'+filter+thin_param
	if combine eq 1 then filenamed7 = 'unit_mntd_'+filter+wd_td
	named7 = dir1+filenamed7
	openr,unit17,named7,/get_lun
	components = [components, 'ntd']
endif
;open file with radiation fields for the first disk
if idisk eq 'yes' then begin
	filenamed = 'unit_md_'+filter+thick_param
	if combine eq 1 then filenamed = 'unit_md_'+filter+wd_d
	named = dir1+filenamed
	openr,unit0,named,/get_lun
	components = [components, 'd']
endif
; inner disk
if swdisk3 eq 'yes' then begin
if idisk3 eq 'yes' then begin
	filenamed3 = 'unit_mdi_'+filter+thick_param
	if combine eq 1 then filenamed3 = 'unit_mdi_'+filter+lmc_d
	named3 = dir1+filenamed3
	openr,unit03,named3,/get_lun
	components = [components, 'di']
endif
endif
;outer disk
if swdisk5 eq 'yes' then begin
if idisk5 eq 'yes' then begin
	filenamed5 = 'unit_mdo_'+filter+thick_param
	if combine eq 1 then filenamed5 = 'unit_mdo_'+filter+wd_d
	named5 = dir1+filenamed5
	openr,unit05,named5,/get_lun
	components = [components, 'do']
endif
endif
if idisk eq 'no' then begin
	filenamed = 'unit_md_b'+thick_param
	if combine eq 1 then filenamed = 'unit_md_b'+lmc_d
	named = dir1+filenamed
	openr,unit0,named,/get_lun
	components = [components, 'b']
endif
;inner disk
if swdisk3 eq 'yes' then begin
if idisk3 eq 'no' then begin
	filenamed3 = 'unit_mdi_b'+thick_param
	if combine eq 1 then filenamed3 = 'unit_mdi_b'+lmc_d
	named3 = dir1+filenamed3
	openr,unit03,named3,/get_lun
	components = [components, 'di']
endif
endif
;outer disk
if swdisk5 eq 'yes' then begin
if idisk5 eq 'no' then begin
	filenamed5 = 'unit_mdo_b'+thick_param
	if combine eq 1 then filenamed5 = 'unit_mdo_b'+wd_d
	named5 = dir1+filenamed5
	openr,unit05,named5,/get_lun
	components = [components, 'do']
endif
endif

;open file with radiation fields for the bulge
if ibulge eq 'yes' then begin
	filenameb = 'unit_mb_'+filter+bulge_param
	if combine eq 1 then filenameb = 'unit_mb_'+filter+lmc_b
	nameb = dir1+filenameb
	openr,unit1,nameb,/get_lun
	components = [components, 'b']
endif

if ibulge eq 'no' then begin
	filenameb = 'unit_mb_b'+bulge_param
	if combine eq 1 then filenameb = 'unit_mb_b'+lmc_b
	nameb = dir1+filenameb
	openr,unit1,nameb,/get_lun
	components = [components, 'b']
endif
if filter eq 'uv09' then print, 'Components found: ', components
;open file for output radiation fields
filename = 'u_'+filter+'_'+model+'_q'+qyear+$
               '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
               '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
               '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
               '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+'.dat'
if combine eq 1 then filename = 'u_'+filter+'_combined_s'+suv+'_o'+sold+'.dat'
name = dir2+filename
openw,unit2,name,/get_lun
;print, 'write ', name

dim_positionsr = 0L
dim_positionsz = 0L
idisk1 = 0L
idisk2 = 0L
ss = ''

;print the input parameters in the header of the output files
printf, unit2, 'face-on central tau B'
printf, unit2, tau
printf, unit2, 'sfr'
printf, unit2, sfr
printf, unit2, 'old'
printf, unit2, newold
printf, unit2, 'bd'
printf, unit2, bd
;read the header of the radiation fields files
readf, unit10, ss
printf, unit2, ss
readf, unit10, tau1
printf, unit2, tau1
readf, unit10, ss
printf, unit2, ss
readf, unit10, hd
printf, unit2, hd
readf, unit10, ss
printf, unit2, ss
readf, unit10, zd
printf, unit2, zd
readf, unit10, ss
printf, unit2, ss
readf, unit10, hdin
printf, unit2, hdin
readf, unit10, ss
printf, unit2, ss
readf, unit10, zdin
printf, unit2, zdin
readf, unit10, ss
printf, unit2, ss
readf, unit10, hdsolar
printf, unit2, hdsolar
readf, unit10, ss
printf, unit2, ss
readf, unit10, zdsolar
printf, unit2, zdsolar
readf, unit10, ss
printf, unit2, ss
readf, unit10, tau2
printf, unit2, tau2
readf, unit10, ss
printf, unit2, ss
readf, unit10, hd1
printf, unit2, hd1
readf, unit10, ss
printf, unit2, ss
readf, unit10, zd1
printf, unit2, zd1
readf, unit10, ss
printf, unit2, ss
readf, unit10, hd1in
printf, unit2, hd1in
readf, unit10, ss
printf, unit2, ss
readf, unit10, zd1in
printf, unit2, zd1in
readf, unit10, ss
printf, unit2, ss
readf, unit10, hd1solar
printf, unit2, hd1solar
readf, unit10, ss
printf, unit2, ss
readf, unit10, zd1solar
printf, unit2, zd1solar

for i = 0,27 do begin
     readf,unit0,ss
endfor
if swdisk3 eq 'yes' then begin
for i = 0,57 do begin;added by JJT 6/10/17
     readf,unit03,ss;added by JJT 6/10/17
endfor;added by JJT 6/10/17
endif
if swdisk5 eq 'yes' then begin
for i = 0,57 do begin;added by JJT 6/10/17
     readf,unit05,ss;added by JJT 6/10/17
endfor;added by JJT 6/10/17
endif
for i = 0,27 do begin
     readf,unit1,ss
endfor
if swdisk4 eq 'yes' then begin
for i = 0,57 do begin;added by JJT 6/10/17
     readf,unit14,ss;added by JJT 6/10/17
endfor;added by JJT 6/10/17
endif
if swdisk6 eq 'yes' then begin
for i = 0,57 do begin;added by JJT 6/10/17
     readf,unit16,ss;added by JJT 6/10/17
endfor;added by JJT 6/10/17
endif
if swdisk7 eq 'yes' then begin
for i = 0,57 do begin;added by JJT 
     readf,unit17,ss;added by JJT
endfor;added by JJT
endif 
readf, unit0, ss
printf, unit2, ss
readf, unit0, pass
printf, unit2, pass
readf, unit0, ss
printf, unit2, ss
readf, unit0, zs
printf,unit2, zs
readf, unit0, ss
printf, unit2, ss
readf, unit0, hsin
printf, unit2, hsin
readf, unit0, ss
printf, unit2, ss
readf, unit0, zsin
printf,unit2, zsin
readf, unit0, ss
printf, unit2, ss
readf, unit0, hssolar
printf,unit2, hssolar
readf, unit0, ss
printf, unit2, ss
readf, unit0, zssolar
printf,unit2, zssolar
readf, unit0, ss
printf, unit2, ss
readf, unit0, rtrun
printf, unit2, rtrun

readf, unit10, ss
printf, unit2, ss
readf, unit10, hs1
printf, unit2, hs1
readf, unit10, ss
printf, unit2, ss
readf, unit10, zs1
printf, unit2, zs1
readf, unit10, ss
printf, unit2, ss
readf, unit10, hs1in
printf, unit2, hs1in
readf, unit10, ss
printf, unit2, ss
readf, unit10, zs1in
printf,unit2, zs1in
readf, unit10, ss
printf, unit2, ss
readf, unit10, hs1solar
printf,unit2, hs1solar
readf, unit10, ss
printf, unit2, ss
readf, unit10, zs1solar
printf,unit2, zs1solar
readf, unit10, ss
printf, unit2, ss
readf, unit10, rtrun1
printf, unit2, rtrun1

readf, unit1, ss
readf, unit1, reff
readf, unit1, ss
readf, unit1, ellipt
readf, unit1, ss
readf, unit1, nsersic

readf, unit10, ss
printf, unit2, ss
readf, unit10, xis0
printf, unit2, xis0

readf, unit0, ss
printf, unit2, ss
readf, unit0, xis1
printf, unit2, xis1

readf, unit10, ss
printf, unit2, ss
readf, unit10, xid0
printf, unit2, xid0
readf, unit10, ss
printf, unit2, ss
readf, unit10, xid1
printf, unit2, xid1
readf, unit10, ss
printf, unit2, ss
readf, unit10, idisk1
printf, unit2, idisk1
readf, unit10, ss
printf, unit2, ss
readf, unit10, idisk2
printf, unit2, idisk2
readf, unit10, ss
printf, unit2, ss
readf, unit10, dim_positionsr
printf, unit2, dim_positionsr
readf, unit10, ss
printf, unit2, ss
readf, unit10, dim_positionsz
printf, unit2, dim_positionsz
readf, unit10, ss
printf, unit2, ss
readf, unit10, ss
printf,unit2, '    r(pc)      z(pc)   urad(erg/pc3/A)        eabs(erg/s/pc3/A)'


for i = 0,13 do begin
     readf,unit0,ss
endfor

for i = 0,13 do begin
     readf,unit1,ss
endfor

dim_positions = dim_positionsr * dim_positionsz

ud = 1.d
ud1 = 1.d
ub = 1.d
uu = 1.d
ed = 1.d
ed1 = 1.d
eb = 1.d
ee = 1.d

ud3 = 1.d
ud4 = 1.d
ed3 = 1.d
ed4 = 1.d

ud5 = 1.d
ud6 = 1.d
ed5 = 1.d
ed6 = 1.d

ud7 = 1.d
ed7 = 1.d

rad = dblarr(dim_positions)
height = dblarr(dim_positions)
rfed = dblarr(dim_positions)
energy_abs = dblarr(dim_positions)

for k = 0L, dim_positions-1 do begin
	; thin main disk, td
	readf, unit10, r, z, ud1, ed1
	ud1 = ud1 * factor_disk1
	ed1 = ed1 * factor_disk1
	uu = ud1
	ee = ed1
	frac = 0.
	if idisk eq 'yes' then begin
		; thick main disk, d
		readf,unit0,r,z,ud,ed
		ud = ud * factor_disk
		ed = ed * factor_disk
		uu = uu + ud
		ee = ee + ed
		if ibulge eq 'yes' then begin
			; bulge, b
			readf,unit1,r,z,ub,eb
			ub = ub * factor_bulge
			eb = eb * factor_bulge
			uu = uu + ub
			ee = ee + eb
		endif
	endif
	; JJT 9/10/17
	;no extra bulge taken into account
	;inner component
	if swdisk4 eq 'yes' then begin
		; thin inner disk, tdi
		readf,unit14,r,z,ud4,ed4
		ud4 = ud4 * factor_disk4
		ed4 = ed4 * factor_disk4
		uu = uu + ud4
		ee = ee + ed4
		frac = 0.
		if idisk3 eq 'yes' then begin
			; thick inner disk, di
			readf,unit03,r,z,ud3,ed3
			ud3 = ud3 * factor_disk3
			ed3 = ed3 * factor_disk3
			uu = uu + ud3
			ee = ee + ed3
		endif
	endif
	;outer component
	if swdisk6 eq 'yes' then begin
		; thin outer disk, tdo
		readf,unit16,r,z,ud6,ed6
		ud6 = ud6 * factor_disk6
		ed6 = ed6 * factor_disk6
		uu = uu + ud6
		ee = ee + ed6
		frac = 0.
		if idisk5 eq 'yes' then begin
			; thick outer disk , do
			readf,unit05,r,z,ud5,ed5
			ud5 = ud5 * factor_disk5
			ed5 = ed5 * factor_disk5
			uu = uu + ud5
			ee = ee + ed5
		endif
	endif
	;nuclear component
	if swdisk7 eq 'yes' then begin
		; nuclear disk, ntd
		readf,unit17,r,z,ud7,ed7
		ud7 = ud7 * factor_disk7
		ed7 = ed7 * factor_disk7
		uu = uu + ud7
		ee = ee + ed7
		frac = 0.
	endif
    
	printf, unit2, r,z,uu,ee
	rad(k) = r
	height(k) = z 
	rfed(k) = uu
	energy_abs(k) = ee
endfor
savename = urad_dir+'/saves/model/'+model+'_rfed_'+filter+'.save'
save, rad, height, rfed, energy_abs, filename = savename
;print, 'saved: '+savename
free_lun, unit10
free_lun, unit0
free_lun, unit1
free_lun, unit2
if (swdisk3 eq 'yes') then free_lun, unit03
if (swdisk4 eq 'yes') then free_lun, unit14
if (swdisk5 eq 'yes') then free_lun, unit05
if (swdisk6 eq 'yes') then free_lun, unit16
if (swdisk7 eq 'yes') then free_lun, unit17

mark11:
end
