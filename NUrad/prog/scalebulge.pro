pro scalebulge,model,qyear,tau,scaabs,nsersic,bulge_disk

; subroutine of 'maps.pro'

; The purpose of this code is correctly scale the bulge according to only the stellar emission (tau = 0), this is what the raw maps output from urad are. This code opens these raw output maps for each component present and sums up the total of each.
; Then the total bulge is divided by the sum of the totals of every other component, and this value is used in the routine 'maps.pro'.
root = '../'
dir=root+'out/'

ss=''
name = root+'indata/geometry.in'
openr,unit,name,/get_lun
readf, unit, ss
readf, unit, tau1
readf, unit, ss
readf, unit, tau2
readf, unit, ss
readf, unit, hd
readf, unit, ss
readf, unit, zd
readf, unit, ss
readf, unit, hdin
readf, unit, ss
readf, unit, zdin
readf, unit, ss
readf, unit, hdsolar
readf, unit, ss
readf, unit, zdsolar
readf, unit, ss
readf, unit, hd1
readf, unit, ss
readf, unit, zd1
readf, unit, ss
readf, unit, hd1in
readf, unit, ss
readf, unit, zd1in
readf, unit, ss
readf, unit, hd1solar
readf, unit, ss
readf, unit, zd1solar
readf, unit, ss
readf, unit, h_bdisc
readf, unit, ss
readf, unit, h_vdisc
readf, unit, ss
readf, unit, h_idisc
readf, unit, ss
readf, unit, h_jdisc
readf, unit, ss
readf, unit, h_kdisc
readf, unit, ss
readf, unit, h_ir34disc
readf, unit, ss
readf, unit, h_ir45disc
readf, unit, ss
readf, unit, h_ir58disc
readf, unit, ss
readf, unit, zs
readf, unit, ss
readf, unit, hsin
readf, unit, ss
readf, unit, zsin
readf, unit, ss
readf, unit, hssolar
readf, unit, ss
readf, unit, zssolar
readf, unit, ss
readf, unit, hs1
readf, unit, ss
readf, unit, zs1
readf, unit, ss
readf, unit, hs1in
readf, unit, ss
readf, unit, zs1in
readf, unit, ss
readf, unit, hs1solar
readf, unit, ss
readf, unit, zs1solar
readf, unit, ss
readf, unit, rtrun
readf, unit, ss
readf, unit, sharp
readf, unit, ss
readf, unit, rtrund
readf, unit, ss
readf, unit, sharpd
readf, unit, ss
readf, unit, rtrun1
readf, unit, ss
readf, unit, sharp1
readf, unit, ss
readf, unit, rtrund1
readf, unit, ss
readf, unit, sharpd1
readf, unit, ss
readf, unit, reff
readf, unit, ss
readf, unit, ellipt
readf, unit, ss
readf, unit, nsersic
readf, unit, ss
readf, unit, xis0
readf, unit, ss
readf, unit, xis1
readf, unit, ss
readf, unit, xid0
readf, unit, ss
readf, unit, xid1
readf, unit, ss
readf, unit, idisc1
readf, unit, ss
readf, unit, idisc2
FOR skip = 0,190 DO BEGIN
        readf, unit, ss
ENDFOR
readf, unit, ss
readf, unit, hstin
readf, unit, ss
readf, unit, hdtin
readf, unit, ss
readf, unit, hs1tin
readf, unit, ss
readf, unit, hd1tin
readf, unit, ss
readf, unit, hs3tin
readf, unit, ss
readf, unit, hd3tin
readf, unit, ss
readf, unit, hs4tin
readf, unit, ss
readf, unit, hd4tin
readf, unit, ss
readf, unit, hs5tin
readf, unit, ss
readf, unit, hd5tin
readf, unit, ss
readf, unit, hs6tin
readf, unit, ss
readf, unit, hd6tin

free_lun, unit
hs = h_bdisc
shd = strcompress(string(fix(hd*1000.)),/remove_all)
szd = strcompress(string(fix(zd*1000.)),/remove_all)
shd1 = strcompress(string(fix(hd1*1000.)),/remove_all)
szd1 = strcompress(string(fix(zd1*1000.)),/remove_all)
shs = strcompress(string(fix(hs*1000.)),/remove_all)
szs = strcompress(string(fix(zs*1000.)),/remove_all)
shs1 = strcompress(string(fix(hs1*1000.)),/remove_all)
szs1 = strcompress(string(fix(zs1*1000.)),/remove_all)
sreff = strcompress(string(fix(reff*1000.)),/remove_all)
sellipt = strcompress(string(fix(ellipt*100.)),/remove_all)
stau = strcompress(string(fix(tau*10.)),/remove_all)
snsersic = strcompress(string(fix(nsersic)),/remove_all)

filename_param = root+'indata/gal_param.in'
name_param = filename_param
openr,unit,name_param,/get_lun
readf, unit, ss
readf, unit, ss
readf, unit, ss
readf, unit, ss
readf, unit, ss
readf, unit, ss
readf, unit, ss
readf, unit, ss
readf, unit, ss
readf, unit, dist_in
readf, unit, ss
readf, unit, inclination
readf, unit, ss
readf, unit, ss
readf, unit, pix_size
readf, unit, ss
readf, unit, nx_b
readf, unit, ss
readf, unit, ny_b
readf, unit, ss
readf, unit, nx_n
readf, unit, ss
readf, unit, ny_n
readf, unit, ss
readf, unit, nx_i
readf, unit, ss
readf, unit, ny_i
readf, unit, ss
readf, unit, nx_m
readf, unit, ss
readf, unit, ny_m
readf, unit, ss
readf, unit, nx_o
readf, unit, ss
readf, unit, ny_o
free_lun, unit

sinclination=strcompress(string(fix(inclination)),/remove_all)


filter = ['uv36','b','v','i','j','k','ir36','ir45','ir58']
dim = n_elements(filter)
bulge_disk = DBLARR(dim)

for ii = 0, dim-1L do begin
	fname_b=dir+'map_mb_'+filter(ii)+'_'+model+'_q'+qyear+'_i'+sinclination+'_t'+stau+'_n'+snsersic+'_'+scaabs+'.dat'
	fname_d=dir+'map_md_'+filter(ii)+'_'+model+'_q'+qyear+'_i'+sinclination+'_t'+stau+'_'+scaabs+'.dat'
	fname_di=dir+'map_mdi_'+filter(ii)+'_'+model+'_q'+qyear+'_i'+sinclination+'_t'+stau+'_'+scaabs+'.dat'
	fname_do=dir+'map_mdo_'+filter(ii)+'_'+model+'_q'+qyear+'_i'+sinclination+'_t'+stau+'_'+scaabs+'.dat'
	fname_td=dir+'map_mtd_'+filter(ii)+'_'+model+'_q'+qyear+'_i'+sinclination+'_t'+stau+'_'+scaabs+'.dat'
	fname_tdi=dir+'map_mtdi_'+filter(ii)+'_'+model+'_q'+qyear+'_i'+sinclination+'_t'+stau+'_'+scaabs+'.dat'
	fname_tdo=dir+'map_mtdo_'+filter(ii)+'_'+model+'_q'+qyear+'_i'+sinclination+'_t'+stau+'_'+scaabs+'.dat'
	
	openr,unit1,fname_b,/get_lun
	openr,unit2,fname_d,/get_lun
	openr,unit3,fname_di,/get_lun
	openr,unit4,fname_do,/get_lun
	openr,unit5,fname_td,/get_lun
	openr,unit6,fname_tdi,/get_lun
	openr,unit7,fname_tdo,/get_lun
	
	zb=dblarr(nx_b,ny_b)
	zbx = dblarr(ny_b)
	for i = 0, nx_b-1L do begin
	    readf,unit1,zbx
	    zb(i,*) = zbx(*)
	endfor
	free_lun,unit1
	print, total(zb)
	
	zd=dblarr(nx_m,ny_m)
	zdi=dblarr(nx_i,ny_i)
	zdo=dblarr(nx_o,ny_o)
	ztd=dblarr(nx_m,ny_m)
	ztdi=dblarr(nx_i,ny_i)
	ztdo=dblarr(nx_o,ny_o)
	
	FOR i = 0, nx_m-1L DO BEGIN
		readf,unit2,zdx
                zd(i,*) = zdx(*)
		readf,unit5,ztdx
                ztd(i,*) = ztdx(*)
	ENDFOR
	FOR i = 0, nx_i-1L DO BEGIN
                readf,unit3,zdix
                zdi(i,*) = zdix(*)
                readf,unit6,ztdix
                ztdi(i,*) = ztdix(*)
        ENDFOR
	FOR i = 0, nx_o-1L DO BEGIN
                readf,unit4,zdox
                zdo(i,*) = zdox(*)
                readf,unit7,ztdox
                ztdo(i,*) = ztdox(*)
        ENDFOR
	
	free_lun,unit2
	free_lun,unit3
	free_lun,unit4
	free_lun,unit5
	free_lun,unit6
	free_lun,unit7
	;
	bulge_disk(ii)=total(zb)/(total(zd)+total(ztd)+total(zdi)+total(ztdi)+total(zdo)+total(ztdo))
	print, bulge_disk
;
endfor ;

end
