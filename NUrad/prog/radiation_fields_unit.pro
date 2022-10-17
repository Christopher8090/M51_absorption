;this program scales the radiation fields corresponding to the unit
;SFR and old
pro radiation_fields_unit,model,qyear,morph,tau,scaabs,nsersic

;input parameters:
; model - the dust model used (e.g. 'wd01' or 'pop00' refers to the model of
;         Weingartner and Draine 2001 or Popescu et al. 2000)
; qyear - year of last updated the extinction emissivities (e.g. '06'
;         or'01')
; morph - the morphological component ('d', 'td', 'b', 'di', 'tdi','do', or 'tdo' )
; tau   - the central face-on opacity in the B band
;scaabs  - can be 'sca' or 'abs'
;nsersic - in the case of bulges the sersic index is required
;          (1,2,4,8)

common dirdef,rootdir
close,/all
rootdir='./'
dir = '../'
dir1 = '/unit/'
dir2 = '/indata/'

common param, dist, factor, zs, hs, zs1, hs1, hsin, hs1in, ztrun, rtrun, ztrun1, rtrun1, xis0, xis1, ll, ll1, S0, S0_opt, ms, $
                            zs3, hs3, zs4, hs4, hs3in, hs4in, rtrun3, rtrun4, xis3, xis4, $
                            zs5, hs5, zs6, hs6, hs5in, hs6in, rtrun5, rtrun6, xis5, xis6, $
			                   hstin, hdtin, hs1tin, hd1tin, hs3tin, hd3tin, hs4tin, hd4tin, $
			                   hs5tin, hd5tin, hs6tin, hd6tin, $
                            rtrund, rtrund1,rtrund3, rtrund4, rtrund5, rtrund6, $
                            zs7, hs7, hs7in, rtrun7, xis7, hs7tin
;define some strings for file names
stau=strcompress(string(round(tau*10)),/remove_all)

factor=3.086 * 1.d+16  ;to convert [pc] in [m]


;truncation of the first disk
ztrun=84000.  ;in [pc]

;truncation of the second disk
ztrun1=84000.  ;in [pc]

ss=''

filename_param='gal_param.in'
name_param = rootdir+dir+dir2+filename_param
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
free_lun, unit
dist = dist_in   ;the distance to M51 in Mpc

ms = 20.         ;standard central face-on SB

filter_opt=['uv36','b','v','i','j','k','ir36','ir45','ir58'];,'nir']
dim_opt = n_elements(filter_opt)
filter_uv = ['uv09','uv13','uv15','uv16','uv20','uv22','uv25','uv28','uv36']
dim_uv = n_elements(filter_uv)
filter = ['uv09','uv13','uv15','uv16','uv20','uv22','uv25','uv28','uv36','b','v','i','j','k','ir36','ir45','ir58'];,'nir']
dim = n_elements(filter)

lambda=[912.,1350.,1500.,1650.,2000.,2200.,2500.,2800.,3650.,4430.,5640.,8090.,12590.,22000.,36000.,45000.,58000.] * 1.d ;in AA

;the luminosity of the young stellar disk in W/Hz corresponding to
;SFR=1Msolar/yr, from PT11
ll1 = [0.344,0.905,0.844,0.863,0.908,0.926,0.843,0.910,1.842,2.271,3.837,5.734,0.931,0.728,0.141,0.141,0.141] * 1.d+21

;the luminosity of the old stellar disk in W/Hz corresponding to old=1, from PT11
ll = [4.771,4.771,9.382,19.54,72.20,64.97,12.58,12.58,12.58] * 1.d+21


;zero point to convert Jy in magnitude
S0 = [460.,249.,460.,554.,710.,728.,759.,954.,1889.,4063.0d0,3636.0d0,2416.0d0,1589.0d0,640.0d0,280.9d0,179.9d0,115.0d0];,640.0d0]

S0_opt=[1889.0d0,4063.0d0,3636.0d0,2416.0d0,1589.0d0,640.0d0,280.9d0,179.9d0,115.0d0];,640.0d0]

if (morph eq 'td' or morph eq 'tdi' or morph eq 'tdo' or morph eq 'ntd') then begin
   dim_filter = dim
   filter = filter
endif

if (morph eq 'd' or morph eq 'di' or morph eq 'do') then begin
   dim_filter = dim_opt
   filter = filter_opt
endif

if morph eq 'b' then begin
	dim_filter = dim_opt
	filter = filter_opt
	snsersic=strcompress(string(round(nsersic)),/remove_all)
	scalebulge, model, qyear, tau, scaabs, nsersic, bulge_disk ; added 21/04/21 CJI
	scaabs=scaabs
endif

hs_pass = dblarr(dim_filter)
ii = 0
ss = ''
idisk1 = 0L
idisk2 = 0L
nsersic=0L

filename = 'geometry.in'
name = rootdir+dir+dir2+filename
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
	readf, unit, h_bdisk
	readf, unit, ss
	readf, unit, h_vdisk
	readf, unit, ss
	readf, unit, h_idisk
	readf, unit, ss
	readf, unit, h_jdisk
	readf, unit, ss
	readf, unit, h_kdisk
	readf, unit, ss
	readf, unit, h_ir36disk
	readf, unit, ss
	readf, unit, h_ir45disk
	readf, unit, ss
	readf, unit, h_ir58disk
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
	readf, unit, idisk1
	readf, unit, ss
	readf, unit, idisk2
	readf, unit, ss
	readf, unit, ss
	readf, unit, tau3
	readf, unit, ss
	readf, unit, tau4
	readf, unit, ss
	readf, unit, hd3
	readf, unit, ss
	readf, unit, zd3
	readf, unit, ss
	readf, unit, hd3in
	readf, unit, ss
	readf, unit, zd3in
	readf, unit, ss
	readf, unit, hd3solar
	readf, unit, ss
	readf, unit, zd3solar
	readf, unit, ss
	readf, unit, hd4
	readf, unit, ss
	readf, unit, zd4
	readf, unit, ss
	readf, unit, hd4in
	readf, unit, ss
	readf, unit, zd4in
	readf, unit, ss
	readf, unit, hd4solar
	readf, unit, ss
	readf, unit, zd4solar
	readf, unit, ss
	readf, unit, h_bdisk3
	readf, unit, ss
	readf, unit, h_vdisk3
	readf, unit, ss
	readf, unit, h_idisk3
	readf, unit, ss
	readf, unit, h_jdisk3
	readf, unit, ss
	readf, unit, h_kdisk3
	readf, unit, ss
	readf, unit, h_ir36disk3
	readf, unit, ss
	readf, unit, h_ir45disk3
	readf, unit, ss
	readf, unit, h_ir58disk3
	readf, unit, ss
	readf, unit, zs3
	readf, unit, ss
	readf, unit, hs3in
	readf, unit, ss
	readf, unit, zs3in
	readf, unit, ss
	readf, unit, hs3solar
	readf, unit, ss
	readf, unit, zs3solar
	readf, unit, ss
	readf, unit, hs4
	readf, unit, ss
	readf, unit, zs4
	readf, unit, ss
	readf, unit, hs4in
	readf, unit, ss
	readf, unit, zs4in
	readf, unit, ss
	readf, unit, hs4solar
	readf, unit, ss
	readf, unit, zs4solar
	readf, unit, ss
	readf, unit, rtrun3
	readf, unit, ss
	readf, unit, sharp3
	readf, unit, ss
	readf, unit, rtrund3
	readf, unit, ss
	readf, unit, sharpd3
	readf, unit, ss
	readf, unit, rtrun4
	readf, unit, ss
	readf, unit, sharp4
	readf, unit, ss
	readf, unit, rtrund4
	readf, unit, ss
	readf, unit, sharpd4
	readf, unit, ss
	readf, unit, xis3
	readf, unit, ss
	readf, unit, xis4
	readf, unit, ss
	readf, unit, xid3
	readf, unit, ss
	readf, unit, xid4
	readf, unit, ss
	readf, unit, idisk3
	readf, unit, ss
	readf, unit, idisk4
	readf, unit, ss
	readf, unit, ss
	readf, unit, tau5
	readf, unit, ss
	readf, unit, tau6
	readf, unit, ss
	readf, unit, hd5
	readf, unit, ss
	readf, unit, zd5
	readf, unit, ss
	readf, unit, hd5in
	readf, unit, ss
	readf, unit, zd5in
	readf, unit, ss
	readf, unit, hd5solar
	readf, unit, ss
	readf, unit, zd5solar
	readf, unit, ss
	readf, unit, hd6
	readf, unit, ss
	readf, unit, zd6
	readf, unit, ss
	readf, unit, hd6in
	readf, unit, ss
	readf, unit, zd6in
	readf, unit, ss
	readf, unit, hd6solar
	readf, unit, ss
	readf, unit, zd6solar
	readf, unit, ss
	readf, unit, h_bdisk5
	readf, unit, ss
	readf, unit, h_vdisk5
	readf, unit, ss
	readf, unit, h_idisk5
	readf, unit, ss
	readf, unit, h_jdisk5
	readf, unit, ss
	readf, unit, h_kdisk5
	readf, unit, ss
	readf, unit, h_ir36disk5
	readf, unit, ss
	readf, unit, h_ir45disk5
	readf, unit, ss
	readf, unit, h_ir58disk5
	readf, unit, ss
	readf, unit, zs5
	readf, unit, ss
	readf, unit, hs5in
	readf, unit, ss
	readf, unit, zs5in
	readf, unit, ss
	readf, unit, hs5solar
	readf, unit, ss
	readf, unit, zs5solar
	readf, unit, ss
	readf, unit, hs6
	readf, unit, ss
	readf, unit, zs6
	readf, unit, ss
	readf, unit, hs6in
	readf, unit, ss
	readf, unit, zs6in
	readf, unit, ss
	readf, unit, hs6solar
	readf, unit, ss
	readf, unit, zs6solar
	readf, unit, ss
	readf, unit, rtrun5
	readf, unit, ss
	readf, unit, sharp5
	readf, unit, ss
	readf, unit, rtrund5
	readf, unit, ss
	readf, unit, sharpd5
	readf, unit, ss
	readf, unit, rtrun6
	readf, unit, ss
	readf, unit, sharp6
	readf, unit, ss
	readf, unit, rtrund6
	readf, unit, ss
	readf, unit, sharpd6
	readf, unit, ss
	readf, unit, xis5
	readf, unit, ss
	readf, unit, xis6
	readf, unit, ss
	readf, unit, xid5
	readf, unit, ss
	readf, unit, xid6
	readf, unit, ss
	readf, unit, idisk5
	readf, unit, ss
	readf, unit, idisk6
	readf, unit, ss
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
	readf, unit, ss
	readf, unit, ss
	readf, unit, hs7
	readf, unit, ss
	readf, unit, zs7
	readf, unit, ss
	readf, unit, hs7in
	readf, unit, ss
	readf, unit, zs7in
	readf, unit, ss
	readf, unit, hs7solar
	readf, unit, ss
	readf, unit, zs7solar
	readf, unit, ss
	readf, unit, rtrun7
	readf, unit, ss
	readf, unit, sharp7
	readf, unit, ss
	readf, unit, xis7
	readf, unit, ss
	readf, unit, idisk7
	readf, unit, ss
	readf, unit, hs7tin
free_lun, unit

hs = h_bdisk
shd=strcompress(string(round(hd*1000.)),/remove_all)
szd=strcompress(string(round(zd*1000.)),/remove_all)
shd1=strcompress(string(round(hd1*1000.)),/remove_all)
szd1=strcompress(string(round(zd1*1000.)),/remove_all)
shs=strcompress(string(round(hs*1000.)),/remove_all)
szs=strcompress(string(round(zs*1000.)),/remove_all)
shs1=strcompress(string(round(hs1*1000.)),/remove_all)
szs1=strcompress(string(round(zs1*1000.)),/remove_all)
sreff=strcompress(string(round(reff*1000.)),/remove_all)
sellipt=strcompress(string(round(ellipt*100.)),/remove_all)

hs=[h_bdisk,h_bdisk,h_vdisk,h_idisk,h_jdisk,h_kdisk,h_kdisk,h_ir36disk,h_ir45disk,h_ir58disk] * 1000.d
zs = zs * 1000.d
hs1 = hs1 * 1000.d
zs1 = zs1 * 1000.d
hsin = hsin * 1000.d
hs1in = hs1in * 1000.d
rtrun = rtrun * 1000.d
rtrun1 = rtrun1 * 1000.d
hstin = hstin * 1000.d
hdtin = hdtin * 1000.d
hs1tin = hs1tin * 1000.d
hd1tin = hd1tin * 1000.d

hs3=[h_bdisk3,h_bdisk3,h_vdisk3,h_idisk3,h_jdisk3,h_kdisk3,h_kdisk3,h_ir36disk3,h_ir45disk3,h_ir58disk3] * 1000.d
zs3 = zs3 * 1000.d
hs4 = hs4 * 1000.d
zs4 = zs4 * 1000.d
hs3in = hs3in * 1000.d
hs4in = hs4in * 1000.d
rtrun3 = rtrun3 * 1000.d
rtrun4 = rtrun4 * 1000.d
hs3tin = hs3tin * 1000.d
hd3tin = hd3tin * 1000.d
hs4tin = hs4tin * 1000.d
hd4tin = hd4tin * 1000.d

hs5=[h_bdisk5,h_bdisk5,h_vdisk5,h_idisk5,h_jdisk5,h_kdisk5,h_kdisk5,h_ir36disk5,h_ir45disk5,h_ir58disk5] * 1000.d
zs5 = zs5 * 1000.d
hs6 = hs6 * 1000.d
zs6 = zs6 * 1000.d
hs5in = hs5in * 1000.d
hs6in = hs6in * 1000.d
rtrun5 = rtrun5 * 1000.d
rtrun6 = rtrun6 * 1000.d
hs5tin = hs5tin * 1000.d
hd5tin = hd5tin * 1000.d
hs6tin = hs6tin * 1000.d
hd6tin = hd6tin * 1000.d

hs7 = hs7 * 1000.d
zs7 = zs7 * 1000.d
hs7in = hs7in * 1000.d
rtrun7 = rtrun7 * 1000.d
hs7tin = hs7tin * 1000.d

;call routine that calculations the scaling factors
conversion_factors_flat,factord_old, factord1_young,factord3_old, factord4_young,factord5_old, factord6_young, factord7_young

factorb_opt = factord_old
if (morph eq 'td') then begin
	ff = factord1_young
endif

if (morph eq 'd') then begin
	ff = factord_old
endif

if morph eq 'b' then begin
	ff = factord_old/bulge_disk
endif

if (morph eq 'tdi') then begin
	ff = factord4_young
endif

if (morph eq 'di') then begin
	ff = factord3_old
endif

if (morph eq 'tdo') then begin
	ff = factord6_young
endif

if (morph eq 'do') then begin
	ff = factord5_old
endif

if (morph eq 'ntd') then begin
	ff = factord7_young
endif

dim_positionsr = 0L
dim_positionsz = 0L

for ii = 0, dim_filter-1 do begin ;loop in wavelength
;open file with input radiation fields
if filter[ii] eq 'nir' then begin
	if ( morph eq 'td' or morph eq 'tdi' or morph eq 'tdo' or morph eq 'ntd') then  filename = $
			'u_m'+morph+'_k_'+model+'_q'+qyear+'_t'+stau+$
			'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
			'_hs1_'+shs1+'_zs1_'+szs1+$
			'_'+scaabs+'.dat'
	if (morph eq 'd' or morph eq 'di' or morph eq 'do') then  filename = $
			'u_m'+morph+'_k_'+model+'_q'+qyear+'_t'+stau+$
			'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
			'_hs'+shs+'_zs'+szs+$
			'_'+scaabs+'.dat'
	if morph eq 'b' then  filename = $
			'u_m'+morph+'_k_'+model+'_q'+qyear+'_t'+stau+$
			'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
			'_reff'+sreff+'_ell'+sellipt+'_'+'n'+snsersic+'_'+scaabs+'.dat'
endif else begin
	if ( morph eq 'td' or morph eq 'tdi' or morph eq 'tdo' or morph eq 'ntd') then  filename = $
			'u_m'+morph+'_'+filter[ii]+'_'+model+'_q'+qyear+$
			'_t'+stau+$
			'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
			'_hs1_'+shs1+'_zs1_'+szs1+$
			'_'+scaabs+'.dat'
	if (morph eq 'd' or morph eq 'di' or morph eq 'do') then  filename = $
			'u_m'+morph+'_'+filter[ii]+'_'+model+'_q'+qyear+$
			'_t'+stau+$
			'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
			'_hs'+shs+'_zs'+szs+$
			'_'+scaabs+'.dat'
	if morph eq 'b' then  filename = $
			'u_m'+morph+'_'+filter[ii]+'_'+model+'_q'+qyear+$
			'_t'+stau+$
			'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
			'_reff'+sreff+'_ell'+sellipt+'_'+'n'+snsersic+'_'+scaabs+'.dat'
endelse 
name = rootdir+dir+'/out/'+filename
print, 'read ', name
openr,unit,name,/get_lun

;open file with output radiation fields
if (morph eq 'td' or morph eq 'tdi' or morph eq 'tdo' or morph eq 'ntd') then begin
	filename1 = 'unit_m'+morph+'_'+filter[ii]+'_'+model+'_q'+qyear+$
		    '_t'+stau+$
		    '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
		    '_hs1_'+shs1+'_zs1_'+szs1+$
		    '_'+scaabs+'.dat'
	name1 = rootdir+dir+dir1+filename1
	openw,unit1,name1,/get_lun
	print, 'write ', name1
	;read the header of the radiation fields files
	readf, unit, ss
	printf, unit1, ss
	readf, unit, tau1
	printf, unit1, tau1
	readf, unit, ss
	printf, unit1, ss
	readf, unit, hd
	printf, unit1, hd
	readf, unit, ss
	printf, unit1, ss
	readf, unit, zd
	printf, unit1, zd
	readf, unit, ss
	printf, unit1, ss
	readf, unit, hdin
	printf, unit1, hdin
	readf, unit, ss
	printf, unit1, ss
	readf, unit, zdin
	printf, unit1, zdin
	readf, unit, ss
	printf, unit1, ss
	readf, unit, hdsolar
	printf, unit1, hdsolar
	readf, unit, ss
	printf, unit1, ss
	readf, unit, zdsolar
	printf, unit1, zdsolar
	readf, unit, ss
	printf, unit1, ss
	readf, unit, tau2
	printf, unit1, tau2
	readf, unit, ss
	printf, unit1, ss
	readf, unit, hd1
	printf, unit1, hd1
	readf, unit, ss
	printf, unit1, ss
	readf, unit, zd1
	printf, unit1, zd1
	readf, unit, ss
	printf, unit1, ss
	readf, unit, hd1in
	printf, unit1, hd1in
	readf, unit, ss
	printf, unit1, ss
	readf, unit, zd1in
	printf, unit1, zd1in
	readf, unit, ss
	printf, unit1, ss
	readf, unit, hd1solar
	printf, unit1, hd1solar
	readf, unit, ss
	printf, unit1, ss
	readf, unit, zd1solar
	printf, unit1, zd1solar
	readf, unit, ss
	printf, unit1, ss
	readf, unit, hs1
	printf, unit1, hs1
	readf, unit, ss
	printf, unit1, ss
	readf, unit, zs1
	printf, unit1, zs1
	readf, unit, ss
	printf, unit1, ss
	readf, unit, hs1in
	printf, unit1, hs1in
	readf, unit, ss
	printf, unit1, ss
	readf, unit, zs1in
	printf, unit1, zs1in
	readf, unit, ss
	printf, unit1, ss
	readf, unit, hs1solar
	printf, unit1, hs1solar
	readf, unit, ss
	printf, unit1, ss
	readf, unit, zs1solar
	printf, unit1, zs1solar
	readf, unit, ss
	printf, unit1, ss
	readf, unit, rtrun1
	printf, unit1, rtrun1
	readf, unit, ss
	printf, unit1, ss
	readf, unit, xis1
	printf, unit1, xis1
	readf, unit, ss
	printf, unit1, ss
	readf, unit, xid0
	printf, unit1, xid0
	readf, unit, ss
	printf, unit1, ss
	readf, unit, xid1
	printf, unit1, xid1
	readf, unit, ss
	printf, unit1, ss
	readf, unit, idisk1
	printf, unit1, idisk1
	readf, unit, ss
	printf, unit1, ss
	readf, unit, idisk2
	printf, unit1, idisk2
endif
if (morph eq 'd' or morph eq 'di' or morph eq 'do') then begin
	filename1 = 'unit_m'+morph+'_'+filter[ii]+'_'+model+'_q'+qyear+'_t'+stau+$
		    '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+'_hs'+shs+'_zs'+szs+'_'+scaabs+'.dat'
	name1 = rootdir+dir+dir1+filename1
	openw,unit1,name1,/get_lun
	print, 'write ', name1

 ;read the header of the radiation fields files
readf, unit, ss
printf, unit1, ss
readf, unit, tau1
printf, unit1, tau1
readf, unit, ss
printf, unit1, ss
readf, unit, hd
printf, unit1, hd
readf, unit, ss
printf, unit1, ss
readf, unit, zd
printf, unit1, zd
readf, unit, ss
printf, unit1, ss
readf, unit, hdin
printf, unit1, hdin
readf, unit, ss
printf, unit1, ss
readf, unit, zdin
printf, unit1, zdin
readf, unit, ss
printf, unit1, ss
readf, unit, hdsolar
printf, unit1, hdsolar
readf, unit, ss
printf, unit1, ss
readf, unit, zdsolar
printf, unit1, zdsolar
readf, unit, ss
printf, unit1, ss
readf, unit, tau2
printf, unit1, tau2
readf, unit, ss
printf, unit1, ss
readf, unit, hd1
printf, unit1, hd1
readf, unit, ss
printf, unit1, ss
readf, unit, zd1
printf, unit1, zd1
readf, unit, ss
printf, unit1, ss
readf, unit, hd1in
printf, unit1, hd1in
readf, unit, ss
printf, unit1, ss
readf, unit, zd1in
printf, unit1, zd1in
readf, unit, ss
printf, unit1, ss
readf, unit, hd1solar
printf, unit1, hd1solar
readf, unit, ss
printf, unit1, ss
readf, unit, zd1solar
printf, unit1, zd1solar
readf, unit, ss
printf, unit1, ss
readf, unit, pass
printf, unit1, pass
readf, unit, ss
printf, unit1, ss
readf, unit, zs
printf, unit1, zs
readf, unit, ss
printf, unit1, ss
readf, unit, hsin
printf, unit1, hsin
readf, unit, ss
printf, unit1, ss
readf, unit, zsin
printf, unit1, zsin
readf, unit, ss
printf, unit1, ss
readf, unit, hssolar
printf, unit1, hssolar
readf, unit, ss
printf, unit1, ss
readf, unit, zssolar
printf, unit1, zssolar
readf, unit, ss
printf, unit1, ss
readf, unit, rtrun
printf, unit1, rtrun
readf, unit, ss
printf, unit1, ss
readf, unit, xis0
printf, unit1, xis0
readf, unit, ss
printf, unit1, ss
readf, unit, xid0
printf, unit1, xid0
readf, unit, ss
printf, unit1, ss
readf, unit, xid1
printf, unit1, xid1
readf, unit, ss
printf, unit1, ss
readf, unit, idisk1
printf, unit1, idisk1
readf, unit, ss
printf, unit1, ss
readf, unit, idisk2
printf, unit1, idisk2

endif
if morph eq 'b' then begin
	filename1 = 'unit_m'+morph+'_'+filter[ii]+'_'+model+'_q'+qyear+$
	'_t'+stau+$
	'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
	'_reff'+sreff+'_ell'+sellipt+'_'+'n'+snsersic+'_'+scaabs+'.dat'
	name1 = rootdir+dir+dir1+filename1
	openw,unit1,name1,/get_lun
	print, 'write ', name1

;read the header of the radiation fields files
readf, unit, ss
printf, unit1, ss
readf, unit, tau1
printf, unit1, tau1
readf, unit, ss
printf, unit1, ss
readf, unit, hd
printf, unit1, hd
readf, unit, ss
printf, unit1, ss
readf, unit, zd
printf, unit1, zd
readf, unit, ss
printf, unit1, ss
readf, unit, hdin
printf, unit1, hdin
readf, unit, ss
printf, unit1, ss
readf, unit, zdin
printf, unit1, zdin
readf, unit, ss
printf, unit1, ss
readf, unit, hdsolar
printf, unit1, hdsolar
readf, unit, ss
printf, unit1, ss
readf, unit, zdsolar
printf, unit1, zdsolar
readf, unit, ss
printf, unit1, ss
readf, unit, tau2
printf, unit1, tau2
readf, unit, ss
printf, unit1, ss
readf, unit, hd1
printf, unit1, hd1
readf, unit, ss
printf, unit1, ss
readf, unit, zd1
printf, unit1, zd1
readf, unit, ss
printf, unit1, ss
readf, unit, hd1in
printf, unit1, hd1in
readf, unit, ss
printf, unit1, ss
readf, unit, zd1in
printf, unit1, zd1in
readf, unit, ss
printf, unit1, ss
readf, unit, hd1solar
printf, unit1, hd1solar
readf, unit, ss
printf, unit1, ss
readf, unit, zd1solar
printf, unit1, zd1solar
readf, unit, ss
printf, unit1, ss
readf, unit, reff
printf, unit1, reff
readf, unit, ss
printf, unit1, ss
readf, unit, ellipt
printf, unit1, ellipt
readf, unit, ss
printf, unit1, ss
readf, unit, nsersic
printf, unit1, nsersic
readf, unit, ss
printf, unit1, ss
readf, unit, xid0
printf, unit1, xid0
readf, unit, ss
printf, unit1, ss
readf, unit, xid1
printf, unit1, xid1
readf, unit, ss
printf, unit1, ss
readf, unit, idisk1
printf, unit1, idisk1
readf, unit, ss
printf, unit1, ss
readf, unit, idisk2
printf, unit1, idisk2
endif
readf, unit, ss
printf, unit1, ss
readf, unit, dim_positionsr
printf, unit1, dim_positionsr
readf, unit, ss
printf, unit1, ss
readf, unit, dim_positionsz
printf, unit1, dim_positionsz
readf, unit, ss
printf, unit1, ss
readf, unit, ss
printf,unit1, '      r(pc)        z(pc)     urad(erg/pc3/A) eabs(erg/s/pc3/A)'

dim_positions = dim_positionsr * dim_positionsz

u = 1.d
e = 1.d
uu = 1.d
ee = 1.d

for k=0L,dim_positions-1 do begin ;loop in position
	readf,unit,r,z,u,e
	uu = u * ff[ii]
	ee = e * ff[ii]
	printf,unit1,r,z,uu,ee	; prints the scaled value to the file
endfor ;end loop in position

free_lun, unit
free_lun, unit1
endfor ;end loop in wavelength 
mark1: 
end


