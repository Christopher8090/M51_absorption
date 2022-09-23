pro maps,model,qyear,morph,tau,scaabs

rootdir='./'
dir = '../'
dir1 = '/out/'
dir2 = '/indata/'

common param, dist, factor, zs, hs, zs1, hs1, hsin, hs1in, ztrun, rtrun, ztrun1, rtrun1, xis0, xis1, ll, ll1, S0, S0_opt, ms, $
                            zs3, hs3, zs4, hs4, hs3in, hs4in, rtrun3, rtrun4, xis3, xis4, $
                            zs5, hs5, zs6, hs6, hs5in, hs6in, rtrun5, rtrun6, xis5, xis6, $
                            hstin, hdtin, hs1tin, hd1tin, hs3tin, hd3tin, hs4tin, hd4tin, $
                            hs5tin, hd5tin, hs6tin, hd6tin, $
                            rtrund, rtrund1,rtrund3, rtrund4, rtrund5, rtrund6, $
                            zs7, hs7, hs7in, rtrun7, xis7, hs7tin

ii = 0
ss = ''
xis0 = 0L
xis1 = 0L
xid0 = 0L
xid1 = 0L
idisk1 = 0L
idisk2 = 0L
nsersic=0L
xis3 = 0L
xis4 = 0L
xid3 = 0L
xid4 = 0L
idisk3 = 0L
idisk4 = 0L

xis5 = 0L
xis6 = 0L
xid5 = 0L
xid6 = 0L
idisk5 = 0L
idisk6 = 0L

xis7 = 0L
idisk7 = 0L

filename_param='gal_param.in'
name_param = rootdir+dir+dir2+filename_param
openr,unit,name_param,/get_lun
print, 'read ', name_param
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
readf, unit, ss
readf, unit, ss
readf, unit, mopt_b
readf, unit, ss
readf, unit, mopt_n
readf, unit, ss
readf, unit, mopt_i
readf, unit, ss
readf, unit, mopt_m
readf, unit, ss
readf, unit, mopt_o
readf, unit, ss
readf, unit, ss
readf, unit, mstep1_b
readf, unit, ss
readf, unit, mlength1_b
readf, unit, ss
readf, unit, mstep2_b
readf, unit, ss
readf, unit, mlength2_b
readf, unit, ss
readf, unit, mstep3_b
readf, unit, ss
readf, unit, ss
readf, unit, mstep1_n
readf, unit, ss
readf, unit, mlength1_n
readf, unit, ss
readf, unit, mstep2_n
readf, unit, ss
readf, unit, mlength2_n
readf, unit, ss
readf, unit, mstep3_n
readf, unit, ss
readf, unit, ss
readf, unit, mstep1_i
readf, unit, ss
readf, unit, mlength1_i
readf, unit, ss
readf, unit, mstep2_i
readf, unit, ss
readf, unit, mlength2_i
readf, unit, ss
readf, unit, mstep3_i
readf, unit, ss
readf, unit, ss
readf, unit, mstep1_m
readf, unit, ss
readf, unit, mlength1_m
readf, unit, ss
readf, unit, mstep2_m
readf, unit, ss
readf, unit, mlength2_m
readf, unit, ss
readf, unit, mstep3_m
readf, unit, ss
readf, unit, ss
readf, unit, mstep1_o
readf, unit, ss
readf, unit, mlength1_o
readf, unit, ss
readf, unit, mstep2_o
readf, unit, ss
readf, unit, mlength2_o
readf, unit, ss
readf, unit, mstep3_o
readf, unit, ss
readf, unit, ss
readf, unit, elm_b
readf, unit, ss
readf, unit, elm_n
readf, unit, ss
readf, unit, elm_i
readf, unit, ss
readf, unit, elm_m
readf, unit, ss
readf, unit, elm_o
free_lun, unit

dist = dist_in   ;the distance to M33 in Mpc

factor=3.086 * 1.d+16  ;to convert pc in meter
ms = 20.         ;standard central face-on SB
;truncation of the first disk
ztrun=84000.  ;in pc

;truncation of the second disk
ztrun1=84000.  ;in pc


filter_opt=['uv36','b','v','i','j','k','ir36','ir45','ir58']
dim_opt = n_elements(filter_opt)
filter_uv = ['uv09','uv13','uv15','uv16','uv20','uv22','uv25','uv28','uv36']
dim_uv = n_elements(filter_uv)
filter = ['uv09','uv13','uv15','uv16','uv20','uv22','uv25','uv28','uv36','b','v','i','j','k','ir36','ir45','ir58']
dim = n_elements(filter)

lambda=[912.,1350.,1500.,1650.,2000.,2200.,2500.,2800.,3650.,4430.,5640.,8090.,12590.,22000.,36000.,45000.,58000.] * 1.d ;in AA
;the luminosity of the young stellar disk in W/Hz corresponding to
;SFR=1Msolar/yr
ll1 = [0.344,0.905,0.844,0.863,0.908,0.926,0.843,0.910,1.842,2.271,3.837,5.734,0.931,0.728,0.141,0.141,0.141] * 1.d+21

;the luminosity of the old stellar disk in W/Hz corresponding to old=1
ll = [4.771,4.771,9.382,19.54,72.20,64.97,12.58,12.58,12.58] * 1.d+21

;zero point to convert Jy in magnitude
S0 = [460.,249.,460.,554.,710.,728.,759.,954.,1889.,4063.0d0,3636.0d0,2416.0d0,1589.0d0,640.0d0,280.9d0,179.9d0,115.0d0]

S0_opt=[1889.0d0,4063.0d0,3636.0d0,2416.0d0,1589.0d0,640.0d0,280.9d0,179.9d0,115.0d0]

filename = 'geometry.in'
name = rootdir+dir+dir2+filename
openr, unit, name, /get_lun
print, 'read ', name
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
readf, unit, sharpd1
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

;read inner component
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

;read outer component
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
;read inner truncation radius
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
;read nuclear component
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
stau=strcompress(string(round(tau*10)),/remove_all)
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
snsersic=strcompress(string(round(nsersic)),/remove_all)
sinclination=strcompress(string(round(inclination)),/remove_all)

if (morph eq 'td' or morph eq 'tdi' or morph eq 'tdo' or morph eq 'ntd') then begin
	dim_filter = dim
	filter = filter
	scaabs = scaabs
endif
if (morph eq 'd' or morph eq 'di' or morph eq 'do') then begin
	dim_filter = dim_opt
	filter = filter_opt
	scaabs = scaabs
endif

if morph eq 'b' then begin
	dim_filter = dim_opt
	filter = filter_opt
	snsersic = strcompress(string(round(nsersic)),/remove_all)
	scalebulge,model,qyear,tau,scaabs,nsersic,bulge_disk	; added 21/04/21 CJI
	scaabs = scaabs
endif

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


rtrund=rtrund * 1000.d
rtrund1=rtrund1 * 1000.d
rtrund3=rtrund3 * 1000.d
rtrund4=rtrund4 * 1000.d
rtrund5=rtrund5 * 1000.d
rtrund6=rtrund6 * 1000.d

hs7 = hs7 * 1000.d
zs7 = zs7 * 1000.d
hs7in = hs7in * 1000.d
rtrun7 = rtrun7 * 1000.d
hs7tin = hs7tin * 1000.d

;call routine that calculations the scaling factors
conversion_factors_flat,factord_old, factord1_young,factord3_old, factord4_young,factord5_old, factord6_young, factord7_young

factorb_opt = factord_old

;define conversion factord* and S0 values to use
if (morph eq 'td') then begin
   ff = factord1_young
   S0_JJT = S0
endif

if (morph eq 'd') then begin
   ff = factord_old
   S0_JJT = S0_opt
endif

if morph eq 'b' then begin
   ff = factord_old/bulge_disk
   S0_JJT = S0_opt
endif

if (morph eq 'tdi') then begin
   ff = factord4_young
   S0_JJT = S0
endif

if (morph eq 'di') then begin
   ff = factord3_old
   S0_JJT = S0_opt
endif

if (morph eq 'tdo') then begin
   ff = factord6_young
   S0_JJT = S0
endif

if (morph eq 'do') then begin
   ff = factord5_old
   S0_JJT = S0_opt
endif

if (morph eq 'ntd') then begin
	ff = factord7_young
	S0_JJT = S0
endif

if (morph eq 'b') then  begin
	ny=nx_b ;nx and ny flipped from in urad.f
	nx=ny_b
	mopt=mopt_b
	mstep1=mstep1_b
	mlength1=mlength1_b
	mstep2=mstep2_b
	mlength2=mlength2_b
	mstep3=mstep3_b	
	elm=elm_b
endif

if (morph eq 'ntd') then  begin
	ny=nx_n ;nx and ny flipped from in urad.f
	nx=ny_n
	mopt=mopt_n
	mstep1=mstep1_n
	mlength1=mlength1_n
	mstep2=mstep2_n
	mlength2=mlength2_n
	mstep3=mstep3_n
	elm=elm_n		
endif

if (morph eq 'tdi' or morph eq 'di') then  begin
	ny=nx_i ;nx and ny flipped from in urad.f
	nx=ny_i
	mopt=mopt_i
	mstep1=mstep1_i
	mlength1=mlength1_i
	mstep2=mstep2_i
	mlength2=mlength2_i
	mstep3=mstep3_i
	elm=elm_i		
endif

if (morph eq 'td' or morph eq 'd') then  begin
	ny=nx_m ;nx and ny flipped from in urad.f
	nx=ny_m
	mopt=mopt_m
	mstep1=mstep1_m
	mlength1=mlength1_m
	mstep2=mstep2_m
	mlength2=mlength2_m
	mstep3=mstep3_m
	elm=elm_m		
endif

if (morph eq 'tdo' or morph eq 'do') then  begin
	ny=nx_o ;nx and ny flipped from in urad.f
	nx=ny_o
	mopt=mopt_o
	mstep1=mstep1_o
	mlength1=mlength1_o
	mstep2=mstep2_o
	mlength2=mlength2_o
	mstep3=mstep3_o
	elm=elm_o		
endif

for ii = 0, dim_filter-1 do begin ;loop in wavelength
if filter[ii] eq 'nir' then begin
	if ( morph eq 'td' or morph eq 'tdi' or morph eq 'tdo' or morph eq 'ntd') then  begin
		filename = $
		'map_m'+morph+'_k_'+model+'_q'+qyear+$
		'_i'+sinclination+'_t'+stau+'_'+scaabs+'.dat'
		filenameout = 'map_m'+morph+'_k_'+model+'_q'+qyear+$
		'_i'+sinclination+'_t'+stau+$
		'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
		'_hs1_'+shs1+'_zs1_'+szs1+$
		'_'+scaabs;+'.fits'  
	endif   
	if (morph eq 'd' or morph eq 'di' or morph eq 'do') then begin
		filename = $
		'map_m'+morph+'_k_'+model+'_q'+qyear+$
		'_i'+sinclination+'_t'+stau+'_'+scaabs+'.dat'
		filenameout = $
		'map_m'+morph+'_k_'+model+'_q'+qyear+$
		'_i'+sinclination+'_t'+stau+$
		'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
		'_hs'+shs+'_zs'+szs+$
		'_'+scaabs;+'.fits'
	endif
	if morph eq 'b' then begin
		filename = $
		'map_m'+morph+'_k_'+model+'_q'+qyear+$
		'_i'+sinclination+'_t'+stau+'_n'+snsersic+'_'+scaabs+'.dat'
		filenameout = $
		'map_m'+morph+'_k_'+model+'_q'+qyear+$
		'_i'+sinclination+'_t'+stau+$
		'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
		'_reff'+sreff+'_ell'+sellipt+'_n'+snsersic+'_'+scaabs;+'.fits' 
	endif
endif else begin
	if (morph eq 'td' or morph eq 'tdi' or morph eq 'tdo' or morph eq 'ntd') then begin
		filename = $
		'map_m'+morph+'_'+filter[ii]+'_'+model+'_q'+qyear+$
		'_i'+sinclination+'_t'+stau+'_'+scaabs+'.dat' 
		filenameout = $
		'map_m'+morph+'_'+filter[ii]+'_'+model+'_q'+qyear+$
		'_i'+sinclination+'_t'+stau+$
		'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
		'_hs1_'+shs1+'_zs1_'+szs1+$
		'_'+scaabs;+'.fits'  
	endif
	if (morph eq 'd' or morph eq 'di' or morph eq 'do') then begin
		filename = $
		'map_m'+morph+'_'+filter[ii]+'_'+model+'_q'+qyear+$
		'_i'+sinclination+'_t'+stau+'_'+scaabs+'.dat'
		filenameout = $
		'map_m'+morph+'_'+filter[ii]+'_'+model+'_q'+qyear+$
		'_i'+sinclination+'_t'+stau+$
		'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
		'_hs'+shs+'_zs'+szs+$
		'_'+scaabs;+'.fits'
	endif
	if morph eq 'b' then  begin
		filename = $
		'map_m'+morph+'_'+filter[ii]+'_'+model+'_q'+qyear+$
		'_i'+sinclination+'_t'+stau+'_n'+snsersic+'_'+scaabs+'.dat'
		filenameout = $
		'map_m'+morph+'_'+filter[ii]+'_'+model+'_q'+qyear+$
		'_i'+sinclination+'_t'+stau+$
		'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
		'_reff'+sreff+'_ell'+sellipt+'_n'+snsersic+'_'+scaabs;+'.fits'
	endif
endelse
name = rootdir+dir+'/out/'+filename
print, 'read ', name
openr,unit,name,/get_lun


;read in counts
z = dblarr(nx,ny)
zx = dblarr(ny)
;interpolate if mask is set
if (mopt eq 1) then begin
	x = findgen(nx)
	y = findgen(ny)
	;length up to which mstep1 applies
	dim1 = mlength1
	;length up to which mstep2 applies
	dim2 = mlength2
	;the number of pixels within mlength1 that are unmasked
	dimin1 = (dim1-1)/mstep1
	;initialise arrays
	yin1 = fltarr(dimin1)
	xin1 = fltarr(dimin1)

	;the number of pixels between mlength2 and mlength1 that are unmasked
	dimin2 = (dim2-dim1)/mstep2
	yin2 = fltarr(dimin2)
	xin2 = fltarr(dimin2)

	dim3 = nx-1
	;the number of pixels between mlength2 and the map limit that are unmasked
	dimin3 = (dim3-dim2)/mstep3
        yin3 = fltarr(dimin3)
        xin3 = fltarr(dimin3)	
	if (morph eq 'td' or morph eq 'd') then begin
		yin3 = fltarr(dimin3+1)
		xin3 = fltarr(dimin3+1)
	endif
	;read in data from ascii (extinction loop.pro)
	;read pixels up to mlength1 in increments of mstep1
	for i = 0, dim1-1 do begin
		readf, unit, zx
		z[i,*] = zx[*]
	endfor

	for i = dim1, dim2-1 do begin
		readf,unit,zx
		z[i,*] = zx[*]
	endfor

	for i = dim2, dim3 do begin
		readf,unit,zx
		z[i,*] = zx[*]
	endfor


	;ellipse to define where the elliptical mask is to have zero value.
	;define ellipse parameters
	ratio=(cos(inclination*!DPI/180.))^(-1.)   ;ratio between major and minor axis
	pos_ang = 90.   
	xc = 0
	yc = (ny/2)-1
	dim = size(z, /dimensions)
	;produce ellipse to be used for annuli
	dist_ellipse, el_map, dim, xc, yc, ratio, pos_ang 
	;make sure only resgion lt the eliptical mask radius "elm" is set to zero
	if (elm ne 0) then begin
		iq = where(el_map lt elm)
		z[iq] = 0
	endif
	for jj = 0,ny-1 do begin
		k = 0L
		for ll = 1, dim1-1, mstep1 do begin
			yin1[k] = z[ll,jj]
			xin1[k] = x[ll]
			k = k+1
		endfor   
		int = interpol(yin1, xin1, x(1:dim1), /spline)
		z(1:dim1,jj) = int
	endfor
	
	for jj = 0,ny-1 do begin
		k = 0L
		for ll = dim1, dim2-mstep2, mstep2 do begin    
			yin2(k) = z(ll,jj)
			xin2(k) = x(ll)
			k = k+1
		endfor   
		int=interpol(yin2, xin2, x(dim1:dim2), /spline)
		z(dim1:dim2,jj) = int
	endfor
	
	for jj = 0,ny-1 do begin
		k = 0L
		for ll = dim2+(mstep3-mstep2), dim3, mstep3 do begin    
			yin3(k) =z(ll,jj)
			xin3(k)= x(ll)
			k = k+1
		endfor   
		int=interpol(yin3,xin3,x(dim2:dim3),/spline)
		z(dim2:dim3,jj) = int
	endfor
	inrb = where (z eq -1,nrb)
	if nrb gt 0 then z(inrb) = 0.

	totalcount1 = mstep1 * total(z[1:mlength1,0:ny-1])
	totalcount2 = mstep2 * total(z[mlength1+1:mlength2,0:ny-1])
	totalcount3 = mstep3 * total(z[mlength2+1:nx-1,0:ny-1])
	z_tot = (totalcount1 + totalcount2 + totalcount3)+ total(z[0:0,0:ny-1]) ;total is only half what the full galaxy would be

endif else begin ;resume as normal if no mask
	for i=0,nx-1 do begin
		readf,unit,zx
		z(i,*) = zx(*)
	endfor

;ellipse to define where the elliptical mask is to have zero value.
   ;define ellipse parameters
   ratio=(cos(inclination*!DPI/180.))^(-1.)   ;ratio between major and minor axis
   pos_ang = 90.   
   xc=0
   yc=(ny/2)-1
   dim=size(z,/dimensions)
   ;produce ellipse to be used for annuli
   dist_ellipse, el_map, dim, xc, yc, ratio, pos_ang 
   ;make sure only resgion lt the eliptical mask radius "elm" is set to zero
   if (elm ne 0) then begin
	   iq=where(el_map lt elm)
	   z(iq)=0
   endif
	z_tot=total(z)
endelse	

;calculate counts to Jy conversion factor
bfactor=1
z_mag=-2.5*ALOG10(z_tot*bfactor) ;z_tot as mag
z_tot_Jy=S0_JJT[ii]*10^(-z_mag/2.5) ;z_tot in Jy
z_count_Jy=z_tot_Jy/z_tot ;Jy/count

;apply calibration factors to z
z_cal=z*z_count_Jy*ff[ii] ;calibrated surface brightness
free_lun,unit
nameout=rootdir+dir+dir1+filenameout+'.fits'
print, 'write', nameout
writefits,nameout,z_cal
endfor
end
