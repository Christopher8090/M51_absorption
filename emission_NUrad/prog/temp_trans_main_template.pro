;this routines calculates the temperature distributions for each grain size and
;composition placed in the radiation fields that are inputed to this routine
;it uses RF having a different spectral shape for old
 
pro temp_trans_main_template,model,qyear,tau,sfr,old,bd,scaabs,check 
;temp_trans_main_template,'wd01','06',6.6,1.0,0.8,0.49,'abs','no' ;call for MW
;do not confuse the keyword "qyear" with "qcomp_year"!!

;input:
; model - dust model, e.g. 'wd01'
; qyear - e.g.'06'
; tau - total B band face-on central opacity
; sfr - star formation rate in Msolar/yr; is only needed for the file names, 
;       since the input radiation fields are already scaled for the right SFR
; old - scaling factor for the optical emission; is only needed for the 
;       file names, since the input radiation fields are already scaled for 
;       the right value
; bd -  bulge-to-disk ratio; is only needed for the file names, since the input
;       radiation fields are scaled for the right BD
; check - a logical switch if check='yes' then the ISRF is considered;
;                          if check='no' then the radiation fields in the 
;                                        galaxy are considered

;subroutines:
;temp_trans.pro
;q_grain.pro
;planck_grain.pro
;u_l_isrf.pro
start_time = systime(/seconds)
common scaling
close,/all

dir = emission_dir+'outdata/'
dir1 = emission_dir+'indata/'
dir2 = emission_dir+'outdata_temp/'
dir3 = urad_dir+'indata/'

;define some strings for file names
suv=strcompress(string(round(sfr*100)),/remove_all)
sbd=strcompress(string(round(bd*100)),/remove_all)
sold=strcompress(string(round(old*100)),/remove_all)
stau=strcompress(string(round(tau*10)),/remove_all)
;stau='86096'

;read file with the geometry of the system
ss = ''
idisk1 = 0L
idisk2 = 0L

filename = 'geometry.in'
name = dir3+filename
openr,unit,name,/get_lun
;print, 'read ', name
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

;read the file with the dust model parameters and grain size distribution
if model eq 'wd01' then namer1=rootdir+dir+'grain_sizeswd01_q'+qyear+'.dat'
if model eq 'lmc1' then namer1=rootdir+dir+'grain_sizeslmc1_q'+qyear+'.dat'
if model eq 'wd01_c60' then namer1=rootdir+dir+'grain_sizeswd01_q'+qyear+'.dat'
if model eq 'wd01_c50' then namer1=rootdir+dir+'grain_sizeswd01_q'+qyear+'.dat'
if model eq 'wd01_c40' then namer1=rootdir+dir+'grain_sizeswd01_q'+qyear+'.dat'
if model eq 'wd01_c30' then namer1=rootdir+dir+'grain_sizeswd01_q'+qyear+'.dat'
openr,unit11,namer1,/get_lun

ss = ' '
dim_size = 0L
dim_comp=0L
dim_Draine=0L
ig_amin=0L
ig_amax=0L
readf,unit11,ss
readf,unit11,dim_comp
comp_array = strarr(dim_comp)
readf,unit11,ss
for ii = 0, dim_comp-1 do begin
	readf, unit11, ss
	comp_array(ii) = ss
endfor

comp = ''
qcomp_year = ''
for j = 0L, dim_comp-1 do begin	; begin loop in composition
	readf,unit11,ss		;
	readf,unit11,comp               ;grain composition
	readf,unit11,ss                 ;
	readf,unit11,qcomp_year         ;model for Q (year when updated)
	readf,unit11,ss                 ;
	readf,unit11,dim_Draine         ;number of grain sizes tabulated by Draine
	readf,unit11,ss                 ;
	readf,unit11,amin               ;minimum grain size in micron
	readf,unit11,ss                 ;
	readf,unit11,amax               ;maximum grain size in micron
	readf,unit11,ss                 ;
	readf,unit11,ig_amin            ;index of amin in Draine sizes
	readf,unit11,ss                 ;
	readf,unit11,ig_amax            ;index of amax in Draine sizes
	readf,unit11,ss                 ;
	readf,unit11,s_grain            ;grain density in g/cm^3
	readf,unit11,ss                 ;
	readf,unit11,dim_size           ;number of grain sizes
	readf,unit11,ss                 ;
	readf,unit11,weight             ;
	readf,unit11,ss                 ;
	readf,unit11,ss                 ;
	readf,unit11,ss                 ;
	readf,unit11,ss                 ;
	a = dblarr(dim_size)           ;grain sizes in cm
	na = dblarr(dim_size)          ;grain size distribution
	
	for ii = 0L,dim_size-1 do begin
		readf,unit11,format='(2(e13.6))',pass1,pass2
		a(ii) = pass1 * 1.D+4 ;transform from cm to micron
		na(ii) = pass2
	endfor
	
	;read the emissivities from Draine
	q_grain,comp,qcomp_year,lambda_q,qq_grain,qq_grain_sc,qq_grain_ph
	ncheck = n_elements(qq_grain(*,0))
	if ncheck ne dim_Draine then begin
	print,'unexpected number of grain sizes; program stops'
	goto, mark1
	endif
	lambda_q = lambda_q*1.E4         ;in Angstroem
	dim_lambda=n_elements(lambda_q)
	
	;truncate the emissivities from Draine to the range of sizes given in the input
	;there is no interpolation or truncation in the lambda from Draine
	qq_abs = dblarr(dim_size,dim_lambda)
	qq_abs  = qq_grain[ig_amin:ig_amax,*]
	
	;call the routine that provides the averaged Planck emissivity
	;read predifined table from Joerg and interpolates according to grains size
	;no interpolation in temperature compared to the sampling in Joerg's table
	;Qp_grain contain the averaged emissivity and TT_Q contains the respective
	;dust temperature 
	planck_grain,comp,qcomp_year,a,tt_q,qp_grain
	
	kkk = 0L ;counter in grain size
	
	if check eq 'yes' then begin
		dim_wave=500L
		a1=Alog10(9.120D2)
		a2=ALOG10(1.d+7)
		a3=(a2-a1)/(dim_wave-1)
		lu_check=10.D^(findgen(dim_wave)*a3+a1)
		ul_check = u_l_isrf(lu_check)
		fact = 3.d-57
		ul_check = 1.d-10 * ul_check/fact
	endif
	
	filter = ['uv09','uv13','uv16','uv20','uv22','uv25','uv28','uv36',$
	           'b','v','i','j','k','ir36','ir45','ir58']
	lu = [912.,1350.,1650.,2000.,2200.,2500.,2800.,3650.,$
		4430.,5640.,8090.,12590.,22000.,34000.,45000.,58000.] * 1.d
	dim_lu = n_elements(lu)
	ul = dblarr(dim_lu)
	
	;read files with radiation fields
	namein = model+'_q'+qyear+$
	              '_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
	              '_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
	              '_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
	              '_reff'+sreff+'_ell'+sellipt+'_'+scaabs+'.dat'
	unitsave = lonarr(dim_lu)
	for ii = 0, dim_lu-1 do begin
		openr,iiunit,dir1+'u_'+filter(ii)+'_'+namein,/get_lun
		unitsave(ii) = iiunit
	endfor 

	dim_positionsr = 0L
	dim_positionsz = 0L
	idisk1 = 0L
	idisk2 = 0L
	;read the header of the radiation fields files
	;the first line contain the input parameters
	unit0 = unitsave(0)
	readf, unit0, ss
	readf, unit0, newtau
	readf, unit0, ss
	readf, unit0, newsfr
	readf, unit0, ss
	readf, unit0, newold
	readf, unit0, ss
	readf, unit0, newbd
	;check that the input parameters read from the file are the same as those
	;inputed to the program (since we have a multiple definition of parameters)
	if newtau ne tau then begin
		print, 'unexpected tau; program stops'
		goto, mark1
	endif
	if newsfr ne sfr then begin
		print, 'unexpected sfr; program stops'
		goto, mark1
	endif
	if newbd ne bd then begin
		print, 'unexpected bd; program stops'
		goto, mark1
	endif
	if newold ne old then begin
		print, 'unexpected old; program stops'
		goto, mark1
	endif

	;read the rest of the header
	readf, unit0, ss
	readf, unit0, tau1
	readf, unit0, ss
	readf, unit0, hd
	readf, unit0, ss
	readf, unit0, zd
	readf, unit0, ss
	readf, unit0, hdin
	readf, unit0, ss
	readf, unit0, zdin
	readf, unit0, ss
	readf, unit0, hdsolar
	readf, unit0, ss
	readf, unit0, zdsolar
	readf, unit0, ss
	readf, unit0, tau2
	readf, unit0, ss
	readf, unit0, hd1
	readf, unit0, ss
	readf, unit0, zd1
	readf, unit0, ss
	readf, unit0, hd1in
	readf, unit0, ss
	readf, unit0, zd1in
	readf, unit0, ss
	readf, unit0, hd1solar
	readf, unit0, ss
	readf, unit0, zd1solar
	readf, unit0, ss
	readf, unit0, hs
	readf, unit0, ss
	readf, unit0, zs
	readf, unit0, ss
	readf, unit0, hsin
	readf, unit0, ss
	readf, unit0, zsin
	readf, unit0, ss
	readf, unit0, hssolar
	readf, unit0, ss
	readf, unit0, zssolar
	readf, unit0, ss
	readf, unit0, rtrun
	readf, unit0, ss
	readf, unit0, hs1
	readf, unit0, ss
	readf, unit0, zs1
	readf, unit0, ss
	readf, unit0, hs1in
	readf, unit0, ss
	readf, unit0, zs1in
	readf, unit0, ss
	readf, unit0, hs1solar
	readf, unit0, ss
	readf, unit0, zs1solar
	readf, unit0, ss
	readf, unit0, rtrun1
	readf, unit0, ss
	readf, unit0, xis0
	readf, unit0, ss
	readf, unit0, xis1
	readf, unit0, ss
	readf, unit0, xid0
	readf, unit0, ss
	readf, unit0, xid1
	readf, unit0, ss
	readf, unit0, idisk1
	readf, unit0, ss
	readf, unit0, idisk2
	readf, unit0, ss
	readf, unit0, dim_positionsr
	readf, unit0, ss
	readf, unit0, dim_positionsz
	readf, unit0, ss
	readf, unit0, ss
	
	dim_positions = dim_positionsr * dim_positionsz
	if check eq 'yes' then dim_positions= 1L
	;skip the header for the rest of the rest of
	;the files with the radiation fields 
	for i = 0,81 do begin
		for ii = 1, dim_lu-1 do begin
			iiunit = unitsave(ii)
			readf, iiunit, ss
		endfor
	endfor

	if check eq 'no' then  begin 
		Name='PT_'+comp+'_'+model+'_'+qyear+$
			'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
			'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
			'_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
			'_reff'+sreff+'_ell'+sellipt+'_'+scaabs+'.dat'
		Name1='energy_absorbed_grain_'+comp+'_'+model+'_'+qyear+$
			'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+$
			'_hd'+shd+'_zd'+szd+'_hd1_'+shd1+'_zd1_'+szd1+$
			'_hs'+shs+'_zs'+szs+'_hs1_'+shs1+'_zs1_'+szs1+$
			'_reff'+sreff+'_ell'+sellipt+'_'+scaabs+'.dat' 
	endif

	if check eq 'yes' then begin
		Name='PT_'+comp+'_'+model+'_'+qyear+'_isrf.dat'
		Name1='energy_absorbed_grain_'+comp+'_'+model+'_'+qyear+'_isrf.dat' 
	endif

	openw, unit10, dir2+Name,/get_lun
	openw, unit200, dir2+Name1,/get_lun
	print, 'Write: '+comp
	printf, unit10, 'central face-on tau B'
	printf, unit200, 'central face-on tau B'
	printf, unit10, tau
	printf, unit200, tau
	printf, unit10, 'sfr'
	printf, unit200, 'sfr'
	printf, unit10, sfr
	printf, unit200, sfr
	printf, unit10, 'old'
	printf, unit200, 'old'
	printf, unit10, old
	printf, unit200, old
	printf, unit10, 'bd'
	printf, unit200, 'bd'
	printf, unit10, bd
	printf, unit200, bd
	printf, unit10, 'tau1'
	printf, unit10, tau1
	printf, unit10, 'hd [pc]'
	printf, unit10, hd
	printf, unit10, 'zd [pc]'
	printf, unit10, zd
	printf, unit10, 'hdin [pc]'
	printf, unit10, hdin
	printf, unit10, 'zdin [pc]'
	printf, unit10, zdin
	printf, unit10, 'hdsolar [pc]'
	printf, unit10, hdsolar
	printf, unit10, 'zdsolar [pc]'
	printf, unit10, zdsolar
	printf, unit10, 'tau2'
	printf, unit10, tau2
	printf, unit10, 'hd1 [pc]'
	printf, unit10, hd1
	printf, unit10, 'zd1 [pc]'
	printf, unit10, zd1
	printf, unit10, 'hd1in [pc]'
	printf, unit10, hd1in
	printf, unit10, 'zd1in [pc]'
	printf, unit10, zd1in
	printf, unit10, 'hd1solar [pc]'
	printf, unit10, hd1solar
	printf, unit10, 'zd1solar [pc]'
	printf, unit10, zd1solar
	printf, unit10, 'hs [pc]'
	printf, unit10, hs
	printf, unit10, 'zs [pc]'
	printf, unit10, zs
	printf, unit10, 'hsin [pc]'
	printf, unit10, hsin
	printf, unit10, 'zsin [pc]'
	printf, unit10, zsin
	printf, unit10, 'hssolar [pc]'
	printf, unit10, hssolar
	printf, unit10, 'zssolar [pc]'
	printf, unit10, zssolar
	printf, unit10, 'hs1 [pc]'
	printf, unit10, hs1
	printf, unit10, 'zs1 [pc]'
	printf, unit10, zs1
	printf, unit10, 'hs1in [pc]'
	printf, unit10, hs1in
	printf, unit10, 'zs1in [pc]'
	printf, unit10, zs1in
	printf, unit10, 'hs1solar [pc]'
	printf, unit10, hs1solar
	printf, unit10, 'zs1solar [pc]'
	printf, unit10, zs1solar
	printf, unit10, 'rtrun [pc]'
	printf, unit10, rtrun
	printf, unit10, 'rtrun1 [pc]'
	printf, unit10, rtrun1
	printf, unit10, 'xis0'
	printf, unit10, xis0
	printf, unit10, 'xis1'
	printf, unit10, xis1
	printf, unit10, 'xid0'
	printf, unit10, xid0
	printf, unit10, 'xid1'
	printf, unit10, xid1
	printf, unit10, 'idisk1'
	printf, unit10, idisk1
	printf, unit10, 'idisk2'
	printf, unit10, idisk2
	printf, unit10, 'number of r positions'
	printf, unit10, dim_positionsr
	printf, unit10, 'number of z positions'
	printf, unit10, dim_positionsz
	printf, unit10
	printf, unit10,  'number of positions in the galaxy:'
	printf, unit10, dim_positions
	printf, unit10
	printf, unit10,  'Number of grain sizes:'
	printf, unit10, dim_size
	printf, unit10
	printf, unit200,  'number of positions in the galaxy:'
	printf, unit200, dim_positions
	printf, unit200
	printf, unit200,  'Number of grain sizes:'
	printf, unit200, dim_size
	printf, unit200
	
	uu = 1.d
	for k=0L,dim_positions-1 do begin ;start loop in position in the galaxy
		for ii = 0,dim_lu-1 do begin
			iiunit = unitsave(ii)
			readf,iiunit,r,z,uu,d
			ul(ii) = uu
		endfor
		if check eq 'yes' then begin
			lu = lu_check
			ul = ul_check
		endif
		printf, unit10, '; R z in pc'
		printf, unit10,  r, z
		printf, unit200, '; R z in pc'
		printf, unit200,  r, z
		printf, unit200
		printf, unit200, 'a      energy absorbed'
		printf, unit200, 'micron  [J/s]'
		
		;call the subroutine that calculates the PT curves
		temp_trans_template,comp,model,qyear,qcomp_year,dim_size,a,lu,ul,lambda_q,qq_abs,$
			tt_q,qp_grain,stau,suv,sold,sbd,dir2,inclbigbang=0,energy_abs,check,scaabs
		;open temporary file for reading
		Name=''
		Name = 'PT_'+comp+'_'+model+'_'+qyear
		if check eq 'no' then begin
			Name = Name+'_t'+stau+'_s'+suv+'_no'+sold+'_bd'+sbd+'_'+scaabs+'.temp'
		endif
		if check eq 'yes' then begin
			Name = Name+'_isrf.temp'
		endif
		openr, unit100, dir2+Name,/get_lun
		ss=''
		Zahl = 0L
		REPEAT BEGIN
		; counting counter clockwise
		; start with largest grain
		readf, unit100, ss
		readf, unit100, kkk
		readf, unit100, ss
		readf, unit100, ss
		readf, unit100, ss
		readf, unit100, aa
		readf, unit100, ss
		readf, unit100, T_min,T_max
		readf, unit100, ss
		readf, unit100, Zahl
		readf, unit100, ss
		readf, unit100, ss
		TT_m=dblarr(Zahl-1)
		prob=dblarr(Zahl-1)
		dT=dblarr(Zahl-1)
		for kk = 0L, Zahl-2 do begin
			readf, unit100, pass1,pass2,pass3
			TT_m(kk)=pass1
			prob(kk)=pass2
			dT(kk) = pass3
		endfor
		printf, unit10
		printf, unit10, kkk
		printf, unit10, '; Temperature of grains ('+comp+')'
		printf, unit10
		printf, unit10, '; grain size in micron'
		printf, unit10,  aa
		printf, unit10, '; T-Interval (Kelvin)'
		printf, unit10, T_min,T_max
		printf, unit10, '; number of points : '
		printf, unit10, Zahl
		printf, unit10
		printf, unit10, '    <T>     dp(T)'
		for kk = 0L, Zahl-2 do printf, unit10, TT_m[kk],prob[kk], dT[kk]
		printf, unit200, aa,energy_abs[kkk]
		ENDREP UNTIL kkk EQ 0
		free_lun, unit100 ;close temporary file
	endfor ;finish loop in galaxy position

	for ii = 0, dim_lu-1 do begin
		iiunit = unitsave(ii)
		free_lun, iiunit
	endfor
	free_lun, unit10
	free_lun,unit200
endfor ;finish loop in composition

mark1:
free_lun,unit11
print, 'DONE: temp_trans_main_template.pro ('+strtrim(ceil(systime(/seconds)-start_time),1)+' s)'
end
