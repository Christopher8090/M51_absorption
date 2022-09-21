pro param

root = '/nfs/d58/vlk/sedmodel/cinman/m51a/'
nurad_indata = root+'NUrad_M51a/indata/'

common param, shd, szd, shd1, szd1, shs, szs, shs1, szs1, sreff, sellipt, n_radial_positions, $
	      n_vertical_positions, bands, wavelengths_AA, model, qyear, scaabs, stau, strsfr, strold, strbd, $
	      tau, sfr, sfr4, sfr6, sfr7, old, old3, old5, bd, ffactor, ffactor4, ffactor6, ffactor7, $
	      f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5, f_bd, $
	      distance_mpc, inclination_deg, sinclination, pixsize_arcsec, snsersic,$
	      nx_b, ny_b, nx_n, ny_n, nx_i, ny_i, nx_m, ny_m, nx_o, ny_o, $
	tau1, tau2, hd, zd, hdin, zdin, hdsolar, zdsolar, hd1, zd1, hd1in, zd1in, hd1solar, zd1solar, $
	h_bdisk, h_vdisk, h_idisk, h_jdisk, h_kdisk, h_ir34disk, h_ir45disk, h_ir58disk, zs, hsin, zsin, $
	hssolar, zssolar, hs1, zs1, hs1in, zs1in, hs1solar, zs1solar, rtruncate, sharp, rtruncated, sharpd, $
	rtruncate1, sharp1, rtruncated1, sharpd1, reff, ellipt, nsersic, xis0, xis1, xid0, xid1, idisk1, idisk2, $
	$
	tau3, tau4, hd3, zd3, hd3in, zd3in, hd3solar, zd3solar, hd4, zd4, hd4in, zd4in, hd4solar, zd4solar, $
        h_bdisk3, h_vdisk3, h_idisk3, h_jdisk3, h_kdisk3, h_ir34disk3, h_ir45disk3, h_ir58disk3, zs3, hs3in, zs3in, $
        hs3solar, zs3solar, hs4, zs4, hs4in, zs4in, hs4solar, zs4solar, rtruncate3, sharp3, rtruncated3, sharpd3, $
        rtruncate4, sharp4, rtruncated4, sharpd4, xis3, xis4, xid3, xid4, idisk3, idisk4, $
	$
        tau5, tau6, hd5, zd5, hd5in, zd5in, hd5solar, zd5solar, hd6, zd6, hd6in, zd6in, hd6solar, zd6solar, $
        h_bdisk5, h_vdisk5, h_idisk5, h_jdisk5, h_kdisk5, h_ir34disk5, h_ir45disk5, h_ir58disk5, zs5, hs5in, zs5in, $
        hs5solar, zs5solar, hs6, zs6, hs6in, zs6in, hs6solar, zs6solar, rtruncate5, sharp5, rtruncated5, sharpd5, $
        rtruncate6, sharp6, rtruncated6, sharpd6, xis5, xis6, xid5, xid6, idisk5, idisk6, $
	hstin, hdtin, hs1tin, hd1tin, hs3tin, hd3tin, hs4tin, hd4tin, hs5tin, hd5tin, hs6tin, hd6tin, $
$
	hii_h1, hii_z1, hii_h1in, hii_z1in, hii_h1solar, hii_rtruncate1, hii_sharp1, hii_xi1, hii_rtin1, $
	hii_h4, hii_z4, hii_h4in, hii_z4in, hii_h4solar, hii_rtruncate4, hii_sharp4, hii_xi4, hii_rtin4, disks



bands=['uv09','uv13','uv15','uv16','uv20','uv22','uv25','uv28','uv36','b','v','i','j','k','ir36','ir45','ir58']
wavelengths_AA=[912.,1350.,1500.,1650.,2000.,2200.,2500.,2800.,3650.,$
                4430.,5640.,8090.,12590.,22000.,36000.,45000.,58000.]       ; [um]
disks = ['tot','tdi','td','tdo','di','d','do','b']

;Read scaling.in
readcol, nurad_indata+'scaling.in', x, y, skipline=0, format="(A,D)"
model = (x[0])[0]
qyear = (x[1])[0]
scaabs = (x[2])[0]
tau = (y[where(x eq "tau")])[0]
sfr = (y[where(x eq "sfr")])[0]
sfr4 = (y[where(x eq "sfr4")])[0]
sfr6 = (y[where(x eq "sfr6")])[0]
sfr7 = (y[where(x eq "sfr7")])[0]
old = (y[where(x eq "old")])[0]
old3 = (y[where(x eq "old3")])[0]
old5 = (y[where(x eq "old5")])[0]
bd = (y[where(x eq "bd")])[0]
ffactor = (y[where(x eq "ffactor")])[0]
ffactor4 = (y[where(x eq "ffactor4")])[0]
ffactor6 = (y[where(x eq "ffactor6")])[0]
ffactor7 = (y[where(x eq "ffactor7")])[0]

stau = strcompress(string(fix(tau*10.)),/remove_all)
strsfr = strcompress(string(fix(sfr*100.)),/remove_all)
strold = strcompress(string(fix(old*100.)),/remove_all)
strbd = strcompress(string(fix(bd*100.)),/remove_all)

;Read wavelength dependent amplitude scaling parameters ff_scaling.in
READCOL, nurad_indata+'ff_scaling.in', wave, f_uv, f_uv4, f_uv6, f_uv7, f_BVIK, f_BVIK3, f_BVIK5, f_bd,$
	 SKIPLINE=1, FORMAT="(A,D,D,D,D,D,D,D,D)"
f_BVIK = f_BVIK[8:16]
f_BVIK3 = f_BVIK3[8:16]
f_BVIK5 = f_BVIK5[8:16]
f_bd = f_bd[8:16]
ss=''
;Read HII_geometry.in
openr, hii_lun, nurad_indata+'wd01_sca/hii_geometry.in', /get_lun
readf, hii_lun, ss
readf, hii_lun, hii_h1
readf, hii_lun, ss
readf, hii_lun, hii_z1
readf, hii_lun, ss
readf, hii_lun, hii_h1in
readf, hii_lun, ss
readf, hii_lun, hii_z1in
readf, hii_lun, ss
readf, hii_lun, hii_h1solar
readf, hii_lun, ss
readf, hii_lun, ss
readf, hii_lun, ss
readf, hii_lun, hii_rtruncate1
readf, hii_lun, ss
readf, hii_lun, hii_sharp1
readf, hii_lun, ss
readf, hii_lun, hii_xi1
readf, hii_lun, ss
readf, hii_lun, hii_rtin1
readf, hii_lun, ss
readf, hii_lun, ss
readf, hii_lun, hii_h4
readf, hii_lun, ss
readf, hii_lun, hii_z4
readf, hii_lun, ss
readf, hii_lun, hii_h4in
readf, hii_lun, ss
readf, hii_lun, hii_z4in
readf, hii_lun, ss
readf, hii_lun, hii_h4solar
readf, hii_lun, ss
readf, hii_lun, hii_rtruncate4
readf, hii_lun, ss
readf, hii_lun, hii_sharp4
readf, hii_lun, ss
readf, hii_lun, hii_xi4
readf, hii_lun, ss
readf, hii_lun, hii_rtin4
free_lun, hii_lun
ss=''
;Read geometry.in
openr, unit, nurad_indata+'wd01_sca/geometry.in', /get_lun
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
readf, unit, rtruncate
readf, unit, ss
readf, unit, sharp
readf, unit, ss
readf, unit, rtruncated
readf, unit, ss
readf, unit, sharpd
readf, unit, ss
readf, unit, rtruncate1
readf, unit, ss
readf, unit, sharp1
readf, unit, ss
readf, unit, rtruncated1
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
readf, unit, rtruncate3
readf, unit, ss
readf, unit, sharp3
readf, unit, ss
readf, unit, rtruncated3
readf, unit, ss
readf, unit, sharpd3
readf, unit, ss
readf, unit, rtruncate4
readf, unit, ss
readf, unit, sharp4
readf, unit, ss
readf, unit, rtruncated4
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
readf, unit, rtruncate5
readf, unit, ss
readf, unit, sharp5
readf, unit, ss
readf, unit, rtruncated5
readf, unit, ss
readf, unit, sharpd5
readf, unit, ss
readf, unit, rtruncate6
readf, unit, ss
readf, unit, sharp6
readf, unit, ss
readf, unit, rtruncated6
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
free_lun, unit

hs = h_bdisk
shd=strcompress(string(fix(hd*1000.)),/remove_all)
szd=strcompress(string(fix(zd*1000.)),/remove_all)
shd1=strcompress(string(fix(hd1*1000.)),/remove_all)
szd1=strcompress(string(fix(zd1*1000.)),/remove_all)
shs=strcompress(string(fix(hs*1000.)),/remove_all)
szs=strcompress(string(fix(zs*1000.)),/remove_all)
shs1=strcompress(string(fix(hs1*1000.)),/remove_all)
szs1=strcompress(string(fix(zs1*1000.)),/remove_all)
sreff=strcompress(string(fix(reff*1000.)),/remove_all)
sellipt=strcompress(string(fix(ellipt*100.)),/remove_all)
snsersic=strcompress(string(fix(nsersic)),/remove_all)

;Read gal_param.in
openr, gplun, nurad_indata+'gal_param.in', /get_lun
	while ss ne 'dist[Mpc]' do readf, gplun, ss
	readf, gplun, distance_mpc
	readf, gplun, ss
	readf, gplun, inclination_deg
	while ss ne 'pixel size[arcsec]' do readf, gplun, ss
	readf, gplun, pixsize_arcsec
	readf, gplun, ss
	readf, gplun, nx_b
	readf, gplun, ss
	readf, gplun, ny_b
	readf, gplun, ss
	readf, gplun, nx_n
	readf, gplun, ss
	readf, gplun, ny_n
	readf, gplun, ss
	readf, gplun, nx_i
	readf, gplun, ss
	readf, gplun, ny_i
	readf, gplun, ss
	readf, gplun, nx_m
	readf, gplun, ss
	readf, gplun, ny_m
	readf, gplun, ss
	readf, gplun, nx_o
	readf, gplun, ss
	readf, gplun, ny_o
free_lun, gplun

sinclination = strcompress(string(fix(inclination_deg)),/remove_all)

;Read urad.in
openr, urad_lun, nurad_indata+'urad.in', /get_lun
	;Skip the first two lines (iurad and icang respectively)
	READF, urad_lun, ss
	READF, urad_lun, ss
	;Read the number of sampled radial positions
	READF, urad_lun, n_radial_positions
	;Skip the sampling positions
	FOR i=0, n_radial_positions-1 DO READF, urad_lun, ss
	;Read the number of sampled vertical positions
	READF, urad_lun, n_vertical_positions
;Close urad.in
FREE_LUN, urad_lun
end
