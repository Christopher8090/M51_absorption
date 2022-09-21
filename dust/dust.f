c       changes 21.11.2012
c       we increase the size of the cylinder where we do the calculations
c       rhmax1 = 5.0 * max(hd,hs)
c       d1=15.0* max(zs, reff)

        Subroutine dust(mask,dzmodel)
c----------------------------------------------------------------------------
c	Computation of the surface brightness of a galaxy that contains dust.
c
c       Developed by Nick Kylafis in 1984.  Re-examined in May 1991 and 
c	some small errors were corrected.
c
c	The following improvements have been made by Manolis Xilouris.
c	February 1995, Double Precision is now available.
c	August 1996, The Hubble profile and the R^1/4 law are used to 
c	describe the bulge.
c
c	The code is given by Nick Kylafis to his collaborators.  They 
c	SHOULD NOT give the code to anybody else.  People interested in
c	the code should contact Nick kylafis (e-mail: kylafis@physics.uoc.gr)
c
c	Heraklion 28 February 2001.
c
c       April 2001:
c       Adapted for use as subroutine in modelfitting. Parameters formerly
c       read from file dust.in now passed into routine via common. 
c       The parameter ellipt was added to common bulge. The model
c       brightness distribution is passed through the array dzmodel.
c
c       18th April 2001: 
c       Only the directions with the array mask=0 are calculated.
c       Directions not calculated are set to -1.0 in dzmodel and in
c       all output files. Additional disk file dust.mask output on
c       unit 11 containing mask array.
c
c       for further changes see file changes.txt
c
c----------------------------------------------------------------------------

	implicit real*8 (a-h, o-z)
        PARAMETER (NSIZE=980000)
        parameter (nsizeu=50000)
c        PARAMETER (NSIZE=1960000)
c        parameter (nsizeu=30000)
        DOUBLE PRECISION DZMODEL(NSIZE)
        double precision u(nsizeu),fa(nsizeu)
        double precision fv(nsizeu,256)
        double precision uct(256),ust(256),ucp(256),usp(256)
        double precision dustplo(2000,2000)
        double precision ru(100),zu(100)
        integer mask(nsize),maskcalc(nsize)
        common/incli/theta_in,tau0
*        common/mincom/dz,dez,npoint,nprm
	common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar
c	common/wei/wt(50),wp(50)
	common/f/fs0,fs1,fsdum,fsurad,dsb1,dsb2
	common/atte/a0,eta0,rtruncate,sharp,athresh,rtrun,sha
	common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
	common/bcons/etabmax,rin,ibulge
        common/ferrers/etaf,af,bf,cf,thetaf,ctf,stf,iferrers
        common/ellexp/eta2,zs2,hsa2,hsb2,thetae,cte,ste,
     >  hda2,hdb2,zd2,tau2,a2,iellexp
	common/g/g,g2,cosalp
	common/integ/k,m,n,km
	common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  eta1arm,a1arm,narm,rup,fillarm
	common/translate/pix2kpc
        common/inout/nx,ny,pixsize
        common/ur/u,fa,fv,ures,rsizeu,zsizeu,albedo,iscat
        common/ur2/ru,zu,nru,nzu,iurad,nru_used,nzu_used
        common/urang/icang,nang,uct,ust,ucp,usp
        common/photom/away,bfactor,bin,ifilter
        common/config/izero,idust,itau,istars,imask,iall,
     >  shiftx,shifty,ifswitch,ifinns1,ifinns2,ifinnd1,ifinnd2

c        write(6,1947)albedo,g
c 1947   format(1h ,'DUST.f:  albedo,g = ',2f8.2)
c


c
        pxsize=pixsize
c
        if(idust .eq. 2 .and. nx .gt. 2000)write(6,2698)
 2698   format(1h ,'warning: nx must be less than 2000 for idust=2')
        if(idust .eq. 2 .and. ny .gt. 2000)write(6,2699)
 2699   format(1h ,'warning: ny must be less than 2000 for idust=2')

c        write(6,890)theta_in
c 890   format('theta_in =',f10.4)
c        write(6,891)eta0
c 891    format('eta0 =',f10.4)
c        write(6,892)tau0
c 892    format('tau =',f10.4)
c        write(6,893)zs
c 893    format('zs =',f10.4)
c        write(6,894)zd
c 894    format('zd =',f10.4)
c        write(6,895)hs
c 895    format('hs =',f10.4)
c        write(6,896)hd
c 896    format('hd =',f10.4)
c        write(6,897)ellipt
c 897    format('ellipt =',f10.4)
c        write(6,898)reff
c 898    format('reff =',f10.4)
c        write(6,899)etab
c 899    format('etab =',f10.4)
c        write(6,794)eta1
c 794    format('eta1 =',f10.4)
c        write(6,795)zs1
c 795    format('zs1 =',f10.4)
c        write(6,796)hs1
c 796    format('hs1 =',f10.4)
c        write(6,7961)zs1t
c 7961   format('zs1t =',f10.4)
c        write(6,797)tau1
c 797    format('tau1 =',f10.4) 
c        write(6,798)zd1
c 798    format('zd1 =',f10.4)
c        write(6,799)hd1
c 799    format('hd1 =',f10.4)
c        write(6,7991)zd1t
c 7991   format('zd1t =',f10.4)
c        write(6,599)etaf
c 599    format('etaf =',f10.4)
c        write(6,598)af
c 598    format('af =',f10.4)
c        write(6,597)bf
c 597    format('bf =',f10.4)
c        write(6,596)cf
c 596    format('cf =',f10.4)
c        write(6,5599)eta2
c 5599   format('eta2 =',f10.4)
c        write(6,5598)hsa2
c 5598   format('hsa2 =',f10.4)
c        write(6,5597)hsb2
c 5597   format('hsb2 =',f10.4)
c        write(6,5596)zs2
c 5596   format('zs2 =',f10.4)
c


c	theta is the inclination angle.
c       For the definitions of eta0, tau0, zs, zd, hs, and hd 
c	see eqs. 12 and 13 of Kylafis & Bahcall (1987).
c	ellipt is the ellipticity of the bulge, reff is the effective
c	radius of the bulge and etab is the central luminosity
c	density of the bulge (Hubble profile).

c	open (unit=1, file='dust.in')
c	open (unit=2, file='config.in')
c	open (unit=3, file='photometry.in')

c       constrain the scattering in the energy density calculation
c       to be indentical with the scattering in the projected image
c       calculation
c
        iscat=izero

	if(idust .ne. 0)open (unit=4, file='dust.plo')
	if(istars .ne. 0)open (unit=7, file='dust.star')
	if(itau .ne. 0)open (unit=8, file='dust.tau')
        if(imask .ne. 0)open (unit=11,file='dust.mask')

       ntot=nx*ny
       if(nsize .lt. ntot)then
       write(6,2762)
 2762  format(1h ,'DUST - aborting - set nsize higher')
       goto 999
       endif  


c	Reading the parameter file 'dust.in'
c	read(1, *) theta_in
c        read(1, *) eta0
c        read(1, *) tau0
c        read(1, *) zs
c        read(1, *) zd
c        read(1, *) hs
c        read(1, *) hd
c	read(1, *) ellipt
c	read(1, *) reff
c	read(1, *) etab
c	read(1, *) eta1
c	read(1, *) zs1
c	read(1, *) hs1

c      do 27765 i=1,100
      do 27765 i=1,nsize
      dzmodel(i)=0.0
27765 continue
c
	ellipt2= ellipt* ellipt

c	Transformation to Model system
	pixsize2=pxsize*pxsize
	bin2=bin*bin
	pix2kpc=4.84d-3*away*pxsize

        write(6,343) pix2kpc
 343    format(1h , 'pix2kpc = ', f10.2)
	hd=hd/pix2kpc
	zd=zd/pix2kpc
        hd1=hd1/pix2kpc
        zd1=zd1/pix2kpc
        hdin=hdin/pix2kpc
        zdin = zdin/pix2kpc
        hdsolar = hdsolar/pix2kpc
        zdsolar = zdsolar/pix2kpc
        hd1in=hd1in/pix2kpc
        zd1in = zd1in/pix2kpc
        hd1solar = hd1solar/pix2kpc
        zd1solar = zd1solar/pix2kpc
	hs=hs/pix2kpc
	zs=zs/pix2kpc
	hs1=hs1/pix2kpc
	zs1=zs1/pix2kpc
        hsin=hsin/pix2kpc
        zsin = zsin/pix2kpc
        hssolar = hssolar/pix2kpc
        zssolar = zssolar/pix2kpc
        hs1in=hs1in/pix2kpc
        zs1in = zs1in/pix2kpc
        hs1solar = hs1solar/pix2kpc
        zs1solar = zs1solar/pix2kpc
        af=af/pix2kpc
        bf=bf/pix2kpc
        cf=cf/pix2kpc
        hsa2=hsa2/pix2kpc
        hsb2=hsb2/pix2kpc
        zs2=zs2/pix2kpc
        
        hda2=hda2/pix2kpc
        hdb2=hdb2/pix2kpc
        zd2=zd2/pix2kpc

	reff=reff/pix2kpc
        rtrun=rtruncate/pix2kpc
        rtrun1=rtruncate1/pix2kpc
c        write(6,3699)rtrun,rtruncate
c 3699     format(1h ,'***** rtrun,rtruncate = ',2e12.3)
        sha=sharp/pix2kpc
        sha1=sharp1/pix2kpc
crjt    convert taus to edge on values. They are input
crjt    as face-on values
c        if(iurad .eq. 2)ru=ru/pix2kpc
c        if(iurad .eq. 2)zu=zu/pix2kpc
        if(iurad .eq. 2)then
         do 790 iru=1,nru
          ru(iru)=ru(iru)/pix2kpc
 790     continue
         do 791 izu=1,nzu
          zu(izu)=zu(izu)/pix2kpc
 791   continue
        endif

	tau0=tau0*hd/zd
        tau1=tau1*hd1/zd1
c       tau for elliptical disk is normalised face on
        tau2=tau2
c       convert face on magitudes per sq arcsec to 
c       volume emissivities (per cubic pixel)
c	eta0=(10.**(-eta0/2.5))*pixsize2/(2.d0*hs*bfactor)/bin2
c        eta1=(10.**(-eta1/2.5))*pixsize2/(2.d0*hs1*bfactor)/bin2
c        etab=(10.**(-etab/2.5))*pixsize2/(5.12d0*reff*bfactor)/bin2
c	etaf=(10.**(-etaf/2.5))*pixsize2/(16.d0*cf*bfactor/15.d0)/bin2

	eta0=(10.**(-eta0/2.5))*pixsize2
c edgeon        eta0=eta0*bin2/(2.d0*hs*bfactor)
        eta0=eta0*bin2/(2.d0*zs*bfactor)

        eta1=(10.**(-eta1/2.5))*pixsize2
c edgeon       eta1=eta1*bin2/(2.d0*hs1*bfactor)
        eta1=eta1*bin2/(2.d0*zs1*bfactor)

        eta2=(10.**(-eta2/2.5))*pixsize2
        eta2=eta2*bin2/(2.d0*zs2*bfactor)

        etab=(10.**(-etab/2.5))*pixsize2
c edgeon       etab=etab*bin2/(5.12d0*reff*bfactor)
c face on

         if(ibulge .gt. 3 .or. ibulge .lt. 0)then
          write(6,2689)ibulge
 2689     format(1h ,'**** INPUT error - ibulge out of bounds')
          goto 999
         endif        



        if(ibulge .eq. 1)etab=etab*bin2/(2.0d0*reff*bfactor*ellipt)
        if(ibulge .eq. 2)then
c        etab=etab*bin2/(5.12d0*reff*bfactor*ellipt)
         factorb=2.0*3.14159/7.66925
         factorb=sqrt(factorb)
         etab=etab*bin2/(factorb*reff*bfactor*ellipt)
        endif
        if(ibulge .eq. 3)then
         if(nsersic .eq. 1)bsersic=1.67835d0
         if(nsersic .eq. 2)bsersic=3.67206d0
         if(nsersic .eq. 3)bsersic=5.67017d0
         if(nsersic .eq. 4)bsersic=7.66925d0
         if(nsersic .eq. 5)bsersic=9.66872d0
         if(nsersic .eq. 6)bsersic=11.6684d0
         if(nsersic .eq. 7)bsersic=13.6681d0
         if(nsersic .eq. 8)bsersic=15.6679d0
         if(nsersic .eq. 9)bsersic=17.6678d0
         if(nsersic .eq. 10)bsersic=19.6677d0
         factorb=2.0*3.14159/bsersic
         factorb=sqrt(factorb)
         etab=etab*bin2/(factorb*reff*bfactor*ellipt)
        endif


c        if(ibulge .eq. 1)etab=etab*bin2*ellipt/(2.0d0*reff*bfactor)
c        if(ibulge .eq. 2)etab=etab*bin2*ellipt/(5.12d0*reff*bfactor)
c
	etaf=(10.**(-etaf/2.5))*pixsize2
        etaf=etaf*bin2*15.0d0/(16.d0*cf*bfactor)




c note that in Nick's original program the normalisation
c of the exponential disks and bulge was done through
c line integrations side-on. The definition is self-consistent
c with the code in function eta for side-on integrations.
c The ferrers bar is normalised face on.
c
c       The anisotropy parameter g (see eq. 5) and the albedos (see
c       definition below eq. 6) are assumed to be correctly specified
c       as input parameters passed through the common blocks g and ur,
c       for the filter specified.
c	The numbers 1,2,3,4,5,6,7,8,9,10,11,12,13,14 15 are for the 
c       uv09,uv15,uv28,uv31,B,V,I,J,K,uv22,uv13,uv16,uv20,uv25 u filters 
c       respectively.
c       uv09 912 A
c       uv15 1500 A
c       uv28 2800 A
c       uv31 3150 A
c       uv22 2200 A
c       uv13 3150 A
c       uv16 1650 A
c       uv20 2000 A
c       uv25 2500 A

        ufactor=(4.0d0*3.1415927d0)/((2.997925*1.0d8)*(1.0d26))
        pxs=(pxsize/3600.0)*3.1415927d0/180.0d0
        ufactor=ufactor/(pxs*pxs)
c        write(6,289)ufactor
c289     format(1h ,'ufactor = ',15e12.6)
c       at this point ufactor should be 1.7833613e-23 for 1" pixels 
c
c       the Jy per zero magnitude in the UV were calculated as
c       monocromatic values from the Vega spectrum of 
c       Dreiling & Bell (ApJ 241 736). The value at 912 A is however
c       arbitrarily set to the value as 1500 A since Vega has
c       almost no flux at 912 A. However, provided the magnitudes
c       per sq arcsec given as the input parameter to this program
c       is calculated from the SFR using this same Jy per zero
c       magnitude factor, the calculated energy densities will
c       be correct and indeoendent of the actual value of
c       the conversion factor in the UV. In the optical, where
c       energy densities are related to observations rather than
c       a model parameter SFR, the absolute value of the conversion
c       factor becomes meaningful.
c
	  if(ifilter.eq.1) then
c       uv09 912
c	  albedo=0.28
c	  g=0.57
c       Jy for zero magnitude
        szero=460.0d0
        c1=1.0d0/(bfactor*szero)
        ufactor=ufactor/c1
	  endif
	  if(ifilter.eq.2) then
c       uv15 1500
c	  albedo=0.39
c	  g=0.60
c       Jy for zero magnitude
        szero=460.0d0
        c1=1.0d0/(bfactor*szero)
        ufactor=ufactor/c1
	  endif
	  if(ifilter.eq.3) then
c       uv28 2800
c	  albedo=0.54
c	  g=0.49
c       Jy for zero magnitude

        szero=954.0d0
        c1=1.0d0/(bfactor*szero)
        ufactor=ufactor/c1
	  endif
	  if(ifilter.eq.4) then
c       uv31 3150
c	  albedo=0.54
c	  g=0.50
c       Jy for zero magnitude
        szero=1160.0d0
        c1=1.0d0/(bfactor*szero)
        ufactor=ufactor/c1
	  endif
	if(ifilter.eq.5) then
c       B-band 4350
c	albedo=0.65 from Nick
c        albedo=0.53 from pop00
c	g=0.58 from Nick
c        g=0.49 from pop00
c       szero is Jy for zero magnitude
        szero=4063.0d0
c       c1 is counts per Jy, if bfactor is counts for zero magnitude
        c1=1.0d0/(bfactor*szero)
c       st this point ufactor is multiplitive factor to obtain 
c       spectral energy densities in Joules per cubic metres per Hz
        ufactor=ufactor/c1
	endif
	  if(ifilter.eq.6) then
c       V-band 5550
c	  albedo=0.60 from Nick
c        albedo=0.52 from pop00
c	  g=0.50 ; from Nick
c        g=0.45 from pop00
c       Jy for zero magnitude
        szero=3636.0d0
        c1=1.0d0/(bfactor*szero)
        ufactor=ufactor/c1
	  endif
	if(ifilter.eq.7) then
c       I-band 8800
c        albedo=0.49
c	g=0.34 from Nick
c        g=0.36 from pop00
c       Jy for zero magnitude
        szero=2416.0d0
        c1=1.0d0/(bfactor*szero)
        ufactor=ufactor/c1
	endif
	  if(ifilter.eq.8) then
c       J-band 12500
c	  albedo=0.37
c	  g=0.16 from Nick
c        g=0.15 from pop00
c       Jy for zero magnitude
        szero=1589.0d0 
        c1=1.0d0/(bfactor*szero)
        ufactor=ufactor/c1
	  endif
	if(ifilter.eq.9) then
c
c       K-band 22000
c	albedo=0.25 from Nick
c        albedo=0.16 from pop00
c	g=0.0 from nick
c        g=0.04 from pop00
c       Jy for zero magnitude
        szero=640.0d0 
        c1=1.0d0/(bfactor*szero)
        ufactor=ufactor/c1
	endif
	if(ifilter.eq.10) then
c       uv22 220 nm
c	   albedo=0.43
c	   g=0.49
c       Jy for zero magnitude
           szero=728.0d0 
           c1=1.0d0/(bfactor*szero)
           ufactor=ufactor/c1
	endif
	if(ifilter.eq.11) then
c       uv13 135 nm
c	   albedo=0.38
c	   g=0.62
c       Jy for zero magnitude
           szero=249.0d0 
           c1=1.0d0/(bfactor*szero)
           ufactor=ufactor/c1
	endif
	if(ifilter.eq.12) then
c       uv16 165 nm
c	   albedo=0.43
c	   g=0.61
c       Jy for zero magnitude
           szero=554.0d0 
           c1=1.0d0/(bfactor*szero)
           ufactor=ufactor/c1
	endif
	if(ifilter.eq.13) then
c       uv20 200 nm
c	   albedo=0.43
c	   g=0.53
c       Jy for zero magnitude
           szero=710.0d0 
           c1=1.0d0/(bfactor*szero)
           ufactor=ufactor/c1
	endif
	if(ifilter.eq.14) then
c       uv25 250 nm
c	   albedo=0.52
c	   g=0.48
c       Jy for zero magnitude
           szero=759.0d0 
           c1=1.0d0/(bfactor*szero)
           ufactor=ufactor/c1
	endif
	if(ifilter.eq.15) then
c       U-band 365 nm
c       Jy for zero magnitude
c          Longair
           szero=188.9
           c1=1.0d0/(bfactor*szero)
           ufactor=ufactor/c1
	endif
c
c       convert ufactor such that spectral energy densities
c       are output in eV per cubic centimetres per Hz (so that
c       the energy densities multiplied by frequency should be of
c       order unity for a galaxy like the Milky Way). 
c
        ufactor=ufactor/1.60219d-13
c        write(6,289)ufactor
c2899    format(1h ,'ufactor = ',15e12.6)
c
	g2=g*g

c       For computational reasons we will take the galaxy to be
c       situated symmetrically inside a cylinder of radius rhomax 
c	and half height d.

	nx2= nx/ 2
	ny2= ny/ 2
	rxx= dfloat(nx2)
	ryy= dfloat(ny)
c
c  find maximum position of unmasked data for ryy
c

      iimax=nx*ny
      ylim=0.0
        do  7828 ii=1,iimax
	iy= int(ii/ nx)+ 1
	iz= ii- int(ii/ nx)* nx

	if(mod(ii,nx) .eq. 0) then
	iy= ii/ nx
	iz= ii- nx*(iy-1)
	endif
        yyy=dfloat(iy)
         if(mask(ii) .eq. 0)then
            if(ylim .le. yyy)then
            ylim=yyy+1
            endif 
         endif
 7828	continue
      ryy=ylim
c      write(6,9329)ryy,ny
c 9329 format(1h ,'ryy = ',e12.3,i6)
c
c	rhomax1= 3.0* max(hd, hs)
        rhomax1= 5.0* max(hd, hs)
c       the first disk is taken to be the thicker of
c       the two disks
c	d1= 7.0* max(zs, reff)
        d1= 15.0* max(zs, reff)
	rhomax= max(rhomax1, ryy)
c rjt 11/01 this is only correct for an edge on system.
c for a face on system rhomax1 should be max(rhomax1,ryy,rxx)
c this might lead to rays being calculated which have yasked,zasked 
c outside the cylinder for face on systems. Hence insert the following
c to ensure the cylinder radius is greater than the major and minor axis
c data radius.
c 

        rhomax=max(rhomax,rxx)
c
	d= max(d1,rxx)
c rjt 11/01 d=max(d1,rxx) is only correct for an edge on system.
c for a face on system d should have nothing to do with rxx and we might
c be calculating an inefficiently long pathlength.
c
c        write(6,388)nx,ny,hd,hs,zs,reff,rxx,ryy,rhomax,d
c 388    format(1h ,
c     >  'DUST: nx,ny,hd,hs,zs,reff,rxx,ryy,rhomax,d = ',2i4,8e12.4)

c	Definition of the step size for the integration in the bulge
c	region.
	dsb1= reff*ellipt/100.0
	dsb2= reff*ellipt/10.0
c        dsb1= reff*ellipt/10.0
c        dsb2= reff*ellipt/10.0

c	CALCULATION OF cos(alpha) FOR WHICH THE HENYEY-GREENSTEIN
c	PHASE FUNCTION IS EQUAL TO 1/2 ITS (MAXIMUM+MINIMUM) VALUE.
c	(SEE EQ. 4).
c	This  is done for later use, when we do integrals of the form
c	of eq. 2.

	cosalp= 0.0d0
	if(g .eq. 0.0d0) goto 9
	fun1=1.0d0+g2
	fun2=(1.0d0-g2)*(1.0d0-g2)
	fun3=(1.0d0+3.0d0*g2)**(2.0d0/3.0d0)
	cosalp=(fun1-fun2/fun3)/(2.0d0*g)
9	continue

c	For the definition of the following quantity see section 3d.
        a0=0.0d0
      	if(idisk1 .eq. 1)a0= tau0/ (2.0d0* hd)
      	if(idisk1 .eq. 3)a0= tau0/ (2.0d0* hd)

crjt
        a1=0.0d0
        if(idisk2 .eq. 1)a1= tau1/ (2.0d0* hd1)
        if(idisk2 .eq. 3)a1= tau1/ (2.0d0* hd1)
c       tau for elliptical disk is normalised face on
        a2=0.0d0
        if(iellexp .eq.1)a2= tau2/ (2.0d0* zd2)
c
c	For the R^1/4 law, we give a constant value to the inner
c	region in a sphere of radius 3 pixels. This value is
c	etabmax.
c 	rin= 3.d0/ reff
c        write(6,3629)rin
 3629   format(1h ,'rin = ',e12.4)
c
c revised by rjt to prevent central spike
c
c       For the R^1/4 law, we give a constant value to the inner
c	region in a sphere of radius reff/3.0 pixels.
c        rin=reff/3.0d0
        rin=0.4d0
        acaprin=rin/reff
        rin=acaprin
        write(6,3629)rin
	etabmax= etab* dexp(-7.67*rin**0.25)/ rin**0.875
c note that etabmax is never used in the program
c concerning the3 pixel limit see also zone2= 3.0+ ds0/reff in FSTEP

         call plot_etabulger
         call plot_etabulgez
         call plot_etadiskr
         call plot_etadiskz
         call plot_etadisk1r
         call plot_etadisk1z
         call plot_taudiskr
         call plot_taudiskz


c       Initialise parameters for spiral arms in second emissivity
c       and dust disk (this only has an effect for idisk2=2). 
c       rup is needed in definition of spiral arms in SETARM,
c       which must be called after 2nd disk parameters
c       eta1 and tau1 have been initialised
c       in main program (above)
        if(idisk2 .eq. 2)then
c       armset must only be called once
        call armset
        endif
c
c	At the center of the galaxy define a Cartesian coordinate
c	system such that the z-axis is perpendicular to the plane
c	of the galaxy and our line  of sight is in the xz-plane.
c	Thus, the spherical coordinates of our line of sight are 
c	theta_in (the  inclination angle of the galaxy) and phi=0.
c	We restrict ourselves to 0 < cos(theta) < 1.  
c	Theta = 90 degrees means edge-on.

	conver=3.14159d0/180.0d0
	theta=theta_in*conver
	ct=dcos(theta)
	st=dsin(theta)
	phi=0.0d0*conver
	cp=dcos(phi)
	sp=dsin(phi)
        thetaf=thetaf*conver
        thetae=thetae*conver
        ctf=dcos(thetaf)
        stf=dsin(thetaf)
        cte=dcos(thetae)
        ste=dsin(thetae)
c        write(6,*)cte,ste

c       PREPARATION FOR CALCULATION OF I sub 1 (see eq. 10)

c       The integral over all solid angles in equation 10 is done with
c       the Gauss method.  k+m directions are used for the integration over
c       mu = cos(theta) and n directions for the phi integration.
c       The larger the values of k,m,n, the more accurate the calculation
c       is.

	k= 3
	m= 5
	n= 8

c	k= 4
c	m= 6
c	n= 10

c	k= 6
c	m= 10
c	n= 16

c       For g>0, k is for the significant part of  the mu integration, and m
c       for the insignificant. For g<0, it is the opposite.
c       k+m, as well as n, should be less than 50
c       n should always be even

	km=k+m

	call direct(ct,st,cp,sp)

c        testsum=0.0
c           write(6,2997)
c 2997      format(1h ,'i,j,phase(i,j),wt(j),wp(i)')
c        do 833 j=1,km
c        do 832 i=1,n
c           write(6,2998)i,j,phase(i,j),wt(j),wp(i)
c 2998      format(1h ,2i6,3f10.2)
c           testsum=testsum+phase(i,j)*wt(j)*wp(i)
c 832       continue
c 833       continue
c        testsum=testsum/4.0d0/3.141592654d0
c        write(6,2996)testsum
c 2996   format(1h ,'testsum = ',f10.3)


c	Our line of sight intersects the projection-plane of the galaxy 
c	at a point (0, yasked, zasked). The coordinates of this point
c	are transformed in the galaxy system (described above), and
c	the integrations are done from this point to the boudaries
c	of the cylinder that the galaxy is in.
c	The projection-plane (the image), consists of nx times ny pixels.

	iimax= nx* ny
c        write(6,9787)iimax
c 9787   format('imax =',i8)
c
c       initialize array of calculated directions

        do 2999 ii=1,iimax
        maskcalc(ii)=mask(ii)
2999    continue
c
c       dust density threshold for calculation of scattering
        zthresh=3.5*zd
        athresh=a(0.0d0,0.0d0,zthresh)
c
	do 6000 ii= 1, iimax
c
        if(mask(ii) .ne. 0)goto 6001
c        write(6,6002)ii
c 6002	format('calculating for ii = ',i6)
c
	iy= int(ii/ nx)+ 1
	iz= ii- int(ii/ nx)* nx

	if(mod(ii,nx) .eq. 0) then
	iy= ii/ nx
        if(ii .eq. 1)write(6,288)
288	format(' ')
c        write(6,289)iy
c289	format(1h ,'calculating row ',i8)
	iz= ii- nx*(iy-1)
	endif

	iynew= iy
	if(iall.eq.1) iynew= iy- ny2
	iznew= iz- nx2

c rjt Jan02 for shiftx and shifty zero make sure rays pass
c through centre of galaxy as appropriate for value of iall 
c this avoids need to input different shifts when nx and/or ny
c change between even and odd. 
c ie we assume that for shifty=0 the lefthand pixel
c on the input imge passes through yasked=0 if iall=0 (half image
c calculated). For iall=1 this is only the case for ny odd
c similarly for x coordinate
crjt	yasked= shifty+ dfloat(iynew)- 0.5
crjt	zasked= shiftx+ dfloat(iznew)- 0.5
        if(iall .ne. 0 .and. iall .ne. 1)then 
           write(6,*) 'abortin - iall must be 0 or 1'
           goto 999
        endif
        if(iall .eq. 0)yasked= shifty+ dfloat(iynew)-1
        if(iall .eq. 1)then 
           if(2*(ny/2) .ne. ny)yasked= shifty+ dfloat(iynew)-1.0
           if(2*(ny/2) .eq. ny)yasked= shifty+ dfloat(iynew)-0.5
        endif
        if(2*(nx/2) .ne. nx)zasked= shiftx+ dfloat(iznew)- 1.0
        if(2*(nx/2) .eq. nx)zasked= shiftx+ dfloat(iznew)- 0.5
c        write(6,6005)ii,iy,iz,yasked,zasked
c 6005	format('ii,iy,iz,yasked,zasked = ',3i6,' ',f6.3,' ',f6.3)
c

c	For a given yasked and zasked in the projection-plane (the image),
c	the surface brightness is calculated and a model image is
c	created.

	x= -zasked*ct
	y= yasked
	z= zasked*st

c        trap for face on case
c        rtest=dsqrt(y*y+x*x)
c        if(rtest .gt. rhomax)write(6,2399)rtest
c 2399   format(1h ,'DUST: ray outside cylinder: r = ',e12.4)

c        if(sqrt(x*x+y*y+z*z) .lt. 2.)write(6,11)x,y,z,yasked,zasked
c 11	format('x y z       =',5f6.3)

	b0=0.0
	b1=0.0

c	calculation of intensity I sub 0 (called b0) (see eq. 9), 
c	unattenuated star intensity (i.e., as if there were no dust)
c	(called bstars) and optical depth along the line of sight 
c	tau(R,z) (called tauRz) for various values of z.
c	Eq. 16 gives tau(R,z) for edge-on galaxies of infinite extend.
c	Here R is equal to yasked and z equal to zasked.

c       fs0 is a small number that determines the integration step
c       (see subroutine bst)
c       edge on
c	fs0=0.15 
c       face on
c        fsdum=0.5

c        fs0=0.05
c        fsdum=0.05

        if(dabs(ct) .gt. 0.173648)fs0=fsdum 
        if(ifswitch .eq. 1)fs0=fsdum
        if(ifswitch .eq. 2)fs0=fsdum
        if(ifswitch .eq. 4)fs0=fsdum
c       The subroutines bound0 and  boundf determine the limits of
c       integration along the line of sight (i.e., the point of the
c       galaxy farthest from the observer and the point closest to
c       the observer)

c        write(6,2046)
c 2046   format(1h ,'call bound0')
	call bound0(x,y,z,x0,y0,z0,si,ct,st,cp,sp,ierr)
        ierr1=ierr
c        write(6,2049)
c2049   format(1h ,'end call bound0')
c        if(y .eq. 0 .and. abs(x) .lt. 1.0 .and. abs(z) .lt. 1.0)
c     >  write(6,*)x,y,z,x0,y0,z0,si
c following is equivalent condition for 90 or 0 inclination:
c        if(y0 .eq. 0 .and. abs(z0) .lt. 1.0)write(6,*)x,y,z,x0,y0,z0,si

	call boundf(x,y,z,xf,yf,zf,sf,ct,st,cp,sp,ierr)
c        if(y .eq. 0 .and. abs(x) .lt. 1.0 .and. abs(z) .lt. 1.0)
c     >  write(6,*)x,y,z,xf,yf,zf,sf
c following is equivalent condition for 90 or 0 inclination:
c        if(y0 .eq. 0 .and. abs(z0) .lt. 1.0)write(6,*)x,y,z,xf,yf,zf,sf
        ierr2=ierr
c        write(6,2047)
c2047   format(1h ,'end call boundf')
c        write(6,*)x,y,z,xf,yf,zf,sf
c
        ierr=ierr1+ierr2
        iextracy=0
        if(ierr .ne. 0)then
c        write(6,2995)yasked,zasked
c 2995   format(1h ,'DUST: extracylinder ray: yasked zasked = ',2e12.4)
        maskcalc(ii)=1
        iextracy=1
        goto 6001
        endif

	if (dabs(z) .le. d) goto 290
	sf=dsqrt((xf-x0)*(xf-x0)+(yf-y0)*(yf-y0)+(zf-z0)*(zf-z0))
290	continue
	stot=si+sf

c        write(6,474)stot,x0,y0,z0,ct,st,cp,sp
c474	format(1h ,'calling br0',8f8.3)
	b0=br0(stot,x0,y0,z0,ct,st,cp,sp)
c        write(6,478)
c 478	format(1h ,'left br0')
c
c        if(b0 .le. 0.0)write(6,2078)iy,iz
c2078    format(1h ,'warnings: zero brightnes calculated: iy,iz = ',2i6)


c        write(6,379)b0
c 379	format(1h ,'b0',f8.3)
c	If istars=0, NO calculation of pure starlight (as if there
c	were no dust) takes place.
	if(istars.ne.0) bstars= bst(stot,x0,y0,z0,ct,st,cp,sp)
c	If itau=0, NO optical depth calculation takes place.
	if(itau.ne.0) tauRz= tau(stot,x0,y0,z0,ct,st,cp,sp)
c	if izero=0, NO scattering is considered.
c        write(6,*)izero
	if (izero .eq. 0) goto 1800

c	This makes the assumption that scattering is not 
c	significant for z > 6*zd .
c        write(6,2890)z,zd
c 2890	format(1h ,'z,zd = ',2f10.2)
c	if (dabs(z/zd) .gt. 6.0) goto 1800
c
c rjt	if (a(x,y,z) .lt. athresh) goto 1800
c
c	CALCULATION OF INTENSITY I sub 1 (called b1) (see eq. 10)

c	In order to calculate I sub 1 one must calculate an integral
c	over all solid angles AT EVERY POINT along the line of sight 
c	in the galaxy.  This takes a lot of computer time and thats 
c	why I take fs0 and fs1 = 0.2.  The smaller these numbers are,
c	the more accurate the calculation is.

c        fs0=0.15
c        fs1=0.2
c        fsdum=0.5
        

        if(dabs(ct) .gt. 0.173648)then
        fs0=fsdum
        fs1=fsdum
        endif
        if(ifswitch .eq. 1 .or. ifswitch .eq. 2)then
           fs0=fsdum
           fs1=fsdum
        endif
        if(ifswitch .eq. 4)then
           fs0=fsdum
           fs1=fsdum
        endif


c
c        write(6,377)x0,y0,z0,ct,st,cp,sp
c 377	format(1h ,'calling br1 with x0,y0,z0,ct,st,cp,sp = ',7e12.4)
	b1= br1(stot,x0,y0,z0,ct,st,cp,sp)
c        write(6,378)
c 378	format(1h ,'exited br1')
1800	continue
c
c	CALCULATION OF I(R,z) (called btotal) (see eq. 20)
c
        call setzero(b0,b1,btotal)
c        write(6,3279)b0,b1,btotal
c 3279   format(1h ,'b0,b1,btotal ',3e15.6)
c
6001   continue
c
c outnocalc is value given to output data for uncalculated masked directions
c these are specified either by the input MASK array or the internal
c array MASKCALC. (MASKCALC are all directions on the data map not
c enclosed by the cylinder).
c
        outnocalc=-1.0 
        if(iextracy .eq. 1)then
           outnocalc=0.0
           iextracy=0
        endif 
        out1=outnocalc
        out2=outnocalc
        out3=outnocalc
        out4=mask(ii)
        if(mask(ii) .eq. 0 .and. maskcalc(ii) .eq.0)then
        out1=btotal
        out2=bstars
        out3=tauRz
        out4=mask(ii)
        endif
        if(abs(out1) .le. 1.0D-50)out1=0.0
        if(abs(out2) .le. 1.0D-50)out2=0.0
        if(abs(out3) .le. 1.0D-50)out3=0.0
	if(idust .eq. 1)write(4, *) out1
        if(idust .eq. 2)dustplo(iy,iz)=out1
c
        dzmodel(ii)=out1
*        write(6,*)ii,dzmodel(ii)
	if(istars .ne. 0)write(7, *) out2
	if(itau .ne. 0)write(8, *) out3
        if(imask .ne. 0)write(11,*) out4
c	write(4, *) btotal
c	write(7, *) bstars
c	write(8, *) tauRz
c
6000	continue
c
c calculate energy densities
c 
        if(iurad .eq. 1 .or. iurad .eq. 2)call urad_cylinder(ufactor)
c
 999    continue
	close (unit=13)
        if(idust .eq. 2)then
         do 923 iy=1,ny
         write(4,*)(dustplo(iy,k),k=1,nx)
 923     continue
        endif
	if(idust .ne. 0)close (unit=4)
	if(istars .ne. 0)close (unit=7)
	if(itau .ne. 0)close (unit=8)
        if(imask .ne. 0)close (unit=11)

	return
	end

c------------------------
        function a(x,y,z)
c------------------------

	implicit real*8 (a-h, o-z)

c	attenuation coefficient due to dust (see eq. 13)

	common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar
	common/atte/a0,eta0,rtruncate,sharp,athresh,rtrun,sha
	common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
	common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  eta1arm,a1arm,narm,rup,fillarm
        common/ellexp/eta2,zs2,hsa2,hsb2,thetae,cte,ste,
     >  hda2,hdb2,zd2,tau2,a2,iellexp
        common/config/izero,idust,itau,istars,imask,iall,
     >  shiftx,shifty,ifswitch,ifinns1,ifinns2,ifinnd1,ifinnd2

	Rho=dsqrt(x*x+y*y)
c       1st disk and all other structure except 2nd disk        
        sigtruncate=sha*rtrun
        rplane=sqrt(x*x+y*y+z*z)
        truncate=0.5*(1.0-derf((rplane-rtrun)/sigtruncate))
        sigtruncate1=sha1*rtrun1
        truncate1=0.5*(1.0-derf((rplane-rtrun1)/sigtruncate1))


c flared disk for the first (old) dust disk
         zd01 = zd
         if (zdin .gt. zd) then
          flared01 = dlog10((zdsolar - zd)/(zdin - zd))
         else 
          flared01 = 0.
         end if
         flared02 = dlog10(hdsolar/hdin)
         flared00 = flared01/flared02
         zd01 = zd + (zdin - zd) * (rho/hdin)**flared00

c flared disk for the second (young) dust disk
         zd11 = zd1
         if (zd1in .gt. zd1) then
          flared11 = dlog10((zd1solar - zd1)/(zd1in - zd1))
         else 
          flared11 = 0.
         end if
         flared12 = dlog10(hd1solar/hd1in)
         flared11 = flared11/flared12
         zd11 = zd1 + (zd1in - zd1) * (rho/hd1in)**flared11

        if(idisk1 .eq. 1) then
         if (Rho .le. hdin) then
          if (ifinnd1 .eq. 1) then
           adum=a0* (zd/zd01)*
     >     dexp(-hdin/hd- dabs(z)/zd01)
          endif
          if (ifinnd1 .eq. 2) then
            adum=a0*(zd/zd01)*(1./hdin)*rho*
     >      dexp(-hdin/hd- dabs(z)/zd01)
          endif
          if (ifinnd1 .eq. 3) then
            adum=(a0/2.)*(zd/zd01)*(rho/hdin + 1.) * 
     >      dexp(-hdin/hd- dabs(z)/zd01)
          endif
         endif
         if (Rho .gt. hdin) then
          adum=a0*(zd/zd01)*dexp(-rho/hd- dabs(z)/zd01)
         endif
        endif

c  sech2 vertical distribution        
        if(idisk1 .eq. 3) then
         sechd01=1.d0/cosh(dabs(z)/zd01)
   	 sechd02=sechd01*sechd01
         if (Rho .le. hdin) then
          if (ifinnd1 .eq. 1) then
           adum=a0*(zd/zd01)*dexp(-hdin/hd)*sechd02
          endif
          if (ifinnd1 .eq. 2) then
            adum=a0*(zd/zd01)*(1./hdin)*rho*
     >      dexp(-hdin/hd)*sechd02
          endif
          if (ifinnd1 .eq. 3) then
            adum=(a0/2.)*(zd/zd01)*(rho/hdin + 1.) * 
     >      dexp(-hdin/hd)*sechd02
          endif
         endif
         if (Rho .gt. hdin) then
          adum=a0*(zd/zd01)*dexp(-rho/hd)*sechd02
         endif
        endif

        adum=adum*truncate
        adum1 = adum

c
c     young dust disk        
c double exponential disk
        if(idisk2 .eq. 1) then
         if (Rho .le. hd1in) then 
          if (ifinnd2 .eq. 1) then
            adum=adum+truncate1*a1*
     >     (zd1/zd11)*dexp(-hd1in/hd1- dabs(z)/zd11)
          endif
          if (ifinnd2 .eq. 2) then
           adum=adum+truncate1*a1*
     >     (zd1/zd11)*(1./hd1in)*rho*dexp(-hd1in/hd1- dabs(z)/zd11)
          endif
          if (ifinnd2 .eq. 3) then 
           adum=adum+truncate1*(a1/2.)*
     >     (zd1/zd11)*(rho/hd1in+1.)*dexp(-hd1in/hd1- dabs(z)/zd11)
           endif
         endif
         if (Rho .gt. hd1in) then
          adum=adum+truncate1*a1*
     >    (zd1/zd11)*dexp(-rho/hd1- dabs(z)/zd11)
         endif
        endif

c spiral arm
        if(idisk2 .eq. 2) then
         adum=adum+truncate1*aarm(rho)*(zd1/zd11)*
     >   dexp(-dabs(z)/zd11)
        endif

c sech2 vertical distribution
        if(idisk2 .eq. 3) then
	 sechd11=1.d0/cosh(dabs(z)/zd11)
	 sechd12=sechd11*sechd11
         if (Rho .le. hd1in) then 
          if (ifinnd2 .eq. 1) then
           adum=adum+truncate1*a1*
     >    (zd1/zd11)*dexp(-hd1in/hd1)*sechd12
          endif
          if (ifinnd2 .eq. 2) then
           adum=adum+truncate1*a1*
     >     (zd1/zd11)*(1./hd1in)*rho*dexp(-hd1in/hd1)*sechd12
           endif
          if (ifinnd2 .eq. 3) then
           adum=adum+truncate1*(a1/2.)*
     >     (zd1/zd11)*(rho/hd1in+1.)*dexp(-hd1in/hd1)*sechd12
          endif
         endif
         if (Rho .gt. hd1in) then
          adum=adum+truncate1*a1*
     >    (zd1/zd11)*dexp(-rho/hd1)*sechd12
         endif
        endif

c       elliptical exponential disk
        if(iellexp .eq. 1)then
        ye=y*cte - x*ste
        xe=x*cte + y*ste
        ye=dabs(ye)
        xe=dabs(xe)
        rhoe2=(ye*ye)/(hda2*hda2)+(xe*xe)/(hdb2*hdb2)
        dummy=0.0d0 -dsqrt(rhoe2) -(dabs(z)/zd2)
        dum=max(-30.0d0,dummy)
        adum = adum + (a2 * truncate * dexp(dum))                
        endif     
        a=adum

	return
	end

c-------------------------------------------------------
        subroutine bound0(x,y,z,x0,y0,z0,si,ct,st,cp,sp,ierr)
c-------------------------------------------------------
c	It computes the lower limit of integration along a line of sight
c       Input:  x,y,z,ct,st,cp,sp,d,rhomax
c       Output: x0,y0,z0, si, ierr
c
	implicit real*8 (a-h, o-z)

	common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar
c
c        write(6,29)x,y,z,ct,st,cp,sp,rhomax
c 29     format(1h ,'bound0: x,y,z,ct,st,cp,sp,rhomax = ',8e12.4)

        ierr=0
c
c        if(st .lt. 0.0)write(6,9)
c 9	format('warning bound0 theta negative: theta = ',f8.3)
	ct1=ct
	if (ct .eq. 0.0d0) ct1=1.0d-6
	z0=-d*ct1/dabs(ct1)
	y0=y+(z0-z)*st*sp/ct1
	x0=x+(z0-z)*st*cp/ct1
	rho0=dsqrt(x0*x0+y0*y0)
c        write(6,299)x0,y0,z0,rho0
c 299    format(1h ,'x0,y0,z0,rho0 = ',4e12.4)
	si=(z-z0)/ct1
	if (rho0 .le. rhomax) goto 100
	alpha=st*st
c next line rjt 11/01 to deal with case ray parallel to
c axis of cylinder but outside radius of cylinder
        if(alpha .eq. 0.0) then
c        write(6,699)x,y,z,ct,st,cp,sp,rhomax
c699    format(1h ,'BOUND0 failed: x,y,z,ct,st,cp,sp,rhomax = ',8e12.4)
        ierr=1
        goto 100
        endif
	beta=2.0*st*(x*cp+y*sp)
	gamma=(x*x+y*y-rhomax*rhomax)
        test=beta*beta-4.0*alpha*gamma

        if(test .lt. 0.0)then
c       ray does not intersect cylinder because it started from outside
cc        write(6,599)x,y,z,ct,st,cp,sp,rhomax
cc 599   format(1h ,'BOUND0 aborting: x,y,z,ct,st,cp,sp,rhomax = ',8e12.4)
        ierr=1
        goto 100
        endif

cc        write(6,799)x,y,z,ct,st,cp,sp,rhomax
cc 799   format(1h ,'BOUND0 cont: x,y,z,ct,st,cp,sp,rhomax = ',8e12.4)

	si=(beta+dsqrt(beta*beta-4.0*alpha*gamma))/(2.0*alpha)
c        si=(beta+dsqrt(test))/(2.0*alpha)
	x0=x-si*st*cp
	y0=y-si*st*sp
	z0=z-si*ct1
100	continue
c        write(6,298)
c 298    format(1h ,'exiting bound0')
c        write(6,378)x,y,z,x0,y0,z0,si,ct,st,cp,sp
c 378    format(1h ,'bound0:',11e12.3)
	return
	end

c-------------------------------------------------------
        subroutine boundf(x,y,z,xf,yf,zf,sf,ct,st,cp,sp,ierr)
c-------------------------------------------------------
c	It computes the upper limit of integration along a line of sight

	implicit real*8 (a-h, o-z)

	common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar

        ierr=0
        
c        write(6,200)
c 200    format(1h ,'entering boundf')
	ct1=ct
	if (ct .eq. 0.0) ct1=1.0d-6
	zf=d*ct1/dabs(ct1)
	yf=y+(zf-z)*st*sp/ct1
	xf=x+(zf-z)*st*cp/ct1
	rhof=dsqrt(xf*xf+yf*yf)
	sf=(zf-z)/ct1
c        write(6,299)xf,yf,zf,rhof,rhomax
c 299    format(1h ,'BOUNDF: xf,yf,zf,rhof,rhomax = ',5f8.1)
	if (rhof .le. rhomax) goto 100
	alpha=st*st
c
c next line rjt 11/01 to deal with case ray parallel to
c axis of cylinder but outside radius of cylinder
        if(alpha .eq. 0.0) then
c        write(6,599)x,y,z,ct,st,cp,sp,rhomax
c 599   format(1h ,'BOUNDF failed: x,y,z,ct,st,cp,sp,rhomax = ',8e12.4)
        ierr=1
        goto 100
        endif
	beta=2.0d0*st*(x*cp+y*sp)
	gamma=(x*x+y*y-rhomax*rhomax)
c
        test=beta*beta-4.0*alpha*gamma
c        test=beta*beta+4.0*alpha*gamma
        if(test .lt. 0.0)then
c       ray does not intersect cylinder because it started from outside
        write(6,699)x,y,z,ct,st,cp,sp,rhomax
 699   format(1h ,'BOUNDF aborting: x,y,z,ct,st,cp,sp,rhomax = ',8e12.4)
        ierr=1
        goto 100
        endif
c rjt	si=(-beta+dsqrt(beta*beta-4.0*alpha*gamma))/(2.0*alpha)
        dum1=-beta/(2.0d0*alpha)
        dum2=dsqrt(test)/(2.0d0*alpha)
c        write(6,269)dum1,dum2
c 269    format(1h ,'dum1,dum2 = ',2f8.1)
        sf=(-beta+dsqrt(test))/(2.0d0*alpha)
c        sf=(-beta-dsqrt(test))/(2.0d0*alpha)
	xf=x+sf*st*cp
	yf=y+sf*st*sp
	zf=z+sf*ct1
100	continue
c        write(6,201)
c 201    format(1h ,'exiting boundf')
	return
	end

c--------------------------------------------
        function br0(s0,x0,y0,z0,ct,st,cp,sp)
c--------------------------------------------

	implicit real*8(a-h, o-z)

	common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar
	common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
	common/f/fs0,fs1,fsdum,fsurad,dsb1,dsb2
	
	fract= 0.5
	esum= 0.0
	esumtot= 0.0
	sum= 0.0
	dstot= 0.0
	x1= x0+ s0* st* cp
	y1= y0+ s0* st* sp
	z1= z0+ s0* ct
        att=a(x1,y1,z1)
c
c	subroutine fstep gives the insegration step
c	depending on the position in the galaxy.
c        write(6,377)
c 377    format(1h ,'BR0: entering fstep')
600	call fstep(ds,x1,y1,z1,ct,s0,fs0,att)
c        write(6,378)
c 378    format(1h ,'BR0: leaving fstep')
	dstot= dstot+ ds
	if(dstot .gt. s0) goto900
	if(ds .eq. 0.) goto 1000
	
	delx= ds* st* cp
	dely= ds* st* sp
	delz= ds* ct
	x2= x1- delx* fract
	y2= y1- dely* fract
	z2= z1- delz* fract
c        if(x0 .eq. 0. .and. y0 .eq. 0.)write(6,28),x2,y2,z2 
c 28         format(1h ,'br0: x2,y2,z2 = ',3f10.3)
	att= a(x2,y2,z2)
	esumtot= esumtot+ att* ds
	esum= esumtot- att* ds* fract

c	The following trick to save time may give wrong results
c	for the problem that we are interested in.
c	if(esum .gt. 7.0) goto 880
	expesu= dexp(-esum)
c        rrr=sqrt(y2*y2+z2*z2)
c        dum=eta(x2,y2,z2)
c        if(rrr .lt. 3.0)write(6,6288)x2,y2,z2,dum
c 6288   format(1h ,4f6.2)
	sum= sum+ ds* eta(x2,y2,z2)* expesu
c       for theta=90
        
c        if(y0 .eq. 5 .and. abs(z0) .lt. 1.0)write(6,28),x2,y2,z2,
c     >eta(x2,y2,z2),expesu,sum 
c28         format(1h ,'br0 (theta=90): x2,y2,z2,sum = ',3f10.3,3e12.3)
c       For theta=0
c        if(y0 .eq. 5 .and. abs(x0) .lt. 1.0)write(6,29),x2,y2,z2,
c     >eta(x2,y2,z2),expesu,sum 
c29      format(1h ,'br0 (theta=0): x2,y2,z2,sum = ',3f10.3,3e12.3)
c        if(rrr .lt. 1.0)write(6,*)x2,ds
880	x1= x1- delx
	y1= y1- dely
	z1= z1- delz
	goto600
900	br0= sum
	goto 1100
1000	br0= 0.0
1100	return
	end

c--------------------------------------------
        function br1(s1,x0,y0,z0,ct,st,cp,sp)
c--------------------------------------------

	implicit real*8(a-h, o-z)

	common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar
	common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
	common/f/fs0,fs1,fsdum,fsurad,dsb1,dsb2
	common/atte/a0,eta0,rtruncate,sharp,athresh,rtrun,sha

c        write(6,5889)fs0,fs1
c 5889   format(1h ,'BR1: fs0,fs1 = ',2e12.4)
	fract= 0.5
	esum= 0.0
	esumtot= 0.0
	sum= 0.0
	dstot= 0.0
	x1= x0+ s1* st* cp
	y1= y0+ s1* st* sp
	z1= z0+ s1* ct
        att=a(x1,y1,z1)
        fs=fs1
600	call fstep(ds,x1,y1,z1,ct,s1,fs,att)
c        write(6,885)ds
c 885    format(1h ,'ds = ',e12.4)
	dstot= dstot+ ds
	if(dstot .gt. s1) goto 900
	if(ds .eq. 0.) goto 1000
	
	delx= ds* st* cp
	dely= ds* st* sp
	delz= ds* ct
	x2= x1- delx* fract
	y2= y1- dely* fract
	z2= z1- delz* fract

	att= a(x2,y2,z2)
	esumtot= esumtot+ att* ds
	esum= esumtot- att* ds* fract
	if(esum .gt. 7.0) goto 880
	expesu= dexp(-esum)
c        write(6,377)
c 377    format(1h ,'calling rdist0')
        if(att .ge. athresh)then
   	sum1= rdist0(x2,y2,z2,ct,st,cp,sp)
c        write(6,378)
c 378    format(1h ,'returning from rdist0')
  	sum= sum+ ds* att* sum1* expesu
        endif
c
c        if(att .lt. athresh)then
c        sum= sum+ ds* eta(x2,y2,z2)* expesu
c        endif
c
880	x1= x1- delx
	y1= y1- dely
	z1= z1- delz
	goto600
900	br1= sum
	goto 1100
1000	br1= 0.0
1100	return
	end

c--------------------------------------------
        function bst(s0,x0,y0,z0,ct,st,cp,sp)
c--------------------------------------------

	implicit real*8(a-h, o-z)

	common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar
	common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
	common/f/fs0,fs1,fsdum,fsurad,dsb1,dsb2

	fract= 0.5
	bstars= 0.0
	dstot= 0.0

c	The unattenuated starlight is calculated bellow.

	x1= x0
	y1= y0
	z1= z0
        att=a(x1,y1,z1)
        fs=fs0
600	call fstep(ds,x1,y1,z1,ct,s0,fs,att)
	dstot= dstot+ ds
	if (dstot .gt. s0) goto 900
	if (ds .eq. 0.) goto 1000
	
	delx= ds* st* cp
	dely= ds* st* sp
	delz= ds* ct
	x2= x1+ delx* fract
	y2= y1+ dely* fract
	z2= z1+ delz* fract

	bstars= bstars+ ds* eta(x2,y2,z2)

	x1= x1+ delx
	y1= y1+ dely
	z1= z1+ delz
	goto 600
900	bst= bstars
	goto 1100
1000	bst= 0.0
1100	return
	end

c-------------------------------------
        subroutine direct(ct,st,cp,sp)
c-------------------------------------
c
	implicit real*8 (a-h, o-z)

c	It calculates the directions that will be used to approximate
c	the integrals over all angles.

c	real*8 cthet1, cthet2, cthet, wt1, wt2, wt, phi, wp

	common/dir/cthet(50),sthet(50),cphi(50),sphi(50),phi(50)
	common/dir1/cosg(50,50),upr(50,50),vpr(50,50),wpr(50,50)
	common/wei/wt(50),wp(50)
	common/g/g,g2,cosalp
	common/integ/k,m,n,km

	dimension cthet1(50), cthet2(50), wt1(50), wt2(50)
	dimension phi1(50), phi2(50), wp1(50), wp2(50)
c
c	direction cosines of the direction of interest in lab frame
c
	u=st*cp
	v=st*sp
	w=ct
c next two line commented rjt 11/01
c since for perfect face on case the 10^-6 trap is avoided if w .ne. 1)
c	if (w .eq. 1.0d0) w=9.9999d-1
c	if (w .eq. -1.0d0) w=-9.9999d-1
	kk=k
	mm=m
	nn=n/2

c	we break the integral over mu into two parts
c        write(6,288)cosalp,mm
c 288    format('cosalp,mm = ',f8.3,i6)
	call gauss(-1.0d0,cosalp,cthet1,wt1,mm,0)
c        write(6,288)cosalp,kk
	call gauss(cosalp,1.0d0,cthet2,wt2,kk,0)
c	we also break the integral over phi into two parts

c        write(6,289)nn
c 289    format('nn = ',i6)
	call gauss(0.d0, 3.1416d0, phi1, wp1, nn, 0)
c        write(6,289)nn
	call gauss(3.1416d0, 6.2832d0, phi2, wp2, nn, 0)

c	combine the two parts into one

	do 100 i=1, nn
	cphi(i)=dcos(phi1(i))
	sphi(i)=dsin(phi1(i))
	wp(i)= wp1(i)
100	continue
	do 110 i=1, nn
	cphi(i+nn)= dcos(phi2(i))
	sphi(i+nn)= dsin(phi2(i))
	wp(i+nn)= wp2(i)
110	continue

	do 150 j=1, m
	cthet(j)=cthet1(j)
	sthet(j)=dsqrt(1.0-cthet(j)*cthet(j))
	wt(j)=wt1(j)
150	continue
	do 170 j=1, k
	cthet(j+m)=cthet2(j)
	sthet(j+m)=dsqrt(1.0-cthet(j+m)*cthet(j+m))
	wt(j+m)=wt2(j)
170	continue

c	compute cosg and convert directions into lab frame

	do 300 j=1, km
	a=cthet(j)
	b=sthet(j)

	do 200 i=1, n
	c=cphi(i)
	d=sphi(i)

	cosg(i,j)=cthet(j)
	if (1.0- dabs(w) .lt. 1.0d-6) goto 180
	upr(i,j)=(b*c*w*u-b*d*v)/dsqrt(1.0-w*w)+a*u
	vpr(i,j)=(b*c*w*v+b*d*u)/dsqrt(1.0-w*w)+a*v
	wpr(i,j)=-b*c*dsqrt(1.0-w*w)+a*w
	goto 190
180	continue
	upr(i,j)=b*c
	vpr(i,j)=b*d
	wpr(i,j)=a*w
190	continue
200	continue
300	continue
	return
	end

c---------------------------
        function eta(x,y,z)
c---------------------------
c	Space distribution of stars (see eq. 12)

	implicit real*8 (a-h, o-z)
	common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  eta1arm,a1arm,narm,rup,fillarm
	common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar
	common/atte/a0,eta0,rtruncate,sharp,athresh,rtrun,sha
	common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
        common/ferrers/etaf,af,bf,cf,thetaf,ctf,stf,iferrers
        common/ellexp/eta2,zs2,hsa2,hsb2,thetae,cte,ste,
     >  hda2,hdb2,zd2,tau2,a2,iellexp
	common/bcons/etabmax,rin,ibulge
	common/translate/pix2kpc
        common/config/izero,idust,itau,istars,imask,iall,
     >  shiftx,shifty,ifswitch,ifinns1,ifinns2,ifinnd1,ifinnd2


c	The DISK(S)
	rho=dsqrt(x*x+y*y)


c
c       overall truncation 
c
c       first disk and all structures except second star and dust disk
        sigtruncate=sha*rtrun
        rplane=sqrt(x*x+y*y+z*z)
c        truncate=0.5*(1.0-derf((x-rtrun)/sigtruncate))
        truncate=0.5*(1.0-derf((rplane-rtrun)/sigtruncate))
c        write(6,*)x,y,z,rplane,sigtruncate,rtrun,truncate
c
c       second disk
        sigtruncate1=sha1*rtrun1
        truncate1=0.5*(1.0-derf((rplane-rtrun1)/sigtruncate1))


c       flared disk for the first (old) stellar disk
         zs01 = zs
         if (zsin .gt. zs) then
          flares01 = dlog10((zssolar - zs)/(zsin - zs))
         else 
          flares01 = 0.
         end if
         flares02 = dlog10(hssolar/hsin)
         flares00 = flares01/flares02
         zs01 = zs + (zsin - zs) * (rho/hsin)**flares00

c       flared disk for the second (young) stellar disk
         zs11 = zs1
         if (zs1in .gt. zs1) then
          flares11 = dlog10((zs1solar - zs1)/(zs1in - zs1))
         else 
          flares11 = 0.
         end if
         flares12 = dlog10(hs1solar/hs1in)
         flares11 = flares11/flares12
         zs11 = zs1 + (zs1in - zs1) * (rho/hs1in)**flares11
c
        if(idisk1 .eq. 0) then etadisk0=0.

c       double exponential disks
        if(idisk1 .eq. 1) then
         if (Rho .le. hsin) then
          if (ifinns1 .eq. 1) then 
           etadisk0 =  eta0 * (zs/zs01) * 
     >      dexp(-hsin/hs-dabs(z)/zs01)
          endif
          if (ifinns1 .eq. 2) then 
           etadisk0 =  eta0 * (zs/zs01) * (1./hsin) * 
     >      rho * dexp(-hsin/hs-dabs(z)/zs01)
          endif
          if (ifinns1 .eq. 3) then 
           etadisk0 =  (eta0/2.) * (zs/zs01) * 
     >      (rho/hsin+1.) *  dexp(-hsin/hs-dabs(z)/zs01)
          endif
         endif
         if (Rho .gt. hsin) then         
          etadisk0= eta0* (zs/zs01) * 
     >     dexp(-rho/hs- dabs(z)/zs01) 
         endif
        endif       
c       constant emissivity to test energy densities
c        etadisk0= eta0*dfloat(idisk1)

c       sech2 vertical distribution
        if(idisk1 .eq. 3) then
         sech01=1.d0/cosh(dabs(z)/zs01)
   	 sech02=sech01*sech01
         if (Rho .le. hsin) then
          if (ifinns1 .eq. 1) then
           etadisk0 = eta0 * (zs/zs01) * 
     >      dexp(-hsin/hs)*sech02
          endif
          if (ifinns1 .eq. 2) then
           etadisk0 = eta0 * (zs/zs01) * (1./hsin) * 
     >      rho * dexp(-hsin/hs)*sech02
          endif          
          if (ifinns1 .eq. 3) then
           etadisk0 = (eta0/2.) * (zs/zs01) * 
     >      (rho/hsin+1.) * dexp(-hsin/hs)*sech02
          endif          
         endif
         if (Rho .gt. hsin) then
           etadisk0= eta0 * (zs/zs01) * 
     >      dexp(-rho/hs)*sech02
         endif
        endif
c next line is constant emissivity within a sphere
c        etadisk0= dfloat(idisk1)*eta0
c next 2 lines is plane disk of uniform emissivity
c        etadisk0=0.0
c        if(abs(z) .lt. zs)etadisk0=dfloat(idisk1)*eta0

        if(idisk2 .eq. 0) then etadisk1=0.

c       double exponential disk
        if(idisk2 .eq. 1) then
         if (Rho .le. hs1in) then 
          if (ifinns2 .eq. 1) then
           etadisk1 = eta1 *
     >      (zs1/zs11)*dexp(-hs1in/hs1- dabs(z)/zs11)
          endif
          if (ifinns2 .eq. 2) then
           etadisk1 =  eta1 *
     >      (zs1/zs11)*(1./hs1in)*rho*dexp(-hs1in/hs1- dabs(z)/zs11)
          endif
          if (ifinns2 .eq. 3) then
           etadisk1 = (eta1/2.) *
     >      (zs1/zs11)*(rho/hs1in+1.)*dexp(-hs1in/hs1- dabs(z)/zs11)
          endif
         endif
         if (Rho .gt. hs1in) then
          etadisk1 = eta1 *
     >     (zs1/zs11)*dexp(-rho/hs1- dabs(z)/zs11)
         endif
        endif

c spiral arm radial distribution
        if(idisk2 .eq. 2) then      
          etadisk1= etaarm(rho) * 
     >     (zs1/zs11) * dexp(-dabs(z)/zs11)
        endif

c   sech2 vertical distribution
        if(idisk2 .eq. 3) then
	 sech11=1.d0/cosh(dabs(z)/zs11)
	 sech12=sech11*sech11
         if (Rho .le. hs1in) then         
          if (ifinns2 .eq. 1) then
           etadisk1 = eta1 *
     >       (zs1/zs11)*dexp(-hs1in/hs1)*sech12
          endif
          if (ifinns2 .eq. 2) then
            etadisk1 = eta1 *
     >       (zs1/zs11)*(1./hs1in)*rho*dexp(-hs1in/hs1)*sech12
          endif
          if (ifinns2 .eq. 3) then
            etadisk1 = (eta1/2.) *
     >       (zs1/zs11)*(rho/hs1in+1.)*dexp(-hs1in/hs1)*sech12
          endif
         endif
         if (Rho .gt. hs1in) then
          etadisk1=eta1*
     >    (zs1/zs11)*dexp(-rho/hs1)*sech12
         endif

        endif             

       
        etadisk0=etadisk0*truncate
        etadisk1=etadisk1*truncate1
c
c	The BULGE 
	acap= ((x*x+ y*y+ (z*z/ellipt2))**0.5)/ reff
c        write(6,6289)x,y,z,reff,ellipt2,acap
c 6289   format(1h ,'x,y,z,reff,ellipt2,acap: ',3f10.4,3f6.2)
c       note - this differs from the definition given in Popescu et al.
c       of acap= (x*x+ y*y+ z*z*a*a/(b*b))**0.5/ reff
c	Modified Hubble profile:
c        acaphubblemax=10000000.0d0
        acaphubblemax=3.0d0
	if(ibulge.eq.1 .and. acap .le. acaphubblemax)etabulge=etab* 
     #               (1.0d0+acap*acap)**(-1.5d0)
c        write(6,6789)etabulge,etab,acap
c 6789   format(1h ,'etabulge,etab,acap: ',3f10.6)
c	R^1/4 law:
	if(ibulge.eq.2) then
          acapdum=max(acap,rin)
c          acapdum=acap
c       rrr=sqrt(x*x+ y*y+ z*z)
c       if(rrr .lt. 3.)write(6,2688)x,y,z,reff,acapdum,etab,acap,rin
c 2688  format(1h ,'ETA: ',8f6.2)
       dum=7.67d0*(acapdum**0.25d0)
       etabulge=dexp(-dum)
c       if(rrr .lt. 3.)write(6,2668)dum,etabulge
c 2678  format(1h ,'ETA: ',f6.2)
c 2668  format(1h ,'ETA: ',2f6.2)
       etabulge=etabulge/(acapdum**0.875d0)
c       if(rrr .lt. 3.)write(6,2678)etabulge
       etabulge=etabulge*etab
c RT 7/8/07       if(acapdum .gt. 3.0d0)etabulge=0.0d0
        if(acapdum .gt. 10.0d0)etabulge=0.0d0
c       if(rrr .lt. 3.)write(6,2678)etabulge
         endif

	if(ibulge.eq.3) then
          acapdum=max(acap,rin)
c          acapdum=acap
c       rrr=sqrt(x*x+ y*y+ z*z)
c       if(rrr .lt. 3.)write(6,2688)x,y,z,reff,acapdum,etab,acap,rin
c 2688  format(1h ,'ETA: ',8f6.2)
       esersic=1.0d0/float(nsersic)
       dum=bsersic*(acapdum**esersic)
       etabulge=dexp(-dum)
c       if(rrr .lt. 3.)write(6,2668)dum,etabulge
c 2678  format(1h ,'ETA: ',f6.2)

c 2668  format(1h ,'ETA: ',2f6.2)
       dsersic=2.0d0*float(nsersic)
       dsersic=(dsersic-1.0d0)/dsersic
       etabulge=etabulge/(acapdum**dsersic)
c       if(rrr .lt. 3.)write(6,2678)etabulge
       etabulge=etabulge*etab
c       if(acapdum .gt. 3.0d0)etabulge=0.0d0       if(acapdum .gt. 10.0d0)etabulge=0.0d0
c       if(acapdum .gt. 10.0d0)etabulge=0.0d0
         if(acapdum .gt. 3.0d0)etabulge=0.0d0 
c        if(acapdum .gt. 10.0d0)etabulge=0.0d0
c       if(rrr .lt. 3.)write(6,2678)etabulge
         endif
       if(ibulge .le. 0)etabulge=0.0d0

        etabulge=etabulge*truncate

c
c	The HALO
c	halocap= (x*x+ y*y+ z*z)**2
c	halo0= eta0*200.
c	rhalo=2.83/pix2kpc
c	halo=halo0*(rhalo*rhalo+ halocap)**(-1.75)
c
c       ferrers bar
c
        etaferrers=0.0d0
        if(iferrers .eq. 1)then
        yf=y*ctf - x*stf
        xf=x*ctf + y*stf
c                write(6,*)x,y,ctf,stf,xf,yf,z
        rhof2=(yf*yf)/(af*af)+(xf*xf)/(bf*bf)+(z*z)/(cf*cf)
        rhof2=min(1.0,rhof2)
c        if(rhof2 .lt. 1.0)write(6,*)rhof2
        etaferrers=(1.0-rhof2)        
c        write(6,*)etaferrers
        etaferrers=etaf*etaferrers*etaferrers
c        if(rhof2 .lt. 1.0)write(6,*)etaferrers
       etaferrers=etaferrers*truncate
        endif      

c       elliptical exponential disk

        etaellexp=0.0d0
        if(iellexp .eq. 1)then
        ye=y*cte - x*ste
        xe=x*cte + y*ste
        ye=dabs(ye)
        xe=dabs(xe)
c                write(6,*)x,y,cte,ste,xe,ye,z
        rhoe2=(ye*ye)/(hsa2*hsa2)+(xe*xe)/(hsb2*hsb2)
        dummy=0.0d0 -dsqrt(rhoe2) -(dabs(z)/zs2)
        dum=max(-30.0d0,dummy)
        etaellexp=eta2 * dexp(dum)                
        etaellexp=etaellexp*truncate
        endif     

c        if(rhof2 .lt. 1.)write(6,*)etadisk0,etadisk1,etabulge,etaferrers
c        write(6,*)etadisk0,etadisk1,etabulge,etaferrers
	eta= etadisk0+etabulge+etadisk1+etaferrers+etaellexp
c	if(ibulge.eq.0) eta= etadisk0+etadisk1
c        write(6,2789)eta,etadisk0,etabulge,truncate
c 2789   format(1h ,'eta,etadisk0,etabulge,truncate: ',4f10.6)
	return
	end

c---------------------------
        function etaarm(rho)
c---------------------------
c	finds mid-plane value for emissivity in second disk
c       for the case that the second disk is prescribed as a
c       set of concentric emitting zones, with emissivity contrast
c       10^6, to mimic spiral arms with no interarm emission
c       The distribution of emission from the arms is defined by the
c       parameters (already set in ARMSET in the common block NEWDISK):
c       NARM: number of arms 
c       FILLARM: linear filling factor of arms along a radius
c       RUP: determines rough position of outer edge of
c            last spiral arm in units of scale length
c       The other input parameters are:
c       RHO: the radius at which the emissivity is to be calculated
c       ETA1: the central face-on brightness of a filled
c       exponential disk with the same volume-integrated emissivity
c       as the set of concentric emission zones, within the volume
c       bounded in radius by the truncation radius RTRUN1
c       RTRUN1: truncation radius of 2nd disk 
c               (used in definition of ETA1, above)
c       ETA1ARM: constant which determines the normalisation of the
c                brightness scaling to the pure exponential disk
c                of central brightness ETA1. ETA1ARM is assigned
c                in a call to the subroutine ARMSET which must be called
c                before ETAARM is called for the first time. The value
c                of ETA1ARM is such that the volume-integrated
c                emission of the set of concentric emission zones
c                out to the truncation radius RTRUN1 of the second disk
c                is the same as the volume-integrated emission
c                out to RTRUN1 of a filled exponential disk with
c                central emissivity ETA1
c
c       note that the routine ARMSET should be run before this
c       is run for the first time 
c
	implicit real*8 (a-h, o-z)
	common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  eta1arm,a1arm,narm,rup,fillarm
c
c       1.0/etainter is the contrast interarm/arm
        etainter=1.0d-6
c
        etaarm0=eta1arm*etainter
        etaarm=etaarm0
c
        spacearm=rup/dfloat(narm)
        width2=spacearm*fillarm/2.0d0
        iarm=0
100     continue
        iarm=iarm+1
        rarm=spacearm*(-0.5d0+dfloat(iarm))
        rarmin=rarm-width2
c        if(iarm .eq. 1)rarmin=0.0d0
        rarmout=rarm+width2
        if(rho .ge. rarmin .and. rho .lt. rarmout)then
        etaarm=eta1arm*dexp(-rarm/hs1)
        goto 200
        endif
        if(iarm .lt. narm)goto 100
200     continue
c
	return
	end
c---------------------------
        function aarm(rho)
c---------------------------
c	finds mid-plane value for dust density in second disk
c       for the case that the second disk is prescribed as a
c       set of concentric emitting zones, with dust column density contrast
c       10^6, to mimic spiral arms with no interarm dust.
c       The distribution of dust in the arms is defined by the
c       parameters (already set in ARMSET in the common block NEWDISK):
c       NARM: number of arms 
c       FILLARM: linear filling factor of arms along a radius
c       RUP: determines rough position of outer edge of
c            last spiral arm in units of scale length
c       The other input parameters are:
c       RHO: the radius at which the emissivity is to be calculated
c       A1: the central dust opacity density of a filled
c       exponential disk with the same total dust mass as the dust
c       as the set of concentric emission zones, within the volume
c       bounded in radius by the truncation radius RTRUN1
c       RTRUN1: truncation radius of 2nd disk 
c               (used in definition of A1, above)
c       A1ARM: constant which determines the normalisation of the
c                brightness scaling to the pure exponential disk
c                of central faceon column density TAU1. A1ARM is assigned
c                in a call to the subroutine ARMSET which must be called
c                before AARM is called for the first time. The value
c                of A1ARM is such that the volume-integrated
c                mass of dust in the set of concentric emission zones
c                out to the truncation radius RTRUN1 of the second disk
c                is the same as the volume-integrated dust mass
c                out to RTRUN1 of a filled exponential disk with
c                central dust opacity density A1
c
c       note that the routine ARMSET should be run before this
c       is run for the first time 
c
	implicit real*8 (a-h, o-z)
	common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  eta1arm,a1arm,narm,rup,fillarm
c
        ainter=1.0e-6
c       
        aarm0=etadisk1*etainter
        aarm=aarm0
c
        spacearm=rup/float(narm)
        width2=spacearm*fillarm/2.0
        iarm=0
100     continue
        iarm=iarm+1
        rarm=spacearm*(-0.5+float(iarm))
        rarmin=rarm-width2
c        if(iarm .eq. 1)rarmin=0.0d0
        rarmout=rarm+width2
        if(rho .ge. rarmin .and. rho .lt. rarmout)then
        aarm=a1arm*dexp(-rarm/hs1)
        goto 200
        endif
        if(iarm .lt. narm)goto 100
200     continue
c
	return
	end
c---------------------------
        subroutine armset
c---------------------------
c       Initializes parameters needed by functions ETAARM and A1ARM
c
c       Input parameters:
c       A1: the central face-on opacity density of a filled
c       exponential disk with the same total dust mass
c       as the set of concentric emission zones, within the volume
c       bounded in radius by the truncation radius RTRUN1
c       RTRUN1: truncation radius of 2nd disk 
c               (used in definition of A1ARM, above)
c
c       Output parameter:
c       A1ARM:   constant which determines the normalisation of the
c                dust mass  scaling to the pure exponential disk
c                of central brightness opacity density A1. A1ARM is assigned
c                in a call to the subroutine ARMSET which must be called
c                before AARM is called for the first time. The value
c                of A1ARM is such that the volume-integrated
c                dust mass of the set of concentric emission zones
c                out to the truncation radius RTRUN1 of the second disk
c                is the same as the volume-integrated dust mass
c                out to RTRUN1 of a filled exponential disk with
c                central dust opacity density A1
c       ETA1ARM: same for emissivity of second dust disk
c       NARM: number of arms 
c       FILLARM: linear filling factor of arms along a radius
c       RUP: determines rough position of outer edge of
c            last spiral arm in units of scale length
c
c       note that the routine ARMSET should be run before this
c       is run for the first time 
c
	implicit real*8 (a-h, o-z)
	common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  eta1arm,a1arm,narm,rup,fillarm
c
c       in the calls to aarm and etaarm used here we use a1 and eta1
c       a1arm and eta1arm, before setting these values at the end of
c       this routine
c
        a1arm=a1
        eta1arm=eta1
c
        pie=3.1415927d0
        eta1in=eta1
        a1in=a1
        if(eta1 .le. 0.0d0)then
        write(6,201)
201     format(1h ,'ARMSET: eta1 not initialised')
        write(6,202)eta1
202     format(1h ,'... aborting ... eta1 = ',15e12.6)
        goto 555
        endif
        if(a1 .le. 0.0 .and. tau1 .gt. 0.0)then
        write(6,301)
301     format(1h ,'ARMSET: eta1 not initialised')
        write(6,302)a1
302     format(1h ,'... aborting ... a1 = ',15e12.6)
        goto 555
        endif
c
        narm=3
        rup=hs1*3.0d0
c        rup=hs1*2.0d0
        fillarm=0.5d0
c        fillarm=0.25d0
c
        ncalc=100001
        armflux=0.0d0
        diskflux=0.0d0
        dr=rtrun1/dfloat(ncalc)
        do 100 i=1,ncalc
        ii=i-1
        rprev=dfloat(ii)*dr
        r=dfloat(i)*dr
        b=etaarm(r)
        bdisk=eta1*dexp(-r/hs1)
        armflux=armflux+pie*b*(r*r-rprev*rprev)
        diskflux=diskflux+pie*bdisk*(r*r-rprev*rprev)
100     continue
c  analytical disk total flux
c        xx=rtrun1/hs1
c        dum1=1.0d0-dexp(-xx)
c        dum1=dum1-xx*dexp(-xx)
c        diskflux=2.0d0*pie*eta1in*dum1
c        diskflux=diskflux*hs1*hs1
        eta1arm=eta1arm*(diskflux/armflux)
c       since dust in second disk has same radial geometry as
c       the emissivity the factor is the same
        a1arm=a1*(diskflux/armflux)
        write(6,790)
790     format('ARMSET:')
        write(6,890)armflux
890     format('armflux  =  ', 15e12.6)
        write(6,891)eta1
891     format('eta1     =  ', 15e12.6)
        write(6,892)diskflux
892     format('diskflux =  ', 15e12.6)
        write(6,893)eta1arm
893     format('eta1arm  =  ', 15e12.6)
        write(6,894)a1
894     format('a1       =  ', 15e12.6)
        write(6,895)a1arm
895     format('a1arm    =  ', 15e12.6)
c
555     continue
c
	return
	end

c---------------------------------------------
        subroutine fstep(ds,x1,y1,z1,ct,s0,fs,att)
c---------------------------------------------
c
c original with ellipt added to common BULGE
c
	implicit real*8(a-h, o-z)

	common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar
	common/bcons/etabmax,rin,ibulge
	common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
	common/f/fs0,fs1,fsdum,fsurad,dsb1,dsb2
        common/config/izero,idust,itau,istars,imask,iall,
     >  shiftx,shifty,ifswitch,ifinns1,ifinns2,ifinnd1,ifinnd2
	common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  eta1arm,a1arm,narm,rup,fillarm
c
c       ifswitch is 0 then we use Nick's integration
c       scheme for inclinations of 80 degrees or more
c       and rt's scheme for inclinations less tha 80 degrees
c
c       if ifswitch is 1 or 2 or 3 or 4 then we use rt's scheme for all
c       inclinations. For ifswitch=1 the step length is
c       determined purely from the geometry of the old disks.
c       For ifswitch=2 the step length is also cconstrained
c       by the need to make it smaller than a certain fraction
c       of a tau. ifswitch = 3 is the same as ifswitch = 1
c       except that the disk scale height is used to determine
c       the step length for all inclinations, whereas for ifswitch = 1
c       the disk scale length is used for inclinations approaching 90. 
c
c
c       maximum fraction of a tau permitted for each step
c       in RT's scheme with ifswitch=2
        taufractmax=0.2
c
	ct1= ct
	if(dabs(ct1) .eq. 0.) ct1= 1.d-6

c       for inclinations greater than 80 degrees use the 
c       original fstep routine from Nick which is quick and 
c       sufficiently accurate for edge on galaxies. 
c       For inclinations less than 80 degrees use the more 
c       accurate fstep needed for the face on galaxies

c       80 degrees
        if(dabs(ct1) .gt. 0.173648)goto 299 

        if(ifswitch .ge. 1)goto 299

c       kylafis scheme for ifswitch=0 and near edge on:

	ds0= fs* hs
	if(dabs(ct1) .gt. zd/hs) ds0= fs* zd/ dabs(ct1)

c        dum=1.0-ct1*ct1
c        st=sqrt(dum)
c        scalestars=hs/st
c        ss=dabs(zs/ct1)
c        if(ss .lt. scalestars)scalestars=ss
c        ds0crit=fs*scalestars
c        if(ds0 .gt. ds0crit)ds0=ds0crit
c        dtau=a(x1,y1,z1)*ds0
c        dtaulim=0.1
c        if(dtau .gt. dtaulim)ds0=ds0*(dtaulim/dtau)      
        
c        write(6,212)s0,ds0
c 212    format(1h ,'s0,ds0 = ',2f12.6) 

	isbins= s0/ ds0
        if(isbins .le. 0)isbins=1
	dsd= s0/ dfloat(isbins)

	dsb= dsd
	zone1= 0.5
	zone2= 3.0+ ds0/reff

        zone1=1.0
        zone2=3.0*reff         

c	zone1= 1.0
c	zone2= 3.0+ ds0/reff

c       	acap= (x1*x1+ y1*y1+ z1*z1/ ellipt2)**0.5/ reff

        acap= (x1*x1+ y1*y1+ z1*z1/ ellipt2)**0.5

c	if(acap .le. zone1) dsb= dsb1
c	if(acap .gt. zone1 .and. acap .le. zone2) dsb= dsb2

        if(acap .lt. zone2)then dsb=0.1*reff


	ds= min(dsd, dsb)


        
c        if(acap .le. zone1)then
c        write(6,*)dsd,dsb,ds
c        endif

c        if(acap .gt. zone2)write(6,*)fs,dsd,dsb,ds

        goto 399
 299    continue

c        
c        write(6,377)x1,y1,z1,ct,s0,fs
c 377    format(1h ,'FSTEP: x1,x2,x3,ct,s0,fs = ',6e12.4)
c        write(6,299)hs,zd,reff,ellipt2,dsb1,dsb2
c 299    format(1h ,'hs,zd,reff,ellipt2,dsb1,dsb2 = ',6e12.4)
c
	ct1= ct
	if(dabs(ct1) .eq. 0.) ct1= 1.d-6
        hscalemin=hs
        vscalemin=zd
        vscalemax=zs
        if(ifswitch .eq. 4)then
          hscalemin=hs1
          vscalemin=zd1
          vscalemax=zs1
        endif
        if(ibulge .le. 0)vscalemax=vscalemax*dabs(ct1)
	ds0= fs*hscalemin
	if(ifswitch .eq. 3)ds0= fs*vscalemin
        dsbulge0=0.1*reff
c	if(dabs(ct1) .gt. vscalemin/hscalemin) ds0= fs*zd/dabs(ct1)
        if(dabs(ct1) .gt. vscalemin/hscalemin) then
          ds0= fs*vscalemin
          if(ibulge .le. 0 .and. dabs(ct1) .gt. 0.0)ds0=ds0/dabs(ct1)
        endif
c        write(6,378)ct1,hscalemin,vscalemin,ds0
c 378    format(1h ,'FSTEP: ct1,hscalemin,vscalemin,ds0 = ',4e12.4)
	isbins= s0/ ds0
        isbinsbulge=s0/dsbulge0
c added rjt 11/01
        if(isbins .le. 1)then
c        write(6,2878)s0,ds0,fs
c 2878   format(1h ,'isbinstrap',3e12.6)
        isbins=1
        ds0=s0
        endif
        if(isbinsbulge .le. 1)then
c        write(6,2878)
c 2878   format(1h ,'isbinstrap')
        isbinsbulge=1
        dsbulge0=s0
        endif
c        write(6,379)isbins
c 379    format(1h ,'FSTEP: isbins = ',i8)
	dsd= s0/ dfloat(isbins)
c added rjt 11/01
c        if(dsd .le. 0.5)dsd=1.0
c    lower limit in pixels everywhere
c        if(dsd .le. 0.01)dsd=0.01
c    upper limit in pixels away from plane
c        zlim=4.0*vscalemax
c        zlimbulge=4.0*reff
        zlim=4.5*vscalemax
        zlimbulge=4.5*reff
        rlim=(x1*x1+ y1*y1)**0.5/ reff
c        rlimit=4.0
        rlimit=3.0
        if(rlim .lt. rlimit)then
           if(zlimbulge .gt. zlim)then
              zdum=(zlimbulge-zlim)*(rlimit-rlim)/rlimit
              zlim=zlim+zdum
           endif
        endif   

c        if(rlim .lt. rlimit)then
c           zlim=max(zlim,zlimbulge)
c           if(dsd .gt. dsbulge0)dsd=dsbulge0
c        endif




        zdiff=dabs(z1)-zlim
        if(dabs(z1).gt.zlim)dsd=1.0+zdiff
	dsb= dsd
	zone1= 0.5
	zone2= 3.0+ ds0/reff
	acap= (x1*x1+ y1*y1+ z1*z1/ ellipt2)**0.5/ reff
	if(acap .le. zone1) dsb= dsb1
	if(acap .gt. zone1 .and. acap .le. zone2) dsb= dsb2
	ds= min(dsd, dsb)
c       turrn off use of bulge to determine integration step
        ds=dsd
c
c       for ifswitch=2 check step not too large a fraction of a tau
c        
        if(ifswitch .eq. 2)then
         test = ds * att
         if(test .gt. taufractmax)then
c          write(6,9568)ds,att,test
c 9568     format(1h ,'FSTEP: ds,att,test = ',3e12.4)
          ds = taufractmax/att
         endif
        endif 
c
c        write(6,199)ds
c 199	format(1h ,'FSTEP: returning with ds = ',e12.4)

 399    continue

        return
	end

c------------------------------------
      SUBROUTINE GAUSS(A,B,XX,W,N,IP)
c------------------------------------

c	It computes the points XX in the interval (A,B)
c	and the corresponding weights W. The integral is then
c	sum over i of f[XX(i)] * W(i), where f is the function
c	to be integrated.

      implicit real*8 (a-h,o-z)

      DIMENSION XX(8),W(8)
      DATA CONV/1.D-14/
      M=(N+1)/2
      K=N-M
      DO 40 I=1,M
      Z=DCOS(3.14159d0*(M-I+.75d0)/(N+.5d0))
10    P1=1.d0
      P2=0.d0
      DO 20 J=1,N
      P3=P2
      P2=P1
20    P1=((2.d0*J-1.d0)*Z*P2-(J-1)*P3)/J
      PP=N*(Z*P1-P2)/(Z*Z-1.d0)
      Z1=Z
      Z=Z1-P1/PP
      IF(DABS(Z-Z1).GT.CONV  )GO TO 10
      K=K+1

C      EVALUATE GAUSSIAN QUADRATURE POINTS AND WEIGHTS ON INTERVAL A-B

      XX(K)=Z
40    W(K)=2.d0/((1.d0-Z*Z)*PP*PP)
      P1=(B-A)*.5d0
      P2=(B+A)*.5d0
C      EVEN
      J=N-M+1
      K=J
      IF(2*M-N.EQ.0) GO TO 31
C      ODD
      XX(K)=XX(K)*P1+P2
      W(K)=P1*W(K)
      IF(N.EQ.1) GO TO 32
      J=J+1
31    DO 30 I=J,N
      P3=XX(I)*P1
      XX(I)=P3+P2
      W(I)=P1*W(I)
      K=K-1
      XX(K)=-P3+P2
      W(K)=W(I)
30    CONTINUE
32    IF(IP.EQ.0) RETURN
      WRITE(6,100) N,A,B
      WRITE(6,101) (I,XX(I),W(I),I=1,N)
      RETURN
100   FORMAT(//,38H GAUSSIAN QUADRATURE POINTS,WEIGHTS,N=,I5,13H  ON INT
     *ERVAL,2D12.4,/,4X,1HI,6X,1HX,19X,1HW)
101   FORMAT(I5,2F20.13)
      END

c--------------------------------------------------
        function phase(i,j)
c--------------------------------------------------
c	Expression for the phase function (see eq. 4)

	implicit real*8 (a-h, o-z)

c	real*8 cthet, phi

	common/dir/cthet(50),sthet(50),cphi(50),sphi(50),phi(50)
	common/dir1/cosg(50,50),upr(50,50),vpr(50,50),wpr(50,50)
	common/g/g,g2,cosalp
c
	phase=(1.0-g2)/(1.0d0+g2-2.0d0*g*cosg(i,j))**1.5d0
	return
	end

c-------------------------------------------------
        function rdist0(xr,yr,zr,ct,st,cp,sp)
c-------------------------------------------------
c
	implicit real*8 (a-h, o-z)

c	It does the integral over all angles of (I0 * phase function)

c	real*8 cthet, phi, wt, wp

	common/f/fs0,fs1,fsdum,fsurad,dsb1,dsb2
	common/dir1/cosg(50,50),upr(50,50),vpr(50,50),wpr(50,50)
	common/wei/wt(50),wp(50)
	common/integ/k,m,n,km
c
	sum1=0.0d0
c        write(6,375)k,m,n,km
c 375    format(1h ,'RDIST0: k,m,n,km = ',4i6)
	do 300 j=1, km
	do 200 i=1, n
	ctp=wpr(i,j)
	stp=dsqrt(1.0d0-ctp*ctp)
c        write(6,379)i,j,ctp,stp
c 379    format(1h ,'RDIST0: i,j,ctp,sto = ',2i6,2e12.4)
	cpp=upr(i,j)/stp
	spp=vpr(i,j)/stp
c        write(6,380)cpp,spp
c 380    format(1h ,'RDIST0: cpp,spp = ',2e12.4)
c        write(6,377)
c 377    format(1h ,'RDIST0: calling bound0')
	call bound0(xr,yr,zr,x0r,y0r,z0r,sir,ctp,stp,cpp,spp,ierr)
c        write(6,378)
c 378    format(1h ,'RDIST0: returning from bound0:')
c        write(6,383)
c 383    format(1h ,'xr,yr,zr,x0r,y0r,z0r,sir,ctp,stp,cpp,spp = ')
c        write(6,384)xr,yr,zr,x0r,y0r,z0r,sir,ctp,stp,cpp,spp
c 384    format(1h ,11e12.4)
c        write(6,381)
c 381    format(1h ,'RDIST0: calling br0')
	sum0=br0(sir,x0r,y0r,z0r,ctp,stp,cpp,spp)
c        write(6,382)
c 382    format(1h ,'RDIST0: leaving br0')
	sum1=sum1+sum0*phase(i,j)*wt(j)*wp(i)
200	continue
300	continue
	rdist0=sum1/4.0d0/3.141592654d0
	return
	end

c---------------------------------------------------
        function tau(s0,x0,y0,z0,ct,st,cp,sp)
c---------------------------------------------------
c
	implicit real*8 (a-h, o-z)

crjt	dimension att0(3000)
crjt        dimension att0(50000)
        dimension att0(400000)
c
	common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar
	common/f/fs0,fs1,fsdum,fsurad,dsb1,dsb2
c
	ct1=ct
	if (dabs(ct) .eq. 0.0d0) ct1=1.0d-6
	ds0=fs0*hd
	if (dabs(ct1) .gt. zd/hd) ds0=fs0*zd/dabs(ct1)
c	this makes the assumptions that zd < zs and hd > hs
	isbins=s0/ds0
	if (isbins .eq. 0) goto 2000
	ds=s0/dfloat(isbins)
c in following three lines rjt replaced 2999 with 99999
	if (isbins .gt. 399999) write(1,100) isbins, s0, ds0, hd, zd
	if (isbins .gt. 399999) write(2,100) isbins, s0, ds0, hd, zd
	if (isbins .gt. 399999) write(6,100) isbins, s0, ds0, hd, zd
100	format(1x,'$$$ increase the dimension of att0 $$$  '
     #  ,i8,4(1pd10.3))
	delx=ds*st*cp
	dely=ds*st*sp
	delz=ds*ct

	fract=0.5d0
	x1=x0+delx*fract
	y1=y0+dely*fract
	z1=z0+delz*fract
	esum=0.0d0
	x2=x1
	y2=y1
	z2=z1
	do 500 is=1,isbins
	att0(is)=a(x2,y2,z2)
	esum=esum+att0(is)*ds
	x2=x2+delx
	y2=y2+dely
	z2=z2+delz
500	continue
	esum=esum-att0(1)*ds*fract
	tau= esum
	goto 2100
2000	continue
	tau= 0.d0
2100	continue
	return
	end
c
c---------------------------------------------------
        subroutine urad_cylinder(ufactor)
c---------------------------------------------------
c
	implicit real*8 (a-h, o-z)
        PARAMETER (NSIZEU=50000)
        double precision u(nsizeu),fa(nsizeu),f1(nsizeu)
        double precision fv(nsizeu,256)
        double precision uct(256),ust(256),ucp(256),usp(256)
        double precision ru(100),zu(100)
c       the dimension 256 is set also in the setuangle routines
	common/f/fs0,fs1,fsdum,fsurad,dsb1,dsb2
	common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar
	common/translate/pix2kpc
        common/ur/u,fa,fv,ures,rsizeu,zsizeu,albedo,iscat
        common/ur2/ru,zu,nru,nzu,iurad,nru_used,nzu_used
        common/urang/icang,nang,uct,ust,ucp,usp
        common/photom/away,bfactor,bin,ifilter
c
c     calculates radiation energy density for 
c     cylindrical symmetry and writes file with values to disk
c
c     resolution of grid in kpc
c
c     ufactor is multiplicative conversion factor such that spectral
c     energy densities are output in file dust.ur (unit 16) in units
c     of eV per cubic centimetres per Hz. Parameters output in other
c     files are for diagnostics only, and are in arbitrary units.
c
c     if IURAD is 1 a regular grid is calculated
c     if IURAD is 2 an irregular grid is calculated at positions
c                   supplied by the user.
c
      fs0_in=fs0
      fs1_in=fs1
      fs0=fsurad
      fs1=fsurad
c
      if(iurad .eq. 1)then 
       write(6,221)
 221   format(1h ,'URAD_CYLINDER (IURAD = 1:')
       write(6,222)iscat
 222   format(1h ,'iscat      = ',i6)
       write(6,23)ures
 23    format(1h ,'ures           = ',f8.2)
       write(6,223)rsizeu,zsizeu
 223   format(1h ,'rsizeu,zsizeu = ',2f8.2)
       write(6,24)icang
 24    format(1h ,'icang  = ',i6)
       write(6,288)ufactor
 288   format(1h ,'using ufactor = ',15e12.6)
      endif

      if(iurad .eq. 2)then 
       write(6,2221)
 2221  format(1h ,'URAD_CYLINDER (IURAD = 2:')
       write(6,222)iscat
       write(6,227)nru,nzu
 227   format(1h ,'nru,nzu = ',2i6)
       write(6,24)icang
       write(6,288)ufactor
      endif


c
c     grid of positions in r,z (both positive) to fill
c     cylinder on which Radtrans calculation has been performed
c

c     number of data map pixels per resolution element URES
c     on which the energy density is to be calculated
c
      if(iurad .eq. 1)then
       ures=ures/pix2kpc
       rsizeu=rsizeu/pix2kpc
       zsizeu=zsizeu/pix2kpc
c
       write(6,224)ures
 224   format(1h ,'ures (data pixels = ',f8.4)
c     dimensions of slice of galaxy in (r,z) within
c     the boundary cylinder for the integrations
       nr1=int(rhomax/ures)
       nz1=int(d/ures)
c     dimensions of slice of galaxy in (r,z) specified
c     by user for the calculation of energy densities
       nr2=int(rsizeu/ures)
       nz2=int(zsizeu/ures)
c     avoid calculations outside 
c     the boundary cylinder for the integrations
       if(nr2 .gt. nr1)then
        nr2=nr1
        write(6,289)
 289    format(1h ,'*** resetting R limit to within the cylinder')
       endif 
       if(nz2 .gt. nz1)then
        nz2=nz1
        write(6,299)
 299    format(1h ,'*** resetting Z limit to within the cylinder')
       endif
c
       nr=nr2
       nz=nz2
c
       write(6,9)nr,nz
 9     format(1h ,'nr nz      = ',2i6)
c
       ntot=nr*nz
       if(nsizeu .lt. ntot)then
        write(6,15)nsizeu,ntot
 15     format(1h ,'URAD_CYLINDER: aborting - set nsizeu higher:
     > nsizeu,ntot = ',2i10)
        goto 999
       endif 
      endif
      if(icang .lt. 1 .or. icang .gt. 3)then
      write(6,17)
 17   format(1h ,'URAD_CYLINDER: aborting - icang out of range')
      goto 999
      endif 

      iicang=icang
c      if(icang .eq. 1)call setuangle1
       if(icang .eq. 1)call setuangle(iicang)
       if(icang .eq. 2)call setuangle(iicang)
       if(icang .eq. 3)call setuangle(iicang)
c
      if(nang .gt. 256)then
      write(6,16)
 16   format(1h ,'URAD_CYLINDER: aborting - set nang lower')
      goto 999
      endif 

      x=0.0
      clight=1.0
      ii=0

c     r z u fa
      open (unit=14,file='dust.urad')
c     r z iang fv
      open (unit=15,file='dust.fvec')
c     u
      open (unit=16,file='dust.ur')
c     fa
      open (unit=17,file='dust.abs')
c     f1
      open (unit=18,file='dust.sca')

      if(iurad .eq. 2)then 
       nr=nru
       nz=nzu
      endif
c
      nru_used=nr
      nzu_used=nz
c
        if(ifilter.eq.1)xlam=912.
        if(ifilter.eq.2)xlam=1500.
        if(ifilter.eq.3)xlam=2800.
        if(ifilter.eq.4)xlam=3150.
        if(ifilter.eq.5)xlam=4430.
        if(ifilter.eq.6)xlam=5640.
        if(ifilter.eq.7)xlam=8090.
        if(ifilter.eq.8)xlam=12590.
        if(ifilter.eq.9)xlam=22000.
        if(ifilter.eq.10)xlam=2200.
        if(ifilter.eq.11)xlam=1350.
        if(ifilter.eq.12)xlam=1650.
        if(ifilter.eq.13)xlam=2000.
        if(ifilter.eq.14)xlam=2500.
        if(ifilter.eq.15)xlam=3650.

        nradtot=nr*nz
        write(14,*)nradtot,nr,nz,ifilter,xlam
c       no header to unit 16 as this is intended 
c       for display using plotresult.pro
c        write(16,*)nradtot,nr,nz,ifilter,xlam
        write(17,*)nradtot,nr,nz,ifilter,xlam
        write(18,*)nradtot,nr,nz,ifilter,xlam
      do 100 ir=1,nr
       if(iurad .eq. 1)y=ures*(float(ir)-1.0)
       if(iurad .eq. 2)y=ru(ir)
c       write(6,257)ir
c 257	format(1h ,'ir = ',i6)
       do 200 iz=1,nz
c       write(6,2557)iz
c 2557	format(1h ,'iz = ',i6)
        if(iurad .eq. 1)z=ures*(float(iz)-1.0)
        if(iurad .eq. 2)z=zu(iz)
        ii=ii+1
        u(ii)=0.0
        fa(ii)=0.0
        f1(ii)=0.0
          do 400 iang=1,nang
c       write(6,3559)iang
c 3559	format(1h ,'iang = ',i6)
           cp=ucp(iang)
           sp=usp(iang)
           ct=uct(iang)
           st=ust(iang)
c integrate distance si from back of galaxy to x,y,z along [ct,st,cp,sp]    
c           write(6,2789)ct,st,cp,sp
c 2789	   format(1h ,'call bound0 with ct,st,cp,sp = ',4f6.2)
           call bound0(x,y,z,x0,y0,z0,si,ct,st,cp,sp,ierr)
c          direct light
           flux0=br0(si,x0,y0,z0,ct,st,cp,sp)
c      if(y .eq. 0.0 .and. z .eq. 0.0)write(6,577)ct,st,y0,z0,si,flux0
c      if(y .eq. ures*(float(nr)-1.0) .and. z .eq. 0.0)write(6,577)
c     >ct,st,y0,z0,si,flux0
 577     format(1h ,5f10.2,e12.4)
           flux0=flux0/dfloat(nang)
           flux=flux0
           flux1=0.0
           if(iscat .ge. 1)then
            flux1=br1(si,x0,y0,z0,ct,st,cp,sp)
            flux1=flux1/dfloat(nang)
            call setzero(flux0,flux1,flux)
c            sumterm= 1.0d0
c            do 4500 iterm= 1, iscat
c	     if(flux0 .ne. 0.)sumterm=sumterm+(albedo*flux1/flux0)**iterm
c4500	    continue
c	    flux=flux0*sumterm
           endif
c
c      if(y .eq. 0.0 .and. z .eq. 0.0)write(6,578)flux0,flux1,flux
c 578  format(1h ,3e12.4)
c      if(y .eq. ures*(float(nr)-1.0) .and. z .eq. 0.0)write(6,*)
c     >flux0,flux1,flux
c           write(6,*)flux0,flux1,flux
           fv(ii,iang)=flux
           u(ii)=u(ii)+flux/clight
c           fa(ii)=fa(ii)+flux*ures*a(x,y,z)
c       absorbed energy in units of ergs/s/A/pc^3
c       found as uout*c(pc/s)*a(x,y,z)*0.001/pix2kpc
      fa(ii)=fa(ii)+(1.0d0-albedo)*
     >(flux/clight)*9.716*a(x,y,z)*1.0d-12/pix2kpc
c           f1(ii)=f1(ii)+flux1*ures*a(x,y,z)
           f1(ii)=f1(ii)+flux1/clight
c     if(y .eq. 0.0 .and. z .eq. 0.0)write(6,579)iang,ii,fv(ii,iang)
c579  format(1h ,2i8,e12.4)
 400      continue
c        write(6,2795)ir,iz,ii,u(ii)
c 2795   format(1h ,'ir,iz,ii,u(ii) = ',3i8,e15.6)
        u(ii)=u(ii)*ufactor    
        fa(ii)=fa(ii)*ufactor
c output y,z in units of 100pc    
        yout=y*pix2kpc*1000.0
        zout=z*pix2kpc*1000.0
c output u in units of erg/pc^3/A
c        factor=1.634d59 
c        factor=1.41106d58
        factor=1.41106d62
        factor=factor/(xlam*xlam)
        uout=u(ii)*factor
        faout=fa(ii)*factor
        write(14,275) yout,zout,uout,faout
 275    format(1h ,4e12.4)
c        write(6,275) yout,zout,uout,faout
c unit 16 in units of eV/cm^3/Hz
        write(16,*) u(ii) 
        write(17,*) fa(ii)
        write(18,*) f1(ii)
 200   continue
 100  continue
c
      close (unit=14)
      close (unit=15)
      close (unit=16)
      close (unit=17)
      close (unit=18)
c
      write(6,892)
 892  format(1h ,'r z u_rad u_abs written to dust.urad')
      write(6,893)
 893  format(1h ,'r z iang f_vec  written to dust.fvec')
      write(6,894)
 894  format(1h ,'u_rad written to dust.ur')
      write(6,895)
 895  format(1h ,'f_abs written to dust.abs')
      write(6,896)
 896  format(1h ,'f_1 written to dust.sca')
c
 999  continue
c
      fs0=fs0_in
      fs1=fs1_in
c
	return
	end
c---------------------------------------------------
        subroutine setuangle(imult)
c---------------------------------------------------
c
	implicit real*8 (a-h, o-z)
        double precision uct(256),ust(256),ucp(256),usp(256)
        double precision tp(3)
        common/urang/icang,nang,uct,ust,ucp,usp
c
c for imult=1 gives 26 vertices of elements of a 3x3 cube 
c             viewed from central element
c for imult=2 gives 98 vertices of elements of a 5x5 cube 
c             viewed from central element
c for imult=3 gives 218 vertices of elements of a 7x7 cube 
c             viewed from central element
c
      if(imult .le. 0 .or. imult .gt.3)then
      write(6,299)imult
 299  format(1h ,'SETUANGLE aborting: imult must be in range 0 to 3')
      goto 999
      endif
c
      nside=1+2*imult
c      nang=26
      nang=6*(nside-2)*(nside-2)
      nang=nang+12*(nside-2)
      nang=nang+8
      write(6,3991)nang
 3991 format(1h ,'SETUANGLE: number of directions = ',i6)
      sideoff=dfloat(nside/2)
      sideoff=sideoff+1.0
c      write(6,2898)imult,nside,sideoff
c 2898 format(1h ,'imult,nside,sideoff = ',2i6,e12.3)
      iang=0
      do 100 i=1,nside
      tp(1)=0.0d0
      if(i .eq. 1 .or. i .eq. nside)tp(1)=1.0d0
      do 200 j=1,nside
      tp(2)=0.0d0
      if(j .eq. 1 .or. j .eq. nside)tp(2)=1.0d0
      do 300 k=1,nside
      tp(3)=0.0d0
      if(k .eq. 1 .or. k .eq. nside)tp(3)=1.0d0
c     only take perimeter points
      if(tp(1).ne. 1.0 .and. tp(2).ne. 1.0 .and. tp(3).ne. 1.0)goto 288
      x=dfloat(i)-sideoff
      y=dfloat(j)-sideoff
      z=dfloat(k)-sideoff
      r=sqrt(x*x+y*y+z*z)
      if(r .le. 0.001)goto 288
      x=x/r
      y=y/r
      z=z/r
c
c test effect of random rotations:
c
      phirot=20.0
      thetrot=20.0
      pie=3.14159
      phirot=phirot*pie/180.0
c     next line added 08/07/13 (transormation had previously been omitted)
      thetrot=thetrot*pie/180.0
c      rotate about z axis by phirot
      xx=x*cos(phirot)-y*sin(phirot)
      yy=x*sin(phirot)+y*cos(phirot)
      zz=z
c     rotate about new yy axis
      xxx=xx*cos(thetrot)-zz*sin(thetrot)
      yyy=yy
      zzz=xx*sin(thetrot)+zz*cos(thetrot)
      x=xxx
      y=yyy
      z=zzz
c
      thet=acos(z)
c     thet in radians
c      write(6,357)thet
c357   format(1h ,'thet = ',f6.2)
c      write(6,279)i,j,k,x,y,z,thet
c 279	format(1h ,'i,j,k,x,y,z,thet = ',3i6,3f7.3,f7.3)
      st=sin(thet)
      ct=cos(thet)
      if(abs(st) .gt. 0.001)then
       cp=y/st
       sp=x/st
      endif
      if(abs(st) .le. 0.001)then
       cp=1.0
       sp=0.0
      endif
      iang=iang+1
c      write(6,*)iang
      uct(iang)=ct
      ust(iang)=st
      ucp(iang)=cp
      usp(iang)=sp
c      write(6,2799)iang,uct(iang),ust(iang),ucp(iang),usp(iang)      
c 2799 format(1h ,i6,4f10.4)
c      write(6,280)iang,i,j,k,x,y,z,uct(iang),
c     >ust(iang),ucp(iang),usp(iang) 
c 280	format(1h ,i7,3i3,3f6.3,f7.3,3f6.3)
c      write(6,29)
c 29   format(1h ,' ')
 288  continue
 300  continue
 200  continue
 100  continue
c
 999  continue
c
      return
      end
c
        Subroutine setzero(b0,b1,btotal)
c
c calculates the total light btotal along current ray as a
c function of input parameters b0 (direct light)
c and b1 (first order scattered light) and izero using the
c geometric series method of Kylafis. 
c
c In cases where there is little or no direct light 
c the geometric series in (b1/b0) the geometric
c expansion is limited to lower terms to avoid
c blowing up. If there is no zeroth order light
c the scattered light is set to the first order term
c
	implicit real*8 (a-h, o-z)
        parameter (nsizeu=50000)
        double precision u(nsizeu),fa(nsizeu)
        double precision fv(nsizeu,256)
        common/ur/u,fa,fv,ures,rsizeu,zsizeu,albedo,iscat
        common/config/izero,idust,itau,istars,imask,iall,
     >  shiftx,shifty,ifswitch,ifinns1,ifinns2,ifinnd1,ifinnd2
c
	sumterm= 1.0d0
c
        ib0flag=0
        if(izero .ge. 1)then
         test01=0.0
         xlim01=0.5
c
         if(b0 .ne. 0.)then
          test01=albedo*(b1/b0)
          if(dabs(test01) .le. xlim01)then
           xizeronew=float(izero)*(1.0-(dabs(test01)/xlim01))
           izeronew=int(xizeronew)
           if(izeronew .lt. 1)izeronew=1
	   do 4500 iterm= 1, izeronew
	    sumterm= sumterm+ (albedo* b1/ b0)** iterm
4500	   continue
c         dabs(test01) .le. xlim01
          endif 
c         if can't calculate High order terms take 1st order
          if(dabs(test01) .gt. xlim01)then
           sumterm=sumterm+(albedo* b1/ b0)
          endif
c        b0 .ne. 0.
         endif
c        if no direct light take first order scattered light 
         if(b0 .eq. 0.)then
          ib0flag=1
         endif
c
c       izero .ge. 1
        endif
c
c        write(6,8537)b0,b1,sumterm,albedo,btotal
c8537   format(1h ,'b0,b1,sumterm,albedo,btotal:',5e12.3)

	if(ib0flag .eq. 0)btotal= b0* sumterm
        if(ib0flag .eq. 1)btotal= b1* albedo
c
      return
      end



      subroutine plot_etabulger
c
c plots a radial profile of bulge emissivity and writes
c to a file
c
	implicit real*8 (a-h, o-z)
	common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
	common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  eta1arm,a1arm,narm,rup,fillarm
	common/atte/a0,eta0,rtruncate,sharp,athresh,rtrun,sha
	common/translate/pix2kpc

      open (unit=15, file='bulgeplotr.dat')
      write(6,221)
 221  format('open bulgeplotr.dat')
      x=0.0
      z=0.0
      eta1_in=eta1
      eta1=0.0
      eta0_in=eta0
      eta0=0.0
c plot from 0 to 11 * reff sampled at 0.01*reff
      write(6,78) reff,pix2kpc
 78   format(1h,'reff,pix2kpc=',2f8.2)
      step=reff*0.01
      write(6,79) step
 79   format(1h,'step',f8.4)
c      nstep=1+(7./pix2kpc)/step
      nstep=1+20. * reff/step
      write(15,*)nstep
      do 234 istep=1,nstep 
       y=step*(float(istep)-1.0)
       etaplot=eta(x,y,z)
       yout=y * pix2kpc
       write(15,*)yout,etaplot
234   continue
      eta1=eta1_in
      eta0=eta0_in
      close (unit=15)
      write(6,222)
 222  format('close bulgeplotr.dat')
c
      return
      end

      subroutine plot_etabulgez
c
c plots a vertical profile of bulge emissivity and writes
c to a file
c
	implicit real*8 (a-h, o-z)
	common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
	common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  eta1arm,a1arm,narm,rup,fillarm
	common/atte/a0,eta0,rtruncate,sharp,athresh,rtrun,sha
	common/translate/pix2kpc

      open (unit=16, file='bulgeplotz.dat')
      write(6,221)
 221  format('open bulgeplotz.dat')
      x=0.0
      y=0.0
      eta1_in=eta1
      eta1=0.0
      eta0_in=eta0
      eta0=0.0
c plot from 0 to 11 * reff sampled at 0.01*reff
      write(6,78) reff,pix2kpc
 78   format(1h,'reff,pix2kpc=',2f8.2)
      step=reff*0.01
      write(6,79) step
 79   format(1h,'step',f8.4)
c      nstep=1+(7./pix2kpc)/step
      nstep=1+10. * reff/step
      write(16,*)nstep
      do 234 istep=1,nstep 
       z=step*(float(istep)-1.0)
       etaplot=eta(x,y,z)
       zout=z * pix2kpc
       write(16,*)zout,etaplot
234   continue
      eta1=eta1_in
      eta0=eta0_in
      close (unit=16)
      write(6,222)
 222  format('close bulgeplotz.dat')
c
      return
      end

cc disk emissivity

      subroutine plot_etadiskr
c
c plots a radial profile of disk emissivity and writes
c to a file
c
	implicit real*8 (a-h, o-z)
	common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
	common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  eta1arm,a1arm,narm,rup,fillarm
	common/atte/a0,eta0,rtruncate,sharp,athresh,rtrun,sha
	common/translate/pix2kpc
        common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar

      open (unit=19, file='diskplotr.dat')
      write(6,221)
 221  format('open diskplotr.dat')
      x=0.0
      z=0.0
      eta1_in=eta1
      eta1=0.0
      etab_in=etab
      etab=0.0
c plot from 0 to 11 * reff sampled at 0.01*reff
      write(6,78) hs,pix2kpc
 78   format(1h,'hs,pix2kpc=',2f8.2)
      step=hs*0.01
      write(6,79) step
 79   format(1h,'step',f8.4)
c      nstep=1+(7./pix2kpc)/step
      nstep=1+20. * hs/step
      write(19,*)nstep
      do 234 istep=1,nstep 
       y=step*(float(istep)-1.0)
       etaplot=eta(x,y,z)
       yout=y * pix2kpc
       write(19,*)yout,etaplot
234   continue
      eta1=eta1_in
      etab=etab_in
      close (unit=19)
      write(6,222)
 222  format('close diskplotr.dat')
c
      return
      end

      subroutine plot_etadiskz
c
c plots a vertical profile of disk emissivity and writes
c to a file
c
	implicit real*8 (a-h, o-z)
	common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
	common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  eta1arm,a1arm,narm,rup,fillarm
	common/atte/a0,eta0,rtruncate,sharp,athresh,rtrun,sha
	common/translate/pix2kpc
        common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar

      open (unit=20, file='diskplotz.dat')
      write(6,221)
 221  format('open diskplotz.dat')
      x=0.0
      y=0.0
      eta1_in=eta1
      eta1=0.0
      etab_in=etab
      etab=0.0
c 
      write(6,78) hs,pix2kpc
 78   format(1h,'hs,pix2kpc=',2f8.2)
      step=hs*0.01
      write(6,79) step
 79   format(1h,'step',f8.4)
c      nstep=1+(7./pix2kpc)/step
      nstep=1+10. * hs/step
      write(20,*)nstep
      do 234 istep=1,nstep 
       z=step*(float(istep)-1.0)
       etaplot=eta(x,y,z)
       zout=z * pix2kpc
       write(20,*)zout,etaplot
234   continue
      eta1=eta1_in
      etab=etab_in
      close (unit=20)
      write(6,222)
 222  format('close diskplotz.dat')
c
      return
      end

c thin disk emissivity
cc disk emissivity

      subroutine plot_etadisk1r
c
c plots a radial profile of thin disk emissivity and writes
c to a file
c
	implicit real*8 (a-h, o-z)
	common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
	common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  eta1arm,a1arm,narm,rup,fillarm
	common/atte/a0,eta0,rtruncate,sharp,athresh,rtrun,sha
	common/translate/pix2kpc
        common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar

      open (unit=23, file='disk1plotr.dat')
      write(6,221)
 221  format('open disk1plotr.dat')
      x=0.0
      z=0.0
      eta0_in=eta0
      eta0=0.0
      etab_in=etab
      etab=0.0
c plot from 0 to 11 * reff sampled at 0.01*reff
      write(6,78) hs1,pix2kpc
 78   format(1h,'hs1,pix2kpc=',2f8.2)
      step=hs1*0.01
      write(6,79) step
 79   format(1h,'step',f8.4)
c      nstep=1+(7./pix2kpc)/step
      nstep=1+20. * hs1/step
      write(23,*)nstep
      do 234 istep=1,nstep 
       y=step*(float(istep)-1.0)
       etaplot=eta(x,y,z)
       yout=y * pix2kpc
       write(23,*)yout,etaplot
234   continue
      eta0=eta0_in
      etab=etab_in
      close (unit=23)
      write(6,222)
 222  format('close disk1plotr.dat')
c
      return
      end

      subroutine plot_etadisk1z
c
c plots a vertical profile of thin disk emissivity and writes
c to a file
c
	implicit real*8 (a-h, o-z)
	common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
	common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  eta1arm,a1arm,narm,rup,fillarm
	common/atte/a0,eta0,rtruncate,sharp,athresh,rtrun,sha
	common/translate/pix2kpc
        common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar

      open (unit=24, file='disk1plotz.dat')
      write(6,221)
 221  format('open disk1plotz.dat')
      x=0.0
      y=0.0
      eta0_in=eta0
      eta0=0.0
      etab_in=etab
      etab=0.0
c 
      write(6,78) zs1,pix2kpc
 78   format(1h,'zs1,pix2kpc=',2f8.2)
      step=zs1*0.01
      write(6,79) step
 79   format(1h,'step',f8.4)
c      nstep=1+(7./pix2kpc)/step
      nstep=1+10. * zs1/step
      write(24,*)nstep
      do 234 istep=1,nstep 
       z=step*(float(istep)-1.0)
       etaplot=eta(x,y,z)
       zout=z * pix2kpc
       write(24,*)zout,etaplot
234   continue
      eta0=eta0_in
      etab=etab_in
      close (unit=24)
      write(6,222)
 222  format('close disk1plotz.dat')
c
      return
      end

cc dust disk opacity

      subroutine plot_taudiskr
c
c plots a radial profile of disk opacity and writes
c to a file
c
	implicit real*8 (a-h, o-z)
	common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
	common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  eta1arm,a1arm,narm,rup,fillarm
	common/atte/a0,eta0,rtruncate,sharp,athresh,rtrun,sha
	common/translate/pix2kpc
        common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar

      open (unit=21, file='disktauplotr.dat')
      write(6,221)
 221  format('open disktauplotr.dat')
      x=0.0
      z=0.0
c      eta1_in=eta1
c      eta1=0.0
c      etab_in=etab
c      etab=0.0
c plot from 0 to 11 * reff sampled at 0.01*reff
      write(6,78) hd,pix2kpc
 78   format(1h,'hd,pix2kpc=',2f8.2)
      step=hd*0.01
      write(6,79) step
 79   format(1h,'step',f8.4)
c      nstep=1+(7./pix2kpc)/step
      nstep=1+20. * hd/step
      write(21,*)nstep
      do 234 istep=1,nstep 
       y=step*(float(istep)-1.0)
       tauplot=a(x,y,z)
c       write(6,80)tauplot
c 80    format(1h,'tauplot=',d10.4)
       yout=y * pix2kpc
       write(21,*)yout,tauplot
234   continue
c      eta1=eta1_in
c      etab=etab_in
      close (unit=21)
      write(6,222)
 222  format('close disktauplotr.dat')
c
      return
      end

      subroutine plot_taudiskz
c
c plots a vertical profile of disk opacity and writes
c to a file
c
	implicit real*8 (a-h, o-z)
	common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
	common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  eta1arm,a1arm,narm,rup,fillarm
	common/atte/a0,eta0,rtruncate,sharp,athresh,rtrun,sha
	common/translate/pix2kpc
        common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar

      open (unit=22, file='disktauplotz.dat')
      write(6,221)
 221  format('open disktauplotz.dat')
      x=0.0
      y=0.0
c      eta1_in=eta1
c      eta1=0.0
c      etab_in=etab
c      etab=0.0
c 
      write(6,78) zd,pix2kpc
 78   format(1h,'zd,pix2kpc=',2f8.2)
      step=zd*0.01
      write(6,79) step
 79   format(1h,'step',f8.4)
c      nstep=1+(7./pix2kpc)/step
      nstep=1+10. * zd/step
      write(22,*)nstep
      do 234 istep=1,nstep 
       z=step*(float(istep)-1.0)
       tauplot=a(x,y,z)
       zout=z * pix2kpc
       write(22,*)zout,tauplot
234   continue
c      eta1=eta1_in
c      etab=etab_in
      close (unit=22)
      write(6,222)
 222  format('close etatauplotz.dat')
c
      return
      end
