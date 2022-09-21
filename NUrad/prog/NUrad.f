        program NUrad
c urad.f updated to include 4 extra disk components - JJT 3 Oct 17
c added inner truncation - JJT 23 Nov 17
c decoupled truncation radii - JJT 15 Mar 18
c added new young stellar component for nuclear disc of M33 - JJT 3 Aug 18
        Implicit real*8 (a-h, o-z)

        parameter (NSIZE=45000000)
        parameter (nsizeu=50000)
        double precision DZMODEL(NSIZE)
        double precision u(nsizeu),fa(nsizeu)
        double precision fv(nsizeu,256)
        double precision uct(256),ust(256),ucp(256),usp(256)
        double precision ru(100),zu(100)
c new param from gal_param ccc
        double precision dist,inclination,pix_size,ellipse,conver
        integer mask_opt,nx_b,ny_b,nx_n,ny_n,nx_i,ny_i,nx_m,ny_m,
     >               nx_o,ny_o,mstep1_b,mstep2_b,mstep3_b,mlength1_b,
     >               mlength2_b,
     >               mask_opt_b,mask_opt_n,mask_opt_i,mask_opt_m,
     >               mask_opt_o,mstep1_n,mstep2_n,mstep3_n,mlength1_n,
     >               mlength2_n,
     >               mstep1_i,mstep2_i,mstep3_i,mlength1_i,mlength2_i,
     >               mstep1_m,mstep2_m,mstep3_m,mlength1_m,mlength2_m,
     >               mstep1_o,mstep2_o,mstep3_o,mlength1_o,mlength2_o,
     >               elm,elm_n,elm_b,elm_i,elm_m,elm_o,iinclination
        integer diag_plots
c
        integer mask(nsize), Narg
        character*80 morph,filter, dustmodel,
     >               filename_dustmodel,sizero,
     >               stau,snsersic,sabssca,path1,path2,filename_uradin,
     >               path,ss,filename_geometry,sdummy,sdiag,
     >               shd,szd,shd1,szd1,shs,szs,shs1,szs1,sreff,sellipt,
     >               sinclination
        character*160 filename_uradd,filename_uradtd,filename_uradb,
     >               filename_image,filename_gal

c      use NUrad_variables
        common/incli/theta,tau0
        common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >  zsin,zssolar,hssolar,zdin,zdsolar,hdsolar,
     >   xis0,xid0,hstin,hdtin
        common/f/fs0,fs1,fsdum,fsurad,dsb1,dsb2
        common/atte/a0,eta0,rtruncate,sharp,athresh,rtrun,sha,
     >   rtruncated,sharpd,rtrund,shad
        common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
        common/ferrers/etaf,af,bf,cf,thetaf,ctf,stf,iferrers
        common/ellexp/eta2,zs2,hsa2,hsb2,thetae,cte,ste,
     >  hda2,hdb2,zd2,tau2,a2,iellexp
        common/bcons/etabmax,rin,ibulge
        common/g/g,g2,cosalp
        common/integ/k,m,n,km
        common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  xis1,xid1,hs1tin,hd1tin,
     >  eta1arm,a1arm,narm,rup,fillarm,
     >   rtruncated1,sharpd1,rtrund1,shad1
        common/translate/pix2kpc
        common/inout/nx,ny,pixsize,nimage
        common/ur/u,fa,fv,ures,rsizeu,zsizeu,albedo,iscat
        common/ur2/ru,zu,nru,nzu,iurad,nru_used,nzu_used
        common/urang/icang,nang,uct,ust,ucp,usp
        common/photom/away,bfactor,bin,ifilter
        common/config/izero,idust,itau,istars,imask,iall,
     >  shiftx,shifty,ifswitch

ccc JJT 3 Oct 17 Common blocks ccc
c extra disks for inner region
c disk
	common/eindisk3/eta3,zs3,hs3,zd3,hd3,a3,tau3,idisk3,
     >  rtruncate3,sharp3,rtrun3,sha3,hd3in,zd3in,hs3in,zs3in,
     >  hd3solar,zd3solar,hs3solar,zs3solar,
     >  xis3,xid3,hs3tin,hd3tin,
     >  eta3arm,a3arm,narm3,rup3,fillarm3,
     >  rtruncated3,sharpd3,rtrund3,shad3
c tdisk
	common/eindisk4/eta4,zs4,hs4,zd4,hd4,a4,tau4,idisk4,
     >  rtruncate4,sharp4,rtrun4,sha4,hd4in,zd4in,hs4in,zs4in,
     >  hd4solar,zd4solar,hs4solar,zs4solar,
     >  xis4,xid4,hs4tin,hd4tin,
     >  eta4arm,a4arm,narm4,rup4,fillarm4,
     >  rtruncated4,sharpd4,rtrund4,shad4

c extra disks for outer regions
c disk
	common/eoutdisk5/eta5,zs5,hs5,zd5,hd5,a5,tau5,idisk5,
     >  rtruncate5,sharp5,rtrun5,sha5,hd5in,zd5in,hs5in,zs5in,
     >  hd5solar,zd5solar,hs5solar,zs5solar,
     >  xis5,xid5,hs5tin,hd5tin,
     >  eta5arm,a5arm,narm5,rup5,fillarm5,
     >  rtruncated5,sharpd5,rtrund5,shad5
c tdisk
	common/eoutdisk6/eta6,zs6,hs6,zd6,hd6,a6,tau6,idisk6,
     >  rtruncate6,sharp6,rtrun6,sha6,hd6in,zd6in,hs6in,zs6in,
     >  hd6solar,zd6solar,hs6solar,zs6solar,
     >  xis6,xid6,hs6tin,hd6tin,
     >  eta6arm,a6arm,narm6,rup6,fillarm6,
     >  rtruncated6,sharpd6,rtrund6,shad6
ccc JJT 3 Oct 17 - Common blocks ccc
ccc JJT 3 Aug 18 - Common blocks ccc
c nulear comp ( young stellar disc
	common/ntdisk7/eta7,zs7,hs7,a7,idisk7,
     >  rtruncate7,sharp7,rtrun7,sha7,hs7in,zs7in,
     >  hs7solar,zs7solar,
     >  xis7,hs7tin,
     >  eta7arm,a7arm,narm7,rup7,fillarm7
ccc
c dust.urad common JJT 19FEB19
            common/dust_urad/yout(nsizeu),zout(nsizeu),
     >  uout(nsizeu),faout(nsizeu),bout(nsizeu)
            common/diag/prof_lim
c
      open(unit=21,file='path.in',status='old')
       read(21,800)path1
       read(21,800)path2
      close(unit=21)
      lenpath1=len_trim(trim(path1))
      lenpath2=len_trim(trim(path2))
c
c      open(unit=1,file='param.in',status='old')
c
c galaxy parameter read in
       filename_gal=
     >      path1(1:lenpath1)//path2(1:lenpath2)//
     >      'indata/'//
     >      'gal_param.in'
       open(unit=41, file=filename_gal,status='old')
       read(41,*) ss
       read(41,*) diag_plots
        print*, diag_plots
       read(41,*) ss
	    read(41,*) prof_lim
       read(41,*) ss
       read(41,*) izero
       read(41,*) ss
       read(41,*) dustmodel
       read(41,*) ss
       read(41,*) dist
       read(41,*) ss
       read(41,*) inclination
       read(41,*) ss
       read(41,*) ss
       print*,ss
       read(41,*) pix_size
       print*,pix_size
       read(41,*) ss
       read(41,*) nx_b
       read(41,*) ss
       read(41,*) ny_b
       read(41,*) ss
       read(41,*) nx_n
       read(41,*) ss
       read(41,*) ny_n
       read(41,*) ss
       read(41,*) nx_i
       read(41,*) ss
       read(41,*) ny_i
       read(41,*) ss
       read(41,*) nx_m
       read(41,*) ss
       read(41,*) ny_m
       read(41,*) ss
       read(41,*) nx_o
       read(41,*) ss
       read(41,*) ny_o
       read(41,*) ss
       read(41,*) ss
       read(41,*) mask_opt_b
       read(41,*) ss
       read(41,*) mask_opt_n
       read(41,*) ss
       read(41,*) mask_opt_i
       read(41,*) ss
       read(41,*) mask_opt_m
       read(41,*) ss
       print*,ss
       read(41,*) mask_opt_o
       print*,mask_opt_o
       read(41,*) ss
       read(41,*) ss
       read(41,*) mstep1_b
       read(41,*) ss
       read(41,*) mlength1_b
       read(41,*) ss
       read(41,*) mstep2_b
       read(41,*) ss
       read(41,*) mlength2_b
       read(41,*) ss
       read(41,*) mstep3_b
       read(41,*) ss
       read(41,*) ss
       read(41,*) mstep1_n
       read(41,*) ss
       read(41,*) mlength1_n
       read(41,*) ss
       read(41,*) mstep2_n
       read(41,*) ss
       read(41,*) mlength2_n
       read(41,*) ss
       read(41,*) mstep3_n
       read(41,*) ss
       read(41,*) ss
       read(41,*) mstep1_i
       read(41,*) ss
       read(41,*) mlength1_i
       read(41,*) ss
       read(41,*) mstep2_i
       read(41,*) ss
       read(41,*) mlength2_i
       read(41,*) ss
       read(41,*) mstep3_i
       read(41,*) ss
       read(41,*) ss
       read(41,*) mstep1_m
       read(41,*) ss
       read(41,*) mlength1_m
       read(41,*) ss
       read(41,*) mstep2_m
       read(41,*) ss
       read(41,*) mlength2_m
       read(41,*) ss
       read(41,*) mstep3_m
       read(41,*) ss
       read(41,*) ss
       read(41,*) mstep1_o
       read(41,*) ss
       read(41,*) mlength1_o
       read(41,*) ss
       read(41,*) mstep2_o
       read(41,*) ss
       read(41,*) mlength2_o
       read(41,*) ss
       print*,ss
       read(41,*) mstep3_o
       read(41,*) ss
       read(41,*) ss
       read(41,*) elm_b
       read(41,*) ss
       read(41,*) elm_n
       read(41,*) ss
       read(41,*) elm_i
       read(41,*) ss
       read(41,*) elm_m
       read(41,*) ss
       read(41,*) elm_o
       close(unit=41)
c former photometry.in
c Distance of galaxy (in Mpc)
c rtruncate, sharp (truncation radius of exponential disk in kpc)
c                  sigma of cutoff = sharp * rtruncate 
c Luminosity scaling factor
c Binning factor
c
c
       away=dist

c
c read in steering parameters for calculation
c
c      READ(1,800)morph

c read in steering parameters for calculation from command line
c thereby avoiding the need for a param.in file.
      Narg=iargc()
      if (Narg < 1) then
      write(*,'(a)')"usage: ./foo op1 op2"
      write(*,'(a)')"op1=morph (b/d/td/di/tdi/do/tdo/ntd)"
c      write(*,'(a)')"op2=izero (0=abs, 10=sca)"
      write(*,'(a)')"op2=filter"
c      write(*,'(a)')"op4=dust model (wd01_q10/pop00_q01)"
c      write(*,'(a)')"op4=diagnostic run (0=no, 1=yes)"
      else      
      call getarg(1,morph)
c      call getarg(2,sizero)
      call getarg(2,filter)
c      call getarg(4,dustmodel)
c      call getarg(5,sdiag)
c      read(sizero,*) izero
c      read(sdiag,*) diag_plots
      endif
      lenmorph=len_trim(trim(morph))
c      write(6,2389)lenmorph
c2389  format(1h ,'lenmorph = ',i10)
      morph=morph(1:lenmorph)
c      if(morph.ne.'d'.and.morph.ne.'b'.and.morph.ne.'td')then
ccc JJT 3 Oct 17 ccc
      if(morph.ne.'d'.and.morph.ne.'b'.and.morph.ne.'td'.and.
     >  morph.ne.'di'.and.morph.ne.'do'.and.
     >  morph.ne.'tdi'.and.morph.ne.'tdo' .and.
     >  morph .ne. 'ntd' )then
ccc JJT 3 Oct 17 ccc
       write(6,8270)
 8270  format(1h ,'URAD: morph incorrectly specified - aborting:')
       write(6,2568)morph
 2568  format(1h ,a10)
       goto 999
      endif



c
c     set up geometry parameters according to Tuffs et al. 2004 which
c     are quoted with reference to stellar scale length in B-band in kpc
c

       filename_geometry=
     >path1(1:lenpath1)//path2(1:lenpath2)//
     >'indata/'//
c    >'geometry.in'
ccc JJT 3 Oct 17 ccc
     >'geometry.in'
ccc JJT 3 Oct 17 ccc 
       print *, filename_geometry
       open(unit=31, file=filename_geometry,status='old')
       read(31,*) ss
       read(31,*) totaltau1
       read(31,*) ss
       read(31,*) totaltau2
       read(31,*) ss
       read(31,*) hd
       read(31,*) ss
       read(31,*) zd
       read(31,*) ss
       read(31,*) hdin
       read(31,*) ss
       read(31,*) zdin
       read(31,*) ss
       read(31,*) hdsolar
       read(31,*) ss
       read(31,*) zdsolar
       read(31,*) ss
       read(31,*) hd1
       read(31,*) ss
       read(31,*) zd1
       read(31,*) ss
       read(31,*) hd1in
       read(31,*) ss
       read(31,*) zd1in
       read(31,*) ss
       read(31,*) hd1solar
       read(31,*) ss
       read(31,*) zd1solar
       read(31,*) ss
       read(31,*) h_bdisk
       read(31,*) ss
       read(31,*) h_vdisk
       read(31,*) ss
       read(31,*) h_idisk
       read(31,*) ss
       read(31,*) h_jdisk
       read(31,*) ss
       read(31,*) h_kdisk
       read(31,*) ss
       read(31,*) h_ir36disk
       read(31,*) ss
       read(31,*) h_ir45disk
       read(31,*) ss
       read(31,*) h_ir58disk
       read(31,*) ss
       read(31,*) zs
       read(31,*) ss
       read(31,*) hsin
       read(31,*) ss
       read(31,*) zsin
       read(31,*) ss
       read(31,*) hssolar
       read(31,*) ss
       read(31,*) zssolar
       read(31,*) ss
       read(31,*) hs1
       read(31,*) ss
       read(31,*) zs1
       read(31,*) ss
       read(31,*) hs1in
       read(31,*) ss
       read(31,*) zs1in
       read(31,*) ss
       read(31,*) hs1solar
       read(31,*) ss
       read(31,*) zs1solar
       read(31,*) ss
       read(31,*) rtruncate
       print*,'test'
       print*,ss
       print*,rtruncate
       read(31,*) ss
       read(31,*) sharp
       read(31,*) ss
       read(31,*) rtruncated
       print*,'test'
       print*,ss
       print*,rtruncated
       read(31,*) ss
       read(31,*) sharpd
       read(31,*) ss
       read(31,*) rtruncate1
       print*,'test'
       print*,ss
       print*,rtruncate1
       read(31,*) ss
       read(31,*) sharp1
       read(31,*) ss
       read(31,*) rtruncated1
       print*,'test'
       print*,ss
       print*,rtruncated1
       read(31,*) ss
       read(31,*) sharpd1
       read(31,*) ss
       read(31,*) reff
       read(31,*) ss
       read(31,*) ellipt
       read(31,*) ss
       read(31,*) nsersic
       read(31,*) ss
       read(31,*) xis0
       read(31,*) ss
       read(31,*) xis1
       read(31,*) ss
       read(31,*) xid0
       read(31,*) ss
       read(31,*) xid1
       read(31,*) ss
       read(31,*) idisk1
       read(31,*) ss
       read(31,*) idisk2
ccc JJT 3 Oct 17 ccc
c Read in parameters for new components
c Inner component
       read(31,*) ss
       read(31,*) ss
       read(31,*) totaltau3
       read(31,*) ss
       read(31,*) totaltau4
       read(31,*) ss
       read(31,*) hd3
       read(31,*) ss
       read(31,*) zd3
       read(31,*) ss
       read(31,*) hd3in
       read(31,*) ss
       read(31,*) zd3in
       read(31,*) ss
       read(31,*) hd3solar
       read(31,*) ss
       read(31,*) zd3solar
       read(31,*) ss
       read(31,*) hd4
       read(31,*) ss
       read(31,*) zd4
       read(31,*) ss
       read(31,*) hd4in
       read(31,*) ss
       read(31,*) zd4in
       read(31,*) ss
       read(31,*) hd4solar
       read(31,*) ss
       read(31,*) zd4solar
       read(31,*) ss
       read(31,*) h_bdisk3
       read(31,*) ss
       read(31,*) h_vdisk3
       read(31,*) ss
       read(31,*) h_idisk3
       read(31,*) ss
       read(31,*) h_jdisk3
       read(31,*) ss
       read(31,*) h_kdisk3
       read(31,*) ss
       read(31,*) h_ir36disk3
       read(31,*) ss
       read(31,*) h_ir45disk3
       read(31,*) ss
       read(31,*) h_ir58disk3
       read(31,*) ss
       read(31,*) zs3
       read(31,*) ss
       read(31,*) hs3in
       read(31,*) ss
       read(31,*) zs3in
       read(31,*) ss
       read(31,*) hs3solar
       read(31,*) ss
       read(31,*) zs3solar
       read(31,*) ss
       read(31,*) hs4
       read(31,*) ss
       read(31,*) zs4
       read(31,*) ss
       read(31,*) hs4in
       read(31,*) ss
       read(31,*) zs4in
       read(31,*) ss
       read(31,*) hs4solar
       read(31,*) ss
       read(31,*) zs4solar
       read(31,*) ss
       read(31,*) rtruncate3
       print*,'test'
       print*,ss
       print*,rtruncate3
       read(31,*) ss
       read(31,*) sharp3
       read(31,*) ss
       read(31,*) rtruncated3
       print*,'test'
       print*,ss
       print*,rtruncated3
       read(31,*) ss
       read(31,*) sharpd3
       read(31,*) ss
       read(31,*) rtruncate4
       print*,'test'
       print*,ss
       print*,rtruncate4
       read(31,*) ss
       read(31,*) sharp4
       read(31,*) ss
       read(31,*) rtruncated4
       print*,'test'
       print*,ss
       print*,rtruncated4
       read(31,*) ss
       read(31,*) sharpd4
       read(31,*) ss
       read(31,*) xis3
       read(31,*) ss
       read(31,*) xis4
       read(31,*) ss
       read(31,*) xid3
       read(31,*) ss
       read(31,*) xid4
       read(31,*) ss
       read(31,*) idisk3
       read(31,*) ss
       read(31,*) idisk4

c Outer component
       read(31,*) ss
       read(31,*) ss
       read(31,*) totaltau5
       read(31,*) ss
       read(31,*) totaltau6
       read(31,*) ss
       read(31,*) hd5
       read(31,*) ss
       read(31,*) zd5
       read(31,*) ss
       read(31,*) hd5in
       read(31,*) ss
       read(31,*) zd5in
       read(31,*) ss
       read(31,*) hd5solar
       read(31,*) ss
       read(31,*) zd5solar
       read(31,*) ss
       read(31,*) hd6
       read(31,*) ss
       read(31,*) zd6
       read(31,*) ss
       read(31,*) hd6in
       read(31,*) ss
       read(31,*) zd6in
       read(31,*) ss
       read(31,*) hd6solar
       read(31,*) ss
       read(31,*) zd6solar
       read(31,*) ss
       read(31,*) h_bdisk5
       read(31,*) ss
       read(31,*) h_vdisk5
       read(31,*) ss
       read(31,*) h_idisk5
       read(31,*) ss
       read(31,*) h_jdisk5
       read(31,*) ss
       read(31,*) h_kdisk5
       read(31,*) ss
       read(31,*) h_ir36disk5
       read(31,*) ss
       read(31,*) h_ir45disk5
       read(31,*) ss
       read(31,*) h_ir58disk5
       read(31,*) ss
       read(31,*) zs5
       read(31,*) ss
       read(31,*) hs5in
       read(31,*) ss
       read(31,*) zs5in
       read(31,*) ss
       read(31,*) hs5solar
       read(31,*) ss
       read(31,*) zs5solar
       read(31,*) ss
       write(6,*) ss
       read(31,*) hs6
       write(6,*) hs6
       read(31,*) ss
       read(31,*) zs6
       read(31,*) ss
       read(31,*) hs6in
       read(31,*) ss
       read(31,*) zs6in
       read(31,*) ss
       read(31,*) hs6solar
       read(31,*) ss
       read(31,*) zs6solar
       read(31,*) ss
       read(31,*) rtruncate5
       print*,'test'
       print*,ss
       print*,rtruncate5
       read(31,*) ss
       read(31,*) sharp5
       read(31,*) ss
       read(31,*) rtruncated5
       print*,'test'
       print*,ss
       print*,rtruncated5
       read(31,*) ss
       read(31,*) sharpd5
       read(31,*) ss
       read(31,*) rtruncate6
       print*,'test'
       print*,ss
       print*,rtruncate6
       read(31,*) ss
       read(31,*) sharp6
       read(31,*) ss
       read(31,*) rtruncated6
       print*,'test'
       print*,ss
       print*,rtruncated6
       read(31,*) ss
       read(31,*) sharpd6
       read(31,*) ss
       read(31,*) xis5
       read(31,*) ss
       read(31,*) xis6
       read(31,*) ss
       read(31,*) xid5
       read(31,*) ss
       read(31,*) xid6
       read(31,*) ss
       read(31,*) idisk5
       read(31,*) ss
       read(31,*) idisk6
cc JJT 8 Oct 17 ccc
c read in inner truncation radii
       read(31,*) ss
       read(31,*) ss
       read(31,*) hstin
       read(31,*) ss
       read(31,*) hdtin
       read(31,*) ss
       read(31,*) hs1tin
       read(31,*) ss
       read(31,*) hd1tin
       read(31,*) ss
       read(31,*) hs3tin
       read(31,*) ss
       read(31,*) hd3tin
       read(31,*) ss
       read(31,*) hs4tin
       read(31,*) ss
       read(31,*) hd4tin
       read(31,*) ss
       read(31,*) hs5tin
       read(31,*) ss
       read(31,*) hd5tin
       read(31,*) ss
       read(31,*) hs6tin
       read(31,*) ss
       read(31,*) hd6tin

ccc JJT 3 Oct 17 ccc    

ccc JJT 3 Aug 18 ccc
c Nuclear component
       read(31,*) ss
       write(6,*) ss
       read(31,*) ss
       write(6,*) ss
       read(31,*) hs7
       write(6,*) hs7
       read(31,*) ss
       read(31,*) zs7
       read(31,*) ss
       read(31,*) hs7in
       read(31,*) ss
       read(31,*) zs7in
       read(31,*) ss
       read(31,*) hs7solar
       read(31,*) ss
       read(31,*) zs7solar
       read(31,*) ss
       read(31,*) rtruncate7
       read(31,*) ss
       read(31,*) sharp7
       read(31,*) ss
       read(31,*) xis7
       read(31,*) ss
       read(31,*) idisk7
       print*,ss
       print*,idisk7
       read(31,*) ss
       read(31,*) hs7tin

       close(unit=31)
c      READ(1,700) totaltau1
c     write(6,3909)totaltau1
c3909 format(1h ,'totaltau1 = ',f10.2) 
c      if(totaltau1.lt.0..or.totaltau1.gt.10.)then
      if(totaltau1.lt.0)then
      write(6,8269)
8269  format(1h ,'URAD: totaltau1 out of range - aborting:')
c      write(6,2567)totaltau1
c2567  format(1h ,f15.2)
       goto 999
      endif
c      READ(1,700) totaltau2
c     write(6,3909)totaltau2
c3909 format(1h ,'totaltau2 = ',f10.2) 
c      if(totaltau2.lt.0..or.totaltau2.gt.10.)then
      if(totaltau2.lt.0)then
      write(6,8272)
 8272 format(1h ,'URAD: totaltau2 out of range - aborting:')
c      write(6,2567)totaltau2
c2567  format(1h ,f15.2)
       goto 999
      endif

ccc JJT 3 Oct 17 ccc
c error message for negative tau values
      if(totaltau3.lt.0)then
      write(6,8269)
c      write(6,2567)totaltau1
c2567  format(1h ,f15.2)
       goto 999
      endif
c      READ(1,700) totaltau2
c     write(6,3909)totaltau2
c3909 format(1h ,'totaltau2 = ',f10.2) 
c      if(totaltau2.lt.0..or.totaltau2.gt.10.)then
      if(totaltau4.lt.0)then
      write(6,8272)
c      write(6,2567)totaltau2
c2567  format(1h ,f15.2)
       goto 999
      endif
      if(totaltau5.lt.0)then
      write(6,8269)
c      write(6,2567)totaltau1
c2567  format(1h ,f15.2)
       goto 999
      endif
c      READ(1,700) totaltau2
c     write(6,3909)totaltau2
c3909 format(1h ,'totaltau2 = ',f10.2) 
c      if(totaltau2.lt.0..or.totaltau2.gt.10.)then
      if(totaltau6.lt.0)then
      write(6,8272)
c      write(6,2567)totaltau2
c2567  format(1h ,f15.2)
       goto 999
      endif
ccc JJT 3 Oct 17 ccc
c      totaltau = totaltau1 + totaltau2
ccc JJT 3 Oct 17 ccc
      totaltau = totaltau1 + totaltau2 + totaltau3 + totaltau4 +
     > totaltau5 + totaltau6
ccc JJT 3 Oct 17 ccc
c     create string stau corresponding to 10*tautotal
      itotaltau=int(totaltau*10.)
      ihd = int(hd *1000.)
      izd = int(zd * 1000.)
      ihd1 = int(hd1 * 1000.)
      izd1 = int(zd1 * 1000.)
      ihs = int(h_bdisk * 1000.)
      izs = int(zs * 1000.)
      ihs1 = int(hs1 * 1000.)
      izs1 = int(zs1 * 1000.)
      ireff = int(reff * 1000.)
      iellipt = int(ellipt * 100.)
      iinclination = int(inclination)

!!! Create strings
10    format(I9)
      write(sdummy,*) itotaltau
      stau = adjustl(sdummy)
      write(sdummy,*) ihd
      shd = adjustl(sdummy)
      write(sdummy,*) izd
      szd = adjustl(sdummy)      
      write(sdummy,*) ihd1
      shd1 = adjustl(sdummy)
      write(sdummy,*) izd1
      szd1 = adjustl(sdummy)
      write(sdummy,*) ihs
      shs = adjustl(sdummy)
      write(sdummy,*) izs
      szs = adjustl(sdummy)
      write(sdummy,*) ihs1
      shs1 = adjustl(sdummy)
      write(sdummy,*) izs1
      szs1 = adjustl(sdummy)
      write(sdummy,*) ireff
      sreff = adjustl(sdummy)
      write(sdummy,*) iellipt
      sellipt =adjustl(sdummy)
      write(sdummy,*) iinclination
      sinclination = adjustl(sdummy)

c     string lengths
      lenstau=len_trim(trim(stau))
      lenshd=len_trim(trim(shd))
      lenszd=len_trim(trim(szd))
      lenshd1=len_trim(trim(shd1))
      lenszd1=len_trim(trim(szd1))
      lenshs=len_trim(trim(shs))
      lenszs=len_trim(trim(szs))
      lenshs1=len_trim(trim(shs1))
      lenszs1=len_trim(trim(szs1))
      lensreff=len_trim(trim(sreff))
      lensellipt=len_trim(trim(sellipt))
      lensinclination=len_trim(trim(sinclination))

c      READ(1,600) izero
      if(izero.lt.0.or.izero.gt.10)then
       write(6,8268)
 8268  format(1h ,'URAD: izero out of range - aborting:')
       write(6,2566)izero
 2566  format(1h ,i10)
       goto 999
      endif
      write(6,601)izero
601   format(1h ,'izero = ',i6)
c
c      READ(1,800) filter
      lenfilter=len_trim(trim(filter))
      filter=filter(1:lenfilter)
c
c      READ(1,800) dustmodel
      lendustmodel=len_trim(trim(dustmodel))
      dustmodel=dustmodel(1:lendustmodel)
c      if(dustmodel.ne.'wd01_q06'.and.dustmodel.ne.'pop00_q01')then
c       write(6,8271)
c 8271  format(1h ,'URAD: dust model incorrectly specified - aborting:')
c       write(6,2569)dustmodel
c 2569  format(1h ,a10)
c       goto 999
c      endif
c
      filename_dustmodel=
     >path1(1:lenpath1)//path2(1:lenpath2)//
     >'indata/'//
     >'dustmodel_'//dustmodel(1:lendustmodel)//'.in'
      write(6,2079)filename_dustmodel
2079  format(1h ,'opening ',a) 
      open(unit=2,file=filename_dustmodel,status='old')
c
600   FORMAT(i10)
700   format(f10.2)  
800   format(a80)
c
c set up filter identifier and check filter recognised
c
      iflag=0
      if(filter.eq.'uv09')then
       ifilter=1
       wave=912.
       iflag=iflag+1
      endif
      if(filter.eq.'uv15')then
       ifilter=2
       wave=1500.
       iflag=iflag+1
      endif
      if(filter.eq.'uv28')then
       ifilter=3
       wave=2800.
       iflag=iflag+1
      endif
      if(filter.eq.'uv31')then
       ifilter=4
       wave=3150.
       iflag=iflag+1
      endif
      if(filter.eq.'b')then
       ifilter=5
       wave=4430.
       iflag=iflag+1
      endif
      if(filter.eq.'v')then
       ifilter=6
       wave=5640.
       iflag=iflag+1
      endif
      if(filter.eq.'i')then
       ifilter=7
       wave=8090.
       iflag=iflag+1
      endif
      if(filter.eq.'j')then
       ifilter=8
       wave=12590.
       iflag=iflag+1
      endif
      if(filter.eq.'k')then
       ifilter=9
       wave=22000.
       iflag=iflag+1
      endif
      if(filter.eq.'uv22')then
       ifilter=10
       wave=2200.
       iflag=iflag+1
      endif
      if(filter.eq.'uv13')then
       ifilter=11
       wave=1350.
       iflag=iflag+1
      endif
      if(filter.eq.'uv16')then
       ifilter=12
       wave=1650.
       iflag=iflag+1
      endif
      if(filter.eq.'uv20')then
       ifilter=13
       wave=2000.
       iflag=iflag+1
      endif
      if(filter.eq.'uv25')then
       ifilter=14
       wave=2500.
       iflag=iflag+1
      endif
      if(filter.eq.'uv36')then
       ifilter=15
       wave=3650.
       iflag=iflag+1
      endif
cJJT 30OCT18
      if(filter.eq.'ir36')then
       ifilter=16
c JJT set correct irac1 wavelength
       wave=36000.
       iflag=iflag+1
      endif
      if(filter.eq.'ir45')then
       ifilter=17
       wave=45000.
       iflag=iflag+1
      endif
      if(filter.eq.'ir58')then
       ifilter=18
       wave=58000.
       iflag=iflag+1
      endif
c
      if(iflag.ne.1)then
       write(6,8267)
 8267  format(1h ,'URAD: filter not found - aborting:')
       write(6,2565)filter
 2565  format(1h ,a10)
       goto 999
      endif
c     
c read in optical constants and set albedo, g factor and tauratio
c for chosen band
c
      write(6,2892)wave
2892  format(1h ,'wave = ',f10.2)
c
      iflag=0
      READ(2,600)noc
      do 2563 i=1,noc
       read(2,*)dum1,dum2,dum3,dum4
c       write(6,6329)dum1,dum2,dum3,dum4
c6329   format(1h ,4f9.2)
        if(dum1 .eq. wave)then
c        tauratio is ratio of tau in specified band to B-band tau
         tauratio=dum2
         albedo=dum3
         g=dum4
         iflag=iflag+1
        endif 
 2563  continue
       if(iflag.ne.1)then
       write(6,8266)
 8266  format(1h ,'URAD: correct albedo,g not found - aborting:')
       write(6,2564)ifilter,albedo,g
 2564  format(1h ,'ifilter,albedo,g = ',i5,2f9.2)
       goto 999
      endif
c
c     opacities
c
c
c       factor1=0.279
c       factor2=0.721

c      B-band tau are input and converted to tau at wavelength wave
       totaltau_in=totaltau
       totaltau=totaltau_in*tauratio
c       tau0=factor1*totaltau
c       tau1=totaltau*factor2
c       tau0_in=factor1*totaltau_in
c       tau1_in=totaltau_in*factor2

      tau0_in = totaltau1
      tau1_in = totaltau2
      tau0 = tau0_in * tauratio
      tau1 = tau1_in * tauratio
ccc JJT 3 Oct 17 ccc
      tau3_in = totaltau3
      tau4_in = totaltau4
      tau3 = tau3_in * tauratio
      tau4 = tau4_in * tauratio

      tau5_in = totaltau5
      tau6_in = totaltau6
      tau5 = tau5_in * tauratio
      tau6 = tau6_in * tauratio
ccc JJT 3 Oct 17 ccc
c
c     benchmark central disk brightnesses in mag per square arcsec
c     and geometry-dependent interpolation scheme
c
      refeta=20.0
      eta0=refeta
      etab=refeta
      eta1=refeta
ccc JJT 3 Oct 17 ccc
      eta3=refeta
      eta4=refeta
      eta5=refeta
      eta6=refeta
ccc JJT 3 Aug 18 ccc
      eta7=refeta
ccc JJT 3 Oct 17 ccc
c      reflambda=4000.
      reflambda=3600.
      ifswitch=1
      if(wave.lt.reflambda.and.morph.eq.'d')then
       write(6,8289)
8289   format(1h ,'URAD: wrong filter morph combination - aborting:')
       write(6,2575)filter,morph
2575   format(1h ,'filter, morph = ',2a10)
       goto 999
      endif
      if(wave.ge.reflambda.and.morph.eq.'d')then
        eta1=50.0
        eta0=refeta
        etab=50.0
ccc JJT 3 Oct 17 ccc
      eta3=50.0
      eta4=50.0
      eta5=50.0
      eta6=50.0
      eta7=50.0 ! JJT 3 Aug 18
ccc JJT 3 Oct 17 ccc
      endif
      if(wave.lt.reflambda.and.morph.eq.'td')then
        eta1=refeta
        eta0=50.0
        etab=50.0
ccc JJT 3 Oct 17 ccc
      eta3=50.0
      eta4=50.0
      eta5=50.0
      eta6=50.0
      eta7=50.0 ! JJT 3 Aug 18
ccc JJT 3 Oct 17 ccc
        ifswitch=4
      endif
      if(wave.ge.reflambda.and.morph.eq.'td')then
        eta1=refeta
        eta0=50.0
        etab=50.0
ccc JJT 3 Oct 17 ccc
      eta3=50.0
      eta4=50.0
      eta5=50.0
      eta6=50.0
      eta7=50.0 ! JJT 3 Aug 18
ccc JJT 3 Oct 17 ccc
        ifswitch=4
      endif
      if(wave.ge.reflambda.and.morph.eq.'b')then
        eta1=50.0
        eta0=50.0
        etab=refeta
ccc JJT 3 Oct 17 ccc
      eta3=50.0
      eta4=50.0
      eta5=50.0
      eta6=50.0
      eta7=50.0 ! JJT 3 Aug 18
ccc JJT 3 Oct 17 ccc
        ifswitch=3
      endif
      if(wave.lt.reflambda.and.morph.eq.'b')then
       write(6,8289)
       write(6,2575)filter,morph
       goto 999
      endif

ccc JJT 3 Oct 17 ccc
ccc inner component
      if(wave.lt.reflambda.and.morph.eq.'di')then
       write(6,8289)
c      format(1h ,'URAD: wrong filter morph combination - aborting:')
       write(6,2575)filter,morph
c      format(1h ,'filter, morph = ',2a10)
       goto 999
      endif
      if(wave.ge.reflambda.and.morph.eq.'di')then
        eta1=50.0
        eta0=50.0
        etab=50.0
        eta3=refeta
        eta4=50.0
        eta5=50.0
        eta6=50.0
        eta7=50.0 ! JJT 3 Aug 18
      endif
      if(wave.lt.reflambda.and.morph.eq.'tdi')then
        eta1=50.0
        eta0=50.0
        etab=50.0
        eta3=50.0
        eta4=refeta
        eta5=50.0
        eta6=50.0
        eta7=50.0 ! JJT 3 Aug 18
        ifswitch=4
      endif
      if(wave.ge.reflambda.and.morph.eq.'tdi')then
        eta1=50.0
        eta0=50.0
        etab=50.0
        eta3=50.0
        eta4=refeta
        eta5=50.0
        eta6=50.0
        eta7=50.0 ! JJT 3 Aug 18
        ifswitch=4
      endif

ccc outer component
      if(wave.lt.reflambda.and.morph.eq.'do')then
       write(6,8289)
c      format(1h ,'URAD: wrong filter morph combination - aborting:')
       write(6,2575)filter,morph
c      format(1h ,'filter, morph = ',2a10)
       goto 999
      endif
      if(wave.ge.reflambda.and.morph.eq.'do')then
        eta1=50.0
        eta0=50.0
        etab=50.0
        eta3=50.0
        eta4=50.0
        eta5=refeta
        eta6=50.0
        eta7=50.0 ! JJT 3 Aug 18
      endif
      if(wave.lt.reflambda.and.morph.eq.'tdo')then
        eta1=50.0
        eta0=50.0
        etab=50.0
        eta3=50.0
        eta4=50.0
        eta5=50.0
        eta6=refeta
        eta7=50.0 ! JJT 3 Aug 18
        ifswitch=4
         write(6,*)'eta6'
         write(6,*) eta6
      endif
      if(wave.ge.reflambda.and.morph.eq.'tdo')then
        eta1=50.0
        eta0=50.0
        etab=50.0
        eta3=50.0
        eta4=50.0
        eta5=50.0
        eta6=refeta
        eta7=50.0 ! JJT 3 Aug 18
        ifswitch=4
         write(6,*)'eta6'
         write(6,*) eta6
      endif


ccc JJT 3 Oct 17 ccc
ccc JJT 3 Aug 18 ccc
      if(wave.le.reflambda.and.morph.eq.'ntd')then
        eta1=50.0
        eta0=50.0
        etab=50.0
        eta3=50.0
        eta4=50.0
        eta5=50.0
        eta6=50.0
        eta7=refeta ! JJT 3 Aug 18
        ifswitch=4
         write(6,*)'eta7'
         write(6,*) eta7
      endif
      if(wave.ge.reflambda.and.morph.eq.'ntd')then
        eta1=50.0
        eta0=50.0
        etab=50.0
        eta3=50.0
        eta4=50.0
        eta5=50.0
        eta6=50.0
        eta7=refeta ! JJT 3 Aug 18
        ifswitch=4
         write(6,*)'eta7'
         write(6,*) eta7
      endif

c

c      fixing the rest of the geometrical parameters


c
      hs=h_bdisk
c
      if(wave.gt.reflambda)then
       if(filter.eq.'v')hs=h_vdisk
       if(filter.eq.'i')hs=h_idisk
       if(filter.eq.'j')hs=h_jdisk
       if(filter.eq.'k')hs=h_kdisk
       if(filter.eq.'ir36')hs=h_ir36disk
       if(filter.eq.'ir45')hs=h_ir45disk
       if(filter.eq.'ir58')hs=h_ir58disk
      endif

ccc JJT 3 Oct 17 ccc
c Setting correct hs for extra old disks
      hs3=h_bdisk3
c
      if(wave.gt.reflambda)then
       if(filter.eq.'v')hs3=h_vdisk3
       if(filter.eq.'i')hs3=h_idisk3
       if(filter.eq.'j')hs3=h_jdisk3
       if(filter.eq.'k')hs3=h_kdisk3
       if(filter.eq.'ir36')hs3=h_ir36disk3
       if(filter.eq.'ir45')hs3=h_ir45disk3
       if(filter.eq.'ir58')hs3=h_ir58disk3
      endif

       hs5=h_bdisk5
c
      if(wave.gt.reflambda)then
       if(filter.eq.'v')hs5=h_vdisk5
       if(filter.eq.'i')hs5=h_idisk5
       if(filter.eq.'j')hs5=h_jdisk5
       if(filter.eq.'k')hs5=h_kdisk5
       if(filter.eq.'ir36')hs5=h_ir36disk5
       if(filter.eq.'ir45')hs5=h_ir45disk5
       if(filter.eq.'ir58')hs5=h_ir58disk5
      endif
ccc JJT 3 Oct 17 ccc
c
c     other parameters
c
c       rtruncate=24.0
c       sharp=0.002
c       rtruncate1 = 24.0
c       sharp1 = 0.002
c       bfactor=1.d-10
c      bfactor to unity for energy density calculations
       bfactor=1.0d0 
       bin=1.d0
c      integration scheme of Nick (near edge-on) and RT (otherwise)
c       ifswitch=1
c
c former config.in
c izero = 0 means pure absorption (no scattering terms)
c       = 1 to include 1st order scattered light (exact calculation)
c       = 2 or higher to include higher order scattered light (approximation)
c idust = 1 (2 means dust.plo written column by column;
c            1 means dust.plo written point by point;
c            0 means no dust.plo file output)
c itau = 0 means No optical depth calculations are taking place
c istars = 0 means No Pure starlight calculations are taking place
c idisk1, idisk2 (old and young disks)
c         =  0 for no emissivity and opacity 
c         =  1 for double exponential emissivity and opacity
c         =  2 spiral arm  - only for the second dust disk and young
c                            stellar disk
c         = 3 exponential in radial and sech2 in vertical
c ibulge= 0 means that NO bulge component is used
c       = 1 means Modified Hubble profile is used
c       = 2 means R^1/4 law is used
c       = 3 means Sersic
c iferrers = 0 means that NO ferrers bar component is used
c          = 1 means that ferrers bar will be calculated
c iellexp  = 0 means that NO elliptical component is used
c          = 1 means that elliptical disk will be calculated
c imask = 0 means no mask file is output
c iall = 1 means that the whole image of the galaxy is created
c      = else half part of the galaxy (folded case?) is created
 
c
c     izero set in input file
c      izero=0
      idust=1
      itau=0
      istars=0
c      idisk1=3
c      idisk2=3
      ibulge=3
c      nsersic=4
      iferrers=0
      iellexp=0
      imask=0
      iall=0
      shiftx=0.0d0
      shifty=0.0d0

c
c set up integration step for dust.f
c      durect light
      fs0=0.15 
c     scattered light
      fs1=0.20
      fsdum=0.5
c     to be used for radiation fields
      fsurad=0.2


c
       theta=inclination
c

c
c     ferrers bar
       etaf=21.0 
       af=6.5
       bf=3.2
       cf=0.4
       thetaf=-43.0
c     elliptical disk
      eta2=28.6
      hsa2=3.8
      hsb2=1.5
      zs2=0.4
      thetae=-43.0 
      hda2=3.8
      hdb2=1.5
      zd2=0.23
      tau2=0.5
c 
c pixel sampling (for the ouput image and also some rays)
c 0.5 kpc resolution is 10.8560 arcsec for N891 at 9.5 Mpc
c Need 48 pixels to reach 24 Kpc which is the end of the dust disk.
c     y is horizontal (= in plane for edge on seen horizontal)
c     x is vertical. These are Nick's conventions.
c       nx=96*5
c       ny=48*5
c       nx=192
c       ny=96

      if(morph .eq. 'd' .or. morph .eq. 'td') then
c        nx=4400
c        ny=2200
c        nx=1500
c        ny=750
        nx=nx_m
        ny=ny_m
c          nx=300
c          ny=150
c         nx=1400
c         ny=700
        mask_opt=mask_opt_m
        mstep1=mstep1_m
        mstep2=mstep2_m
        mstep3=mstep3_m
        mlength1=mlength1_m
        mlength2=mlength2_m
        elm=elm_m
      endif

ccc JJT 3 Oct 17 ccc
ccc fits file dimensions ccc
      if(morph .eq. 'di' .or. morph .eq. 'tdi') then
c        nx=4400
c        ny=2200
        nx=nx_i !will give a 2kpc radius
        ny=ny_i !at d=859kpc
c        nx=3000
c        ny=1500
c          nx=300
c          ny=150
c         nx=1400
c         ny=700
        mask_opt=mask_opt_i
        mstep1=mstep1_i
        mstep2=mstep2_i
        mstep3=mstep3_i
        mlength1=mlength1_i
        mlength2=mlength2_i
        elm=elm_i
      endif

      if(morph .eq. 'do' .or. morph .eq. 'tdo') then
c        nx=4400
c        ny=2200
c        nx=1500
c        ny=750
        nx=nx_o
        ny=ny_o
c          nx=300
c          ny=150
c         nx=1400
c         ny=700
        mask_opt=mask_opt_o
        mstep1=mstep1_o
        mstep2=mstep2_o
        mstep3=mstep3_o
        mlength1=mlength1_o
        mlength2=mlength2_o
        elm=elm_o
      endif
ccc JJT 3 Aug 18ccc 
      if(morph .eq. 'ntd') then
c        nx=4400
c        ny=2200
c        nx=1500
c        ny=750
        nx=nx_n 
        ny=ny_n 
c          nx=300
c          ny=150
c         nx=1400
c         ny=700
        mask_opt=mask_opt_n
        mstep1=mstep1_n
        mstep2=mstep2_n
        mstep3=mstep3_n
        mlength1=mlength1_n
        mlength2=mlength2_n
        elm=elm_n
      endif

ccc JJT 3 Oct 17 ccc
      if(morph .eq. 'b') then
c        nx=4400
c        ny=2200
c         nx=1500
c         ny=750
        nx=nx_b
        ny=ny_b
        mask_opt=mask_opt_b
        mstep1=mstep1_b
        mstep2=mstep2_b
        mstep3=mstep3_b
        mlength1=mlength1_b
        mlength2=mlength2_b
        elm=elm_b
      endif

        print*,nx
        print*,ny

c pix size in arcsec
c       pixsize=1.4
c        pixsize=4.8
        pixsize=pix_size
c      pixsize=5.1455
c      pixsize=0.24
c       pixsize=10.8560/5.
c       pixsize=10.8560/2.0
c
c      do j=1,nx*ny
c         mask(j)=0
c      enddo
      ii=0
c JJT pixels between next mask strip | could set this for different radii
c large towards inner larger outer
      icont=1
c       icont=10
c      icont=5
      do 105 j=1,ny
      do 106 i=1,nx
      ii=ii+1

cif mask opt=1 then apply mask as described below
       if (mask_opt .eq. 1) then
         if (ii .le. mlength1*nx) then
            icont=mstep1      
         endif
         if (ii .gt. mlength1*nx .and. ii .le. mlength2*nx) then
            icont=mstep2
         endif
         if (ii .gt. mlength2*nx .and. ii .le. ny*nx) then
            icont=mstep3
         endif
       endif   

      mask(ii)=1
      jtest=(j-1)/icont
c     calculate in vertical stripes separated by icont pixels
      if((jtest*icont) .eq. j-1)mask(ii)=0
	   conver=3.14159d0/180.0d0 !deg to rad
c     add eliptical mask
      ellipse=sqrt(((i-nx/2)/(dcos(theta*conver)))**2+j**2)
      if(ellipse .lt. elm) mask(ii)=1   
c      if(i .eq. 5 .and. j .eq. 2)then 
c         mask(ii)=1
c         write(6,107)i,j,ii
c 107     format('direction masked with i,j,ii = ',3i6)
c      endif
c      if(j .ge. 84 .and. j .le. 88)then 
c         mask(ii)=0
c         write(6,107)i,j,ii
c      endif
 106    continue
 105    continue
c
c      call printout
      write(6,982)
982   format(1h ,'URAD input parameters:')
      write(6,983)izero
983   format(1h ,'izero = ',i6)
      write(6,984)filter
984   format(1h ,'filter = ',a6)
      write(6,985)dustmodel
985   format(1h ,'dustmodel = ',a10)
c
c      write(6,8267)ifilter,albedo,g
c
c

c     create string snsersic corresponding to nsersic
      insersic = int(nsersic)
      write(sdummy,*) insersic
      snsersic = adjustl(sdummy)
      lensnsersic=len_trim(trim(snsersic))
      snsersic=snsersic(1:lensnsersic)
c      istart=0
c      do 1277 i=1,lensnsersic
c       if(snsersic(i:i).ne.' '.and.istart.eq.0)istart=i
c 1277 continue
c      snsersic=snsersic(istart:lensnsersic)
c      print*, 'test2=',snsersic      
c      lensnsersic=len_trim(trim(snsersic))
c      snsersic=snsersic(1:lensnsersic)
c      print*, 'test3=',snsersic
c      STOP 'NSERSIC TEST'

c      read to common block the steering
c      parameters from file urad.in in directory path_urad 
c
      filename_uradin=
     >path1(1:lenpath1)//path2(1:lenpath2)//
     >'indata/urad.in'
      filename_uradin=trim(filename_uradin)
c
      call readuradpar(filename_uradin)
c
c create energy density file with correct name
c
      if(izero.eq.0)sabssca='abs'
      if(izero.ne.0)sabssca='sca'
      lensabssca=len_trim(trim(sabssca))
      sabssca=sabssca(1:lensabssca)

c      path=path1(1:lenpath1)//
c     >               path2(1:lenpath2)//
c     >              'urad_'//filter(1:lenfilter)//
c     >              '/outdata/'   
c      filename_urad=
c     >              'urad_m'//morph(1:lenmorph)//
c     >              '_'//filter(1:lenfilter)//
c     >              '_'//dustmodel(1:lendustmodel)//
c     >              '_t'//stau(1:lenstau)//
c     >              '_'//sabssca(1:lensabssca)//
c     >              '.dat'
c      write(6,9826)path, filename_urad
c9826  format(1h ,'URAD: writing ',a, a) 
c
c      open(unit=12,file=path//filename_urad)
ccc JJT 3 Oct 17 ccc
      if (morph .eq. 'd' .or. morph .eq. 'di' .or. morph .eq. 'do') then
ccc JJT 3 Oct 17 ccc
       filename_uradd=path1(1:lenpath1)//
     >               path2(1:lenpath2)//
     >              '/out/'//
     >              'u_m'//morph(1:lenmorph)//
     >              '_'//filter(1:lenfilter)//
     >              '_'//dustmodel(1:lendustmodel)//
     >              '_t'//stau(1:lenstau)//
     >              '_hd'//shd(1:lenshd)//
     >              '_zd'//szd(1:lenszd)//
     >              '_hd1_'//shd1(1:lenshd1)//
     >              '_zd1_'//szd1(1:lenszd1)//
     >              '_hs'//shs(1:lenshs)//
     >              '_zs'//szs(1:lenszs)//
     >              '_'//sabssca(1:lensabssca)//
     >              '.dat'
       write(6,9826)filename_uradd
       open(unit=12,file=filename_uradd)
c
c     note that tau1 in dust SED generation corresponds to tau0 in radtrans
c     there are also iferences in other parameters
c     note too that whereas for dust.f lengths are input in kpc,
c     lengths are input in pc to the sed generation program
        write(12,701)'tau1'
        write(12,702)tau0_in
        write(12,701)'hd [pc]'
        write(12,702)hd*1000.
        write(12,701)'zd [pc]'
        write(12,702)zd*1000.
        write(12,701)'hdin [pc]'
        write(12,702) hdin * 1000.
        write(12,701) 'zdin [pc]'
        write(12,702) zdin * 1000.
        write(12,7016) 'hdsolar [pc]'
        write(12,702) hdsolar * 1000.
        write(12,7016) 'zdsolar [pc]'
        write(12,702) zdsolar * 1000.
        write(12,701)'tau2'
        write(12,702)tau1_in
        write(12,701)'hd1 [pc]'
        write(12,702)hd1*1000.
        write(12,701)'zd1 [pc]'
        write(12,702)zd1*1000.
        write(12,701)'hd1in [pc]'
        write(12,702) hd1in * 1000.
        write(12,701) 'zd1in [pc]'
        write(12,702) zd1in * 1000.
        write(12,7016) 'hd1solar [pc]'
        write(12,702) hd1solar * 1000.
        write(12,7016) 'zd1solar [pc]'
        write(12,702) zd1solar * 1000.
        write(12,701) 'hs [pc]'
        write(12,702) hs * 1000.
        write(12,701) 'zs [pc]'
        write(12,702) zs * 1000.
        write(12,701)'hsin [pc]'
        write(12,702) hsin * 1000.
        write(12,701) 'zsin [pc]'
        write(12,702) zsin * 1000.
        write(12,7016) 'hssolar [pc]'
        write(12,702) hssolar * 1000.
        write(12,7016) 'zssolar [pc]'
        write(12,702) zssolar * 1000.
        write(12,7016) 'rtruncate [pc]'
        write(12,702) rtruncate * 1000.
        write(12,701) 'xis0'
        write(12,702) xis0
        write(12,701) 'xid0'
        write(12,702) xid0
        write(12,701) 'xid1'
        write(12,702) xid1      
        write(12,701) 'idisk1'
        write(12,703) idisk1
        write(12,701) 'idisk2'
        write(12,703) idisk2
      endif
ccc JJT 3 Oct 17 ccc
              if (morph.eq.'td'.or.morph.eq.'tdi'.or.morph.eq.'tdo'
     >               .or.morph.eq.'ntd') then  !update 3 Aug 18 JJT
ccc JJT 3 Oct 17 ccc
       filename_uradtd=path1(1:lenpath1)//
     >               path2(1:lenpath2)//
     >              '/out/'//
     >              'u_m'//morph(1:lenmorph)//
     >              '_'//filter(1:lenfilter)//
     >              '_'//dustmodel(1:lendustmodel)//
     >              '_t'//stau(1:lenstau)//
     >              '_hd'//shd(1:lenshd)//
     >              '_zd'//szd(1:lenszd)//
     >              '_hd1_'//shd1(1:lenshd1)//
     >              '_zd1_'//szd1(1:lenszd1)//
     >              '_hs1_'//shs1(1:lenshs1)//
     >              '_zs1_'//szs1(1:lenszs1)//
     >              '_'//sabssca(1:lensabssca)//
     >              '.dat'
       write(6,9826)filename_uradtd
       open(unit=12,file=filename_uradtd)

c     note that tau1 in dust SED generation corresponds to tau0 in radtrans
c     there are also iferences in other parameters
c     note too that whereas for dust.f lengths are input in kpc,
c     lengths are input in pc to the sed generation program
       write(12,701)'tau1'
       write(12,702)tau0_in
       write(12,701)'hd [pc]'
       write(12,702)hd*1000.
       write(12,701)'zd [pc]'
       write(12,702)zd*1000.
       write(12,701)'hdin [pc]'
       write(12,702) hdin * 1000.
       write(12,701) 'zdin [pc]'
       write(12,702) zdin * 1000.
       write(12,7016) 'hdsolar [pc]'
       write(12,702) hdsolar * 1000.
       write(12,7016) 'zdsolar [pc]'
       write(12,702) zdsolar * 1000.
       write(12,701)'tau2'
       write(12,702)tau1_in
       write(12,701)'hd1 [pc]'
       write(12,702)hd1*1000.
       write(12,701)'zd1 [pc]'
       write(12,702)zd1*1000.
       write(12,701)'hd1in [pc]'
       write(12,702) hd1in * 1000.
       write(12,701) 'zd1in [pc]'
       write(12,702) zd1in * 1000.
       write(12,7016) 'hd1solar [pc]'
       write(12,702) hd1solar * 1000.
       write(12,7016) 'zd1solar [pc]'
       write(12,702) zd1solar * 1000.
       write(12,701) 'hs1 [pc]'
       write(12,702) hs1 * 1000.
       write(12,701) 'zs1 [pc]'
       write(12,702) zs1 * 1000.
       write(12,701)'hs1in [pc]'
       write(12,702) hs1in * 1000.
       write(12,701) 'zs1in [pc]'
       write(12,702) zs1in * 1000.
       write(12,7016) 'hs1solar [pc]'
       write(12,702) hs1solar * 1000.
       write(12,7016) 'zs1solar [pc]'
       write(12,702) zs1solar * 1000.
       write(12,7016) 'rtruncate1 [pc]'
       write(12,702) rtruncate1 * 1000.
       write(12,701) 'xis1'
       write(12,702) xis1
       write(12,701) 'xid0'
       write(12,702) xid0
       write(12,701) 'xid1'
       write(12,702) xid1      
       write(12,701) 'idisk1'
       write(12,703) idisk1
       write(12,701) 'idisk2'
       write(12,703) idisk2
      endif
      if (morph .eq. 'b') then
       filename_uradb=path1(1:lenpath1)//
     >               path2(1:lenpath2)//
     >              '/out/'//
     >              'u_m'//morph(1:lenmorph)//
     >              '_'//filter(1:lenfilter)//
     >              '_'//dustmodel(1:lendustmodel)//
     >              '_t'//stau(1:lenstau)//
     >              '_hd'//shd(1:lenshd)//
     >              '_zd'//szd(1:lenszd)//
     >              '_hd1_'//shd1(1:lenshd1)//
     >              '_zd1_'//szd1(1:lenszd1)//
     >              '_reff'//sreff(1:lensreff)//
     >              '_ell'//sellipt(1:lensellipt)//
     >              '_n'//snsersic(1:lensnsersic)//
     >              '_'//sabssca(1:lensabssca)//
     >              '.dat'
       write(6,9826)filename_uradb
       open(unit=12,file=filename_uradb)
c
c     note that tau1 in dust SED generation corresponds to tau0 in radtrans
c     there are also iferences in other parameters
c     note too that whereas for dust.f lengths are input in kpc,
c     lengths are input in pc to the sed generation program
       write(12,701)'tau1'
       write(12,702)tau0_in
       write(12,701)'hd [pc]'
       write(12,702)hd*1000.
       write(12,701)'zd [pc]'
       write(12,702)zd*1000.
       write(12,701)'hdin [pc]'
       write(12,702) hdin * 1000.
       write(12,701) 'zdin [pc]'
       write(12,702) zdin * 1000.
       write(12,7016) 'hdsolar [pc]'
       write(12,702) hdsolar * 1000.
       write(12,7016) 'zdsolar [pc]'
       write(12,702) zdsolar * 1000.
       write(12,701)'tau2'
       write(12,702)tau1_in
       write(12,701)'hd1 [pc]'
       write(12,702)hd1*1000.
       write(12,701)'zd1 [pc]'
       write(12,702)zd1*1000.
       write(12,701)'hd1in [pc]'
       write(12,702) hd1in * 1000.
       write(12,701) 'zd1in [pc]'
       write(12,702) zd1in * 1000.
       write(12,7016) 'hd1solar [pc]'
       write(12,702) hd1solar * 1000.
       write(12,7016) 'zd1solar [pc]'
       write(12,702) zd1solar * 1000.
       write(12,701) 'reff [pc]'
       write(12,702) reff * 1000.
       write(12,701) 'ellipt'
       write(12,702) ellipt
       write(12,701) 'nsersic'
       write(12,703) nsersic
       write(12,701) 'xid0'
       write(12,702) xid0
       write(12,701) 'xid1'
       write(12,702) xid1      
       write(12,701) 'idisk1'
       write(12,703) idisk1
       write(12,701) 'idisk2'
       write(12,703) idisk2
      endif
  
701   format(a10)
7011  format(a14)
7016  format(a15)
70111 format(a64)
702   format(f10.2)
722   format(f10.4)
9826  format(1h ,'URAD: writing ',a) 


c
      call NUdust(mask,dzmodel,morph(1:lenmorph),filter(1:lenfilter),
     >   diag_plots)
c
      write(12,7011)'nposition in r'
      write(12,703)nru_used
703   format(i10)
      write(12,7011)'nposition in z'
      write(12,703)nzu_used
      write(12,701)''
      write(12,70111)'r(pc)       z(pc)       urad(erg/pc3/A)  
     >   eabs(erg/s/pc3/A)'
c
!      open(unit=20,file='dust.urad')
!      read(20,*)nradtot,ndum1,ndum2,ndum3,xdum
      nradtot = nru * nzu
!      print*, shape(yout)
!      STOP
      do 750 i=1,nradtot
!       read(20,275)yout,zout,uout,faout
c       write(6,2688)yout,zout,uout,faout
c2688   format(1h ,4f10.2)
275    format(1h ,4e12.4)
       write(12,275)yout(i),zout(i),uout(i),faout(i)
c       write(6,2688)yout,zout,uout,faout
750   continue

c create image file with correct name
c
ccc JJT 4 Oct 17 ccc
c commented out what appears to be unused old code
c for aborting if there is a theta value ne 0
c      if (theta.ne.0.)then
c       write(6,6666)
c 6666  format(1h,'URAD: inclination incorrectly specified - aborting')
c       write(6,6667) theta
c 6667  format(f10.2)
c      endif
ccc JJT 4 Oct 17 ccc
c      if (morph .eq. 'd' .or. morph .eq. 'td') then
ccc JJT 3 Oct 17 ccc
c di - inner disk, do- outer disk, d - original disk, n - nuclear disk
      if (morph .eq. 'd' .or. morph .eq. 'di' .or. morph .eq. 'do' .or.
     >    morph.eq.'td'.or.morph.eq.'tdi'.or.morph.eq.'tdo' .or.
     >    morph .eq. 'ntd')then !JJT 3 Aug 18
ccc JJT 3 Oct 17 ccc
       filename_image=path1(1:lenpath1)//
     >               path2(1:lenpath2)//
     >               '/out/'//
     >               'map_m'//morph(1:lenmorph)//
     >               '_'//filter(1:lenfilter)//
     >               '_'//dustmodel(1:lendustmodel)//
     >               '_i'//sinclination(1:lensinclination)//
     >               '_t'//stau(1:lenstau)//
     >               '_'//sabssca(1:lensabssca)//
     >               '.dat'
      endif

      if (morph .eq. 'b') then
       filename_image=path1(1:lenpath1)//
     >               path2(1:lenpath2)//
     >               '/out/'//
     >               'map_m'//morph(1:lenmorph)//
     >               '_'//filter(1:lenfilter)//
     >               '_'//dustmodel(1:lendustmodel)//
     >               '_i'//sinclination(1:lensinclination)//
     >               '_t'//stau(1:lenstau)//
     >               '_n'//snsersic(1:lensnsersic)//
     >               '_'//sabssca(1:lensabssca)//
     >               '.dat'
      endif

      write(6,6668)filename_image

 6668 format(1h ,'URAD: writing ',a) 

      open(unit=14,file=filename_image)
c      open(unit=22,file='dust.plo')
      nimage = nx*ny
      do 777 i=1,nimage
c         read(22,*) bout
         write(14,*) dzmodel(i)
 777  continue

c
c      call printout
c
c      close(unit=1)
      close(unit=2)
      close(unit=12)
      close(unit=20)
      close(unit=14)
!      close(unit=22)
c
 999  continue
c       write(6,8287)ifilter,albedo,g
c 8287  format(1h ,'URAD: ifilter,albedo,g = ',i6,2f8.2)
      stop
      end

        subroutine printout

        Implicit real*8 (a-h, o-z)
c        PARAMETER (NSIZE=9680000)
c        PARAMETER (NSIZE=1125000)
        PARAMETER (NSIZE=45000000)
c        PARAMETER (NSIZE=980000)
        parameter (nsizeu=50000)
        DOUBLE PRECISION DZMODEL(NSIZE)
        double precision u(nsizeu),fa(nsizeu)
        double precision fv(nsizeu,256)
        double precision uct(256),ust(256),ucp(256),usp(256)
        integer mask(nsize)
        common/incli/theta,tau0
     	  common/param/zs,zd,hs,hd,d,rhomax,hdin,hsin,
     >   zsin,zssolar,hssolar,zdin,zdsolar,hdsolar,
     >   xis0,xid0,hstin,hdtin
        common/f/fs0,fs1,fsdum,fsurad,dsb1,dsb2
        common/atte/a0,eta0,rtruncate,sharp,athresh,rtrun,sha,
     >   rtruncated,sharpd,rtrund,shad
        common/bulge/etab,reff,ellipt2,ellipt,nsersic,bsersic
        common/ferrers/etaf,af,bf,cf,thetaf,ctf,stf,iferrers
        common/ellexp/eta2,zs2,hsa2,hsb2,thetae,cte,ste,
     >  hda2,hdb2,zd2,tau2,a2,iellexp
        common/bcons/etabmax,rin,ibulge
        common/g/g,g2,cosalp
        common/integ/k,m,n,km
        common/newdisk/eta1,zs1,hs1,zd1,hd1,a1,tau1,idisk1,idisk2,
     >  rtruncate1,sharp1,rtrun1,sha1,hd1in,zd1in,hs1in,zs1in,
     >  hd1solar,zd1solar,hs1solar,zs1solar,
     >  xis1,xid1,hs1tin,hd1tin,
     >  eta1arm,a1arm,narm,rup,fillarm,
     >   rtruncated1,sharpd1,rtrund1,shad1
        common/translate/pix2kpc
        common/inout/nx,ny,pixsize,nimage
        common/ur/u,fa,fv,ures,rsizeu,zsizeu,albedo,iscat
        common/ur2/ru,zu,nru,nzu,iurad,nru_used,nzu_used
        common/urang/icang,nang,uct,ust,ucp,usp
        common/photom/away,bfactor,bin,ifilter
        common/config/izero,idust,itau,istars,imask,iall,
     >            shiftx,shifty,ifswitch
	     common/eindisk3/eta3,zs3,hs3,zd3,hd3,a3,tau3,idisk3,
     >  rtruncate3,sharp3,rtrun3,sha3,hd3in,zd3in,hs3in,zs3in,
     >  hd3solar,zd3solar,hs3solar,zs3solar,
     >  xis3,xid3,hs3tin,hd3tin,
     >  eta3arm,a3arm,narm3,rup3,fillarm3,
     >   rtruncated3,sharpd3,rtrund3,shad3
c tdisk
	     common/eindisk4/eta4,zs4,hs4,zd4,hd4,a4,tau4,idisk4,
     >  rtruncate4,sharp4,rtrun4,sha4,hd4in,zd4in,hs4in,zs4in,
     >  hd4solar,zd4solar,hs4solar,zs4solar,
     >  xis4,xid4,hs4tin,hd4tin,
     >  eta4arm,a4arm,narm4,rup4,fillarm4,
     >   rtruncated4,sharpd4,rtrund4,shad4

c extra disks for outer regions
c disk
	     common/eoutdisk5/eta5,zs5,hs5,zd5,hd5,a5,tau5,idisk5,
     >  rtruncate5,sharp5,rtrun5,sha5,hd5in,zd5in,hs5in,zs5in,
     >  hd5solar,zd5solar,hs5solar,zs5solar,
     >  xis5,xid5,hs5tin,hd5tin,
     >  eta5arm,a5arm,narm5,rup5,fillarm5,
     >   rtruncated5,sharpd5,rtrund5,shad5
c tdisk
	     common/eoutdisk6/eta6,zs6,hs6,zd6,hd6,a6,tau6,idisk6,
     >  rtruncate6,sharp6,rtrun6,sha6,hd6in,zd6in,hs6in,zs6in,
     >  hd6solar,zd6solar,hs6solar,zs6solar,
     >  xis6,xid6,hs6tin,hd6tin,
     >  eta6arm,a6arm,narm6,rup6,fillarm6,
     >   rtruncated6,sharpd6,rtrund6,shad6
ccc JJT 3 Aug 18 - Common blocks ccc
c nulear comp ( young stellar disc
	common/ntdisk7/eta7,zs7,hs7,a7,idisk7,
     >  rtruncate7,sharp7,rtrun7,sha7,hs7in,zs7in,
     >  hs7solar,zs7solar,
     >  xis7,hs7tin,
     >  eta7arm,a7arm,narm7,rup7,fillarm7
c
c
        write(6,890)theta
 890   format('theta =',f10.4)
        write(6,891)eta0
 891    format('eta0 =',f14.4)
        write(6,892)tau0
 892    format('tau =',f10.4)
        write(6,893)zs
 893    format('zs =',f10.4)
        write(6,894)zd
 894    format('zd =',f10.4)
        write(6,895)hs
 895    format('hs =',f10.4)
        write(6,896)hd
 896    format('hd =',f10.4)
        write(6,2966)hsin
2966    format('hsin=',f10.4)
        write(6,8966)hdin
 8966   format('hdin=',f10.4)
        write(6,8967) zd1t
 8967   format('zd1t=',f10.4)        
        write(6,897)ellipt
 897    format('ellipt =',f10.4)
        write(6,898)reff
 898    format('reff =',f10.4)
        write(6,899)etab
 899    format('etab =',f14.4)
        write(6,8999)nsersic
 8999   format('nsersic =',i10)
        write(6,794)eta1
 794    format('eta1 =',f10.4)
        write(6,795)zs1
 795    format('zs1 =',f10.4)
        write(6,796)hs1
 796    format('hs1 =',f10.4)
        write(6,798)zd1
 798    format('zd1 =',f10.4)
        write(6,799)hd1
 799    format('hd1 =',f10.4)
        write(6,7999)hd1in
 7999   format('hd1in=',f10.4)
        write(6,7996)zd1t
 7996   format('zd1t=',f10.4)
        write(6,797)tau1
 797    format('tau1 =',f10.4) 
        write(6,599)etaf
 599    format('etaf =',f10.4)
        write(6,598)af
 598    format('af =',f10.4)
        write(6,597)bf
 597    format('bf =',f10.4)
        write(6,596)cf
 596    format('cf =',f10.4)
        write(6,196)thetaf
 196    format('thetaf =',f10.4)
        write(6,5599)eta2
 5599   format('eta2 =',f10.4)
        write(6,5598)hsa2
 5598   format('hsa2 =',f10.4)
        write(6,5597)hsb2
 5597   format('hsb2 =',f10.4)
        write(6,5596)zs2
 5596   format('zs2 =',f10.4)
        write(6,197)thetae
 197    format('thetae =',f10.4)
        write(6,6599)tau2
 6599   format('tau2 =',f10.4)
        write(6,6598)hda2
 6598   format('hda2 =',f10.4)
        write(6,6597)hdb2
 6597   format('hdb2 =',f10.4)
        write(6,6596)zd2
 6596   format('zd2  =',f10.4)
        write(6,699)rtrun
 699    format('rtrun =',f10.4)
        write(6,199)sha
 199    format('sha =',f10.4)
        write(6,966)rhomax
 966    format('rhomax =',f10.4)
        write(6,967)d
 967    format('d =',f10.4)
        write(6,968)pixsize
 968    format('pixsize =',f10.4)
        write(6,969)ny,nx
 969    format('ny nz =',2i6)
        write(6,975)ifilter,albedo,g
 975    format('ifilter, albedo,g =',i6,2f10.4)
        write(6,976)idisk1,idisk2,ibulge,iferrers,iellexp,izero
 976    format('idisk1,idisk2,ibulge,iferrers,iellexp,izero =',6i4)
        write(6,977)iurad,ures,icang,rsizeu,zsizeu
 977    format('iurad,ures,icang,rsizeu,zsizeu =',i6,f10.4,i6,2f10.4)
c
      return
      end
c
       subroutine readuradpar(filename_uradin)
c
c      reads parameters from file urad.in 
c      in a directory given by path
c
        Implicit real*8 (a-h, o-z)
        parameter (nsizeu=50000)
        double precision u(nsizeu),fa(nsizeu)
        double precision fv(nsizeu,256)
        double precision uct(256),ust(256),ucp(256),usp(256)
        double precision ru(100),zu(100)
        character*80 filename_uradin
        common/ur/u,fa,fv,ures,rsizeu,zsizeu,albedo,iscat
        common/ur2/ru,zu,nru,nzu,iurad,nru_used,nzu_used
        common/urang/icang,nang,uct,ust,ucp,usp
c
      write(6,2079)filename_uradin
2079  format(1h ,'READURADPAR: opening ',a) 
c
	open (unit=13, file=filename_uradin)
c
c       Reading the energy density steering file 'urad.in'
        read(13, *) iurad
        write(6,388)iurad
388     format(1h ,'iurad = ',i6)
        if(iurad .eq. 1)then
        write(6,209)iurad
 209     format(1h ,'READURADPAR: iurad = ',i6)
 	 read(13, *) ures
 	 read(13, *) icang
 	 read(13, *) rsizeu
 	 read(13, *) zsizeu
        write(6,201)ures
 201    format(1h ,'ures   = ',f10.2)
        write(6,202)icang
 202    format(1h ,'icang  = ',i6)
        write(6,203)rsizeu,zsizeu
 203    format(1h ,'rsizeu, zsizeu = ',2f10.2)
        endif
        if(iurad .eq. 2)then
 	 read(13, *) icang
c        write(6,2898)icang
c2898    format(1h ,'icang = ',i6)
         read(13, *) nru
c        write(6,2899)nru
c2899    format(1h ,'nru = ',i6)
         if(nru .gt. 100)then
          write(6,268)nru
268      format(1h ,'READURADPAR: INPUT error - nru out of bounds')
          goto 999
         endif        
         do 2988 i=1,nru
         read(13,*)ru(i)
c        write(6,*)ru(i)
 2988    continue
         read(13, *) nzu
c         write(6,289)nzu
c289      format(1h ,'nzu = ',i6)
         if(nzu .gt. 100)then
          write(6,269)nzu
269      format(1h ,'READURADPAR: INPUT error - nrz out of bounds')
          goto 999
         endif     
         do 2989 i=1,nzu
          read(13,*)zu(i)
c         write(6,*)zu(i)
 2989    continue 
        endif
c
999   continue
c
      close(unit=13)
c
      return
      end

