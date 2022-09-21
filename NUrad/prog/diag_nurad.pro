PRO diag_nurad
COMPILE_OPT IDL2						;set compile options
start
ss=''								;make a dummy string
;read gal_param
OPENR,gal_param_lun, "../indata/gal_param.in",  /GET_LUN	;open gal_param.in
WHILE ss NE 'prof_lim[kpc]' DO BEGIN				;find prof_lim[kpc], the radial extent of the galaxy
	READF, gal_param_lun, ss
ENDWHILE
READF, gal_param_lun, prof_lim					;read prof_lim
FREE_LUN, gal_param_lun						;close gal_param.in

indir='../out/profiles/'					;directory containing .dat diagnostic files, e.g. 'diskplotr_d_uv09.dat'.
outdir='../figures/nurad_diag_plots/'					;directory to output plots.

;arrays to make filenames
filter_arr=['b','v','i','j','k','ir36','ir45','ir58','uv09','uv13','uv15','uv16','uv20','uv22','uv25','uv28','uv36'];array of filters, looped through with j.
disk_arr=['bulge','disk','disk1','disk2','disk3','disk4','disk5','disk6','disk7']				    ;array of disks, looped through with k.
emiss_opac=['','tau']											       	    ;array of emissivity ('') and opacity ('tau'), looped through with l.
r_z=['r','z']								 		;array of radial ('r') and vertical ('z'), looped through with m.
r_z_in='r'
filter_in='b'
FOR i = 0,8 DO BEGIN
	disk_type = disk_arr[i]
	FOR p = 0,0 DO BEGIN
		emiss_type = emiss_opac[p]
		IF r_z_in EQ 'r' THEN m=0
		IF r_z_in EQ 'z' THEN m=1

		filename = indir+disk_type+emiss_type+'plot'+r_z[m]+'_'+filter_in+'.dat'
		file_exist=FILE_TEST(filename)	;check whether the file exists.
		IF file_exist EQ 0 THEN BEGIN									
			PRINT, "File not found: "+filename
			GOTO, skip;abort
		ENDIF
		READCOL, filename, x, y, SKIPLINE=1, FORMAT="(D,D)";read the data file, skip the first line (which is the number of lines in the file) and format the data as double-precision
		nan_index=WHERE(~FINITE(y))									;find any NaN values 
		IF N_ELEMENTS(nan_index) EQ N_ELEMENTS(y) THEN BEGIN
			PRINT, 'No data found in file: '+filename
			GOTO, skip	;if all of the y values are NaNs, skip the file
		ENDIF
		IF nan_index NE -1 THEN y[nan_index]=0.D		;set the NaN values to zero
		
   	outfile=outdir+emiss_type+disk_type+'_'+filter_in+'_'+r_z[m]+'.ps'	;output filename
	IF disk_type EQ 'bulge' THEN morph='b'
	IF disk_type EQ 'disk' AND emiss_type EQ '' THEN morph='d'
	IF disk_type EQ 'disk' AND emiss_type EQ 'tau' THEN morph='d'
	IF disk_type EQ 'disk1' AND emiss_type EQ 'tau' THEN morph='d'
        IF disk_type EQ 'disk1' AND emiss_type EQ '' THEN morph='td'
        IF disk_type EQ 'disk2' AND emiss_type EQ 'tau' THEN morph='td'
        IF disk_type EQ 'disk3' THEN morph='di'
        IF disk_type EQ 'disk4' THEN morph='tdi'
        IF disk_type EQ 'disk5' THEN morph='do'
        IF disk_type EQ 'disk6' THEN morph='tdo'
        IF disk_type EQ 'disk7' THEN morph='ntd'

	IF emiss_type EQ '' THEN emiss_opac2='emissivity'
        IF emiss_type EQ 'tau' THEN emiss_opac2='opacity'


       	outfile=outdir+morph+'_'+emiss_opac2+'_'+r_z[m]+'_'+filter_in+'.ps'	;output filename

		;Make the plot
		PRINT, 'Plotting: ', outfile
	
	r_in = 0
	r_out = 5
		aspect=cgPSWINDOW()			;set the aspect ratio using Coyote graphics
		thisDevice=!D.NAME			;set the device name
		SET_PLOT, 'PS', /COPY			;set the plot type to postscript and copy the default colour table
		DEVICE, FILENAME=outfile,$
                 	XSIZE=aspect.xsize, $
                       	YSIZE=aspect.ysize, $
                       	XOFFSET=aspect.xoffset, $
                       	YOFFSET=aspect.yoffset, $
                       	COLOR=1,$
                       	$av
                       	ENCAPSULATED=encapsulated, $
                       	Inches=aspect.inches, $
                       	bits=8, $
                       	set_character_size=[180,200], $
                       	set_font='HELVETICA', $
                       	/LANDSCAPE
                       	!p.thick=6
                       	!x.thick=3
                       	!y.thick=3
                       	!p.charthick=3
                       	!p.charsize=3

			;y=(100/MAX(y))*y			;scale the y-values as a percentage
			IF emiss_type EQ '' AND r_z[m] EQ 'r' THEN BEGIN
                       		CGPLOT, x, y,$
				XTITLE='Radius (kpc)',$
				YTITLE='Emissivity (percentage of maximum)',$
				;TITLE=morph+', '+filter_in+' band',$
				XRANGE=[0,1.2], yrange=[1e-15,1e-9], /ylog;1.2*prof_lim], $
			;	YRANGE=[-5,110];, /YLOG
			ENDIF
                        IF emiss_type EQ '' AND r_z[m] EQ 'z' THEN BEGIN
                        	CGPLOT, x, y,$
                                XTITLE='Height (kpc)',$
                                YTITLE='Emissivity (percentage of maximum)',$
                                TITLE=morph+', '+filter_in+' band',$
				XRANGE=[0,20];, $;1.2*prof_lim], $
			;	YRANGE=[-5,110]
                        ENDIF
                        IF emiss_type EQ 'tau' AND r_z[m] EQ 'r' THEN BEGIN
                        	CGPLOT, x, y,$
                                XTITLE='Radius (kpc)',$
                                YTITLE='Opacity (percentage of maximum)',$
                                ;TITLE=morph+', '+filter_in+' band',$S
				XRANGE=[r_in,r_out];, $;1.2*prof_lim], $
			;	YRANGE=[-5,110];, /YLOG
                        ENDIF
                        IF emiss_type EQ 'tau' AND r_z[m] EQ 'z' THEN BEGIN
                                CGPLOT, x, y,$
                                XTITLE='Height (kpc)',$
                                YTITLE='Opacity (percentage of maximum)',$
                                TITLE=morph+', '+filter_in+' band',$
				XRANGE=[0,10];, $;1.2*prof_lim], $
			;	YRANGE=[-5,110]
                        ENDIF                                                                                      

		        DEVICE, /CLOSE_FILE
                	cgFIXPS, outfile
skip:
ENDFOR
ENDFOR
abort:
END
