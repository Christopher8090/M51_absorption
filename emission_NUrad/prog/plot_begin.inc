;IDL include file.
;This is a standard preamble to make .ps plots look nice.
;Requires plot_name to be defined as a string with no file extension
;Requires the following line after 'DEVICE, /CLOSE_FILE'
;cgFixPS, plot_name
;For use with Coyote graphics library.


;plot_name='myplot'

aspect = cgPSWINDOW()
thisDevice = !D.NAME
PRINT, 'Plotting: '+plot_name+'.ps'
SET_PLOT, 'PS', /Copy
DEVICE, FILENAME=plot_name+'.ps', $
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
        !p.thick=3
        !x.thick=3
        !y.thick=3
        !p.charthick=2
        !p.charsize=1.5

;cgPlot, x, y

;DEVICE, /CLOSE_FILE
;cgFixPS, plot_name+'.ps'
