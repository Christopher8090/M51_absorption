;IDL include file
;end a plot started by 'plot_begin.inc'
;should be preceeded by:

;plot_name='myplot'
;@plot_begin.inc
;cgPlot, x, y

DEVICE, /CLOSE_FILE
cgFixPS, plot_name+'.ps'
