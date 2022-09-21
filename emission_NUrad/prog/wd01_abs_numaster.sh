#!/bin/bash

hostname
df /d1/vlk
source ~/.bashrc

model="'wd01'"
qyear="'06'"
tau=25.4
sfr=7.9
sfr4=5.8
sfr6=1.1
sfr7=0.
old=2.6
old3=1.4
old5=0.8
bd=0.01
scaabs="'abs'"
nsersic=4.
swdisk3="'yes'"
swdisk4="'yes'"
swdisk5="'yes'"
swdisk6="'yes'"
swdisk7="'no'"
ffactor=0.1
ffactor4=0.25
ffactor6=0.
ffactor7=0.
swheatr="'no'"
check="'no'"
dr_req=100.
dz_req=25.

PATH=$PATH:/d1/vlk/iso/idl/bin
cd /nfs/d58/vlk/sedmodel/cinman/m51a/emission_NUrad_M51a/prog
/d1/vlk/iso/idl/bin/idl_astro <<eof
start
@init
read_scaling,$model,$qyear,$scaabs,$tau,$nsersic,$sfr,$sfr4,$sfr6,$sfr7,$old,$old3,$old5,$bd,$ffactor,$ffactor4,$ffactor6,$ffactor7,$f_uv,$f_uv4,$f_uv6,$f_uv7,$f_BVIK,$f_BVIK3,$f_BVIK5,$f_bd
radiation_fields_template,$model,$qyear,$tau,$sfr,$sfr4,$sfr6,$sfr7,$old,$old3,$old5,$bd,$scaabs,$nsersic,$swdisk3,$swdisk4,$swdisk5,$swdisk6,$swdisk7
temp_trans_main_template,$model,$qyear,$tau,$sfr,$old,$bd,$scaabs,$check
luminosity_ec_template,$model,$qyear,$tau,$sfr,$old,$bd,$scaabs,$check
integrated_luminosity_template,$model,$qyear,$tau,$sfr,$sfr4,$sfr6,$sfr7,$old,$old3,$old5,$bd,$scaabs,$dr_req,$dz_req,$swdisk3,$swdisk4,$swdisk5,$swdisk6,$swheatr
local_luminosity,$model,$qyear,$tau,$sfr,$sfr4,$sfr6,$sfr7,$old,$bd,$ffactor,$ffactor4,$ffactor6,$ffactor7,$scaabs
local_emissivity_flat,$model,$qyear,$tau,$sfr,$sfr4,$sfr6,$sfr7,$old,$bd,$ffactor,$ffactor4,$ffactor6,$ffactor7,$scaabs
plot_total_luminosity_template,$model,$qyear,$tau,$sfr,$old,$bd,$ffactor,$scaabs,$swheatr
eof
