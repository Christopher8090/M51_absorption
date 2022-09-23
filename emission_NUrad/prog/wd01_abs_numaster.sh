#!/bin/tcsh

hostname

set model="'wd01'"
set qyear="'06'"
set tau=25.4
set sfr=7.9
set sfr4=5.8
set sfr6=1.1
set sfr7=0.
set old=2.6
set old3=1.4
set old5=0.8
set bd=0.01
set scaabs="'abs'"
set nsersic=4.
set swdisk3="'yes'"
set swdisk4="'yes'"
set swdisk5="'yes'"
set swdisk6="'yes'"
set swdisk7="'no'"
set ffactor=0.1
set ffactor4=0.25
set ffactor6=0.
set ffactor7=0.
set swheatr="'no'"
set check="'no'"
set dr_req=100.
set dz_req=25.

cd /net/triangulum/work/cjinman/M51_absorption/emission_NUrad/prog
/applications/idl/current/idl/bin/idl <<eof
read_scaling
@init
radiation_fields_template,$model,$qyear,$tau,$sfr,$sfr4,$sfr6,$sfr7,$old,$old3,$old5,$bd,$scaabs,$nsersic,$swdisk3,$swdisk4,$swdisk5,$swdisk6,$swdisk7
temp_trans_main_template,$model,$qyear,$tau,$sfr,$old,$bd,$scaabs,$check
luminosity_ec_template,$model,$qyear,$tau,$sfr,$old,$bd,$scaabs,$check
integrated_luminosity_template,$model,$qyear,$tau,$sfr,$sfr4,$sfr6,$sfr7,$old,$old3,$old5,$bd,$scaabs,$dr_req,$dz_req,$swdisk3,$swdisk4,$swdisk5,$swdisk6,$swheatr
local_luminosity,$model,$qyear,$tau,$sfr,$sfr4,$sfr6,$sfr7,$old,$bd,$scaabs
local_emissivity_flat,$model,$qyear,$tau,$sfr,$sfr4,$sfr6,$sfr7,$old,$bd,$scaabs
plot_total_luminosity_template,$model,$qyear,$tau,$sfr,$old,$bd,$ffactor,$scaabs,$swheatr
eof
