#!/bin/tcsh

hostname

set swdisk3="'yes'"
set swdisk4="'yes'"
set swdisk5="'yes'"
set swdisk6="'yes'"
set swdisk7="'no'"
set swheatr="'no'"
set check="'no'"
set dr_req=100.
set dz_req=25.

cd /net/triangulum/work/cjinman/M51_absorption/emission_NUrad/prog
/applications/idl/current/idl/bin/idl <<eof
read_scaling
@init
radiation_fields_template,$swdisk3,$swdisk4,$swdisk5,$swdisk6,$swdisk7
temp_trans_main_template,$check
luminosity_ec_template,$check
integrated_luminosity_template,$dr_req,$dz_req,$swdisk3,$swdisk4,$swdisk5,$swdisk6,$swheatr
local_luminosity
local_emissivity_flat
plot_total_luminosity_template,$swheatr
eof
