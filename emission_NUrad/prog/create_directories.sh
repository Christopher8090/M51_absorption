#!/bin/tcsh

if ( ! -d "../outdata_intlum/") then
	mkdir ../outdata_intlum/
	echo "create directory: ../outdata_intlum/"
endif

if ( ! -d "../outdata_lum/") then
        mkdir ../outdata_lum/
        echo "create directory: ../outdata_lum/"
endif

if ( ! -d "../outdata_temp/") then
        mkdir ../outdata_temp/
        echo "create directory: ../outdata_temp/"
endif

if ( ! -d "../indata/") then
        mkdir ../indata/
        echo "create directory: ../indata/"
endif

if ( ! -d "../figures/") then
	mkdir ../figures/
	echo "create directory: ../figures/"
endif
