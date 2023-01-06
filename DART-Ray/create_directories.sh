#!/bin/tcsh

if ( ! -d "./MAPS/") then
	mkdir ./MAPS/
	echo "create directory: ./MAPS/"
endif

if ( ! -d "./MAPS/M51a/") then
        mkdir ./MAPS/M51a/
        echo "create directory: ./MAPS/M51a/"
endif

if ( ! -d "./RUNS/") then
        mkdir ./RUNS/
        echo "create directory: ./RUNS/"
endif

if ( ! -d "./RUNS/GALAXY_M51a/") then
        mkdir ./RUNS/GALAXY_M51a/
        echo "create directory: ./RUNS/GALAXY_M51a/"
endif

