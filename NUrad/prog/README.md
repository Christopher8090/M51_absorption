# NUrad

These codes calculate the radiation fields according to the files in `/NUrad/indata/`. The order of operation is as follow:

1. Ensure the `NUrad` program is properly compiled, run `Makefile` to compile this program.
2. Run `Nurad` for each morphological component and each wavelength. All components are calculated by running: `master.sh`.
3. Now that the uncalibrated radiation fields have been calculated, they can be calibrated by running `radition_fields_unit.pro`. The surface brightness maps are made by running `maps.pro`. This must be done for each morphological omponent which can be done by simply running `twod_all.pro`.
4. To produce the azimuthally averaged surface brightness profiles, run `twod_phot.pro`. This produces `.save` files in /saves/model/ and also plots in /figures/
