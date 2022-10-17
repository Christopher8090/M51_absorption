## geometry.in
The contents of `geometry.in` define the geometry of the:
- Stellar disk 
- Bulge
- Thin stellar disk
- Dust disk
- Thin dust disk

geometry.in contains the parameters for different morphological components which are indicated by the following numbers:
- Inner disks: 3/4
- Main disks: *blank*/1
- Outer disks: 5/6
- Nuclear disk: 7

The parameters for the disk components (stellar, thin stellar, dust and thin dust disks) are:
- tau: The central face-on dust opacity in the B band (only applicable to dust disks). This is effectively the amplitude parameter for the dust.
- hs/hd: Scale-length of the disk (hs for stellar disks, hd for dust disks).
- h_b/v/i/j/k/ir36/ir45/ir58: Scale-length of the old stellar disk in the B/V/I/J/K/ir36/ir45/i58 bands. The scale-length of the U band is equal to h_b.
- hsin/hdin: The inner radius of the disk (hsin for stellar disks, hdin for dust disks).
- hstin/hdtin: The inner truncation radius of the disk (hstin for stellar disks, hdtin for dust disks).
- rtruncate/rtruncated: The outer truncation radious of the disk (rtruncate for stellar disks, rtruncated for dust disks).
- hssolar/hdsolar: Radius at which the scale-height of the disk can be flared.
- zs/zd: Scale-height of the disk.
- zsin/zdin: Scale-height of the dist at R_e.
- zssolar/zdsolar: Scale-height of the disk at hssolar/hdsolar.
- sharp/sharpd: Smoothing factor at rtruncate/rtruncated.
- xis/xid: Parameter describing the linear slope from R=0 to hstin/hdtin.
- idisk: Boolean flag to include the component.

Due to the vast number of parameters it may be difficult to peice together the above explanations. To help elucidate the naming conventions a few examples are given:
- The central face-on dust opacity in the B band (tau) for the outer dust disk (5) would be 'tau5', for the thin dust disk (6) it would be 'tau6'.
- The inner truncation radius (hstin) for the inner thin stellar disk (4) would be 'hs4tin'.
- The scale-length (hs) of the main thin stellar disk (1) would be 'hs1'.


## scaling.in
scaling.in contains information about the model being run, including the dust model, whether the calculation includes scattering or not, as well as all the amplitude scaling factors.

## ff_scaling.in
ff_scaling.in contains a table of the wavelength dependent amplitude scaling factors used to fit the UV/optical/NIR profiles. N.B., the coeffeicients for uv22 (f_uv, f_uv4, f_uv6, f_uv7) and the B band (f_BVIK, f_BVIK3, f_BVIK5, f_bd) should all be kept equal to 1, to scale these wavelengths change the values in `scaling.in`.
