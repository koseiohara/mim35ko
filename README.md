# MIM v0.35.2-ko

This version is based on MIM v0.35 provided by Yuki Kanno (https://github.com/mim-proj/mim/tags) and modified by Kosei Ohara.  
It is easier to read namelist in this version than the original one.
Sample namelists are in this package.
Terms related to diabatic heating are decomposed into TTSWR (short wave radiation), TTLWR (long wave radiation), LRGHR (large scale condensation), CNVHR (convection), and VDFHR (vertical diffuction) and contribusions of each parameter are outputted.

## Features
- Namelist is read from files
- NaN is detected before output
- Program do not stop by too many warns
- Usage of Makefile is slightly changed
- q, qe, qgz in zonal and vint are decomposed into 5 parameters
- qz in gmean is decomposed into 5 parameters

## Samples of Namelists and Control Files
Their samples are for Japanese Reanalysis for Three Quarters of a Century (JRA3Q).

### input.nml
- INPUT\_U\_FILENAME
 &rarr; anl\_p125\_ugrd
- INPUT\_V\_FILENAME
 &rarr; anl\_p125\_vgrd
- INPUT\_T\_FILENAME
 &rarr; anl\_p125\_tmp
- INPUT\_PS\_FILENAME
 &rarr; anl\_surf125 (PRES)
- INPUT\_Z\_FILENAME
 &rarr; anl\_p125\_hgt
- INPUT\_OMEGA\_FILENAME
 &rarr; anl\_p125\_vvel
- INPUT\_TTSWR\_FILENAME
 &rarr; fcst\_phyp125\_ttswr
- INPUT\_TTLWR\_FILENAME
 &rarr; fcst\_phyp125\_ttlwr
- INPUT\_LRGHR\_FILENAME
 &rarr; fcst\_phyp125\_lrghr
- INPUT\_CNVHR\_FILENAME
 &rarr; fcst\_phyp125\_cnvhr
- INPUT\_VDFHR\_FILENAME
 &rarr; fcst\_phyp125\_vdfhr
- INPUT\_TOPO\_FILENAME
 &rarr; LL125\_surf
