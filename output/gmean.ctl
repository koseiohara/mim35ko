dset ^gmean.grd
title MIM global mean

undef  -9.99E33

options little_endian yrev

xdef 1 linear 0.0  2.5 
ydef 1 linear 0.0  2.5
zdef 1 levels 1000 
tdef 77432 linear 00Z01JAN1949 6hr

VARS 2
az 1 99  GM(=Global Mean) Available Potential Energy [J/m^2]
qz 1 99  GM Available Zonal Diabatic Heating [W/m^2]
ENDVARS

