dset ^JRA3Q_1975_2022_ALL_GMEAN_FEWER.grd
title MIM global mean

undef  9.99E20

options little_endian yrev

xdef 1 linear 0.0  2.5 
ydef 1 linear 0.0  2.5
zdef 1 levels 1000 
tdef 70128 linear 00Z01JAN1975 6hr

VARS 7
az       1 99  GM(=Global Mean) Available Potential Energy [J/m^2]
qz       1 99  GM Available Zonal Diabatic Heating [W/m^2]
ttswr_qz 1 99  GM Available Zonal Diabatic Heating [W/m^2] by short wave radiation
ttlwr_qz 1 99  GM Available Zonal Diabatic Heating [W/m^2] by long wave radiation
lrghr_qz 1 99  GM Available Zonal Diabatic Heating [W/m^2] by large scale condensation
cnvhr_qz 1 99  GM Available Zonal Diabatic Heating [W/m^2] by convection
vdfhr_qz 1 99  GM Available Zonal Diabatic Heating [W/m^2] by vertical diffusion
ENDVARS

