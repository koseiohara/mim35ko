dset ^JRA3Q_1975_2023_ALL_VINT_FEWER.grd
title MIM vertical integral

undef 9.999e+20

options little_endian yrev

xdef 1 linear 1 1
ydef 145 linear -90.0 1.25
zdef 1 levels 1000
tdef 71588 linear 00Z01JAN1975 6hr

VARS 36
kz           1 99  VI (=Vertically Integrated) Zonal Kinetic Energy [J/m^2]
ke           1 99  VI Eddy Kinetic Energy [J/m^2]
az           1 99  VI Zonal Available Potential Energy (not for gmean) [J/m^2]
ae           1 99  VI Eddy Available Potential Energy [J/m^2]
c_az_kz      1 99  VI C(Az->Kz) [W/m^2]
c_kz_pe_u    1 99  VI C(Kz->Ae) by U (i.e. form drag) [W/m^2]
c_kz_ae_v    1 99  VI C(Kz->Ae) by V [W/m^2]
c_kz_ae      1 99  VI C(Kz->Ae) [W/m^2]
c_ae_ke_u    1 99  VI C(Ae->Ke) by U [W/m^2]
c_ae_ke_v    1 99  VI C(Ae->Ke) by V [W/m^2]
c_ae_ke      1 99  VI C(Ae->Ke) [W/m^2]
c_kz_ke_uy   1 99  VI C(Kz->Ke) by u * DFy [W/m^2]
c_kz_ke_uz   1 99  VI C(Kz->Ke) by u * DFz [W/m^2]
c_kz_ke_vy   1 99  VI C(Kz->Ke) by u * DGy [W/m^2]
c_kz_ke_vz   1 99  VI C(Kz->Ke) by u * DGz [W/m^2]
c_kz_ke_tan  1 99  VI C(Kz->Ke) by tan [W/m^2]
c_kz_ke      1 99  VI C(Kz->Ke) [W/m^2]
c_kz_w       1 99  VI C(Kz->W) = C(Kz->Ke) + C(Kz->Ae) [W/m^2]
q            1 99  Diabatic Heating [W/m^2]
ttswr        1 99  Diabatic Heating [W/m^2] by short wave radiation
ttlwr        1 99  Diabatic Heating [W/m^2] by long wave radiation
lrghr        1 99  Diabatic Heating [W/m^2] by large scale condensation
cnvhr        1 99  Diabatic Heating [W/m^2] by convection
vdfhr        1 99  Diabatic Heating [W/m^2] by vertical diffusion
qgz          1 99  Zonal Diabatic Heating (+Ground State) [W/m^2]
ttswr_gz     1 99  Zonal Diabatic Heating (+Ground State) [W/m^2] by short wave radiation
ttlwr_gz     1 99  Zonal Diabatic Heating (+Ground State) [W/m^2] by long wave radiation
lrghr_gz     1 99  Zonal Diabatic Heating (+Ground State) [W/m^2] by large scale condensation
cnvhr_gz     1 99  Zonal Diabatic Heating (+Ground State) [W/m^2] by convection
vdfhr_gz     1 99  Zonal Diabatic Heating (+Ground State) [W/m^2] by vertical diffusion
qe           1 99  Eddy Available Diabatic Heating [W/m^2]
ttswr_qe     1 99  Eddy Available Diabatic Heating [W/m^2] by short wave radiation
ttlwr_qe     1 99  Eddy Available Diabatic Heating [W/m^2] by long wave radiation
lrghr_qe     1 99  Eddy Available Diabatic Heating [W/m^2] by large scale condensation
cnvhr_qe     1 99  Eddy Available Diabatic Heating [W/m^2] by convection
vdfhr_qe     1 99  Eddy Available Diabatic Heating [W/m^2] by vertical diffusion
ENDVARS
