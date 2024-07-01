dset ^JRA55_2000_2000_ALL_ZONAL_FEWER.grd
title MIM

undef 9.999e+20

options little_endian yrev

xdef 1 linear 1 1
ydef 145 linear -90.0 1.25
zdef 37 levels
1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 225, 200, 175, 150, 125, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1
tdef 1464 linear 00Z01JAN2000 6hr

vars 59
u           37 99  Zonal Wind [m/s]
v           37 99  Meridional Wind [m/s]
pt          37 99  Potential Temperature [K]
t           37 99  Temperature Dagger derived from pt [K]
st          37 99  Mass Streamfunction [kg/s]
w           37 99  Vertical Velocity derived from st [m/s]
z           37 99  Geopotential Height [m]
epy         37 99  Meridional Component of EP Flux [kg/s^2]
depy        37 99  EP Flux Divergence dut to epy [m/s^2]
epz_form    37 99  Vertical Component of EP Flux (Form drag) [kg/s^2]
depz_form   37 99  EP Flux Divergence dut to epz_form [m/s^2]
epz_w       37 99  Vertical Component of EP Flux (W) [kg/s^2]
depz_w      37 99  EP Flux Divergence dut to epz_w [m/s^2]
epz_ut      37 99  Vertical Component of EP Flux (u'T') [kg/s^2]
depz_ut     37 99  EP Flux Divergence dut to epz_ut [m/s^2]
epz         37 99  Vertical Component of EP Flux (Total) [kg/s^2]
depz        37 99  EP Flux Divergence dut to epz [m/s^2]
divf        37 99  EP Flux Divergence [m/s^2]
gy          37 99  Meridional Comp. of G (Meridional Momentum) Flux [kg/s^2]
dgy         37 99  G Flux Divergence due to gy [m/s^2]
gz          37 99  Vertical Component of G (Meridional Momentum) Flux [kg/s^2]
dgz         37 99  G Flux Divergence due to gz [m/s^2]
uux         37 99  Zonal mean (u'u') [m^2/s^2]
c_az_kz     37 99  C(Az->Kz) [J/(kg s)] or [W/kg]
c_kz_ae_u   37 99  C(Kz->Ae) by U (i.e. form drag) [J/(kg s)] or [W/kg]
c_kz_ae_v   37 99  C(Kz->Ae) by V [J/(kg s)] or [W/kg]
c_kz_ae     37 99  C(Kz->Ae) [J/(kg s)] or [W/kg]
c_ae_ke_u   37 99  C(Ae->Ke) by U [J/(kg s)] or [W/kg]
c_ae_ke_v   37 99  C(Ae->Ke) by V [J/(kg s)] or [W/kg]
c_ae_ke     37 99  C(Ae->Ke) [J/(kg s)] or [W/kg]
c_kz_ke_uy  37 99  C(Kz->Ke) by u * DFy [J/(kg s)] or [W/kg]
c_kz_ke_uz  37 99  C(Kz->Ke) by u * DFz^uw [J/(kg s)] or [W/kg]
c_kz_ke_vy  37 99  C(Kz->Ke) by v * DGy [J/(kg s)] or [W/kg]
c_kz_ke_vz  37 99  C(Kz->Ke) by v * DGz [J/(kg s)] or [W/kg]
c_kz_ke_tan 37 99  C(Kz->Ke) by tan [J/(kg s)] or [W/kg]
c_kz_ke     37 99  C(Kz->Ke) [J/(kg s)] or [W/kg]
c_kz_w      37 99  C(Kz->W) = C(Kz->Ke) + C(Kz->Ae) [J/(kg s)] or [W/kg]
q           37 99  Diabatic Heating (q/cp -> dT/dt) [J/(kg s)]
ttswr       37 99  Diabatic Heating (q/cp -> dT/dt) [J/(kg s)] by short wave radiation
ttlwr       37 99  Diabatic Heating (q/cp -> dT/dt) [J/(kg s)] by long wave radiation
lrghr       37 99  Diabatic Heating (q/cp -> dT/dt) [J/(kg s)] by large scale condensation
cnvhr       37 99  Diabatic Heating (q/cp -> dT/dt) [J/(kg s)] by convection
vdfhr       37 99  Diabatic Heating (q/cp -> dT/dt) [J/(kg s)] by vertical diffusion
qgz         37 99  Zonal Diabatic Heating (+Ground State) [J/(kg s)]
ttswr_gz    37 99  Zonal Diabatic Heating (+Ground State) [J/(kg s)] by short wave radiation
ttlwr_gz    37 99  Zonal Diabatic Heating (+Ground State) [J/(kg s)] by long wave radiation
lrghr_gz    37 99  Zonal Diabatic Heating (+Ground State) [J/(kg s)] by large scale condensation
cnvhr_gz    37 99  Zonal Diabatic Heating (+Ground State) [J/(kg s)] by convection
vdfhr_gz    37 99  Zonal Diabatic Heating (+Ground State) [J/(kg s)] by vertical diffusion
qe          37 99  Eddy Available Diabatic Heating [J/(kg s)]
ttswr_qe    37 99  Eddy Available Diabatic Heating [J/(kg s)] by short wave radiation
ttlwr_qe    37 99  Eddy Available Diabatic Heating [J/(kg s)] by long wave radiation
lrghr_qe    37 99  Eddy Available Diabatic Heating [J/(kg s)] by large scale condensation
cnvhr_qe    37 99  Eddy Available Diabatic Heating [J/(kg s)] by convection
vdfhr_qe    37 99  Eddy Available Diabatic Heating [J/(kg s)] by vertical diffusion
kz          37 99  Zonal Kinetic Energy [m^2/s^2] or [J/kg]
ke          37 99  Eddy Kinetic Energy [m^2/s^2] or [J/kg]
pz          37 99  Potential Energy (including groud state) [m^2/s^2] or [J/kg]
ae_total    37 99  Eddy Available Potential Energy (including surface effect)
ENDVARS
