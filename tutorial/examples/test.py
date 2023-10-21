TairC = 25
RH = 0.5
hPa2cm = 1.0197
dEauPure = (999.83952 + TairC * (16.952577 + TairC *
                                 (- 0.0079905127 + TairC * (- 0.000046241757 + TairC *
                                                            (0.00000010584601 + TairC * (
                                                                - 0.00000000028103006)))))) / (
                       1 + 0.016887236 * TairC)
siPhi = (30 - TairC) / (91 + TairC)
siEnne = 0
mu = pow(10, (- 0.114 + (siPhi * (1.1 + 43.1 * pow(siEnne, 1.25)))))
mu = mu / (
            24 * 60 * 60) / 100 / 1000;  # //mPa s to hPa d, 1.11837e-10 hPa d for pure water at 293.15K
mu = mu * hPa2cm  # hPa d to cmh2o d

# number of vascular bundles
VascBundle_leaf = 32
VascBundle_stem = 52
VascBundle_root = 1  # valid for all root type

# radius of xylem type^4 * number per bundle
rad_x_l_1 = (0.0015 ** 4) * 2;
rad_x_l_2 = (0.0005 ** 4) * 2
rad_x_s_1 = (0.0017 ** 4) * 3;
rad_x_s_2 = (0.0008 ** 4) * 1
rad_x_r0_1 = (0.0015 ** 4) * 4
rad_x_r12_1 = (0.00041 ** 4) * 4;
rad_x_r12_2 = (0.00087 ** 4) * 1
rad_x_r3_1 = (0.00068 ** 4) * 1

# axial conductivity [cm^3/day]
kz_l = VascBundle_leaf * (rad_x_l_1 + rad_x_l_2) * np.pi / (mu * 8)
kz_s = VascBundle_stem * (rad_x_s_1 + rad_x_s_2) * np.pi / (mu * 8)
kz_r0 = VascBundle_root * rad_x_r0_1 * np.pi / (mu * 8)
kz_r1 = VascBundle_root * (rad_x_r12_1 + rad_x_r12_2) * np.pi / (mu * 8)
kz_r2 = VascBundle_root * (rad_x_r12_1 + rad_x_r12_2) * np.pi / (mu * 8)
kz_r3 = VascBundle_root * rad_x_r3_1 * np.pi / (mu * 8)  # 4.32e-1
