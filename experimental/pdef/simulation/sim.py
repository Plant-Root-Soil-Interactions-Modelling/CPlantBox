import sys;
CPBdir = "../../.."
sys.path.append(CPBdir + "/src");
sys.path.append(CPBdir);
sys.path.append("../../..");
sys.path.append("..");
sys.path.append(CPBdir + "/src/python_modules");
sys.path.append("../build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../modules/")  # python wrappers
sys.path.append("../modules/functional/")  # python wrappers
sys.path.append("../../experimental/photosynthesis/")

from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import visualisation.vtk_plot as vp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


MaxSimtime = 28

""" Parameters """
def sinusoidal(t):
    return (np.sin(np.pi*t*2)+1)/2

# Environment
coefhours = sinusoidal(MaxSimtime)
Tmax = 24
Tmin = 18
TairC = Tmin + (Tmax - Tmin) * coefhours
hPa2cm = 1.0197
siPhi = (30 - TairC) / (91 + TairC)
siEnne=0
mu =  pow(10, (- 0.114 + (siPhi * (1.1 + 43.1 * pow(siEnne, 1.25) ))))
mu = mu /(24*60*60)/100/1000; #//mPa s to hPa d, 1.11837e-10 hPa d for pure water at 293.15K
mu = mu * hPa2cm #hPa d to cmh2o d
#g/ml
dEauPure = (999.83952 + TairC * (16.952577 + TairC *
        (- 0.0079905127 + TairC * (- 0.000046241757 + TairC *
        (0.00000010584601 + TairC * (- 0.00000000028103006)))))) /  (1 + 0.016887236 *
                                                                     TairC)/1000
molarMassWater = 18 #g/mol
mmol2cm3 = molarMassWater/dEauPure /1000#cm3/mmol

a = 1.185 #mean stem radi from all P levels

kx_stem = (np.pi*(a*0.017)**4)/(8*mu)

kr_s = np.array([[0, 0], [1e4, 0]])
kz_s = np.array([[0, kx_stem], [1e4, kx_stem]])
kr_l_temp = 25

s2d = 1/3600/24 #[d/s]
MPa2cm = hPa2cm * 1000 #[cm/hPa * hPa/MPa]
m22cm2 = 10000 #cm2/m2


#25 mmol/m2/s/MPa
kr_l_temp2 = kr_l_temp*mmol2cm3 /m22cm2/s2d/MPa2cm #cm2/cm3/d

kr_l = np.array([[0, 3.83e-6 * hPa2cm], [1e4, 1000 * hPa2cm]])
kz_l = np.array([[0, kx_stem/10], [1e4, kx_stem/10]])

kr0 = kr_l
kz0 = kz_s

'''Root conductivities with exponential growth'''
kr1 = np.array([
    [0.10748958, 0.000106343],
    [0.161379976, 0.000106008],
    [0.215367912, 0.000106236],
    [0.26945374, 0.000106033],
    [0.323637817, 0.000106347],
    [0.3779205, 0.00010643],
    [0.432302149, 0.000106263],
    [0.650825663, 0.000106385],
    [0.87096191, 0.000106457],
    [1.092734874, 0.000105761],
    [1.541289585, 0.000113566],
    [1.768122056, 0.00011429],
    [1.996692724, 4.00E-05],
    [2.227028433, 3.99E-05],
    [2.459156657, 4.00E-05],
    [2.693105512, 3.99E-05],
    [2.928903785, 4.00E-05],
    [3.166580949, 3.68E-05],
    [3.406167186, 3.68E-05],
    [4.633843597, 3.89E-05],
    [7.252043149, 4.10E-05],
    [10.12236428, 4.29E-05]
])

kz1 = np.array([
    [0.10748958, 0.001428701],
    [0.161379976, 0.001411575],
    [0.215367912, 0.001503535],
    [0.26945374, 0.001474202],
    [0.323637817, 0.001372298],
    [0.3779205, 0.001462689],
    [0.432302149, 0.001406009],
    [0.650825663, 0.001467545],
    [0.87096191, 0.001407742],
    [1.092734874, 0.00151302],
    [1.541289585, 0.001336126],
    [1.768122056, 0.001350818],
    [1.996692724, 0.001348176],
    [2.227028433, 0.00195872],
    [2.693105512, 0.161388966],
    [2.928903785, 0.163310872],
    [3.166580949, 0.171753617],
    [3.406167186, 0.173627719],
    [4.633843597, 0.161798441],
    [7.252043149, 0.243831044],
    [10.12236428, 0.333556853]
])

#seminal
kr4 = np.array([
    [0.12909831, 0.000112594],
    [0.193822343, 0.000112978],
    [0.258663523, 0.000113248],
    [0.323622275, 0.000112912],
    [0.388699028, 0.000113171],
    [0.453894209, 0.000113439],
    [0.519208251, 0.000113637],
    [0.781661749, 0.000121677],
    [1.046052191, 0.00012169],
    [1.312408379, 0.000122518],
    [1.580759763, 0.00012298],
    [1.85113646, 0.000122788],
    [2.123569273, 0.000123039],
    [2.398089714, 0.000123636],
    [2.674730025, 0.0001233],
    [2.953523201, 4.22E-05],
    [3.234503012, 4.23E-05],
    [3.517704031, 4.24E-05],
    [3.803161655, 4.25E-05],
    [4.090912136, 4.26E-05],
    [5.565389475, 3.93E-05],
    [8.709928112, 4.16E-05],
    [12.15727256, 4.31E-05]
])

kz4 = np.array([
    [0.12909831, 0.001043448],
    [0.193822343, 0.001084252],
    [0.258663523, 0.001052552],
    [0.323622275, 0.001103842],
    [0.388699028, 0.00108313],
    [0.453894209, 0.001067554],
    [0.519208251, 0.001006955],
    [0.781661749, 0.000973576],
    [1.046052191, 0.000976876],
    [1.312408379, 0.000985499],
    [1.580759763, 0.000907123],
    [1.85113646, 0.000865667],
    [2.123569273, 0.000801025],
    [2.398089714, 0.000818666],
    [2.674730025, 0.000870364],
    [2.953523201, 0.000843755],
    [3.234503012, 0.000833276],
    [3.517704031, 0.000911009],
    [3.803161655, 0.000875389],
    [4.090912136, 0.000898246],
    [5.565389475, 0.112154361],
    [8.709928112, 0.101482174],
    [12.15727256, 0.150518048]
])

kr2 = np.array([
    [0.243904673, 0.000118148],
    [0.375243914, 0.000117727],
    [0.513747967, 0.000117599],
    [0.660243877, 0.000118013],
    [0.815710704, 0.000117724],
    [0.981319231, 0.000117804],
    [1.158485533, 0.000117637],
    [2.024739266, 0.000117638],
    [3.349677885, 0.000117924],
    [6.282752364, 3.75E-05]
])

kz2 = np.array([
    [0.243904673, 4.48E-05],
    [0.375243914, 4.63E-05],
    [0.513747967, 3.62E-05],
    [0.660243877, 4.87E-05],
    [0.815710704, 3.97E-05],
    [0.981319231, 4.29E-05],
    [1.158485533, 4.66E-05],
    [2.024739266, 4.22E-05],
    [3.349677885, 3.86E-05],
    [6.282752364, 4.21E-05]
])

kr3 = np.array([
    [0.296274974, 0.000129391],
    [0.530053444, 0.000129699],
    [0.92489003, 0.000129233],
    [2.834286141, 0.000129567]
])

kz3 = np.array([
    [0.296274974, 1.49E-05],
    [0.530053444,1.60E-05],
    [0.92489003, 1.58E-05]])

# shoot-born
kr5 = np.array([
    [0.107343906, 0.000133676],
    [0.161161267, 0.000133827],
    [0.215076036, 0.000133532],
    [0.269088565, 0.000133534],
    [0.323199209, 0.000134758],
    [0.377408326, 0.000134415],
    [0.431716275, 0.000135286],
    [0.649943636, 0.000135067],
    [0.869781546, 0.000140551],
    [1.091253953, 0.000141006],
    [1.314385345, 0.00014217],
    [1.539200764, 0.000142592],
    [1.765725822, 0.000147261],
    [1.993986721, 0.000148067],
    [2.22401027, 0.000148443],
    [2.455823903, 0.00013956],
    [2.689455702, 0.000144219],
    [2.924934411, 5.28E-05],
    [3.162289464, 5.28E-05],
    [3.401551004, 5.31E-05],
    [4.627563616, 5.08E-05],
    [7.242214873, 5.24E-05],
    [10.10864602, 5.67E-05]
])

kz5 = np.array([
    [0.107343906, 0.006101211],
    [0.161161267, 0.006128716],
    [0.215076036, 0.006138274],
    [0.269088565, 0.006037813],
    [0.323199209, 0.006311727],
    [0.377408326, 0.005902104],
    [0.431716275, 0.005960428],
    [0.649943636, 0.005785935],
    [0.869781546, 0.005770295],
    [1.091253953, 0.006697144],
    [1.314385345, 0.007411226],
    [1.539200764, 0.006981372],
    [1.765725822, 0.006388351],
    [1.993986721, 0.00680202],
    [2.22401027, 0.006945852],
    [2.689455702, 1.427732171],
    [2.924934411, 1.516871436],
    [3.162289464, 1.607222779],
    [3.401551004, 1.691691634],
    [4.627563616, 2.344891587],
    [7.242214873, 4.094611567],
    [10.10864602, 6.608610426]
])

dt = 1.
steps = round(MaxSimtime / dt)  # steps
runs = 100

krs_P0,krs_P1,krs_P2,krs_P3, count  = [], [],[], [],[]
krshoot_P0,krshoot_P1,krshoot_P2,krshoot_P3 = [], [],[], []
krplant_P0,krplant_P1,krplant_P2,krplant_P3 = [], [],[], []
jc_P0,jc_P1,jc_P2,jc_P3  = [], [],[], []
eswp_P0,eswp_P1,eswp_P2,eswp_P3  = [], [],[], []
eawp_P0,eawp_P1,eawp_P2,eawp_P3  = [], [],[], []
cwp_P0,cwp_P1,cwp_P2,cwp_P3  = [], [],[], []
vol_P0_, vol_P1_, vol_P2_, vol_P3_ = [], [], [], []
len_P0_, len_P1_, len_P2_, len_P3_ = [], [], [], []
sur_P0_, sur_P1_, sur_P2_, sur_P3_ = [], [], [], []

for i in range(0,runs):
    allRS = []
    for i in range(0, 4):
        name = 'P'+str(i)
        rs = pb.MappedPlant()
        if name == 'P0':
            for p in rs.getOrganRandomParameter(pb.leaf):
                p.la,  p.lmax = 38.41053981, 38.41053981
                p.areaMax = 54.45388021  # cm2, area reached when length = lmax
                NLeaf = 100  
                phi = np.array([-90,-80, -45, 0., 45, 90]) / 180. * np.pi    
                l = np.array([38.41053981,1 ,1, 0.3, 1, 38.41053981]) #distance from leaf center
                p.tropismT = 1 # 6: Anti-gravitropism to gravitropism
                p.tropismS = 0.05
                p.tropismAge = 5 #< age at which tropism switch occures, only used if p.tropismT = 6
                p.createLeafRadialGeometry(phi,l,NLeaf)
            for p in rs.getOrganRandomParameter(pb.stem):
                p.r = 1

        if name == 'P1':
            for p in rs.getOrganRandomParameter(pb.leaf):
                p.lb =  0 # length of leaf stem
                p.la,  p.lmax = 42.60617256, 42.60617256
                p.areaMax = 66.69532685  # cm2, area reached when length = lmax
                NLeaf = 100  
                phi = np.array([-90,-80, -45, 0., 45, 90]) / 180. * np.pi
                l = np.array([42.60617256,1 ,1, 0.3, 1, 42.60617256]) #distance from leaf center
                p.tropismT = 1 # 6: Anti-gravitropism to gravitropism
                #p.tropismN = 5
                p.tropismS = 0.05
                p.tropismAge = 5 #< age at which tropism switch occures, only used if p.tropismT = 6v
                p.createLeafRadialGeometry(phi, l, NLeaf)
            for p in rs.getOrganRandomParameter(pb.stem):
                p.r= 1

        if name == 'P2':
            for p in rs.getOrganRandomParameter(pb.leaf):
                p.lb =  0 # length of leaf stem
                p.la,  p.lmax = 52.23664394, 52.23664394
                p.areaMax = 80.68274258  # cm2, area reached when length = lmax
                NLeaf = 100  
                phi = np.array([-90,-80, -45, 0., 45, 90]) / 180. * np.pi
                l = np.array([52.23664394,1 ,1, 0.3, 1, 52.23664394]) #distance from leaf center
                p.tropismT = 1 # 6: Anti-gravitropism to gravitropism
                #p.tropismN = 5
                p.tropismS = 0.05
                p.tropismAge = 5 #< age at which tropism switch occures, only used if p.tropismT = 6
                p.createLeafRadialGeometry(phi, l, NLeaf)
            for p in rs.getOrganRandomParameter(pb.stem):
                p.r= 1

        if name == 'P3':
            for p in rs.getOrganRandomParameter(pb.leaf):
                p.lb =  0 # length of leaf stem
                p.la,  p.lmax = 49.12433414, 49.12433414
                p.areaMax = 71.95670914  # cm2, area reached when length = lmax
                NLeaf = 100  
                phi = np.array([-90,-80, -45, 0., 45, 90]) / 180. * np.pi
                l = np.array([49.12433414,1 ,1, 0.3, 1, 49.12433414]) #distance from leaf center
                p.tropismT = 1 # 6: Anti-gravitropism to gravitropism
                p.tropismN = 5
                p.tropismS = 0.05
                p.tropismAge = 5 #< age at which tropism switch occures, only used if p.tropismT = 6
                p.createLeafRadialGeometry(phi, l, NLeaf)
            for p in rs.getOrganRandomParameter(pb.stem):
                p.r= 1

        rs.setGeometry(pb.SDF_PlantBox(1.e6, 1.e6, 1.e6))  # not allowed to grow upwards out of soil
        p_s = np.linspace(-500, -2000, 30001)  #  
        soil_index = lambda x, y, z : 0
        p_s = -200  # static soil pressure [cm]
        p_a = -15000 #static air pressure
        rs.setSoilGrid(soil_index)
        rs.readParameters( name + ".xml")
        rs.initialize()  
        allRS.append(rs)
    ö = 0
    for rs in allRS:
        ö +=1
        krs_, krshoot_, krplant_, suf_, jc_, ti, eswp_,eawp_, cwp_, vol_, len_, sur_ = [], [], [], [], [], [], [], [],[],[], [], []
        k_soil = []
        simtime = 0
        r = XylemFluxPython(rs)
        r.setKrTables([[ kr1[:, 1], kr2[:, 1], kr3[:, 1], kr4[:, 1], kr5[:, 1]],[kr_s[:, 1],kr_s[:, 1]],[kr_l[:, 1]]],
                [[kr1[:, 0], kr2[:, 0], kr3[:, 0], kr4[:, 0], kr5[:, 0]],[kr_s[:, 0],kr_s[:, 0]],[kr_l[:, 0]]])
        r.setKxTables([[kz1[:, 1], kz2[:, 1], kz3[:, 1], kz4[:, 1], kz5[:, 1]],[kz_s[:, 1],kz_s[:, 1]],[kz_l[:, 1]]],
                [[kz1[:, 0], kz2[:, 0], kz3[:, 0], kz4[:, 0], kz5[:, 0]],[kz_s[:, 0],kz_s[:, 0]],[kz_l[:, 0]]])
        for j in range(steps):
            r.rs.simulate(dt, False)
            simtime += dt

            """ set up xylem parameters """
            r.airPressure = p_a
            suf = r.get_suf(j)
            print('SUF:', suf)
            suf_.append(suf)
            krs,krshoot,krplant, jc, eswp, eawp,cwp = r.get_krs(j, plant = True)
            print("P nr", ö,"simtime",int(simtime))
            print("\tkrs",np.around(krs,3))
            print("\tsum of SUF",np.around( np.sum(suf),2), "summed positive",np.around( np.sum(suf[suf >= 0]),2) )
            vol = np.sum(np.array(rs.getParameter('volume')))
            len = np.sum(np.array(rs.getParameter('length')))
            sur = np.sum(np.array(rs.getParameter('surface')))
            print('Volume of root:',vol)
            print('Length of root:',len)
            print('Surface of root:',sur)

            krs_.append(krs)
            krshoot_.append(krshoot)
            krplant_.append(krplant)
            jc_.append(jc)
            eswp_.append(eswp)
            eawp_.append(eawp)
            cwp_.append(cwp)
            ti.append(simtime)
            suf_.append(suf)
            vol_.append(vol)
            len_.append(len)
            sur_.append(sur)

        """ plot """
        if ö == 1:
            count.extend(ti)
            krs_P0.extend(krs_)
            krshoot_P0.extend(krshoot_)
            krplant_P0.extend(krplant_)
            jc_P0.extend(jc_)
            eswp_P0.extend(eswp_)
            eawp_P0.extend(eawp_)
            cwp_P0.extend(cwp_)
            vol_P0_.extend(vol_)
            len_P0_.extend(len_)
            sur_P0_.extend(sur_)
        elif ö == 2:
            krs_P1.extend(krs_)
            krshoot_P1.extend(krshoot_)
            krplant_P1.extend(krplant_)
            jc_P1.extend(jc_)
            eswp_P1.extend(eswp_)
            eawp_P1.extend(eawp_)
            cwp_P1.extend(cwp_)
            vol_P1_.extend(vol_)
            len_P1_.extend(len_)
            sur_P1_.extend(sur_)
        elif ö == 3:
            krs_P2.extend(krs_)
            krshoot_P2.extend(krshoot_)
            krplant_P2.extend(krplant_)
            jc_P2.extend(jc_)
            eswp_P2.extend(eswp_)
            eawp_P2.extend(eawp_)
            cwp_P2.extend(cwp_)
            vol_P2_.extend(vol_)
            len_P2_.extend(len_)
            sur_P2_.extend(sur_)
        elif ö == 4:
            krs_P3.extend(krs_)
            krshoot_P3.extend(krshoot_)
            krplant_P3.extend(krplant_)
            jc_P3.extend(jc_)
            eswp_P3.extend(eswp_)
            eawp_P3.extend(eawp_)
            cwp_P3.extend(cwp_)
            vol_P3_.extend(vol_)
            len_P3_.extend(len_)
            sur_P3_.extend(sur_)

def getFigdata(pdFinal, name, legend, ti):
    krs_final_mean = pdFinal.groupby('day').mean()
    krs_final_std = pdFinal.groupby('day').std()
    ti = ti
    
    fig = plt.figure()
    plt.plot(ti, np.array(krs_final_mean.krs_P0),label = 'P0' )
    plt.fill_between(ti, krs_final_mean.krs_P0 - krs_final_std.krs_P0 ,krs_final_mean.krs_P0 + krs_final_std.krs_P0, alpha = 0.1)
    plt.plot(ti, np.array(krs_final_mean.krs_P1),label = 'P1' )
    plt.fill_between(ti, krs_final_mean.krs_P1 - krs_final_std.krs_P1 ,krs_final_mean.krs_P1 + krs_final_std.krs_P1, alpha = 0.1)
    plt.plot(ti, np.array(krs_final_mean.krs_P2),label = 'P2')
    plt.fill_between(ti, krs_final_mean.krs_P2 - krs_final_std.krs_P2 ,krs_final_mean.krs_P2 + krs_final_std.krs_P2, alpha = 0.1)
    plt.plot(ti, np.array(krs_final_mean.krs_P3),label = 'P3')
    plt.fill_between(ti, krs_final_mean.krs_P3 - krs_final_std.krs_P3 ,krs_final_mean.krs_P3 + krs_final_std.krs_P3, alpha = 0.1)
    plt.ylabel(legend)
    plt.xlabel('time (d)')
    plt.xlim(0,29)
    plt.legend()
    plt.show()
    fig.savefig(name+'.png', dpi=fig.dpi)
   
krs_final = pd.DataFrame({'day': count, 'krs_P0': krs_P0,'krs_P1': krs_P1,'krs_P2': krs_P2,'krs_P3': krs_P3})   
vol_final = pd.DataFrame({'day':count, 'vol_P0': vol_P0_,'vol_P1': vol_P1_,'vol_P2': vol_P2_,'vol_P3': vol_P3_})
len_final = pd.DataFrame({'day':count, 'len_P0': len_P0_,'len_P1': len_P1_,'len_P2': len_P2_,'len_P3': len_P3_})
sur_final = pd.DataFrame({'day':count, 'sur_P0': sur_P0_,'sur_P1': sur_P1_,'sur_P2': sur_P2_,'sur_P3': sur_P3_})

krs_final.to_csv('krs_final.csv')
vol_final.to_csv('vol_final.csv')
len_final.to_csv('len_final.csv')
sur_final.to_csv('sur_final.csv')

krs_final_m = krs_final.groupby('day').mean()
krs_final_sd = krs_final.groupby('day').std()
vol_final_m = vol_final.groupby('day').mean()
vol_final_m.to_csv('VOL_final.csv', index=False)

vol_final_sd = vol_final.groupby('day').std()
len_final_m = len_final.groupby('day').mean()
len_final_sd = len_final.groupby('day').std()
sur_final_m = sur_final.groupby('day').mean()
sur_final_sd = sur_final.groupby('day').std()

P0_krs_vol = krs_final_m.krs_P0/vol_final_m.vol_P0
P1_krs_vol = krs_final_m.krs_P1/vol_final_m.vol_P1
P2_krs_vol = krs_final_m.krs_P2/vol_final_m.vol_P2
P3_krs_vol = krs_final_m.krs_P3/vol_final_m.vol_P3

P0_krs_len = krs_final_m.krs_P0/len_final_m.len_P0
P1_krs_len = krs_final_m.krs_P1/len_final_m.len_P1
P2_krs_len = krs_final_m.krs_P2/len_final_m.len_P2
P3_krs_len = krs_final_m.krs_P3/len_final_m.len_P3

P0_krs_vol_sd = krs_final_sd.krs_P0/vol_final_sd.vol_P0
P1_krs_vol_sd = krs_final_sd.krs_P1/vol_final_sd.vol_P1
P2_krs_vol_sd = krs_final_sd.krs_P2/vol_final_sd.vol_P2
P3_krs_vol_sd = krs_final_sd.krs_P3/vol_final_sd.vol_P3

fig = plt.figure()
plt.plot(np.array(vol_final_m.index), np.array(vol_final_m.vol_P0), label = 'P0')
plt.fill_between(np.array(vol_final_m.index), np.array(vol_final_m.vol_P0) - np.array(vol_final_sd.vol_P0) ,np.array(vol_final_m.vol_P0) + np.array(vol_final_sd.vol_P0), alpha = 0.1)
plt.plot(np.array(vol_final_m.index), np.array(vol_final_m.vol_P1), label = 'P1')
plt.fill_between(np.array(vol_final_m.index), np.array(vol_final_m.vol_P1) - np.array(vol_final_sd.vol_P1) ,np.array(vol_final_m.vol_P1) + np.array(vol_final_sd.vol_P1), alpha = 0.1)
plt.plot(np.array(vol_final_m.index), np.array(vol_final_m.vol_P2), label = 'P2')
plt.fill_between(np.array(vol_final_m.index), np.array(vol_final_m.vol_P2) - np.array(vol_final_sd.vol_P2) ,np.array(vol_final_m.vol_P2) + np.array(vol_final_sd.vol_P2), alpha = 0.1)
plt.plot(np.array(vol_final_m.index), np.array(vol_final_m.vol_P3), label = 'P3')
plt.fill_between(np.array(vol_final_m.index), np.array(vol_final_m.vol_P3) - np.array(vol_final_sd.vol_P3) ,np.array(vol_final_m.vol_P3) + np.array(vol_final_sd.vol_P3), alpha = 0.1)
plt.legend()

plt.savefig('vol_final_m' +'.png',dpi=fig.dpi)

fig = plt.figure()
plt.plot(np.array(len_final_m.index), np.array(len_final_m.len_P0), label = 'P0')
plt.fill_between(np.array(len_final_m.index), np.array(len_final_m.len_P0) - np.array(len_final_sd.len_P0) ,np.array(len_final_m.len_P0) + np.array(len_final_sd.len_P0), alpha = 0.1)
plt.plot(np.array(len_final_m.index), np.array(len_final_m.len_P1), label = 'P1')
plt.fill_between(np.array(len_final_m.index), np.array(len_final_m.len_P1) - np.array(len_final_sd.len_P1) ,np.array(len_final_m.len_P1) + np.array(len_final_sd.len_P1), alpha = 0.1)
plt.plot(np.array(len_final_m.index), np.array(len_final_m.len_P2), label = 'P2')
plt.fill_between(np.array(len_final_m.index), np.array(len_final_m.len_P2) - np.array(len_final_sd.len_P2) ,np.array(len_final_m.len_P2) + np.array(len_final_sd.len_P2), alpha = 0.1)
plt.plot(np.array(len_final_m.index), np.array(len_final_m.len_P3), label = 'P3')
plt.fill_between(np.array(len_final_m.index), np.array(len_final_m.len_P3) - np.array(len_final_sd.len_P3) ,np.array(len_final_m.len_P3) + np.array(len_final_sd.len_P3), alpha = 0.1)
plt.legend()

plt.savefig('len_final_m' +'.png',dpi=fig.dpi)

fig = plt.figure()
plt.plot(np.array(sur_final_m.index), np.array(sur_final_m.sur_P0), label = 'P0')
plt.fill_between(np.array(sur_final_m.index), np.array(sur_final_m.sur_P0) - np.array(sur_final_sd.sur_P0) ,np.array(sur_final_m.sur_P0) + np.array(sur_final_sd.sur_P0), alpha = 0.1)
plt.plot(np.array(sur_final_m.index), np.array(sur_final_m.sur_P1), label = 'P1')
plt.fill_between(np.array(sur_final_m.index), np.array(sur_final_m.sur_P1) - np.array(sur_final_sd.sur_P1) ,np.array(sur_final_m.sur_P1) + np.array(sur_final_sd.sur_P1), alpha = 0.1)
plt.plot(np.array(sur_final_m.index), np.array(sur_final_m.sur_P2), label = 'P2')
plt.fill_between(np.array(sur_final_m.index), np.array(sur_final_m.sur_P2) - np.array(sur_final_sd.sur_P2) ,np.array(sur_final_m.sur_P2) + np.array(sur_final_sd.sur_P2), alpha = 0.1)
plt.plot(np.array(sur_final_m.index), np.array(sur_final_m.sur_P3), label = 'P3')
plt.fill_between(np.array(sur_final_m.index), np.array(sur_final_m.sur_P3) - np.array(sur_final_sd.sur_P3) ,np.array(sur_final_m.sur_P3) + np.array(sur_final_sd.sur_P3), alpha = 0.1)
plt.legend()

plt.savefig('sur_final_m' +'.png',dpi=fig.dpi)

krs_final_m.to_csv('KRS_results.csv', index=False)

getFigdata(krs_final,"Krs", ' Krs (cm$^2$ d$^{-1}$ )', ti)

