""" water movement within the root (static soil) """
#directoryN = "/7to14dry/"

import sys; 
CPBdir = "../.."
sys.path.append(CPBdir+"/src");
sys.path.append(CPBdir);
sys.path.append("../../..");sys.path.append(".."); 
sys.path.append(CPBdir+"/src/python_modules");
sys.path.append("../build-cmake/cpp/python_binding/") # dumux python binding
sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../modules/") # python wrappers 

#import matplotlib
#matplotlib.use('AGG') #otherwise "import matplotlib.pyplot" hangs

#from rosi_richards import RichardsSP  # C++ part (Dumux binding)
#from richards import RichardsWrapper  # Python part
from functional.phloem_flux import PhloemFluxPython  # Python hybrid solver
#from Leuning_speedup import Leuning #about 0.7 for both
#from photosynthesis_cython import Leuning
import plantbox as pb
#from plantbox import Photosynthesis as ph
#import vtk_plot as vp
import math
import os
import numpy as np
#import vtk_plot as vp
#import matplotlib.pyplot as plt
from datetime import datetime, timedelta
#from kr_kx import *  



isCluster = (os.environ['HOME'] == '/home/m.giraud')
#if isCluster:
 #   def print(*args, **kwargs):
  #      """ custom print() function.
   #         for cluster: can get output even if program stop
    #        unexpectedly (i.e., without creating the outputfile)
     #   """
      #  # Adding new arguments to the print function signature
      #  # is probably a bad idea.
   #     # Instead consider testing if custom argument keywords
   #     # are present in kwargs
   #     if 'sep' in kwargs:
   #         sep = kwargs['sep']
   #     else:
   #         sep = ' '
   #     home_dir = os.getcwd()
   #     dir_name =  "/results"+directoryN
   #     dir_name2 = home_dir + dir_name
   #     name2 = dir_name2 + 'prints.txt'
   #     with open(name2, 'a') as log:
   #         for arg in args: log.write(str(arg) + sep)
   #         log.write('\n')

#       qr, qs, alpha, n, ks (in [cm/d])
#vgSoil = [0.059, 0.45, 0.00644, 1.503, 1]
def theta2H(vg,theta):#(-) to cm
    thetar =vg[0]# 0.059
    thetas = vg[1]#0.445
    alpha = vg[2]#0.00644
    n = vg[3]#1.503
    nrev = 1/(1-1/n)
    H =-(((( (thetas - thetar)/(theta - thetar))**nrev) - 1)**(1/n))/alpha
    return(H)#cm

def sinusoidal(t):
    return (np.sin(np.pi*t*2)+1)/2 #( (t%1) < 0.5)#
#https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html
#https://cds.climate.copernicus.eu/cdsapp#!/dataset/projections-cmip6?tab=form
#https://cds.climate.copernicus.eu/toolbox-editor/141198/getdata
#https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
def qair2rh(qair, press,TK):
    #pressPa = press *100
    #eaPa = qair* pressPa/(0.622+0.378*qair)
    #ea = eaPa/100
    T0 = 273.16
    RH = 26.3 * press * qair  /(exp((17.67*(TK - T0))/(TK- 29.65)))
    return RH



def div0(a, b, c):        
    return np.divide(a, b, out=np.full(len(a), c), where=b!=0)
    
def div0f(a, b, c):    
    if b != 0:
        return a/b
    else:
        return a/c
        



def setKrKx_xylem(TairC, RH,r,kr_l): #inC
    #mg/cm3
    hPa2cm = 1.0197
    dEauPure = (999.83952 + TairC * (16.952577 + TairC * 
        (- 0.0079905127 + TairC * (- 0.000046241757 + TairC * 
        (0.00000010584601 + TairC * (- 0.00000000028103006)))))) /  (1 + 0.016887236 * TairC)
    siPhi = (30 - TairC) / (91 + TairC)
    siEnne=0
    mu =  pow(10, (- 0.114 + (siPhi * (1.1 + 43.1 * pow(siEnne, 1.25) )))) 
    mu = mu /(24*60*60)/100/1000; #//mPa s to hPa d, 1.11837e-10 hPa d for pure water at 293.15K
    mu = mu * hPa2cm #hPa d to cmh2o d 

    #number of vascular bundles
    VascBundle_leaf = 32
    VascBundle_stem = 52
    VascBundle_root = 1 #valid for all root type
            
    #radius of xylem type^4 * number per bundle
    rad_x_l_1   = (0.0015 **4) * 2; rad_x_l_2   = (0.0005 **4) * 2   
    rad_x_s_1   = (0.0017 **4) * 3; rad_x_s_2   = (0.0008 **4) * 1     
    rad_x_r0_1  = (0.0015 **4) * 4    
    rad_x_r12_1 = (0.00041**4) * 4; rad_x_r12_2 = (0.00087**4) * 1
    rad_x_r3_1  = (0.00068**4) * 1      

    # axial conductivity [cm^3/day]  
    betaXylX =1# 0.1      
    kz_l  = VascBundle_leaf *(rad_x_l_1 + rad_x_l_2)    *np.pi /(mu * 8)  * betaXylX
    kz_s  = VascBundle_stem *(rad_x_s_1 + rad_x_s_2)    *np.pi /(mu * 8) * betaXylX 
    kz_r0 = VascBundle_root * rad_x_r0_1                *np.pi /(mu * 8) * betaXylX  
    kz_r1 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*np.pi /(mu * 8)  * betaXylX
    kz_r2 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*np.pi /(mu * 8)  * betaXylX 
    kz_r3 = VascBundle_root * rad_x_r3_1                *np.pi /(mu * 8) * betaXylX

    #radial conductivity [1/day],0.00014 #
    betaXyl = 1#0.1#0.1
    #kr_l  = 3.83e-5 * hPa2cm * betaXyl# init: 3.83e-4 cm/d/hPa
    kr_s  = 0.#1.e-20  * hPa2cm # set to almost 0
    kr_r0 =6.37e-5 * hPa2cm * betaXyl
    kr_r1 =7.9e-5  * hPa2cm * betaXyl
    kr_r2 =7.9e-5  * hPa2cm * betaXyl
    kr_r3 =6.8e-5  * hPa2cm * betaXyl
    l_kr = 0.8 #cm
    # r.setKr([[kr_r0],[kr_s],[kr_l]]) 
    #r.setKr_meso([kr_l]) 
    # r.setKx([[kz_r0],[kz_s],[kz_l]])
    
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
        # kr of tap root at age 0

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
    #too high
    #m2/d
    kr_l = np.array([[0, 3.83e-6 * hPa2cm], [1e4, 1000 * hPa2cm]])
    kz_l = np.array([[0, kx_stem/10], [1e4, kx_stem/10]])

    kr0 = kr_l
    kz0 = kz_s
    '''Artificial Shoot kx'''

    # artificial shoot
    kr0 = np.array([[0, 1.e-12], [1e4, 1.e-12]])
    kz0 = np.array([[0, 1.], [1e4, 1.]])	

    # kz0 = np.array([[0., 0.33355685296285653], [1.e-20, 0.33355685296285653]])		#
    # kz0 = np.array([[0, 0.356832], [1e20, 0.356832]])		# kx of tap root at age 1e20

    # tap root
    # kr0 = np.array([[0, 1.14048e-3], [2, 1.08864e-3], [4, 1.0368e-3], [6, 9.8486e-4], [8, 9.3312e-4], [10, 8.8992e-4], [12, 8.47584e-4], [14, 8.06112e-4], [16, 7.67232e-4], [18, 7.3008e-4], [20, 6.9552e-4], [22, 6.61824e-4], [24, 6.29856e-4], [26, 5.99616e-4], [28, 5.7024e-4], [30, 5.42592e-4], [32, 5.16672e-4], [1e20, 5.16672e-4]])
    # kz0 = np.array([[0, 0.067392], [2, 0.074736], [4, 0.082944], [6, 0.092448], [8, 0.101952], [10, 0.113184], [12, 0.126144], [14, 0.139968], [16, 0.154656], [18, 0.171936], [20, 0.190944], [22, 0.21168], [24, 0.235008], [26, 0.260928], [28, 0.28944], [30, 0.321408], [32, 0.356832], [1e20, 0.356832]])


    '''Root conductivities with exponential growth'''
    kr1 = np.array([[0.08243227637739019,0.0001063427809106925],[0.12403083800012145,0.00010600819000027001],[0.1658885823638555,0.0001062360127185785],[0.20800875944065764,0.000106033016896955],[0.2503946807170684,0.000106346827879547],[0.2930497207564021,0.00010642950047885449],[0.33597731881095994,0.0001062634438556405],[0.5104844350538453,0.00010638469493007],[0.6896456188450584,0.0001064569258243585],[0.8737159341936453,0.0001057606079808335],[1.2577145170642563,0.000113565939216522],[1.4582710891214414,0.0001142901954040455],[1.6649995920099006,4.00076833588706e-05],[1.8782919901321284,3.994978345837505e-05],[2.098578814953127,3.99623454358767e-05],[2.3263343967827597,3.9885567944048e-05],[2.562083014935793,4.0027559345494755e-05],[2.8064061667829177,3.676500658534805e-05],[3.0599512092164862,3.6763439236325344e-05],[4.494773247525064,3.890549566596665e-05],[8.838909102521768,4.10210723656745e-05],[24.94832020698389,4.29179747391326e-05]])
    kz1 = np.array([[0.08243227637739019,0.001428701491936845],[0.12403083800012145,0.00141157523117865],[0.1658885823638555,0.00150353452632553],[0.20800875944065764,0.0014742021745241351],[0.2503946807170684,0.00137229824465936],[0.2930497207564021,0.001462688585623995],[0.33597731881095994,0.001406009043608815],[0.5104844350538453,0.0014675452059425699],[0.6896456188450584,0.0014077415905594049],[0.8737159341936453,0.0015130204053065551],[1.2577145170642563,0.00133612567641224],[1.4582710891214414,0.001350817954244865],[1.6649995920099006,0.0013481762863067248],[1.8782919901321284,0.001958719658258785],[2.3263343967827597,0.16138896628702398],[2.562083014935793,0.1633108719691035],[2.8064061667829177,0.1717536169909495],[3.0599512092164862,0.1736277194342965],[4.494773247525064,0.1617984411846205],[8.838909102521768,0.243831044057843],[24.94832020698389,0.33355685296285653]])

    #seminal
    kr4 =np.array([[0.112816849,0.000112594],[0.169748659,0.000112978],[0.227035186,0.000113248],[0.284680878,0.000112912],[0.342690268,0.000113171],[0.401067975,0.000113439],[0.459818704,0.000113637],[0.698649219,0.000121677],[0.943849293,0.00012169],[1.195768006,0.000122518],[1.454783947,0.00012298],[1.721308633,0.000122788],[1.995790444,0.000123039],[2.278719163,0.000123636],[2.570631231,0.0001233],[2.872115875,4.21906e-05],[3.183822263,4.2315e-05],[3.506467924,4.23972e-05],[3.840848695,4.24768e-05],[4.18785055,4.25832e-05],[6.151548613,3.92503e-05],[12.09693483,4.16085e-05],[34.14428184,4.30515e-05]])
    kz4 = np.array([[0.11281684914998105,0.00104344808387714],[0.16974865860241395,0.001084251914759505],[0.2270351856664305,0.0010525523998907351],[0.2846808782552067,0.001103841878199095],[0.34269026847065615,0.00108312954572133],[0.40106797474159234,0.0010675539603493649],[0.4598187040302041,0.0010069550887520464],[0.6986492188960015,0.0009735756104467685],[0.9438492926240302,0.0009768755278932784],[1.1957680059275344,0.0009854992034376783],[1.4547839467929418,0.0009071226233668881],[1.7213086327468925,0.000865667148399697],[1.9957904441216752,0.0008010248002760815],[2.278719162705088,0.0008186655439126395],[2.570631231148137,0.000870364391740444],[2.8721158750002798,0.000843754941907991],[3.183822262929019,0.0008332761382412236],[3.5064679238316003,0.000911009348285576],[3.840848695261385,0.000875389202048173],[4.187850550141607,0.000898246121518689],[6.151548613165449,0.11215436059298149],[12.096934825678332,0.10148217370333051],[34.144281839931345,0.150518047844499]])

    # l-type lateral
    # kr2 = np.array([[0, 4.11264e-3], [1, 3.888e-3], [2, 3.672e-3], [3, 3.47328e-3], [4, 3.2832e-3], [5, 3.10176e-3], [6, 2.92896e-3], [7, 2.77344e-3], [8, 2.61792e-3], [9, 2.47968e-3], [10, 2.34144e-3], [11, 2.21184e-3], [12, 2.09952e-3], [13, 1.97856e-3], [14, 1.86624e-3], [15, 1.76256e-3], [16, 1.66752e-3], [17, 1.58112e-3], [1e20, 1.58112e-3]])
    # kz2 = np.array([[0, 4.06944e-4], [1, 5.00256e-4], [2, 6.15168e-4], [3, 7.56e-4], [4, 9.3312e-4], [5, 1.14048e-3], [6, 1.40832e-3], [7, 1.728e-3], [8, 2.12544e-3], [9, 2.60928e-3], [10, 3.21408e-3], [11, 3.94848e-3], [12, 4.85568e-3], [13, 5.97024e-3], [14, 7.344e-3], [15, 8.9856e-3], [16, 0.0110592], [17, 0.0136512], [1e20, 0.0136512]])

    kr2 = np.array([[0.17216300608705393,0.00011814790856374001],[0.26033859445570695,0.00011772680594141399],[0.34997169639406706,0.00011759868155066849],[0.44111130735821585,0.000118013417528984],[0.5338089356304391,0.00011772378367250949],[0.6281187771458294,0.0001178036806349285],[0.7240979057925947,0.000117636862682687],[1.125961665200367,0.000117637966391041],[1.560007591489408,0.00011792394054094201],[2.0318423759545428,3.7481355372504104e-05],[2.548678905856449,3.748182033219565e-05],[3.120015182978816,3.437261303083105e-05],[3.75871480095179,3.433395216767585e-05],[4.482812706744385,3.43365166096648e-05],[5.318722392441198,3.439823298059135e-05],[6.30739370680824,3.4359288049365e-05],[7.51742960190358,3.44139257267471e-05],[9.077437193392988,3.43931832083891e-05],[11.276144402855369,3.4370715979269455e-05],[15.03485920380716,3.4363613145709196e-05]])
    kz2 = np.array([[0.17216300608705393,4.4768072266751145e-05],[0.26033859445570695,4.6349054075485354e-05],[0.34997169639406706,3.6227416994815655e-05],[0.44111130735821585,4.8713327058618e-05],[0.5338089356304391,3.9706371946602345e-05],[0.6281187771458294,4.2895381334739246e-05],[0.7240979057925947,4.65580350149153e-05],[1.125961665200367,4.22104082674339e-05],[1.560007591489408,3.8594389322433053e-05],[2.0318423759545428,4.2146643130030445e-05],[3.120015182978816,0.004316757905359855],[3.75871480095179,0.00431239331466016],[4.482812706744385,0.004323225988052145],[5.318722392441198,0.00430169788625176],[6.30739370680824,0.004345754295000255],[7.51742960190358,0.00434487892228969],[9.077437193392988,0.0043265316331582305],[11.276144402855369,0.004304888540330575],[15.03485920380716,0.0042966613356701655]])

    # s-type lateral
    # kr3 = np.array([[0, 4.11264e-3], [1, 3.888e-3], [2, 3.672e-3], [3, 3.47328e-3], [4, 3.2832e-3], [5, 3.10176e-3], [6, 2.92896e-3], [7, 2.77344e-3], [8, 2.61792e-3], [9, 2.47968e-3], [10, 2.34144e-3], [11, 2.21184e-3], [12, 2.09952e-3], [13, 1.97856e-3], [14, 1.86624e-3], [15, 1.76256e-3], [16, 1.66752e-3], [17, 1.58112e-3], [1e20, 1.58112e-3]])
    # kz3 = np.array([[0, 4.06944e-4], [1, 5.00256e-4], [2, 6.15168e-4], [3, 7.56e-4], [4, 9.3312e-4], [5, 1.14048e-3], [6, 1.40832e-3], [7, 1.728e-3], [8, 2.12544e-3], [9, 2.60928e-3], [10, 3.21408e-3], [11, 3.94848e-3], [12, 4.85568e-3], [13, 5.97024e-3], [14, 7.344e-3], [15, 8.9856e-3], [16, 0.0110592], [17, 0.0136512], [1e20, 0.0136512]])

    kr3 = np.array([[0.200254122123127,0.000129390622450947],[0.3039865472077396,0.000129699289462529],[0.41028052389123926,0.000129232506285035],[0.51926577126878,0.000129567186814999],[0.631082118606905,0.0001292076875711025],[0.7458805840656142,0.0001295009114293915],[0.8638246012425924,0.0001296203509850465],[1.3708437430993439,0.0001073382759855375],[1.9456546377424146,3.365280024629185e-05],[2.6092251168618796,3.379687849238945e-05],[3.3940624132730615,3.10777016348636e-05],[4.3546256324811665,3.110714184120505e-05],[5.5930070062437025,2.70109885520665e-05],[7.338407521862991,2.694561559503295e-05],[10.32218941124481,2.70282984028366e-05]])
    kz3 = np.array([[0.200254122123127,1.4876873204742652e-05],[0.3039865472077396,1.60334194419736e-05],[0.41028052389123926,1.5787513976139952e-05],[0.8638246012425924,0.000776201512075445],[1.3708437430993439,0.0007679983531071785],[1.9456546377424146,0.0007745252668304815],[2.6092251168618796,0.0007652181272701905],[3.3940624132730615,0.0007729805938091705],[4.3546256324811665,0.0007682011183304045],[5.5930070062437025,0.0007532090253400926],[7.338407521862991,0.0007783360746834],[10.32218941124481,0.0007613256253168725]])

    # shoot-born
    kr5 = np.array([[0.08233900727871345,0.000133675979105574],[0.1238905016540107,0.000133826931499557],[0.16570088551454124,0.00013353191679336398],[0.207773405155146,0.000133533966573953],[0.25011136831553993,0.000134758249008305],[0.2927181457408428,0.0001344153385482685],[0.33559717279196777,0.00013528554438529348],[0.5099068406304184,0.0001350674852631835],[0.6888653101103879,0.0001405508737968255],[0.8727273566453477,0.00014100645240813547],[1.0617692914353087,0.0001421701293576315],[1.2562914591971874,0.00014259222363723702],[1.456621108845681,0.0001472611740842825],[1.663115706011982,0.00014806700032824903],[1.8761667716052333,0.00014844271074813552],[2.0962043499598675,0.0001395596875498485],[2.323702234698392,0.000144218683466018],[2.5591841119327503,5.28436517340678e-05],[2.803230821090675,5.2815110144092294e-05],[3.056488986603928,5.30836636709535e-05],[4.489687576378076,5.08139842681165e-05],[8.828908200914931,5.23839696968556e-05],[24.920092097298415,5.6739924599214845e-05]])
    kz5 = np.array([[0.08233900727871345,0.0061012108489030895],[0.1238905016540107,0.0061287157563131],[0.16570088551454124,0.0061382740130904445],[0.207773405155146,0.00603781321781165],[0.25011136831553993,0.00631172718713725],[0.2927181457408428,0.005902104241609545],[0.33559717279196777,0.00596042829982248],[0.5099068406304184,0.0057859354283172106],[0.6888653101103879,0.00577029495359218],[0.8727273566453477,0.00669714417581805],[1.0617692914353087,0.007411225692761311],[1.2562914591971874,0.0069813721132007],[1.456621108845681,0.0063883507696567456],[1.663115706011982,0.006802020444944466],[1.8761667716052333,0.006945851992287625],[2.323702234698392,1.427732171171995],[2.5591841119327503,1.51687143593099],[2.803230821090675,1.60722277926901],[3.056488986603928,1.69169163420242],[4.489687576378076,2.34489158735256],[8.828908200914931,4.09461156681657],[24.920092097298415,6.608610426314295]])
    # kr5 = kr1  
    # kz5 = kz1

    '''Root conductivities with linear growth from Adrien'''

    kr1 = np.array([[0.227272727,0.000106343],[0.340909091,0.000106008],[0.454545455,0.000106236],[0.568181818,0.000106033],[0.681818182,0.000106347],[0.795454545,0.00010643],[0.909090909,0.000106263],[1.363636364,0.000106385],[1.818181818,0.000106457],[2.272727273,0.000105761],[3.181818182,0.000113566],[3.636363636,0.00011429],[4.090909091,4.00E-05],[4.545454545,3.99E-05],[5.,	4.00E-05],[5.454545455,3.99E-05],[5.909090909,4.00E-05],[6.363636364,3.68E-05],[6.818181818,3.68E-05],[9.090909091,3.89E-05],[13.63636364,4.10E-05],[18.18181818,4.29E-05]])

    kr2= np.array([[0.3125,0.000112594],[0.46875,0.000112978],[0.625,0.000113248],[0.78125,0.000112912],[0.9375,0.000113171],[1.09375,0.000113439],[1.25,0.000113637],[1.875,0.000121677],[2.5,0.00012169],[3.125,	0.000122518],[3.75,0.00012298],[4.375,0.000122788],[5,0.000123039],[5.625,0.000123636],[6.25,0.0001233],[6.875,4.22E-05],[7.5,4.23E-05],[8.125,4.24E-05],[8.75,4.25E-05],[9.375,4.26E-05],[12.5,3.93E-05],[18.75,4.16E-05],[25,4.31E-05]])

    kr3 = np.array([[2,0.000118148],[3,0.000117727],[4,0.000117599],[5,0.000118013],[6,0.000117724],[7,0.000117804],[8,0.000117637],[12,0.000117638],[16,0.000117924],[20,3.75E-05],[24,3.75E-05],[28,3.44E-05],[32,3.43E-05],[36,3.43E-05],[40,3.44E-05],[44,3.44E-05],[48,3.44E-05],[52,3.44E-05],[56,3.44E-05],[60,3.44E-05]])

    kr4 = np.array([[3.649635036,0.000129391],[5.474452555,0.000129699],[7.299270073,0.000129233],[9.124087591,0.000129567],[10.94890511,0.000129208],[12.77372263,0.000129501],[14.59854015,0.00012962],[21.89781022,0.000107338],[29.19708029,3.37E-05],[36.49635036,3.38E-05],[43.79562044,3.11E-05],[51.09489051,3.11E-05],[58.39416058,2.70E-05],[65.69343066,2.69E-05],[72.99270073,2.70E-05]])

    kr5 = np.array([[0.185185185,0.000133676],[0.277777778,0.000133827],[0.37037037,0.000133532],[0.462962963,0.000133534],[0.555555556,0.000134758],[0.648148148,0.000134415],[0.740740741,0.000135286],[1.111111111,0.000135067],[1.481481481,0.000140551],[1.851851852,0.000141006],[2.222222222,0.00014217],[2.592592593,0.000142592],[2.962962963,.000147261],[3.333333333,0.000148067],[3.703703704,0.000148443],[4.074074074,0.00013956],[4.444444444,0.000144219],[4.814814815,5.28E-05],[5.185185185,5.28E-05],[5.555555556,5.31E-05],[7.407407407,5.08E-05],[11.11111111,5.24E-05],[14.81481481,5.67E-05]])


    kz1 = np.array([[0.227272727,0.001428701],[0.340909091,0.001411575],[0.454545455,0.001503535],[0.568181818,0.001474202],[0.681818182,0.001372298],[0.795454545,0.001462689],[0.909090909,0.001406009],[1.363636364,0.001467545],[1.818181818,0.001407742],[2.272727273,0.00151302],[3.181818182,0.001336126],[3.636363636,0.001350818],[4.090909091,0.001348176],[4.545454545,0.00195872],[5.454545455,0.161388966],[5.909090909,0.163310872],[6.363636364,0.171753617],[6.818181818,0.173627719],[9.090909091,0.161798441],[13.63636364,0.243831044],[18.18181818,0.333556853]])

    kz2 = np.array([[0.3125,0.001043448],[0.46875,0.001084252],[0.625,0.001052552],[0.78125,0.001103842],[0.9375,0.00108313],[1.09375,0.001067554],[1.25,0.001006955],[1.875,0.000973576],[2.5,0.000976876],[3.125,0.000985499],[3.75,0.000907123],[4.375,0.000865667],[5,0.000801025],[5.625,0.000818666],[6.25,0.000870364],[6.875,0.000843755],[7.5,0.000833276],[8.125,0.000911009],[8.75,0.000875389],[9.375,0.000898246],[12.5,0.112154361],[18.75,0.101482174],[25,0.150518048]])

    kz3 = np.array([[2,4.48E-05],[3,4.63E-05],[4,3.62E-05],[5,4.87E-05],[6,3.97E-05],[7,4.29E-05],[8,4.66E-05],[12,4.22E-05],[16,3.86E-05],[20,4.21E-05],[28,0.004316758],[32,0.004312393],[36,0.004323226],[40,0.004301698],[44,0.004345754],[48,0.004344879],[52,0.004326532],[56,0.004304889],[60,0.004296661]])

    kz4 = np.array([[3.649635036,1.49E-05],[5.474452555,1.60E-05],[7.299270073,1.58E-05],[14.59854015,0.000776202],[21.89781022,0.000767998],[29.19708029,0.000774525],[36.49635036,0.000765218],[43.79562044,0.000772981],[51.09489051,0.000768201],[58.39416058,0.000753209],[65.69343066,0.000778336],[72.99270073,0.000761326]])

    kz5 = np.array([[0.185185185,0.006101211],[0.277777778,0.006128716],[0.37037037,0.006138274],[0.462962963,0.006037813],[0.555555556,0.006311727],[0.648148148,0.005902104],[0.740740741,0.005960428],[1.111111111,0.005785935],[1.481481481,0.005770295],[1.851851852,0.006697144],[2.222222222,0.007411226],[2.592592593,0.006981372],[2.962962963,0.006388351],[3.333333333,0.00680202],[3.703703704,0.006945852],[4.444444444,1.427732171],[4.814814815,1.516871436],[5.185185185,1.607222779],[5.555555556,1.691691634],[7.407407407,2.344891587],[11.11111111,4.094611567],[14.81481481,6.608610426]])




    r.setKrTables([[ kr1[:, 1], kr2[:, 1], kr3[:, 1], kr4[:, 1], kr5[:, 1]],[kr_s[:, 1],kr_s[:, 1]],[kr_l[:, 1]]], [[kr1[:, 0], kr2[:, 0], kr3[:, 0], kr4[:, 0], kr5[:, 0]],[kr_s[:, 0],kr_s[:, 0]],[kr_l[:, 0]]])
    r.setKxTables([[kz1[:, 1], kz2[:, 1], kz3[:, 1], kz4[:, 1], kz5[:, 1]],[kz_s[:, 1],kz_s[:, 1]],[kz_l[:, 1]]],[[kz1[:, 0], kz2[:, 0], kz3[:, 0], kz4[:, 0], kz5[:, 0]],[kz_s[:, 0],kz_s[:, 0]],[kz_l[:, 0]]])
    
    Rgaz=8.314 #J K-1 mol-1 = cm^3*MPa/K/mol
    rho_h2o = dEauPure/1000#g/cm3
    Mh2o = 18.05 #g/mol
    MPa2hPa = 10000
    hPa2cm = 1/0.9806806
    #log(-) * (cm^3*MPa/K/mol) * (K) *(g/cm3)/ (g/mol) * (hPa/MPa) * (cm/hPa) =  cm                      
    #p_a = np.log(RH) * Rgaz * rho_h2o * (TairC + 273.15)/Mh2o * MPa2hPa * hPa2cm
    #done withint solve photosynthesis
    #r.psi_air = p_a #*MPa2hPa #used only with xylem
    return r

    
def setKrKx_phloem(r): #inC

    #number of vascular bundles
    VascBundle_leaf = 32
    VascBundle_stem = 52
    VascBundle_root = 1 #valid for all root type
            
    #numPerBundle
    numL = 18
    numS = 21
    numr0 = 33
    numr1 = 25
    numr2 = 25
    numr3 = 1
        
    #radius of phloem type^4 * number per bundle
    rad_s_l   = numL* (0.00025 **4)# * 2; rad_x_l_2   = (0.0005 **4) * 2   
    rad_s_s   = numS *(0.00019 **4) #* 3; rad_x_s_2   = (0.0008 **4) * 1     
    rad_s_r0  = numr0 *(0.00039 **4) #* 4    
    rad_s_r12 = numr1*(0.00035**4) #* 4; rad_x_r12_2 = (0.00087**4) * 1
    rad_s_r3  = numr3 *(0.00068**4) #* 1      

    # axial conductivity [cm^3/day] , mu is added later as it evolves with CST  
    beta = 0.9 #Thompson 2003a
    kz_l   = VascBundle_leaf * rad_s_l   * np.pi /8 * beta  
    kz_s   = VascBundle_stem * rad_s_s   * np.pi /8 * beta
    kz_r0  = VascBundle_root * rad_s_r0  * np.pi /8 * beta
    kz_r12 = VascBundle_root * rad_s_r12 * np.pi /8 * beta
    kz_r3  = VascBundle_root * rad_s_r3  * np.pi /8 * beta
    
    #print([[kz_r0,kz_r12,kz_r12,kz_r3],[kz_s,kz_s ],[kz_l]])
    #raise Exception
    #radial conductivity [1/day],
    kr_l  = 0.#3.83e-4 * hPa2cm# init: 3.83e-4 cm/d/hPa
    kr_s  = 0.#1.e-20  * hPa2cm # set to almost 0
    # kr_r0 = 1e-1
    # kr_r1 = 1e-1
    # kr_r2 = 1e-1
    # kr_r3 = 1e-1
    kr_r0 = 5e-2
    kr_r1 = 5e-2
    kr_r2 = 5e-2
    kr_r3 = 5e-2
    l_kr = 0.8 #0.8 #cm
    
    r.setKr_st([[kr_r0,kr_r1 ,kr_r2 ,kr_r0],[kr_s,kr_s ],[kr_l]] , kr_length_= l_kr)
    r.setKx_st([[kz_r0,kz_r12,kz_r12,kz_r0],[kz_s,kz_s ],[kz_l]])
    
    a_ST = [[0.00039,0.00035,0.00035,0.00039 ],[0.00019,0.00019],[0.00025]]
    Across_s_l   = numL*VascBundle_leaf *(a_ST[2][0]**2)*np.pi# (0.00025 **2)# * 2; rad_x_l_2   = (0.0005 **4) * 2   
    Across_s_s   = numS *VascBundle_stem * (a_ST[1][0]**2)*np.pi#(0.00019 **2) #* 3; rad_x_s_2   = (0.0008 **4) * 1     
    Across_s_r0  = numr0 *VascBundle_root * (a_ST[0][0]**2)*np.pi#(0.00039 **2) #* 4    
    Across_s_r12 = numr1*VascBundle_root * (a_ST[0][1]**2)*np.pi#(0.00035**2) #* 4; rad_x_r12_2 = (0.00087**4) * 1
    Across_s_r3  =  numr3 *VascBundle_root *(a_ST[0][2]**2)*np.pi# (0.00068**2) #* 1    
    
    Perimeter_s_l   = numL*VascBundle_leaf *(a_ST[2][0])* 2 * np.pi# (0.00025 **2)# * 2; rad_x_l_2   = (0.0005 **4) * 2   
    Perimeter_s_s   = numS *VascBundle_stem * (a_ST[1][0])* 2 * np.pi#(0.00019 **2) #* 3; rad_x_s_2   = (0.0008 **4) * 1     
    Perimeter_s_r0  = numr0 *VascBundle_root * (a_ST[0][0])* 2 * np.pi#(0.00039 **2) #* 4    
    Perimeter_s_r12 = numr1*VascBundle_root * (a_ST[0][1])* 2 * np.pi#(0.00035**2) #* 4; rad_x_r12_2 = (0.00087**4) * 1
    Perimeter_s_r3  =  numr3 *VascBundle_root *(a_ST[0][2])* 2 * np.pi# (0.00068**2) #* 1  
    #print(a_ST[2][0],a_ST[1][0],a_ST[0][0],a_ST[0][1],a_ST[0][2])
    #r.a_ST = a_ST #to check for water equilibrium assumption
    #tot surface/np.pi of sieve tube  (np.pi added after)
    #r.a_ST_eqs = [[rad_s_r0,rad_s_r12,rad_s_r12,rad_s_r0],[rad_s_s,rad_s_s],[rad_s_l]]
    r.setAcross_st([[Across_s_r0,Across_s_r12,Across_s_r12,Across_s_r0],[Across_s_s,Across_s_s],[Across_s_l]])
    return r
    
""" Parameters """
def launchUQR(directoryN,simInit,simStartSim, condition,csChoise,spellDuration):
    def write_file_array(name, data):
        name2 = 'results'+ directoryN+ name+ '.txt'
        with open(name2, 'a') as log:
            log.write(','.join([num for num in map(str, data)])  +'\n')

    def write_file_float(name, data):
        name2 = 'results' + directoryN+  name+ '.txt'
        with open(name2, 'a') as log:
            log.write(repr( data)  +'\n')
            
    def weather(simDuration, hp):
        vgSoil = [0.059, 0.45, 0.00644, 1.503, 1]
        loam = [0.08, 0.43, 0.04, 1.6, 50]
        Qnigh = 0; Qday = 960e-6 #458*2.1
        if ((condition == "wet") or (simDuration <= simStartSim) or (simDuration > simStartSim +7)):
            Tnigh = 15.8; Tday = 22
            #Tnigh = 13; Tday = 20.7
            #specificHumidity = 0.0097
            RHday = 0.60; RHnigh = 0.88
            Pair = 1010.00 #hPa
            thetaInit = 40/100
            if csChoise == "yesCS":
                cs = 390e-6#50e-6
            elif csChoise == "noCS": 
                cs = 350e-6
            else:
                print("csChoise",csChoise)
                raise Exception("csChoise not recognised")
        elif condition == "dry":
            Tnigh = 20.7; Tday = 30.27
            #Tnigh = 15.34; Tday = 23.31
            #specificHumidity = 0.0097# 0.0111
            RHday = 0.7; RHnigh =0.7#0.44; RHnigh = 0.78
            Pair = 1070.00 #hPa
            thetaInit = 28/100   
            #cs = 1230e-6
            if csChoise == "yesCS":
                cs = 1230e-6
            elif csChoise == "noCS": 
                cs = 350e-6
            else:
                print("csChoise",csChoise)
                raise Exception("csChoise not recognised")
        else:
            print("condition",condition)
            raise Exception("condition not recognised")

        coefhours = sinusoidal(simDuration)
        RH_ = RHnigh + (RHday - RHnigh) * coefhours
        TairC_ = Tnigh + (Tday - Tnigh) * coefhours
        Q_ = Qnigh + (Qday - Qnigh) * coefhours
         #co2 paartial pressure at leaf surface (mol mol-1)
        #390, 1231
        #RH = 0.5 # relative humidity
        es =  6.112 * np.exp((17.67 * TairC_)/(TairC_ + 243.5))
        ea = es*RH_#qair2ea(specificHumidity,  Pair)
        assert ea < es
        #RH = ea/es
        assert ((RH_ > 0) and(RH_ < 1))
        bl_thickness = 1/1000 #m
        diffusivity= 2.5e-5#m2/sfor 25°C
        rbl =bl_thickness/diffusivity #s/m 13
        #cs = 350e-6
        Kcanopymean = 1e-1 # m2/s
        meanCanopyL = (2/3) * hp /2
        rcanopy = meanCanopyL/Kcanopymean
        windSpeed = 2 #m/s
        zmzh = 2 #m
        karman = 0.41 #[-]
        
        rair = 1
        if hp > 0:
            rair = np.log((zmzh - (2/3)*hp)/(0.123*hp)) * np.log((zmzh - (2/3)*hp)/(0.1*hp)) / (karman*karman*windSpeed)
            #print()
            #raise Exception
            

        pmean = theta2H(vgSoil, thetaInit)

        weatherVar = {'TairC' : TairC_,'TairK' : TairC_ + 273.15,'Pair':Pair,"es":es,
                        'Qlight': Q_,'rbl':rbl,'rcanopy':rcanopy,'rair':rair,"ea":ea,
                        'cs':cs, 'RH':RH_, 'p_mean':pmean, 'vg':loam}
        print("Env variables at", round(simDuration//1),"d",round((simDuration%1)*24),"hrs :\n", weatherVar)
        return weatherVar

    weatherInit = weather(0,0)
    simDuration = simInit # [day] init simtime
    #spellDuration = 5
    simMax = 25#simStartSim+ spellDuration
    depth = 60
    dt = 1/24 #10min
    verbose = True

    # plant system 
    pl = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
    #pl2 = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
    path = CPBdir+"/modelparameter/structural/plant/"
    name = "Triticum_aestivum_adapted_2023"#"Triticum_aestivum_adapted_2021"#
    print("check file", path + name + ".xml")

    pl.readParameters(path + name + ".xml")
    #pl2.readParameters(path + name + ".xml")



    #raise Exception
    sdf = pb.SDF_PlantBox(np.inf, np.inf, depth )

    pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil
    #pl2.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil


    pl.initialize(verbose = True)#, stochastic = False)
    pl.simulate(simDuration, False)#, "outputpm15.txt")
    #pl2.initialize(verbose = True)#, stochastic = False)
    #pl2.simulate(simDuration, False)#, "outputpm15.txt")
    ot_ =np.array(pl.organTypes)
    segments_ = np.array(pl.segLength())
    print(len(segments_), sum(segments_))
    print(len(segments_[np.where(ot_ ==2)]), sum(segments_[np.where(ot_ ==2)]))
    print(len(segments_[np.where(ot_ ==3)]), sum(segments_[np.where(ot_ ==3)]))
    print(len(segments_[np.where(ot_ ==4)]), sum(segments_[np.where(ot_ ==4)]))
    print(np.array(pl.leafBladeSurface)[np.where(ot_ ==4)],
        len(pl.leafBladeSurface))
    #raise Exception
    """ Coupling to soil """



    min_b = [-3./2, -12./2, -41.]#distance between wheat plants
    max_b = [3./2, 12./2, 0.]
    rez = 0.5
    cell_number = [int(6*rez), int(24*rez), int(40*rez)]#1cm3? 
    layers = depth; soilvolume = (depth / layers) * 3 * 12
    k_soil = []
    initial = weatherInit["p_mean"]#mean matric potential [cm] pressure head
    if False:
        s = RichardsWrapper(RichardsSP())
        s.initialize()
        periodic = True
        s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
        s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
        s.setTopBC("noFlux")
        s.setBotBC("constantPressure")#("freeflow")
        s.setParameter("Newton.EnableChop", "True")
        s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
        s.setVGParameters([weatherInit['vg']])
        s.initializeProblem()
        sx = s.getSolutionHead()  # inital condition, solverbase.py
        picker = lambda x, y, z: s.pick([x, y, z])    
    else:
        p_mean = initial
        p_bot = p_mean + depth/2
        p_top = initial - depth/2
        sx = np.linspace(p_top, p_bot, depth)
        picker = lambda x,y,z : max(int(np.floor(-z)),-1) 
    sx_static_bu = sx    
    pl.setSoilGrid(picker)  # maps segment
    #pl2.setSoilGrid(picker)  # maps segment


    """ Parameters phloem and photosynthesis """
    r = PhloemFluxPython(pl,psiXylInit = min(sx),ciInit = weatherInit["cs"]*0.5) #XylemFluxPython(pl)#
    #r2 = PhloemFluxPython(#pl2,psiXylInit = min(sx),ciInit = weatherInit["cs"]*0.5) #XylemFluxPython(pl)#

    r = setKrKx_phloem(r)
    print("oldciEq ",r.oldciEq)
    r.oldciEq = True
    
    print("oldciEq2 ",r.oldciEq)
    r.Rd_ref = 0 #to avoid error (C < 0 in meso, mention this in paper)
    r.g0 = 8e-3
    r.VcmaxrefChl1 =1.28#/2
    r.VcmaxrefChl2 = 8.33#/2
    r.a1 = 0.6/0.4#0.7/0.3#0.6/0.4 #ci/(cs - ci) for ci = 0.6*cs
    r.a3 = 1.5
    r.alpha = 0.4#0.2#/2
    r.theta = 0.6#0.9#/2
    r.k_meso = 1e-3#1e-4
    r.setKrm2([[2e-5]])
    r.setKrm1([[10e-2]])#([[2.5e-2]])
    r.setRhoSucrose([[0.51],[0.65],[0.56]])#0.51
    #([[14.4,9.0,0,14.4],[5.,5.],[15.]])
    rootFact = 2
    r.setRmax_st([[2.4*rootFact,1.5*rootFact,0.6*rootFact,2.4*rootFact],[2.,2.],[8.]])#6.0#*6 for roots, *1 for stem, *24/14*1.5 for leaves
    #r.setRmax_st([[12,9.0,6.0,12],[5.,5.],[15.]])
    r.KMrm = 0.1#VERY IMPORTANT TO KEEP IT HIGH
    #r.exud_k = np.array([2.4e-4])#*10#*(1e-1)
    #r.k_gr = 1#0
    r.sameVolume_meso_st = False
    r.sameVolume_meso_seg = True
    r.withInitVal =True
    r.initValST = 0.#0.6#0.0
    r.initValMeso = 0.#0.9#0.0
    r.beta_loading = 0.6
    r.Vmaxloading = 0.05 #mmol/d, needed mean loading rate:  0.3788921068507634
    r.Mloading = 0.2
    r.Gr_Y = 0.8
    r.CSTimin = 0.4
    r.surfMeso=0.0025
    r.leafGrowthZone = 2 # cm
    r.StemGrowthPerPhytomer = True # 
    r.psi_osmo_proto = -10000*1.0197 #schopfer2006
    r.fwr = 0

    r.cs = weatherInit["cs"]

    #r.r_forPhloem(24/14*1.5, 4)
    #r.r_forPhloem(24/14, 3)
    #r.r_forPhloem(6, 2) #because roots will have high C_ST less often
    r.expression = 6
    r.update_viscosity = True
    r.solver = 1
    r.atol = 1e-10
    r.rtol = 1e-6
    #r.doNewtonRaphson = False;r.doOldEq = False
    SPAD= 41.0
    chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
    r.Chl = np.array( [chl_]) 
    r.Csoil = 1e-4


    """ for post processing """
    structSumInit = 0
    orgs_all = r.plant.getOrgans(-1, True)

    for org in orgs_all:
        if org.organType() < 2:
            raise Exception("ot < 3")
        structSumInit += org.orgVolume(-1,False) * r.rhoSucrose_f(int(org.getParameter("subType")),
                                                                    org.organType())

    AnSum = 0
    Nt = len(r.plant.nodes) 
    Q_Rmbu      = np.array([0.])

    Q_Grbu      = np.array([0.])
    Q_Gr4bu      = np.array([0.])
    Q_Gr3bu      = np.array([0.])
    Q_Gr2bu      = np.array([0.])
    Q_Grbuth      = np.array([0.])
    Q_Gr4buth      = np.array([0.])
    Q_Gr3buth      = np.array([0.])
    Q_Gr2buth      = np.array([0.])

    Q_Exudbu    = np.array([0.])
    Q_Rmmaxbu   = np.array([0.])
    Q_Grmaxbu   = np.array([0.])
    Q_Exudmaxbu = np.array([0.])
    Q_STbu      = np.array([0.])
    Q_Parbu      = np.array([0.])
    Q_mesobu    = np.array([0.])
    volSegbu =  np.array([0.])
    NOrg = r.plant.getNumberOfOrgans()
    delta_ls_bu = np.full(NOrg, 0.)
    delta_ls_max = 0


    Ntbu = 1
    Q_in  = 0
    Q_out = 0



    orgs_all = r.plant.getOrgans(-1, True)
    volOrgbu = np.array([org.orgVolume(-1,False) for org in orgs_all]) 
    volOrg = np.array([org.orgVolume(-1,False) for org in orgs_all]) 
    ot_orgs = np.array([org.organType() for org in orgs_all])
    #len_orgs = np.array([org.getLength(False) for org in orgs_all])
    st_orgs = np.array([org.getParameter("subType") for org in orgs_all])

    volOrgini = np.array([org.orgVolume(-1,False) for org in orgs_all])

    volOrgi_th = 0.
    lenOrgbu = np.array([org.getLength(False) for org in orgs_all]) 
    lenOrg = np.array([org.getLength(False) for org in orgs_all]) 
    lenOrgi_th = 0.
    Orgidsbu = np.array([org.getId() for org in orgs_all])
    Orgids = np.array([org.getId() for org in orgs_all]) #true:realized
    #raise Exception
    ö=0


    volOrgini_type =  np.array([sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(2, True)]),
                                sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)]), 
                                sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(4, True)])]) 

    volOrgini2 = np.array([org.orgVolume(-1,False) for org in r.plant.getOrgans(2, True)])
    volOrgini3 = np.array([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)])
    volOrgini4 = np.array([org.orgVolume(-1,False) for org in r.plant.getOrgans(4, True)])
    sucOrgini_type =  np.array([sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(2, True)])*r.rhoSucrose_f(0,2),
                                sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)])*r.rhoSucrose_f(0,3), 
                                sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(4, True)])*r.rhoSucrose_f(0,4)]) 
    typeOrg_unit =  np.array([org.organType() for org in r.plant.getOrgans(-1, True)])
    idOrg_unit =  np.array([org.getId() for org in r.plant.getOrgans(-1, True)])
    sucOrg_unit =  np.array([org.orgVolume(-1,False)*r.rhoSucrose_f(int(org.getParameter("subType")),org.organType()) for org in r.plant.getOrgans(-1, True)])
    sucOrgini_unit = sucOrg_unit
    #print(volOrgini2_type, volOrgini_type)
    sucOrg_type = sucOrgini_type

    volOrg_type = volOrgini_type
    orgs_roots = r.plant.getOrgans(2, True)
    orgs_ln = np.array([])
    orgs_st = np.array([])
    for org_ in orgs_roots:
        if org_.getParameter("subType") == 2 :
            orgs_ln= np.append(orgs_ln,len(org_.param().ln) ) 
            orgs_st= np.append(orgs_st,org_.getParameter("subType")  ) 


    orgs_all = r.plant.getOrgans(-1, True)
    volOrgbu = np.array([np.full( (len(org.getNodeIds())-1),org.getId())  for org in orgs_all], dtype=object) 


    beginning = datetime.now()
    #1h for 1d when dxMin = 0.3

    
    Q_ST_init = np.array([])
    Q_meso_init  = np.array([])
    Q_Gr4bu =Q_Gr3bu=Q_Gr2bu=[0]
    deltasucorgbu = np.array([])
    AnSum = 0
    EvSum = 0
    Andt= 0
    Evdt = 0
    Nt = len(r.plant.nodes) 
    Q_Rmbu      = np.array([0.])

    Q_Grbu      = np.array([0.])
    Q_Gr4bu      = np.array([0.])
    Q_Gr3bu      = np.array([0.])
    Q_Gr2bu      = np.array([0.])
    Q_Grbuth      = np.array([0.])
    Q_Gr4buth      = np.array([0.])
    Q_Gr3buth      = np.array([0.])
    Q_Gr2buth      = np.array([0.])

    Q_Exudbu    = np.array([0.])
    Q_Rmmaxbu   = np.array([0.])
    Q_Grmaxbu   = np.array([0.])
    Q_Exudmaxbu = np.array([0.])
    Q_STbu      = np.array([0.])
    Q_Parbu      = np.array([0.])
    Q_mesobu    = np.array([0.])


    def resistance2conductance(resistance,r):
        resistance = resistance* (1/100) #[s/m] * [m/cm] = [s/cm]
        resistance = resistance * r.R_ph * weatherX["TairK"] / r.Patm # [s/cm] * [K] * [hPa cm3 K−1 mmol−1] * [hPa] = [s] * [cm2 mmol−1]
        resistance = resistance * (1000) * (1/10000)# [s cm2 mmol−1] * [mmol/mol] * [m2/cm2] = [s m2 mol−1]
        return 1/resistance
    

    dynamic_soil_first = True
    dtVTP = 10
    dtPrint = 10
    rootLength=0
    leafLength=0
    while simDuration < simMax: 

        print('simDuration:',simDuration )

        ot_4phloem = r.plant.organTypes # np.insert(,0,2)
        ot_4phloem.insert(0,2)#first node 
        ot_4phloem = np.array(ot_4phloem)
        
        hp = max([tempnode[2] for tempnode in r.get_nodes()]) /100 #maxnode canopy [m]
        #print([tempnode[2] for tempnode in r.get_nodes()], hp)

        weatherX = weather(simDuration, hp)
        r.Patm = weatherX["Pair"]
        ##resistances
        r.g_bl = resistance2conductance(weatherX["rbl"],r) / r.a2_bl
        r.g_canopy = resistance2conductance(weatherX["rcanopy"],r) / r.a2_canopy
        r.g_air = resistance2conductance(weatherX["rair"],r) / r.a2_air

        r.Qlight = weatherX["Qlight"] #; TairC = weatherX["TairC"] ; text = "night"


        r = setKrKx_xylem(weatherX["TairC"], weatherX["RH"],r)
        
        dynamic_soil = ((simDuration > simStartSim) and (simDuration <= simStartSim +spellDuration))
        if dynamic_soil and dynamic_soil_first:
            dynamic_soil_first = False
            s = RichardsWrapper(RichardsSP())
            s.initialize()
            periodic = True
            s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
            s.setHomogeneousIC(weatherX["p_mean"], True)  # cm pressure head, equilibrium
            s.setTopBC("noFlux")
            s.setBotBC("fleeFlow")
            s.setVGParameters([weatherX['vg']])
            s.initializeProblem()
            s.setCriticalPressure(-15000)
            sx = s.getSolutionHead()  # inital condition, solverbase.py
            picker = lambda x, y, z: s.pick([x, y, z])    
            pl.setSoilGrid(picker)  # maps segment
            dtWater = 1/24/6
        if not dynamic_soil:
            sx = sx_static_bu
            picker = lambda x,y,z : max(int(np.floor(-z)),-1) 
            pl.setSoilGrid(picker)  # maps segment
            dtWater = dt#1/24
        
        r.es = weatherX["es"]
        
        if dynamic_soil:
            waterDt = 1/2/24
            
            
        dtWatertot = 0
        Evdt =0
        Andt = 0
        print(dtWatertot , dt,Andt)
        while dtWatertot < dt:
            dtWatertot += dtWater
            print("#### IN WATER ####",dtWatertot , dtWater)

            r.solve_photosynthesis(sim_time_ = simDuration, sxx_=sx, cells_ = True,ea_ = weatherX["ea"],
                verbose_ = False, doLog_ = False,TairC_= weatherX["TairC"] )

            #trans = np.array(r.outputFlux)
            """ dumux """   
            fluxesSoil = r.soilFluxes(simDuration, r.psiXyl, sx, approx=False)
            doSmall = False
            #minSoilVal
            if dynamic_soil:
                if doSmall:
                    dtdiff =0 
                    while(dtdiff < dt):
                        dtdiff += dt/100
                        s.setSource(fluxesSoil.copy())  # richards.py 
                        s.solve(dt/100)
                        sx = s.getSolutionHead()  # richards.py    
                        min_sx, min_rx, max_sx, max_rx = np.min(sx), np.min(r.psiXyl), np.max(sx), np.max(r.psiXyl)
                        n = round((simDuration- simInit)/(simMax-simInit) * 100.)

                        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g} days {:g}"
                                .format(min_sx, max_sx, min_rx, max_rx, s.simTime, r.psiXyl[0]))
                else:
                    s.setSource(fluxesSoil.copy())  # richards.py 
                    s.solve(dtWater)
                    sx = s.getSolutionHead()  # richards.py  
                    if False:
                        try:
                            s.solve(dt)
                            sx = s.getSolutionHead()  # richards.py    
                        except:
                            print("DUMUXFAILURE")
                            s = RichardsWrapper(RichardsSP())
                            s.initialize()
                            periodic = True
                            s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
                            s.setHomogeneousIC(np.mean(sx), True)  # cm pressure head, equilibrium
                            s.setTopBC("noFlux")
                            s.setBotBC("fleeFlow")
                            s.setVGParameters([weatherX['vg']])
                            s.initializeProblem()
                            #s.setCriticalPressure(-15000)
                            sx = s.getSolutionHead()  # inital condition, solverbase.py
                    min_sx, min_rx, max_sx, max_rx = np.min(sx), np.min(r.psiXyl), np.max(sx), np.max(r.psiXyl)
                    n = round((simDuration- simInit)/(simMax-simInit) * 100.)

                    print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g} days {:g}"
                            .format(min_sx, max_sx, min_rx, max_rx, s.simTime, r.psiXyl[0]))
                    if False:# min_sx < -20000:
                        s = RichardsWrapper(RichardsSP())
                        s.initialize()
                        periodic = True
                        s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
                        s.setHomogeneousIC(min(np.mean(sx),), True)  # cm pressure head, equilibrium
                        s.setTopBC("noFlux")
                        s.setBotBC("fleeFlow")
                        s.setVGParameters([weatherX['vg']])
                        s.initializeProblem()
                        #s.setCriticalPressure(-15000)
                        sx = s.getSolutionHead()  # inital condition, solverbase.py



            #print(r.Ag4Phloem)
            AnSum += np.sum(r.Ag4Phloem)*dtWater
            Andt += np.array(r.Ag4Phloem)*dtWater
            EvSum += np.sum(r.Ev)*dtWater
            Evdt += np.sum(r.Ev)*dtWater
        r.Ag4Phloem = Andt/dt
        startphloem= simDuration
        endphloem = startphloem + dt
        stepphloem = 1
        #filename = "results/"++"/pmincpb_" + str(simDuration) + ".txt" 

        errLeuning = sum(r.outputFlux)
        fluxes = np.array(r.outputFlux)
        # save_stdout = sys.stdout
        # sys.stdout = open('trash', 'w')
        verbose_phloem = True
        #deltaVol = r.calcDeltaVolOrgNode(dt,1)
        #print(deltaVol)
        #print("start pm")
        #raise Exception
        filename = "results"+ directoryN +"inPM_"+str(ö)+".txt"
        print("startpm", hp)
        
        try:
            r.startPM(startphloem, endphloem, stepphloem, ( weatherX["TairC"]  +273.15) , verbose_phloem, filename)
        except:
            #print("startpm WITH troubleshoot NONO", hp)
            #r.doTroubleshooting = True
            r.startPM(startphloem, endphloem, stepphloem, ( weatherX["TairC"]  +273.15) , verbose_phloem, filename)
            
        try:
            os.remove(filename)
        except OSError:
            pass    
        
        if r.withInitVal and (len(Q_ST_init) ==0) :
            Q_ST_init = np.array(r.Q_init[0:Nt])
            Q_meso_init = np.array(r.Q_init[Nt:(Nt*2)])

        #print("start pm_done")
        # sys.stdout = save_stdout

        Q_ST    = np.array(r.Q_out[0:Nt])
        Q_meso  = np.array(r.Q_out[Nt:(Nt*2)])
        Q_Rm    = np.array(r.Q_out[(Nt*2):(Nt*3)])
        Q_Exud  = np.array(r.Q_out[(Nt*3):(Nt*4)])
        Q_Gr    = np.array(r.Q_out[(Nt*4):(Nt*5)])
        Q_Gr4       = Q_Gr[np.where(ot_4phloem==4)[0]]#Q_Gr4     - Q_Gr4bu
        Q_Gr3       = Q_Gr[np.where(ot_4phloem==3)[0]]#Q_Gr3     - Q_Gr3bu
        Q_Gr2       = Q_Gr[np.where(ot_4phloem==2)[0]]#Q_G#r2     - Q_G#r2bu


        C_ST    = np.array(r.C_ST)
        Q_Par   = np.array(r.Q_out[(Nt*8):(Nt*9)])
        Fl      = np.array(r.Fl)
        volST   = np.array(r.vol_ST)
        volMeso   = np.array(r.vol_Meso)
        C_Par   = Q_Par/volST
        C_meso  = Q_meso/volMeso
        Q_in   += sum(np.array(r.AgPhl)*dt)
        Q_out   = Q_Rm + Q_Exud + Q_Gr
        error   = sum(Q_ST +Q_Par+ Q_meso + Q_out )- Q_in - sum(Q_ST_init)  - sum(Q_meso_init)
        # Q_ST_dot    = np.array(r.Q_out_dot[0:Nt])
        # Q_meso_dot  = np.array(r.Q_out_dot[Nt:(Nt*2)])
        # Q_Rm_dot    = np.array(r.Q_out_dot[(Nt*2):(Nt*3)])
        # Q_Exud_dot  = np.array(r.Q_out_dot[(Nt*3):(Nt*4)])
        # Q_Gr_dot    = np.array(r.Q_out_dot[(Nt*4):(Nt*5)])

        #delta_ls_max += sum(np.array(#r2.rmaxSeg(dt, r.k_gr)) * dt)
        delta_ls_max_i = np.array(r.delta_ls_org_imax)
        delta_ls_max = np.array(r.delta_ls_org_max)
        delta_ls_i = np.array(r.delta_ls_org_i)
        delta_ls = np.array(r.delta_ls_org)

        Q_Rmmax       = np.array(r.Q_out[(Nt*5):(Nt*6)])
        Q_Grmax       = np.array(r.Q_out[(Nt*6):(Nt*7)])
        Q_Exudmax     = np.array(r.Q_out[(Nt*7):(Nt*8)])

        Q_ST_i        = Q_ST      - Q_STbu
        Q_Par_i       = Q_out     - Q_Parbu
        Q_Rm_i        = Q_Rm      - Q_Rmbu
        Q_Gr_i        = Q_Gr      - Q_Grbu
        Q_Gr4_i       = Q_Gr_i[np.where(ot_4phloem==4)[0]]#Q_Gr4     - Q_Gr4bu
        Q_Gr3_i       = Q_Gr_i[np.where(ot_4phloem==3)[0]]#Q_Gr3     - Q_Gr3bu
        Q_Gr2_i       = Q_Gr_i[np.where(ot_4phloem==2)[0]]#Q_G#r2     - Q_G#r2bu

        #Q_Gr_ith        = Q_Grth      - Q_Grbuth
        #Q_Gr4_ith       = Q_Gr4th     - Q_Gr4buth
        #Q_Gr3_ith       = Q_Gr3th     - Q_Gr3buth
        #Q_G#r2_ith       = Q_G#r2th     - Q_G#r2buth

        Q_Exud_i      = Q_Exud    - Q_Exudbu
        Q_meso_i      = Q_meso    - Q_mesobu

        Q_Rmmax_i     = Q_Rmmax   - Q_Rmmaxbu
        Q_Grmax_i     = Q_Grmax   - Q_Grmaxbu
        Q_Exudmax_i   = Q_Exudmax - Q_Exudmaxbu

        Q_out_i       = Q_Rm_i    + Q_Exud_i      + Q_Gr_i
        Q_outmax_i    = Q_Rmmax_i + Q_Exudmax_i   + Q_Grmax_i


        orgs = r.plant.getOrgans(-1, True)
        id_orgs = np.array([org.getId() for org in orgs])
        orgs_all = r.plant.getOrgans(-1, True)
        ot_orgs_all = np.array([org.organType() for org in orgs_all])
        volOrgi_th = 0#volOrg2 - volOrgini2

        volOrg2 = np.array([org.orgVolume(-1,False) for org in r.plant.getOrgans(2, True)])
        volOrg3 = np.array([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)])
        volOrg4 = np.array([org.orgVolume(-1,False) for org in r.plant.getOrgans(4, True)])

        volOrg_typei = volOrg_type - volOrgini_type
        #volOrg2_typei = volOrg2_type - volOrgini2_type

        JW_ST = np.array(r.JW_ST)
        length_ST = np.array(r.plant.segLength())
        #0.0001037
        Lp = 0.005#0.004320 #cm d-1 hPa, assumed resistance between xylem and phloem
        Rgaz =  83.14 #hPa cm3 K-1 mmol-1
        a_STs = np.array(r.a_ST)#np.array([a_ST[ot][st] for ot, st in ])
        #RhatFhat =   (weatherX["TairC"] + 273.15) * C_ST[1:] * Rgaz * (2/a_STs) * length_ST * Lp /np.abs(JW_ST) /   ((a_STs**2)*np.pi)  
        RhatFhat =   (weatherX["TairC"] + 273.15) * C_ST[1:] * Rgaz * 2 /a_STs* length_ST * Lp /np.abs(JW_ST) *   (25*a_STs*a_STs*np.pi) 
        C_ST_ = C_ST[1:]
        ids = np.where(RhatFhat ==  min(RhatFhat))
        
        rootLength = sum([mr.getLength(False) for mr in r.plant.getOrgans(2, False)])
        leafLength = sum([mr.getLength(False) for mr in r.plant.getOrgans(4, False)])
        if False:#min(RhatFhat) < 1) :
            #print()
            C_ST_ = C_ST[1:]
            ids = np.where(RhatFhat ==  min(RhatFhat))
            print(min(RhatFhat))
            print(C_ST_[ids] , Rgaz  )
            print(a_STs[ids]  )
            print(length_ST[ids]   )
            print(JW_ST[ids]) 
            print( RhatFhat[ids],(weatherX["TairC"] + 273.15)  )
            print("issue RhatFhat")
            #raise Exception

        if verbose :
            print("\n\n\n\t\tat ", int(np.floor(simDuration)),"d", int((simDuration%1)*24),"h",  round(r.Qlight *1e6),"mumol m-2 s-1")
            print("Error in Suc_balance:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(error, div0f(error,Q_in, 1.)))
            #print("Error in growth:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(errorGri, relErrorGri))
            print("Error in photos:\n\tabs (cm3/day) {:5.2e}".format(errLeuning))
            print("water fluxes (cm3/day):\n\ttrans {:5.2e}\tminExud {:5.2e}\tmaxExud {:5.2e}".format(sum(fluxesSoil.values()), min(fluxesSoil.values()), max(fluxesSoil.values())))
            print("C_ST (mmol ml-1):\n\tmean {:.2e}\tmin  {:5.2e} at {:d} segs \tmax  {:5.2e}".format(np.mean(C_ST), min(C_ST), len(np.where(C_ST == min(C_ST) )[0]), max(C_ST)))        
            print("C_me (mmol ml-1):\n\tmean {:.2e}\tmin  {:5.2e}\tmax  {:5.2e}".format(np.mean(C_meso), min(C_meso), max(C_meso)))        
            print('Q_X (mmol Suc): \n\tST   {:.2e}\tmeso {:5.2e}\tin   {:5.2e}'.format(sum(Q_ST), sum(Q_meso), Q_in))
            print('\tRm   {:.2e}\tGr   {:.2e}\tExud {:5.2e}'.format(sum(Q_Rm), sum(Q_Gr), sum(Q_Exud)))
            print("aggregated sink satisfaction at last time step (%) :\n\ttot  {:5.1f}\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(
                sum(Q_out_i)/sum(Q_outmax_i)*100,sum(Q_Rm_i)/sum(Q_Rmmax_i)*100, 
                 div0f(sum(Q_Gr_i), sum(Q_Grmax_i), 1.)*100,div0f(sum(Q_Exud_i),sum(Q_Exudmax_i), 1.)*100))
            print("aggregated sink repartition at last time step (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rm_i)/sum(Q_out_i)*100, 
                 sum(Q_Gr_i)/sum(Q_out_i)*100,sum(Q_Exud_i)/sum(Q_out_i)*100))
            print("aggregated sink repartition (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rm)/sum(Q_out)*100, 
                 sum(Q_Gr)/sum(Q_out)*100,sum(Q_Exud)/sum(Q_out)*100))
            print("aggregated sink repartition for max (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rmmax_i)/sum(Q_outmax_i)*100, 
                 sum(Q_Grmax_i)/sum(Q_outmax_i)*100,sum(Q_Exudmax_i)/sum(Q_outmax_i)*100))
            print("abs val for max :\n\tRm   {:5.5f}\tGr   {:5.5f}\tExud {:5.5f}".format(sum(Q_Rmmax_i), 
                 sum(Q_Grmax_i),sum(Q_Exudmax_i)))
            print("Q_Par {:5.2e}, C_Par {:5.2e}".format(sum(Q_Par), np.mean(C_Par)))
            print("growth (cm)\ttot {:5.2e}\ti {:5.2e}".format(sum(delta_ls), sum(delta_ls_i)))      
            print("max growth (cm)\ttot {:5.2e}\ti {:5.2e}".format(sum(delta_ls_max), sum(delta_ls_max_i)))     
            print("amount Suc (cm)\tAn {:5.2e}\tGr {:5.2e}\tRGr {:5.2e}\tRm {:5.2e}\tExud {:5.2e}".format(AnSum, sum(Q_Gr)*r.Gr_Y,sum(Q_Gr)*(1-r.Gr_Y), sum(Q_Rm), sum(Q_Exud))) 
            #print("growth (cm3)\n\tobs {:5.2e}\ttotth {:5.2e}\ttotobs {:5.2e}".format(sum(volOrgi_obs), volOrgi_th,volOrg_obstot ))  
            print("growth (cm3) per type\n\ttotobs", volOrg_typei)# , volOrg2_typei)       
            print("sucOrg obs (mmol)\t th (mmol)\t", sucOrg_type - sucOrgini_type)#, sucOrg2_type - sucOrgini2_type)       
            print("Grobs (mmol) root\tstem\tleaf\t", sum(Q_Gr2bu)*r.Gr_Y,sum(Q_Gr3bu)*r.Gr_Y, sum(Q_Gr4bu)*r.Gr_Y)# , gr4ith) 
            print("RhatFhat ", min(RhatFhat))#,C_ST_[ids],a_STs[ids], length_ST[ids], JW_ST[ids]  )
            print("root length {:5.2e}\tlef length {:5.2e}".format(rootLength,leafLength))#,C_ST_[ids],a_STs[ids], length_ST[ids], JW_ST[ids]  )
        #raise Exception     
        if min(C_ST) < 0.0:
            print("min(C_ST) < 0.015", min(C_ST),np.mean(C_ST),max(C_ST))
            raise Exception
        write_file_array("RhatFhat", RhatFhat)
        write_file_array("a_STs", a_STs)
        write_file_array("JW_ST", JW_ST)#cm3/d
        write_file_array("length_ST", length_ST)
        #NB: delta_ls_max inexact as does not account for growth of organxs which are not here at the beginning
        #print([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)],[org.orgVolume(-1,False) for org in #r2.plant.getOrgans(3, True)])
        #print([org.getParameter("lengthTh") for org in r.plant.getOrgans(3, True)],[org.getParameter("lengthTh") for org in #r2.plant.getOrgans(3, True)])
        # print([org.getId() for org in r.plant.getOrgans(3, True)],[org.orgVolume(-1,False) for org in #r2.plant.getOrgans(3, True)])
        # print(len(r.plant.getOrgans(2, True)),len(r.plant.getOrgans(3, True)),len(r.plant.getOrgans(4, True)))
        #dist2tip = np.array(r.plant.dist2tips)
        ana = pb.SegmentAnalyser(r.plant.mappedSegments())

        cutoff = 1e-15 #is get value too small, makes paraview crash
        C_ST_p = C_ST
        C_ST_p[abs(C_ST_p) < cutoff] = 0
        fluxes_p = fluxes
        fluxes_p[abs(fluxes_p) < cutoff] = 0
        #print(fluxes,fluxes_p)
        Q_Exud_i_p = Q_Exud_i
        Q_Exud_i_p[abs(Q_Exud_i_p) < cutoff] = 0
        Q_Rm_i_p = Q_Rm_i
        Q_Rm_i_p[abs(Q_Rm_i_p) < cutoff] = 0
        Q_Gr_i_p = Q_Gr_i
        Q_Gr_i_p[abs(Q_Gr_i_p) < cutoff] = 0

        Q_Exudmax_i_p = Q_Exudmax_i
        Q_Exudmax_i_p[abs(Q_Exudmax_i_p) < cutoff] = 0
        Q_Rmmax_i_p = Q_Rmmax_i
        Q_Rmmax_i_p[abs(Q_Rmmax_i_p) < cutoff] = 0
        Q_Grmax_i_p = Q_Grmax_i
        Q_Grmax_i_p[abs(Q_Grmax_i_p) < cutoff] = 0


        C_Exud_i_p = Q_Exud_i/volST
        C_Exud_i_p[abs(C_Exud_i_p ) < cutoff] = 0
        C_Rm_i_p = Q_Rm_i/volST
        C_Rm_i_p[abs(C_Rm_i_p) < cutoff] = 0
        C_Gr_i_p = Q_Gr_i/volST
        C_Gr_i_p[abs(C_Gr_i_p) < cutoff] = 0

        C_Exudmax_i_p = Q_Exudmax_i/volST
        C_Exudmax_i_p[abs(C_Exudmax_i_p) < cutoff] = 0
        C_Rmmax_i_p = Q_Rmmax_i/volST
        C_Rmmax_i_p[abs(C_Rmmax_i_p) < cutoff] = 0
        C_Grmax_i_p = Q_Grmax_i/volST
        C_Grmax_i_p[abs(C_Grmax_i_p) < cutoff] = 0

        psiXyl_p = np.array(r.psiXyl)
        psiXyl_p[abs(psiXyl_p) < cutoff] = 0
        psi_p_symplasm_p = np.array(r.psi_p_symplasm)
        psi_p_symplasm_p[abs(psi_p_symplasm_p) < cutoff] = 0
        
        
        ana.addData("CST", C_ST_p)
        #do as konz or div per vol or surf?
        #ana.addData("Q_Exud", Q_Exud)  # cut off for vizualisation
        ana.addData("fluxes", fluxes_p)
        ana.addData("Fpsi", np.array(r.Fpsi))
        ana.addData("Flen", np.array(r.Flen))
        ana.addData("GrowthZone", np.array(r.GrowthZone))
        ana.addData("GrowthZoneLat", np.array(r.GrowthZoneLat))

        ana.addData("QExud", Q_Exud_i_p)  # cut off for vizualisation
        ana.addData("QRm", Q_Rm_i_p)  # cut off for vizualisation
        ana.addData("QGr", Q_Gr_i_p)  # cut off for vizualisation
        ana.addData("QExudmax", Q_Exudmax_i_p)  # cut off for vizualisation
        ana.addData("QRmmax", Q_Rmmax_i_p)  # cut off for vizualisation
        ana.addData("QGrmax", Q_Grmax_i_p)  # cut off for vizualisation

        ana.addData("CExud", C_Exud_i_p)  # cut off for vizualisation
        ana.addData("CRm", C_Rm_i_p)  # cut off for vizualisation
        ana.addData("CGr", C_Gr_i_p)  # cut off for vizualisation
        ana.addData("CExudmax", C_Exudmax_i_p)  # cut off for vizualisation
        ana.addData("CRmmax", C_Rmmax_i_p)  # cut off for vizualisation
        ana.addData("CGrmax", C_Grmax_i_p)  # cut off for vizualisation

        ana.addData("psi_Xyl",psiXyl_p)
        ana.addData("psi_p_symplasm_p",psi_p_symplasm_p)
        if True:
            ana.write("results"+directoryN+"plot_"+str(int(simStartSim))+str(condition)+"at"+ str(ö) +".vtp", 
                      ["organType", "subType",
                       "CST", "fluxes","psi_Xyl",
                       "QExud", "QGr", "QRm",
                       "CExud", "CGr", "CRm",
                       "QExudmax", "QGrmax", "QRmmax",
                       "CExudmax", "CGrmax", "CRmmax",
                       "Fpsi","Flen","GrowthZone","GrowthZoneLat",
                      "psi_p_symplasm"]) 

    #     fluxes_p = np.insert(fluxes_p,0,0)# "[sucrose]",
        dtVTP += dt
        print(simDuration,dtVTP)
        if((simDuration > 0) and (dtVTP >=0.29/24)):
            dtVTP = 0
            vp.plot_plant(r.plant,p_name = [ "fluxes","Exud","Gr","Rm","xylem pressure (cm)","sucrose concentration (mmol/cm3)"],
                                vals =[ fluxes_p, Q_Exud_i_p, Q_Gr_i_p, Q_Rm_i_p, psiXyl_p, C_ST_p], 
                                filename = "results"+ directoryN +"plotpsi_"+str(int(simStartSim))+str(condition)+"at"+ str(ö), 
                          range_ = [0,5000])
            #vp.plot_plant(r.plant,p_name = [ "fluxes","Exud","Gr","Rm","sucrose concentration (mmol/cm3)"],
             #                   vals =[ fluxes_p, Q_Exud_i_p, Q_Gr_i_p, Q_Rm_i_p, C_ST_p], 
              #                  filename = "results"+ directoryN +"plotsuc_"+str(int(simStartSim))+str(condition)+"at"+ str(ö), range_ = [0,3])   
        ö +=1
        #raise Exception()
        r_ST_ref = np.array(r.r_ST_ref)
        r_ST = np.array(r.r_ST)

        #with open("results"+ directoryN +"CWGr_max_15pm.txt",  'a') as data: 
         #     data.write(str(r.deltaSucOrgNode))

        ots = np.concatenate((np.array([0]), r.get_organ_types()))#per node
        #ots_org = np.concatenate((np.array([0]), r.get_organ_types()))#per org
        #write_file_array("id_org", ot_orgs_all)
        #write_file_array("ots_org", ot_orgs_all)
        #sucOrg_unit =  np.array([org.orgVolume(-1,False)*r.rhoSucrose[org.organType()] for org in r.plant.getOrgans(-1, True)])
        #typeOrg_unit =  np.array([org.organType() for org in r.plant.getOrgans(-1, True)])
        
        dtPrint += dt

        if dtPrint >= 0.9/24:
            dtPrint = 0
            rl_ = ana.distribution("length", 0., -depth, layers, True)                   
            rl_ = np.array(rl_)/ soilvolume  # convert to density
            write_file_array("RLD", rl_)
            write_file_array("Fpsi", r.Fpsi)
            #write_file_array("deltasucorgbu_type", typeOrg_unit)
            #write_file_array("deltasucorgbu_phloem", deltasucorgbu)
            #write_file_array("deltasucorgbu_plant",  (sucOrg_unit - sucOrgini_unit)[idOrg_unit.argsort()])
            #write_file_array("deltasucorgbu_plantid",  idOrg_unit[idOrg_unit.argsort()])

            #write_file_array("volOrg2", volOrg2-volOrgini2)
            #write_file_array("volOrg3", volOrg3-volOrgini3)
            #write_file_array("volOrg4", volOrg4-volOrgini4)

            write_file_array("leafBladeSurface", np.array(r.plant.leafBladeSurface))
            write_file_array("fw", r.fw)
            write_file_array("gco2", r.gco2)
            write_file_array("pvd", r.PVD)
            write_file_array("pg", r.pg)
            write_file_float("ea", r.ea)
            write_file_float("es", r.es)
            write_file_float("psi_air", r.psi_air)
            write_file_float("RH", weatherX["RH"])
            write_file_array("EAL", r.EAL)
            write_file_array("hrelL", r.hrelL)
            write_file_array("ci", r.ci)
            write_file_array("k_stomatas", r.k_stomatas)

            #write_file_array("length_blade", length_blade)
            write_file_array("ots", ots)
            write_file_array("soilWatPot", sx)
            write_file_array("fluxes", fluxes)#cm3 day-1
            write_file_array("volST", volST)
            write_file_array("volOrg",  volOrg) #with epsilon
            write_file_array("Fl", Fl)
            write_file_array("AgPhl", np.array(r.AgPhl))
            write_file_array("Q_ST_dot", Q_ST_i/dt)
            write_file_array("Q_meso_dot", Q_meso_i/dt)
            write_file_array("Q_Rm_dot", Q_Rm_i/dt)
            write_file_array("Q_Exud_dot", Q_Exud_i/dt)
            write_file_array("Q_Gr_dot", Q_Gr_i/dt)
            write_file_array("Q_Rmmax_dot", Q_Rmmax_i/dt)
            write_file_array("Q_Exudmax_dot", Q_Exudmax_i/dt)
            write_file_array("Q_Grmax_dot", Q_Grmax_i/dt)

            #write_file_array("Q_Par", Q_Par)
            #write_file_array("C_Par", C_Par)
            write_file_array("Q_ST", Q_ST)
            write_file_array("C_ST", C_ST)
            write_file_array("C_meso", C_meso)
            write_file_array("Q_meso", Q_meso)
            write_file_array("Q_Rm", Q_Rm)
            write_file_array("Q_Exud", Q_Exud)
            write_file_array("Q_Gr", Q_Gr)
            write_file_array("psiXyl", r.psiXyl)
            write_file_array("psi_p_symplasm", r.psi_p_symplasm)
            write_file_float("trans", EvSum)
            write_file_float("transrate",Evdt)
            write_file_array("Anrate",Andt)
            
            leavesSegs = np.where(ots[1:] ==4)
            fluxes_leaves = fluxes[leavesSegs]
            if (min(r.Ev) < 0) or (min(r.Jw) < 0) or (min(fluxes_leaves)<0):
                print("leaf looses water", min(r.Ev),min(r.Jw), min(fluxes_leaves))
                raise Exception
            write_file_array("Q_Grmax", Q_Grmax)
            write_file_array("Q_Rmmax", Q_Rmmax)
            write_file_array("Q_Exudmax", r.Q_Exudmax)


            write_file_array("nodes_X",np.array([tempnode[0] for tempnode in r.get_nodes()]))
            write_file_array("nodes_Y", np.array([tempnode[1] for tempnode in r.get_nodes()]))
            write_file_array("nodes_Z", np.array([tempnode[2] for tempnode in r.get_nodes()]))


            write_file_array("ratio_Gr ", div0(Q_Gr_i,Q_Grmax_i, np.nan))
            write_file_array("ratio_Rm ", Q_Rm_i/Q_Rmmax_i)
            write_file_array("ratio_Exud ", div0(Q_Exud_i,Q_Exudmax_i, np.nan))
            write_file_array("satisfaction ", Q_out_i/Q_outmax_i)

            write_file_array("r_ST_ref", r_ST_ref)
            write_file_array("r_ST", r_ST)
            write_file_array("mu", r_ST/r_ST_ref)#hPa d

            write_file_array("id_orgs", id_orgs)
            write_file_array("ot_orgs", ot_orgs)
            write_file_array("ot_orgs_all", ot_orgs)#with arg to small to be represented
            write_file_array("st_orgs", st_orgs)
            write_file_array("len_orgs", lenOrg)
            #write_file_array("delta_ls", delta_ls)
            #write_file_array("Q_Ag", r.AgPhl)
            #write_file_array("delta_ls_max", delta_ls_max)
            #write_file_array("delta_ls_i", delta_ls_i)
            #write_file_array("delta_ls_max_i", delta_ls_max_i)
            write_file_float("time", simDuration)
            write_file_float("computeTime", datetime.now() - beginning)
            organTypes = r.get_organ_types()
            subTypes =r.get_subtypes()
            nodes_organtype = np.concatenate(([1], organTypes))#seg_Indx = node_y_Indx - 1. just need to add organ type for seed (set to type 'root')     #np.zeros(len(nods))#-1)
            nodes_subtype = np.concatenate(([1], subTypes))
            write_file_array("organTypes", organTypes)
            write_file_array("subTypes", subTypes)

            write_file_array("nodesOrganTypes", nodes_organtype)
            write_file_array("nodesSubTypes", nodes_subtype)
        
        Ntbu = Nt
        #Ntbu2 = Nt2
        orgs_all = r.plant.getOrgans(-1, True)
        NOrgbu = len(orgs_all)
        orgradibu= np.array([org.getParameter("a") for org in orgs_all])
        volOrgbu = np.array([org.orgVolume(-1,False) for org in orgs_all])
        lenOrgbu = np.array([org.getLength(False) for org in orgs_all])
        Orgidsbu = np.array([org.getId() for org in orgs_all]) 

        verbose_simulate = False
        r.plant.simulate(dt, verbose_simulate)#, "outputpm15.txt") #time span in days /(60*60*24)
        
        sucOrg_type =  np.array([sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(2, True)])*r.rhoSucrose_f(0,2),
                                sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)])*r.rhoSucrose_f(0,3), 
                                sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(4, True)])*r.rhoSucrose_f(0,4)]) 

        sucOrg_unit =  np.array([org.orgVolume(-1,False)*r.rhoSucrose_f(int(org.getParameter("subType")),org.organType()) for org in r.plant.getOrgans(-1, True)])
        typeOrg_unit =  np.array([org.organType() for org in r.plant.getOrgans(-1, True)])
        idOrg_unit =  np.array([org.getId() for org in r.plant.getOrgans(-1, True)])

        deltasucorgbu = np.array(r.delta_suc_org)   

        #write_file_array("volOrgth", volOrg2_type)
        write_file_array("volOrgobs", volOrg_type)
        #write_file_array("sucOrgth", sucOrg2_type)
        write_file_array("sucOrgobs", sucOrg_type)

        #orgs_all2 = #r2.plant.getOrgans(-1, True)
        #volOrgi_th2 =  sum(np.array([org.orgVolume(-1,False) for org in orgs_all2]) )#true:realized

        orgs = r.plant.getOrgans(-1, True)
        orgs_all = r.plant.getOrgans(-1, True)

        Nt = len(r.plant.nodes)
        NOrg = len(orgs_all)#r.plant.getNumberOfOrgans()#att! not same number
        sucOrgini_unit = np.concatenate((sucOrgini_unit, np.full(NOrg - NOrgbu, 0.)))

        orgradi= np.array([org.getParameter("a") for org in orgs_all])
        orgradibu = np.concatenate((orgradibu, np.full(NOrg - NOrgbu, 0.)))
        #Orgidsbu = Orgids
        Orgids = np.array([org.getId() for org in orgs_all]) 

        Orgidsbu = np.concatenate((Orgidsbu, np.full(NOrg - NOrgbu, 0.)))

        volOrg = np.array([org.orgVolume(-1,False) for org in orgs_all])
        ot_orgs = np.array([org.organType() for org in orgs_all])
        st_orgs = np.array([org.getParameter("subType") for org in orgs_all])
        stroot = st_orgs[np.where(ot_orgs==2)]
        # errorst = sum(np.where(stroot > 3))
        # if errorst > 0:
            # raise Exception

        volOrgbu     = np.concatenate((volOrgbu, np.full(NOrg - NOrgbu, 0.)))#len(volOrg) - len(volOrgbu), 0.)))
        #volOrgi_th = sum(Q_Gr_i)*r.Gr_Y/r.rhoSucrose


        #not sur I get the organs always in same order
        lenOrg = np.array([org.getLength(False) for org in orgs_all]) #true:realized
        lenOrgbu     = np.concatenate((lenOrgbu, np.full(NOrg - NOrgbu, 0.)))#len(volOrg) - len(volOrgbu), 0.)))
        lenOrgi_th = np.concatenate((delta_ls_i, np.full(NOrg - NOrgbu, 0.)))

        #delta_ls_max =  np.concatenate((delta_ls_max, np.full(Nt2 - Ntbu2, 0.)))
        Q_Rmbu       =   np.concatenate((Q_Rm, np.full(Nt - Ntbu, 0.)))
        Q_Grbu       =   np.concatenate((Q_Gr, np.full(Nt - Ntbu, 0.))) 
        Q_Exudbu     =   np.concatenate((Q_Exud, np.full(Nt - Ntbu, 0.))) 

        Q_Rmmaxbu    =   np.concatenate((Q_Rmmax, np.full(Nt - Ntbu, 0.)))
        Q_Grmaxbu    =   np.concatenate((Q_Grmax, np.full(Nt - Ntbu, 0.))) 
        Q_Exudmaxbu  =   np.concatenate((Q_Exudmax, np.full(Nt - Ntbu, 0.))) 

        Q_STbu       =   np.concatenate((Q_ST, np.full(Nt - Ntbu, 0.)))
        Q_mesobu     =   np.concatenate((Q_meso, np.full(Nt - Ntbu, 0.)))
        delta_ls_bu  =   np.concatenate((delta_ls, np.full(NOrg - NOrgbu, 0.)))

        Q_Gr4bu= Q_Gr4
        Q_Gr3bu =Q_Gr3
        Q_Gr2bu= Q_Gr2
        ##
        #       CHECKS
        ##
        highNeed1 = sum(Q_Rm_i[np.greater(Q_Rm_i,Q_Rmmax_i)] - Q_Rmmax_i[np.greater(Q_Rm_i,Q_Rmmax_i)])
        highNeed2 = sum(Q_Gr_i[np.greater(Q_Gr_i,Q_Grmax_i)] - Q_Grmax_i[np.greater(Q_Gr_i,Q_Grmax_i)])
        highNeed3 = sum(Q_Exud_i[np.greater(Q_Exud_i,Q_Exudmax_i)] - Q_Exudmax_i[np.greater(Q_Exud_i,Q_Exudmax_i)])

        if highNeed1 > abs(error):
            print(np.where([np.greater(Q_Rm_i,Q_Rmmax_i)]))
        #assert div0f(error,Q_in,1.) < 1e-3, "balance error > 0.1%"

        if (len(Orgids) != len(Orgidsbu)):
            print(len(Orgids), len(Orgidsbu), NOrg, NOrgbu, Nt, Ntbu)
        orderId  = Orgids - Orgidsbu

        lenOrgi_obs = lenOrg - lenOrgbu
        lenOrgi_obs2 = sum(lenOrg) - sum(lenOrgbu)
        errorleni = lenOrgi_obs -  lenOrgi_th 
        errorleni2 = lenOrgi_obs2 -  sum(lenOrgi_th) 

        simDuration += dt
        #raise Exception
        # if  abs(errorleni2) > 1e-3:
            # print(lenOrgi_obs)
            # print(lenOrgi_th)
            # print(errorleni)
            # print(lenOrgi_obs2,sum(lenOrgi_th) )
            # print(errorleni2)
            # print( "len error")
            # raise Exception


        # volOrgi_obs = volOrg - volOrgbu
        # volOrgi_obs2 = sum(volOrg) - sum(volOrgbu)
        # errorGri = volOrgi_obs2 -  volOrgi_th 
        # relErrorGri = div0f(errorGri,volOrgi_th ,1.)
        # assert abs(errorGri) < 1e-3, "vol errorend"
    structSum = 0
    orgs_all = r.plant.getOrgans(-1, True)

    for org in orgs_all: #cm3 * mmol Suc /cm3 = mmol Suc 
        structSum += org.orgVolume(-1,False) * r.rhoSucrose_f(int(org.getParameter("subType")),org.organType())


    GrowthSum = structSum - structSumInit
    print("simDuration", simDuration, "d")
    end = datetime.now()
    print(end - beginning)

if __name__ == '__main__':
    simInit = sys.argv[1]
    simStartSim = sys.argv[2]
    condition = sys.argv[3]
    csChoise = "noCS"#sys.argv[4]
    spellDuration = 7
    spellEnd = str(int(simStartSim) + spellDuration)
    simEnd = str(28)
    directoryN = "/"+os.path.basename(__file__)[:-3]+"/"+simInit+"to"+simStartSim+"to"+simEnd+condition+"lowerRez/"

    simInit = float(simInit)
    simStartSim = float(simStartSim)

    main_dir=os.environ['PWD']#dir of the file
    results_dir = main_dir +"/results"+directoryN
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    else:
        import shutil
        shutil.rmtree(results_dir)
        os.makedirs(results_dir)


    launchUQR(directoryN,simInit,simStartSim, condition,csChoise,spellDuration)