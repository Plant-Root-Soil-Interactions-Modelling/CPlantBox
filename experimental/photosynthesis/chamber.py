import base64
import os
import sys
import warnings
# warnings.filterwarnings("ignore", category=DeprecationWarning)
# warnings.filterwarnings("ignore", category=FutureWarning)
# warnings.filterwarnings("ignore", category=UserWarning)
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import constants
import pandas as pd
import sys; sys.path.append("../.."); sys.path.append("../../src")
import datetime
import time
# import json
# from enum import Enum
# # progress bar
# from tqdm import tqdm

# if the runtime folder is not the python script folder, switch there
if os.path.dirname(os.path.abspath(__file__)) != os.getcwd():
  os.chdir(os.path.dirname(os.path.abspath(__file__)))


class Soil :
  def __init__(self) :
    self.source = pd.read_csv("./Soil_parameter.csv")
  #enddef
  def get_pressure_head(self, depth) :
    hydroprop = self.source["Hyprop"]
    measuredepth = self.source["depth"]
    return np.interp(depth, measuredepth, hydroprop)
  #enddef
  def soil_matric_potential(self, depth) :
    return self.get_pressure_head(depth)
  #enddef
#endclass


class Weather :
  def __init__(self, Lat, Long, Time ) :
    self.Lat = Lat
    self.Long = Long
    self.Time = Time
    self.Day = self.Time.day
    with open ("SE_EC_001.1711883710395.csv", "r") as f :
      # get line 93
      for i in range(92) :
        f.readline()
      # get the column data
      line = f.readline()[1:]
      line = [c.strip() for c in line.split(",")]
      self.columns = line
    # endwith
    self.data = pd.read_csv("SE_EC_001.1711883710395.csv", skiprows = 93, names = self.columns)
    for col in self.data.columns:
      if "QualityFlag" in col:
        self.data[col] = self.data[col].str.split("_", expand=True)[1]
        # convert the values to numeric
        self.data[col] = pd.to_numeric(self.data[col], errors='coerce')
      # endif
    # endfor
    # exchange "noData" with NaN
    self.data = self.data.replace("noData", np.nan)
    # convert the time to datetime
    self.data["Time"] = pd.to_datetime(self.data["Time"], format="%Y-%m-%dT%H:%M:%S%z")
    self.data = self.data.set_index("Time")
    self.quality_flags = self.data.filter(like="QualityFlag")
    self.data = self.data.drop(columns = self.quality_flags.columns)
    self.data = self.data.drop(columns="feature")
    self.data = self.data.apply(pd.to_numeric, errors='coerce')
    self.data.sort_index(inplace=True)
    self.data.index = pd.to_datetime(self.data.index, format="%Y-%m-%dT%H:%M:%S%z", utc=True)
    self.quality_flags.index = self.data.index
    self.data = self.data.asfreq("10T")
    self.quality_flags = self.quality_flags.asfreq("10T")
  #enddef
  def fill_nans(self) :
    self.data = self.data.interpolate(method="time")
  #enddef
  def __call__(self, tp, column) :
    # get the closest time point
    timest = pd.Timestamp(tp, tz="UTC")
    timed = self.data.index.get_indexer([timest], method="nearest")
    return self.data[column][timed].values[0]
  #enddef
  def print_columns(self) :
    print(self.data.columns)
  #enddef
  def relative_humidity(self, time) :
    absolute_humidity = self.__call__(time, "AirHumidity_2m_Avg10min [g*m-3]")
    temperature = self.__call__(time, "AirTemperature_2m_Avg30min [°C]")
    # calculate the relative humidity
    es = 6.1078 * 10 ** (7.5 * temperature / (237.3 + temperature))
    return absolute_humidity / es
  #enddef
  def wind_speed(self, time) :
    return self.__call__(time, "WindSpeed_10m_Avg10min [m/s]")
  #enddef
  def wind_direction(self, time) :
    return self.__call__(time, "WindDirection_2.52m_Avg10min [°N]")
  #enddef
  def air_pressure(self, time) :
    return self.__call__(time, "AirPressure_1m_Avg10min [mbar]")
  #enddef
  def radiation_watts(self, time) :
    return self.__call__(time, "RadiationGlobal_Avg10min [W*m-2]")
  #enddef
  def radiation_lux(self, time) :
    irradiance = self.__call__(time, "RadiationGlobal_Avg10min [W*m-2]")
    Eqf = lambda l, I: 10 * constants.pi * I * l * constants.nano / (constants.h * constants.c * constants.Avogadro * constants.micro) # [umol/m2s = muE]
    # conversion to UE Lux
    lux = Eqf(555, irradiance) # [lux]
    lux = max(1e-7, lux)
    return lux
  #enddef
  def radiation_par(self, time) :
    photo_radiation = self.__call__(time, "RadiationPhotosyntheticActive_2m_Avg10min [umol*m-2*s-1]")
    photo_radiation = max(1e-7, photo_radiation)
    # umol to mol
    return photo_radiation / 1e6
  def precipitation(self, time) :
    return self.__call__(time, "Precipitation_Avg10min [mm]")
  #enddef
  def temperature(self, time) :
    return self.__call__(time, "AirTemperature_2m_Avg30min [°C]")
  #enddef
  def temperature_k(self, time) :
    return self.__call__(time, "AirTemperature_2m_Avg30min [°C]") + 273.15
  #enddef
  def saturated_vapour_pressure(self, time) :
    temperature = self.__call__(time, "AirTemperature_2m_Avg30min [°C]")
    e0 = 6.1078
    L = 2.5 * 10 ** 6
    R0 = 461.5
    return e0 * np.exp(L/R0 * (1.0/273.15 - 1.0/(temperature + 273.15)))
  #enddef
  def actual_vapour_pressure(self, time) :
    relative_humidity = self.relative_humidity(time)
    saturated_vapour_pressure = self.saturated_vapour_pressure(time)
    return relative_humidity * saturated_vapour_pressure
  #enddef
  def mean_metric_potential(self, time) :
    soil_moisture = self.__call__(time, "SoilMoisture_10cm_Avg10min [m3*m-3]")
    # water content 0% -> -10kPa, 100% -> -1500kPa
    return -10.0 + soil_moisture * (-1500.0 + 10.0)
  #enddef
  def molar_fraction_co2(self, time) :
    co2_molar = self.__call__(time, "AirConcentration_CO2_2m_Avg30min [mmol*m-3]")
    air_pressure = self.air_pressure(time)
    air_molar = (air_pressure * 100.0) / (8.314 * self.temperature_k(time)) * 1000.0
    return co2_molar / air_molar
  #enddef
  def get_time(self) :
    return self.data.index
  #enddef
  def get_columns(self) :
    return self.data.columns
  #enddef
  def rbl(self, time) :
    leaf_thickness = 0.0001 # m
    diffusivity = 2.5 * 10 ** -5
    return leaf_thickness / diffusivity
  #enddef
  def rcanopy(self, time) :
    # resistivity to water vapour flow in the canopy
    eddy_covariance_canopy = 0.1 # m2 s-1
    return 1.0 / eddy_covariance_canopy
  #enddef
  def generator_t(self, property, start_t, dt, end_time) :
    t = start_t
    while t < end_time :
      p = self.__call__(t, property)
      yield t, p
      t += dt
    #endwhile
  #enddef
#endclass

soil = Soil()

# PREAMBLE FOR SYNAVIS+CPLANTBOX COUPLING ###############################

# if the runtime folder is not the python script folder, switch there
if os.path.dirname(os.path.abspath(__file__)) != os.getcwd():
  os.chdir(os.path.dirname(os.path.abspath(__file__)))

# # add Synavis (some parent directory to this) to path
# sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "unix"))
# import PySynavis as syn
# syn.SetGlobalLogVerbosity(syn.LogVerbosity.LogError)
# pylog = syn.Logger()
# pylog.setidentity("Chamber Experiment")
# # make a stand-in for print
# def logfun(*args) :
  # pylog.log(" ".join([str(a) for a in args]))
# #enddef

# sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "modules"))
# import signalling_server as ss

cplantbox_dir = "../../"

# # check if we have a CPLANTBOX_PATH in the environment
# if "CPLANTBOX_PATH" in os.environ:
  # sys.path.append(os.environ["CPLANTBOX_PATH"])
  # sys.path.append(os.environ["CPLANTBOX_PATH"] + "/src")
  # cplantbox_dir = os.environ["CPLANTBOX_PATH"]
# else:
  # # request input for the path to the cplantbox directory
  # print("Please select the path to the cplantbox directory")
  # folder_selected = input("Path: ")
  # if not os.path.isdir(folder_selected):
    # print("The path does not exist")
    # sys.exit(1)
  # # append the path to the sys path
  # sys.path.append(folder_selected)
  # # write the export line to .bashrc
  # with open(os.path.expanduser("~/.bashrc"), "a") as f:
    # f.write("export CPLANTBOX_PATH=" + folder_selected + "\n")
  # cplantbox_dir = folder_selected
  # print("CPLANTBOX_PATH exported to ~/.bashrc")
  # print("Please restart the terminal to update the environment variables")
  # print("The path has been added to the sys path for this session")
# #endif

import plantbox as pb
from functional.xylem_flux import XylemFluxPython
from functional.photosynthesis_cpp import PhotosynthesisPython as Photosynthesis


# model parameters
simtime = (1 / 24) * 2 # days
depth = 40 # cm
# dt = 0.1 hours to days
dt = 0.2 / 24.0
spacing = 25


class Signal :
  def __init__(self, args, callback) :
    self.args = args
    self.callback = callback
  #enddef
  def __call__(self, msg) :
    self.callback(msg)
  #enddef
  # equal operator
  def __eq__(self, args:dict) :
    return all(self.args[key] == args[key] for key in args if key in self.args)
    # any formulation checks if the key is in the args and if the value is different
    #return not any(self.args[key] != args[key] if key in self.args else True for key in args)
  #enddef
  def __str__(self) :
    return str(self.args)
  #enddef
#endclass

message_buffer = []
signal_relay = []

def reset_message() :
  global message_buffer
  message_buffer = []
#enddef

def get_message() :
  global message_buffer
  while len(message_buffer) == 0 :
    time.sleep(0)
  #endwhile
  message = message_buffer.pop(0)
  return message
#enddef

def message_callback(msg) :
  global message_buffer, signal_relay, pylog
  try :
    msg = json.loads(msg)
    if "type" in msg :
      pylog.log("Message received with parameters " + ", ".join([key + ": " + str(msg[key])[0:10] for key in msg]))
      # ignore certain types
      if msg["type"] == "parameter" :
        return
      handler = next((h for h in signal_relay if h == msg), None)
      logfun("Handler is ", str(handler), " and not any of ", ",".join([str(h) for h in signal_relay]))
      if not handler is None :
        handler(msg)
      else :
        pylog.log("Appended to message buffer (now " + str(len(message_buffer)) + " messages)")
        message_buffer.append(msg)
      #endif
    #endif
  except :
    pylog.log("Skipping message of size " + str(len(msg)) + " for failing to convert to JSON")
  #endtry
#enddef

# END OF PREAMBLE ########################################################

hours_dt = 0.5

start_time = datetime.datetime(2016, 4, 20, 12, 13, 48, 0)
#weather = Weather(55.7, 13.2, start_time)
end_time = start_time + datetime.timedelta(days = simtime)
timerange_numeric = np.arange(0, simtime, dt)
delta_time = datetime.timedelta(hours = hours_dt)
timerange = [start_time]
while timerange[-1] < end_time:
  timerange.append(timerange[-1] + delta_time)

parameter_file = os.path.join(cplantbox_dir, "modelparameter", "structural", "plant", "Triticum_aestivum_adapted_2023.xml")

# SETUP SYNAVIS ##########################################################

# dataconnector = syn.DataConnector()
# dataconnector.Initialize()
# dataconnector.SetMessageCallback(message_callback)
# dataconnector.SetConfig({
    # "SignallingIP": "172.20.16.1", #ss.get_interface_ip("ib0"),
    # "SignallingPort": 8080
  # })
# dataconnector.SetTakeFirstStep(True)
# dataconnector.StartSignalling()
# dataconnector.SetRetryOnErrorResponse(True)
# dataconnector.LockUntilConnected(500)
# dataconnector.SendJSON({"type":"delete"})


# dataconnector.SendJSON({"type":"command", "name":"cam", "camera": "scene"})
# dataconnector.SendJSON({"type":"console", "command":"t.maxFPS 18"})

# dataconnector.SendJSON({"type":"spawnmeter", "number": 20, "calibrate": False})
# m_ = {"type":"placeplant", 
                        # "number": 1,
                        # "rule": "square",
                        # "mpi_world_size": 1,
                        # "mpi_rank": 0,
                        # "spacing": 10,
                      # }
# dataconnector.SendJSON(m_)

# SETUP CPLANTBOX ########################################################

# logfun("Setting up plantbox")

# plot three
# Sowing: 2015-10-26
sowing_time = datetime.datetime(2015, 10, 26, 6, 0, 0, 1)
# Emergence: 2015-11-01
emergence_time = datetime.datetime(2015, 2, 1, 6, 0, 0, 1)#(2015, 11, 1, 6, 0, 0, 1)
# Tasseling: -
# Flowering: 2016-06-03
flowering_time = datetime.datetime(2016, 3, 1, 6, 0, 0, 1)

# time from sowing until start
time_from_sowing = start_time - sowing_time

plant_relative_scaling = 10.0

plant = pb.MappedPlant()
#seednum = 1
plant.readParameters(parameter_file)
#sp = plant.getOrganRandomParameter(1)[0]
#sp.seedPos = pb.Vector3d(0,0,0)
sdf = pb.SDF_PlantBox(np.inf, np.inf, 40)
plant.setGeometry(sdf)
plant.initialize(False, True)
r = Photosynthesis(plant, -500, 360/2)
r.setKr([[1.728e-4], [0], [3.83e-5]])  
r.setKx([[4.32e-1]])
# copy of parameters used in dumux-CPB photosynthesis study (C3, pseudo-wheat)
r.g0 = 8e-6
r.VcmaxrefChl1 =1.28
r.VcmaxrefChl2 = 8.33
SPAD= 41.0
chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
r.Chl = np.array( [chl_]) 
r.a1 = 0.6/0.4
r.a3 = 1.5
r.alpha = 0.4 
r.theta = 0.6 
r.fwr = 1e-16
r.fw_cutoff = 0.04072 
r.gm =0.03
r.sh =5e-4
#r.maxErrLim = 1/100
r.maxLoop = 10000
r.minLoop=900

#logfun("Pre-simulating until chamber measurement")
# pre-simulating until the day of the gas chamber measurement
delta_time_days = time_from_sowing.days
#delta_time_days = 30
plant.simulate(delta_time_days,False)

has_measured = False
measurement = np.array([])
def ReceiveMeasurement(msg) :
  global has_measured, measurement
  if "i" in msg and msg["type"] == "mm" :
    #logfun("Received measurement of light from UE")
    has_measured = True
    mlight = base64.b64decode(msg["i"])
    mlight = np.frombuffer(mlight, dtype=np.float32)
    light = np.nan_to_num(mlight)
    measurement = light
  #endif
#enddef

#logfun("Starting simulation")

# start the measurement
for t in timerange:
  plant.simulate(hours_dt, False)
  vis = pb.PlantVisualiser(plant)
  vis.SetLeafResolution(20)
  vis.SetGeometryResolution(6)
  vis.SetComputeMidlineInLeaf(False)
  #vis.SetLeafMinimumWidth(1.0)
  #vis.SetRightPenalty(0.5)
  #vis.SetVerbose(False)
  time_string_ue = t.strftime("%Y.%m.%d-%H:%M:%S")
  # get organs
  organs = plant.getOrgans()
  slot = 0
  time.sleep(0)
  leaf_points = 0
  leaf_amount = 0
  # send all stem or leaf organs
  for organ in organs :
    if organ.organType() != 2 :
      vis.ResetGeometry()
      vis.ComputeGeometryForOrgan(organ.getId())
      points,triangles,normals = np.array(vis.GetGeometry())*plant_relative_scaling, np.array(vis.GetGeometryIndices()), np.array(vis.GetGeometryNormals())
      if organ.organType() == pb.leaf.value :
        leaf_points += len(points)
        leaf_amount += 1
      #endif
      # texture = np.array(vis.GetGeometryTextureCoordinates())
      # message = {"type":"do",
      # "p": base64.b64encode(points.astype("float64")).decode("utf-8"),
      # "n": base64.b64encode(normals.astype("float64")).decode("utf-8"),
      # "i": base64.b64encode(triangles.astype("int32")).decode("utf-8"),
      # "t": base64.b64encode(texture.astype("float64")).decode("utf-8"),
      # "l": 0,
      # "s": slot,
      # "o": int(organ.organType()),
      # }
      # dataconnector.SendJSON(message)
      # time.sleep(0.05)
      slot += 1
    #endif
  #endfor
  # #logfun("Sent " + str(leaf_amount) + " leaf organs with " + str(leaf_points) + " points")
  # surrounding variables
  intensity = 120000.0      # from callibration (UE! CHECK WITH EXPERIMENT!)
  radiation_par = 1549.310181  # from experiment!  micromole m-2 s-1
  leaf_nodes = r.get_nodes_index(4)[:-1]
  if len(leaf_nodes) == 0 :
    #logfun("No leaf nodes found")
    exit(-1)
  plant_nodes = np.array(plant.getNodes())
  # dataconnector.SendJSON({"type":"resetlights"})
  # dataconnector.SendJSON({"type":"t", "s":time_string_ue})
  # #pylog.logjson({"type":"t", "s":time_string_ue})
  # dataconnector.SendJSON({
    # "type": "parameter",
    # "object": "SunSky_C_1.DirectionalLight",
    # "property": "Intensity",
    # "value": intensity
  # })
  # time.sleep(0.2)
  # dataconnector.SendJSON({
    # "type":"calibrate",
    # "flux": radiation_par
  # })
  leaf_nodes = plant_nodes[leaf_nodes]*plant_relative_scaling
  # message = {"type":"mms", "l":0, "p": base64.b64encode(leaf_nodes.tobytes()).decode("utf-8")}
  # signal_relay.append(Signal({"type":"mm"}, ReceiveMeasurement))
  # logfun("Registered handler for measurement")
  # dataconnector.SendJSON(message)
  # logfun("Sent measurement request, waiting now...")
  #while not has_measured :
  #  time.sleep(0)
  has_measured = False
  # logfun("Received measurement, starting photosynthesis")
  # simulation data
  p_s_input = -200. #soil.soil_matric_potential(0.1) # soil matric potential, lookup
  print('p_s_input',p_s_input)
  #RH_input = weather.relative_humidity(t) # relative humidity, lookup
  RH_input = 0.85 # from average historic data for april
  Tair_input = 19.19866943 # air temperature, lookup
  P_s_input = 101.7164764 # air pressure, lookup Pa
  cs_input = 0.0 # CO2 molar fraction, lookup
  es = 0.61078 * np.exp(17.27 * Tair_input / (Tair_input + 237.3))  # FAO56
  ea = es * RH_input
  Patm = 1070.00  # hPa 101.3 * ((293.0 - 0.0065 * depth) / 293.0) ** 5.26 # lookup air pressure at time
  r.Qlight = intensity
  r.solve_photosynthesis(sim_time_ = hours_dt, sxx_= [p_s_input], cells_= True, ea_ = ea, 
  es_ = es, verbose_ = False, doLog_ = False, TairC_ = Tair_input, outputDir_="./")
  fluxes = r.outputFlux # cm3/day
  idsC3 = [np.where(np.array(r.ci)> 0)[0]]
  print('idsC3',idsC3)
  print(r.fw)#np.min(np.array(r.fw[idsC3])),np.max(np.array(r.fw[idsC3])))
  print('rx',r.rx)
  raise Exception
