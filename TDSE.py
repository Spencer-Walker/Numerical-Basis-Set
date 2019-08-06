import os 
import json
import h5py
import numpy as np

def potential(c0,zc,a,b,c,r,cap_eta,cap_present,gobbler):
  V =  np.zeros((2,len(r)))
  V[0,:] += -c0/r
  V[0,:] += -zc*np.exp(-c*r)/r
  for i in range(len(a)):
    V[0,:] += -a[i]*np.exp(-b[i]*r)
  if cap_present == 1 :
    iabs = int(round(gobbler*Rmax/dr)) - 1
    V[1,iabs:] += -cap_eta*np.sin(np.pi*(r[iabs:]/Rmax-gobbler)/(2.0*(1.0-gobbler)))**2.0
  return V

def electric_field(lam,I,shape,frequency_shift,t,num_cycles,t_start,custom_envalope_phase,cep):
  if not isinstance(t,float):
    E =  np.zeros(len(t))
  conEV = 27.2110
  w = 1239.8006/(lam*conEV)
  conWpCM2 = 3.5101e16
  E0 = np.sqrt(I/conWpCM2)
  if shape == 'sin2': 
    if frequency_shift == 1:
      mu = 4.0*np.arcsin(np.exp(-0.25))**2.0
    else:
      mu = 1.0
    wa = 2.0*w/(1.0 + np.sqrt(1.0 + mu/num_cycles**2.0))
    T0 = num_cycles*2.0*np.pi/wa
    t_cep = T0/2.0 + t_start
    if custom_envalope_phase == 0:
      t_cep = T0/2.0 + t_start
      cep  = 0.0
    t_end = t_start + T0
    if not isinstance(t,float):
      for i in range(len(t)):
        if t_start <= t[i] < t_end:
          E[i] =  -((E0*np.sin((np.pi*(t[i]-t_start))/T0)* \
            (T0*wa*np.cos(cep + (t[i] - t_cep)*wa)* \
              np.sin((np.pi*(t[i] - t_start))/T0) + \
                2.0*np.pi*np.cos((np.pi*(t[i] - t_start))/T0)* \
                  np.sin(cep + (t[i] - t_cep)*wa)))/(T0*wa))
        else:
          E[i] = 0.0
    else: 
      if t_start <= t < t_end:
        E =  -((E0*np.sin((np.pi*(t-t_start))/T0)* \
          (T0*wa*np.cos(cep + (t - t_cep)*wa)* \
            np.sin((np.pi*(t - t_start))/T0) + \
              2.0*np.pi*np.cos((np.pi*(t - t_start))/T0)* \
                np.sin(cep + (t - t_cep)*wa)))/(T0*wa))
      else:
        E = 0.0
  elif shape == 'gaussian': 
    if frequency_shift == 1:
      mu = 8.0*np.log(2.0)/np.pi**2.0
    else:
      mu = 1.0
    wa = 2.0*w/(1.0 + np.sqrt(1.0 + mu/num_cycles**2.0))
    T0 = num_cycles*2.0*np.pi/wa
    t_cep = T0/2.0 + t_start
    if custom_envalope_phase == 0:
      t_cep = 5.0*T0 + t_start
      cep  = 0.0
    t_end = t_start + 10.0*T0
    if not isinstance(t,float):
      for i in range(len(t)):
        if t_start <= t[i] < t_end:
          E[i] = -((E0*(T0**2.0*wa*np.cos(cep + (t[i] - t_cep )*wa) + 8.0* \
            (-t[i] + t_cep )*np.log(2.0)*np.sin(cep +  (t[i] - t_cep )*wa)))/ \
              (2.0**((4.0*(-t[i] + t_cep)**2.0)/T0**2.0)*T0**2.0*wa))
        else:
          E[i] = 0.0
    else:
      if t_start <= t < t_end:
        E = -((E0*(T0**2.0*wa*np.cos(cep + (t - t_cep )*wa) + 8.0* \
          (-t + t_cep )*np.log(2.0)*np.sin(cep +  (t - t_cep )*wa)))/ \
            (2.0**((4.0*(-t + t_cep)**2.0)/T0**2.0)*T0**2.0*wa))
      else:
        E = 0.0
  return E

# Recursive function that traverses the json file and 
# converts it to hdf5
def json_to_hdf5(key,val,group):
  if isinstance(val,(dict)):
    grp = group.create_group(str(key))
    data = val
    for key,val in data.items():
      json_to_hdf5(key,val,grp)
  elif isinstance(val,unicode):
    dset = group.create_dataset(str(key),(1,),dtype = "S"+str(len(val)))
    dset[0] =  np.string_(val) 
  elif isinstance(val,int):
    dset = group.create_dataset(str(key),(1,),dtype = "i4")
    dset[0] =  val 
  elif isinstance(val,float):
    dset = group.create_dataset(str(key),(1,),dtype = "f8")
    dset[0] =  val
  elif isinstance(val,list):
    if all(isinstance(x,int) for x in val):
      dset = group.create_dataset(str(key),(len(val),),dtype = "i4")
      dset[:] =  val[:] 
    elif all(isinstance(x,float) for x in val):
      dset = group.create_dataset(str(key),(len(val),),dtype = "f8")
      dset[:] =  val[:]      
    elif isinstance(val[0],dict):
      for i in range(len(val)):
        grp = group.create_group(str(key)+str(i))
        for k,v in val[i].items():
          json_to_hdf5(k,v,grp)
    else:
      print("Innaproprate datatype present in dict")
  else:
    print("Innaproprate datatype present in input.json")

# Gets the current working directory
cwd = os.getcwd()
# If an old version of the parameters.h5 file exits remove it
# since we are are reading in a new input.json
os.system('rm parameters.h5')
# Create a new parameters.h5 file to write to
params = h5py.File('parameters.h5','w')

# Read the json file and write it to hdf5.
with open(cwd+"/input.json", 'r') as f:
  data = json.load(f)
  for key,val in data.items():
    json_to_hdf5(key,val,params)
Rmax = data["EPS"]["R_max"]
dr = data["EPS"]["delta_x"]
num_points = int(round(Rmax/dr))
r = np.zeros((num_points,))
for i in range(num_points):
  r[i] = dr*(i+1)
a = data["EPS"]["nuclei"]["a"]
b = data["EPS"]["nuclei"]["b"]
c = data["EPS"]["nuclei"]["c"]
c0 = data["EPS"]["nuclei"]["c0"]
zc = data["EPS"]["nuclei"]["Zc"]
cap_present = data["EPS"]["cap_present"]
cap_eta = data["EPS"]["cap_eta"]
gobbler = data["EPS"]["gobbler"]
V = potential(c0,zc,a,b,c,r,cap_eta,cap_present,gobbler)
dset = params["EPS"].create_dataset("V",(2,num_points),dtype = "f8")
dset[:,:] = V[:,:]
dset = params["EPS"].create_dataset("r",(num_points,),dtype = "f8")
dset[:] = r[:]
conEV = 27.2110
dt = np.zeros((len(data["laser"]["pulse"]),))
L = np.zeros((len(data["laser"]["pulse"]),))
t_start = np.zeros((len(data["laser"]["pulse"]),))
t_end = np.zeros((len(data["laser"]["pulse"]),))
for i in range(len(data["laser"]["pulse"])):
  if data["laser"]["pulse"][i]["pulse_shape"] == 'sin2': 
    if data["laser"]["frequency_shift"] == 1:
      mu = 4.0*np.arcsin(np.exp(-0.25))**2.0
    else:
      mu = 1.0
  elif data["laser"]["pulse"][i]["pulse_shape"] == 'gaussian': 
    if data["laser"]["frequency_shift"] == 1:
      mu = 8.0*np.log(2.0)/np.pi**2.0
    else:
      mu = 1.0

  w = 1239.8006/(data["laser"]["pulse"][i]["wavelength"]*conEV)
  wa = 2.0*w/(1.0 + np.sqrt(1.0 + mu/data["laser"]["pulse"][i]["cycles"]**2.0))
  T0 = data["laser"]["pulse"][i]["cycles"]*2.0*np.pi/wa
  dt[i] = data["laser"]["pulse"][i]["relative_delay"]
  if data["laser"]["pulse"][i]["pulse_shape"] == 'sin2':
    L[i] = T0
  elif data["laser"]["pulse"][i]["pulse_shape"] == 'gaussian':
    L[i] = 2.0*data["laser"]["pulse"][i]["gaussian_length"]*T0
  if i == 0:
    Lfirst = L[0]
    Llast  = L[0]
    dtfirst = dt[0]
    dtlast = dt[0]
  for j in range(i+1):
    if L[i]/2.0-dt[i] > L[j]/2.0-dt[j]:
      Lfirst = L[i]
      dtfirst = dt[i]
    if L[i]/2.0+dt[i] > L[j]/2.0+dt[j]:
      Llast = L[i]
      dtlast = dt[i]

max_time = 0.5*(Lfirst+Llast) + dtlast - dtfirst 
for i in range(len(data["laser"]["pulse"])):
  t_start[i] = (Lfirst/2.0 - dtfirst) - (L[i]/2.0-dt[i])
  t_end[i] = t_start[i] + L[i]
num_steps = int(np.ceil(max_time/data["TDSE"]["delta_t"]))+1
t = np.zeros((num_steps,))
dset = params["TDSE"].create_dataset("t",(num_steps,),dtype="f8")
for i in range(num_steps):
  t[i] = i*data["TDSE"]["delta_t"]
dset[:] = t[:]
E = np.zeros(t.shape)
for i in range(len(data["laser"]["pulse"])):
  lam = data["laser"]["pulse"][i]["wavelength"]
  I = data["laser"]["pulse"][i]["intensity"]
  shape = data["laser"]["pulse"][i]["pulse_shape"]
  frequency_shift = data["laser"]["frequency_shift"]
  num_cycles = data["laser"]["pulse"][i]["cycles"]
  custom_envalope_phase = data["laser"]["pulse"][i]["custom_envalope_phase"]
  cep = data["laser"]["pulse"][i]["cep"]
  Ei = electric_field(lam,I,shape,frequency_shift,t,num_cycles,t_start[i],custom_envalope_phase,cep)
  dset = params["laser"]["pulse"+str(i)].create_dataset("E"+str(i),(num_steps,),dtype = "f8")
  dset[:] = Ei[:]
  E += Ei
dset = params["laser"].create_dataset("E",(num_steps,),dtype = "f8")
dset[:] = E[:]
dset = params["laser"].create_dataset("E_length",(1,),dtype = "i4")
dset[0] = num_steps
params.close()


# Compute basis
if data["EPS"]["compute"] == 1 :
  if data["mpi"]["assign_host"] == 1:
    print("basisf90")
    os.system("mpirun -np " + str(data["mpi"]["np"]) + \
      " --host " + str(data["mpi"]["host"]) + " ./basisf90")
  else:
    print("basisf90")
    os.system("mpirun -np " + str(data["mpi"]["np"]) + \
      " ./basisf90")

# Compute operators  
if data["operators"]["compute"] == 1:
  # If manually assigning host (done when running on multiple nodes)
  if data["mpi"]["assign_host"] == 1:
    print("generateH0f90")
    os.system("mpirun -np " + str(min(data["mpi"]["np"],data["TDSE"]["l_max"])) + \
      " --host " + str(data["mpi"]["host"]) + \
        " ./generateH0f90")
    if data["TDSE"]["l_max"] > 0:
      print("generateH1f90")
      os.system("mpirun -np " + str(min(data["mpi"]["np"],data["TDSE"]["l_max"]))+ \
        " --host " + str(data["mpi"]["host"]) + \
          " ./generateH1f90")
      print("generateDipoleAccelerationf90")
      os.system("mpirun -np " + str(min(data["mpi"]["np"],data["TDSE"]["l_max"])) + \
        " --host " + str(data["mpi"]["host"]) + \
          " ./generateDipoleAccelerationf90")
  else:
    print("generateH0f90")
    os.system("mpirun -np " + str(min(data["mpi"]["np"],data["TDSE"]["l_max"])) + \
      " ./generateH0f90")
    if data["TDSE"]["l_max"] > 0:
      print("generateH1f90")
      os.system("mpirun -np " + str(min(data["mpi"]["np"],data["TDSE"]["l_max"])) + \
        " ./generateH1f90")
      print("generateDipoleAccelerationf90")
      os.system("mpirun -np " + str(min(data["mpi"]["np"],data["TDSE"]["l_max"])) + \
        " ./generateDipoleAccelerationf90")

# Propagate 
if data["TDSE"]["propagate"] == 1:
  # If running on multiple nodes 
  if data["mpi"]["assign_host"] == 1:
    print("propagatef90")
    os.system("mpirun -np " + str(data["mpi"]["np"]) + \
      " --host " + str(data["mpi"]["host"]) + " ./propagatef90")
  else:
    print("propagatef90")
    os.system("mpirun -np " + str(data["mpi"]["np"]) + \
      " ./propagatef90")
