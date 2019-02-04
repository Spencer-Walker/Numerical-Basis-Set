import os
import shutil

# create folders for each lmax 
for lmax in range(20,51,2):
  path = os.getcwd() + "/nmax500lmax" + str(lmax)

  try:
    os.mkdir(path)
  except OSError:
    print ("Creation of the directory %s failed" % path)
  else:
    print ("Successfully created the directory %s" % path)
  
  # create all of the needed files
  for file_name in ["makefile", "schrodinger1Df90.F90",\
    "simulation_parametersf90.F90","generateH0f90.F90",\
    "generateH0f90.F90","generateDipoleAccelerationf90.F90",\
    "propagatef90.F90"]:
    
    try:
      shutil.copyfile(os.getcwd()+"/" + file_name,path +"/" + file_name)
    except OSError:
      print ("Creation of %s failed" % file_name)
    else: 
      print ("Successfully created %s" % file_name)
