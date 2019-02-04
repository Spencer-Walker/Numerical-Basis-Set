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
  #end try

  # create all of the needed files
  for file_name in ["makefile", "schrodinger1Df90.F90",\
    "generateH0f90.F90","generateH0f90.F90",\
    "generateDipoleAccelerationf90.F90",\
    "propagatef90.F90","submit_sim"]:
    
    try:
      shutil.copyfile(os.getcwd()+"/" + file_name,path +"/" + file_name)
    except OSError:
      print ("Creation of %s failed" % file_name)
    else: 
      print ("Successfully created %s" % file_name)
    #end try

    file_in = open(os.getcwd()+"/simulation_parametersf90.F90",'r')
    file_out = open(path+"/simulation_parametersf90.F90",'w')

    for line in file_in:
      if "l_max" in line:
        file_out.write("  integer,   parameter :: l_max = "+str(lmax)+"\n")
      else:
        file_out.write(line)
      #end if 
    #end for 
  #end for 
#end for 