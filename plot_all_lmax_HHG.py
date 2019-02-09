import os
import shutil
# create folders for each lmax 
home = os.getcwd()
for lmax in range(20,51,2):
  path = os.getcwd() + "/nmax500lmax" + str(lmax)

  file_name = 'plot_hhg.py'

  try:
    shutil.copyfile(os.getcwd()+"/" + file_name,path +"/" + file_name)
  except OSError:
    print ("Creation of %s failed" % file_name)
  else: 
    print ("Successfully created %s" % file_name)
  #end try
  os.chdir(path)
  os.system('pwd') 
  os.system('python plot_hhg.py')
  os.chdir(home)
#end for 