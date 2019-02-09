import os 

path = os.getcwd()

for lmax in range(20,51,2):
  os.chdir(path+'/nmax500lmax'+str(lmax))
  os.system('pwd') 
  os.system('sbatch submit_sim')
