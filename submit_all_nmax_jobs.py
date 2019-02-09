import os 

path = os.getcwd()

for nmax in range(400,500,10):
  os.chdir(path+"/nmax"+str(nmax) + "lmax50")
  os.system('pwd') 
  os.system('sbatch submit_sim')
