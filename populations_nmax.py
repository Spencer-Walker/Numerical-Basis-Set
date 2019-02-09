import os
path = os.getcwd() 
pop2 = []
pop4 = []
pop20 = []
for nmax in range(400,500,10):
  os.chdir(path+ "/nmax"+str(nmax) + "lmax50")
  file_in = open(os.getcwd()+"/H_rho.output",'r')
  data = []
  for line in file_in:
    if not(('Vec' in line) or ('type' in line) or ('Process' in line)):
      data = data + [float(line)]
    #end if 
  #end for 


  low_index = 0

  pop = []
  print('nmax' + str(nmax))
  for i in range(1,21):
    pop = pop + [[i,sum(data[low_index:low_index+i])]]
    low_index = low_index + i
  #end for 

  pop2 = pop2 + [pop[1][1]]
  pop4 = pop4 + [pop[3][1]]
  pop20 = pop20 + [pop[19][1]]
print('pop2')
print(pop2)
print('pop4')
print(pop4)
print('pop20')
print(pop20)

