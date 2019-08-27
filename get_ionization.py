import io
import os
from glob import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def find_between_r( s, first, last ):
  try:
    start = s.rindex( first ) + len( first )
    end = s.rindex( last, start )
    return s[start:end]
  except ValueError:
    return ""

root = os.getcwd()

x = np.zeros(3001)
y = np.zeros(3001)

for i in range(0,3001):
  print i
  os.chdir(str(i))
  if glob("H_rho.output") != []:
    file_name = glob('slurm*')
    if file_name == []:
      file_name = glob('run*')
    file_in = open(file_name[0])
    for line in file_in:
      last = line
    x[i] = i - 1500
    y[i] = find_between_r(last,',',')')
  else:
    x[i] = i - 1500
    y[i] = 0.0
  os.chdir(root)

np.savetxt("relative_delay.txt", x, fmt="%s")
np.savetxt("ionization.txt",y, fmt="%s")

fig = plt.figure(figsize=(24, 18), dpi=80)
plt.plot(x,y)

plt.grid(True, which='both')
plt.tight_layout()
fig.savefig(os.getcwd()+'/ionization.png')
plt.clf()
plt.close(fig)

energy = 0.1875
fig = plt.figure(figsize=(24, 18), dpi=80)
data = y
data = data * np.blackman(data.shape[0])
padd2 = 2**np.ceil(np.log2(data.shape[0] * 4))
paddT = 3000*padd2 / data.shape[0]
dH = 2 * np.pi / paddT 
if np.max(data) > 1e-19:
  data = np.absolute(
    np.fft.fft(
      np.lib.pad(
        data, (int(np.floor((padd2 - data.shape[0]) / 2)),
          int(np.ceil((padd2 - data.shape[0]) / 2))),
        'constant',
        constant_values=(0.0, 0.0))))
data = data**2.0
plt.semilogy(
  np.arange(data.shape[0]) * dH,
  data)
plt.legend()
x_min = 0
x_max = 1.5
plt.xlim([x_min, x_max])

xticks = np.zeros(30)
for i in range(0,10):
  xticks[0+3*i] = i*0.1875
  xticks[1+3*i] = i*0.1875 + .0625
  xticks[2+3*i] = (i+1)*0.1875-.0625

plt.xticks(xticks)
plt.ylim([1e-12*data.max(), 1*data.max()])
plt.grid(True, which='both')
plt.tight_layout()
fig.savefig(os.getcwd()+'/fft.png')
plt.clf()
plt.close(fig)
