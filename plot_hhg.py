import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np 
file_in = open("H_dipoleAcceleration.output",'r')

array = []
for line in file_in:
  if ("Vec" not in line) and ("type" not in line) and ("Process" not in line): 
    row = line.split()
    array = array + [row[0]]
#end for 

dipole_acceleration = []
for num in array:
  dipole_acceleration = dipole_acceleration + [float(num)]
#end for 

print "Plotting HHG Spectrum"
fig = plt.figure(figsize=(24, 18), dpi=80)
energy =   0.0570
data = np.array(dipole_acceleration)
data = data
data = data * np.blackman(data.shape[0])
padd2 = 2**np.ceil(np.log2(data.shape[0] * 4))
paddT = 2.58E+02*padd2 / data.shape[0]
dH = 2 * np.pi / paddT / energy
if np.max(data) > 1e-19:
  data = np.absolute(
    np.fft.fft(
      np.lib.pad(
        data, (int(np.floor((padd2 - data.shape[0]) / 2)),
          int(np.ceil((padd2 - data.shape[0]) / 2))),
        'constant',
        constant_values=(0.0, 0.0))))
data /= data.max()
data = data**2.0
plt.semilogy(
  np.arange(data.shape[0]) * dH,
  data)
plt.ylabel("HHG Spectrum (a.u.)")
plt.title("HHG Spectrum")
plt.legend()
x_min = 0
x_max = 50
plt.xticks(np.arange(x_min + 1, x_max + 1, 2.0))
plt.xlim([x_min, x_max])
plt.ylim([1e-16, 1])
plt.grid(True, which='both')
plt.tight_layout()
fig.savefig(os.getcwd()+'/HHG_Spectrum.png')
plt.clf()
plt.close(fig)

