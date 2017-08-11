import numpy as np
import matplotlib.pyplot as plt
import stray_field_calc as sfc
import format_plot as fp

height = 50
Ms = 1
t = 1
Mst = Ms*t
pi = np.pi

def dipole_field(m, h, res, size):
	dlen = int(size/res)
	bex = np.zeros((3, dlen, dlen))
	r = np.array([0, 0, h])
	for j in range(0, dlen):
		r[1] = res*j-(size/2)
		for i in range(0, dlen):
			r[0] = res*i-(size/2)
			for k in range(0, 3):
				bex[k][j, i] = (1e-3) * ( (3*r[k]*(m[0]*r[0]+m[1]*r[1]+m[2]*r[2]) / ((r[0]**2 + r[1]**2 + h**2)**(5/2))) -
					(m[k] / ((r[0]**2 + r[1]**2 + h**2)**(3/2))) )
	return bex

baspath = "/Users/alec/UCSB/LTSPM/cofeb_analysis/ta/stray_field_sim/"
filenames = ["mbx","mby","mbz"]
numfiles = len(filenames)
m = []
for i in range(0,numfiles):
    m.append(np.transpose(np.loadtxt(baspath+filenames[i]+".dat")))

scd, vcd, meff, hk, b = sfc.stray_field_calc(m[0],m[1],m[2],Mst,2500,height)

flen = len(b[0])
fres = 2500/flen
fx = np.arange(-flen*fres/2, flen*fres/2, fres)

md = [Ms * (fres**3), 0, 0]

bexact = dipole_field(md, height, fres, 2500)

plt.close('all')

fig1, ax1 = plt.subplots()
plt.plot(fx, bexact[0][int(flen/2), :], 'r-')
plt.plot(fx, b[0][int(flen/2), :], 'b-')
fp.format_plot(plt, 650, 400, 0, 50)

fig1, ax1 = plt.subplots()
plt.plot(fx, bexact[1][int(flen/2), :], 'r-')
plt.plot(fx, b[1][int(flen/2), :], 'b-')
fp.format_plot(plt, 650, 400, 0, 450)

# fig1, ax1 = plt.subplots()
# plt.imshow(bexact[1])
# plt.imshow(b[1])
# fp.format_plot(plt, 650, 400, 0, 450)

plt.show()