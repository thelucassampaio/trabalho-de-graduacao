import requests
import argparse, sys
from galpy.orbit import Orbit
from galpy.potential import plotRotcurve
from galpy.potential import SpiralArmsPotential
from galpy.potential.JunqueiraPotential import JunqueiraPotential
####### Potencial do Junqueira está em

####### C:\Users\lucas\Anaconda3\Lib\site-packages\galpy\potential

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

from galpy.potential import MWPotential2014
from galpy.potential import MiyamotoNagaiPotential
from galpy.potential import LogarithmicHaloPotential

# ep = evaluatePotentials(MiyamotoNagaiPotential, 1., 0.)
# print(SpiralArmsPotential)
# print(JunqueiraPotential)

sap = SpiralArmsPotential(N = 4, omega = 25)
lhp = LogarithmicHaloPotential()
# plotRotcurve(lhp + sap, Rrange=[0.01,10.],grid=1001,yrange=[0.,1.2], phi = 0.)
mp = MiyamotoNagaiPotential(a=0.5,b=0.0375,normalize=1.)
jp = JunqueiraPotential()

# plotRotcurve(lhp, Rrange=[0.01,10.],grid=1001,yrange=[0.,0.4], phi = 0.1)
# plotRotcurve(jp, Rrange=[0.01,10.], yrange=[0.,3.], phi = 0.1)
# plt.show()


raios = np.linspace(0.1, 12., 100)
phis = np.linspace(0.0001, 2*np.pi, 100)


frj = jp._Rforce(raios, 0.)
frsap = sap.Rforce(raios, 0.)
frj2 = jp._R2deriv(raios, 0.)


ep = jp(raios, 0., phi = np.pi/2)
esap = sap(raios, 0., phi = np.pi/2)


o= Orbit([20.,30.,2.,-10.,20.,50.],radec=True,ro=8.,vo=220.)
ts= np.linspace(0,100,10000)
print(o.vxvv)
# o.integrate(ts, lhp+sap, method='rk6_c')
# o.plot()

#
plt.figure(1)
plt.plot(raios, frj)
plt.xlabel("Raio (kpc)")
plt.ylabel("Força em R")
plt.xlim(2, 12)
plt.ylim(-2.3, 3.3)

plt.figure(2)
plt.plot(raios, frj2)
plt.xlabel("Raio (kpc)")
plt.ylabel("Derivada segunda em R")
plt.xlim(2, 12)
plt.ylim(-10, 20.3)
plt.grid()
#
# plt.plot(raios, esap, 'r')
#
plt.show()

# plt.xlabel("Phi (rad)")
# plt.ylabel("Potencial (?)")
# print(ep)
# generate 2 2d grids for the x & y bounds
y, x = np.meshgrid(np.linspace(-10, 10, 100), np.linspace(-10, 10, 100))

# print(z.shape, x.shape, y.shape)
z = [[0 for i in range(len(raios))] for j in range(len(phis))]

# for i in range(len(x)):
#     for j in range(len(y)):
#         # z[i][j] = " Raio: " + str(raios[i]) + " Phi: " + str(phis[j])
#         z[i][j] = jp(np.sqrt(x[i]**2 + y[j]**2), 0., phi = np.arctan(y[j]/x[i]))
z = jp(np.sqrt(x**2 + y**2), 0., phi = -np.arctan(y/x))
# + mp(np.sqrt(x**2 + y**2), 0., phi = -np.arctan(y/x))


# z = z[:-1][:-1]
z = np.array(z)
# print(z.shape, x.shape, y.shape)
# x and y are bounds, so z should be the value *inside* those bounds.
# Therefore, remove the last value from the z array.
# z = z[:-1, :-1]
z_min, z_max = -np.abs(z).max(), np.abs(z).max()
#
fig, ax = plt.subplots()


# fig = plt.figure()
# ax = plt.axes(projection='3d')


#
# ax.plot3D(x, y, z, 'gray')
c = ax.pcolormesh(x, y, z, cmap='RdBu', vmin=z_min, vmax=z_max)
plt.xlabel("x (kpc)")
plt.ylabel("y (kpc)")
# ax.set_title('pcolormesh')
# set the limits of the plot to the limits of the data
ax.axis([x.min(), x.max(), y.min(), y.max()])
fig.colorbar(c, ax=ax)

plt.show()
# -1.3733506513947895

# SpiralArmsPotential.__init__()



# SpiralArmsPotential.__init__()
