from galpy.orbit import Orbit
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import time

from galpy.potential import MWPotential2014
from galpy.potential import LogarithmicHaloPotential
from galpy.potential import SpiralArmsPotential
from galpy.potential import plotRotcurve
from galpy.potential.JunqueiraPotential import JunqueiraPotential

start = time.time()
c= SkyCoord(ra=20.*u.deg,dec=30.*u.deg,distance=2.*u.kpc,
                pm_ra_cosdec=-5.*u.mas/u.yr,pm_dec=10.*u.mas/u.yr,
                radial_velocity=100.*u.km/u.s)
o= Orbit(c)
o2 = Orbit(c)


lhp = LogarithmicHaloPotential()
sap = SpiralArmsPotential(N = 4, omega = 25)
jp = JunqueiraPotential()

# plotRotcurve(MWPotential2014 + jp, Rrange=[0.01,10.],grid=1001,yrange=[0.,1.2], phi = 0.)
# plotRotcurve(MWPotential2014+sap, Rrange=[0.01,10.],grid=1001,yrange=[0.,1.2], phi = 0.)
# plotRotcurve(lhp+jp, Rrange=[0.01,10.],yrange=[0.,2.], phi = 0.)

ts= np.linspace(0.,1000.,200001)

# o.integrate(ts,lhp+jp, method='rk6_c', numcores = 4)
o.integrate(ts,MWPotential2014+jp, method='rk6_c')
# o.plot(d1='x',d2='y')
# plot([o.R()],[o.z()],'ro')
# plt.show()
# o2.integrate(ts,MWPotential2014+jp, method='rk6_c', numcores = 4)
o2.integrate(ts,lhp+sap, method='odeint')
o2.plot(d1='x',d2='y')

end = time.time()
# plot([o.R()],[o.z()],'ro')

print(end-start, "s")
plt.savefig('Orbita.png')
#plt.show()
