import numpy as np
import pylab as pl

from atmosphere import AtmosphereConstant
from surface import SurfaceFixedTemperature

# Simple calculation examples to bypass complex code
surface_input_simple = { "theta_surf" : 305., "ra" : 50. }
surface_simple = SurfaceFixedTemperature(surface_input_simple)

atmosphere_input_simple = { "theta" : 300. }
atmosphere_simple = AtmosphereConstant(atmosphere_input_simple)

surface_simple.state(atmosphere_simple)
print("Calculated wtheta = {0}".format(surface_simple.wtheta))

ra_array = np.linspace(10., 50, 100)
wtheta_array = np.empty(ra_array.size)

for n in range(ra_array.size):
    surface_simple.ra = ra_array[n]
    surface_simple.state(atmosphere_simple)
    wtheta_array[n] = surface_simple.wtheta

pl.figure()
pl.plot(ra_array, wtheta_array)
pl.xlabel('ra')
pl.ylabel('wtheta')
