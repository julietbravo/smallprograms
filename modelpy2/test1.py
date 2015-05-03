import pylab as pl
import numpy as np

from model import *

# Test case 1: Box model for atmosphere
model_input = {}
model_input["runtime"] = 7200.
model_input["dt"] = 60.

atmosphere_input = {}
atmosphere_input["h"] = 100.
atmosphere_input["theta"] = 300.
atmosphere_input["dtheta"] = 1.

surface_input = {}
surface_input["wtheta"] = 0.1

model1 = Model( "Box with fixed flux",
                model_input,
                AtmosphereBox(atmosphere_input),
                SurfaceFixedFlux(surface_input) )

model1.run()

# Test case 2: Mixed-layer model for atmosphere
atmosphere_input2 = atmosphere_input.copy()
atmosphere_input2["gamma_theta"] = 0.006
atmosphere_input2["beta"] = 0.2

surface_input2 = surface_input.copy()
surface_input2["wtheta"] = 0.1

model2 = Model( "Growing mixed layer with fixed flux",
                model_input,
                AtmosphereMixedLayer(atmosphere_input2),
                SurfaceFixedFlux(surface_input2) )
model2.run()

# Test case 3: Mixed-layer model with a interacting surface
surface_input3 = {}
surface_input3["theta_surf"] = 305.
surface_input3["ra"] = 50.

model3 = Model( "Mixed layer with fixed temperature",
                model_input,
                AtmosphereMixedLayer(atmosphere_input2),
                SurfaceFixedTemperature(surface_input3) )
model3.run()

models = [ model1, model2, model3 ]

# Plot the output of the three cases
pl.figure()
for m in models:
    t = np.asarray(m.output["t"])
    theta = np.asarray(m.atmosphere.output["theta"])
    pl.plot(t, theta, label=m.label)
pl.xlabel('t')
pl.ylabel('theta')
pl.legend(loc=0, frameon=False)

