import abc
import pylab as pl
import numpy as np

class Model(object):
    def __init__(self, label, input, atmosphere, surface):
        self.label = label

        self.runtime = input["runtime"]
        self.dt = input["dt"]
        self.t = 0.

        self.atmosphere = atmosphere
        self.surface = surface

        self.output = {}
        self.output["t"] = []

        print("Atmosphere Type: {0}".format(type(self.atmosphere).__name__))
        print("Surface Type: {0}".format(type(self.surface).__name__))

    def run(self):
        nt = int(self.runtime / self.dt) + 1

        # Time loop.
        for n in range(nt):
            # Calculate state
            self.surface.state(self.atmosphere)

            # Calculate tendency
            self.atmosphere.tendency(self.surface)

            # Integrate in time
            self.atmosphere.integrate(self.dt)

            # Save output
            self.output["t"].append(self.t)
            self.atmosphere.save_output()

            # Increate time
            self.t += self.dt

class Atmosphere(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, atmosphere_input):
        self.theta = atmosphere_input["theta"]

        self.output = {}
        self.output["theta"] = []

    @abc.abstractmethod
    def tendency(self, surface):
        """Calculate the tendencies"""
        return

    @abc.abstractmethod
    def integrate(self, dt):
        """Time integrate the variables"""
        return

    def save_output(self):
        self.output["theta"].append(self.theta)

class AtmosphereConstant(Atmosphere):
    def __init__(self, atmosphere_input):
        Atmosphere.__init__(self, atmosphere_input)

    def tendency(self, surface):
        return

    def integrate(self, dt):
        return

class AtmosphereBox(Atmosphere):
    def __init__(self, atmosphere_input):
        Atmosphere.__init__(self, atmosphere_input)
        self.h = atmosphere_input["h"]

    def tendency(self, surface):
        self.theta_tend = surface.wtheta / self.h

    def integrate(self, dt):
        self.theta += dt * self.theta_tend
        self.theta_tend = 0.

class AtmosphereMixedLayer(Atmosphere):
    def __init__(self, atmosphere_input):
        Atmosphere.__init__(self, atmosphere_input)
        self.h = atmosphere_input["h"]
        self.dtheta = atmosphere_input["dtheta"]
        self.gamma_theta = atmosphere_input["gamma_theta"]
        self.beta = atmosphere_input["beta"]

    def tendency(self, surface):
        self.h_tend = self.beta * surface.wtheta / self.dtheta
        self.theta_tend = (1. + self.beta) * surface.wtheta / self.h
        self.dtheta_tend = self.gamma_theta * self.h_tend - self.theta_tend

    def integrate(self, dt):
        self.h += dt * self.h_tend
        self.theta += dt * self.theta_tend
        self.dtheta += dt * self.dtheta_tend

        self.h_tend = 0.
        self.theta_tend = 0.
        self.dtheta_tend = 0.

class Surface(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, surface_input):
        return

    @abc.abstractmethod
    def state(self, atmosphere):
        """Update the value of the surface flux"""
        return

class SurfaceFixedFlux(Surface):
    def __init__(self, surface_input):
        Surface.__init__(self, surface_input)
        self.wtheta = surface_input["wtheta"]

    def state(self, atmosphere):
        return

class SurfaceFixedTemperature(Surface):
    def __init__(self, surface_input):
        Surface.__init__(self, surface_input)
        self.theta_surf = surface_input["theta_surf"]
        self.ra = surface_input["ra"]

    def state(self, atmosphere):
        self.wtheta = (self.theta_surf - atmosphere.theta) / self.ra

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

# Simple calculation examples to bypass complex code
surface_input_simple = { "theta_surf" : 305., "ra" : 50. }
surface_simple = SurfaceFixedTemperature(surface_input_simple)

atmosphere_input_simple = { "theta" : 300. }
atmosphere_simple = AtmosphereConstant(atmosphere_input_simple)

surface_simple.state(atmosphere_simple)
print("Calculated wtheta = {0}".format(surface_simple.wtheta))
