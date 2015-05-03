import abc

class Model(object):
    def __init__(self, input, atmosphere, surface):
        self.runtime = input["runtime"]
        self.dt = input["dt"]

        self.atmosphere = atmosphere
        self.surface = surface

        print("Atmosphere Type: ", type(self.atmosphere))

    def run(self):
        nt = int(self.runtime / self.dt) + 1
        t = 0.

        # Time loop.
        for n in range(nt):
            # Calculate state
            self.surface.state(self.atmosphere)

            # Calculate tendency
            self.atmosphere.tendency(self.surface)

            # Integrate in time
            self.atmosphere.integrate(self.dt)

            # Save output
            self.atmosphere.output(t)

            # Increate time
            t += self.dt

class Atmosphere(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, atmosphere_input):
        self.theta = atmosphere_input["theta"]

    @abc.abstractmethod
    def tendency(self, surface):
        """Calculate the tendencies"""
        return

    @abc.abstractmethod
    def integrate(self, dt):
        """Time integrate the variables"""
        return


    def output(self, t):
        print("Saving output at t = {0}".format(t))
        print("theta = {0}".format(self.theta))

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

class SurfaceFixedTemp(Surface):
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

model = Model( model_input,
               AtmosphereBox(atmosphere_input),
               SurfaceFixedFlux(surface_input) )

model.run()

# Test case 2: Mixed-layer model for atmosphere
atmosphere_input2 = atmosphere_input.copy()
atmosphere_input2["gamma_theta"] = 0.006
atmosphere_input2["beta"] = 0.2

surface_input2 = surface_input.copy()
surface_input2["wtheta"] = 0.1

model2 = Model( model_input,
                AtmosphereMixedLayer(atmosphere_input2),
                SurfaceFixedFlux(surface_input2) )
model2.run()

# Test case 3: Mixed-layer model with a interacting surface
surface_input3 = {}
surface_input3["theta_surf"] = 303.
surface_input3["ra"] = 50.

model3 = Model( model_input,
                AtmosphereMixedLayer(atmosphere_input2),
                SurfaceFixedTemp(surface_input3) )
model3.run()
