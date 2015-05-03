class Model(object):
    def __init__(self, atmosphere, surface=None):
        self.runtime = 7200.
        self.dt = 60.

        self.atmosphere = atmosphere
        self.surface = surface

        print("Atmosphere Type: ", type(self.atmosphere))

    def run(self):
        nt = int(self.runtime / self.dt) + 1
        t = 0.

        # Time loop.
        for n in range(nt):
            # Calculate state
            # Calculate tendency
            self.atmosphere.tendency(self.surface)

            # Integrate in time
            self.atmosphere.integrate(self.dt)

            # Save output
            self.atmosphere.output(t)

            # Increate time
            t += self.dt

class Atmosphere(object):
    def __init__(self, atmosphere_input):
        self.theta = atmosphere_input["theta"]

    def integrate(self, dt):
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
    def __init__(self, surface_input):
        self.wtheta = surface_input["wtheta"]

class SurfaceFixedFlux(Surface):
    def __init__(self, surface_input):
        Surface.__init__(self, surface_input)

# Test case 1: Box model for atmosphere
atmosphere_input = {}
atmosphere_input["h"] = 100.
atmosphere_input["theta"] = 300.
atmosphere_input["dtheta"] = 1.

surface_input = {}
surface_input["wtheta"] = 0.1

model = Model( AtmosphereBox(atmosphere_input),
               SurfaceFixedFlux(surface_input) )

model.run()

# Test case 2: Mixed-layer model for atmosphere
atmosphere_input2 = atmosphere_input.copy()
atmosphere_input2["gamma_theta"] = 0.006
atmosphere_input2["beta"] = 0.2

surface_input2 = surface_input.copy()
surface_input2["wtheta"] = 0.1

model2 = Model( AtmosphereMixedLayer(atmosphere_input2),
                SurfaceFixedFlux(surface_input2) )
model2.run()
