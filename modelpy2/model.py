class Model(object):
    def __init__(self, atmosphere):
        self.runtime = 7200.
        self.dt = 60.

        self.atmosphere = atmosphere

        print("Atmosphere Type: ", type(self.atmosphere))

    def run(self):
        nt = int(self.runtime / self.dt) + 1
        t = 0.

        # Time loop.
        for n in range(nt):
            # Calculate state
            # Calculate tendency
            self.atmosphere.tendency()

            # Integrate in time
            self.atmosphere.integrate(self.dt)

            # Save output
            self.atmosphere.output(t)

            # Increate time
            t += self.dt

class AtmosphereInput(object):
    def __init__(self):
        self.h = None
        self.theta = None
        self.dtheta = None

class Atmosphere(object):
    def __init__(self, atmosphere_input):
        self.h = atmosphere_input.h
        self.theta = atmosphere_input.theta
        self.dtheta = atmosphere_input.dtheta

    def tendency(self):
        return

    def integrate(self, dt):
        return

    def output(self, t):
        print("Saving output at t = {0}".format(t))
        print("h = {0}, theta = {1}, dtheta = {2}".format(self.h, self.theta, self.dtheta))

class AtmosphereBox(Atmosphere):
    def __init__(self, atmosphere_input):
        Atmosphere.__init__(self, atmosphere_input)

class AtmosphereMixedLayer(Atmosphere):
    def __init__(self, atmosphere_input):
        Atmosphere.__init__(self, atmosphere_input)
        self.gamma_theta = 0.006
        self.wtheta = 0.1
        self.beta = 0.2

    def tendency(self):
        self.h_tend = self.beta * self.wtheta / self.dtheta
        self.theta_tend = (1. + self.beta) * self.wtheta / self.h
        self.dtheta_tend = self.gamma_theta * self.h_tend - self.theta_tend

    def integrate(self, dt):
        self.h += dt * self.h_tend
        self.theta += dt * self.theta_tend
        self.dtheta += dt * self.dtheta_tend

        self.h_tend = 0.
        self.theta_tend = 0.
        self.dtheta_tend = 0.

# Set there the input options
atmosphere_input = AtmosphereInput()
atmosphere_input.h = 100.
atmosphere_input.theta = 300.
atmosphere_input.dtheta = 1.

model = Model( AtmosphereBox( atmosphere_input) )
model.run()

model = Model( AtmosphereMixedLayer( atmosphere_input) )
model.run()
