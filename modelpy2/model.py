class Model(object):
    def __init__(self, atmosphere):
        self.runtime = 7200.
        self.dt = 360.

        self.atmosphere = atmosphere

        print("Atmosphere Type: ", type(self.atmosphere))

    def run(self):
        nt = int(self.runtime / self.dt) + 1
        t = 0.

        # Time loop.
        for n in range(nt):
            # Calculate state
            # Calculate tendency
            # Integrate tendency

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

    def output(self, t):
        print("Saving output at t = {0}".format(t))

class AtmosphereBox(Atmosphere):
    def __init__(self, atmosphere_input):
        Atmosphere.__init__(self, atmosphere_input)

# Set there the input options
atmosphere_input = AtmosphereInput()
atmosphere_input.h = 1000.
atmosphere_input.theta = 300.
atmosphere_input.dtheta = 1.

model = Model( AtmosphereBox( atmosphere_input) )

model.run()
