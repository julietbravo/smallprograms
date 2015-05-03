from atmosphere import *
from surface import *

class Model:
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

