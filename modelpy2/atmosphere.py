class AtmosphereConstant:
    def __init__(self, atmosphere_input):
        self.theta = atmosphere_input["theta"]

        self.output = {}
        self.output["theta"] = []

    def tendency(self, surface):
        return

    def integrate(self, dt):
        return

    def save_output(self):
        self.output["theta"].append(self.theta)

class AtmosphereBox:
    def __init__(self, atmosphere_input):
        self.theta = atmosphere_input["theta"]
        self.h = atmosphere_input["h"]

        self.output = {}
        self.output["theta"] = []
        self.output["h"] = []

    def tendency(self, surface):
        self.theta_tend = surface.wtheta / self.h

    def integrate(self, dt):
        self.theta += dt * self.theta_tend
        self.theta_tend = 0.

    def save_output(self):
        self.output["theta"].append(self.theta)

class AtmosphereMixedLayer:
    def __init__(self, atmosphere_input):
        self.theta = atmosphere_input["theta"]
        self.h = atmosphere_input["h"]
        self.dtheta = atmosphere_input["dtheta"]
        self.gamma_theta = atmosphere_input["gamma_theta"]
        self.beta = atmosphere_input["beta"]

        self.output = {}
        self.output["theta"] = []
        self.output["h"] = []
        self.output["dtheta"] = []

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

    def save_output(self):
        self.output["theta"].append(self.theta)
        self.output["h"].append(self.h)
        self.output["dtheta"].append(self.dtheta)

