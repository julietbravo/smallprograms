class SurfaceFixedFlux:
    def __init__(self, surface_input):
        self.wtheta = surface_input["wtheta"]

    def state(self, atmosphere):
        return

class SurfaceFixedTemperature:
    def __init__(self, surface_input):
        self.theta_surf = surface_input["theta_surf"]
        self.ra = surface_input["ra"]

    def state(self, atmosphere):
        self.wtheta = (self.theta_surf - atmosphere.theta) / self.ra

