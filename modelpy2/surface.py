import abc

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

