from atmosphere import AtmosphereConstant
from surface import SurfaceFixedTemperature

# Simple calculation examples to bypass complex code
surface_input_simple = { "theta_surf" : 305., "ra" : 50. }
surface_simple = SurfaceFixedTemperature(surface_input_simple)

atmosphere_input_simple = { "theta" : 300. }
atmosphere_simple = AtmosphereConstant(atmosphere_input_simple)

surface_simple.state(atmosphere_simple)
print("Calculated wtheta = {0}".format(surface_simple.wtheta))
