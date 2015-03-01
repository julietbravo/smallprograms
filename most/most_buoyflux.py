# !/usr/bin/python
import numpy as np
import pylab as pl

# Input parameters.
#B0 = 0.1 * (9.81/300.)
#u0 = 1.
B0  = -2.68e-2 * (9.81/300.)
u0  = 10.
z0m = 0.1
z0h = 0.1
zsl = 10

# Constants.
kappa = 0.4

# Integrated flux gradient relationships following Businger-Dyer.
def psim(zeta):
  if(zeta <= 0):
    x     = (1. - 16. * zeta) ** (0.25)
    psim  = 3.14159265 / 2. - 2. * np.arctan(x) + np.log( (1.+x) ** 2. * (1. + x ** 2.) / 8.)
  else:
    psim  = -2./3. * (zeta - 5./0.35)*np.exp(-0.35 * zeta) - zeta - (10./3.) / 0.35
  return psim

def psih(zeta):
  if(zeta <= 0):
    x     = (1. - 16. * zeta) ** (0.25)
    psih  = 2. * np.log( (1. + x ** 2.) / 2. )
  else:
    psih  = -2./3. * (zeta - 5./0.35)*np.exp(-0.35 * zeta) - (1. + (2./3.) * zeta) ** (1.5) - (10./3.) / 0.35 + 1.
  return psih

# Integrated flux gradient relationships following Wilson.
def psimw(zeta):
  if(zeta <= 0):
    x     = (1. + 3.6 * abs(zeta) ** (2./3.)) ** (-0.5)
    psimw = 3. * np.log( (1. + 1. / x) / 2.)
  else:
    psimw  = -2./3. * (zeta - 5./0.35)*np.exp(-0.35 * zeta) - zeta - (10./3.) / 0.35
  return psimw
  
def psihw(zeta):
  if(zeta <= 0):
    x     = (1. + 7.9 * abs(zeta) ** (2./3.)) ** (-0.5)
    psihw  = 3. * np.log( (1. + 1. / x) / 2.)
  else:
    psihw  = -2./3. * (zeta - 5./0.35)*np.exp(-0.35 * zeta) - (1. + (2./3.) * zeta) ** (1.5) - (10./3.) / 0.35 + 1.
  return psihw

def fm(L):
  return kappa / (np.log(zsl/z0m) - psim(zsl/L) + psim(z0m/L))

def fh(L):
  return kappa / (np.log(zsl/z0h) - psih(zsl/L) + psih(z0h/L))

def fmw(L):
  return kappa / (np.log(zsl/z0m) - psimw(zsl/L) + psimw(z0m/L))

def fhw(L):
  return kappa / (np.log(zsl/z0h) - psihw(zsl/L) + psihw(z0h/L))

def eval_bd(L, gamma):
  return zsl/L * fm(L)**3 + gamma

def eval_w(L, gamma):
  return zsl/L * fmw(L)**3 + gamma

def create_zL(nzL):
  zL_tmp = np.zeros(nzL)
  zL = np.zeros(nzL)

  # Calculate the non-streched part between -10 to 10 z/L with 75% of the points.
  dzL = 20. / (3./4.*nzL-1)
  zL_tmp[0] = -10.

  for n in range(1, 3*nzL/4):
    zL_tmp[n] = zL_tmp[n-1] + dzL

  # Stretch the remainder of the z/L values far down for free convection.
  zLend = 1.e4 - 10.

  # Find stretching that ends up at the correct value using geometric progression.
  r  = 1.01
  r0 = 1.e9
  while (abs( (r-r0)/r0 ) > 1.e-10):
    r0 = r
    r  = ( 1. - (zLend/dzL)*(1.-r) )**(4./nzL)

  for n in range(3*nzL/4, nzL):
    zL_tmp[n] = zL_tmp[n-1] + dzL
    dzL *= r

  # Calculate the final array and delete the temporary array.
  for n in range(nzL):
    zL[n] = -zL_tmp[nzL-n-1];

  return zL

zL = create_zL(100000)
L  = zsl / zL

gamma0 = zsl*kappa*B0 / u0**3

eval0_bd = np.zeros(L.size)
eval0_w  = np.zeros(L.size)
for i in range(L.size):
  eval0_bd[i] = eval_bd(L[i], gamma0)
  eval0_w [i] = eval_w (L[i], gamma0)

if (max(eval0_bd)*min(eval0_bd) > 0):
  zL0_bd = np.nan
else:
  zL0_bd = np.interp(0., eval0_bd, zL)
db_bd  = zL0_bd * fm(zsl/zL0_bd)**2 * u0**2 / (kappa * zsl * fh(zsl/zL0_bd) )

if (max(eval0_w)*min(eval0_w) > 0):
  zL0_w = np.nan
else:
  zL0_w = np.interp(0., eval0_w, zL)
db_w  = zL0_w * fm(zsl/zL0_w)**2 * u0**2 / (kappa * zsl * fh(zsl/zL0_w) )

ustar_bd = u0 * fm (zsl/zL0_bd)
ustar_w  = u0 * fmw(zsl/zL0_bd)
ustar_n  = u0 * fm (np.inf)

print('BD: z/L = {0}, db = {1}, B0 = {2}, ustar = {3}'.format(zL0_bd, db_bd, B0, ustar_bd))
print('W:  z/L = {0}, db = {1}, B0 = {2}, ustar = {3}'.format(zL0_w , db_w , B0, ustar_w ))
print("ustar_neutral = {0}".format(ustar_n))

#pl.close('all')

pl.figure(1)
pl.plot(zL, eval0_bd)
pl.plot(zL, eval0_w )
pl.xlabel('z/L')
pl.ylabel('eval')

"""
dzL = zL[1] - zL[0]
deval0_bd = np.gradient(eval0_bd, dzL)
deval0_w  = np.gradient(eval0_w , dzL)

pl.figure()
pl.plot(zL, deval0_bd)
pl.plot(zL, deval0_w )
pl.xlabel('z/L')
pl.ylabel('deval/dzL')
"""
