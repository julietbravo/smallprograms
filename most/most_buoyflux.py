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

def eval_bd(L):
  return zsl/L * fm(L)**3

def eval_w(L):
  return zsl/L * fmw(L)**3

zL = np.linspace(-100., 10., 1e5)
L  = zsl / zL

# Evaluate the function (this has to be done only once for a range of Ri).
eval0_bd = np.zeros(L.size)
eval0_w  = np.zeros(L.size)
for i in range(L.size):
  eval0_bd[i] = eval_bd(L[i])
  eval0_w [i] = eval_w (L[i])

# Find the value that matches with the value of Ri.
Ri = -zsl*kappa*B0 / u0**3

if (max(eval0_bd-Ri)*min(eval0_bd-Ri) > 0):
  zL0_bd = np.nan
else:
  zL0_bd = np.interp(0., eval0_bd-Ri, zL)
db_bd  = zL0_bd * fm(zsl/zL0_bd)**2 * u0**2 / (kappa * zsl * fh(zsl/zL0_bd) )

if (max(eval0_w-Ri)*min(eval0_w-Ri) > 0):
  zL0_w = np.nan
else:
  zL0_w = np.interp(0., eval0_w-Ri, zL)
db_w  = zL0_w * fm(zsl/zL0_w)**2 * u0**2 / (kappa * zsl * fh(zsl/zL0_w) )

ustar_bd = u0 * fm (zsl/zL0_bd)
ustar_w  = u0 * fmw(zsl/zL0_bd)
ustar_n  = u0 * fm (np.inf)

print('BD: z/L = {0}, db = {1}, B0 = {2}, ustar = {3}'.format(zL0_bd, db_bd, B0, ustar_bd))
print('W:  z/L = {0}, db = {1}, B0 = {2}, ustar = {3}'.format(zL0_w , db_w , B0, ustar_w ))
print("ustar_neutral = {0}".format(ustar_n))

pl.close('all')

pl.figure()
pl.plot(zL, eval0_bd)
pl.plot(zL, eval0_w )
pl.xlabel('z/L')
pl.ylabel('eval')

dzL = zL[1] - zL[0]
deval0_bd = np.gradient(eval0_bd, dzL)
deval0_w  = np.gradient(eval0_w , dzL)

pl.figure()
pl.plot(zL, deval0_bd)
pl.plot(zL, deval0_w )
pl.xlabel('z/L')
pl.ylabel('deval/dzL')
