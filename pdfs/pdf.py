import numpy as np
import struct
import sys
from pylab import *

nx = 1024
ny = 1024

nbins = int(sys.argv[1])
if(len(sys.argv) > 2):
    alphaval = float(sys.argv[2])
else:
    alphaval = 1.

dx = 2./1024

# class pdf:
#   def __init__(self, pdfx, pdfy):
#     self.pdfx = pdfx
#     self.pdfy = pdfy
#   def p(self, val):
#     i = 0
#     while(True):
#       if(pdfx[i] > val):
#         return pdfy[i]

fin = open("w.xy.00030.0000250", "rb")
raw = fin.read(nx*ny*8)
w = np.array(struct.unpack('<{}d'.format(nx*ny), raw))
del(raw)

wavg = mean(w)
wstd = std(w)

tmp = w.reshape((ny,nx))
w2d = np.zeros((ny, nx+4))
w2d[:,2:nx+2] = tmp
w2d[:,0:2] = tmp[:,nx-2:nx]
w2d[:,nx+2::] = tmp[:,0:2]
del(tmp)

wgrad = (1.*w2d[:,0:-4] - 27.*w2d[:,1:-3] + 27.*w2d[:,2:-2] - 1.*w2d[:,3:-1]) / (24.*dx)
del(w2d)
wgrad = wgrad.reshape(nx*ny)

wgradavg = mean(wgrad)
wgradstd = std(wgrad)

#close('all')

# plot histogram
# figure()
# hist(wstat, nbins, normed=1, facecolor='#eeaaaa', histtype='stepfilled', alpha=0.4, label=r'w')

# create PDF
n, bins = np.histogram(w, nbins, density=1)
wpdfx = zeros(size(n))
wpdfy = zeros(size(n))
for k in range(size(n)):
  wpdfx[k] = 0.5*(bins[k]+bins[k+1])
  wpdfy[k] = n[k]

wbin = wpdfx[1] - wpdfx[0]

figure(1)
plot(wpdfx / wstd, wpdfx**0.*wpdfy / sum(wpdfx**0.*wpdfy*wbin), 'k+-', label='n = 0', alpha=alphaval)
plot(wpdfx / wstd, wpdfx**2.*wpdfy / sum(wpdfx**2.*wpdfy*wbin), 'r+-', label='n = 2', alpha=alphaval)
plot(wpdfx / wstd, wpdfx**4.*wpdfy / sum(wpdfx**4.*wpdfy*wbin), 'g+-', label='n = 4', alpha=alphaval)
plot(wpdfx / wstd, wpdfx**6.*wpdfy / sum(wpdfx**6.*wpdfy*wbin), 'b+-', label='n = 6', alpha=alphaval)
plot(wpdfx / wstd, wpdfx**8.*wpdfy / sum(wpdfx**8.*wpdfy*wbin), 'm+-', label='n = 8', alpha=alphaval)
legend(loc=0, frameon=False)

# create PDF
n, bins = np.histogram(wgrad, nbins, density=1)
wgradpdfx = zeros(size(n))
wgradpdfy = zeros(size(n))
for k in range(size(n)):
  wgradpdfx[k] = 0.5*(bins[k]+bins[k+1])
  wgradpdfy[k] = n[k]

wgbin = wgradpdfx[1] - wgradpdfx[0]

figure(2)
plot(wgradpdfx * (dx / wstd), wgradpdfx**0.*wgradpdfy / sum(wgradpdfx**0.*wgradpdfy*wgbin), 'k+-', label='n = 0', alpha=alphaval)
plot(wgradpdfx * (dx / wstd), wgradpdfx**2.*wgradpdfy / sum(wgradpdfx**2.*wgradpdfy*wgbin), 'r+-', label='n = 2', alpha=alphaval)
plot(wgradpdfx * (dx / wstd), wgradpdfx**4.*wgradpdfy / sum(wgradpdfx**4.*wgradpdfy*wgbin), 'g+-', label='n = 4', alpha=alphaval)
plot(wgradpdfx * (dx / wstd), wgradpdfx**6.*wgradpdfy / sum(wgradpdfx**6.*wgradpdfy*wgbin), 'b+-', label='n = 6', alpha=alphaval)
plot(wgradpdfx * (dx / wstd), wgradpdfx**8.*wgradpdfy / sum(wgradpdfx**8.*wgradpdfy*wgbin), 'm+-', label='n = 8', alpha=alphaval)
legend(loc=0, frameon=False)

