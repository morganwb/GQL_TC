"""
Usage:
  tc_eigenvalue_3d.py [--re=<reynolds> --eta=<eta> --m=<initial_m> --k=<k> --mu=<mu> --ar=<Gamma> --alpha=<alpha>]

Options:
  --re=<reynolds>  Reynolds number for simulation [default: 80]
  --eta=<eta>      Eta - ratio of R1/R2 [default: 0.692520775623]
  --m=<initial_m>  m mode [default: 0]
  --k=<k>          k mode [default: 1]
  --mu=<mu>        mu = Omega2/Omega1 [default: 0]
  --ar=<Gamma>     Aspect ratio (height/width) [default: 3]
  --alpha=<alpha>  z wavenumber [default: 3.13]
"""

import numpy as np
from dedalus import public as de
import logging
from docopt import docopt
from eigentools import Eigenproblem
logger = logging.getLogger(__name__)

nr = 32

args=docopt(__doc__)
Re1=float(args['--re'])
eta=np.float(args['--eta'])
mu = np.float(args['--mu'])
Gamma = float(args['--ar'])
alpha = float(args['--alpha'])
m = int(args['--m'])
k = int(args['--k'])

"""
delta = R2 - R1
mu = Omega2/Omega1
eta = R1/R2

scale [L] = delta
scale [T] = delta/(R1 Omega1)
scale [V] = R1 Omega1

Default parameters from Barenghi (1991, J. Comp. Phys.).
"""
#derived parameters
R1 = eta/(1. - eta)
R2 = 1./(1-eta)
Omega1 = 1/R1
Omega2 = mu*Omega1
nu = R1 * Omega1/Re1
midpoint = R1 + (R2-R1) / 2
if alpha:
    Lz = 2*np.pi/alpha
else:
    Lz = Gamma 

#Taylor Number
#Ta = ((1+eta)**4/(64*eta**2) ) * ( (R2-R1)**2 * (R2+R1)**2 * (Omega1-Omega2)**2 ) / nu**2 #Grossman lohse
#Ta1 = (2*Omega1**2*(R2-R1)**4*eta**2 ) / (nu**2 *(11-eta**2)  )
A = (1/eta - 1.)*(mu-eta**2)/(1-eta**2)
B = eta*(1-mu)/((1-eta)*(1-eta**2))
#Ta = -4*A*Omega1*Re1**2
logger.info("-4*A = {}".format(-4*A))
Ta = 64/9. * (Re1*eta/(1-eta))**2 * (1-mu)*(1-4*mu)
#Ta = -4*Omega1**2 * R1**4 * Re1**2 * (1-mu)*(1-mu/eta**2)/(1-eta**2)**2
Ro_inv = (2 * Omega2 * (R2-R1) ) /  (np.abs(Omega1-Omega2)*R1 )

logger.info("Re:{:.3e}, eta:{:.4e}, mu:{:.4e}".format(Re1,eta,mu))
logger.info("Taylor Number:{:.2e}, Ro^(-1):{:.2e}".format(Ta,Ro_inv))

logger.info("Lz set to {:.6e}".format(Lz))

variables = ['u','ur','v','vr','w','wr','p']

#domain
r_basis = de.Chebyshev('r', nr, interval=[R1, R2])

bases = [r_basis]
domain = de.Domain(bases) 

#problem
problem = de.EVP(domain, eigenvalue='sigma', variables=variables)

#params into equations
problem.parameters['eta']=eta
problem.parameters['mu']=mu
problem.parameters['Lz']=Lz
problem.parameters['nu']=nu
problem.parameters['R1']=R1
problem.parameters['R2']=R2
problem.parameters['Omega1']=Omega1
problem.parameters['Omega2']=Omega2

problem.parameters['pi']=np.pi
problem.parameters['kz'] = 2*np.pi/Lz * k
problem.parameters['m'] = m

#Substitutions

"""
this implements the cylindrical del operators. 
NB: ASSUMES THE EQUATION SET IS PREMULTIPLIED BY A POWER OF r (SEE BELOW)!!!

Lap_s --> scalar laplacian
Lap_r --> r component of vector laplacian
Lap_t --> theta component of vector laplacian
Lap_z --> z component of vector laplacian

"""

problem.substitutions['A'] = '(1/eta - 1.)*(mu-eta**2)/(1-eta**2)'
problem.substitutions['B'] = 'eta*(1-mu)/((1-eta)*(1-eta**2))'

problem.substitutions['v0'] = 'A*r + B/r'       #background profile? forcing instead of forcing the boundaries
problem.substitutions['dv0dr'] = 'A - B/(r*r)'  #d/dr of background forcing

problem.substitutions['dtheta(f)'] = '1j*m*f'
problem.substitutions['dz(f)'] = '1j*kz*f'
problem.substitutions['dt(f)'] = 'sigma*f'

# assume pre-multiplication by r*r
problem.substitutions['Lap_s(f, f_r)'] = "r*r*dr(f_r) + r*f_r + dtheta(dtheta(f)) + r*r*dz(dz(f))"
problem.substitutions['Lap_r'] = "Lap_s(u, ur) - u - 2*dtheta(v)"
problem.substitutions['Lap_t'] = "Lap_s(v, vr) - v + 2*dtheta(u)"
problem.substitutions['Lap_z'] = "Lap_s(w, wr)"

# momentum equations
problem.add_equation("r*r*dt(u) - nu*Lap_r - 2*r*v0*v + r*v0*dtheta(u) + r*r*dr(p) = 0")
problem.add_equation("r*r*dt(v) - nu*Lap_t + r*r*dv0dr*u + r*v0*u + r*v0*dtheta(v) + r*dtheta(p)  = 0")
problem.add_equation("r*r*dt(w) - nu*Lap_z + r*r*dz(p) + r*v0*dtheta(w) = 0.")

#continuity
problem.add_equation("r*ur + u + dtheta(v) + r*dz(w) = 0")

#Auxillilary equations
problem.add_equation("ur - dr(u) = 0")
problem.add_equation("vr - dr(v) = 0")
problem.add_equation("wr - dr(w) = 0")

#Boundary Conditions
problem.add_bc("left(u) = 0")
problem.add_bc("right(u) = 0")
problem.add_bc("left(v) = 0")
problem.add_bc("right(v) = 0")
problem.add_bc("left(w) = 0")
problem.add_bc("right(w) = 0")


ep = Eigenproblem(problem)

growth, index, freq = ep.growth_rate({})

logger.info("Growth rate = {:16.15e}; frequency = {:16.15e}".format(growth, freq[0]))

#ep.spectrum(spectype='hires')
