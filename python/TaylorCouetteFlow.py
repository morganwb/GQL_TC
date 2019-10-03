#!/usr/bin/env python
# coding: utf-8

# # Axisymmetric Taylor-Couette flow in Dedalus
# 
# 
# ![Taylor Couette Flow](http://upload.wikimedia.org/wikipedia/commons/3/3d/CouetteTaylorSystem.svg "Taylor Couette Flow")
# (image: wikipedia)
# 
# Taylor-Couette flow is characterized by three dimensionless numbers:
# 
# $\eta = R1/R2$, the ratio of the inner cylidner radius $R_1$ to the outer cylinder radius $R_2$
# 
# $\mu = \Omega_2/\Omega_1$, the ratio of the OUTER cylinder rotation rate $\Omega_2$ to the inner rate $\Omega_1$
# 
# $\mathrm{Re} = \Omega_1 R_1 \delta/\nu$, the Reynolds numbers, where $\delta = R_2 - R_1$, the gap width between the cylinders
# 
# 
# We non dimensionalize the flow in terms of 
# 
# $[\mathrm{L}] = \delta = R_2 - R_1$ 
# 
# $[\mathrm{V}] = R_1 \Omega_1$ 
# 
# $[\mathrm{M}] = \rho \delta^3$
# 
# And choose $\delta = 1$, $R_1 \Omega_1 = 1$, and $\rho = 1$.

import numpy as np
import time
import h5py

import dedalus.public as de 
from dedalus.extras import flow_tools

import logging
root = logging.root
for h in root.handlers:
    h.setLevel("INFO")
    
logger = logging.getLogger(__name__)

eta = 1./1.444 # R1/R2
alpha = 3.13   # vertical wavenumber
Re = 800.       # in units of R1*Omega1*delta/nu
mu = 0.        # Omega2/Omega1

# computed quantitites
omega_in = 1.
omega_out = mu * omega_in
r_in = eta/(1. - eta)
r_out = 1./(1. - eta)
height = 2.*np.pi/alpha 
v_l = 1. # by default, we set v_l to 1.
v_r = omega_out*r_out

r_basis = de.Chebyshev('r', 128, interval=(r_in, r_out), dealias=3/2)
z_basis = de.Fourier('z', 128, interval=(0., height), dealias=3/2)
domain = de.Domain([z_basis, r_basis], grid_dtype=np.float64)

TC = de.IVP(domain, variables=['p', 'u', 'v', 'w', 'ur', 'vr', 'wr'], ncc_cutoff=1e-8)
TC.meta[:]['r']['dirichlet'] = True
TC.parameters['nu'] = 1./Re
TC.parameters['v_l'] = v_l
TC.parameters['v_r'] = v_r
mu = TC.parameters['v_r']/TC.parameters['v_l'] * eta

TC.add_equation("r*ur + u + r*dz(w) = 0")
TC.add_equation("r*r*dt(u) - r*r*nu*dr(ur) - r*nu*ur - r*r*nu*dz(dz(u)) + nu*u + r*r*dr(p) = -r*r*u*ur - r*r*w*dz(u) + r*v*v")
TC.add_equation("r*r*dt(v) - r*r*nu*dr(vr) - r*nu*vr - r*r*nu*dz(dz(v)) + nu*v  = -r*r*u*vr - r*r*w*dz(v) - r*u*v")
TC.add_equation("r*dt(w) - r*nu*dr(wr) - nu*wr - r*nu*dz(dz(w)) + r*dz(p) = -r*u*wr - r*w*dz(w)")
TC.add_equation("ur - dr(u) = 0")
TC.add_equation("vr - dr(v) = 0")
TC.add_equation("wr - dr(w) = 0")

r = domain.grid(1, scales=domain.dealias)
z = domain.grid(0, scales=domain.dealias)
 
p_analytic = (eta/(1-eta**2))**2 * (-1./(2*r**2*(1-eta)**2) -2*np.log(r) +0.5*r**2 * (1.-eta)**2)
v_analytic = eta/(1-eta**2) * ((1. - mu)/(r*(1-eta)) - r * (1.-eta) * (1 - mu/eta**2)) 

# boundary conditions
TC.add_bc("left(u) = 0")
TC.add_bc("left(v) = v_l")
TC.add_bc("left(w) = 0")
TC.add_bc("right(u) = 0", condition="nz != 0")
TC.add_bc("right(v) = v_r")
TC.add_bc("right(w) = 0")
TC.add_bc("left(p) = 0", condition="nz == 0")

dt = max_dt = 1.
omega1 = TC.parameters['v_l']/r_in
period = 2*np.pi/omega1

ts = de.timesteppers.RK443
IVP = TC.build_solver(ts)
IVP.stop_sim_time = 15.*period
IVP.stop_wall_time = np.inf
IVP.stop_iteration = np.inf

p = IVP.state['p']
u = IVP.state['u']
v = IVP.state['v']
w = IVP.state['w']
ur = IVP.state['ur']
vr = IVP.state['vr']
wr = IVP.state['wr']

phi = domain.new_field(name='phi')

for f in [phi,p,u,v,w,ur,vr,wr]:
    f.set_scales(domain.dealias, keep_data=False)

v['g'] = v_analytic
v.differentiate('r',out=vr)


def filter_field(field,frac=0.5):
    dom = field.domain
    local_slice = dom.dist.coeff_layout.slices(scales=dom.dealias)
    coeff = []
    for i in range(dom.dim)[::-1]:
        coeff.append(np.linspace(0,1,dom.global_coeff_shape[i],endpoint=False))
    cc = np.meshgrid(*coeff)
    
    field_filter = np.zeros(dom.local_coeff_shape,dtype='bool')
    for i in range(dom.dim):
        field_filter = field_filter | (cc[i][local_slice] > frac)
    field['c'][field_filter] = 0j


# incompressible perturbation, arbitrary vorticity
# u = -dz(phi)
# w = dr(phi) + phi/r
phi['g'] = 1e-3* np.random.randn(*v['g'].shape)*np.sin(np.pi*(r - r_in))
filter_field(phi)
phi.differentiate('z',out=u)
u['g'] *= -1
phi.differentiate('r',out=w)
w['g'] += phi['g']/r

u.differentiate('r',out=ur)
w.differentiate('r',out=wr)

CFL = flow_tools.CFL(IVP, initial_dt=1e-3, cadence=5, safety=0.3,
                     max_change=1.5, min_change=0.5)
CFL.add_velocities(('u', 'w'))

analysis1 = IVP.evaluator.add_file_handler("scalar_data", iter=10)
analysis1.add_task("integ(0.5 * (u*u + v*v + w*w))", name="total kinetic energy")
analysis1.add_task("integ(0.5 * (u*u + w*w))", name="meridional kinetic energy")
analysis1.add_task("integ((u*u)**0.5)", name='u_rms')
analysis1.add_task("integ((w*w)**0.5)", name='w_rms')

# Snapshots every half an inner rotation period
analysis2 = IVP.evaluator.add_file_handler('snapshots',sim_dt=0.5*period, max_size=2**30)
analysis2.add_system(IVP.state, layout='g')

# radial profiles every 100 timestpes
analysis3 = IVP.evaluator.add_file_handler("radial_profiles", iter=100)
analysis3.add_task("integ(r*v, 'z')", name='Angular Momentum')

dt = CFL.compute_dt()
# Main loop
start_time = time.time()

while IVP.ok:
    IVP.step(dt)
    if IVP.iteration % 10 == 0:
        logger.info('Iteration: %i, Time: %e, dt: %e' %(IVP.iteration, IVP.sim_time, dt))
    dt = CFL.compute_dt()


end_time = time.time()

# Print statistics
logger.info('Total time: %f' %(end_time-start_time))
logger.info('Iterations: %i' %IVP.iteration)
logger.info('Average timestep: %f' %(IVP.sim_time/IVP.iteration))

