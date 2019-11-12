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
from configparser import ConfigParser
import argparse
from pathlib import Path
import dedalus.public as de 
from dedalus.extras import flow_tools
from dedalus.tools  import post
import logging
root = logging.root
for h in root.handlers:
    h.setLevel("INFO")
    
logger = logging.getLogger(__name__)

# Parses .cfg filename passed to script
parser = argparse.ArgumentParser(description='Passes filename')
parser.add_argument('filename', metavar='Rc', type=str, help='Provides config file to parser')
args = parser.parse_args()
filename = Path(vars(args)['filename'])
outbase = Path("data")

# Parse .cfg file to set global parameters for script
config = ConfigParser()
config.read(str(filename))
alpha = config.getfloat('parameters','alpha')
eta = eval(config.get('parameters','eta'))
Re = config.getfloat('parameters','Re')
mu = config.getfloat('parameters','mu')

nx = config.getint('parameters','nx')
nz = config.getint('parameters','nz')
ntheta = config.getint('parameters','ntheta')

if ntheta == 0:
	threeD = False
else:
	threeD = True

# computed quantitites
omega_in = 1.
omega_out = mu * omega_in
r_in = eta/(1. - eta)
r_out = 1./(1. - eta)
height = 2.*np.pi/alpha 
v_l = 1. # by default, we set v_l to 1.
v_r = omega_out*r_out


r_basis = de.Chebyshev('r', nx, interval=(r_in, r_out), dealias=3/2)
z_basis = de.Fourier('z', nz, interval=(0., height), dealias=3/2)

if threeD:
    theta_basis = de.Fourier('theta', ntheta, interval=(0., 2*np.pi), dealias=3/2)
    domain = de.Domain([z_basis, theta_basis, r_basis], grid_dtype=np.float64)
else:
    domain = de.Domain([z_basis, r_basis], grid_dtype=np.float64)

TC = de.IVP(domain, variables=['p', 'u', 'v', 'w', 'ur', 'vr', 'wr'], ncc_cutoff=1e-8)
TC.meta[:]['r']['dirichlet'] = True
TC.parameters['nu'] = 1./Re
TC.parameters['v_l'] = v_l
TC.parameters['v_r'] = v_r
TC.parameters['eta'] = eta
mu = TC.parameters['v_r']/TC.parameters['v_l'] * eta
TC.parameters['mu'] = mu

# adds substitutions

TC.substitutions['A'] = '(1/eta - 1.)*(mu-eta**2)/(1-eta**2)'
TC.substitutions['B'] = 'eta*(1-mu)/((1-eta)*(1-eta**2))'
TC.substitutions['v0'] = 'A*r + B/r'
TC.substitutions['dv0dr'] = 'A - B/(r*r)'
TC.substitutions['DivU'] = 'dr(r*u) + dtheta(v) + r*dz(w)'

if threeD:
    TC.substitutions['Lap_s(f, f_r)'] = "r*r*dr(f_r) + r*f_r + dtheta(dtheta(f)) + r*r*dz(dz(f))"
    TC.substitutions['Lap_r'] = "Lap_s(u, ur) - u - 2*dtheta(v)"
    TC.substitutions['Lap_t'] = "Lap_s(v, vr) - v + 2*dtheta(u)"
    TC.substitutions['Lap_z'] = "Lap_s(w, wr)"
    TC.substitutions['UdotGrad_s(f, f_r)'] = "r*r*u*f_r + r*v*dtheta(f) + r*r*w*dz(f)"
    TC.substitutions['UdotGrad_r'] = "UdotGrad_s(u, ur) - r*v*v"
    TC.substitutions['UdotGrad_t'] = "UdotGrad_s(v, vr) + r*u*v"
    TC.substitutions['UdotGrad_z'] = "UdotGrad_s(w, wr)"
else:
    TC.substitutions['Lap_s(f, f_r)'] = "r*dr(f_r) + f_r + r*dz(dz(f))"
    TC.substitutions['Lap_r'] = "r*Lap_s(u, ur) - u"
    TC.substitutions['Lap_t'] = "r*Lap_s(v, vr) - v"
    TC.substitutions['Lap_z'] = "Lap_s(w,wr)"
    TC.substitutions['UdotGrad_s(f, f_r)'] = "r*u*f_r + r*w*dz(f)"
    TC.substitutions['UdotGrad_r'] = "r*UdotGrad_s(u, ur) - r*v*v"
    TC.substitutions['UdotGrad_t'] = "r*UdotGrad_s(v, vr) + r*u*v"
    TC.substitutions['UdotGrad_z'] = "UdotGrad_s(w, wr)"

# adds different equations to TC object depending on whether solving 2D or 3D equations

if threeD:
    TC.add_equation("r*ur + u + dtheta(v) + r*dz(w) = 0")
else:
    TC.add_equation("r*ur + u + r*dz(w) = 0")

r_mom = "r*r*dt(u) - nu*Lap_r - 2*r*v0*v"
if threeD:
    r_mom += "+ r*v0*dtheta(u)"
r_mom += "+ r*r*dr(p) = r*v0*v0 - UdotGrad_r"
TC.add_equation(r_mom)

theta_mom = "r*r*dt(v) - nu*Lap_t + r*r*dv0dr*u + r*v0*u"
if threeD:
    theta_mom += "+ r*v0*dtheta(v) + r*dtheta(p)"
theta_mom += " = -UdotGrad_t"
TC.add_equation(theta_mom)

if threeD:
    TC.add_equation("r*r*dt(w) - nu*Lap_z + r*r*dz(p) + r*v0*dtheta(w) = -UdotGrad_z")
else:
    TC.add_equation("  r*dt(w) - nu*Lap_z +   r*dz(p)                  = -UdotGrad_z")

TC.add_equation("ur - dr(u) = 0")
TC.add_equation("vr - dr(v) = 0")
TC.add_equation("wr - dr(w) = 0")

if threeD:
    r_index = 2
else:
    r_index = 1

r = domain.grid(r_index, scales=domain.dealias)
z = domain.grid(0, scales=domain.dealias)

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
IVP.stop_sim_time = 10.*period
IVP.stop_wall_time = np.inf
IVP.stop_iteration = np.inf

p = IVP.state['p']
u = IVP.state['u']
v = IVP.state['v']
w = IVP.state['w']
ur = IVP.state['ur']
vr = IVP.state['vr']
wr = IVP.state['wr']

for f in [p,u,v,w,ur,vr,wr]:
    f.set_scales(domain.dealias, keep_data=False)

def filter_field(field,frac=0.5):
    field.require_coeff_space()
    dom = field.domain                                                                                                                                                      
    local_slice = dom.dist.coeff_layout.slices(scales=dom.dealias)                                                                                                          
    coeff = []                                                                                                                                                              
    for n in dom.global_coeff_shape:
        coeff.append(np.linspace(0,1,n,endpoint=False))                                                                                             
    cc = np.meshgrid(*coeff,indexing='ij')
    field_filter = np.zeros(dom.local_coeff_shape,dtype='bool')                                                                                                             
    for i in range(dom.dim):                                                                                                                                                
        field_filter = field_filter | (cc[i][local_slice] > frac)                                                                                                           
    field['c'][field_filter] = 0j


# incompressible perturbation, arbitrary vorticity
# u = -dz(phi)
# w = dr(phi) + phi/r
if threeD:
    A_r = domain.new_field(name='A_r')
    A_theta = domain.new_field(name='A_r')
    A_z = domain.new_field(name='A_r')
    rr = domain.new_field()
    for f in [A_r, A_theta, A_z,rr]:
        f.set_scales(domain.dealias, keep_data=False)

    rr['g']= r
    A_r['g'] = 1e-3* np.random.randn(*v['g'].shape)*np.sin(np.pi*(r - r_in))
    filter_field(A_r)
    A_theta['g'] = 1e-3* np.random.randn(*v['g'].shape)*np.sin(np.pi*(r - r_in))
    filter_field(A_theta)
    A_z['g'] = 1e-3* np.random.randn(*v['g'].shape)*np.sin(np.pi*(r - r_in))
    filter_field(A_z)
    
    u['g'] = A_z.differentiate('theta')['g']/r - A_theta.differentiate('z')['g']
    v['g'] = A_theta.differentiate('z')['g'] - A_z.differentiate('r')['g']
    scratch = (rr * A_theta).evaluate()
    w['g'] = (scratch.differentiate('r')['g']- A_r.differentiate('theta')['g'])/r
    
else:
    phi = domain.new_field(name='phi')
    phi.set_scales(domain.dealias, keep_data=False)
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
flow = flow_tools.GlobalFlowProperty(IVP, cadence=10)
flow.add_property("abs(DivU)", name='divu')

if threeD:
    CFL.add_velocities(('u', 'v', 'w'))
else:
    CFL.add_velocities(('u', 'w'))

analysis_tasks = []

analysis1 = IVP.evaluator.add_file_handler("scalar_data", iter=10)
analysis1.add_task("integ(0.5 * (u*u + v*v + w*w))", name="total kinetic energy")
analysis1.add_task("integ(0.5 * (u*u + w*w))", name="meridional kinetic energy")
analysis1.add_task("integ((u*u)**0.5)", name='u_rms')
analysis1.add_task("integ((w*w)**0.5)", name='w_rms')
analysis_tasks.append(analysis1)

# Snapshots every half an inner rotation period
analysis2 = IVP.evaluator.add_file_handler('snapshots',sim_dt=0.5*period, max_size=2**30)
analysis2.add_system(IVP.state, layout='g')
analysis_tasks.append(analysis2)

# radial profiles every 100 timestpes
analysis3 = IVP.evaluator.add_file_handler("radial_profiles", iter=100)
analysis3.add_task("integ(r*v, 'z')", name='Angular Momentum')
analysis_tasks.append(analysis3)

dt = CFL.compute_dt()
# Main loop
start_time = time.time()

while IVP.ok:
    IVP.step(dt)
    if IVP.iteration % 10 == 0:
        logger.info('Iteration: %i, Time: %e, dt: %e' %(IVP.iteration, IVP.sim_time, dt))
        logger.info('Max |divu| = {}'.format(flow.max('divu')))
    dt = CFL.compute_dt()


end_time = time.time()

# Print statistics
logger.info('Total time: %f' %(end_time-start_time))
logger.info('Iterations: %i' %IVP.iteration)
logger.info('Average timestep: %f' %(IVP.sim_time/IVP.iteration))

for task in analysis_tasks:
    logger.info(task.base_path)
    post.merge_analysis(task.base_path)
    
