import numpy as np
from python.natconst import *
from scipy.integrate import odeint, ode

"""
Solving the Rankine Hugionot conditions
"""
def get_RHcondition(rho1=None, u1=None, P1=None, gamma=3.0/2.0):
    """
    Given the initial conditions of the gas
    Return the post shock conditions
    """
    # solve for velocity after the shock front
    # Quadratic equation solver: ax^2 + bx + c = 0
    print u1, P1, rho1
    gam1 = (gamma/(gamma-1.0))
    aa = (1.0/2.0) - gam1
    bb = (gam1*u1) + gam1*P1/(u1*rho1)
    cc = -(1.0/2.0 * u1*u1 + gam1 * P1/rho1)
    
    discrim = bb*bb - 4.*aa*cc
    
    # check the discriminant
    if discrim < 0:
        print 'This has no solution'
        raise ValueError
    elif discrim == 0:
        u2 = (-bb + np.sqrt(bb*bb - 4.0*aa*cc))/(2.0*aa)
    else:
        u21 = (-bb + np.sqrt(bb*bb - 4.0*aa*cc))/(2.0*aa)
        u22 = (-bb - np.sqrt(bb*bb - 4.0*aa*cc))/(2.0*aa)
        # Solve the pressure
        P21 = rho1*u1*u1 + P1 - rho1*u1*u21
        P22 = rho1*u1*u1 + P1 - rho1*u1*u22
        if (P21 > P1):
            u2 = u21
            P2 = P21
        elif (P22 > P1):
            u2 = u22
            P2 = P22
        else:
            print 'No shock solutions'
            print u1, P1, rho1
            print u21, u22
            raise SystemError
    """"""
    
    # Solve rho2
    rho2 = rho1 * u1/u2
    print 'Density: %2.5e -->  %2.5e'%(rho1/(2.0*mp), rho2/(2.0*mp))
    print 'Pressure: %2.5e --> %2.5e'%(P1, P2)
    print 'Velocity: %2.5e --> %2.5e'%(u1, u2)
    
    return P2, rho2, u2
""""""

"""
HD Solver
"""
def solveHD(x=None, gas=None, dust=None, numden = None, mass = None, mugas=2.8, v0 = None, grid=None, t0 = None, tdust=None, vdust=None, dv = 0.0, dT = 0.0, gamma=3.0/2.0, abserr=1e-3, telerr = 1e-6):
    """
    Call the solver with initial conditions v0 and t0
    """
    w0 = [v0, t0]
    if gas is None:
        p = [mass, mugas, gamma]
    else:
        w0 = w0 + [a*v0 for a in gas.numden]
        p = [gas]
        if dust is not None:
            p = p + [dust]
            if tdust is None:
                w0 = w0 + [dust.numden[0], v0, t0, dust.size[0]]
            else:
                w0 = w0 + [dust.numden[0], vdust, tdust, dust.size[0]]
    """"""
    # call the solver
    from shock1d import vectorfield
    #
    # Iterator
    #
    for ixrange in xrange(x.shape[0]):
        vode = ode(vectorfield).set_integrator('lsoda', atol=abserr,
            rtol=telerr, nsteps=1e4)
        vode.set_initial_value(w0, -1.).set_f_params(p)
        wsol = []
        wsol.append(w0)
        #
        # Integrate with increasing dt starting from 5e-3
        #
        iiter = 1
        dt = (x[ixrange+1]-x[ixrange])
        print '%e'%(dt)
        while vode.successful():
            dt = 1e-3
            vode.integrate(vode.t+dt)
            print iiter, vode.successful(), dt, vode.t
            print vode.y
            print 
            if not vode.successful():
                raise SystemError
            iiter += 1
            w0 = vode.y
            tnow = vode.t
            vode = ode(vectorfield).set_integrator('lsoda', atol=abserr,
                                               rtol=telerr, nsteps=1e4)
            vode.set_initial_value(w0, tnow).set_f_params(p)
        """"""
        print iiter, dt
        print '%4.8e %4.8e'%(x[ixrange], vode.t)
        raise SystemError
        dt = (x[ixrange+1]-x[ixrange])
        #
        # Update the w0 and x0
        #
        w0 = vode.y
    """"""
    return wsol
""""""
"""
HDrt Solver
"""
def solveHDrt(x=None, gas=None, dust=None, numden = None, mass = None, mugas=2.8, v0 = None, t0 = None, tdust=None, vdust=None, Jrad=None, gamma=3.0/2.0, abserr=1e-5, telerr = 1e-5):
    """
    Call the solver with initial conditions v0 and t0
    """
    w0 = [v0, t0]
    if gas is None:
        p = [mass, mugas, gamma, Jrad]
    else:
        w0 = w0 + [a*v0 for a in gas.numden]
        p = [gas, Jrad]
        if dust is not None:
            p = p + [dust]
            if tdust is None:
                w0 = w0 + [dust.numden[0], v0, t0, dust.size[0]]
            else:
                w0 = w0 + [dust.numden[0], vdust, tdust, dust.size[0]]
            """"""
    """"""
    # call the solver
    from shock1d import vectorfieldrt
    #wsol = odeint(vectorfieldrt, w0, x, args=(p,),
                  #atol=abserr, rtol=telerr, mxstep=50000)
    vode = ode(vectorfieldrt).set_integrator('lsoda', method='bdf',
        with_jacobian=False, atol=abserr, rtol=telerr, nsteps=30000)
    vode.set_initial_value(w0, x[0]).set_f_params(p)
    dx = x[1:]-x[:-1]
    wsol = []
    wsol.append(w0)
    #
    # Iterator
    #
    indx = 0
    while (vode.successful()) and (vode.t < x[-1]):
        dt = np.minimum(1e4, dx[indx]/1e3)
        vode.integrate(vode.t+dt)
        if vode.t >= x[indx+1]:
            wsol.append(vode.y)
            indx += 1
    """"""
    return wsol
""""""
