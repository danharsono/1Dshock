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
def solveHD(x=None, gas=None, dust=None, numden = None, mass = None, mugas=2.8, v0 = None, grid=None, t0 = None, tdust=None, vdust=None, dv = 0.0, dT = 0.0, gamma=3.0/2.0, abserr=1e-12, telerr = 1e-6):
    """
    Call the solver with initial conditions v0 and t0
    """
    w0 = [v0, t0]
    w0 = w0 + [a*v0 for a in gas.numden]
    p = [gas, dust, False]
    w0 = w0 + [dust.numden[0], v0, t0, dust.size[0]]
    #
    # call the solver
    #
    from shock1d import vectorfield
    #
    # Iterator
    #
    wsol = []
    for ixrange in xrange(x.shape[0]-1):
        if x[ixrange] == 0.0:
            print
            print 'This is the shock front'
            print x[ixrange]
            vshock = w0[0]
            #
            # Solve the shock front
            #
            print '####################################################'
            print '  Condition before shock      '
            print '  velocity:   %2.2f  km/s'%(vshock*1e-5)
            print '  Pre shock T   :   %d     K   '%(w0[1])
            if dust.nspecs is not None:
                print '  Dust Temperature: %d     K   '%(w0[7])
                print '  Dust density: %2.2e     '%(w0[6])
            print '####################################################'
            print
            P1 = (kk/(gas.mugas*mp) * gas._sumRho() *
                  w0[1])
            P2, rho2, u2 = get_RHcondition(rho1 = gas._sumRho(),
                u1 = w0[0], P1=P1)
            t2 = P2*(gas.mugas*mp)/(kk*rho2)
            #
            # Need to check the species fractions
            #
            numdentot = sum(gas.numden)
            gas.specfrac = gas.numden/numdentot
            gas._updateRho(rho=rho2)
            #
            # Set the new inputs and then integrate
            #
            dt = (x[ixrange+1]-x[ixrange])/1e3
            vode = ode(vectorfield).set_integrator('vode', atol=1e-4,
                rtol=1e-4, order=15, method='bdf',nsteps=1e6,
                first_step = dt*1e-9, with_jacobian=True)
            w0 = [u2, t2]
            w0 = w0 + [a*u2 for a in gas.numden]
            w0 = w0 + [dust.numden[0], dust.vel[0], dust.temp[0], dust.size[0]]
            p = [gas,dust, False]
            vode.set_initial_value(w0, x[ixrange]).set_f_params(p)
            wsol1 = [vode.t]+w0
            wsol.append(wsol1)
            #
            # Integrate this
            #
            iiter = 1
            while vode.successful() and (vode.t<x[ixrange+1]):
                vode.integrate(vode.t+dt)
                if not vode.successful():
                    raise SystemError
                #
                # Create the new ode solver from this point
                #
                w0 = vode.y
                tnow = vode.t
                vode = ode(vectorfield).set_integrator('vode', atol=abserr,
                    rtol=telerr, order=15, method='bdf',nsteps=1e6,
                    first_step = np.minimum(2e-8, dt/1e9),
                    with_jacobian=True)
                #
                # Update the dust and gas
                #
                gas._updateGas(allns=[w0[2]/w0[0], w0[3]/w0[0], w0[4]/w0[0]])
                dust._updateDust(allns=[w0[5]], size=w0[8])
                p = [gas, dust, False]
                vode.set_initial_value(w0, tnow).set_f_params(p)
                iiter += 1
                """"""
            """"""
            print
            print 'FINISH shock front'
            print
        else:
            dt = (x[ixrange+1]-x[ixrange])/1e5
            vode = ode(vectorfield).set_integrator('vode', atol=abserr,
                rtol=telerr, order=5, method='bdf',nsteps=5e3,
                first_step = np.minimum(dt, 2e-4), max_step=dt/5.,
                with_jacobian=True)
            vode.set_initial_value(w0, x[ixrange]).set_f_params(p)
            wsol1 = [vode.t]+w0
            wsol.append(wsol1)
            #
            # Integrate with increasing dt starting from 5e-3
            #
            iiter = 1
            while vode.successful() and (vode.t<x[ixrange+1]):
                dt *= np.float(iiter)**(1.0/np.float(iiter))
                if (vode.t + dt) > x[ixrange+1]:
                    dt = np.minimum((x[ixrange+1]-vode.t)+1e-3, dt)
                vode.integrate(vode.t+dt)
#                print iiter, ixrange, vode.successful()
#                print '%e  %e'%(dt, vode.t), (vode.t < x[ixrange+1])
#                print vode.y
#                print
                if not vode.successful():
                    raise SystemError
                    """
                        Switch to lsoda
                    """
                    vode = ode(vectorfield).set_integrator('lsoda', atol=abserr,
                        rtol=telerr, nsteps=1e4, max_order_ns=10,
                        first_step = np.minimum(dt, 2e-4), max_step=dt/1e2)
                    vode.set_initial_value(w0, tnow).set_f_params(p)
                    vode.integrate(vode.t+dt)
                """"""
                #
                # Create the new ode solver from this point
                #
                w0 = vode.y
                tnow = vode.t
                #
                # Check for NAN
                #
                if len(np.isnan(w0).nonzero()[0]) > 1:
                    raise SystemError
                vode = ode(vectorfield).set_integrator('vode', atol=abserr,
                    rtol=telerr, order=5, method='bdf',nsteps=5e3,
                    first_step = np.minimum(dt, 2e-4), max_step=dt/5.,
                    with_jacobian=True)
                #
                # Update the dust and gas
                #
                gas._updateGas(allns=[w0[2]/w0[0], w0[3]/w0[0], w0[4]/w0[0]])
                dust._updateDust(allns=[w0[5]], size=w0[8])
                p = [gas, dust, False]
                vode.set_initial_value(w0, tnow).set_f_params(p)
                iiter += 1
            """"""
        """"""
        #
        # Update the w0 and x0
        #
        w0 = [a for a in vode.y]
    """"""
    wsol1 = [vode.t]+w0
    wsol.append(wsol1)
    return np.array(wsol)
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
