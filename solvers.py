import numpy as np
from python.natconst import *
from scipy.integrate import odeint, ode
from progressbar import ProgressBar, Percentage, Bar
from copy import deepcopy
"""
    ##################################################################
    ##################################################################
"""
"""
    Solving the Rankine Hugionot conditions 2
"""
def get_RHcondition2(u1=None, t1 = None, gas=None):
    """
        Given the initial conditions of the gas
        Return the post shock conditions
    """
    #
    # Find the mach number
    #
    gam = gas._getGamma()
    M1 = gas._sumRho()*u1*u1/(gam*sum(gas.numden)*kk*t1)
    #
    # Solve the densities, velocities and temperatures
    #
    rho2 = gas._sumRho()*( (gam+1.0)*M1/( (gam-1.0)*M1+2.) )
    u2 = u1 * (gas._sumRho()/rho2)
    t2 = t1 * ( (gam - 1.) * M1 + 2.)*(2.*gam*M1 - (gam-1.))/( (gam + 1.)*
        (gam+1.)*M1)
    print 'Mach: %d'%(np.sqrt(M1))
    print 'Density: %2.5e -->  %2.5e'%(gas._sumRho()/(2.0*mp), rho2/(2.0*mp))
    print 'Temperature: %2.5e --> %2.5e'%(t1, t2)
    print 'Velocity: %2.5e --> %2.5e'%(u1, u2)
    return rho2, u2, t2, M1
""""""

"""
Solving the Rankine Hugionot conditions
"""
def get_RHcondition(rho1=None, u1=None, P1=None, gamma=3.0/2.0):
    """
    Given the initial conditions of the gas
    Return the post shock conditions
    """
    #
    # solve for velocity after the shock front
    # Quadratic equation solver: ax^2 + bx + c = 0
    # print u1, P1, rho1
    #
    gam1 = (gamma/(gamma-1.0))
    aa = (1.0/2.0) - gam1
    bb = (gam1*u1) + gam1*P1/(u1*rho1)
    cc = -(1.0/2.0 * u1*u1 + gam1 * P1/rho1)
    discrim = bb*bb - 4.*aa*cc
    #
    # check the discriminant
    #
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
def solveHD(x=None, gas=None, dust=None, numden = None, mass = None, mugas=2.8, v0 = None, grid=None, t0 = None, tdust=None, vdust=None, dv = 0.0, dT = 0.0, gamma=3.0/2.0, abserr=1e-10, telerr = 1e-10):
    """
        Call the solver with initial conditions v0 and t0
    """
    debug=False # turn this on for debugging
    #
    # TODO: inputs are generated from gas and dust properties
    # INPUTS: w0 and p
    #
    w0 = [v0, t0]
    w0 = w0 + [a*v0 for a in gas.numden]
    p = [gas, dust, debug]
    w0 = w0 + [dust.numden[0], v0, t0, dust.size[0]]
    #
    # Starting VALUES
    #
    print w0
    print ['%8.5e'%(a) for a in gas.numden]
    print
    #
    # call the solver
    #
    from shock1d import vectorfield
    #
    # Iterator
    #
    wsol = []
    #
    # Progress bar
    #
    pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=50.0).start()
    for ixrange in xrange(x.shape[0]-1):
        pbar.update((np.float(ixrange)/np.float(x.shape[0]))*50.0)
        dtnow = (x[ixrange+1]-x[ixrange]) # current step to next point
        if x[ixrange] == 0.0:
            debug = False
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
                print '  Dust Temperature: %d     K   '%(w0[8])
                print '  Dust density: %2.2e     '%(w0[7])
                print '  Dust size:  %2.3e       '%(w0[9])
            print '####################################################'
            print gas.specfrac
            print gas._sumRho(), ['%2.4e'%(a) for a in gas.numden]
            print
            rho2, u2, t2, M1 = get_RHcondition2(u1=vshock, t1=w0[1], gas=gas)
            #
            # Need to check the species fractions
            #
            numdentot = sum(gas.numden)
            gas.specfrac = [a /numdentot for a in gas.numden]
            gas._updateRho(rho=rho2, M1=M1)
            print gas.specfrac
            print rho2, gas._sumRho(), ['%2.4e'%(a) for a in gas.numden]
            #
            # After the shock
            #
            print '####################################################'
            print '  Condition after shock      '
            print '  Velocity      :   %2.2f km/s'%(u2*1e-5)
            print '  Temperature   :   %d     K   '%(t2)
            if dust.nspecs is not None:
                print '  Dust Temperature: %d     K   '%(w0[8])
                print '  Dust density: %2.2e     '%(w0[7])
                print '  Dust size:  %2.3e       '%(w0[9]), dust.size[0]
            print '####################################################'
            print
            #
            # Set the new inputs and then integrate
            #
            dt = (x[ixrange+1]-x[ixrange])/5e1
            vode = ode(vectorfield).set_integrator('vode', atol=abserr,
                rtol=telerr, order=5, method='bdf',nsteps=1e6,
                first_step = dt*1e-9, with_jacobian=True)
            w0 = [u2, t2]
            w0 = w0 + [a*u2 for a in gas.numden]
            w0 = w0 + [dust.numden[0], dust.vel[0],
                dust.temp[0], dust.size[0]]
            p = [gas,dust, debug]
            vode.set_initial_value(w0, x[ixrange]).set_f_params(p)
            wsol1 = [vode.t]+w0
            wsol.append(wsol1)
            #
            # Integrate this
            #
            iiter = 1
            while vode.successful() and (vode.t<x[ixrange+1]):
                if (vode.t + dt) > x[ixrange+1]:
                    dt = np.minimum((x[ixrange+1]-vode.t)+1e-30, dt)
                vode.integrate(vode.t+dt)
                if not vode.successful():
                    raise SystemError
                #
                # Create the new ode solver from this point
                #
                w0 = deepcopy(vode.y)
                tnow = vode.t
                vode = ode(vectorfield).set_integrator('vode', atol=abserr,
                    rtol=telerr, order=5, method='bdf',nsteps=1e6,
                    first_step = np.minimum(1e-15, dt/1e12),
                    with_jacobian=True)
                #
                # Update the dust and gas
                #
                gas._updateGas(allns=[w0[2+ispec]/w0[0] for ispec
                    in xrange(gas.nspecs)], tgas=w0[1], vgas=w0[0])
                dust._updateDust(allns=[w0[gas.nspecs+2]],
                   size=w0[gas.nspecs+5], tdust=w0[gas.nspecs+4],
                   vdust=w0[gas.nspecs+3])
                p = [gas, dust, debug]
                vode.set_initial_value(w0, tnow).set_f_params(p)
                iiter += 1
                """"""
            """"""
            print
            print 'FINISH shock front'
            print vode.y
            print
        else:
            """
                Update with adaptive stepping
            """
            dt = dtnow/1e1
            if dtnow < 500.0:
                dt = dtnow + 1e-15
            #
            # Define the steps criteria
            #
            maxstep = np.maximum(dt/1e2, 1e-5) # maximum stepping
            minstep = np.abs(x[ixrange]*1e-10) # numerical stability
            minstep = np.maximum(minstep, 1e-12)
            if x[ixrange] > 0:
                debug = False
            #
            # Setup the vode
            #
            vode = ode(vectorfield).set_integrator('vode', atol=abserr,
                rtol=telerr, order=5, method='bdf',nsteps=1e5,
                first_step = minstep, max_step=maxstep,
                with_jacobian=True)
            vode.set_initial_value(w0, x[ixrange]).set_f_params(p)
            wsol1 = [vode.t]+w0
            #
            # Save the solution
            #
            wsol.append(wsol1)
            #
            # Integrate with increasing dt starting from 5e-3
            #
            iiter = 1 # iteration stepping
            while vode.successful() and (vode.t<x[ixrange+1]):
                if (dust.numden[0] < 0.0):
                    print iiter, ixrange, vode.successful()
                    print '%e  %e'%(dt, vode.t), (vode.t < x[ixrange+1])
                    print vode.y
                    print 'NEGATIVE DUST'
                    raise SystemError
                """"""
                vode.integrate(vode.t+dt)
                ierror = 0 # iteration for error keeping
                while not vode.successful():
                    """
                        Change the step size again
                    """
                    dt = dtnow/np.float(2e2+1e1*ierror)
                    dt = np.minimum(dt, minstep)
                    p = [gas, dust, True]
                    vode.set_initial_value(w0, tnow).set_f_params(p)
                    vode.integrate(vode.t+dt)
                    ierror+=1
                    print 'Warning! %d'%(ierror)
                    print 'Error in VODE and decreasing step size'
                    print
                    if ierror > 10: raise SystemError
                """"""
                #
                # Update the gas and dust first and then create a new input
                # Create the new ode solver from this point
                #
                w0 = vode.y
                tnow = vode.t
                #
                # Check for NAN
                #
                if len(np.isnan(w0).nonzero()[0]) > 1:
                    raise SystemError
                #
                # Setup the vode
                #
                vode = ode(vectorfield).set_integrator('vode', atol=abserr,
                    rtol=telerr, order=5, method='bdf',nsteps=1e5,
                    first_step = minstep, max_step=maxstep,
                    with_jacobian=True)
                #
                # Update the dust and gas
                #
                gas._updateGas(allns=[w0[2+ispec]/w0[0] for ispec
                    in xrange(gas.nspecs)], tgas=w0[1], vgas=w0[0])
                dust._updateDust(allns=[w0[gas.nspecs+2]],
                    size=w0[gas.nspecs+5], tdust=w0[gas.nspecs+4],
                    vdust=w0[gas.nspecs+3])
                p = [gas, dust, debug]
                vode.set_initial_value(w0, tnow).set_f_params(p)
                iiter += 1
                #
                # Update the step size
                #
                if iiter > 1:
                    dt *= np.float(iiter)
                    #
                    # Define the steps criteria
                    #
                    maxstep = dt/5e1 # maximum stepping
                    if (x[ixrange+1]-vode.t) < 1e1:
                        dt = (x[ixrange+1]-vode.t) + 1e-15
                        maxstep = dt/5e1
            """"""
        """"""
        #
        # Uncomment here to debug
        #
#        print vode.y[1], '%2.5e'%(gas._calculateR(t=vode.y[1])),
#        print '%2.6e'%(gas._calculateFreeEnergy(t=vode.y[1]))
        #
        # Update the w0 and x0
        #
        w0 = [a for a in vode.y]
    """"""
    wsol1 = [vode.t]+w0
    wsol.append(wsol1)
    return np.array(wsol), vshock
""""""
"""
HDrt Solver
"""
def solveHDrt(x=None, gas=None, dust=None, numden = None, mass = None, mugas=2.8, v0 = None, t0 = None, tdust=None, vdust=None, Jrad=None, gamma=3.0/2.0, abserr=1e-12, telerr = 1e-10):
    """
        Call the solver with initial conditions v0 and t0
    """
    debug   = False # turn this on for debugging
    w0      = [v0, t0]
    w0      = w0 + [a*v0 for a in gas.numden]
    p       = [gas, dust, Jrad, debug]
    w0      = w0 + [dust.numden[0], v0, t0, dust.size[0]]
    #
    # Starting VALUES
    #
    print w0
    print ['%8.5e'%(a) for a in gas.numden]
    print
    #
    # call the solver
    #
    from shock1d import vectorfieldrt
    #
    # Iterator
    #
    wsol = []
    #
    # Progress bar
    #
    pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=50.0).start()
    for ixrange in xrange(x.shape[0]-1):
        pbar.update((np.float(ixrange)/np.float(x.shape[0]))*50.0)
        dtnow = (x[ixrange+1]-x[ixrange]) # current step to next point
        if x[ixrange] == 0.0:
            debug = False
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
                print '  Dust Temperature: %d     K   '%(w0[8])
                print '  Dust density: %2.2e     '%(w0[7])
                print '  Dust size:  %2.3e       '%(w0[9])
            print '####################################################'
            print gas.specfrac
            print gas._sumRho(), ['%2.4e'%(a) for a in gas.numden]
            print
            rho2, u2, t2, M1 = get_RHcondition2(u1=vshock, t1=w0[1], gas=gas)
            #
            # Need to check the species fractions
            #
            numdentot = sum(gas.numden)
            gas.specfrac = [a /numdentot for a in gas.numden]
            gas._updateRho(rho=rho2, M1=M1)
            print gas.specfrac
            print rho2, gas._sumRho(), ['%2.4e'%(a) for a in gas.numden]
            #
            # After the shock
            #
            print '####################################################'
            print '  Condition after shock      '
            print '  Velocity      :   %2.2f km/s'%(u2*1e-5)
            print '  Temperature   :   %d     K   '%(t2)
            if dust.nspecs is not None:
                print '  Dust Temperature: %d     K   '%(w0[8])
                print '  Dust density: %2.2e     '%(w0[7])
                print '  Dust size:  %2.3e       '%(w0[9]), dust.size[0]
            print '####################################################'
            print
            #
            # Set the new inputs and then integrate
            #
            dtnow = (x[ixrange+1]-x[ixrange])
            dt = dtnow/5e1
            vode = ode(vectorfieldrt).set_integrator('vode', atol=abserr,
                rtol=telerr, order=5, method='bdf',nsteps=1e6,
                first_step = dt*1e-9, with_jacobian=True)
            w0 = [u2, t2]
            w0 = w0 + [a*u2 for a in gas.numden]
            w0 = w0 + [dust.numden[0], dust.vel[0],
                dust.temp[0], dust.size[0]]
            p = [gas,dust, Jrad, debug]
            vode.set_initial_value(w0, x[ixrange]).set_f_params(p)
            wsol1 = [vode.t]+w0
            wsol.append(wsol1)
            #
            # Integrate this
            #
            iiter = 1
            while vode.successful() and (vode.t<x[ixrange+1]):
                if (vode.t + dt) > x[ixrange+1]:
                    dt = np.minimum((x[ixrange+1]-vode.t)+1e-10, dt)
                vode.integrate(vode.t+dt)
                if not vode.successful():
                    raise SystemError
                #
                # Create the new ode solver from this point
                #
                w0 = deepcopy(vode.y)
                tnow = vode.t
                vode = ode(vectorfieldrt).set_integrator('vode', atol=abserr,
                    rtol=telerr, order=5, method='bdf',nsteps=1e6,
                    first_step = np.minimum(1e-15, dt/1e12),
                    with_jacobian=True)
                #
                # Update the dust and gas
                #
                gas._updateGas(allns=[w0[2+ispec]/w0[0] for ispec
                    in xrange(gas.nspecs)], tgas=w0[1], vgas=w0[0])
                dust._updateDust(allns=[w0[gas.nspecs+2]],
                    size=w0[gas.nspecs+5], tdust=w0[gas.nspecs+4],
                    vdust=w0[gas.nspecs+3])
                p = [gas, dust, Jrad, debug]
                vode.set_initial_value(w0, tnow).set_f_params(p)
                iiter += 1
                dt = np.minimum(dt*np.float(iiter), dtnow/3.)
                if iiter > 1000:
                    raise SystemError
                """"""
            """"""
            print
            print 'FINISH shock front'
            print vode.y
            print
        else:
            dt = dtnow/5e1
            if dtnow < 1e1:
                dt = dtnow
            #
            # Define the steps criteria
            #
            maxstep = np.maximum(dt/1e2, 1e-5)
            minstep = np.abs(x[ixrange]*1e-10) # numerical stability
            minstep = np.maximum(minstep, 1e-12)
            #
            # Setup the vode
            #
            vode = ode(vectorfieldrt).set_integrator('vode', atol=abserr,
                rtol=telerr, order=5, method='bdf',nsteps=1e6,
                first_step = minstep, max_step=maxstep,
                with_jacobian=True)
            vode.set_initial_value(w0, x[ixrange]).set_f_params(p)
            wsol1 = [vode.t]+w0
            wsol.append(wsol1)
            #
            # Integrate with increasing dt starting from 5e-3
            #
            iiter = 1 # iteration stepping
            while vode.successful() and (vode.t<x[ixrange+1]):
                if (dust.numden[0] < 0.0):
                    print iiter, ixrange, vode.successful()
                    print '%e  %e'%(dt, vode.t), (vode.t < x[ixrange+1])
                    print vode.y
                    print 'NEGATIVE DUST'
                    raise SystemError
                """"""
                vode.integrate(vode.t+dt)
                tnow = vode.t
                ierror = 0 # iteration for error keeping
                while not vode.successful():
                    """
                    Change the step size again
                    """
                    dt = dtnow/np.float(3e2+1e1*ierror)
                    dt = np.minimum(dt, minstep)
                    p = [gas, dust, Jrad, True]
                    vode.set_initial_value(w0, tnow).set_f_params(p)
                    vode.integrate(vode.t+dt)
                    ierror+=1
                    print 'Warning! %d'%(ierror), ixrange
                    print 'Error in VODE and decreasing step size'
                    print
                    if ierror > 10: raise SystemError
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
                #
                # Setup the vode
                #
                vode = ode(vectorfieldrt).set_integrator('vode', atol=abserr,
                    rtol=telerr, order=5, method='bdf',nsteps=1e5,
                    first_step = minstep, max_step=maxstep,
                    with_jacobian=True)
                #
                # Update the dust and gas
                #
                gas._updateGas(allns=[w0[2+ispec]/w0[0] for ispec
                    in xrange(gas.nspecs)], tgas=w0[1], vgas=w0[0])
                dust._updateDust(allns=[w0[gas.nspecs+2]],
                    size=w0[gas.nspecs+5], tdust=w0[gas.nspecs+4],
                    vdust=w0[gas.nspecs+3])
                p = [gas, dust, Jrad, debug]
                vode.set_initial_value(w0, tnow).set_f_params(p)
                iiter += 1
                #
                # Update the step size
                #
                if iiter > 1:
                    dt *= np.float(iiter)
                    #
                    # Define the steps criteria
                    #
                    maxstep = dt/1e2 # maximum stepping
                    if (x[ixrange+1]-vode.t) < 1e1:
                        dt = (x[ixrange+1]-vode.t) + 1e-15
                        maxstep = dt/1e2
                """"""
            """"""
        #
        # Update the w0 and x0
        #
        w0 = [a for a in vode.y]
    """"""
    wsol1 = [vode.t]+w0
    wsol.append(wsol1)
    return np.array(wsol), vshock
""""""
