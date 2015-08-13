import numpy as np
from python.natconst import *
from scipy.integrate import odeint, ode
from progressbar import ProgressBar, Percentage, Bar
from copy import deepcopy
"""
    ##################################################################
    ##################################################################
"""
#
# Global absolute error
#
glbabserr = [1.0, 1e-2, 0.1, 0.1, 0.1, 0.1, 1e-12, 1.0, 1e-2, 1e-8]
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
def solveHD(x=None, gas=None, dust=None, v0 = None, t0 = None, haveJ = False, Jrad = None, abserr=1e-10, telerr = 1e-10):
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
    if haveJ:
        p  = [gas, dust, Jrad, True, debug]
    else:
        p  = [gas, dust, None, False, debug]
    w0 = w0 + [dust.numden[0], v0, t0, dust.size[0]]
    #
    # Starting VALUES
    #
    print ['%8.5e'%(a) for a  in w0]
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
    gasKappa = [] # need to store kappa
    #
    # Progress bar
    #
    pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=50.0).start()
    for ixrange in xrange(x.shape[0]-1):
        pbar.update((np.float(ixrange)/np.float(x.shape[0]))*50.0)
        dtnow = (x[ixrange+1]-x[ixrange]) # current step to next point
        if x[ixrange] == 0.0:
            debug = False
            vshock = w0[0]
            #
            # Solve the shock front
            #
            rho2, u2, t2, M1 = get_RHcondition2(u1=vshock, t1=w0[1], gas=gas)
            #
            # Need to check the species fractions
            #
            numdentot = sum(gas.numden)
            gas.specfrac = [a /numdentot for a in gas.numden]
            gas._updateRho(rho=rho2, M1=M1)
            #
            # Set the new inputs and then integrate
            #
            dt = (x[ixrange+1]-x[ixrange])/1e1
            vode = ode(vectorfield).set_integrator('vode', atol=glbabserr,
                rtol=telerr, order=5, method='bdf',nsteps=1e3,
                first_step = dt*1e-3, with_jacobian=True)
            w0 = [u2, t2]
            w0 = w0 + [a*u2 for a in gas.numden]
            w0 = w0 + [dust.numden[0], dust.vel[0],
                dust.temp[0], dust.size[0]]
            if haveJ:
                p  = [gas, dust, Jrad, True, debug]
            else:
                p  = [gas, dust, None, False, debug]
            vode.set_initial_value(w0, x[ixrange]).set_f_params(p)
            wsol1 = [vode.t]+w0
            #
            # Save the solution
            #
            wsol.append(wsol1)
            gasKappa.append(gas._getKap(w0[1], destroyed=dust.destroyed))
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
                    rtol=telerr, order=5, method='bdf',nsteps=1e4,
                    first_step = np.maximum(1e-15, dt/1e3),
                    with_jacobian=True)
                #
                # Update the dust and gas
                #
                gas._updateGas(allns=[w0[2+ispec]/w0[0] for ispec
                    in xrange(gas.nspecs)], tgas=w0[1], vgas=w0[0])
                dust._updateDust(allns=[w0[gas.nspecs+2]],
                    size=w0[gas.nspecs+5], tdust=w0[gas.nspecs+4],
                    vdust=w0[gas.nspecs+3])
                if haveJ:
                    p  = [gas, dust, Jrad, True, debug]
                else:
                    p  = [gas, dust, None, False, debug]
                vode.set_initial_value(w0, tnow).set_f_params(p)
                iiter += 1
                """"""
            """"""
        else:
            """
                Update with adaptive stepping
            """
            dt = dtnow
            #
            # Setup the vode
            #
            vode = ode(vectorfield).set_integrator('vode', atol=glbabserr,
                rtol=telerr, order=5, method='bdf',nsteps=1e5,
                with_jacobian=True, min_step=dt*1e-5, first_step=dt*1e-8)
            tnow = x[ixrange]
            vode.set_initial_value(w0, x[ixrange]).set_f_params(p)
            wsol1 = [vode.t]+w0
            #
            # Save the solution
            #
            wsol.append(wsol1)
            gasKappa.append(gas._getKap(w0[1], destroyed=dust.destroyed))
            #
            # Integrate with increasing dt starting from 5e-3
            #
            iiter = 1 # iteration stepping
            while vode.successful() and (vode.t<x[ixrange+1]):
                vode.integrate(vode.t+dt)
                ierror = 0 # iteration for error keeping
                while not vode.successful():
                    """
                        Change the step size again
                    """
                    print '%2.5e  %2.5e'%(vode.t, dt)
                    dt = dt/np.float(10.+5.*ierror)
                    if haveJ:
                        p  = [gas, dust, Jrad, True, debug]
                    else:
                        p  = [gas, dust, None, False, debug]
                    vode = ode(vectorfield).set_integrator('vode',
                        atol=glbabserr, rtol=telerr, order=5, method='bdf',
                        nsteps=1e5, with_jacobian=True, min_step=1e-5,
                        first_step=1e-3)
                    vode.set_initial_value(w0, tnow).set_f_params(p)
                    vode.integrate(vode.t+dt)
                    ierror+=1
                    print 'Warning! %d'%(ierror), vode.successful()
                    print 'Error in VODE and decreasing step size'
                    print '%2.3e  %2.3e  %2.3e %2.3e %2.3e %2.3e'%(
                        vode.t, tnow, dt, dtnow, x[ixrange], x[ixrange+1])
                    print '%2.5e  %2.5e'%(vode.t, dt)
                    print vode.y
                    print
                    if ierror > 15: raise SystemError
                """"""
                #
                # Update the gas and dust first and then create a new input
                # Create the new ode solver from this point
                #
                w0      = deepcopy(vode.y)
                tnow    = vode.t
                #
                # Check for NAN
                #
                if len(np.isnan(w0).nonzero()[0]) > 1:
                    raise SystemError
                #
                # Setup the vode
                #
                vode = ode(vectorfield).set_integrator('vode', atol=glbabserr,
                    rtol=telerr, order=5, method='bdf',nsteps=1e5,
                    with_jacobian=True)
                #
                # Update the dust and gas
                #
                gas._updateGas(allns=[w0[2+ispec]/w0[0] for ispec
                    in xrange(gas.nspecs)], tgas=w0[1], vgas=w0[0])
                dust._updateDust(allns=[w0[gas.nspecs+2]],
                    size=w0[gas.nspecs+5], tdust=w0[gas.nspecs+4],
                    vdust=w0[gas.nspecs+3])
                #
                # Check if dust is destroyed
                #
                if (not dust.destroyed) and (w0[-2] > 2e3):
                    dust.destroyed = True
                """"""
                if haveJ:
                    p  = [gas, dust, Jrad, True, debug]
                else:
                    p  = [gas, dust, None, False, debug]
                vode.set_initial_value(w0, tnow).set_f_params(p)
                iiter += 1
                #
                # Change the step size
                #
                dt = (x[ixrange+1] - vode.t) + 1e-15
            """"""
        """"""
        #
        # Update the w0 and x0
        #
        w0 = [a for a in vode.y]
        #
        # DEBUG HERE
        #
#        print '%2.5e'%(x[ixrange]), ['%2.3e'%(a) for a in w0], gas._sumRho(), dust._sumRho()
    """"""
    wsol1 = [vode.t]+w0
    #
    # Save the solution
    #
    wsol.append(wsol1)
    gasKappa.append(gas._getKap(w0[1], destroyed=dust.destroyed))
    return np.array(wsol), vshock, gasKappa
""""""
