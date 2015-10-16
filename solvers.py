import numpy as np
from python.natconst import *
from scipy.integrate import odeint, ode
from progressbar import ProgressBar, Percentage, Bar
from copy import deepcopy
#
# call the solver
#
#from shock1d import vectorfield
import pyximport; pyximport.install(setup_args={
    "include_dirs":np.get_include()},
    reload_support = True)
from cshock1d import solve
"""
    ##################################################################
    ##################################################################
"""
#
# Global absolute error
#
glbabserr = [1e-2, 1e-5, 1e-4, 1e-4, 1e-4, 1e-4, 1e-25, 1e-2, 1e-5, 1e-20]
def createInputs(v0, t0, gas, dust):
    """
    Function to create input parameters
    """
    #
    # Input structure: vgas, tgas, nspec, ndust1, vdust1, tdust1, adust1,...
    #
    w = [v0, t0] # initial gas velocity and temperature
    w += [a for a in gas.numden]
    for idust in range(dust.nspecs):
        w += [dust.numden[idust]]
        w += [v0]
        w += [t0]
        w += [dust.size[idust]]
    """"""
    return w
""""""
"""
    Solving the Rankine Hugionot conditions 2
"""
def get_RHcondition2(u1=None, par=None, gas=None):
    """
        Given the initial conditions of the gas
        Return the post shock conditions
    """
    #
    # Find the mach number
    #
    rhogas  = sum([a*b for (a,b) in zip(par[2:gas.nspecs+2]/par[0], gas.mass)])
    Ngas    = sum(par[2:gas.nspecs+2]/par[0])
    mbar    = rhogas/Ngas
    gam = sum([a * (b * 2. - 2.) for (a,b) in
        zip(par[2:gas.nspecs+2], gas.gamma)])
    gam = (gam/sum(par[2:gas.nspecs+2]) + 2.)/(gam/sum(par[2:gas.nspecs+2]))
    M1 = mbar/gam * u1*u1 /(kk*par[1])
    #
    # Solve the densities, velocities and temperatures
    #
    rho2 = rhogas*( (gam+1.0)*M1/( (gam-1.0)*M1+2.) )
    u2 = u1 * ((gam-1.)*M1+2.)/( (gam+1.)*M1)
    t2 = par[1] * (2.*gam*M1 - gam +1.)*( (gam-1.)*M1 + 2.)/( (gam + 1.)*
        (gam + 1.) * M1)
    #print
    #print 'Mach: %d'%(np.sqrt(M1)), par[1], u1, mbar, gam
    #print 'Density: %2.5e -->  %2.5e'%(rhogas/(2.0*mp), rho2/(2.0*mp))
    #print 'Temperature: %2.5e --> %2.5e'%(par[1], t2)
    #print 'Velocity: %2.5e --> %2.5e'%(u1, u2)
    return rho2, u2, t2, M1
""""""
def checkDust(w, gas, dust):
    if (not dust.destroyed and (w[1] > 1.4e3)):
        dust.destroyed=True
        #
        # Add the components
        #
        Rhogas = sum([a*b for (a,b) in zip(w[2:2+gas.nspecs+1], gas.mass)])
        for ispec in xrange(gas.nspecs):
            w[2+ispec]  += (gas.dustfrac[ispec] * dust.dustfrac*(
                1.-dust.chonfrac)*Rhogas/gas.mass[ispec])
"""
HD Solver
"""
def solveHD(x=None, gas=None, dust=None, v0 = None, t0 = None, haveJ = False, Jrad = None, abserr=1e-10, telerr = 1e-10):
    """
        Call the solver with initial conditions v0 and t0
    """
    debug=False # turn this on for debugging
    #
    # INPUTS: w0 and p
    #
    w0 = createInputs(v0,t0,gas,dust)
    if haveJ:
        p  = [gas, dust, Jrad, True, debug]
    else:
        p  = [gas, dust, np.zeros((20,2)), False, debug]
    #    #
    #    # Starting VALUES
    #    #
    #    print 'Starting VALUES: '
    #    print ['%8.5e'%(a) for a  in w0]
    #    print
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
            w0 = [a for a in vode.y]
            rho2, u2, t2, M1 = get_RHcondition2(u1=vshock, par = w0, gas=gas)
            #
            # Need to check the species fractions
            #
            gam = sum([a * (b * 2. - 2.) for (a,b) in
                zip(w0[2:gas.nspecs+2], gas.gamma)])
            gam = (gam/sum(w0[2:gas.nspecs+2]) + 2.)/(gam/
                sum(w0[2:gas.nspecs+2]))
            #
            # The new inputs
            # Calculate the new number densities
            #
            for iinps in xrange(gas.nspecs+2):
                if iinps    == 0:
                    w0[iinps] = u2
                elif iinps  == 1:
                    w0[iinps] = t2
                else:
                    w0[iinps] *= ((gam+1.)*M1 / ( (gam-1.)*M1 + 2.))
                """"""
            """"""
            #
            # Set the new inputs and then integrate
            #
            dt = (x[ixrange+1]-x[ixrange])/1e1
            vode = ode(solve).set_integrator('vode', atol=glbabserr,
                rtol=telerr, order=5, method='bdf',nsteps=1e3,
                first_step = dt*1e-3, with_jacobian=True)
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
                #
                # Check if dust is destroyed
                #
                checkDust(w0, gas, dust)
                tnow = vode.t
                vode = ode(solve).set_integrator('vode', atol=abserr,
                    rtol=telerr, order=5, method='bdf',nsteps=1e4,
                    first_step = np.maximum(1e-15, dt/1e3),
                    with_jacobian=True)
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
            vode = ode(solve).set_integrator('vode', atol=glbabserr,
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
                    vode = ode(solve).set_integrator('vode',
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
                #
                # Check if dust is destroyed
                #
                checkDust(w0, gas, dust)
                tnow    = vode.t
                #
                # Check for NAN
                #
                if len(np.isnan(w0).nonzero()[0]) > 1:
                    raise SystemError
                #
                # Setup the vode
                #
                vode = ode(solve).set_integrator('vode', atol=glbabserr,
                    rtol=telerr, order=5, method='bdf',nsteps=1e5,
                    with_jacobian=True)
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
        # Solver debug
        #
        #if ixrange < 2:
        #    print '%d  %2.5e'%(ixrange, x[ixrange]), solve(vode.t, vode.y, p)
        #if x[ixrange] >0:
        #    print '%d  %2.5e'%(ixrange, x[ixrange]), solve(vode.t, vode.y, p)
        #    raise SystemExit
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
