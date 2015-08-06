import numpy as np
from python.natconst import *
from solvers import get_RHcondition, solveHD, solveHDrt
from oneDRT import calc_tau, calcJrad, calc_tauall
from gasspec import gasSpecs
from dustspec import dustSpecs
import matplotlib.pyplot as plt
from progressbar import ProgressBar, Percentage, Bar


"""
The main part of the shock code
"""

def shock_main(numden=1e14, rhogas=1e-9, nspecs=None, ndust=None, adust=300e-4, mass=2.0*mp, v0=6e5, t0=300., sizex=5.0, numpoints=1e5, mugas=2.8, dustfrac=0.005, mdust=3.3, ncpu=3, niter=1, restart = False):
    """
    The main part of the shock code
    Here:
        - Solve the pre shock conditions
        - Solve the post shock conditions
        - Merge them and make sure they are consistent
        - iterates the solutions
    """
    """
    Calculate the number densities of the different species
    use the rhogas
    """
    gas = gasSpecs(rho=rhogas, nspecs=nspecs)
    print
    print 'Using gas Species: %d species'%(nspecs)    
    """
    Generate the dust number densities
    """
    dust = dustSpecs(dustfrac=dustfrac, nspecs=ndust, 
        mdust=mdust, gas=gas, size=adust)
    print 'Using dust Speces: %d sizes'%(ndust)
    print 'Initial size: %1.2f micron'%(dust.size[0]*1e4)
    print
    #
    # Check whether there are Nones in the input
    #
    print 
    print '####################################################'
    print '  Input parameters       '
    print '  Shock velocity:   %2.2f  km/s'%(v0/1e5)
    print '  Pre shock T   :   %d     K   '%(t0)
    print '  Pre-shock n   : 1E%2.1f  cm-3'%(np.log10(gas._sumRho()))
    print '  Dust densities: %2.3e   cm-3'%(dust.numden[0])
    print '####################################################'
    print 
    #
    # Solve the whole shock
    # Pre-shock part
    #
    xpre = np.logspace(np.log10(sizex), -8, numpoints)
    xpre = -xpre[:-1]
    #
    # Add the shock front
    #
    xpre = np.concatenate((xpre, np.array([0.0])))
    #
    # Add the post shock
    #
    xpost = -xpre[:-1]
    xpre = np.concatenate((xpre, xpost[::-1]))
    dx = xpre[1:]-xpre[:-1]
    #
    # Solve this
    #
    if not restart:
        wpre, vshock = solveHD(x=xpre, gas=gas, dust=dust, v0=v0, t0=t0,
            abserr=1e-6, telerr=1e-8)
    sol0 = wpre
    print
    print '####################################################'
    print '  Condition at the end      '
    print '  velocity:   %2.2f  km/s'%(sol0[-2][1]/1e5)
    print '  Gas temperature   :   %d     K   '%(sol0[-2][2])
    if ndust is not None:
        print '  Dust Temperature: %d     K   '%(sol0[-2][9])
        print '  Dust Densities: %1.3e   cm^_3'%(sol0[-1][7])
        print '  Dust size     :   %1.2f microns'%(sol0[-2][10])
    print '####################################################'
    print
    print wpre[-1,:]
    #
    Tpre    = sol0[0,-2]
    Tpost   = np.minimum(sol0[-1,2], sol0[-1,-2])
    """
    Initialize the radiative transfer
      - Calculate the Jrad
      - Caluclate tau
    """
    print
    print 'Solving radiation field...'
    tau, srcall = calc_tauall(sol=sol0, gas=gas,
        dust=dust)
    #
    # Vectorize method as of Jul 2015
    #
    Jrad, Frad = calcJrad(Tpre=Tpre, Tpost=Tpost, srcall=srcall, tau=tau)
    print 'Tpre: %d  --  Tpost: %d'%(Tpre, Tpost)
    print 'Frad: %2.3e -- %2.3e'%(Frad[0], Frad[-1])
    print
    Jrad = np.array([sol0[:,0],Jrad[:]])
    """
    Start solving the HD equations with radiative transfer
    iterate this such that Tpost change a bit
    """
    for iiter in xrange(niter):
        print
        print ' Iteration: %d'%(iiter+1)
        print 'Tpost: %8.1f'%(Tpost)
        #
        # Reset the gas and dust conditions
        #
        gas     = None
        dust    = None
        gas = gasSpecs(rho=rhogas, nspecs=nspecs)
        print
        print 'Using gas Species: %d species'%(nspecs)
        """
            Generate the dust number densities
            """
        dust = dustSpecs(dustfrac=dustfrac, nspecs=ndust,
            mdust=mdust, gas=gas, size=adust)
        print 'Using dust Speces: %d sizes'%(ndust)
        print 'Initial size: %1.2f micron'%(dust.size[0]*1e4)
        print
        #
        # Check whether there are Nones in the input
        #
        print
        print '####################################################'
        print '  Input parameters       '
        print '  Shock velocity:   %2.2f  km/s'%(v0/1e5)
        print '  Pre shock T   :   %d     K   '%(t0)
        print '  Pre-shock n   : 1E%2.1f  cm-3'%(np.log10(gas._sumRho()))
        print '  Dust densities: %2.3e   cm-3'%(dust._sumRho())
        print '####################################################'
        print
        #
        # Solve this
        #
        wpre, vshock = solveHDrt(x=xpre, gas=gas, dust=dust, v0=v0, t0=t0,
                 Jrad=Jrad, abserr=1e-6, telerr=1e-8)
        sol0 = wpre
        print
        print '####################################################'
        print '  Condition at the end      '
        print '  velocity:   %2.2f  km/s'%(sol0[-2][1]/1e5)
        print '  Pre shock T   :   %d     K   '%(sol0[-2][2])
        if ndust is not None:
            print '  Dust Temperature: %d     K   '%(sol0[-2][9])
            print '  Dust Densities: %1.3e   cm^_3'%(sol0[-1][7])
            print '  Dust size     :   %1.2f microns'%(sol0[-2][10])
        print '####################################################'
        print
        """
            update the radiative transfer
            - Calculate the Jrad
            - Caluclate tau
        """
        Tpre=sol0[0, -2]
        oldtpost = Tpost
        #
        # Check the Frad and calculate Tpost accordingly
        #
        if np.abs(Frad[0]) > np.abs(Frad[-1]):
            corrFrad = Frad[0]
        else:
            corrFrad = Frad[-1]
        changeTpost = np.abs(np.power(np.abs(corrFrad/ss), 0.25))/2.
        if corrFrad > 0.0:
            Tpost += changeTpost
        else:
            Tpost -= changeTpost
        """"""
        delT = np.abs(oldtpost - Tpost)/Tpost
        print 'Frad: %2.3e -- %2.3e'%(Frad[0], Frad[-1])
        print 'Iter: %d -- Tpost: %d %2.3e %2.3e'%(iiter, Tpost, changeTpost,
            delT)
        print
        """"""
        print
        print 'Solving radiation field...'
        print
        from joblib import Parallel, delayed
        if nspecs is None:
            taumax, tau, dtau, srcall = calc_tau(sol0[:,0], numden, mass,
                             sol0[:, 2])
            Jrad = [Parallel(n_jobs=ncpu)(delayed(calcJrad)(Tpre, Tpost,
                tau[ix], tau, dtau,taumax,temps=sol0[:,2]) for ix in
                xrange(dtau.shape[0]))]
        else:
            if ndust is None: # no dust
                taumax, tau, dtau, srcall = calc_tauall(sol=sol0, gas=gas,
                    dust=dust)
                Jrad = calcJrad(Tpre=Tpre, Tpost=Tpost, srcall=srcall, tau=tau)
            else: # with dust
                tau, srcall = calc_tauall(sol=sol0, gas=gas,
                    dust=dust)
                #
                # Vectorize method as of Jul 2015
                #
                Jrad, Frad = calcJrad(Tpre=Tpre, Tpost=Tpost, srcall=srcall,
                    tau=tau)
            """"""
        Jrad = np.array([sol0[:,0],Jrad[:]])
        """"""
        if delT < 1e-4:
            iiter = niter+1
            break
        if Tpost > 1e5: break
    """"""
    print
    print '#### DONE ###'
    print
    return sol0, [tau, Jrad], vshock, Frad
""""""
    
