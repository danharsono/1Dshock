import numpy as np
from python.natconst import *
from solvers import solveHD
from oneDRT import calc_tauall, calcJrad
from gasspec import gasSpecs
from dustspec import dustSpecs
from progressbar import ProgressBar, Percentage, Bar

"""
The main part of the shock code
"""
def shock_main(numden=1e14, rhogas=1e-9, nspecs=None, ndust=None, adust=300e-4, v0=6e5, t0=300., sizex=5.0, numpoints=1e5, dustfrac=0.005, mdust=3.3, ncpu=3, niter=5, restart = False, Tpost0=1300.0):
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
    print '  Pre-shock n   : 1E%2.1f  cm-3'%(np.log10(sum(gas.numden)))
    print '  Dust densities: 1E%2.1f   cm-3'%(np.log10(dust.numden[0]))
    print '####################################################'
    print 
    #
    # Solve the whole shock
    # Pre-shock part
    #
    xpre = np.logspace(-3, np.log10(sizex), numpoints)
    xpre = xpre[::-1]
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
        wpre, vshock, gasKap= solveHD(x=xpre, gas=gas, dust=dust, v0=v0, t0=t0,
            telerr=1e-5)
    sol0 = wpre
    print
    print '####################################################'
    print '  Condition at the end      '
    print '  velocity:   %2.2f  km/s'%(sol0[-2][1]/1e5)
    print '  Gas temperature   :   %d     K   '%(sol0[-2][2])
    print '  Postshock n   : 1E%2.1f  cm-3'%(np.log10(sum(sol0[-2][
        3:gas.nspecs+3])))
    if ndust is not None:
        dumid = 3+gas.nspecs
        print '  Dust Temperature: %d     K   '%(sol0[-2][dumid+2])
        print '  Dust Densities: 1E%2.1f   cm^-3'%(np.log10(sol0[-2][dumid+0]))
        print '  Dust size     :   %1.4f microns'%(sol0[-2][dumid+3])
    print '####################################################'
    print
    #
    Tpre    = sol0[0,-2]
    Tpost   = Tpost0
    """
    Initialize the radiative transfer
      - Calculate the Jrad
      - Caluclate tau
    """
    #
    # Get the opacities
    #
    gas._getOpacs()
    print
    print 'Solving radiation field...'
    print
    tau, srcall = calc_tauall(sol=sol0, gas=gas, dust=dust, gasKap=gasKap)
    print 'TAU is done...'
    #
    # Vectorize method as of Jul 2015
    #
    Jrad, Frad = calcJrad(Tpre=Tpre, Tpost=Tpost, srcall=srcall,
        tau=tau)
    #
    # Assume a radiative mean intensidies
    #
    Jrad = np.zeros(sol0[:,0].shape[0])
    Jrad[:numpoints] = ss*np.power(Tpre, 4.)/np.pi
    Jrad[numpoints:] = ss*np.power(Tpost,4.)/np.pi
    Jrad = np.array([sol0[:,0],Jrad[:]])
    corrFrad = 0.0
    """
    Start solving the HD equations with radiative transfer
    iterate this such that Tpost change a bit
    """
    for iiter in xrange(niter):
        print
        print ' Iteration: %d, Tpost: %8.2f, Frad: %2.5e'%(iiter+1,
            Tpost, corrFrad)
        #
        # Old temperature solution
        #
        oldT = sol0[:,2]
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
        # Get opacities
        #
        gas._getOpacs()
        #
        # Check whether there are Nones in the input
        #
        print
        print '####################################################'
        print '  Input parameters       '
        print '  Shock velocity:   %2.2f  km/s'%(v0/1e5)
        print '  Pre shock T   :   %d     K   '%(t0)
        print '  Pre-shock n   : 1E%2.1f  cm-3'%(np.log10(sum(gas.numden)))
        print '  Dust densities: 1E%2.1f   cm-3'%(np.log10(dust.numden[0]))
        print '####################################################'
        print
        #
        # Solve this
        #
        wpre, vshock, gasKap = solveHD(x=xpre, gas=gas, dust=dust, v0=v0, t0=t0,
            haveJ=True, Jrad=Jrad, abserr=1e-8, telerr=1e-8)
        sol0 = wpre
        print
        print '####################################################'
        print '  Condition at the end      '
        print '  velocity:   %2.2f  km/s'%(sol0[-2][1]/1e5)
        print '  Gas temperature   :   %d     K   '%(sol0[-2][2])
        print '  Postshock n   : 1E%2.1f  cm-3'%(np.log10(sum(sol0[-2][
            3:gas.nspecs+3])))
        if ndust is not None:
            dumid = 3+gas.nspecs
            print '  Dust Temperature: %d     K   '%(sol0[-2][dumid+2])
            print '  Dust Densities: 1E%2.1f   cm^-3'%(np.log10(sol0[-2][
                dumid+0]))
            print '  Dust size     :   %1.4f microns'%(sol0[-2][dumid+3])
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
        if Frad[0] < 0.0:
            corrFrad = Frad[0]
        else:
            corrFrad = Frad[-5]
        changeTpost = np.minimum(np.abs(np.power(np.abs(corrFrad/ss), 0.25))
            /100., 5.0)
        if corrFrad > 0.0:
            Tpost += changeTpost
        else:
            Tpost -= changeTpost
        """"""
        #
        # Check the temperature at boundary
        #
        Tpost1   = np.minimum(sol0[-1,2], sol0[-1,-2])
        Tpost = (Tpost + Tpost1)/2.
        delT = np.abs(oldtpost - Tpost)
        print
        print 'Solving radiation field...'
        print
        tau, srcall = calc_tauall(sol=sol0, gas=gas, dust=dust, gasKap=gasKap)
        #
        # Vectorize method as of Jul 2015
        #
        Jrad, Frad = calcJrad(Tpre=Tpre, Tpost=Tpost, srcall=srcall,
            tau=tau)
        """"""
        if iiter < 2:
            """
            #
            # Assume a radiative mean intensidies
            #
            """
            Tpost = Tpost0
            Jrad = np.zeros(sol0[:,0].shape[0])
            Jrad[:numpoints] = ss*np.power(Tpre, 4.)/np.pi
            Jrad[numpoints:] = ss*np.power(Tpost,4.)/np.pi
        """"""
        Jrad = np.array([sol0[:,0],Jrad[:]])
        #
        # Calculate the change in temperature
        #
        Tchange = np.sqrt(np.power(np.abs(oldT - sol0[:,2]),2.0).sum())
        delT = np.maximum(delT, Tchange/np.float(xpre.shape[0]))
        print 'Frad: %2.5e -- %2.5e'%(Frad[0], Frad[-3:].mean())
        print 'Iteration: %d -- delT: %2.4f'%(iiter+1, delT)
        """"""
        #
        # Check for convergence
        #
        if (delT < 1e-2):
            iiter = iiter+1
            break
        """"""
    """"""
    print
    print '#### DONE ###'
    print
    return sol0, [tau, Jrad], vshock, Frad
""""""
    
