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

def shock_main(numden=1e14, rhogas=1e-9, nspecs=None, ndust=None, adust=300e-4, mass=2.0*mp, v0=6e5, t0=300., sizex=5.0, numpoints=1e5, mugas=2.8, dustfrac=0.005, mdust=3.3, ncpu=3, niter=1):
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
    print '  Dust densities: %2.3e   cm-3'%(dust._sumRho())
    print '####################################################'
    print 
    #
    # Solve the whole shock
    # Pre-shock part
    #
    xpre = np.logspace(-10, np.log10(sizex), numpoints)
    xpre = -xpre[::-1]
    xpre = xpre[:-1]
    #
    # Add the shock front
    #
    xpre = np.concatenate((xpre, np.array([0.0])))
    #
    # Add the post shock
    #
    xpost = -xpre[20:-1]
    xpre = np.concatenate((xpre, xpost[::-1]))
    #
    # Solve this
    #
    wpre, vshock = solveHD(x=xpre, gas=gas, dust=dust, v0=v0, t0=t0)
    sol0 = wpre
    print
    print '####################################################'
    print '  Condition at the end      '
    print '  velocity:   %2.2f  km/s'%(sol0[-2][1]/1e5)
    print '  Pre shock T   :   %d     K   '%(sol0[-2][2])
    if ndust is not None:
        print '  Dust Temperature: %d     K   '%(sol0[-2][9])
    print '####################################################'
    print
    print wpre[-1,:]
    #
    Tpre=sol0[0,2]
    Tpost=1500.0
    #"""
    #Initialize the radiative transfer
      #- Calculate the Jrad
      #- Caluclate tau
    #"""
    print 
    print 'Solving radiation field...'
    print
    from joblib import Parallel, delayed
    if nspecs is None:
        taumax, tau, dtau = calc_tau(sol0[:,0], numden, mass, 
            sol0[:, 2])    
        Jrad = [Parallel(n_jobs=ncpu)(delayed(calcJrad)(Tpre, Tpost, 
            tau[ix], tau, dtau,taumax,temps=sol0[:,2]) for ix in 
            xrange(dtau.shape[0]))]
    else:
        if ndust is None: # no dust
            taumax, tau, dtau = calc_tau(sol0[:,0], sum(gas.numden), 
                gas.mugas*mp, sol0[:, 2])
            Jrad = [Parallel(n_jobs=ncpu)(delayed(calcJrad)(Tpre, 
                Tpost, tau[ix], tau, dtau, taumax, temps=sol0[:,2])
                for ix in xrange(dtau.shape[0]))]
        else: # with dust
            taumax, tau, dtau = calc_tauall(sol=sol0, gas=gas,
                dust=dust)
            Jrad = [Parallel(n_jobs=ncpu)(delayed(calcJrad)(Tpre, 
                Tpost, tau[ix], tau, dtau, taumax, sol=sol0, gas=gas, 
                ix=ix) for ix in xrange(dtau.shape[0]))]
    Jrad = np.array(Jrad)
    Jrad = np.array([sol0[:,0],Jrad[0,:]])
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
        gas = gasSpecs(rho=rhogas, nspecs=nspecs)
        if ndust is not None:
            """
            Generate the dust number densities
            """
            dust = dustSpecs(dustfrac=dustfrac, nspecs=ndust, 
                mdust=mdust, gas=gas, size=adust)
        """"""
        # Solve the pre shock conditions
        if ndust is None:
            wpre = solveHDrt(x=xpre, gas=gas, v0=v0, t0=t0, Jrad=Jrad)
        else:
            wpre = solveHDrt(x=xpre, gas=gas, dust=dust, 
                v0=v0, t0=t0, Jrad=Jrad)
        """"""
        sol1 = []
        for ix in xrange(len(xpre)):
            vel = wpre[ix][0]
            if gas.nspecs > 1:
                if ndust is not None:
                    sol1.append([ xpre[ix],wpre[ix][0], 
                        wpre[ix][1], wpre[ix][2]/vel, 
                        wpre[ix][3]/vel, wpre[ix][4]/vel,
                        wpre[ix][5], wpre[ix][6], wpre[ix][7],
                        wpre[ix][8]
                        ])
                else:
                    sol1.append([ xpre[ix],wpre[ix][0], 
                        wpre[ix][1], wpre[ix][2]/vel, 
                        wpre[ix][3]/vel, wpre[ix][4]/vel])
            else:
                sol1.append([ xpre[ix],wpre[ix][0], wpre[ix][1], 
                    wpre[ix][2]])
        """"""        
        print
        print '####################################################'
        print '  Condition before shock      '
        print '  velocity:   %2.2f  km/s'%(sol1[-2][1]/1e5)
        print '  Pre shock T   :   %d     K   '%(sol1[-2][2])
        print '  Pre-shock n   : 1E%2.1f  cm-3'%(np.log10(gas._sumRho()))
        if ndust is not None:
            print '  Dust Temperature: %d     K   '%(sol1[-2][8])
            print '  Dust densities: 1E%2.1f cm-3  '%(
                np.log10(dust._sumRho()))
        print '####################################################'
        print
        #
        # If the velocity is 0
        # or the temperature  is 0
        # exit
        #
        if (sol1[-2][1] < 0.5e5) or (sol1[-2][2] < 5.):
            print ' ##!!!!!!!!!!!!!!!!!!!!!@@'
            print 'ERROR HERE: Conditions before the shock are not satisfied'
            print
            raise SystemExit
        """"""
        #
        #Shock front
        #
        if sol1[-1][1] < 4.0e5: sol1[-1][1] = 4.0e5
        if nspecs is None:
            P1 = (kk/(mugas*mp) * numden * 2.0*mp * sol1[-1][2])
            P2, rho2, u2 = get_RHcondition(rho1 = mass*numden, 
                u1 = sol1[-1][1], P1=P1)
            t2 = P2*(mugas*mp)/(kk*rho2)
        else:
            P1 = (kk/(gas.mugas*mp) * gas._sumRho() * 
                sol1[-1][2])
            P2, rho2, u2 = get_RHcondition(rho1 = gas._sumRho(), 
                u1 = sol1[-1][1], P1=P1)
            t2 = P2*(gas.mugas*mp)/(kk*rho2)
            gas._updateRho(rho=rho2)
        """"""
        #raise SystemExit
        print 
        print 'Solving post shock....'
        print
        #
        # Solve the post shock
        #
        if nspecs is None:
            wpost = solveHDrt(x=xpost, numden=rho2/(2.0*mp), mass=mass, 
                v0=u2, t0=t2)
        else:
            if ndust is None:
                wpost = solveHD(x=xpost, gas=gas, v0=u2, t0=t2)
            else:
                wpost = solveHD(x=xpost, gas=gas, dust=dust, 
                    v0=u2, t0=t2, vdust=wpre[-1][6], tdust=wpre[-1][7])
            """"""
        """"""
        #
        # Save the values
        #
        """"""
        for ix in xrange(len(xpost)):
            vel = wpost[ix][0]
            if gas.nspecs > 1:
                if ndust is not None:
                    sol1.append([ xpost[ix],wpost[ix][0], 
                        wpost[ix][1], wpost[ix][2]/vel, 
                        wpost[ix][3]/vel, wpost[ix][4]/vel,
                        wpost[ix][5], wpost[ix][6], wpost[ix][7],
                        wpost[ix][8]
                        ])
                else:
                    sol1.append([ xpost[ix],wpost[ix][0], 
                        wpost[ix][1], wpost[ix][2]/vel, 
                        wpost[ix][3]/vel, wpost[ix][4]/vel])
                
            else:
                sol1.append([xpost[ix],wpost[ix][0], wpost[ix][1], 
                    wpost[ix][2]])
        """"""
        sol1 = np.array(sol1)
        """
            
        update the radiative transfer
        - Calculate the Jrad
        - Caluclate tau
        """
        sol0 = np.array(sol1)
        # Update Tpost
        # Use bisection method
        oldtpost = Tpost
        Tpost = 0.5*(Tpost + sol0[-1,2])
        delT = np.abs(oldtpost - Tpost)/Tpost
        """"""
        print
        print '####################################################'
        print '  Condition at the end      '
        print '  velocity:   %2.2f  km/s'%(sol0[-2][1]/1e5)
        print '  Pre shock T   :   %d     K   '%(sol0[-2][2])
        if ndust is not None:
            print '  Dust Temperature: %d     K   '%(sol0[-2][8])
        print '  Temperature change:   %d  per   '%(delT*1e2)
        print '####################################################'
        print
        print 'Solving radiation field...'
        print
        if nspecs is None:
            taumax, tau, dtau = calc_tau(sol0[:,0], numden, mass, 
                sol0[:, 2])    
            Jrad = [Parallel(n_jobs=ncpu)(delayed(calcJrad)(Tpre, Tpost, 
                tau[ix], tau, dtau,taumax,temps=sol0[:,2]) for ix in 
                xrange(dtau.shape[0]))]
        else:
            if ndust is None: # no dust
                taumax, tau, dtau = calc_tau(sol0[:,0], sum(gas.numden), 
                    gas.mugas*mp, sol0[:, 2])
                Jrad = [Parallel(n_jobs=ncpu)(delayed(calcJrad)(Tpre, 
                    Tpost, tau[ix], tau, dtau, taumax, temps=sol0[:,2])
                    for ix in xrange(dtau.shape[0]))]
            else: # with dust
                taumax, tau, dtau = calc_tauall(sol=sol0, gas=gas,
                    dust=dust)
                Jrad = [Parallel(n_jobs=ncpu)(delayed(calcJrad)(Tpre, 
                    Tpost, tau[ix], tau, dtau, taumax, sol=sol0, gas=gas, 
                    ix=ix) for ix in xrange(dtau.shape[0]))]            
        Jrad = np.array(Jrad)
        Jrad = np.array([sol0[:,0],Jrad[0,:]])
        """"""
        if delT < 1e-4: iiter = niter+1
        if Tpost > 1e5: break
    """"""
    print
    print '#### DONE ###'
    print
    #raise SystemExit
    return sol0, [tau, Jrad], vshock
""""""
    
