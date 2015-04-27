import numpy as np
from python.natconst import *
from solvers import get_RHcondition, solveHD, solveHDrt
from oneDRT import calc_tau, calcJrad, calc_tauall
from gasspec import gasSpecs
from dustspec import dustSpecs

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
    if nspecs is not None:
        """ 
        Calculate the number densities of the different species
        use the rhogas
        """
        gas = gasSpecs(rho=rhogas, nspecs=nspecs)
        print
        print 'Using gas Species: %d species'%(nspecs)    
        if ndust is not None:
            """
            Generate the dust number densities
            """
            dust = dustSpecs(dustfrac=dustfrac, nspecs=ndust, 
                mdust=mdust, gas=gas, size=adust)
            print 'Using dust Speces: %d sizes'%(ndust)
            print 'Initial size: %1.2f micron'%(dust.size[0]*1e4)
        print
    """"""
    # Check whether there are Nones in the input
    print 
    print '####################################################'
    print '  Input parameters       '
    print '  Shock velocity:   %2.2f  km/s'%(v0/1e5)
    print '  Pre shock T   :   %d     K   '%(t0)
    print '  Pre-shock n   : 1E%2.1f  cm-3'%(np.log10(gas._sumRho()))
    print '####################################################'
    print 
    #
    # Solve the pre shock conditions
    #
    # The pre shock grid
    xpre = np.logspace(np.log10(sizex), np.log10(2.0/3.0 * sizex), 
        np.int(numpoints*0.1))
    xpre = xpre[::-1]
    xpre = np.concatenate([np.logspace(-5, np.log10(2.0/3.0 * sizex),
        np.int(numpoints*0.9)+2), xpre[:-1]])
    xpre = -xpre[::-1]
    xpre = xpre[:-1]
    if nspecs is None:
        wpre = solveHD(x=xpre, numden=numden, mass=mass, v0=v0, t0=t0)
    else:
        if ndust is None:
            wpre = solveHD(x=xpre, gas=gas, v0=v0, t0=t0)
        else:
            wpre = solveHD(x=xpre, gas=gas, dust=dust, v0=v0, t0=t0)
    sol0 = []
    for ix in xrange(len(xpre)):
        vel = wpre[ix][0]
        if gas.nspecs > 1:
            if ndust is not None:
                sol0.append([ xpre[ix],wpre[ix][0], wpre[ix][1], 
                    wpre[ix][2]/vel, wpre[ix][3]/vel, wpre[ix][4]/vel,
                    wpre[ix][5], wpre[ix][6], wpre[ix][7], wpre[ix][8]
                    ])
            else:
                sol0.append([ xpre[ix],wpre[ix][0], wpre[ix][1], 
                    wpre[ix][2]/vel, wpre[ix][3]/vel, wpre[ix][4]/vel])
            
        else:
            sol0.append([ xpre[ix],wpre[ix][0], wpre[ix][1], 
                wpre[ix][2]])
    """"""
    print
    print '####################################################'
    print '  Condition before shock      '
    print '  velocity:   %2.2f  km/s'%(sol0[-2][1]/1e5)
    print '  Pre shock T   :   %d     K   '%(sol0[-2][2])
    if ndust is not None:
        print '  Dust Temperature: %d     K   '%(sol0[-2][8])
        print '  Dust density: %2.2e     '%(sol0[-2][7])
    print '####################################################'
    print 
    raise SystemExit
    #
    # Check few things here
    #
    #print xpre[-1], wpre[-1][5], wpre[-1][6], wpre[-1][7], dust.numden[0]
    #raise SystemExit
    #
    #Shock front
    #
    if nspecs is None:
        P1 = (kk/(mugas*mp) * numden * 2.0*mp * sol0[-2][2])
        P2, rho2, u2 = get_RHcondition(rho1 = mass*numden, 
            u1 = sol0[-1][1], P1=P1)
        t2 = P2*(mugas*mp)/(kk*rho2)
    else:
        P1 = (kk/(gas.mugas*mp) * gas._sumRho() * 
              sol0[-2][2])
        P2, rho2, u2 = get_RHcondition(rho1 = gas._sumRho(), 
            u1 = sol0[-1][1], P1=P1)
        t2 = P2*(gas.mugas*mp)/(kk*rho2)
        gas._updateRho(rho=rho2)
    """"""
    #
    # Update dust
    #
    dust._updateDust(allns=[sol0[-1][6]], size=sol0[-1][9])
    print u2, rho2, t2
    print sol0[-1][6], sol0[-1][9]
    if True:
        u2 = 2.1e5
        rho2 = 3.2e-9
        t2 = 2850.
        ndust = 2.5e-5
        dsize = 0.03
        gas._updateRho(rho=rho2)
        dust._updateDust(allns=[ndust], size=dsize)
    #
    # Solve the post shock
    #
    #print 
    #print 'Solving post shock....'
    #print
    #xpost = np.logspace(-5, np.log10(sizex), np.int(numpoints))
    #if nspecs is None:
        #wpost = solveHD(x=xpost, numden=rho2/(2.0*mp), mass=mass, 
            #v0=u2, t0=t2)
    #else:
        #if ndust is None:
            #wpost = solveHD(x=xpost, gas=gas, v0=u2, t0=t2)
        #else:
            #wpost = solveHD(x=xpost, gas=gas, dust=dust, 
                #v0=u2, t0=t2, tdust=wpre[-1][7], vdust=wpre[-1][6])
    #for ix in xrange(len(xpost)):
        #vel = wpost[ix][0]
        #if gas.nspecs > 1:
            #if ndust is not None:
                #sol0.append([ xpost[ix],wpost[ix][0], wpost[ix][1], 
                    #wpost[ix][2]/vel, wpost[ix][3]/vel, wpost[ix][4]/vel,
                    #wpost[ix][5], wpost[ix][6], wpost[ix][7], 
                    #wpost[ix][8]
                    #])
            #else:
                #sol0.append([ xpost[ix],wpost[ix][0], wpost[ix][1], 
                    #wpost[ix][2]/vel, wpost[ix][3]/vel, wpost[ix][4]/vel])
            
        #else:
            #sol0.append([xpost[ix],wpost[ix][0], wpost[ix][1], 
                #wpost[ix][2]])
    #""""""
    #
    #
    #
    sol0 = np.array(sol0)
    print
    print '####################################################'
    print '  Condition at the end      '
    print '  velocity:   %2.2f  km/s'%(sol0[-2][1]/1e5)
    print '  Pre shock T   :   %d     K   '%(sol0[-2][2])
    if ndust is not None:
        print '  Dust Temperature: %d     K   '%(sol0[-2][8])
    print '####################################################'
    print 
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
    return sol0, [tau, Jrad]
""""""
    
