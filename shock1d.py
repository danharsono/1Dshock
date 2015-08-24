import numpy as np
from python.natconst import *
import time
"""
Decipher the input parameters from the gas and dust parameters
"""
def decipherW(w, gas, dust):
    x1, x2, n1, n2, n3, n4, nd, vd, td, ad = w
    if (ad < 0.0): ad = 1e-10
    w = [x1,x2,n1/x1,n2/x1,n3/x1,n4/x1,nd,vd,td,ad]
    return w
""""""
"""
The vector field of the matrix to solve
"""
def vectorfield(x,w, p):
    """
    w   :  the variables to solve -> x, y 
    x   :  position with respect to the shock
    p   :  parameters
    """
    gas, dust, Jrad, haveJ, debug = p
    if haveJ:
        from oneDRT import getEps
        #
        # Interpolate to find the current Jrad
        #
        curJ = np.interp(x, Jrad[0], Jrad[1])
    else:
        curJ = None
    """"""
    w = decipherW(w, gas, dust)
    #
    # Calculate the constants first
    #
    Rhogas  = sum([a*b for (a,b) in zip(w[2:gas.nspecs+2], gas.mass)])
    Ngas    = sum(w[2:gas.nspecs+2])
    #
    # Dust changes
    #
    mdust       = np.zeros(dust.nspecs)
    fdrag       = np.zeros(dust.nspecs)
    dxa         = np.zeros(dust.nspecs)
    dxtd        = np.zeros(dust.nspecs)
    Cpd         = np.zeros(dust.nspecs)
    dustLoss    = 0.0
    for idust in xrange(dust.nspecs):
        nd = w[2+gas.nspecs + idust*4 + 0]
        vd = w[2+gas.nspecs + idust*4 + 1]
        td = w[2+gas.nspecs + idust*4 + 2]
        ad = w[2+gas.nspecs + idust*4 + 3]
        mdust[idust]    = (4.0/3.0)*np.pi*dust.mdust*np.power(ad, 3.)
        fdrag[idust]    = dust._calculateFdrag(w=w, gas=gas, idust=idust,
            rhoGas=Rhogas, NGas = Ngas)
        dxa[idust]      = dust._calcDxa(w=w, idust=idust,
            rhoGas=Rhogas, NGas = Ngas, gas=gas, Jrad=curJ)
        dxtd[idust]     = dust._calcDxTd(w=w, idust=idust,
            rhoGas=Rhogas, NGas=Ngas, gas=gas, Jrad=curJ)
        Cpd[idust]      = dust._getCpd(td = w[-2])
        if dxa[idust] > 0.0:
            dustLoss    += 0.0
        else:
            dustLoss    += -4.0 *np.pi * (nd * vd * dust.mdust *
                ad*ad * dxa[idust])
        """"""
    """"""
    #
    # Calculate the variable here
    #
    varA = Rhogas*w[0] - (kk*w[1]/w[0])*Ngas
    varB = Ngas*kk
    #
    # Variable C now has dust properties
    # Calculate the gas part
    #
    RatesMat = (gas._calculateR(vars=w)).sum(axis=0)
    dumGas = sum([a*(w[0]*b + (kk*w[1]/w[0])) for (a,b) in zip(RatesMat,
        gas.mass)])
    #
    # Add the dust evaporation mass
    #
    dumGas += dustLoss/(gas.mass[3])*(w[0]*gas.mass[3] + (kk*w[1]/w[0]))
    #
    # Get the dust part
    #
    #    dumDust = nd * fdrag + 4.0*np.pi*w[-1]*w[-1]*dust.mdust*(
    #        w[gas.nspecs+2]*w[gas.nspecs+2])*dxa
    dumDust = [w[2+gas.nspecs+id*4 +0] * (fdrag[id] +
        4.0*np.pi * np.power(w[2+gas.nspecs+id*4 + 3],2) * dust.mdust *
        np.power(w[2+gas.nspecs+id*4+1], 2.) * dxa[id])
        for id in xrange(dust.nspecs)]
    varC = -(dumGas + sum(dumDust))
    #
    # Calculate the variables for energy equation
    #
    varD = w[0]*w[0]*Rhogas
    varE = kk*w[0]*sum([a*b for (a,b) in zip(gas.gamma, w[2:gas.nspecs+2])])
    #
    # Calculate the variable F with the dust
    #
    dumGas = sum([(a*(b*kk*w[1] + 0.5*w[0]*w[0]*c)) for (a,b,c,) in
        zip(RatesMat, gas.gamma, gas.mass)])
    dumGas += dustLoss/(gas.mass[3]) * (gas.gamma[3]*kk*w[1] +
        0.5*w[0]*w[0]*gas.mass[3])
    #
    # Dust part for variable F
    #
    #    dumDust = vd*nd * (fdrag + 2*np.pi*w[-1]*w[-1]*dust.mdust*vd*vd*dxa)
    #    dumDust += nd * mdust * vd * Cpd * dxtd
    #    dumDust += (4.0 * np.pi * w[-1]*w[-1] * nd * Cpd * w[-3] * w[-2] *
    #        dust.mdust * dxa)
    dumid = gas.nspecs+2
    dumDust = [ w[dumid+id*4+1] * w[dumid+id*4+0] *
        (fdrag[id] + 2*np.pi*np.power(w[dumid+id*4+3],2.)*dust.mdust*
        np.power(w[dumid+id*4+1], 2.)*dxa[id]) + w[dumid+id*4+0] * (
        mdust[id] * w[dumid+id*4+1] * Cpd[id] * dxtd[id]) + (4.0 * np.pi * (
        np.power(w[dumid+id*4+3],2.) * w[dumid+id*4+0] * Cpd[id] *
        w[dumid+id*4+1] * w[dumid+id*4+2]*dust.mdust*dxa[id]))
        for id in xrange(dust.nspecs)]
    varF = - (dumGas + sum(dumDust)) + gas._calculateFreeEnergy(w=w)
    radiation = 0.0
    if haveJ:
        radiation += (4.0*np.pi*Rhogas * gas._getKap(x2,
            destroyed=dust.destroyed) *
            (curJ - (ss*np.power(w[1], 4.)/np.pi)) +
            sum([4.0* np.pi * np.pi * w[dumid+4*id+0] *
                np.power(w[dumid+4*id+3], 2.) * dust._getEps(size =
                w[dumid+4*id+3]) * (curJ-(ss*np.power(w[dumid+4*id+2], 4.)
                /np.pi)) for id in xrange(dust.nspecs)]))
    varF += radiation
    #
    # The RHS matrix
    # Check few numbers to make sure these are not too large
    # or too small
    #
    f1 = (varC*varE - varB*varF)/(varA*varE - varB*varD)
    f2 = (varA*varF - varD*varC)/(varA*varE - varB*varD)
    #
    # Integration is over space dx
    # The change of n * v
    # changes in chemical species
    #
    fchem = [a for a in RatesMat]
    #
    # Add the dust evaporation
    #
    fchem = [a + b*dustLoss/c for (a,b,c) in zip(fchem, gas.dustfrac,
        gas.mass)]
    #
    #
    # DUST from here on
    #
    if dust.nspecs == 0:
        fdust = np.zeros(4)
    else:
        fdust = np.zeros(dust.nspecs*4)
        for idust in xrange(dust.nspecs):
            nd = w[2+gas.nspecs + idust*4 + 0]
            vd = w[2+gas.nspecs + idust*4 + 1]
            td = w[2+gas.nspecs + idust*4 + 2]
            ad = w[2+gas.nspecs + idust*4 + 3]
            fdust[idust*4+1] = fdrag[idust]/(mdust[idust]*vd)
            fdust[idust*4+0] = - (nd/vd) * fdust[idust*4+1]
            fdust[idust*4+2] = dxtd[idust]
            fdust[idust*4+3] = dxa[idust]
            if fdust[idust*4+0] < 1e-25: fdust[idust*4+0] = 0.0
        """"""
        #
        # Crash and NAN handles
        #
        if np.isnan(fdust).any():
            print 'NAN is found!'
            print 'Tgas and vgas: %e, %e'%(w[1], w[0]*1e-5)
            print 'ndust, vdust and Tdust: ', w[2+gas.nspecs:]
            print fdust
            print
            print fdrag, mdust
            print w
            raise SystemExit
    """"""
    """
        Limit the changes with respect to values to 1e-15
        f1 -> vgas
        f2 -> tgas
        f3 -- f6 -> number densities
        rest -> dust
    """
    f = [f1,f2]+fchem+[a for a in fdust]
    for (a,b) in zip(f,w):
        if np.abs(a/(b+1e-30)) < 1e-25: a = 0.0
    return f
""""""
