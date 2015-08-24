import numpy as np
from python.natconst import *
import time

"""
Set global variables here
"""
isstart = 0
def isstart():
    isstart = 1
""""""

def set_oldx(a):
    global oldx
    oldx = a
""""""

tau = 1e-5
def set_tau(a=None, b=None):
    global tau
    tau = a + b
""""""

def getTau():
    return tau

"""
Define functions here
"""
def get_dtau(alpha, b):
    ds = b-oldx
    return alpha*ds
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
    x1, x2, n1, n2, n3, n4, nd, vd, td, ad = w
    if (ad < 0.0): ad = 1e-10
    w = [x1,x2,n1/x1,n2/x1,n3/x1,n4/x1,nd,vd,td,ad]
    #
    # Calculate the constants first
    #
    Rhogas  = sum([a*b for (a,b) in zip(w[2:gas.nspecs+2], gas.mass)])
    Ngas    = sum(w[2:gas.nspecs+2])
    #
    # Dust changes
    #
    mdust   = (4.0/3.0)*np.pi*dust.mdust*np.power(w[-1], 3.)
    fdrag   = dust._calculateFdrag(w=w, rhoGas=Rhogas, NGas = Ngas)
    dxa     = dust._calcDxa(w=w, rhoGas=Rhogas, NGas = Ngas, gas=gas, Jrad=curJ)
    dxtd    = dust._calcDxTd(w=w, rhoGas=Rhogas, NGas=Ngas, gas=gas, Jrad=curJ)
    Cpd     = dust._getCpd(td = w[-2])
    #
    # Calculate the variable here
    #
    varA = Rhogas*x1 - (kk*x2/x1)*Ngas
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
    dustLoss = -4.0 *np.pi * nd * vd * dust.mdust * ad*ad * dxa
    if dxa > 0.0:
        dustLoss = 0.0
    dumGas += dustLoss/(gas.mass[3])*(w[0]*gas.mass[3] + (kk*w[1]/w[0]))
    #
    # Get the dust part
    #
    dumDust = nd * fdrag + 4.0*np.pi*w[-1]*w[-1]*dust.mdust*(
        w[gas.nspecs+2]*w[gas.nspecs+2])*dxa
    varC = -(dumGas + dumDust)
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
    dumGas += dustLoss/(gas.mass[3]) * (gas.gamma[3]*kk*w[1] + 0.5*w[0]*w[0]*gas.mass[3])
    #
    # Dust part for variable F
    #
    dumDust = vd*nd * (fdrag + 2*np.pi*w[-1]*w[-1]*dust.mdust*vd*vd*dxa)
    dumDust += nd * mdust * vd * Cpd * dxtd
    dumDust += (4.0 * np.pi * w[-1]*w[-1] * nd * Cpd * w[-3] * w[-2] *
        dust.mdust * dxa)
    varF = - (dumGas + dumDust) + gas._calculateFreeEnergy(w=w)
    radiation = 0.0
    if haveJ:
        radiation += (4.0*np.pi*Rhogas * gas._getKap(x2,
            destroyed=dust.destroyed) *
            (curJ - (ss*np.power(w[1], 4.)/np.pi)) +
            4.0* np.pi * np.pi * nd  * ad * ad * dust._getEps(size = ad) *
            (curJ-(ss*np.power(w[-2],4.)/np.pi)) )
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
    #
    f3 = RatesMat[0] # (s cm^-3)^-1
    f4 = RatesMat[1]
    f5 = 0.0
    f6 = dustLoss/gas.mass[3]
    #
    #
    # DUST from here on
    #
    if dust.nspecs == 0:
        f7 = 0.0
        f8 = 0.0
        f9 = 0.0
        f10 = 0.0
    else:
        f8 = fdrag/(mdust*vd) # Dust velocities change
        f7 = -(nd/vd)*f8 # This is for number of dust
        if f7 < 1e-25: f7 = 0.0
        f9 = dxtd
        f10 = dxa
        #
        # Crash and NAN handles
        #
        if np.isnan(f7) or np.isnan(f8) or np.isnan(f9) or np.isnan(f10):
            print 'NAN is found!'
            print 'Densities: %2.5e  %2.5e'%(gas._sumRho(), dust._sumRho())
            print 'Masses: ', gas.mass, dust.mass
            print
            print 'Tgas and vgas: %e, %e'%(x1*1e-5, x2)
            print 'ndust, vdust and Tdust: %2.5e  %2.5e  %d'%(nd, vd, td)
            print '%e %e  %e  %e'%(f6, f7, f8, f9)
            print
            print fdrag, mdust, nd/vd
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
    if np.abs(f1/x1) < 1e-25: f1 = 0.0
    if np.abs(f2/x2) < 1e-25: f2 = 0.0
    if np.abs(f3/n1) < 1e-25: f3 = 0.0
    if np.abs(f4/n2) < 1e-25: f4 = 0.0
    if np.abs(f5/n3) < 1e-25: f5 = 0.0
    if (n4 > 1e-25 and (np.abs(f6/n4) < 1e-25)): f6 = 0.0
    if np.abs(f7/nd) < 1e-25: f7 = 0.0
    if np.abs(f8/vd) < 1e-25: f8 = 0.0
    if np.abs(f9/td) < 1e-25: f9 = 0.0
    if np.abs(f10/ad) < 1e-25: f10 = 0.0
    f = np.array([f1, f2, f3, f4, f5, f6, f7, f8, f9, f10])
    return f
""""""
