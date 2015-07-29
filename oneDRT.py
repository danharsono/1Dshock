import numpy as np
from python.natconst import *
from scipy.special import expn
from python.my_header import *
#import multiprocessing

"""
The 1D radiative transfer module: Plane parallel
"""

"""
Get the planck mean opacities
"""
def getKap(temp=None):
    if (temp < 1000):
        return 5.0
    elif (temp < 0):
        print 'Negative Temperature!'
        return 1.0
    elif ( (temp >1000.) and (temp < 3000.)):
        return 0.5
    else:
        return 0.05
""""""

def calc_tau(x, numden, mass, temps):
    """
    Calculate the optical depth
    dtau = alpha ds = rho kapp ds
    """
    x = np.array(x)
    dtau = np.zeros(x.shape[0])
    tau = np.zeros(x.shape[0])
    dx = np.zeros(x.shape[0]+1)
    dx[0] = 0.5*(x[1]-x[0])
    dx[1:-1] = x[1:]-x[:-1]
    
    # calculate this
    kapp = np.array([getKap(a) for a in temps])
    dtau = numden*mass*kapp*dx[1:]
    tau[-1] = 0.0
    for ix in xrange(dtau.shape[0]-2,-1,-1):
        tau[ix] = tau[ix+1]+ dtau[ix+1]
    taumax = tau.max()
    
    return taumax,tau, dtau
""""""

def getEps(adust):
    """
    Return the absorption efficiencies
    """
    try:
        return 0.8*np.minimum(np.ones(adust.shape), adust/(2.0e-4))
    except StandardError:
        return 0.8*np.minimum(1.0, adust/(2.0e-4))
""""""

def calc_tauall(sol=None, gas=None, dust=None):
    """
    Calculate the optical depth
    dtau = alpha ds = rho kapp ds
    """
    nspec = gas.nspecs
    idust   = 2+nspec # where dust information starts
    #
    # Dust: ndust, vdust, tdust, adust
    #
    x = np.array(sol[:,0])
    dtau = np.zeros(x.shape[0])
    tau = np.zeros(x.shape[0])
    dx = np.zeros(x.shape[0]+1)
    dx[0] = 0.5*(x[1]-x[0])
    dx[1:-1] = x[1:]-x[:-1]
    #
    # calculate this
    # number densities are in terms of n * v
    #
    kapp = np.array([getKap(a) for a in sol[:,2]]) # gas
    kapd = sol[:,idust+1]*np.pi*(sol[:,idust+4]*sol[:,idust+4])*getEps(
        sol[:,idust+4])
    gasrho = np.sum(sol[:,3:3+nspec]*gas.mass, axis=1)/sol[:,1]
    dtau = (gasrho*kapp + kapd)*dx[1:]
    tau[-1] = 0.0
    for ix in xrange(dtau.shape[0]-2,-1,-1):
        tau[ix] = tau[ix+1]+ dtau[ix+1]
    taumax = tau.max()
    #
    # calculate the source function
    #
    srcall = np.zeros(x.shape[0])
    fmt = '%d %2.2e %2.2e %2.2e %2.2e %d %d %2.2e %2.3e'
    for ix in xrange(x.shape[0]):
        bottom = gasrho[ix]*kapp[ix] + kapd[ix]
        top1    = gasrho[ix]*kapp[ix]*(ss/np.pi)*np.power(sol[ix,2],4.)
        top2    = kapd[ix] * (ss/np.pi)*np.power(sol[ix,idust+3], 4.)
        srcall[ix] = (top1+top2)/bottom
#        print fmt%(ix, gasrho[ix], sol[ix,idust+1], bottom, top1+top2,
#                   sol[ix,2], sol[ix, idust+3], srcall[ix],
#                   np.power(sol[ix,2],4.))
#        if ix == 10: raise SystemExit
#    raise SystemExit
#
#    srcall = (gasrho*kapp*(ss/np.pi)*np.power(sol[:,2],4.) + kapd *
#              np.power(sol[:,idust+3], 4.) )/(gasrho*kapp + kapd)
#    semilogy(x,srcall,'ko')
#    show()
#    close()
#    raise SystemExit
    return tau, srcall
""""""

def addJrad(tau, arrtau, temps, dtau, Jrad,ind):
    delttau = np.abs(arrtau[ind]-tau)
    if delttau < 1e-10:
        Jrad+= 10.0
    else:
        Jrad += 0.5*( (ss/np.pi)*temps[ind]**(4.0)*
            expn(1,delttau)*dtau[ind])
    return Jrad
""""""

def calcJrad(Tpre=None, Tpost=None, srcall=None, tau=None):
    """ 
    Calculate the mean itensity  and radiative flux
    """
    Ipre = (ss/np.pi)*np.power(Tpre, 4.)
    Ipost = (ss/np.pi)*np.power(Tpost, 4.)
    taumax = tau.max()
    #
    # Calculate Jrad
    #
    Jrad = np.zeros(tau.shape[0])
    Frad = np.zeros(tau.shape[0])
    for ix in xrange(tau.shape[0]):
        Jrad[ix] = srcall[ix]
        #
        # Add the other terms
        #
        Jrad[ix] += 0.5 * (Ipre - srcall[0])*expn(2.0, taumax-tau[ix])
        Frad[ix] = 2.*np.pi*expn(3., taumax-tau[ix]) * (Ipre - srcall[0])
        Jrad[ix] += 0.5 * (Ipost - srcall[-1])*expn(2.0, tau[ix])
        Frad[ix] -= 2.*np.pi*expn(3., tau[ix]) * (Ipost - srcall[-1])
        #
        # Add the derivatives and contribution from all the other cells
        #
        if ix == 0:
            Jrad[ix] += -0.5 * ((srcall[:-1]-srcall[1:])*expn(
                2., tau[ix]-tau[:-1])).sum()
        elif ix == tau.shape[0]-1:
            Jrad[ix] += 0.5 * ((srcall[:-1]-srcall[1:])*expn(
                2., np.abs(tau[ix]-tau[:-1]))).sum()
        else:
            """
                Need to split the positive and negative
            """
            dumleft     = 0.5 * ((srcall[:ix]-srcall[1:ix+1]) *expn(
                2., np.abs(tau[ix]-tau[:ix]))).sum()
            dumright    = -0.5 * ((srcall[ix:-1]-srcall[ix+1:])*expn(
                2., tau[ix]-tau[ix:-1])).sum()
            Jrad[ix] += dumleft + dumright
        """"""
        #
        # Add all the radiative fluxes
        #
        Frad[ix] += 2.*np.pi * ( (srcall[:-1]-srcall[1:]) *
            expn(3., np.abs(tau[ix] - tau[:-1])) ).sum()
    """"""
    return Jrad, Frad
""""""
