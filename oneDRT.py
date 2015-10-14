import numpy as np
from python.natconst import *
from scipy.special import expn
from python.my_header import *
#import multiprocessing as mp
from joblib import Parallel, delayed
ncpu = 1

"""
The 1D radiative transfer module: Plane parallel
"""
def getEps(adust):
    """
    Return the absorption efficiencies
    """
    try:
        return 0.8*np.minimum(np.ones(adust.shape), adust/(2.0e-4))
    except StandardError:
        return 0.8*np.minimum(1.0, adust/(2.0e-4))
""""""
def calc_tauall(sol=None, gas=None, dust=None, gasKap=None):
    """
    Calculate the optical depth
    dtau = alpha ds = rho kapp ds
    """
    nspec   = gas.nspecs
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
    #    kapp = np.array([getKap(a) for a in sol[:,2]]) # gas
    kapp = gasKap
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
    return tau, srcall
""""""
def getJrad(ix, Ipre, Ipost, taumax, src, tau):
    """
    Return radiation field at each grid
    """
    Jrad = 0.0
    Frad = 0.0
    Jrad = src[ix]
    Frad = 2. * src[ix]
    #
    # Add the other terms
    #
    Jrad += 0.5 * (Ipre - src[0])*expn(2.0, taumax-tau[ix])
    Frad = 2.*expn(3., taumax-tau[ix]) * (Ipre - src[0])
    Jrad += 0.5 * (Ipost - src[-1])*expn(2.0, tau[ix])
    Frad -= 2.*expn(3., tau[ix]) * (Ipost - src[-1])
    #
    # Use the derivative of the source functions
    #
    if ix == 0:
        Jrad += -0.5 * ((src[:-1]-src[1:])*expn(
            2., tau[ix]-tau[:-1])).sum()
        Frad += 2.0 * ((src[:-1]-src[1:])*expn(
            3., tau[ix]-tau[:-1])).sum()
    elif ix == tau.shape[0]-1:
        Jrad += 0.5 * ((src[:-1]-src[1:])*expn(
            2., np.abs(tau[ix]-tau[:-1]))).sum()
        Frad += 2. * ((src[:-1]-src[1:])*expn(
            3., np.abs(tau[ix]-tau[:-1]))).sum()
    else:
        """
            Need to split the positive and negative
        """
        dumleft     = 0.5 * ((src[:ix]-src[1:ix+1]) *expn(
            2., np.abs(tau[ix]-tau[:ix]))).sum()
        dumright    = -0.5 * ((src[ix:-1]-src[ix+1:])*expn(
            2., tau[ix]-tau[ix:-1])).sum()
        Jrad += dumleft + dumright
        #
        # Frad
        #
        dumleft     = 2. * ((src[:ix]-src[1:ix+1]) *expn(
            3., np.abs(tau[ix]-tau[:ix]))).sum()
        dumright    = 2. * ((src[ix:-1]-src[ix+1:])*expn(
            3., tau[ix]-tau[ix:-1])).sum()
        Frad += dumleft + dumright
    """"""
    return Jrad,Frad
""""""

def calcJrad(Tpre=None, Tpost=None, srcall=None, tau=None, ncpu=3):
    """ 
    Calculate the mean itensity  and radiative flux
    """
    Ipre = (ss/np.pi)*np.power(Tpre, 4.)
    Ipost = (ss/np.pi)*np.power(Tpost, 4.)
    taumax = tau.max()
    #
    # Empty arrays
    #
    Jrad = np.empty(tau.shape[0])
    Frad = np.empty(tau.shape[0])
    #
    # Calculate Jrad
    #
    results = np.array(Parallel(n_jobs=ncpu, backend='multiprocessing')(
        delayed(getJrad)(ix,Ipre, Ipost, taumax, srcall, tau)
        for ix in xrange(tau.shape[0])))
    Jrad[:]    = results[:,0]
    Frad[:]    = results[:,1]
    return Jrad, Frad
""""""
