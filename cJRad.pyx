"""
Cython version of the radiation calculations
"""
from __future__ import division
from numpy cimport ndarray, dtype; from scipy.special import expn
import numpy as np; cimport numpy as np
cimport cython; from cpython cimport bool
from cython.parallel import prange
from joblib import Parallel, delayed
import threading
#
# Types
#
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t
#
# Constants here and externs from math
#
cdef extern from "math.h":
    double sqrt(double m)
    double pow(double m, double n)
    double exp(double a)
    double fabs(double a)
    double M_PI
    double fmin(double, double)
cdef double kk  = 1.3806488e-16
cdef double ss  = 5.670367e-5
#
# Functions HERE
#
"""
EPS
"""
cdef double getEPS(double ad):
    return 0.8*fmin(1.0, ad/(2.0e-4))
""""""
"""
The 1D radiative transfer module: Plane parallel
"""
@cython.boundscheck(False) # turns off bounds-checking
cdef np.ndarray [DTYPE_t, ndim=2] calctaus(int nspecs, ndarray[DTYPE_t, ndim=1] x, ndarray[DTYPE_t, ndim=1] kapp, ndarray[DTYPE_t, ndim=2] sol, ndarray[DTYPE_t, ndim=1] gmasses):
    """
    Calculate the optical depth
    dtau = alpha ds = rho kapp ds
    """
    #
    # Variables
    #
    cdef int idust = 2+nspecs # where dust information starts
    cdef double taumax, bottom, top1, top2
    cdef int ix
    #
    # Dust: ndust, vdust, tdust, adust
    #
    cdef ndarray[DTYPE_t, ndim=1] dtau  = np.zeros(x.shape[0], dtype=DTYPE)
    cdef ndarray[DTYPE_t, ndim=1] tau   = np.zeros(x.shape[0], dtype=DTYPE)
    cdef ndarray[DTYPE_t, ndim=1] dx    = np.zeros(x.shape[0]+1, dtype=DTYPE)
    dx[0] = 0.5*(x[1]-x[0]) # assign the first
    dx[1:-1] = x[1:]-x[:-1] # assing the last
    #
    # calculate this
    # number densities are in terms of n * v
    #
    cdef ndarray[DTYPE_t, ndim=1] kapd  = np.zeros(x.shape[0], dtype=DTYPE)
    for ix in xrange(x.shape[0]):
        kapd[ix]    = (sol[ix, <unsigned int> (idust+1)]*M_PI*
            (sol[ix,<unsigned int> (idust+4)]*
            sol[ix,<unsigned int> (idust+4)])*
            getEPS(sol[ix,<unsigned int> (idust+4)]))
    cdef ndarray[DTYPE_t, ndim=1] gasrho    = (np.sum(sol[:,3:3+nspecs]*
        gmasses, axis=1)/sol[:,1])
    for ix in range(x.shape[0]):
        dtau[ix]    = (gasrho[ix]*kapp[ix] + kapd[ix])*dx[<unsigned int> (1+ix)]
    tau[-1]     = 0.0
    for ix in range(dtau.shape[0]-2,-1,-1):
        tau[ix] = tau[ix+1]+ dtau[ix+1]
    taumax = tau.max()
    #
    # calculate the source function
    #
    cdef ndarray[DTYPE_t, ndim=1] srcall = np.zeros(x.shape[0], dtype=DTYPE)
    for ix in range(x.shape[0]):
        bottom  = gasrho[ix]*kapp[ix] + kapd[ix]
        top1    = gasrho[ix]*kapp[ix]*(ss/M_PI)*pow(sol[ix,2],4.)
        top2    = (kapd[ix] * (ss/M_PI)*pow(
            sol[ix,<unsigned int> (idust+3)], 4.) )
        srcall[ix]  = (top1+top2)/bottom
    #
    # arrays
    #
    cdef ndarray[DTYPE_t, ndim=2] stuffs = np.zeros((x.shape[0], 2), dtype=DTYPE)
    stuffs[:,0]     = tau
    stuffs[:,1]     = srcall
    return stuffs
""""""
#
# This is to call from outside of c
#
cpdef np.ndarray [DTYPE_t, ndim=2] calc_tauall(list p, ndarray[DTYPE_t, ndim=1] gasKap):
    """
    Call this from python as a buffer
    """
    #
    # variables
    #
    cdef ndarray[DTYPE_t, ndim=1] x     = np.array(p[0][:,0])
    return calctaus(<unsigned int> p[1].nspecs, x, gasKap, np.array(p[0]),
        np.array(p[1].mass))
""""""
#
# Radiaton field
#
#def getJrad(ix, ntau, Ipre, Ipost, taumax, src, tau):
cpdef getJrad(int ix, int ntau, double Ipre, double Ipost, double taumax, ndarray[DTYPE_t, ndim=1] src, ndarray[DTYPE_t, ndim=1] tau):
    """
    Return radiation field at each grid
    """
    #
    # Variables
    #
    cdef int ii
    cdef double Jrad, Frad, dumleft, dumright
    #
    #
    #
    Jrad = 0.0 # Reset
    Frad = 0.0
    Jrad += src[ix]
    Frad += 2. * src[ix]
    #
    # Add the other terms
    #
    Jrad += 0.5 * (Ipre - src[0])*expn(2., taumax-tau[ix])
    Frad += 2. * expn(3., taumax-tau[ix]) * (Ipre - src[0])
    Jrad += (0.5 * (Ipost - src[<unsigned int> (ntau-1)])*
        expn(2., tau[ix]) )
    Frad -= 2. * expn(3., tau[ix]) * (
        Ipost - src[<unsigned int> (ntau-1)])
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
    #if ix == 0:
    #    Jrad    += -0.5 * ((src[ii]-src[<unsigned int> (1+ii)]) *
    #        expn(2., tau[ix]-tau[ii]))
    #    Frad    += 2.0 * ((src[ii]-src[<unsigned int> (1+ii)]) *
    #        expn(3., tau[ix]-tau[ii]))
    #elif ix == ntau-1:
    #    for ii in range(ntau-1):
    #        Jrad += 0.5 * ((src[ii]-src[<unsigned int> (1+ii)]) *
    #            expn(2.,fabs(tau[ix]-tau[ii])))
    #        Frad += 2. * ((src[ii]-src[<unsigned int> (1+ii)]) *
    #            expn(3., fabs(tau[ix]-tau[ii])))
    #else:
    #    """
    #        Need to split the positive and negative
    #    """
    #    dumleft     = 0.0
    #    dumright    = 0.0
    #    #
    #    # Jrad
    #    #
    #    for ii in range(0,ix,1):
    #        dumleft     += 0.5 * ((src[ii]-src[<unsigned int> (1+ii)]) *
    #            expn(2., fabs(tau[ix]-tau[ii])))
    #    for ii in range(ix,ntau-1,1):
    #        dumright    += -0.5 * ((src[ii]-src[<unsigned int> (1+ii)]) *
    #            expn(2., tau[ix]-tau[ii]))
    #    Jrad += dumleft + dumright
    #    #
    #    # Frad
    #    #
    #    for ii in range(0,ix,1):
    #        dumleft     += 2. * ((src[ii]-src[<unsigned int> (1+ii)]) *
    #            expn(3., fabs(tau[ix]-tau[ii])))
    #    for ii in range(ix,ntau-1,1):
    #        dumright    += 2. * ((src[ii]-src[<unsigned int> (1+ii)]) *
    #            expn(3., tau[ix]-tau[ii]))
    #    Frad += dumleft + dumright
    #""""""
    return Jrad,Frad
""""""
#
# Call to outside
#
#cpdef calcJrad(double Tpre, double Tpost, ndarray[DTYPE_t, ndim=1] srcall, ndarray[DTYPE_t, ndim=1] tau, int ncpu = 2):
def calcJrad(Tpre, Tpost, srcall, tau, ncpu=2):
    """ 
    Calculate the mean itensity  and radiative flux
    """
    #
    # Variables
    #
    cdef int ip, ntau
    cdef double Ipre, Ipost, taumax, out1, out2
    Ipre    = (ss/M_PI)*pow(Tpre, 4.)
    Ipost   = (ss/M_PI)*pow(Tpost, 4.)
    taumax  = tau.max()
    ntau    = <unsigned int> tau.shape[0]
    #
    # Empty arrays
    #
    cdef ndarray [DTYPE_t, ndim=1] Jrad = np.zeros(tau.shape[0], dtype=DTYPE)
    cdef ndarray [DTYPE_t, ndim=1] Frad = np.zeros(tau.shape[0], dtype=DTYPE)
    #
    # Calculate Jrad
    #
    results = np.array(Parallel(n_jobs=ncpu, backend='multiprocessing')(
        delayed(getJrad)(ix, ntau, Ipre, Ipost, taumax, srcall, tau)
        for ix in xrange(tau.shape[0])))
    Jrad[:]    = results[:,0]
    Frad[:]    = results[:,1]
    return Jrad, Frad
""""""


