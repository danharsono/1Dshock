import numpy as np
from python.my_header import *
from python.natconst import *; from oneDRT import getKap, getEps
from scipy.special import expn

"""
    Calculate the radiative flux
"""

def calcJrad(sol):
    """
        Test the radiative flux
    """
    #
    # Read the grid
    #
    grid    = np.empty(sol['grid'].shape, sol['grid'].dtype)
    sol['grid'].read_direct(grid)
    #
    # Gas data
    #
    gmasses = np.array([mp, 2.0*mp, 4.0*mp, 44.0*mp])
    gas     = np.empty(sol['gas'].shape, sol['gas'].dtype)
    sol['gas'].read_direct(gas)
    #
    # Dust data
    #
    dust    = np.empty(sol['dust'].shape, sol['dust'].dtype)
    sol['dust'].read_direct(dust)
    #
    # Calculate the optical depth
    #
    rhodust = dust[:,2]*3.3 # n*m
    rhogas  = (gas[:,2:]*gmasses[np.newaxis, :]).sum(axis=1)
    dtau    = np.zeros(grid.shape[0])
    tau     = np.zeros(grid.shape[0])
    dx      = np.zeros(grid.shape[0]+1)
    dx[0]   = 0.5*(grid[0]+grid[1])
    dx[-1]  = 0.5*(grid[-2]+grid[-1])
    dx[1:-1]= grid[1:]-grid[:-1]
    #
    # Kappa planck
    #
    kapp    = np.array([getKap(a) for a in gas[:,1]])
    galp    = rhogas*kapp # gas alpha
    deps    = np.array([getEps(a) for a in dust[:,2]])
    dalp    = dust[:,2]*np.pi*(dust[:,3]*dust[:,3])*deps # dust alpha
    dtau    = (galp+dalp)*dx[1:]
    tau     = np.zeros(grid.shape[0])
    tau[-1] = 0.
    for ix in xrange(dtau.shape[0]-2, -1, -1):
        tau[ix] = tau[ix+1] + dtau[ix+1]
    #
    # Get the taumax
    #
    taumax = tau.max()
    #
    # Calculate the radiation field
    #
    Tpre    = gas[0,1]
#    Tpost   = np.minimum(gas[-1,1], dust[-1,2])
    Tpost   = 1700.
    Ipre    = (ss/np.pi)*Tpre**(4.)
    Ipost   = (ss/np.pi)*Tpost**(4.)
    #
    # Source functions
    #
    srcall = (galp*(ss/np.pi)*np.power(gas[:,1],4.)+
        dalp*(ss/np.pi)*np.power(dust[:,1], 4.))/(rhogas*kapp+dalp)
    #
    # Calculate Jrad
    #
    Jrad = np.zeros(grid.shape[0])
    Frad = np.zeros(grid.shape[0])
    for ix in xrange(grid.shape[0]):
        Jrad[ix] = srcall[ix]
        Frad[ix] = 2. * srcall[ix]
        #
        # Add the other terms
        #
        Jrad[ix] += 0.5 * (Ipre - srcall[0])*expn(2.0, taumax-tau[ix])
        Frad[ix] = 2.*expn(3., taumax-tau[ix]) * (Ipre - srcall[0])
        Jrad[ix] += 0.5 * (Ipost - srcall[-1])*expn(2.0, tau[ix])
        Frad[ix] -= 2.*expn(3., tau[ix]) * (Ipost - srcall[-1])
        #
        # Use the derivative of the source functions
        #
        if ix == 0:
            Jrad[ix] += -0.5 * ((srcall[:-1]-srcall[1:])*expn(
                2., tau[ix]-tau[:-1])).sum()
            Frad[ix] += 2.0 * ((srcall[:-1]-srcall[1:])*expn(
                3., tau[ix]-tau[:-1])).sum()
        elif ix == grid.shape[0]-1:
            Jrad[ix] += 0.5 * ((srcall[:-1]-srcall[1:])*expn(
                2., np.abs(tau[ix]-tau[:-1]))).sum()
            Frad[ix] += 2. * ((srcall[:-1]-srcall[1:])*expn(
                3., np.abs(tau[ix]-tau[:-1]))).sum()
        else:
            """
                Need to split the positive and negative
            """
            dumleft     = 0.5 * ((srcall[:ix]-srcall[1:ix+1]) *expn(
                2., np.abs(tau[ix]-tau[:ix]))).sum()
            dumright    = -0.5 * ((srcall[ix:-1]-srcall[ix+1:])*expn(
                2., tau[ix]-tau[ix:-1])).sum()
            Jrad[ix] += dumleft + dumright
            #
            # Frad
            #
            dumleft     = 2. * ((srcall[:ix]-srcall[1:ix+1]) *expn(
                3., np.abs(tau[ix]-tau[:ix]))).sum()
            dumright    = 2. * ((srcall[ix:-1]-srcall[ix+1:])*expn(
                3., tau[ix]-tau[ix:-1])).sum()
            Frad[ix] += dumleft + dumright
        """"""
    """"""
    semilogy(grid, Jrad,'r-')
    semilogy(grid,(ss/np.pi)*np.power(gas[:,1],4.),'b-.')
    semilogy(grid,srcall,'k--')
    show()
    close()

    print
    print 'Tpre: %d  Tpost: %d'%(Tpre, Tpost)
    print 'Frad: %2.3e -- %2.3e'%(Frad[0], Frad[-1])

    semilogy(grid, np.abs(Frad),'r-')
    ylim(1e-10, 1e10)
    show()
    close()
    raise SystemExit

""""""

