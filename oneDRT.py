import numpy as np
from python.natconst import *
from scipy.special import expn
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
    x = np.array(sol[:,0])
    dtau = np.zeros(x.shape[0])
    tau = np.zeros(x.shape[0])
    dx = np.zeros(x.shape[0]+1)
    dx[0] = 0.5*(x[1]-x[0])
    dx[1:-1] = x[1:]-x[:-1]
    
    # calculate this
    kapp = np.array([getKap(a) for a in sol[:,2]]) # gas
    kapd = sol[:,6]*np.pi*(sol[:,9]*sol[:,9])*getEps(sol[:,9])
    gasrho = np.sum(sol[:,[3,4,5,6]]*gas.mass, axis=1)
    dtau = (gasrho*kapp + kapd)*dx[1:]
    tau[-1] = 0.0
    for ix in xrange(dtau.shape[0]-2,-1,-1):
        tau[ix] = tau[ix+1]+ dtau[ix+1]
    taumax = tau.max()
    return taumax,tau, dtau
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

def calcJrad(Tpre, Tpost, tau, arrtau, dtau, taumax, ix=None, temps=None, sol=None, gas=None):
    """ 
    Calculate the mean itensity
    """
    Ipre = (ss/np.pi)*Tpre**(4.0)
    Ipost = (ss/np.pi)*Tpost**(4.0)
    if tau < 1e-5: tau=1e-5
    if expn(2.0,tau) < 1e-25:
        Jrad = 0.5*Ipre*expn(2.0, taumax-tau)
    else:
        Jrad = 0.5*Ipre*expn(2.0, taumax-tau) + 0.5*Ipost*expn(2.0,tau)
    """"""
    # Add the all intermediate terms
    #for ix in xrange(dtau.shape[0]-2, -1, -1):
        #delttau = np.abs(arrtau[ix]-tau)
        #if sol is None:
            #if delttau < 1e-25:
                #delttau = 1e-25
            #Jrad += 0.5*( (ss/np.pi)*temps[ix]**(4.0)*
                #expn(1,delttau)*dtau[ix])
        #else:
            #numdens = sol[ix,[3,4,5]]
            #rhogas = sum([a*b for (a,b) in zip(numdens,     
                #gas.mass)])
            #Srcg = (rhogas*getKap(sol[ix,2])*(ss/np.pi)*
                #sol[ix,2]**(4.))
            #Srcd = (sol[ix,6]*np.pi*sol[ix,9]*sol[ix,9]*
                #getEps(sol[ix,9])*(ss/np.pi)*sol[ix,8]**(4.))
            #Src = (Srcg + Srcd)/(rhogas*getKap(sol[ix,2]) +
                #np.pi*sol[ix,6]*sol[ix,9]*sol[ix,9]*
                #getEps(sol[ix,9]))
            ##Src = Srcg/(rhogas*getKap(sol[ix,2]))
            #if delttau < 1e-25:
                #delttau = 1e-25
            #addjrad = 0.5*(Src * expn(1,delttau)*dtau[ix])
            #if addjrad < 1e-20: addjrad=0.0
            #Jrad += addjrad
        #""""""
    """"""
    return Jrad
""""""
