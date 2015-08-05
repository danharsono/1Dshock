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
    gas, dust, debug = p
    x1, x2, n1, n2, n3, n4, nd, vd, td, ad = w
    #
    # Calculate the variable here
    #
    varA = gas._sumRho()*x1 - (kk*x2/x1)*sum(gas.numden)
    varB = sum(gas.numden)*kk
    #
    # Variable C now has dust properties
    #
    varC = -gas._calculateVarC(vel=x1, t=x2, veld=vd, 
        td=td, dust=dust)
    #
    # Calculate the variables for energy equation
    #
    #    varD = (gas._sumEnergy()*kk*x2 +
    #        gas._sumRho()*x1*x1*1.5 +
    #        (kk*gas._sumRho()/(gas._getMugas()*mp)) * x2)
    #    varE = (kk*gas._sumEnergy() +
    #        (kk*gas._sumRho()/(gas._getMugas()*mp)))*x1
    varD = x1*x1*gas._sumRho()
    varE = kk*x1*gas._sumGammas()
    #
    # Calculate the variable F with the dust
    #
    varF = (-gas._calculateVarF(vel=x1, t=x2, veld=vd, 
        td=td, dust=dust) + 
        gas._calculateFreeEnergy(t=x2))
    #
    # Drag forces
    #
    fdrag = dust._calculateFdrag(vd=vd, vg=x1, Tg=x2,
        Td=td, gas=gas)
    #
    # Dust mass
    #
    mdust = (4.0/3.0)*np.pi*ad*ad*ad*dust.mdust
    #
    # Dust temperature 
    #
    Dtd = dust._calcDxTd(vd=vd, vg=x1, Tg=x2, Td=td, 
        gas=gas)
    #
    # Dust size change
    #
    dxa = dust._calcDxa(vd=vd, vg=x1, Tg=x2, Td=td,
         gas=gas)
    #
    # The RHS matrix
    # Check few numbers to make sure these are not too large
    # or too small
    #
    f1 = (varC*varE - varB*varF)/(varA*varE - varB*varD)
    f2 = (varA*varF - varD*varC)/(varA*varE - varB*varD)
    f3 = -2.0*gas._calculateR(t=x2)
    f4 = gas._calculateR(t=x2)
    f5 = 0.0
    f6 = -(nd*vd*4.0*np.pi*dust._sumRho()*ad*ad)/(gas.mass[3])*dxa
    #
    # DUST from here on
    #
    if dust.nspecs == 0:
        f7 = 0.0
        f8 = 0.0
        f9 = 0.0
        f10 = 0.0
    else:
        f7 = -(nd/vd)*(fdrag/(mdust*vd)) # This is for number of dust
        f8 = fdrag/(mdust*vd) # Dust velocities change
        f9 = Dtd
        f10 = dxa
#        if np.abs(f10) > 0.0:
#            f9 = 0.0
#            f10 = dxa
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
    if np.abs(f1/x1) < 1e-15: f1 = 0.0
    if np.abs(f2/x2) < 1e-15: f2 = 0.0
    if np.abs(f3/n1) < 1e-15: f3 = 0.0
    if np.abs(f4/n2) < 1e-15: f4 = 0.0
    if np.abs(f5/n3) < 1e-15: f5 = 0.0
    if (n4 > 1e-15 and (np.abs(f6/n4) < 1e-15)): f6 = 0.0
    if np.abs(f7/nd) < 1e-15: f7 = 0.0
    if np.abs(f8/vd) < 1e-15: f8 = 0.0
    if np.abs(f9/td) < 1e-15: f9 = 0.0
    if np.abs(f10/ad) < 1e-15: f10 = 0.0
    f = np.array([f1, f2, f3, f4, f5, f6, f7, f8, f9, f10])
    if nd < 0.0:
        print 'Tgas and vgas: %e, %e'%(x1*1e-5, x2)
        print 'ndust, vdust and Tdust: %2.5e  %2.5e  %d'%(nd, vd, td)
        print '%e %e  %e  %e'%(f6, f7, f8, f9)
        print
        print fdrag, mdust, vd, nd/vd
        print
        print f
        print w
        print
        print 'DUST is negative!'
        raise SystemExit
    if debug or np.abs(f2) > 1e4:
        print 'Densities: %2.5e  %2.5e'%(gas._sumRho(), dust._sumRho())
        print 'Masses: ', gas.mass, dust.mass
        print varA, varB, varC, varD, varE, varF
        print
        print 'Tgas and vgas: %e, %e'%(x1*1e-5, x2)
        print 'ndust, vdust and Tdust: %2.5e  %2.5e  %e'%(nd, vd, td)
        print 'Fdrag: %2.5e'%(fdrag)
        print 'Dtd: %2.5e'%(Dtd)
        print 'Dxa: %2.5e'%(dxa)
        print '%e %e  %e  %e'%(f6, f7, f8, f9)
        print
        print fdrag, mdust, nd/vd
        print f
        print w
        print
        time.sleep(0.05)
#        raise SystemExit
    """"""
    return f
""""""

"""
The vector field of the matrix to solve
with radiative transfer
"""
def vectorfieldrt(x, w, p):
    """
    w   :  the variables to solve -> x, y 
    x   :  position with respect to the shock
    p   :  parameters
    """
    from oneDRT import getKap,getEps
    gas, dust, Jrad, debug = p
    x1, x2, n1, n2, n3, n4, nd, vd, td, ad = w
    #
    # Interpolate to find the current Jrad
    #
    curJ = np.interp(x, Jrad[0], Jrad[1])
    #
    # Calculate the variable here
    #
    varA = 2.0*gas._sumRho()*x1 - (kk*x2/x1)*sum(gas.numden)
    varB = sum(gas.numden)*kk
    #
    # Variable C now has dust properties
    #
    varC = -gas._calculateVarC(vel=x1, t=x2, veld=vd, 
        td=td, dust=dust, Jrad=curJ)
    #
    # Calculate the variables for energy equation
    #
    varD = (gas._sumEnergy()*kk*x2 + 
        gas._sumRho()*x1*x1*1.5 + 
        (kk*gas._sumRho()/(gas.mugas*mp)) * x2)
    varE = (kk*gas._sumEnergy() +
        (kk*gas._sumRho()/(gas.mugas*mp)))*x1
    #
    # Calculate the variable F with the dust
    #
    varF = ( -gas._calculateVarF(vel=x1, t=x2, veld=vd, 
            td=td, dust=dust, Jrad=curJ) + 
            gas._calculateFreeEnergy(t=x2) )
    if curJ is not None:
        radiation = (4.0*np.pi*gas._sumRho()*getKap(temp=x2) *
             (curJ - (ss*np.power(x2, 4.)/np.pi)) +
             4.0*np.pi*nd*ad*ad*getEps(ad)*
             (curJ-(ss*np.power(td,4.)/np.pi)) )
        varF += (4.0*np.pi*gas._sumRho()*getKap(temp=x2) *
            (curJ - (ss*x2**(4.)/np.pi)) + 
            4.0*np.pi*nd*ad*ad*getEps(ad)*
            (curJ-(ss*td**(4.)/np.pi)) )
    #
    # Drag forces
    #
    fdrag = dust._calculateFdrag(vd=vd, vg=x1, Tg=x2,
        Td=td, gas=gas)
    #
    # Dust mass
    #
    mdust = (4.0/3.0)*np.pi*ad*ad*ad*dust.mdust
    #
    # Dust temperature 
    #
    Dtd = dust._calcDxTd(vd=vd, vg=x1, Tg=x2, Td=td, 
        gas=gas, Jrad=curJ)
    #
    # Dust size change
    #
    dxa = dust._calcDxa(vd=vd, vg=x1, Tg=x2, Td=td,
         gas=gas, Jrad=curJ)
    #
    # The RHS matrix
    #
    # Check few numbers to make sure these are not too large
    #
    # or too small
    f1 = (varC*varE - varB*varF)/(varA*varE - varB*varD)
    f2 = (varA*varF - varD*varC)/(varA*varE - varB*varD)
    f3 = -2.0*gas._calculateR(t=x2)
    f4 = gas._calculateR(t=x2)
    f5 = 0.0
    f6 = -(nd*vd*4.0*np.pi*dust._sumRho()*ad*ad)/(gas.mass[3])*dxa
#    print x, curJ, radiation, x1, x2, vd, td, fdrag, dxa, Dtd
#    raise SystemExit
    #
    # DUST from here on
    #
    if dust.nspecs == 0:
        f7 = 0.0
        f8 = 0.0
        f9 = 0.0
        f10 = 0.0
    else:
        f7 = -(nd/vd)*(fdrag/(mdust*vd)) # This is for number of dust
        f8 = fdrag/(mdust*vd) # Dust velocities change
        f9 = Dtd
        f10 = dxa
        if np.abs(f10) > 0.0:
            f9 = 0.0
            f10 = dxa
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
        Check if f7 ~ nd
    """
    if (np.abs(f7) >= np.abs(nd)):
        print 'Tgas and vgas: %e, %e'%(x1*1e-5, x2)
        print 'ndust, vdust and Tdust: %2.5e  %2.5e  %d'%(nd, vd, td)
        print '%e %e  %e  %e'%(f6, f7, f8, f9)
        print
        print fdrag, mdust, vd, nd/vd
        print w
        print
        raise SystemExit
    """
        Check if gas number densities change by quite a lot
    """
    if np.abs(f1/x1) < 1e-15: f1 = 0.0
    if np.abs(f2/x2) < 1e-15: f2 = 0.0
    if np.abs(f3/n1) < 1e-15: f3 = 0.0
    if np.abs(f4/n2) < 1e-15: f4 = 0.0
    if np.abs(f5/n3) < 1e-15: f5 = 0.0
    if (n4 > 1e-15 and (np.abs(f6/n4) < 1e-15)): f6 = 0.0
    if np.abs(f7/nd) < 1e-15: f7 = 0.0
    if np.abs(f8/vd) < 1e-15: f8 = 0.0
    if np.abs(f9/td) < 1e-15: f9 = 0.0
    if np.abs(f10/ad) < 1e-15: f10 = 0.0
    """
        Limit the changes with respect to original value
    """
    f = np.array([f1, f2, f3, f4, f5, f6, f7, f8, f9, f10])
    """ Do something below here """
    if nd < 1e-19:
        print 'Tgas and vgas: %e, %e'%(x1*1e-5, x2)
        print 'ndust, vdust and Tdust: %2.5e  %2.5e  %d'%(nd, vd, td)
        print '%e %e  %e  %e'%(f6, f7, f8, f9)
        print
        print fdrag, mdust, vd, nd/vd
        print
        print f
        print w
        print
        print 'DUST is negative!'
        raise SystemExit
    if debug:
        print 'Densities: %2.5e  %2.5e'%(gas._sumRho(), dust._sumRho())
        print 'Masses: ', gas.mass, dust.mass
        print varA, varB, varC, varD, varE, varF
        print
        print 'Tgas and vgas: %e, %e'%(x1*1e-5, x2)
        print 'ndust, vdust and Tdust: %2.5e  %2.5e  %e'%(nd, vd, td)
        print 'Fdrag: %2.5e'%(fdrag)
        print 'Dtd: %2.5e'%(Dtd)
        print 'Dxa: %2.5e'%(dxa)
        print '%e %e  %e  %e'%(f6, f7, f8, f9)
        print
        print fdrag, mdust, nd/vd
        print f
        print w
        print
        raise SystemExit
        time.sleep(0.05)
    """"""
    return f
""""""
