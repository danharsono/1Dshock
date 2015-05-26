import numpy as np
from python.natconst import *

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
    x1, x2, n1, n2, n3, nd, vd, td, ad = w
    gas, dust, debug = p
    #
    # Calculate the variable here
    #
    varA = 2.0*gas._sumRho()*x1 - (kk*x2/x1)*sum(gas.numden)
    varB = sum(gas.numden)*kk
    #
    # Variable C now has dust properties
    #
    varC = -gas._calculateVarC(vel=x1, t=x2, veld=vd, 
        td=td, dust=dust)
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
    #
    # Check few numbers to make sure these are not too large
    #
    # or too small
    f1 = (varC*varE - varB*varF)/(varA*varE - varB*varD)
    f2 = (varA*varF - varD*varC)/(varA*varE - varB*varD)
    f3 = -2.0*gas._calculateR(t=x2)
    f4 = gas._calculateR(t=x2)
    f5 = 0.0
    #
    # DUST from here on
    #
    if dust.nspecs == 0:
        f8 = 0.0
        f9 = 0.0
        f7 = 0.0
        f6 = 0.0
    else:
        f6 = -(nd/vd)*(fdrag/(mdust*vd)) # This is for number of dust
#        f6 = 0.0
        f7 = fdrag/(mdust*vd) # Dust velocities change
        if np.abs(f7) > (vd * 1e-5):
            if dust.nspecs == 0:
                f7 = 0.0
                raise SystemExit
                #
                # This should not be here anymore
                #
        if np.isnan(f7):
            if dust.nspecs == 0:
                f7 = 0.0
                raise SystemExit
                #
                # This should not be here anymore
                #
            else:
                print 'NAN here in dust density change'
                raise SystemExit
        """"""
        f8 = Dtd
        f9 = dxa
    """"""
    f = np.array([f1, f2, f3, f4, f5, f6, f7, f8, f9])
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
    if debug:
        print 'Densities: %2.5e  %2.5e'%(gas._sumRho(), dust._sumRho())
        print 'Masses: ', gas.mass, dust.mass
        print
        print 'Tgas and vgas: %e, %e'%(x1*1e-5, x2)
        print 'ndust, vdust and Tdust: %2.5e  %2.5e  %d'%(nd, vd, td)
        print '%e %e  %e  %e'%(f6, f7, f8, f9)
        print
        print fdrag, mdust, nd/vd
        print f
        print w
        print
#        raise SystemExit
    """"""
    f[(np.abs(f/w)<1e-15)] = 0.0
    #
    # Need to make sure that the change in v and t is not too much
    #
    if (np.abs(f1) > x1*1e5) or (np.abs(f2) > x2*1e5):
        f[0] = 0.0
        f[1] = 0.0
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
    if (len(w) > 3) and (len(p) == 2):
        x1, x2, n1, n2, n3= w
        if (x2 > 1e8):
            print 
            print '%2.3e,  %2.3e,  %2.3e'%(x, x1,x2)
            print '%d, %2.3e, %2.3e'%(idx, grid[idx], grid.max())
            print 'tau: %e  taumax: %2.1e'%(tau[idx], taumax)
            print '%2.4e'%(varF)
            print 'Jrad: %2.3e'%(calcJrad(Tpre, Tpost, tau[idx], 
                    tau,dtau, temps, taumax))
            print 'Too high temperatures!!'
            raise SystemExit
        """"""
        gas = p[0]
        Jrad = p[1]
        #
        # Interpolate to find the current Jrad
        #
        curJ = np.interp(x, Jrad[0], Jrad[1])
        # update the gas velocities
        gas._updateGas(allns=[n1/x1, n2/x1, n3/x1])
        #
        # Calculate the different variables
        #
        varA = 2.0*gas._sumRho()*x1 - (kk*x2/x1)*sum(gas.numden)
        varB = sum(gas.numden)*kk
        varC = -gas._calculateVarC(vel=x1, t=x2)
        varD = (gas._sumEnergy()*kk*x2 + 
            gas._sumRho()*x1*x1*1.5 + 
            (kk*gas._sumRho()/(gas.mugas*mp)) * x2)
        varE = (kk*gas._sumEnergy() +
            (kk*gas._sumRho()/(gas.mugas*mp)))*x1
        varF = (-gas._calculateVarF(vel=x1, t=x2) + 
            gas._calculateFreeEnergy(t=x2) +
            4.0*np.pi*gas._sumRho()*getKap(temp=x2)*(
            curJ - (ss*x2**(4.)/np.pi)))
        f = [(varC*varE - varB*varF)/(varA*varE - varB*varD),
            (varA*varF - varD*varC)/(varA*varE - varB*varD), 
            -2.0*gas._calculateR(t=x2), gas._calculateR(t=x2), 0.0]
    else:
        """
        This means that there is dust
        """
        x1, x2, n1, n2, n3, nd, vd, td, ad = w
        if (x2 < 5.0): x2 = 25.0
        if (td < 5.0): td = 25.0
        if (x1 < 1e3): x1 = 1e4
        if (nd < 1e-11): nd = 1e-8
        gas, Jrad, dust = p
        #
        # Interpolate to find the current Jrad
        #
        curJ = np.interp(x, Jrad[0], Jrad[1])
        #curJ = None
        # update the gas velocities
        gas._updateGas(allns=[n1/x1, n2/x1, n3/x1])
        dust._updateDust(allns=[nd], size=ad)
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
            varF += (4.0*np.pi*gas._sumRho()*getKap(temp=x2) *
                (curJ - (ss*x2**(4.)/np.pi)) + 
                4.0*np.pi*nd*ad*ad*getEps(ad)*
                (curJ-(ss*td**(4.)/np.pi)) )
        #
        # Drag forces
        #
        fdrag = dust._calculateFdrag(vd=vd, vg=x1, Tg=x2,
            Td=td, rhogas=gas._sumRho())
        #
        # Dust mass
        #
        mdust = (4.0/3.0)*np.pi*ad*ad*ad*dust._sumRho()
        #
        # Dust temperature 
        #
        Dtd = dust._calcDxTd(vd=vd, vg=x1, Tg=x2, Td=td, 
            rhogas=gas._sumRho(), Jrad=curJ)
        #
        # Dust size change
        #
        dxa = dust._calcDxa(vd=vd, vg=x1, Tg=x2, Td=td,
             rhogas=gas._sumRho(), Jrad=curJ)
        #
        # The RHS matrix
        #
        f = [ (varC*varE - varB*varF)/(varA*varE - varB*varD),
             (varA*varF - varD*varC)/(varA*varE - varB*varD), 
             -2.0*gas._calculateR(t=x2), 
             gas._calculateR(t=x2),
             0.0, 
             -nd/vd*(fdrag/(mdust*vd)), 
             fdrag/(mdust*vd),
             Dtd, 
             dxa
             ]
        """ Do something below here """
        #print '%1.2e %1.2e, %.1f, %1.2e, %1.2e, %.1f, %1.2e'%(
            #x, x1,x2,nd,vd,td,ad)
        #print varA, varB, varD, varE
        #print varC, varF, Dtd, '%2.3f'%(np.abs(x1-vd))
        #print f
        #print
        #raise SystemExit
        if (x1 > 1e6) or (x2 > 1e4):
            print '%1.2e %1.2e, %.1f, %1.2e, %1.2e, %.1f, %1.2e'%(
                x, x1,x2,nd,vd,td,ad)
            print varA, varB, varD, varE
            print varC, varF, Dtd, '%2.3f'%(np.abs(x1-vd))
            print f
            print            
        if (x1 > 1e9) or (x2 > 1e8):
            print
            print '%1.2e %1.2e, %.1f, %1.2e, %1.2e, %.1f, %1.2e'%(
                x, x1,x2,nd,vd,td,ad)
            print varA, varB, varD, varE
            print varC, varF, Dtd, '%2.3f'%(np.abs(x1-vd))
            print f
            print
            raise SystemExit
    """"""
    """ Check to make sure the values are ok """
    if (x2 < 0):
        print
        print '%2.3e,  %2.3e,  %2.3e'%(x, x1,x2)
        print '%2.4e'%(varF)
        print '%2.4e'%(curJ)
        print 'Negative temperatures!!'
        raise SystemExit
    return f
""""""
