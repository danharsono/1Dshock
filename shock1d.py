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
    if (len(w) > 3) and (len(p) == 1):
        x1, x2, n1, n2, n3 = w
        gas = p[0]
        # update the gas velocities
        gas._updateGas(allns=[n1/x1, n2/x1, n3/x1])
        varA = 2.0*gas._sumRho()*x1 - (kk*x2/x1)*sum(gas.numden)
        varB = sum(gas.numden)*kk
        varC = -gas._calculateVarC(vel=x1, t=x2)
        varD = (gas._sumEnergy()*kk*x2 + 
            gas._sumRho()*x1*x1*1.5 + 
            (kk*gas._sumRho()/(gas.mugas*mp)) * x2)
        varE = (kk*gas._sumEnergy() +
            (kk*gas._sumRho()/(gas.mugas*mp)))*x1
        varF = (-gas._calculateVarF(vel=x1, t=x2) + 
                gas._calculateFreeEnergy(t=x2))
        #print x2, varF, gas._calculateFreeEnergy(t=x2)
        #raise SystemExit
        f = [(varC*varE - varB*varF)/(varA*varE - varB*varD),
            (varA*varF - varD*varC)/(varA*varE - varB*varD), 
            -2.0*gas._calculateR(t=x2), gas._calculateR(t=x2), 0.0]
    elif (len(w) > 6) and (len(p) == 2):
        """
        This means that there is dust
        """
        x1, x2, n1, n2, n3, nd, vd, td, ad = w
        gas, dust = p
        nd = np.abs(nd)
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
        mdust = (4.0/3.0)*np.pi*ad*ad*ad*dust._sumRho()
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
        #
        # Check few numbers to make sure these are not too large
        # or too small
        f1 = (varC*varE - varB*varF)/(varA*varE - varB*varD)
        f2 = (varA*varF - varD*varC)/(varA*varE - varB*varD)
        f3 = -2.0*gas._calculateR(t=x2)
        f4 = gas._calculateR(t=x2)
        f5 = 0.0
        f6 = -(nd/vd)*(fdrag/(mdust*vd)) # This is for number of dust
        f7 = fdrag/(mdust*vd) # Dust velocities
        f8 = Dtd
        f9 = dxa
        f = [f1, f2, f3, f4, f5, f6, f7, f8, f9]
        #
        #
        #
        #print '%2.2e %2.2e %d %2.2e %2.2e %2.2e %d'%(
            #x, x1, x2, n1, nd, vd, td)
        ###print '%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e'%(
            ###x1, x2, td, nd, f7, f6)
        ###
        #print '%2.2e %2.2e %2.2e'%(varA, varB, varC)
        #print '%2.2e %2.2e %2.2e'%(varD, varE, varF)
        #print 'Fdrag: %2.3e'%(fdrag)
        #print 'Dtd: %2.3e'%(Dtd)
        #print 'Dvd: %2.3e'%(f7)
        #print 'Dxa: %2.3e'%(dxa)
        #print '%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e'%(
            #f1, f2, f6, f7, f8, f9)
        #print
        if x2 < 10.0 or (td < 10.0):
            print 'Negative temperatures!'
            raise SystemExit
        #if nd < 1e-15:
            #print 'Negative dust...'
            #raise SystemExit
        if x > -8.9e9:
            raise SystemExit
    else:
        x1, x2, n1 = w
        gas = p[0]
        gas._updateGas(allns=[n1])
        #varA = 2.0*n1*gas.mass[0]*x1 - (kk*x2/x1)*n1
        varA = 2.0*gas._sumRho()*x1 - (kk*x2/x1)*n1
        varB = n1*kk
        varC = -gas._calculateVarC(vel=x1, t=x2)
        varD = (gas._sumEnergy()*kk*x2 + 
            gas._sumRho()*x1*x1*1.5 + 
            (kk*gas._sumRho()/(gas.mugas*mp)) * x2)
        varE = (kk*gas._sumEnergy() +
            (kk*gas._sumRho()/(gas.mugas*mp)))*x1
        varF = 0.0
        f = [(varC*varE - varB*varF)/(varA*varE - varB*varD),
            (varA*varF - varD*varC)/(varA*varE - varB*varD), 
            0.0]
    """"""
    # 
    # Update the gas properties
    #
    #gas._updateGas(allns=[n1])
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
