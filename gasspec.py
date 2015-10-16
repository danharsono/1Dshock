from python.natconst import *
import numpy as np; from copy import deepcopy
from scipy.integrate import ode

"""
The gas components for the shock code
"""

def solveChem(x,w,p):
    """
        Solve the chemistry
    """
    gas, t = p
    #
    # calculate the changes
    #
    dx = gas._calculateChem(t=t)
    return dx
""""""

class gasSpecs():
    """ 
    Constructor
    Destructor is not required in python
    """
    def __init__(self, rho = 1e-9, vgas=6.5e5, tgas = 300., nspecs=3):
        self.rho = rho
        self.vgas = vgas
        self.tgas = tgas
        self.logT   = np.zeros(100)
        self.logK   = np.zeros(100)
        if nspecs > 0:
            self.nspecs = nspecs
        else:
            print 
            print 'Error initializing gas species'
            print 'Nspecs: %d'%(nspecs)
            print 'Nspecs must be > 0'
            print
            raise SystemExit
        """"""
        if nspecs == 1:
            """
            Only use He
            """
            self.nspecs = nspecs
            self.specfrac = [1.]
            self.mass = [4.0*mp]
            self.numden = [self.rho/a for a in self.mass]
            self.gamma = [5.0/2.0]
            self.mugas = 2.4
            self.dustfrac = [0]
        elif nspecs == 3:
            """
            Use H, H2 and He
            """
            self.nspecs =3
            self.mass = [mp, 2.0*mp, 4.0*mp]
            self.specfrac = [0.01, 0.74, 0.25]
            mu = 0.71 + 0.27/4.
            self.mugas = mu**(-1.0)
            n0 = self.rho/sum(self.specfrac)
            self.numden = [a*n0/b for (a,b) in zip(self.specfrac,
                self.mass)]
            self.gamma = [5.0/2.0, 7.0/2.0, 5.0/2.0]
            self.dustfrac = [0, 0, 0]
        elif nspecs == 4:
            """
            Add SiO
            """
            self.nspecs=4
            self.mass = [mp, 2.0*mp, 4.0*mp, 44.0*mp]
            self.gamma = [5.0/2.0, 7.0/2.0, 5.0/2.0, 7.0/2.]
            self.dustfrac = [0., 0., 0., 1.]
            #
            # get columns
            #
            self.mugas = 1./(0.83*(0.5 + 0.2*0.25))
            self.specfrac = [1e-4, 0.83, 0.16656,0.0]
            n0 = self.rho/sum(self.specfrac)
            self.numden = [a*n0/b for (a,b) in zip(self.specfrac,
                self.mass)]
            print ['%2.5e'%(a) for a in self.numden]
            #
            # Evolve this until steady state
            #
#            temp = self._getSteady()
#            self.numden = [a for a in temp]
    """"""
    def getInps(self):
        """
        Create a list of the input parameters for the gas
        """
        return [self.vgas, self.tgas]+[a for a in self.numden]
    """"""
    def _getOpacs(self):
        """ 
        Get gas planck opacities
        """
        data = np.genfromtxt('planckOpacs.data')
        self.logT = np.log10(data[:,0])
        self.logK = np.log10(data[:,1])
    """"""
    def _getKap(self, t, destroyed=False):
        """
        Get opacities
        """
        if destroyed:
            return 0.5
        else:
            try:
                return 10.0**(np.interp(t, self.logT, self.logK))
            except AttributeError:
                return 0.0
    """"""
    def _calculateR(self, vars=None):
        """
        Calculate the reactions and update the values
        This rate is the creation of gas species i 
        in terms of number of reactions per unit time 
        per unit volume
        """
        t = vars[1]
        rate1 = 8.72e-33*(t/300.0)**(-0.6) # cm^6 per sec
        rate2 = 1.83e-31*(t/300.0)**(-1) # cm^6 per sec
        rate3 = 1.50e-9*np.exp(-46350./t)
        rate4 = 3.75e-8*(t/300.0)**(-0.5)*np.exp(-53280./t)
        #
        # Create the reaction matrix
        #
        temprates = np.zeros((self.nspecs, self.nspecs), dtype=np.float64)
        if self.nspecs > 1:
            totrate = (vars[2]*vars[2]*(rate1*vars[3] +
                rate2*vars[2]) - vars[3]*(rate3*vars[3]+rate4*vars[2]) )
        else:
            totrate = 0.0
        temprates[0,0] = -2. * totrate
        temprates[0,1] = totrate
        return temprates
    """"""
    def _calculateFreeEnergy(self, w=None):
        """
        Calculate the net energy
        H + H + energy -> H2 
        H2 + energy -> H + H
        """
        normrate = self._calculateR(vars=w)
        onev = 1.6021772e-12
        return sum(normrate.sum(axis=0))*(-4.48*onev)
    """"""
    def _getChemInput(self):
        """
        Create the input for the chemistry
        """
        return [a for a in self.numden]
    """"""
    def _calculateChem(self, t):
        """
            Calculate the rates
        """
        #
        # ToDO: match the species and rates
        # Assume H, H2, He, SiO
        #
        rate = self._calculateR(vars=self.getInps())
        return rate.sum(axis=0)
    """"""
    def _getSteady(self):
        """
        Solve the rates to get steady state solution and abundances at 
        fixed temperatures
        """
        #
        # INPUTS: w0 and p
        #
        p = [self, self.tgas]
        w0 = self._getChemInput()
        dt = 1e5
        t = 0.0
        #
        # Setup the vode
        #
        vode = ode(solveChem).set_integrator('vode', atol=1e-5,
            rtol=1e-5, order=5, method='bdf',nsteps=1e5,
            first_step = 1e1, max_step=1e5,
            with_jacobian=True)
        #
        # stop it when this occurs
        #
        reltoll = 1e-4
        reached = True
        while(reached):
            vode.set_initial_value(w0, t).set_f_params(p)
            #
            # store old values
            #
            oldy = vode.y
            #
            # Integrate
            #
            vode.integrate(vode.t+dt)
            t = vode.t
            #
            # Set the new values
            #
            self._updateGas(allns=deepcopy(vode.y), tgas=self.tgas)
            w0 = self._getChemInput()
            p = [self, self.tgas]
            #
            # Check changes
            #
            dy = max(np.abs(vode.y-oldy)/(vode.y+1e-15))
            dx = max(self._calculateChem(t=self.tgas))
            if np.isnan(dx): raise SystemExit
            if (dy < reltoll) or (dx < reltoll):
                reached = False
        """"""
        return w0
""""""