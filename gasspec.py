from python.natconst import *; from numpy import loadtxt
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
#
# Read in the species
#
def readSpecs(infile='species.dat'):
    """
    Read the species.dat
    return it
    """
    dum = loadtxt(infile, dtype={
        'names': ['id', 'species', 'frac', 'mass', 'gam', 'dustfrac'],
        'formats': ['i', 'S5', 'd', 'd', 'd', 'd']})
    specs = np.array([a for a in dum])
    return specs
""""""
#
# Indexer
#
def getIndx(a, specs):
    if a == 'PHOTON':
        return 'PH'
    elif a == 'CRP':
        return 'CRP'
    elif a == 'CRPHOTON':
        return 'CRPH'
    elif a == 'ACC': # accrete
        return 'ACC'
    else:
        try:
            return (a== specs['species']).nonzero()[0][0]
        except IndexError:
            return -99
    """"""
""""""
#
# Read in the rates: 2 body reaction only
#
def readRates(specs, infile='rates.dat'):
    """
    Read the rates.dat file and create
    """
    fin     = open('rates.dat','r')
    rates   = []
    #
    # This has to be 7 species
    #
    for line in fin:
        columns = line.strip().split()
        ncol    = len(columns)
        rate    = {}
        if ncol == 11:
            rate['R1']      = getIndx(columns[1], specs)
            rate['R2']      = getIndx(columns[2], specs)
            rate['R3']      = getIndx(columns[3], specs)
            rate['P1']      = getIndx(columns[4], specs)
            rate['P2']      = getIndx(columns[5], specs)
            rate['P3']      = getIndx(columns[6], specs)
            rate['P4']      = getIndx(columns[7], specs)
            rate['alp']     = float(columns[8])
            rate['beta']    = float(columns[9])
            rate['gamma']   = float(columns[10])
        elif ncol == 10:
            rate['R1']      = getIndx(columns[1], specs)
            rate['R2']      = getIndx(columns[2], specs)
            rate['R3']      = getIndx(columns[3], specs)
            rate['P1']      = getIndx(columns[4], specs)
            rate['P2']      = getIndx(columns[5], specs)
            rate['P3']      = getIndx(columns[6], specs)
            rate['P4']      = -99
            rate['alp']     = float(columns[7])
            rate['beta']    = float(columns[8])
            rate['gamma']   = float(columns[9])
        elif ncol == 9:
            rate['R1']      = getIndx(columns[1], specs)
            rate['R2']      = getIndx(columns[2], specs)
            rate['R3']      = -99
            rate['P1']      = getIndx(columns[3], specs)
            rate['P2']      = getIndx(columns[4], specs)
            rate['P3']      = getIndx(columns[5], specs)
            rate['P4']      = -99
            rate['alp']     = float(columns[6])
            rate['beta']    = float(columns[7])
            rate['gamma']   = float(columns[8])
        elif ncol == 8:
            rate['R1']      = getIndx(columns[1], specs)
            rate['R2']      = getIndx(columns[2], specs)
            rate['R3']      = -99
            rate['P1']      = getIndx(columns[3], specs)
            rate['P2']      = getIndx(columns[4], specs)
            rate['P3']      = -99
            rate['P4']      = -99
            rate['alp']     = float(columns[5])
            rate['beta']    = float(columns[6])
            rate['gamma']   = float(columns[7])
        elif ncol == 7: # ICE
            rate['R1']      = getIndx(columns[1], specs)
            rate['R2']      = getIndx(columns[2], specs)
            rate['R3']      = -99
            rate['P1']      = getIndx(columns[3], specs)
            rate['P2']      = -99
            rate['P3']      = -99
            rate['P4']      = -99
            rate['alp']     = float(columns[4])
            rate['beta']    = float(columns[5])
            rate['gamma']   = float(columns[6])
        else:
            if columns[0] == '9999':
                break
            else:
                print 'UNKNOWN text format! %d \n'%(ncol)
                print line, columns[0]
                print 'Check input file\n'
                raise SystemExit
        """"""
        rates.append(rate)
    """"""
    fin.close()
    del rate
    #
    # Get the indices and rates directly
    #
    rates1 = np.zeros((len(rates),11))
    for irate in range(len(rates)):
        """
        Options for rate reactions
        """
        rates1[irate,0]         = 0.0
        if rates[irate]['R2'] == 'CRP':
            rates1[irate, 0]    = 1.
            rates1[irate, 2]    = -99
        elif rates[irate]['R2'] == 'CRPH':
            rates1[irate, 0]    = 2.
            rates1[irate, 2]    = -99
        elif rates[irate]['R2'] == 'PH':
            rates1[irate, 0]    = 3.
            rates1[irate, 2]    = -99
        elif rates[irate]['R2'] == 'ACC':
            rates1[irate, 0]    = 4.
            rates1[irate, 2]    = -99
        else:
            rates1[irate, 2]    = rates[irate]['R2']
        """"""
        rates1[irate, 1]    = rates[irate]['R1']
        rates1[irate, 3]    = rates[irate]['R3']
        rates1[irate, 4]    = rates[irate]['P1']
        try:
            rates1[irate, 5]    = rates[irate]['P2']
        except ValueError:
            rates1[irate, 5]    = -99
        rates1[irate, 6]    = rates[irate]['P3']
        rates1[irate, 7]    = rates[irate]['P4']
        rates1[irate, 8]    = rates[irate]['alp']
        rates1[irate, 9]    = rates[irate]['beta']
        rates1[irate, 10]   = rates[irate]['gamma']
    """"""
    rates = []
    del rates
    return np.array(rates1)
""""""
#
# Gas
#
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
        #
        # Updated version: read in species.dat
        #
        dum     = readSpecs()
        nspecs  = len(dum)
        #
        # Put it into the object
        #
        self.nspecs     = nspecs
        self.mass       = [a*mp for a in dum['mass']]
        self.gamma      = [a/2. for a in dum['gam']]
        self.dustfrac   = [a for a in dum['dustfrac']]
        #
        # get columns
        #
        self.specfrac   = [a for a in dum['frac']]
        n0              = self.rho/sum(self.specfrac)
        self.numden     = [a*n0/b for (a,b) in zip(
            self.specfrac, self.mass)]
        self.specs      = dum
        del dum
        #
        # Read the rates
        #
        self.rates      = readRates(self.specs)
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
        try:
            kap = np.power(10., np.interp(np.log10(t), self.logT, self.logK))
        except AttributeError:
            kap =0.0
        if destroyed and t < 1.4e3:
            kap =0.5
        return kap
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