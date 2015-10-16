from python.natconst import *
import numpy as np
import warnings

"""
The gas components for the shock code
"""
class dustSpecs():
    """ 
    Constructor
    Destructor is not required in python
    """
    def __init__(self, dustfrac = 0.01, nspecs=1, mdust=3.3, size=100e-4,
            gas=None):
        self.dustfrac   = dustfrac
        self.chonfrac   = 0.75
        self.rho        = gas.rho * dustfrac
        self.nspecs     = nspecs
        self.mdust      = 3.3
        self.destroyed  = False
        if nspecs == 1:
            self.nspecs = nspecs
            self.size = [size]
            self.mass = [self.mdust*((4./3.)*np.pi*size*size*size)]
            self.numden = [self.rho*self.chonfrac/a for a in self.mass]
            self.temp = [10.0]
            self.vel = [1e5]
        elif nspecs == 0:
            self.nspecs = nspecs
            self.size = [size]
            self.mass = [mdust*((4.0/3.)*np.pi*size*size*size)]
            self.numden = [1e-30]
            self.temp = [10.0]
            self.vel = [1e5]
        else:
            """
            Multiple dust species! 
            Not yet supported
            """
            print 
            print 'Multiple dust species not supported'
            print 'Use ndust = 1'
            print 
            raise SystemExit
    """"""
    def _sumRho(self):
        return sum([a*b for (a,b) in zip(self.numden, self.mass)])
    """"""
    def _sumEnergy(self):
        return sum([a*b*c for (a,b,c) in 
            zip(self.gamma,self.numden,self.mass)])
    """"""
    def _calculateFdrag(self, w=None, idust = None, gas=None, rhoGas=None, NGas=None):
        """
        Calculate the drag forces
        This depends on 
        vdust, vgas, Tgas, Tdust, gas density
        """
        #
        # Get the proper dust properties from idust
        #
        nd = w[2+gas.nspecs + idust*4 + 0]
        vd = w[2+gas.nspecs + idust*4 + 1]
        td = w[2+gas.nspecs + idust*4 + 2]
        ad = w[2+gas.nspecs + idust*4 + 3]
        #
        # veldif
        #
        vg = w[0]
        vdg = (vd - vg)
        #
        # The s variable
        #
        from scipy.special import erf
        Fdrag = 0.0
        mbar = rhoGas/NGas
        s = np.abs(vdg)/np.sqrt(2.*kk*w[1]/mbar)
        #
        # Calculate the drag coefficients
        #
        if s == 0.0:
            return 0.0
        Cd = ( (2./(3.*s))*np.sqrt(np.pi*td/w[1]) + (
            (2.0*s*s+1.)/(np.sqrt(np.pi)*s*s*s) * np.exp(-s*s)) + (
            (4.0*s**(4.) + 4.*s*s - 1.)/(2.*s**(4.)) * erf(s)) )
        #
        # Calculate the drag
        #
        Fdrag1 = (-np.pi*ad*ad*rhoGas * Cd/2. * (np.abs(vdg)*vdg) )
        if np.isnan(Cd) or np.isnan(Fdrag1):
            Fdrag += 0.0
        else:
            Fdrag += Fdrag1
        """"""
        return Fdrag
    """"""
    def _calcqdust(self, Tg=None, vdg=None, Td=None, mass=None, gamma=5.0/3.0):
        """
        Calculate the dust heating rate
        qd = rhogas CHd ( Trec - Td) |vg - vd|
        """
        from scipy.special import erfinv, erf
        s = np.abs(vdg)/np.sqrt(2.*kk*Tg/mass)
        #
        # Calculate the heating rate
        #
        if np.abs(vdg) <= 1e0:
            qd = (gamma+1.)/(gamma-1.) * (Tg - Td) * (np.sqrt(
                (kk*kk*kk*Tg)/(8.*np.pi*mass*mass*mass)))
        else:
            #
            # Calculate the recovery temperature
            #
            fact = ((2.*gamma)/(gamma-1.) + 2.*s*s - ( 1./(
                0.5 + s*s + (s/np.sqrt(np.pi)) * np.exp(-s*s) * (1./erf(s)))) )
            Trec = Tg*(gamma-1.)/(gamma+1.) * fact
            if s == 0.0:
                Trec = Tg
            #
            # Calculate the heat transfer coefficient
            #
            CHd = (gamma+1.)/(gamma-1.) * (kk/(8.0*mass*s*s)) * (
                (s/np.sqrt(np.pi))*np.exp(-s*s) + ((0.5 + s*s) * 
                erf(s)) )
            """"""
            #
            # Calculate the heating rate
            #
            qd = CHd*(Trec - Td)*np.abs(vdg)
        return qd
    """"""
    def _fevap(self, Td=None):
        """
        Return the evaporation fraction
        """
        dT = 1e2
        if Td < 1950.0:
            return 0.0
        elif Td > 2050.0:
            return 1.0
        else:
            return ( 3./(dT*dT) * (Td - 2e3 + dT/2.)**(2.) -
                2./(dT*dT*dT) * (Td - 2e3 + dT/2.)**(3.))
    """"""
    def _calcDxa(self, w=None, idust = None, rhoGas=None, NGas=None, Jrad=None, gas=None, debug=False):
        """
        Calculate the change of the dust size
        dx ad = - f/(rho * H * vd) * (qd + radiation)
        This also depends on the recovery temperature 
        and the heat transfer coefficient
        -> needs the differential velocities
        """
        #
        # Get the proper dust properties from idust
        #
        nd = w[2+gas.nspecs + idust*4 + 0]
        vd = w[2+gas.nspecs + idust*4 + 1]
        td = w[2+gas.nspecs + idust*4 + 2]
        ad = w[2+gas.nspecs + idust*4 + 3]
        Hevap = 1.1e11
        #
        # veldif
        #
        vg = w[0]
        vdg = (vd - vg)
        #
        # The s variable
        #
        qd = 0.0
        mbar = rhoGas/NGas
        #
        # Calculate the heating rate
        # The function returns the middle part of the equation
        # qd = rhogas * qd * |vg - vd|
        #
        gam = sum([a * (b * 2. - 2.) for (a,b) in
            zip(w[2:gas.nspecs+2], gas.gamma)])
        gam = (gam/sum(w[2:gas.nspecs+2]) + 2.)/(gam/sum(w[2:gas.nspecs+2]))
        tempqd = self._calcqdust(Tg=w[1], vdg=vdg, Td=td, gamma=gam,
            mass=mbar)
        if np.isnan(tempqd):
            qd += 0.0
        else:
            qd += tempqd*rhoGas
        """"""
        #
        # Calculate the evaporation fraction
        #
        fevap = self._fevap(Td=td)
        #
        # Calculate the dxa
        #
        if Jrad is None:
            netheat = qd
        else:
            netheat = (qd + self._getEps(size=ad)*( np.pi*Jrad - ss*
                np.power(td, 4.)) )
        dxa = -fevap/(self.mdust*Hevap*vd) * netheat
        if netheat < 0.0:
            dxa = 0.0
        if debug:
            print '%2.3e %2.3e %2.3e'%(w[-3], w[-2], w[-1])
            print '%2.3e %2.3e %2.3e %2.3e'%(fevap, qd, netheat, dxa)
            raise SystemExit
        return dxa
    """"""
    def _getEps(self, size=None):
        """
        Return the absorption efficiencies
        """
        return 0.8*np.minimum(1.0, size/(2.0e-4))
    """"""
    def _calcDxTd(self, w=None, idust=None, rhoGas=None, NGas=None, gas=None,  Jrad=None):
        """
        Calculate the rate of change of the dust temperature
        """
        #
        # Get the proper dust properties from idust
        #
        nd = w[2+gas.nspecs + idust*4 + 0]
        vd = w[2+gas.nspecs + idust*4 + 1]
        td = w[2+gas.nspecs + idust*4 + 2]
        ad = w[2+gas.nspecs + idust*4 + 3]
        #
        # veldif
        #
        vg  = w[0]
        vdg = (vd - vg)
        #
        # Temperatures
        #
        Tg = w[1]
        Td = td
        #
        # The s variable
        # This has to be done per gas species
        #
        qd      = 0.0
        mbar    = rhoGas/ NGas
        s       = np.abs(vdg)/np.sqrt(2.*kk*w[1]/mbar)
        #
        # Calculate the heating rate
        # The function returns the middle part of the equation
        # qd = rhogas * qd * |vg - vd|
        #
        gam = sum([a * (b * 2. - 2.) for (a,b) in
            zip(w[2:gas.nspecs+2], gas.gamma)])
        gam = (gam/sum(w[2:gas.nspecs+2]) + 2.)/(gam/sum(w[2:gas.nspecs+2]))
        tempqd = self._calcqdust(Tg=w[1], vdg=vdg, Td=td, gamma=gam,
            mass=mbar)
        if np.isnan(tempqd):
            qd += 0.0
        else:
            qd += tempqd*rhoGas
        """"""
        #
        # Calculate the heat capacity
        #
        Cpd = self._getCpd(td = Td)
        #
        # Calculate the evaporation fraction
        #
        fevap = self._fevap(Td=Td)
        #
        # Calculate the rate of change
        #
        if Jrad is None:
            netheat = qd
        else:
            netheat = qd + self._getEps(size=ad)*( np.pi*Jrad -
                ss*np.power(Td, 4.))
        dxtd = ( (3. * (1.0 - fevap))/( vd * Cpd * self.mdust * ad) * netheat)
        return dxtd
    """"""
    def _getCpd(self, td = None):
        """ 
        Calculate and return the heat capacity
        """
        #
        # Calculate the heat capacity
        #
        if (td > 1400.0) and (td < 1820.0):
            Cpd = 1e7 + 5e9/(1820.-1400.0)
        else:
            Cpd = 1e7
        return Cpd
    """"""
""""""

