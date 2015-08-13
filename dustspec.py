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
        self.gas = gas # gas properties
        self.rho = self.gas.rho * dustfrac
        self.mdust = 3.3
        self.destroyed = False
        if nspecs == 1:
            self.nspecs = nspecs
            self.size = [size]
            self.mass = [self.mdust*((4./3.)*np.pi*size*size*size)]
            self.numden = [self.rho*0.75/a for a in self.mass]
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
    def _updateRho(self, rho=None):
        """
        Update rho
        """
        self.rho = rho
        self.mass = [self.mdust*((4.0/3.)*np.pi*self.size*
            self.size*self.size)]
        self.numden=[self.rho/a for a in self.mass]

    """"""
    def _updateDust(self, allns=None, size=None, tdust=None, vdust=None):
        """
        Update dust densities
        """
        self.numden = allns
        if size < 1e-20:
            size = 1e-20
        self.size = [size]
        self.mass = [self.mdust*((4.0/3.)*np.pi*size*
            size*size)]
        self.vel = [vdust]
        self.temp=[tdust]
    """"""
    def _calculateFdrag(self, w=None, rhoGas=None, NGas=None):
        """
        Calculate the drag forces
        This depends on 
        vdust, vgas, Tgas, Tdust, gas density
        """
        #
        # veldif
        #
        vd = w[-3]
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
        Cd = ( (2./(3.*s))*np.sqrt(np.pi*w[-2]/w[1]) + (
            (2.0*s*s+1.)/(np.sqrt(np.pi)*s*s*s) * np.exp(-s*s)) + (
            (4.0*s**(4.) + 4.*s*s - 1)/(2.*s**(4.)) * erf(s)) )
        #
        # Calculate the drag
        #
        Fdrag1 = (-np.pi*w[-1]**(2.)*rhoGas * Cd/2. * (np.abs(vdg)*vdg) )
        if np.isnan(Cd) or np.isnan(Fdrag1):
            Fdrag += 0.0
        else:
            Fdrag += Fdrag1
        """"""
        return Fdrag
    """"""
    def _calcqdust(self, Tg=None, s=None, Td=None, mass=None, gamma=5.0/3.0):
        """
        Calculate the dust heating rate
        qd = rhogas CHd ( Trec - Td) |vg - vd|
        """
        from scipy.special import erfinv, erf
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
        qd = CHd*(Trec - Td)
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
    def _calcDxa(self, w=None, rhoGas=None, NGas=None, Jrad=None, gas=None, debug=False):
        """
        Calculate the change of the dust size
        dx ad = - f/(rho * H * vd) * (qd + radiation)
        This also depends on the recovery temperature 
        and the heat transfer coefficient
        -> needs the differential velocities
        """
        Hevap = 1.1e11
        #
        # veldif
        #
        vd = w[-3]
        vg = w[0]
        vdg = (vd - vg)
        #
        # The s variable
        #
        qd = 0.0
        mbar = rhoGas/NGas
        s = np.abs(vdg)/np.sqrt(2.*kk*w[1]/mbar)
        #
        # Calculate the heating rate
        # The function returns the middle part of the equation
        # qd = rhogas * qd * |vg - vd|
        #
        gam = sum([a * (b * 2. - 2.) for (a,b) in
            zip(w[2:gas.nspecs+2], gas.gamma)])
        gam = (gam/sum(w[2:gas.nspecs+2]) + 2.)/(gam/sum(w[2:gas.nspecs+2]))
        if np.abs(vdg) <= 1e2:
            tempqd = (gam+1.)/(gam-1.) * (w[1] - w[-2]) * (np.sqrt(
                (kk*kk*kk*w[1])/(8.*np.pi*mbar*mbar*mbar)))
        else:
            tempqd = ( self._calcqdust(Tg=w[1], s=s, Td=w[-2],
                gamma=gam, mass=mbar) * np.abs(vg - vd))
        if np.isnan(tempqd):
            qd += 0.0
        else:
            qd += tempqd*rhoGas
        """"""
        #
        # Calculate the evaporation fraction
        #
        fevap = self._fevap(Td=w[-2])
        #
        # Calculate the dxa
        #
        if Jrad is None:
            Jrad = ss*np.power(w[-2],4.)/(np.pi)
        netheat = (qd + self._getEps(size=w[-1])*( np.pi*Jrad - ss*
            np.power(w[-2], 4.)) )
        if netheat < 0.0:
            fevap = 0.0
        dxa = -fevap/(self.mdust*Hevap*w[-3]) * (qd +
            self._getEps(size=w[-1])*( np.pi*Jrad - ss*np.power(w[-2],4.) ))
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
    def _calcDxTd(self, w=None, rhoGas=None, NGas=None, gas=None,  Jrad=None):
        """
        Calculate the rate of change of the dust temperature
        """
        #
        # veldif
        #
        vd  = w[-3]
        vg  = w[0]
        vdg = (vd - vg)
        #
        # Temperatures
        #
        Tg = w[1]
        Td = w[-2]
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
        if np.abs(vdg) <= 1e2:
            tempqd = (gam+1.)/(gam-1.) * (w[1] - w[-2]) * (np.sqrt(
                (kk*kk*kk*w[1])/(8.*np.pi*mbar*mbar*mbar)))
        else:
            tempqd = ( self._calcqdust(Tg=w[1], s=s, Td=w[-2],
                gamma=gam, mass=mbar) * np.abs(vg - vd))
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
        if qd < 0.0:
            fevap = 0.0
        #
        # Calculate the rate of change
        #
        if Jrad is None:
            Jrad = ss*Td**(4.)/(np.pi)
        dxtd = ( (3. * (1.0 - fevap))/( vd * Cpd *
            self.mdust * w[-1]) * (qd +
            self._getEps(size=w[-1])*( np.pi*Jrad - ss*np.power(Td, 4.))) )
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

