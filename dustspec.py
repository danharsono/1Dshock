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
        if nspecs == 1:
            self.nspecs = nspecs
            self.size = [size]
            self.mass = [self.mdust*((4./3.)*np.pi*size*size*size)]
            self.numden = [self.rho/a for a in self.mass]
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
        if size < 1e-5:
            size = 1e-5
        self.size = [size]
        self.mass = [self.mdust*((4.0/3.)*np.pi*size*
            size*size)]
        self.vel = [vdust]
        self.temp=[tdust]
    """"""
    def _calculateFdrag(self, vd=None, vg=None, Tg=None, Td=None, gas=None):
        """
        Calculate the drag forces
        This depends on 
        vdust, vgas, Tgas, Tdust, gas density
        """
        #
        # veldif
        #
        vdg = (vd - vg)
        if (Tg < 0.0) or (Td < 0):
            print vd, vg, Tg, Td
            print 'Negative temperatures! in Drag'
            raise SystemExit
        """"""
        #
        # The s variable
        #
        from scipy.special import erf
        Fdrag = 0.0
        for ispec in xrange(gas.nspecs):
            s = np.abs(vdg)/np.sqrt(2.*kk*Tg/(gas.mass[ispec]))
            #
            # Calculate the drag coefficients
            #
            gam = gas.gamma[ispec]
            mass = gas.mass[ispec]
            if s > 1e3:
                Cd = 2.
            elif np.abs(vdg) < 1e-3:
                if (vdg == 0.0):
                    Cd = 0.0 # If it is 0.0 then set this to 0
                else:
                    Cd = (2.0/(3.*s))*np.sqrt(np.pi*Td/Tg)
            else:
                Cd = ( (2./(3.*s))*np.sqrt(np.pi*Td/Tg) + (
                    (2.0*s*s+1.)/(np.sqrt(np.pi)*s*s*s) * np.exp(-s*s)) + (
                    (4.0*s**(4.) + 4.*s*s - 1)/(2.*s**(4.)) * erf(s)) )
            #
            # Calculate the drag
            #
            Fdrag += (-np.pi*self.size[0]**(2.)*gas.numden[ispec]*
                 gas.mass[ispec]*Cd/2. *(np.abs(vdg)*vdg) )
            if np.isnan(Cd) or np.isnan(Fdrag):
                print ispec, vdg, s, Cd, Fdrag
                print 'NAN is found in Drag force!'
                raise SystemExit
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
        if Td < 1950.0:
            return 0.0
        elif Td > 2050.0:
            return 1.0
        else:
            return (Td - 1950.0)/100.0
    """"""
    def _calcDxa(self, vd=None, vg=None, Tg=None, Td=None, gas=None, Jrad=None):
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
        vdg = (vd - vg)
        #
        # The s variable
        #
        qd = 0.0
        for ispec in xrange(gas.nspecs):
            s = np.abs(vdg)/np.sqrt(2.*kk*Tg/(gas.mass[ispec]))
            #
            # Calculate the heating rate
            # The function returns the middle part of the equation
            # qd = rhogas * qd * |vg - vd|
            #
            gam = gas.gamma[ispec]
            mass = gas.mass[ispec]
            if np.abs(s) > 5e1:
                tempqd = (1./8.)*gas.mass[ispec] * gas.numden[ispec]*np.abs(
                    vdg)**(3.)
            elif np.abs(vdg) < 1e-3:
                tempqd = (gam+1.)/(gam-1.) * (Tg - Td) * (np.sqrt(
                    (kk*kk*kk*Tg)/(8.*np.pi*mass*mass*mass)) *
                    gas.numden[ispec]*mass)
            else:
                tempqd = ( self._calcqdust(Tg=Tg, s=s, Td=Td,
                    gamma=gam, mass=mass) *
                    gas.numden[ispec]*gas.mass[ispec] * np.abs(vg - vd))
            qd += tempqd
        """"""
        #
        # Calculate the evaporation fraction
        #
        fevap = self._fevap(Td=Td)
        #
        # Calculate the dxa
        #
        if Jrad is None:
            Jrad = ss*Td**(4.)/(np.pi)
        netheat = (qd + self._getEps()*( np.pi*Jrad - ss*Td**(4.)) )
        if netheat < 0.0:
            fevap = 0.0
        dxa = -fevap/(self.mdust*Hevap*vd) * (qd +
            self._getEps()*( np.pi*Jrad - ss*Td**(4.)) )
        return dxa
    """"""
    def _getEps(self):
        """
        Return the absorption efficiencies
        """
        return 0.8*np.minimum(1.0, self.size[0]/(2.0e-4))
    """"""
    def _calcDxTd(self, vd=None, vg=None, Tg=None, Td=None, gas=None,  Jrad=None):
        """
        Calculate the rate of change of the dust temperature
        """
        #
        # veldif
        #
        vdg = (vd - vg)
        #
        # The s variable
        # This has to be done per gas species
        #
        qd = 0.0
        for ispec in xrange(gas.nspecs):
            s = np.abs(vdg)/np.sqrt(2.*kk*Tg/(gas.mass[ispec]))
            #
            # Calculate the heating rate
            # The function returns the middle part of the equation
            # qd = rhogas * qd * |vg - vd|
            #
            gam = gas.gamma[ispec]
            mass = gas.mass[ispec]
            if np.abs(vdg) > 1e4:
                tempqd = (1./8.)*gas.mass[ispec] * gas.numden[ispec]*np.abs(
                  vdg)**(3.)
            if vdg < 1e-3:
                tempqd = (gam+1.)/(gam-1.) * (Tg - Td) * (np.sqrt(
                    (kk*kk*kk*Tg)/(8.*np.pi*mass*mass*mass)) *
                    gas.numden[ispec]*mass)
            else:
                tempqd = ( self._calcqdust(Tg=Tg, s=s, Td=Td,
                    gamma=gam, mass=mass) *
                    gas.numden[ispec]*gas.mass[ispec] * np.abs(vg - vd))
            qd += tempqd
        """"""
        #
        # Calculate the heat capacity
        #
        if (Td > 1400.0) and (Td < 1820.0):
            Cpd = 2.19e7
        else:
            Cpd = 1e7
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
            self.mdust * self.size[0]) * (qd +
            self._getEps()*( np.pi*Jrad -
            ss*Td**(4.)) ) )
        return dxtd
    """"""
""""""
