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
    def __init__(self, dustfrac = 0.01, nspecs=1, mdust=3.3, size=100e-4,  gas=None):
        self.gas = gas # gas properties
        self.rho = self.gas.rho * dustfrac
        self.mdust = 3.3
        if nspecs == 1:
            self.nspecs = nspecs
            self.size = [size]
            self.mass = [mdust*((4.0/3.)*np.pi*size*size*size)]
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
        n0 = self.rho/sum(self.specfrac)
        self.numden=[a*n0/b for (a,b) in zip(self.specfrac,
            self.mass)]
        self.mass = [self.mdust*((4.0/3.)*np.pi*self.size*
            self.size*self.size)]
    """"""
    def _updateDust(self, allns=None, size=None):
        """
        Update dust densities
        """
        self.numden = allns
        self.size = [size]
        self.mass = [self.mdust*((4.0/3.)*np.pi*size*
            size*size)]
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
        if vdg == 0.0:
            return 0.0
        #
        # The s variable
        #
        from scipy.special import erf
        Fdrag = 0.0
        for ispec in xrange(gas.nspecs):
            s = np.abs(vdg)/np.sqrt(2.*kk*Tg/(gas.mass[ispec]))
            #
            # Calculate the drag coefficient
            #
            Cd = ( (2./(3.*s))*np.sqrt(np.pi*Td/Tg) + (
                (2.0*s*s+1.)/(np.sqrt(np.pi)*s*s*s) * np.exp(-s*s)) + (
                (4.0*s**(4.) + 4.*s*s - 1)/(2.*s**(4.)) * erf(s)) )
            #
            # Calculate the drag
            #
            Fdrag = (-np.pi*self.size[0]**(2.)*gas.numden[ispec]*gas.mass[ispec]
                *Cd/2. *(np.abs(vdg)*vdg) )
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
        if s < 1e-25:
            Trec = Tg
            CHd = ((gamma+1.0)/(gamma-1.)) * (kk/(4.0*np.sqrt(np.pi)*
                mass*s))
        elif s > 1e2:
            Trec = 2.0*(gamma-1.0)/(gamma+1) * s*s*Tg
            CHd = ((gamma+1.)/(gamma-1.))*(kk/(8.0*mass))
        else:
            bottom = ((2.*gamma)/(gamma-1.) * 2.*s*s - 0.5 + 
                (2./np.sqrt(np.pi)) * np.exp(-s*s) * erfinv(s))
            Trec = (Tg*(gamma-1.)/(gamma+1.) / bottom)
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
        if (vdg ==0) and (Tg == Td):
            return 0.0
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
            if s < 1e-10:
                tempqd = ( gas.numden[ispec]*0.25*np.sqrt(
                    (8*kk*Tg)/(np.pi*gas.mass[ispec])) * 0.5 * ( (
                    gas.gamma[ispec]+1.0)/(gas.gamma[ispec]-1.0) ) * 
                    kk * (Tg-Td))
            elif s > 10.:
                tempqd = (1.0/8.0)*gas.numden[ispec]*gas.mass[ispec]*(
                    np.abs(vdg)**(3.) )
            else:
                tempqd = self._calcqdust(Tg=Tg, s=s, Td=Td, 
                    gamma=gas.gamma[ispec], mass=gas.mass[ispec])
                tempqd = gas.numden[ispec]*gas.mass[ispec] * (
                    tempqd * np.abs(vg - vd) )
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
            Jrad = 0.0
        dxa = -fevap/(self._sumRho()*Hevap*vd) * (qd +
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
        if (vdg ==0) and (Tg == Td):
            return 0.0
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
            if s < 1e-25:
                tempqd = ( gas.numden[ispec]*0.25*np.sqrt(
                    (8*kk*Tg)/(np.pi*gas.mass[ispec])) * 0.5 * ( (
                    gas.gamma[ispec]+1.0)/(gas.gamma[ispec]-1.0) ) * 
                    kk * (Tg-Td))
            elif s > 1e0:
                tempqd = (1.0/8.0)*gas.numden[ispec]*gas.mass[ispec]*(
                    np.abs(vdg)**(3.) )
            else:
                tempqd = self._calcqdust(Tg=Tg, s=s, Td=Td, 
                    gamma=gas.gamma[ispec], mass=gas.mass[ispec])
                tempqd = gas.numden[ispec]*gas.mass[ispec] * (
                    tempqd * np.abs(vg - vd) )
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
        #
        # Calculate the rate of change
        #
        if Jrad is None:
            Jrad = 0.0
        dxtd = ( (3. * (1.0 - fevap))/( vd * Cpd * 
            self._sumRho() * self.size[0]) * (qd + 
            self._getEps()*( np.pi*Jrad -
            ss*Td**(4.)) ) )
        if dxtd < 0.0: dxtd = 0
        return dxtd
    """"""
""""""
