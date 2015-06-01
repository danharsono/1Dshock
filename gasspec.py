from python.natconst import *
import numpy as np

"""
The gas components for the shock code
"""
class gasSpecs():
    """ 
    Constructor
    Destructor is not required in python
    """
    def __init__(self, rho = 1e-9, vgas=6.5e5, tgas = 300., nspecs=3):
        self.rho = rho
        self.vgas = vgas
        self.tgas = tgas
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
        elif nspecs == 3:
            """
            Use H, H2 and He
            """
            self.nspecs =3
            self.mass = [mp, 2.0*mp, 4.0*mp]
            self.specfrac = [0.01, 0.74, 0.25]
            mu = 0.01 + 0.74/2.0 + 0.25/4.0
            self.mugas = mu**(-1.0)
            n0 = self.rho/sum(self.specfrac)
            self.numden = [a*n0/b for (a,b) in zip(self.specfrac,
                self.mass)]
            self.gamma = [5.0/2.0, 7.0/2.0, 5.0/2.0]
        elif nspecs == 4:
            """
            Add SiO
            """
            self.nspecs=4
            self.mass = [mp, 2.0*mp, 4.0*mp, 44.0*mp]
            self.specfrac = [0.01, 0.74, 0.25, 0.0]
            mu = 0.01 + 0.65/2.0 + 0.34/4.0
            self.mugas = mu**(-1.0)
            n0 = self.rho/sum(self.specfrac)
            self.numden = [a*n0/b for (a,b)
                in zip(self.specfrac, self.mass)]
            self.gamma = [5.0/2.0, 7.0/2.0, 5.0/2.0, 7.0/2.]

    """"""
    def _sumRho(self):
        return sum([a*b for (a,b) in zip(self.numden, self.mass)])
    """"""
    def _getMugas(self):
        curfrac = [(a*b)/self._sumRho() for (a,b) in
               zip(self.numden, self.mass)]
        mu = [a/(b/mp) for (a,b) in zip(curfrac, self.mass)]
        return sum(mu)**(-1.)
    """"""
    def _sumEnergy(self):
        return sum([a*b*c for (a,b,c) in 
            zip(self.gamma,self.numden, self.mass)])
    """"""
    def _updateRho(self, rho=None):
        """
        Update rho
        """
        self.rho = rho
        n0 = self.rho/sum(self.specfrac)
        self.numden=[a*n0/b for (a,b) in zip(self.specfrac,
                self.mass)]
    """"""
    def _getGamma(self):
        """
            Return the total specific heats
        """
        top = sum([a*b for (a,b) in zip(self.gamma, self.numden)])
        bottom = sum([(a-1.)*b for (a,b) in zip(self.gamma, self.numden)])
        return top/bottom
    """"""
    def _updateGas(self, allns=None, tgas=None, vgas=None):
        """
        Update gas densities
        """
        self.numden = allns
        self.tgas = tgas
        self.vgas = vgas
    """"""
    def _calculateR(self, t=None):
        """
        Calculate the reactions and update the values
        This rate is the creation of gas species i 
        in terms of number of reactions per unit time 
        per unit volume
        """
        rate1 = 8.72e-33*(t/300.0)**(-0.6) # cm^6 per sec
        rate2 = 1.83e-31*(t/300.0)**(-0.6)
        rate3 = 1.50e-9*np.exp(-46350./t)**(1.0)
        rate4 = 3.75e-8*(t/300.0)**(-0.6)*np.exp(-53280./t)
        if self.nspecs > 1:
            totrate = (self.numden[0]**(2.)*(rate1*self.numden[1] + 
                rate2*self.numden[0]) - self.numden[1]*(rate3*
                self.numden[1]+rate4*self.numden[0]) )
            return totrate
        else:
            return 0.0
    """"""
    def _calculateVarC(self, vel=None, t=None, veld=None, td=None, dust=None, Jrad=None):
        """
        Calculate the right hand side of the momentum equation
        """
        normrate = self._calculateR(t=t)
        if self.nspecs > 1:
            first = 2.0*normrate*(vel*self.mass[0]+(kk*t/vel))
            second = -normrate*(vel*self.mass[1]+(kk*t/vel))
            third = 0.0
            #
            # This is now dust stuff
            #
            if dust is not None:
                if Jrad is None:
                    dxa = dust._calcDxa(vd=veld, vg=vel, Tg=t, Td=td, gas=self)
                    fdrag = dust._calculateFdrag(vd=veld,
                         vg=vel, Tg=t, Td=td, gas=self)
                else:
                    dxa = dust._calcDxa(vd=veld, vg=vel, Tg=t,
                        Td=td,rhogas=self._sumRho(), Jrad=Jrad)
                    fdrag = dust._calculateFdrag(vd=veld,
                         vg=vel, Tg=t, Td=td, gas=self, Jrad = Jrad)
                """"""
                fourth = ( (dust.numden[0]*dust.vel[0]*
                    4.0*np.pi*dust._sumRho()*dust.size[0]*dust.size[0])/
                    (self.mass[3]) )*dxa
                fifth = - (dust.numden[0]* fdrag +
                    4.0*np.pi*dust.size[0]*dust.size[0]*
                    dust._sumRho() * veld*veld*dxa )
                if np.abs(fourth) > 1e2:
                    print '%2.5e %2.5e %2.5e %2.5e'%(vel, t, veld, td)
                    print '%2.5e %2.5e %2.5e %2.5e'%(dust.numden[0],
                         dust.vel[0], dust._sumRho(), dust.size[0])
                    print '%2.5e %2.5e %2.5e %2.5e %2.5e'%(first,
                        second, third, fourth, fifth)
                    top = (dust.numden[0]*dust.vel[0]*
                           4.0*np.pi*dust._sumRho()*dust.size[0]*dust.size[0])
                    print '%2.5e %2.5e'%(dxa, self.mass[3])
                    print
            else:
                fourth = 0.0
                fifth = 0.0
            return first+second+third+fourth+fifth
        else:
            return 0.0
        
    """"""
    def _calculateVarF(self, vel=None, t=None, veld=None, td=None, dust=None, Jrad=None):
        """
        Calculate the right hand side of the energy equation
        """
        normrate = self._calculateR(t=t)
        if self.nspecs > 1:
            first = -2.0*normrate*(self.mass[0]*self.gamma[0]*kk*t +
                0.5*vel*vel*self.mass[0])
            second = normrate*(self.mass[1]*self.gamma[1]*kk*t +
                    0.5*vel*vel*self.mass[1])
            third = 0.0
            #
            # Dust stuffs below
            #
            if dust is not None:
                if Jrad is None:
                    dxa = dust._calcDxa(vd=veld, vg=vel, Tg=t, Td=td, gas=self)
                    fdrag = dust._calculateFdrag(vd=veld,
                         vg=vel, Tg=t, Td=td, gas=self)
                    dxtd = dust._calcDxTd(vd=veld, vg=vel, Tg=t,
                        Td=td, gas=self)
                else:
                    dxa = dust._calcDxa(vd=veld, vg=vel, Tg=t,
                        Td=td,rhogas=self._sumRho(), Jrad=Jrad)
                    fdrag = dust._calculateFdrag(vd=veld,
                        vg=vel, Tg=t, Td=td, gas=self, Jrad = Jrad)
                    dxtd = dust._calcDxTd(vd=veld, vg=vel, Tg=t,
                      Td=td, gas=self, Jrad = Jrad)
                """"""
                fourth = ((dust.numden[0]*dust.vel[0]*4.0*np.pi*
                    dust._sumRho()*dust.size[0]*dust.size[0]*dxa)*(
                    self.mass[3]*self.gamma[3]*kk*t +
                    0.5*vel*vel*self.mass[3]))
                fifth = veld*dust.numden[0]*(3./2. * fdrag +
                    2.0*np.pi*dust.size[0]*dust.size[0]*
                    dust._sumRho()*veld*veld*dxa)
                dustmass = (4.0/3.0)*np.pi*dust._sumRho()* dust.size[0]**(3.)
                #
                # Calculate the heat capacity
                #
                if (td > 1400.0) and (td < 1820.0):
                    Cpd = 2.19e7
                else:
                    Cpd = 1e7
                """"""
                sixth = (dust.numden[0]*dustmass*veld*Cpd*dxtd)
                seventh = (4.0*np.pi*dust.size[0]**(2.)*dust.numden[0]*
                   Cpd*veld*td*dust._sumRho()*dxa )
            else:
                fourth = 0.0
                fifth = 0.0
                sixth = 0.0
                seventh= 0.0
            return first+second+third+fourth+fifth+sixth+seventh
        else:
            return 0.0
        
    """"""
    def _calculateFreeEnergy(self, t=None):
        """
        Calculate the net energy
        H + H + energy -> H2 
        H2 + energy -> H + H
        """
        normrate = self._calculateR(t=t)
        onev = 1.6021772e-12
        return 2*normrate*4.48*onev - normrate*4.48*onev
    """"""
    
""""""