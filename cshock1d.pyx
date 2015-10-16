"""
Cython implementation for the 1D shock code
"""
from __future__ import division
from scipy.special import erf, erfinv
from numpy cimport ndarray, dtype
import numpy as np; cimport numpy as np
cimport cython; from cpython cimport bool
#
# Types
#
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t
#
# Constants here and externs from math
#
cdef extern from "math.h":
    double sqrt(double m)
    double pow(double m, double n)
    double exp(double a)
    double fabs(double a)
    double M_PI
    bint isnan(double x)
    double fmin(double, double)
cdef double kk  = 1.3806488e-16
cdef double ss  = 5.670367e-5
#
# Functions HERE
#
"""
Decipher the input parameters from the gas and dust parameters
"""
cdef decipherW(ndarray[DTYPE_t, ndim=1] w, int nspecs, int ndust):
    """
        Variables
    """
    cdef int ig, id
    cdef ndarray[DTYPE_t, ndim=1] w1 = np.zeros(w.shape[0], dtype=DTYPE)
    #
    # First copy
    #
    for ig in range(w.shape[0]):
        w1[ig]  = w[ig]
    #
    # modify gas
    #
    for ig in range(nspecs):
        w1[2+ig]    = w[2+ig]
    #
    # modify dust: limit size
    #
    for id in range(ndust):
        if (w1[<unsigned int> (2+nspecs+4*id+3)] < 0.0):
            w1[2+nspecs+4*id+3] = 1e-30
    return w1
""""""
cdef double calcQdust(double Tg, double vdg, double Td, double mass, double gamma):
    """
    Calculate the dust heating rate
    qd = rhogas CHd ( Trec - Td) |vg - vd|
    """
    #
    # Variables
    #
    cdef double s, qd, fact, Trec, CHd
    #
    #
    #
    s = fabs(vdg)/sqrt(2.*kk*Tg/mass)
    #
    # Calculate the heating rate
    #
    qd = 0.0
    if fabs(vdg) <= 1e0:
        qd += (gamma+1.)/(gamma-1.) * (Tg - Td) * (sqrt(
            (pow(kk,3.)*Tg)/(8.*M_PI*pow(mass,3.))))
    else:
        """
        #
        # Calculate the recovery temperature
        #
        """
        fact = ((2.*gamma)/(gamma-1.) + 2.*s*s - ( 1./(
            0.5 + s*s + (s/sqrt(M_PI)) * exp(-s*s) * (1./erf(s)))) )
        Trec = Tg*(gamma-1.)/(gamma+1.) * fact
        if s == 0.0:
            Trec = Tg
        #
        # Calculate the heat transfer coefficient
        #
        CHd = (gamma+1.)/(gamma-1.) * (kk/(8.0*mass*s*s)) * (
            (s/sqrt(M_PI))*exp(-s*s) + ((0.5 + s*s) *
            erf(s)) )
        """"""
        #
        # Calculate the heating rate
        #
        qd += CHd*(Trec - Td)*fabs(vdg)
    return qd
""""""
cdef double getCpd(double Td):
    """
    Calculate and return the heat capacity
    """
    #
    # Calculate the heat capacity
    #
    cdef double Cpd = 0.0
    if (Td > 1400.0) and (Td < 1820.0):
        Cpd += 1e7 + 5e9/(1820.-1400.0)
    else:
        Cpd += 1e7
    return Cpd
""""""
cdef double getfEvap(double Td):
    """
    Return the evaporation fraction
    """
    cdef double dT = 1e2
    if Td < 1950.0:
        return 0.0
    elif Td > 2050.0:
        return 1.0
    else:
        return ( 3./(dT*dT) * (Td - 2e3 + dT/2.)**(2.) -
            2./(dT*dT*dT) * (Td - 2e3 + dT/2.)**(3.))
""""""
cdef double getEPS(double ad):
    return 0.8*fmin(1.0, ad/(2.0e-4))
""""""
cdef double calcFDrag(double mbar, double nd, double vd, double td, double ad, double vg, double Tg, double Rhogas):
    """
    Calculate the drag forces
    This depends on 
    vdust, vgas, Tgas, Tdust, gas density
    """
    #
    # Variables
    #
    cdef double vdg, Fdrag, Cd, s
    #
    # veldif
    #
    vdg = (vd - vg)
    #
    # The s variable
    #
    Fdrag = 0.0
    s = fabs(vdg)/sqrt(2.* kk * Tg/ mbar)
    #
    # Calculate the drag coefficients
    #
    if s == 0.0:
        return 0.0
    else:
        Cd = ( (2./(3.*s)) * sqrt(M_PI*td/Tg) + (
            (2.0*s*s+1.)/(sqrt(M_PI)* s * s * s) * np.exp(-s*s)) + (
            (4.0* pow(s,4.) + 4.*s*s - 1.)/(2.* pow(s,4.)) * erf(s)) )
    #
    # Calculate the drag
    #
    Fdrag1 = (-M_PI*ad*ad*Rhogas * Cd/2. * (fabs(vdg)*vdg) )
    if isnan(Cd) or isnan(Fdrag1):
        Fdrag += 0.0
    else:
        Fdrag += Fdrag1
    """"""
    return Fdrag
""""""
cdef double calcDXa(double nd, double vd, double td, double ad, double vg, double gam, double Tg, double Jrad, double rhod, double Rhogas, double Ngas):
    """
        Calculate the change of the dust size
        dx ad = - f/(rho * H * vd) * (qd + radiation)
        This also depends on the recovery temperature 
        and the heat transfer coefficient
        -> needs the differential velocities
    """
    #
    # Variables
    #
    cdef double vdg, qd, mbar, fevap, netheat, dxa
    cdef double Hevap = 1.1e11
    #
    # veldif
    #
    vdg = (vd - vg)
    #
    # The s variable
    #
    qd      = 0.0
    mbar    = Rhogas/Ngas
    #
    # Calculate the heating rate
    # The function returns the middle part of the equation
    # qd = rhogas * qd * |vg - vd|
    #
    tempqd = calcQdust(Tg, vdg, td, mbar, gam)
    if isnan(tempqd):
        qd += 0.0
    else:
        qd += tempqd*Rhogas
    """"""
    #
    # Calculate the evaporation fraction
    #
    fevap = getfEvap(td)
    #
    # Calculate the dxa
    #
    netheat     = 0.0
    if Jrad == 0.0:
        netheat += qd
    else:
        netheat = (qd + getEPS(ad)*( M_PI*Jrad - ss*pow(td, 4.)) )
    dxa = -fevap/(rhod*Hevap*vd) * netheat
    if netheat < 0.0:
        dxa = 0.0
    return dxa
""""""
cdef double calcDxTd(double nd, double vd, double td, double ad, double vg, double Tg, double Rhogas, double Ngas, double gam, double Jrad, double rhod):
    """
    Calculate the rate of change of the dust temperature
    """
    #
    # variables
    #
    cdef double vdg, qd, mbar, s, tempqd, Cpd, fevap, netheat, dxtd
    #
    # veldif
    #
    vdg = (vd - vg)
    #
    # The s variable
    # This has to be done per gas species
    #
    qd      = 0.0
    mbar    = Rhogas/ Ngas
    s       = fabs(vdg)/sqrt(2.* kk * Tg / mbar)
    #
    # Calculate the heating rate
    # The function returns the middle part of the equation
    # qd = rhogas * qd * |vg - vd|
    #
    tempqd = calcQdust(Tg, vdg, td, mbar, gam)
    if isnan(tempqd):
        qd += 0.0
    else:
        qd += tempqd*Rhogas
    """"""
    #
    # Calculate the heat capacity
    #
    Cpd = getCpd(td)
    #
    # Calculate the evaporation fraction
    #
    fevap = getfEvap(td)
    #
    # Calculate the rate of change
    #
    netheat = 0.0
    if Jrad == 0.0:
        netheat += qd
    else:
        netheat += qd + getEPS(ad)*(M_PI*Jrad - ss*pow(td, 4.))
    """"""
    dxtd    = 0.0
    dxtd    += ( (3. * (1.0 - fevap))/( vd * Cpd * rhod * ad) * netheat)
    return dxtd
""""""
#
# Reaction rate
#
cdef double calcR(double Tg, double nH, double nH2):
    """
    Reactions
    """
    cdef double rate1, rate2, rate3, rate4
    cdef double totrate
    #
    rate1 = 8.72e-33*pow((Tg/300.0), -0.6) # cm^6 per sec
    rate2 = 1.83e-31*pow((Tg/300.0), -1) # cm^6 per sec
    rate3 = 1.50e-9*exp(-46350./Tg)
    if (-53280/Tg > -40.0):
        rate4   = 3.75e-8*pow(Tg/300.0, -0.5)*exp(-53280./Tg)
    else:
        rate4   = 0.0
    totrate = (nH * nH *(rate1*nH2 +
        rate2*nH) - nH2*(rate3*nH2+rate4*nH) )
    return totrate
""""""
#
# Chemistry of gas -> move to another file later along with some above
#
cdef calculateR(double Tg, int nspecs, double nH, double nH2):
    """
    Calculate the reactions and update the values
    This rate is the creation of gas species i 
    in terms of number of reactions per unit time 
    per unit volume
    """
    #
    # Variables
    #
    cdef ndarray[DTYPE_t, ndim=2] temprates = np.zeros(
        (nspecs, nspecs), dtype=DTYPE)
    #
    # Create the reaction matrix
    #
    if nspecs > 1:
        totrate = calcR(Tg, nH,nH2)
    else:
        totrate = 0.0
    temprates[0,0] = -2. * totrate
    temprates[0,1] = totrate
    return temprates
""""""
cdef double calculateFreeEnergy(double Tg, int nspecs, double nH, double nH2):
    """
    Calculate the net energy
    H + H + energy -> H2 
    H2 + energy -> H + H
    """
    #
    # Variables
    #
    cdef double onev = 1.6021772e-12
    #
    totrate = calcR(Tg, nH, nH2)
    return <double> (totrate*4.48*onev)
""""""
cdef double gasKap(double Tg, int destroyed, ndarray[DTYPE_t, ndim=1] Tkap,
    ndarray[DTYPE_t, ndim=1] Kaps):
    """
        Averaged gas opacities
    """
    cdef double kap
    kap     = <double> pow(10.0, np.interp(Tg, Tkap, Kaps))
    if (destroyed == 1 and (Tg < 1.4e3)):
        kap = 0.5
    return kap
""""""
"""
The vector field of the matrix to solve
"""
@cython.boundscheck(False)
@cython.cdivision
cdef vectorfield(double x, np.ndarray[DTYPE_t, ndim=1] w, np.ndarray[DTYPE_t, ndim=2] Jrad, int haveJ, int nspecs, int ndust, ndarray[DTYPE_t, ndim=1] gmasses, ndarray[DTYPE_t, ndim=1] gammas, double rhod, int destroyed, ndarray[DTYPE_t, ndim=1] Tkap, ndarray[DTYPE_t, ndim=1] Kaps, ndarray[DTYPE_t, ndim=1] gdustFrac):
    """
    w   :  the variables to solve -> x, y 
    x   :  position with respect to the shock
    p   :  parameters
    """
    #
    # Variables here
    #
    cdef int ii, idust, ig, dumid
    cdef double curJ
    cdef double Rhogas, Ngas, vgas, Tgas, gam, gam1
    cdef double dustLoss, nd, vd, td, ad
    cdef double varA, varB, varC, varD, varE, varF
    cdef double dumGas, dumDust, radiation
    cdef double f1, f2
    dumid   = nspecs+2
    cdef ndarray[DTYPE_t, ndim=1] YDOT = np.zeros(
        <unsigned int> (2 + nspecs+ 4*ndust), dtype=DTYPE)
    #
    # Rest of functions
    #
    curJ = 0.0
    if haveJ == 1:
        curJ += np.interp(x, Jrad[0,:], Jrad[1,:])
    """"""
    #
    # Calculate the gas constants first
    #
    vgas    = w[0]
    Tgas    = w[1]
    Rhogas  = 0.0
    Ngas    = 0.0
    for ii in range(nspecs):
        Rhogas  += w[<unsigned int> (2+ii)] * gmasses[<unsigned int> ii]
        Ngas    += w[<unsigned int> (2+ii)]
    #
    # Gamma
    #
    gam1    = 0.0
    gam     = 0.
    for ii in xrange(nspecs):
        gam1 += w[<unsigned int> (2+ii)] * (gammas[ii] *2. - 2.)
        gam  += (gam1/Ngas + 2.)/(gam1/Ngas)
    #
    # Dust changes
    #
    cdef np.ndarray[DTYPE_t, ndim=1] mdust  = np.zeros(ndust, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] fdrag  = np.zeros(ndust, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] dxa    = np.zeros(ndust, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] dxtd   = np.zeros(ndust, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] Cpd    = np.zeros(ndust, dtype=DTYPE)
    dustLoss    = 0.0
    #
    # Dust position is 2 + gas species + dusttype*4
    # since each dust takes 4 positions
    #
    for idust in range(ndust):
        nd  = w[<unsigned int> (2+nspecs + idust*4 + 0)]
        vd  = w[<unsigned int> (2+nspecs + idust*4 + 1)]
        td  = w[<unsigned int> (2+nspecs + idust*4 + 2)]
        ad  = w[<unsigned int> (2+nspecs + idust*4 + 3)]
        mdust[idust]    = (4.0/3.0)*M_PI*rhod*pow(ad, 3.)
        fdrag[idust]    = calcFDrag(Rhogas/Ngas, nd, vd, td, ad, vgas, Tgas,
            Rhogas)
        dxa[idust]      = calcDXa(nd, vd, td, ad, vgas, gam, Tgas, curJ, rhod,
            Rhogas, Ngas)
        dxtd[idust]     = calcDxTd(nd, vd, td, ad, vgas, Tgas, Rhogas,
            Ngas, gam, curJ, rhod)
        Cpd[idust]      = getCpd(td)
        if dxa[idust] > 0.0:
            dustLoss    += 0.0
        else:
            dustLoss    += -4.0 * M_PI * (nd * vd * rhod * ad*ad * dxa[idust])
        """"""
    """"""
    #
    # Calculate the variable here
    #
    varA    = Rhogas*w[0] - (kk*w[1]/w[0])*Ngas
    varB    = Ngas*kk
    varC    = 0.0 # reset variable C
    #
    # Variable C now has dust properties
    # Calculate the gas part
    #
    cdef ndarray[DTYPE_t, ndim=1] RatesMat = (calculateR(Tgas, nspecs,
        w[2], w[3])).sum(axis=0)
    dumGas  = 0.0 # Reset gas component
    for ii in range(nspecs):
        dumGas += RatesMat[ii]*(w[0]*gmasses[ii] + (kk*w[1]/w[0]))
    #
    # Add the dust evaporation mass
    #
    dumGas += dustLoss/(gmasses[3])*(w[0]*gmasses[3] + (kk*w[1]/w[0]))
    #
    # Get the dust part
    #
    dumDust     = 0. # reset the dust part
    for idust in range(ndust):
        dumDust     += (w[<unsigned int> (2+nspecs+idust*4 +0)] *
            (fdrag[idust] + 4.0*M_PI *
            pow(w[<unsigned int> (2+nspecs+idust*4 + 3)],2) * rhod *
            pow(w[<unsigned int> (2+nspecs+idust*4+1)], 2.) * dxa[idust]) )
    varC    += -(dumGas + dumDust)
    #
    # Calculate the variables for energy equation
    #
    varD    = w[0]*w[0]*Rhogas
    varE    = 0.0
    for ii in range(nspecs):
        varE += gammas[ii] * w[2+ii]
    varE    *= kk*w[0]
    varF    = 0.0 # reset variable F
    #
    # Calculate the variable F with the dust
    #
    dumGas  = 0.0 # reset gas
    for ii in range(nspecs):
        dumGas  += (RatesMat[ii]*(gammas[ii]*kk*w[1] +
            0.5*w[0]*w[0]*gmasses[ii]))
    dumGas += dustLoss/(gmasses[3]) * (gammas[3]*kk*w[1] +
        0.5*w[0]*w[0]*gmasses[3])
    #
    # Dust part for variable F
    #
    dumDust = 0.0 # reset dust
    for idust in range(ndust):
        dumDust     += ( w[<unsigned int> (dumid + idust*4+1)] *
            w[<unsigned int> (dumid+idust*4+0)] *
            (fdrag[idust] + 2.*M_PI*
            pow(w[<unsigned int> (dumid+idust*4+3)], 2.) * rhod *
            pow(w[<unsigned int> (dumid+idust*4+1)], 2.) *dxa[idust]) +
            w[<unsigned int> (dumid+idust*4+0)] * (
            mdust[idust] * w[<unsigned int> (dumid+idust*4+1)] * Cpd[idust] *
            dxtd[idust]) + (4.0 * M_PI * (
            pow(w[<unsigned int> (dumid+idust*4+3)],2.) *
            w[<unsigned int> (dumid+idust*4+0)] * Cpd[idust] *
            w[<unsigned int> (dumid+idust*4+1)] *
            w[<unsigned int> (dumid+idust*4+2)] * rhod *dxa[idust])) )
    varF    += - (dumGas + dumDust) + calculateFreeEnergy(
        Tgas, nspecs, w[2], w[3])
    #
    # Radiation part
    #
    radiation = 0.0
    if haveJ == 1:
        radiation   += (4.0 * M_PI * Rhogas *
            gasKap(w[1], destroyed, Tkap, Kaps) *
            (curJ - (ss*pow(w[1], 4.)/M_PI)))
        for idust in range(ndust):
            radiation   += (4.0* M_PI * M_PI *
                w[<unsigned int> (dumid+4*idust+0)] *
                pow(w[<unsigned int> (dumid+4*idust+3)], 2.) *
                getEPS(w[<unsigned int> (dumid+4*idust+3)]) *
                (curJ- (ss * pow(w[<unsigned int> (dumid+4*idust+2)], 4.)/
                M_PI)))
    varF += radiation # add the radiation
    #
    # The RHS matrix
    # Check few numbers to make sure these are not too large
    # or too small
    #
    YDOT[0]     = (varC*varE - varB*varF)/(varA*varE - varB*varD)
    YDOT[1]     = (varA*varF - varD*varC)/(varA*varE - varB*varD)
    #
    # Integration is over space dx
    # The change of n * v
    # changes in chemical species
    # Add the dust evaporation
    #
    for ii in range(nspecs):
        YDOT[<unsigned int> (2+ii)]     += RatesMat[ii]
        YDOT[<unsigned int> (2+ii)]     += (gdustFrac[ii]*dustLoss/gmasses[ii])
        YDOT[<unsigned int> (2+ii)]     -= (w[2+ii] * YDOT[0])
        YDOT[<unsigned int> (2+ii)]     /= w[0]
    #
    #
    # DUST from here on
    #
    if ndust > 0:
        for idust in xrange(ndust):
            nd  = w[<unsigned int> (2+nspecs + idust*4 + 0)]
            vd  = w[<unsigned int> (2+nspecs + idust*4 + 1)]
            td  = w[<unsigned int> (2+nspecs + idust*4 + 2)]
            ad  = w[<unsigned int> (2+nspecs + idust*4 + 3)]
            YDOT[<unsigned int> (dumid + idust*4 + 1)]  += (
                fdrag[idust]/(mdust[idust]*vd))
            YDOT[<unsigned int> (dumid + idust*4 + 0)]  += (
                - (nd/vd) * YDOT[<unsigned int> (dumid + idust*4 + 1)])
            YDOT[<unsigned int> (dumid + idust*4 + 2)]  += dxtd[idust]
            YDOT[<unsigned int> (dumid + idust*4 + 3)]  += dxa[idust]
            #
            # Crash and BURN handles
            #
            if (isnan(YDOT[<unsigned int> (dumid + idust*4 + 1)])):
                print 'NAN is found!'
                print 'Tgas and vgas: %e, %e'%(w[1], w[0]*1e-5)
                print 'ndust, vdust and Tdust: ', w[2+nspecs:]
                print
                print fdrag, mdust
                print w
                raise SystemExit
            if (isnan(YDOT[<unsigned int> (dumid + idust*4 + 0)])):
                print 'NAN is found!'
                print 'Tgas and vgas: %e, %e'%(w[1], w[0]*1e-5)
                print 'ndust, vdust and Tdust: ', w[2+nspecs:]
                print
                print fdrag, mdust
                print w
                raise SystemExit
            if (isnan(YDOT[<unsigned int> (dumid + idust*4 + 2)])):
                print 'NAN is found!'
                print 'Tgas and vgas: %e, %e'%(w[1], w[0]*1e-5)
                print 'ndust, vdust and Tdust: ', w[2+nspecs:]
                print
                print fdrag, mdust
                print w
                raise SystemExit
            if (isnan(YDOT[<unsigned int> (dumid + idust*4 + 3)])):
                print 'NAN is found!'
                print 'Tgas and vgas: %e, %e'%(w[1], w[0]*1e-5)
                print 'ndust, vdust and Tdust: ', w[2+nspecs:]
                print
                print fdrag, mdust
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
    return YDOT
""""""
#
# Function to be called from outside
#
cpdef solve(double x, np.ndarray[DTYPE_t, ndim=1] w, list p):
    """
    Here is the function to pass to the cython function
    """
    #
    # Variables
    #
    cdef int nspecs, ndust, haveJ, destroyed
    cdef ndarray[DTYPE_t, ndim=1] w0    = np.zeros(w.shape[0], dtype=DTYPE)
    #
    # Rescale the input
    #
    nspecs  = <unsigned int> (p[0].nspecs)
    ndust   = <unsigned int> (p[1].nspecs)
    w0  = decipherW(w, nspecs, ndust)
    #
    # return YDOT
    #
    haveJ = 0
    if p[3]:
        haveJ = 1
    destroyed = 0
    if p[1].destroyed:
        destroyed = 1
    return (vectorfield(x, w0, p[2], haveJ, nspecs, ndust,
        np.array(p[0].mass), np.array(p[0].gamma),
        p[1].mdust, destroyed, p[0].logT, p[0].logK,
        np.array(p[0].dustfrac)))
""""""

