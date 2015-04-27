from python.my_header import *
from scipy.integrate import odeint
from scipy.special import expn

"""
Global variables here
"""
isshock = 0
def set_shock_occur():
    global isshock
    isshock = 1
""""""

def setOldT(a):
    global oldT
    oldT = a
""""""

def setdX(a,b):
    global dx, oldx
    oldx = a
    dx = b
""""""

"""
Get the planck mean opacities
"""
def getKap(temp=None):
    if (temp < 1000):
        return 5.0
    elif ( (temp >1000.) and (temp < 3000.)):
        return 0.5
    else:
        return 0.05

"""
Calculate the mean intensity: solve the radiative transfer
"""
def getJrad(temp1=None, temp2=None, alpha=None, step=None):
    preconst = ss/np.pi
    first = preconst*expn(getTau(), 2)
    
"""
Solving the Rankine Hugionot conditions
"""
def get_RHcondition(rho1=None, u1=None, P1=None, gamma=3.0/2.0):
    """
    Given the initial conditions of the gas
    Return the post shock conditions
    """
    # solve for velocity after the shock front
    # Quadratic equation solver: ax^2 + bx + c = 0
    gam1 = (gamma/(gamma-1.0))
    aa = (1.0/2.0) - gam1
    bb = (gam1*u1) + gam1*P1/(u1*rho1)
    cc = -(1.0/2.0 * u1*u1 + gam1 * P1/rho1)
    
    discrim = bb*bb - 4.*aa*cc
    
    # check the discriminant
    if discrim < 0:
        print 'This has no solution'
        raise ValueError
    elif discrim == 0:
        u2 = (-bb + np.sqrt(bb*bb - 4.0*aa*cc))/(2.0*aa)
    else:
        u21 = (-bb + np.sqrt(bb*bb - 4.0*aa*cc))/(2.0*aa)
        u22 = (-bb - np.sqrt(bb*bb - 4.0*aa*cc))/(2.0*aa)
        if (u21 < u1):
            u2 = u21
        else:
            u2 = u22
    """"""
    
    # Solve the pressure
    P2 = rho1*u1*u1 + P1 - rho1*u1*u2
    # Solve rho2
    rho2 = rho1 * u1/u2
    print 'Density: %2.5e -->  %2.5e'%(rho1, rho2)
    print 'Pressure: %2.5e --> %2.5e'%(P1, P2)
    print 'Velocity: %2.5e --> %2.5e'%(u1, u2)
    print 'Temperature: %2.5e --> %2.5e'%(u1, u2)    
    
    return P2, rho2, u2
""""""

"""
The vector field of the matrix to solve
with radiative transfer
"""
def vectorfieldrt(w, x, p):
    """
    w   :  the variables to solve -> x, y 
    x   :  position with respect to the shock
    p   :  parameters
    """
    
    x1, x2 = w
    numden, mass, mugas, gamma = p
    varA = 2.0*numden*mass
    varB = numden*kk
    varC = 0.
    varD = (gamma*numden*mass*x2 + numden*mass + (numden*kk/
        mugas) * x2)
    varE = (gamma*numden*kk*mass*x1 + (numden*kk/mugas)*x1)
    """
    Set the global variables
    """
    setOldT(x2)
    setdX(x, x-oldx)
    varF = -4.0*nppi*getKap(x2) * (getJrad(oldT, x2, dx) - 
            (ss/np.pi)*x2**(4.))
    
    f = [ (varA*varF - varD*varC)/( varA*varE - varB*varD ),
         (varC*varE - varB*varA)/(varA*varE - varB*varD) ]
    return f
""""""

""" 
Vector field without radiative transfer
This needs to be run before with radiative transfer
to estimate the post shock
"""
def vectorfield(w, x, p):
    """
    w   :  the variables to solve -> x, y 
    x   :  position with respect to the shock
    p   :  parameters
    """
    
    x1, x2 = w
    numden, mass, mugas, gamma = p
    print x, x1, x2
    if x == 0:
        """ # put in the shock """
        P1 = (kk/(mugas*mp) * numden * 2.0*mp * x2)
        P2, rho2, u2 = get_RHcondition(rho1 = 2.0*mp*numden, 
            u1 =x1, P1=P1)
        numden = rho2/(2.0*mp)
        print 'Here'
        p = [numden, mass, mugas,gamma]
        x1 = u2
        x2 = P2*(mugas*mp)/(kk*numden * 2.0*mp)        
    varA = 2.0*numden*mass
    varB = numden*kk
    varC = 0.
    varD = (gamma*numden*mass*x2 + numden*mass + (numden*kk/
        mugas) * x2)
    varE = (gamma*numden*kk*mass*x1 + (numden*kk/mugas)*x1)
    varF = 0.0
    f = [ (varA*varF - varD*varC)/( varA*varE - varB*varD ),
         (varC*varE - varB*varA)/(varA*varE - varB*varD) ]
    return f
""""""

"""
Solver and parameters
"""
def solveHD(x=None, numden = None, mass = None, mugas=2.8, v0 = None, t0 = None, dv = 0.0, dT = 0.0, gamma=3.0/2.0, abserr=1e-8, telerr = 1e-6, sizex = 10.0, numpoints = 500, rt=False):
    """
    Call the solver with initial conditions v0 and t0
    """
    x0 = [np.abs(sizex) * (np.float(i)/(numpoints - 1)) - sizex 
         for i in xrange(numpoints)]
    x1 = [np.abs(sizex) * 
         (np.float(i)/(numpoints - 1)) for i in xrange(numpoints)]
    x =  x0[:-1] + [0.0] + x1[1:]
    p = [numden, mass, mugas, gamma]
    w0 = [v0, t0]
    
    # call the solver
    wsol = odeint(vectorfield, w0, x, args=(p,),
                  atol=abserr, rtol=telerr)
    
    return x, wsol
""""""

"""
###########################################################
"""
""" parameters """
setdX(-5.0, 0.0)
numden = 1e14
#mass = 2.0*mp
v0 = 6e5
t0 = 300
mugas = 2.8

x, shocksol = solveHD(numden=numden, mass=2.0*mp, v0=v0, t0=t0,               
    sizex = 5.0, numpoints=3)
solutions = []
for x1, w1 in zip(x, shocksol):
    solutions.append([x1, w1[0], w1[1]])
""""""
## put in the shock
#P1 = (kk/(mugas*mp) * numden * 2.0*mp * solutions[-1][2])
#P2, rho2, u2 = get_RHcondition(rho1 = 2.0*mp*numden, 
    #u1 = solutions[-1][1], P1=P1)
#numden = rho2/(2.0*mp)
#v0 = u2
#t0 = P2*(mugas*mp)/(kk*numden * 2.0*mp)
x, shocksol = solveHD(numden=numden, mass=2.0*mp, v0=v0, t0=t0,               
    sizex = 5.0, numpoints=20)
for x1, w1 in zip(x, shocksol):
    solutions.append([x1, w1[0], w1[1]])
""""""
solutions = np.array(solutions)

#print
#print 'Tpre = %8.2f, Tpost = %8.2f'%(solutions[1][2], solutions[-2][2])
#print

#"""
#Perform the runs with radiative transfer
#"""


#raise SystemError


fig, ax0 = subplots(1,1, figsize=[aaonecol, aaonecol])
subplots_adjust(left=0.15, right=0.8, top=0.98, bottom=0.12)

ax0.plot(solutions[:,0], solutions[:,1]*1e-5, 'b-', lw=1.2)
ax0 = fig_labs(ax=ax0, xlab=r'\boldmath$x$',
    ylab=r'\boldmath$\upsilon_{\rm gas}$', fontsize=8,
    xlim=[-5.0,5.0], xform=r'\boldmath$%1.1f$', yform=r'\boldmath$%1.2f$')

ax1 = ax0.twinx()
ax1.plot(solutions[:,0], solutions[:,2], 'g-', lw=1.2)
ax1 = fig_labs(ax=ax1, xlab=r'\boldmath$x$',
    ylab=r'\boldmath$T_{\rm gas}$', fontsize=8,
    xlim=[-5.0,5.0], xform=r'\boldmath$%1.1f$', yform=r'\boldmath$%d$')

for t1 in ax1.get_yticklabels():
    t1.set_fontsize(8)
    t1.set_color('g')

#leg = ax0.legend(loc='upper right')
#leg.draw_frame(False)

fig.savefig('shocktest.pdf', dpi=500)
close()
os.system('okular shocktest.pdf')



