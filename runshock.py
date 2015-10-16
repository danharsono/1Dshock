from python.my_header import *
import main
import fileio

"""
This is where the shock code is being run
"""
ndust = 1
ngas = 4
vshock=6.5e5
tgas = 300.
#
# SOLVE
#
solutions, Jrad, vshock, Frad = main.shock_main(sizex=3e10, numpoints=3e3,
    nspecs=ngas, ndust=ndust, v0=vshock, niter=2, ncpu=3,
    t0=tgas, restart=False, Tpost0=1100.)
#print vshock
#
# Write out the solutions
#
fileio.writeOut(solutions=solutions, Jrad=Jrad, Frad=Frad,
    vshock=vshock, ngas=ngas, ndust=ndust, fname='shockout')
#
# Read the output
#
sols = fileio.readOut(fname='shockout')
#
# position is relative to the shock front
#
x0 = sols['grid']/sols['gas'].attrs['vshock']
x0 = x0/3600.0
#
# Plot these
#
xlims = [-15, 15]
fig, ax0 = subplots(1,1, figsize=[aaonecol, aaonecol])
subplots_adjust(left=0.15, right=0.8, top=0.98, bottom=0.12)

ax0.plot(x0, sols['gas'][:,0]*1e-5, 'b-', lw=1.2, alpha=0.6,
    label=r'\boldmath$\upsilon_{\rm gas}$')
if ndust is not None:
    ax0.plot(x0, sols['dust'][:,0]*1e-5, 'b-.', lw=1.2,
        label=r'\boldmath$\upsilon_{\rm dust}$')
ax0 = fig_labs(ax=ax0, xlab=r'\textbf{\boldmath$t$ [h]}',
    ylab=r'\boldmath$\upsilon$', fontsize=8,
    xlim=xlims, ylim=[0.0,10.0], yminloc=0.5, xminloc=1.0,
    xform=r'\boldmath$%1.1f$', yform=r'\boldmath$%d$')
for t1 in ax0.get_yticklabels():
    t1.set_fontsize(8)
    t1.set_color('b')
ax0.grid(lw=0.5, color='0.6', ls=':', alpha=0.6)

ax1 = ax0.twinx()
ax1.plot(x0, sols['gas'][:,1], 'g-', lw=1.2, alpha=0.6,
    label=r'\boldmath$T_{\rm gas}$')
if ndust is not None:
    ax1.plot(x0, sols['dust'][:,1], 'r-', lw=1.2, alpha=0.6,
        label=r'\boldmath$T_{\rm dust}$')
ax1 = fig_labs(ax=ax1, xlab=r'\textbf{\boldmath$t$ [h]}',
    ylab=r'\boldmath$T_{\rm gas}$', fontsize=8, xlim=xlims,
    ylim=[0, 3500.0],
    xform=r'\boldmath$%1.1f$', yform=r'\boldmath$%d$')

for t1 in ax1.get_yticklabels():
    t1.set_fontsize(8)
    t1.set_color('g')

leg = ax0.legend(loc='upper left', fontsize=7)
leg.draw_frame(False)
leg1 = ax1.legend(loc='upper right', fontsize=7)
leg1.draw_frame(False)

fig.savefig('shocktest.pdf', dpi=500)
close()
os.system('open shocktest.pdf')
#
# Number densities and chemistry
#
fig, ax0 = subplots(1,1, figsize=[aaonecol, aaonecol])
subplots_adjust(left=0.15, right=0.85, top=0.98, bottom=0.12)

ax0.semilogy(x0, sols['gas'][:,2], 'k-', lw=1.2,
    label=r'\boldmath${\rm H}$')
ax0.semilogy(x0, sols['gas'][:,3], 'b-', lw=1.2,
    label=r'\boldmath${\rm H_{2}}$')
ax0.semilogy(x0, sols['gas'][:,4], 'g-', lw=1.2,
    label=r'\boldmath${\rm He}$')
if ngas == 4:
    ax0.semilogy(x0, sols['gas'][:,5], 'b--', lw=1.2,
        label=r'\boldmath${\rm SiO}$')
if ndust is not None:
    ax0.semilogy(x0, sols['dust'][:,2]*1e15, 'r--', lw=1.2,
    label=r'\boldmath${\rm d}$')

ax0 = fig_labs(ax=ax0, xlab=r'\textbf{\boldmath$t$ [h]}',
    ylab=r'\boldmath$n$', fontsize=8, xminloc=1.0, 
    xlim=xlims, xform=r'\boldmath$%1.1f$', yform=1)

leg = ax0.legend(loc='upper left', fontsize=6)
leg.draw_frame(False)
ax0.grid(lw=0.5, color='0.6', ls=':', alpha=0.6)

fig.savefig('numden.pdf', dpi=500)
close()
os.system('open numden.pdf')
#
# Radiation field and optical depth
#
fig, ax0 = subplots(1,1, figsize=[aaonecol, aaonecol])
subplots_adjust(left=0.15, right=0.85, top=0.97, bottom=0.11)

ax0.semilogy(x0, sols['radiation'][:,0], 'b-', lw=1.2)
ax0 = fig_labs(ax=ax0, xlab=r'\textbf{\boldmath$t$ [h]}',
    ylab=r'\boldmath$\tau$', fontsize=8, xminloc=1.0,
    xlim=xlims, xform=r'\boldmath$%1.1f$', yform=1)
for t1 in ax0.get_yticklabels():
    t1.set_fontsize(8)
    t1.set_color('b')
ax0.grid(lw=0.5, color='0.6', ls=':', alpha=0.6)

ax1 = ax0.twinx()
ax1.plot(x0, np.log10(sols['radiation'][:,1]), 'g-', lw=1.2, ms=5)
ax1.plot(x0, np.log10(np.abs(sols['radiation'][:,2])), 'r--', lw=1.2, ms=5)
ax1 = fig_labs(ax=ax1, xlab=r'\textbf{\boldmath$t$ [h]}',
    ylab=r'\boldmath$J_{\rm rad}$', fontsize=8, xlim=xlims,
    xform=r'\boldmath$%1.1f$', yform=r'\boldmath$%d$',xminloc=1.0,
    ylim=[-5,20.], yminloc=2.0)

for t1 in ax1.get_yticklabels():
    t1.set_fontsize(8)
    t1.set_color('g')

#leg = ax0.legend(loc='upper right')
#leg.draw_frame(False)

fig.savefig('radiation.pdf', dpi=500)
close()
os.system('open radiation.pdf')

