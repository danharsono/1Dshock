from python.my_header import *
import main 
"""
This is where the shock code is being run
"""
ndust = 1
ngas = 3
vshock=6.5e5
solutions, Jrad = main.shock_main(sizex=3e10, numpoints=8e3, nspecs=ngas,
    ndust=ndust, v0=vshock, niter=0, ncpu=2)

#x0 = solutions[:,0]/solutions[:,1]
x0 = solutions[:,0]/vshock
x0 = x0/3600.0

#raise SystemExit

fig, ax0 = subplots(1,1, figsize=[aaonecol, aaonecol])
subplots_adjust(left=0.15, right=0.8, top=0.98, bottom=0.12)

ax0.plot(x0, solutions[:,1]*1e-5, 'b-', lw=1.2, alpha=0.6)
if ndust is not None:
    ax0.plot(x0, solutions[:,7]*1e-5, 'b--', lw=1.2)
ax0 = fig_labs(ax=ax0, xlab=r'\textbf{\boldmath$t$ [h]}',
    ylab=r'\boldmath$\upsilon_{\rm gas}$', fontsize=8,
    xlim=[-15,10], ylim=[0.0,10.0], yminloc=0.5, xminloc=1.0,
    xform=r'\boldmath$%1.1f$', yform=r'\boldmath$%d$')
for t1 in ax0.get_yticklabels():
    t1.set_fontsize(8)
    t1.set_color('b')
ax0.grid(lw=0.5, color='0.6', ls=':', alpha=0.6)
f
ax1 = ax0.twinx()
ax1.plot(x0, solutions[:,2], 'g-', lw=1.2, alpha=0.6)
if ndust is not None:
    ax1.plot(x0, solutions[:,8], 'r-', lw=1.2, alpha=0.6)
ax1 = fig_labs(ax=ax1, xlab=r'\textbf{\boldmath$t$ [h]}',
    ylab=r'\boldmath$T_{\rm gas}$', fontsize=8, xlim=[-15, 10],
    ylim=[0, 3500.0],
    xform=r'\boldmath$%1.1f$', yform=r'\boldmath$%d$')

for t1 in ax1.get_yticklabels():
    t1.set_fontsize(8)
    t1.set_color('g')

#leg = ax0.legend(loc='upper right')
#leg.draw_frame(False)

fig.savefig('shocktest.pdf', dpi=500)
close()
#os.system('okular shocktest.pdf')


fig, ax0 = subplots(1,1, figsize=[aaonecol, aaonecol])
subplots_adjust(left=0.15, right=0.85, top=0.98, bottom=0.12)

ax0.semilogy(x0, solutions[:,3], 'k-', lw=1.2, 
    label=r'\boldmath${\rm H}$')
ax0.semilogy(x0, solutions[:,4], 'b-', lw=1.2, 
    label=r'\boldmath${\rm H_{2}}$')
ax0.semilogy(x0, solutions[:,5], 'g-', lw=1.2, 
    label=r'\boldmath${\rm He}$')

if ndust is not None:
    ax0.semilogy(x0, solutions[:,6]*1e8, 'k--', lw=1.2,
    label=r'\boldmath${\rm d}$')

ax0 = fig_labs(ax=ax0, xlab=r'\textbf{\boldmath$t$ [h]}',
    ylab=r'\boldmath$n$', fontsize=8, xminloc=1.0, 
    xlim=[-15,10], xform=r'\boldmath$%1.1f$', yform=1)

leg = ax0.legend(loc='upper left')
leg.draw_frame(False)
ax0.grid(lw=0.5, color='0.6', ls=':', alpha=0.6)

fig.savefig('numden.pdf', dpi=500)
close()
#os.system('okular numden.pdf')


fig, ax0 = subplots(1,1, figsize=[aaonecol, aaonecol])
subplots_adjust(left=0.15, right=0.85, top=0.97, bottom=0.11)

ax0.semilogy(x0, Jrad[0], 'b-', lw=1.2)
ax0 = fig_labs(ax=ax0, xlab=r'\textbf{\boldmath$t$ [h]}',
    ylab=r'\boldmath$\tau$', fontsize=8, xminloc=1.0,
    xlim=[-15,10], xform=r'\boldmath$%1.1f$', yform=1)
for t1 in ax0.get_yticklabels():
    t1.set_fontsize(8)
    t1.set_color('b')
ax0.grid(lw=0.5, color='0.6', ls=':', alpha=0.6)

ax1 = ax0.twinx()
ax1.plot(x0, np.log10(Jrad[1][1]), 'go', lw=1.2, ms=5)
ax1 = fig_labs(ax=ax1, xlab=r'\textbf{\boldmath$t$ [h]}',
    ylab=r'\boldmath$J_{\rm rad}$', fontsize=8, xlim=[-15, 10], 
    xform=r'\boldmath$%1.1f$', yform=r'\boldmath$%d$',xminloc=1.0,
    ylim=[-5,20.], yminloc=2.0)

for t1 in ax1.get_yticklabels():
    t1.set_fontsize(8)
    t1.set_color('g')

#leg = ax0.legend(loc='upper right')
#leg.draw_frame(False)

fig.savefig('radiation.pdf', dpi=500)
close()
#os.system('okular radiation.pdf')

