import dustspec, gasspec
import numpy as np
from python.natconst import *
from python.my_header import *

gas = gasspec.gasSpecs()
dust = dustspec.dustSpecs(gas=gas)

s = np.logspace(-9, -1, 500)
Tg = 3.000405e2
Td = 3.000405e2
vd = 6.49999e5
vg = np.linspace(6.3e5, 6.5e5, 500)

qall = np.zeros(s.shape)
dtall = np.zeros(s.shape)
for ias in xrange(s.shape[0]):
    qall[ias]  = dust._calcqdust(Tg=Tg, Td=Td, s=s[ias], mass=1.*mp)
    dtall[ias] = dust._calcDxTd(vd=vd, vg=vg[ias], Tg=Tg, Td=Td, gas=gas)
    print qall[ias], dtall[ias]

semilogx(s,qall, 'ko')
show()
close()

plot(vg*1e-5,dtall, 'ko')
show()
close()

