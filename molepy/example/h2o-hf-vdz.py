#!/usr/bin/env python
from molepy.coord.coord import Mole
from molepy.scf.rhf import RHF
from molepy.lib import mout

""" example h2o/cc-pVDZ

1. define geometry with Mole([list of atomic coordinates])
   supply list in the format as shown below.
2. molog() writes an output file, name the output file using
   set_name('NAME').
3. call RHF(geometry,'basis set name')
   RHF returns total energy, orbital energies and coefficients.
4. Basis sets:
   'vdz' : cc-pVDZ
   'vtz' : cc-pVTZ, etc...
   check molepy/basis/readme.txt for basis sets.

"""

h2o = Mole([
    ('O',   -0.317069806115,    3.11002116482,    0.205500963755),
    ('H',    1.15917268226 ,    2.08363883745,    0.214906081107),
    ('H',   -0.177295060282,    4.17574360331,   -1.21625722816 )])

out = mout.molog()
out.set_name('h2o-hf-vdz.log')

E,orbe,orbc = RHF(h2o,'vdz')
print 'HF energy ',E
