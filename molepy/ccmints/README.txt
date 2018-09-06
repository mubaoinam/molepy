
Molepy evaluates one-electron and two-electron integrals 
using the Obara-Saika scheme over Cartesian GTOs. 
The integrals are then transformed to ones over spherical harmonic GTOs.
This can be turned off by setting 

hfparam.Spheri = False

in molepy/lib/mole_param.py.

The transformation to spherical harmonics utilizes pre-computed
coefficients for angular-momentum quantum number less than 4.


Molepy uses the codes in ccmints but you can choose to turn it off by setting

hfparam.ccmints = False

in molepy/lib/mole_param.py.

Evaluation of Boys function with ccmints uses pre-computed gamma values, this
speeds up the calculation to some extend.
