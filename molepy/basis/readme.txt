Available basis sets and there name which molepy recognizes,

cc-pVDZ	                vdz
aug-cc-pVDZ		avdz
cc-pVTZ			vtz
aug-cc-pVTZ		avtz
def2-SVP		def2-svp
def2-SVPD		def2-svpd
def2-TZVP		def2-tzvp
def2-TZVPD		def2-tzvpd
STO-3G			sto3g

Incase your favourite basis set is not available,
generate one to use in molepy with the script 'genbasis.py'.

Just run,

genbasis.py basisdata basisname

'basisdata' should be a file containing basis functions
in GAMESS-US format. This can be obtained easily from 
https://bse.pnl.gov/bse

'basisname' should be the name of the basis set you would
want to use in molepy.
