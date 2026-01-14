source leaprc.gaff2
loadoff NP_correct.lib
loadamberparams NP_charged.frcmod
loadamberparams NP_neut.frcmod
NP = loadpdb NP_dry.pdb
charge NP
quit

