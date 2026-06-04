# EOM-DLPNO-CC

This repository will be where I put work relating to my dissertation research. The core of this research surrounds implementing general-order DLPNO-CC(n), DLPNO-CC(n)(n - 1), DLPNO-CCn, EOM-DLPNO-CC(n), and EOM-DLPNO-CCn, with excitations, properties, ionization energies, electron attachment energies, and spin-flip energies from EOM-DLPNO-CC.

## General-order DLPNO-CC

To solve the problem of general-order DLPNO-CC, a general-order CC method based on the work of Kállay can be used. This can then be converted to the DLPNO formalism. Work on the standard DLPNO-CC(n) and DLPNO-CC(n)(n - 1) has yielded good results. To perform the conversion, the equations will first need to be transformed into the density-fitted forms.
Then, we can use localization tensors to convert the terms into full DLPNO. For DLPNO-CCn, the cutoffs will be treated in a similar way to DLPNO-CC(n)(n - 1) with the main difference being that perturbative corrections will use CCn instead of CC(n)(n - 1). DLPNO-CC2 will end up replacing LPNO-MP2 in the EOM-DLPNO-CCn calculations. It may also end up
replacing LPNO-MP2 in the regular DLPNO-CCn as well, but this will need to be compared.

## DLPNO-CC Properties

Due to the fact that the energy derivatives from DLPNO-CC may not exist, only non-derived properties may be considered. This would include things like multipole moments, but we would be unable to get good gradients. For this, we would need the lambda equations.

## EOM-DLPNO-CC

EOM-DLPNO-CC(n)(n - 1) does not really exist, so only EOM-DLPNO-CC(n) and EOM-DLPNO-CCn will be developed. Since DLPNO-CCSD, CCSD(T), CCSDT, CCSDT(Q), and CCSDTQ have been developed already, some combination of LPNO-CC2, DLPNO-CC3, DLPNO-CC4, EOM-LPNO-CC2, EOM-DLPNO-CCSD, EOM-DLPNO-CC3, EOM-DLPNO-CCSDT, EOM-DLPNO-CC4, and EOM-DLPNO-CCSDTQ
will likely be implemented by hand, while the rest of the methods will be implemented using antisymmetrized diagrams. Note that LPNO-CC2 would probably be better than DLPNO-CC2 since we already have all of the singles amplitudes that we would need. 
