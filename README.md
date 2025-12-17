# HLINOP

- New version with bug fix thanks to Richard Hoppe, 11 July 2023; see history.txt
- hbop_auto.f added 17 December 2025; small changes to hbop.f

A set of codes for computing hydrogen line profiles and opacities developed by Paul Barklem and Nik Piskunov.  This code was first presented at IAU 210, Modelling Stellar Atmospheres, in Uppsala June 2002 (Barklem & Piskunov 2003). There a need for codes to compute these opacities for atellar atmosphere applications. Codes such as HLINOP and BALMER9 have been provided by Kurucz and Peterson (http://kurucz.harvard.edu/) in the 1970's based on the broadening theories of the time.  Since then, understanding of H line broadening has improved and there is a need to update these codes. Three codes have been produced, namely:

HLINPROF: For detailed, accurate calculation of lower Balmer line profiles, suitable for detailed analysis of Balmer lines.  This code depends on binary files containing tables and are thus provided in little- and big-endian versions.  There is also a version (_nofiles) which attempts to avoid the need for external files, but is not as well tested as the default version.

HLINOP: For approximate, quick calculation of any line of neutral H, suitable for model atmosphere calculations.  Based on the original Kurucz and Peterson version.

HBOP: A wrapper around hlinop.f to implement the occupation probability formalism of Daeppen, Anderson and Milhalas (1987) and thus account for the merging of lines into the continuum.  As bound-free and bound-bound opacities are merged, the distinction between line and continuous opacity becomes blurred, so one should be careful not to include things twice (including Rayleigh scattering off Lyman alpha). For *efficient* implementation you may have to adapt this, in particular so that populations are not recomputed.

HBOP_AUTO: This code is a wrapper to HBOP, which will give the total bound-x (bound-bound and bound-free) neutral H opacity at a given wavelength, from the occupation probability formalism, WITHOUT the need for a line list.  This wrapper will simply produce a line list including all lines up to a highest principal quantum number NMAX, in the relevant spectral region (i.e. at the interrogation wavelength).  

The main codes are in the root directory.  An example of their application to the SYNTH spectral synthesis program by Nik Piskunov is also provided.

Views, comments, bug reports and suggested improvements all welcome: paul.barklem@physics.uu.se

Please cite the DOI if you make use of the code in research:
[![DOI](https://zenodo.org/badge/21607/barklem/hlinop.svg)](https://zenodo.org/badge/latestdoi/21607/barklem/hlinop)

References:

- Barklem, P. S., and N. Piskunov. “Hydrogen Balmer Lines as Probes of Stellar Atmospheres,” Vol. 210, 2003. http://adsabs.harvard.edu/abs/2003IAUS..210P.E28B.

- Dappen, Werner, Lawrence Anderson, and Dimitri Mihalas. “Statistical Mechanics of Partially Ionized Stellar Plasma - The Planck-Larkin Partition Function, Polarization Shifts, and Simulations of Optical Spectra.” The Astrophysical Journal 319 (August 1, 1987): 195–206. doi:10.1086/165446.





