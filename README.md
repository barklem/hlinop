# hlinop

A set of codes for computing H line profiles and opacities has been developed by myself and Nik Piskunov, and was first presented at IAU 210, Modelling Stellar Atmospheres, in Uppsala June 2002. There is clearly a need for freely available code to compute these opacities. Codes such as HLINOP and BALMER9 have been provided by Kurucz and Peterson (http://cfaku5.harvard.edu/) in the 1970's based on the broadening theories of the time. Since then, understanding of H line broadening has improved significantly and so there is a need to update these codes. Three codes have been produced, namely:

HLINPROF: For detailed, accurate calculation of lower Balmer line profiles, suitable for detailed analysis of Balmer lines.  This code depends on binary files containing tables and are thus provided in little- and big-endian versions.  There is also a version (_nofiles) which attempts to avoid the need for external files, but is not as well tested as the default version.

HLINOP: For approximate, quick calculation of any line of neutral H, suitable for model atmosphere calculations.

HBOP:A wrapper around hlinop.f to implement the occupation probability formalism of Daeppen, Anderson and Milhalas and thus can account for the merging of lines into the continuum.  As bound-free and bound-bound opacities are merged, this code computes both. The distinction between line and continuous opacity becomes blurred, so be careful that you are not including things twice. That goes for Rayleigh scattering off Lyman alpha also.

You can download the codes here.  An example of their application to the SYNTH spectral synthesis program by Nik Piskunov is also provided.  For *fast* implementation you may have to adapt this. For example, in many applications it will be stupid to recompute all the populations for the same model every time you change wavelength.

Views, comments, bug reports and suggested improvements all welcome

History:

19 Nov 2002: A small bug in HLINOP which allowed the wavelength to go negative if the REFLECT switch was turned on was fixed. Thanks to Kjell Eriksson for reporting it.

4 Jun 2003: I fixed the example files so they actually work;-), the old line file referenced the wrong model (one that is on my machine but not the tarred one!). I also have added compiler options for Intel's Fortran compiler under Linux to the makefile.

26 Aug 2003: Two minor bugs fixed which meant the resonance and Stark broadening in the Lyman alpha red wing was incorrect (bugs meant they were not computed at all).

27 April 2005: Small bug in HLINPROF, meant that log(0) occurred if using the self-broadening grids and interrogating a wavelength beyond the range of the profile grids. Now fixed thanks to Ulrike Heiter.

9 Jul 2008: SYNTH code updated so that it warns, rather than crashes, if one tries to feed it a spherical model. Thanks to Nik.

2 Mar 2011: Small bug for case where XHOLT4000 goes to zero fixed. Thanks to Kjell.

2 Mar 2011:  hbop code added.



