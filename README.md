# hlinop

A set of codes for computing H line profiles and opacities has been developed by myself and Nik Piskunov, and was first presented at IAU 210, Modelling Stellar Atmospheres, in Uppsala June 2002. There is clearly a need for freely available code to compute these opacities. Codes such as HLINOP and BALMER9 have been provided by Kurucz and Peterson (http://cfaku5.harvard.edu/) in the 1970's based on the broadening theories of the time. Since then, understanding of H line broadening has improved significantly and so there is a need to update these codes. Three codes have been produced, namely:

HLINPROF: For detailed, accurate calculation of lower Balmer line profiles, suitable for detailed analysis of Balmer lines.  This code depends on binary files containing tables and are thus provided in little- and big-endian versions.  There is also a version (_nofiles) which attempts to avoid the need for external files, but is not as well tested as the default version.

HLINOP: For approximate, quick calculation of any line of neutral H, suitable for model atmosphere calculations.  Based on the original Kurucz and Peterson version.

HBOP:A wrapper around hlinop.f to implement the occupation probability formalism of Daeppen, Anderson and Milhalas and thus can account for the merging of lines into the continuum.  As bound-free and bound-bound opacities are merged, this code computes both. The distinction between line and continuous opacity becomes blurred, so be careful that you are not including things twice. That goes for Rayleigh scattering off Lyman alpha also.

You can download the codes here.  An example of their application to the SYNTH spectral synthesis program by Nik Piskunov is also provided.  For *fast* implementation you may have to adapt this. For example, in many applications it will be stupid to recompute all the populations for the same model every time you change wavelength.


Views, comments, bug reports and suggested improvements all welcome.

