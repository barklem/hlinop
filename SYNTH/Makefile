# HP compiler options
CFLAG = +O3 +Onolimit +U77
LFLAG = +O3 +Onolimit +U77
#############################
# Sun options
#CFLAG = -O3
#LFLAG = -O3
FF = f90
#############################
# LINUX Intel compiler options
#CFLAG = -O3 -w95 -Vaxlib
#LFLAG = -O3 -Vaxlib
#FF = ifc
#############################
# Mac OSX, intel fortran  options
CFLAG = -O3 
LFLAG = -O3 
FF = ifort
#############################
# Directory for executable
T=.
VERSION = iau210
# Use these for big-endian machines like Sun, HP
#HFILE     = $(PWD)/stehle.dat.SUN
#HSELFFILE = $(PWD)/bpo_self.grid.SUN
# Use these for little-endian machines like Intel PC, DEC, IBM
HFILE     = $(PWD)/stehle.dat.DEC
HSELFFILE = $(PWD)/bpo_self.grid.DEC
# VCS file is ascii
HVCSFILE  = $(PWD)/vcsbalmer.dat
EXECNAME = synth_hlinop


$(T)/$(EXECNAME): synth.o transfs.o hlinop.o hlinprof.o contins.o \
                  Makefile
	$(FF) $(LFLAG) -o $(T)/$(EXECNAME) synth.o transfs.o hlinprof.o \
                  hlinop.o contins.o $(LIBS)
#
synth.o: synth.f COMMONS.SYN SIZES.SYN Makefile
	echo "      PARAMETER (VERSION='$(VERSION)')" >  VERSION
	echo "      PARAMETER (VERDATE='`date +%d-%B-%y`')" >> VERSION
	$(FF) -c $(CFLAG) synth.f
	rm VERSION
transfs.o: transfs.f COMMONS.SYN SIZES.SYN Makefile
	$(FF) -c $(CFLAG) transfs.f
# the below ${VAR:start:finish} works for bash, unsure about other shells	
hlinprof.o: hlinprof.f $(HFILE) $(HSELFFILE) $(HVCSFILE) Makefile
	echo "      PARAMETER (HFILE="            >  DATA.FILES
	echo "     ;'$(HFILE)')"                 >> DATA.FILES
	echo "      PARAMETER (HSELFFILE="        >> DATA.FILES
	echo "     ;'$(HSELFFILE)')"             >> DATA.FILES
	echo "      PARAMETER (HVCSFILE="         >> DATA.FILES
	echo "     ;'$(HVCSFILE)')"              >> DATA.FILES
	$(FF) -c $(CFLAG) -132 hlinprof.f
	rm DATA.FILES
hlinop.o: hlinop.f Makefile
	$(FF) -c $(CFLAG) hlinop.f
contins.o: contins.f COMMONS.SYN SIZES.SYN Makefile
	$(FF) -c $(CFLAG) contins.f
