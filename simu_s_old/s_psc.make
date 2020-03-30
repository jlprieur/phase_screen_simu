##############################################################################
# Makefile for s_psc and s_pscorg 
# JLP
# Version 06-11-91
##############################################################################
OBJ=  s_psc.o s_pscorg.o
midir=$(JLPLIB)/midas
mylib=$(JLPLIB)/jlp
smgolib=$(JLPLIB)/smgo
#JLIB=$(mylib)/jlpacc.a $(mylib)/jlputil.a \
#     $(mylib)/newplot0.a $(mylib)/jlp_splot.a 
JLIB=$(mylib)/jlpacc.a $(mylib)/jlputil.a
LIBC= /lib/libc.a
XLIB= -lX11
#XLIB= -lXaw -lXmu -lXt -lX11
F77LIB=/usr/lib/libF77.a /usr/lib/libI77.a /usr/lib/libU77.a 
#MIDLIB=$(nagdir)/tbllib.a $(midir)/stlib.a $(nagdir)/mathlib.a \
#       $(midir)/oslib.a
MIDLIB=$(midir)/stlib.a $(midir)/udiolib.a $(midir)/oslib.a
#MATHLIB=/usr1/midas/91MAY/lib/mathlib.a
#MATHLIB=$(JLPLIB)/math/mymath.a
#MATHLIB= fft_nag.o nagfft.o jlp_short.o
MATHLIB= 
F77=f77
INC=.
FFLAGS=-g
#

.SUFFIXES:
.SUFFIXES:  .o .for .exe $(SUFFIXES) 

.for.exe:
	runs esoext1 -I . -f $*.for
	$(F77) -c $(FFLAGS) $*.f
	rm $*.f
	$(F77) $(FFLAGS) -o $(EXEC)/$*.exe $*.o \
	$(JLIB) $(MIDLIB) $(MATHLIB) $(XLIB) $(LIBC)

.for.o:
	runs esoext1 -I . -f $*.for
	$(F77) -c $(FFLAGS) $*.f
	rm $*.f
	$(F77) $(FFLAGS) -o $(EXEC)/$*.exe $*.o \
	$(JLIB) $(MIDLIB) $(MATHLIB) $(XLIB) $(LIBC)

all: $(OBJ) 
	@sleep 1

s_psc.exe: s_psc.for

s_pscorg.exe: s_pscorg.for

clean:
	rm -f $(OBJ)
