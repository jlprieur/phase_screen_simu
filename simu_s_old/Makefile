######################################################################
# Makefile for s_ programs 
#
# JLP
# Version 19-07-2008
######################################################################
include $(JLPSRC)/jlp_make.mk
mylib=$(JLPLIB)/jlp

.SUFFIXES: 
.SUFFIXES: .o .c .for .exe $(SUFFIXES)

.for.o:
	runs esoext1 -f $*.for
	$(F77) $(FFLAGS) -c $*.f

.for.exe:
	runs esoext1 -f $*.for
	$(F77) $(FFLAGS) -c $*.f
	$(F77) -g -o $(EXEC)/$*.exe $*.o \
	$(JLIB) $(MATHLIB) $(MIDLIB) $(F77LIB) -lm
	rm $*.o

.c.o:
	cc -c $(CFLAGS) $*.c

.o.exe:
	cc -o $(EXEC)/$*.exe $*.o $(SIMU_LIB) \
	$(JLIB) $(MIDLIB) $(F77LIB) $(XLIB) -lm

.c.exe:
	cc -c $(CFLAGS) $*.c
	cc -o $(EXEC)/$*.exe $*.o $(SIMU_LIB) \
	$(JLIB) $(MIDLIB) $(F77LIB) $(XLIB) -lm

all: $(OUT)

s_screen.exe: s_screen.for

s_pscorg.exe: s_pscorg.for $(SIMU_LIB) 
	runs esoext1 -f $*.for
	$(F77) $(FFLAGS) -c $*.f
	$(F77) -g -o $(EXEC)/$*.exe $*.o ../hrsa/jlp_cover_mask.o \
	$(JLIB) $(MATHLIB) $(MIDLIB) $(F77LIB) -lm

s_fourn1.o: s_fourn1.for

clear:
	rm -f *.o
