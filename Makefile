######################################################################
# Makefile for decode...
#
# JLP
# Version 19-07-2008
######################################################################
include $(JLPSRC)/jlp_make.mk
mylib=$(JLPLIB)/jlp
SIMU_LIB=fourn.o kolmo.o jlp_poidev.o bessel_j.o 

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

all: $(SIMU_LIB)

kolmo.o : kolmo.c

fstructure.o : fstructure.c
	cc $(CFLAGS) -c $(INC) $*.c 

fstructure.exe : fstructure.c 
	cc $(CFLAGS) -c $(INC) -DMAIN_PROGRAM $*.c 
	cc $(CFLAGS) $(INC) -o $(EXEC)/$*.exe $*.o jlp_poidev.o \
	$(JLIB) $(MATHLIB) $(MIDLIB) $(F77LIB) -lm
	rm $*.o

s_simu1.exe : s_simu1.c 
	cc $(CFLAGS) -c $(INC) $*.c 
	cc $(CFLAGS) $(INC) -o $(EXEC)/$*.exe $*.o \
	../hrsa/jlp_cover_mask.o fourn.o jlp_poidev.o fstructure.o \
	$(JLIB) $(MATHLIB) $(MIDLIB) $(F77LIB) -lm

bessel_j.o : bessel_j.c

fourn.o : fourn.c

jlp_poidev.o : jlp_poidev.c

simu2.exe : simu2.c jlp_poidev.o 

simu3.exe : simu3.c kolmo.o jlp_poidev.o 

kolmo.exe : kolmo.c 
	cc $(CFLAGS) -c $(INC) -DMAIN_PROGRAM $*.c 
	cc $(CFLAGS) $(INC) -o $(EXEC)/$*.exe \
	$(SIMU_LIB) $(JLIB) $(MIDLIB) $(F77LIB) -lm

create_object.exe : create_object.c
	cc $(CFLAGS) $(INC) $*.c -o $(EXEC)/$*.exe \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

object_binary.exe : object_binary.c
	cc $(CFLAGS) $(INC) $*.c -o $(EXEC)/$*.exe \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

s_simu.exe: s_simu.for fstructure.c $(SIMU_LIB) 
	runs esoext1 -f $*.for
	$(F77) $(FFLAGS) -c $*.f
	$(F77) -g -o $(EXEC)/$*.exe $*.o $(SIMU_LIB) \
	fstructure.o ../hrsa/jlp_cover_mask.o \
	$(JLIB) $(MATHLIB) $(MIDLIB) $(F77LIB) -lm

clear:
	rm -f *.o

