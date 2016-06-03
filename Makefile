PROG =	radsim

SRCS =	array_lib.f90 dsd.f90 dsd_melt.f90 gases.f90 math_lib.f90 \
	mrgrnk.f90 optical_melt.f90 optical_sphere.f90 optics_lib.f90 \
	zeff.f90

OBJS =	array_lib.o dsd.o dsd_melt.o gases.o math_lib.o mrgrnk.o \
	optical_melt.o optical_sphere.o optics_lib.o \
	zeff.o

LIBS =	

F90 = gfortran
F90FLAGS = -fPIC -O
LDFLAGS = 

F2PY = f2py
F2PY_FLAGS = --quiet

all: $(PROG)

$(PROG): $(OBJS) $(PROG).pyf $(PROG).so

$(PROG).pyf: $(OBJS)
	$(F2PY) $(F2PY_FLAGS) -m $(PROG) -h $(PROG).pyf *.f90

$(PROG).so: $(OBJS) $(PROG).pyf
	$(F2PY) $(F2PY_FLAGS) -c $(PROG).pyf *.f90

clean:
	rm -f $(PROG).so $(PROG).pyf $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

array_lib.o: mrgrnk.o
dsd.o: array_lib.o math_lib.o
math_lib.o: array_lib.o mrgrnk.o
optical_sphere.o: math_lib.o optics_lib.o
zeff.o: math_lib.o
