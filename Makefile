PROG =	radsim

SRCS =	array_lib.f90 dsd.f90 gases.f90 math_lib.f90 \
	m_mrgrnk.f90 optical_sphere.f90 optics_lib.f90 \
	zeff.f90

OBJS =	array_lib.o dsd.o gases.o math_lib.o m_mrgrnk.o \
	optical_sphere.o optics_lib.o \
	zeff.o

LIB_SRC = dsd.f90 math_lib.f90 optical_sphere.f90 zeff.f90 gases.f90

LIBS =	

F90 = gfortran
F90FLAGS = -fPIC -O
LDFLAGS = 

F2PY = f2py3
F2PY_FLAGS = --quiet

# See https://stackoverflow.com/questions/43618725/makefile-for-f2py-in-python3
EXT_SUFFIX := $(shell python3-config --extension-suffix)

all: $(PROG)

$(PROG): $(OBJS) $(PROG).pyf $(PROG).so

$(PROG).pyf: $(OBJS)
	$(F2PY) $(F2PY_FLAGS) -m $(PROG) -h $(PROG).pyf $(LIB_SRC)

$(PROG).so: $(OBJS) $(PROG).pyf
	$(F2PY) $(F2PY_FLAGS) -c $(PROG).pyf $(SRCS)
	mv $(PROG).*.so $(PROG).so

clean:
	rm -f $(PROG).so $(OBJS) *.mod quickbeam.pyc

clean-all:
	rm -f $(PROG).so $(PROG).pyf $(OBJS) *.mod quickbeam.pyc

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

array_lib.o: m_mrgrnk.o
dsd.o: array_lib.o math_lib.o
math_lib.o: array_lib.o m_mrgrnk.o
optical_sphere.o: math_lib.o optics_lib.o
zeff.o: math_lib.o
