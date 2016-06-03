PROG =	driver

SRCS =	array_lib.f90 atmos_lib.f90 driver.f90 dsd.f90 dsd_melt.f90 \
	format_input.f90 gases.f90 load_hydrometeor_classes.f90 \
	load_mie_table.f90 math_lib.f90 mrgrnk.f90 optical_melt.f90 \
	optical_sphere.f90 optics_lib.f90 radar_simulator.f90 \
	radar_simulator_types.f90 zeff.f90

OBJS =	array_lib.o atmos_lib.o driver.o dsd.o dsd_melt.o format_input.o \
	gases.o load_hydrometeor_classes.o load_mie_table.o math_lib.o \
	mrgrnk.o optical_melt.o optical_sphere.o optics_lib.o \
	radar_simulator.o radar_simulator_types.o zeff.o

LIBS =	

F90 = ifort
F90FLAGS = -O
LDFLAGS = 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

array_lib.o: mrgrnk.o
driver.o: array_lib.o atmos_lib.o format_input.o radar_simulator_types.o
dsd.o: array_lib.o math_lib.o
format_input.o: array_lib.o
load_hydrometeor_classes.o: radar_simulator_types.o
load_mie_table.o: array_lib.o radar_simulator_types.o
math_lib.o: array_lib.o mrgrnk.o
optical_sphere.o: math_lib.o optics_lib.o
radar_simulator.o: array_lib.o math_lib.o optics_lib.o \
	radar_simulator_types.o
zeff.o: math_lib.o
