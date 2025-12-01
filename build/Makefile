FC = gfortran
FCFLAGS = -Ofast -I../source `nf-config --fflags` -fPIC
LDFLAGS = -L`nc-config --libdir` `nf-config --flibs`

all: model

model : ../source/model.f90 module_parameters.o module_types.o \
        module_physics.o module_output.o
	$(FC) $(FCFLAGS) $^ -o $@ $(LDFLAGS)

module_parameters.o : ../source/module_parameters.f90
	$(FC) $(FCFLAGS) -c $< -o $@

module_types.o : ../source/module_types.F90 module_parameters.o
	$(FC) $(FCFLAGS) -c $< -o $@

module_physics.o : ../source/module_physics.f90 module_types.o
	$(FC) $(FCFLAGS) -c $< -o $@

module_output.o : ../source/module_output.F90 module_physics.o
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	rm -f *.o *.mod model *.nc

