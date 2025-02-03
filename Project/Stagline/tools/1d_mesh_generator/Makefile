# Uncomment the two following lines to compile using gfortran
FC=gfortran
FFLAGS= -fdefault-double-8 -fdefault-real-8 

# Uncomment the two following lines to compile using ifort
#FC=ifort
#FFLAGS = -i4 -r8 

all:
	@echo "Compiling...";
	@$(FC) $(FFLAGS) -o mesh_generator mesh_generator.F90;

clean:
	@echo "Deleting executable file";
	@rm -f  mesh_generator;


