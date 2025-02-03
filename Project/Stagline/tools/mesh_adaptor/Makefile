# Uncomment the two following lines to compile using gfortran
FC=gfortran
FFLAGS= -fdefault-double-8 -fdefault-real-8 

# Uncomment the two following lines to compile using ifort
#FC=ifort
#FFLAGS = -i4 -r8 

all:
	@$(FC) $(FFLAGS) -o adapt adapt.f90;

clean:
	@rm -f  adapt *.mod;


