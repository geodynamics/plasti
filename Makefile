
## compilers
FORT  = gfortran
F90FLAGS = -O3

# g95-specific flags
#FORT  = g95
#F90FLAGS = -O3 -fsloppy-char

## libraries
LIBS = -llapack -lblas

all: plasti_oly meshg_oly

#####
#####  PLASTI
#####
## object files to link (includes any header and module files)
PLAS_OBJS = SRC/plasti_oly.o SRC/thermal_oly.o
## Link all files into main program
plasti_oly: $(PLAS_OBJS)
	$(FORT) $(LINKFLAGS) $(PLAS_OBJS) -o plasti_oly $(LIBS)
## compile object files
SRC/plasti_oly.o: SRC/plasti_oly.f 
	$(FORT) $(F90FLAGS) -c SRC/plasti_oly.f -o SRC/plasti_oly.o
SRC/thermal_oly.o: SRC/thermal_oly.f
	$(FORT) $(F90FLAGS) -c SRC/thermal_oly.f -o SRC/thermal_oly.o

## clean
clean: 
	rm -f $(PLAS_OBJS) *.mod $(MESH_OBJS) $(PLAS2DX_OBJS) plasti_oly meshg_oly

#####
##### MESHG
#####
## object files to link
MESH_OBJS = SRC/meshg_oly.o SRC/erfc.o
## Link files into main program
meshg_oly: $(MESH_OBJS)
	$(FORT) $(LINKFLAGS) $(MESH_OBJS) -o meshg_oly
## compile object files
SRC/meshg_oly.o: SRC/meshg_oly.f
	$(FORT) $(F90FLAGS) -c -o SRC/meshg_oly.o SRC/meshg_oly.f
SRC/erfc.o: SRC/erfc.f
	$(FORT) $(F90FLAGS) -c -o SRC/erfc.o SRC/erfc.f

