#
#
## Commandes shell
RM = rm -f

#
#
## Compilateur
FC = gfortran
LD = $(FC)

#
#
## Options compilateur
#FFLAGS  = -march=native -O3 -ffree-form -fbacktrace -g
FFLAGS  = -O2 #-ffixed-form -fbacktrace -g -fcheck=all -Wall -ffpe-trap=zero,overflow,invalid
LDFLAGS = -O3
#LDFLAGS = 

#
#
## Executable
EXEC = exe

#
###
# objects Fortran
FOBJS = \
	numerics.o \
	algebre.o \
	vit_et_topo.o \
	advection.o \
	cond.o \
	ecriture.o \
	flux.o \
	calcul.o \
	correction.o \
	main.o

FMODS = \
	numerics.mod \
	mod_algebre.mod \
	mod_vit_top.mod \
	mod_advection.mod \
	mod_cond.mod \
	mod_dat.mod  \
	mod_flux.mod \
	mod_calcul.mod \
	mod_correction.mod 
# sources Fortran
FSRCS = \
	numerics.f90 \
	algebre.f90 \
	vit_et_topo.f90 \
	advection.f90 \
	cond.f90 \
	ecriture.f90 \
	flux.f90 \
	calcul.f90 \
	correction.f90 \
	main.f90 
#
%.o: %.f90
	$(FC) $(FFLAGS) $(FSRCS) -c $(FOBJS)
#
#
## Dependances compilation
all: $(FOBJS)
	$(LD) $(LDFLAGS) $(FOBJS) -o $(EXEC)

#
#
# nettoyage
clean:
	   $(RM) $(FOBJS) $(FMODS)

#
numerics.o: numerics.f90
algebre.o: algebre.f90 
vit_et_topo.o: vit_et_topo.f90 numerics.o
avection.o: advection.f90 vit_et_topo.o numerics.o
cond.o: cond.f90 numerics.o vit_et_topo.o
ecriture.o: ecriture.f90 numerics.o cond.o
flux.o: flux.f90 cond.o numerics.o vit_et_topo.o
calcul.o: calcul.f90 numerics.o cond.o flux.o vit_et_topo.o
correction.o: correction.f90 numerics.o
main.o: main.f90 flux.o cond.o advection.o ecriture.o correction.o numerics.o 
