.SUFFIXES: .o .f

### adjust to site requirements ########################################

NPROCS  = 1
PRGDIR  = .

COMPILER = INTEL
#COMPILER = PGI
#COMPILER = G95

#if ($COMPILER == "IFORT") THEN
#	MYLIBS = /data/M3HOME/lib/

MPIPATH     = /sw/spack-levante/openmpi-4.1.2-yfwe6t/lib
NETCDFPATH  = /sw/spack-levante/netcdf-c-4.8.1-2k3cmu/lib
NETCDFFPATH = /sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g/lib
PNETCDFPATH = /sw/spack-levante/parallel-netcdf-1.12.2-mc24h4/lib
NETCDFINC   = /sw/spack-levante/netcdf-c-4.8.1-2k3cmu/include
NETCDFFINC  = /sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g/include
IOAPIPATH   = /work/gg0305/g260102/ioapi/ioapi-3.2/Linux2_x86_64ifortmpi/
IOAPIINC    = /work/gg0305/g260102/ioapi/ioapi-3.2/ioapi/fixed_src/

	FC = mpif90
	FLINKER = mpif90

#DEBUGGING FLAGS
#	FFLAGS = -g -cpp -march=core-avx2 -O3 -fp-model source -assume byterecl -132 -qopenmp -check all -traceback 

#PRODUCTION FLAGS
# 	FFLAGS = -g -cpp -fp-model source -assume byterecl -132 -qopenmp -check all -traceback

LFLAGS  = -L$(MPIPATH) -L/usr/lib64 -L${NETCDFPATH} -L${NETCDFFPATH} -L${PNETCDFPATH} -L${IOAPIPATH} -qopenmp# -check all #-L${IFPATH}
LIBS    = -lnetcdf -lnetcdff -lpnetcdf -lioapi #-lmpi #-libimf -libsvml -libintlc -libifport -libifcore
IFLAGS  = -I$(NETCDFINC) -I${NETCDFFINC} -I${IOAPIINC} #-I$(MPIINC)

PROGRAM = polpar_main 
PROGEXT = .x

PROGHOME = .
#SRCDIR  = $(PROGHOME)/SOURCE
SRCDIR  = $(PROGHOME)

OBJS =  cparam.o		   \
	cutil.o			  \
	cout.o			 \
	cin.o                   \
	cdata.o		       \
	polpar_main.o         \
        polpar_ncrw.o        \
        polpar_light.o      \
        polpar_boun.o      \
        partpart.o        \
	chem.o		 \
	airsea.o	\
	icechem.o      \
	bioacc_mehg.o \
        bioacc_hg.o  \
	sediment.o  \
	strom.o

EXEPROG = $(PRGDIR)/$(PROGRAM)$(PROGEXT)

$(EXEPROG) : $(OBJS)
	$(FLINKER) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

#include /usr/local/include/netcdf.inc

run : $(EXEPROG)
	-git describe --tags > which_tag.txt
	-git rev-list --tags --max-count=1 > which_commit.txt
	mpirun -np $(NPROCS) -machinefile MACHINES.USE  $(EXEPROG) 

clean : 
	-rm $(EXEPROG) *.mod *.o

#depend :
#	cd $(SRCDIR); \
#	grep -i '^[ ]*use ' *.f | sed -f dep_script | sort -u > use.dep;

%.o: $(SRCDIR)/%.f
	$(FC) $(FFLAGS) $(IFLAGS) -o $@ -c $<

