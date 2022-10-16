MEX=/opt/apps/matlab/r2017b/bin/mex
#location of mex compiler on ronaldo
FORTRAN = ifort
FFLAGS = -c -fpic -qopenmp -O3
MEXLIBDIR = "/opt/apps/intel/17.0u4/parallel_studio_xe_2017/"
#for ronaldo
MEXLIB = -lifcore -lirc -lsvml -limf -lgfortran -liomp5

OBJS=\
bhcommon_oldtree.o\
bhfmm2dpart.o\
bhfmm2drouts.o\
d2tstrcr_omp.o\
d2mtreeplot.o\
l2dterms.o\
laprouts2d.o\
bhrouts.o\
prini.o\
hkrand.o\
dlaran.o

all: bhcommon_oldtree bhfmm2dpart bhfmm2drouts d2tstrcr_omp d2mtreeplot l2dterms laprouts2d bhrouts prini hkrand dlaran mex

mex: stokesDLPnewfmm.F $(OBJS)
	$(MEX) -O stokesDLPnewfmm.F $(OBJS) -L$(MEXLIBDIR) $(MEXLIB)

bhfmm2dpart_dr: bhfmm2dpart_dr.f
	$(FORTRAN) $(FFLAGS) bhfmm2dpart_dr.f

bhcommon_oldtree: bhcommon_oldtree.f
	$(FORTRAN) $(FFLAGS) bhcommon_oldtree.f

bhfmm2dpart: bhfmm2dpart.f
	$(FORTRAN) $(FFLAGS) bhfmm2dpart.f

bhfmm2drouts: bhfmm2drouts.f
	$(FORTRAN) $(FFLAGS) bhfmm2drouts.f

d2tstrcr_omp: d2tstrcr_omp.f
	$(FORTRAN) $(FFLAGS) d2tstrcr_omp.f

d2mtreeplot: d2mtreeplot.f
	$(FORTRAN) $(FFLAGS) d2mtreeplot.f

l2dterms: l2dterms.f
	$(FORTRAN) $(FFLAGS) l2dterms.f

laprouts2d: laprouts2d.f
	$(FORTRAN) $(FFLAGS) laprouts2d.f

bhrouts: bhrouts.f
	$(FORTRAN) $(FFLAGS) bhrouts.f

prini: prini.f
	$(FORTRAN) $(FFLAGS) prini.f

hkrand: hkrand.f
	$(FORTRAN) $(FFLAGS) hkrand.f

dlaran: dlaran.f
	$(FORTRAN) $(FFLAGS) dlaran.f

second-r8: second-r8.f
	$(FORTRAN) $(FFLAGS) second-r8.f

clean:
	rm *.o
	cp *.mex* ../src
	rm *.mex*
