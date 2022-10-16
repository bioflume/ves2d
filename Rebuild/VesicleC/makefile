OBJS = GMRES.o DatOps.o Misc.o Contour.o Vesicle.o StokesStructs.o Stokes.o Boundary.o TestTop.o VesXML.o
CC = gcc
DEBUG = -g -pg
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)
DFLAGS = -L /usr/lib64/atlas/ 
REFS = -llapack -lm -lfftw3 -lxml2 -latlas -lcblas
INC = -I /usr/include/libxml2 -I /usr/include/atlas -I /usr/lib64/atlas/

Ves2DC: $(OBJS)
	$(CC) $(LFLAGS) $(DFLAGS) $(OBJS) -o Ves2DC $(REFS)

GMRES.o: GMRES.h GMRES.c
	$(CC) $(CFLAGS) GMRES.c $(INC)

DatOps.o: DatOps.h DatOps.c
	$(CC) $(CFLAGS) DatOps.c $(INC)

Misc.o: Misc.h Misc.c
	$(CC) $(CFLAGS) Misc.c $(INC)

Contour.o: Contour.h Contour.c
	$(CC) $(CFLAGS) Contour.c $(INC)

Vesicle.o: Vesicle.h Vesicle.c
	$(CC) $(CFLAGS) Vesicle.c $(INC)

StokesStucts.o: StokesStructs.h StokesStructs.c
	$(CC) $(CFLAGS) StokesStructs.c $(INC)

Stokes.o: Stokes.h Stokes.c 
	$(CC) $(CFLAGS) Stokes.c $(INC)

Boundary.o: Boundary.h Boundary.c
	$(CC) $(CFLAGS) Boundary.c $(INC)

TestTop.o: TestTop.c
	$(CC) $(CFLAGS) TestTop.c $(INC)

VesXML.o: VesXML.h VesXML.c
	$(CC) $(CFLAGS) VesXML.c $(INC)

clean:
	\rm *.o