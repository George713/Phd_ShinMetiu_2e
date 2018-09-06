FC  		= f95
FFLAGS		= -fopenmp -lfftw3 -funroll-loops -O3 -mcmodel=medium
LDFLAGS		= -L/home/ki96zib/apps/fftw/lib
# LDFLAGS		=
SRC			= ./src/
BIN			= ./bin/

OBJS		=	$(BIN)modules.o		$(BIN)subroutines.o\
				$(BIN)input.o		$(BIN)potential.o\
				$(BIN)adiabatic.o	$(BIN)nuclear_wf.o\
				$(BIN)prop1d.o		$(BIN)prop3d.o\
				$(BIN)main.o

control:	$(OBJS)
		$(FC) $(BIN)*.o ${LDFLAGS} ${FFLAGS} -o control

$(BIN)modules.o:	 $(SRC)modules.f90
	${FC} $(SRC)modules.f90 ${FFLAGS} -c -o $@ -J$(BIN)

$(BIN)subroutines.o:	 $(SRC)subroutines.f90
	${FC} $(SRC)subroutines.f90 ${FFLAGS} -c -o $@ -J$(BIN)

$(BIN)input.o:	 $(SRC)input.f90
	${FC} $(SRC)input.f90 ${FFLAGS} -c -o $@ -J$(BIN)

$(BIN)potential.o:	$(SRC)potential.f90
	${FC} $(SRC)potential.f90 ${FFLAGS} -c -o $@ -J$(BIN)	

$(BIN)adiabatic.o:	$(SRC)adiabatic.f90
	${FC} $(SRC)adiabatic.f90 ${FFLAGS} -c -o $@ -J$(BIN)

$(BIN)nuclear_wf.o:	$(SRC)nuclear_wf.f90
	${FC} $(SRC)nuclear_wf.f90 ${FFLAGS} -c -o $@ -J$(BIN)

$(BIN)prop1d.o:	$(SRC)prop1d.f90
	${FC} $(SRC)prop1d.f90 ${FFLAGS} -c -o $@ -J$(BIN)

$(BIN)prop3d.o:	$(SRC)prop3d.f90
	${FC} $(SRC)prop3d.f90 ${FFLAGS} -c -o $@ -J$(BIN)

$(BIN)main.o:	 $(SRC)main.f90
	${FC} $(SRC)main.f90 ${FFLAGS} -c -o $@ -J$(BIN)

clean:
	rm -f $(BIN)*.o 
	rm -f $(BIN)*.mod
	rm -f ./out/*
	rm -f control
	rm -f time*
	rm -f slurm*
	rm -f std.out
