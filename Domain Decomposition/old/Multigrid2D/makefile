f90=ifort
flags= -qopenmp -heap-arrays -O3 -mkl=sequential 
src= ./src/main.f90 
lib=./lib/modules.f90
program=RUN

.PHONY: clean

all: modules.o
	$(f90) $(src) modules.o $(flags) -o $(program)
modules.o:
	$(f90) -c $(lib) $(flags) -o modules.o
clean:
	rm -rf *.mod modules.o $(program) ./out/*.vtk ./out/*.plt
