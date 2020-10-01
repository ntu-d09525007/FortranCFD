all:mod.o	
	ifort MAIN.F90 mod.o -qopenmp -heap-arrays -O3 -o run
mod.o:
	ifort -c  DATA_MOD.F90 -qopenmp -heap-arrays -o mod.o
clean:
	rm -rf *.mod *.o *.PLT *.vtk *.TXT *.DAT run nohup.*
