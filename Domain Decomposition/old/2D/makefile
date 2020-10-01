all:
	cd ./lib && ifort ../src/main.f90 modules.o -qopenmp -O3 -o ../run
mod:
	cd ./lib && ifort -c modules.f90 -qopenmp
clean:
	rm -rf ./run ./lib/*.mod ./lib/*.o 
clear:
	rm -rf ./run ./lib/*.mod ./lib/*.o ./out/*.dat ./out/*.plt ./out/*.vtk
cout:
	rm -rf ./out/*.dat ./out/*.plt ./out/*.vtk

