f90=mpiifort
flags= -qopenmp -heap-arrays -O3 -mkl=sequential -module $(moddir)
program=RUN

src= ./src/main.f90 
lib=./lib/modules.f90

moddir=./mod/
plotdir=./lib/plot/

packages=modules.o ir_precision.o lib_base64.o lib_pack_data.o lib_vtk_io.o

.PHONY: clean
.PHONY: clear

all:$(packages)
	$(f90) $(src) $(packages) -I$(moddir) $(flags) -o $(program) 
modules.o:
	$(f90) -c $(lib) $(flags) -o modules.o
###########################################
lib_vtk_io.o: lib_pack_data.o lib_base64.o ir_precision.o
	$(f90) -cpp -c $(plotdir)$*.f90 $^ -O3 -module $(moddir) 
lib_base64.o: ir_precision.o lib_pack_data.o
	$(f90) -cpp -c $(plotdir)$*.f90 $^ -O3 -module $(moddir)
lib_pack_data.o: ir_precision.o
	$(f90) -cpp -c $(plotdir)$*.f90 $^ -O3 -module $(moddir)
ir_precision.o:
	$(f90) -cpp -c $(plotdir)$*.f90 -O3 -module $(moddir)
##########################################
clean:
	rm -rf ./mod/* $(program) $(packages)
clear:
	rm -rf ./mod/* $(program) ./out/* $(packages) 
