all:
	ifort Shallow2D.f90 -qopenmp -heap-arrays -O3 -o run
single:
	ifort Shallow2D.f90	-heap-arrays -O3 -o run
clean:
	rm Output_* cpu\ time.txt Breach\ Width.plt run
