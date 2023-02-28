CC=g++
CFLAGS=-I
DEPS= SpatialLinearSystemConstruction.h QuasiEuler.h Geometry.h

%.o: %.cpp SpatialLinearSystemConsruction.h QuasiEuler.h Geometry.h
	$(CC) -c -o $@ $< $(CFLAGS)

hw2: hw2.o
	$(CC) -o hw2 hw2.o 
