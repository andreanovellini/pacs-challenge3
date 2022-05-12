CXX=g++
CC=$(CXX)
CXXFLAGS=-Wall -std=c++20

all: main_matrix

main_matrix: main_matrix.o mmio.o

main_matrix.o: main_matrix.cpp Matrix.hpp mmio.c 

clean: 
	$(RM) *.o

distclean: clean
	$(RM) main_matrix

run:
	./main_matrix lnsp_131.mtx
