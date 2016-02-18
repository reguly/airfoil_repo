all: airfoil_seq
airfoil_seq: airfoil_seq.cpp Makefile
	gcc airfoil_seq.cpp -o airfoil_seq -O3 -fopenmp
