CC = gfortran #compiler
CFLAGS = -O3 -ftree-vectorize # runtime
CCPACK = -llapack -lblas # sudo apt-get install liblapack-dev

## executables
von_karman_swirling_flow.x: von_karman_swirling_flow.o module_subroutines.o module_algebra.o
	$(CC) $(CFLAGS) von_karman_swirling_flow.o module_subroutines.o module_algebra.o -o von_karman_swirling_flow.x $(CCPACK)

## object files
von_karman_swirling_flow.o: von_karman_swirling_flow.f90
	$(CC) $(CFLAGS) -c von_karman_swirling_flow.f90 -o von_karman_swirling_flow.o
module_subroutines.o: module_subroutines.f90
	$(CC) $(CFLAGS) -c module_subroutines.f90 -o module_subroutines.o
module_algebra.o: module_algebra.f90
	$(CC) $(CFLAGS) -c module_algebra.f90 -o module_algebra.o
