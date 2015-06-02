a.out: codicino.o print_matrix.o delta.o print_array.o jacobi.o fock.o combinatorics.o
	gfortran -o a.out -L${ACML}/lib -lacml codicino.o print_matrix.o delta.o print_array.o jacobi.o fock.o combinatorics.o

codicino.o: codicino.f90
	gfortran -O3 -c -L${ACML}/lib -lacml codicino.f90

print_matrix.o: print_matrix.f90
	gfortran -O3 -c print_matrix.f90

print_array.o: print_array.f90
	gfortran -O3 -c print_array.f90

jacobi.o: jacobi.f90
	gfortran -O3 -c jacobi.f90

delta.o: delta.f90
	gfortran -O3 -c delta.f90

combinatorics.o: combinatorics.f90
	gfortran -O3 -c combinatorics.f90

fock.o: fock.f90
	gfortran -O3 -c fock.f90

clean:
	rm -rf *.o
