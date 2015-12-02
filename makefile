all:parallel
serial:*
	g++ -g *.cpp
	./a.out


parallel:*
	mpic++ -DUSE_MPI *.cpp
	mpirun -n 5 ./a.out
