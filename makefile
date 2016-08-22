all:parallel
serial:*
	g++ *.cpp
	./a.out


parallel:*
	mpic++ -DUSE_MPI *.cpp
	mpirun -f machinefile -n 4 ./a.out


clean:
	rm *.txt
