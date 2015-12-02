#ifdef MPI_HELPER_H
#else
#define MPI_HELPER_H

#ifdef USE_MPI
//be careful of, this block should be the first code in main file.
#define OMPI_IMPORTS
#include "mpi.h"

class MPIHelper{
#define TAG 99
	int id;
	int numProcess;
	MPI_Status status;
	public:
		MPIHelper(int &argc,char ** &argv){
			MPI_Init(&argc,&argv);
			MPI_Comm_rank(MPI_COMM_WORLD,&id);
			MPI_Comm_size(MPI_COMM_WORLD,&numProcess);
		}
		bool isMaster()const{return id==0;}
		int getID()const{return id;}
		int getNumProcesses()const{return numProcess;}
		void send(const double *buf,int count,int dest,int tag=TAG){
			MPI_Send(buf,count,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
		}
		void recv(double *buf,int count,int src,int tag=TAG){
			MPI_Recv(buf,count,MPI_DOUBLE,src,tag,MPI_COMM_WORLD,&status);
		}
		const MPI_Status *getMPIStatus()const{
			return &status;
		}
		~MPIHelper(){
			MPI_Finalize();
		}
};
#endif
#endif

