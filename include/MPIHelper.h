//Define USE_MPI to use mpi.
//without USE_MPI, mpi is not used.
#ifdef MPI_HELPER_H
#else
#define MPI_HELPER_H

#ifdef USE_MPI
//be careful of, this block should be the first code in main file.
#define OMPI_IMPORTS
#include "mpi.h"
#endif

class MPIHelper{
#define TAG 99
	int id;
	int numProcess;
#ifdef USE_MPI
	MPI_Status status;
	char name[MPI_MAX_PROCESSOR_NAME];
#endif
	public:
		MPIHelper(int &argc,char ** &argv){
#ifdef USE_MPI
			int len;
			MPI_Init(&argc,&argv);
			MPI_Comm_rank(MPI_COMM_WORLD,&id);
			MPI_Comm_size(MPI_COMM_WORLD,&numProcess);
			MPI_Get_processor_name(name,&len);
#else
		id=0;
		numProcess=1;
#endif
		}
		const char *getName()const{
#ifdef USE_MPI
			return name;
#else
			return "localhost";
#endif
		}
		bool isMaster()const{return id==0;}
		int getID()const{return id;}
		int getNumProcesses()const{return numProcess;}
		void send(double *buf,int count,int dest,int tag=TAG){
#ifdef USE_MPI
			MPI_Send(buf,count,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
#endif
		}
		void recv(double *buf,int count,int src,int tag=TAG){
#ifdef USE_MPI
			MPI_Recv(buf,count,MPI_DOUBLE,src,tag,MPI_COMM_WORLD,&status);
#endif
		}
/*
		const MPI_Status *getMPIStatus()const{
			return &status;
		}
*/
		~MPIHelper(){
#ifdef USE_MPI
			MPI_Finalize();
#endif
		}
};
#endif

