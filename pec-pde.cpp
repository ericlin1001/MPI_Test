#include "include/template.h"
#include "include/BasicDE.h"
#include "pecFunction.h"
#include<iostream>
#include<fstream>
using namespace std;
DefFunction(MyF1,-100,100,0)
	double res=0.0;
	for(int i=0;i<size;i++){
		double x=xs[i];
		res+=x*x;
	}
	return res;
EndDef

class ParallelDE:public EA
{
	private:
		//about function:f
		Function *f;
		vector<vector<double> >range;
		int numDim;
		//algorithm related parameters.
		int PopSize;
		double F,CR;
		//
		vector<vector<double> >x;//x,trail x.
		vector<vector<double> >tmpX;
		vector<double>fx;
		vector<double>tmpFx;
		vector<double>tx;
		//
		int bestI;
		MPIHelper*mpi;
	private:
		void updateX(){
			//main process
			vector<vector<double> >txs;
			vector<double> ftxs;
			txs.resize(PopSize);

			RandomPermutation perm(PopSize);
			for(int i=0;i<PopSize;i++){
				perm.generate();
				int a=bestI; int b=perm.next(); int c=perm.next();
				int randDim=rand()%numDim;
				for(int j=0;j<numDim;j++){
					if(j==randDim||drand()<CR){
						tx[j]=x[a][j]+F*(x[b][j]-x[c][j]);
						if(tx[j]<range[j][0] || tx[j]>range[j][1]){
							tx[j]=drand(range[j][0],range[j][1]);
						}
					}else{
						tx[j]=x[i][j];
					}
				}
				txs[i]=tx;
			}
			evaluatePopulation(txs,ftxs);
			for(int i=0;i<PopSize;i++){
				vector<double>&tx=txs[i];
				double &ftx=ftxs[i];
				if(ftx<fx[i]){
					x[i]=tx;
					fx[i]=ftx;
					if(ftx<fx[bestI]){
						bestI=i;
					}
				}
			}
		}
	private:
		double getBestFx()const{
			return fx[bestI];
		}
		void update(int maxGeneration){
#define SaveData s.add(getBestFx());
			Save s(f->getName(),"Generation","F");
			SaveData;
			for(int g=1;g<=maxGeneration;g++){
				updateX();
				SaveData;
				//Tic::tic("update one X");
			}
		}
	public:
		ParallelDE(MPIHelper *h):mpi(h){
		}
		void initParam(SearchParam *param){
			PopSize=param->getInt("PopSize");
			F=param->getDouble("F");
			CR=param->getDouble("CR");
			param->getBiVector("Range",range);
			setName(param->getString("Name"));
		}
		void evaluatePopulation(vector<vector<double> >&xs,vector<double>&fx){
			fx.resize(xs.size());
			int numSlaves=mpi->getNumProcesses()-1;
			if(mpi->isMaster()){
				ASSERT(numSlaves>0);
				for(int i=0;i<xs.size();i++){
					mpi->send(&xs[i][0],xs[i].size(),i%numSlaves+1);
					//fx[i]=f->evaluate(&xs[i][0]);
				}
				for(int i=0;i<xs.size();i++){
					mpi->recv(&fx[i],1,i%numSlaves+1);
				}
			}else{
				for(int i=0;i<xs.size();i++){
					if((i%numSlaves+1)==mpi->getID()){
						mpi->recv(&xs[i][0],xs[i].size(),0);
					}
				}
				for(int i=0;i<xs.size();i++){
					if((i%numSlaves+1)==mpi->getID()){
						fx[i]=f->evaluate(&xs[i][0]);
					}
				}
				for(int i=0;i<xs.size();i++){
					if((i%numSlaves+1)==mpi->getID()){
						mpi->send(&fx[i],1,0);
					}
				}
			}
		}
		void calulateBestI(){
			for(int i=0;i<PopSize;i++){
				if(fx[i]<fx[bestI]){
					bestI=i;
				}
			}
		}
		virtual double getMin(Function *f,int MaxFEs,vector<double>&out_x,double &out_fx){
			//if(mpi->isMaster()){
				this->f=f;
				numDim=f->getNumDim();
				//population initializing....
				//allocate space.
				tx.resize(numDim);
				x.resize(PopSize);
				fx.resize(PopSize);
				for(int i=0;i<PopSize;i++){
					x[i].resize(numDim);
				}
				bestI=0;
				for(int i=0;i<PopSize;i++){
					for(int d=0;d<numDim;d++){
						x[i][d]=drand(range[d][0],range[d][1]);
					} 
				}
				evaluatePopulation(x,fx);
				calulateBestI();
				//update, main process.
				update(MaxFEs/PopSize-1);
				out_x=x[bestI];
				out_fx=fx[bestI];
				return out_fx;
				/*
			}else{
				//only evaluate the f(x).
				return 0;
			}
			*/
		}
};
int main(int argc,char *argv[]){
	//srand(time(NULL));
	MPIHelper mpi(argc,argv);
//	Save s("title","x","y");
//	s.add(mpi.getID());
//	s.add(100);
//	ofstream out;
//	out.open("abc.txt");
//	out<<"hello from id:"<<mpi.getID()<<endl;
//	out.close();
	Trace(mpi.getName());
	Trace(mpi.getID());
}
int odl_main(int argc,char *argv[]){
	//srand(time(NULL));
	MPIHelper mpi(argc,argv);
//	SearchParam param("DE.json");
	SearchParam param("MyF1.json");
	ParallelDE *de=new ParallelDE(&mpi);
	de->initParam(&param);
	if(mpi.isMaster()){
	cout<<"Runing "<<de->getName()<<" "<<endl;
	}
	Tic::tic("begin");
	vector<double>x;
	double fx;
	cout<<"FunName(MyBestF,Optima)"<<endl;
	//Function*f=new PECFunction(param.getInt("NumDim"));
	Function*f=new MyF1(param.getInt("NumDim"));
	//Trace(de->getMin(f,param.getInt("MaxFEs"),x,fx));
	Trace(de->getMax(f,param.getInt("MaxFEs"),x,fx));
	if(mpi.isMaster()){
		printf("%s(%g,%g)",f->getName(),fx,f->getFBest());
	}
	Tic::tic("end");
	delete de;
	return 0;
}
