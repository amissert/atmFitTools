#include "histoCompare.cxx"

using namespace std;


int main(int argc, char* argv[]){
 
 // get par file
 TString parfilename = argv[1];
 cout<<"Parameter file"<<parfilename.Data()<<endl;

 // make new histogram comparison object
 histoCompare* hc= new histoCompare(parfilename.Data());

 // run the mcmc
 hc->runMCMC(-1); 

 //
 return 0;

}


/*
int rootrun(const char* parfilename){

 // make new histogram comparison object
 histoCompare* hc= new histoCompare(parfilename);

 // run the mcmc
 hc->runMCMC(-1); 

 //
 return 0;
}

*/



