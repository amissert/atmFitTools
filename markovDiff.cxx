#ifndef MARKOVDIFF_CXX
#define MARKOVDIFF_CXX

#include "markovDiff.h"

using namespace std;



void markovDiff::setouttree(){


  // open output file
  outfile = new TFile(outFileName.Data(),"RECREATE");
  

  //use same file

  diffpars = new TTree("MCMCdiff","MCMCdiff");
  diffpars->Branch("npars",&npars,"npars/I");
  diffpars->Branch("par1",par1,"par1[500]/F");
  diffpars->Branch("par2",par1,"par2[500]/F");
  diffpars->Branch("parindex",parindex,"parindex[500]/I");
  diffpars->Branch("pardiff",pardiff,"pardiff[500]/F");

  return;
}

void markovDiff::setaddresses(){
 
  // only turn on the important mcmc branches
  mcmcpars->SetBranchStatus("*",0);
  mcmcpars->SetBranchStatus("npars",1);
  mcmcpars->SetBranchStatus("par",1);
  mcmcpars->SetBranchStatus("parindex",1);

  // set branch address
  mcmcpars->SetBranchAddress("par",par);
  mcmcpars->SetBranchAddress("parindex",parindex);

  // read in number of pars
  int nparstmp;
  mcmcpars->SetBranchAddress("npars",&nparstmp);
  mcmcpars->GetEntry(0);
  npars = nparstmp;

  cout<<"Number of mcmc pars: "<<npars<<endl;
  return; 
}

//////////////////////////////////////////////////////////////////
void markovDiff::fillDiffPars(int npairs){

  // set branch addresses of input tree
  setaddresses();

  // setup new output tree
  setouttree();

  int npts = mcmcpars->GetEntries() - nburn;

  // Loop over mcmc points and select random pairs of points after
  // some burn-in.  For these pairs, find the difference between the parameters and 
  // save this difference to the new difference tree
  mcmcpars->GetEntry(nburn); //< start somewhere
  for (int i=0; i<npairs; i++){
    if ((i%500)==0) cout<<i<<endl;
    // get a pair of randoms
  //  int ipoint = randy->Integer(npts); 
    int fpoint = randy->Integer(npts); 
   // mcmcpars->GetEntry(ipoint+nburn);
    for (int ipar=0; ipar<npars; ipar++){
      par1[ipar] = par[ipar];
    }
    mcmcpars->GetEntry(fpoint+nburn);
    for (int ipar=0; ipar<npars; ipar++){
      par2[ipar] = par[ipar];
      pardiff[ipar] = par2[ipar]-par1[ipar];
    }
    diffpars->Fill(); 
  }

  diffpars->Write();
  outfile->Close();

  return;
}


//////////////////////////////////////////////////////////////////
void markovDiff::fillDiffPars2(int npairs){

//  TRandom2* randy = new TRandom2(npairs);

  // set branch addresses of input tree
  setaddresses();

  // setup new output tree
  setouttree();

  int npts = mcmcpars->GetEntries() - nburn;

  // steps between difference evaluation
  int nbetween = 50;

  if (npars>(npts/nbetween)) npars = (npts/nbetween) - 1;

  // Loop over mcmc points and select random pairs of points after
  // some burn-in.  For these pairs, find the difference between the parameters and 
  // save this difference to the new difference tree
  int ichoose = nburn;
  for (int i=0; i<npairs; i++){

    if ((i%100)==0) cout<<i<<endl;
    // get a pair of randoms
    
    mcmcpars->GetEntry(ichoose);
    for (int ipar=0; ipar<npars; ipar++){
      par1[ipar] = par[ipar];
    }
    mcmcpars->GetEntry(ichoose+nbetween);
    for (int ipar=0; ipar<npars; ipar++){
      par2[ipar] = par[ipar];
      pardiff[ipar] = par2[ipar]-par1[ipar];
    }
    diffpars->Fill(); 

    ichoose+=nbetween;
  }

  diffpars->Write();
  outfile->Close();

  return;
}

//////////////////////////////////////////////////////////
// constructor
markovDiff::markovDiff(const char* fname, int burn, int chainflg){

  // read in single mcmc path
  if (chainflg==0){
    mcmcfile = new TFile(fname);
    mcmcpars = (TTree*) mcmcfile->Get("MCMCpath");
    nburn = burn; 
  }

  // read in many mcmc files
  else{
    TChain *mcmcChain = new TChain("MCMCpath");
    mcmcChain->Add(fname);
    mcmcpars = (TTree*)mcmcChain;
    nburn = burn;
  }

  outFileName = "mcmcdiff.root";
  return;
}



#endif


