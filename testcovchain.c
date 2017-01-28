{
 gROOT->ProcessLine(".L makeCov.cxx+");

 TChain *ch = new TChain("MCMCpath");
 ch->Add("/nfs/data41/t2k/amissert/atmos/head/atmFitTools/run/results/demcmc*2.root");
 TTree* tr = (TTree*)ch;
 makeCov *maker = new makeCov("shimpars.dat");
 maker->setParTree(tr);
 maker->nburn = 40000;
 maker->buildMatrix();
}
