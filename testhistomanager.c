{
gROOT->ProcessLine(".L histoManager.C++");
gROOT->ProcessLine(".L hSplines.C+");
//gROOT->ProcessLine(".x ~/style.c");
gROOT->ProcessLine(".L atmFitPars.C+");

int nbin=3;
int ncomp=7;
int nsamp=3;
int natt=1;

//atmFitPars* fitpars = new atmFitpars(nbin,ncomp,nsamp,natt,0);
//TChain chdat("h1");
//TChain chmc("h1");
//chdat.Add("/nfs/data40/t2k/amissert/skdata/cosmicsMultiGeV/14a/data/patmue.run073010.004*");
//chmc.Add("/nfs/data40/t2k/amissert/skdata/cosmicsMultiGeV/14a/mc/stopmu_mc_1390_for_atmnu_test_0649*001*fQ.root");
TFile fdata("nominal2.root");
TFile fmc("nominal1.root");
atmFitPars* fitpars = new atmFitPars(nbin,ncomp,nsamp,natt,1);
 
TTree* trdata = (TTree*)fdata.Get("h1");
TTree* trmc   = (TTree*)fmc.Get("h1");
//histoManager* hm = new histoManager(3,3,7,"test2"); 
histoManager* hm = new histoManager("factoryOut_factorytest.root",3,3,7,1); 
hm->readSplinesFromFile("splineOut_debug.root");
hm->setFitPars(fitpars);
//hSplines* hs = hm->getSplines(0,0,0,0);

//hm->addAttribute(1);
//hm->addAttribute(2);
//hm->setDataTree(trdata);
//hm->setMCTree(trmc);
//hm->init();
//hm->readFromFile("hManager_atmos_emuratio");
//hm->fillHistos();
//hm->saveToFile();
}
