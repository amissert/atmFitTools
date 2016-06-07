#ifndef SPLINEFACTORY_C
#define SPLINEFACTORY_C

#include "splineFactory.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//run spline factory using parameters from file
void splineFactory::runSplineFactory(){
  

  
  //make them splines
  buildTheSplines();

  ////////////////////
  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void splineFactory::resetModHistos(){
  for (int ibin=0;ibin<nBin;ibin++){
    for (int isamp=0;isamp<nSamp;isamp++){
      for (int icomp=0;icomp<nComp;icomp++){
        for (int iatt=0;iatt<nAtt;iatt++){
          for (int ipt=0;ipt<NPTSMAX;ipt++){
            hMC[isamp][ibin][icomp][iatt][ipt]->Reset();
          }
        }
      }
    }
  }
  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void splineFactory::setMCTree(TTree* tr){
  mcTree = tr;
  mcEvt = new fqProcessedEvent(mcTree);
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void splineFactory::debugtest(){
  //create new factory
  splineFactory* sfactory = new splineFactory(3,3,7,1,1,"debugtest");
  //set histogram manager and create template histograms
  sfactory->makeManagerFromFile("factoryOut_factorytest.root");
  

  return;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fills branches of output tree with information required to build a spline 
void splineFactory::fillBranches(int isamp,int ibin,int icomp,int iatt,int isyst){

   // these numbers identify which histogram bin the spline is for
   nsample = isamp;
   nbin = ibin;
   ncomponent=icomp;
   nattribute=iatt;
   nsystpar = isyst;

   // loop over the evaluation points for this bin to create array of bin weights
   nhistobins = hMC[isamp][ibin][icomp][iatt][0]->GetNbinsX();
   npoints = NPTSMAX;
   for (int ipt=0;ipt<NPTSMAX;ipt++){
     for (int jhistobin=0;jhistobin<=nhistobins;jhistobin++){
       if (hManager->getHistogram(isamp,ibin,icomp,iatt)->GetBinContent(jhistobin)==0){
         binWeight[ipt][jhistobin] = 1.;
       }
       else{
         binWeight[ipt][jhistobin] =
           ((double)hMC[isamp][ibin][icomp][iatt][ipt]->GetBinContent(jhistobin))/
           (double)hManager->getHistogram(isamp,ibin,icomp,iatt)->GetBinContent(jhistobin);
       }
     }
   }

   // fill array of all systematic parameter points used to fill histograms
   double sigvals[5] = {-5.,-2.,0.,2.,5.};
   for (int jpt=0;jpt<NPTSMAX;jpt++){
     incrementSystPars(isyst,sigvals[jpt]);
     systParValues[jpt]=fitPars->sysPar[isyst];
   }

   // a spline can now be created to interpolate binWeight[ibin][isyst] as a function of
   // systParValues[isyst]

   return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Build splines using event-by-event reweighting for various error 
void splineFactory::buildTheSplines(){

  // file setup
  TString fname = nameTag.Data();
  fname.Append("_splineOut.root");
  fout = new TFile(foutName.Data(),"recreate");

  //setup the output tree
  splineTree = new TTree("splinePars","spinePars");
  splineTree->Branch("nbin",&nbin,"nbin/I");
  splineTree->Branch("nhistobins",&nhistobins,"nhistobins/I");
  splineTree->Branch("ncomponent",&ncomponent,"ncomponent/I");
  splineTree->Branch("nattribute",&nattribute,"nattribute/I");
  splineTree->Branch("nsample",&nsample,"nsample/I");
  splineTree->Branch("nsystpar",&nsystpar,"nsystpar/I");
  splineTree->Branch("nsyspartot",&nSyst,"nsyspartot/I");
  splineTree->Branch("npoints",&npoints,"npoints/I");
  splineTree->Branch("systParValues",systParValues,Form("systParValues[%d]/D",NPTSMAX));
  splineTree->Branch("binWeight",binWeight,Form("binWeight[%d][%d]/D",NPTSMAX,NHBINSMAX));
  

  //setup systematic deviations (in sigma)
  double sigvals[5] = {-4.,-2.,0.,2.,4.};

  cout<<"creating spines"<<endl; 

  //loop over each systematic parameter specified in the parameter file
  for (int isyst=0;isyst<fitPars->nSysPars;isyst++){
    resetModHistos(); //< sets all bin contents to zero
    cout<<"splineFactory: Filling histograms while modifying systematic parameter: "<<isyst<<endl;
    //loop over the MC events
    for (int iev=0;iev<mcTree->GetEntries();iev++){
      mcTree->GetEvent(iev);  //< read event
      //fillAttributes();  //fill all attributes of this event
      //loop over the values of this parameter
      for (int ipt=0;ipt<NPTSMAX;ipt++){
        incrementSystPars(isyst,sigvals[ipt]); //< set parameter value
      //  getEvtWeight(isyst); //< get new new weight for this event 
        getEvtWeight(mcEvt,isyst,fitPars->sysPar[isyst]);
        //loop over fiTQun attributes and fill histograms
        for (int iatt=0;iatt<nAtt;iatt++){
          hMC[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]
            ->Fill(mcEvt->attribute[iatt],eventWeight);
        } //end attribute loop
      }  //end point of evaluation loop
    }//end event loop
    //loop over all histograms to create spline information
    for (int kbin=0;kbin<nBin;kbin++){
      for (int kcomp=0;kcomp<nComp;kcomp++){
        for (int ksamp=0;ksamp<nSamp;ksamp++){
          for (int katt=0;katt<nAtt;katt++){
            //get spline information for this bin 
            fillBranches(ksamp,kbin,kcomp,katt,isyst);
            splineTree->Fill();
          }
        }
      }
    }
  } 
 
  //write output tree with spline parameters
  splineTree->Write();

  //close output file
  fout->Close();

  ///////
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void splineFactory::incrementSystPars(int isyspar, double nsig){

  //reset parameters to defaults
  fitPars->resetDefaults();
  
  //adjust systematic parameter in question
  double adjustment = nsig*fitPars->sysParUnc[isyspar];
  double value = fitPars->sysPar[isyspar];
  fitPars->setSysParameter(isyspar,(value+adjustment));

  
  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Initialize the histogram arrays needed to builde the splines
void splineFactory::setupHistos(){
  cout<<"called setupHistos()"<<endl;
  TH1D* htemplate;
  TString hname;
  for (int ipt=0;ipt<NPTSMAX;ipt++){
    for (int isamp=0;isamp<nSamp;isamp++){
      for (int ibin=0;ibin<nBin;ibin++){
        for (int icomp=0;icomp<nComp;icomp++){
          for (int iatt=0;iatt<nAtt;iatt++){
            htemplate = hManager->getHistogram(isamp,ibin,icomp,iatt);
            hname = htemplate->GetName();
            hname.Append(Form("_%d%d%d%d%d",ipt,isamp,ibin,icomp,iatt));
            cout<<"setting histogram template: "<<hname.Data()<<endl;
            hMC[isamp][ibin][icomp][iatt][ipt]=(TH1D*)htemplate->Clone(hname.Data());
            hMC[isamp][ibin][icomp][iatt][ipt]->Reset();
          }
        }
      }
    }
  }
  return;
}
 

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void splineFactory::makeManagerFromFile(const char* fname){
  hManager = new histoManager(fname,nSamp,nBin,nComp,nAtt);
  //setupHistos();
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Fill array of fiTQun outputs to be fit
void splineFactory::fillAttributes(){


  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void splineFactory::fillHistograms(int ipt, int isyst){
  //fills histograms after applying systematic error parameter
  
  //reset bin contents
  for (int jbin=0;jbin<nBin;jbin++){
    for (int jsamp=0;jsamp<nSamp;jsamp++){
      for (int jatt=0;jatt<nAtt;jatt++){
        for (int jcomp=0;jcomp<nComp;jcomp++){
          hMC[jsamp][jbin][jcomp][jatt][ipt]->Reset();
        }
      }
    }
  }

  //fill new bin contents
  for (int iev=0;iev<mcTree->GetEntries();iev++){
     mcTree->GetEvent(iev);
     fillAttributes();
    // getEvtWeight(isyst);
     getEvtWeight(mcEvt,isyst,fitPars->sysPar[isyst]);
     for (int jatt=0;jatt<nAtt;jatt++){
      hMC[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][jatt][ipt]->Fill(attribute[jatt],eventWeight);
    }
  }

  return;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Re-weight the event given the value "value" of additional systematic parameter "ipar"
//Character keys determine the type of systematic error parameters to use
double splineFactory::getEvtWeight(fqProcessedEvent* mcevent,int ipar,double value){;;;
  
  //first get the original weight of this event:
  double ww = mcevent->evtweight;

  //get the additional weight from systematic parameters
  double systweight = getSystWeight(sysParType.Data(),mcevent,ipar,value);

  //apply additional weight
  ww *= systweight;

  //set internal variable
  eventWeight = ww; 

  //////////
  return ww;


}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Construct from parameter file
splineFactory::splineFactory(const char*  parfile){

  //set file name and read in pars
  parFileName = parfile;
  cout<<"splineFactory: Parameter File: "<<parfile<<endl;
  runpars = new sharedPars(parFileName.Data());
  runpars->readParsFromFile();

  //setup systematic parameters object
  fitPars = new atmFitPars(parfile);

  //fix internal variables
  nSamp=runpars->nSamples;
  nBin=runpars->nFVBins;
  nComp=runpars->nComponents;
  nAtt=runpars->nAttributes;
  nSyst=runpars->nSysPars;
  nameTag=runpars->globalRootName.Data();
  foutName = runpars->splineFactoryOutput.Data(); //< name of output file with spline parameters
  sysParType = runpars->sysParType; 

  //create histogram manager from file with prefilled histograms
  makeManagerFromFile(runpars->hFactoryOutput.Data());

  //initialize histogram arrays using the hManager as a template
  setupHistos();

  
  //get pointer to MC tree
  TChain *chmc = new TChain("h1");
  chmc->Add(runpars->hFactoryMCFiles.Data());
  setMCTree((TTree*)chmc);

  return;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Original constructor
splineFactory::splineFactory(int isamp, int ibin, int icomp, int iatt, int isyst, const char* name){
  nameTag = name;
  nSamp = isamp;
  nBin  = ibin;
  nComp = icomp;
  nAtt  = iatt;
  nSyst = isyst;
  foutName = "splineFactoryOutput.root";
}


#endif
