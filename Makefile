

ROOTCFLAGS = `root-config --cflags`
ROOTLIBS   = `root-config --libs` -lTreePlayer -lMinuit

CXXFLAGS += -I. -w

CXX = g++

INCS = calcResErr.h atmFitPars.h covBANFF.h covBase.h covXsec.h fqEvent.h fqProcessedEvent.h histoCompare.h histoFactory.h histoManager.h mcmcReader.h preProcess.h shared.h splineFactory.h splineParReader.h t2kHistoFactory.h t2kPreProcess.h t2kSplineFactory.h t2kfqEvent.h t2kfqReader.h visRing.h movieTools.h

OBS = movieTools.o calcEnu.o markovDiff.o mcmcReader.o globalRandom.o fqProcessedEvent.o FVCalculators.o TH2FV.o ThrowParms.o atmFitPars.o contInt.o covBANFF.o covBase.o covXsec.o getSystWeight.o hSplines.o histoCompare.o histoFactory.o histoManager.o histoTransforms.o keyread.o likelihood.o makeCov.o  markovTools.o  mcmcReader.o  preProcess.o  sharedPars.o splineFactory.o splineParReader.o visRing.o fqEvent.o masktools.o calcResErr.o

XLOBS = markovParManager.o t2kHistoFactory.o t2kPreProcess.o t2kSplineFactory.o t2kfqEvent.o t2kfqReader.o

SRCS = $(wildcard *.cxx)
#OBS = $(SRCS:.cxx=.o)
MOREOBS = fqEvent.o

printobs: 
	echo $(OBS)

#fqProcessedEvent.o : fqProcessedEvent.h
#	$(CXX) -c $(CXXFLAGS) $(ROOTCFLAGS) -o $@ $< 

#fqEvent.o : fqEvent.h
#	$(CXX) -c $(CXXFLAGS) $(ROOTCFLAGS) -o $@ $< 

%.o: %.cxx $(INCS)
	$(CXX) -c $(CXXFLAGS) $(ROOTCFLAGS) -o $@ $<

runMCMCFit: runMCMCFit.o $(OBS) $(MOREOBS) 
	$(CXX) -o ./bin/$@ $(CXXFLAGS) $(ROOTCFLAGS) -L. $^ $(ROOTLIBS)

runCalcDiff: runCalcDiff.o $(OBS) $(MOREOBS) 
	$(CXX) -o ./bin/$@ $(CXXFLAGS) $(ROOTCFLAGS) -L. $^ $(ROOTLIBS)

runDEMCMCFit: runDEMCMCFit.o  $(OBS) $(MOREOBS) 
	$(CXX) -o ./bin/$@ $(CXXFLAGS) $(ROOTCFLAGS) -L. $^ $(ROOTLIBS)

all: runMCMCFit runCalcDiff runDEMCMCFit

clean:
	rm -rf ./*.o
	rm -rf ./*.so
	rm -rf ./*.d
	rm -rf ./bin/runMCMCFit
	rm -rf ./bin/runDEMCMCFit
	rm -rf ./bin/runCalcDiff
