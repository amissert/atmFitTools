#ifndef SELECTOR_CXX
#define SELECTOR_CXX

#include "eventSelector.h"


eventSelector::eventSelector(){

}

int eventSelector::selectNuMu(double momentum, double pidlike, double nring, double wall, double towall){

  if (momentum<100.) return 0;
  if (pidlike>0.) return 0;
  if (nring>1.) return 0;


  return 1;
}


#endif
