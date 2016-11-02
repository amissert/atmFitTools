#ifndef SELECTOR_H
#define SELECTOR_H


#include <iostream>

using namespace std;

class eventSelector{

 public:

 eventSelector();

 int selectNuMu(int nhitac,double momentum, double pidlike, double nring, double wall=0, double towall=0);

};







#ifdef CINTMODE
#include "eventSelector.cxx"
#endif




#endif
