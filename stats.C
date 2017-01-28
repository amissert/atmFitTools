#ifndef STATS_C
#define STATS_C

#include "TMath.h"
#include <iostream>


//////////////////////////////////////////////
float arraymean(float arr[], int n){
  
  float mean = 0.;
  float N = 0.;

  for (int i=0; i<n; i++){
    mean += arr[i];
    N++;
  }

  return mean/N;
}



///////////////////////////////////////////////
float arrayvar(float arr[], int n, float mean){

  float var = 0.;
  float N = 0.;

  for (int i=0; i<n; i++){
    var += ((arr[i]-mean)*(arr[i]-mean));
    N++;
  }

  return var/N;
}



//////////////////////////////////////////////
double arraymeanD(double arr[], int n){
  
  double mean = 0.;
  double N = 0.;

  for (int i=0; i<n; i++){
    mean += arr[i];
    N++;
  }
//  cout<<"N: "<<N;
  return mean/N;
}



///////////////////////////////////////////////
double arrayvarD(double arr[], int n){

  double var = 0.;
  double N = 0.;
  double mean = arraymeanD(arr,n);

  for (int i=0; i<n; i++){
    var += ((arr[i]-mean)*(arr[i]-mean));
    N++;
  }

  return var/N;
}



#endif
