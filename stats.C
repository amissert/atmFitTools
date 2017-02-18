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
float arrayvar(float arr[], int n){

  float mean = arraymean(arr,n);
  float var = 0.;
  float N = 0.;

  for (int i=0; i<n; i++){
    var += ((arr[i]-mean)*(arr[i]-mean));
    N++;
  }

  return var/N;
}


///////////////////////////////////////////////
float arraycov(float arrx[], float arry[], int n){

  // intermediate vars 
  float meanx = arraymean(arrx,n);
  float meany = arraymean(arry,n);
  float N = (float)n;
  float cov=0.;

  if ((N-1) <= 0.) return 0.;
  for (int i=0; i<n; i++){
    cov += (arrx[i]-meanx)*(arry[i]-meany)/(N-1.); 
  }

  return cov;
}


///////////////////////////////////////////////
float arraycor(float arrx[], float arry[], int n){

  // intermediate vars 
  float meanx = arraymean(arrx,n);
  float meany = arraymean(arry,n);
  float varx = TMath::Sqrt(arrayvar(arrx,n));
  float vary = TMath::Sqrt(arrayvar(arry,n));
  
  // checks
  if (varx*vary == 0) return 0.;
  float N = (float)n;
  if ((N-1) <= 0.) return 0.;

  // calc
  float cov=0.;
  for (int i=0; i<n; i++){
    cov += (arrx[i]-meanx)*(arry[i]-meany)/(N-1.); 
  }

  return cov/(varx*vary);
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
