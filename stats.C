#ifndef STATS_C
#define STATS_C

#include "TMath.h"
#include "shared.h"

#include <iostream>

/*
//////////////////////////////////////////////
void fillnormal(float arr[], int n, float  mean, float sig){
   for (int i=0; i<n; i++){
      arr[i] = randy->Gaus(mean,sig);
   }
   return;
}


//////////////////////////////////////////////
void fillnormal(float arrx[], float arry[], int n
                ,float  mean, float sig){
   for (int i=0; i<n; i++){
      float val = randy->Gaus(mean,sig);
      arrx[i] = val;
      arry[i] = val;
   }
   return;
}
*/

const double ZEROISH = 1e-20;


//////////////////////////////////////////////
float arraymean(float arr[], int n){
  
  float mean = 0.;
  float N = (float)n;

  for (int i=0; i<n; i++){
    mean += arr[i];
  }

  return mean/N;
}



//////////////////////////////////////////////
float arraymin(float arr[], int n){
  
  float arrmin = 1e9;

  for (int i=0; i<n; i++){
    if (arr[i]<arrmin){
      arrmin = arr[i];
    }
  }

  return arrmin;
}


//////////////////////////////////////////////
double arrayminD(double arr[], int n){
  
  double arrmin = 1e9;

  for (int i=0; i<n; i++){
    if (arr[i]<arrmin){
      arrmin = arr[i];
    }
  }

  return arrmin;
}



//////////////////////////////////////////////
float arraymax(float arr[], int n){
  
  float arrmax = -1e9;

  for (int i=0; i<n; i++){
    if (arr[i]>arrmax){
      arrmax = arr[i];
    }
  }

  return arrmax;
}



//////////////////////////////////////////////
double arraymaxD(double arr[], int n){
  
  double arrmax = -1e9;

  for (int i=0; i<n; i++){
    if (arr[i]>arrmax){
      arrmax = arr[i];
    }
  }

  return arrmax;
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

  if (var<ZEROISH) return 0.;
  return var/(N-1.);
}


///////////////////////////////////////////////
float arraycov(float arrx[], float arry[], int n){

  // intermediate vars 
  float meanx = arraymean(arrx,n);
  float meany = arraymean(arry,n);
  float N = (float)n;
  float cov=0.;

  if ((N) <= 0.) return 0.;
  for (int i=0; i<n; i++){
    cov += (arrx[i]-meanx)*(arry[i]-meany)/(N); 
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
  if ((N) <= 0.) return 0.;

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
  if (var<ZEROISH) return 0.;
  return var/(N-1.);
}

///////////////////////////////////////////////
double arraycovD(double arrx[], double arry[], int n){

  // intermediate vars 
  double meanx = arraymeanD(arrx,n);
  double meany = arraymeanD(arry,n);
  double N = (double)n;
  double cov=0.;

  if ((N-1) <= 0.) return 0.;
  for (int i=0; i<n; i++){
    cov += (arrx[i]-meanx)*(arry[i]-meany)/(N); 
  }

  return cov;
}


///////////////////////////////////////////////
double arraycorD(double arrx[], double arry[], int n){

  // intermediate vars 
  double meanx = arraymeanD(arrx,n);
  double meany = arraymeanD(arry,n);
  double varx = TMath::Sqrt(arrayvarD(arrx,n));
  double vary = TMath::Sqrt(arrayvarD(arry,n));
  
  // checks
  if (varx*vary <= ZEROISH) return 0.;
  double N = (double)n;
  if ((N-1) <= 0.) return 0.;

  // calc
//  double cov=0.;
//  for (int i=0; i<n; i++){
//    cov += (arrx[i]-meanx)*(arry[i]-meany)/(N); 
//  }
  double cov = arraycovD(arrx,arry,n);
  return cov/(varx*vary);
}



#endif
