#include "TROOT.h"

#ifndef BINS_H
#define BINS_H
class Bins
{
public:

      int getBins(double  *Bins,int nBins,double MinBin=15.0,double MaxBin=1000.,bool log=false);
      int getBin(int nBins,double  *Bins,double value,double*x0=0,double*x1=0);
      void getBins_int( int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog=true);
private:
};
#endif
