// ------------------------------------------------------------
//  
//    QGLikelihoodCalculator - Class
//    for the computation of the QG likelihood.
//    Needs files provided by having run the
//    Ntp1Finalizer_QG on QCD samples.
//
// ------------------------------------------------------------

#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"

#ifndef QGL
#define QGL

class QGLikelihoodCalculator {

 public:

  QGLikelihoodCalculator( const char * nchargedFileName="ncharged_fitresults.root",const char * nneutralFileName="nneutral_fitresults.root",const char * PtDFileName="PtD_fitresults.root");
  virtual ~QGLikelihoodCalculator() {};

  float computeQGLikelihood( float pt, int nCharged, int nNeutral, float ptD );

 private:

  TFile *nchargedFile, *nneutralFile,*PtDFile;
  double x;
  double *nChargedPar,*nNeutralPar,*PtDPar;

  float Interpolate(float pt, TGraph *g);
  //Fit functinos 
  
   inline double gammadistr_(double* x, double* par);
   inline double functionPtD_(double * x ,double*par);

};

#endif
