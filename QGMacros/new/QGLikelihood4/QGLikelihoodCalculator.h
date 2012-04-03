// ------------------------------------------------------------
//  
//    QGLikelihoodCalculator - Class
//    for the computation of the QG likelihood.
//    Needs files provided by having run the
//    double fit procedure. output.txt?
//
// ------------------------------------------------------------

#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph2D.h"
#include <map>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>
using namespace std;


class QGLikelihoodCalculator {

 public:

  QGLikelihoodCalculator( const std::string& fileNameNC , const std::string &fileNameNN,const std::string& fileNamePTD);
   ~QGLikelihoodCalculator();

  float computeQGLikelihoodPU( float pt, float rhoPF, int nCharged, int nNeutral, float ptD, float rmsCand=-1. );


 private:

    FILE *fr_;
    int ReadBinTxt(const char*fileName,int parameter,TH2F *quark,TH2F *gluon);
    int ReadParTxt( const char *fileName,std::map< pair<int,int>, double*> *parq,std::map< pair<int,int>, double*> *parg ,int nPar=2 );
    char GetChar(FILE *fr);
    TH2F* ptD0_q,*ptD0_g;
    TH2F* ptD1_q,*ptD1_g;
    TH2F* ptD2_q,*ptD2_g;
    std::map< std::pair<int,int>,double* > *parqNC,*parqNN,*pargNC,*pargNN,*parqPTD,*pargPTD;

};

