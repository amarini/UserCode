#include "QGLikelihoodCalculator.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include "TMath.h"
#include "TF1.h"
#include <map>
#include "ReadTxt.C"
#include <stdio.h>
using namespace std;

//#define DEBUG


// constructor:

QGLikelihoodCalculator::QGLikelihoodCalculator( const std::string& fileNameNC , const std::string &fileNameNN,const std::string& fileNamePTD) {
	parqNC=new map< pair<int,int>, double* >;
	pargNC=new map< pair<int,int>, double* >;
	parqNN=new map< pair<int,int>, double* >;
	pargNN=new map< pair<int,int>, double* >;
	parqPTD=new map< pair<int,int>, double* >;
	pargPTD=new map< pair<int,int>, double* >;
	
	ReadParTxt(fileNameNC.c_str(),parqNC,pargNC);
	ReadParTxt(fileNameNN.c_str(),parqNN,pargNN);
	ReadParTxt(fileNamePTD.c_str(),parqPTD,pargPTD,3);
	

}

// ADD map destructor
QGLikelihoodCalculator::~QGLikelihoodCalculator()
{
}

inline double gammadistr_(double* x, double* par)
{
        return TMath::Exp( - x[0] *par[0]/par[1] ) * TMath::Power(x[0],par[0]-1) * TMath::Power(par[1]/par[0],-par[0])/TMath::Gamma(par[0]) ;
}

//half gamma+ offset
inline double functionPtD_(double * x ,double*par)
{
        if((x[0]-par[0])<0)return 0;
        return TMath::Exp( - (x[0]-par[0]) *par[1]/par[2] ) * TMath::Power((x[0]-par[0]),par[1]-1) * TMath::Power(par[2]/par[1],-par[1])/TMath::Gamma(par[1]) ;
}


float QGLikelihoodCalculator::computeQGLikelihoodPU( float pt, float rhoPF, int nCharged, int nNeutral, float ptD, float rmsCand ) {
//plots are already normalized to unity.
double Q=1;
double G=1;


TF1 *pol3=new TF1("pol3","[0]+[1]*TMath::Log(x)+[2]*TMath::Log(x)*TMath::Log(x)+[3]*TMath::Log(x)*TMath::Log(x)*TMath::Log(x)",20,3500);//LOG! ->PT
TF1 *pol1=new TF1("pol1","[1]+[0]*x",0,20); //NOT LOG -> RHO

double *par=new double[5];
double *x=new double[5];
double a,b;
x[0]=nCharged;
//NC Q
pol3->SetParameters( (*parqNC)[pair<int,int>(0,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*parqNC)[pair<int,int>(0,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[0]=pol1->Eval(rhoPF);
pol3->SetParameters( (*parqNC)[pair<int,int>(1,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*parqNC)[pair<int,int>(1,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[1]=pol1->Eval(rhoPF);
Q*=gammadistr_(x,par);
//NC G
pol3->SetParameters( (*pargNC)[pair<int,int>(0,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*pargNC)[pair<int,int>(0,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[0]=pol1->Eval(rhoPF);
pol3->SetParameters( (*pargNC)[pair<int,int>(1,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*pargNC)[pair<int,int>(1,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[1]=pol1->Eval(rhoPF);
G*=gammadistr_(x,par);

//fprintf(stderr,"NC %.2lf -> %.3lf %.3lf <-",Q/(Q+G),pol3->GetParameter(0),(((*pargNC)[pair<int,int>(1,1)])[0]) ); //DEBUG
x[0]=nNeutral;
//NN Q
pol3->SetParameters( (*parqNN)[pair<int,int>(0,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*parqNN)[pair<int,int>(0,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[0]=pol1->Eval(rhoPF);
pol3->SetParameters( (*parqNN)[pair<int,int>(1,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*parqNN)[pair<int,int>(1,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[1]=pol1->Eval(rhoPF);
Q*=gammadistr_(x,par);
//NN G
pol3->SetParameters( (*pargNN)[pair<int,int>(0,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*pargNN)[pair<int,int>(0,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[0]=pol1->Eval(rhoPF);
pol3->SetParameters( (*pargNN)[pair<int,int>(1,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*pargNN)[pair<int,int>(1,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[1]=pol1->Eval(rhoPF);
G*=gammadistr_(x,par);

//fprintf(stderr,"NCN %.2lf \n",Q/(Q+G));//DEBUG
x[0]=ptD;
//PtD Q
pol3->SetParameters( (*parqPTD)[pair<int,int>(0,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*parqPTD)[pair<int,int>(0,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[0]=pol1->Eval(rhoPF);
pol3->SetParameters( (*parqPTD)[pair<int,int>(1,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*parqPTD)[pair<int,int>(1,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[1]=pol1->Eval(rhoPF);
pol3->SetParameters( (*parqPTD)[pair<int,int>(2,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*parqPTD)[pair<int,int>(2,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[2]=pol1->Eval(rhoPF);
Q*=functionPtD_(x,par);
//PtD G
pol3->SetParameters( (*pargPTD)[pair<int,int>(0,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*pargPTD)[pair<int,int>(0,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[0]=pol1->Eval(rhoPF);
pol3->SetParameters( (*pargPTD)[pair<int,int>(1,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*pargPTD)[pair<int,int>(1,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[1]=pol1->Eval(rhoPF);
pol3->SetParameters( (*pargPTD)[pair<int,int>(2,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*pargPTD)[pair<int,int>(2,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[2]=pol1->Eval(rhoPF);
G*=functionPtD_(x,par);

delete[] par;
delete[] x;
delete pol3;
delete pol1;
if(Q==0)return 0;
return float(Q/(Q+G));
}

