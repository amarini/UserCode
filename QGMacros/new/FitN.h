#include "TH1F.h"
#include "TF1.h"
#include "TGraphErrors.h"

/*********************************************** 
 *Author: Andrea Carlo Marini                  *
 *Email: andrea.carlo.marini@cern.ch           *
 *Description: Fit excluding points            *
 *                                             *
 *                                             *
 ***********************************************/

class FitN
{
public:
	FitN();
	~FitN();
	//Prepare the TGraphError from a TH1F
        //if Error > E and x in range
	int SetTH1F(TH1F*a);
	int SetTF1(TF1*b);
	//Fit the TGraph with the TF1
	double Fit();
	double Fit(TGraphErrors*a);
	//Minimum Probability
	int SetP(double x);
	//Minimum Number of Points
	int SetN(int a);
	//Minimum Error of each points
	int SetE(double x);
	int SetRange(double x,double y);
private:
	TH1F* h;
	TGraphErrors *G;
	TF1 *f;
	double p;//=0.05;
	double E;
	int N;
	double xmin,xmax;
};
