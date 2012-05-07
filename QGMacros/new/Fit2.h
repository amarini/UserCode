/******************************************
 *Author: Andrea Carlo marini             *
 *email: andrea.carlo.marini@cern.ch      *
 *Description: class to fit two histos    *
 *     with a pol1 with the same angular  *
 *     coefficient                        *
 *                                        *
 ******************************************/
#include <vector>
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1F.h"
class Fit2
{
public:
	Fit2();
	~Fit2();
	std::vector<float> *Fit();
	void SetTGRAPH(TGraphErrors*h,int n=0);
	void SetTH1F(TH1F*h,int n=0);
private:
	TGraphErrors *a;
	TGraphErrors *b;
	TF1  *f;

	//static float M;// Maximum allowed value for an Histo
};
double myfunc(double*x,double*par);
