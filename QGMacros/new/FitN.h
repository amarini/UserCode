#include "TH1F.h"
#include "TF1.h"
#include "TGraphErrors.h"

class FitN
{
public:
	FitN();
	~FitN();
	int SetTH1F(TH1F*a);
	int SetTF1(TF1*b);
	double Fit();
	double Fit(TGraphErrors*a);
	int SetP(double x);
	int SetN(int a);
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
