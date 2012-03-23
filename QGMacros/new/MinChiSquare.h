#include "TH1F.h"
#include "TGraph.h"

#ifndef MINCHISQUARE_H
#define MINCHISQUARE_H
class MinChiSquare
{
public:
	MinChiSquare();	
	~MinChiSquare();	
	static double GetMinChiSquare(TH1F*Q,TH1F*G,TH1F*T,double x=0.01,const char *options="CHI2/NDF WW",TGraph*graph=NULL);
	static double GetError(TGraph*graph);
private:
};

#endif
