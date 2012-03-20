#include "TH1F.h"
#include "TCanvas.h"
#include <vector>
#ifndef PLOT_H
#define PLOT_H

class Plot
{
public:
	Plot();	
	~Plot();
	int Draw();
	int AddHisto(TH1F*);
	int SetDefaultColors();
	int SetDrawOptions(const char *);
private:
 	TCanvas *c;
	std::vector<TH1F*> histos;
	double Max;
	double Min;
	char DrawOption[255];
};
#endif
