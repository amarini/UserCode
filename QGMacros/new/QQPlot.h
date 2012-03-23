#include "TH1F.h"
#include "TGraph.h"

class QQPlot
{
public:
	QQPlot();
	~QQPlot();
	int AddHisto(TH1F*h,int n=0);	
	TGraph* GetQQPlot();
private:
	TH1F *h1,*h2;
	TGraph* QQ;
};
