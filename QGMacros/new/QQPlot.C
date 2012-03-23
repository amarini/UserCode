
#include "QQPlot.h"

QQPlot::QQPlot()
{
h1=NULL;
h2=NULL;
QQ=NULL;
}

QQPlot::~QQPlot()
{
if(h1)delete h1;
if(h2)delete h2;
if(QQ)delete QQ;
}

int QQPlot::AddHisto(TH1F*h,int n)
{
if(n<=0) 
	{
	if(h1==NULL) n=1;
	else if(h2==NULL) n=2;
	else return 1;
	}
if(n==1)h1=(TH1F*)h->Clone();
else if(n==2) h2=(TH1F*)h->Clone();
else return 1;
return 0; 
}

TGraph* QQPlot::GetQQPlot()
{
if(h1->GetNbinsX()!=h2->GetNbinsX())return NULL; // they must have the same number of bins
if(QQ) delete QQ;
QQ=new TGraph();
int count=0;
QQ->SetName("QQ");
h1->Scale(1./h1->Integral("width"));
h2->Scale(1./h2->Integral("width"));
for(int i=1;i<=h1->GetNbinsX();i++)
	{
	QQ->SetPoint(count,h1->Integral(1,i,"width"),h2->Integral(1,i,"width"));	
	count++;
	}
return QQ;
}
