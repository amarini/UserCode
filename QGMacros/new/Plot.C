#include "Plot.h"
#include "TLegend.h"
#include <stdio.h>
using namespace std;
Plot::Plot()
{
sprintf(DrawOption,"HIST SAME");
}
Plot::~Plot()
{}
//=========================================
int Plot::Draw()
{
if(histos.size()==0)return 1;
c=new TCanvas("c1","c1",800,800);
histos.at(0)->SetMaximum(Max*1.1);
histos.at(0)->SetMinimum(Min*0.9);
histos.at(0)->Draw("AXIS");
histos.at(0)->Draw("AXIS X+ Y+ SAME");
for(int i=0;i<histos.size();i++)	
	histos.at(i)->Draw(DrawOption);

TLegend *L=new TLegend(0.3,0.4,0.5,0.6,"");
for(int i=0;i<histos.size();i++)	
	L->AddEntry(histos.at(i)->GetName(),(i==0)?"Quark":"Gluon","F");

return 0;
}
//=========================================
int Plot::AddHisto(TH1F*a)
{
TH1F*b=(TH1F*)a->Clone();
b->Scale(1./b->Integral());
if(histos.size()==0)
	{
	Max=b->GetMaximum();
	Min=b->GetMinimum();
	}
else
	{
	Max=(Max>b->GetMaximum())?Max:b->GetMaximum();
	Min=(Min<b->GetMinimum())?Min:b->GetMinimum();
	}
histos.push_back(b);
return 0;
}
//=========================================
int Plot::SetDefaultColors()
{
	switch(histos.size())
	{
	case 2:
		histos.at(1)->SetLineColor(kRed+2);
		histos.at(1)->SetLineWidth(2);
		histos.at(1)->SetFillColor(kRed+2);
		histos.at(1)->SetFillStyle(3005);
	case 1:
		histos.at(0)->SetLineColor(kBlue+2);
		histos.at(0)->SetLineWidth(2);
		histos.at(0)->SetFillColor(kBlue+2);
		histos.at(0)->SetFillStyle(3004);
	case 0:
		;
	}
	return 0;
}
int Plot::SetDrawOptions(const char *a)
{
sprintf(DrawOption,a);
}
