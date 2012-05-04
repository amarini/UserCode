#include "Fit2.h"
#include <stdio.h>
#include "TAxis.h"
#include "TCanvas.h"
#ifndef NULL
	#define NULL 0
#endif

#define MINERROR 0.001
#define M 100

//constructor
Fit2::Fit2(){
//	f=new TF1("pol1","[0]+[1]*x",0,100);
//	f=NULL;
}
//destructor
Fit2::~Fit2()
{
//	if(f!=NULL)delete f;
}

//Set The Histograms
void Fit2::SetTGRAPH(TGraphErrors*h,int n)
{
	switch (n)
	{
	case 0: a=h;break;
	case 1: b=h;break;
	default: break;
	}
return ;
}

double myfunc(double*x,double*par)
	{ return (x[0]<M*.9)?(par[0]+x[0]*par[1]):(par[2]+par[1]*(x[0]-M)); }	

std::vector<float> *Fit2::Fit()
{
//M=a->GetXaxis()->GetXmax()+20.;
char str[1023];
//sprintf(str,"(x<%f)?([0]+x*[1]):([2] + (x-%f)*[1])",M,M);
//f=new TF1("2pol1",str,3,0,1000);
TGraphErrors *R=new TGraphErrors();int count=0;
f=new TF1("2pol1",myfunc,0,1000,3);

for(int i=0;i<a->GetN();i++)
	{
	double x,y,ex,ey;
	a->GetPoint(i,x,y);
	if(a->GetErrorY(i)/y>MINERROR){
		R->SetPoint(count,x,y);	
		R->SetPointError(count,a->GetErrorX(i),a->GetErrorY(i));
		count++;
		}
	}

for(int i=0;i<b->GetN();i++)
	{
	double x,y,ex,ey;
	b->GetPoint(i,x,y);
	if(b->GetErrorY(i)/y>MINERROR){
		R->SetPoint(count,x+M,y);	
		R->SetPointError(count,b->GetErrorX(i),b->GetErrorY(i));
		count++;
		}
	}
fprintf(stderr,"Total Numebr of points=%d\n",count);
R->Fit("2pol1","W");
std::vector<float> *v=new std::vector<float>;
v->push_back(f->GetParameter(0));
v->push_back(f->GetParameter(1));
v->push_back(f->GetParameter(2));
//DEBUG
//TCanvas *C=new TCanvas();
//R->Draw("A P");
//f->DrawCopy("SAME");
//END DEBUG
return v;
}
//Create a TGraphErrors which olds the points
void Fit2::SetTH1F(TH1F*h,int n)
{
TGraphErrors*G=new TGraphErrors();
for(int i=0; i<=h->GetNbinsX();i++)
	{
	G->SetPoint(i,h->GetBinCenter(i),h->GetBinContent(i) );
	G->SetPointError(i,h->GetBinWidth(i)/2.0,h->GetBinError(i)   );
	}
return SetTGRAPH(G,n);
}

