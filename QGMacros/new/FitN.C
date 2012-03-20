#include "FitN.h"
#include "TMath.h"
#include <stdio.h>
#define MAXPAR 10

FitN::FitN()
{
p=0.05;
N=4;
E=0.001;
xmin=0;xmax=1000;
G=NULL;
h=NULL;
f=NULL;
}

FitN::~FitN()
{}

int FitN::SetTH1F(TH1F*a)
	{
	if(a==NULL)return -1;
	h=a;
	G=new TGraphErrors();
	int k=0;
	for(int i=1;i<=h->GetNbinsX();i++)
		if((h->GetBinError(i)>E)&&(h->GetBinCenter(i)>xmin)&&(h->GetBinCenter(i)<xmax)){G->SetPoint(k,h->GetBinCenter(i),h->GetBinContent(i));G->SetPointError(k,h->GetBinWidth(i)/2.0,h->GetBinError(i));k++;}
	return k;
	}
int FitN::SetTF1(TF1*b)
	{
	if(b==NULL)return 1;
	f=b;return 0;
	}
double FitN::Fit()
	{
	return Fit(G);
	}
double FitN::Fit(TGraphErrors*a)
	{
//	if(f==NULL)return -1;
	double par[MAXPAR];
	f->SetRange(xmin,xmax);	
	G->Fit(f->GetName(),"N Q R");
	G->Fit(f->GetName(),"N Q M R");
	for(int i=0;i<f->GetNpar();++i)par[i]=f->GetParameter(i);
	double R=a->Chisquare(f);
	double P=TMath::Gamma( (a->GetN()-f->GetNpar())/2.0, R/2.0);
	if(TMath::ChisquareQuantile(p,a->GetN()-f->GetNpar()) < R )
	if(G->GetN()>N)//we do not want low n of bins
		{
		for(int i=0;i<a->GetN();i++)
			{
			TGraphErrors *tmp=(TGraphErrors*)a->Clone("tmp");
			tmp->RemovePoint(i);
			double r;
			tmp->Fit(f->GetName(),"N Q R");
			tmp->Fit(f->GetName(),"N Q M R");
			r=TMath::Gamma( (tmp->GetN()-f->GetNpar())/2.0, (tmp->Chisquare(f))/2.0  );
			if(r<P){
				for(int j=0;j<f->GetNpar();++j)par[j]=f->GetParameter(j);
				P=r;
				Fit(tmp);
				}
			else {
				for(int j=0;j<f->GetNpar();++j)f->SetParameter(j,par[j]);
				}
			delete tmp;
			}
		}
	return P;
	}
int FitN::SetP(double x){p=x;return 0;}
int FitN::SetN(int a){N=a;return 0;}
int FitN::SetE(double x){E=x;return 0;}
int FitN::SetRange(double x,double y){ xmin=x; xmax=y;return 0;}
