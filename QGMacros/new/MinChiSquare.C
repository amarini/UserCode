#include <stdlib.h>
#include <stdio.h>
#include "MinChiSquare.h"
#include "math.h"

MinChiSquare::MinChiSquare(){}
MinChiSquare::~MinChiSquare(){}

double MinChiSquare::GetMinChiSquare(TH1F*Q,TH1F*G,TH1F*T,double x,const char *options,TGraph*graph)
{
Q->Scale(1./Q->Integral());
G->Scale(1./G->Integral());
double R=-1;
double amin=-1;
double ea=-1;
for(double alpha=0;alpha<=1.0;alpha+=x)
	{
	TH1F*mixed=(TH1F*)Q->Clone("mixed");
	mixed->Scale(alpha);
	mixed->Add(G,1-alpha);
	double r=mixed->Chi2Test(T,options);
	if((R<0)||(r<R)){R=r; amin=alpha;}
	if(graph!=NULL){graph->SetPoint( graph->GetN(), alpha,r);}	
	}

return amin;
}


double MinChiSquare::GetError(TGraph*graph)
	{
	graph->Sort();
	if(graph->GetN()<=0)return -1;
	//find min
	int iMin=0;double alpha;
	double Min;graph->GetPoint(iMin,alpha,Min);//initializing
	for(int i=0;i<graph->GetN();i++)
		{
		double r;
		graph->GetPoint(i,alpha,r);
		if(r<=Min){Min=r;iMin=i;}
		}
	for(int i=1;;i++)
		{
		double r2,r3;
		double alpha2,alpha3;
		int boundL=0,boundR=0;
		if(iMin+i<graph->GetN())graph->GetPoint(iMin+i,alpha2,r2);
			else boundR=1;
		if(iMin-i>=0)graph->GetPoint(iMin-i,alpha3,r3);
			else boundL=1;
		switch( boundL | (boundR<<1))
		{
		case 0: //no bound
			//cout<<"NO BOUND"<<endl;
			if((r2-Min)>1 && (r3-Min)>1)return fabs(alpha2-alpha3)/2.0;
			break;
		case 1:	//boundL
			//cout<<"BOUND L"<<endl;
			if((r2-Min)>1)return fabs(alpha2-alpha);
			break;
			
		case 2: //boundR
			//cout<<"BOUND R"<<endl;
			if((r3-Min)>1)return fabs(alpha3-alpha);
			break;
		case 3: return 1.0; //double bound
		}//switch	
		}//for
	//
	}

