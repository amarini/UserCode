

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TAxis.h"
#include <map>
#include <vector>
#include <algorithm>

#include "TLegend.h"

using namespace std;


TGraph *FindMinimum(TGraph2D*g, double *returnMin=NULL){
	TGraph *R=new TGraph();
	R->SetName("min");
	R->SetTitle("Min Point");
	
	if(g==NULL)printf("ERROR\n");
	
	Double_t *x=g->GetX(), *y=g->GetY(), *z=g->GetZ();
	
	Double_t min=-1;float x0,y0;int k0;	
	for(int i=0;i<g->GetN();i++)
		{
		if((min<0) || (min> z[i]))
			{
			k0=i;
			min=z[i];
			x0=x[i];
			y0=y[i];
			}
		}
	
	R->SetPoint(0,x0,y0);
		R->SetMarkerStyle(20);
	if(returnMin!=NULL)(*returnMin)=min;
	return R;
}

TGraph *FindContour(TGraph2D *g,double min,double Delta=1.0){
	map< pair<float,float>,int> A;
	TGraph *R=new TGraph();
	R->SetName("cont");
	R->SetTitle("Cont Point");
	Double_t *x=g->GetX(),*y=g->GetY(),*z=g->GetZ();
	for(int i=0;i<g->GetN();i++)
		{
		if( (z[i] < min+Delta) && ( A[ pair<float,float>(x[i],y[i]) ] ==0 ) )
			{
			A[pair<float,float>(x[i],y[i])]++;
			R->SetPoint(R->GetN(),x[i],y[i]);
			}
		}
	return R;
	
}

void FindDrawMin(TGraph*g,double &xMin,double&xMax,double&yMin,double&yMax)
{
	Double_t *x=g->GetX(),*y=g->GetY();
	xMin=x[0];xMax=x[0];yMin=y[0];yMax=y[0];
	for(int i=0;i<g->GetN();i++)
		{
		xMax = (xMax < x[i])?x[i]:xMax;
		yMax = (yMax < y[i])?y[i]:yMax;

		xMin = (xMin > x[i])?x[i]:xMin;
		yMin = (yMin > y[i])?y[i]:yMin;
		}
	return;
}

void InsertPoint(TGraph*g,int i,double x0, double y0)
	{
	if(g->GetN()<i)return;
	Double_t *x=g->GetX(),*y=g->GetY();
	int N=g->GetN();
	for(int k=N;k>i;k--)
		{
		g->SetPoint(k,x[k-1],y[k-1]);
		}
	g->SetPoint(i,x0,y0);
	}

int RemovePoints(TGraph *t)
{
	//prova a rimuovere punti finche puoi
	float area=0;
	Double_t *xt=t->GetX(),*yt=t->GetY();
	for(int i=0;i<t->GetN();i++)
		{
		area=t->Integral();
		double x0=xt[i],y0=yt[i];	
		t->RemovePoint(i);
		if( (area >= t->Integral()) || (area <0.0001)) InsertPoint(t,i,x0,y0);
		else i=0;
		}
}

int ReorderGraph(TGraph*g){
	float area=0;
	TGraph *t=new TGraph();
	Double_t *x=g->GetX(),*y=g->GetY();
	
	t->SetPoint(t->GetN(),x[0],y[0]);
	t->SetPoint(t->GetN(),x[1],y[1]);
	t->SetPoint(t->GetN(),x[2],y[2]);

	for(int i=3;i<g->GetN();i++)
		{
		area=t->Integral();
		int k0=-1;
		for(int k=0;k<t->GetN();k++)
			{
			InsertPoint(t,k,x[i],y[i]);
			if(t->Integral()>=area){ area=t->Integral(); k0=k; }
			t->RemovePoint(k);
			}
		//trovo il punto che massimizza l'area e lo inserisco li - se esiste	
		if(k0>=0)InsertPoint(t,k0,x[i],y[i]);
		//provo a togliere i vari punti
		RemovePoints(t);
		//
		}
	//copy to the original
	Double_t *xt=t->GetX(),*yt=t->GetY();
	for(int i=0;i<t->GetN();i++)
		g->SetPoint(i,xt[i],yt[i]);
	for(int i=g->GetN();i>=t->GetN();i--)
		g->RemovePoint(i);
	}


int DrawParameter(const char *varName="QGLHisto",float PtMin=80,float PtMax=120,float RhoMin=0,float RhoMax=15, float EtaMin=0,float EtaMax=2.0,float Delta=4.7){

//output_QGLHisto_pt120_250_rho0_15_eta0_2g2_g.root   output_QGLMLP_pt120_250_rho0_15_eta0_2g2_q.root
//output_QGLHisto_pt120_250_rho0_15_eta0_2g2_q.root 

float Delta2=Delta*4;

TFile  *fq =TFile::Open( Form("Results/output_%s_pt%.0f_%.0f_rho%.0f_%.0f_eta%.0f_%.0fg2_q.root",varName,PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax) );
TFile  *fg =TFile::Open( Form("Results/output_%s_pt%.0f_%.0f_rho%.0f_%.0f_eta%.0f_%.0fg2_g.root",varName,PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax) );


TGraph2D *gq=(TGraph2D*)fq->Get("g2_q");
TGraph2D *gg=(TGraph2D*)fg->Get("g2_g");

TGraph *a=new TGraph(); a->SetName("axis");
a->SetPoint(0,-10,-10);
a->SetPoint(1,10,10);


double min;
TGraph *min_q=FindMinimum(gq,&min);min_q->SetName("min_q"); TGraph *cont_q=FindContour(gq,min,Delta); cont_q->SetName("cont_q");TGraph *cont2_q=FindContour(gq,min,Delta2); cont2_q->SetName("cont2_q");
TGraph *min_g=FindMinimum(gg,&min);min_g->SetName("min_g"); TGraph *cont_g=FindContour(gg,min,Delta); cont_g->SetName("cont_g");TGraph *cont2_g=FindContour(gg,min,Delta2); cont2_g->SetName("cont2_g");


//Style
min_q->SetMarkerColor(kBlue-4);
min_g->SetMarkerColor(kRed+2);

cont_q->SetFillColor(kBlue-7);
cont_g->SetFillColor(kRed-7);

cont_q->SetMarkerStyle(20);cont_q->SetMarkerColor(kBlue-7);
cont_g->SetMarkerStyle(20);cont_g->SetMarkerColor(kRed-7);

cont2_q->SetFillColor(kBlue-9);
cont2_g->SetFillColor(kRed-9);

cont2_q->SetMarkerStyle(20);cont2_q->SetMarkerColor(kBlue-9);
cont2_g->SetMarkerStyle(20);cont2_g->SetMarkerColor(kRed-9);
//ReorderGraph(cont_q);
//ReorderGraph(cont_g);

//cont_q->SetFillStyle(3004);
//cont_g->SetFillStyle(3005);

//Line
TGraph *l1=new TGraph();l1->SetName("l1");
TGraph *l2=new TGraph();l2->SetName("l2");

l1->SetPoint(0,-5,0);l1->SetPoint(1,5,0);
l2->SetPoint(0,1,-5);l2->SetPoint(1,1,5);


//Find Ranges
double xMin,xMax,yMin,yMax;
double xMin2,xMax2,yMin2,yMax2;
FindDrawMin(cont2_q,xMin,xMax,yMin,yMax);
FindDrawMin(cont2_g,xMin2,xMax2,yMin2,yMax2);
xMin=(xMin<xMin2)?xMin:xMin2;
yMin=(yMin<yMin2)?yMin:yMin2;
xMax=(xMax>xMax2)?xMax:xMax2;
yMax=(yMax>yMax2)?yMax:yMax2;

xMin2=1;yMin2=0; xMax2=1;yMax2=0;
xMin=(xMin<xMin2)?xMin:xMin2;yMin=(yMin<yMin2)?yMin:yMin2;xMax=(xMax>xMax2)?xMax:xMax2;yMax=(yMax>yMax2)?yMax:yMax2;

xMin+=-.05;
xMax+=.05;
yMin+=-.05;
yMax+=.05;

//Draw
TCanvas *c=new TCanvas("c1","c1",800,800);
a->Draw("AP");
//a->GetXaxis()->SetRangeUser(0,2);
//a->GetYaxis()->SetRangeUser(-1,1);
a->GetXaxis()->SetRangeUser(xMin,xMax);
a->GetYaxis()->SetRangeUser(yMin,yMax);

printf("COUNT=%d\n",cont_q->GetN());
cont2_q->Draw("P SAME ");
cont2_g->Draw("P SAME ");


cont_q->Draw("P SAME ");
cont_g->Draw("P SAME ");

l1->Draw("L SAME");
l2->Draw("L SAME");

min_q->Draw("P SAME");
min_g->Draw("P SAME");


TLegend *L=new TLegend(0.7,0.7,0.89,.89,Form("#splitline{%s %.0f<P_{T}[GeV]<%.0f}{%.0f<#rho<%.0f  %.0f<#eta<%.0f}",varName,PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax));
L->SetBorderSize(0);
L->SetFillStyle(0);
L->AddEntry(min_q,"Quark","LP");
L->AddEntry(min_g,"Gluon","LP");
L->AddEntry(cont_q,Form("(min+%.0f)",Delta),"F");
L->AddEntry(cont_g,Form("(min+%.0f)",Delta),"F");
L->AddEntry(cont2_q,Form("(min+%.0f)",Delta2),"F");
L->AddEntry(cont2_g,Form("(min+%.0f)",Delta2),"F");
L->Draw();

c->SaveAs(Form("Results/Plot/Parameters_%s_pt%.0f_%.0f_rho%.0f_%.0f_eta%.0f_%.0f.pdf",varName,PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax));

}
