

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCut.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TCanvas.h"

#include <string>
#include <cstdio>

#define MAXP 10

class ComputeSystMap
{

public:
	ComputeSystMap();
	~ComputeSystMap();
	void InitHisto(int nBins,float xmin,float xmax);
	void SetTree(TTree *T,int type){if(type==0)t_data=T;else if(type>0)t_mc=T;};
	int FillHisto(int type, const char *varName="QGL2012");//0 data; 1 mc
	void SetSelection(const char *selection);
	void ScaleMCtoData();
	void SetOverFlow(int i=1){overflow=i;return;}
	int GetOverFlow(){return overflow;}
	void SetNpars(int i=1){if(i<=MAXP)nPar=i;else {fprintf(stderr,"Parameters set to %d",MAXP);nPar=MAXP;}return;}
	int GetNpars(){return nPar;}
	void SetActivePar(int i=1){if(i>=MAXP){fprintf(stderr,"Nothing can be done for Par %d",i);return ;};active[i]=1;return ;}
	void UnsetActivePar(int i=1){if(i>=MAXP){fprintf(stderr,"Nothing can be done for Par %d",i);}else active[i]=0; return; }
	void MinimizePars(int nCycle=4,int depth=8,float delta=0.1,const char *varName="QGL2012");
	void Print();
	void DrawHist();
	//--- Access to sturctures
	TH1F* GetHistData(){return h_data;}
	TH1F* GetHistMC(){return h_mc;}
	TTree* GetTree(int type){if(type==0)return t_data;if(type>0)return t_mc; return NULL;}
	TEventList* GetEventList(){return elist;}
	float GetPar(int i){return pars[i];} //ADD CHECK TODO
	void SetPar(int i,float x){pars[i]=x;return;}
protected:
	TH1F* h_data;
	TH1F* h_mc;
	int overflow;
	int nPar;
	int nBins_;float xmin_;float xmax_;
	double pars[MAXP];
	short active[MAXP]; //char = 1byte, short 2, int 4, long 4, float 4, double 8;
	TEventList *elist;
	std::string Sel;
	TTree *t_data;
	TTree *t_mc;
	void PrintCount();
	int iPrint;

};

ComputeSystMap::ComputeSystMap() 
	{
	iPrint=0;
	h_data=NULL;
	h_mc=NULL;
	overflow=0;
	elist=NULL;
	for(int i=0;i<MAXP;i++) active[i]=0;
	for(int i=0;i<MAXP;i++) pars[i]=0;
	}
ComputeSystMap::~ComputeSystMap()
	{
	//if(h_data!=NULL) delete h_data;
	//if(h_mc!=NULL) delete h_mc;
	}
void ComputeSystMap::InitHisto(int nBins,float xmin, float xmax){
	nBins_=nBins;
	xmin_=xmin;
	xmax_=xmax;
	if(h_data!=NULL) delete h_data;
	if(h_mc!=NULL) delete h_mc;
	h_data=new TH1F("h_data","Data;L;events",nBins,xmin,xmax);
	h_mc=new TH1F("h_mc","MC;L;events",nBins,xmin,xmax);
	
	h_data->Sumw2();
	h_mc->Sumw2();
	return;
}

void ComputeSystMap::SetSelection(const char *selection)
	{
//	elist=new TEventList("elist","elist");
	//t_mc->Draw(">>elist",selection);
	//elist=(TEventList*)gDirectory->Get("elist");
	//t->SetEventList(elist);
	Sel=selection;
	}
int ComputeSystMap::FillHisto(int type,const char *varName)//0 data; 1 mc
	{
	using namespace std;
	TH1F *h=NULL;
	if(type==0) h=h_data;
	if(type>0) h=h_mc;

	TTree *t;		
	if(type==0) t=t_data;
	if(type>0) t=t_mc;

	if(h==NULL)return -1;
	if(t==NULL)return -1;
	
	string target;
	string var;
	if(type==0) target="h_data";
	if(type>0) target="h_mc";


	if(type==0) var=varName;
	if(type>0){var=varName;
		if(nPar>0 && active[0]) { var += Form("+( %f )",pars[0]);}
		if(nPar>1 && active[1]) { var += Form("+( %f*(%s-0.5) )",pars[1],varName);} //Not TMath -> Speed
		for(int i=2;i<nPar;i++) if(active[i]) { var += Form("+( %f*TMath::Power( %s-0.5,%d) )",pars[i],varName,i);}
		}

	//no need of selection: EVENT LIST
	//fprintf(stderr,"DEBUG: Going to Draw: ---%s--- with selection ---%s---\n",Form("%s>>%s",var.c_str(),target.c_str()),Sel.c_str() );
	int R=t->Draw(   Form("%s>>%s",var.c_str(),target.c_str()),Sel.c_str() ,"E");
	if(type==0)h_data=(TH1F*)gDirectory->Get("h_data");
	if(type>0)h_mc=(TH1F*)gDirectory->Get("h_mc");
//	fprintf(stderr,"DEBUG: %d entries\n",R);
		
	return 0;
	}
void ComputeSystMap::ScaleMCtoData(){
	if(!overflow)
		h_mc->Scale(h_data->Integral()/h_mc->Integral());
	else 
		h_mc->Scale(h_data->Integral(0,h_data->GetNbinsX()+1)/h_mc->Integral(0, h_mc->GetNbinsX()+1) );
	}
void ComputeSystMap::MinimizePars(int nCycle/*4*/,int depth/*8*/,float delta/*=0.2*/,const char *varName/*="QGL2012"*/){
for(int c=0;c<nCycle;c++){
	
	const char opt[]="CHI2 WW";
	FillHisto(0,varName);	
	FillHisto(1,varName);	
	ScaleMCtoData();
	double min=h_data->Chi2Test(h_mc,opt);
	for(int i=0;i<nPar;i++)
	{
	float d= delta*2;	
	if(!active[i])continue;
	float parmin=pars[i];
	for(int k=0;k<depth;k++)
		{
		d*=0.5;	
		bool change=0;
		parmin=pars[i];
		pars[i]+=d;	
			if(pars[i]==parmin) return; //precision check
		FillHisto(1,varName);
		ScaleMCtoData();
		PrintCount(); //
		fprintf(stderr,"min=%f parmin=%f\n",min,parmin);
		if(min>h_data->Chi2Test(h_mc,opt) ){
							min=h_data->Chi2Test(h_mc,opt);
							parmin=pars[i]; 
							change=1;
							fprintf(stderr,"change +\n");
							}
		else{pars[i]=parmin-d;
			FillHisto(1,varName);
			ScaleMCtoData();
			if(min>h_data->Chi2Test(h_mc,opt) ){min=h_data->Chi2Test(h_mc,opt);
								parmin=pars[i];change=1;
							fprintf(stderr,"change -\n");
								}
    		     }

		if(change)k--;
		else pars[i]=parmin;
		printf("k=%d\n",k);	
		}
	}
	}
}
void ComputeSystMap::Print()
{
for(int i=0;i<nPar;i++){if(!active[i])continue;printf("par %d: %f\n",i,pars[i]);}	
return;
}
void ComputeSystMap::DrawHist(){
	h_mc->SetLineColor(kRed);
	h_data->SetMarkerStyle(20);

	TCanvas *c=new TCanvas("c1","c1",800,800);
	h_mc->Draw("HIST");
	h_data->Draw("P SAME");
	return;
}
void ComputeSystMap::PrintCount()
	{
	if(iPrint<20){fprintf(stderr,".");}
	else{iPrint=0; fprintf(stderr,"\r");fprintf(stderr,"                    \r");}
	}
