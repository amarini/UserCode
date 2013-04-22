

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMath.h"
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cstdio>


using namespace std;

class Analyzer{
public:
	Analyzer(){
		PtMin=30;PtMax=80;RhoMin=0;RhoMax=15;EtaMin=0;EtaMax=2.0;
		t_mc=NULL;t_data=NULL;
		varName="QGLHisto";
		nBins=30;xMin=0;xMax=1.000001;
		nstep=10;stp0=.01;stp1=0.01;
		lmin=0;lmax=1.0;
		opt="CHI2 WW";
		aMin=0.5;aMax=1.3;bMin=0.5;bMax=1.5;
		}
	string varName;//QGL HISTO
	int nstep;
	float stp0;
	float stp1;
	
	
	float function(float x0, float a ,float b,float min=0,float max=1);
	void ComputeMin();
	void Loop(TChain *t,int type); //type|=4 : compute lmin,lmax; type|=1 data type |=2 mc
	void FillHisto(TH1F *h,string var);
	pair<float, float> MinG(TGraph2D *g);
	void SpanMin();
	void CreateHisto();
	void SetTrees(TChain *mc,TChain*data){t_data=data;t_mc=mc;}
	void Reset(TH1F *h);
	pair<float,float> SmearDoubleMin(float a0_q,float b0_q , float a0_g,float b0_g,int type); //type = 0 Q, 1 G
	void ComputeDoubleMin();
	void ComputeMinFast();

//private:
	TChain *t_mc;
	TChain *t_data;
	float PtMin;
	float PtMax;
	float RhoMin;
	float RhoMax;
	float EtaMin;
	float EtaMax;
	
	float lmin;
	float lmax;
	float a_q,a_g,b_q,b_g;
	float alpha,beta;

	float aMin;
	float aMax;
	float bMin;
	float bMax;
	
	vector<float> alphaFast;
	vector<float> betaFast;
	vector<TH1F*> h_mcFast;

	Int_t nBins;
	string opt;
	Double_t xMin,xMax;
	
//---- internal not modify
	TH1F*h_mc,*h_data;
	TGraph2D *g2;
	
	map<string,float> treeVar;
	map<string,int> treeVarInt;
};


float Analyzer::function(float x0, float a ,float b,float min,float max)
{
using namespace TMath;
float x=(x0-min)/(max-min); 
if(x<0)x=0;
if(x>1)x=1;
float x1= (TanH( a* ATanH(2*x-1)+b )/2+.5 ) ;
return x1*(max-min)+min;
}

void Analyzer::Loop(TChain *t,int type){ //type|=4 : compute lmin,lmax; type|=1 data type |=2 mc

		treeVar["betaStarJet0"]		=-999;t->SetBranchAddress("betaStarJet0",	&treeVar["betaStarJet0"]	);
		treeVarInt["nvertex"]		=-999;t->SetBranchAddress("nvertex",		&treeVarInt["nvertex"]		);
		treeVar["deltaPhi_jet"]		=-999;t->SetBranchAddress("deltaPhi_jet",	&treeVar["deltaPhi_jet"]	);
		treeVar["mZ"]			=-999;t->SetBranchAddress("mZ",			&treeVar["mZ"]			);
		treeVarInt["nPFCand_QC_ptCutJet"]=-999;t->SetBranchAddress("nPFCand_QC_ptCutJet",&treeVarInt["nPFCand_QC_ptCutJet"]			);
		treeVar["ptD_QCJet0"]			=-999;t->SetBranchAddress("ptD_QCJet0",			&treeVar["ptD_QCJet0"]			);
		treeVar["axis1_QCJet0"]			=-999;t->SetBranchAddress("axis1_QCJet0",			&treeVar["axis1_QCJet0"]			);
		treeVar["axis2_QCJet0"]			=-999;t->SetBranchAddress("axis2_QCJet0",			&treeVar["axis2_QCJet0"]			);
		treeVar["ptZ"]			=-999;t->SetBranchAddress("ptZ",		&treeVar["ptZ"]			);
		treeVar["ptJet0"]		=-999;t->SetBranchAddress("ptJet0",		&treeVar["ptJet0"]		);
		treeVar["rhoPF"]		=-999;t->SetBranchAddress("rhoPF",		&treeVar["rhoPF"]		);
		treeVar["etaJet"]		=-999;t->SetBranchAddress("etaJet0",		&treeVar["etaJet"]		);
		treeVarInt["pdgIdPartJet0"]	=-999;if(type!=1)t->SetBranchAddress("pdgIdPartJet0",	&treeVarInt["pdgIdPartJet0"]); //only mc 

		treeVar["QGLHisto"]		=-999;t->SetBranchAddress("QGLHisto",		&treeVar["QGLHisto"]		);
		treeVar["QGLMLP"]		=-999;t->SetBranchAddress("QGLMLP",		&treeVar["QGLMLP"]		);
		treeVar["QGLHistoFwd"]		=-999;if(type==1)t->SetBranchAddress("QGLHistoFwd",	&treeVar["QGLHistoFwd"]		); //only data
		treeVar["QGLMLPFwd"]		=-999;if(type==1)t->SetBranchAddress("QGLMLPFwd",	&treeVar["QGLMLPFwd"]		); //only data
		
		treeVar["PUReWeight"]		=1   ;if(type!=1)t->SetBranchAddress("PUReWeight",	&treeVar["PUReWeight"]		); //only mc
		treeVar["eventWeight"]		=1   ;if(type!=1)t->SetBranchAddress("eventWeight",	&treeVar["eventWeight"]		); //only mc
		
		if(type&4) {lmin=1.0;lmax=0;} //reset lmin-lmax

		for(int i=0;i<t->GetEntries();i++)
			{
			t->GetEntry(i);
			if( treeVarInt["nPFCand_QC_ptCutJet"] <=0 )continue;
			if((treeVar["ptJet0"]<PtMin)||(treeVar["ptJet0"]>PtMax)||(abs(treeVar["etaJet0"])<EtaMin)||(abs(treeVar["etaJet0"])>EtaMax)|| (treeVar["rhoPF"]<RhoMin)||(treeVar["rhoPF"]>RhoMax))continue;
			//Z
			if( (treeVar["mZ"]<70) || (treeVar["mZ"]>110)||(treeVar["betaStarJet0"] > 0.2 * TMath::Log( treeVarInt["nvertex"] - 0.67)) || (treeVar["deltaPhi_jet"]<3.1415-0.5) ) continue;
			if( (treeVar["axis1_QCJet0"] <=0 ) || (treeVar["axis2_QCJet0"] <=0) )continue;
			if( (treeVar["ptD_QCJet0"]<=0)) continue;
			//---------------------------
		
			if(type&1){
				string var=varName;
				if(EtaMin>2.5)var+="Fwd"; //only in data fwd
				FillHisto(h_data,var);
				}	
			if(type&2){
				//mc
				FillHisto(h_mc,varName);
				}
			if(type&4){
				string var=varName;
				if(lmin>treeVar[var]) lmin=treeVar[var];
				if(lmax<treeVar[var]) lmax=treeVar[var];
				}
			if(type&8){
				alpha=1;beta=0;
				if( treeVarInt["pdgIdPartJet0"] ==21) {alpha=a_g;beta=b_g;}
				if( fabs(treeVarInt["pdgIdPartJet0"]) < 5) {alpha=a_q;beta=b_q;}
				if( fabs(treeVarInt["pdgIdPartJet0"])== 0) {alpha=1;beta=0;}

				FillHisto(h_mc,varName);
				}
			if(type&16){
				for(int j=0;j<alphaFast.size();j++)
					{
					alpha=alphaFast[j];
					beta=betaFast[j];
					FillHisto(h_mcFast[j],varName);
					}
				}
			if(type&32){ //TODO
				for(int j=0;j<alphaFast.size();j++)
					{
					alpha=1;beta=0;
					if( treeVarInt["pdgIdPartJet0"] ==21) {alpha=a_g;beta=b_g;}
					if( fabs(treeVarInt["pdgIdPartJet0"]) < 5) {alpha=a_q;beta=b_q;}
					if( fabs(treeVarInt["pdgIdPartJet0"])== 0) {alpha=1;beta=0;}
					}
				}
			}
}

void Analyzer::FillHisto(TH1F *h,string var){
	h->Fill(  function(treeVar[var],alpha,beta,lmin,lmax ), treeVar["eventWeight"]*treeVar["PUReWeight"] );
	}

void Analyzer::CreateHisto()
	{	
	h_mc=new TH1F("hmc","hmc",nBins,xMin,xMax);
	h_data=new TH1F("hmc","hmc",nBins,xMin,xMax);
	h_mc->Sumw2();
	h_data->Sumw2();
	}
void Analyzer::Reset(TH1F *h)
	{
	for(int i=0;i<=h->GetNbinsX()+1;i++) {h->SetBinContent(i,0);h->SetBinError(i,0);}
	h->SetEntries(0);
	//h->Sumw2();
	}

pair<float, float> Analyzer::MinG(TGraph2D *g){
	pair<float,float> R(-99,-99);
	if(g==NULL) return R;
	if(g->GetN()==0)return R;
	double *x1,*y1,*z1;
	x1=g->GetX();
	y1=g->GetY();
	z1=g->GetZ();
	
	float a=z1[0];R=pair<float,float>(x1[0],y1[0]);
	for(int i=0;i<g->GetN();i++){if((z1[i]<a)||(a<0)){a=z1[i];    R=pair<float,float>(x1[i],y1[i]); }}
	
	return R;
}

void Analyzer::ComputeMin(){
	g2=new TGraph2D(); //TODO
	
	alpha=1.0;beta=0;
	Loop(t_data,1);
	if(varName=="QGLMLP")
		Loop(t_mc,4);
	//scan
	alpha=1.0;beta=0;
	for(float ai=0.7; ai<=1.1; ai+=0.02)
		{
		Reset(h_mc);	
		alpha=ai;
		Loop(t_mc,2);
		h_mc->Scale(h_data->Integral()/h_mc->Integral());
		g2->SetPoint(g2->GetN(),alpha,beta, h_data->Chi2Test(h_mc,opt.c_str())  );	
		}
	alpha=1.0;beta=0;
	for(float bi=-0.5; bi<=0.5; bi+=0.01)
		{
		Reset(h_mc);	
		beta=bi;
		Loop(t_mc,2);
		h_mc->Scale(h_data->Integral()/h_mc->Integral());
		g2->SetPoint(g2->GetN(),alpha,beta, h_data->Chi2Test(h_mc,opt.c_str())  );	
		}
	
	//Find min0;min1
	float min0=1,min1=0;
	pair<float,float> R=MinG(g2);
	min0=R.first;min1=R.second;

	for(int i=-nstep;i<=nstep;i++)
	for(int j=-nstep;j<=nstep;j++)
        	{
        	alpha=min0+i*stp0;
        	beta=min1+j*stp1;
		Loop(t_mc,2);
		h_mc->Scale(h_data->Integral()/h_mc->Integral());
		g2->SetPoint(g2->GetN(),alpha,beta, h_data->Chi2Test(h_mc,opt.c_str())  );
		}
	
	R=MinG(g2);
	printf("a=%.3f;b=%.3f;lmin=%.3f;lmax=%.3f;break;\n",R.first,R.second,lmin,lmax);
	return;
	
}

void Analyzer::SpanMin(){
	vector<pair<float,float> > PtBins;
	vector<pair<float,float> > EtaBins;
	vector<pair<float,float> > RhoBins;
	
		PtBins.push_back(  pair<float,float>(30,50) );
		PtBins.push_back(  pair<float,float>(50,80) );
		PtBins.push_back(  pair<float,float>(80,120) );
		PtBins.push_back(  pair<float,float>(120,250) );
		
		RhoBins.push_back(  pair<float,float>(0,15) );
		RhoBins.push_back(  pair<float,float>(15,40) );
		
		EtaBins.push_back(  pair<float,float>(0,2) );
		EtaBins.push_back(  pair<float,float>(3,4.7) );
	
	for ( int e=0; e< EtaBins.size();e++)
	for ( int p=0; p< PtBins.size();p++)
	for ( int r=0; r< RhoBins.size();r++)
		{
		PtMin=PtBins[p].first;PtMax=PtBins[p].second;
		EtaMin=EtaBins[e].first;EtaMax=EtaBins[e].second;
		RhoMin=RhoBins[r].first;RhoMax=RhoBins[r].second;
	
		if(g2!=NULL)delete g2;	
		int t;
		if(varName=="QGLMLP")t=2;
		if(varName=="QGLHisto")t=3;
		int bin=(p+1)+(r+1)*10+(e+1)*100 + t*1000;
		printf("case %d:",bin);
		ComputeMinFast(); //Be Fast!
		}
		
	}
void Analyzer::ComputeDoubleMin(){
	nstep=5; //otherwise too slow?
	pair<float,float> R_q,R_g;
	R_q=SmearDoubleMin(1,0,1,0,0); //
	R_g=SmearDoubleMin(R_q.first,R_q.second,1,0,1); //
	R_q=SmearDoubleMin(R_q.first,R_q.second,R_g.first,R_g.second,0); //
	R_g=SmearDoubleMin(R_q.first,R_q.second,R_g.first,R_g.second,1); //

	printf("a_q=%.3f;b_q=%.3f;a_g=%.3f;b_g=%.3f;lmin=%.3f;lmax=%.3f;break;\n",R_q.first,R_q.second,R_g.first,R_g.second,lmin,lmax);
	}
pair<float,float> Analyzer::SmearDoubleMin(float a0_q,float b0_q , float a0_g,float b0_g,int type){ //type = 0 Q, 1 G
	TGraph2D *g2_q=new TGraph2D(); 
	TGraph2D *g2_g=new TGraph2D(); 
	
	alpha=1.0;beta=0;
	Loop(t_data,1);
	if(varName=="QGLMLP")
		Loop(t_mc,4);
	//scan
	a_q=a0_q;b_q=b0_q;
	a_g=a0_g;b_g=b0_g;

	alpha=1.0;beta=0;
	for(float ai=0.7; ai<=1.1; ai+=0.02)
		{
		Reset(h_mc);	
		if(type==0)a_q=ai;
		if(type==1)a_g=ai;
		Loop(t_mc,8);
		h_mc->Scale(h_data->Integral()/h_mc->Integral());

		if(type==0)g2_q->SetPoint(g2_q->GetN(),a_q,b_q, h_data->Chi2Test(h_mc,opt.c_str())  );	
		if(type==1)g2_g->SetPoint(g2_g->GetN(),a_g,b_g, h_data->Chi2Test(h_mc,opt.c_str())  );	
		}
	alpha=1.0;beta=0;
	a_q=a0_q;b_q=b0_q;
	a_g=a0_g;b_g=b0_g;

	for(float bi=-0.5; bi<=0.5; bi+=0.01)
		{
		Reset(h_mc);	
		if(type==0)b_q=bi;
		if(type==1)b_g=bi;
		Loop(t_mc,8);
		h_mc->Scale(h_data->Integral()/h_mc->Integral());
		if(type==0)g2_q->SetPoint(g2_q->GetN(),a_q,b_q, h_data->Chi2Test(h_mc,opt.c_str())  );	
		if(type==1)g2_g->SetPoint(g2_g->GetN(),a_g,b_g, h_data->Chi2Test(h_mc,opt.c_str())  );	
		}
	
	//Find min0;min1
	float min0=1,min1=0;
	pair<float,float> R;
		if(type==0)R=MinG(g2_q);
		if(type==0)R=MinG(g2_g);
	min0=R.first;min1=R.second;

	for(int i=-nstep;i<=nstep;i++)
	for(int j=-nstep;j<=nstep;j++)
        	{
        	if(type==0)a_q=min0+i*stp0;
        	if(type==1)a_g=min0+i*stp0;
        	if(type==0)b_q=min1+j*stp1;
        	if(type==1)b_g=min1+j*stp1;
		Loop(t_mc,8);
		h_mc->Scale(h_data->Integral()/h_mc->Integral());
		if(type==0)g2_q->SetPoint(g2_q->GetN(),a_q,b_q, h_data->Chi2Test(h_mc,opt.c_str())  );	
		if(type==1)g2_g->SetPoint(g2_g->GetN(),a_g,b_g, h_data->Chi2Test(h_mc,opt.c_str())  );	
		}
	
		if(type==0)R=MinG(g2_q);
		if(type==0)R=MinG(g2_g);
	//SAME ON G,& REDO
	//printf("a=%.3f;b=%.3f;lmin=%.3f;lmax=%.3f;break;\n",R.first,R.second,lmin,lmax);
	return R;
}
void Analyzer::ComputeMinFast(){
	g2=new TGraph2D(); //TODO
	
	alpha=1.0;beta=0;
	Loop(t_data,1);
	if(varName=="QGLMLP")
		Loop(t_mc,4);
	//scan
	alpha=1.0;beta=0;
	for(float ai=aMin; ai<=aMax; ai+=0.02)
		{
		Reset(h_mc);	
		alphaFast.push_back(ai);	
		betaFast.push_back(0);	
		h_mcFast.push_back(new TH1F( Form("hmc_%d",int(h_mcFast.size())),"hmc",nBins,xMin,xMax) );
		}
	alpha=1.0;beta=0;
	for(float bi=bMin; bi<=bMax; bi+=0.01)
		{
		Reset(h_mc);	
		alphaFast.push_back(1.0);	
		betaFast.push_back(bi);	
		h_mcFast.push_back(new TH1F( Form("hmc_%d",int(h_mcFast.size())),"hmc",nBins,xMin,xMax) );
		}
	
	Loop(t_mc,16);
		for(int i=0 ;i<alphaFast.size();i++)
			{
			h_mcFast[i]->Scale(h_data->Integral()/h_mcFast[i]->Integral());
			g2->SetPoint(g2->GetN(),alphaFast[i],betaFast[i], h_data->Chi2Test(h_mcFast[i],opt.c_str())  );
			}
	//Find min0;min1
	float min0=1,min1=0;
	pair<float,float> R=MinG(g2);
	min0=R.first;min1=R.second;
		
	alphaFast.clear();
	betaFast.clear();
	for(int i=0;i<h_mcFast.size();i++){delete h_mcFast[i]; };
	h_mcFast.clear();//	

	delete g2;
	g2=new TGraph2D();

	for(int i=-nstep;i<=nstep;i++)
	for(int j=-nstep;j<=nstep;j++)
        	{
        	alphaFast.push_back(min0+i*stp0 );
        	betaFast.push_back(min1+j*stp1 );
		h_mcFast.push_back(new TH1F( Form("hmc_%d",int(h_mcFast.size())),"hmc",nBins,xMin,xMax) );
		//g2->SetPoint(g2->GetN(),alpha,beta, h_data->Chi2Test(h_mc,opt.c_str())  );
		}
	
	Loop(t_mc,16);
		for(int i=0 ;i<alphaFast.size();i++)
			{
			h_mcFast[i]->Scale(h_data->Integral()/h_mcFast[i]->Integral());
			g2->SetPoint(g2->GetN(),alphaFast[i],betaFast[i], h_data->Chi2Test(h_mcFast[i],opt.c_str())  );
			}
	R=MinG(g2);
	printf("a=%.3f;b=%.3f;lmin=%.3f;lmax=%.3f;break;\n",R.first,R.second,lmin,lmax);
	return;
}

int ComputeDoubleMin(){
	TChain *mc=new TChain("tree_passedEvents");
	TChain *data=new TChain("tree_passedEvents");
		data->Add("/Users/andreamarini/Documents/QGDiscriminator/ZJet/ZJet_DoubleMu*.root");
		data->Add("/Users/andreamarini/Documents/QGDiscriminator/ZJet/ZJet_DoubleE*.root");

		mc->Add("/Users/andreamarini/Documents/QGDiscriminator/ZJet/ZJet_DYJetsToLL_M-50*.root");
	Analyzer A;
	A.nstep=10;
	//A.varName="QGLHisto";
	A.varName="QGLMLP";
	A.CreateHisto();
	A.SetTrees(mc,data);
		freopen("/dev/null","w",stderr);
	A.ComputeMinFast(); //A.ComputDoubleMin;
		/*
		A.alpha=1;
		A.beta=0;
		A.Loop(mc,2)
		A.Loop(data,1)
		*/
	A.SpanMin();
	return 0;
	}

