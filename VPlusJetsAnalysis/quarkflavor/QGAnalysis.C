#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "THStack.h"

#include "TMath.h"

#include <iostream>
#include <cstdio>
#include <cstdlib>

//Roofit
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"

#ifdef STANDALONE
#include "src/ReadParameters.C"
#endif

using namespace std;
using namespace RooFit;
int CHID2=-4;

//TFile
TFile *fdata_4;
TFile *fdata_1;
TFile *fdata_2;
TFile *fDY_4;
TFile *fTT_4;
TFile *fWJ_4;
TFile *fWW_4;
TFile *fWZ_4;
TFile *fZZ_4;

TFile *fDY_1;
TFile *fTT_1;
TFile *fWJ_1;
TFile *fWW_1;
TFile *fWZ_1;
TFile *fZZ_1;

TFile *fDY_2;
TFile *fTT_2;
TFile *fWJ_2;
TFile *fWW_2;
TFile *fWZ_2;
TFile *fZZ_2;
//----
TFile *fout;

const bool DataBKG=true;

map<string,TH1F*> histos;


int QGFit(map<string,TH1F*> &histos,map<string,vector<float>* > &frac,float lumi,int WriteResults=1);

int QGAnalysis(float lumi=18.7,const char *OutFile="")
{ // Main function: open files, call fit, call subfunctions, writes results.
fout=TFile::Open( OutFile ,"RECREATE");
//Take llPt
TH1F* llPt_data=(TH1F*)fdata_4->Get("llPt")->Clone("llPt_data");
	llPt_data->Add((TH1F*)fdata_1->Get("llPt"));
TH1F* llPt_DY=(TH1F*)fDY_4->Get("llPt")->Clone("llPt_DY");
	llPt_DY->Add((TH1F*)fDY_1->Get("llPt"));
TH1F* llPt_TT=(TH1F*)fTT_4->Get("llPt")->Clone("llPt_TT");
	llPt_TT->Add((TH1F*)fTT_1->Get("llPt"));
TH1F* llPt_WJ=(TH1F*)fWJ_4->Get("llPt")->Clone("llPt_WJ");
	llPt_WJ->Add((TH1F*)fWJ_1->Get("llPt"));
TH1F* llPt_WW=(TH1F*)fWW_4->Get("llPt")->Clone("llPt_WW");
	llPt_WW->Add((TH1F*)fWW_1->Get("llPt"));
TH1F* llPt_WZ=(TH1F*)fWZ_4->Get("llPt")->Clone("llPt_WZ");
	llPt_WZ->Add((TH1F*)fWZ_1->Get("llPt"));
TH1F* llPt_ZZ=(TH1F*)fZZ_4->Get("llPt")->Clone("llPt_ZZ");
	llPt_ZZ->Add((TH1F*)fZZ_1->Get("llPt"));

//Scale MC to luminosity
cout<<"LUMI="<<lumi<<"1/fb"<<endl;
llPt_DY->Scale(lumi);
llPt_TT->Scale(lumi);
llPt_WJ->Scale(lumi);
llPt_WW->Scale(lumi);
llPt_WZ->Scale(lumi);
llPt_ZZ->Scale(lumi);

if(!DataBKG){
//Subtract BKG from data
cout<< "BKG Subtraction from MC"<<endl;
llPt_data->Add(llPt_TT,-1.0);
llPt_data->Add(llPt_WJ,-1.0);
llPt_data->Add(llPt_WW,-1.0);
llPt_data->Add(llPt_WZ,-1.0);
llPt_data->Add(llPt_ZZ,-1.0);
}else {
cout<< "BKG Subtraction from Data"<<endl;
TH1F* llPt_data_2=(TH1F*)fdata_2->Get("llPt")->Clone("llPt_data_2");
 llPt_data->Add(llPt_data_2,-1.0);
 llPt_data->Add(llPt_WJ,-1.0);
 llPt_data->Add(llPt_WW,-1.0);
 llPt_data->Add(llPt_WZ,-1.0);
 llPt_data->Add(llPt_ZZ,-1.0);
}

histos["llPt_DY"]=llPt_DY;
histos["llPt_TT"]=llPt_TT;
histos["llPt_WJ"]=llPt_WJ;
histos["llPt_WW"]=llPt_WW;
histos["llPt_WZ"]=llPt_WZ;
histos["llPt_ZZ"]=llPt_ZZ;
histos["llPt_data"]=llPt_data;

for(int bin=0; bin<=llPt_data->GetNbinsX()+1;bin++)
	{
	// Import data Histograms
	TH1F* qgl_llPt_data =(TH1F*) fdata_4->Get(  Form("qgl_llPt_bin%d",bin ) ) ;
		qgl_llPt_data->Add( (TH1F*)fdata_1->Get(Form("qgl_llPt_bin%d",bin)));
	TH1F* btag_llPt_data=(TH1F*) fdata_4->Get(  Form("btag_llPt_bin%d",bin ) );
		qgl_llPt_data->Add((TH1F*)fdata_1->Get(Form("qgl_llPt_bin%d",bin)));
	//BKG Sub
	if(DataBKG){
		qgl_llPt_data->Add((TH1F*)fdata_2->Get(Form("qgl_llPt_bin%d",bin)));
		btag_llPt_data->Add((TH1F*)fdata_2->Get(Form("btag_llPt_bin%d",bin)));
	}

	// Import MC histogram for flavor
	   //QGL
	   TH1F* qgl_llPt_q=(TH1F*) fDY_4->Get(  Form("qgl_llPt_bin%d_flavor%d",bin,1 ) )->Clone( Form("qgl_llPt_bin%d_q",bin) )  ;
	  		qgl_llPt_q->Add((TH1F*) fDY_4->Get(  Form("qgl_llPt_bin%d_flavor%d",bin,2) ) );
	  		qgl_llPt_q->Add((TH1F*) fDY_4->Get(  Form("qgl_llPt_bin%d_flavor%d",bin,3) ) );
	   TH1F* qgl_llPt_g=(TH1F*) fDY_4->Get(  Form("qgl_llPt_bin%d_flavor%d",bin,21 ))->Clone( Form("qgl_llPt_bin%d_g",bin) ); 
	   TH1F* qgl_llPt_c=(TH1F*) fDY_4->Get(  Form("qgl_llPt_bin%d_flavor%d",bin,4 ))->Clone( Form("qgl_llPt_bin%d_c",bin) ); 
	   TH1F* qgl_llPt_b=(TH1F*) fDY_4->Get(  Form("qgl_llPt_bin%d_flavor%d",bin,5 ))->Clone( Form("qgl_llPt_bin%d_b",bin) ); 
	   //BTAG
	   TH1F* btag_llPt_q=(TH1F*) fDY_4->Get(  Form("btag_llPt_bin%d_flavor%d",bin,1 ))->Clone( Form("btag_llPt_bin%d_q",bin) ); 
	  		btag_llPt_q->Add((TH1F*) fDY_4->Get(  Form("btag_llPt_bin%d_flavor%d",bin,2) ) );
	  		btag_llPt_q->Add((TH1F*) fDY_4->Get(  Form("btag_llPt_bin%d_flavor%d",bin,3) ) );
	   TH1F* btag_llPt_g=(TH1F*) fDY_4->Get(  Form("btag_llPt_bin%d_flavor%d",bin,21 ))->Clone( Form("btag_llPt_bin%d_g",bin) ); 
	   TH1F* btag_llPt_c=(TH1F*) fDY_4->Get(  Form("btag_llPt_bin%d_flavor%d",bin,4  ))->Clone( Form("btag_llPt_bin%d_c",bin) ); 
	   TH1F* btag_llPt_b=(TH1F*) fDY_4->Get(  Form("btag_llPt_bin%d_flavor%d",bin,5  ))->Clone( Form("btag_llPt_bin%d_b",bin) ); 
	
	//Electrons

	
	int nbinsIntegration=5; //template done for mc in +- nbI
	cout << "BIN INTEGRATION="<<nbinsIntegration<<endl;
	  for(int i=-nbinsIntegration;i<=nbinsIntegration;i++)
		{
		if(bin+i<0)continue;
		if(bin+i>llPt_data->GetNbinsX()+1)continue;

			qgl_llPt_q->Add((TH1F*) fDY_1->Get(  Form("qgl_llPt_bin%d_flavor%d",bin+i,1) ) );
			qgl_llPt_q->Add((TH1F*) fDY_1->Get(  Form("qgl_llPt_bin%d_flavor%d",bin+i,2) ) );
	  		qgl_llPt_q->Add((TH1F*) fDY_1->Get(  Form("qgl_llPt_bin%d_flavor%d",bin+i,3) ) );
	  		qgl_llPt_c->Add((TH1F*) fDY_1->Get(  Form("qgl_llPt_bin%d_flavor%d",bin+i,4) ) );
	  		qgl_llPt_b->Add((TH1F*) fDY_1->Get(  Form("qgl_llPt_bin%d_flavor%d",bin+i,5) ) );
	  		qgl_llPt_g->Add((TH1F*) fDY_1->Get(  Form("qgl_llPt_bin%d_flavor%d",bin+i,21) ) );
	
			btag_llPt_q->Add((TH1F*) fDY_1->Get(  Form("btag_llPt_bin%d_flavor%d",bin+i,1) ) );
			btag_llPt_q->Add((TH1F*) fDY_1->Get(  Form("btag_llPt_bin%d_flavor%d",bin+i,2) ) );
	  		btag_llPt_q->Add((TH1F*) fDY_1->Get(  Form("btag_llPt_bin%d_flavor%d",bin+i,3) ) );
	  		btag_llPt_c->Add((TH1F*) fDY_1->Get(  Form("btag_llPt_bin%d_flavor%d",bin+i,4) ) );
	  		btag_llPt_b->Add((TH1F*) fDY_1->Get(  Form("btag_llPt_bin%d_flavor%d",bin+i,5) ) );
	  		btag_llPt_g->Add((TH1F*) fDY_1->Get(  Form("btag_llPt_bin%d_flavor%d",bin+i,21) ) );

		if(i==0)continue; //0 already exists -- for muons

			qgl_llPt_q->Add((TH1F*) fDY_4->Get(  Form("qgl_llPt_bin%d_flavor%d",bin+i,1) ) );
			qgl_llPt_q->Add((TH1F*) fDY_4->Get(  Form("qgl_llPt_bin%d_flavor%d",bin+i,2) ) );
	  		qgl_llPt_q->Add((TH1F*) fDY_4->Get(  Form("qgl_llPt_bin%d_flavor%d",bin+i,3) ) );
	  		qgl_llPt_c->Add((TH1F*) fDY_4->Get(  Form("qgl_llPt_bin%d_flavor%d",bin+i,4) ) );
	  		qgl_llPt_b->Add((TH1F*) fDY_4->Get(  Form("qgl_llPt_bin%d_flavor%d",bin+i,5) ) );
	  		qgl_llPt_g->Add((TH1F*) fDY_4->Get(  Form("qgl_llPt_bin%d_flavor%d",bin+i,21) ) );
	
			btag_llPt_q->Add((TH1F*) fDY_4->Get(  Form("btag_llPt_bin%d_flavor%d",bin+i,1) ) );
			btag_llPt_q->Add((TH1F*) fDY_4->Get(  Form("btag_llPt_bin%d_flavor%d",bin+i,2) ) );
	  		btag_llPt_q->Add((TH1F*) fDY_4->Get(  Form("btag_llPt_bin%d_flavor%d",bin+i,3) ) );
	  		btag_llPt_c->Add((TH1F*) fDY_4->Get(  Form("btag_llPt_bin%d_flavor%d",bin+i,4) ) );
	  		btag_llPt_b->Add((TH1F*) fDY_4->Get(  Form("btag_llPt_bin%d_flavor%d",bin+i,5) ) );
	  		btag_llPt_g->Add((TH1F*) fDY_4->Get(  Form("btag_llPt_bin%d_flavor%d",bin+i,21) ) );
	
		}

	//Normalize
	qgl_llPt_q->Scale(1./qgl_llPt_q->Integral());
	qgl_llPt_g->Scale(1./qgl_llPt_g->Integral());
	qgl_llPt_c->Scale(1./qgl_llPt_c->Integral());
	qgl_llPt_b->Scale(1./qgl_llPt_b->Integral());

	btag_llPt_q->Scale(1./btag_llPt_q->Integral());
	btag_llPt_g->Scale(1./btag_llPt_g->Integral());
	btag_llPt_c->Scale(1./btag_llPt_c->Integral());
	btag_llPt_b->Scale(1./btag_llPt_b->Integral());
	
	qgl_llPt_data->Scale(1./qgl_llPt_data->Integral());
	btag_llPt_data->Scale(1./btag_llPt_data->Integral());
	
	histos[Form("qgl_llPt_bin%d_flavor_q",bin)]=qgl_llPt_q;
	histos[Form("qgl_llPt_bin%d_flavor_g",bin)]=qgl_llPt_g;
	histos[Form("qgl_llPt_bin%d_flavor_c",bin)]=qgl_llPt_c;
	histos[Form("qgl_llPt_bin%d_flavor_b",bin)]=qgl_llPt_b;
	histos[Form("qgl_llPt_bin%d_data",bin)]=qgl_llPt_data;

	histos[Form("btag_llPt_bin%d_flavor_q",bin)]=btag_llPt_q;
	histos[Form("btag_llPt_bin%d_flavor_g",bin)]=btag_llPt_g;
	histos[Form("btag_llPt_bin%d_flavor_c",bin)]=btag_llPt_c;
	histos[Form("btag_llPt_bin%d_flavor_b",bin)]=btag_llPt_b;
	histos[Form("btag_llPt_bin%d_data",bin)]=btag_llPt_data;

	cout<<"TODO: BKG Sub in Fit"<<endl;
	}//loop on llPt_bin

vector<float> *q_frac=new vector<float>;
vector<float> *g_frac=new vector<float>;
vector<float> *c_frac=new vector<float>;
vector<float> *b_frac=new vector<float>;

vector<float> *q_err =new vector<float>;
vector<float> *g_err =new vector<float>;
vector<float> *c_err =new vector<float>;
vector<float> *b_err =new vector<float>;

map<string,vector<float>* > frac;

frac["q_frac"]=q_frac;
frac["g_frac"]=g_frac;
frac["c_frac"]=c_frac;
frac["b_frac"]=b_frac;
frac["q_err"]=q_err;
frac["g_err"]=g_err;
frac["c_err"]=c_err;
frac["b_err"]=b_err;
//----------------------------------------------------------FIT------------------------
 QGFit(histos,frac,lumi);
//-------------------------------------------------------------------------------------


//	if(WriteResults)
	{	
	TCanvas *c=new TCanvas("llPtComp","Composition",800,800);
	c->cd();
	llPt_data->SetMarkerStyle(20);llPt_data->SetMarkerColor(kBlack);llPt_data->SetMarkerSize(1.0);llPt_data->Draw("P");
	TH1F* llPt_q=(TH1F*)llPt_data->Clone("llPt_q");
	TH1F* llPt_g=(TH1F*)llPt_data->Clone("llPt_g");
	TH1F* llPt_c=(TH1F*)llPt_data->Clone("llPt_c");
	TH1F* llPt_b=(TH1F*)llPt_data->Clone("llPt_b");
	for(int bin=0;bin<=llPt_data->GetNbinsX()+1;bin++)
		{
		llPt_q->SetBinContent(bin,llPt_q->GetBinContent(bin)* (*q_frac)[bin]);
			llPt_q->SetBinError  (bin,llPt_q->GetBinContent(bin)*(*q_err)[bin]);
		llPt_g->SetBinContent(bin,llPt_g->GetBinContent(bin)*(*g_frac)[bin]);
			llPt_g->SetBinError  (bin,llPt_g->GetBinContent(bin)*(*g_err)[bin]);
		llPt_c->SetBinContent(bin,llPt_c->GetBinContent(bin)*(*c_frac)[bin]);
			llPt_c->SetBinError  (bin,llPt_c->GetBinContent(bin)*(*c_err)[bin]);
		llPt_b->SetBinContent(bin,llPt_b->GetBinContent(bin)*(*b_frac)[bin]);
			llPt_q->SetBinError  (bin,llPt_q->GetBinContent(bin)*(*b_err)[bin]);
		}
	THStack *S=new THStack("stack","stack");
		llPt_b->SetMarkerStyle(24);llPt_b->SetMarkerColor(kGreen+2) ;llPt_b->SetMarkerSize(0.8);S->Add(llPt_b);
		llPt_c->SetMarkerStyle(24);llPt_c->SetMarkerColor(kYellow+2);llPt_c->SetMarkerSize(0.8);S->Add(llPt_c);
		llPt_g->SetMarkerStyle(24);llPt_g->SetMarkerColor(kRed+2)   ;llPt_g->SetMarkerSize(0.8);S->Add(llPt_g);
		llPt_q->SetMarkerStyle(24);llPt_q->SetMarkerColor(kBlue+2)  ;llPt_q->SetMarkerSize(0.8);S->Add(llPt_q);
	S->Draw("SAME");
	c->Write();
	TCanvas *c2=new TCanvas("llPtComp_ratio","Composition Ratio",800,800);
		TH1F* r_q=(TH1F*)llPt_q->Clone("r_q");	
		TH1F* r_g=(TH1F*)llPt_g->Clone("r_g");	
		TH1F* r_c=(TH1F*)llPt_g->Clone("r_c");	
		TH1F* r_b=(TH1F*)llPt_g->Clone("r_b");	
	for(int bin=0;bin<=llPt_data->GetNbinsX()+1;bin++)
		{
		r_q->SetBinContent(bin,(*q_frac)[bin]);
		r_g->SetBinContent(bin,(*g_frac)[bin]);
		r_b->SetBinContent(bin,(*b_frac)[bin]);
		r_c->SetBinContent(bin,(*c_frac)[bin]);

		r_q->SetBinError(bin,(*q_err)[bin]);
		r_g->SetBinError(bin,(*g_err)[bin]);
		r_b->SetBinError(bin,(*b_err)[bin]);
		r_c->SetBinError(bin,(*c_err)[bin]);

		}
		r_q->Draw("P");
		r_g->Draw("P SAME");
		r_c->Draw("P SAME");
		r_b->Draw("P SAME");
	c2->Write();
	//Write Histos in the OutFile
	llPt_data->Write();
	llPt_q->Write();
	llPt_g->Write();
	llPt_c->Write();
	llPt_b->Write();

	r_q->Write();
	r_g->Write();
	r_c->Write();
	r_b->Write();

	//Get MC and Write in the Analysis file
	TH1F* llPt_DY_u =(TH1F*) fDY_4->Get(      Form("llPt_flavor%d",0 )) ->Clone( "llPt_DY_u");
		llPt_DY_u->Add( (TH1F*)fDY_1->Get(Form("llPt_flavor%d",0)));

	TH1F* llPt_DY_q =(TH1F*) fDY_4->Get(      Form("llPt_flavor%d",1 ) ) ->Clone( "llPt_DY_q");
		llPt_DY_q->Add( (TH1F*)fDY_4->Get(Form("llPt_flavor%d",2)));
		llPt_DY_q->Add( (TH1F*)fDY_4->Get(Form("llPt_flavor%d",3)));
		llPt_DY_q->Add( (TH1F*)fDY_1->Get(Form("llPt_flavor%d",1)));
		llPt_DY_q->Add( (TH1F*)fDY_1->Get(Form("llPt_flavor%d",2)));
		llPt_DY_q->Add( (TH1F*)fDY_1->Get(Form("llPt_flavor%d",3)));
	TH1F* llPt_DY_g =(TH1F*) fDY_4->Get(      Form("llPt_flavor%d",21 )) ->Clone( "llPt_DY_g");
		llPt_DY_g->Add( (TH1F*)fDY_1->Get(Form("llPt_flavor%d",21)));
	TH1F* llPt_DY_c =(TH1F*) fDY_4->Get(      Form("llPt_flavor%d",4 ) ) ->Clone( "llPt_DY_c");
		llPt_DY_c->Add( (TH1F*)fDY_1->Get(Form("llPt_flavor%d",4)));
	TH1F* llPt_DY_b =(TH1F*) fDY_4->Get(      Form("llPt_flavor%d",5 ) ) ->Clone( "llPt_DY_b");
		llPt_DY_b->Add( (TH1F*)fDY_1->Get(Form("llPt_flavor%d",5)));
	
	llPt_DY_u->Scale(lumi);llPt_DY_u->Write();
	llPt_DY_q->Scale(lumi);llPt_DY_q->Write();
	llPt_DY_g->Scale(lumi);llPt_DY_g->Write();
	llPt_DY_c->Scale(lumi);llPt_DY_c->Write();
	llPt_DY_b->Scale(lumi);llPt_DY_b->Write();
	}


}//end QGAnalysis


int QGFit(map<string,TH1F*> &histos,map<string,vector<float>* > &frac,float lumi,int WriteResults)
{
TH1F* llPt_DY=histos["llPt_DY"];
TH1F* llPt_TT=histos["llPt_TT"];
TH1F* llPt_WJ=histos["llPt_WJ"];
TH1F* llPt_WW=histos["llPt_WW"];
TH1F* llPt_WZ=histos["llPt_WZ"];
TH1F* llPt_ZZ=histos["llPt_ZZ"];
TH1F* llPt_data=histos["llPt_data"];

vector<float> *q_frac = frac["q_frac"];
vector<float> *g_frac = frac["g_frac"];
vector<float> *c_frac = frac["c_frac"];
vector<float> *b_frac = frac["b_frac"];
vector<float> *q_err = frac["q_err"];
vector<float> *g_err = frac["g_err"];
vector<float> *c_err = frac["c_err"];
vector<float> *b_err = frac["b_err"];

q_frac->clear();q_frac->resize(llPt_data->GetNbinsX()+1);
g_frac->clear();g_frac->resize(llPt_data->GetNbinsX()+1);
c_frac->clear();c_frac->resize(llPt_data->GetNbinsX()+1);
b_frac->clear();b_frac->resize(llPt_data->GetNbinsX()+1);

q_err->clear(); q_err->resize(llPt_data->GetNbinsX()+1);
g_err->clear(); g_err->resize(llPt_data->GetNbinsX()+1);
c_err->clear(); c_err->resize(llPt_data->GetNbinsX()+1);
b_err->clear(); b_err->resize(llPt_data->GetNbinsX()+1);

//create a variable for the likelihood x
RooRealVar l("l","l",0,1.00001) ;
//create a variable for the btag b
RooRealVar b("b","b",0,1.00001) ;



for(int bin=0; bin<=llPt_data->GetNbinsX()+1;bin++)
{
	TH1F* qgl_llPt_q=histos[Form("qgl_llPt_bin%d_flavor_q",bin)];
	TH1F* qgl_llPt_g=histos[Form("qgl_llPt_bin%d_flavor_g",bin)];
	TH1F* qgl_llPt_c=histos[Form("qgl_llPt_bin%d_flavor_c",bin)];
	TH1F* qgl_llPt_b=histos[Form("qgl_llPt_bin%d_flavor_b",bin)];
	TH1F* qgl_llPt_data=histos[Form("qgl_llPt_bin%d_data",bin)];

	TH1F* btag_llPt_q=histos[Form("btag_llPt_bin%d_flavor_q",bin)];
	TH1F* btag_llPt_g=histos[Form("btag_llPt_bin%d_flavor_g",bin)];
	TH1F* btag_llPt_c=histos[Form("btag_llPt_bin%d_flavor_c",bin)];
	TH1F* btag_llPt_b=histos[Form("btag_llPt_bin%d_flavor_b",bin)];
	TH1F* btag_llPt_data=histos[Form("btag_llPt_bin%d_data",bin)];

	//This is a check that all is working written in the out-file	
	if(WriteResults){	
		TCanvas *c_comp=new TCanvas(Form("c_comp_%d",bin),"c_comp");
			qgl_llPt_data->Draw("P");
			qgl_llPt_q->SetLineColor(kBlue);qgl_llPt_q->Draw("HIST SAME");
			qgl_llPt_g->SetLineColor(kRed) ;qgl_llPt_g->Draw("HIST SAME");
			qgl_llPt_q->SetLineColor(kYellow+2);qgl_llPt_c->Draw("HIST SAME");
			qgl_llPt_q->SetLineColor(kGreen+2);qgl_llPt_b->Draw("HIST SAME");
		c_comp->Write();
		TCanvas *c_comp2=new TCanvas(Form("c_comp2_%d",bin),"c_comp2");
			btag_llPt_data->Draw("P");
			btag_llPt_q->SetLineColor(kBlue);btag_llPt_q->Draw("HIST SAME");
			btag_llPt_g->SetLineColor(kRed) ;btag_llPt_g->Draw("HIST SAME");
			btag_llPt_q->SetLineColor(kYellow+2);btag_llPt_c->Draw("HIST SAME");
			btag_llPt_q->SetLineColor(kGreen+2);btag_llPt_b->Draw("HIST SAME");
		c_comp2->Write();
	} //end WriteResults

	if( (qgl_llPt_data->Integral()==0) || (qgl_llPt_q->Integral()==0) || (qgl_llPt_g->Integral()==0) || (qgl_llPt_c->Integral()==0) || (qgl_llPt_b->Integral()==0) )
		{
		(*q_frac)[bin]=0 ;
	        (*g_frac)[bin]=0 ;
	        (*c_frac)[bin]=0 ;
	        (*b_frac)[bin]=0 ;
		
		(*q_err)[bin]= 0 ;
	        (*g_err)[bin]= 0 ;
	        (*c_err)[bin]= 0 ;
	        (*b_err)[bin]= 0 ;
		continue;
		}

	RooDataHist qgl_data (Form("qgl_data_bin%d" ,bin),"qgl_data"   ,l,  qgl_llPt_data   ); 
		qgl_data.weightError(RooAbsData::SumW2) ;
	RooDataHist btag_data(Form("btag_data_bin%d",bin),"btag_data",b,  btag_llPt_data  ); 
		btag_data.weightError(RooAbsData::SumW2) ;

	//Import MC histo in RooFit	
	RooDataHist qgl_q(Form("qgl_DY_q_bin%d",bin)  ,"qgl_DY_q"   ,l,  qgl_llPt_q ); 
		qgl_q.weightError(RooAbsData::SumW2) ;
	RooDataHist btag_q(Form("btag_DY_q_bin%d",bin),"btag_DY_q"  ,b,  btag_llPt_q); 
		btag_q.weightError(RooAbsData::SumW2) ;

	RooDataHist qgl_g(Form("qgl_DY_g_bin%d",bin)  ,"qgl_DY_g"   ,l,  qgl_llPt_g ); 
		qgl_g.weightError(RooAbsData::SumW2) ;
	RooDataHist btag_g(Form("btag_DY_g_bin%d",bin),"btag_DY_g"  ,b,  btag_llPt_g); 
		btag_g.weightError(RooAbsData::SumW2) ;

	RooDataHist qgl_c(Form("qgl_DY_c_bin%d",bin)  ,"qgl_DY_c"   ,l,  qgl_llPt_c ); 
		qgl_c.weightError(RooAbsData::SumW2) ;
	RooDataHist btag_c(Form("btag_DY_c_bin%d",bin),"btag_DY_c"  ,b,  btag_llPt_c); 
		btag_c.weightError(RooAbsData::SumW2) ;

	RooDataHist qgl_b(Form("qgl_DY_b_bin%d",bin)  ,"qgl_DY_b"   ,l,  qgl_llPt_b ); 
		qgl_b.weightError(RooAbsData::SumW2) ;
	RooDataHist btag_b(Form("btag_DY_b_bin%d",bin),"btag_DY_b"  ,b,  btag_llPt_b); 
		btag_b.weightError(RooAbsData::SumW2) ;


	// Fractions
	RooRealVar c0("c0","quark fraction",.56,0.,1.) ;
	RooRealVar c1("c1","gluon fraction",0.4,0.,1.) ;
	RooRealVar c2("c2","c fraction"    ,0.02,0.,.1) ;
	RooRealVar c3("c3","b fraction"    ,0.02,0.,.1) ;

	//Type: 0=QGL 1=BTAG 2=...
	RooCategory type("type","type") ;
  		type.defineType("qgl")  ;
  		type.defineType("btag") ;

	//Construct combined dataset
	RooDataHist combData("combData","combined data", RooArgList(l,b) ,Index(type), Import("qgl",qgl_data),Import("btag",btag_data)) ;
	
	//Create Model QGL
	RooHistPdf qgl_q_histPdf(Form("qgl_q_histPdf_bin%d",bin),"qgl_q_histPdf",l,qgl_q,0) ;
	RooHistPdf qgl_g_histPdf(Form("qgl_g_histPdf_bin%d",bin),"qgl_g_histPdf",l,qgl_g,0) ;
	RooHistPdf qgl_c_histPdf(Form("qgl_c_histPdf_bin%d",bin),"qgl_c_histPdf",l,qgl_c,0) ;
	RooHistPdf qgl_b_histPdf(Form("qgl_b_histPdf_bin%d",bin),"qgl_b_histPdf",l,qgl_b,0) ;
	RooAddPdf modelQGL(Form("modelQGL_bin%d",bin),"modelQGL",RooArgList(qgl_q_histPdf,qgl_g_histPdf,qgl_c_histPdf,qgl_b_histPdf),RooArgList(c0,c1,c2)) ;

	//Create Model BTAG
	RooHistPdf btag_q_histPdf(Form("btag_q_histPdf_bin%d",bin),"btag_q_histPdf",b,btag_q,0) ;
	RooHistPdf btag_g_histPdf(Form("btag_g_histPdf_bin%d",bin),"btag_g_histPdf",b,btag_g,0) ;
	RooHistPdf btag_c_histPdf(Form("btag_c_histPdf_bin%d",bin),"btag_c_histPdf",b,btag_c,0) ;
	RooHistPdf btag_b_histPdf(Form("btag_b_histPdf_bin%d",bin),"btag_b_histPdf",b,btag_b,0) ;
	//not c3 so sum=1
	RooAddPdf modelBTAG(Form("modelBTAG_bin%d",bin),"modelBTAG",RooArgList(btag_q_histPdf,btag_g_histPdf,btag_c_histPdf,btag_b_histPdf),RooArgList(c0,c1,c2)) ;

	//Create Combined Model	
	RooSimultaneous simPdf(Form("simPdf_bin%d",bin),"simultaneous pdf",type) ;
		simPdf.addPdf(modelQGL,"qgl");
		simPdf.addPdf(modelBTAG,"btag");

	//FIT	
	RooFitResult* r=simPdf.fitTo(combData,Save(),SumW2Error(kTRUE));

	//Save Results
	RooPlot* frame1 =l.frame();
	RooPlot* frame2 =b.frame();

	qgl_data.plotOn (frame1,DataError(RooAbsData::SumW2));
	modelQGL.plotOn (frame1,DataError(RooAbsData::SumW2));
	btag_data.plotOn(frame2,DataError(RooAbsData::SumW2));
	modelBTAG.plotOn(frame2,DataError(RooAbsData::SumW2));
	
	//Save Fraction in Vectors
			//q_frac.push_back( r->floatParsFinal().at( 0 ) );
        		//g_frac.push_back( r->floatParsFinal().at( 1 ) );
        		//c_frac.push_back( r->floatParsFinal().at( 2 ) );
        		//b_frac.push_back( r->floatParsFinal().at( 3 ) );
			//
			//q_err.push_back( TMath::Sqrt(r->reducedCovarianceMatrix( RooArgList("c0","c0") ) ) );
        		//g_err.push_back( TMath::Sqrt(r->reducedCovarianceMatrix( RooArgList("c1","c1") ) ) );
        		//c_err.push_back( TMath::Sqrt(r->reducedCovarianceMatrix( RooArgList("c2","c2") ) ) );
        		//b_err.push_back( TMath::Sqrt(r->reducedCovarianceMatrix( RooArgList("c3","c3") ) ) );
	cout<<"NOT USING RooFitResults but direct variables"<<endl;
	
	(*q_frac)[bin]= c0.getVal() ;
        (*g_frac)[bin]= c1.getVal() ;
        (*c_frac)[bin]= c2.getVal() ;
        (*b_frac)[bin]= 1-c0.getVal()-c1.getVal()-c2.getVal() ;
	
	(*q_err)[bin]= c0.getError() ;
        (*g_err)[bin]= c1.getError() ;
        (*c_err)[bin]= c2.getError() ;
        (*b_err)[bin]= c2.getError() ; //TODO

	//Save Histo in OutFile
	if(WriteResults){
		TCanvas *c_frame=new TCanvas(Form("c_bin%d",bin),"c");
		c_frame->Divide(2);
		c_frame->cd(1); frame1->Draw();
		c_frame->cd(2); frame2->Draw();
		fout->cd();
		c_frame->Write();
	}
	
	} //end bin llPt
}

#ifdef STANDALONE
int main(int argc,char**argv)
{

string configFile;

if(argc<2) configFile="data/config.ini";
else configFile=argv[1];
string configFile2;
if(argc>=3){configFile2=argv[2];}
else configFile2="";

Read A(configFile.c_str());
//sscanf(A.ReadParFromMultFile(configFile2.c_str(),"CHID2"),"%d",&CHID2 );
string DirOut=A.ReadParFromMultFile(configFile2.c_str(),"OUTDIR"); DirOut+="/";  
float lumi; sscanf(A.ReadParFromMultFile(configFile2.c_str(),"LUMI"),"%f",&lumi);

fdata_4=TFile::Open(Form("%sQG_DoubleMu_4.root",DirOut.c_str()));
fdata_1=TFile::Open(Form("%sQG_DoubleE_1.root",DirOut.c_str()));
fdata_2=TFile::Open(Form("%sQG_MuEG_2.root",DirOut.c_str()));

CHID2=-4;
fDY_4=TFile::Open(Form("%sQG_DY_%d.root",DirOut.c_str(),-CHID2));
fTT_4=TFile::Open(Form("%sQG_TT_%d.root",DirOut.c_str(),-CHID2));
fWJ_4=TFile::Open(Form("%sQG_WJ_%d.root",DirOut.c_str(),-CHID2));
fWW_4=TFile::Open(Form("%sQG_WW_%d.root",DirOut.c_str(),-CHID2));
fWZ_4=TFile::Open(Form("%sQG_WZ_%d.root",DirOut.c_str(),-CHID2));
fZZ_4=TFile::Open(Form("%sQG_ZZ_%d.root",DirOut.c_str(),-CHID2));

CHID2=-2;
fDY_2=TFile::Open(Form("%sQG_DY_%d.root",DirOut.c_str(),-CHID2));
fTT_2=TFile::Open(Form("%sQG_TT_%d.root",DirOut.c_str(),-CHID2));
fWJ_2=TFile::Open(Form("%sQG_WJ_%d.root",DirOut.c_str(),-CHID2));
fWW_2=TFile::Open(Form("%sQG_WW_%d.root",DirOut.c_str(),-CHID2));
fWZ_2=TFile::Open(Form("%sQG_WZ_%d.root",DirOut.c_str(),-CHID2));
fZZ_2=TFile::Open(Form("%sQG_ZZ_%d.root",DirOut.c_str(),-CHID2));

CHID2=-1;
fDY_1=TFile::Open(Form("%sQG_DY_%d.root",DirOut.c_str(),-CHID2));
fTT_1=TFile::Open(Form("%sQG_TT_%d.root",DirOut.c_str(),-CHID2));
fWJ_1=TFile::Open(Form("%sQG_WJ_%d.root",DirOut.c_str(),-CHID2));
fWW_1=TFile::Open(Form("%sQG_WW_%d.root",DirOut.c_str(),-CHID2));
fWZ_1=TFile::Open(Form("%sQG_WZ_%d.root",DirOut.c_str(),-CHID2));
fZZ_1=TFile::Open(Form("%sQG_ZZ_%d.root",DirOut.c_str(),-CHID2));


QGAnalysis(lumi,Form("%sQG_Analysis.root",DirOut.c_str()) );

return 0;
}
#endif
