#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH1.h"
#include "TCanvas.h"
#include <stdlib.h>
#include <stdio.h>
#include "TH2F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TSystem.h"
#include "TRandom.h"
#include "RooUnfold/RooUnfoldResponse.h"
#ifdef __CINT__
gSystem->Load("libRooUnfold.so");
#endif

TH1F* MergeHistos(TH1F* a,TH1F*b,const char*name="merged")
{
if(a->GetNbinsX()!=b->GetNbinsX())return NULL;
float aMin=a->GetBinLowEdge(1);
float aMax=a->GetBinLowEdge(a->GetNbinsX()+1);
float bMin=b->GetBinLowEdge(1);
float bMax=b->GetBinLowEdge(b->GetNbinsX()+1);
//fprintf(stderr,"%d: %.5f %.5f\n",a->GetNbinsX(),aMin,aMax);
TH1F*R=new TH1F(name,name,(a->GetNbinsX())*2,aMin,aMax+(aMax-aMin));
R->Sumw2();
for(int i=1;i<=a->GetNbinsX();i++)
	{
	R->SetBinContent(i,a->GetBinContent(i));	
	R->SetBinError(i,a->GetBinError(i));	
	}
for(int i=1;i<=b->GetNbinsX();i++)
	{
	R->SetBinContent(i+a->GetNbinsX(),b->GetBinContent(i));	
	R->SetBinError(i+a->GetNbinsX(),b->GetBinError(i));	
	}
return R;
}

int DivideHisto(TH1F* a,TH1F**b, TH1F**c)
{
if(  (a->GetNbinsX()&1)!=0 ) return 1; //not event number of bins
char name[1023];
sprintf(name,"%s_1",a->GetName());
*b=new TH1F(name,name,a->GetNbinsX()/2,a->GetBinLowEdge(1),a->GetBinLowEdge(a->GetNbinsX()+1)/2);
sprintf(name,"%s_2",a->GetName());
*c=new TH1F(name,name,a->GetNbinsX()/2,a->GetBinLowEdge(1),a->GetBinLowEdge(a->GetNbinsX()+1)/2);
//TODO	
for(int i=1;i<=a->GetNbinsX()/2;i++){(*b)->SetBinContent(i,a->GetBinContent(i));(*b)->SetBinError(i,a->GetBinError(i));}
for(int i=1;i<=a->GetNbinsX()/2;i++){(*c)->SetBinContent(i,a->GetBinContent(i+a->GetNbinsX()/2));(*c)->SetBinError(i,a->GetBinError(i+a->GetNbinsX()/2));}
return 0;
}


int Unfold(const char*fileName1,
	  const char*fileName2,
  	  const char *varName="ptD",
	  const char *range="(100,0,1)",
	  float PtMin=100,float PtMax=110,
	  float RhoMin=12,float RhoMax=13,
	  float alpha1=0.85,float alpha2=0.55,
	  float isMC=true,
	  const char*MCFile="",
	  const char*jet1="Jet0",
	  const char*jet2="Jet0",
	  const char*treeName="tree_passedEvents")
{

 //Some stuff
        gROOT->SetStyle("Plain");
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
//
TFile *f1=TFile::Open(fileName1);
TFile *f2=TFile::Open(fileName2);
if( (f1==NULL) || (f2==NULL)){fprintf(stderr,"FILES DOES NOT EXIST!\n");return 1;}

//Create the Unfold Matrix - merge of the histograms
//use alpha1 alpha2 
	int nBins=100;
	float xMin,xMax;
	sscanf(range,"(%d,%f,%f)",&nBins,&xMin,&xMax);
	TH1F *aux=new TH1F("aux","aux",nBins*2,xMin,xMax+(xMax-xMin));
	RooUnfoldResponse *R=new RooUnfoldResponse(aux,aux,"Unfold");
	for(int i=1;i<=nBins;i++)R->Fill(aux->GetBinCenter(i)      ,aux->GetBinCenter(i)     ,alpha1); //RECO - GEN (PhDi - QG)
	for(int i=1;i<=nBins;i++)R->Fill(aux->GetBinCenter(nBins+i),aux->GetBinCenter(nBins+i),1-alpha2);
	for(int i=1;i<=nBins;i++)R->Fill(aux->GetBinCenter(nBins+i),aux->GetBinCenter(i)      ,alpha2);//TODO: tocheck
	for(int i=1;i<=nBins;i++)R->Fill(aux->GetBinCenter(i)      ,aux->GetBinCenter(nBins+i),1-alpha1);
	
	//I do not want it normalized to 1 in row of Q/G
	for(int i=1;i<=nBins;i++)R->Miss(aux->GetBinCenter(i),1-alpha1-alpha2);
	for(int i=1;i<=nBins;i++)R->Miss(aux->GetBinCenter(nBins+i),alpha1+alpha2-1);

	TCanvas *c1=new TCanvas();
	TH2F*Matrix=((TH2F*)R->Hresponse());
	Matrix->GetXaxis()->SetTitle("RECO");
	Matrix->GetYaxis()->SetTitle("Q/G");
	Matrix->Draw("BOX");
//Getting the histograms
	TTree *t1=(TTree*)f1->Get(treeName);
	TTree *t2=(TTree*)f2->Get(treeName);
	char name[1023];
	char selection[1023];
	TH1F *h1=new TH1F("var1","var1",nBins,xMin,xMax);h1->Sumw2();
	TH1F *h2=new TH1F("var2","var2",nBins,xMin,xMax);h2->Sumw2();
	sprintf(name,"%s%s>>var1",varName,jet1);
	if(isMC)sprintf(selection,"eventWeight*( %f < pt%s && pt%s<%f && %f<rhoPF && rhoPF<%f)",PtMin,jet1,jet1,PtMax,RhoMin,RhoMax);
	else sprintf(selection,"( %f < pt%s && pt%s<%f && %f<rhoPF && rhoPF<%f)",PtMin,jet1,jet1,PtMax,RhoMin,RhoMax);
	t1->Draw(name,selection,"goff");h1->Scale(1./h1->Integral());
	fprintf(stderr,"1:===%s===%s===\n",name,selection);

	sprintf(name,"%s%s>>var2",varName,jet2);
	if(isMC)sprintf(selection,"eventWeight*( %f < pt%s && pt%s<%f && %f<rhoPF && rhoPF <%f&& passedID_FULL && !btagged) ",PtMin,jet2,jet2,PtMax,RhoMin,RhoMax);
	else sprintf(selection,"( %f < pt%s && pt%s<%f && %f<rhoPF && rhoPF <%f&& passedID_FULL && !btagged) ",PtMin,jet2,jet2,PtMax,RhoMin,RhoMax);
	t2->Draw(name,selection,"goff");h2->Scale(1./h2->Integral());
	fprintf(stderr,"2:===%s===%s===\n",name,selection);
//Unfold the distributions
	TH1F* H=MergeHistos(h1,h2,"merged");
	//regularisation
	for(int i=1;i<H->GetNbinsX();i++)if((H->GetBinContent(i)==0)&&(H->GetBinContent(i-1)==0) && (H->GetBinContent(i+1)==0))H->SetBinError(i,0.0001);
	RooUnfoldSvd *Unfold=new RooUnfoldSvd(R,H,int(nBins));
	//RooUnfoldInvert *Unfold=new RooUnfoldInvert(R,H);
//Print the output
	TFile *mc;
	if(MCFile[0]!='\0')mc=TFile::Open(MCFile);
	TTree *t_mc;
	if(MCFile[0]!='\0')t_mc=(TTree*)mc->Get(treeName);
	TH1F *hq=new TH1F("quark","quark",nBins,xMin,xMax);hq->Sumw2();
	sprintf(name,"%s%s>>quark",varName,jet1);
	sprintf(selection,"eventWeight*( %f < pt%s && pt%s<%f && %f<rhoPF && rhoPF<%f && abs(pdgIdPart%s) <5) ",PtMin,jet1,jet1,PtMax,RhoMin,RhoMax,jet1);

	if(isMC){t1->Draw(name,selection,"goff");hq->Scale(1./hq->Integral());}
	else if(MCFile[0]!='\0'){t_mc->Draw(name,selection,"goff");hq->Scale(1./hq->Integral());}
	TH1F *hg=new TH1F("gluon","gluon",nBins,xMin,xMax);hg->Sumw2();

	sprintf(name,"%s%s>>gluon",varName,jet1);
	sprintf(selection,"eventWeight*( %f < pt%s && pt%s<%f && %f<rhoPF&&rhoPF<%f && abs(pdgIdPart%s) ==21) ",PtMin,jet1,jet1,PtMax,RhoMin,RhoMax,jet1);

fprintf(stderr,"=%s==%s=\n",name,selection);//DEBUG

	if(isMC){t1->Draw(name,selection,"goff");hg->Scale(1./hg->Integral());}
	else if(MCFile[0]!='\0'){t_mc->Draw(name,selection,"goff");hg->Scale(1./hg->Integral());}
	

	TCanvas *c=new TCanvas();
	TH1F* unf=Unfold->Hreco(3);
	unf->SetLineColor(kBlue+2);

	TH1F*allTrue=MergeHistos(hg,hq,"true-g-q");
	allTrue->SetLineColor(kRed+2);

	fprintf(stderr,"DEBUG:%.3f == 2.0?\n",unf->Integral());

	allTrue->Draw("");
	unf->Draw("HIST SAME"); //unfold
	H->SetMarkerColor(kMagenta+2);
	H->SetMarkerStyle(20);
	H->Draw("P SAME");
	
	TH1F *Inv=MergeHistos(h2,h1,"merged2");Inv->Draw("HIST SAME");

	char LegendName[1023];
	sprintf(LegendName,"%.0f<P_{T}<%.0f & %.0f<#rho<%.0f ",PtMin,PtMax,RhoMin,RhoMax);
	TLegend *L=new TLegend(0.45,0.75,0.55,.89,LegendName);
	L->AddEntry("merged","merged");
	L->AddEntry("true-g-q","g-q");
	L->AddEntry("Unfold","Unfolded");
	L->Draw();
	
//DEBUG
	TCanvas *c2=new TCanvas();
	hq->SetLineWidth(2);
	hq->SetLineColor(kBlue+2);
	hq->SetFillColor(kBlue-3);
	hq->SetFillStyle(3004);
	hq->Draw("HIST");
	hg->SetLineWidth(2);
	hg->SetLineColor(kRed+2);
	hg->SetFillColor(kRed-3);
	hg->SetFillStyle(3005);
	hg->Draw("HIST SAME");

	h1->SetMarkerColor(kMagenta+2);
	h1->SetMarkerStyle(20);
	h1->Draw("P SAME");
	h2->SetMarkerColor(kGreen+2);
	h2->SetMarkerStyle(29);
	h2->Draw("P SAME");
	
	TH1F *unf_q=0,*unf_g=0;
	DivideHisto(unf,&unf_q,&unf_g); 
		fprintf(stderr,"\033[01;31m Scaling\033[00m unfold histos: %.3f - %.3f\n",1./unf_q->Integral(),1./unf_g->Integral()); 
		unf_q->Scale(1./unf_q->Integral()); unf_g->Scale(1./unf_g->Integral());
	unf_q->SetMarkerStyle(29);
	unf_g->SetMarkerStyle(30);
	unf_q->SetMarkerColor(kBlue);
	unf_g->SetMarkerColor(kRed);
	unf_q->Draw("P SAME");
	unf_g->Draw("P SAME");
	
	sprintf(LegendName,"%.0f<P_{T}<%.0f & %.0f<#rho<%.0f ",PtMin,PtMax,RhoMin,RhoMax);
	TLegend *L2=new TLegend(0.7,.7,.89,.89,LegendName);
	L2->AddEntry("var1","dijet","FP");
	L2->AddEntry("var2","phjet","FP");
	L2->AddEntry("quark","quark","F");
	L2->AddEntry("gluon","gluon","F");
	L2->Draw();

//DEBUG	

//Draw log|di|
//TCanvas *c3=new TCanvas();
//c3->SetLogy();
//TMatrixD *U,*V,*A;
//TVectorD *b=new TVectorD(),*SS,*d;
//TDecompSVD svd;
//A=(TMatrixD*)R->Mresponse().Clone();
//svd.SetMatrix(*A);
//U=(TMatrixD*)svd.GetU().Clone();
//V=(TMatrixD*)svd.GetV().Clone();
//SS=(TVectorD*)svd.GetSig().Clone();
//b->ResizeTo(H->GetNbinsX());
//for(int i=1;i<=H->GetNbinsX();i++)(*b)[i-1]=H->GetBinContent(i);
//U->Transpose(*U);
//d=(TVectorD*)b->Clone();
//(*d)*=(*U);
//U->Transpose(*U);//transpose back
//for(int i=0;i<H->GetNbinsX();i++)printf("%.3lf ",(*d)[i]);printf("\n");
//for(int i=0;i<H->GetNbinsX();i++)(*d)[i]=fabs((*d)[i]);
//TH1F*hd=new TH1F(*d);
//hd->Draw("HIST");
//
////Draw Eingen functions
//TCanvas *c5=new TCanvas("eigen","eigen");
//
//c5->Divide(6,5);
//for(int i=0;i<30;i++)
//        {
//        c5->cd(i+1);
//        TVectorD z;z.ResizeTo(H->GetNbinsX());
//        z[i]=1.0;
//        z*=(*U);
//        sprintf(name,"eigen%d",i+1);
//        TH1F *h=new TH1F(z);
//        h->SetName(name);
//        h->Draw("HIST");
//        }

return 0;
}

