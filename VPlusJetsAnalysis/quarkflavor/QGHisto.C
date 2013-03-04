#include "TH1F.h"
#include "TFile.h"
#include "THStack.h"
#include <string>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"

using namespace std;
bool LogScale=true;

TFile *fIn;


int Plot(string plotDir ){

//Get Histograms
TH1F*	llPt_data	= (TH1F*)fIn->Get("llPt_data");
TH1F*	llPt_q		= (TH1F*)fIn->Get("llPt_q");
TH1F*	llPt_g		= (TH1F*)fIn->Get("llPt_g");
TH1F*	llPt_c		= (TH1F*)fIn->Get("llPt_c");
TH1F*	llPt_b		= (TH1F*)fIn->Get("llPt_b");
TH1F*	r_q		= (TH1F*)fIn->Get("r_q");
TH1F*	r_g		= (TH1F*)fIn->Get("r_g");
TH1F*	r_c		= (TH1F*)fIn->Get("r_c");
TH1F*	r_b		= (TH1F*)fIn->Get("r_b");
TH1F*	llPt_DY_u	= (TH1F*)fIn->Get("llPt_DY_u");
TH1F*	llPt_DY_q	= (TH1F*)fIn->Get("llPt_DY_q");
TH1F*	llPt_DY_g	= (TH1F*)fIn->Get("llPt_DY_g");
TH1F*	llPt_DY_c	= (TH1F*)fIn->Get("llPt_DY_c");
TH1F*	llPt_DY_b	= (TH1F*)fIn->Get("llPt_DY_b");

TH1F*	llPt_DY		= (TH1F*) llPt_DY_q->Clone("llPt_DY");
		llPt_DY->Add(llPt_DY_g);
		llPt_DY->Add(llPt_DY_c);
		llPt_DY->Add(llPt_DY_b);
		llPt_DY->Add(llPt_DY_u);
	
TH1F* r_DY_q		= (TH1F*) llPt_DY_q->Clone("r_DY_q"); r_DY_q->Divide(llPt_DY);
TH1F* r_DY_g		= (TH1F*) llPt_DY_g->Clone("r_DY_g"); r_DY_g->Divide(llPt_DY);
TH1F* r_DY_c		= (TH1F*) llPt_DY_c->Clone("r_DY_c"); r_DY_c->Divide(llPt_DY);
TH1F* r_DY_b		= (TH1F*) llPt_DY_b->Clone("r_DY_b"); r_DY_b->Divide(llPt_DY);
TH1F* r_DY_u		= (TH1F*) llPt_DY_u->Clone("r_DY_u"); r_DY_u->Divide(llPt_DY);
	
	//Set Style
	llPt_data->SetMarkerStyle(20);llPt_data->SetMarkerColor(kBlack   );llPt_data->SetLineColor(kBlack   );llPt_data->SetFillColor(kBlack   ); llPt_data->SetMarkerSize(1  );llPt_data->SetFillStyle(3004);
	llPt_q   ->SetMarkerStyle(20);llPt_q   ->SetMarkerColor(kBlue+2  );llPt_q   ->SetLineColor(kBlue+2  );llPt_q   ->SetFillColor(kBlue+2  ); llPt_q   ->SetMarkerSize(0.8);llPt_q   ->SetFillStyle(3004);
	llPt_g   ->SetMarkerStyle(20);llPt_g   ->SetMarkerColor(kRed+2   );llPt_g   ->SetLineColor(kRed+2   );llPt_g   ->SetFillColor(kRed+2   ); llPt_g   ->SetMarkerSize(0.8);llPt_g   ->SetFillStyle(3005);
	llPt_c   ->SetMarkerStyle(20);llPt_c   ->SetMarkerColor(kYellow+2);llPt_c   ->SetLineColor(kYellow+2);llPt_c   ->SetFillColor(kYellow+2); llPt_c   ->SetMarkerSize(0.8);llPt_c   ->SetFillStyle(3004);
	llPt_b   ->SetMarkerStyle(20);llPt_b   ->SetMarkerColor(kGreen+2 );llPt_b   ->SetLineColor(kGreen+2 );llPt_b   ->SetFillColor(kGreen+2 ); llPt_b   ->SetMarkerSize(0.8);llPt_b   ->SetFillStyle(3005);
	
	llPt_DY_q   ->SetMarkerStyle(24);llPt_DY_q   ->SetMarkerColor(kBlue+4  );llPt_DY_q   ->SetLineColor(kBlue+4  ); llPt_DY_q   ->SetFillColor(kBlue+4  );llPt_DY_q   ->SetMarkerSize(0.8);llPt_DY_q   ->SetFillStyle(3003);
	llPt_DY_g   ->SetMarkerStyle(24);llPt_DY_g   ->SetMarkerColor(kRed+4   );llPt_DY_g   ->SetLineColor(kRed+4   ); llPt_DY_g   ->SetFillColor(kRed+4   );llPt_DY_g   ->SetMarkerSize(0.8);llPt_DY_g   ->SetFillStyle(3003);
	llPt_DY_c   ->SetMarkerStyle(24);llPt_DY_c   ->SetMarkerColor(kOrange  );llPt_DY_c   ->SetLineColor(kOrange  ); llPt_DY_c   ->SetFillColor(kOrange  );llPt_DY_c   ->SetMarkerSize(0.8);llPt_DY_c   ->SetFillStyle(3003);
	llPt_DY_b   ->SetMarkerStyle(24);llPt_DY_b   ->SetMarkerColor(kGreen+4 );llPt_DY_b   ->SetLineColor(kGreen+4 ); llPt_DY_b   ->SetFillColor(kGreen+4 );llPt_DY_b   ->SetMarkerSize(0.8);llPt_DY_b   ->SetFillStyle(3003);
	llPt_DY_u   ->SetMarkerStyle(24);llPt_DY_u   ->SetMarkerColor(kGray+2  );llPt_DY_u   ->SetLineColor(kGray+2  ); llPt_DY_u   ->SetFillColor(kGray+2  );llPt_DY_u   ->SetMarkerSize(0.8);llPt_DY_u   ->SetFillStyle(3003);
	llPt_DY     ->SetMarkerStyle(24);llPt_DY     ->SetMarkerColor(kBlack   );llPt_DY     ->SetLineColor(kBlack   ); llPt_DY    ->SetFillColor(kBlack    );llPt_DY     ->SetMarkerSize(0.8);llPt_DY     ->SetFillStyle(3003);
	
	
	r_q   ->SetMarkerStyle(20);r_q   ->SetMarkerColor(kBlue+2  );r_q   ->SetLineColor(kBlue+2  ); r_q   ->SetFillColor(kBlue+2  );r_q   ->SetMarkerSize(0.8);r_q   ->SetFillStyle(3004);
	r_g   ->SetMarkerStyle(20);r_g   ->SetMarkerColor(kRed+2   );r_g   ->SetLineColor(kRed+2   ); r_g   ->SetFillColor(kRed+2   );r_g   ->SetMarkerSize(0.8);r_g   ->SetFillStyle(3005);
	r_c   ->SetMarkerStyle(20);r_c   ->SetMarkerColor(kYellow+2);r_c   ->SetLineColor(kYellow+2); r_c   ->SetFillColor(kYellow+2);r_c   ->SetMarkerSize(0.8);r_c   ->SetFillStyle(3004);
	r_b   ->SetMarkerStyle(20);r_b   ->SetMarkerColor(kGreen+2 );r_b   ->SetLineColor(kGreen+2 ); r_b   ->SetFillColor(kGreen+2 );r_b   ->SetMarkerSize(0.8);r_b   ->SetFillStyle(3005);
	
	r_DY_q   ->SetMarkerStyle(24);r_DY_q   ->SetMarkerColor(kBlue+4  );r_DY_q   ->SetLineColor(kBlue+4  ); r_DY_q   ->SetFillColor(kBlue+4  );r_DY_q   ->SetMarkerSize(0.8);r_DY_q   ->SetFillStyle(3003);
	r_DY_g   ->SetMarkerStyle(24);r_DY_g   ->SetMarkerColor(kRed+4   );r_DY_g   ->SetLineColor(kRed+4   ); r_DY_g   ->SetFillColor(kRed+4   );r_DY_g   ->SetMarkerSize(0.8);r_DY_g   ->SetFillStyle(3003);
	r_DY_c   ->SetMarkerStyle(24);r_DY_c   ->SetMarkerColor(kOrange  );r_DY_c   ->SetLineColor(kOrange  ); r_DY_c   ->SetFillColor(kOrange  );r_DY_c   ->SetMarkerSize(0.8);r_DY_c   ->SetFillStyle(3003);
	r_DY_b   ->SetMarkerStyle(24);r_DY_b   ->SetMarkerColor(kGreen+4 );r_DY_b   ->SetLineColor(kGreen+4 ); r_DY_b   ->SetFillColor(kGreen+4 );r_DY_b   ->SetMarkerSize(0.8);r_DY_b   ->SetFillStyle(3003);
	r_DY_u   ->SetMarkerStyle(24);r_DY_u   ->SetMarkerColor(kGray+2  );r_DY_u   ->SetLineColor(kGray+2  ); r_DY_u   ->SetFillColor(kGray+2  );r_DY_u   ->SetMarkerSize(0.8);r_DY_u   ->SetFillStyle(3003);

	//Add to a Stack
	THStack *llPtData=new THStack("data","data");
		llPtData->Add(llPt_b);
		llPtData->Add(llPt_c);
		llPtData->Add(llPt_g);
		llPtData->Add(llPt_q);
	THStack *llPtMC=new THStack("mc","mc");
		llPtMC->Add(llPt_DY_b);
		llPtMC->Add(llPt_DY_c);
		llPtMC->Add(llPt_DY_g);
		llPtMC->Add(llPt_DY_q);
		llPtMC->Add(llPt_DY_u);

TCanvas *c=new TCanvas("c1","c1",800,1000);

	TPad*    upperPad = new TPad("upperPad", "upperPad",.005, .25, .995, .995);
        TPad*    lowerPad = new TPad("lowerPad", "lowerPad",.005, .005, .995, .2475);
        upperPad->Draw();
        lowerPad->Draw();
        upperPad->cd();
        if(LogScale)upperPad->SetLogy();
	llPtData->Draw("P E4");
	llPtMC->Draw("P E4 SAME");

	//printf("going to draw legend\n");	
	TLegend *L=new TLegend(0.75,0.45,0.89,.89,"");
		L->SetBorderSize(0);
		L->SetFillColor(0);
		L->SetFillStyle(0);
		L->AddEntry(llPt_data,"Data","PF");
		L->AddEntry(llPt_q   ,"Data (q)","PF");
		L->AddEntry(llPt_g   ,"Data (g)","PF");
		L->AddEntry(llPt_c   ,"Data (c)","PF");
		L->AddEntry(llPt_b   ,"Data (b)","PF");
		L->AddEntry(llPt_DY  ,"DY","PF");
		L->AddEntry(llPt_DY_u,"DY (u)","PF");
		L->AddEntry(llPt_DY_q,"DY (q)","PF");
		L->AddEntry(llPt_DY_g,"DY (g)","PF");
		L->AddEntry(llPt_DY_c,"DY (c)","PF");
		L->AddEntry(llPt_DY_b,"DY (b)","PF");
	L->Draw();
	//printf("going to draw latex\n");	
	TLatex *latex=new TLatex();
		latex->SetNDC();
		latex->SetTextFont(63);
		latex->SetTextSize(28);
		latex->SetTextAlign(21);
		latex->DrawLatex(0.5,.91,"CMS, Work in Progress, #sqrt{s}=8 TeV,L=18.7fb^{-1}");
		
	lowerPad->cd();
		r_q->GetXaxis()->SetTitle("P_{T}^{ll}");
		r_q->GetYaxis()->SetTitle("Fractions to all");
		r_q->GetYaxis()->SetRangeUser(0,1.15);
		
		r_q->Draw("P E4 ");
		r_g->Draw("P E4 SAME");
		r_b->Draw("P E4 SAME");
		r_c->Draw("P E4 SAME");

		r_DY_u->Draw("P E4 SAME");
		r_DY_q->Draw("P E4 SAME");
		r_DY_g->Draw("P E4 SAME");
		r_DY_b->Draw("P E4 SAME");
		r_DY_c->Draw("P E4 SAME");
	lowerPad->SetGridy(2);
	

c->SaveAs(Form("%sQG_Analysis.pdf",plotDir.c_str()));
}


int OpenFiles(string outDir, string plotDir){
        gROOT->SetStyle("Plain");
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
fIn=TFile::Open(Form("%sQG_Analysis.root",outDir.c_str()));
Plot(plotDir);
}

#ifdef STANDALONE
#include "src/ReadParameters.C"
int main(int argc,char **argv)
{
//Read A("data/config.ini");
string configFile;

if(argc<2) configFile="data/config.ini";
else configFile=argv[1];

Read A(configFile.c_str());

string DirOut=A.ReadParameter("OUTDIR"); DirOut+="/";
string PlotDir=A.ReadParameter("PLOTDIR"); PlotDir+="/";

OpenFiles(DirOut,PlotDir);
}
#endif 
