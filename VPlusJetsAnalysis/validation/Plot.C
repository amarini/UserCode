#include "TChain.h"
#include "TH1F.h"
#include <vector>
#include "TTree.h"
#include "TLorentzVector.h"
#include <iostream>

using namespace std;

int Loop(TChain *data, TH1F* h_data,int type);
int Plot(const char *dir="./")
{
TChain *mc=new TChain("accepted/events");
TChain *data=new TChain("accepted/events");
mc->Add( (string(dir)+"TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2.root").c_str());
data->Add(  (string(dir)+"MuEG*.root").c_str());

//data->Draw("llM>>h(100,60,120)","lepChId[0]*lepChId[1]==-2")
//mc->Draw("llM>>h2(100,60,120)","PUWeight*18.74*(lepChId[0]*lepChId[1]==-2)","HIST SAME")

TH1F *h_data=new TH1F("h_data","h_data",30,60,120);
TH1F *h_mc=new TH1F("h_mc","h_mc",30,60,120);
Loop(data,h_data,0);
Loop(mc,h_mc,1);

h_data->SetMarkerStyle(20);
h_mc->SetLineColor(kRed);

h_data->Draw("P");
h_mc->Draw("HIST SAME");
return 0;
}


int Loop(TChain *t, TH1F* h_t,int type)
{
float llM;
double PUWeight;
float lumi=18.74;
vector<int> *lepChId=NULL;
vector<float> *jetPt=NULL;
vector<float> *jetEta=NULL;
vector<float> *jetPhi=NULL;
vector<float> *jetE=NULL;
vector<float> *lepPt=NULL;
vector<float> *lepEta=NULL;
vector<float> *lepPhi=NULL;
vector<float> *lepE=NULL;

t->SetBranchAddress("llM",&llM);
t->SetBranchAddress("lepChId",&lepChId);
t->SetBranchAddress("lepPt",&lepPt);
t->SetBranchAddress("lepEta",&lepEta);
t->SetBranchAddress("lepPhi",&lepPhi);
t->SetBranchAddress("lepE",&lepE);
t->SetBranchAddress("jetPt",&jetPt);
t->SetBranchAddress("jetEta",&jetEta);
t->SetBranchAddress("jetPhi",&jetPhi);
t->SetBranchAddress("jetE",&jetE);

if(type!=0){
t->SetBranchAddress("PUWeight",&PUWeight);
}

for(long i=0;i<t->GetEntries() ;i++)
	{
	t->GetEntry(i);	
	if(lepPt->size()<2) continue;
	//	cout<<"Cut 1"<<endl;
	if( (*lepChId)[0]*(*lepChId)[1] != -2) continue;
	//	cout<<"Cut 2"<<endl;
	if( (*lepPt)[0] <20 ) continue;
	//	cout<<"Cut 3"<<endl;
	if( (*lepPt)[1] <20 ) continue;
	//	cout<<"Cut 4"<<endl;
	TLorentzVector l1,l2,j;
	l1.SetPtEtaPhiE( (*lepPt)[0],
		 	 (*lepEta)[0],
		 	 (*lepPhi)[0],
		 	 (*lepE)[0]
			);
	l2.SetPtEtaPhiE( (*lepPt)[1],
		 	 (*lepEta)[1],
		 	 (*lepPhi)[1],
		 	 (*lepE)[1]
			);
	bool foundJet=false;
	for(int k=0;k<jetPt->size();k++)
		{
		if( (*jetPt)[k] <50 ) continue;
			j.SetPtEtaPhiE( (*jetPt)[k],(*jetEta)[k],(*jetPhi)[k],(*jetE)[k]);
		float R1=l1.DeltaR(j);
		float R2=l2.DeltaR(j);
		float R=TMath::Min(R1,R2);
		if(R>0.4 && j.Pt()>50 )foundJet=true;

		}
	if(!foundJet)continue;
//		cout<<"Cut 5"<<endl;
	if(type==0)h_t->Fill(llM,1);
	//if(type!=0)h_t->Fill(llM,PUWeight*lumi);
	if(type!=0)h_t->Fill(double(llM),double(PUWeight*lumi));
//		cout<<"llM "<<llM<<"Weight= "<<lumi*PUWeight<<endl;

	}
}

