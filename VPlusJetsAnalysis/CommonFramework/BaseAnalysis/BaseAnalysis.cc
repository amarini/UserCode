
#include "BaseAnalysis.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include <cstdio>
#include <stdlib.h>
#include <iostream>

int BaseAnalysis::SetConfigFile(string a,string b)
{
configFileBase=a;
configFileNew=b;
}

int BaseAnalysis::OpenFiles()
{
t=new TChain("accepted/events");
for(int i=0; i < int(bookedFiles.size()) ; i++ )
	{
	t->Add( ReadMult(configFileBase.c_str(),configFileNew.c_str(),bookedFiles[i].c_str()) );	
	}
	return 0;
}

int BaseAnalysis::BookFile(string a)
{
bookedFiles.push_back(a);
	return 0;
}

int BaseAnalysis::CreateHistos()
{return 0;}

int BaseAnalysis::Smear(){
	if(SmearType==0) return 0;
}

int BaseAnalysis::Loop(int type,int j ,int NJobs)
{
//all variables 
t->SetBranchAddress("nVtx",&nVtx);
t->SetBranchAddress("nLeptons",&nLeptons);
t->SetBranchAddress("nPhotons",&nPhotons);
t->SetBranchAddress("nJets",&nJets);
t->SetBranchAddress("rho",&rho);
t->SetBranchAddress("rho25",&rho25);
t->SetBranchAddress("llM",&llM);
t->SetBranchAddress("llPt",&llPt);
t->SetBranchAddress("llPhi",&llPhi);
t->SetBranchAddress("llY",&llY);
t->SetBranchAddress("llDPhi",&llDPhi);
t->SetBranchAddress("llEta",&llEta);
t->SetBranchAddress("lepPt",&lepPt);
t->SetBranchAddress("lepEta",&lepEta);
t->SetBranchAddress("lepPhi",&lepPhi);
t->SetBranchAddress("lepE",&lepE);
t->SetBranchAddress("lepChId",&lepChId);
t->SetBranchAddress("jetPt",&jetPt);
t->SetBranchAddress("jetEta",&jetEta);
t->SetBranchAddress("jetPhi",&jetPhi);
t->SetBranchAddress("jetE",&jetE);
t->SetBranchAddress("jetUNC",&jetUNC);
t->SetBranchAddress("jetQGL",&jetQGL);
t->SetBranchAddress("jetBtag",&jetBtag);
t->SetBranchAddress("jetBeta",&jetBeta);
jetVeto=new vector<int>;
t->SetBranchAddress("photonPt",&photonPt);
t->SetBranchAddress("photonEta",&photonEta);
t->SetBranchAddress("photonPhi",&photonPhi);
t->SetBranchAddress("photonE",&photonE);
if(type>0)t->SetBranchAddress("PUWeight",&PUWeight);

if(type>0)t->SetBranchAddress("llMGEN",&llMGEN);
if(type>0)t->SetBranchAddress("llPtGEN",&llPtGEN);
if(type>0)t->SetBranchAddress("lepPtGEN",&lepPtGEN);
if(type>0)t->SetBranchAddress("lepEtaGEN",&lepEtaGEN);
if(type>0)t->SetBranchAddress("lepPhiGEN",&lepPhiGEN);
if(type>0)t->SetBranchAddress("lepEGEN",&lepEGEN);
if(type>0)t->SetBranchAddress("lepChIdGEN",&lepChIdGEN);
if(type>0)t->SetBranchAddress("jetPtGEN",&jetPtGEN);
if(type>0)t->SetBranchAddress("jetEtaGEN",&jetEtaGEN);
if(type>0)t->SetBranchAddress("jetPhiGEN",&jetPhiGEN);
if(type>0)t->SetBranchAddress("jetEGEN",&jetEGEN);

//t->SetBranchStatus("*",0); 
//ACTIVATE BRANCHES

if(debug>1)cout<<"Beginning loop "<<endl;
unsigned long long EntryMin=0;
unsigned long long EntryMax=t->GetEntries();

if((j>=0)  && (NJobs>0) && (j<NJobs)) //split in jobs
	{
	unsigned long long EntriesPerJob= EntryMax/NJobs;
	unsigned long long r= EntryMax % NJobs;
	if(r>0) EntriesPerJob++;
	
	EntryMin= j*EntriesPerJob;
	if(EntryMax>(j+1)*EntriesPerJob)
		EntryMax=(j+1)*EntriesPerJob;
	
	}
for(unsigned long long iEntry=EntryMin;iEntry<EntryMax ;iEntry++)
        {
        t->GetEntry(iEntry);
        if(debug>1)cout<<"Getting Entry "<<iEntry<<endl;
	Smear();
	
	FillHistosGEN();

	if(!BaseSelection())continue;

        //Triggers??

        if(debug>1)cout<<"Passed Selection"<<endl;
                //redo vetos
        jetVeto->clear();
        jetVeto->resize(jetPt->size());
        TLorentzVector j,l1,l2;
	if(lepPt->size()>=2){
        	l1.SetPtEtaPhiE( (*lepPt)[0],(*lepEta)[0],(*lepPhi)[0],(*lepE)[0]);
        	l2.SetPtEtaPhiE( (*lepPt)[1],(*lepEta)[1],(*lepPhi)[1],(*lepE)[1]);
	}
        for(int k=0;k<jetPt->size();k++)
                {
                if( (*jetPt)[k] <JetPtCut) continue;
                        j.SetPtEtaPhiE( (*jetPt)[k],(*jetEta)[k],(*jetPhi)[k],(*jetE)[k]);
                float R1=l1.DeltaR(j);
                float R2=l2.DeltaR(j);
                (*jetVeto)[k]=0;
                if(R1<JetDRCut) (*jetVeto)[k]=1;
                if(R2<JetDRCut) (*jetVeto)[k]|=2;

                }

        //int nJetsVeto=0;
        for(int iJ=0;iJ<jetVeto->size() && (*jetPt)[iJ]>=JetPtCut;iJ++)if( !( ( (*jetVeto)[iJ]&JetSel)  ) ) nJetsVeto++;
        if(debug>1)cout<<"nJetsVeto "<<nJetsVeto<<endl;
        if(debug>1)cout<<"jetSize Pt "<<jetPt->size()<<" Veto "<<jetVeto->size()<<endl;
        int jet0=-1,jet1=-1,jet2=-1;
                if(nJetsVeto>0)for(int iJ=0;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&JetSel)  ) ){jet0=iJ; break;}
                if(nJetsVeto>1)for(int iJ=jet0+1;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&JetSel)  ) ){jet1=iJ; break;}
                if(nJetsVeto>2)for(int iJ=jet1+1;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&JetSel)  ) ){jet2=iJ; break;}
        if( (jet0>=0) && (*jetPt)[jet0]< JetPtCut) jet0=-1;
        if( (jet1>=0) && (*jetPt)[jet1]< JetPtCut) jet1=-1;
        if( (jet2>=0) && (*jetPt)[jet2]< JetPtCut) jet2=-1;

	jet[0]=jet0;
	jet[1]=jet1;
	jet[2]=jet2;

        //double weight=1;
	weight=1;
        if(type>0)weight=PUWeight;
      if( jet0 < 0 )continue; //cut on the first jet

        //if((photonPt->size()>0)&&((*photonPt)[0]>130 )){
        //        //find the first isolated photon
        //        int pho0=0;
        //        //if(type==0){for(pho0=0;pho0<photonPt->size();pho0++) if( ((*photonIsoFPRNeutral)[pho0]+(*photonIsoFPRCharged)[pho0]+(*photonIsoFPRPhoton)[pho0])/(*photonPt)[pho0] < .5) {break;} }
        //        TLorentzVector g;
        //        if( (pho0<photonPt->size()) && ((*photonPt)[pho0]>130)){
        //        g.SetPtEtaPhiE( (*photonPt)[pho0],(*photonEta)[pho0],(*photonPhi)[pho0],(*photonE)[pho0]);
        //        histos["llgM"]->Fill((l1+l2+g).M(),weight);
        //        histos["l1gM"]->Fill((l1+g).M(),weight);
        //        histos["l2gM"]->Fill((l2+g).M(),weight);
        //        }
        //}
	
	

        //----------------------Jets with BStar: requirements already applied  E jet0 + veto
        for(int k=0;k<jetPt->size();k++)
                {
                if( (*jetPt)[k] <JetPtCut) continue;
                if( ( 1. - (*jetBeta)[k] >= 0.2 * TMath::Log( nVtx - 0.67))) //betaStar=1-beta
                        (*jetVeto)[k]|=8; //1= lept1, 2=lept2, 4= gamma, 8=betaStar
                }

        int jet0_BS=-1,jet1_BS=-1,jet2_BS=-1;
                if(nJetsVeto>0)for(int iJ=0;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&JetSel) || ((*jetVeto)[iJ]&8)) ) {jet0_BS=iJ; break;}
                if(nJetsVeto>1)for(int iJ=jet0_BS+1;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&JetSel) || ((*jetVeto)[iJ]&8) ) ){jet1_BS=iJ; break;}
                if(nJetsVeto>2)for(int iJ=jet1_BS+1;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&JetSel)  || ((*jetVeto)[iJ]&8)) ){jet2_BS=iJ; break;}
       if( (jet0_BS>=0) && (*jetPt)[jet0_BS]< JetPtCut) jet0_BS=-1;
        if( (jet1_BS>=0) && (*jetPt)[jet1_BS]< JetPtCut) jet1_BS=-1;
        if( (jet2_BS>=0) && (*jetPt)[jet2_BS]< JetPtCut) jet2_BS=-1;

	jet_BS[0]=jet0_BS;
	jet_BS[1]=jet1_BS;
	jet_BS[2]=jet2_BS;

        //llM selection
	


        bool isZlead=false;
                if( (jet1>=0) &&(llPt>(*jetPt)[jet1]))isZlead=true;
	
	FillHistos();
	
        } //END LOOP OVER ENTRIES

        jetVeto->clear();
        delete jetVeto;


}

int BaseAnalysis::Write(string outputFile){
	TFile *out=TFile::Open(outputFile.c_str(),"RECREATE");
	out->cd();
	for(map<string,TH1F*>::iterator it=histos.begin();it!=histos.end();it++) it->second->Write();
	for(map<string,TH2F*>::iterator it=histos2.begin();it!=histos2.end();it++) it->second->Write();
	for(map<string,RooUnfoldResponse*>::iterator it=response.begin();it!=response.end();it++) it->second->Write();
}

int BaseAnalysis::Sumw2(){
	for(map<string,TH1F*>::iterator it=histos.begin();it!=histos.end();it++) it->second->Sumw2();
	for(map<string,TH2F*>::iterator it=histos2.begin();it!=histos2.end();it++) it->second->Sumw2();
}

//constructor
BaseAnalysis::BaseAnalysis(){ JetSel=3;extraLabel="";}
//destructor
BaseAnalysis::~BaseAnalysis(){}
