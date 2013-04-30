#include "QGDumper.h"
#include "math.h"
#include "TLorentzVector.h"

int QGDumper::Smear(){
	if(SmearType==0) return 0;
	if(SmearType==1){ //JES UP
			for(int i=0;i<(jetPt->size());i++)
				(*jetPt)[i]*=(1+(*jetUNC)[i] );
			return 0;
			}
	if(SmearType==2){ //JES DN
			for(int i=0;i<(jetPt->size());i++)
				(*jetPt)[i]*=(1-(*jetUNC)[i] );
			return 0;
			}
}

bool QGDumper::BaseSelection()
{
if(lepPt->size()<2)return false; //2 leptons
if( (*lepChId)[0]*(*lepChId)[1]!=CHID2)return false; //two muons

return true;
}

int QGDumper::CreateHistos()
{

printf("Creating files: %s\n",outFileName.c_str());
outFile=TFile::Open(outFileName.c_str(),"RECREATE");
outFile->cd();
t_=new TTree(treeName.c_str(),treeName.c_str());

}


int QGDumper::SetActive(){
//-----
printf("Create Branches: llM\n");
t_->Branch("llM",		&llM,		"llM/F");
printf("Done\n");
t_->Branch("llPt",		&llPt,		"llPt/F");
t_->Branch("llPhi",		&llPhi,		"llPhi/F");
t_->Branch("eventWeight",	&eventWeight,	"eventWeight/F");
t_->Branch("PUWeight",		&PUWeight,	"PUWeight/F");
//--- this must be fixed in non vector style
PUWeight=1;eventWeight=1;
out_pdgIdPartJet0=0;
t_->Branch("QGLMLP",		&out_QGLMLP,	"QGLMLP/F");
t_->Branch("QGLHisto",		&out_QGLHisto,	"QGLHisto/F");
t_->Branch("PUReWeight",	&out_PUReWeight,"PUReWeight/F");
t_->Branch("ptJet0",		&out_ptJet0,	"ptJet0/F");
t_->Branch("etaJet0",		&out_etaJet0,	"etaJet0/F");
t_->Branch("phiJet0",		&out_phiJet0,	"phiJet0/F");
t_->Branch("pdgIdPartJet0",	&out_pdgIdPartJet0,"pdgIdPartJet0/I");

}

int QGDumper::FillHistos(){
	
	//Final Selection
	TLorentzVector j0,ll;
		ll.SetPtEtaPhiM(llPt,llEta,llPhi,llM);
	
	//llM selection
	if (!(fabs(llM-90)<llMCut)) return 0; //70-110
	// 1jet selection
	if (!(jet_BS[0]>=0)) return 0;
		j0.SetPtEtaPhiE((*jetPt)[jet_BS[0]],(*jetEta)[jet_BS[0]],(*jetPhi)[jet_BS[0]],(*jetE)[jet_BS[0]]);
	//if (!()) return 0;
	//if (!()) return 0;
	// jet balance - phi
	if (!( fabs(ll.DeltaPhi(j0)) > 3.1415 -0.5 )) return 0;
	//anti b-tag
	if (!((*jetBtag)[jet_BS[0]]<0.75)) return 0;
	//
	out_QGLMLP	=(*jetQGMLP)[jet_BS[0]];
	out_QGLHisto	=(*jetQGL)[jet_BS[0]];
	out_PUReWeight	=PUWeight/eventWeight;
	out_ptJet0	=(*jetPt)[jet_BS[0]];
	out_etaJet0	=(*jetEta)[jet_BS[0]];
	out_phiJet0	=(*jetPhi)[jet_BS[0]];
	out_pdgIdPartJet0=(*jetPdgId)[jet_BS[0]];

	t_->Fill();
	
	return 0;
}
int QGDumper::Write(string outFileName){
	outFile->cd();
	t_->Write();
	outFile->Close();

}

//constructor
QGDumper::QGDumper():BaseAnalysis(){t_=NULL;}
//destructor
QGDumper::~QGDumper(){}

