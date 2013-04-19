#include "ZAnalysis.h"
#include "math.h"

int ZAnalysis::Smear(){
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

bool ZAnalysis::BaseSelection()
{
if(lepPt->size()<2)return false; //2 leptons
if( (*lepChId)[0]*(*lepChId)[1]!=CHID2)return false; //two muons
        //Triggers??
return true;
}

bool ZAnalysis::BaseSelectionGEN()
{
if(lepPtGEN->size()<2)return false; //2 leptons
int CHID2GEN=0;
	if(CHID2==-4) CHID2GEN=-13*13;
	if(CHID2==-1) CHID2GEN=-11*11;
	if(CHID2==-2) CHID2GEN=-11*13;
if( (*lepChIdGEN)[0]*(*lepChIdGEN)[1]!=CHID2GEN)return false; //two muons
        //Triggers??
return true;
}

int ZAnalysis::FillMiss()
{
response[ ("llPt"+extraLabel).c_str()]->Miss(llPtGEN,weight);
return 0;
}

int ZAnalysis::FillHistosGEN()
{
	//Gen Selection
	selGen=BaseSelectionGEN() && ( fabs(llMGEN-91)<llMCut) ;
	//llMCut GEN
	if( selGen && ( (*jetPtGEN)[0] > JetPtCut))
		selGen=true;
		else selGen=false;
	if(!selGen) return 0; //---------------------------------------------
	if(!BaseSelection()){//RECO
		FillMiss();
	}
	
return 0;	
}
int ZAnalysis::FillHistos()
{
//FILL MISS AT EACH SELECTION FAILED
//before llM cut 
histos[ ("nVtx" + extraLabel).c_str() ]->Fill(nVtx,weight);
histos[ ("nLeptons" + extraLabel).c_str() ]->Fill(nLeptons,weight);
histos[ ("nPhotons" + extraLabel).c_str() ]->Fill(nPhotons,weight);
histos[ ("nJets" + extraLabel).c_str() ]->Fill(nJetsVeto,weight);
histos[ ("rho" + extraLabel).c_str() ]->Fill(rho,weight);
histos[ ("rho25" + extraLabel).c_str() ]->Fill(rho25,weight);
histos[ ("llM" + extraLabel).c_str() ]->Fill(llM,weight);

if(!(fabs(llM-91)<llMCut))
	{
	if(selGen)FillMiss();
	return 0;} //---------------------------------------------
//after llM cut


histos[ ("llY" + extraLabel).c_str() ]->Fill(llY,weight);
histos[ ("llPhi" + extraLabel).c_str() ]->Fill(llPhi,weight);

histos[ ("lep0Pt" + extraLabel).c_str() ]->Fill( (*lepPt)[0],weight);
histos[ ("lep1Pt" + extraLabel).c_str() ]->Fill( (*lepPt)[1],weight);

histos[ ("lep0Eta" + extraLabel).c_str() ]->Fill( (*lepEta)[0],weight);
histos[ ("lep1Eta" + extraLabel).c_str() ]->Fill( (*lepEta)[1],weight);

histos[ ("lep0Phi" + extraLabel).c_str() ]->Fill( (*lepPhi)[0],weight);
histos[ ("lep1Phi" + extraLabel).c_str() ]->Fill( (*lepPhi)[1],weight);


if(jet_BS[0]>=0){
	histos[ ("jet0Pt" + extraLabel).c_str() ]->Fill( (*jetPt)[jet_BS[0]],weight);
	histos[ ("jet0Eta" + extraLabel).c_str() ]->Fill( (*jetEta)[jet_BS[0]],weight);
	histos[ ("jet0Phi" + extraLabel).c_str() ]->Fill( (*jetPhi)[jet_BS[0]],weight);
	histos[ ("jet0QGL" + extraLabel).c_str() ]->Fill( (*jetQGL)[jet_BS[0]],weight);
	histos[ ("jet0Btag" + extraLabel).c_str() ]->Fill( (*jetBtag)[jet_BS[0]],weight);
	
	histos[ ("llPtBetaStar" + extraLabel).c_str() ]->Fill( llPt,weight);
	histos[ ("llPt" + extraLabel).c_str() ]->Fill(llPt,weight);
	//fill Response
	{
	if(selGen){
		response[ ("llPt"+extraLabel).c_str()]->Fill(llPt,llPtGEN,weight);
		}
	else response[ ("llPt"+extraLabel).c_str()]->Fake(llPt,weight);
	}
}else{
	if(selGen) FillMiss();
	}

if(jet_BS[1]>=0){
	histos[ ("jet1Pt" + extraLabel).c_str() ]->Fill( (*jetPt)[jet_BS[1]],weight);
	histos[ ("jet1Eta" + extraLabel).c_str() ]->Fill( (*jetEta)[jet_BS[1]],weight);
	histos[ ("jet1Phi" + extraLabel).c_str() ]->Fill( (*jetPhi)[jet_BS[1]],weight);
	histos[ ("jet1QGL" + extraLabel).c_str() ]->Fill( (*jetQGL)[jet_BS[1]],weight);
	histos[ ("jet1Btag" + extraLabel).c_str() ]->Fill( (*jetBtag)[jet_BS[1]],weight);
}

}


int ZAnalysis::OpenFiles()
{
}

int ZAnalysis::CreateHistos()
{
	for(map<string,TH1F*>::iterator it=histos.begin();it!=histos.end();it++){delete it->second;}
histos.clear();
	for(map<string,TH2F*>::iterator it=histos2.begin();it!=histos2.end();it++){delete it->second;}
histos2.clear();

vector<string> labels;
labels.push_back("");
labels.push_back("_JESUP");
labels.push_back("_JESDN");

for(int i=0;i< int(labels.size());i++)
{

histos[("nVtx"+labels[i]).c_str()] 		= new TH1F("nVtx","nVtx;N_{Vtx};events",50,0,50);
histos[("nLeptons"+labels[i]).c_str()] 	= new TH1F("nLeptons","nLeptons;N_{lep};events",10,0,10);
histos[("nPhotons"+labels[i]).c_str()] 	= new TH1F("nPhotons","nPhotons;N_{#gamma};events",10,0,10);
histos[("nJets"+labels[i]).c_str()] 	= new TH1F("nJets","nJets;N_{jets};events",10,0,10);

//rho
histos[("rho"+labels[i]).c_str()]		= new TH1F("rho","rho;#rho;events",50,0,50);
histos[("rho25"+labels[i]).c_str()]		= new TH1F("rho25","rho25;#rho(2.5); events",50,0,50);

//DiBoson
histos[("llM"+labels[i]).c_str()]		= new TH1F("llM","llM;M^{ll};events",50,40,150);
histos[("llPt"+labels[i]).c_str()]		= new TH1F("llPt","llPt;P_{T}^{ll};events",100,0,450);
	response[("llPt"+labels[i]).c_str()] = new RooUnfoldResponse(histos["llPt"],histos["llPt"],("R_llPt"+labels[i]).c_str());
histos[("llY"+labels[i]).c_str()]		= new TH1F("llY","llY;Y^{ll};events",50,-5,5);
histos[("llPhi"+labels[i]).c_str()]		= new TH1F("llPhi","llPhi;#phi^{ll};events",50,-3.1416,3.1416);

//photons?
histos[("llgM"+labels[i]).c_str()]		= new TH1F("llgM","llgM;M^{ll#gamma};events/10GeV",200,0,2000); //todo
histos[("l1gM"+labels[i]).c_str()]		= new TH1F("l1gM","l1gM;M^{l1#gamma};events/10GeV",200,0,2000); //todo
histos[("l2gM"+labels[i]).c_str()]		= new TH1F("l2gM","l2gM;M^{l2#gamma};events/10GeV",200,0,2000); //todo

//lepton
histos[("lep0Pt"+labels[i]).c_str()]	= new TH1F("lep0Pt","lep0Pt;P_{T}^{1st lep};events",50,20,150);
histos[("lep0Eta"+labels[i]).c_str()]	= new TH1F("lep0Eta","lep0Eta;#eta^{1st lep};events",50,-5,5);
histos[("lep0Phi"+labels[i]).c_str()]	= new TH1F("lep0Phi","lep0Phi;#phi^{1st lep};events",50,-3.1416,3.1416);

histos[("lep1Pt"+labels[i]).c_str()]	= new TH1F("lep1Pt","lep1Pt;#P_{T}^{2nd lep};events",50,20,150);
histos[("lep1Eta"+labels[i]).c_str()]	= new TH1F("lep1Eta","lep1Eta;#eta^{2nd lep};events",50,-5,5);
histos[("lep1Phi"+labels[i]).c_str()]	= new TH1F("lep1Phi","lep1Phi;#phi^{2nd lep};events",50,-3.1416,3.1416);

//jets
histos[("jet0Pt"+labels[i]).c_str()]	= new TH1F("jet0Pt","jet0Pt;P_{T}^{1st jet};events",50,50,150);
histos[("jet0Eta"+labels[i]).c_str()]	= new TH1F("jet0Eta","jet0Eta;#eta^{1st jet};events",50,-5,5);
histos[("jet0Phi"+labels[i]).c_str()]	= new TH1F("jet0Phi","jet0Phi;#phi^{1st jet};events",50,-3.1416,3.1416);
histos[("jet0QGL"+labels[i]).c_str()]	= new TH1F("jet0QGL","jet0QGL;QGL^{1st jet};events",50,-1.001,1.001);
histos[("jet0Btag"+labels[i]).c_str()]	= new TH1F("jet0Btag","jet0Btag;Btag^{1st jet};events",50,-1.001,1.001);

histos[("jet1Pt"+labels[i]).c_str()]	= new TH1F("jet1Pt","jet1Pt;P_{T}^{2nd jet};events",50,50,150);
histos[("jet1Eta"+labels[i]).c_str()]	= new TH1F("jet1Eta","jet1Eta;#eta^{2nd jet};events",50,-5,5);
histos[("jet1Phi"+labels[i]).c_str()]	= new TH1F("jet1Phi","jet1Phi;#phi^{2nd jet};events",50,-3.1416,3.1416);
histos[("jet1QGL"+labels[i]).c_str()]	= new TH1F("jet1QGL","jet1QGL;QGL^{1st jet};events",50,-1.001,1.001);
histos[("jet1Btag"+labels[i]).c_str()]	= new TH1F("jet1Btag","jet1Btag;Btag^{2nd jet};events",50,-1.001,1.001);

// a bit of variables
histos[("jetLLDPhi0"+labels[i]).c_str()]	= new TH1F("JetLLDPhi0","JetLLDPhi0;#Delta#phi(Z,j_{1});events",50,0,3.1416); //todo
histos[("Sum3j"+labels[i]).c_str()]	= new TH1F("Sum3j","Sum3j;#sum_{ij#elem 1..3}d#phi_{ij};events",50,0,3.1416*2); //todo

//Introduce some cuts:
histos[("llPt_betaStar"+labels[i]).c_str()]		= new TH1F("llPt_betaStar","llPt;P_{T}^{ll} (betaStar on Jets);events",100,0,450);
histos[("jetLLDPhi0_PtZ_50"+labels[i]).c_str()]	= new TH1F("JetLLDPhi0_PtZ_50","JetLLDPhi0;#Delta#phi(Z,j_{1}) [P_{T}^{ll}>50 GeV];events",50,0,3.1416); //todo
histos[("llPt_nJets_3"+labels[i]).c_str()]		= new TH1F("llPt_nJets_3","llPt;P_{T}^{ll} (N_{jets} #geq 3);events",100,0,450); //todo
//histos["llPt_Analysis"]		= new TH1F("llPt_Analysis","llPt;P_{T}^{ll};events/[10GeV]",400,0,4000);

}
Sumw2();
}

//constructor
ZAnalysis::ZAnalysis():BaseAnalysis(){}
//destructor
ZAnalysis::~ZAnalysis(){}

