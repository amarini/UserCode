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
if(lepPtGEN==NULL)return false;
if(lepPtGEN->size()<2)return false; //2 leptons
if(lepChIdGEN->size()<2)return false; //2 leptons
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
#ifdef USE_ROOUNFOLD
response[ ("llPt"+extraLabel).c_str()]->Miss(llPtGEN,weight);
#endif
return 0;
}

int ZAnalysis::FillHistosGEN()
{
	//Gen Selection
	selGen=BaseSelectionGEN() && ( fabs(llMGEN-91)<llMCut) && (fabs( llYGEN) < ZYCut) ;
	//llMCut GEN
	if( selGen && ( (*jetPtGEN)[0] > JetPtCut))
		selGen=true;
		else selGen=false;
	if(!selGen) return 0; //---------------------------------------------
	if(!BaseSelection()){//RECO
		FillMiss();
	}
	
	//in any case fill histo gen
	if(extraLabel=="")
		histos[ "llPtGEN" ]->Fill(llPtGEN,weight);
	//fill always the Response matrix
	histos["R_llPt_GEN"+extraLabel]->Fill(llPtGEN,weight);    	
	
return 0;	
}
int ZAnalysis::FillHistos()
{
//FILL MISS AT EACH SELECTION FAILED
//before llM cut 
//debug=2;
if(debug>1)printf("Filling Histo\n");
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
	if(debug>1)printf("--> HERE C\n");
	histos[ ("jet0Pt" + extraLabel).c_str() ]->Fill( (*jetPt)[jet_BS[0]],weight);
	if(debug>1)printf("--> HERE B\n");
	histos[ ("jet0Eta" + extraLabel).c_str() ]->Fill( (*jetEta)[jet_BS[0]],weight);
	if(debug>1)printf("--> HERE A\n");
	histos[ ("jet0Phi" + extraLabel).c_str() ]->Fill( (*jetPhi)[jet_BS[0]],weight);
	if(debug>1)printf("--> HERE D\n");
	histos[ ("jet0QGL" + extraLabel).c_str() ]->Fill( (*jetQGL)[jet_BS[0]],weight);
	if(debug>1)printf("--> HERE E\n");
	histos[ ("jet0Btag" + extraLabel).c_str() ]->Fill( (*jetBtag)[jet_BS[0]],weight);
	if(debug>1)printf("--> HERE F\n");
	
	//histos[ ("llPtBetaStar" + extraLabel).c_str() ]->Fill( llPt,weight);
	if(debug>1)printf("--> HERE G\n");


	if( fabs(llY)<ZYCut)
		histos[ ("llPt" + extraLabel).c_str() ]->Fill(llPt,weight);
	else {if(selGen) FillMiss();}
	//fill Response
	{
	if(selGen){
		#ifdef USE_ROOUNFOLD
		response[ ("llPt"+extraLabel).c_str()]->Fill(llPt,llPtGEN,weight);
		#endif
                histos2["R_llPt_MATRIX"+extraLabel]->Fill(llPt,llPtGEN,weight);	
		}
	else{
		#ifdef USE_ROOUNFOLD
		 response[ ("llPt"+extraLabel).c_str()]->Fake(llPt,weight);
		#endif
	    }
                histos["R_llPt_RECO"+extraLabel]->Fill(llPt,weight); //fill regardless to selGen
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

if(debug>1)printf("Done Filling Histo\n");
}


int ZAnalysis::CreateHistos()
{

if(debug>0)printf("Creating %s -- %d \n",("nVtx"+extraLabel).c_str(),SmearType);

histos[("nVtx"+extraLabel).c_str()] 		= new TH1F(("nVtx"+extraLabel).c_str(),"nVtx;N_{Vtx};events",50,0,50);
histos[("nLeptons"+extraLabel).c_str()] 	= new TH1F(("nLeptons"+extraLabel).c_str(),"nLeptons;N_{lep};events",10,0,10);
histos[("nPhotons"+extraLabel).c_str()] 	= new TH1F(("nPhotons"+extraLabel).c_str(),"nPhotons;N_{#gamma};events",10,0,10);
histos[("nJets"+extraLabel).c_str()] 	= new TH1F(("nJets"+extraLabel).c_str(),"nJets;N_{jets};events",10,0,10);

//rho
histos[("rho"+extraLabel).c_str()]		= new TH1F(("rho"+extraLabel).c_str(),"rho;#rho;events",50,0,50);
histos[("rho25"+extraLabel).c_str()]		= new TH1F(("rho25"+extraLabel).c_str(),"rho25;#rho(2.5); events",50,0,50);

//DiBoson
histos[("llM"+extraLabel).c_str()]		= new TH1F(("llM"+extraLabel).c_str(),"llM;M^{ll};events",50,40,150);
histos[("llPt"+extraLabel).c_str()]		= new TH1F(("llPt"+extraLabel).c_str(),"llPt;P_{T}^{ll};events",100,0,450);
if(extraLabel=="") //no syst on the GEN LEVEL
	histos["llPtGEN"]		= new TH1F("llPtGEN","llPt;P_{T}^{ll};events",100,0,450);
if(debug>2)printf("Creating Response\n");
#ifdef USE_ROOUNFOLD
	response[("llPt"+extraLabel).c_str()] = new RooUnfoldResponse(histos["llPt"+extraLabel],histos["llPt"+extraLabel],("R_llPt"+extraLabel).c_str());
#endif
	histos["R_llPt_GEN"+extraLabel] 	=new TH1F(("R_llPt_GEN"+extraLabel).c_str(),"Response Matrix GEN",100,0,450);
	histos["R_llPt_RECO"+extraLabel]	=new TH1F(("R_llPt_RECO"+extraLabel).c_str(),"Response Matrix RECO",100,0,450);
	histos2["R_llPt_MATRIX"+extraLabel]	=new TH2F(("R_llPt_MATRIX"+extraLabel).c_str(),"Repsonse Matrix;RECO;GEN;events",100,0,450,100,0,450);
if(debug>2) printf("done 1 \n");
histos[("llY"+extraLabel).c_str()]		= new TH1F(("llY"+extraLabel).c_str(),"llY;Y^{ll};events",50,-5,5);
histos[("llPhi"+extraLabel).c_str()]		= new TH1F(("llPhi"+extraLabel).c_str(),"llPhi;#phi^{ll};events",50,-3.1416,3.1416);

//photons?
histos[("llgM"+extraLabel).c_str()]		= new TH1F(("llgM"+extraLabel).c_str(),"llgM;M^{ll#gamma};events/10GeV",200,0,2000); //todo
histos[("l1gM"+extraLabel).c_str()]		= new TH1F(("l1gM"+extraLabel).c_str(),"l1gM;M^{l1#gamma};events/10GeV",200,0,2000); //todo
histos[("l2gM"+extraLabel).c_str()]		= new TH1F(("l2gM"+extraLabel).c_str(),"l2gM;M^{l2#gamma};events/10GeV",200,0,2000); //todo

//lepton
histos[("lep0Pt"+extraLabel).c_str()]	= new TH1F(("lep0Pt"+extraLabel).c_str(),"lep0Pt;P_{T}^{1st lep};events",50,20,150);
histos[("lep0Eta"+extraLabel).c_str()]	= new TH1F(("lep0Eta"+extraLabel).c_str(),"lep0Eta;#eta^{1st lep};events",50,-5,5);
histos[("lep0Phi"+extraLabel).c_str()]	= new TH1F(("lep0Phi"+extraLabel).c_str(),"lep0Phi;#phi^{1st lep};events",50,-3.1416,3.1416);

histos[("lep1Pt"+extraLabel).c_str()]	= new TH1F(("lep1Pt"+extraLabel).c_str(),"lep1Pt;#P_{T}^{2nd lep};events",50,20,150);
histos[("lep1Eta"+extraLabel).c_str()]	= new TH1F(("lep1Eta"+extraLabel).c_str(),"lep1Eta;#eta^{2nd lep};events",50,-5,5);
histos[("lep1Phi"+extraLabel).c_str()]	= new TH1F(("lep1Phi"+extraLabel).c_str(),"lep1Phi;#phi^{2nd lep};events",50,-3.1416,3.1416);

//jets
histos[("jet0Pt"+extraLabel).c_str()]	= new TH1F(("jet0Pt"+extraLabel).c_str(),"jet0Pt;P_{T}^{1st jet};events",50,50,150);
histos[("jet0Eta"+extraLabel).c_str()]	= new TH1F(("jet0Eta"+extraLabel).c_str(),"jet0Eta;#eta^{1st jet};events",50,-5,5);
histos[("jet0Phi"+extraLabel).c_str()]	= new TH1F(("jet0Phi"+extraLabel).c_str(),"jet0Phi;#phi^{1st jet};events",50,-3.1416,3.1416);
histos[("jet0QGL"+extraLabel).c_str()]	= new TH1F(("jet0QGL"+extraLabel).c_str(),"jet0QGL;QGL^{1st jet};events",50,-1.001,1.001);
histos[("jet0Btag"+extraLabel).c_str()]	= new TH1F(("jet0Btag"+extraLabel).c_str(),"jet0Btag;Btag^{1st jet};events",50,-1.001,1.001);

histos[("jet1Pt"+extraLabel).c_str()]	= new TH1F(("jet1Pt"+extraLabel).c_str(),"jet1Pt;P_{T}^{2nd jet};events",50,50,150);
histos[("jet1Eta"+extraLabel).c_str()]	= new TH1F(("jet1Eta"+extraLabel).c_str(),"jet1Eta;#eta^{2nd jet};events",50,-5,5);
histos[("jet1Phi"+extraLabel).c_str()]	= new TH1F(("jet1Phi"+extraLabel).c_str(),"jet1Phi;#phi^{2nd jet};events",50,-3.1416,3.1416);
histos[("jet1QGL"+extraLabel).c_str()]	= new TH1F(("jet1QGL"+extraLabel).c_str(),"jet1QGL;QGL^{1st jet};events",50,-1.001,1.001);
histos[("jet1Btag"+extraLabel).c_str()]	= new TH1F(("jet1Btag"+extraLabel).c_str(),"jet1Btag;Btag^{2nd jet};events",50,-1.001,1.001);

// a bit of variables
histos[("jetLLDPhi0"+extraLabel).c_str()]	= new TH1F(("JetLLDPhi0"+extraLabel).c_str(),"JetLLDPhi0;#Delta#phi(Z,j_{1});events",50,0,3.1416); //todo
histos[("Sum3j"+extraLabel).c_str()]	= new TH1F(("Sum3j"+extraLabel).c_str(),"Sum3j;#sum_{ij#elem 1..3}d#phi_{ij};events",50,0,3.1416*2); //todo

//Introduce some cuts:
histos[("llPt_betaStar"+extraLabel).c_str()]		= new TH1F(("llPt_betaStar"+extraLabel).c_str(),"llPt;P_{T}^{ll} (betaStar on Jets);events",100,0,450);
histos[("jetLLDPhi0_PtZ_50"+extraLabel).c_str()]	= new TH1F(("JetLLDPhi0_PtZ_50"+extraLabel).c_str(),"JetLLDPhi0;#Delta#phi(Z,j_{1}) [P_{T}^{ll}>50 GeV];events",50,0,3.1416); //todo
histos[("llPt_nJets_3"+extraLabel).c_str()]		= new TH1F(("llPt_nJets_3"+extraLabel).c_str(),"llPt;P_{T}^{ll} (N_{jets} #geq 3);events",100,0,450); //todo
//histos["llPt_Analysis"]		= new TH1F(("llPt_Analysis"+extraLabel).c_str(),"llPt;P_{T}^{ll};events/[10GeV]",400,0,4000);

//}
if(debug>2) printf("Sumw2 \n");
Sumw2();
if(debug>2) printf("Done\n");
}

//constructor
ZAnalysis::ZAnalysis():BaseAnalysis(){}
//destructor
ZAnalysis::~ZAnalysis(){}

