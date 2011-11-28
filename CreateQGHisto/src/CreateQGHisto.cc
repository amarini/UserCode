// -*- C++ -*-
//
// Package:    CreateQGHisto
// Class:      CreateQGHisto
// 
/**\class CreateQGHisto CreateQGHisto.cc CreateQGHisto/CreateQGHisto/src/CreateQGHisto.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrea Carlo Marini,32 2-C24,+41227676319,
//         Created:  Sun Nov 27 16:37:55 CET 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "libraries.h"
#include <vector>
#include <map>
#include <algorithm>
	class tris {   
			public:
			char a;
			double b ; 
			double c;
			tris(){}
			tris(char x, double y, double z){a=x; b=y; c=z;}
			};
	const inline bool operator<(tris x, tris y){ 
					if(x.a!=y.a)return x.a<y.a;
					if(x.b!=y.b)return x.b<y.b;
					return x.c<y.c;				
					}

const double ptBins[]={15,30,45,70};
const int nBins=sizeof(ptBins)/sizeof(double);
//
// class declaration
//

class CreateQGHisto : public edm::EDAnalyzer {
   public:
      explicit CreateQGHisto(const edm::ParameterSet&);
      ~CreateQGHisto();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      int evtNum_;
      // ---- run number ------------------------------------------------  
      int runNum_;
      // ---- lumi section ----------------------------------------------
      int lumi_;
      // ---- flag to identify real data --------------------------------
      int isRealData_;
      // ---- TFile Service ---------------------------------------------
      edm::Service<TFileService> fTFileService_;
      // ---- Tree ------------------------------------------------------
      TTree *tree_;
      // ---- Histos ----------------------------------------------------
      // ---- Variables
      int nCharged_,nNeutral_;
      double ptD_;
      int nCharged_nopu_,
		nNeutral_nopu_;
      double ptD_nopu_;
      //jet pt
      double pt_;
      double pt_nopu_;
      double eta_,phi_,e_;
      double eta_nopu_,phi_nopu_,e_nopu_;
      int jetNumber_;
      int jetNumber_nopu_;
      // ---- intime pu (MC truth)
      int pu_;
      // ----- rho
      double rho_,rho_nopu_;
      // ----- npv
      int npv_;
	//
	int pdgid_;
      //
	bool isNOPU_;
        bool isNormal_;
       // clear variables
	void ClearVariables();
      // histograms
	//low bin value of pt and rho
      std::map<tris,TH1D*> histos;

};
void CreateQGHisto::ClearVariables(){ pt_=0;
				pt_nopu_=0;
				nCharged_=0;
				nCharged_nopu_=0;
				nNeutral_=0;
				nNeutral_nopu_=0;
				ptD_=0;
				ptD_nopu_=0;
				isNOPU_=false;
				isNormal_=false;
				pdgid_=0;
				jetNumber_=-1;
				jetNumber_nopu_=-1;
				}

//
// constructors and destructor
//
CreateQGHisto::CreateQGHisto(const edm::ParameterSet& iConfig){}


CreateQGHisto::~CreateQGHisto(){}


//
// member functions
//

// ------------ method called for each event  ------------
void
CreateQGHisto::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;
isRealData_ = iEvent.isRealData() ? 1:0;
// ---- PU ----------------------------------------------------------
    Handle<vector<PileupSummaryInfo> > pileupInfo;
    iEvent.getByLabel("addPileupInfo", pileupInfo);
    vector<PileupSummaryInfo>::const_iterator PUI;
    pu_ = 0; //in-time pileup
    for(PUI = pileupInfo->begin(); PUI != pileupInfo->end(); ++PUI) {
      if (PUI->getBunchCrossing() == 0)
        pu_ += PUI->getPU_NumInteractions();
    }// PUI loop
  //---- Rho ------------------------------------------------------------
  Handle<double> rho;
  iEvent.getByLabel(edm::InputTag("kt6PFJets","rho"),rho);
  rho_=*rho;
  //---- Rho nopu --------------------------------------------------------
  Handle<double> rho_nopu;
  iEvent.getByLabel(edm::InputTag("kt6PFJets_nopu","rho"),rho_nopu);
  rho_nopu_=*rho_nopu;
  //---- reco vertices block --------------------------------------------
  edm::Handle<VertexCollection> vertices_;
  iEvent.getByLabel("offlinePrimaryVertices", vertices_);
  //const reco::Vertex *primVtx = &(*(vertices_.product()))[0];
  npv_=vertices_->size();
  //---- RUN
  runNum_ = iEvent.id().run();
  evtNum_ = iEvent.id().event();
  lumi_ = iEvent.luminosityBlock();

  //---- JETS

  edm::Handle < vector<pat::Jet> >  patJets;
  iEvent.getByLabel("patJetsAK5PF",patJets);
 	//take only the ns leading jets
  const reco::Jet *jet = 0;
	for(int i=0; i< int(patJets->size()) && i<2; i++)
		{
		ClearVariables();
		jet = (const reco::Jet*) ( & ((*patJets)[i]) );
  		const pat::Jet *patJet = dynamic_cast<const pat::Jet*>(&*jet);
		isNormal_=true;
		jetNumber_=i;
		pt_=jet->pt();
		phi_=jet->phi();
		eta_=jet->eta();
		e_=jet->energy();
		nCharged_=patJet->chargedMultiplicity();
		nNeutral_=patJet->neutralMultiplicity();
		double ptS=0,ptSquare=0;
		for(long int j=0;j<jet->nConstituents();++j)
        		   {
        		     ptS+=jet->getJetConstituentsQuick()[j]->pt();
        		     ptSquare+= TMath::Power(jet->getJetConstituentsQuick()[j]->pt(),2);
        		   }    
		ptD_=ptSquare/TMath::Power(ptS,2);
		 if(patJet->genParton())
                	{
                	pdgid_= (  (patJet->genParton())->pdgId());
               		}
            		else 	pdgid_=0;
		tree_->Fill();	
		}
  // ----- JETS NOPU
  edm::Handle < vector<pat::Jet> >  patJets_nopu;
  iEvent.getByLabel("patJetsAK5PFNOPU",patJets_nopu);
	for(int i=0; i<int( patJets->size()) && i<2; i++)
		{
		ClearVariables();
		jet = (const reco::Jet*) ( & ((*patJets_nopu)[i]) );
  		const pat::Jet *patJet = dynamic_cast<const pat::Jet*>(&*jet);
		isNOPU_=true;
		jetNumber_nopu_=i;
		pt_nopu_=jet->pt();
		phi_nopu_=jet->phi();
		eta_nopu_=jet->eta();
		e_nopu_=jet->energy();
		nCharged_nopu_=patJet->chargedMultiplicity();
		nNeutral_nopu_=patJet->neutralMultiplicity();
		double ptS=0,ptSquare=0;
		for(long int j=0;j<jet->nConstituents();++j)
        		   {
        		     ptS+=jet->getJetConstituentsQuick()[j]->pt();
        		     ptSquare+= TMath::Power(jet->getJetConstituentsQuick()[j]->pt(),2);
        		   }    
		ptD_nopu_=ptSquare/TMath::Power(ptS,2);
		 if(patJet->genParton())
                	{
                	pdgid_= (  (patJet->genParton())->pdgId());
               		}
            		else 	pdgid_=0;
		tree_->Fill();	
		}

//	%s_pt%.0f_%.0f_rho%d
}


// ------------ method called once each job just before starting event loop  ------------
void 
CreateQGHisto::beginJob()
{
tree_ = fTFileService_->make<TTree>("jetTree","jetTree");
tree_->Branch("pu",&pu_,"pu/I");
tree_->Branch("runNum",&runNum_,"runNum/I");
tree_->Branch("evtNum",&evtNum_,"evtNum/I");
tree_->Branch("lumi",&lumi_,"lumi/I");
tree_->Branch("nCharged",&nCharged_,"nCharged/I");
tree_->Branch("nCharged_nopu",&nCharged_nopu_,"nCharged_nopu/I");
tree_->Branch("nNeutral",&nNeutral_,"nNeutral/I");
tree_->Branch("nNeutral_nopu",&nNeutral_nopu_,"nNeutral_nopu/I");
tree_->Branch("npv",&npv_,"npv/I");
tree_->Branch("pdgid",&pdgid_,"pdgid/I");
tree_->Branch("jetNumber",&jetNumber_,"jetNumber/I");
tree_->Branch("jetNumber_nopu",&jetNumber_nopu_,"jetNumber_nopu/I");
//DOUBL
tree_->Branch("ptD",&ptD_,"ptD/D");
tree_->Branch("ptD_nopu",&ptD_nopu_,"ptD_nopu/D");
tree_->Branch("rho",&rho_,"rho/D");
tree_->Branch("rho_nopu",&rho_nopu_,"rho_nopu/D");
tree_->Branch("eta",&eta_,"eta/D");
tree_->Branch("phi",&phi_,"phi/D");
tree_->Branch("e",&e_,"e/D");
tree_->Branch("eta_nopu",&eta_nopu_,"eta_nopu/D");
tree_->Branch("phi_nopu",&phi_nopu_,"phi_nopu/D");
tree_->Branch("e_nopu",&e_nopu_,"e_nopu/D");
//bool
tree_->Branch("isNOPU",&isNOPU_,"isNOPU/B");
tree_->Branch("isNormal",&isNormal_,"isNormal/B");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CreateQGHisto::endJob() {}

// ------------ method called when starting to processes a run  ------------
void 
CreateQGHisto::beginRun(edm::Run const&, edm::EventSetup const&){}

// ------------ method called when ending the processing of a run  ------------
void 
CreateQGHisto::endRun(edm::Run const&, edm::EventSetup const&){}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CreateQGHisto::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CreateQGHisto::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CreateQGHisto::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CreateQGHisto);
