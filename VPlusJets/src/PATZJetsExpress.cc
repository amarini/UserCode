// -*- C++ -*-
//
// Package:    ZJetsExpress
// Class:      ZJetsExpress
// 
/**\class ZJetsExpress ZJetsExpress.cc ElectroWeakAnalysis/VPlusJets/src/ZJetsExpress.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  A. Marini, K. Kousouris,  K. Theofilatos
//         Created:  Mon Oct 31 07:52:10 CDT 2011
// $Id: PATZJetsExpress.cc,v 1.40 2013/02/20 14:25:17 webermat Exp $
//
//

// system include files
#include <memory>

#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>


#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1I.h"
#include "TH1F.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h" 
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h" 
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h" 
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRawData/interface/EcalRawDataCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "RecoCaloTools/Navigation/interface/EcalPreshowerNavigator.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

// Added photon stuff
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"

//Marco Isolation
#include "PFIsolation/SuperClusterFootprintRemoval/interface/SuperClusterFootprintRemoval.h"

#include "amarini/VPlusJets/interface/CiCPhotonID.h"

//
// class declaration
//
using namespace edm;
using namespace std;
using namespace reco;

class PATZJetsExpress : public edm::EDAnalyzer {
   public:
      explicit PATZJetsExpress(const edm::ParameterSet&);
      ~PATZJetsExpress();

   private:
      virtual void beginJob();
      virtual void beginRun(edm::Run const &, edm::EventSetup const& iSetup);
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob();
      bool checkTriggerName(string,std::vector<string>); //checks if string belongs to any of the vector<string>
      inline double getEffectiveAreaForMuons(const double& eta) const;
      inline double getEffectiveAreaForElectrons(const double& eta) const;
      // ---- method that builds the tree -------------------------------
      void buildTree();
      // ---- method that re-initializes the tree branches --------------
      void clearTree();
      // ---- structures for temporary object storage -------------------
      struct PARTICLE {
        // ---- momentum 4-vector ---------------------------------------
        TLorentzVector p4; 
        // ---- standard PF isolation uncorrected for PU ----------------
        float isoPFUnc;
        // ---- modified isolation Delta beta corrected -----------------
        float isoPFDb; 
        // ---- modified isolation rho corrected ------------------------
        float isoPFRho; 
        // ---- charge id (+/-1: electrons, +/-2: muons) ---------------- 
        int chid; 
        // ---- tight id ------------------------------------------------
        int id;
        // --- flag bit
        int bit;
        // --- float paerameters
        std::vector<float> parameters;
        
        // --- MP foot print removal isolation
        float isoFPRNeutral;
        float isoFPRCharged;
        float isoFPRPhoton;

        enum photonParameters {  passconv, pfIsoCH, pfIsoNH, pfIsoP, pfCic03P, pfCic03N, pfCic04P, pfCic04N , pfCic03Cg, pfCic04Cg, pfCic03Cb, pfCic04Cb, sieie, sieip, etaw, phiw,  r9, lR, s4, e25, sceta, ESEff,
        	hcalTowerSumEtConeDR04,
        	ecalRecHitSumEtConeDR04,
        	nTrkSolidConeDR04,
        	trkSumPtSolidConeDR04,
        	nTrkHollowConeDR04,
        	trkSumPtHollowConeDR04,
        	oldHadronicOverEm,
        	hadronicOverEm2012, nPhotonParameters } ;

        // --- sigma IetaIeta
        float sigmaIEtaIEta;
        float hadronicOverEm;
	// --- for electrons r9, for muons chi2/ndof for global track
	float r9_or_chi2ndof;
        
        // --- triggerObjectMatches
        pat::TriggerObjectStandAloneCollection triObjMatchF1;
        pat::TriggerObjectStandAloneCollection triObjMatchF2;
        pat::TriggerObjectStandAloneCollection triObjMatchF3mu;
        pat::TriggerObjectStandAloneCollection triObjMatchF3ele;
        pat::TriggerObjectStandAloneCollection triObjMatchF4;
        pat::TriggerObjectStandAloneCollection triObjMatchF5;
        pat::TriggerObjectStandAloneCollection triObjMatchF6;
      };
      struct GENPARTICLE {
        // ---- momentum 4-vector ---------------------------------------
        TLorentzVector p4;  
        // ---- pdgid ---------------------------------------------------
        int pdgId; 
        // ---- motherid ---------------------------------------------------
        int motherId; 
      };
      struct GENPARTICLEPHO {
        // ---- momentum 4-vector ---------------------------------------
        TLorentzVector p4;  
        // ---- pdgid ---------------------------------------------------
        int pdgId; 
        // ---- motherid ---------------------------------------------------
        int motherId; 
	float isoPtDR03;
	float isoEDR03;
	float isoPtDR04;
	float isoEDR04;
	float isoPtDR05;
	float isoEDR05;
      };
      struct JET {
        // ---- momentum 4-vector ---------------------------------------
        TLorentzVector p4; 
        // ---- tight id ------------------------------------------------
        int   id; 
        // ---- jet area (needed for pu subtraction) --------------------
        float area;
        // ---- track pt fraction associated with the PV ----------------
        float beta;
        // ---- btagger -------------------------------------------------
        float btag;
        // ---- btag extra info -----------------------------------------
        float taginfoNvtx;
        float taginfoNtracks;
        float taginfoVtxMass;
        // ---- MC flavour (0 in case of data) --------------------------
        int mcflavour;
        // ---- jet energy correction factor ----------------------------
        float jec; 
        // ---- jet energy uncertainty ----------------------------------
        float unc;
        // ---- charged hadron energy fraction --------------------------
        float chf;
        // ---- neutral hadron energy fraction --------------------------
        float nhf;
        // ---- photon energy fraction ---------------------------------- 
        float phf;
        // ---- muon energy fraction ------------------------------------ 
        float muf;
        // ---- electron energy fraction --------------------------------
        float elf;
        // ---- qgl ----------------------------------------------------
        float qgl;
        // ---- rms ----------------------------------------------------
        float rms;
        // ---- veto ---------------------------------------------------
        int veto;
      };
      class GENJET : public TLorentzVector {
         public:
            GENJET(float x=0., float y=0., float z=0., float t=0.) : TLorentzVector( x, y, z, t ){veto=0;}
            GENJET(TLorentzVector &v):TLorentzVector(v){veto=0;id=0;nparton=0;}
            virtual ~GENJET(){}
            int veto;
	    int id;
	    int nparton;
      };

      vector<float> *ComputeQGVariables(edm::View<pat::Jet>::const_iterator & jet,const Event& iEvent,int index);
      // ---- sorting rules ---------------------------------------------
      static bool lepSortingRule(PARTICLE x, PARTICLE y)                {return x.p4.Pt() > y.p4.Pt();}
      static bool lepSortingRuleGEN(GENPARTICLE x, GENPARTICLE y)       {return x.p4.Pt() > y.p4.Pt();}
      static bool phoSortingRuleGEN(GENPARTICLEPHO x, GENPARTICLEPHO y)       {return x.p4.Pt() > y.p4.Pt();}
      //static bool jetSortingRule(JET x, JET y)                      {return x.p4.Pt() > y.p4.Pt();}
      static bool p4SortingRule(TLorentzVector x, TLorentzVector y) {return x.Pt() > y.Pt();}
      // ---------- member data -----------------------------------------
      edm::Service<TFileService> fTFileService;
      TTree *myTree_;
      // ---- histogram to record the number of events ------------------
      TH1I  *hRecoLeptons_,*hRecoPhotons_,*hGenLeptons_,*hEvents_,*hWEvents_;
      TH1F  *hMuMuMass_,*hElElMass_,*hElMuMass_;
      TH1F  *hZMuMuMass_,*hZElElMass_;
      TH1F  *hMuMuMassWeighted_,*hElElMassWeighted_,*hElMuMassWeighted_;
      TH1F  *hElElEBMass_,*hLepLepMass_;
      TH1F  *hTriggerNames_,*hTriggerPass_;
      TH1F  *hGenPhotonPt_,*hGenPhotonMatchedPt_;
      TH1F  *hGenPhotonEta_,*hGenPhotonMatchedEta_;
      TH1F  *hGenMuonPt_,*hGenMuonMatchedPt_;
      TH1F  *hGenMuonEta_,*hGenMuonMatchedEta_;
      TH1F  *hGenElectronPt_,*hGenElectronMatchedPt_;
      TH1F  *hGenElectronEta_,*hGenElectronMatchedEta_;
      TH1F* hXSec_;
      // ---- simulated in-time pile-up ---------------------------------
      TH1D  *mcPU_;
      // ---- flag to set the JEC uncertainty object --------------------
      //bool mIsJECuncSet;
      // ---- jet energy corrector object -------------------------------
      //const JetCorrector *mJEC;
      // ---- jet energy uncertainty object -----------------------------
      //JetCorrectionUncertainty *mJECunc;
      // ---- trigger ---------------------------------------------------
      std::string   processName_;
      std::vector<std::string> triggerNames_,triggerNamesFull_;
      std::vector<std::string> triggerFamily1_;
      std::vector<std::string> triggerFamily2_;
      std::vector<std::string> triggerFamily3_;
      std::vector<std::string> triggerFamily4_;
      std::vector<std::string> triggerFamily5_;
      std::vector<std::string> triggerFamily6_;
      std::vector<std::string> triggerFamily7_;
      std::vector<std::string> triggerFamily8_;
      std::vector<std::string> prescaleDontAsk_;
      std::vector<unsigned int> triggerIndex_;
      edm::InputTag triggerResultsTag_;
      edm::InputTag triggerEventTag_;
      edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
      edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
      HLTConfigProvider hltConfig_;
      // ---- configurable parameters -----------------------------------
      bool          mIsMC,mOnlyMC,mReducedPh;
      int           mMinNjets,mGENType, mDressedRadius;
      double        mMinJetPt,mMaxJetEta,mMinLepPt,mMaxLepEta,mMaxCombRelIso03,mMaxCombRelIso04,mJetLepIsoR,mJetPhoIsoR,mMinLLMass,mMinPhoPt,mMinPhoPtId,mMaxPhoEta;
      edm::InputTag pfIsoValEleCH03Name,pfIsoValEleNH03Name,pfIsoValEleG03Name;
      //int      mGENCrossCleaning;
      //string        mJECserviceMC, mJECserviceDATA, mPayloadName;
      edm::InputTag mJetsName,mSrcRho,mSrcRho25;
      // ---- tree variables --------------------------------------------
      // ---- event number ----------------------------------------------
      ULong64_t eventNum_;
      // ---- run number ------------------------------------------------  
      int runNum_; 
      // ---- lumi section ----------------------------------------------
      int lumi_;
      // ---- flag to identify real data --------------------------------
      int isRealData_;
      // ---- trigger decisions -----------------------------------------
      std::vector<int> *fired_;
      int isTriggered_;
      // ---- trigger Matching -----------------------------------------
      int isTriggerMatchedFamily1_;
      int isTriggerMatchedFamily2_;
      int isTriggerMatchedFamily3_;
      int isTriggerMatchedFamily4_;
      int isTriggerMatchedFamily5_;
      int isTriggerMatchedFamily6_;
      // ---- L1 prescale -----------------------------------------------
      std::vector<int> *prescaleL1_;
      // ---- HLT prescale -----------------------------------------------
      std::vector<int> *prescaleHLT_;
      // ---- dilepton mass ---------------------------------------------
      float llM_,llMGEN_;
      // ---- dilepton rapidity -----------------------------------------  
      float llY_,llYGEN_;
      // ---- dilepton pseudorapidity -----------------------------------------  
      float llEta_,llEtaGEN_;
      // ---- dilepton pt -----------------------------------------------
      float llPt_,llPtGEN_;
      // ---- dilepton phi ----------------------------------------------
      float llPhi_,llPhiGEN_;
      // ---- dphi between the two leptons ------------------------------
      float llDPhi_,llDPhiGEN_; 
      // ---- dphi between jets and dilepton ----------------------------
      vector<float> *jetllDPhi_,*jetllDPhiGEN_;
      // ---- lepton kinematics -----------------------------------------
      vector<float> *lepPt_,*lepEta_,*lepPhi_,*lepE_,*lepPtGEN_,*lepEtaGEN_,*lepPhiGEN_,*lepEGEN_;
      vector<float> *lepSigmaIEtaIEta_,*lepHadronicOverEm_;
      vector<float> *lepR9orChi2ndof_;
      // ---- lepton properties ----------------------------------------- 
      vector<int>   *lepChId_,*lepId_,*lepChIdGEN_,*lepMatchedGEN_;
      vector<float> *lepPFIsoUnc_,*lepPFIsoDBCor_,*lepPFIsoRhoCor_,*lepMatchedDRGEN_;
      // ---- number of leptons -----------------------------------------
      int nLeptons_,nLeptonsGEN_;
      // ---- photon variables
      int nPhotons_;
      vector<float>* photonPt_;
      vector<float>* photonE_;
      vector<float>* photonEta_;
      vector<float>* photonPhi_;
      vector<int>* photonBit_;
      vector<float>* photonPassConversionVeto_;
      vector<float>* photonPfIsoChargedHad_;
      vector<float>* photonPfIsoNeutralHad_;
      vector<float>* photonPfIsoPhoton_;
      vector<float>* photonPfIsoPhotons03ForCic_;
      vector<float>* photonPfIsoNeutrals03ForCic_;
      vector<float>* photonPfIsoCharged03ForCicVtx0_;
      vector<float>* photonPfIsoCharged03BadForCic_;
      vector<float>* photonPfIsoPhotons04ForCic_;
      vector<float>* photonPfIsoNeutrals04ForCic_;
      vector<float>* photonPfIsoCharged04ForCicVtx0_;
      vector<float>* photonPfIsoCharged04BadForCic_;
      vector<float>* photonid_sieie_;
      vector<float>* photonid_sieip_;
      vector<float>* photonid_etawidth_;
      vector<float>* photonid_phiwidth_;
      vector<float>* photonid_r9_;
      vector<float>* photonid_lambdaRatio_;
      vector<float>* photonid_s4Ratio_;
      vector<float>* photonid_e25_;
      vector<float>* photonid_sceta_;
      vector<float>* photonid_ESEffSigmaRR_;
      vector<float>* photonid_hadronicOverEm_;
      vector<float>* photonid_hadronicOverEm2012_;
      vector<float>* photonhcalTowerSumEtConeDR04_;
      vector<float>* photonecalRecHitSumEtConeDR04_;
      vector<float>* photonnTrkSolidConeDR04_;
      vector<float>* photontrkSumPtSolidConeDR04_;
      vector<float>* photonnTrkHollowConeDR04_;
      vector<float>* photontrkSumPtHollowConeDR04_;
	
      //PhotonIso
      vector<float>* photonIsoFPRCharged_; //with FOOT PRINT REMOVAL FROM MP
      vector<float>* photonIsoFPRNeutral_;
      vector<float>* photonIsoFPRPhoton_;
  
      vector<float> *jetPhotonDPhi_;
//       vector<float> *photonPar_;
      int nPhotonsGEN_;
      float photonEGEN_;
      float photonPtGEN_;
      float photonEtaGEN_;
      float photonPhiGEN_;
      float photonIsoPtDR03GEN_;
      float photonIsoEDR03GEN_;
      float photonIsoPtDR04GEN_;
      float photonIsoEDR04GEN_;
      float photonIsoPtDR05GEN_;
      float photonIsoEDR05GEN_;
      int photonMotherIdGEN_;
      float photonRECODRGEN_;
      // ---- VBParton variables
      float VBPartonM_;
      float VBPartonE_;
      float VBPartonPt_;
      float VBPartonEta_;
      float VBPartonPhi_;
      int VBPartonDM_; // decay mode
      // ---- jet kinematics --------------------------------------------
      vector<float> *jetPt_,*jetEta_,*jetPhi_,*jetE_,*jetPtGEN_,*jetEtaGEN_,*jetPhiGEN_,*jetEGEN_;
      //---- forward jets------  save two leading jets
      vector<float> *fwjetPt_,*fwjetEta_,*fwjetPhi_,*fwjetE_;
      int fwnJets_;
      vector<int> *jetIdGEN_,*jetNpartonsGEN_;
      vector<int> *jetVetoGEN_;
      // ---- other jet properties --------------------------------------
      vector<float> *jetBeta_,*jetBtag_,*jetTagInfoNVtx_,*jetTagInfoNTracks_,*jetTagInfoVtxMass_,*jetArea_,*jetJEC_,*jetUNC_,*jetQGL_,*jetRMS_;
      vector<int> *jetVeto_, *jetMCFlavour_;
      // ---- QG ----
      vector<float> *QGVars_;
      // ---- number of jets --------------------------------------------
      int nJets_,nJetsGEN_;
      // int nRJets_;
      // ---- pf met ----------------------------------------------------
      float pfmet_;
      // ---- pf sumEt --------------------------------------------------
      float pfSumEt_;
      // ---- sum of the pt of all status 3 parton: madgraph cut? ---------------------- 
      float HTParSum_;
      // ---- pf met phi ------------------------------------------------
      float pfmetPhi_;
      // ---- pt of the hadronic recoil ---------------------------------
      float pfhadPt_; 
      // ---- pf pt density ---------------------------------------------
      float rho_;
      float rho25_;
      // ---- reconstructed vertices' prperties -------------------------
      vector<float> *vtxZ_,*vtxNdof_;
      // ---- number of good reconstructed vertices ---------------------
      int   nVtx_;
      // ---- number of simulated pu interactions -----------------------
      int   puINT_,puOOT_,puTrueINT_,puTrueOOT_;
      // ---- MC weight
      float mcWeight_;

      // PF isolation calculator for photon
      PFIsolationEstimator isolator;

      CiCPhotonID* cicPhotonId;

      // new H/E calculator for photon
      ElectronHcalHelper::Configuration hcalCfg;
      ElectronHcalHelper *hcalHelper;
	    //xSec variables
	    //---- xSec information
	    TH1F*uXSec_; 

      std::vector<float> getESHits(double X, double Y, double Z, std::map<DetId, EcalRecHit> rechits_map, const CaloGeometry& geometry, CaloSubdetectorTopology *topology_p, int row=0);
      std::vector<float> getESShape(std::vector<float> ESHits0);
  
};
//
// class implemetation
//
// ---- constructor -----------------------------------------------------
PATZJetsExpress::PATZJetsExpress(const ParameterSet& iConfig)
{
  mMinNjets          = iConfig.getParameter<int>                       ("minNjets");   
  mJetLepIsoR        = iConfig.getParameter<double>                    ("jetLepIsoRadius");
  mJetPhoIsoR        = iConfig.getParameter<double>                    ("jetLepIsoRadius");
  mMinJetPt          = iConfig.getParameter<double>                    ("minJetPt");
  mMaxJetEta         = iConfig.getParameter<double>                    ("maxJetEta");
  mMinLepPt          = iConfig.getParameter<double>                    ("minLepPt");
  mMaxLepEta         = iConfig.getParameter<double>                    ("maxLepEta");
  mMinPhoPt          = iConfig.getParameter<double>                    ("minPhoPt");
  mMinPhoPtId       = iConfig.getParameter<double>                    ("minPhoPtId");
  mMaxPhoEta         = iConfig.getParameter<double>                    ("maxPhoEta");
  mMinLLMass         = iConfig.getParameter<double>                    ("minLLMass");
  mMaxCombRelIso03   = iConfig.getParameter<double>                    ("maxCombRelIso03");
  mMaxCombRelIso04   = iConfig.getParameter<double>                    ("maxCombRelIso04");
  mJetsName          = iConfig.getParameter<edm::InputTag>             ("jets");
  mSrcRho            = iConfig.getParameter<edm::InputTag>             ("srcRho");
  mSrcRho25          = iConfig.getParameter<edm::InputTag>             ("srcRho25");
  //mJECserviceDATA    = iConfig.getParameter<std::string>               ("jecServiceDATA");
  //mJECserviceMC      = iConfig.getParameter<std::string>               ("jecServiceMC");
  //mPayloadName       = iConfig.getParameter<std::string>               ("payload");
  processName_       = iConfig.getParameter<std::string>               ("processName");
  triggerNames_      = iConfig.getParameter<std::vector<std::string> > ("triggerName");
  triggerFamily1_    = iConfig.getParameter<std::vector<std::string> > ("triggerFamily1");
  triggerFamily2_    = iConfig.getParameter<std::vector<std::string> > ("triggerFamily2");
  triggerFamily3_    = iConfig.getParameter<std::vector<std::string> > ("triggerFamily3");
  triggerFamily4_    = iConfig.getParameter<std::vector<std::string> > ("triggerFamily4");
  triggerFamily5_    = iConfig.getParameter<std::vector<std::string> > ("triggerFamily5");
  triggerFamily6_    = iConfig.getParameter<std::vector<std::string> > ("triggerFamily6");
  triggerFamily7_    = iConfig.getParameter<std::vector<std::string> > ("triggerFamily7");
  triggerFamily8_    = iConfig.getParameter<std::vector<std::string> > ("triggerFamily8");
  prescaleDontAsk_   = iConfig.getParameter<std::vector<std::string> > ("prescaleDontAsk");
  triggerResultsTag_ = iConfig.getParameter<edm::InputTag>             ("triggerResults");
  triggerEventTag_   = iConfig.getParameter<edm::InputTag>             ("triggerEvent");   
  mGENType           = iConfig.getParameter<int>                       ("GENType");
  mDressedRadius     = iConfig.getParameter<double>                    ("dressedRadius");
  mOnlyMC	     = iConfig.getParameter<bool>		       ("OnlyMC");
  mReducedPh         = iConfig.getParameter<bool>		       ("ReducedPh");
  pfIsoValEleCH03Name= iConfig.getParameter<edm::InputTag>("pfIsoValEleCH03");
  pfIsoValEleNH03Name= iConfig.getParameter<edm::InputTag>("pfIsoValEleNH03");
  pfIsoValEleG03Name = iConfig.getParameter<edm::InputTag>("pfIsoValEleG03");

  //mGENCrossCleaning  = iConfig.getParameter<int>                       ("GENCrossCleaning"); //0: Nothing 1: Leptons 2: photons (Bit) 4:

  // initializing the PF isolation estimator for photon isolation 2012
  isolator.initializePhotonIsolation(kTRUE);
  isolator.setConeSize(0.3);

  cicPhotonId = new CiCPhotonID(iConfig); // chiara

  // initializing the new H/E calculator for photon 2012 H/E
  hcalCfg.hOverEConeSize = 0.15;
  hcalCfg.useTowers = true;
  hcalCfg.hcalTowers = edm::InputTag("towerMaker");
  hcalCfg.hOverEPtMin = 0;
  hcalHelper = new ElectronHcalHelper(hcalCfg);

}
// ---- destructor ------------------------------------------------------
PATZJetsExpress::~PATZJetsExpress()
{
  delete cicPhotonId;
  delete hcalHelper;
}
// ---
bool PATZJetsExpress::checkTriggerName(string aString,std::vector<string> aFamily)
{
  bool result(false);	
  for(unsigned int i=0;i<aFamily.size();i++) // checks if any of the aFamily strings contains aString
  {
    size_t found = aFamily[i].find(aString);
    if(found!=string::npos)result=true;
  }
  return result;
} 


std::vector<float> PATZJetsExpress::getESHits(double X, double Y, double Z, std::map<DetId, EcalRecHit> rechits_map, const CaloGeometry& geometry, CaloSubdetectorTopology *topology_p, int row) {
  std::vector<float> esHits;

  const GlobalPoint point(X,Y,Z);

  const CaloSubdetectorGeometry *geometry_p ;
  geometry_p = geometry.getSubdetectorGeometry (DetId::Ecal,EcalPreshower) ;

  DetId esId1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 1);
  DetId esId2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 2);
  ESDetId esDetId1 = (esId1 == DetId(0)) ? ESDetId(0) : ESDetId(esId1);
  ESDetId esDetId2 = (esId2 == DetId(0)) ? ESDetId(0) : ESDetId(esId2);

  std::map<DetId, EcalRecHit>::iterator it;
  ESDetId next;
  ESDetId strip1;
  ESDetId strip2;

  strip1 = esDetId1;
  strip2 = esDetId2;
    
  EcalPreshowerNavigator theESNav1(strip1, topology_p);
  theESNav1.setHome(strip1);
    
  EcalPreshowerNavigator theESNav2(strip2, topology_p);
  theESNav2.setHome(strip2);

  if (row == 1) {
    if (strip1 != ESDetId(0)) strip1 = theESNav1.north();
    if (strip2 != ESDetId(0)) strip2 = theESNav2.east();
  } else if (row == -1) {
    if (strip1 != ESDetId(0)) strip1 = theESNav1.south();
    if (strip2 != ESDetId(0)) strip2 = theESNav2.west();
  }

  // Plane 1 
  if (strip1 == ESDetId(0)) {
    for (int i=0; i<31; ++i) esHits.push_back(0);
  } else {

    it = rechits_map.find(strip1);
    if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
    else esHits.push_back(0);
    //cout<<"center : "<<strip1<<" "<<it->second.energy()<<endl;      

    // east road 
    for (int i=0; i<15; ++i) {
      next = theESNav1.east();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"est "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"east "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // west road 
    theESNav1.setHome(strip1);
    theESNav1.home();
    for (int i=0; i<15; ++i) {
      next = theESNav1.west();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"west "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"west "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }

  if (strip2 == ESDetId(0)) {
    for (int i=0; i<31; ++i) esHits.push_back(0);
  } else {

    it = rechits_map.find(strip2);
    if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
    else esHits.push_back(0);
    //cout<<"center : "<<strip2<<" "<<it->second.energy()<<endl;      

    // north road 
    for (int i=0; i<15; ++i) {
      next = theESNav2.north();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"north "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;  
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"north "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // south road 
    theESNav2.setHome(strip2);
    theESNav2.home();
    for (int i=0; i<15; ++i) {
      next = theESNav2.south();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"south "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"south "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }

  return esHits;
}

std::vector<float> PATZJetsExpress::getESShape(std::vector<float> ESHits0)
{
  std::vector<float> esShape;
    
  const int nBIN = 21;
  float esRH_F[nBIN];
  float esRH_R[nBIN];
  for (int idx=0; idx<nBIN; idx++) {
    esRH_F[idx] = 0.;
    esRH_R[idx] = 0.;
  }

  for(int ibin=0; ibin<((nBIN+1)/2); ibin++) {
    if (ibin==0) {
      esRH_F[(nBIN-1)/2] = ESHits0[ibin];
      esRH_R[(nBIN-1)/2] = ESHits0[ibin+31];
    } else {
      esRH_F[(nBIN-1)/2+ibin] = ESHits0[ibin];
      esRH_F[(nBIN-1)/2-ibin] = ESHits0[ibin+15];
      esRH_R[(nBIN-1)/2+ibin] = ESHits0[ibin+31];
      esRH_R[(nBIN-1)/2-ibin] = ESHits0[ibin+31+15];
    }
  } 

  // ---- Effective Energy Deposit Width ---- //
  double EffWidthSigmaXX = 0.;
  double EffWidthSigmaYY = 0.;
  double totalEnergyXX   = 0.;
  double totalEnergyYY   = 0.;
  double EffStatsXX      = 0.;
  double EffStatsYY      = 0.;
  for (int id_X=0; id_X<21; id_X++) {
    totalEnergyXX  += esRH_F[id_X];
    EffStatsXX     += esRH_F[id_X]*(id_X-10)*(id_X-10);
    totalEnergyYY  += esRH_R[id_X];
    EffStatsYY     += esRH_R[id_X]*(id_X-10)*(id_X-10);
  }
  EffWidthSigmaXX  = (totalEnergyXX>0.)  ? sqrt(fabs(EffStatsXX  / totalEnergyXX))   : 0.;
  EffWidthSigmaYY  = (totalEnergyYY>0.)  ? sqrt(fabs(EffStatsYY  / totalEnergyYY))   : 0.;

  esShape.push_back(EffWidthSigmaXX);
  esShape.push_back(EffWidthSigmaYY);
    
  return esShape;
}

// ---- method called once each job just before starting event loop -----
void PATZJetsExpress::beginJob()
{
  // ---- create the objects --------------------------------------------
  hXSec_          	 = fTFileService->make<TH1F>("XSec", "XSec",20,-.5,19.5);hXSec_->Sumw2();
  hRecoLeptons_          = fTFileService->make<TH1I>("RecoLeptons", "RecoLeptons",6,0,6);hRecoLeptons_->Sumw2();
  hRecoPhotons_          = fTFileService->make<TH1I>("RecoPhotons", "RecoPhotons",6,0,6);hRecoPhotons_->Sumw2();  // chiara
  hGenLeptons_           = fTFileService->make<TH1I>("GenLeptons", "GenLeptons",6,0,6);hGenLeptons_->Sumw2();
  hEvents_               = fTFileService->make<TH1I>("Events", "Events",1,0,1);hEvents_->Sumw2();
  hMuMuMass_             = fTFileService->make<TH1F>("MuMuMass", "MuMuMass",40,71,111);hMuMuMass_->Sumw2();
  hElElMass_             = fTFileService->make<TH1F>("ElElMass", "ElElMass",40,71,111);hElElMass_->Sumw2();
  hElMuMass_             = fTFileService->make<TH1F>("ElMuMass", "ElMuMass",40,71,111);hElMuMass_->Sumw2();
  hElElEBMass_           = fTFileService->make<TH1F>("ElElEBMass", "ElElEBMass",40,71,111);hElElEBMass_->Sumw2();
  hLepLepMass_           = fTFileService->make<TH1F>("LepLepMass", "LepLepMass",980,40,2000);hLepLepMass_->Sumw2();
  hTriggerNames_         = fTFileService->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  hTriggerNames_ ->SetBit(TH1::kCanRebin);
  for(unsigned i=0;i<triggerNames_.size();i++)
    hTriggerNames_->Fill(triggerNames_[i].c_str(),1);
  hTriggerPass_          = fTFileService->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  hTriggerPass_  ->SetBit(TH1::kCanRebin);
  mcPU_                  = fTFileService->make<TH1D>("mcPU", "mcPU",40,0,40);

  hWEvents_              = fTFileService->make<TH1I>("WEvents", "Weighted Events",1,0,1);hWEvents_->Sumw2();

  hMuMuMassWeighted_     = fTFileService->make<TH1F>("MuMuMassWeighted", "MuMuMassWeighted",40,71,111);hMuMuMassWeighted_->Sumw2();
  hElElMassWeighted_     = fTFileService->make<TH1F>("ElElMassWeighted", "ElElMassWeighted",40,71,111);hElElMassWeighted_->Sumw2();
  hElMuMassWeighted_     = fTFileService->make<TH1F>("ElMuMassWeighted", "ElMuMassWeighted",40,71,111);hElMuMassWeighted_->Sumw2();

  hGenPhotonPt_          = fTFileService->make<TH1F>("hGenPhotonPt", "hGenPhotonPt;photon p_{T} [GeV];events",300,50,350);hGenPhotonPt_->Sumw2();
  hGenPhotonEta_         = fTFileService->make<TH1F>("hGenPhotonEta","hGenPhotonPt;photon p_{T} [GeV];events",300,-3.0,3.0);hGenPhotonEta_->Sumw2();
  hGenPhotonMatchedPt_   = fTFileService->make<TH1F>("hGenPhotonMatchedPt", "hGenPhotonMatchedPt;photon p_{T} [GeV];events",300,50,350);
  hGenPhotonMatchedEta_  = fTFileService->make<TH1F>("hGenPhotonMatchedEta","hGenPhotonMatchedPt;photon p_{T} [GeV];events",300,-3.0,3.0);
  hGenPhotonMatchedPt_->Sumw2();
  hGenPhotonMatchedEta_->Sumw2();

  hGenMuonPt_            = fTFileService->make<TH1F>("hGenMuonPt", "hGenMuonPt;muon p_{T} [GeV];events",300,0,300);hGenMuonPt_->Sumw2();
  hGenMuonEta_           = fTFileService->make<TH1F>("hGenMuonEta","hGenMuonPt;muon p_{T} [GeV];events",300,-3.0,3.0);hGenMuonEta_->Sumw2();
  hGenMuonMatchedPt_     = fTFileService->make<TH1F>("hGenMuonMatchedPt", "hGenMuonMatchedPt;muon p_{T} [GeV];events",300,0,300);
  hGenMuonMatchedEta_    = fTFileService->make<TH1F>("hGenMuonMatchedEta","hGenMuonMatchedPt;muon p_{T} [GeV];events",300,-3.0,3.0);
  hGenMuonMatchedPt_->Sumw2();
  hGenMuonMatchedEta_->Sumw2();

  hGenElectronPt_            = fTFileService->make<TH1F>("hGenElectronPt", "hGenElectronPt;muon p_{T} [GeV];events",300,0,300);hGenElectronPt_->Sumw2();
  hGenElectronEta_           = fTFileService->make<TH1F>("hGenElectronEta","hGenElectronPt;muon p_{T} [GeV];events",300,-3.0,3.0);hGenElectronEta_->Sumw2();
  hGenElectronMatchedPt_     = fTFileService->make<TH1F>("hGenElectronMatchedPt", "hGenElectronMatchedPt;muon p_{T} [GeV];events",300,0,300);
  hGenElectronMatchedEta_    = fTFileService->make<TH1F>("hGenElectronMatchedEta","hGenElectronMatchedPt;muon p_{T} [GeV];events",300,-3.0,3.0);
  hGenElectronMatchedPt_->Sumw2();
  hGenElectronMatchedEta_->Sumw2();

  hZMuMuMass_             = fTFileService->make<TH1F>("ZMuMuMass", "ZMuMuMass",40,71,111);hZMuMuMass_->Sumw2();
  hZElElMass_             = fTFileService->make<TH1F>("ZElElMass", "ZElElMass",40,71,111);hZElElMass_->Sumw2();


  

  myTree_                = fTFileService->make<TTree>("events", "events");
  // ---- build the tree ------------------------------------------------
  buildTree();
  // ---- set the jec uncertainty flag ----------------------------------
  //mIsJECuncSet = false; 
}
// ---- method called everytime there is a new run ----------------------
void PATZJetsExpress::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  if(!mOnlyMC){
    if (triggerNames_.size() > 0) {
      bool changed(true);
      if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
	if (changed) {
	  triggerNamesFull_.clear();
	  // check if trigger names in (new) config
	  cout<<"New trigger menu found !!!"<<endl;
	  triggerIndex_.clear(); 
	  const unsigned int n(hltConfig_.size());
	  for(unsigned itrig=0;itrig<triggerNames_.size();itrig++) {
	    bool found(false);
	    for(unsigned iName=0;iName<n;iName++) {
	      std::string ss(hltConfig_.triggerName(iName));
	      if (int(ss.find(triggerNames_[itrig])) > -1) {
		triggerNamesFull_.push_back(ss);
		found = true;
		continue;
	      } 
	    }
	    if (!found) {
	      triggerNamesFull_.push_back("");
	    }
	    triggerIndex_.push_back(hltConfig_.triggerIndex(triggerNamesFull_[itrig]));
	    cout<<triggerNames_[itrig]<<" "<<triggerNamesFull_[itrig]<<" "<<triggerIndex_[itrig]<<" ";  
	    if (triggerIndex_[itrig] >= n)
	      cout<<"does not exist in the current menu"<<endl;
	    else
	      cout<<"exists"<<endl;
	  }// trigger names loop
	}
      } 
      else {
	cout << "ProcessedTreeProducer::analyze:"
	     << " config extraction failure with process name "
	     << processName_ << endl;
      }
    }
  }
  //isRealData_ Not set yet
  if(iRun.run()<100){
  	Handle<GenRunInfoProduct> genruninfo;  
  	iRun.getByLabel("generator",genruninfo);
  	cout<<"in begin Run  intXS/extXSLO/extXSNLO "<<genruninfo->internalXSec().value()<<"/"<<genruninfo->externalXSecLO().value()<<"/"<<genruninfo->externalXSecNLO().value()<<endl;	
  	hXSec_->Fill(0. ,genruninfo->internalXSec().value()/pow(genruninfo->internalXSec().error(),2)    );
  	hXSec_->Fill(1 ,1./pow(genruninfo->internalXSec().error(),2) );
  	
  	hXSec_->Fill(2 ,genruninfo->externalXSecLO().value()/pow(genruninfo->externalXSecLO().error(),2)  );
  	hXSec_->Fill(3 ,1./pow(genruninfo->externalXSecLO().error(),2)  );
  	
  	hXSec_->Fill(4 ,genruninfo->externalXSecNLO().value()/pow(genruninfo->externalXSecNLO().error(),2) );
  	hXSec_->Fill(5 ,1./pow(genruninfo->externalXSecNLO().error(),2) );
  	
  	hXSec_->Fill(6 , 1 );
  	
  	hXSec_->Fill(7 ,pow(genruninfo->internalXSec().value(),1) );
  	hXSec_->Fill(6 ,pow(genruninfo->externalXSecLO().value(),1)  );
  	hXSec_->Fill(9 ,pow(genruninfo->externalXSecNLO().value(),1) );
  	
  	hXSec_->Fill(10 ,pow(genruninfo->internalXSec().value(),2) );
  	hXSec_->Fill(11 ,pow(genruninfo->externalXSecLO().value(),2)  );
  	hXSec_->Fill(12 ,pow(genruninfo->externalXSecNLO().value(),2) );
	}
}
// ---- event loop ------------------------------------------------------
void PATZJetsExpress::analyze(const Event& iEvent, const EventSetup& iSetup)
{
  if(!mOnlyMC){
    // setup for new H/E calculator for 2012
    hcalHelper->checkSetup(iSetup);
    hcalHelper->readEvent(iEvent);
  }
  // ---- event counter -------------------------------------------------
  hEvents_->Fill(0.5);  
  // ---- initialize the tree branches ----------------------------------
  clearTree();
  isRealData_ = iEvent.isRealData() ? 1:0;

  int nJets_lepVeto=0;
  int nJets_phoVeto=0;

  int nJets_lepVetoGEN=0;
  int nJets_phoVetoGEN=0;


  // ----  MC truth block -----------------------------------------------
  vector<GENPARTICLE>      myGenLeptons;
  //vector<TLorentzVector> myGenJets;  
  vector<GENJET> myGenJets;  
  vector<GENPARTICLEPHO> myGenPhotons;  
  TLorentzVector VBParton(0,0,0,0);
  if (!isRealData_) {
    // ---- PU ----------------------------------------------------------
    if(!mOnlyMC){
      Handle<vector<PileupSummaryInfo> > pileupInfo;
      iEvent.getByLabel("addPileupInfo", pileupInfo);
      vector<PileupSummaryInfo>::const_iterator PUI;
      puOOT_ = 0;
      puINT_ = 0;
      puTrueINT_= 0;
      puTrueOOT_= 0;
      for(PUI = pileupInfo->begin(); PUI != pileupInfo->end(); ++PUI) {
	if (PUI->getBunchCrossing() == 0){
	  puINT_ += PUI->getPU_NumInteractions();     
      	  puTrueINT_+=PUI->getTrueNumInteractions();
	  }
	else{
	  puOOT_ += PUI->getPU_NumInteractions();
      	  puTrueOOT_+=PUI->getTrueNumInteractions();
	  }

      }// PUI loop
      mcPU_->Fill(puINT_);
    }
    // --- MC weight
    Handle<GenEventInfoProduct> geninfo;  
    if(!isRealData_)
	{
    	iEvent.getByLabel("generator",geninfo);
    	mcWeight_ = geninfo->weight();
    	hWEvents_->Fill(0.5,mcWeight_);
	}
    // --- Gen Jets
    Handle<GenJetCollection> genjets;
    if(!isRealData_)
    	iEvent.getByLabel("ak5GenJets",genjets);
    Handle<GenParticleCollection> gen;
    if(!isRealData_)
    	iEvent.getByLabel("genParticles", gen);
    GenParticleCollection::const_iterator i_gen;
    GenParticleCollection::const_iterator j_gen;
    GenJetCollection::const_iterator i_genjet;

    int VBPartonDM=0;
    // ---- loop over the gen particles ---------------------------------
    for(i_gen = gen->begin(); i_gen != gen->end(); i_gen++) {
    // save MC Vector Boson partons momenta and their decay mode, when applicable
    if( (i_gen->pdgId() ==23 || i_gen->pdgId()==22) && i_gen->status()==3) {
      float Px = i_gen->p4().Px();
      float Py = i_gen->p4().Py();
      float Pz = i_gen->p4().Pz();
      float E  = i_gen->p4().E();
      VBParton.SetPxPyPzE(Px,Py,Pz,E);
      const GenParticle* gen_Dau;
      for(int kk = 0; kk < int(i_gen-> numberOfDaughters()); ++kk) {
	gen_Dau = static_cast<const GenParticle*>(i_gen->daughter(kk)); // find daughter
	if(gen_Dau->pdgId()==23) continue;
	VBPartonDM = abs(gen_Dau->pdgId());
      }
      if(VBParton.Pt()>0) {
	VBPartonDM_  = VBPartonDM;
	VBPartonE_   = VBParton.E();
	VBPartonPt_  = VBParton.Pt();
	VBPartonEta_ = VBParton.Eta();
	VBPartonPhi_ = VBParton.Phi();
       	VBPartonM_   = VBParton.M();

        if(VBPartonDM_==13)hZMuMuMass_->Fill(VBPartonM_,mcWeight_);
        if(VBPartonDM_==11)hZElElMass_->Fill(VBPartonM_,mcWeight_);
      }
    }
 
      // ---- consider only final state particles -----------------------
      if (i_gen->status() == ((mGENType==2)?3:1)) { 
        // ---- consider only electron and muon flavors -----------------
        if (abs(i_gen->pdgId()) == 11 || abs(i_gen->pdgId()) == 13) {
            GENPARTICLE aGenLepton;
            TLorentzVector lepP4GEN(i_gen->p4().Px(),i_gen->p4().Py(),i_gen->p4().Pz(),i_gen->p4().E());
            TLorentzVector lepP4BARE(i_gen->p4().Px(),i_gen->p4().Py(),i_gen->p4().Pz(),i_gen->p4().E()); //same as before but this will not be modified
		//look for the GEN type: 0 GEN - 1 JEF  - 2 PAR
		switch (mGENType)
		{
		case 1:
			//--- search for status 1 photons in DR=0.1 cone around the bare lepton
			for(j_gen=gen->begin();j_gen !=gen->end();j_gen++){
			  if(j_gen->status() == 1 && abs(j_gen->pdgId()) ==22){
			    TLorentzVector phoP4GEN(j_gen->p4().Px(),j_gen->p4().Py(),j_gen->p4().Pz(),j_gen->p4().E());
			    float DR = phoP4GEN.DeltaR(lepP4BARE);
			    if(DR<mDressedRadius) lepP4GEN+=phoP4GEN;
			  }//if: st=1, pdgid=22
			}//for: gen particle (2nd loop)
			break;
		case 2: 
		case 0:
		default: break;
		}
            aGenLepton.p4    = lepP4GEN; 
            aGenLepton.pdgId = i_gen->pdgId();
           // ---- apply geometric and kinematic acceptance -------------
          if ((i_gen->pt() > mMinLepPt) && (fabs(i_gen->eta())) < mMaxLepEta) { // N.B. It is moved after JEF PAR GEN sel
            myGenLeptons.push_back(aGenLepton);
          }// geometrica and kinematics acc
        }
	
        // ---- consider only photons -----------------
        if (abs(i_gen->pdgId()) == 22) {
           // ---- apply geometric and kinematic acceptance -------------
          if ((i_gen->pt() > mMinPhoPtId) && (fabs(i_gen->eta())) < mMaxPhoEta) {
            GENPARTICLEPHO aGenPhoton;
            TLorentzVector phoP4GEN(i_gen->p4().Px(),i_gen->p4().Py(),i_gen->p4().Pz(),i_gen->p4().E());
	    //cout<<"evt "<<iEvent.id().event()<<" ph pt/eta "<<i_gen->pt() <<"/"<<i_gen->eta()<<endl;
	    //--- search for stable particles around photon in DR=0.3,0.4 and 0.5 cone around the photon
	    //--- for particle isolation, don't check for invisible particles
	    //correct for the fact that also the photon itself IS among the list of particles looped over
	    //subtract the momentum vector and the energy
	    float isoPxDR03=-i_gen->p4().Px();
	    float isoPxDR04=-i_gen->p4().Px();
	    float isoPxDR05=-i_gen->p4().Px();
	    float isoPyDR03=-i_gen->p4().Py();
	    float isoPyDR04=-i_gen->p4().Py();
	    float isoPyDR05=-i_gen->p4().Py();
	    float isoEDR03=-i_gen->p4().E();
	    float isoEDR04=-i_gen->p4().E();
	    float isoEDR05=-i_gen->p4().E();
	    for(j_gen=gen->begin();j_gen !=gen->end();j_gen++){
	      //don't use neutrinos
	      if(j_gen->status() == 1 && abs(j_gen->pdgId()) !=12 && abs(j_gen->pdgId()) !=14 &&  abs(j_gen->pdgId()) !=16){
		TLorentzVector partP4GEN(j_gen->p4().Px(),j_gen->p4().Py(),j_gen->p4().Pz(),j_gen->p4().E());
		float DR = partP4GEN.DeltaR(phoP4GEN);
		if(DR<0.5){
		  isoPxDR05+=j_gen->p4().Px();
		  isoPyDR05+=j_gen->p4().Py();
		  isoEDR05+=j_gen->p4().E();
		  if(DR<0.4){
		    isoPxDR04+=j_gen->p4().Px();
		    isoPyDR04+=j_gen->p4().Py();
		    isoEDR04+=j_gen->p4().E();
		    if(DR<0.3){
		      isoPxDR03+=j_gen->p4().Px();
		      isoPyDR03+=j_gen->p4().Py();
		      isoEDR03+=j_gen->p4().E();
		    }
		  }
		}
	      }
	    }

            aGenPhoton.pdgId = i_gen->pdgId();
            aGenPhoton.p4    = phoP4GEN;
	    aGenPhoton.isoEDR03=isoEDR03;
	    aGenPhoton.isoEDR04=isoEDR04;
	    aGenPhoton.isoEDR05=isoEDR05;
	    aGenPhoton.isoPtDR03=sqrt(isoPxDR03*isoPxDR03+isoPyDR03*isoPyDR03);
	    aGenPhoton.isoPtDR04=sqrt(isoPxDR04*isoPxDR04+isoPyDR04*isoPyDR04);
	    aGenPhoton.isoPtDR05=sqrt(isoPxDR05*isoPxDR05+isoPyDR05*isoPyDR05);
	    // find mother id
	    for(int kk = 0; kk < int(i_gen-> numberOfMothers()); ++kk) {
	      const GenParticle * gen_Moth = static_cast<const GenParticle*>(i_gen->mother(kk)); // find mother
	      aGenPhoton.motherId = gen_Moth->pdgId();
	    }
	    if(aGenPhoton.isoPtDR03/phoP4GEN.Pt()<0.25){
	      myGenPhotons.push_back(aGenPhoton);
	    }
          }//eta and pt cutoff
        }
      }//gen status IF
	//HTParSum
      if(i_gen->status()== 3){
		switch(abs(i_gen->pdgId())){
		case 1:case 2: case 3: case 4: case 5:case 6: case 21: //quark udscbt and gluons?
			HTParSum_+=i_gen->pt();
			break;
		default: break;
		}
		}
    }
    hGenLeptons_->Fill(int(myGenLeptons.size()));
    // ---- sort the genLeptons -----------------------------------------
    sort(myGenLeptons.begin(),myGenLeptons.end(),lepSortingRuleGEN);
    // ---- sort the genPhotons -----------------------------------------
    sort(myGenPhotons.begin(),myGenPhotons.end(),phoSortingRuleGEN);
    // ---- genjets -----------------------------------------------------
    for(i_genjet = genjets->begin(); i_genjet != genjets->end(); i_genjet++) {
      // ---- genlepton - genjet cross cleaning -------------------------
      //bool isISO(true);
      int isVETO=0;
      //if(mGENCrossCleaning&1)
      {
      for(unsigned l=0;l<myGenLeptons.size();l++) { 
        // ---- genjet vs 2 leading genlepton cleaning ------------------
        if (l >= 2) continue; 
        if (deltaR(i_genjet->eta(),i_genjet->phi(),myGenLeptons[l].p4.Eta(),myGenLeptons[l].p4.Phi()) < mJetLepIsoR) {
          isVETO |= (1<<l);
	  if (i_genjet->pt()>=mMinJetPt && (fabs(i_genjet->eta()) <= mMaxJetEta)){
	    //cout<<"DR is in GEN, lep "<<deltaR(i_genjet->eta(),i_genjet->phi(),myGenLeptons[l].p4.Eta(),myGenLeptons[l].p4.Phi())<<"/"<<l<<" pt gen/lep "<<i_genjet->pt()<<"/"<<myGenLeptons[l].p4.Pt()<<" ID "<<myGenLeptons[l].pdgId<<"/"<<i_genjet->energy()<<"/"<<i_genjet->emEnergy()<<"/"<<i_genjet->hadEnergy()<<endl;
	    //if(abs(myGenLeptons[l].pdgId)==13){
	    //for(unsigned int i=0;i<i_genjet->getGenConstituents().size();i++){
	    //cout<<"particle "<<i<<" ID/pt "<<i_genjet->getGenConstituents()[i]->pdgId()<<"/"<<i_genjet->getGenConstituents()[i]->pt()<<endl;
	    //}
	    //}
	    nJets_lepVetoGEN+=1;
	  }
          continue;
        }
      }
      }
      //if(mGENCrossCleaning&2)
      {
	// ---- genjet vs leading leading photon cleaning ------------------
	//if(myGenPhotons.size()>0){
	for(unsigned l=0;l<myGenPhotons.size();l++) { 
	  if (deltaR(i_genjet->eta(),i_genjet->phi(),myGenPhotons[0].p4.Eta(),myGenPhotons[0].p4.Phi()) < mJetPhoIsoR) { 
	    isVETO|= (1<<2);
	    if (i_genjet->pt()>=mMinJetPt && (fabs(i_genjet->eta()) <= mMaxJetEta)){
	      nJets_phoVetoGEN+=1;
	    }
	    continue;
	  }
	}
      }
      //if (!isISO) continue;
      // ---- preselection on genjets -----------------------------------
      if ((i_genjet->pt() < mMinJetPt) || (fabs(i_genjet->eta()) > mMaxJetEta)) continue;
      //for(unsigned int i=0;i<i_genjet->getGenConstituents().size();i++){
	    //cout<<"particle "<<i<<" ID/pt "<<i_genjet->getGenConstituents()[i]->pdgId()<<"/"<<i_genjet->getGenConstituents()[i]->pt()<<endl;
	    //}
      //}
      GENJET aGenJet(i_genjet->p4().Px(),i_genjet->p4().Py(),i_genjet->p4().Pz(),i_genjet->p4().E());
      aGenJet.veto=isVETO;  
      //---- try first B hadron code

      //---- gen jet matching
      float DR_parton=0.5;
      float pt_parton=0;
      for(i_gen = gen->begin(); i_gen != gen->end(); i_gen++) {
	if(i_gen->status()==3){
	  TLorentzVector i_gen_lv(i_gen->p4().Px(),i_gen->p4().Py(),i_gen->p4().Pz(),i_gen->p4().E());
	  //take out the incoming protons before the collision - otherwise code crashes
	  if(i_gen_lv.Pt()!=0){
	    if(aGenJet.DeltaR(i_gen_lv) <DR_parton) 
	      {
		if(pt_parton< i_gen->pt()) //keep the leading object inside the cone
		  {aGenJet.id=i_gen->pdgId(); 
		    pt_parton=i_gen->pt();
		    DR_parton=aGenJet.DeltaR(i_gen_lv) ;
		  }
		aGenJet.nparton++; //count how many partons are inside the cone
	      }
	  }
	}
      }
      myGenJets.push_back(aGenJet);  
    }// genjet loop
  }// if MC

  //define lepton, photon and jet member to continue the MC only loop
  //with less effort
  vector<PARTICLE> myLeptons;
  vector<PARTICLE> myPhotons;
  //vector<PARTICLE> myFSRphotons;
  vector<JET> myJets;
  vector<JET> myFwJets;
  //vector<JET> myRJets;

  TLorentzVector mypfmetP4(0,0,0,0);

  if(!mOnlyMC){
    //---- Rho ------------------------------------------------------------
    Handle<double> rho;
    iEvent.getByLabel(mSrcRho,rho);
    //---- Rho25 ------------------------------------------------------------
    Handle<double> rho25;
    iEvent.getByLabel(mSrcRho25,rho25);
    rho_        = *rho; 
    rho25_      = *rho25;

    //---- reco vertices block --------------------------------------------
    edm::Handle<VertexCollection> vertices_;
    iEvent.getByLabel("offlinePrimaryVertices", vertices_);
    const reco::Vertex *primVtx = &(*(vertices_.product()))[0];
    for(VertexCollection::const_iterator i_vtx = vertices_->begin(); i_vtx != vertices_->end(); ++i_vtx) {  
      if (!i_vtx->isFake() && (fabs(i_vtx->z()) < 24) && (i_vtx->ndof() >= 4)) {
	vtxZ_   ->push_back(i_vtx->z());
	vtxNdof_->push_back(i_vtx->ndof());
      }
    }
    //---- beam spot for conversion-safe electron veto --------------------
    Handle<reco::BeamSpot> bsHandle;
    iEvent.getByLabel("offlineBeamSpot", bsHandle);
    const reco::BeamSpot &beamspot = *bsHandle.product();
    
    //---- conversions for conversion-safe electron veto -------------------
    Handle<reco::ConversionCollection> hConversions;
    iEvent.getByLabel("allConversions", hConversions);
    
    //---- leptons block --------------------------------------------------
    //  Handle<View<GsfElectron> > electrons_;
    Handle<GsfElectronCollection> electrons_;
    iEvent.getByLabel("gsfElectrons", electrons_);
    GsfElectronRef::key_type eleIndex = 0;
    // Fetch pf iso maps
    Handle<ValueMap<double> > pfIsoValEleCH03;
    iEvent.getByLabel(pfIsoValEleCH03Name, pfIsoValEleCH03);
    Handle<ValueMap<double> > pfIsoValEleNH03;
    iEvent.getByLabel(pfIsoValEleNH03Name, pfIsoValEleNH03);
    Handle<ValueMap<double> > pfIsoValEleG03;
    iEvent.getByLabel(pfIsoValEleG03Name, pfIsoValEleG03);

    Handle<View<Muon> > muons_;
    iEvent.getByLabel("muons",muons_);
    
    // chiara
    //----- preshower rechits block for gamma id --------------------------------------------- 
    Handle<ESRecHitCollection> ecalhitses_;
    iEvent.getByLabel("reducedEcalRecHitsES", ecalhitses_);
    std::map<DetId,EcalRecHit> rechits_map_;
    if (ecalhitses_.isValid()) {
      EcalRecHitCollection::const_iterator it;
      for (it = ecalhitses_->begin(); it != ecalhitses_->end(); ++it) {
	// remove bad ES rechits
	if (it->recoFlag()==1 || it->recoFlag()==14 || (it->recoFlag()<=10 && it->recoFlag()>=5)) continue;
	//Make the map of DetID, EcalRecHit pairs
	rechits_map_.insert(std::make_pair(it->id(), *it));
      }
    }
    
    // chiara
    // geometry and topology for gamma id
    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloGeometry* geometry = geoHandle.product();

    // chiara
    edm::ESHandle<CaloTopology> pTopology;
    iSetup.get<CaloTopologyRecord>().get(pTopology);
    //    const CaloTopology *topology = pTopology.product();
    CaloSubdetectorTopology *topology_p = new EcalPreshowerTopology(geoHandle);

    //---- tracks block -----------------------
    Handle<TrackCollection> tracks_;
    iEvent.getByLabel("generalTracks",tracks_); 
    
    // ---- loop over muons -----------------------------------------------
    for(View<Muon>::const_iterator i_mu = muons_->begin(); i_mu != muons_->end(); ++i_mu) {
      //---- apply kinematic and geometric acceptance 
      if ((i_mu->pt() < mMinLepPt) || (fabs(i_mu->eta()) > mMaxLepEta))  continue;
      
      //---- apply tight definition for muon selection
      //---- https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
      if (!muon::isTightMuon(*i_mu,*primVtx))          continue;
      
      //---- pfi isolation not corrected, for efficiency studies look at
      //---- https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs
      //---- https://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=128982 
      float muonIsoPFUnc = (i_mu->pfIsolationR04().sumChargedHadronPt + i_mu->pfIsolationR04().sumNeutralHadronEt + i_mu->pfIsolationR04().sumPhotonEt)/i_mu->pt();
      
      // Muon pf isolation rho corrected
      float muonIsoPFRho = (i_mu->pfIsolationR04().sumChargedHadronPt
			  + max(0.,i_mu->pfIsolationR04().sumNeutralHadronEt + i_mu->pfIsolationR04().sumPhotonEt - (*rho25)*getEffectiveAreaForMuons(i_mu->eta())))/i_mu->pt();
      
      // Muon pf isolation DB corrected
      float muonIsoPFdb  = (i_mu->pfIsolationR04().sumChargedHadronPt 
			    + max(0., i_mu->pfIsolationR04().sumNeutralHadronEt + i_mu->pfIsolationR04().sumPhotonEt - 0.5*i_mu->pfIsolationR04().sumPUPt))/i_mu->pt();    


      if (muonIsoPFdb > mMaxCombRelIso04)                                 continue;

      pat::Muon PAT_muonF1, PAT_muonF3, PAT_muonF5;
      Handle<View<pat::Muon> > PATmuonsF1_, PATmuonsF3_, PATmuonsF5_;
      iEvent.getByLabel("patMuonsTriggerMatchHLTFamily1",PATmuonsF1_);
      iEvent.getByLabel("patMuonElectronTriggerMatchHLTFamily3mu",PATmuonsF3_);
      iEvent.getByLabel("patMuonsTriggerMatchHLTFamily5",PATmuonsF5_);
      for(View<pat::Muon>::const_iterator PATi_muF1 = PATmuonsF1_->begin(); PATi_muF1 != PATmuonsF1_->end(); ++PATi_muF1){
      	edm::Ptr<reco::Candidate> originalRef_muF1 = PATi_muF1->originalObjectRef();
      	const reco::Muon *originalMuF1 = dynamic_cast<const reco::Muon *>(originalRef_muF1.get());
    		if(originalMuF1->pt() == i_mu->pt()) {
	  			PAT_muonF1 = *PATi_muF1; 
	  			break;
				}
      }
      for(View<pat::Muon>::const_iterator PATi_muF3 = PATmuonsF3_->begin(); PATi_muF3 != PATmuonsF3_->end(); ++PATi_muF3){
      	edm::Ptr<reco::Candidate> originalRef_muF3 = PATi_muF3->originalObjectRef();
      	const reco::Muon *originalMuF3 = dynamic_cast<const reco::Muon *>(originalRef_muF3.get());
    		if(originalMuF3->pt() == i_mu->pt()) {
	  			PAT_muonF3 = *PATi_muF3; 
	  			break;
				}
      }
      for(View<pat::Muon>::const_iterator PATi_muF5 = PATmuonsF5_->begin(); PATi_muF5 != PATmuonsF5_->end(); ++PATi_muF5){
      	edm::Ptr<reco::Candidate> originalRef_muF5 = PATi_muF5->originalObjectRef();
      	const reco::Muon *originalMuF5 = dynamic_cast<const reco::Muon *>(originalRef_muF5.get());
    		if(originalMuF5->pt() == i_mu->pt()) {
	  			PAT_muonF5 = *PATi_muF5; 
	  			break;
				}
      }

      
      PARTICLE aLepton;
      TLorentzVector lepP4(i_mu->p4().Px(),i_mu->p4().Py(),i_mu->p4().Pz(),i_mu->p4().E());
      aLepton.p4   = lepP4;
      aLepton.chid = 2*i_mu->charge();
      aLepton.id   = 1; // all muons are tight
      aLepton.isoPFUnc = muonIsoPFUnc;
      aLepton.isoPFDb  = muonIsoPFdb;
      aLepton.isoPFRho = muonIsoPFRho;
      aLepton.sigmaIEtaIEta  = -1;
      aLepton.hadronicOverEm = -1;
      aLepton.triObjMatchF1   = PAT_muonF1.triggerObjectMatches();
      aLepton.triObjMatchF3mu = PAT_muonF3.triggerObjectMatches();
      aLepton.triObjMatchF5   = PAT_muonF5.triggerObjectMatches();
      aLepton.r9_or_chi2ndof  = i_mu->globalTrack()->normalizedChi2();
      myLeptons.push_back(aLepton);
    }// muon loop
    // ---- loop over electrons -------------------------------------------
    //  for(View<GsfElectron>::const_iterator i_el = electrons_->begin();i_el != electrons_->end(); ++i_el) {
    for(GsfElectronCollection::const_iterator i_el = electrons_->begin();i_el != electrons_->end(); ++i_el) {
      float elPt                           = i_el->p4().Pt();
      float elEta                          = i_el->p4().Eta();
      float sigmaIetaIeta                  = i_el->sigmaIetaIeta();
      float hadronicOverEm                 = i_el->hadronicOverEm();
      float deltaPhiSuperClusterTrackAtVtx = i_el->deltaPhiSuperClusterTrackAtVtx();
      float deltaEtaSuperClusterTrackAtVtx = i_el->deltaEtaSuperClusterTrackAtVtx();
      float r9=i_el->r9();
      bool  isMedium(false);

      float epDifference 			 = fabs( 1./i_el->ecalEnergy() - i_el->eSuperClusterOverP()/i_el->ecalEnergy()  );
      
      if ((elPt < mMinLepPt) || (fabs(elEta) > mMaxLepEta)) continue;
      // ---- use WP90 as default preselection, store also WP80 subset (https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification)
      if (i_el->isEB()) {
	if (sigmaIetaIeta                  > 0.01)          continue; 
	if (deltaPhiSuperClusterTrackAtVtx > 0.15)          continue; //2011: 0.8
	if (deltaEtaSuperClusterTrackAtVtx > 0.007)         continue;
	if (hadronicOverEm                 > 0.12)          continue; //2011: 0.15 	
	if (epDifference                   > 0.05)	    continue; //2011: Not implemented
	if (sigmaIetaIeta < 0.01) {  // WP80 subset: 2012 Medium
	  if (deltaPhiSuperClusterTrackAtVtx    < 0.06 
	      && deltaEtaSuperClusterTrackAtVtx < 0.004
	      && hadronicOverEm                 < 0.12)
	    isMedium = true;
	}
      }// if EE
      if (i_el->isEE()) {
	if (sigmaIetaIeta                  > 0.03)          continue;
	if (deltaPhiSuperClusterTrackAtVtx > 0.1)           continue; //2011:.7
	if (deltaEtaSuperClusterTrackAtVtx > 0.009)         continue;
	if (hadronicOverEm                 > 0.10)          continue; //2011:.15
	if ( epDifference                  > 0.05)          continue; //2011: Not implemented
	if (sigmaIetaIeta<0.03) {  // WP80 subset
	  if (deltaPhiSuperClusterTrackAtVtx    < 0.03
	      && deltaEtaSuperClusterTrackAtVtx < 0.007
	      && hadronicOverEm                 < 0.10) 
	    isMedium = true;
	} 
      }
      if(ConversionTools::hasMatchedConversion(*i_el, hConversions, beamspot.position())) continue;
      if(i_el->gsfTrack()->trackerExpectedHitsInner().numberOfHits() > 1)                 continue;


      GsfElectronRef electronRef(electrons_, eleIndex++);
      float sumChargedHadronPt = (*pfIsoValEleCH03)[electronRef];
      float sumNeutralHadronEt = (*pfIsoValEleNH03)[electronRef];
      float sumPhotonEt        = (*pfIsoValEleG03)[electronRef];

      float electronIsoPFUnc   = (sumChargedHadronPt + sumNeutralHadronEt + sumPhotonEt)/i_el->pt();
      float electronIsoPFRho   = (sumChargedHadronPt + std::max(sumNeutralHadronEt + sumPhotonEt -  (*rho) * getEffectiveAreaForElectrons(i_el->eta()), 0.))/i_el->pt();     

      if(electronIsoPFRho > mMaxCombRelIso03)              continue;

      pat::Electron PAT_electronF2, PAT_electronF3, PAT_electronF6;
      Handle<View<pat::Electron> > PATelectronsF2_, PATelectronsF3_, PATelectronsF6_;
      iEvent.getByLabel("patElectronsTriggerMatchHLTFamily2",PATelectronsF2_);
      iEvent.getByLabel("patMuonElectronTriggerMatchHLTFamily3ele",PATelectronsF3_);
      iEvent.getByLabel("patElectronsTriggerMatchHLTFamily6",PATelectronsF6_);
      for(View<pat::Electron>::const_iterator PATi_elF2 = PATelectronsF2_->begin(); PATi_elF2 != PATelectronsF2_->end(); ++PATi_elF2){
      	edm::Ptr<reco::Candidate> originalRef_elF2 = PATi_elF2->originalObjectRef();
      	const reco::GsfElectron *originalElF2 = dynamic_cast<const reco::GsfElectron *>(originalRef_elF2.get());
    		if(originalElF2->pt() == i_el->pt()){
          PAT_electronF2 = *PATi_elF2;
          break;
        }
      }
      for(View<pat::Electron>::const_iterator PATi_elF3 = PATelectronsF3_->begin(); PATi_elF3 != PATelectronsF3_->end(); ++PATi_elF3){
      	edm::Ptr<reco::Candidate> originalRef_elF3 = PATi_elF3->originalObjectRef();
      	const reco::GsfElectron *originalElF3 = dynamic_cast<const reco::GsfElectron *>(originalRef_elF3.get());
    		if(originalElF3->pt() == i_el->pt()){
          PAT_electronF3 = *PATi_elF3;
          break;
        }
      }
      for(View<pat::Electron>::const_iterator PATi_elF6 = PATelectronsF6_->begin(); PATi_elF6 != PATelectronsF6_->end(); ++PATi_elF6){
      	edm::Ptr<reco::Candidate> originalRef_elF6 = PATi_elF6->originalObjectRef();
      	const reco::GsfElectron *originalElF6 = dynamic_cast<const reco::GsfElectron *>(originalRef_elF6.get());
    		if(originalElF6->pt() == i_el->pt()){
          PAT_electronF6 = *PATi_elF6;
          break;
        }
      }
      

      PARTICLE aLepton;
      TLorentzVector lepP4(i_el->p4().Px(),i_el->p4().Py(),i_el->p4().Pz(),i_el->p4().E());
      aLepton.p4   = lepP4;
      aLepton.chid = i_el->charge();
      aLepton.id   = 0;
      if (isMedium) {
	aLepton.id = 1;
      }

      aLepton.isoPFUnc = electronIsoPFUnc;
      aLepton.isoPFDb  = -1;
      aLepton.isoPFRho = electronIsoPFRho;
      aLepton.r9_or_chi2ndof=r9;
      aLepton.sigmaIEtaIEta = sigmaIetaIeta;
      aLepton.hadronicOverEm = hadronicOverEm;
      aLepton.triObjMatchF2 = PAT_electronF2.triggerObjectMatches();
      aLepton.triObjMatchF3ele = PAT_electronF3.triggerObjectMatches();
      aLepton.triObjMatchF6 = PAT_electronF6.triggerObjectMatches();
      myLeptons.push_back(aLepton);
    } // electrons loop
    hRecoLeptons_->Fill(int(myLeptons.size()));
    // ---- sort the leptons according to their pt ------------------------
    sort(myLeptons.begin(),myLeptons.end(),lepSortingRule); 
  

   Handle<PFCandidateCollection>  pfCandidates;
   iEvent.getByLabel("particleFlow", pfCandidates);
  
    //---- photons block --------------------------------------------------

    if(true) 
      {
	Handle<reco::PhotonCollection> photons_;
	iEvent.getByLabel("photons",photons_);
	
	std::vector<std::string> photonIDCollectionTags_;
	photonIDCollectionTags_.push_back("PhotonCutBasedIDLoose");
	photonIDCollectionTags_.push_back("PhotonCutBasedIDLooseEM");
	photonIDCollectionTags_.push_back("PhotonCutBasedIDTight");
	
	const int nPhoIDC = photonIDCollectionTags_.size();
	std::vector< const edm::ValueMap<Bool_t>* > phoIds;
	
	for(int j=0; j<nPhoIDC; j++) {
	  edm::Handle<edm::ValueMap<Bool_t> > phoIDValueMap_;
	  iEvent.getByLabel("PhotonIDProd",photonIDCollectionTags_[j], phoIDValueMap_);
	  phoIds.push_back( phoIDValueMap_.product() );
	}
	
	//---- loop over the photon collection ---------------------------------
	
	// pass the collection to the ID calculator // chiara
	cicPhotonId->configure(vertices_, tracks_, electrons_, pfCandidates, rho_);
	EcalClusterLazyTools lazyTools(iEvent, iSetup, edm::InputTag("reducedEcalRecHitsEB"),
				       edm::InputTag("reducedEcalRecHitsEE"));
	int ipho = 0;
	for(reco::PhotonCollection::const_iterator it = photons_->begin();it != photons_->end(); it++) {
	  
	  //for low pt photons we check isolation in order to reduce rate of saved events
	  bool passPhotCaloId = false;

	  //---- don't bother if it has pt less than 15 GeV ----------------------
	  if(it->pt() < 15) continue;
	  if(abs(it->eta()) > mMaxPhoEta) continue;
	  if(it->isEBEEGap())continue;
	  if(ConversionTools::hasMatchedPromptElectron(it->superCluster(), electrons_, hConversions, beamspot.position())) continue;

//  	  if(it->isEB() && it->hadronicOverEm()>0.15)continue; // on-line requirement 
//  	  if(it->isEB() && it->sigmaIetaIeta()>0.024)continue; // on-line requirement 
//  	  if(it->isEE() && it->hadronicOverEm()>0.10)continue; // on-line requirement 
//  	  if(it->isEE() && it->sigmaIetaIeta()>0.040)continue; // on-line requirement 
	  
	  reco::PhotonRef phoRef(photons_,ipho++);
	  int photonID=0;
	  
	  TLorentzVector aPhoton(it->p4().Px(),it->p4().Py(),it->p4().Pz(),it->p4().E());
	  //   photonBit |= (it->isEB()          << 0);
	  
	  std::map<TString,UChar_t> idPairs;
	  for(int k=0; k<nPhoIDC; k++) {
	    idPairs[ TString(photonIDCollectionTags_[k].c_str()) ] = (*phoIds[k])[phoRef];
	    if(photonIDCollectionTags_[k] == "PhotonCutBasedIDLoose")photonID |= (*phoIds[k])[phoRef] <<1;
	    if(photonIDCollectionTags_[k] == "PhotonCutBasedIDTight")photonID |= (*phoIds[k])[phoRef] <<2;
	  }// for id
	  
	  float hcalTowerSumEtConeDR03            = it->hcalTowerSumEtConeDR03(); // hcalTowerSumEtConeDR03
	  float ecalRecHitSumEtConeDR03           = it->ecalRecHitSumEtConeDR03(); // ecalRecHitSumEtConeDR03
	  float trkSumPtHollowConeDR03            = it->trkSumPtHollowConeDR03();
	  
	  
	  float hcalTowerSumEtConeDR04            = it->hcalTowerSumEtConeDR04(); // hcalTowerSumEtConeDR04
	  float ecalRecHitSumEtConeDR04           = it->ecalRecHitSumEtConeDR04(); // ecalRecHitSumEtConeDR04
	  float trkSumPtHollowConeDR04            = it->trkSumPtHollowConeDR04();
	  float nTrkSolidConeDR04                 = it->nTrkSolidConeDR04();
	  float nTrkHollowConeDR04                = it->nTrkHollowConeDR04();
	  float trkSumPtSolidConeDR04             = it->trkSumPtSolidConeDR04();
	  
	  
	  float sigmaIetaIeta                     = it->sigmaIetaIeta();
	  //	  float phoHasConvTrks                    = it->hasConversionTracks();
	  // 	  float r9                                = it->r9();
	  float hadronicOverEm                    = it->hadronicOverEm();
	  float hadronicOverEm2012                = -1; // to be computed later
	  float sigmaIphiIphi                     = -1; // to be computed later
	  float sigmaIetaIphi                     = -1; // to be computed later
	  
	  float gammaPt = aPhoton.Pt();
	  bool  isTriggerISO = false;
	  
	  
	  // --- calculate new H/E for 2012
	  std::vector<CaloTowerDetId> hcalTowersBehindClusters = hcalHelper->hcalTowersBehindClusters(*(it->superCluster()));
	  float hcalDepth1 = hcalHelper->hcalESumDepth1BehindClusters(hcalTowersBehindClusters);
	  float hcalDepth2 = hcalHelper->hcalESumDepth2BehindClusters(hcalTowersBehindClusters);
	  hadronicOverEm2012 = (hcalDepth1 + hcalDepth2)/it->superCluster()->energy();
	  //in principle in CMSSW_5_2_X and later we could use the direct accessor, not available in 4_X
	  //hadronicOverEm2012 = it->hadTowOverEm();
	  
	  
	  // --- get iphi-iphi
	  EcalClusterLazyTools lazyTool(iEvent, iSetup, InputTag("reducedEcalRecHitsEB"), InputTag("reducedEcalRecHitsEE"));
	  
	  // Next few lines get sigma-iphi-iphi
	  const reco::CaloClusterPtr  seed = it->superCluster()->seed();
	  if (seed.isAvailable()) {
	    // Get sigma-iphi-iphi
	    std::vector<float> vCov = lazyTool.covariances(*seed);
	    sigmaIphiIphi  = sqrt(vCov[2]);   // This is Sqrt(covPhiPhi)
	    sigmaIetaIphi  = sqrt(vCov[1]);   // This is Sqrt(covEtaPhi)
	    // sigmaIetaIeta  = sqrt(vCov[0]);   // This is Sqrt(covEtaEta)  // chiara: controlla se e' uguale a quello sopra
	  }
	  
	  // -- conversion-safe electron veto
	  float photonPassConv = !ConversionTools::hasMatchedPromptElectron(it->superCluster(), electrons_, hConversions, beamspot.position());

	  // --- PF isolation ---
	  unsigned int ivtx = 0;
	  VertexRef myVtxRef(vertices_, ivtx);
	 
	  //MP ISO
	  SuperClusterFootprintRemoval remover(iEvent,edm::ParameterSet(),iSetup); 
	  float photonIsoFPRCharged = remover.PFIsolation("charged",it->superCluster(),0);
	  float photonIsoFPRNeutral = remover.PFIsolation("neutral",it->superCluster());
	  float photonIsoFPRPhoton = remover.PFIsolation("photon",it->superCluster());

	  // chiara: ma va fatto cosi'? c'era gia'
	  isolator.fGetIsolation((&*it),pfCandidates.product(), myVtxRef, vertices_);
	  float photonPfIsoCH = isolator.getIsolationCharged();
	  float photonPfIsoNH = isolator.getIsolationNeutral();
	  float photonPfIsoP   = isolator.getIsolationPhoton();

	  // chiara PF isolation for MVA ID and CICs - init
	  float photonPfIsoP03Cic = cicPhotonId->pfEcalIso(phoRef, 0.3, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);  
	  float photonPfIsoN03Cic = cicPhotonId->pfHcalIso(phoRef, 0.3, 0.0, reco::PFCandidate::h0);	                          
	  float photonPfIsoP04Cic = cicPhotonId->pfEcalIso(phoRef, 0.4, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma); 
	  float photonPfIsoN04Cic = cicPhotonId->pfHcalIso(phoRef, 0.4, 0.0, reco::PFCandidate::h0);	                          

	  std::vector<float> vPhotonPfIsoCharged03ForCic;
	  std::vector<float> vPhotonPfIsoCharged04ForCic;
	  vPhotonPfIsoCharged03ForCic = cicPhotonId->pfTkIsoWithVertex(phoRef, 0.3, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h); 
	  vPhotonPfIsoCharged04ForCic = cicPhotonId->pfTkIsoWithVertex(phoRef, 0.4, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h); 
	  if (vPhotonPfIsoCharged03ForCic.size() != vertices_->size()) cout << " problem " << endl;
	  if (vPhotonPfIsoCharged04ForCic.size() != vertices_->size()) cout << " problem " << endl;
	  float photonPfIsoC03Vtx0 = vPhotonPfIsoCharged03ForCic[0];  
	  float photonPfIsoC04Vtx0 = vPhotonPfIsoCharged04ForCic[0];  

	  float photonPfIsoC03Bad=-999.;   
	  float photonPfIsoC04Bad=-999.;   
	  for(unsigned int ivtx=0; ivtx<vertices_->size(); ivtx++) {
	    photonPfIsoC03Bad=vPhotonPfIsoCharged03ForCic[ivtx]>photonPfIsoC03Bad?vPhotonPfIsoCharged03ForCic[ivtx]:photonPfIsoC03Bad;
	    photonPfIsoC04Bad=vPhotonPfIsoCharged04ForCic[ivtx]>photonPfIsoC04Bad?vPhotonPfIsoCharged04ForCic[ivtx]:photonPfIsoC04Bad;
	  }
  
	  // chiara: variables for gammaID MVA
	  float my_photonid_sieie_       = it->sigmaIetaIeta();                
	  float my_photonid_sieip_       = sigmaIetaIphi;                      
	  float my_photonid_etawidth_    = it->superCluster()->etaWidth();     
	  float my_photonid_phiwidth_    = it->superCluster()->phiWidth();     
	  float my_photonid_r9_          = it->e3x3()/it->superCluster()->rawEnergy();  // chiara: e' uguale a R9 calcolato sopra?
	  float lambdaMinus     = (sigmaIetaIeta + sigmaIphiIphi - sqrt(pow(sigmaIetaIeta - sigmaIphiIphi, 2) + 4*pow(sigmaIetaIphi, 2)));
	  float lambdaPlus      = (sigmaIetaIeta + sigmaIphiIphi + sqrt(pow(sigmaIetaIeta - sigmaIphiIphi, 2) + 4*pow(sigmaIetaIphi, 2)));
	  float my_photonid_lambdaRatio_ = lambdaMinus/lambdaPlus;   
	  float e4phot  = lazyTools.e2x2(*(it->superCluster()->seed()));  
	  float e25phot = lazyTools.e5x5(*(it->superCluster()->seed()));   // chiara: controlla che venga come dal metodo it->e5x5()
	  float my_photonid_s4ratio_ = e4phot/it->e5x5();                      
	  float my_photonid_sceta_   = it->superCluster()->position().eta();   
  
	  float pid_esXwidth = 0.;
	  float pid_esYwidth = 0.;
	  if (ecalhitses_.isValid() && (fabs(it->superCluster()->eta()) > 1.6 && fabs(it->superCluster()->eta()) < 3)) {
	    std::vector<float> phoESHits0 = getESHits(it->superCluster()->x(), it->superCluster()->y(), it->superCluster()->z(), rechits_map_, *geometry, topology_p, 0);
	    std::vector<float> phoESShape = getESShape(phoESHits0);
	    pid_esXwidth = phoESShape[0];
	    pid_esYwidth = phoESShape[1];
	  }

	  float rr2=pid_esXwidth*pid_esXwidth+pid_esYwidth*pid_esYwidth;
	  float my_photonid_ESEffSigmaRR = 0.0; 
	  if(rr2>0. && rr2<999999.) my_photonid_ESEffSigmaRR = sqrt(rr2);  
	  // chiara PF isolation for ID and CICs - end



	  // --- https://twiki.cern.ch/twiki/bin/viewauth/CMS/QCDPhotonsTrigger2011
	  // chiara: vale ancora? su questi tagliano
	  if(ecalRecHitSumEtConeDR03 < 6.0 + 0.012*gammaPt) // requirement of _IsoVL_ type photon triggers 
	    if(hcalTowerSumEtConeDR03  < 4.0 + 0.005*gammaPt)
	      if(trkSumPtHollowConeDR03  < 4.0 + 0.002*gammaPt)isTriggerISO=true;
	  
	  // --- https://twiki.cern.ch/twiki/bin/viewauth/CMS/Vgamma2011PhotonID
	  // --- recommended photon isolation + id in one step
	  bool isVgamma2011 = false;
	  float Rho25 = *rho25;
	  float ATr = 0.0167 ;if(it->isEE()) ATr = 0.032;
	  float AEc = 0.183  ;if(it->isEE()) AEc = 0.090;
	  float AHc = 0.062  ;if(it->isEE()) AHc = 0.180;
	  float sigmaIetaIetaMax = 0.011 ;if(it->isEE()) sigmaIetaIetaMax = 0.03;
	  if(trkSumPtHollowConeDR04  < 2.0 + 0.001*gammaPt + ATr*Rho25)
	    if(ecalRecHitSumEtConeDR04 < 4.2 + 0.006*gammaPt + AEc*Rho25) 
	      if(hcalTowerSumEtConeDR04  < 2.2 + 0.0025*gammaPt + AHc*Rho25)
		if(sigmaIetaIeta<sigmaIetaIetaMax)
		  if(!it->hasPixelSeed())
		    if(it->isEE() || (it->isEB() && sigmaIetaIeta > 0.001 && sigmaIphiIphi > 0.001)) // additional EB spike cleaning
		      if(hadronicOverEm<0.05) isVgamma2011 = true;
	  photonID |= isVgamma2011 << 3;
	  
	  // --- online isolation + Vgamma2011 id
	  if(isTriggerISO) 
	    if(sigmaIetaIeta<sigmaIetaIetaMax)
	      if(!it->hasPixelSeed())
		if(it->isEE() || (it->isEB() && sigmaIetaIeta > 0.001 && sigmaIphiIphi > 0.001)) // additional EB spike cleaning
		  if(hadronicOverEm<0.05) 
		    photonID |= 1 << 4;
	  
	  // --- Vgamma2011 photon id w/o isolation
	  if(sigmaIetaIeta<sigmaIetaIetaMax) 
	    if(!it->hasPixelSeed())
	      if(it->isEE() || (it->isEB() && sigmaIetaIeta > 0.001 && sigmaIphiIphi > 0.001)) // additional EB spike cleaning
		if(hadronicOverEm<0.05) 
		  photonID |= 1 << 5;
	  
	  
	  // photon near masked region
	  float gammaEta = aPhoton.Eta();
	  float gammaPhi = aPhoton.Phi();
	  bool mask_0  = ( gammaPhi>=-2.72 && gammaPhi<=-2.61 && gammaEta>=-1.33 && gammaEta<=-1.25 );
	  bool mask_1  = ( gammaPhi<=-3.05 && gammaEta<=-1.40 );
	  bool mask_2  = ( gammaPhi<=-3.05 && gammaEta>=1.04 && gammaEta<=1.15 );
	  bool mask_3  = ( gammaPhi>=-2.64 && gammaPhi<=-2.52 && gammaEta>=-0.20 && gammaEta<=-0.08 );
	  bool mask_4  = ( gammaPhi>=-2.64 && gammaPhi<=-2.52 && gammaEta>=1.38 && gammaEta<=1.46 );
	  bool mask_5  = ( gammaPhi>=-2.03 && gammaPhi<=-1.91 && gammaEta>=-0.47 && gammaEta<=-0.3 );
	  bool mask_6  = ( gammaPhi>=-1.25 && gammaPhi<=-1.12 && gammaEta>=-1.25 && gammaEta<=-1.13 );
	  bool mask_7  = ( gammaPhi>=-0.81 && gammaPhi<=-0.69 && gammaEta>=-0.82 && gammaEta<=-0.69 );
	  bool mask_8  = ( gammaPhi>=-0.55 && gammaPhi<=-0.42 && gammaEta>=1.21 && gammaEta<=1.33 );
	  bool mask_9  = ( gammaPhi>=-0.29 && gammaPhi<=-0.16  && gammaEta>=1.38 );
	  bool mask_10 = ( gammaPhi>=0.07 && gammaPhi<=0.19 && gammaEta>=-0.29 && gammaEta<=-0.16 );
	  bool mask_11 = ( gammaPhi>=0.15 && gammaPhi<=0.28 && gammaEta>=-0.37 && gammaEta<=-0.25 );
	  bool mask_12 = ( gammaPhi>=0.69 && gammaPhi<=0.79 && gammaEta>=-0.20 && gammaEta<=-0.07 );
	  bool mask_13 = ( gammaPhi>=0.86 && gammaPhi<=0.97 && gammaEta>=-0.11 && gammaEta<=0.02 );
	  bool mask_14 = ( gammaPhi>=0.60 && gammaPhi<=0.70 && gammaEta>=0.96 && gammaEta<=1.06 );
	  bool mask_15 = ( gammaPhi>=1.74 && gammaPhi<=1.84 && gammaEta>=0.08 && gammaEta<=0.19 );
	  bool mask_16 = ( gammaPhi>=1.65 && gammaPhi<=1.75 && gammaEta>=0.87 && gammaEta<=0.97 );
	  bool mask_17 = ( gammaPhi>=1.99 && gammaPhi<=2.10 && gammaEta>=-0.97 && gammaEta<=-0.87 );
	  bool mask_18 = ( gammaPhi>=2.95 && gammaPhi<=3.05 && gammaEta>=-0.98 && gammaEta<=-0.87 );
	  bool mask_19 = ( gammaPhi>=2.78 && gammaPhi<=2.89 && gammaEta>=0.86 && gammaEta<=0.98);
	  bool mask_20 = ( gammaPhi>=2.69 && gammaPhi<=2.81 && gammaEta>= 1.39);
	  
	  bool isMasked = mask_0 || mask_1 || mask_2 || mask_3 || mask_4 || mask_5 || mask_6 || mask_7 || mask_8 || mask_9 || mask_10 || mask_11;
	  isMasked      = isMasked || mask_12 || mask_13 || mask_14 || mask_15 || mask_16 || mask_17 || mask_18 || mask_19 || mask_20;      
	  
	  
	  
	  // chiara: da sistemare ?
	  int photonBit=0;
	  photonBit |= (it->isEB()          << 0);
	  photonBit |= (it->isEE()          << 1);
	  photonBit |= (it->isEBEtaGap()    << 2);
	  photonBit |= (it->isEBPhiGap()    << 3);
	  photonBit |= (it->isEERingGap()   << 4);
	  photonBit |= (it->isEEDeeGap()    << 5);
	  photonBit |= (it->isEBEEGap()     << 6);
	  photonBit |= (it->hasPixelSeed()  << 7);
	  photonBit |= (isTriggerISO        << 8);
	  photonBit |= (isMasked            << 9);

	  pat::Photon PAT_photon;
    Handle<View<pat::Photon> > PATphotons_;
    iEvent.getByLabel("patPhotonsTriggerMatchHLTFamily4",PATphotons_);
    for(View<pat::Photon>::const_iterator PATi_ph = PATphotons_->begin(); PATi_ph != PATphotons_->end(); ++PATi_ph){
    	edm::Ptr<reco::Candidate> originalRef_ph = PATi_ph->originalObjectRef();
    	const reco::Photon *originalPh = dynamic_cast<const reco::Photon *>(originalRef_ph.get());
    	if(originalPh->pt() == it->pt()){ 
    		PAT_photon = *PATi_ph;
    		break;
      }
    }
	  

	  if(it->isEB()){
	    //if H/E<0.15 and sigmaietaieta<0.024 -> ID is passed
	    if(sigmaIetaIeta<0.024 && hadronicOverEm<0.15){
	      passPhotCaloId =true;
	    }
	  }else if (it->isEE()){
	    //if H/E<0.10 and sigmaietaieta<0.040 -> ID is passed
	    if(sigmaIetaIeta<0.040 && hadronicOverEm<0.10){
	      passPhotCaloId = true;
	    }
	  }

	  // chiara
	  PARTICLE gamma; // define pseudo lepton half of the p4 of the real photon
	  gamma.p4    = aPhoton;          // chiara ok
	  gamma.chid  = 0;
	  gamma.id    = photonID;
// 	  gamma.iso   = isTriggerISO;    // this bool 0/1
// 	  gamma.isoPF = 0; 
	  gamma.bit   = photonBit;
	  // chiara
	  gamma.parameters.resize(PARTICLE::nPhotonParameters);
	  gamma.parameters[PARTICLE::passconv]  = photonPassConv;    
	  gamma.parameters[PARTICLE::pfIsoCH]   = photonPfIsoCH;
	  gamma.parameters[PARTICLE::pfIsoNH]  = photonPfIsoNH;
	  gamma.parameters[PARTICLE::pfIsoP]    = photonPfIsoP;
	  gamma.parameters[PARTICLE::pfCic03P]  = photonPfIsoP03Cic;
	  gamma.parameters[PARTICLE::pfCic03N]  = photonPfIsoN03Cic;
	  gamma.parameters[PARTICLE::pfCic04P]  = photonPfIsoP04Cic;
	  gamma.parameters[PARTICLE::pfCic04N]  = photonPfIsoN04Cic;
	  gamma.parameters[PARTICLE::pfCic03Cg] = photonPfIsoC03Vtx0;
	  gamma.parameters[PARTICLE::pfCic04Cg] = photonPfIsoC04Vtx0;
	  gamma.parameters[PARTICLE::pfCic03Cb] = photonPfIsoC03Bad;
	  gamma.parameters[PARTICLE::pfCic04Cb] = photonPfIsoC04Bad;
	  gamma.parameters[PARTICLE::sieie] = my_photonid_sieie_;
	  gamma.parameters[PARTICLE::sieip] = my_photonid_sieip_;
	  gamma.parameters[PARTICLE::etaw]  = my_photonid_etawidth_;
	  gamma.parameters[PARTICLE::phiw]  = my_photonid_phiwidth_;
	  gamma.parameters[PARTICLE::r9]    = my_photonid_r9_;
	  gamma.parameters[PARTICLE::lR]    = my_photonid_lambdaRatio_;
	  gamma.parameters[PARTICLE::s4]    = my_photonid_s4ratio_;
	  gamma.parameters[PARTICLE::e25]    = e25phot;
	  gamma.parameters[PARTICLE::sceta] = my_photonid_sceta_;
	  gamma.parameters[PARTICLE::ESEff] = my_photonid_ESEffSigmaRR;
	  gamma.parameters[PARTICLE::hcalTowerSumEtConeDR04]=hcalTowerSumEtConeDR04;
	  gamma.parameters[PARTICLE::ecalRecHitSumEtConeDR04]=ecalRecHitSumEtConeDR04;
	  gamma.parameters[PARTICLE::nTrkSolidConeDR04]=nTrkSolidConeDR04;
	  gamma.parameters[PARTICLE::trkSumPtSolidConeDR04]=trkSumPtSolidConeDR04;     
	  gamma.parameters[PARTICLE::nTrkHollowConeDR04]=nTrkHollowConeDR04;        
	  gamma.parameters[PARTICLE::trkSumPtHollowConeDR04]=trkSumPtHollowConeDR04;    
	  //	  gamma.parameters[PARTICLE::phoHasConvTrks]=phoHasConvTrks;            
	  gamma.parameters[PARTICLE::oldHadronicOverEm]=hadronicOverEm;             
	  gamma.parameters[PARTICLE::hadronicOverEm2012]=hadronicOverEm2012;         
    gamma.triObjMatchF4 = PAT_photon.triggerObjectMatches();
	
	  gamma.isoFPRCharged=photonIsoFPRCharged ;
	  gamma.isoFPRNeutral=photonIsoFPRNeutral ;
	  gamma.isoFPRPhoton=photonIsoFPRPhoton ;

	  // --- save FSR photon candidates and gamma+jets in seperate paths
	  
//	  myFSRphotons.push_back(gamma);                          // FSR photons are *not* used in the photon+jet DR cone rejection
	  
	  //---- consider single event interpretation, exclusive di-lepton/photon interpretation (myLeptons.size()==0)  SINGLEEVENT:
	  //if(it->pt() > mMinPhoPt && myLeptons.size()==0 && isTriggerISO) myPhotons.push_back(gamma);    //note: hard photons imply later a DR cone rejection wrt the jets
	  //if(it->pt() > mMinPhoPt && isTriggerISO) myPhotons.push_back(gamma);    //note: hard photons imply later a DR cone rejection wrt the jets
	  //check first for harder photon cutoff without isolation requirement
	  if(it->pt() > mMinPhoPt ){ myPhotons.push_back(gamma);  }  //note: hard photons imply later a DR cone rejection wrt the jets
	  else if(it->pt() > mMinPhoPtId && passPhotCaloId ){//if photon not hard enough - check in addition for ID (CaloBased)
	    myPhotons.push_back(gamma); 
	  }
	}
	hRecoPhotons_->Fill(int(myPhotons.size()));   // chiara
      }
    sort(myPhotons.begin(),myPhotons.end(),lepSortingRule);
    //   sort(myFSRphotons.begin(),myFSRphotons.end(),lepSortingRule);
    //---- jets block -----------------------------------------------------
    edm::Handle<edm::View<pat::Jet> > jets;
    iEvent.getByLabel(mJetsName,jets);
    edm::View<pat::Jet> pat_jets = *jets;
    
    //Handle<PFJetCollection> jets_;
    //iEvent.getByLabel(mJetsName,jets_);

   //---- get the jet energy corrector -----------------------------------
    //if(isRealData_)mJEC = JetCorrector::getJetCorrector(mJECserviceDATA,iSetup);
    //if(!isRealData_)mJEC = JetCorrector::getJetCorrector(mJECserviceMC,iSetup);
    //---- the JEC uncertainty only needs to be set-up once ---------------
    //if (!mIsJECuncSet) {
    //edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
    //iSetup.get<JetCorrectionsRecord>().get(mPayloadName,JetCorParColl); 
    //JetCorrectorParameters const& JetCorPar = (*JetCorParColl)["Uncertainty"];
    //mJECunc = new JetCorrectionUncertainty(JetCorPar);
    //mIsJECuncSet = true;
    //}
    // ---- jets loop -----------------------------------------------------
    for(edm::View<pat::Jet>::const_iterator i_jet = pat_jets.begin();i_jet != pat_jets.end(); ++i_jet) {
      int index= i_jet - pat_jets.begin();  //index QG
      TLorentzVector jetP4(i_jet->px(),i_jet->py(),i_jet->pz(),i_jet->energy());
      int jetIsDuplicate(0);
      bool jetIsInAcceptance(true);
      bool jetIsIDed(true);
      bool jetIsInFWAcceptance(true);
      
      //----- remove the leptons ------------------------------------------
      for(unsigned int i_lep = 0; i_lep < myLeptons.size(); i_lep++) {
	if(i_lep>=2)continue;
	float DR = myLeptons[i_lep].p4.DeltaR(jetP4);
	if (DR < mJetLepIsoR) {
	  jetIsDuplicate |= 1<<i_lep; 
	   //set veto counter on jets which would survive the final selection
	  if (jetP4.Pt()>=mMinJetPt && (fabs(i_jet->eta()) <= mMaxJetEta)){
	    nJets_lepVeto+=1;
	  }
	}
	
      }// lepton loop 
      
      // chiara: questo va cambiato e vanno rimossi tutti
      // paolo: Keeping the overlap between photon andn jets will be removed later after the photonID requirements are applied
      //----- remove the leading photon ------------------------------------ (reminder nPhotons>0 only IF nLeptons==0)
	     bool leadPhotIdedFound=false;
             for(unsigned int i_pho = 0; i_pho < myPhotons.size(); i_pho++) { 
	       //rejection only wrt the leading photon (ided) in the longrun
	       //so far save more photon candidates - check with all those
		if(leadPhotIdedFound)continue;
       		float DR = myPhotons[i_pho].p4.DeltaR(jetP4);
		//if(myPhotons[i_pho].id)leadPhotIdedFound=true;
       		if (DR < mJetPhoIsoR /*&& (myPhotons[i_pho].id & (1<<3) )*/) {
       		  jetIsDuplicate |= 1<<2;
		  //set veto counter on jets which would survive the final selection
		  if (jetP4.Pt()>mMinJetPt && (fabs(i_jet->eta()) < mMaxJetEta)){
		    nJets_phoVeto+=1;
		  }
		}       		
       	      }// photon loop
      
      // ---- get the jec and the uncertainty -----------------------------    
      //int index = i_jet - jets_->begin();
      //edm::RefToBase<reco::Jet> jetRef(edm::Ref<PFJetCollection>(jets_,index));
      //double jec = mJEC->correction(*i_jet,jetRef,iEvent,iSetup);
      // ---- only keep jets within the kinematic acceptance --------------
      if ((i_jet->pt() < 15) || (fabs(i_jet->eta()) > 5.0)) continue; // apply first a hardcode preselection
      if ((i_jet->pt() < mMinJetPt) || (fabs(i_jet->eta()) > mMaxJetEta)) jetIsInAcceptance = false;
      if ((i_jet->pt() < mMinJetPt) || (fabs(i_jet->eta()) < mMaxJetEta)) jetIsInFWAcceptance = false;
      if(jetIsInFWAcceptance && jetIsInAcceptance){
	cout<<"jet is FW and not FW, obviously a contradiction"<<endl;
      }
      //mJECunc->setJetEta(i_jet->eta());
      // ---- the unc is a function of the corrected pt -------------------
      //mJECunc->setJetPt(jec * i_jet->pt());
      //double unc = mJECunc->getUncertainty(true);
      float jec  = 1./i_jet->jecFactor(0);
      float unc  = i_jet->userFloat("jecUnc");
      float btag = i_jet->bDiscriminator("combinedSecondaryVertexBJetTags");
     
      // ---- Extra B-Tag info: vertex properties -----
      const reco::SecondaryVertexTagInfo *svTagInfo = i_jet->tagInfoSecondaryVertex(); //get's the SV info from the input Jet      
      float taginfoNvtx    = svTagInfo->nVertices();
      float taginfoNtracks = svTagInfo->nVertexTracks();
      float taginfoVtxMass = -999.;
      if(taginfoNvtx) 	//pick the first secondary vertex (i.e. the best one)
	taginfoVtxMass = svTagInfo->secondaryVertex(0).p4().M();
      
      // ---- Gen-Jet match ---------------------------
      // will be 0 in case of data
      int mcflavour = i_jet->partonFlavour();

      // ---------------------------------------------
      float beta = i_jet->userFloat("beta");
      // ---- keep only jets that pass the tight id -----------------------
      float chf = i_jet->chargedHadronEnergyFraction();
      float nhf = i_jet->neutralHadronEnergyFraction() + i_jet->HFHadronEnergyFraction();
      float phf = i_jet->photonEnergy()/(i_jet->jecFactor(0) * i_jet->energy());
      float elf = i_jet->electronEnergy()/(i_jet->jecFactor(0) * i_jet->energy());
      float muf = i_jet->muonEnergy()/(i_jet->jecFactor(0) * i_jet->energy());
      int chm   = i_jet->chargedHadronMultiplicity();
      int npr   = i_jet->chargedMultiplicity() + i_jet->neutralMultiplicity();
      bool id = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(i_jet->eta())<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && chf>0 && chm>0) || fabs(i_jet->eta())>2.4));
      float rms=TMath::Sqrt( TMath::Power(i_jet->userFloat("axis1"),2) + TMath::Power(i_jet->userFloat("axis2"),2) );
      if (!id) jetIsIDed = false;
      // ---- jet vertex association --------------------------------------
      // ---- get the vector of tracks ------------------------------------ 
      //reco::TrackRefVector vTrks(i_jet->getTrackRefs());
      //float sumTrkPt(0.0),sumTrkPtBeta(0.0),sumTrkPtBetaStar(0.0),beta(-1.0),betaStar(-1.0);
      // ---- loop over the tracks of the jet -----------------------------
      /*
	for(reco::TrackRefVector::const_iterator i_trk = vTrks.begin(); i_trk != vTrks.end(); i_trk++) {
	sumTrkPt += (*i_trk)->pt();
	// ---- loop over all vertices ------------------------------------
	for(unsigned i_vtx = 0;i_vtx < vertices_->size();i_vtx++) {
        // ---- loop over the tracks associated with the vertex ---------
        if ((*vertices_)[i_vtx].isFake() || (*vertices_)[i_vtx].ndof() < 4) continue; 
        for(reco::Vertex::trackRef_iterator i_vtxTrk = (*vertices_)[i_vtx].tracks_begin(); i_vtxTrk != (*vertices_)[i_vtx].tracks_end(); ++i_vtxTrk) {
	// ---- match the jet track to the track from the vertex ------
	reco::TrackRef trkRef(i_vtxTrk->castTo<reco::TrackRef>());
	// ---- check if the tracks match -----------------------------
	if (trkRef == (*i_trk)) {
	if (i_vtx == 0) {
	sumTrkPtBeta += (*i_trk)->pt();
	}
	else {
	sumTrkPtBetaStar += (*i_trk)->pt();
	}   
	continue;
	}
        }
	}// vertices loop 
	}// jet tracks loop
	if (sumTrkPt > 0) {
	beta     = sumTrkPtBeta/sumTrkPt;
	betaStar = sumTrkPtBetaStar/sumTrkPt;
	}
      */
      JET aJet; 
      aJet.p4       = jetP4;
      aJet.jec      = jec;
      aJet.unc      = unc;
      aJet.area     = i_jet->jetArea();
      aJet.chf      = chf;
      aJet.nhf      = nhf;
      aJet.phf      = phf;
      aJet.elf      = elf;
      aJet.muf      = muf;
      aJet.id       = 0;
      aJet.veto	    =jetIsDuplicate;
      if (id) {
	aJet.id     = 1;
      }
      aJet.beta     = beta;
      aJet.btag           = btag;
      aJet.taginfoNvtx    = taginfoNvtx;      
      aJet.taginfoNtracks = taginfoNtracks;
      aJet.taginfoVtxMass = taginfoVtxMass;
      aJet.mcflavour      = mcflavour;
      //if(jetIsDuplicate){aJet.p4       = jetP4; myRJets.push_back(aJet);}  // store the uncorrected jet (this is virtually the matched lepton in DR) 
    Handle<double> rhoQG;
    iEvent.getByLabel("kt6PFJetsIsoQG",rhoQG);
//
      if( jetIsInAcceptance && jetIsIDed){ //store QG Variable for myJets
	vector<float> *QGvars=ComputeQGVariables(i_jet,iEvent,index);
	{
	  int i;
	  for( i=0;i<int(QGvars->size());i++)
	    QGVars_->push_back(QGvars->at(i));
	  if(QGVars_->size()<10) QGVars_->push_back(rhoQG);i++;
	  for(;i<9;i++)
	    QGVars_->push_back(-1.0);
	  QGVars_->push_back(-99.);
	}
	aJet.qgl=QGvars->at(1);
	QGvars->clear();
	delete QGvars;
      }
      aJet.rms=rms;
    if( jetIsInAcceptance && jetIsIDed)myJets.push_back(aJet);
    if( jetIsInFWAcceptance && jetIsIDed)myFwJets.push_back(aJet);
    }// jet loop
      

    // ---- sort jets according to their corrected pt ---------------------
    //sort(myJets.begin(),myJets.end(),jetSortingRule);    
    //sort(myRJets.begin(),myRJets.end(),jetSortingRule);    
    // ---- MET block -----------------------------------------------------
    Handle<View<PFMET> > pfmetCol_;
    iEvent.getByLabel("pfMet", pfmetCol_);
    float pfmetPx = (pfmetCol_->front()).px();
    float pfmetPy = (pfmetCol_->front()).py();
    pfmet_      = (pfmetCol_->front()).pt();
    pfmetPhi_   = (pfmetCol_->front()).phi();
    pfSumEt_    = (pfmetCol_->front()).sumEt();
    mypfmetP4.SetPxPyPzE(pfmetPx,pfmetPy,0,sqrt(pfmetPx * pfmetPx + pfmetPy * pfmetPy));
  }
  // ---- define di-lepton pair (if exist)
  TLorentzVector llP4(0,0,0,0);
  if(myLeptons.size()>1) llP4 = myLeptons[0].p4 + myLeptons[1].p4;
  
  TLorentzVector llP4GEN(0,0,0,0); 

  if(myGenLeptons.size()>1) llP4GEN = myGenLeptons[0].p4 + myGenLeptons[1].p4;  
  // ---- counters ------------------------------------------------------
  nVtx_        = int(vtxZ_->size());
  nLeptons_    = int(myLeptons.size()); 
  nPhotons_    = int(myPhotons.size());
  nPhotonsGEN_ = int(myGenPhotons.size());
  nJets_       = int(myJets.size());
  fwnJets_       = int(myFwJets.size());
 // nRJets_      = int(myRJets.size());
  nLeptonsGEN_ = int(myGenLeptons.size()); 
  nJetsGEN_    = int(myGenJets.size()); 
  // ---- Gen To Reco Matching for leptons ------------------------------------------------------
  for(unsigned gg=0; gg < myGenLeptons.size(); gg++) {                               
    lepMatchedDRGEN_->push_back(999);
    lepMatchedGEN_->push_back(999);
    for(unsigned rr = 0; rr < myLeptons.size(); rr++) {
      float DR = myGenLeptons[gg].p4.DeltaR(myLeptons[rr].p4);
      bool isSameFlavor = (abs(myGenLeptons[gg].pdgId) == 11 && abs(myLeptons[rr].chid)==1);
      isSameFlavor = isSameFlavor || (abs(myGenLeptons[gg].pdgId) == 13 && abs(myLeptons[rr].chid)==2);
      if(DR<lepMatchedDRGEN_->at(gg) && isSameFlavor) {
	lepMatchedGEN_->at(gg) = rr;                              //store recolepton matched index
	lepMatchedDRGEN_->at(gg) = DR;                           //store DR with matched GEN lepton
      }
    }
  }
  // ---- plot some inclusive di-lepton spectra (prior jet requirement)
  if(nLeptons_>1) {
    int dileptonId = myLeptons[0].chid*myLeptons[1].chid;
    hLepLepMass_->Fill(llP4.M());
    if(dileptonId==-4)hMuMuMass_->Fill(llP4.M());
    if(dileptonId==-1)hElElMass_->Fill(llP4.M());
    if(dileptonId==-2)hElMuMass_->Fill(llP4.M());
    if(dileptonId==-1 && abs(myLeptons[0].p4.Eta())<1.4 && abs(myLeptons[1].p4.Eta())<1.4)hElElEBMass_->Fill(llP4.M());
    if(dileptonId==-4)hMuMuMassWeighted_->Fill(llP4.M(),mcWeight_);
    if(dileptonId==-1)hElElMassWeighted_->Fill(llP4.M(),mcWeight_);
    if(dileptonId==-2)hElMuMassWeighted_->Fill(llP4.M(),mcWeight_);
  }
  if(nPhotonsGEN_>0) {
    float pt = myGenPhotons[0].p4.Pt();
    float eta = myGenPhotons[0].p4.Eta();
    float DR=999;
    if(myPhotons.size()>0)DR = myPhotons[0].p4.DeltaR(myGenPhotons[0].p4);
    hGenPhotonPt_->Fill(pt);
    hGenPhotonEta_->Fill(eta);
    if(DR<0.2) {
      hGenPhotonMatchedPt_->Fill(pt);
      hGenPhotonMatchedEta_->Fill(eta);
    }  
  }
  
  for(int rr=0; rr< int(lepMatchedDRGEN_->size()); rr++) {
    if(rr>=2)continue; // do this only for first 2 leptons
    float pt = myGenLeptons[rr].p4.Pt();
    float eta = myGenLeptons[rr].p4.Eta();
    int pdgId =  myGenLeptons[rr].pdgId;
    float DR = lepMatchedDRGEN_->at(rr);
    if(abs(pdgId)==13) {
      hGenMuonPt_->Fill(pt);
      hGenMuonEta_->Fill(eta);
    }
    if(abs(pdgId)==13 && DR<0.2) {
      hGenMuonMatchedPt_->Fill(pt);
      hGenMuonMatchedEta_->Fill(eta);
    }
    if(abs(pdgId)==11) {
      hGenElectronPt_->Fill(pt);
      hGenElectronEta_->Fill(eta);
    }
    if(abs(pdgId)==11 && DR<0.2) {
      hGenElectronMatchedPt_->Fill(pt);
      hGenElectronMatchedEta_->Fill(eta);
    }
  }

  //we only want to veto on the leading photon - so maximally 1 jet will be removed:
  //check if for n photons we have n vetoed jets - then nJets_phoVeto=1
  //if we have less vetoed jets - put nJets_phoVeto=0 
  //as the "real" photon candidate might be the one which doesn't veto a jet
  
  if(nPhotons_< nJets_phoVeto){
    nJets_phoVeto=0;
  }else{
    nJets_phoVeto=1;
  }

  if(nJets_phoVetoGEN>1){
    nJets_phoVetoGEN=1;
  }

  // ---- keep only selected events -------------------------------------
  bool selectionRECO = false;
  if(!mOnlyMC){
    selectionRECO = ((nVtx_ > 0) && (nLeptons_ > 1) && ((nJets_-nJets_lepVeto)>=mMinNjets)/*(nJets_ >= mMinNjets)*/ && llP4.M()>mMinLLMass);
  }

  selectionRECO = selectionRECO || ((nVtx_ > 0) &&  nPhotons_>0 &&  ((nJets_ -nJets_phoVeto)>=mMinNjets)/*(nJets_ >= mMinNjets)*/); // add photon logic for RECO
  bool selection(selectionRECO);
  bool selectionGEN(false);
  if (!isRealData_) {
    selectionGEN = ((nLeptonsGEN_ > 1) && ((nJetsGEN_ - nJets_lepVetoGEN) >= mMinNjets) && llP4GEN.M()>mMinLLMass); 
    selectionGEN = (selectionGEN || ((nPhotonsGEN_ > 0) && (nJetsGEN_ - nJets_phoVetoGEN)>= mMinNjets));      // add photon logic for GEN 
    selection +=  selectionGEN;
  }
  if (selection) {
    eventNum_   = iEvent.id().event();
    runNum_     = iEvent.id().run();
    lumi_       = iEvent.luminosityBlock();
    isRealData_ = isRealData_; // just pass through 
    if (selectionRECO) {
      // ---- hadronic recoil vector --------------------------------------
      TLorentzVector pfmetP4=mypfmetP4;
      TLorentzVector pfhadP4(0,0,0,0);
      if(myLeptons.size()>1) { 
	llM_                   = llP4.M();
	llPt_                  = llP4.Pt();
	llPhi_                 = llP4.Phi();
	if(llPt_>0)llY_        = llP4.Rapidity();
	if(llPt_>0)llEta_      = llP4.Eta();
	llDPhi_     = fabs(myLeptons[0].p4.DeltaPhi(myLeptons[1].p4));
        pfhadP4     = -pfmetP4 - llP4;      
      }
      TLorentzVector lepP4(0,0,0,0); 
      for(unsigned l = 0; l < myLeptons.size(); l++) {
        lepP4 += myLeptons[l].p4;
        lepPt_          ->push_back(myLeptons[l].p4.Pt());
        lepEta_         ->push_back(myLeptons[l].p4.Eta());
        lepPhi_         ->push_back(myLeptons[l].p4.Phi());
        lepE_           ->push_back(myLeptons[l].p4.Energy());
        lepPFIsoUnc_    ->push_back(myLeptons[l].isoPFUnc);
        lepPFIsoDBCor_  ->push_back(myLeptons[l].isoPFDb);
        lepPFIsoRhoCor_ ->push_back(myLeptons[l].isoPFRho);
        lepId_          ->push_back(myLeptons[l].id);
        lepChId_        ->push_back(myLeptons[l].chid);
	lepR9orChi2ndof_->push_back(myLeptons[l].r9_or_chi2ndof);

        lepHadronicOverEm_     ->push_back(myLeptons[l].hadronicOverEm);
        lepSigmaIEtaIEta_      ->push_back(myLeptons[l].sigmaIEtaIEta);
      }      
      pfhadPt_    = pfhadP4.Pt();
      // ---- store photon variables ------------------------------------
      for(unsigned gg = 0; gg < myPhotons.size(); gg++) {
	// myTree_->Branch("mPhotonj1"        ,"vector<float>"        ,&mPhotonj1_);
	// myTree_->Branch("ptPhotonj1"       ,"vector<float>"        ,&ptPhotonj1_);
	// myTree_->Branch("jetPhotonDPhi"    ,"vector<float>"        ,&jetPhotonDPhi_);  // chiara: questi erano gia' vector... 
        photonE_   -> push_back(myPhotons[gg].p4.Energy());
        photonPt_  -> push_back(myPhotons[gg].p4.Pt());
        photonEta_ -> push_back(myPhotons[gg].p4.Eta());
        photonPhi_ -> push_back(myPhotons[gg].p4.Phi());
	photonPassConversionVeto_ -> push_back(myPhotons[gg].parameters[PARTICLE::passconv]);   
	photonPfIsoChargedHad_    -> push_back(myPhotons[gg].parameters[PARTICLE::pfIsoCH]);
	photonPfIsoNeutralHad_    -> push_back(myPhotons[gg].parameters[PARTICLE::pfIsoNH]);
	photonPfIsoPhoton_        -> push_back(myPhotons[gg].parameters[PARTICLE::pfIsoP]);
	photonPfIsoPhotons03ForCic_  -> push_back(myPhotons[gg].parameters[PARTICLE::pfCic03P]);
	photonPfIsoNeutrals03ForCic_ -> push_back(myPhotons[gg].parameters[PARTICLE::pfCic03N]);
	photonPfIsoPhotons04ForCic_  -> push_back(myPhotons[gg].parameters[PARTICLE::pfCic04P]);
	photonPfIsoNeutrals04ForCic_ -> push_back(myPhotons[gg].parameters[PARTICLE::pfCic04N]);
	photonPfIsoCharged03ForCicVtx0_ -> push_back(myPhotons[gg].parameters[PARTICLE::pfCic03Cg]);
	photonPfIsoCharged04ForCicVtx0_ -> push_back(myPhotons[gg].parameters[PARTICLE::pfCic04Cg]);
	photonPfIsoCharged03BadForCic_ -> push_back(myPhotons[gg].parameters[PARTICLE::pfCic03Cb]);
	photonPfIsoCharged04BadForCic_ -> push_back(myPhotons[gg].parameters[PARTICLE::pfCic04Cb]); 
	photonid_sieie_     -> push_back(myPhotons[gg].parameters[PARTICLE::sieie]);
	photonid_sieip_     -> push_back(myPhotons[gg].parameters[PARTICLE::sieip]);
	photonid_etawidth_  -> push_back(myPhotons[gg].parameters[PARTICLE::etaw]);
	photonid_phiwidth_   -> push_back(myPhotons[gg].parameters[PARTICLE::phiw]);
	photonid_r9_           -> push_back(myPhotons[gg].parameters[PARTICLE::r9]);
	photonid_lambdaRatio_  -> push_back(myPhotons[gg].parameters[PARTICLE::lR]);
	photonid_s4Ratio_      -> push_back(myPhotons[gg].parameters[PARTICLE::s4]);
	photonid_e25_      -> push_back(myPhotons[gg].parameters[PARTICLE::e25]);
	photonid_sceta_        -> push_back(myPhotons[gg].parameters[PARTICLE::sceta]);
	photonid_ESEffSigmaRR_ ->push_back(myPhotons[gg].parameters[PARTICLE::ESEff]);
	photonid_hadronicOverEm_->push_back(myPhotons[gg].parameters[PARTICLE::oldHadronicOverEm]);
	photonid_hadronicOverEm2012_->push_back(myPhotons[gg].parameters[PARTICLE::hadronicOverEm2012]);
	photonhcalTowerSumEtConeDR04_->push_back(myPhotons[gg].parameters[PARTICLE::hcalTowerSumEtConeDR04]);
	photonecalRecHitSumEtConeDR04_->push_back(myPhotons[gg].parameters[PARTICLE::ecalRecHitSumEtConeDR04]);
	photonnTrkSolidConeDR04_->push_back(myPhotons[gg].parameters[PARTICLE::nTrkSolidConeDR04]);
	photontrkSumPtSolidConeDR04_->push_back(myPhotons[gg].parameters[PARTICLE::trkSumPtSolidConeDR04]);
	photonnTrkHollowConeDR04_->push_back(myPhotons[gg].parameters[PARTICLE::nTrkHollowConeDR04]);
	photontrkSumPtHollowConeDR04_->push_back(myPhotons[gg].parameters[PARTICLE::trkSumPtHollowConeDR04]);
	photonBit_->push_back(myPhotons[gg].bit);
	photonIsoFPRCharged_->push_back(myPhotons[gg].isoFPRCharged);
	photonIsoFPRNeutral_->push_back(myPhotons[gg].isoFPRNeutral);
	photonIsoFPRPhoton_->push_back(myPhotons[gg].isoFPRPhoton);
      }



      vector<TLorentzVector> allP4;
      if(myLeptons.size()>1)allP4.push_back(llP4);
      if(myPhotons.size()>0)allP4.push_back(myPhotons[0].p4);
      double prod(1.0),sum(0.0);
      for(unsigned j = 0; j < myJets.size(); j++) {
        prod *= myJets[j].p4.Pt();
        sum  += myJets[j].p4.Pt();
        allP4.push_back(myJets[j].p4); 
        if(nLeptons_ > 1) jetllDPhi_     ->push_back(fabs(llP4.DeltaPhi(myJets[j].p4)));
        if(nPhotons_ > 0) jetPhotonDPhi_ ->push_back(fabs(myPhotons[0].p4.DeltaPhi(myJets[j].p4)));
        jetPt_       ->push_back(myJets[j].p4.Pt()); 
        jetEta_      ->push_back(myJets[j].p4.Eta()); 
        jetPhi_      ->push_back(myJets[j].p4.Phi()); 
        jetE_        ->push_back(myJets[j].p4.Energy()); 
        jetArea_     ->push_back(myJets[j].area);
        jetBeta_     ->push_back(myJets[j].beta);
        jetBtag_     ->push_back(myJets[j].btag);
	jetTagInfoNVtx_    ->push_back(myJets[j].taginfoNvtx);
	jetTagInfoNTracks_ ->push_back(myJets[j].taginfoNtracks);
	jetTagInfoVtxMass_ ->push_back(myJets[j].taginfoVtxMass);
	jetMCFlavour_ ->push_back(myJets[j].mcflavour);
        jetJEC_      ->push_back(myJets[j].jec);
        jetUNC_      ->push_back(myJets[j].unc);
        jetQGL_     ->push_back(myJets[j].qgl);
        jetRMS_     ->push_back(myJets[j].rms);
        jetVeto_     ->push_back(myJets[j].veto);
      }
      for(unsigned j = 0; j < myFwJets.size(); j++) {
	//fill only leading fw jets
	if(j<2){
	  fwjetPt_       ->push_back(myFwJets[j].p4.Pt()); 
	  fwjetEta_      ->push_back(myFwJets[j].p4.Eta()); 
	  fwjetPhi_      ->push_back(myFwJets[j].p4.Phi()); 
	  fwjetE_        ->push_back(myFwJets[j].p4.Energy());
	} 
      }
      

      sort(allP4.begin(),allP4.end(),p4SortingRule);

      // ----  Trigger block: Bother for trigger info only if event is selected (saves time)-------------
      iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
      if (!triggerResultsHandle_.isValid()) {
        cout << "ProcessedTreeProducer::analyze: Error in getting TriggerResults product from Event!" << endl;
        return;
      }
      iEvent.getByLabel(triggerEventTag_,triggerEventHandle_);
      if (!triggerEventHandle_.isValid()) {
        cout << "ProcessedTreeProducer::analyze: Error in getting TriggerEvent product from Event!" << endl;
        return;
      }
      // sanity check
      assert(triggerResultsHandle_->size() == hltConfig_.size());
      //------ loop over all trigger names ---------
      for(unsigned itrig=0;itrig<triggerNames_.size();itrig++) {
        bool accept(false);
        int preL1(-1);
        int preHLT(-1);
        int tmpFired(-1); 
        
    
        if (triggerIndex_[itrig] < hltConfig_.size()) {
          accept = triggerResultsHandle_->accept(triggerIndex_[itrig]);
          // --- check if your trigger bit is in the list which we don't ask for prescale (emu paths)
	  bool doCheckForPrescale = true;
          string reducedTriggerName = "";
          string reducedTriggerName2 = "";
	  int arraySize = int(triggerNamesFull_[itrig].size());
	  if(arraySize-1>0)
	    {reducedTriggerName=triggerNamesFull_[itrig].substr(0,triggerNamesFull_[itrig].size()-1); // remove last char from the str
	    reducedTriggerName2=triggerNamesFull_[itrig].substr(0,triggerNamesFull_[itrig].size()-2); // remove 2 last chars from the str
	    }
	  for(int nn = 0; nn<int(prescaleDontAsk_.size()); nn++) {
	    if(reducedTriggerName==prescaleDontAsk_[nn] || reducedTriggerName2==prescaleDontAsk_[nn]){
	      doCheckForPrescale=false;
	    }
	  }
	  //if(!doCheckForPrescale)cout << "skipping to check = " << triggerNamesFull_[itrig] << endl;

          if (triggerNamesFull_[itrig] != "" && doCheckForPrescale ) {
            const std::pair<int,int> prescales(hltConfig_.prescaleValues(iEvent,iSetup,triggerNamesFull_[itrig]));
            preL1  = prescales.first;
            preHLT = prescales.second;
          }  
          if (!accept) {
            tmpFired = 0;
          }
          else { 
            std::string ss(triggerNames_[itrig]); 
            hTriggerPass_->Fill((ss.erase(ss.find("v")-1,ss.find("v"))).c_str(),1);
            tmpFired = 1;
            // save trigger bit (0001) if family1 has fired, (0100) if family 3 has triggered
            isTriggered_ |= checkTriggerName(triggerNames_[itrig],triggerFamily1_) << 0; // if true 0001
            isTriggered_ |= checkTriggerName(triggerNames_[itrig],triggerFamily2_) << 1; // if true 0010
            isTriggered_ |= checkTriggerName(triggerNames_[itrig],triggerFamily3_) << 2; // if true 0100
            isTriggered_ |= checkTriggerName(triggerNames_[itrig],triggerFamily4_) << 3; // if true 1000
            isTriggered_ |= checkTriggerName(triggerNames_[itrig],triggerFamily5_) << 4; // if true 1000
            isTriggered_ |= checkTriggerName(triggerNames_[itrig],triggerFamily6_) << 5; // if true 1000
            isTriggered_ |= checkTriggerName(triggerNames_[itrig],triggerFamily7_) << 6; // if true 1000
            isTriggered_ |= checkTriggerName(triggerNames_[itrig],triggerFamily8_) << 7; // if true 1000
//          std::cout << "f1 " << checkTriggerName(triggerNames_[itrig],triggerFamily1_) << " " <<triggerNames_[itrig] << " " << isTriggered_ << std::endl;
//          std::cout << "f2 " << checkTriggerName(triggerNames_[itrig],triggerFamily2_) << " " <<triggerNames_[itrig] << " " << isTriggered_ << std::endl;
//          std::cout << "f3 " << checkTriggerName(triggerNames_[itrig],triggerFamily3_) << " " <<triggerNames_[itrig] << " " << isTriggered_ << std::endl;
//          std::cout << "f4 " << checkTriggerName(triggerNames_[itrig],triggerFamily4_) << " " <<triggerNames_[itrig] << " " << isTriggered_ << std::endl;
          }
        }
        fired_      ->push_back(tmpFired);
        prescaleL1_ ->push_back(preL1);
        prescaleHLT_->push_back(preHLT);
      }

      //------ trigger matching ---------
      if(myLeptons.size()>1){
        if(myLeptons[0].chid*myLeptons[1].chid == -4){
        	if(myLeptons[0].triObjMatchF1.size() && myLeptons[1].triObjMatchF1.size()){
            isTriggerMatchedFamily1_ = 1;
          } 
        }
        else if(myLeptons[0].chid*myLeptons[1].chid == -1){
        	if(myLeptons[0].triObjMatchF2.size() && myLeptons[1].triObjMatchF2.size()){
        		isTriggerMatchedFamily2_ = 1;
          }
        }
        else if(myLeptons[0].chid*myLeptons[1].chid == -2){
        	if((myLeptons[0].triObjMatchF3mu.size() && myLeptons[1].triObjMatchF3ele.size()) || (myLeptons[0].triObjMatchF3ele.size() && myLeptons[1].triObjMatchF3mu.size())){
        		isTriggerMatchedFamily3_ = 1;
          }
        }
        if((abs(myLeptons[0].chid) == 2 && myLeptons[0].triObjMatchF5.size()) || 
           (abs(myLeptons[1].chid) == 2 && myLeptons[1].triObjMatchF5.size()) || 
           ((abs(myLeptons[0].chid) == 2 && myLeptons[0].triObjMatchF5.size()) &&
            (abs(myLeptons[1].chid) == 2 && myLeptons[1].triObjMatchF5.size()))){
            isTriggerMatchedFamily5_ = 1; 
        }   
        if((abs(myLeptons[0].chid) == 1 && myLeptons[0].triObjMatchF6.size()) || 
           (abs(myLeptons[1].chid) == 1 && myLeptons[1].triObjMatchF6.size()) || 
           ((abs(myLeptons[0].chid) == 1 && myLeptons[0].triObjMatchF6.size()) &&
            (abs(myLeptons[1].chid) == 1 && myLeptons[1].triObjMatchF6.size()))){
            isTriggerMatchedFamily6_ = 1;
        }    
      }

      if(myPhotons.size()>0){
        if(myPhotons[0].triObjMatchF4.size()){
          isTriggerMatchedFamily4_ = 1;
        }  
      }	
    }// if selection GEN 
    if (selectionGEN) {

      if(llP4GEN.Pt()>0) {
       	llMGEN_                      = llP4GEN.M();
        llPtGEN_                     = llP4GEN.Pt();
        llPhiGEN_                    = llP4GEN.Phi();
        if(llPtGEN_>0)llYGEN_        = llP4GEN.Rapidity();
        if(llPtGEN_>0)llEtaGEN_      = llP4GEN.Eta();
        llDPhiGEN_                   = fabs(myGenLeptons[0].p4.DeltaPhi(myGenLeptons[1].p4));
      }
      TLorentzVector lepP4GEN(0,0,0,0); 

      for(unsigned l = 0; l < myGenLeptons.size(); l++) {
        lepP4GEN += myGenLeptons[l].p4;
        lepPtGEN_     ->push_back(myGenLeptons[l].p4.Pt());
        lepEtaGEN_    ->push_back(myGenLeptons[l].p4.Eta());
        lepPhiGEN_    ->push_back(myGenLeptons[l].p4.Phi());
        lepEGEN_      ->push_back(myGenLeptons[l].p4.Energy());
        lepChIdGEN_   ->push_back(myGenLeptons[l].pdgId);
      }   

      if(myGenPhotons.size()>0) {     
        photonPtGEN_     = myGenPhotons[0].p4.Pt();
        photonEtaGEN_    = myGenPhotons[0].p4.Eta();
        photonPhiGEN_    = myGenPhotons[0].p4.Phi();
        photonEGEN_      = myGenPhotons[0].p4.Energy();
	photonIsoPtDR03GEN_= myGenPhotons[0].isoPtDR03;
	photonIsoEDR03GEN_= myGenPhotons[0].isoEDR03;
	photonIsoPtDR04GEN_= myGenPhotons[0].isoPtDR04;
	photonIsoEDR04GEN_= myGenPhotons[0].isoEDR05;
	photonIsoPtDR05GEN_= myGenPhotons[0].isoPtDR05;
	photonIsoEDR05GEN_= myGenPhotons[0].isoEDR05;
	photonMotherIdGEN_ = myGenPhotons[0].motherId;
        if(myPhotons.size()>0)photonRECODRGEN_ = myPhotons[0].p4.DeltaR(myGenPhotons[0].p4); // GEN TO RECO matching

      }   

      vector<TLorentzVector> allP4GEN;
      allP4GEN.push_back(llP4GEN);
      double prod(1.0),sum(0.0);
      for(unsigned j = 0; j < myGenJets.size(); j++) {
        prod *= myGenJets[j].Pt();
        sum  += myGenJets[j].Pt();
        allP4GEN.push_back(myGenJets[j]); 
        jetllDPhiGEN_   ->push_back(fabs(llP4GEN.DeltaPhi(myGenJets[j])));
        jetPtGEN_       ->push_back(myGenJets[j].Pt()); 
        jetEtaGEN_      ->push_back(myGenJets[j].Eta()); 
        jetPhiGEN_      ->push_back(myGenJets[j].Phi()); 
        jetEGEN_        ->push_back(myGenJets[j].Energy()); 
        jetVetoGEN_     ->push_back(myGenJets[j].veto); 
	jetIdGEN_	->push_back(myGenJets[j].id);
	jetNpartonsGEN_->push_back(myGenJets[j].nparton);
      }
      sort(allP4GEN.begin(),allP4GEN.end(),p4SortingRule);
    }// if selection GEN 
    myTree_->Fill();
  }// if selectionRECO || selectionGEN
}
// ---- QG
vector<float> *PATZJetsExpress::ComputeQGVariables(edm::View<pat::Jet>::const_iterator &jet,const Event& iEvent,int index)
{
vector<float> *Vars=new vector<float>;
//const reco::PFJet& recoJet=dynamic_cast<const reco::PFJet&> (*( jet->originalObject()));
// Loop over candidates
float SumPt=0.;
float SumPt2=0.;
for(int i=0; i<jet->nConstituents();i++){
	SumPt+=jet->getJetConstituentsQuick()[i]->pt();
	SumPt2+=pow(jet->getJetConstituentsQuick()[i]->pt(),2);
	}

Handle<ValueMap<float> > qglMap;
Handle<ValueMap<float> > qglMapMLP;
   iEvent.getByLabel("QGTagger","qgLikelihood",qglMap);
   iEvent.getByLabel("QGTagger","qgMLP",qglMapMLP);

   Handle<edm::View<pat::Jet> > pfjetsakt5;
   iEvent.getByLabel(mJetsName, pfjetsakt5);
   
  edm::RefToBase<reco::Jet> jetRef(edm::Ref<edm::View<pat::Jet> >(pfjetsakt5,index)); 
if(jet->genParton() != NULL)Vars->push_back(  jet->genParton()->pdgId() ); 	//0 PART INFO
else Vars->push_back( 0 );
//Vars->push_back( jet->partonFlavour () );
Vars->push_back(  (*qglMap)[jetRef]  ); 		//1 QGL
Vars->push_back(  (*qglMapMLP)[jetRef]);
Vars->push_back(	jet->userFloat("ptD_QC") ); 	//6
Vars->push_back(	jet->userFloat("axis1_QC") ); 	//7
Vars->push_back(	jet->userFloat("axis2_QC") ); 	//8
Vars->push_back(	jet->userFloat("nNeutral_ptCut") ) ; 	//9
Vars->push_back(	jet->userFloat("nChg_QC") ) ; 	//9
Vars->push_back(-99);
return Vars;
}

// ---- method called once each job just after ending the event loop  ---
void PATZJetsExpress::endJob() 
{
}
// ---- method for tree building ----------------------------------------
void PATZJetsExpress::buildTree()
{
  QGVars_	     = new std::vector<float>();
  fired_             = new std::vector<int>();
  prescaleL1_        = new std::vector<int>();
  prescaleHLT_       = new std::vector<int>();
  lepPt_             = new std::vector<float>();
  lepEta_            = new std::vector<float>();
  lepPhi_            = new std::vector<float>();
  lepE_              = new std::vector<float>(); 
  lepPFIsoUnc_       = new std::vector<float>();
  lepPFIsoDBCor_     = new std::vector<float>();
  lepPFIsoRhoCor_    = new std::vector<float>();
  lepChId_           = new std::vector<int>();
  lepId_             = new std::vector<int>();
  lepSigmaIEtaIEta_  = new std::vector<float>();
  lepHadronicOverEm_ = new std::vector<float>();
  lepR9orChi2ndof_   = new std::vector<float>();
  fwjetPt_           = new std::vector<float>(); 
  fwjetEta_          = new std::vector<float>();
  fwjetPhi_          = new std::vector<float>();
  fwjetE_            = new std::vector<float>();
  jetPt_             = new std::vector<float>(); 
  jetEta_            = new std::vector<float>();
  jetPhi_            = new std::vector<float>();
  jetE_              = new std::vector<float>();
  jetArea_           = new std::vector<float>();
  jetBeta_           = new std::vector<float>();
  jetQGL_            = new std::vector<float>();
  jetRMS_            = new std::vector<float>();
  jetBtag_           = new std::vector<float>();
  jetTagInfoVtxMass_ = new std::vector<float>();
  jetTagInfoNTracks_ = new std::vector<float>();
  jetTagInfoNVtx_    = new std::vector<float>();
  jetMCFlavour_      = new std::vector<int>();
  jetJEC_            = new std::vector<float>();
  jetUNC_            = new std::vector<float>();
  jetVeto_           = new std::vector<int>(); 
  vtxZ_              = new std::vector<float>();
  vtxNdof_           = new std::vector<float>();
  lepPtGEN_          = new std::vector<float>();
  lepEtaGEN_         = new std::vector<float>();
  lepPhiGEN_         = new std::vector<float>();
  lepEGEN_           = new std::vector<float>(); 
  lepChIdGEN_        = new std::vector<int>();
  lepMatchedGEN_     = new std::vector<int>();
  lepMatchedDRGEN_   = new std::vector<float>();
  jetPtGEN_          = new std::vector<float>(); 
  jetEtaGEN_         = new std::vector<float>();
  jetPhiGEN_         = new std::vector<float>();
  jetEGEN_           = new std::vector<float>();
  jetVetoGEN_          = new std::vector<int>(); 
  jetIdGEN_          = new std::vector<int>(); 
  jetNpartonsGEN_    = new std::vector<int>(); 
  jetllDPhiGEN_      = new std::vector<float>();
  jetllDPhi_         = new std::vector<float>();
  jetPhotonDPhi_     = new std::vector<float>();

  photonPt_= new std::vector<float>();
  photonE_= new std::vector<float>();
  photonEta_= new std::vector<float>();
  photonPhi_= new std::vector<float>();
  photonBit_= new std::vector<int>();
  photonPassConversionVeto_= new std::vector<float>();
  photonPfIsoChargedHad_= new std::vector<float>();
  photonPfIsoNeutralHad_= new std::vector<float>();
  photonPfIsoPhoton_= new std::vector<float>();
  photonPfIsoPhotons03ForCic_= new std::vector<float>();
  photonPfIsoNeutrals03ForCic_= new std::vector<float>();
  photonPfIsoCharged03ForCicVtx0_= new std::vector<float>();
  photonPfIsoCharged03BadForCic_= new std::vector<float>();
  photonPfIsoPhotons04ForCic_= new std::vector<float>();
  photonPfIsoNeutrals04ForCic_= new std::vector<float>();
  photonPfIsoCharged04ForCicVtx0_= new std::vector<float>();
  photonPfIsoCharged04BadForCic_= new std::vector<float>();
  photonid_sieie_= new std::vector<float>();
  photonid_sieip_= new std::vector<float>();
  photonid_etawidth_= new std::vector<float>();
  photonid_phiwidth_= new std::vector<float>();
  photonid_r9_= new std::vector<float>();
  photonid_lambdaRatio_= new std::vector<float>();
  photonid_s4Ratio_= new std::vector<float>();
  photonid_e25_= new std::vector<float>();
  photonid_sceta_= new std::vector<float>();
  photonid_ESEffSigmaRR_= new std::vector<float>();
  photonid_hadronicOverEm_= new std::vector<float>();
  photonid_hadronicOverEm2012_= new std::vector<float>();
  photonhcalTowerSumEtConeDR04_= new std::vector<float>();
  photonecalRecHitSumEtConeDR04_= new std::vector<float>();
  photonnTrkSolidConeDR04_= new std::vector<float>();
  photontrkSumPtSolidConeDR04_= new std::vector<float>();
  photonnTrkHollowConeDR04_= new std::vector<float>();
  photontrkSumPtHollowConeDR04_= new std::vector<float>();
  photonIsoFPRCharged_ = new std::vector<float>();
  photonIsoFPRNeutral_ = new std::vector<float>();
  photonIsoFPRPhoton_ = new std::vector<float>();

  //  photonPar_         = new std::vector<float>();
//  FSRphotonPar_      = new std::vector<float>();
  // ---- global event variables ----------------------------------------
  myTree_->Branch("isRealData"       ,&isRealData_        ,"isRealData/I");
  myTree_->Branch("eventNum"         ,&eventNum_          ,"eventNum/l"); // 64 bit unsigned integer (ULong64_t)
  myTree_->Branch("runNum"           ,&runNum_            ,"runNum/I");
  myTree_->Branch("lumi"             ,&lumi_              ,"lumi/I");
  myTree_->Branch("nVtx"             ,&nVtx_              ,"nVtx/I");
  myTree_->Branch("nLeptons"         ,&nLeptons_          ,"nLeptons/I");
  myTree_->Branch("nPhotons"         ,&nPhotons_          ,"nPhotons/I");
  myTree_->Branch("nJets"            ,&nJets_             ,"nJets/I"); 
  myTree_->Branch("fwnJets"          ,&fwnJets_           ,"fwnJets/I"); 
  myTree_->Branch("rho"              ,&rho_               ,"rho/F");
  myTree_->Branch("rho25"            ,&rho25_             ,"rho25/F");
  // ---- met variables -------------------------------------------------
  myTree_->Branch("pfmet"            ,&pfmet_             ,"pfmet/F");
  myTree_->Branch("pfmetPhi"         ,&pfmetPhi_          ,"pfmetPhi/F");
  myTree_->Branch("pfhadPt"          ,&pfhadPt_           ,"pfhadPt/F");
  myTree_->Branch("pfSumEt"          ,&pfSumEt_           ,"pfSumEt/F");  
  // ---- dilepton variables --------------------------------------------
  myTree_->Branch("llM"              ,&llM_               ,"llM/F");
  myTree_->Branch("llPt"             ,&llPt_              ,"llPt/F");
  myTree_->Branch("llPhi"            ,&llPhi_             ,"llPhi/F");
  myTree_->Branch("llDPhi"           ,&llDPhi_            ,"llDPhi/F");
  myTree_->Branch("llY"              ,&llY_               ,"llY/F");
  myTree_->Branch("llEta"            ,&llEta_             ,"llEta/F");
  // ---- photon variables ----------------------------------------------
  // chiara
  myTree_->Branch("photonPt"         ,"vector<float>"     ,&photonPt_);
  myTree_->Branch("photonE"          ,"vector<float>"     ,&photonE_);
  myTree_->Branch("photonEta"        ,"vector<float>"     ,&photonEta_);
  myTree_->Branch("photonPhi"        ,"vector<float>"     ,&photonPhi_);
  if(!mReducedPh){
    myTree_->Branch("photonPassConversionVeto","vector<float>"   ,&photonPassConversionVeto_);
    myTree_->Branch("photonPfIsoChargedHad"   ,"vector<float>"   ,&photonPfIsoChargedHad_);  // chiara: li lasciamo? 
    myTree_->Branch("photonPfIsoNeutralHad"   ,"vector<float>"   ,&photonPfIsoNeutralHad_);  // chiara: li lasciamo? 
    myTree_->Branch("photonPfIsoPhoton"       ,"vector<float>"   ,&photonPfIsoPhoton_);      // chiara: li lasciamo? 
    myTree_->Branch("photonPfIsoPhoton03ForCic",      "vector<float>"   ,&photonPfIsoPhotons03ForCic_);      
    myTree_->Branch("photonPfIsoNeutrals03ForCic",    "vector<float>"   ,&photonPfIsoNeutrals03ForCic_);      
    myTree_->Branch("photonPfIsoCharged03ForCicVtx0", "vector<float>"   ,&photonPfIsoCharged03ForCicVtx0_);      
    myTree_->Branch("photonPfIsoCharged03BadForCic",  "vector<float>"   ,&photonPfIsoCharged03BadForCic_);      
    myTree_->Branch("photonPfIsoPhoton04ForCic",      "vector<float>"   ,&photonPfIsoPhotons04ForCic_);      
    myTree_->Branch("photonPfIsoNeutrals04ForCic",    "vector<float>"   ,&photonPfIsoNeutrals04ForCic_);      
    myTree_->Branch("photonPfIsoCharged04ForCicVtx0", "vector<float>"   ,&photonPfIsoCharged04ForCicVtx0_);      
    myTree_->Branch("photonPfIsoCharged04BadForCic",  "vector<float>"   ,&photonPfIsoCharged04BadForCic_);     
    myTree_->Branch("photonid_sieie",                 "vector<float>"   ,&photonid_sieie_);
    myTree_->Branch("photonid_sieip",                 "vector<float>"   ,&photonid_sieip_);
    myTree_->Branch("photonid_etawidth",              "vector<float>"   ,&photonid_etawidth_);
    myTree_->Branch("photonid_phiwidth",              "vector<float>"   ,&photonid_phiwidth_);  
    myTree_->Branch("photonid_r9",                    "vector<float>"   ,&photonid_r9_);
    myTree_->Branch("photonid_lambdaRatio",           "vector<float>"   ,&photonid_lambdaRatio_);
    myTree_->Branch("photonid_s4Ratio",               "vector<float>"   ,&photonid_s4Ratio_);
    myTree_->Branch("photonid_e25",               "vector<float>"   ,&photonid_e25_);
    myTree_->Branch("photonid_sceta",                 "vector<float>"   ,&photonid_sceta_);
    myTree_->Branch("photonid_ESEffSigmaRR",          "vector<float>"   ,&photonid_ESEffSigmaRR_);
    myTree_->Branch("photonid_hadronicOverEm",          "vector<float>"   ,&photonid_hadronicOverEm_);
    myTree_->Branch("photonid_hadronicOverEm2012",          "vector<float>"   ,&photonid_hadronicOverEm2012_);
    myTree_->Branch("photonhcalTowerSumEtConeDR04",          "vector<float>"   ,&photonhcalTowerSumEtConeDR04_);
    myTree_->Branch("photonecalRecHitSumEtConeDR04",          "vector<float>"   ,&photonecalRecHitSumEtConeDR04_);
    myTree_->Branch("photonnTrkSolidConeDR04",          "vector<float>"   ,&photonnTrkSolidConeDR04_);
    myTree_->Branch("photontrkSumPtSolidConeDR04",          "vector<float>"   ,&photontrkSumPtSolidConeDR04_);
    myTree_->Branch("photonnTrkHollowConeDR04",          "vector<float>"   ,&photonnTrkHollowConeDR04_);
    myTree_->Branch("photontrkSumPtHollowConeDR04",          "vector<float>"   ,&photontrkSumPtHollowConeDR04_);
    myTree_->Branch("photonIsoFPRCharged",          "vector<float>"   ,&photonIsoFPRCharged_);
    myTree_->Branch("photonIsoFPRNeutral",          "vector<float>"   ,&photonIsoFPRNeutral_);
    myTree_->Branch("photonIsoFPRPhoton",          "vector<float>"   ,&photonIsoFPRPhoton_);
    myTree_->Branch("photonBit","vector<int>",&photonBit_ );
  }
  // ---- trigger variables ---------------------------------------------
  myTree_->Branch("fired"            ,"vector<int>"       ,&fired_);
  myTree_->Branch("prescaleL1"       ,"vector<int>"       ,&prescaleL1_);
  myTree_->Branch("prescaleHLT"      ,"vector<int>"       ,&prescaleHLT_);
  myTree_->Branch("isTriggered"      ,&isTriggered_        ,"isTriggered/I");
  // ---- trigger matching variables ---------------------------------------------
  myTree_->Branch("isTriggerMatchedFamily1"          ,&isTriggerMatchedFamily1_           ,"isTriggerMatchedFamily1/I");
  myTree_->Branch("isTriggerMatchedFamily2"          ,&isTriggerMatchedFamily2_           ,"isTriggerMatchedFamily2/I");
  myTree_->Branch("isTriggerMatchedFamily3"          ,&isTriggerMatchedFamily3_           ,"isTriggerMatchedFamily3/I");
  myTree_->Branch("isTriggerMatchedFamily4"          ,&isTriggerMatchedFamily4_           ,"isTriggerMatchedFamily4/I");
  myTree_->Branch("isTriggerMatchedFamily5"          ,&isTriggerMatchedFamily5_           ,"isTriggerMatchedFamily5/I");
  myTree_->Branch("isTriggerMatchedFamily6"          ,&isTriggerMatchedFamily6_           ,"isTriggerMatchedFamily6/I");
  // ---- lepton variables ----------------------------------------------
  myTree_->Branch("lepPt"            ,"vector<float>"     ,&lepPt_);
  myTree_->Branch("lepEta"           ,"vector<float>"     ,&lepEta_);
  myTree_->Branch("lepPhi"           ,"vector<float>"     ,&lepPhi_);
  myTree_->Branch("lepE"             ,"vector<float>"     ,&lepE_);
  myTree_->Branch("lepPFIsoUnc"      ,"vector<float>"     ,&lepPFIsoUnc_);
  myTree_->Branch("lepPFIsoDBCor"    ,"vector<float>"     ,&lepPFIsoDBCor_);
  myTree_->Branch("lepPFIsoRhoCor"   ,"vector<float>"     ,&lepPFIsoRhoCor_);
  myTree_->Branch("lepChId"          ,"vector<int>"       ,&lepChId_);
  myTree_->Branch("lepR9orChi2ndof" ,"vector<float>"     ,&lepR9orChi2ndof_);
  myTree_->Branch("lepId"            ,"vector<int>"       ,&lepId_);
  // ---- jet variables -------------------------------------------------
  myTree_->Branch("jetVeto"          ,"vector<int>"       ,&jetVeto_);
  myTree_->Branch("jetPt"            ,"vector<float>"     ,&jetPt_);
  myTree_->Branch("jetEta"           ,"vector<float>"     ,&jetEta_);
  myTree_->Branch("jetPhi"           ,"vector<float>"     ,&jetPhi_);
  myTree_->Branch("jetE"             ,"vector<float>"     ,&jetE_);
  myTree_->Branch("jetArea"          ,"vector<float>"     ,&jetArea_);
  myTree_->Branch("jetBeta"          ,"vector<float>"     ,&jetBeta_);
  myTree_->Branch("jetQGL"           ,"vector<float>"     ,&jetQGL_);
  myTree_->Branch("jetRMS"           ,"vector<float>"     ,&jetRMS_);
  myTree_->Branch("jetBtag"          ,"vector<float>"     ,&jetBtag_);
  myTree_->Branch("jetTagInfoNVtx"   ,"vector<float>"     ,&jetTagInfoNVtx_);
  myTree_->Branch("jetTagInfoNTracks","vector<float>"     ,&jetTagInfoNTracks_);
  myTree_->Branch("jetTagInfoVtxMass","vector<float>"     ,&jetTagInfoVtxMass_);
  myTree_->Branch("jetMCFlavour"     ,"vector<int>"       ,&jetMCFlavour_);
  myTree_->Branch("jetJEC"           ,"vector<float>"     ,&jetJEC_);
  myTree_->Branch("jetUNC"           ,"vector<float>"     ,&jetUNC_);
  myTree_->Branch("jetllDPhi"        ,"vector<float>"     ,&jetllDPhi_);
  //-----forward jets - two leading forward jets
  myTree_->Branch("fwjetPt"          ,"vector<float>"     ,&fwjetPt_);
  myTree_->Branch("fwjetEta"         ,"vector<float>"     ,&fwjetEta_);
  myTree_->Branch("fwjetPhi"         ,"vector<float>"     ,&fwjetPhi_);
  myTree_->Branch("fwjetE"           ,"vector<float>"     ,&fwjetE_);
  // ---- vertex variables ----------------------------------------------
  myTree_->Branch("vtxZ"             ,"vector<float>"     ,&vtxZ_);
  myTree_->Branch("vtxNdof"          ,"vector<float>"     ,&vtxNdof_);
  // ---- gen variables ----------------------------------------------
  myTree_->Branch("puINT"            ,&puINT_             ,"puINT/I");
  myTree_->Branch("puOOT"            ,&puOOT_             ,"puOOT/I");
  myTree_->Branch("puTrueINT"        ,&puTrueINT_         ,"puTrueINT/I");
  myTree_->Branch("puTrueOOT"        ,&puTrueOOT_         ,"puTrueOOT/I");
  myTree_->Branch("nLeptonsGEN"      ,&nLeptonsGEN_       ,"nLeptonsGEN/I");
  myTree_->Branch("nJetsGEN"         ,&nJetsGEN_          ,"nJetsGEN/I");
  myTree_->Branch("llMGEN"           ,&llMGEN_            ,"llMGEN/F");
  myTree_->Branch("llPtGEN"          ,&llPtGEN_           ,"llPtGEN/F");
  myTree_->Branch("llPhiGEN"         ,&llPhiGEN_          ,"llPhiGEN/F");
  myTree_->Branch("llDPhiGEN"        ,&llDPhiGEN_         ,"llDPhiGEN/F");
  myTree_->Branch("llYGEN"           ,&llYGEN_            ,"llYGEN/F");
  myTree_->Branch("llEtaGEN"         ,&llEtaGEN_          ,"llEtaGEN/F");
  myTree_->Branch("lepPtGEN"         ,"vector<float>"     ,&lepPtGEN_);
  myTree_->Branch("lepEtaGEN"        ,"vector<float>"     ,&lepEtaGEN_);
  myTree_->Branch("lepPhiGEN"        ,"vector<float>"     ,&lepPhiGEN_);
  myTree_->Branch("lepEGEN"          ,"vector<float>"     ,&lepEGEN_);
  myTree_->Branch("lepChIdGEN"       ,"vector<int>"       ,&lepChIdGEN_);
  myTree_->Branch("lepMatchedDRGEN"  ,"vector<float>"     ,&lepMatchedDRGEN_);
  myTree_->Branch("lepMatchedGEN"    ,"vector<int>"       ,&lepMatchedGEN_);
  myTree_->Branch("jetPtGEN"         ,"vector<float>"     ,&jetPtGEN_);
  myTree_->Branch("jetEtaGEN"        ,"vector<float>"     ,&jetEtaGEN_);
  myTree_->Branch("jetPhiGEN"        ,"vector<float>"     ,&jetPhiGEN_);
  myTree_->Branch("jetEGEN"          ,"vector<float>"     ,&jetEGEN_);
  myTree_->Branch("jetVetoGEN"         ,"vector<int>"     ,&jetVetoGEN_);
  myTree_->Branch("jetllDPhiGEN"     ,"vector<float>"     ,&jetllDPhiGEN_);
  myTree_->Branch("jetIdGEN"         ,"vector<int>"     ,&jetIdGEN_);
  myTree_->Branch("jetNpartonsGEN"   ,"vector<int>"     ,&jetNpartonsGEN_);
  myTree_->Branch("HTParSum"         ,&HTParSum_          ,"HTParSum/F");  
  myTree_->Branch("mcWeight"         ,&mcWeight_          ,"mcWeight/F");
  myTree_->Branch("nPhotonsGEN"      ,&nPhotonsGEN_       ,"nPhotonsGEN/I");
  myTree_->Branch("photonPtGEN"      ,&photonPtGEN_       ,"photonPtGEN/F");
  myTree_->Branch("photonEGEN"       ,&photonEGEN_        ,"photonEGEN/F");
  myTree_->Branch("photonEtaGEN"     ,&photonEtaGEN_      ,"photonEtaGEN/F");
  myTree_->Branch("photonPhiGEN"     ,&photonPhiGEN_      ,"photonPhiGEN/F");
  myTree_->Branch("photonIsoPtDR03GEN"    ,&photonIsoPtDR03GEN_     ,"photonIsoPtDR03GEN/F");
  myTree_->Branch("photonIsoEDR03GEN"     ,&photonIsoEDR03GEN_     ,"photonIsoEDR03GEN/F");
  myTree_->Branch("photonIsoPtDR04GEN"    ,&photonIsoPtDR04GEN_     ,"photonIsoPtDR04GEN/F");
  myTree_->Branch("photonIsoEDR04GEN"     ,&photonIsoEDR04GEN_     ,"photonIsoEDR04GEN/F");
  myTree_->Branch("photonIsoPtDR05GEN"    ,&photonIsoPtDR05GEN_     ,"photonIsoPtDR05GEN/F");
  myTree_->Branch("photonIsoEDR05GEN"     ,&photonIsoEDR05GEN_     ,"photonIsoEDR05GEN/F");
  myTree_->Branch("photonMotherIdGEN",&photonMotherIdGEN_ ,"photonMotherIdGEN/I");
  myTree_->Branch("photonRECODRGEN"  ,&photonRECODRGEN_   ,"photonRECODRGEN/F");
  myTree_->Branch("VBPartonDM"       ,&VBPartonDM_        ,"VBPartonDM/I");
  myTree_->Branch("VBPartonM"        ,&VBPartonM_         ,"VBPartonM/F");
  myTree_->Branch("VBPartonE"        ,&VBPartonE_         ,"VBPartonE/F");
  myTree_->Branch("VBPartonPt"       ,&VBPartonPt_        ,"VBPartonPt/F");
  myTree_->Branch("VBPartonEta"      ,&VBPartonEta_       ,"VBPartonEta/F");
  myTree_->Branch("VBPartonPhi"      ,&VBPartonPhi_       ,"VBPartonPhi/F");
  myTree_->Branch("QGVars"	     ,"vector<float>"	  ,&QGVars_);
  myTree_->Branch("lepSigmaIEtaIEta" ,"vector<float>"	  ,&lepSigmaIEtaIEta_);
  myTree_->Branch("lepHadronicOverEm","vector<float>"	  ,&lepHadronicOverEm_);
}
// ---- method for tree initialization ----------------------------------
void PATZJetsExpress::clearTree()
{
   lepSigmaIEtaIEta_->clear();
   lepHadronicOverEm_->clear();

  QGVars_->clear();
  isRealData_        = -999;
  eventNum_          = 0; //ULong64_t can't be negative
  runNum_            = -999;
  lumi_              = -999;
  nVtx_              = -999;
  nLeptons_          = -999;
  nPhotonsGEN_       = -999;
  nPhotons_          = -999;
  nJets_             = -999;
  fwnJets_           = -999;
 // nRJets_            = -999;
  rho_               = -999;
  rho25_             = -999;
  pfmet_             = -999;
  pfmetPhi_          = -999;
  pfhadPt_           = -999;
  pfSumEt_           = -999;
  HTParSum_          = -999;
  llM_               = -999;
  llPt_              = -999; 
  llPhi_             = -999;
  llDPhi_            = -999;
  llY_               = -999;
  llEta_             = -999;

  //Photon variables
  photonE_           ->clear();
  photonPt_          ->clear();
  photonEta_         ->clear();
  photonPhi_         ->clear();
  photonBit_         ->clear();
  photonPassConversionVeto_ ->clear();
  photonPfIsoChargedHad_   ->clear();
  photonPfIsoNeutralHad_   ->clear();
  photonPfIsoPhoton_       ->clear();
  photonPfIsoPhotons03ForCic_ ->clear(); 
  photonPfIsoNeutrals03ForCic_ ->clear();
  photonPfIsoPhotons04ForCic_  ->clear();
  photonPfIsoNeutrals04ForCic_ ->clear();
  // chiara
  photonPfIsoCharged03ForCicVtx0_ ->clear();
  photonPfIsoCharged04ForCicVtx0_ ->clear();
  photonPfIsoCharged03BadForCic_  ->clear();
  photonPfIsoCharged04BadForCic_  ->clear();
  photonid_sieie_                 ->clear();
  photonid_sieip_                 ->clear();
  photonid_etawidth_              ->clear();
  photonid_phiwidth_              ->clear();
  photonid_r9_                    ->clear();
  photonid_lambdaRatio_           ->clear();
  photonid_s4Ratio_               ->clear();
  photonid_e25_                   ->clear();
  photonid_sceta_                 ->clear();
  photonid_ESEffSigmaRR_          ->clear();
  photonid_hadronicOverEm_        ->clear();
  photonid_hadronicOverEm2012_    ->clear();
  photonhcalTowerSumEtConeDR04_   ->clear();
  photonecalRecHitSumEtConeDR04_  ->clear();
  photonnTrkSolidConeDR04_        ->clear();
  photontrkSumPtSolidConeDR04_    ->clear();
  photonnTrkHollowConeDR04_       ->clear();
  photontrkSumPtHollowConeDR04_   ->clear();
  photonIsoFPRCharged_            ->clear();
  photonIsoFPRNeutral_            ->clear();
  photonIsoFPRPhoton_             ->clear();


  isTriggered_       =    0; // please keep this 0
  isTriggerMatchedFamily1_      =    0;
  isTriggerMatchedFamily2_      =    0;
  isTriggerMatchedFamily3_      =    0;
  isTriggerMatchedFamily4_      =    0;
  isTriggerMatchedFamily5_      =    0;
  isTriggerMatchedFamily6_      =    0;
  jetPhotonDPhi_     ->clear();
  //  photonPar_         ->clear();
  //FSRphotonPar_      ->clear();
  fired_             ->clear();
  prescaleL1_        ->clear();
  prescaleHLT_       ->clear();
  lepPt_             ->clear();
  lepEta_            ->clear();
  lepPhi_            ->clear();
  lepE_              ->clear();
  lepPFIsoUnc_       ->clear();
  lepPFIsoDBCor_     ->clear();
  lepPFIsoRhoCor_    ->clear();
  lepChId_           ->clear();
  lepR9orChi2ndof_   ->clear();
  lepMatchedDRGEN_   ->clear();
  lepMatchedGEN_     ->clear();
  lepId_             ->clear();
  jetPt_             ->clear();
  jetEta_            ->clear();
  jetPhi_            ->clear();
  jetE_              ->clear();
  jetArea_           ->clear();
  jetBeta_           ->clear();
  jetQGL_            ->clear();
  jetRMS_            ->clear();
  jetBtag_           ->clear();
  jetTagInfoVtxMass_ ->clear();
  jetTagInfoNTracks_ ->clear();
  jetTagInfoNVtx_    ->clear();
  jetMCFlavour_      ->clear();
  jetJEC_            ->clear();
  jetUNC_            ->clear();
  jetllDPhi_         ->clear();
  jetVeto_           ->clear();
  fwjetPt_           ->clear();
  fwjetEta_          ->clear();
  fwjetPhi_          ->clear();
  fwjetE_            ->clear();
  vtxZ_              ->clear();
  vtxNdof_           ->clear();
  puINT_             = -999;
  puOOT_             = -999;
  puTrueINT_         = -999;
  puTrueOOT_         = -999;
  isRealData_        = -999;
  nLeptonsGEN_       = -999;
  nJetsGEN_          = -999;
  llMGEN_            = -999;
  llPtGEN_           = -999; 
  llPhiGEN_          = -999;
  llDPhiGEN_         = -999;
  llYGEN_            = -999;
  llEtaGEN_          = -999;
  lepPtGEN_          ->clear();
  lepEtaGEN_         ->clear();
  lepPhiGEN_         ->clear();
  lepEGEN_           ->clear();
  lepChIdGEN_        ->clear();
  jetPtGEN_          ->clear();
  jetEtaGEN_         ->clear();
  jetPhiGEN_         ->clear();
  jetEGEN_           ->clear();
  jetVetoGEN_        ->clear();
  jetllDPhiGEN_      ->clear();
  jetIdGEN_          ->clear();
  jetNpartonsGEN_    ->clear();
  photonEGEN_        = -999;
  photonPtGEN_       = -999;
  photonEtaGEN_      = -999;
  photonPhiGEN_      = -999;
  photonMotherIdGEN_ = -999;
  photonIsoEDR03GEN_ = -999; 
  photonIsoPtDR03GEN_= -999;
  photonIsoEDR04GEN_ = -999; 
  photonIsoPtDR04GEN_= -999; 
  photonIsoEDR05GEN_ = -999; 
  photonIsoPtDR05GEN_= -999; 
  photonRECODRGEN_   = +999; // please keep this positive (will cut offline to <0.2 for matched)
  VBPartonDM_        = -999;
  VBPartonM_         = -999;
  VBPartonE_         = -999;
  VBPartonPt_        = -999;
  VBPartonEta_       = -999;
  VBPartonPhi_       = -999;
  mcWeight_          = -999;
}

double PATZJetsExpress::getEffectiveAreaForMuons(const double& eta) const {
  double abseta = abs(eta);
  // values for 0.4 cone
  
  if(abseta <= 1.)                        return 0.674;
  else if(abseta > 1.  && abseta <= 1.5)  return 0.565;
  else if(abseta > 1.5 && abseta <= 2.0)  return 0.442;
  else if(abseta > 2.0 && abseta <= 2.2)  return 0.515;
  else if(abseta > 2.2 && abseta <=  2.3) return 0.821;
  else if(abseta > 2.3 && abseta <= 2.4)  return 0.660;
  else                                    return 9999;
}

double PATZJetsExpress::getEffectiveAreaForElectrons(const double& eta) const {
  double abseta = abs(eta);
    // values for 0.3 cone
  
  if(abseta <= 1.)                            return 0.13;
  else if(abseta > 1.    && abseta <= 1.479)  return 0.14;
  else if(abseta > 1.479 && abseta <= 2.0)    return 0.07;
  else if(abseta > 2.    && abseta <= 2.2)    return 0.09;
  else if(abseta > 2.2   && abseta <=  2.3)   return 0.11;
  else if(abseta > 2.3   && abseta <= 2.4)    return 0.11;
  else if(abseta > 2.4)                       return 0.14;
  else                                        return 9999;
}

// ---- define this as a plug-in ----------------------------------------
DEFINE_FWK_MODULE(PATZJetsExpress);
