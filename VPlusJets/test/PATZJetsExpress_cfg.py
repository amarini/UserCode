from PhysicsTools.PatAlgos.patTemplate_cfg import *

isMC=False;

process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
process.load('RecoJets.JetProducers.TrackJetParameters_cfi')
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
if(isMC):
	process.load("RecoJets.Configuration.GenJetParticles_cff")
	process.load("RecoJets.Configuration.GenJetParticles_cff")
	process.load('RecoJets.Configuration.RecoGenJets_cff')

from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

# ---- access the global tag (needed for the JEC) -----------------------
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
if(isMC):
	process.GlobalTag.globaltag = 'START53_V16::All'    # MC 53Y release
else:
	#process.GlobalTag.globaltag = 'FT_53_V6C_AN3::All' # DATA 2012B: multijet parked dataset re-reco with CMSSW_5_3_6 
	#process.GlobalTag.globaltag = 'FT53_V10A_AN3:All'  # DATA 2012Cv1: multijet parked dataset re-reco with CMSSW_5_3_6
	#process.GlobalTag.globaltag = 'FT_P_V42C_AN3::All' # DATA 2012Cv2: multijet parked dataset re-reco with CMSSW_5_3_6
	process.GlobalTag.globaltag = 'FT_P_V42_AN3::All'  # DATA 2012D: multijet parked dataset re-reco with CMSSW_5_3_6

##--------- good primary vertices ---------------
process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
    src          = cms.InputTag('offlinePrimaryVertices'),
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) )
)

from amarini.VPlusJets.hggPhotonIDCuts_cfi import *

##--------- remove cleaning --------------------
removeCleaning(process)
##--------- jets -------------------------------
process.patJets.embedPFCandidates = False
process.patJets.embedCaloTowers = False
process.patJets.addTagInfos = True

# ---- format the message service ---------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
# ---- load geometry package --------------------------------------------
#process.load("Configuration.StandardSequences.Geometry_cff")
# ---- maximum number of events to run over -----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
# ---- define the source ------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	#'file:/scratch0/webermat/DYJets_MadGraph_START_53_V7A.root'
	#'file:/scratch0/webermat/WJets_MadGraph_START_53_V7A.root'
#	'/store/relval/CMSSW_5_3_6-START53_V14/RelValProdTTbar/AODSIM/v2/00000/76ED0FA6-1E2A-E211-B8F1-001A92971B72.root'
	#'/store/relval/CMSSW_5_3_6-START53_V14/RelValH130GGgluonfusion/GEN-SIM-RECO/v2/00000/202DD4DB-F929-E211-8F53-001A92810AF2.root'
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/752/504D95A3-789B-E111-9B6C-003048D3C944.root', 
    )
)
## 
process.qglAK5PF   = cms.EDProducer("QuarkGluonTagger",
         # jets     = cms.InputTag("ak5PFJets"),
   	 # jets            = cms.InputTag('extendedPatJets'),
   	  jets            = cms.InputTag('jetExtender','extendedPatJets'),
          rho      = cms.InputTag('kt6PFJets','rho'),
          jec      = cms.string('ak5PFL1FastL2L3'),
	  isPatJet = cms.bool(True),
          #jec      = cms.string('ak5PFL1FastL2L3Residual'),
)


##--------- remove MC matching -----------------
if not isMC:
	removeMCMatching(process)

addPfMET(process, 'PF')

if(isMC):
	Corrections=cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])
else:
	Corrections=cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'])

#BTagDiscriminators= ['jetBProbabilityBJetTags','jetProbabilityBJetTags','trackCountingHighPurBJetTags','trackCountingHighEffBJetTags','simpleSecondaryVertexHighEffBJetTags','simpleSecondaryVertexHighPurBJetTags','combinedSecondaryVertexBJetTags','combinedSecondaryVertexMVABJetTags','softMuonBJetTags','softMuonByPtBJetTags','softMuonByIP3dBJetTags','simpleSecondaryVertexNegativeHighEffBJetTags','simpleSecondaryVertexNegativeHighPurBJetTags']
BTagDiscriminators= ['jetBProbabilityBJetTags','jetProbabilityBJetTags','trackCountingHighPurBJetTags','trackCountingHighEffBJetTags','simpleSecondaryVertexHighEffBJetTags','simpleSecondaryVertexHighPurBJetTags','combinedSecondaryVertexBJetTags','combinedSecondaryVertexMVABJetTags']

switchJetCollection(process,cms.InputTag('ak5PFJets'),
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', Corrections),
		 genJetCollection=cms.InputTag('ak5GenJets'),
                 doType1MET   = False,
                 doJetID      = True,
		btagdiscriminators = BTagDiscriminators 
                 )

process.selectedPatJets.cut        = "pt > 10 && abs(eta) < 4.7"

##--------- keep only jet and MET PAT objects ---
#removeAllPATObjectsBut(process,["Jets","METs"])
if(isMC):
	process.patJets.addGenPartonMatch               = cms.bool(True)
	process.patJets.embedGenPartonMatch             = cms.bool(True)
#process.patJets.genPartonMatch  = cms.InputTag('patJetPartonMatch');

process.jetExtender = cms.EDProducer("JetExtendedProducer",
    jets    = cms.InputTag('selectedPatJets'),
    result  = cms.string('extendedPatJets'),
    payload = cms.string('AK5PF')
)

# ---- define the output file -------------------------------------------
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("PATZJetsExpress.root"),
    closeFileFast = cms.untracked.bool(True)
)


# ---- Gen-Jet flavour matching -----------------------------
process.load("PhysicsTools.JetMCAlgos.SelectPartons_cff")

process.jetPartonAssociationAK5PF = cms.EDProducer("JetPartonMatcher",
						   jets = cms.InputTag("ak5PFJets"),
						   coneSizeToAssociate = cms.double(0.3),
						   partons = cms.InputTag("myPartons")
						   )

process.jetFlavourAssociationAK5PF = cms.EDProducer("JetFlavourIdentifier",
						    srcByReference = cms.InputTag("patJetPartonAssociationAK5PF"),
						    physicsDefinition = cms.bool(False)
						    )

process.GenJetFlavourMatching = cms.Sequence(process.myPartons*process.jetPartonAssociationAK5PF*process.jetFlavourAssociationAK5PF)


# ---- Recompute electron's PF iso deposits -----------------------------
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.phoIsoSequence = setupPFPhotonIso(process, 'photons')

process.options = cms.untracked.PSet(
	wantSummary = cms.untracked.bool(False)
	)

# ---- ZJetsExpress analyzer --------------------------------------------
process.accepted = cms.EDAnalyzer('PATZJetsExpress',
    jets            = cms.InputTag('jetExtender','extendedPatJets'),
    srcRho          = cms.InputTag('kt6PFJets','rho'),
    srcRho25        = cms.InputTag('kt6PFJetsCentralNeutral','rho'),
    pfIsoValEleCH03 = cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
    pfIsoValEleNH03 = cms.InputTag('elPFIsoValueNeutral03PFIdPFIso'),			  
    pfIsoValEleG03  = cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),                                    
    minNjets        = cms.int32(1),
    jetLepIsoRadius = cms.double(0.4),
    jetLepPhoRadius = cms.double(0.4),
    minJetPt        = cms.double(30),
    maxJetEta       = cms.double(2.5),
    minPhoPt        = cms.double(50),
    minPhoPtId      = cms.double(30),
    maxPhoEta       = cms.double(3.0),
    minLepPt        = cms.double(20),
    maxLepEta       = cms.double(2.4),
    maxCombRelIso03 = cms.double(0.15),
    maxCombRelIso04 = cms.double(0.12),
    minLLMass       = cms.double(15),
    OnlyMC	    = cms.bool(False), 
    dressedRadius   = cms.double(0.1),
   # GENCrossCleaning= cms.int32(1), #OBSOLETE
    GENType	    = cms.int32(0), #
    processName     = cms.string('HLT'),
    triggerName     = cms.vstring('HLT_DoubleMu7_v',
                                  'HLT_Mu13_Mu8_v',
                                  'HLT_Mu17_Mu8_v',
                                  'HLT_Mu17_TkMu8_v', #end 2011
                                  'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v',
                                  'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v', # end 2011
                                  'HLT_Mu17_Ele8_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v',
                                  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v',
				  'HLT_DoubleMu5_Ele8_CaloIdT_TrkIdVL_v',
				  'HLT_DoubleMu8_Ele8_CaloIdT_TrkIdVL_v',
				  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
				  'HLT_Mu30_Ele30_CaloIdL_v',
				  'HLT_Mu7_Ele7_CaloIdT_CaloIsoVL_v',
				  'HLT_Mu8_DoubleEle8_CaloIdT_TrkIdVL_v',
				  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
				  'HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Ele8_CaloIdL_TrkIdVL_v',
				  'HLT_Photon30_v',
				  'HLT_Photon30_CaloIdVL_v',
				  'HLT_Photon50_CaloIdVL_IsoL_v',
				  'HLT_Photon75_CaloIdVL_v',
				  'HLT_Photon90_CaloIdVL_v',
				  'HLT_Photon90_CaloIdVL_IsoL_v',
				  'HLT_Photon125_v',
				  'HLT_Photon135_v',
                                  'HLT_Mu15_v',
                                  'HLT_Mu24_v',
                                  'HLT_Mu30_v',
                                  'HLT_Mu40_v',
                                  'HLT_Mu40_eta2p1_v',
                                  'HLT_IsoMu17_v',
                                  'HLT_IsoMu20_v',
                                  'HLT_IsoMu24_v',
                                  'HLT_IsoMu24_eta2p1_v', # end 2011
                                  'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v',
                                  'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v',
                                  'HLT_Ele52_CaloIdVT_TrkIdT_v',
                                  'HLT_Ele65_CaloIdVT_TrkIdT_v',
                                  'HLT_Ele80_CaloIdVT_TrkIdT_v', # end of 2011
                                  'HLT_Ele27_WP80_v'
         	                 ),                     
    triggerResults  = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    triggerFamily1  = cms.vstring('HLT_DoubleMu7_v',
                                  'HLT_Mu13_Mu8_v',
                                  'HLT_Mu17_Mu8_v',
                                  'HLT_Mu17_TkMu8_v'),
    triggerFamily2  = cms.vstring('HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v',
                                  'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v'),
    triggerFamily3  = cms.vstring('HLT_Mu17_Ele8_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v',
                                  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v',
				  'HLT_DoubleMu5_Ele8_CaloIdT_TrkIdVL_v',
				  'HLT_DoubleMu8_Ele8_CaloIdT_TrkIdVL_v',
				  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
				  'HLT_Mu30_Ele30_CaloIdL_v',
				  'HLT_Mu7_Ele7_CaloIdT_CaloIsoVL_v',
				  'HLT_Mu8_DoubleEle8_CaloIdT_TrkIdVL_v',
				  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
				  'HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Ele8_CaloIdL_TrkIdVL_v'),
    triggerFamily4  = cms.vstring('HLT_Photon30_CaloIdVL_v',
				  'HLT_Photon30_v',
				  'HLT_Photon50_CaloIdVL_IsoL_v',
                                  'HLT_Photon75_CaloIdVL_v',
                                  'HLT_Photon90_CaloIdVL_v',
                                  'HLT_Photon90_CaloIdVL_IsoL_v',
                                  'HLT_Photon125_v',
                                  'HLT_Photon135_v'),
    triggerFamily5  = cms.vstring('HLT_Mu15_v',
                                  'HLT_Mu24_v',
                                  'HLT_Mu30_v',
                                  'HLT_Mu40_v',
                                  'HLT_Mu40_eta2p1_v',
                                  'HLT_IsoMu17_v',
                                  'HLT_IsoMu20_v',
                                  'HLT_IsoMu24_v',
                                  'HLT_IsoMu24_eta2p1_v', # end 2011
                                  'HLT_DoubleMu7_v',
                                  'HLT_Mu13_Mu8_v',
                                  'HLT_Mu17_Mu8_v',
                                  'HLT_Mu17_TkMu8_v'), # end 2011
    triggerFamily6  = cms.vstring('HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v',
                                  'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v',
                                  'HLT_Ele52_CaloIdVT_TrkIdT_v',
                                  'HLT_Ele65_CaloIdVT_TrkIdT_v',
                                  'HLT_Ele80_CaloIdVT_TrkIdT_v', # end of 2011
                                  'HLT_Ele27_WP80_v',
                                  'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v',
                                  'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v'), # end 2011
    triggerFamily7  = cms.vstring('HLT_SingleE_v'), ##TO BE CHECKED!
    triggerFamily8  = cms.vstring([]), ##TO BE CHECKED!

				  hggPhotonIDConfiguration = cms.PSet(hggPhotonIDCuts),
    prescaleDontAsk = cms.vstring('HLT_Mu17_Ele8_CaloIdL_v', # don't ask for L1 prescales for these bits
                                  'HLT_Mu8_Ele17_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v',
                                  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v',
				  'HLT_Ele27_WP80_v',
			          'HLT_Photon30_CaloIdVL_v',
				  'HLT_Photon50_CaloIdVL_IsoL_v',
                                  'HLT_Photon75_CaloIdVL_v',
                                  'HLT_Photon90_CaloIdVL_v',
                                  'HLT_Photon90_CaloIdVL_IsoL_v',
				  'HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Ele8_CaloIdL_TrkIdVL_v',
				  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
				  'HLT_DoubleMu8_Ele8_CaloIdT_TrkIdVL_v',
				  'HLT_Mu8_DoubleEle8_CaloIdT_TrkIdVL_v'),
)
# ---- duplicate the analyzer with different name -----------------------
process.rejected = process.accepted.clone()
# ---- filter the required HLT bits -------------------------------------
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring(
                                     'HLT_DoubleMu7_v*',
                                     'HLT_Mu13_Mu8_v*',
                                     'HLT_Mu17_Mu8_v*',
                                     'HLT_Mu17_TkMu8_v*', #end 2011
                                     'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*',
                                     'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*', # end 2011
                                     'HLT_Mu17_Ele8_CaloIdL_v*',
                                     'HLT_Mu8_Ele17_CaloIdL_v*',
                                     'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*',
                                     'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*',
                                     'HLT_Mu15_v2*',
                                     'HLT_Mu24_v*',
                                     'HLT_Mu30_v*',
                                     'HLT_Mu40_v*',
                                     'HLT_Mu40_eta2p1_v*',
                                     'HLT_IsoMu17_v*',
                                     'HLT_IsoMu20_v*',
                                     'HLT_IsoMu24_v*',
                                     'HLT_IsoMu24_eta2p1_v*', # end 2011
                                     'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*',
                                     'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*',
                                     'HLT_Ele52_CaloIdVT_TrkIdT_v*',
                                     'HLT_Ele65_CaloIdVT_TrkIdT_v*',
                                     'HLT_Ele80_CaloIdVT_TrkIdT_v*', # end of 2011
                                     'HLT_Ele27_WP80_v*'
),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)
############# turn-on the fastjet area calculation needed for the L1Fastjet ##############
process.kt6PFJets.doRhoFastjet = True
process.ak5PFJets.doAreaFastjet = True
############# turn-on the fastjet area in |eta|<2.5 needed for the photonISO #############
process.kt6PFJets25 = process.kt6PFJets.clone( doRhoFastjet = True )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)

############# Trigger Matching #############
pathTriggerMuons = '(path("HLT_DoubleMu7_v*", 1, 1) || path("HLT_Mu13_Mu8_v*", 1, 1) || path("HLT_Mu17_Mu8_v*", 1, 1) || path("HLT_Mu17_TkMu8_v*", 1, 1) || path("HLT_Mu15_v2*", 1, 1) || path("HLT_Mu24_v*", 1, 1) || path("HLT_Mu30_v*", 1, 1) || path("HLT_Mu40_v*", 1, 1) || path("HLT_Mu40_eta2p1_v*", 1, 1) || path("HLT_IsoMu17_v*", 1, 1) || path("HLT_IsoMu20_v*", 1, 1) || path("HLT_IsoMu24_v*", 1, 1) || path("HLT_IsoMu24_eta2p1_v*", 1, 1))' # selecting the trigger objects

pathTriggerElectrons = '(path("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*", 1, 1) || path("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", 1, 1)  || path("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*", 1, 1) || path("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*", 1, 1) || path("HLT_Ele52_CaloIdVT_TrkIdT_v*", 1, 1) || path("HLT_Ele65_CaloIdVT_TrkIdT_v*", 1, 1) || path("HLT_Ele80_CaloIdVT_TrkIdT_v*", 1, 1) || path("HLT_Ele27_WP80_v*", 1, 1))' 

pathTriggerPhotons = '(path("HLT_Photon30_v*", 1, 1) ||path("HLT_Photon30_CaloIdVL_v*", 1, 1) || path("HLT_Photon50_CaloIdVL_IsoL_v*", 1, 1) || path("HLT_Photon75_CaloIdVL_v*", 1, 1) || path("HLT_Photon90_CaloIdVL_v*", 1, 1) || path("HLT_Photon90_CaloIdVL_IsoL_v*", 1, 1) || path("HLT_Photon125_v*", 1, 1) || path("HLT_Photon135_v*", 1, 1))'

process.muonTriggerMatchHLTMuons = cms.EDProducer(
    "PATTriggerMatcherDRLessByR", 
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string(pathTriggerMuons),
    maxDPtRel             = cms.double(0.5),
    maxDeltaR             = cms.double(0.5),
    # maxDeltaEta         = cms.double(0.2),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(True)
)

process.electronTriggerMatchHLTElectrons = cms.EDProducer(
    "PATTriggerMatcherDRLessByR",
    src                   = cms.InputTag('cleanPatElectrons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string(pathTriggerElectrons),
    maxDPtRel             = cms.double(0.5),
    maxDeltaR             = cms.double(0.5),
    # maxDeltaEta         = cms.double(0.2),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(True)
)

process.photonTriggerMatchHLTPhotons = cms.EDProducer(
    "PATTriggerMatcherDRLessByR", 
    src                   = cms.InputTag('cleanPatPhotons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string(pathTriggerPhotons),
    maxDPtRel             = cms.double(1.0),
    maxDeltaR             = cms.double(0.2),
    # maxDeltaEta         = cms.double(0.2),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(True)
)

from PhysicsTools.PatAlgos.tools.coreTools import removeCleaning
removeCleaning(process)

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger(process) 

#switchOnTriggerMatching(process, triggerMatchers = ['muonTriggerMatchHLTMuons', 'electronTriggerMatchHLTElectrons', 'photonTriggerMatchHLTPhotons'])
switchOnTriggerMatchEmbedding(process , triggerMatchers = ['muonTriggerMatchHLTMuons', 'electronTriggerMatchHLTElectrons', 'photonTriggerMatchHLTPhotons'])

# Switch to selected PAT objects in the trigger matching
removeCleaningFromTriggerMatching(process)

process.muonTriggerMatchHLTMuons.src         = cms.InputTag('patMuons')
process.selectedPatMuonsTriggerMatch.src     = cms.InputTag('patMuons')
process.selectedPatMuonsTriggerMatch.matches = cms.VInputTag('muonTriggerMatchHLTMuons')

process.electronTriggerMatchHLTElectrons.src     = cms.InputTag('patElectrons')
process.selectedPatElectronsTriggerMatch.src     = cms.InputTag('patElectrons')
process.selectedPatElectronsTriggerMatch.matches = cms.VInputTag('electronTriggerMatchHLTElectrons')

process.photonTriggerMatchHLTPhotons.src      = cms.InputTag('patPhotons')
process.selectedPatPhotonsTriggerMatch.src    = cms.InputTag('patPhotons')
process.selectedPatPhotonsTriggerMatch.matches = cms.VInputTag('photonTriggerMatchHLTPhotons')

process.patTriggerSequence = cms.Sequence(process.patTrigger + process.muonTriggerMatchHLTMuons + process.selectedPatMuonsTriggerMatch + process.electronTriggerMatchHLTElectrons + process.selectedPatElectronsTriggerMatch + process.photonTriggerMatchHLTPhotons + process.selectedPatPhotonsTriggerMatch)


del process.outpath

# ---- save all events for any trigger ---------------------
process.p = cms.Path(process.pfParticleSelectionSequence
		     + process.eleIsoSequence
		     + process.phoIsoSequence
		     + process.kt6PFJets
		     + process.ak5PFJets
		     + process.kt6PFJets25
		     + process.goodOfflinePrimaryVertices)

if(isMC):
	process.p += process.genParticlesForJets
process.tail = cms.Sequence(process.patDefaultSequence + process.jetExtender + process.qglAK5PF + process.accepted)

process.p += process.tail



