from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
process.load('RecoJets.JetProducers.TrackJetParameters_cfi')
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load('RecoJets.Configuration.RecoGenJets_cff')

from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

# ---- access the global tag (needed for the JEC) -----------------------
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'GR_R_42_V23::All'
#process.GlobalTag.globaltag = 'FT_R_52_V8D::All'
process.GlobalTag.globaltag = 'START53_V7A::All'

##--------- good primary vertices ---------------
process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
    src          = cms.InputTag('offlinePrimaryVertices'),
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) )
)

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
process.load("Configuration.Geometry.GeometryIdeal_cff")
# ---- maximum number of events to run over -----------------------------
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
# ---- define the source ------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
 "file:/scratch0/webermat/test_DYJetsToLL_madgraph_Z2star_8TeV_ptZ_100_Summer12_START53_V7A-v2.root"
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
#removeMCMatching(process)
addPfMET(process, 'PF')
switchJetCollection(process,cms.InputTag('ak5PFJets'),
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'])),
		 genJetCollection=cms.InputTag('ak5GenJets'),
                 doType1MET   = False,
                 doJetID      = True,
                 )

process.selectedPatJets.cut        = "pt > 10 && abs(eta) < 4.7"

##--------- keep only jet and MET PAT objects ---
#removeAllPATObjectsBut(process,["Jets","METs"])

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
    fileName = cms.string("DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10-START53_V7A-v2_PREFSR.root"),
    closeFileFast = cms.untracked.bool(True)
)
# ---- ZJetsExpress analyzer --------------------------------------------
process.accepted = cms.EDAnalyzer('PATZJetsExpress',
    jets            = cms.InputTag('jetExtender','extendedPatJets'),
    srcRho          = cms.InputTag('kt6PFJets','rho'),
    srcRho25        = cms.InputTag('kt6PFJets25','rho'),
    minNjets        = cms.int32(1),
    dressedRadius   = cms.double(0.15),                              
    jetLepIsoRadius = cms.double(0.4),
    jetLepPhoRadius = cms.double(0.4),
    minJetPt        = cms.double(30),
    maxJetEta       = cms.double(2.5),
    minPhoPt        = cms.double(130),                          
    maxPhoEta       = cms.double(3.0),
    minLepPt        = cms.double(20),
    maxLepEta       = cms.double(2.4),
    maxCombRelIso03 = cms.double(0.15),
    minLLMass       = cms.double(40),
    OnlyMC	    = cms.bool(False), #NOT IMPL
    #GENCrossCleaning= cms.int32(1),#2 for photons, 1 for leptons
    GENType	    = cms.int32(0),#0 GEN - 1 JEF  - 2 PAR
    processName     = cms.string('HLT'),
    triggerName     = cms.vstring('HLT_DoubleMu6_v', ### IMPORTANT always put _v in the end of each bit 
                                  'HLT_DoubleMu7_v',
                                  'HLT_Mu13_Mu8_v',
                                  'HLT_Mu17_Mu8_v',
                                  'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v',
                                  'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v',
                                  'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
                                  'HLT_Mu17_Ele8_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v',
                                  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v',
				  #'HLT_Photon20_CaloIdVL_IsoL_v',
				  #'HLT_Photon30_CaloIdVL_v',
				  #'HLT_Photon30_CaloIdVL_IsoL_v', // we start offline at 70 GeV
				  'HLT_Photon50_CaloIdVL_IsoL_v',
				  'HLT_Photon75_CaloIdVL_v',
				  'HLT_Photon90_CaloIdVL_v',
				  'HLT_Photon90_CaloIdVL_IsoL_v',
				  'HLT_Photon125_v',
				  'HLT_Photon135_v',
         	                 ),
    triggerResults  = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    triggerFamily1  = cms.vstring('HLT_DoubleMu6_v','HLT_DoubleMu7_v','HLT_Mu13_Mu8_v','HLT_Mu17_Mu8_v'),
    triggerFamily2  = cms.vstring('HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v',
                                  'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v',
                                  'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v'),
    triggerFamily3  = cms.vstring('HLT_Mu17_Ele8_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v',
                                  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v'),
    triggerFamily4  = cms.vstring('HLT_Photon50_CaloIdVL_IsoL_v',
                                  'HLT_Photon75_CaloIdVL_v',
                                  'HLT_Photon90_CaloIdVL_v',
                                  'HLT_Photon90_CaloIdVL_IsoL_v',
                                  'HLT_Photon125_v',
                                  'HLT_Photon135_v'
                                  ),
    triggerFamily5  = cms.vstring('HLT_IsoMu20_WCandPt80_v',
                                  'HLT_IsoMu24_v',
                                  'HLT_IsoMu30_v',
                                  'HLT_Mu40_v',
                                  'HLT_Mu5_v',
                                  'HLT_Mu24_v',
                                  'HLT_Mu30_v'),                                 
    triggerFamily6  = cms.vstring('HLT_Ele27_WP80_CentralPFJet80_v',
                                  'HLT_Ele27_WP80_WCandPt80_v',
                                  'HLT_Ele27_WP80_v',
                                  'HLT_Ele80_CaloIdVT_GsfTrkIdT_v',
                                  'HLT_Ele90_CaloIdVT_GsfTrkIdT_v',
                                  'HLT_Ele30_CaloIdVT_TrkIdT_v',
                                  'HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v'),
    triggerFamily7  = cms.vstring('HLT_PFJet320_v',
                                  'HLT_PFJet200_v',
                                  'HLT_Jet370_NoJetID_v'),                                
    triggerFamily8  = cms.vstring(),  
    prescaleDontAsk = cms.vstring('HLT_Mu17_Ele8_CaloIdL_v', # don't ask for L1 prescales for these bits
                                  'HLT_Mu8_Ele17_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v',
                                  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v',
                                  'HLT_Photon50_CaloIdVL_IsoL_v',
                                  'HLT_Photon75_CaloIdVL_v',
                                  'HLT_Photon90_CaloIdVL_v',
                                  'HLT_Photon90_CaloIdVL_IsoL_v',
                                  ),


)
# ---- duplicate the analyzer with different name -----------------------
process.rejected = process.accepted.clone()
# ---- filter the required HLT bits -------------------------------------
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring(
                                     'HLT_DoubleMu6_v*', # di-muon triggers
                                     'HLT_DoubleMu7_v*',
                                     'HLT_Mu13_Mu8_v*',
                                     'HLT_Mu17_Mu8*', 
                                     'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*', #di-electron trigger
                                     'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*',
                                     'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
				     'HLT_Mu17_Ele8_CaloIdL_v*', #di-emu trigger (for data-driven ttbar estimation)
                                     'HLT_Mu8_Ele17_CaloIdL_v*',
                                     'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*',
                         	     'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*'
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

del process.outpath

# ---- save all events for any trigger ---------------------
#GEN
#process.p = cms.Path(process.genParticlesForJets+process.ak5GenJets+  process.accepted)
#RECO, data and MC
process.p = cms.Path(process.kt6PFJets + process.ak5PFJets + process.kt6PFJets25 + process.goodOfflinePrimaryVertices 
	+ process.genParticlesForJets
	+  process.patDefaultSequence  + process.jetExtender + process.qglAK5PF +  process.accepted)

# ---- first sequence: events that pass the trigger ---------------------
#process.s1 = cms.Sequence(process.hltFilter + process.kt6PFJets + process.ak5PFJets + process.accepted)
#process.p1 = cms.Path(process.s1)
# ---- second sequence: events that fail the trigger --------------------
#process.s2 = cms.Sequence(~process.hltFilter + process.kt6PFJets + process.ak5PFJets + process.rejected)
#process.p2 = cms.Path(process.s2)



# ---- schedule ---------------------------------------------------------
#process.schedule = cms.Schedule(process.p1,process.p2)
