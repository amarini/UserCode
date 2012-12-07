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
process.GlobalTag.globaltag = 'START52_V9::All'

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
process.MessageLogger.cerr.FwkReport.reportEvery = 10
# ---- load geometry package --------------------------------------------
process.load("Configuration.StandardSequences.Geometry_cff")
# ---- maximum number of events to run over -----------------------------
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
# ---- define the source ------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#	'/store/mc/Summer12/DY1JetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph/DQM/PU_S7_START52_V9-v1/0000/6E20576C-9CA1-E111-90B8-20CF3027A5B7.root',
#'/store/mc/Summer12/DY1JetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph/DQM/PU_S7_START52_V9-v1/0000/B25B2ACE-1196-E111-B0DD-20CF305B0509.root',
#'/store/mc/Summer12/DY1JetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph/DQM/PU_S7_START52_V9-v1/0000/364A10D7-9B96-E111-910A-E0CB4E1A117B.root'

#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/752/504D95A3-789B-E111-9B6C-003048D3C944.root', 
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/774/0CDC3936-889B-E111-9F82-001D09F25041.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/806/C8EE0F38-C89B-E111-9623-001D09F29619.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/812/A66C70B4-D09B-E111-9589-001D09F24664.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/818/5AF2B989-E69B-E111-9ADF-0019B9F72F97.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/822/0E10C6A2-EF9B-E111-B28A-001D09F291D7.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/828/0A3345B7-009C-E111-9EAB-001D09F25041.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/829/B0B54F12-049C-E111-871F-0030486780B4.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/830/1635A18E-059C-E111-B6BE-003048D2BED6.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/833/DA7007CB-1C9C-E111-B053-001D09F242EF.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/834/0281DB15-5C9C-E111-A58E-001D09F26509.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/835/7A99502C-4F9C-E111-9D71-0025901D624A.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/836/3E413843-3F9C-E111-8435-003048D37560.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/840/DE717A5B-1E9C-E111-B03A-001D09F251FE.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/842/9CE9122E-1E9C-E111-8E93-00237DDBE41A.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/844/7C90F03A-2F9C-E111-8DA4-BCAEC53296F8.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/849/2A223F7B-389C-E111-AA16-001D09F2905B.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/852/2EA28439-389C-E111-A2A4-001D09F28EA3.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/865/0A28DF77-429C-E111-B75E-0025901D5DB8.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/870/8ED00431-469C-E111-A16F-001D09F27067.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/871/2E0CDCD0-489C-E111-9A00-001D09F2305C.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/875/3CC38E9E-499C-E111-BD9B-003048D3C932.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/877/ECD954C2-4D9C-E111-94FF-BCAEC518FF44.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/878/F87C3D36-B09C-E111-8CDE-0025901D6288.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/895/A25B9129-639C-E111-A056-003048F1183E.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/898/440D88FF-F19C-E111-8173-0019B9F581C9.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/913/84BC6806-D59C-E111-98AB-001D09F25479.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/917/26695F97-C29C-E111-A66A-5404A63886AF.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/919/A8CF7C3B-DC9C-E111-B2B7-001D09F28F25.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/922/68762909-D39C-E111-BA15-00215AEDFD98.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/925/64F05EB7-DD9C-E111-A2D4-0025B32036D2.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/928/10D325A5-C29C-E111-BFA0-003048D2BC38.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/928/1E55A719-009D-E111-A80F-001D09F29619.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/928/A24670EC-D59C-E111-A453-5404A63886D6.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/928/F44F623C-D49C-E111-9331-003048D3C90E.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/935/6C232E04-9B9C-E111-945A-00237DDC5C24.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/938/DE593E7D-AF9C-E111-99CC-0025901D626C.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/942/6EEA2E74-B59D-E111-A05D-002481E0D7EC.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/944/9C48C67E-B59C-E111-A3D6-0025901D6272.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/948/9AF7F90E-B89C-E111-9A69-0030486730C6.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/950/EA96A3D5-D59C-E111-971B-E0CB4E5536AE.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/955/7047695F-D79C-E111-A9F9-002481E0D958.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/958/0CD07946-E39C-E111-9EEF-001D09F2437B.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/959/C43A38F1-E39C-E111-85CD-5404A63886BE.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/962/86010CFE-DE9C-E111-831F-003048678098.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/983/660A543F-EA9C-E111-A4CA-001D09F23A20.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/984/10BDF483-E99C-E111-A34B-0015C5FDE067.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/987/2685AB5A-F19C-E111-99BC-001D09F27067.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/992/B26B342F-F49C-E111-9753-001D09F2512C.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/996/5C231B02-F79C-E111-AE6F-001D09F2512C.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/998/00F7A237-679D-E111-A8EB-003048D3C932.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/998/12C002D3-599D-E111-8662-003048D2BC30.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/998/1A5B38AB-499D-E111-BB3F-BCAEC5364C4C.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/998/30597AAA-499D-E111-9D8F-001D09F24FEC.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/998/6EF1FB00-469D-E111-AB23-003048D2BEAA.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/998/EEB9C87E-4E9D-E111-B4A2-003048D37524.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/998/FC46DA74-B59D-E111-B3FE-002481E0D7C0.root',
#'/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/999/3C2C8C48-B59D-E111-9371-0025901D624E.root'



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
		btagdiscriminators = ['jetBProbabilityBJetTags','jetProbabilityBJetTags','trackCountingHighPurBJetTags','trackCountingHighEffBJetTags','simpleSecondaryVertexHighEffBJetTags','simpleSecondaryVertexHighPurBJetTags','combinedSecondaryVertexBJetTags','combinedSecondaryVertexMVABJetTags','softMuonBJetTags','softMuonByPtBJetTags','softMuonByIP3dBJetTags','simpleSecondaryVertexNegativeHighEffBJetTags','simpleSecondaryVertexNegativeHighPurBJetTags']
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
    fileName = cms.string("PATZJetsExpress.root"),
    closeFileFast = cms.untracked.bool(True)
)
# ---- ZJetsExpress analyzer --------------------------------------------
process.accepted = cms.EDAnalyzer('PATZJetsExpress',
    jets            = cms.InputTag('jetExtender','extendedPatJets'),
    srcRho          = cms.InputTag('kt6PFJets','rho'),
    srcRho25        = cms.InputTag('kt6PFJets25','rho'),
    minNjets        = cms.int32(1),
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
   # GENCrossCleaning= cms.int32(1),
    GENType	    = cms.int32(1),
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
                                  'HLT_Photon135_v'),
  
    triggerFamily5  = cms.vstring('HLT_SingleMu_v'), ##NAMES TO BE CHECKED!
    triggerFamily6  = cms.vstring('HLT_SingleE_v'), ##TO BE CHECKED!
    triggerFamily7  = cms.vstring('HLT_SingleE_v'), ##TO BE CHECKED!
    triggerFamily8  = cms.vstring([]), ##TO BE CHECKED!

    prescaleDontAsk = cms.vstring('HLT_Mu17_Ele8_CaloIdL_v', # don't ask for L1 prescales for these bits
                                  'HLT_Mu8_Ele17_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v',
                                  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v'),

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
