from PhysicsTools.PatAlgos.patTemplate_cfg import *

isMC=True;

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
process.GlobalTag.globaltag = 'GR_P_V41_AN1::All'
#process.GlobalTag.globaltag = 'GR_R_42_V23::All'
#process.GlobalTag.globaltag = 'FT_R_52_V8D::All'
#process.GlobalTag.globaltag = 'START52_V9::All'

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
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
# ---- load geometry package --------------------------------------------
process.load("Configuration.StandardSequences.Geometry_cff")
# ---- maximum number of events to run over -----------------------------
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
# ---- define the source ------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#	'/store/relval/CMSSW_5_3_6-START53_V14/RelValProdTTbar/AODSIM/v2/00000/76ED0FA6-1E2A-E211-B8F1-001A92971B72.root'
	'/store/relval/CMSSW_5_3_6-START53_V14/RelValH130GGgluonfusion/GEN-SIM-RECO/v2/00000/202DD4DB-F929-E211-8F53-001A92810AF2.root'
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/752/504D95A3-789B-E111-9B6C-003048D3C944.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/774/0CDC3936-889B-E111-9F82-001D09F25041.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/806/C8EE0F38-C89B-E111-9623-001D09F29619.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/812/A66C70B4-D09B-E111-9589-001D09F24664.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/818/5AF2B989-E69B-E111-9ADF-0019B9F72F97.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/822/0E10C6A2-EF9B-E111-B28A-001D09F291D7.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/828/0A3345B7-009C-E111-9EAB-001D09F25041.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/829/B0B54F12-049C-E111-871F-0030486780B4.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/830/1635A18E-059C-E111-B6BE-003048D2BED6.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/833/DA7007CB-1C9C-E111-B053-001D09F242EF.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/834/0281DB15-5C9C-E111-A58E-001D09F26509.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/835/7A99502C-4F9C-E111-9D71-0025901D624A.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/836/3E413843-3F9C-E111-8435-003048D37560.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/840/DE717A5B-1E9C-E111-B03A-001D09F251FE.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/842/9CE9122E-1E9C-E111-8E93-00237DDBE41A.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/844/7C90F03A-2F9C-E111-8DA4-BCAEC53296F8.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/849/2A223F7B-389C-E111-AA16-001D09F2905B.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/852/2EA28439-389C-E111-A2A4-001D09F28EA3.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/865/0A28DF77-429C-E111-B75E-0025901D5DB8.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/870/8ED00431-469C-E111-A16F-001D09F27067.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/871/2E0CDCD0-489C-E111-9A00-001D09F2305C.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/875/3CC38E9E-499C-E111-BD9B-003048D3C932.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/877/ECD954C2-4D9C-E111-94FF-BCAEC518FF44.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/878/F87C3D36-B09C-E111-8CDE-0025901D6288.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/895/A25B9129-639C-E111-A056-003048F1183E.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/898/440D88FF-F19C-E111-8173-0019B9F581C9.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/913/84BC6806-D59C-E111-98AB-001D09F25479.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/917/26695F97-C29C-E111-A66A-5404A63886AF.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/919/A8CF7C3B-DC9C-E111-B2B7-001D09F28F25.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/922/68762909-D39C-E111-BA15-00215AEDFD98.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/925/64F05EB7-DD9C-E111-A2D4-0025B32036D2.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/928/10D325A5-C29C-E111-BFA0-003048D2BC38.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/928/1E55A719-009D-E111-A80F-001D09F29619.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/928/A24670EC-D59C-E111-A453-5404A63886D6.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/928/F44F623C-D49C-E111-9331-003048D3C90E.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/935/6C232E04-9B9C-E111-945A-00237DDC5C24.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/938/DE593E7D-AF9C-E111-99CC-0025901D626C.root', 
# '/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/193/942/6EEA2E74-B59D-E111-A05D-002481E0D7EC.root'	
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
	wantSummary = cms.untracked.bool(True)
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
    minPhoPt        = cms.double(130),
    maxPhoEta       = cms.double(3.0),
    minLepPt        = cms.double(20),
    maxLepEta       = cms.double(2.4),
    maxCombRelIso03 = cms.double(0.15),
    maxCombRelIso04 = cms.double(0.12),
    minLLMass       = cms.double(40),
    OnlyMC	    = cms.bool(False), 
    dressedRadius   = cms.double(0.1),
   # GENCrossCleaning= cms.int32(1), #OBSOLETE
    GENType	    = cms.int32(1), #
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

				  hggPhotonIDConfiguration = cms.PSet(hggPhotonIDCuts),
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
process.p = cms.Path(process.pfParticleSelectionSequence
		     + process.eleIsoSequence
		     + process.phoIsoSequence
		     + process.kt6PFJets
		     + process.ak5PFJets
		     + process.kt6PFJets25
		     + process.goodOfflinePrimaryVertices)

if(isMC):
	process.p += process.genParticlesForJets
process.tail = cms.Sequence(process.patDefaultSequence  + process.jetExtender + process.qglAK5PF +  process.accepted)

process.p += process.tail



