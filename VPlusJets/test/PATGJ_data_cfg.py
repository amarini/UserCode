from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
process.load('RecoJets.JetProducers.TrackJetParameters_cfi')
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
#process.load("RecoJets.Configuration.GenJetParticles_cff")
#process.load("RecoJets.Configuration.GenJetParticles_cff")
#process.load('RecoJets.Configuration.RecoGenJets_cff')

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

###################### MC
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/004CA74B-3198-E111-9786-003048FF9AC6.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/02A11927-2D98-E111-902E-002354EF3BD0.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/04A8C6D7-4D98-E111-ABAD-002618943978.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/0EFD3E1B-3D98-E111-AB78-002618FDA211.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/1210DFA2-3098-E111-8D18-00261894382D.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/1294263F-2898-E111-930C-003048678FD6.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/12AD5991-4D98-E111-880F-002618943860.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/18409BDE-4D98-E111-B43A-001A92810AE0.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/1A8127B9-4E98-E111-9076-0026189438C9.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/1AB69670-5C98-E111-81DC-003048678F9C.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/20B397A7-3598-E111-B4CF-002618943950.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/228DCD00-4D98-E111-B83E-003048D15E14.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/22F2F7E8-3998-E111-B04B-0026189438DD.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/24D83FE6-7A98-E111-97BA-0026189438A5.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/267CFD91-8498-E111-9D16-002618943963.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/2879A188-3198-E111-AFD5-0026189438D3.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/2A9304D4-6198-E111-824A-00304867BFA8.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/2C007199-2598-E111-85DE-002618943947.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/2E74EAE7-5A98-E111-8650-002618FDA248.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/3093779C-7A98-E111-A28A-002618943842.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/343A7770-7998-E111-8CD0-00304867C1BC.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/3440F927-5898-E111-8B67-0018F3D096A6.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/3A07C1E3-5798-E111-BACB-001A928116BC.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/409D4126-4E98-E111-A6C9-003048678B3C.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/445A4B3A-5B98-E111-8874-003048FFD728.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/4678FC9E-2A98-E111-B915-003048FFCB6A.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/4AF02F72-3498-E111-B89C-003048678A76.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/56A6B554-3798-E111-BB98-001A928116CE.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/58356EDE-4F98-E111-87B8-00261894385D.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/625F2B0E-2A98-E111-B10C-003048FFCB84.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/6CA86817-3B98-E111-B991-003048FFD720.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/6E868516-3098-E111-848A-00248C55CC97.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/70C0B391-4098-E111-B3C2-001A92971B7E.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/725B10AC-B898-E111-99C3-003048FFD796.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/72B9B5BB-5A98-E111-9BF6-003048678B04.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/762F5EF8-5A98-E111-9147-0026189438FD.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/7C68AA01-3598-E111-99D5-003048678FF8.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/7E97FC0F-5398-E111-95EC-001A92971B8C.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/8031749B-5E98-E111-9E2D-00261894388A.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/824A1D54-7D98-E111-8ED7-002618FDA279.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/864613E5-4B98-E111-971E-003048FFD7C2.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/868DA08D-5B98-E111-B44D-003048FFCBF0.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/8696787E-F298-E111-BA4A-001A928116BC.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/88C2B77D-5998-E111-B703-00261894389C.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/8ED0E125-5598-E111-AE8C-003048678B16.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/926FFEE7-5098-E111-9E70-002618FDA265.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/948C7DFB-5398-E111-A49C-00261894383C.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/9C205367-2B98-E111-B60E-002618943810.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/9E17136B-7898-E111-8627-001A92810AA2.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/9E8821C7-3898-E111-9F15-003048FFD7A2.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/A2D4AF50-3598-E111-A741-0026189438F8.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/AAA4D11B-5298-E111-B590-0026189438CC.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/ACE6BFFF-5498-E111-9606-002618943977.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/AE2DC452-5098-E111-8609-002618943922.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/B0BC7200-4F98-E111-8D34-002618943831.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/B0CE1A25-3498-E111-98F1-002618FDA211.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/B270074A-3398-E111-B27B-0026189438AC.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/B6F8102F-3698-E111-9F2D-00261894383B.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/BA037BAF-3298-E111-A133-00261894398C.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/BA8B8D9C-3798-E111-BCED-0026189438C0.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/C06AF954-2E98-E111-B1CC-0026189437FE.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/C26262C6-3698-E111-A58F-003048678D6C.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/C6E2021E-6098-E111-B996-002618943934.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/CA6214CA-2898-E111-B04A-003048678BF4.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/CC2E2912-5898-E111-84D7-003048678BE6.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/DC238E4A-5B98-E111-B4D7-003048678E8A.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/DC8D5AE9-7B98-E111-8E27-002618FDA262.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/E0F70543-4F98-E111-856D-002618943966.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/E2A79357-5A98-E111-9482-00261894392B.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/EA9DB2BC-2D98-E111-83B8-002618943902.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/EADD3D24-3298-E111-A0A2-001BFCDBD160.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/ECDCF32E-2F98-E111-945B-003048678B20.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/F4692886-7798-E111-AD06-0026189438B4.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/F6A7736A-4E98-E111-B899-002618943923.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/FA42324D-2C98-E111-978A-00304867920C.root',
#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/FCEC475A-3098-E111-A0A1-0026189438D3.root'


###2011
#    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/06634156-5FA4-E011-BEA2-003048D293B4.root',
#    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/0660AC04-B59C-E011-BEC7-E0CB4E19F9B4.root',
#    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/062DC08D-AF9C-E011-B713-90E6BA0D09B9.root',
#    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/061080FF-A99C-E011-90CD-0030487CDB24.root',
#    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/04F65B00-A89C-E011-978E-E0CB4E5536EF.root',
#    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/04B76CEA-AF9C-E011-8526-E0CB4E1A116A.root',
#    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/047CEEC1-57A4-E011-86DC-0030487C6A1E.root'

#    '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/06B46752-04CD-E011-B8B9-0015178C4994.root'
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/06883EDE-09CD-E011-899C-00A0D1EEE660.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/066F1232-05CD-E011-A14B-00A0D1EE8F34.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/065D3A98-0BCD-E011-8325-00A0D1EE95CC.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/0648740D-FECC-E011-A2C5-0024E876636C.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/063052DF-FCCC-E011-976C-00A0D1EE8A20.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/04920741-0DCD-E011-84FC-0024E8769B05.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/04437AFF-07CD-E011-8675-001D09676BAB.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/02A456DD-F8CC-E011-B0DF-00151796D680.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/02A3A00C-FECC-E011-82CC-0015178C0100.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/027A27F3-FFCC-E011-8210-00266CF97FF4.root'
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
removeMCMatching(process)
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

#process.patJets.addGenPartonMatch               = cms.bool(True)
#process.patJets.embedGenPartonMatch             = cms.bool(True)
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
    maxCombRelIso03 = cms.double(0.12),
    minLLMass       = cms.double(40),
    OnlyMC	    = cms.bool(False), #NOT IMPL
    #GENCrossCleaning= cms.int32(2), #1 leptons 2 Photons 3 both
    GENType	    = cms.int32(0),
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
	+  process.patDefaultSequence  + process.jetExtender + process.qglAK5PF +  process.accepted)

# ---- first sequence: events that pass the trigger ---------------------
#process.s1 = cms.Sequence(process.hltFilter + process.kt6PFJets + process.ak5PFJets + process.accepted)
#process.p1 = cms.Path(process.s1)
# ---- second sequence: events that fail the trigger --------------------
#process.s2 = cms.Sequence(~process.hltFilter + process.kt6PFJets + process.ak5PFJets + process.rejected)
#process.p2 = cms.Path(process.s2)



# ---- schedule ---------------------------------------------------------
#process.schedule = cms.Schedule(process.p1,process.p2)
