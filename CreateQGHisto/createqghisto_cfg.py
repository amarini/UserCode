import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")                               
# ---- access the global tag (needed for the JEC) -----------------------
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_42_V19::All'
# ---- load the reco jet configurations --------------------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')      
process.load('RecoJets.Configuration.RecoJets_cff')
# ---- load the JEC services --------------------------------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')    
# ---- format the message service ---------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000    
# ---- maximum number of events to run over -----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(30000))                                             
# ---- define the source ------------------------------------------------
process.source = cms.Source("PoolSource",                  
    fileNames = cms.untracked.vstring(
       # 'file:///tmp/amarini/0E578C40-D05C-E011-B8F6-0017A4771000.root'
	'rfio:///castor/cern.ch/cms/store/mc/Fall11/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6/AODSIM/RecoTest_PU_S5_START44_V4-v1/0000/F474926B-AACF-E011-B064-0018F3D0960C.root/castor/cern.ch/cms/store/mc/Fall11/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6/AODSIM/RecoTest_PU_S5_START44_V4-v1/0000/F474926B-AACF-E011-B064-0018F3D0960C.root'
    )
)
# ---- define the output file -------------------------------------------
process.TFileService = cms.Service("TFileService",         
    fileName = cms.string("jets.root"),
    closeFileFast = cms.untracked.bool(True)               
)
# ---- ZJetsExpress analyzer --------------------------------------------                                              
process.accepted = cms.EDAnalyzer('CreateQGHisto')
############# turn-on the fastjet area calculation needed for the L1Fastjet ##############                             
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.Rho_EtaMax = cms.double(5.0)             
process.ak5PFJets.doAreaFastjet = True                     
process.ak5PFJets.Rho_EtaMax = cms.double(5.0)             
############# slimming the PFJet collection by raising the pt cut #################
process.ak5PFJets.jetPtMin = cms.double(15.0)              
##############
#PAT

from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.trigTools import *
from PhysicsTools.PatAlgos.tools.metTools import *

process.load("PhysicsTools.PatAlgos.patSequences_cff")

addJetCollection(
                process,
                cms.InputTag('ak5PFJets'),
                'AK5', 'PF' ,
                doJTA        = True,
                doBTagging   = True,
                jetCorrLabel = ('AK5PF',cms.vstring(['L2Relative','L3Absolute'])),
                doType1MET   = True,
                doL1Cleaning = True,
                doL1Counters = False,
                genJetCollection=cms.InputTag("ak5GenJets"),
                doJetID      = True,
                jetIdLabel   = 'ak5',
                outputModules    = []
                )
addJetCollection(
                process,
                cms.InputTag('ak5PFJetsNoPU'),
                'AK5', 'PFNOPU' ,
                doJTA        = True,
                doBTagging   = True,
                jetCorrLabel = ('AK5PF',cms.vstring(['L2Relative','L3Absolute'])),
                doType1MET   = True,
                doL1Cleaning = True,
                doL1Counters = False,
                genJetCollection=cms.InputTag("ak5GenJets"),
                doJetID      = True,
                jetIdLabel   = 'ak5noPU',
                outputModules    = []
                )
process.patJets.addTagInfos = cms.bool(True)
process.patJets.addGenPartonMatch               = cms.bool(True)
process.patJets.embedGenPartonMatch             = cms.bool(True)

#switchJetCollection(process,
#                    cms.InputTag('ak5PFJets'),
#                    doJTA            = True,
#                    doBTagging       = True,
#                    jetCorrLabel     = ('AK5PF', ['L2Relative', 'L3Absolute']),
#                    doType1MET       = False,
#                    genJetCollection = cms.InputTag("ak5GenJets"),
#                    doJetID      = True,
#                    jetIdLabel   = "ak5",
#                    outputModules        = []
#                    )
#
#
    
#process.s1 = cms.Sequence(process.kt6PFJets + process.ak5PFJets + process.accepted)                                   
process.s1 = cms.Sequence(
	#process.genParticlesForJets
        #+ process.ak5GenJets
         process.makePatJets 
	+process.accepted )                
process.p1 = cms.Path(process.s1)


