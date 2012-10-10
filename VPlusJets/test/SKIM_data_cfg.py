import FWCore.ParameterSet.Config as cms
process = cms.Process("FILTER")
# ---- access the global tag (needed for the JEC) -----------------------
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'GR_R_42_V19::All'
# ---- load the JEC services --------------------------------------------
#process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
# ---- format the message service ---------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
# ---- maximum number of events to run over -----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
# ---- define the source ------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring([
	'/store/data/Run2012C/DoubleMu/AOD/PromptReco-v1/000/197/770/7EF0FF77-77C3-E111-941A-5404A63886EE.root',
	'/store/data/Run2012C/DoubleMu/AOD/PromptReco-v1/000/197/772/8E02C4D0-79C3-E111-9773-001D09F28EA3.root',
	'/store/data/Run2012C/DoubleMu/AOD/PromptReco-v1/000/197/774/3C5FCCD6-79C3-E111-8249-001D09F25479.root'
])
)
## ---- define the output file -------------------------------------------
#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string("ZJetsFilter.root"),
#    closeFileFast = cms.untracked.bool(True)
#)
# ---- filter the required HLT bits -------------------------------------
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring('HLT_DoubleMu6_v*','HLT_DoubleMu7_v*','HLT_DoubleMu8_v*','HLT_Mu13_Mu8_v*'),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)
# ---- filter the required leptons -------------------------------------
process.muonsFilter = cms.EDFilter('MuonRefSelector',
	src	=	cms.InputTag("muons"),
	cut	= 	cms.string("abs(eta)<2.5 && (isGlobalMuon || isTrackerMuon) && pt>15"),
)
process.electronsFilter = cms.EDFilter('GsfElectronRefSelector',
	src	=	cms.InputTag("gsElectrons"),
	cut	= 	cms.string("abs(eta)<2.5 && pt>15"),
)
process.MuMuFilter = cms.EDFilter('CandViewCountFilter',
	src	=	cms.InputTag("muonsFilter"),
	checkCharge = 	cms.bool(False),
	minNumber =	cms.uint32(2)
)
process.ElElFilter = cms.EDFilter('CandViewCountFilter',
	src	=	cms.InputTag("electronsFilter"),
	checkCharge = 	cms.bool(False),
	minNumber =	cms.uint32(2)
)
process.ElMuFilter = cms.EDProducer('CandViewShallowCloneCombiner',
	decay	=	cms.string("electronsFilter muonsFilter"),
	checkCharge = 	cms.bool(False),
	cut	= 	cms.string("min(daughter(0).pt,daughter(1).pt)>15")
	#minNumber =	cms.uint32(2)
)
# ---- output Module
process.Out = cms.OutputModule("PoolOutputModule",
	fileName = cms.untracked.string ("ZJetsFilter.root"),
	SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring( 'p1', 'p2','p3' )
    )
)

# ---- first sequence: events that pass the trigger ---------------------
process.s1 = cms.Sequence(process.hltFilter * process.muonsFilter * process.MuMuFilter)
process.p1 = cms.Path(process.s1)

process.s2 = cms.Sequence(process.hltFilter * process.electronsFilter * process.ElElFilter)
process.p2 = cms.Path(process.s2)

process.s3 = cms.Sequence(process.hltFilter * process.muonsFilter * process.electronsFilter * process.ElMuFilter)
process.p3 = cms.Path(process.s3)
# ---- end path
process.end = cms.EndPath(process.Out)
# ---- schedule ---------------------------------------------------------
#process.schedule = cms.Schedule(process.p1,process.p2,process.p3)
