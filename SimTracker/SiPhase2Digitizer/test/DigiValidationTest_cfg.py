import FWCore.ParameterSet.Config as cms

process = cms.Process("digiTest")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.MessageLogger = cms.Service("MessageLogger",
    debugModules = cms.untracked.vstring('siPixelRawData'),
    destinations = cms.untracked.vstring("cout"),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('ERROR')
    )
)
process.source = cms.Source("PoolSource",
    fileNames =  cms.untracked.vstring(
       #'file:/afs/cern.ch/work/d/dutta/public/Digitizer/CMSSW_6_1_2_SLHC5/src/prod/SingleMu_TK_BE.root'
	'file:/user/aolbrech/TrackTrigger_6XX/CMSSW_6_2_0_SLHC4/src/SuchandraFiles/SingleMu_TK_BE5D_500Events.root'
    )
)
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS161_V15::All', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS162_V1::All', '') #Global tag for 62Y UPgrade Studies (MC) -- SWGuideFrontierConditions Twiki
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3','') #Global tag used in SingleMuonFlat python file

#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('./DigiValidation.root')
)
process.analysis = cms.EDAnalyzer("DigiValidation",
    Verbosity = cms.untracked.bool(False),
    src = cms.InputTag("simSiPixelDigis"),
    simG4 = cms.InputTag("g4SimHits")                              
)
#process.digi_step = cms.Sequence(process.siPixelRawData*process.siPixelDigis)
process.p = cms.Path(process.analysis)
