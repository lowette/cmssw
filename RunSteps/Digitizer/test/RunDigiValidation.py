import os, sys
import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("digiTest")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

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
    	fileNames = cms.untracked.vstring('file:' + os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/DIGI.root')
)

process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS161_V15::All', '')

process.TFileService = cms.Service("TFileService",
    	fileName = cms.string(os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/DigiValidation.root')
)

process.analysis = cms.EDAnalyzer("DigiValidationTest",
    	Verbosity = cms.untracked.bool(False),
   	src = cms.InputTag("simSiPixelDigis"),
	simG4 = cms.InputTag("g4SimHits")
)

process.p = cms.Path(process.analysis)
