import FWCore.ParameterSet.Config as cms

process = cms.Process("cluTest")

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(-1)
)

process.MessageLogger = cms.Service("MessageLogger",
    	debugModules = cms.untracked.vstring('siPixelClusters'),
    	destinations = cms.untracked.vstring('cout'),
	cout = cms.untracked.PSet(
    		threshold = cms.untracked.string('ERROR')
    	)
)

process.source = cms.Source("PoolSource",
  	fileNames = cms.untracked.vstring(    
		'file:../../Output/CLUSTER.root',
  	)
)

process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_P_V28::All"

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('../../Output/ClusterValidation.root')
)

process.analysis = cms.EDAnalyzer("ClusterValidationTest",
    Verbosity = cms.untracked.bool(True),
    src = cms.InputTag("siPixelClusters"),
)

process.p = cms.Path(process.analysis)
