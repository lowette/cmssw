import os, sys
import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag

input_file = os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/CLUSTER.root'
output_file = os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/ClusterValidation.root'

for i in range(2, len(sys.argv)):
        if (sys.argv[i] == '_output' and len(sys.argv) > i + 1 and sys.argv[i+1][0] != '_'):
                output_file = sys.argv[i+1]
                i += 1
        elif (sys.argv[i] == '_input' and len(sys.argv) > i + 1 and sys.argv[i+1][0] != '_'):
                input_file = sys.argv[i+1]
                i += 1
                
print 'Input file: ' + input_file
print 'Output file: ' + output_file

process = cms.Process("cluTest")

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
	debugModules = cms.untracked.vstring('siPixelClusters'),
	destinations = cms.untracked.vstring("cout"),
	cout = cms.untracked.PSet(
		threshold = cms.untracked.string('ERROR')
	)
)

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring('file:' + input_file)
)

process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS161_V15::All', '')

process.TFileService = cms.Service("TFileService",
            fileName = cms.string(output_file)
)

process.analysis = cms.EDAnalyzer("ClusterValidationTest",
		Verbosity = cms.untracked.bool(False),
		src = cms.InputTag("siPixelClusters")
)

process.p = cms.Path(process.analysis)

