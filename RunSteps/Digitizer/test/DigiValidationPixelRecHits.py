# Imports
import os, sys
import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag

# Default parameters
input_file = 'DIGI.root'
output_file = 'DigiValidationRecHits.root'

# Look for updates in the parameters using the program's input
for i in range(2, len(sys.argv)):
    if (sys.argv[i] == '_output' and len(sys.argv) > i + 1 and sys.argv[i+1][0] != '_'):
        output_file = sys.argv[i+1]
        i += 1
    elif (sys.argv[i] == '_input' and len(sys.argv) > i + 1 and sys.argv[i+1][0] != '_'):
        input_file = sys.argv[i+1]
        i += 1

# Greetings
print '------------------------------------------------------------'
print '-- Running the DigiValidation step with the following arguments:'
print '-- Input file: ' + input_file
print '-- Output file: ' + output_file
print '------------------------------------------------------------'

# Create a new CMS process
process = cms.Process('digiTest')

# Import all the necessary files
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Number of events (-1 = all)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input file
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring('file:' + input_file)
)

# TAG
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('file:' + output_file)
)

# DEBUG
process.MessageLogger = cms.Service('MessageLogger',
    	debugModules = cms.untracked.vstring('siPixelRawData'),
    	destinations = cms.untracked.vstring('cout'),
   	cout = cms.untracked.PSet(
        	threshold = cms.untracked.string('ERROR')
    	)
)

# Producers for recHits
process.siPixelClusters = cms.EDProducer("SiPixelClusterProducer",
    src = cms.InputTag("siPixelDigis"),
    ChannelThreshold = cms.int32(1000),
    maxNumberOfClusters = cms.int32(-1),
    SplitClusters = cms.bool(False),
    MissCalibrate = cms.untracked.bool(True),
    VCaltoElectronGain = cms.int32(65),
    VCaltoElectronOffset = cms.int32(-414),
    payloadType = cms.string('Offline'),
    SeedThreshold = cms.int32(1000),
    ClusterThreshold = cms.double(4000.0)
)

process.siPixelDigis = cms.EDProducer("SiPixelRawToDigi",
    Timing = cms.untracked.bool(False),
    IncludeErrors = cms.bool(False),
    UseCablingTree = cms.untracked.bool(True),
    InputLabel = cms.InputTag("source"),
    CheckPixelOrder = cms.bool(False)
)

process.siPixelRecHits= cms.EDProducer("SiPixelRecHitConverter",
    VerboseLevel = cms.untracked.int32(0),
    src = cms.InputTag("siPixelClusters"),
    CPE = cms.string('PixelCPEGeneric')
)

# Analyzer
process.analysis = cms.EDAnalyzer('RunStepsDigiValidation',
    Verbosity = cms.untracked.bool(False),
   	src = cms.InputTag('simSiPixelDigis'),
	simG4 = cms.InputTag('g4SimHits')
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
    ,outputCommands = cms.untracked.vstring('drop *',
      "keep *_simSiPixelDigis_*_*",
      "keep *_siPixelClusters_*_*",
      "keep *_siPixelRecHits_*_*")      
)

#RecoLocalTrackerRECO = cms.PSet(
  #  outputCommands = cms.untracked.vstring(
    #'keep DetIdedmEDCollection_siStripDigis_*_*',
#    'keep DetIdedmEDCollection_siPixelDigis_*_*',
  #  'keep *_siPixelClusters_*_*', 
#    'keep *_siStripClusters_*_*',
#    'keep *_clusterSummaryProducer_*_*')
#)

# Processes to run
process.pixeltrackerlocalreco = cms.Sequence(process.siPixelClusters*process.siPixelRecHits)
process.trackerlocalreco = cms.Sequence(process.pixeltrackerlocalreco)

process.p = cms.Path(process.trackerlocalreco)
process.e = cms.EndPath(process.out)
#process.p = cms.Path(process.analysis)
