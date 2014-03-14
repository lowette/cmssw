# Imports
import os, sys
import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag

# Default parameters
input_file = os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/DIGI.root'
output_file = os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/DigiValidation.root'

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
    input = cms.untracked.int32(-1)
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

# Analyzer
process.analysis = cms.EDAnalyzer('RunStepsDigiValidation',
    Verbosity = cms.untracked.bool(False),
   	src = cms.InputTag('simSiPixelDigis'),
	simG4 = cms.InputTag('g4SimHits')
)

# Processes to run
process.p = cms.Path(process.analysis)
