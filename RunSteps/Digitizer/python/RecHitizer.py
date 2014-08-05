# Imports
import os, sys
import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_phase2_BE5D

# Default parameters
input_file = '../../Output/CLUSTER.root' #os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/CLUSTER.root'
output_file = '../../Output/RECHIT.root' #os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/RECHIT.root'

# Look for updates in the parameters using the program's input
#for i in range(2, len(sys.argv)):
#    if (sys.argv[i] == '_output' and len(sys.argv) > i + 1 and sys.argv[i+1][0] != '_'):
#        output_file = sys.argv[i+1]
#        i += 1
#    elif (sys.argv[i] == '_input' and len(sys.argv) > i + 1 and sys.argv[i+1][0] != '_'):
#        input_file = sys.argv[i+1]
#        i += 1

# Greetings
print '------------------------------------------------------------'
print '-- Running the Clusterizer step with the following arguments:'
print '-- Input file: ' + input_file
print '-- Output file: ' + output_file
print '------------------------------------------------------------'

# Create a new CMS process
process = cms.Process('RecHitizerAnnik')

# Import all the necessary files
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('RunSteps.Clusterizer.Configuration')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')


# Number of events (-1 = all)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

# Input file
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring('file:' + input_file)
)

# Options
process.options = cms.untracked.PSet()

# Metadata (info)
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 0.1 $'),
    annotation = cms.untracked.string('RunSteps RecHitizer'),
    name = cms.untracked.string('Applications')
)

# TAG
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

# Output
process.FEVTDEBUGoutput = cms.OutputModule('PoolOutputModule',
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('file:' + output_file),
)
process.FEVTDEBUGoutput.outputCommands.extend([
  'keep *_siPixelRecHits_*_*'
])

# Steps
process.recHitizer_step = cms.Path(process.siPixelRecHits)

#process.endjob_step = cms.EndPath(process.endOfProcess)

process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Processes to run
process.schedule = cms.Schedule(process.recHitizer_step, process.FEVTDEBUGoutput_step)

process = cust_phase2_BE5D(process)
