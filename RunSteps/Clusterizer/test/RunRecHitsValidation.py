# Imports
import os, sys
import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_phase2_BE5D

# Default parameters
#input_file = os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/CLURecHit.root'
input_file = os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/CLUSTER.root'
output_file = os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/RecHitValidation.root'

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
print '-- Running the RecHitsValidation step with the following arguments:'
print '-- Input file: ' + input_file
print '-- Output file: ' + output_file
print '------------------------------------------------------------'

# Create a new CMS process
process = cms.Process('recHitTest')

# Import all the necessary files
#process.load('RunSteps.Clusterizer.Configuration')                                #Our personal clusterizer configuration file
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
#process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi')               #Needed to avoid siPixelRecHits attribute error
process.load('RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff')           #Needed for producing these PixelCPE (Cluster Parameter Estimator)

#Suggestions from Ryo  --> Works with these changes!
process.PixelCPEGenericESProducer.Upgrade = cms.bool(True)
process.PixelCPEGenericESProducer.UseErrorsFromTemplates = cms.bool(False)
process.PixelCPEGenericESProducer.LoadTemplatesFromDB = cms.bool(False)
process.PixelCPEGenericESProducer.TruncatePixelCharge = cms.bool(False)
process.PixelCPEGenericESProducer.IrradiationBiasCorrection = cms.bool(False)
process.PixelCPEGenericESProducer.DoCosmics = cms.bool(False)
#process.siPixelRecHits.CPE = cms.string('PixelCPEGeneric')                        #Needed??

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

# Output
process.FEVTDEBUGoutput = cms.OutputModule('PoolOutputModule',
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('file:' + output_file),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEM-SIM-DIGI-CLU-REC')
    )
)

process.FEVTDEBUGoutput.outputCommands.extend([
  'keep *_siPixelRecHits_*_*'
])


# DEBUG
process.MessageLogger = cms.Service('MessageLogger',
	debugModules = cms.untracked.vstring('siPixelClusters'),
	destinations = cms.untracked.vstring('cout'),
	cout = cms.untracked.PSet(
		threshold = cms.untracked.string('ERROR')
	)
)

# Analyzer
process.analysis = cms.EDAnalyzer('RunStepsRecHitsValidation',
    rechit = cms.InputTag('siPixelRecHits'),
    clu = cms.InputTag('siPixelClusters')
)

process.clusterizer_step = cms.Path(cms.Sequence(process.siPixelRecHits*process.analysis));

process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Processes to run
process.schedule = cms.Schedule(process.clusterizer_step)

# Processes to run
#process.p = cms.Path(process.analysis)
