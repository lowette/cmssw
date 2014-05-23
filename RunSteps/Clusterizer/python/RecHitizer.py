# Imports
import os, sys
import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_phase2_BE5D

# Default parameters
input_file = os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/CLUSTER.root'
output_file = os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/CLURechHit.root'

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
print '-- Running the Clusterizer step with the following arguments:'
print '-- Input file: ' + input_file
print '-- Output file: ' + output_file
print '------------------------------------------------------------'

# Create a new CMS process
process = cms.Process('RecHitAnnik')

# Import all the necessary files
process.load('RunSteps.Clusterizer.Configuration')                                #Our personal clusterizer configuration file
process.load('Configuration.EventContent.EventContent_cff')                       #Needed to avoid FEVTDEBUGEventContent error
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')       #Needed to avoid CSCGeometryESModule error
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')     #Without this, no siPixelRecHits produced
#process.load('Configuration.StandardSequences.EndOfProcess_cff')                 #Not needed when the step is removed from the sequence
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

# Number of events (-1 = all)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input file
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring('file:' + input_file)
)

# Options
process.options = cms.untracked.PSet()

# Metadata (info)
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 0.5 $'),
    annotation = cms.untracked.string('RunSteps Clusterizer'),
    name = cms.untracked.string('Applications')
)

#Full trig and Time report
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)  #Set this to true for full Trig and TimeReport
)

# TAG
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

# Output
process.FEVTDEBUGoutput = cms.OutputModule('PoolOutputModule',
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('file:' + output_file),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEM-SIM-DIGI-CLU')
    )
)
process.FEVTDEBUGoutput.outputCommands.extend([
  'keep *_siPixelRecHits_*_*'
]) 

# Steps
process.clusterizer_step = cms.Path(cms.Sequence(process.siPixelClusters*process.siPixelRecHits));
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Processes to run
process.schedule = cms.Schedule(process.clusterizer_step, process.FEVTDEBUGoutput_step)

process = cust_phase2_BE5D(process)
