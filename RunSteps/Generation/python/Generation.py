# Imports
import os, sys
import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_phase2_BE5D

# Default parameters
n = 10
output_file = os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/GEN_SIM.root'

# Look for updates in the parameters using the program's input
for i in range(2, len(sys.argv)):
    if (sys.argv[i] == '_n' and len(sys.argv) > i + 1 and sys.argv[i+1][0] != '_'):
        n = int(sys.argv[i+1])
        i += 1
    if (sys.argv[i] == '_output' and len(sys.argv) > i + 1 and sys.argv[i+1][0] != '_'):
        output_file = sys.argv[i+1]
        i += 1

# Greetings
print '------------------------------------------------------------'
print '-- Running the Generation step with the following arguments:'
print '-- Number of events: ' + str(n)
print '-- Output file: ' + output_file
print '------------------------------------------------------------'

# Create a new CMS process
process = cms.Process('SIM')

# Import all the necessary files
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Number of events to generate
process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(n)
)

# Input file (none as it will be created)
process.source = cms.Source('EmptySource')

# Options
process.options = cms.untracked.PSet()

# Metadata (info)
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 0.1 $'),
    annotation = cms.untracked.string('RunSteps Generator'),
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
	dataset = cms.untracked.PSet(
    	filterName = cms.untracked.string(''),
    	dataTier = cms.untracked.string('GEN-SIM')
	),
	SelectEvents = cms.untracked.PSet(
		SelectEvents = cms.vstring('generation_step')
	)
)

# Generator options
process.genstepfilter.triggerConditions=cms.vstring('generation_step')

process.generator = cms.EDProducer('FlatRandomPtGunProducer',
    PGunParameters = cms.PSet(
        MaxPt = cms.double(200.),
        MinPt = cms.double(1.),
        PartID = cms.vint32(-13),
        MaxEta = cms.double(2.5),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(-2.5),
        MinPhi = cms.double(-3.14159265359)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('One mu pt 1 to 200'),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
)

# Steps
process.generation_step = cms.Path(process.pgen)

process.simulation_step = cms.Path(process.psim)

process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)

process.endjob_step = cms.EndPath(process.endOfProcess)

process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Processes to run
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)

for path in process.paths:
	getattr(process, path)._seq = process.generator * getattr(process, path)._seq

process = cust_phase2_BE5D(process)
