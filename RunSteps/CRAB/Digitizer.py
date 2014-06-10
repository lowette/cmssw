# Imports
import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_phase2_BE5D

# Create a new CMS process
process = cms.Process('DIGI')

# Import all the necessary files
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('CRAB.Configuration')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load('SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Number of events (-1 = all)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input file
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('/store/relval/CMSSW_6_2_0_SLHC7/RelValFourMuPt1_200_UPG2023_BE5D/GEN-SIM/DES19_62_V8_UPG2023-v2/00000/EC4C06C3-2890-E311-BD11-003048FEB966.root')
   #fileNames = cms.untracked.vstring('/RelValFourMuPt1_200_UPG2023_BE5D/CMSSW_6_2_0_SLHC7-DES19_62_V8_UPG2023-v2/GEN-SIM')/RelValFourMuPt1_200_UPG2023_BE5D/CMSSW_6_2_0_SLHC7-DES19_62_V8_UPG2023-v2/GEN-SIM')
)

# Options
process.options = cms.untracked.PSet()

# Metadata (info)
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 0.1 $'),
    annotation = cms.untracked.string('RunSteps Digitizer'),
    name = cms.untracked.string('Applications')
)

# TAG
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

# Output
process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('DigiPU.root'),
    dataset = cms.untracked.PSet(        
	filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
    )
)

# Steps
process.digitisation_step = cms.Path(process.pdigi)

process.L1simulation_step = cms.Path(process.SimL1Emulator)

process.digi2raw_step = cms.Path(process.DigiToRaw)

process.L1TrackTrigger_step = cms.Path(process.TrackTriggerClustersStubs)

process.L1TTAssociator_step = cms.Path(process.TrackTriggerAssociatorClustersStubs)

process.endjob_step = cms.EndPath(process.endOfProcess)

process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Processes to run
process.schedule = cms.Schedule(process.digitisation_step, process.L1simulation_step, process.digi2raw_step, process.L1TrackTrigger_step, process.L1TTAssociator_step, process.endjob_step, process.FEVTDEBUGoutput_step)

process = cust_phase2_BE5D(process)
