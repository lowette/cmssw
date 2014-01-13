import os, sys
import FWCore.ParameterSet.Config as cms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_phase2_BE5D 
from Configuration.AlCa.GlobalTag import GlobalTag

input_file = os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/DIGI.root'
output_file = os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/L1.root'

for i in range(2, len(sys.argv)):
        if (sys.argv[i] == '_output' and len(sys.argv) > i + 1 and sys.argv[i+1][0] != '_'):
                output_file = sys.argv[i+1]
                i += 1
        elif (sys.argv[i] == '_input' and len(sys.argv) > i + 1 and sys.argv[i+1][0] != '_'):
                input_file = sys.argv[i+1]
                i += 1
                
print 'Input file: ' + input_file
print 'Output file: ' + output_file

process = cms.Process('RECO')

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:' + input_file)
)

process.options = cms.untracked.PSet()

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 0.1 $'),
    annotation = cms.untracked.string('RunSteps RECO'),
    name = cms.untracked.string('Applications')
)

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('file:' + output_file),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
    )
)

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.L1simulation_step = cms.Path(process.SimL1Emulator)

process.digi2raw_step = cms.Path(process.DigiToRaw)

process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)

process.endjob_step = cms.EndPath(process.endOfProcess)

process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

process.schedule = cms.Schedule(process.L1simulation_step, process.digi2raw_step, process.L1TrackTrigger_step, process.endjob_step, process.FEVTDEBUGoutput_step)

process = cust_phase2_BE5D(process)

