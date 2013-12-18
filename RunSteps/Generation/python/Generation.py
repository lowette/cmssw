import os, sys
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(1000)
)

process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet()

process.configurationMetadata = cms.untracked.PSet(
    	version = cms.untracked.string('$Revision: 1.20 $'),
    	annotation = cms.untracked.string('FourMuPt_1_200_cfi nevts:10'),
    	name = cms.untracked.string('Applications')
)

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    	splitLevel = cms.untracked.int32(0),
    	eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    	outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    	fileName = cms.untracked.string('file:' + os.path.dirname(os.path.realpath(sys.argv[1])) + '/../../Output/GEN_SIM.root'),
    	dataset = cms.untracked.PSet(
        	filterName = cms.untracked.string(''),
        	dataTier = cms.untracked.string('GEN-SIM')
    	),
    	SelectEvents = cms.untracked.PSet(
    		SelectEvents = cms.vstring('generation_step')
    	)
)

process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.generator = cms.EDProducer("FlatRandomPtGunProducer",
    	PGunParameters = cms.PSet(
        	MaxPt = cms.double(200.0),
        	MinPt = cms.double(0.9),
        	PartID = cms.vint32(-13),
       	 	MaxEta = cms.double(2.5),
        	MaxPhi = cms.double(3.14159265359),
        	MinEta = cms.double(-2.5),
        	MinPhi = cms.double(-3.14159265359)
    	),
    	Verbosity = cms.untracked.int32(0),
    	psethack = cms.string('Four mu pt 1 to 200'),
    	AddAntiParticle = cms.bool(True),
    	firstRun = cms.untracked.uint32(1)
)

process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)

for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_phase2_BE5D 

process = cust_phase2_BE5D(process)
