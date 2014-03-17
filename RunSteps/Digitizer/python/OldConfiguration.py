# Imports
import FWCore.ParameterSet.Config as cms
from SimGeneral.MixingModule.mixObjects_cfi import theMixObjects
from SimGeneral.MixingModule.mixPoolSource_cfi import *
from SimGeneral.MixingModule.aliases_cfi import *
from SimGeneral.MixingModule.pixelDigitizer_cfi import *
from SimGeneral.MixingModule.stripDigitizer_cfi import *
from SimGeneral.MixingModule.ecalDigitizer_cfi import *
from SimGeneral.MixingModule.hcalDigitizer_cfi import *
from SimGeneral.MixingModule.castorDigitizer_cfi import *
from SimGeneral.MixingModule.trackingTruthProducer_cfi import *
from SimCalorimetry.Configuration.SimCalorimetry_cff import *
from SimMuon.Configuration.SimMuon_cff import *
from SimGeneral.PileupInformation.AddPileupSummary_cfi import *

# Digitizers
doAllDigi = cms.Sequence(calDigi + muonDigi)
pdigi = cms.Sequence(cms.SequencePlaceholder('randomEngineStateProducer') * cms.SequencePlaceholder('mix') * doAllDigi * addPileupInfo)

# Digitizers for the mixing
theDigitizers = cms.PSet(
    pixel = cms.PSet(pixelDigitizer),
    strip = cms.PSet(stripDigitizer),
    ecal = cms.PSet(ecalDigitizer),
    hcal = cms.PSet(hcalDigitizer),
    castor  = cms.PSet(castorDigitizer),
    mergedtruth = cms.PSet(trackingParticles)
)

# Mixing module
mix = cms.EDProducer("MixingModule",
    digitizers = cms.PSet(theDigitizers),
    LabelPlayback = cms.string(''),
    maxBunch = cms.int32(3),
    minBunch = cms.int32(-5),
    bunchspace = cms.int32(25),
    mixProdStep1 = cms.bool(False),
    mixProdStep2 = cms.bool(False),
    playback = cms.untracked.bool(False),
    useCurrentProcessOnly = cms.bool(False),
    mixObjects = cms.PSet(theMixObjects)
)
