# Imports
import FWCore.ParameterSet.Config as cms
from SimGeneral.MixingModule.mixObjects_cfi import theMixObjects
from SimGeneral.MixingModule.mixPoolSource_cfi import *
from SimGeneral.MixingModule.aliases_cfi import *
from RunSteps.Digitizer.RunStepsDigitizer_cfi import *
from SimGeneral.MixingModule.stripDigitizer_cfi import *
from SimGeneral.MixingModule.ecalDigitizer_cfi import *
from SimGeneral.MixingModule.hcalDigitizer_cfi import *
from SimGeneral.MixingModule.castorDigitizer_cfi import *
from SimGeneral.MixingModule.trackingTruthProducer_cfi import *

# Digitizers to use
theDigitizers = cms.PSet(
    pixel = cms.PSet(runStepsDigitizer_cfi),
    strip = cms.PSet(stripDigitizer),
    ecal = cms.PSet(ecalDigitizer),
    hcal = cms.PSet(hcalDigitizer),
    castor  = cms.PSet(castorDigitizer),
    mergedtruth = cms.PSet(trackingParticles)
)

theDigitizersValid = cms.PSet(
    pixel = cms.PSet(runStepsDigitizer_cfi),
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
    minBunch = cms.int32(-5), ## in terms of 25 ns
    bunchspace = cms.int32(25),
    mixProdStep1 = cms.bool(False),
    mixProdStep2 = cms.bool(False),
    playback = cms.untracked.bool(False),
    useCurrentProcessOnly = cms.bool(False),
    mixObjects = cms.PSet(theMixObjects)
)
