import FWCore.ParameterSet.Config as cms
from CondTools.SiPixel.SiPixelGainCalibrationService_cfi import *

# Clusterizer options
siPixelClusters = cms.EDProducer('RunStepsClusterizer',
    SiPixelGainCalibrationServiceParameters,
    src = cms.InputTag('simSiPixelDigis'),
    DigiSearchWidthX = cms.int32(2),
    DigiSearchWidthY = cms.int32(2),
    MinimumWidthX = cms.int32(1),
    MaximumWidthX = cms.int32(5),
    MinimumWidthY = cms.int32(1),
    MaximumWidthY = cms.int32(5),
    MinimumNumberOfDigis = cms.int32(0),
    MaximumNumberOfDigis = cms.int32(25),
    SplitIfOverThresholds = cms.bool(True)
)
