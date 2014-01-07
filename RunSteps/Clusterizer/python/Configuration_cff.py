import FWCore.ParameterSet.Config as cms
from CondTools.SiPixel.SiPixelGainCalibrationService_cfi import *

siPixelClusters = cms.EDProducer("SiPhase2Clusterizer",
    SiPixelGainCalibrationServiceParameters,
    src = cms.InputTag("simSiPixelDigis"),
)
