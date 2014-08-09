import FWCore.ParameterSet.Config as cms
from CondTools.SiPixel.SiPixelGainCalibrationService_cfi import *

# Clusterizer options
siPixelClusters = cms.EDProducer('RunStepsClusterizer',
    SiPixelGainCalibrationServiceParameters,
    src = cms.InputTag('simSiPixelDigis'),
    algorithm = cms.string('WeightedMeans2D'), # WeightedMeans2D or AdjacentHits or RealClusterizer1D
    clusterSimLink = cms.bool(True)
)



