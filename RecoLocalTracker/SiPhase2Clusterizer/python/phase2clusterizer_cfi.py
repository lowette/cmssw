import FWCore.ParameterSet.Config as cms

# Clusterizer options
siPhase2Clusters = cms.EDProducer('SiPhase2Clusterizer',
    src = cms.InputTag('simSiPixelDigis'),
    clusterSimLink = cms.bool(True)
)



