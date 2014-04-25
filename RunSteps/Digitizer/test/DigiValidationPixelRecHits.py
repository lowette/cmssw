# Imports
import os, sys
import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag

# Default parameters
input_file = 'DIGI.root'
output_file = 'DigiValidationRecHits.root'

# Look for updates in the parameters using the program's input
for i in range(2, len(sys.argv)):
    if (sys.argv[i] == '_output' and len(sys.argv) > i + 1 and sys.argv[i+1][0] != '_'):
        output_file = sys.argv[i+1]
        i += 1
    elif (sys.argv[i] == '_input' and len(sys.argv) > i + 1 and sys.argv[i+1][0] != '_'):
        input_file = sys.argv[i+1]
        i += 1

# Greetings
print '------------------------------------------------------------'
print '-- Running the DigiValidation step with the following arguments:'
print '-- Input file: ' + input_file
print '-- Output file: ' + output_file
print '------------------------------------------------------------'

# Create a new CMS process
process = cms.Process('digiTest')

# Import all the necessary files
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Number of events (-1 = all)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input file
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring('file:' + input_file)
)

# TAG
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
#process.GlobalTag.globaltag = cms.string('DES23_62_V1::All')

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('file:' + output_file)
)

# DEBUG
process.MessageLogger = cms.Service('MessageLogger',
    	debugModules = cms.untracked.vstring('siPixelRawData'),
    	destinations = cms.untracked.vstring('cout'),
   	cout = cms.untracked.PSet(
        	threshold = cms.untracked.string('ERROR')
    	)
)

# Producers for recHits
process.siPixelClusters = cms.EDProducer("SiPixelClusterProducer",
    src = cms.InputTag("simSiPixelDigis"),
    ChannelThreshold = cms.int32(1000),
    maxNumberOfClusters = cms.int32(-1),
    SplitClusters = cms.bool(False),
    MissCalibrate = cms.untracked.bool(False),
    VCaltoElectronGain = cms.int32(65),
    VCaltoElectronOffset = cms.int32(-414),
    payloadType = cms.string('Offline'),
    SeedThreshold = cms.int32(1000),
    ClusterThreshold = cms.double(4000.0),
	EnforcePairedClusters=cms.bool(False)	
)

#process.siPixelDigis = cms.EDProducer("SiPixelRawToDigi",
#    Timing = cms.untracked.bool(False),
#    IncludeErrors = cms.bool(False),
#    UseCablingTree = cms.untracked.bool(True),
#    InputLabel = cms.InputTag("source"),
#    CheckPixelOrder = cms.bool(False)
#)

process.siPixelRecHits= cms.EDProducer("SiPixelRecHitConverter",
    VerboseLevel = cms.untracked.int32(0),
    src = cms.InputTag("siPixelClusters"),
    CPE = cms.string('PixelCPEGeneric'),
	EnforcePairedClusters=cms.bool(False)
)

# Analyzer
process.analysis = cms.EDAnalyzer('RunStepsDigiValidation',
    Verbosity = cms.untracked.bool(False),
   	src = cms.InputTag('simSiPixelDigis'),
	simG4 = cms.InputTag('g4SimHits')
)

#process.FEVTDEBUGoutput.outputCommands.extend([
#  'keep *_siPixelRecHits_*_*'
#])
 
process.out = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('file:myOutputFile.root'),
    outputCommands = cms.untracked.vstring( ('drop *',
        'drop *',
        'drop *',
        'keep  FEDRawDataCollection_rawDataCollector_*_*',
        'keep  FEDRawDataCollection_source_*_*',
        'keep  FEDRawDataCollection_rawDataCollector_*_*',
        'keep  FEDRawDataCollection_source_*_*',
        'drop *_hlt*_*_*',
        'keep *_hltL1GtObjectMap_*_*',
        'keep FEDRawDataCollection_rawDataCollector_*_*',
        'keep FEDRawDataCollection_source_*_*',
        'keep edmTriggerResults_*_*_*',
        'keep triggerTriggerEvent_*_*_*',
        'keep *_g4SimHits_*_*',
        'keep edmHepMCProduct_source_*_*',
        'keep *_allTrackMCMatch_*_*',
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*',
        'keep CSCDetIdCSCComparatorDigiMuonDigiCollection_simMuonCSCDigis_*_*',
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*',
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*',
        'keep EBSrFlagsSorted_simEcalDigis_*_*',
        'keep EESrFlagsSorted_simEcalDigis_*_*',
        'keep CrossingFramePlaybackInfoExtended_*_*_*',
        'keep PileupSummaryInfos_*_*_*',
        'keep LHERunInfoProduct_*_*_*',
        'keep LHEEventProduct_*_*_*',
        'keep GenRunInfoProduct_generator_*_*',
        'keep GenEventInfoProduct_generator_*_*',
        'keep edmHepMCProduct_generator_*_*',
        'keep GenFilterInfo_*_*_*',
        'keep *_genParticles_*_*',
        'keep recoGenJets_*_*_*',
        'keep *_genParticle_*_*',
        'keep recoGenMETs_*_*_*',
        'keep FEDRawDataCollection_source_*_*',
        'keep FEDRawDataCollection_rawDataCollector_*_*',
        'keep *_MEtoEDMConverter_*_*',
        'keep *_randomEngineStateProducer_*_*',
        'keep DetIdedmEDCollection_siStripDigis_*_*',
        'keep DetIdedmEDCollection_siPixelDigis_*_*',
        'keep *_siPixelClusters_*_*',
        'keep *_siStripClusters_*_*',
        'keep *_clusterSummaryProducer_*_*',
        'keep *_dt1DRecHits_*_*',
        'keep *_dt4DSegments_*_*',
        'keep *_dt1DCosmicRecHits_*_*',
        'keep *_dt4DCosmicSegments_*_*',
        'keep *_csc2DRecHits_*_*',
        'keep *_cscSegments_*_*',
        'keep *_rpcRecHits_*_*',
        'keep *_hbhereco_*_*',
        'keep *_hfreco_*_*',
        'keep *_horeco_*_*',
        'keep *_hbheUpgradeReco_*_*',
        'keep *_hfUpgradeReco_*_*',
        'keep HBHERecHitsSorted_hbherecoMB_*_*',
        'keep HORecHitsSorted_horecoMB_*_*',
        'keep HFRecHitsSorted_hfrecoMB_*_*',
        'keep ZDCDataFramesSorted_*Digis_*_*',
        'keep ZDCRecHitsSorted_*_*_*',
        'keep *_reducedHcalRecHits_*_*',
        'keep *_castorreco_*_*',
        'keep HcalUnpackerReport_*_*_*',
        'keep *_ecalPreshowerRecHit_*_*',
        'keep *_ecalRecHit_*_*',
        'keep *_ecalCompactTrigPrim_*_*',
        'keep *_ecalTPSkim_*_*',
        'keep *_selectDigi_*_*',
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*',
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*',
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*',
        'keep *_hybridSuperClusters_*_*',
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*',
        'keep *_multi5x5SuperClusters_*_*',
        'keep recoSuperClusters_multi5x5SuperClusters_*_*',
        'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*',
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*',
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*',
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*',
        'keep *_particleFlowSuperClusterECAL_*_*',
        'drop recoClusterShapes_*_*_*',
        'drop recoBasicClustersToOnerecoClusterShapesAssociation_*_*_*',
        'drop recoBasicClusters_multi5x5BasicClusters_multi5x5BarrelBasicClusters_*',
        'drop recoSuperClusters_multi5x5SuperClusters_multi5x5BarrelSuperClusters_*',
        'keep *_CkfElectronCandidates_*_*',
        'keep *_GsfGlobalElectronTest_*_*',
        'keep *_electronMergedSeeds_*_*',
        'keep recoGsfTracks_electronGsfTracks_*_*',
        'keep recoGsfTrackExtras_electronGsfTracks_*_*',
        'keep recoTrackExtras_electronGsfTracks_*_*',
        'keep TrackingRecHitsOwned_electronGsfTracks_*_*',
        'keep recoTracks_generalTracks_*_*',
        'keep recoTrackExtras_generalTracks_*_*',
        'keep TrackingRecHitsOwned_extraFromSeeds_*_*',
        'keep uints_extraFromSeeds_*_*',
        'keep TrackingRecHitsOwned_generalTracks_*_*',
        'keep recoTracks_beamhaloTracks_*_*',
        'keep recoTrackExtras_beamhaloTracks_*_*',
        'keep TrackingRecHitsOwned_beamhaloTracks_*_*',
        'keep recoTracks_regionalCosmicTracks_*_*',
        'keep recoTrackExtras_regionalCosmicTracks_*_*',
        'keep TrackingRecHitsOwned_regionalCosmicTracks_*_*',
        'keep recoTracks_rsWithMaterialTracks_*_*',
        'keep recoTrackExtras_rsWithMaterialTracks_*_*',
        'keep TrackingRecHitsOwned_rsWithMaterialTracks_*_*',
        'keep recoTracks_conversionStepTracks_*_*',
        'keep recoTrackExtras_conversionStepTracks_*_*',
        'keep TrackingRecHitsOwned_conversionStepTracks_*_*',
        'keep *_ctfPixelLess_*_*',
        'keep *_dedxTruncated40_*_*',
        'keep *_dedxDiscrimASmi_*_*',
        'keep *_dedxHarmonic2_*_*',
        'keep *_trackExtrapolator_*_*',
        'keep floatedmValueMap_generalTracks_*_*',
        'keep *_kt6CaloJets_*_*',
        'keep *_ak5CaloJets_*_*',
        'keep double_kt6PFJets_rho_*',
        'keep *_ak5PFJets_*_*',
        'keep *_ak8PFJets_*_*',
        'keep *_ak5PFJetsCHS_*_*',
        'keep *_ak8PFJetsCHS_*_*',
        'keep *_ca8PFJetsCHS_*_*',
        'keep *_ca8PFJetsCHSPruned_*_*',
        'keep *_cmsTopTagPFJetsCHS_*_*',
        'keep *_JetPlusTrackZSPCorJetAntiKt5_*_*',
        'keep *_ak5TrackJets_*_*',
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*',
        'keep *_caloTowers_*_*',
        'keep *_towerMaker_*_*',
        'keep *_CastorTowerReco_*_*',
        'keep *_kt4JetTracksAssociatorAtVertex_*_*',
        'keep *_kt4JetTracksAssociatorAtCaloFace_*_*',
        'keep *_kt4JetExtender_*_*',
        'keep *_ak5JetTracksAssociatorAtVertex_*_*',
        'keep *_ak5JetTracksAssociatorAtVertexPF_*_*',
        'keep *_ak5JetTracksAssociatorAtCaloFace_*_*',
        'keep *_ak5JetTracksAssociatorExplicit_*_*',
        'keep *_ak5JetExtender_*_*',
        'keep *_ak7JetTracksAssociatorAtVertex*_*_*',
        'keep *_ak7JetTracksAssociatorAtCaloFace*_*_*',
        'keep *_ak7JetExtender_*_*',
        'keep *_ak5JetID_*_*',
        'keep *_ak7JetID_*_*',
        'keep *_ic5JetID_*_*',
        'keep *_kt4JetID_*_*',
        'keep *_kt6JetID_*_*',
        'keep *_ak7BasicJets_*_*',
        'keep *_ak7CastorJetID_*_*',
        'keep double_kt6CaloJetsCentral_*_*',
        'keep double_kt6PFJetsCentralChargedPileUp_*_*',
        'keep double_kt6PFJetsCentralNeutral_*_*',
        'keep double_kt6PFJetsCentralNeutralTight_*_*',
        'keep *_fixedGridRho*_*_*',
        'keep recoCaloMETs_met_*_*',
        'keep recoCaloMETs_metNoHF_*_*',
        'keep recoCaloMETs_metHO_*_*',
        'keep recoCaloMETs_corMetGlobalMuons_*_*',
        'keep recoMETs_tcMet_*_*',
        'keep recoMETs_tcMetWithPFclusters_*_*',
        'keep recoPFMETs_pfMet_*_*',
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*',
        'keep recoMuonMETCorrectionDataedmValueMap_muonTCMETValueMapProducer_*_*',
        'keep recoHcalNoiseRBXs_hcalnoise_*_*',
        'keep HcalNoiseSummary_hcalnoise_*_*',
        'keep *HaloData_*_*_*',
        'keep *BeamHaloSummary_BeamHaloSummary_*_*',
        'keep *_MuonSeed_*_*',
        'keep *_ancientMuonSeed_*_*',
        'keep *_mergedStandAloneMuonSeeds_*_*',
        'keep TrackingRecHitsOwned_globalMuons_*_*',
        'keep TrackingRecHitsOwned_tevMuons_*_*',
        'keep recoCaloMuons_calomuons_*_*',
        'keep *_CosmicMuonSeed_*_*',
        'keep recoTrackExtras_cosmicMuons_*_*',
        'keep TrackingRecHitsOwned_cosmicMuons_*_*',
        'keep recoTrackExtras_globalCosmicMuons_*_*',
        'keep TrackingRecHitsOwned_globalCosmicMuons_*_*',
        'keep recoTrackExtras_cosmicMuons1Leg_*_*',
        'keep TrackingRecHitsOwned_cosmicMuons1Leg_*_*',
        'keep recoTrackExtras_globalCosmicMuons1Leg_*_*',
        'keep TrackingRecHitsOwned_globalCosmicMuons1Leg_*_*',
        'keep recoTracks_cosmicsVetoTracks_*_*',
        'keep *_SETMuonSeed_*_*',
        'keep recoTracks_standAloneSETMuons_*_*',
        'keep recoTrackExtras_standAloneSETMuons_*_*',
        'keep TrackingRecHitsOwned_standAloneSETMuons_*_*',
        'keep recoTracks_globalSETMuons_*_*',
        'keep recoTrackExtras_globalSETMuons_*_*',
        'keep TrackingRecHitsOwned_globalSETMuons_*_*',
        'keep recoMuons_muonsWithSET_*_*',
        'keep *_muons_*_*',
        'keep *_*_muons_*',
        'drop *_muons_muons1stStep2muonsMap_*',
        'drop recoIsoDepositedmValueMap_muons_*_*',
        'drop doubleedmValueMap_muons_muPFIso*_*',
        'keep recoTracks_standAloneMuons_*_*',
        'keep recoTrackExtras_standAloneMuons_*_*',
        'keep TrackingRecHitsOwned_standAloneMuons_*_*',
        'keep recoTracks_globalMuons_*_*',
        'keep recoTrackExtras_globalMuons_*_*',
        'keep recoTracks_tevMuons_*_*',
        'keep recoTrackExtras_tevMuons_*_*',
        'keep recoTracks_generalTracks_*_*',
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*',
        'keep recoTracks_cosmicMuons_*_*',
        'keep recoTracks_globalCosmicMuons_*_*',
        'keep recoMuons_muonsFromCosmics_*_*',
        'keep recoTracks_cosmicMuons1Leg_*_*',
        'keep recoTracks_globalCosmicMuons1Leg_*_*',
        'keep recoMuons_muonsFromCosmics1Leg_*_*',
        'keep recoTracks_refittedStandAloneMuons_*_*',
        'keep recoTrackExtras_refittedStandAloneMuons_*_*',
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*',
        'keep *_muIsoDepositTk_*_*',
        'keep *_muIsoDepositCalByAssociatorTowers_*_*',
        'keep *_muIsoDepositCalByAssociatorHits_*_*',
        'keep *_muIsoDepositJets_*_*',
        'keep *_muGlobalIsoDepositCtfTk_*_*',
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*',
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*',
        'keep *_muGlobalIsoDepositJets_*_*',
        'keep *_impactParameterTagInfos_*_*',
        'keep *_trackCountingHighEffBJetTags_*_*',
        'keep *_trackCountingHighPurBJetTags_*_*',
        'keep *_jetProbabilityBJetTags_*_*',
        'keep *_jetBProbabilityBJetTags_*_*',
        'keep *_secondaryVertexTagInfos_*_*',
        'keep *_ghostTrackVertexTagInfos_*_*',
        'keep *_simpleSecondaryVertexBJetTags_*_*',
        'keep *_simpleSecondaryVertexHighEffBJetTags_*_*',
        'keep *_simpleSecondaryVertexHighPurBJetTags_*_*',
        'keep *_combinedSecondaryVertexBJetTags_*_*',
        'keep *_combinedSecondaryVertexMVABJetTags_*_*',
        'keep *_ghostTrackBJetTags_*_*',
        'keep *_softPFMuonsTagInfos_*_*',
        'keep *_softPFElectronsTagInfos_*_*',
        'keep *_softPFElectronBJetTags_*_*',
        'keep *_softPFMuonBJetTags_*_*',
        'keep *_softMuonTagInfos_*_*',
        'keep *_softMuonBJetTags_*_*',
        'keep *_softMuonByIP3dBJetTags_*_*',
        'keep *_softMuonByPtBJetTags_*_*',
        'keep *_ak5PFJetsRecoTauPiZeros_*_*',
        'keep *_hpsPFTauProducer_*_*',
        'keep *_hpsPFTauDiscrimination*_*_*',
        'keep  *_offlinePrimaryVertices__*',
        'keep  *_offlinePrimaryVerticesWithBS_*_*',
        'keep  *_offlinePrimaryVerticesFromCosmicTracks_*_*',
        'keep  *_nuclearInteractionMaker_*_*',
        'keep *_generalV0Candidates_*_*',
        'keep recoGsfElectronCores_gsfElectronCores_*_*',
        'keep recoGsfElectrons_gsfElectrons_*_*',
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*',
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*',
        'keep floatedmValueMap_eidRobustLoose_*_*',
        'keep floatedmValueMap_eidRobustTight_*_*',
        'keep floatedmValueMap_eidRobustHighEnergy_*_*',
        'keep floatedmValueMap_eidLoose_*_*',
        'keep floatedmValueMap_eidTight_*_*',
        'keep *_gedPhotonCore_*_*',
        'keep *_gedPhotons_*_*',
        'keep recoPhotons_photons_*_*',
        'keep recoPhotonCores_photonCore_*_*',
        'keep recoConversions_conversions_*_*',
        'drop *_conversions_uncleanedConversions_*',
        'keep recoConversions_allConversions_*_*',
        'keep recoTracks_ckfOutInTracksFromConversions_*_*',
        'keep recoTracks_ckfInOutTracksFromConversions_*_*',
        'keep recoTrackExtras_ckfOutInTracksFromConversions_*_*',
        'keep recoTrackExtras_ckfInOutTracksFromConversions_*_*',
        'keep TrackingRecHitsOwned_ckfOutInTracksFromConversions_*_*',
        'keep TrackingRecHitsOwned_ckfInOutTracksFromConversions_*_*',
        'keep recoConversions_uncleanedOnlyAllConversions_*_*',
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*',
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*',
        'keep recoTrackExtras_uncleanedOnlyCkfOutInTracksFromConversions_*_*',
        'keep recoTrackExtras_uncleanedOnlyCkfInOutTracksFromConversions_*_*',
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfOutInTracksFromConversions_*_*',
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfInOutTracksFromConversions_*_*',
        'keep *_PhotonIDProd_*_*',
        'keep *_hfRecoEcalCandidate_*_*',
        'keep *_hfEMClusters_*_*',
        'keep *_gedGsfElectronCores_*_*',
        'keep *_gedGsfElectrons_*_*',
        'keep *_pixelTracks_*_*',
        'keep *_pixelVertices_*_*',
        'drop CaloTowersSorted_towerMakerPF_*_*',
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*',
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*',
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*',
        'keep recoPFRecHits_particleFlowClusterHFEM_Cleaned_*',
        'keep recoPFRecHits_particleFlowClusterHFHAD_Cleaned_*',
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*',
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*',
        'keep recoPFRecHits_particleFlowRecHitHCAL_Cleaned_*',
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*',
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*',
        'keep recoPFClusters_particleFlowClusterECAL_*_*',
        'keep recoPFClusters_particleFlowClusterHCAL_*_*',
        'keep recoPFClusters_particleFlowClusterHO_*_*',
        'keep recoPFClusters_particleFlowClusterPS_*_*',
        'keep recoPFBlocks_particleFlowBlock_*_*',
        'keep recoPFCandidates_particleFlowEGamma_*_*',
        'keep recoCaloClusters_particleFlowEGamma_*_*',
        'keep recoSuperClusters_particleFlowEGamma_*_*',
        'keep recoPFCandidates_particleFlow_*_*',
        'keep recoPFCandidates_particleFlowTmp_electrons_*',
        'keep recoPFCandidates_particleFlowTmp_*_*',
        'drop recoPFCandidates_particleFlowTmp__*',
        'keep recoPFDisplacedVertexs_particleFlowDisplacedVertex_*_*',
        'keep *_pfElectronTranslator_*_*',
        'keep *_pfPhotonTranslator_*_*',
        'keep *_particleFlow_electrons_*',
        'keep *_particleFlow_photons_*',
        'keep *_trackerDrivenElectronSeeds_preid_*',
        'keep *_particleFlowPtrs_*_*',
        'keep *_particleFlowTmpPtrs_*_*',
        'keep *_offlineBeamSpot_*_*',
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*',
        'keep *_l1GtRecord_*_*',
        'keep *_l1GtTriggerMenuLite_*_*',
        'keep *_conditionsInEdm_*_*',
        'keep *_l1extraParticles_*_*',
        'keep *_l1L1GtObjectMap_*_*',
        'keep L1MuGMTReadoutCollection_gtDigis_*_*',
        'keep L1GctEmCand*_gctDigis_*_*',
        'keep L1GctJetCand*_gctDigis_*_*',
        'keep L1GctEtHad*_gctDigis_*_*',
        'keep L1GctEtMiss*_gctDigis_*_*',
        'keep L1GctEtTotal*_gctDigis_*_*',
        'keep L1GctHtMiss*_gctDigis_*_*',
        'keep L1GctJetCounts*_gctDigis_*_*',
        'keep L1GctHFRingEtSums*_gctDigis_*_*',
        'keep L1GctHFBitCounts*_gctDigis_*_*',
        'keep LumiDetails_lumiProducer_*_*',
        'keep LumiSummary_lumiProducer_*_*',
        'drop *_hlt*_*_*',
        'keep *_hltL1GtObjectMap_*_*',
        'keep edmTriggerResults_*_*_*',
        'keep triggerTriggerEvent_*_*_*',
        'keep LHERunInfoProduct_*_*_*',
        'keep LHEEventProduct_*_*_*',
        'keep GenRunInfoProduct_generator_*_*',
        'keep GenEventInfoProduct_generator_*_*',
        'keep edmHepMCProduct_generator_*_*',
        'keep GenFilterInfo_*_*_*',
        'keep *_genParticles_*_*',
        'keep recoGenMETs_*_*_*',
        'keep *_kt4GenJets_*_*',
        'keep *_kt6GenJets_*_*',
        'keep *_ak5GenJets_*_*',
        'keep *_ak7GenJets_*_*',
        'keep *_iterativeCone5GenJets_*_*',
        'keep *_genParticle_*_*',
        'keep edmHepMCProduct_source_*_*',
        'keep SimTracks_g4SimHits_*_*',
        'keep SimVertexs_g4SimHits_*_*',
        'keep *_allTrackMCMatch_*_*',
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*',
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*',
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*',
        'keep PileupSummaryInfos_*_*_*',
        'keep L1AcceptBunchCrossings_*_*_*',
        'keep L1TriggerScalerss_*_*_*',
        'keep Level1TriggerScalerss_*_*_*',
        'keep LumiScalerss_*_*_*',
        'keep BeamSpotOnlines_*_*_*',
        'keep DcsStatuss_*_*_*',
        'keep *_logErrorHarvester_*_*',
        'keep *_pfIsolatedElectronsEI_*_*',
        'keep *_pfIsolatedMuonsEI_*_*',
        'keep recoPFJets_pfJetsEI_*_*',
        'keep *_pfJetTrackAssociatorEI_*_*',
        'keep *_impactParameterTagInfosEI_*_*',
        'keep *_secondaryVertexTagInfosEI_*_*',
        'keep *_combinedSecondaryVertexBJetTagsEI_*_*',
        'keep recoPFTaus_pfTausEI_*_*',
        'keep recoPFTauDiscriminator_pfTausDiscrimination*_*_*',
        'keep *_pfMetEI_*_*',
        'keep *_simCscTriggerPrimitiveDigis_*_*',
        'keep *_simDtTriggerPrimitiveDigis_*_*',
        'keep *_simRpcTriggerDigis_*_*',
        'keep *_simRctDigis_*_*',
        'keep *_simCsctfDigis_*_*',
        'keep *_simCsctfTrackDigis_*_*',
        'keep *_simDttfDigis_*_*',
        'keep *_simGctDigis_*_*',
        'keep *_simGmtDigis_*_*',
        'keep *_simGtDigis_*_*',
        'keep *_cscTriggerPrimitiveDigis_*_*',
        'keep *_dtTriggerPrimitiveDigis_*_*',
        'keep *_rpcTriggerDigis_*_*',
        'keep *_rctDigis_*_*',
        'keep *_csctfDigis_*_*',
        'keep *_csctfTrackDigis_*_*',
        'keep *_dttfDigis_*_*',
        'keep *_gctDigis_*_*',
        'keep *_gmtDigis_*_*',
        'keep *_gtDigis_*_*',
        'keep *_gtEvmDigis_*_*',
        'keep *_l1GtRecord_*_*',
        'keep *_l1GtTriggerMenuLite_*_*',
        'keep *_conditionsInEdm_*_*',
        'keep *_l1extraParticles_*_*',
        'keep *_l1L1GtObjectMap_*_*',
        'keep LumiDetails_lumiProducer_*_*',
        'keep LumiSummary_lumiProducer_*_*',
        'drop *_trackingtruthprod_*_*',
        'drop *_electrontruth_*_*',
        'keep *_mix_MergedTrackTruth_*',
        'keep CrossingFramePlaybackInfoExtended_*_*_*',
        'keep *_simSiPixelDigis_*_*',
        'keep *_simSiStripDigis_*_*',
        'drop *_mix_simSiPixelDigis*_*',
        'drop *_mix_simSiStripDigis*_*',
        'keep *_allTrackMCMatch_*_*',
        'keep *_trackingParticleRecoTrackAsssociation_*_*',
        'keep *_assoc2secStepTk_*_*',
        'keep *_assoc2thStepTk_*_*',
        'keep *_assoc2GsfTracks_*_*',
        'keep *_assocOutInConversionTracks_*_*',
        'keep *_assocInOutConversionTracks_*_*',
        'keep *_simMuonCSCDigis_*_*',
        'keep *_simMuonDTDigis_*_*',
        'keep *_simMuonRPCDigis_*_*',
        'keep *_simEcalDigis_*_*',
        'keep *_simEcalPreshowerDigis_*_*',
        'keep *_simEcalTriggerPrimitiveDigis_*_*',
        'keep *_simHcalDigis_*_*',
        'keep ZDCDataFramesSorted_simHcalUnsuppressedDigis_*_*',
        'drop ZDCDataFramesSorted_mix_simHcalUnsuppressedDigis*_*',
        'keep *_simHcalTriggerPrimitiveDigis_*_*',
        'keep *_siPixelRecHits_*_*',
        'keep *_simMuonCSCDigis_*_*',
        'keep *_simMuonRPCDigis_*_*',
        'keep *_simHcalUnsuppressedDigis_*_*',
        'keep *_TTClustersFromPixelDigis_*_*',
        'keep *_TTStubsFromPixelDigis_*_*',
        'keep *_TTTracksFromPixelDigis_*_*',
        'keep *_TTClusterAssociatorFromPixelDigis_*_*',
        'keep *_TTStubAssociatorFromPixelDigis_*_*',
        'keep *_TTTrackAssociatorFromPixelDigis_*_*',
        'drop PixelDigiSimLinkedmDetSetVector_mix_*_*',
        'drop PixelDigiedmDetSetVector_mix_*_*',
        'keep *_simSiPixelDigis_*_*' ) )  
)


#cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('myOutputFile.root')
#    ,outputCommands = cms.untracked.vstring(
#      "keep *_simSiPixelDigis_*_*",
#      "keep *_siPixelClusters_*_*",
#      "keep *_siPixelRecHits_*_*")      
#)

#RecoLocalTrackerRECO = cms.PSet(
  #  outputCommands = cms.untracked.vstring(
    #'keep DetIdedmEDCollection_siStripDigis_*_*',
#    'keep DetIdedmEDCollection_siPixelDigis_*_*',
  #  'keep *_siPixelClusters_*_*', 
#    'keep *_siStripClusters_*_*',
#    'keep *_clusterSummaryProducer_*_*')
#)

# Processes to run
process.pixeltrackerlocalreco = cms.Sequence(process.siPixelClusters*process.siPixelRecHits)
process.trackerlocalreco = cms.Sequence(process.pixeltrackerlocalreco)

process.p = cms.Path(process.trackerlocalreco)
process.e = cms.EndPath(process.out)
#process.p = cms.Path(process.analysis)

 # Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_phase2_BE5D 

#call to customisation function cust_phase2_BE5D imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_phase2_BE5D(process)
