import FWCore.ParameterSet.Config as cms

phase2PixelDigitizer = cms.PSet(
    accumulatorType = cms.string("SiPhase2Digitizer"),
    hitsProducer = cms.string('g4SimHits'),
    makeDigiSimLinks = cms.untracked.bool(True),
#D.B.:the noise should be a function of strip capacitance, roughly: ReadoutNoiseInElec=500+(64*Cdet[pF]) ~= 500+(64*1.5[cm])
    ReadoutNoiseInElec = cms.double(1000.0),#D.B.:Fill readout noise, including all readout chain, relevant for smearing
    DeltaProductionCut = cms.double(0.03),
    ROUList = cms.vstring(
        'TrackerHitsPixelBarrelLowTof', 
        'TrackerHitsPixelBarrelHighTof', 
        'TrackerHitsPixelEndcapLowTof', 
        'TrackerHitsPixelEndcapHighTof'),
    ThresholdInElectrons_Barrel = cms.double(4759.8), #D.B.: this should correspond to a threshold of 530mV    
    ThresholdInElectrons_Endcap = cms.double(4759.8),
    AddThresholdSmearing = cms.bool(True),
    ThresholdSmearing_Barrel = cms.double(204.0),#D.B.: changed (~5mV peakToPeak --> 1.76mV rms) (was 210.0)
    ThresholdSmearing_Endcap = cms.double(204.0),#D.B.: changed (~5mV peakToPeak --> 1.76mV rms) (was 245.0)
    NoiseInElectrons = cms.double(300),	         # 30% of the readout noise (should be changed in future)
    DigitalReadout           = cms.bool(True), # Flag to decide analog or digital readout 
    ElectronPerAdc = cms.double(135.0),	#D.B.:used for misscalibration
    TofUpperCut = cms.double(12.5),
    AdcFullScale = cms.int32(255),
    TofLowerCut = cms.double(-12.5),
    AddNoisyPixels = cms.bool(True),
    Alpha2Order = cms.bool(True),			#D.B.: second order effect, does not switch off magnetic field as described
    AddNoise = cms.bool(True),
    AddXTalk = cms.bool(True),			#D.B.
    InterstripCoupling = cms.double(0.08),	#D.B.
    SigmaZero = cms.double(0.00037),  		#D.B.: 3.7um spread for 300um-thick sensor, renormalized in digitizerAlgo
    SigmaCoeff = cms.double(1.80),  		#D.B.: to be confirmed with simulations in CMSSW_6.X
    ClusterWidth = cms.double(3),		#D.B.: this is used as number of sigmas for charge collection (3=+-3sigmas)
    
    AddInefficiency = cms.bool(False),
    Inefficiency_DB = cms.bool(True),				

    LorentzAngle_DB = cms.bool(True),			
    TanLorentzAnglePerTesla_Endcap = cms.double(0.106),	#D.B.:this I have not checked yet
    TanLorentzAnglePerTesla_Barrel = cms.double(0.106),	#D.B.:this I have not checked yet

    KillModules = cms.bool(True),
    DeadModules_DB = cms.bool(True),

    EfficiencyFactors_Barrel = cms.vdouble(0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999 ),
    EfficiencyFactors_Endcap = cms.vdouble(0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999 ),#Efficiencies kept as Side2Disk1,Side1Disk1 and so on
    DeadModules = cms.VPSet()
)
