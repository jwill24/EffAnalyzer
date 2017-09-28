import FWCore.ParameterSet.Config as cms

treeProducer = cms.EDAnalyzer('TreeProducer',

    # general parameters
    isData = cms.bool(True),
    sqrtS = cms.double(13.e3),
    outputFilename = cms.string('output.root'),

    # input collections
    vertexLabel = cms.InputTag('offlineSlimmedPrimaryVertices'),
    photonsLabel = cms.InputTag('photons'),

    # "tight" single/double photon selection
    minPtSinglePhoton = cms.double(50.),
    minR9SinglePhoton = cms.double(0.94),
    maxEtaSinglePhoton = cms.double(2.5),
    minMassDiPhoton = cms.double(500.),

    # totem RP information extraction
    totemRPTracksLabel = cms.InputTag('totemRPLocalTrackFitter'),
    totemRPUVPatternsLabel = cms.InputTag('totemRPUVPatternFinder'),
    totemRPClustersLabel = cms.InputTag('totemRPClusterProducer'),
    totemRPRecHitsLabel = cms.InputTag('totemRPRecHitProducer'),
    useXiInterpolation = cms.bool(True),
    xiInterpolationFile = cms.FileInPath('EffAnalyzer/TreeProducer/data/ctpps_optics_9mar2017.root'),
    fillNumLUTFile = cms.FileInPath('EffAnalyzer/TreeProducer/data/fill_run_lut_v2.dat'),
    alignmentLUTFile = cms.FileInPath('EffAnalyzer/TreeProducer/data/alignment_collection_v2.out'),

    # generator-level input collections
    maxGenLevelDeltaR = cms.double(5.0),
)
