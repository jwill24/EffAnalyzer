import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:/afs/cern.ch/work/j/juwillia/CMSSW_8_0_25/src/EffAnalyzer/TreeProducer/test/040C2700-A388-E611-AF56-00259073E390.root'
#        'file:/afs/cern.ch/work/j/juwillia/040C2700-A388-E611-AF56-00259073E390.root'
#        'file:/afs/cern.ch/work/j/juwillia/00101B51-B397-E611-889D-0025907B4E6C.root'
#        'file:/afs/cern.ch/work/j/juwillia/029C9BE0-C997-E611-9426-F04DA275BFEC.root',
        'file:/afs/cern.ch/work/j/juwillia/081016B7-6A9E-E711-B5E2-0CC47A7C361E.root',
#        'file:/afs/cern.ch/work/j/juwillia/00948452-C997-E611-8B51-001E67E6F8DC.root'
#'/store/group/phys_pps/diphoton/DoubleEG/lforthom-microAOD-ctpps_Run2016C-23Sep2016_v3/170303_022624/0000/myMicroAODOutputFile_64.root',
    )
)

process.load('EffAnalyzer.TreeProducer.TreeProducer_cfi')

process.p = cms.Path(
    process.treeProducer
)
