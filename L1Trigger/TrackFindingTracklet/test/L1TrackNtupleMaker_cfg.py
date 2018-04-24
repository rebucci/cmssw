############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("L1TrackNtuple")

## Configured for D17 (T5) Geometry ONLY


############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17_cff')

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')



############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

Source_Files = cms.untracked.vstring(
)

process.source = cms.Source("PoolSource", fileNames = Source_Files)

process.TFileService = cms.Service("TFileService", fileName = cms.string('output.root'), closeFileFast = cms.untracked.bool(True))


############################################################
# L1 tracking
############################################################

# remake stubs ?
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *
process.load("SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff") #emulator

from SimTracker.TrackTriggerAssociation.TTClusterAssociation_cfi import *
TTClusterAssociatorFromPixelDigis.digiSimLinks = cms.InputTag("simSiPixelDigis","Tracker")

process.TTClusterStub      = cms.Path(process.TrackTriggerClustersStubs)
process.TTClusterStubTruth = cms.Path(process.TrackTriggerAssociatorClustersStubs)

from L1Trigger.TrackFindingTracklet.Tracklet_cfi import *

#### floating-point version
#
#process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff")
#if GEOMETRY == "D10": 
#    TTTracksFromTracklet.trackerGeometry = cms.untracked.string("flat")
#TTTracksFromTracklet.asciiFileName = cms.untracked.string("evlist.txt")
#
## run only the tracking (no MC truth associators)
#process.TTTracks = cms.Path(process.L1TrackletTracks)
#
## run the tracking AND MC truth associators)
#process.TTTracksWithTruth = cms.Path(process.L1TrackletTracksWithAssociators)


### emulation instead
process.load("L1Trigger.TrackFindingTracklet.L1TrackletEmulationTracks_cff")
process.TTTracksEmulation          = cms.Path(process.L1TrackletEmulationTracks)
process.TTTracksEmulationWithTruth = cms.Path(process.L1TrackletEmulationTracksWithAssociators)
#TTTracksFromTrackletEmulation.asciiFileName = cms.untracked.string("evlist.txt")


############################################################
# Define the track ntuple process, MyProcess is the (unsigned) PDGID corresponding to the process which is run
# e.g. single electron/positron = 11
#      single pion+/pion- = 211
#      single muon+/muon- = 13
#      pions in jets = 6
#      taus = 15
#      all TPs = 1
############################################################

process.L1TrackNtuple = cms.EDAnalyzer('L1TrackNtupleMaker',
        MyProcess        = cms.int32(1),
        DebugMode        = cms.bool(False),  # printout lots of debug statements
        SaveStubs        = cms.bool(True),   # save some info for *all* stubs
        L1Tk_nPar        = cms.int32(4),     # use 4 or 5-parameter L1 track fit ??
        L1Tk_minNStub    = cms.int32(4),     # L1 tracks with >= 4 stubs
        TP_minNStub      = cms.int32(4),     # require TP to have >= X number of stubs associated with it
        TP_minNStubLayer = cms.int32(4),     # require TP to have stubs in >= X layers/disks
        TP_minPt         = cms.double(2.0),  # only save TPs with pt > X GeV
        TP_maxEta        = cms.double(2.4),  # only save TPs with |eta| < X
        TP_maxZ0         = cms.double(30.0), # only save TPs with |z0| < X cm
        #L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),           # floating point
        L1TrackInputTag  = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"),  # emulator
        MCTruthTrackInputTag     = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"), ## MCTruth input
        # other input collections
        L1StubInputTag           = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted"),
        MCTruthClusterInputTag   = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
        MCTruthStubInputTag      = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
        TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
        TrackingVertexInputTag   = cms.InputTag("mix", "MergedTrackTruth"),
        GenJetInputTag           = cms.InputTag("ak4GenJets")
        )


# ### Output definition to save stubs
# process.FEVTDEBUGHLToutput = cms.OutputModule(
#     "PoolOutputModule",
#     dataset = cms.untracked.PSet(
#         dataTier   = cms.untracked.string('GEN-SIM-DIGI-RAW'),
#         filterName = cms.untracked.string('')
#     ),
#     fileName       = cms.untracked.string('stubs_932_TTBar_PU0.root'),
#     outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
#     splitLevel     = cms.untracked.int32(0)
# )
# process.FEVTDEBUGHLToutput.outputCommands.append(
#     'keep Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTClusterAssociationMap_TTClusterAssociatorFromPixelDigis_*_*')


### Processes
process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.ana  = cms.Path(process.L1TrackNtuple)
process.d    = cms.Path(process.dump)
#process.o    = cms.EndPath(process.FEVTDEBUGHLToutput)


##### PROCESS SCHEDULER #####

# To re-run the stub making
#process.schedule = cms.Schedule(process.TTClusterStub,process.TTClusterStubTruth,process.TTTracksEmulationWithTruth,process.ana,process.o) 

# If cluster/stub associators not available
#process.schedule = cms.Schedule(process.TTClusterStubTruth,process.TTTracksEmulationWithTruth,process.ana,process.o)
process.schedule = cms.Schedule(process.TTClusterStubTruth,process.TTTracksEmulationWithTruth,process.ana)

# To run only tracking + track associator
#process.schedule = cms.Schedule(process.TTTracksEmulationWithTruth,process.ana) 
#does not work with files straight from DAS