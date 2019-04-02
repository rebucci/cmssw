#########
#
# Example script to run the python extractor on MC events
# for the skimmed geometry
# 
# Usage: cmsRun SLHC_Extr.py
#
# Author: S.Viret (s.viret@ipnl.in2p3.fr)
# Date  : 24/05/2017
#
# Script tested with release CMSSW_10_4_0
#
# Script, extractors, and tools updated by R.Bucci (rbucci@nd.edu). 2018.
# 
#########



############################################################
# define process
############################################################

import FWCore.ParameterSet.Config as cms

process = cms.Process("MIBextractor")

# Select geometry
flat=False

# Select stub windows
STUBWINDOWS = "Tab"
    # Tab  		 stub  windows in CMSSW_9_4_0 and 10_0_0
    # Tight    tight windows in CMSSW_10_2_X
    # Loose    loose windows in CMSSW_10_2_X
    # Zero     everything is zero!



############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Global tag for PromptReco
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')



############################################################
# input
############################################################

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# input files
Source_Files = cms.untracked.vstring(
        #'/store/relval/CMSSW_10_4_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/103X_upgrade2023_realistic_v2_2023D21noPU-v1/20000/F4344045-AEDE-4240-B7B1-27D2CF96C34E.root',
)

process.source = cms.Source("PoolSource", 
                            fileNames = Source_Files, 
                            duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
                            )

# name of output file
OUTPUT_NAME="extracted.root"


############################################################
# Extractor
############################################################

process.load("StubWindowsPkg.RecoExtractor.MIB_extractor_cff")


############################################################
# Customization

process.MIBextraction.doMC       = True     # adds tracking particle information
process.MIBextraction.doSTUB     = True     # adds official stubs and clusters
process.MIBextraction.doPixel    = True     # activate digis or not (not for official)
process.MIBextraction.doMatch    = True     # adds digi matching to MC (used only if doPixel is true)

process.MIBextraction.doL1TRK    = True     # adds patterns and L1 track info; default false
process.MIBextraction.doBANK 	   = False    # default false
process.MIBextraction.getCoords  = False    # default false
process.MIBextraction.fullInfo   = True     # default false

process.MIBextraction.TP_hitTracker = True        # Require tracking particle to hit in the tracker
process.MIBextraction.TP_minPt      = 2.0         # minimum pT of tracking particles  (default=0.0)
process.MIBextraction.TP_maxEta     = 2.4         # minimum eta of tracking particles (default=5.5)
process.MIBextraction.TP_maxR       = 10000000.0  # maximum R of tracking particles   (default=10000000.0)



############################################################
# Geometry
if flat:
	print 'You choose the flat geometry'
	process.MIBextraction.flatBarrel = True
	process.load('L1Trigger.TrackTrigger.TkOnlyFlatGeom_cff') # Special config file for TkOnly geometry
else:
	print 'You choose the tilted geometry'
	process.MIBextraction.flatBarrel = False
	process.load('L1Trigger.TrackTrigger.TkOnlyTiltedGeom_cff') # Special config file for TkOnly geometry


############################################################
# OUTPUT
process.MIBextraction.extractedRootFile=cms.string(OUTPUT_NAME)



############################################################
# Import/Load Processes
############################################################

process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *

# Stub Windows
if   STUBWINDOWS == "Tab": 
    print "using stub windows for CMSSW_9_4_0 and 10_0_0"
    TTStubAlgorithm_official_Phase2TrackerDigi_.BarrelCut = cms.vdouble( 0, 2.0, 2.0, 3.5, 4.5, 5.5, 6.5)
    TTStubAlgorithm_official_Phase2TrackerDigi_.TiltedBarrelCutSet = cms.VPSet(
        cms.PSet( TiltedCut = cms.vdouble( 0 ) ),
        cms.PSet( TiltedCut = cms.vdouble( 0, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2., 2., 1.5, 1.5, 1., 1.) ),
        cms.PSet( TiltedCut = cms.vdouble( 0, 3., 3., 3., 3., 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2, 2) ),
        cms.PSet( TiltedCut = cms.vdouble( 0, 4.5, 4.5, 4, 4, 4, 4, 3.5, 3.5, 3.5, 3, 3, 3) ),
    )
    TTStubAlgorithm_official_Phase2TrackerDigi_.EndcapCutSet = cms.VPSet(
        cms.PSet( EndcapCut = cms.vdouble( 0 ) ),
        cms.PSet( EndcapCut = cms.vdouble( 0, 1, 1.5, 1.5, 2, 2, 2.5, 3, 3, 3.5, 4, 2.5, 3, 3.5, 4.5, 5.5) ),
        cms.PSet( EndcapCut = cms.vdouble( 0, 1, 1.5, 1.5, 2, 2, 2, 2.5, 3, 3, 3, 2, 3, 4, 5, 5.5) ),
        cms.PSet( EndcapCut = cms.vdouble( 0, 1.5, 1.5, 2, 2, 2.5, 2.5, 2.5, 3.5, 2.5, 5, 5.5, 6) ),
        cms.PSet( EndcapCut = cms.vdouble( 0, 1.0, 1.5, 1.5, 2, 2, 2, 2, 3, 3, 6, 6, 6.5) ),
        cms.PSet( EndcapCut = cms.vdouble( 0, 1.0, 1.5, 1.5, 1.5, 2, 2, 2, 3, 3, 6, 6, 6.5) ),
    )
elif STUBWINDOWS == "Tight": 
    print "using the tight stub windows for CMSSW_10_2_X"
    TTStubAlgorithm_official_Phase2TrackerDigi_.BarrelCut = cms.vdouble( 0, 2, 2.5, 3.5, 4.5, 5.5, 7)
    TTStubAlgorithm_official_Phase2TrackerDigi_.TiltedBarrelCutSet = cms.VPSet(
                cms.PSet( TiltedCut = cms.vdouble( 0 ) ),
                cms.PSet( TiltedCut = cms.vdouble( 0, 3, 3, 2.5, 3, 3, 2.5, 2.5, 2, 1.5, 1.5, 1, 1) ),
                cms.PSet( TiltedCut = cms.vdouble( 0, 3.5, 3, 3, 3, 3, 2.5, 2.5, 3, 3, 2.5, 2.5, 2.5) ),
                cms.PSet( TiltedCut = cms.vdouble( 0, 4, 4, 4, 3.5, 3.5, 3.5, 3.5, 3, 3, 3, 3, 3) ),
    )
    TTStubAlgorithm_official_Phase2TrackerDigi_.EndcapCutSet = cms.VPSet(
                cms.PSet( EndcapCut = cms.vdouble( 0 ) ),
                cms.PSet( EndcapCut = cms.vdouble( 0, 1, 2.5, 2.5, 3, 2.5, 3, 3.5, 4, 4, 4.5, 3.5, 4, 4.5, 5, 5.5) ),
                cms.PSet( EndcapCut = cms.vdouble( 0, 0.5, 2.5, 2.5, 3, 2.5, 3, 3, 3.5, 3.5, 4, 3.5, 3.5, 4, 4.5, 5) ),
                cms.PSet( EndcapCut = cms.vdouble( 0, 1, 3, 3, 2.5, 3.5, 3.5, 3.5, 4, 3.5, 3.5, 4, 4.5) ),
                cms.PSet( EndcapCut = cms.vdouble( 0, 1, 2.5, 3, 2.5, 3.5, 3, 3, 3.5, 3.5, 3.5, 4, 4) ),
                cms.PSet( EndcapCut = cms.vdouble( 0, 0.5, 1.5, 3, 2.5, 3.5, 3, 3, 3.5, 4, 3.5, 4, 3.5) ),
    )
elif STUBWINDOWS == "Loose": 
    print "using the loose stub windows for CMSSW_10_2_X"
    TTStubAlgorithm_official_Phase2TrackerDigi_.BarrelCut = cms.vdouble( 0, 2.0, 3, 4.5, 6, 6.5, 7.0)
    TTStubAlgorithm_official_Phase2TrackerDigi_.TiltedBarrelCutSet = cms.VPSet(
                cms.PSet( TiltedCut = cms.vdouble( 0 ) ),
                cms.PSet( TiltedCut = cms.vdouble( 0, 3, 3., 2.5, 3., 3., 2.5, 2.5, 2., 1.5, 1.5, 1, 1) ),
                cms.PSet( TiltedCut = cms.vdouble( 0, 4., 4, 4, 4, 4., 4., 4.5, 5, 4., 3.5, 3.5, 3) ),
                cms.PSet( TiltedCut = cms.vdouble( 0, 5, 5, 5, 5, 5, 5, 5.5, 5, 5, 5.5, 5.5, 5.5) ),
    )
    TTStubAlgorithm_official_Phase2TrackerDigi_.EndcapCutSet = cms.VPSet(
                cms.PSet( EndcapCut = cms.vdouble( 0 ) ),
                cms.PSet( EndcapCut = cms.vdouble( 0, 1., 2.5, 2.5, 3.5, 5.5, 5.5, 6, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 7, 7) ),
                cms.PSet( EndcapCut = cms.vdouble( 0, 0.5, 2.5, 2.5, 3, 5, 6, 6, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 7, 7) ),
                cms.PSet( EndcapCut = cms.vdouble( 0, 1, 3., 4.5, 6., 6.5, 6.5, 6.5, 7, 7, 7, 7, 7) ),
                cms.PSet( EndcapCut = cms.vdouble( 0, 1., 2.5, 3.5, 6., 6.5, 6.5, 6.5, 6.5, 7, 7, 7, 7) ),
                cms.PSet( EndcapCut = cms.vdouble( 0, 0.5, 1.5, 3., 4.5, 6.5, 6.5, 7, 7, 7, 7, 7, 7) ),
    )
elif STUBWINDOWS == "Zero": 
    print "If there's nothing wrong with me...maybe there's something wrong with the universe."
    TTStubAlgorithm_official_Phase2TrackerDigi_.BarrelCut = cms.vdouble( 0, 0, 0, 0, 0, 0, 0)
    TTStubAlgorithm_official_Phase2TrackerDigi_.TiltedBarrelCutSet = cms.VPSet(
                cms.PSet( TiltedCut = cms.vdouble( 0 ) ),
                cms.PSet( TiltedCut = cms.vdouble( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) ),
                cms.PSet( TiltedCut = cms.vdouble( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) ),
                cms.PSet( TiltedCut = cms.vdouble( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) ),
    )
    TTStubAlgorithm_official_Phase2TrackerDigi_.EndcapCutSet = cms.VPSet(
                cms.PSet( EndcapCut = cms.vdouble( 0 ) ),
                cms.PSet( EndcapCut = cms.vdouble( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) ),
                cms.PSet( EndcapCut = cms.vdouble( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) ),
                cms.PSet( EndcapCut = cms.vdouble( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) ),
                cms.PSet( EndcapCut = cms.vdouble( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) ),
                cms.PSet( EndcapCut = cms.vdouble( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) ),
    )


process.load("SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff")
from SimTracker.TrackTriggerAssociation.TTClusterAssociation_cfi import *
TTClusterAssociatorFromPixelDigis.digiSimLinks = cms.InputTag("simSiPixelDigis","Tracker")

process.TTClusterStub      = cms.Path(process.TrackTriggerClustersStubs)
process.TTClusterStubTruth = cms.Path(process.TrackTriggerAssociatorClustersStubs)

process.TTTracksTruth = cms.Path(process.TrackTriggerAssociatorTracks)
process.TTCSTTruth    = cms.Path(process.TrackTriggerAssociatorComplete) # both the TTAssociatorTracks and the TTAssociatorClustersStubs


from L1Trigger.TrackFindingTracklet.Tracklet_cfi import *
process.load("L1Trigger.TrackFindingTracklet.L1TrackletEmulationTracks_cff")
process.TTTracksEmulation          = cms.Path(process.L1TrackletEmulationTracks)
process.TTTracksEmulationWithTruth = cms.Path(process.L1TrackletEmulationTracksWithAssociators)

from L1Trigger.TrackTrigger.TTCluster_cfi import * # attempted fix


############################################################
# Processes
############################################################

#process.MessageLogger = cms.Service("MessageLogger",
#		destinations = cms.untracked.vstring('logDebug'),
#			logDebug   = cms.untracked.PSet(threshold = cms.untracked.string('DEBUG') )
#)

process.dump = cms.EDAnalyzer("EventContentAnalyzer")
#process.Tracer = cms.Service("Tracer")

process.p    = cms.Path(process.MIBextraction)
#process.d    = cms.Path(process.dump)

# Don't re-run stub making (will use to the stub windows original to the release area)
# process.schedule = cms.Schedule(process.TTClusterStubTruth,process.TTTracksEmulationWithTruth,process.p) 

# Re-run stub making
process.schedule = cms.Schedule(process.TTClusterStub,process.TTClusterStubTruth,process.TTTracksEmulationWithTruth,process.p) 
