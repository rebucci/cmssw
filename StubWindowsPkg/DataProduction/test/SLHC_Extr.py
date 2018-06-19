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
# Script tested with release CMSSW_10_0_0_pre1 (works either for Tilted of Flat geometries)
#
#########



############################################################
# define process
############################################################

import FWCore.ParameterSet.Config as cms

process = cms.Process("MIBextractor")

# Select geometry
flat=False



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
    input = cms.untracked.int32(1)
)

# input files
Source_Files = cms.untracked.vstring(
		##### EXAMPLES #####
		#'file:PGun_example.root',       
		#'file:PU_sample.root', 
		#'file:TT_example.root',       
		#'file:QCD_example.root', 

		##### RELVALS #####
		'/store/relval/CMSSW_10_0_0_pre1/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_94X_upgrade2023_realistic_v2_2023D17PU200-v1/10000/52A9842F-C9CF-E711-84DE-0242AC130002.root',  #TTBar PU200
		#'/store/relval/CMSSW_10_0_0_pre1/RelValSingleMuPt10/GEN-SIM-DIGI-RAW/94X_upgrade2023_realistic_v2_2023D17noPU-v2/10000/4A774661-DECE-E711-9826-0CC47A4C8F18.root',  #SingleMu PU0
)

process.source = cms.Source("PoolSource", fileNames = Source_Files, duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)

# name of output file
OUTPUT_NAME="TTBar200"  #output file will be "extracted_'OUTPUT_NAME'.root"



############################################################
# Extractor
############################################################

process.load("StubWindowsPkg.RecoExtractor.MIB_extractor_cff")


############################################################
# Customization

process.MIBextraction.doMC       = True
process.MIBextraction.doSTUB     = True
process.MIBextraction.doPixel    = True
process.MIBextraction.doMatch    = True

process.MIBextraction.doL1TRK	 = True
process.MIBextraction.doBANK 	 = False
process.MIBextraction.getCoords  = False
process.MIBextraction.fullInfo   = True


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
process.MIBextraction.extractedRootFile=cms.string('extracted_'+OUTPUT_NAME+'.root')

# Note: preferred method uses TFileService. It's unclear if the current method is compatible with lobster.
#process.TFileService = cms.Service("TFileService", fileName = cms.string('output.root'), closeFileFast = cms.untracked.bool(True))
# This will require changes to MIBextractor (remove the output options) as well as all the extractors (add the TFileService output to each) and possibly more. Pending review.



############################################################
# Import/Load Processes
############################################################

############################################################
# General
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load("SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff")


############################################################
# Clusters and Stubs
from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *

from SimTracker.TrackTriggerAssociation.TTClusterAssociation_cfi import *
TTClusterAssociatorFromPixelDigis.digiSimLinks = cms.InputTag("simSiPixelDigis","Tracker")

process.TTClusterStubTruth = cms.Path(process.TrackTriggerAssociatorClustersStubs) #clusters and stubs


############################################################
# L1 Tracks
# from SimTracker.TrackTriggerAssociation.TTTrackAssociation_cfi import *
# TTTrackAssociatorFromPixelDigis.TTTracks = cms.VInputTag(cms.InputTag("TTTracksFromPhase2TrackerDigis", "Level1TTTracks") )

# process.TTTracksTruth         = cms.Path(process.TrackTriggerAssociatorTracks) #ONLY tracks with truth
# process.TTTracksCompleteTruth = cms.Path(process.TrackTriggerAssociatorComplete) #stubs, clusters, and tracks with trutha

from L1Trigger.TrackFindingTracklet.Tracklet_cfi import *
process.load("L1Trigger.TrackFindingTracklet.L1TrackletEmulationTracks_cff")
process.TTTracksEmulation          = cms.Path(process.L1TrackletEmulationTracks)
process.TTTracksEmulationWithTruth = cms.Path(process.L1TrackletEmulationTracksWithAssociators)


############################################################
# Processes
############################################################

process.dump = cms.EDAnalyzer("EventContentAnalyzer")
process.p    = cms.Path(process.MIBextraction)
process.d    = cms.Path(process.dump)

# process.TTClusterStubTruth     -  clusters and stubs
# process.TTTracksTruth          -  l1 tracks
# process.TTTracksCompleteTruth  -  l1 tracks, clusters, stubs
# process.TTTracksEmulationWithTruth - l1 tracks from emulator

process.schedule = cms.Schedule(process.TTClusterStubTruth,process.TTTracksEmulationWithTruth,process.p) 
#process.schedule = cms.Schedule(process.TTTracksCompleteTruth,process.p) 
#process.schedule = cms.Schedule(process.TTClusterStubTruth,process.p) 
