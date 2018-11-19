# MIB_extractor_cfi.py

############################################################
# Define the MIB Extraction Process
############################################################
import FWCore.ParameterSet.Config as cms


MIBextraction = cms.EDAnalyzer("RecoExtractor",    
  
    ## Setup 
    extractedRootFile = cms.string('extracted.root'), #name of the output ROOTfile                         
    inputRootFile     = cms.string('default.root'),   #name of the input ROOTfile, if you start from already extracted file
    flatBarrel        = cms.untracked.bool(True), #flat geometry (true) is default; but the extr.py's can override this
    fillTree          = cms.untracked.bool(True), #start from a RECO file (True) or an extracted ROOTuple (False)
    n_events          = cms.untracked.int32(10),  # How many events you want to analyze (only if fillTree=False)
    skip_events       = cms.untracked.int32(0),   # How many events you want to skip (only if fillTree=False)

    ## MC                      
    doMC             = cms.untracked.bool(False),          # Extract the MC information (MC tree)
    GenParticles     = cms.InputTag("genParticles", ""),
    TrkParticles     = cms.InputTag("mix" , "MergedTrackTruth"),
    SimHits          = cms.InputTag("g4SimHits"),
    TP_hitTracker    = cms.untracked.bool(False),  # only save TPs that have hits in the tracker
    TP_minPt         = cms.untracked.double(0.0),        # default at zero
    TP_maxEta        = cms.untracked.double(5.5),        # default taken from original StubExtractor.cc call of clearTP
    TP_maxR          = cms.untracked.double(10000000.0), # default taken from original StubExtractor.cc call of clearTP

    ## Stubs
    doSTUB                = cms.untracked.bool(False),     # Extract the official STUB information (TkStub tree)
    TTClusters            = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterInclusive"),
    TTClustersAssociators = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterInclusive"),
    TTStubs               = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
    TTStubsAssociators    = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),

    CLUS_container   = cms.string( "TTClustersFromPhase2TrackerDigis" ),
    STUB_container   = cms.string( "TTStubsFromPhase2TrackerDigis" ),
    CLUS_name        = cms.string( "ClusterInclusive" ),
    STUB_name        = cms.string( "StubAccepted" ),

    ## L1 Tracks
    doL1TRK          = cms.untracked.bool(False),          # Extract the official L1track information
    L1track_tag	     = cms.InputTag( "TTTracksFromTrackletEmulation", "Level1TTTracks"), 
                                 
    ## Pixel information                              
    doPixel          = cms.untracked.bool(False),          # Extract the Tracker information (Pixel tree)
    digi_tag         = cms.InputTag( "mix","Tracker" ),    # The collection where to find the pixel info
    digisimlink_tag  = cms.InputTag( "simSiPixelDigis","Tracker" ),    # The collection where to find the pixel info
    PUinfo_tag       = cms.InputTag( "addPileupInfo","" ), # The collection where to find the PU stuff
    doMatch          = cms.untracked.bool(False),          # Add the simtrack index to each digi                  

    ## Additional Options                            
    doBANK           = cms.untracked.bool(False),          # Stub container creation for bank (need official stubs)
    doL1TT           = cms.untracked.bool(False),          # Extract the cluster/stub information
    getCoords        = cms.untracked.bool(False),          # Extract the tracker coordinates
    fullInfo         = cms.untracked.bool(False),          # Extract the ALL tracker coordinates or just modids
    
    # Analysis Settings
    analysisSettings = cms.untracked.vstring() # Format is "STRING VALUE" where STRING is the name of the cut, and VALUE the value of the cut. See demo scripts for usage

)
