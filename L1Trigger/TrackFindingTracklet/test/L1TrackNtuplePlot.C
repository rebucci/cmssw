// -----------------------------------------------------------------------------------------------------------
// Basic example script for making tracking performance plots using the ntuples produced by L1TrackNtupleMaker.cc
// By Louise Skinnari, June 2013
// -----------------------------------------------------------------------------------------------------------

#include "TROOT.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "THStack.h"
#include "TMath.h"
#include <TError.h>

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void SetPlotStyle();
void mySmallText(Double_t x,Double_t y,Color_t color,char *text);

double getIntervalContainingFractionOfEntries( TH1* histogram, double interval );
void makeResidualIntervalPlot( TString type, TString dir, TString variable, TH1F* h_68, TH1F* h_90, TH1F* h_99, double minY, double maxY, TString htitle );


// -----------------------------------------------------------------------------------------------------------
// Main script
// -----------------------------------------------------------------------------------------------------------


void L1TrackNtuplePlot(TString type, int TP_select_pdgid=0, int TP_select_eventid=0, float TP_minPt=2.0, float TP_maxPt=7000.0, float TP_maxEta=2.4) {

  // type:              this is the input file you want to process (minus ".root" extension)
  // TP_select_pdgid:   if non-zero, only select TPs with a given PDG ID
  // TP_select_eventid: if zero, only look at TPs from primary interaction, else, include TPs from pileup
  // TP_minPt:          only look at TPs with pt > X GeV
  // TP_maxPt:          only look at TPs with pt < X GeV
  // TP_maxEta:         only look at TPs with |eta| < X


  gROOT->SetBatch();
  gErrorIgnoreLevel = kWarning;

  SetPlotStyle();


  // ---------------------------------------------------------------------------------------------------------
  // Define Input Options

  int L1Tk_minNstub     = 4;
  float L1Tk_maxChi2    = 999999;
  float L1Tk_maxChi2dof = 999999.;

  bool doLooseMatch      = false;    //MC truth matching

  bool doDetailedPlots   = false;    //Full set of plots for all pT and eta ranges
  bool doResolutionPlots = false;    //Advanced resolution plots including confidence 
  bool doGausFit         = false;    //Do gaussian fit for resolution vs eta/pt plots
  


  // Jet Counting
  const int njets = 3;     // number of (closest) jets counted in jet tracking

  // Some Counters for Integrated Efficiencies
  int n_all_eta2p5    = 0;
  int n_all_eta1p75   = 0;
  int n_all_eta1p0    = 0;
  int n_match_eta2p5  = 0;
  int n_match_eta1p75 = 0;
  int n_match_eta1p0  = 0;

  // Counters for Total Track Rates
  int ntrk_pt2  = 0;
  int ntrk_pt3  = 0;
  int ntrk_pt10 = 0;
  int ntp_pt2   = 0;
  int ntp_pt3   = 0;
  int ntp_pt10  = 0;


  // ---------------------------------------------------------------------------------------------------------
  // Read Ntuples
  TString srcDIR = "/hadoop/store/user/rbucci/Datasets_Emulator_20180420/";

  TChain* tree = new TChain("L1TrackNtuple/eventTree");
  tree->Add(srcDIR+type+".root");

  if (tree->GetEntries() == 0) {
    cout << "File doesn't exist or is empty, returning..." << endl;
    return;
  }


  // ---------------------------------------------------------------------------------------------------------
  // Define Leafs & Branches

  // tracking particles
  vector<float>* tp_pt;
  vector<float>* tp_eta;
  vector<float>* tp_phi;
  vector<float>* tp_dxy;
  vector<float>* tp_z0;
  vector<float>* tp_d0;
  vector<int>*   tp_pdgid;
  vector<int>*   tp_nmatch;
  vector<int>*   tp_nstub;
  vector<int>*   tp_eventid;
  vector<int>*   tp_nloosematch;

  // *L1 track* properties, for tracking particles matched to a L1 track
  vector<float>* matchtrk_pt;
  vector<float>* matchtrk_eta;
  vector<float>* matchtrk_phi;
  vector<float>* matchtrk_d0;
  vector<float>* matchtrk_z0;
  vector<float>* matchtrk_chi2;
  vector<int>*   matchtrk_nstub;

  // *L1 track* properties if m_tp_nloosematch > 0
  vector<float>* loosematchtrk_pt;
  vector<float>* loosematchtrk_eta;
  vector<float>* loosematchtrk_phi;
  vector<float>* loosematchtrk_d0;
  vector<float>* loosematchtrk_z0;
  vector<float>* loosematchtrk_chi2; 
  vector<int>*   loosematchtrk_nstub;

  // all L1 tracks
  vector<float>* trk_pt;
  vector<float>* trk_eta;
  vector<float>* trk_phi;
  vector<float>* trk_d0;
  vector<float>* trk_z0;
  vector<float>* trk_chi2;
  vector<int>*   trk_nstub;
  vector<int>*   trk_injet;
  vector<int>*   trk_injet_highpt;

  vector<int>*   trk_combinatoric;
  vector<int>*   trk_genuine;
  vector<int>*   trk_loose;
  vector<int>*   trk_unknown;
  vector<int>*   trk_fake;

  // Jets
  vector<float>* jet_eta;
  vector<float>* jet_pt;
  vector<float>* jet_phi;
  vector<float>* jet_tp_sumpt;
  vector<float>* jet_trk_sumpt;
  vector<float>* jet_matchtrk_sumpt;  
  vector<float>* jet_loosematchtrk_sumpt;
  

  // TBranches
  TBranch* b_tp_pt;
  TBranch* b_tp_eta;
  TBranch* b_tp_phi;
  TBranch* b_tp_dxy;
  TBranch* b_tp_z0;
  TBranch* b_tp_d0;
  TBranch* b_tp_pdgid;
  TBranch* b_tp_nmatch;
  TBranch* b_tp_nstub;
  TBranch* b_tp_eventid;
  TBranch* b_tp_nloosematch;

  TBranch* b_matchtrk_pt;
  TBranch* b_matchtrk_eta;
  TBranch* b_matchtrk_phi;
  TBranch* b_matchtrk_d0;
  TBranch* b_matchtrk_z0;
  TBranch* b_matchtrk_chi2;
  TBranch* b_matchtrk_nstub;

  // *L1 track* properties if m_tp_nloosematch > 0
  TBranch* b_loosematchtrk_pt;
  TBranch* b_loosematchtrk_eta;
  TBranch* b_loosematchtrk_phi;
  TBranch* b_loosematchtrk_d0;
  TBranch* b_loosematchtrk_z0;
  TBranch* b_loosematchtrk_chi2; 
  TBranch* b_loosematchtrk_nstub;

  TBranch* b_trk_pt;
  TBranch* b_trk_eta;
  TBranch* b_trk_phi;
  TBranch* b_trk_d0;
  TBranch* b_trk_z0;
  TBranch* b_trk_chi2;
  TBranch* b_trk_nstub;

  TBranch* b_trk_combinatoric;
  TBranch* b_trk_genuine;
  TBranch* b_trk_loose;
  TBranch* b_trk_unknown;
  TBranch* b_trk_fake;

  TBranch* b_jet_eta;
  TBranch* b_jet_pt;
  TBranch* b_jet_phi;
  TBranch* b_jet_tp_sumpt;
  TBranch* b_jet_trk_sumpt;
  TBranch* b_jet_matchtrk_sumpt;
  TBranch* b_jet_loosematchtrk_sumpt;
  

  // Initialize
  tp_pt  = 0;
  tp_eta = 0;
  tp_phi = 0;
  tp_dxy = 0;
  tp_z0  = 0;
  tp_d0  = 0;
  tp_pdgid   = 0;
  tp_nmatch  = 0;
  tp_nstub   = 0;
  tp_eventid = 0;
  tp_nloosematch = 0;

  matchtrk_pt  = 0;
  matchtrk_eta = 0;
  matchtrk_phi = 0;
  matchtrk_d0  = 0;
  matchtrk_z0  = 0;
  matchtrk_chi2  = 0;
  matchtrk_nstub = 0;

  loosematchtrk_pt  = 0;
  loosematchtrk_eta = 0;
  loosematchtrk_phi = 0;
  loosematchtrk_d0  = 0;
  loosematchtrk_z0  = 0;
  loosematchtrk_chi2   = 0; 
  loosematchtrk_nstub  = 0;

  trk_pt    = 0;
  trk_eta   = 0;
  trk_phi   = 0;
  trk_d0    = 0;
  trk_z0    = 0;
  trk_chi2  = 0;
  trk_nstub = 0;

  trk_combinatoric = 0;
  trk_genuine = 0;
  trk_loose   = 0;
  trk_unknown = 0;
  trk_fake    = 0;

  jet_eta    = 0;
  jet_pt     = 0;
  jet_phi    = 0;
  jet_tp_sumpt       = 0;
  jet_trk_sumpt      = 0;
  jet_matchtrk_sumpt = 0;
  jet_loosematchtrk_sumpt = 0;

  // Set Branch Addresses
  tree->SetBranchAddress("tp_pt",     &tp_pt,     &b_tp_pt);
  tree->SetBranchAddress("tp_eta",    &tp_eta,    &b_tp_eta);
  tree->SetBranchAddress("tp_phi",    &tp_phi,    &b_tp_phi);
  tree->SetBranchAddress("tp_dxy",    &tp_dxy,    &b_tp_dxy);
  tree->SetBranchAddress("tp_z0",     &tp_z0,     &b_tp_z0);
  tree->SetBranchAddress("tp_d0",     &tp_d0,     &b_tp_d0);
  tree->SetBranchAddress("tp_pdgid",  &tp_pdgid,  &b_tp_pdgid);
  tree->SetBranchAddress("tp_nstub",  &tp_nstub,  &b_tp_nstub);
  tree->SetBranchAddress("tp_eventid",&tp_eventid,&b_tp_eventid);
  if (doLooseMatch) tree->SetBranchAddress("tp_nloosematch",&tp_nmatch, &b_tp_nmatch);
  else              tree->SetBranchAddress("tp_nmatch",     &tp_nmatch, &b_tp_nmatch);

  if (doLooseMatch) {
    tree->SetBranchAddress("loosematchtrk_pt",    &matchtrk_pt,    &b_matchtrk_pt);
    tree->SetBranchAddress("loosematchtrk_eta",   &matchtrk_eta,   &b_matchtrk_eta);
    tree->SetBranchAddress("loosematchtrk_phi",   &matchtrk_phi,   &b_matchtrk_phi);
    tree->SetBranchAddress("loosematchtrk_d0",    &matchtrk_d0,    &b_matchtrk_d0);
    tree->SetBranchAddress("loosematchtrk_z0",    &matchtrk_z0,    &b_matchtrk_z0);
    tree->SetBranchAddress("loosematchtrk_chi2",  &matchtrk_chi2,  &b_matchtrk_chi2);
    tree->SetBranchAddress("loosematchtrk_nstub", &matchtrk_nstub, &b_matchtrk_nstub);
  }
  else {
    tree->SetBranchAddress("matchtrk_pt",    &matchtrk_pt,    &b_matchtrk_pt);
    tree->SetBranchAddress("matchtrk_eta",   &matchtrk_eta,   &b_matchtrk_eta);
    tree->SetBranchAddress("matchtrk_phi",   &matchtrk_phi,   &b_matchtrk_phi);
    tree->SetBranchAddress("matchtrk_d0",    &matchtrk_d0,    &b_matchtrk_d0);
    tree->SetBranchAddress("matchtrk_z0",    &matchtrk_z0,    &b_matchtrk_z0);
    tree->SetBranchAddress("matchtrk_chi2",  &matchtrk_chi2,  &b_matchtrk_chi2);
    tree->SetBranchAddress("matchtrk_nstub", &matchtrk_nstub, &b_matchtrk_nstub);
  }

  tree->SetBranchAddress("trk_pt",    &trk_pt,   &b_trk_pt);
  tree->SetBranchAddress("trk_eta",   &trk_eta,  &b_trk_eta);
  tree->SetBranchAddress("trk_phi",   &trk_phi,  &b_trk_phi);
  tree->SetBranchAddress("trk_d0",    &trk_d0,   &b_trk_d0);
  tree->SetBranchAddress("trk_z0",    &trk_z0,   &b_trk_z0);
  tree->SetBranchAddress("trk_chi2",  &trk_chi2, &b_trk_chi2);
  tree->SetBranchAddress("trk_nstub", &trk_nstub,&b_trk_nstub);

  tree->SetBranchAddress("jet_eta", &jet_eta,  &b_jet_eta);
  tree->SetBranchAddress("jet_pt",  &jet_pt,   &b_jet_pt);
  tree->SetBranchAddress("jet_phi", &jet_phi,  &b_jet_phi);
  tree->SetBranchAddress("jet_tp_sumpt",       &jet_tp_sumpt,       &b_jet_tp_sumpt);
  tree->SetBranchAddress("jet_trk_sumpt",      &jet_trk_sumpt,      &b_jet_trk_sumpt);
  if (doLooseMatch) tree->SetBranchAddress("jet_loosematchtrk_sumpt", &jet_matchtrk_sumpt, &b_jet_matchtrk_sumpt);
  else              tree->SetBranchAddress("jet_matchtrk_sumpt",      &jet_matchtrk_sumpt, &b_jet_matchtrk_sumpt);


  

  // ---------------------------------------------------------------------------------------------------------
  // HISTOGRAMS
  // ---------------------------------------------------------------------------------------------------------

  /////////////////////////////////////////////////
  // NOTATION:                                   //
  // 'C' - Central eta range, |eta|<0.8          //
  // 'I' - Intermediate eta range, 0.8<|eta|<1.6 //
  // 'F' - Forward eta range, |eta|>1.6          //
  //                                             //
  // 'L' - Low pt range,  pt<8 GeV               //
  // 'H' - High pt range, pt>8 GeV               //
  /////////////////////////////////////////////////


  // ---------------------------------------------------------------------------------------------------------
  // Tracking Particles
  TH1F* h_tp_p     = new TH1F("tp_p",    ";Tracking particle p [GeV]; Tracking particles",     200, 0.0, 200.0);
  TH1F* h_tp_pt    = new TH1F("tp_pt",   ";Tracking particle p_{T} [GeV]; Tracking particles", 100, 0.0, 100.0);
  TH1F* h_tp_pt_L  = new TH1F("tp_pt_L", ";Tracking particle p_{T} [GeV]; Tracking particles",  80, 0.0,   8.0);
  TH1F* h_tp_pt_LC = new TH1F("tp_pt_LC",";Tracking particle p_{T} [GeV]; Tracking particles",  80, 0.0,   8.0);
  TH1F* h_tp_pt_H  = new TH1F("tp_pt_H", ";Tracking particle p_{T} [GeV]; Tracking particles",  92, 8.0, 100.0);

  TH1F* h_tp_eta   = new TH1F("tp_eta",   ";Tracking particle #eta; Tracking particles", 50, -2.5, 2.5);
  TH1F* h_tp_eta_L = new TH1F("tp_eta_L", ";Tracking particle #eta; Tracking particles", 50, -2.5, 2.5);
  TH1F* h_tp_eta_H = new TH1F("tp_eta_H", ";Tracking particle #eta; Tracking particles", 50, -2.5, 2.5);
  TH1F* h_tp_eta_23= new TH1F("tp_eta_23",";Tracking particle #eta; Tracking particles", 50, -2.5, 2.5);
  TH1F* h_tp_eta_35= new TH1F("tp_eta_35",";Tracking particle #eta; Tracking particles", 50, -2.5, 2.5);
  TH1F* h_tp_eta_5 = new TH1F("tp_eta_5", ";Tracking particle #eta; Tracking particles", 50, -2.5, 2.5);

  TH1F* h_tp_phi  = new TH1F("tp_phi",  ";Tracking particle #phi [rad];   Tracking particles",  64, -3.20, 3.20);
  TH1F* h_tp_d0   = new TH1F("tp_d0",   ";Tracking particle d_{0} [cm];   Tracking particles", 100, -0.02, 0.02);
  TH1F* h_tp_absd0= new TH1F("tp_absd0",";Tracking particle |d_{0}| [cm]; Tracking particles",  50,  0.00, 2.00);
  TH1F* h_tp_z0   = new TH1F("tp_z0",   ";Tracking particle z_{0} [cm]; Tracking particles",    50, -25.0, 25.0);
  TH1F* h_tp_z0_L = new TH1F("tp_z0_L", ";Tracking particle z_{0} [cm]; Tracking particles",    50, -25.0, 25.0);
  TH1F* h_tp_z0_H = new TH1F("tp_z0_H", ";Tracking particle z_{0} [cm]; Tracking particles",    50, -25.0, 25.0);

  // Matched TPs
  TH1F* h_match_tp_p     = new TH1F("match_tp_p",     ";Tracking particle p [GeV]; Tracking particles",    200, 0.0, 200.0);
  TH1F* h_match_tp_pt    = new TH1F("match_tp_pt",    ";Tracking particle p_{T} [GeV]; Tracking particles",100, 0.0, 100.0);
  TH1F* h_match_tp_pt_L  = new TH1F("match_tp_pt_L",  ";Tracking particle p_{T} [GeV]; Tracking particles", 80, 0.0,   8.0);
  TH1F* h_match_tp_pt_LC = new TH1F("match_tp_pt_LC", ";Tracking particle p_{T} [GeV]; Tracking particles", 80, 0.0,   8.0);
  TH1F* h_match_tp_pt_H  = new TH1F("match_tp_pt_H",  ";Tracking particle p_{T} [GeV]; Tracking particles", 92, 8.0, 100.0);

  TH1F* h_match_tp_eta   = new TH1F("match_tp_eta",   ";Tracking particle #eta; Tracking particles", 50, -2.5, 2.5);
  TH1F* h_match_tp_eta_L = new TH1F("match_tp_eta_L", ";Tracking particle #eta; Tracking particles", 50, -2.5, 2.5);
  TH1F* h_match_tp_eta_H = new TH1F("match_tp_eta_H", ";Tracking particle #eta; Tracking particles", 50, -2.5, 2.5);
  TH1F* h_match_tp_eta_23= new TH1F("match_tp_eta_23",";Tracking particle #eta; Tracking particles", 50, -2.5, 2.5);
  TH1F* h_match_tp_eta_35= new TH1F("match_tp_eta_35",";Tracking particle #eta; Tracking particles", 50, -2.5, 2.5);
  TH1F* h_match_tp_eta_5 = new TH1F("match_tp_eta_5", ";Tracking particle #eta; Tracking particles", 50, -2.5, 2.5);

  TH1F* h_match_tp_phi  = new TH1F("match_tp_phi",  ";Tracking particle #phi [rad]; Tracking particles",  64, -3.20, 3.20);
  TH1F* h_match_tp_d0   = new TH1F("match_tp_d0",   ";Tracking particle d_{0} [cm]; Tracking particles", 100, -0.02, 0.02);
  TH1F* h_match_tp_absd0= new TH1F("match_tp_absd0",";Tracking particle d_{0} [cm]; Tracking particles",  50,  0.00, 2.00);
  TH1F* h_match_tp_z0   = new TH1F("match_tp_z0",   ";Tracking particle z_{0} [cm]; Tracking particles",  50, -25.0, 25.0);
  TH1F* h_match_tp_z0_L = new TH1F("match_tp_z0_L", ";Tracking particle z_{0} [cm]; Tracking particles",  50, -25.0, 25.0);
  TH1F* h_match_tp_z0_H = new TH1F("match_tp_z0_H", ";Tracking particle z_{0} [cm]; Tracking particles",  50, -25.0, 25.0);


  // ---------------------------------------------------------------------------------------------------------
  // L1 Tracks
  TH1F* h_trk_p     = new TH1F("trk_p",     ";L1 Track p [GeV]; L1 Tracks",     200, 0.0, 200.0);
  TH1F* h_trk_pt    = new TH1F("trk_pt",    ";L1 Track p_{T} [GeV]; L1 Tracks", 200, 0.0, 100.0);
  TH1F* h_trk_pt_L  = new TH1F("trk_pt_L",  ";L1 Track p_{T} [GeV]; L1 Tracks", 160, 0.0,   8.0);
  TH1F* h_trk_pt_LC = new TH1F("trk_pt_LC", ";L1 Track p_{T} [GeV]; L1 Tracks", 160, 0.0,   8.0);
  TH1F* h_trk_pt_H  = new TH1F("trk_pt_H",  ";L1 Track p_{T} [GeV]; L1 Tracks", 184, 8.0, 100.0);

  TH1F* h_trk_eta   = new TH1F("trk_eta",   ";L1 Track #eta; L1 Tracks", 50, -2.5, 2.5);
  TH1F* h_trk_eta_L = new TH1F("trk_eta_L", ";L1 Track #eta; L1 Tracks", 50, -2.5, 2.5);
  TH1F* h_trk_eta_H = new TH1F("trk_eta_H", ";L1 Track #eta; L1 Tracks", 50, -2.5, 2.5);
  TH1F* h_trk_eta_23= new TH1F("trk_eta_23",";L1 Track #eta; L1 Tracks / 0.1", 50, -2.5, 2.5);
  TH1F* h_trk_eta_35= new TH1F("trk_eta_35",";L1 Track #eta; L1 Tracks / 0.1", 50, -2.5, 2.5);
  TH1F* h_trk_eta_5 = new TH1F("trk_eta_5", ";L1 Track #eta; L1 Tracks / 0.1", 50, -2.5, 2.5);

  TH1F* h_trk_phi  = new TH1F("trk_phi",  ";L1 Track #phi [rad]; L1 Tracks",  64, -3.20, 3.20);
  TH1F* h_trk_d0   = new TH1F("trk_d0",   ";L1 Track d_{0} [cm]; L1 Tracks", 100, -0.02, 0.02);
  TH1F* h_trk_absd0= new TH1F("trk_absd0",";L1 Track d_{0} [cm]; L1 Tracks",  50,  0.00, 2.00);
  TH1F* h_trk_z0   = new TH1F("trk_z0",   ";L1 Track z_{0} [cm]; L1 Tracks",  50, -25.0, 25.0);
  TH1F* h_trk_z0_L = new TH1F("trk_z0_L", ";L1 Track z_{0} [cm]; L1 Tracks",  50, -25.0, 25.0);
  TH1F* h_trk_z0_H = new TH1F("trk_z0_H", ";L1 Track z_{0} [cm]; L1 Tracks",  50, -25.0, 25.0);

  // Matched Tracks
  TH1F* h_matchtrk_p     = new TH1F("matchtrk_p",     ";L1 Track p [GeV]; L1 Tracks",     200, 0.0, 200.0);
  TH1F* h_matchtrk_pt    = new TH1F("matchtrk_pt",    ";L1 Track p_{T} [GeV]; L1 Tracks", 200, 0.0, 100.0);
  TH1F* h_matchtrk_pt_L  = new TH1F("matchtrk_pt_L",  ";L1 Track p_{T} [GeV]; L1 Tracks", 160, 0.0,   8.0);
  TH1F* h_matchtrk_pt_LC = new TH1F("matchtrk_pt_LC", ";L1 Track p_{T} [GeV]; L1 Tracks", 160, 0.0,   8.0);
  TH1F* h_matchtrk_pt_H  = new TH1F("matchtrk_pt_H",  ";L1 Track p_{T} [GeV]; L1 Tracks", 184, 8.0, 100.0);

  TH1F* h_matchtrk_eta   = new TH1F("matchtrk_eta",   ";L1 Track #eta; L1 Tracks", 50, -2.5, 2.5);
  TH1F* h_matchtrk_eta_L = new TH1F("matchtrk_eta_L", ";L1 Track #eta; L1 Tracks", 50, -2.5, 2.5);
  TH1F* h_matchtrk_eta_H = new TH1F("matchtrk_eta_H", ";L1 Track #eta; L1 Tracks", 50, -2.5, 2.5);
  TH1F* h_matchtrk_eta_23= new TH1F("matchtrk_eta_23",";L1 Track #eta; L1 Tracks", 50, -2.5, 2.5);
  TH1F* h_matchtrk_eta_35= new TH1F("matchtrk_eta_35",";L1 Track #eta; L1 Tracks", 50, -2.5, 2.5);
  TH1F* h_matchtrk_eta_5 = new TH1F("matchtrk_eta_5", ";L1 Track #eta; L1 Tracks", 50, -2.5, 2.5);

  TH1F* h_matchtrk_phi  = new TH1F("matchtrk_phi",  ";L1 Track #phi [rad]; L1 Tracks",  64, -3.20, 3.20);
  TH1F* h_matchtrk_d0   = new TH1F("matchtrk_d0",   ";L1 Track d_{0} [cm]; L1 Tracks", 100, -0.02, 0.02);
  TH1F* h_matchtrk_absd0= new TH1F("matchtrk_absd0",";L1 Track d_{0} [cm]; L1 Tracks",  50,  0.00, 2.00);
  TH1F* h_matchtrk_z0   = new TH1F("matchtrk_z0",   ";L1 Track z_{0} [cm]; L1 Tracks",  50, -25.0, 25.0);
  TH1F* h_matchtrk_z0_L = new TH1F("matchtrk_z0_L", ";L1 Track z_{0} [cm]; L1 Tracks",  50, -25.0, 25.0);
  TH1F* h_matchtrk_z0_H = new TH1F("matchtrk_z0_H", ";L1 Track z_{0} [cm]; L1 Tracks",  50, -25.0, 25.0);


  // ---------------------------------------------------------------------------------------------------------
  // Jet Plots
  TH1F* h_jet_p    = new TH1F("jet_p",    ";Jet p [GeV]; GenJets",     200, 0.0, 200.0);
  TH1F* h_jet_pt   = new TH1F("jet_pt",   ";Jet p_{T} [GeV]; GenJets", 200, 0.0, 100.0);
  TH1F* h_jet_eta  = new TH1F("jet_eta",  ";Jet #eta; GenJets",         50, -2.5, 2.5);
  TH1F* h_jet_phi  = new TH1F("jet_phi",  ";Jet #phi [rad]; GenJets",   64, -3.2, 3.2);

  //-----Efficiency around Jet Axis -----//

  //----- Tracking Particles
  
  // TPs
  TH1F* h_deta_tp   = new TH1F("deta_tp",   ";Tracking particle distance to jet axis in #eta; Tracking particles", 250, -2.5, 2.5);
  TH1F* h_deta_tp_L = new TH1F("deta_tp_L", ";Tracking particle distance to jet axis in #eta; Tracking particles", 250, -2.5, 2.5);
  TH1F* h_deta_tp_H = new TH1F("deta_tp_H", ";Tracking particle distance to jet axis in #eta; Tracking particles", 250, -2.5, 2.5);

  TH1F* h_dphi_tp   = new TH1F("dphi_tp",   ";Tracking particle distance to jet axis in #phi [rad]; Tracking particles", 320, -3.2, 3.2);
  TH1F* h_dphi_tp_L = new TH1F("dphi_tp_L", ";Tracking particle distance to jet axis in #phi [rad]; Tracking particles", 320, -3.2, 3.2);
  TH1F* h_dphi_tp_H = new TH1F("dphi_tp_H", ";Tracking particle distance to jet axis in #phi [rad]; Tracking particles", 320, -3.2, 3.2);

  TH1F* h_dR_tp   = new TH1F("dR_tp",   ";Tracking particle distance to jet axis in dR; Tracking particles", 100, 0.0, 1.0);
  TH1F* h_dR_tp_L = new TH1F("dR_tp_L", ";Tracking particle distance to jet axis in dR; Tracking particles", 100, 0.0, 1.0);
  TH1F* h_dR_tp_H = new TH1F("dR_tp_H", ";Tracking particle distance to jet axis in dR; Tracking particles", 100, 0.0, 1.0);

  TH1F* h_dpt_tp   = new TH1F("dpt_tp",   ";Tracking Particle in dR<0.4 around jet axis, p_{T} [GeV]; Tracking particles", 300, 0.0, 100.0);
  TH1F* h_dpt_tp_L = new TH1F("dpt_tp_L", ";Tracking Particle in dR<0.4 around jet axis, p_{T} [GeV]; Tracking particles", 240, 0.0, 8.0);
  TH1F* h_dpt_tp_H = new TH1F("dpt_tp_H", ";Tracking Particle in dR<0.4 around jet axis, p_{T} [GeV]; Tracking particles", 276, 8.0, 100.0);

  TH1F* h_dpt_jet   = new TH1F("dpt_jet",   ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 500, 0.0, 1000.0);
  TH1F* h_dpt_jet_C = new TH1F("dpt_jet_C", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 500, 0.0, 1000.0);
  TH1F* h_dpt_jet_I = new TH1F("dpt_jet_I", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 500, 0.0, 1000.0);
  TH1F* h_dpt_jet_F = new TH1F("dpt_jet_F", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 500, 0.0, 1000.0);

  TH1F* h_dpt_jet_L  = new TH1F("dpt_jet_L",  ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles",  150, 25.0, 100.0);
  TH1F* h_dpt_jet_LC = new TH1F("dpt_jet_LC", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles",  150, 25.0, 100.0);
  TH1F* h_dpt_jet_LI = new TH1F("dpt_jet_LI", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles",  150, 25.0, 100.0);
  TH1F* h_dpt_jet_LF = new TH1F("dpt_jet_LF", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles",  150, 25.0, 100.0);

  TH1F* h_dpt_jet_H  = new TH1F("dpt_jet_H",  ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles",  225, 100.0, 1000.0);
  TH1F* h_dpt_jet_HC = new TH1F("dpt_jet_HC", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles",  225, 100.0, 1000.0);
  TH1F* h_dpt_jet_HI = new TH1F("dpt_jet_HI", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles",  225, 100.0, 1000.0);
  TH1F* h_dpt_jet_HF = new TH1F("dpt_jet_HF", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles",  225, 100.0, 1000.0);

  TH1F* h_dpt_jet_full  = new TH1F("dpt_jet_full",  ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 875, 0.0, 7000.0);
  TH1F* h_dpt_jet_full_C = new TH1F("dpt_jet_fullC", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 875, 0.0, 7000.0);
  TH1F* h_dpt_jet_full_I = new TH1F("dpt_jet_fullI", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 875, 0.0, 7000.0);
  TH1F* h_dpt_jet_full_F = new TH1F("dpt_jet_fullF", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 875, 0.0, 7000.0);


  // Matched TPs around Jet Axis
  TH1F* h_deta_match_tp   = new TH1F("deta_match_tp",   ";Tracking particle distance to jet axis in #eta; Tracking particles", 250, -2.5, 2.5);
  TH1F* h_deta_match_tp_L = new TH1F("deta_match_tp_L", ";Tracking particle distance to jet axis in #eta; Tracking particles", 250, -2.5, 2.5);
  TH1F* h_deta_match_tp_H = new TH1F("deta_match_tp_H", ";Tracking particle distance to jet axis in #eta; Tracking particles", 250, -2.5, 2.5);

  TH1F* h_dphi_match_tp   = new TH1F("dphi_match_tp",   ";Tracking particle distance to jet axis in #phi [rad]; Tracking particles", 320,  -3.2,  3.2);
  TH1F* h_dphi_match_tp_L = new TH1F("dphi_match_tp_L", ";Tracking particle distance to jet axis in #phi [rad]; Tracking particles", 320,  -3.2,  3.2);
  TH1F* h_dphi_match_tp_H = new TH1F("dphi_match_tp_H", ";Tracking particle distance to jet axis in #phi [rad]; Tracking particles", 320,  -3.2,  3.2);

  TH1F* h_dR_match_tp   = new TH1F("dR_match_tp",   ";Tracking particle distance to jet axis in dR; Tracking particles", 100, 0.0, 1.0);
  TH1F* h_dR_match_tp_L = new TH1F("dR_match_tp_L", ";Tracking particle distance to jet axis in dR; Tracking particles", 100, 0.0, 1.0);
  TH1F* h_dR_match_tp_H = new TH1F("dR_match_tp_H", ";Tracking particle distance to jet axis in dR; Tracking particles", 100, 0.0, 1.0);

  TH1F* h_dpt_match_tp   = new TH1F("dpt_match_tp",   ";Tracking Particle in dR<0.4 around jet axis, p_{T} [GeV]; Tracking particles", 300, 0.0, 100.0);
  TH1F* h_dpt_match_tp_L = new TH1F("dpt_match_tp_L", ";Tracking Particle in dR<0.4 around jet axis, p_{T} [GeV]; Tracking particles", 240, 0.0, 8.0);
  TH1F* h_dpt_match_tp_H = new TH1F("dpt_match_tp_H", ";Tracking Particle in dR<0.4 around jet axis, p_{T} [GeV]; Tracking particles", 276, 8.0, 100.0);

  TH1F* h_dpt_match_jet   = new TH1F("dpt_match_jet",   ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 500, 0.0, 1000.0);
  TH1F* h_dpt_match_jet_C = new TH1F("dpt_match_jet_C", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 500, 0.0, 1000.0);
  TH1F* h_dpt_match_jet_I = new TH1F("dpt_match_jet_I", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 500, 0.0, 1000.0);
  TH1F* h_dpt_match_jet_F = new TH1F("dpt_match_jet_F", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 500, 0.0, 1000.0);

  TH1F* h_dpt_match_jet_L  = new TH1F("dpt_match_jet_L",  ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 150, 25.0, 100.0);
  TH1F* h_dpt_match_jet_LC = new TH1F("dpt_match_jet_LC", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 150, 25.0, 100.0);
  TH1F* h_dpt_match_jet_LI = new TH1F("dpt_match_jet_LI", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 150, 25.0, 100.0);
  TH1F* h_dpt_match_jet_LF = new TH1F("dpt_match_jet_LF", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 150, 25.0, 100.0);

  TH1F* h_dpt_match_jet_H  = new TH1F("dpt_match_jet_H",  ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 225, 100.0, 1000.0);
  TH1F* h_dpt_match_jet_HC = new TH1F("dpt_match_jet_HC", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 225, 100.0, 1000.0);
  TH1F* h_dpt_match_jet_HI = new TH1F("dpt_match_jet_HI", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 225, 100.0, 1000.0);
  TH1F* h_dpt_match_jet_HF = new TH1F("dpt_match_jet_HF", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 225, 100.0, 1000.0);

  TH1F* h_dpt_match_jet_full  = new TH1F("dpt_match_jet_full",  ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 875, 0.0, 7000.0);
  TH1F* h_dpt_match_jet_full_C = new TH1F("dpt_match_jet_fullC", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 875, 0.0, 7000.0);
  TH1F* h_dpt_match_jet_full_I = new TH1F("dpt_match_jet_fullI", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 875, 0.0, 7000.0);
  TH1F* h_dpt_match_jet_full_F = new TH1F("dpt_match_jet_fullF", ";Tracking Particle in dR<0.4 around jet axis, Jet p_{T} [GeV]; Tracking particles", 875, 0.0, 7000.0);


  //----- L1 Tracks

  // L1 Tracks
  TH1F* h_deta_trk   = new TH1F("deta_trk",   ";L1 Track distance to jet axis in #eta; L1 Tracks", 250, -2.5, 2.5);
  TH1F* h_deta_trk_L = new TH1F("deta_trk_L", ";L1 Track distance to jet axis in #eta; L1 Tracks", 250, -2.5, 2.5);
  TH1F* h_deta_trk_H = new TH1F("deta_trk_H", ";L1 Track distance to jet axis in #eta; L1 Tracks", 250, -2.5, 2.5);

  TH1F* h_dphi_trk   = new TH1F("dphi_trk",   ";L1 Track distance to jet axis in #phi [rad]; L1 Tracks", 320, -3.2, 3.2);
  TH1F* h_dphi_trk_L = new TH1F("dphi_trk_L", ";L1 Track distance to jet axis in #phi [rad]; L1 Tracks", 320, -3.2, 3.2);
  TH1F* h_dphi_trk_H = new TH1F("dphi_trk_H", ";L1 Track distance to jet axis in #phi [rad]; L1 Tracks", 320, -3.2, 3.2);

  TH1F* h_dR_trk   = new TH1F("dR_trk",   ";L1 Tracks distance to jet axis in dR; Tracking particles", 100, 0.0, 1.0);
  TH1F* h_dR_trk_L = new TH1F("dR_trk_L", ";L1 Tracks distance to jet axis in dR; Tracking particles", 100, 0.0, 1.0);
  TH1F* h_dR_trk_H = new TH1F("dR_trk_H", ";L1 Tracks distance to jet axis in dR; Tracking particles", 100, 0.0, 1.0);

  TH1F* h_dpt_trk   = new TH1F("dpt_trk",   ";L1 Tracks in dR<0.4 around jet axis, p_{T} [GeV]; L1 Tracks", 300, 0.0, 100.0);
  TH1F* h_dpt_trk_L = new TH1F("dpt_trk_L", ";L1 Tracks in dR<0.4 around jet axis, p_{T} [GeV]; L1 Tracks", 240, 0.0,   8.0);
  TH1F* h_dpt_trk_H = new TH1F("dpt_trk_H", ";L1 Tracks in dR<0.4 around jet axis, p_{T} [GeV]; L1 Tracks", 276, 8.0, 100.0);

  TH1F* h_dpt_trk_jet      = new TH1F("dpt_trk_jet",     ";L1 Tracks in dR<0.4 around jet axis, Jet p_{T} [GeV]; L1 Tracks",  500,   0.0, 1000.0);
  TH1F* h_dpt_trk_jet_L    = new TH1F("dpt_jet_trk_L",   ";L1 Tracks in dR<0.4 around jet axis, Jet p_{T} [GeV]; L1 Tracks",  150,  25.0,  100.0);
  TH1F* h_dpt_trk_jet_H    = new TH1F("dpt_jet_trk_H",   ";L1 Tracks in dR<0.4 around jet axis, Jet p_{T} [GeV]; L1 Tracks",  225, 100.0, 1000.0);
  TH1F* h_dpt_trk_jet_full = new TH1F("dpt_trk_jet_full",";L1 Tracks in dR<0.4 around jet axis, Jet p_{T} [GeV]; L1 Tracks", 875,   0.0, 7000.0);


  // L1 Matched Tracks
  TH1F* h_deta_matchtrk   = new TH1F("deta_matchtrk",   ";L1 Track distance to jet axis in #eta; L1 Tracks", 250, -2.5, 2.5);
  TH1F* h_deta_matchtrk_L = new TH1F("deta_matchtrk_L", ";L1 Track distance to jet axis in #eta; L1 Tracks", 250, -2.5, 2.5);
  TH1F* h_deta_matchtrk_H = new TH1F("deta_matchtrk_H", ";L1 Track distance to jet axis in #eta; L1 Tracks", 250, -2.5, 2.5);

  TH1F* h_dphi_matchtrk   = new TH1F("dphi_matchtrk",   ";L1 Track distance to jet axis in #phi [rad]; L1 Tracks", 320,  -3.2,  3.2);
  TH1F* h_dphi_matchtrk_L = new TH1F("dphi_matchtrk_L", ";L1 Track distance to jet axis in #phi [rad]; L1 Tracks", 320,  -3.2,  3.2);
  TH1F* h_dphi_matchtrk_H = new TH1F("dphi_matchtrk_H", ";L1 Track distance to jet axis in #phi [rad]; L1 Tracks", 320,  -3.2,  3.2);

  TH1F* h_dR_matchtrk   = new TH1F("dR_matchtrk",   ";L1 Tracks distance to jet axis in dR; Tracking particles", 100, 0.0, 1.0);
  TH1F* h_dR_matchtrk_L = new TH1F("dR_matchtrk_L", ";L1 Tracks distance to jet axis in dR; Tracking particles", 100, 0.0, 1.0);
  TH1F* h_dR_matchtrk_H = new TH1F("dR_matchtrk_H", ";L1 Tracks distance to jet axis in dR; Tracking particles", 100, 0.0, 1.0);

  TH1F* h_dpt_matchtrk   = new TH1F("dpt_matchtrk",   ";L1 Tracks in dR<0.4 around jet axis, p_{T} [GeV]; L1 Tracks", 300, 0.0, 100.0);
  TH1F* h_dpt_matchtrk_L = new TH1F("dpt_matchtrk_L", ";L1 Tracks in dR<0.4 around jet axis, p_{T} [GeV]; L1 Tracks", 240, 0.0,   8.0);
  TH1F* h_dpt_matchtrk_H = new TH1F("dpt_matchtrk_H", ";L1 Tracks in dR<0.4 around jet axis, p_{T} [GeV]; L1 Tracks", 276, 8.0, 100.0);

  TH1F* h_dpt_matchtrk_jet      = new TH1F("dpt_matchtrk_jet",     ";L1 Tracks in dR<0.4 around jet axis, Jet p_{T} [GeV]; L1 Tracks",  500,   0.0, 1000.0);
  TH1F* h_dpt_matchtrk_jet_L    = new TH1F("dpt_matchtrk_jet_L",   ";L1 Tracks in dR<0.4 around jet axis, Jet p_{T} [GeV]; L1 Tracks",  150,  25.0,  100.0);
  TH1F* h_dpt_matchtrk_jet_H    = new TH1F("dpt_matchtrk_jet_H",   ";L1 Tracks in dR<0.4 around jet axis, Jet p_{T} [GeV]; L1 Tracks",  225, 100.0, 1000.0);
  TH1F* h_dpt_matchtrk_jet_full = new TH1F("dpt_matchtrk_jet_full",";L1 Tracks in dR<0.4 around jet axis, Jet p_{T} [GeV]; L1 Tracks", 875,   0.0, 7000.0);


  // ----- Jet Profile Plots ----- //
  TProfile* h_jet_profile_TP = new TProfile("jet_profile_TP", ";Jet p_{T} [GeV]; <Tracking Particles per Jet>", 100, 0.0, 1000.0);
  TProfile* h_jet_profile_L1 = new TProfile("jet_profile_L1", ";Jet p_{T} [GeV]; <L1 Tracks per Jet>",             100, 0.0, 1000.0);
  TProfile* h_jet_profile_MTP = new TProfile("jet_profile_MTP", ";Jet p_{T} [GeV]; <Matched Tracking Particles per Jet>", 100, 0.0, 1000.0);
  TProfile* h_jet_profile_ML1 = new TProfile("jet_profile_ML1", ";Jet p_{T} [GeV]; <Matched L1 Tracks per Jet>",          100, 0.0, 1000.0); // duplicate of MTP
  // TProfile* h_jet_profile_looseMTP = new TProfile("jet_profile_looseMTP", ";Jet p_{T} [GeV]; <Loose Matched Tracking Particles per Jet>", 100, 0.0, 1000.0); 
  // TProfile* h_jet_profile_looseML1 = new TProfile("jet_profile_looseML1", ";Jet p_{T} [GeV]; <Loose Matched L1 Tracks  per Jet>", 100, 0.0, 1000.0); // duplicate of looseMTP


  //-----Sum pT in Jets-----//
  TH1F* h_jet_tp_sumpt_vspt   = new TH1F("jet_tp_sumpt_vspt",   ";sum(TP p_{T}) [GeV]; Gen jets", 20,  0.0, 200.0);
  TH1F* h_jet_tp_sumpt_vseta  = new TH1F("jet_tp_sumpt_vseta",  ";Gen jet #eta; Gen jets",        24, -2.4, 2.4);
  TH1F* h_jet_tp_sumpt_vsphi  = new TH1F("jet_tp_sumpt_vsphi",  ";Gen jet #phi; Gen jets",        24, -3.2, 3.2);


  TH1F* h_jet_trk_sumpt_vspt  = new TH1F("jet_trk_sumpt_vspt",  ";sum(TP p_{T}) [GeV]; Gen jets", 20,  0.0, 200.0);
  TH1F* h_jet_trk_sumpt_vseta = new TH1F("jet_trk_sumpt_vseta", ";Gen jet #eta; Gen jets",        24, -2.4, 2.4);
  TH1F* h_jet_trk_sumpt_vsphi = new TH1F("jet_trk_sumpt_vsphi", ";Gen jet #phi; Gen jets",        24, -3.2, 3.2);
  
  TH1F* h_jet_matchtrk_sumpt_vspt  = new TH1F("jet_matchtrk_sumpt_vspt",  ";sum(TP p_{T}) [GeV]; Gen jets", 20,  0.0, 200.0);
  TH1F* h_jet_matchtrk_sumpt_vseta = new TH1F("jet_matchtrk_sumpt_vseta", ";Gen jet #eta; Gen jets",        24, -2.4, 2.4);
  TH1F* h_jet_matchtrk_sumpt_vsphi = new TH1F("jet_matchtrk_sumpt_vsphi", ";Gen jet #phi; Gen jets",        24, -3.2, 3.2);

  h_jet_tp_sumpt_vspt  ->Sumw2();
  h_jet_tp_sumpt_vseta ->Sumw2();
  h_jet_tp_sumpt_vsphi ->Sumw2();

  h_jet_trk_sumpt_vspt  ->Sumw2();
  h_jet_trk_sumpt_vseta ->Sumw2();
  h_jet_trk_sumpt_vsphi ->Sumw2();

  h_jet_matchtrk_sumpt_vspt  ->Sumw2();
  h_jet_matchtrk_sumpt_vseta ->Sumw2();
  h_jet_matchtrk_sumpt_vsphi ->Sumw2();


  // ---------------------------------------------------------------------------------------------------------
  // Chi2 and Chi2/dof Histograms
  // Note: last bin is an overflow bin
  TH1F* h_match_trk_chi2     = new TH1F("match_trk_chi2",     ";#chi^{2}; L1 Tracks", 100, 0, 100);
  TH1F* h_match_trk_chi2_C_L = new TH1F("match_trk_chi2_C_L", ";#chi^{2}; L1 Tracks", 100, 0, 100);
  TH1F* h_match_trk_chi2_I_L = new TH1F("match_trk_chi2_I_L", ";#chi^{2}; L1 Tracks", 100, 0, 100);
  TH1F* h_match_trk_chi2_F_L = new TH1F("match_trk_chi2_F_L", ";#chi^{2}; L1 Tracks", 100, 0, 100);
  TH1F* h_match_trk_chi2_C_H = new TH1F("match_trk_chi2_C_H", ";#chi^{2}; L1 Tracks", 100, 0, 100);
  TH1F* h_match_trk_chi2_I_H = new TH1F("match_trk_chi2_I_H", ";#chi^{2}; L1 Tracks", 100, 0, 100);
  TH1F* h_match_trk_chi2_F_H = new TH1F("match_trk_chi2_F_H", ";#chi^{2}; L1 Tracks", 100, 0, 100);

  TH1F* h_match_trk_chi2_dof     = new TH1F("match_trk_chi2_dof",     ";#chi^{2} / D.O.F.; L1 Tracks", 100, 0, 20);
  TH1F* h_match_trk_chi2_dof_C_L = new TH1F("match_trk_chi2_dof_C_L", ";#chi^{2} / D.O.F.; L1 Tracks", 100, 0, 20);
  TH1F* h_match_trk_chi2_dof_I_L = new TH1F("match_trk_chi2_dof_I_L", ";#chi^{2} / D.O.F.; L1 Tracks", 100, 0, 20);
  TH1F* h_match_trk_chi2_dof_F_L = new TH1F("match_trk_chi2_dof_F_L", ";#chi^{2} / D.O.F.; L1 Tracks", 100, 0, 20);
  TH1F* h_match_trk_chi2_dof_C_H = new TH1F("match_trk_chi2_dof_C_H", ";#chi^{2} / D.O.F.; L1 Tracks", 100, 0, 20);
  TH1F* h_match_trk_chi2_dof_I_H = new TH1F("match_trk_chi2_dof_I_H", ";#chi^{2} / D.O.F.; L1 Tracks", 100, 0, 20);
  TH1F* h_match_trk_chi2_dof_F_H = new TH1F("match_trk_chi2_dof_F_H", ";#chi^{2} / D.O.F.; L1 Tracks", 100, 0, 20);


  // ---------------------------------------------------------------------------------------------------------
  // Total Track Rates
  TH1F* h_trk_vspt = new TH1F("trk_vspt", ";Track p_{T} [GeV]; Tracking particles", 40, 0, 20);
  TH1F* h_tp_vspt  = new TH1F("tp_vspt",  ";TP p_{T} [GeV]; Tracking particles",    40, 0, 20);


  // ---------------------------------------------------------------------------------------------------------
  // Stubs
  TH1F* h_trk_nstub   = new TH1F("trk_nstub",   ";Number of L1 Track stubs; L1 Tracks", 15, 0, 15);
  TH1F* h_trk_nstub_C = new TH1F("trk_nstub_C", ";Number of L1 Track stubs; L1 Tracks", 15, 0, 15);
  TH1F* h_trk_nstub_I = new TH1F("trk_nstub_I", ";Number of L1 Track stubs; L1 Tracks", 15, 0, 15);
  TH1F* h_trk_nstub_F = new TH1F("trk_nstub_F", ";Number of L1 Track stubs; L1 Tracks", 15, 0, 15);

  TH1F* h_tp_nstub   = new TH1F("tp_nstub",   ";Number of TP stubs; Tracking Particles", 15, 0, 15);
  TH1F* h_tp_nstub_C = new TH1F("tp_nstub_C", ";Number of TP stubs; Tracking Particles", 15, 0, 15);
  TH1F* h_tp_nstub_I = new TH1F("tp_nstub_I", ";Number of TP stubs; Tracking Particles", 15, 0, 15);
  TH1F* h_tp_nstub_F = new TH1F("tp_nstub_F", ";Number of TP stubs; Tracking Particles", 15, 0, 15);

  TH1F* h_match_trk_nstub   = new TH1F("match_trk_nstub",   ";Number of Matched stubs; L1 Tracks", 15, 0, 15);
  TH1F* h_match_trk_nstub_C = new TH1F("match_trk_nstub_C", ";Number of Matched stubs; L1 Tracks", 15, 0, 15);
  TH1F* h_match_trk_nstub_I = new TH1F("match_trk_nstub_I", ";Number of Matched stubs; L1 Tracks", 15, 0, 15);
  TH1F* h_match_trk_nstub_F = new TH1F("match_trk_nstub_F", ";Number of Matched stubs; L1 Tracks", 15, 0, 15);


  // ---------------------------------------------------------------------------------------------------------
  // Number of Tracks Per Event
  TH1F* h_ntrk_pt2  = new TH1F("ntrk_pt2",  ";# tracks (p_{T} > 2 GeV) / event; Events",  300, 0, 300.0);
  TH1F* h_ntrk_pt3  = new TH1F("ntrk_pt3",  ";# tracks (p_{T} > 3 GeV) / event; Events",  300, 0, 300.0);
  TH1F* h_ntrk_pt10 = new TH1F("ntrk_pt10", ";# tracks (p_{T} > 10 GeV) / event; Events", 100, 0, 100.0);



  // ---------------------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------------------
  // RESOLUTION

  // resolution vs. pt histograms

  // ----------------------------------------------
  // for new versions of resolution vs pt/eta plots
  unsigned int nBinsPtRes = 500;
  double maxPtRes = 30.;

  unsigned int nBinsPtRelRes = 1000;
  double maxPtRelRes = 10.;

  unsigned int nBinsEtaRes = 500;
  double maxEtaRes = 0.1;

  unsigned int nBinsPhiRes = 500;
  double maxPhiRes = 0.2;

  unsigned int nBinsZ0Res = 100;
  double maxZ0Res = 4.0;
  // ----------------------------------------------
  
  const int nRANGE = 20;
  TString ptrange[nRANGE] = {"0-5","5-10", "10-15","15-20","20-25","25-30","30-35","35-40","40-45","45-50","50-55",
             "55-60","60-65","65-70","70-75","75-80","80-85","85-90","90-95","95-100"};

  const int nRANGE_L = 12;
  TString ptrange_L[nRANGE] = {"2-2.5", "2.5-3", "3-3.5", "3.5-4", "4-4.5", "4.5-5", "5-5.5", "5.5-6", "6-6.5", "6.5-7", "7-7.5", "7.5-8" };

  TH1F* h_absResVsPt_pt[nRANGE];
  TH1F* h_absResVsPt_ptRel[nRANGE];
  TH1F* h_absResVsPt_z0[nRANGE];
  TH1F* h_absResVsPt_phi[nRANGE];
  TH1F* h_absResVsPt_eta[nRANGE];
  TH1F* h_absResVsPt_d0[nRANGE];

  TH1F* h_absResVsPt_pt_L[nRANGE_L];
  TH1F* h_absResVsPt_ptRel_L[nRANGE_L];
  TH1F* h_absResVsPt_z0_L[nRANGE_L];
  TH1F* h_absResVsPt_phi_L[nRANGE_L];
  TH1F* h_absResVsPt_eta_L[nRANGE_L];
  TH1F* h_absResVsPt_d0_L[nRANGE_L];

  for (int i=0; i<nRANGE; i++) {
    h_absResVsPt_pt[i]    = new TH1F("absResVsPt_pt_"+ptrange[i],   ";p_{T} residual (L1 - sim) [GeV]; L1 tracks / 0.1",   nBinsPtRes,    0, maxPtRes);
    h_absResVsPt_ptRel[i] = new TH1F("absResVsPt_ptRel_"+ptrange[i],";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02",nBinsPtRelRes, 0, maxPtRelRes);
    h_absResVsPt_z0[i]    = new TH1F("absResVsPt_z0_"+ptrange[i],   ";z_{0} residual (L1 - sim) [GeV]; L1 tracks / 0.1",   nBinsZ0Res,    0, maxZ0Res);
    h_absResVsPt_phi[i]   = new TH1F("absResVsPt_phi_"+ptrange[i],  ";#phi residual (L1 - sim) [GeV]; L1 tracks / 0.1",    nBinsPhiRes,   0, maxPhiRes);
    h_absResVsPt_eta[i]   = new TH1F("absResVsPt_eta_"+ptrange[i],  ";#eta residual (L1 - sim) [GeV]; L1 tracks / 0.1",    nBinsEtaRes,   0, maxEtaRes);
    h_absResVsPt_d0[i]    = new TH1F("absResVsPt_d0_"+ptrange[i],   ";d_{0}residual (L1 - sim) [GeV]; L1 tracks / 0.1",    100,           0, 0.02);
  }

  for (int i=0; i<nRANGE_L; i++) {
    h_absResVsPt_pt_L[i]    = new TH1F("absResVsPt_L_pt_"+ptrange_L[i],   ";p_{T} residual (L1 - sim) [GeV]; L1 tracks / 0.1",   nBinsPtRes,    0, maxPtRes);
    h_absResVsPt_ptRel_L[i] = new TH1F("absResVsPt_L_ptRel_"+ptrange_L[i],";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02",nBinsPtRelRes, 0, maxPtRelRes);
    h_absResVsPt_z0_L[i]    = new TH1F("absResVsPt_L_z0_"+ptrange_L[i],   ";z_{0} residual (L1 - sim) [GeV]; L1 tracks / 0.1",   nBinsZ0Res,    0, maxZ0Res);
    h_absResVsPt_phi_L[i]   = new TH1F("absResVsPt_L_phi_"+ptrange_L[i],  ";#phi residual (L1 - sim) [GeV]; L1 tracks / 0.1",    nBinsPhiRes,   0, maxPhiRes);
    h_absResVsPt_eta_L[i]   = new TH1F("absResVsPt_L_eta_"+ptrange_L[i],  ";#eta residual (L1 - sim) [GeV]; L1 tracks / 0.1",    nBinsEtaRes,   0, maxEtaRes);
    h_absResVsPt_d0_L[i]    = new TH1F("absResVsPt_L_d0_"+ptrange_L[i],   ";d_{0}residual (L1 - sim) [GeV]; L1 tracks / 0.1",    100,           0, 0.02);
  }

  // resolution vs. eta histograms
  /*
  const int nETARANGE = 24;
  TString etarange[nETARANGE] = {"0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0",
         "1.1","1.2","1.3","1.4","1.5","1.6","1.7","1.8","1.9","2.0",
         "2.1","2.2","2.3","2.4"};
  */
  const int nETARANGE = 12;
  TString etarange[nETARANGE] = {"0.2","0.4","0.6","0.8","1.0",
         "1.2","1.4","1.6","1.8","2.0",
         "2.2","2.4"};

  TH1F* h_absResVsEta_eta[nETARANGE];
  TH1F* h_absResVsEta_z0[nETARANGE];
  TH1F* h_absResVsEta_phi[nETARANGE];
  TH1F* h_absResVsEta_ptRel[nETARANGE];
  TH1F* h_absResVsEta_d0[nETARANGE];

  TH1F* h_absResVsEta_eta_L[nETARANGE];
  TH1F* h_absResVsEta_z0_L[nETARANGE];
  TH1F* h_absResVsEta_phi_L[nETARANGE];
  TH1F* h_absResVsEta_ptRel_L[nETARANGE];
  TH1F* h_absResVsEta_d0_L[nETARANGE];

  TH1F* h_absResVsEta_eta_H[nETARANGE];
  TH1F* h_absResVsEta_z0_H[nETARANGE];
  TH1F* h_absResVsEta_phi_H[nETARANGE];
  TH1F* h_absResVsEta_ptRel_H[nETARANGE];
  TH1F* h_absResVsEta_d0_H[nETARANGE];

  for (int i=0; i<nETARANGE; i++) {
    h_absResVsEta_eta[i]   = new TH1F("absResVsEta_eta_"+etarange[i],  ";#eta residual (L1 - sim) [GeV]; L1 tracks / 0.1",     nBinsEtaRes,   0, maxEtaRes);
    h_absResVsEta_z0[i]    = new TH1F("absResVsEta_z0_"+etarange[i],   ";|z_{0} residual (L1 - sim)| [cm]; L1 tracks / 0.01",  nBinsZ0Res,    0, maxZ0Res);
    h_absResVsEta_phi[i]   = new TH1F("absResVsEta_phi_"+etarange[i],  ";#phi residual (L1 - sim) [GeV]; L1 tracks / 0.1",     nBinsPhiRes,   0, maxPhiRes);    
    h_absResVsEta_ptRel[i] = new TH1F("absResVsEta_ptRel_"+etarange[i],";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02", nBinsPtRelRes, 0, maxPtRelRes);
    h_absResVsEta_d0[i]    = new TH1F("absResVsEta_d0_"+etarange[i],   ";d_{0}residual (L1 - sim) [GeV]; L1 tracks / 0.1",     100,           0, 0.02);

    h_absResVsEta_eta_L[i]   = new TH1F("absResVsEta_eta_L_"+etarange[i],  ";#eta residual (L1 - sim) [GeV]; L1 tracks / 0.1",     nBinsEtaRes,   0, maxEtaRes);
    h_absResVsEta_z0_L[i]    = new TH1F("absResVsEta_z0_L_"+etarange[i],   ";|z_{0} residual (L1 - sim)| [cm]; L1 tracks / 0.01",  nBinsZ0Res,    0, maxZ0Res);
    h_absResVsEta_phi_L[i]   = new TH1F("absResVsEta_phi_L_"+etarange[i],  ";#phi residual (L1 - sim) [GeV]; L1 tracks / 0.1",     nBinsPhiRes,   0, maxPhiRes);    
    h_absResVsEta_ptRel_L[i] = new TH1F("absResVsEta_ptRel_L_"+etarange[i],";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02", nBinsPtRelRes, 0, maxPtRelRes);
    h_absResVsEta_d0_L[i]    = new TH1F("absResVsEta_d0_L_"+etarange[i],   ";d_{0}residual (L1 - sim) [GeV]; L1 tracks / 0.1",     100,           0, 0.02);

    h_absResVsEta_eta_H[i]   = new TH1F("absResVsEta_eta_H_"+etarange[i],  ";#eta residual (L1 - sim) [GeV]; L1 tracks / 0.1",     nBinsEtaRes,   0, maxEtaRes);
    h_absResVsEta_z0_H[i]    = new TH1F("absResVsEta_z0_H_"+etarange[i],   ";|z_{0} residual (L1 - sim)| [cm]; L1 tracks / 0.01",  nBinsZ0Res,    0, maxZ0Res);
    h_absResVsEta_phi_H[i]   = new TH1F("absResVsEta_phi_H_"+etarange[i],  ";#phi residual (L1 - sim) [GeV]; L1 tracks / 0.1",     nBinsPhiRes,   0, maxPhiRes);    
    h_absResVsEta_ptRel_H[i] = new TH1F("absResVsEta_ptRel_H_"+etarange[i],";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02", nBinsPtRelRes, 0, maxPtRelRes);
    h_absResVsEta_d0_H[i]    = new TH1F("absResVsEta_d0_H_"+etarange[i],   ";d_{0}residual (L1 - sim) [GeV]; L1 tracks / 0.1",     100,           0, 0.02);
  }

  // resolution vs phi 

  const int nPHIRANGE = 32;
  TString phirange[nPHIRANGE] = {"-3.0","-2.8","-2.6","-2.4","-2.2",
           "-2.0","-1.8","-1.6","-1.4","-1.2",
           "-1.0","-0.8","-0.6","-0.4","-0.2",
           "0.0","0.2","0.4","0.6","0.8",
           "1.0","1.2","1.4","1.6","1.8",
           "2.0","2.2","2.4","2.6","2.8",
           "3.0","3.2"};

  TH1F* h_absResVsPhi_pt[nPHIRANGE];
  TH1F* h_absResVsPhi_ptRel[nPHIRANGE];

  for (int i=0; i<nPHIRANGE; i++) {
    h_absResVsPhi_pt[i]    = new TH1F("absResVsPt_pt_"+phirange[i],   ";p_{T} residual (L1 - sim) [GeV]; L1 tracks / 0.1",   nBinsPtRes,    0, maxPtRes);
    h_absResVsPhi_ptRel[i] = new TH1F("absResVsPt_ptRel_"+phirange[i],";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02",nBinsPtRelRes, 0, maxPtRelRes);
  }


  // ----------------------------------------------------------------------------------------------------------------
  // resolution histograms
  TH1F* h_res_pt    = new TH1F("res_pt",    ";p_{T} residual (L1 - sim) [GeV]; L1 tracks / 0.05",   200,-5.0,   5.0);
  TH1F* h_res_ptRel = new TH1F("res_ptRel", ";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.01", 200,-1.0,   1.0);
  TH1F* h_res_eta   = new TH1F("res_eta",   ";#eta residual (L1 - sim); L1 tracks / 0.0002",        100,-0.01,  0.01);
  TH1F* h_res_phi   = new TH1F("res_phi",   ";#phi residual (L1 - sim) [rad]; L1 tracks / 0.0001",  100,-0.005, 0.005);

  TH1F* h_res_z0    = new TH1F("res_z0",    ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100, -1.0, 1.0);
  TH1F* h_res_z0_C  = new TH1F("res_z0_C",  ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100, -1.0, 1.0);
  TH1F* h_res_z0_I  = new TH1F("res_z0_I",  ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100, -1.0, 1.0);
  TH1F* h_res_z0_F  = new TH1F("res_z0_F",  ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100, -1.0, 1.0);
  TH1F* h_res_z0_L  = new TH1F("res_z0_L",  ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100, -1.0, 1.0);
  TH1F* h_res_z0_H  = new TH1F("res_z0_H",  ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100, -1.0, 1.0);
  
  TH1F* h_res_z0_C_L = new TH1F("res_z0_C_L", ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100,(-1)*1.0, 1.0);
  TH1F* h_res_z0_I_L = new TH1F("res_z0_I_L", ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100,(-1)*1.0, 1.0);
  TH1F* h_res_z0_F_L = new TH1F("res_z0_F_L", ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100,(-1)*1.0, 1.0);
  TH1F* h_res_z0_C_H = new TH1F("res_z0_C_H", ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100,(-1)*1.0, 1.0);
  TH1F* h_res_z0_I_H = new TH1F("res_z0_I_H", ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100,(-1)*1.0, 1.0);
  TH1F* h_res_z0_F_H = new TH1F("res_z0_F_H", ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100,(-1)*1.0, 1.0);

  TH1F* h_res_d0 = new TH1F("res_d0", ";d_{0} residual (L1 - sim) [cm]; L1 tracks / 0.0002 cm", 200,-0.02,0.02);
  TH1F* h_res_d0_C = new TH1F("res_d0_C", ";d_{0} residual (L1 - sim) [cm]; L1 tracks / 0.0001 cm", 200,-0.05,0.05);
  TH1F* h_res_d0_I = new TH1F("res_d0_I", ";d_{0} residual (L1 - sim) [cm]; L1 tracks / 0.0001 cm", 200,-0.05,0.05);
  TH1F* h_res_d0_F = new TH1F("res_d0_F", ";d_{0} residual (L1 - sim) [cm]; L1 tracks / 0.0001 cm", 200,-0.05,0.05);
  TH1F* h_res_d0_L = new TH1F("res_d0_L", ";d_{0} residual (L1 - sim) [cm]; L1 tracks / 0.0001 cm", 200,-0.05,0.05);
  TH1F* h_res_d0_H = new TH1F("res_d0_H", ";d_{0} residual (L1 - sim) [cm]; L1 tracks / 0.0001 cm", 200,-0.05,0.05);

  TH1F* h_res_d0_C_L = new TH1F("res_d0_C_L", ";d_{0} residual (L1 - sim) [cm]; L1 tracks / 0.0001 cm", 200,-0.05,0.05);
  TH1F* h_res_d0_I_L = new TH1F("res_d0_I_L", ";d_{0} residual (L1 - sim) [cm]; L1 tracks / 0.0001 cm", 200,-0.05,0.05);
  TH1F* h_res_d0_F_L = new TH1F("res_d0_F_L", ";d_{0} residual (L1 - sim) [cm]; L1 tracks / 0.0001 cm", 200,-0.05,0.05);
  TH1F* h_res_d0_C_H = new TH1F("res_d0_C_H", ";d_{0} residual (L1 - sim) [cm]; L1 tracks / 0.0001 cm", 200,-0.05,0.05);
  TH1F* h_res_d0_I_H = new TH1F("res_d0_I_H", ";d_{0} residual (L1 - sim) [cm]; L1 tracks / 0.0001 cm", 200,-0.05,0.05);
  TH1F* h_res_d0_F_H = new TH1F("res_d0_F_H", ";d_{0} residual (L1 - sim) [cm]; L1 tracks / 0.0001 cm", 200,-0.05,0.05);

  // ----------------------------------------------------------------------------------------------------------------
  // more resolution vs pt 

  TH1F* h_resVsPt_pt[nRANGE];
  TH1F* h_resVsPt_pt_C[nRANGE];
  TH1F* h_resVsPt_pt_I[nRANGE];
  TH1F* h_resVsPt_pt_F[nRANGE];

  TH1F* h_resVsPt_ptRel[nRANGE];
  TH1F* h_resVsPt_ptRel_C[nRANGE];
  TH1F* h_resVsPt_ptRel_I[nRANGE];
  TH1F* h_resVsPt_ptRel_F[nRANGE];

  TH1F* h_resVsPt_z0[nRANGE];
  TH1F* h_resVsPt_z0_C[nRANGE];
  TH1F* h_resVsPt_z0_I[nRANGE];
  TH1F* h_resVsPt_z0_F[nRANGE];

  TH1F* h_resVsPt_phi[nRANGE];
  TH1F* h_resVsPt_phi_C[nRANGE];
  TH1F* h_resVsPt_phi_I[nRANGE];
  TH1F* h_resVsPt_phi_F[nRANGE];

  TH1F* h_resVsPt_eta[nRANGE];
  TH1F* h_resVsPt_d0[nRANGE];

  for (int i=0; i<nRANGE; i++) {
    h_resVsPt_pt[i]   = new TH1F("resVsPt_pt_"+ptrange[i],   ";p_{T} residual (L1 - sim) [GeV]; L1 tracks / 0.1", 100, -5.0, 5.0);
    h_resVsPt_pt_C[i] = new TH1F("resVsPt_pt_C_"+ptrange[i], ";p_{T} residual (L1 - sim) [GeV]; L1 tracks / 0.1", 100, -5.0, 5.0);
    h_resVsPt_pt_I[i] = new TH1F("resVsPt_pt_I_"+ptrange[i], ";p_{T} residual (L1 - sim) [GeV]; L1 tracks / 0.1", 100, -5.0, 5.0);
    h_resVsPt_pt_F[i] = new TH1F("resVsPt_pt_F_"+ptrange[i], ";p_{T} residual (L1 - sim) [GeV]; L1 tracks / 0.1", 100, -5.0, 5.0);

    h_resVsPt_ptRel[i]   = new TH1F("resVsPt_ptRel_"+ptrange[i],   ";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02", 300, -0.15, 0.15);
    h_resVsPt_ptRel_C[i] = new TH1F("resVsPt_ptRel_c_"+ptrange[i], ";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02", 300, -0.15, 0.15);
    h_resVsPt_ptRel_I[i] = new TH1F("resVsPt_ptRel_I_"+ptrange[i], ";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02", 300, -0.15, 0.15);
    h_resVsPt_ptRel_F[i] = new TH1F("resVsPt_ptRel_F_"+ptrange[i], ";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02", 300, -0.15, 0.15);

    h_resVsPt_z0[i]   = new TH1F("resVsPt_z0_"+ptrange[i],   ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100, -1, 1);
    h_resVsPt_z0_C[i] = new TH1F("resVsPt_z0_C_"+ptrange[i], ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100, -1, 1);
    h_resVsPt_z0_I[i] = new TH1F("resVsPt_z0_I_"+ptrange[i], ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100, -1, 1);
    h_resVsPt_z0_F[i] = new TH1F("resVsPt_z0_F_"+ptrange[i], ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100, -1, 1);

    h_resVsPt_phi[i]   = new TH1F("resVsPt_phi_"+ptrange[i],   ";#phi residual (L1 - sim) [rad]; L1 tracks / 0.0001", 100, -0.005, 0.005);
    h_resVsPt_phi_C[i] = new TH1F("resVsPt_phi_C_"+ptrange[i], ";#phi residual (L1 - sim) [rad]; L1 tracks / 0.0001", 100, -0.005, 0.005);
    h_resVsPt_phi_I[i] = new TH1F("resVsPt_phi_I_"+ptrange[i], ";#phi residual (L1 - sim) [rad]; L1 tracks / 0.0001", 100, -0.005, 0.005);
    h_resVsPt_phi_F[i] = new TH1F("resVsPt_phi_F_"+ptrange[i], ";#phi residual (L1 - sim) [rad]; L1 tracks / 0.0001", 100, -0.005, 0.005);

    h_resVsPt_eta[i] = new TH1F("resVsPt_eta_"+ptrange[i], ";#eta residual (L1 - sim); L1 tracks / 0.0002", 100, -0.01, 0.01);

    h_resVsPt_d0[i]  = new TH1F("resVsPt_d0_"+ptrange[i],   ";d_{0} residual (L1 - sim) [cm]; L1 tracks / 0.0004",   100, -0.02, 0.02);
  }

  // ----------------------------------------------------------------------------------------------------------------
  // more resolution vs eta

  TH1F* h_resVsEta_eta[nETARANGE];
  TH1F* h_resVsEta_eta_L[nETARANGE];
  TH1F* h_resVsEta_eta_H[nETARANGE];

  TH1F* h_resVsEta_z0[nETARANGE];
  TH1F* h_resVsEta_z0_L[nETARANGE];
  TH1F* h_resVsEta_z0_H[nETARANGE];

  TH1F* h_resVsEta_phi[nETARANGE];
  TH1F* h_resVsEta_phi_L[nETARANGE];
  TH1F* h_resVsEta_phi_H[nETARANGE];

  TH1F* h_resVsEta_ptRel[nETARANGE];
  TH1F* h_resVsEta_ptRel_L[nETARANGE];
  TH1F* h_resVsEta_ptRel_H[nETARANGE];

  TH1F* h_resVsEta_d0[nETARANGE];
  TH1F* h_resVsEta_d0_L[nETARANGE];
  TH1F* h_resVsEta_d0_H[nETARANGE];

  for (int i=0; i<nETARANGE; i++) {
    h_resVsEta_eta[i]   = new TH1F("resVsEta_eta_"+etarange[i],   ";#eta residual (L1 - sim); L1 tracks / 0.0002", 100,-0.01,0.01);
    h_resVsEta_eta_L[i] = new TH1F("resVsEta_eta_L_"+etarange[i], ";#eta residual (L1 - sim); L1 tracks / 0.0002", 100,-0.01,0.01);
    h_resVsEta_eta_H[i] = new TH1F("resVsEta_eta_H_"+etarange[i], ";#eta residual (L1 - sim); L1 tracks / 0.0002", 100,-0.01,0.01);

    h_resVsEta_z0[i]   = new TH1F("resVsEta2_z0_"+etarange[i],   ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.01",   100,-1, 1);
    h_resVsEta_z0_L[i] = new TH1F("resVsEta2_z0_L_"+etarange[i], ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.01",   100,-1, 1);
    h_resVsEta_z0_H[i] = new TH1F("resVsEta2_z0_H_"+etarange[i], ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.01",   100,-1, 1);

    h_resVsEta_phi[i]   = new TH1F("resVsEta2_phi_"+etarange[i],   ";#phi residual (L1 - sim) [rad]; L1 tracks / 0.0001", 100,-0.005,0.005);
    h_resVsEta_phi_L[i] = new TH1F("resVsEta2_phi_L_"+etarange[i], ";#phi residual (L1 - sim) [rad]; L1 tracks / 0.0001", 100,-0.005,0.005);
    h_resVsEta_phi_H[i] = new TH1F("resVsEta2_phi_H_"+etarange[i], ";#phi residual (L1 - sim) [rad]; L1 tracks / 0.0001", 100,-0.005,0.005);

    h_resVsEta_ptRel[i]   = new TH1F("resVsEta2_ptRel_"+etarange[i],  ";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.01", 100,-0.5,0.5);
    h_resVsEta_ptRel_L[i] = new TH1F("resVsEta2_ptRel_L_"+etarange[i],";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02", 100,-0.1,0.1);
    h_resVsEta_ptRel_H[i] = new TH1F("resVsEta2_ptRel_H_"+etarange[i],";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02", 100,-0.25,0.25);

    h_resVsEta_d0[i]    = new TH1F("resVsEta_d0_"+etarange[i],   ";d_{0} residual (L1 - sim) [cm]; L1 tracks / 0.004",    100,-0.02, 0.02);
    h_resVsEta_d0_L[i]  = new TH1F("resVsEta_d0_L_"+etarange[i], ";d_{0} residual (L1 - sim) [cm]; L1 tracks / 0.004",    100,-0.02, 0.02);
    h_resVsEta_d0_H[i]  = new TH1F("resVsEta_d0_H_"+etarange[i], ";d_{0} residual (L1 - sim) [cm]; L1 tracks / 0.004",    100,-0.02, 0.02);

  }



  // ---------------------------------------------------------------------------------------------------------
  //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
  // ---------------------------------------------------------------------------------------------------------

  int nevt = tree->GetEntries();
  cout << "Number of Events = " << nevt << endl;
  if (doLooseMatch) cout << "Using Loose Matching" << endl;
  else              cout << "Using Standard Matching" << endl;


  // ---------------------------------------------------------------------------------------------------------
  // EVENT LOOP
  for (int i=0; i<nevt; i++) {
    tree->GetEntry(i,0);

    // -------------------------------------------------------------------------------------------------------
    // Jet Tracking

    // ----- Jet Preliminaries ----- //

    // Constrained Jet Parameters
    float jet_njet_pt[njets]            = {0};
    float jet_njet_eta[njets]           = {0};
    float jet_njet_phi[njets]           = {0};
    float jet_njet_p[njets]             = {0};

    // Jet Variables
    float a_jet_p[]                     = {0};
    float a_tp_p[]                      = {0};
    float a_trk_p[]                     = {0};
    float a_match_tp_p[]                = {0};
    float a_matchtrk_p[]                = {0};

    // Jet Profile Plot Variables
    float a_tp_per_jet_min[]            = {0};
    float a_trk_per_jet_min[]           = {0};
    float a_match_tp_per_jet_min[]      = {0};
    float a_matchtrk_per_jet_min[]      = {0};
    

    // ----- Jet Loop with njet Constraint ----- //
    //  For Jet Profile and Scatter Plots
    for (int jj=0; jj<(int)jet_pt->size(); jj++) {
      if (jet_pt->at(jj) < 0.2) continue;
      if (jet_pt->at(jj) > TP_maxPt) continue;
      if (fabs(jet_eta->at(jj)) > TP_maxEta) continue;

      // Fill General Jet Histograms
      h_jet_pt  ->Fill(jet_pt->at(jj));
      h_jet_eta ->Fill(jet_eta->at(jj));
      h_jet_phi ->Fill(jet_phi->at(jj));

      // Define Total Momentum and the Histogram
      a_jet_p[jj] = jet_pt->at(jj) * cosh(jet_eta->at(jj));
      h_jet_p->Fill(a_jet_p[jj]);

      // Store Constrained Values
      if(jj<njets) {
        jet_njet_pt[jj]  = jet_pt->at(jj);     //store jet pT for jet profile plot
        jet_njet_eta[jj] = jet_eta->at(jj);    //store jet eta
        jet_njet_phi[jj] = jet_phi->at(jj);    //store jet phi
        jet_njet_p[jj]   = a_jet_p[jj];        //store jet total momentum
      }
    }


    // ----- Sum pT in Jets ----- //
    for (int ij=0; ij<(int)jet_tp_sumpt->size(); ij++) {
      float fraction           = 0;
      float fractionMatch      = 0;
      if (jet_tp_sumpt->at(ij) > 0) {
        fraction      = jet_trk_sumpt->at(ij)/jet_tp_sumpt->at(ij);
        fractionMatch = jet_matchtrk_sumpt->at(ij)/jet_tp_sumpt->at(ij);
      }

      h_jet_tp_sumpt_vspt      ->Fill(jet_tp_sumpt->at(ij),1.0);
      h_jet_trk_sumpt_vspt     ->Fill(jet_tp_sumpt->at(ij),fraction);
      h_jet_matchtrk_sumpt_vspt->Fill(jet_tp_sumpt->at(ij),fractionMatch);

      h_jet_tp_sumpt_vseta      ->Fill(jet_eta->at(ij),1.0);
      h_jet_trk_sumpt_vseta     ->Fill(jet_eta->at(ij),fraction);
      h_jet_matchtrk_sumpt_vseta->Fill(jet_eta->at(ij),fractionMatch);

      h_jet_tp_sumpt_vsphi      ->Fill(jet_phi->at(ij),1.0);
      h_jet_trk_sumpt_vsphi     ->Fill(jet_phi->at(ij),fraction);
      h_jet_matchtrk_sumpt_vsphi->Fill(jet_phi->at(ij),fractionMatch);
    }


    // -------------------------------------------------------------------------------------------------------
    // L1 Tracks

    // Initialize Total Rates Variables
    int ntrkevt_pt2  = 0;
    int ntrkevt_pt3  = 0;
    int ntrkevt_pt10 = 0;

    // ---------- L1 Loop ---------- //
    for (int ik=0; ik<(int)trk_pt->size(); ik++) {
      // For Total Rates
      if (trk_pt->at(ik) > 2.0 && fabs(trk_eta->at(ik))<2.5) {
        ntrk_pt2++;
        ntrkevt_pt2++;
        h_trk_vspt->Fill(trk_pt->at(ik));
      }
      if (trk_pt->at(ik) > 3.0 && fabs(trk_eta->at(ik))<2.5) {
         ntrk_pt3++;
         ntrkevt_pt3++;
      }
      if (trk_pt->at(ik) > 10.0 && fabs(trk_eta->at(ik))<2.5) {
         ntrk_pt10++;
         ntrkevt_pt10++;
      }

      // For L1 Track Efficiency
      if (trk_pt->at(ik) < 0.2) continue;
      if (trk_pt->at(ik) > TP_maxPt) continue;
      if (fabs(trk_eta->at(ik)) > TP_maxEta) continue;

      a_trk_p[ik] = trk_pt->at(ik) * cosh(trk_eta->at(ik));
      h_trk_p->Fill(a_trk_p[ik]);

      h_trk_pt->Fill(trk_pt->at(ik));

      if (trk_pt->at(ik) < 8.0) h_trk_pt_L->Fill(trk_pt->at(ik));
      else h_trk_pt_H->Fill(trk_pt->at(ik));

      if (trk_pt->at(ik) < 8.0 && fabs(trk_eta->at(ik))<1.0) h_trk_pt_LC->Fill(trk_pt->at(ik));
      ///////////

      if (trk_pt->at(ik) > TP_minPt) {
        h_trk_eta  ->Fill(trk_eta->at(ik));
        h_trk_phi  ->Fill(trk_phi->at(ik));
        h_trk_d0   ->Fill(trk_d0->at(ik));
        h_trk_absd0->Fill(fabs(trk_d0->at(ik)));
        h_trk_z0   ->Fill(trk_z0->at(ik));

        if      (trk_pt->at(ik) < 3.0) h_trk_eta_23->Fill(trk_eta->at(ik));
        else if (trk_pt->at(ik) < 5.0) h_trk_eta_35->Fill(trk_eta->at(ik));
        else h_trk_eta_5->Fill(trk_eta->at(ik));

        if (trk_pt->at(ik) < 8.0) {
          h_trk_eta_L->Fill(trk_eta->at(ik));
          h_trk_z0_L ->Fill(trk_z0->at(ik));
        }
        else {
          h_trk_eta_H->Fill(trk_eta->at(ik));
          h_trk_z0_H ->Fill(trk_z0->at(ik));
        }


        // Distance to Jet Axis - L1 Tracks
        float deta_trk    = 0;
        float dphi_trk    = 0;
        float dR_trk      = 0;
        float dpt_trk     = 0;
        float dpt_trk_jet = 0;

        float deta_trk_min = 10000.;
        float dphi_trk_min = 10000.;
        float dR_trk_min   = 10000.;

        int jetID_trk_min  = 0;

        for (int jj=0; jj<(int)jet_pt->size(); jj++) {

          // dR Calculation
          deta_trk = trk_eta->at(ik) - jet_eta->at(jj);
          dphi_trk = trk_phi->at(ik) - jet_phi->at(jj);
          dpt_trk  = trk_pt ->at(ik);
          while (dphi_trk > 3.14159) dphi_trk = fabs(2*3.14159 - dphi_trk);
          dR_trk = sqrt(deta_trk*deta_trk + dphi_trk*dphi_trk);

          if (dR_trk < dR_trk_min) {
            deta_trk_min = deta_trk;
            dphi_trk_min = dphi_trk;
            dR_trk_min   = dR_trk;
            dpt_trk_jet  = jet_pt->at(jj);
            jetID_trk_min = jj;
          }

        }


        // Fill Histograms
        h_dR_trk   ->Fill(dR_trk_min);
        h_deta_trk ->Fill(deta_trk_min);
        h_dphi_trk ->Fill(dphi_trk_min);
        if (dR_trk_min < 0.4) {
          h_dpt_trk          ->Fill(dpt_trk);
          h_dpt_trk_jet      ->Fill(dpt_trk_jet);
          h_dpt_trk_jet_full ->Fill(dpt_trk_jet);
          if (jetID_trk_min<njets) {
            a_trk_per_jet_min[jetID_trk_min]++; //for jet profile plot
          }
        }

        if (trk_pt->at(ik) < 8.0) {
          h_dR_trk_L   ->Fill(dR_trk_min);
          h_deta_trk_L ->Fill(deta_trk_min);
          h_dphi_trk_L ->Fill(dphi_trk_min);
          if (dR_trk_min < 0.4) {
            h_dpt_trk_L     ->Fill(dpt_trk);
            h_dpt_trk_jet_L ->Fill(dpt_trk_jet);
          }
        }

        else {
          h_dR_trk_H   ->Fill(dR_trk_min);
          h_deta_trk_H ->Fill(deta_trk_min);
          h_dphi_trk_H ->Fill(dphi_trk_min);
          if (dR_trk_min < 0.4) {
            h_dpt_trk_H     ->Fill(dpt_trk);
            h_dpt_trk_jet_H ->Fill(dpt_trk_jet);
          }
        }

      }

      // ----- Fill L1 Track nstub Histograms
      if (trk_pt->at(ik) > TP_minPt){
        h_trk_nstub->Fill(trk_nstub->at(ik));

          //Central eta
          if (fabs(trk_eta->at(ik)) < 0.8) {
            h_trk_nstub_C->Fill(trk_nstub->at(ik));
          }

          //Intermediate eta
          else if (fabs(trk_eta->at(ik)) < 1.6 && fabs(trk_eta->at(ik)) >= 0.8) {
            h_trk_nstub_I->Fill(trk_nstub->at(ik));
          }

          //Forward eta
          else if (fabs(trk_eta->at(ik)) >= 1.6) {
            h_trk_nstub_F->Fill(trk_nstub->at(ik));
          }
      }


    }

    // Fill Total Track Rates Histograms
    h_ntrk_pt2  ->Fill(ntrkevt_pt2);
    h_ntrk_pt3  ->Fill(ntrkevt_pt3);
    h_ntrk_pt10 ->Fill(ntrkevt_pt10);



    

    // -------------------------------------------------------------------------------------------------------
    // L1 Matched Track Loop
    for (int im=0; im<(int)matchtrk_pt->size(); im++) {
      
      // -----------------------------------------------------------------------------------------------------
      // Parameters

      // Mins and Maxes
      if (matchtrk_pt->at(im) < 0.2) continue;
      if (matchtrk_pt->at(im) > TP_maxPt) continue;
      if (fabs(matchtrk_eta->at(im)) > TP_maxEta) continue;

      // Use only tracks with min X stubs
      if (matchtrk_nstub->at(im) < L1Tk_minNstub) continue;


      // -----------------------------------------------------------------------------------------------------
      // Proceed with Matched L1s

      // Total Momentum
      a_matchtrk_p[im] = matchtrk_pt->at(im) * cosh(matchtrk_eta->at(im)); //ML1 total momentum
      h_matchtrk_p->Fill(a_matchtrk_p[im]);

      // Fill some histograms
      h_matchtrk_pt->Fill(matchtrk_pt->at(im));
      if (matchtrk_pt->at(im) < 8.0) h_matchtrk_pt_L->Fill(matchtrk_pt->at(im));
      else h_matchtrk_pt_H->Fill(matchtrk_pt->at(im));
      if (matchtrk_pt->at(im) < 8.0 && fabs(matchtrk_eta->at(im))<1.0) h_matchtrk_pt_LC->Fill(matchtrk_pt->at(im));

      if (matchtrk_pt->at(im) > TP_minPt) {
        h_matchtrk_eta  ->Fill(matchtrk_eta->at(im));
        h_matchtrk_phi  ->Fill(matchtrk_phi->at(im));
        h_matchtrk_d0   ->Fill(matchtrk_d0->at(im));
        h_matchtrk_absd0->Fill(fabs(matchtrk_d0->at(im)));
        h_matchtrk_z0   ->Fill(matchtrk_z0->at(im));

        if (matchtrk_pt->at(im) < 8.0) h_matchtrk_eta_L->Fill(matchtrk_eta->at(im));
        else h_matchtrk_eta_H ->Fill(matchtrk_eta->at(im));

        // Distance to Jet Axis - L1 Matched Tracks
        float deta_matchtrk    = 0;
        float dphi_matchtrk    = 0;
        float dR_matchtrk      = 0;
        float dpt_matchtrk     = 0;
        float dpt_matchtrk_jet = 0;

        float deta_matchtrk_min = 10000.;
        float dphi_matchtrk_min = 10000.;
        float dR_matchtrk_min   = 10000.;

        int jetID_matchtrk_min  = 0;

        for (int jj=0; jj<(int)jet_pt->size(); jj++) {

          // dR Calculation
          deta_matchtrk = matchtrk_eta->at(im) - jet_eta->at(jj);
          dphi_matchtrk = matchtrk_phi->at(im) - jet_phi->at(jj);
          dpt_matchtrk  = matchtrk_pt ->at(im);
          while (dphi_matchtrk > 3.14159) dphi_matchtrk = fabs(2*3.14159 - dphi_matchtrk);
          dR_matchtrk = sqrt(deta_matchtrk*deta_matchtrk + dphi_matchtrk*dphi_matchtrk);

          if (dR_matchtrk < dR_matchtrk_min) {
            deta_matchtrk_min = deta_matchtrk;
            dphi_matchtrk_min = dphi_matchtrk;
            dR_matchtrk_min   = dR_matchtrk;
            dpt_matchtrk_jet  = jet_pt->at(jj);
            jetID_matchtrk_min = jj;
          }

        }


        // Fill Histograms
        h_dR_matchtrk   ->Fill(dR_matchtrk_min);
        h_deta_matchtrk ->Fill(deta_matchtrk_min);
        h_dphi_matchtrk ->Fill(dphi_matchtrk_min);
        if (dR_matchtrk_min < 0.4) {
          h_dpt_matchtrk          ->Fill(dpt_matchtrk);
          h_dpt_matchtrk_jet      ->Fill(dpt_matchtrk_jet);
          h_dpt_matchtrk_jet_full ->Fill(dpt_matchtrk_jet);
          if (jetID_matchtrk_min<njets) {
            a_matchtrk_per_jet_min[jetID_matchtrk_min]++; //for jet profile plot
          }
        }

        if (matchtrk_pt->at(im) < 8.0) {
          h_dR_matchtrk_L   ->Fill(dR_matchtrk_min);
          h_deta_matchtrk_L ->Fill(deta_matchtrk_min);
          h_dphi_matchtrk_L ->Fill(dphi_matchtrk_min);
          if (dR_matchtrk_min < 0.4) {
            h_dpt_matchtrk_L     ->Fill(dpt_matchtrk);
            h_dpt_matchtrk_jet_L ->Fill(dpt_matchtrk_jet);
          }
        }

        else {
          h_dR_matchtrk_H   ->Fill(dR_matchtrk_min);
          h_deta_matchtrk_H ->Fill(deta_matchtrk_min);
          h_dphi_matchtrk_H ->Fill(dphi_matchtrk_min);
          if (dR_matchtrk_min < 0.4) {
            h_dpt_matchtrk_H     ->Fill(dpt_matchtrk);
            h_dpt_matchtrk_jet_H ->Fill(dpt_matchtrk_jet);
          }
        }

      }

    }

    // ----- Fill L1 Jet Profile Plot Histograms
    // plot # of tracks per jet
    for(int jjet = 0; jjet<njets; jjet++) {

      if(a_trk_per_jet_min[jjet]>0) {
        h_jet_profile_L1 ->Fill(jet_njet_pt[jjet],a_trk_per_jet_min[jjet]);
      }
      if(a_matchtrk_per_jet_min[jjet]>0) {
        h_jet_profile_ML1->Fill(jet_njet_pt[jjet],a_matchtrk_per_jet_min[jjet]);
      }
    }



    // -------------------------------------------------------------------------------------------------------
    // Tracking Particles

    // ---------- Tracking Particle Loop ---------- //
    for (int it=0; it<(int)tp_pt->size(); it++) {

      // For Total Track Rates
      if (fabs(tp_dxy->at(it)) < 1 && fabs(tp_eta->at(it)) < TP_maxEta) { //the stub requirements are applied when making the ntuples
        if (tp_pt->at(it) > 2.0) {
          ntp_pt2++;
          h_tp_vspt->Fill(tp_pt->at(it));
        }
        if (tp_pt->at(it) > 3.0) ntp_pt3++;
        if (tp_pt->at(it) > 10.0) ntp_pt10++;
      }

      // Cut on PDG ID at plot stage?
      if (TP_select_pdgid != 0) {
        if (abs(tp_pdgid->at(it)) != abs(TP_select_pdgid)) continue;
      }

      // cut on event ID (eventid=0 means the TP is from the primary interaction, so *not* selecting only eventid=0 means including stuff from pileup)
      if (TP_select_eventid == 0 && tp_eventid->at(it) != 0) continue;


      //-----Tracking Particles-----//
      if (tp_dxy->at(it) > 1) continue;
      if (tp_pt->at(it) < 0.2) continue;
      if (tp_pt->at(it) > TP_maxPt) continue;
      if (fabs(tp_eta->at(it)) > TP_maxEta) continue;

      a_tp_p[it] = tp_pt->at(it) * cosh(tp_eta->at(it)); //TP total momentum
      h_tp_p->Fill(a_tp_p[it]);

      h_tp_pt->Fill(tp_pt->at(it));
      if (tp_pt->at(it) < 8.0) h_tp_pt_L->Fill(tp_pt->at(it));
      else h_tp_pt_H->Fill(tp_pt->at(it));
      if (tp_pt->at(it) < 8.0 && fabs(tp_eta->at(it))<1.0) h_tp_pt_LC->Fill(tp_pt->at(it));

      if (tp_pt->at(it) > TP_minPt) {
        if (fabs(tp_eta->at(it)) < 1.0) n_all_eta1p0++;
        else if (fabs(tp_eta->at(it)) < 1.75) n_all_eta1p75++;
        else n_all_eta2p5++;

        h_tp_eta  ->Fill(tp_eta->at(it));
        h_tp_phi  ->Fill(tp_phi->at(it));
        h_tp_z0   ->Fill(tp_z0->at(it));
        h_tp_d0   ->Fill(tp_d0->at(it));
        h_tp_absd0->Fill(fabs(tp_d0->at(it)));

        if      (tp_pt->at(it) < 3.0) h_tp_eta_23->Fill(tp_eta->at(it));
        else if (tp_pt->at(it) < 5.0) h_tp_eta_35->Fill(tp_eta->at(it));
        else h_tp_eta_5->Fill(tp_eta->at(it));

        if (tp_pt->at(it) < 8.0) {
          h_tp_eta_L->Fill(tp_eta->at(it));
          h_tp_z0_L ->Fill(tp_z0->at(it));
        }
        else {
          h_tp_eta_H->Fill(tp_eta->at(it));
          h_tp_z0_H ->Fill(tp_z0->at(it));
        }

        // Distance to Jet Axis - TP
        float deta_tp = 0;
        float dphi_tp = 0;
        float dR_tp   = 0;
        float dpt_tp  = 0;
        float dpt_jet = 0;

        float deta_tp_min = 10000.;
        float dphi_tp_min = 10000.;
        float dR_tp_min   = 10000.;

        int jetID_tp_min  = 0;

        for (int jj=0; jj<(int)jet_pt->size(); jj++) {

          // dR Calculation
          deta_tp = tp_eta->at(it) - jet_eta->at(jj);
          dphi_tp = tp_phi->at(it) - jet_phi->at(jj);
          dpt_tp  = tp_pt ->at(it);
          while (dphi_tp > 3.14159) dphi_tp = fabs(2*3.14159 - dphi_tp);
          dR_tp = sqrt(deta_tp*deta_tp + dphi_tp*dphi_tp);

          if (dR_tp < dR_tp_min) {
            deta_tp_min = deta_tp;
            dphi_tp_min = dphi_tp;
            dR_tp_min   = dR_tp;
            dpt_jet     = jet_pt->at(jj);
            jetID_tp_min = jj;
          }
        }


        // Fill Histograms
        h_dR_tp  ->Fill(dR_tp_min);
        h_deta_tp->Fill(deta_tp_min);
        h_dphi_tp->Fill(dphi_tp_min);

        if (dR_tp_min < 0.4) {
          h_dpt_tp       ->Fill(dpt_tp);
          h_dpt_jet      ->Fill(dpt_jet);
          h_dpt_jet_full ->Fill(dpt_jet);
          if (jetID_tp_min<njets) {
            a_tp_per_jet_min[jetID_tp_min]++; //for jet profile plot
          }

            //central eta
            if (fabs(tp_eta->at(it)) < 0.8) {
              h_dpt_jet_C      ->Fill(dpt_jet);
              h_dpt_jet_full_C ->Fill(dpt_jet);
            }
            //intermediate eta
            else if (fabs(tp_eta->at(it)) < 1.6 && fabs(tp_eta->at(it)) >= 0.8) {
              h_dpt_jet_I      ->Fill(dpt_jet);
              h_dpt_jet_full_I ->Fill(dpt_jet);
            }
            //forward eta
            else if (fabs(tp_eta->at(it)) >= 1.6) {
              h_dpt_jet_F      ->Fill(dpt_jet);
              h_dpt_jet_full_F ->Fill(dpt_jet);
            }

        }

        if (tp_pt->at(it) < 8.0) {
          h_dR_tp_L  ->Fill(dR_tp_min);
          h_deta_tp_L->Fill(deta_tp_min);
          h_dphi_tp_L->Fill(dphi_tp_min);
          if (dR_tp_min < 0.4) {
            h_dpt_tp_L  ->Fill(dpt_tp);
            h_dpt_jet_L ->Fill(dpt_jet);
              //central eta
              if (fabs(tp_eta->at(it)) < 0.8) {
                h_dpt_jet_LC     ->Fill(dpt_jet);
              }
              //intermediate eta
              else if (fabs(tp_eta->at(it)) < 1.6 && fabs(tp_eta->at(it)) >= 0.8) {
                h_dpt_jet_LI     ->Fill(dpt_jet);
              }
              //forward eta
              else if (fabs(tp_eta->at(it)) >= 1.6) {
                h_dpt_jet_LF     ->Fill(dpt_jet);
              }
          }
        }

        else {
          h_dR_tp_H  ->Fill(dR_tp_min);
          h_deta_tp_H->Fill(deta_tp_min);
          h_dphi_tp_H->Fill(dphi_tp_min);
          if (dR_tp_min < 0.4) {
            h_dpt_tp_H  ->Fill(dpt_tp);
            h_dpt_jet_H ->Fill(dpt_jet);
              //central eta
              if (fabs(tp_eta->at(it)) < 0.8) {
                h_dpt_jet_HC     ->Fill(dpt_jet);
              }
              //intermediate eta
              else if (fabs(tp_eta->at(it)) < 1.6 && fabs(tp_eta->at(it)) >= 0.8) {
                h_dpt_jet_HI     ->Fill(dpt_jet);
              }
              //forward eta
              else if (fabs(tp_eta->at(it)) >= 1.6) {
                h_dpt_jet_HF     ->Fill(dpt_jet);
              }
          }

        }

      }

      // -----------------------------------------------------------------------------------------------------
      // Was the tracking particle matched to a L1 track? (also covers unique matches)
      if (tp_nmatch->at(it) < 1) continue;

      // -----------------------------------------------------------------------------------------------------
      // Use only tracks with min X stubs
      if (matchtrk_nstub->at(it) < L1Tk_minNstub) continue;

      // -----------------------------------------------------------------------------------------------------
      // fill chi2 & chi2/dof histograms before making chi2 cut
      float chi2    = matchtrk_chi2->at(it);
      int ndof      = 2*matchtrk_nstub->at(it)-4;
      float chi2dof = (float)chi2/ndof;
      if (chi2    > 100) chi2    = 99.9; //for overflow bin
      if (chi2dof > 20)  chi2dof = 19.99; //for overflow bin

      if (tp_pt->at(it) > TP_minPt) {

        h_match_trk_chi2    ->Fill(chi2);
        h_match_trk_chi2_dof->Fill(chi2dof);

        // central eta
        if (fabs(tp_eta->at(it)) < 0.8) {
          if (tp_pt->at(it) < 8.0) {
            h_match_trk_chi2_C_L    ->Fill(chi2);
            h_match_trk_chi2_dof_C_L->Fill(chi2dof);
          }
          else {
            h_match_trk_chi2_C_H    ->Fill(chi2);
            h_match_trk_chi2_dof_C_H->Fill(chi2dof);
          }
        }

        // intermediate eta
        else if (fabs(tp_eta->at(it)) < 1.6 && fabs(tp_eta->at(it)) >= 0.8) {
          if (tp_pt->at(it) < 8.0) {
            h_match_trk_chi2_I_L    ->Fill(chi2);
            h_match_trk_chi2_dof_I_L->Fill(chi2dof);
          }
          else {
            h_match_trk_chi2_I_H    ->Fill(chi2);
            h_match_trk_chi2_dof_I_H->Fill(chi2dof);
          }
        }

        // forward eta
        else if (fabs(tp_eta->at(it)) >= 1.6) {
          if (tp_pt->at(it) < 8.0) {
            h_match_trk_chi2_F_L    ->Fill(chi2);
            h_match_trk_chi2_dof_F_L->Fill(chi2dof);
          }
          else {
            h_match_trk_chi2_F_H    ->Fill(chi2);
            h_match_trk_chi2_dof_F_H->Fill(chi2dof);
          }
        }

      }

      // -----------------------------------------------------------------------------------------------------
      // Cut on chi2?
      if (matchtrk_chi2->at(it)      > L1Tk_maxChi2)    continue;
      if (matchtrk_chi2->at(it)/ndof > L1Tk_maxChi2dof) continue;



      // -----------------------------------------------------------------------------------------------------
      // Matched TP

      // Total momentum for matched TPs
      vector<float>* match_tp_p;
      TBranch* b_match_tp_p;
      match_tp_p  = 0;

      a_match_tp_p[it] = tp_pt->at(it) * cosh(tp_eta->at(it)); //MTP total momentum
      h_match_tp_p->Fill(a_match_tp_p[it]);

      //-----Matched Tracking Particles-----//
      h_match_tp_pt->Fill(tp_pt->at(it));
      if (tp_pt->at(it) < 8.0) h_match_tp_pt_L->Fill(tp_pt->at(it));
      else h_match_tp_pt_H->Fill(tp_pt->at(it));
      if (tp_pt->at(it) < 8.0 && fabs(tp_eta->at(it))<1.0) h_match_tp_pt_LC->Fill(tp_pt->at(it));

      if (tp_pt->at(it) > TP_minPt) {
        h_match_tp_eta  ->Fill(tp_eta->at(it));
        h_match_tp_phi  ->Fill(tp_phi->at(it));
        h_match_tp_z0   ->Fill(tp_z0->at(it));
        h_match_tp_d0   ->Fill(tp_d0->at(it));
        h_match_tp_absd0->Fill(fabs(tp_d0->at(it)));

        if (fabs(tp_eta->at(it)) < 1.0) n_match_eta1p0++;
        else if (fabs(tp_eta->at(it)) < 1.75) n_match_eta1p75++;
        else n_match_eta2p5++;

        if      (tp_pt->at(it) < 3.0) h_match_tp_eta_23->Fill(tp_eta->at(it));
        else if (tp_pt->at(it) < 5.0) h_match_tp_eta_35->Fill(tp_eta->at(it));
        else h_match_tp_eta_5->Fill(tp_eta->at(it));

        if (tp_pt->at(it) < 8.0) {
          h_match_tp_eta_L->Fill(tp_eta->at(it));
          h_match_tp_z0_L ->Fill(tp_z0->at(it));
        }
        else {
          h_match_tp_eta_H->Fill(tp_eta->at(it));
          h_match_tp_z0_H ->Fill(tp_z0->at(it));
        }

        // Distance to Jet Axis - Matched TP
        float deta_match_tp = 0;
        float dphi_match_tp = 0;
        float dR_match_tp   = 0;
        float dpt_match_tp  = 0;
        float dpt_match_jet = 0;

        float deta_match_tp_min = 10000.;
        float dphi_match_tp_min = 10000.;
        float dR_match_tp_min   = 10000.;

        int jetID_match_tp_min  = 0;

        for (int jj=0; jj<(int)jet_pt->size(); jj++) {

          // dR Calculation
          deta_match_tp = tp_eta->at(it) - jet_eta->at(jj);
          dphi_match_tp = tp_phi->at(it) - jet_phi->at(jj);
          dpt_match_tp  = tp_pt->at(it);
          while (dphi_match_tp > 3.14159) dphi_match_tp = fabs(2*3.14159 - dphi_match_tp);
          dR_match_tp = sqrt(deta_match_tp*deta_match_tp + dphi_match_tp*dphi_match_tp);

          if (dR_match_tp < dR_match_tp_min) {
            deta_match_tp_min = deta_match_tp;
            dphi_match_tp_min = dphi_match_tp;
            dR_match_tp_min   = dR_match_tp;
            dpt_match_jet     = jet_pt->at(jj);
            jetID_match_tp_min = jj;
          }

        }


        // Fill Histograms
        h_dR_match_tp  ->Fill(dR_match_tp_min);
        h_deta_match_tp->Fill(deta_match_tp_min);
        h_dphi_match_tp->Fill(dphi_match_tp_min);

        if (dR_match_tp_min < 0.4) {
          h_dpt_match_tp       ->Fill(dpt_match_tp);
          h_dpt_match_jet      ->Fill(dpt_match_jet);
          h_dpt_match_jet_full ->Fill(dpt_match_jet);
          if (jetID_match_tp_min<njets) {
            a_match_tp_per_jet_min[jetID_match_tp_min]++; //for jet profile plot
          }

            //central eta
            if (fabs(tp_eta->at(it)) < 0.8) {
              h_dpt_match_jet_C      ->Fill(dpt_match_jet);
              h_dpt_match_jet_full_C ->Fill(dpt_match_jet);
            }
            //intermediate eta
            else if (fabs(tp_eta->at(it)) < 1.6 && fabs(tp_eta->at(it)) >= 0.8) {
              h_dpt_match_jet_I      ->Fill(dpt_match_jet);
              h_dpt_match_jet_full_I ->Fill(dpt_match_jet);
            }
            //forward eta
            else if (fabs(tp_eta->at(it)) >= 1.6) {
              h_dpt_match_jet_F      ->Fill(dpt_match_jet);
              h_dpt_match_jet_full_F ->Fill(dpt_match_jet);
            }

        }

        if (tp_pt->at(it) < 8.0) {
          h_dR_match_tp_L  ->Fill(dR_match_tp_min);
          h_deta_match_tp_L->Fill(deta_match_tp_min);
          h_dphi_match_tp_L->Fill(dphi_match_tp_min);
          if (dR_match_tp_min < 0.4) {
            h_dpt_match_tp_L  ->Fill(dpt_match_tp);
            h_dpt_match_jet_L ->Fill(dpt_match_jet);
            //central eta
              if (fabs(tp_eta->at(it)) < 0.8) {
                h_dpt_match_jet_LC ->Fill(dpt_match_jet);
              }
            //intermediate eta
              else if (fabs(tp_eta->at(it)) < 1.6 && fabs(tp_eta->at(it)) >= 0.8) {
                h_dpt_match_jet_LI ->Fill(dpt_match_jet);
              }
            //forward eta
              else if (fabs(tp_eta->at(it)) >= 1.6) {
                h_dpt_match_jet_LF ->Fill(dpt_match_jet);
              }
          }

        }

        else {
          h_dR_match_tp_H  ->Fill(dR_match_tp_min);
          h_deta_match_tp_H->Fill(deta_match_tp_min);
          h_dphi_match_tp_H->Fill(dphi_match_tp_min);
          if (dR_match_tp_min < 0.4) {
            h_dpt_match_tp_H  ->Fill(dpt_match_tp);
            h_dpt_match_jet_H ->Fill(dpt_match_jet);
            //central eta
              if (fabs(tp_eta->at(it)) < 0.8) {
                h_dpt_match_jet_HC ->Fill(dpt_match_jet);
              }
            //intermediate eta
              else if (fabs(tp_eta->at(it)) < 1.6 && fabs(tp_eta->at(it)) >= 0.8) {
                h_dpt_match_jet_HI ->Fill(dpt_match_jet);
              }
            //forward eta
              else if (fabs(tp_eta->at(it)) >= 1.6) {
                h_dpt_match_jet_HF ->Fill(dpt_match_jet);
              }
          }
        }

      }


      // -----------------------------------------------------------------------------------------------------
      // For the following, only consider TPs with pT > TP_minPt
      if (tp_pt->at(it) < TP_minPt) continue;


      // ----- Fill TP nstub Histograms
      h_tp_nstub->Fill(tp_nstub->at(it));
        //Central eta
        if (fabs(tp_eta->at(it)) < 0.8) {
          h_tp_nstub_C->Fill(tp_nstub->at(it));
        }
        //Intermediate eta
        else if (fabs(tp_eta->at(it)) < 1.6 && fabs(tp_eta->at(it)) >= 0.8) {
          h_tp_nstub_I->Fill(tp_nstub->at(it));
        }
        //Forward eta
        else if (fabs(tp_eta->at(it)) >= 1.6) {
          h_tp_nstub_F->Fill(tp_nstub->at(it));
        }


      // ----- Fill Matched TP nstub Histograms
      h_match_trk_nstub->Fill(matchtrk_nstub->at(it));
        //Central eta
        if (fabs(tp_eta->at(it)) < 0.8) {
          h_match_trk_nstub_C->Fill(matchtrk_nstub->at(it));
        }
        //Intermediate eta
        else if (fabs(tp_eta->at(it)) < 1.6 && fabs(tp_eta->at(it)) >= 0.8) {
          h_match_trk_nstub_I->Fill(matchtrk_nstub->at(it));
        }
        //Forward eta
        else if (fabs(tp_eta->at(it)) >= 1.6) {
          h_match_trk_nstub_F->Fill(matchtrk_nstub->at(it));
        }



	    // ----------------------------------------------------------------------------------------------------------------
      // Fill Resolution Histograms
      
      // Define values
      h_res_pt   ->Fill(matchtrk_pt->at(it)  - tp_pt->at(it));
      h_res_ptRel->Fill((matchtrk_pt->at(it) - tp_pt->at(it))/tp_pt->at(it));
      h_res_eta  ->Fill(matchtrk_eta->at(it) - tp_eta->at(it));
      h_res_phi  ->Fill(matchtrk_phi->at(it) - tp_phi->at(it));
      h_res_z0   ->Fill(matchtrk_z0->at(it)  - tp_z0->at(it));
      if (matchtrk_d0->at(it) < 999.) h_res_d0->Fill(matchtrk_d0->at(it) - tp_d0->at(it));

      // ----- Fill z0 resolution histograms
      if      (fabs(tp_eta->at(it)) < 0.8)  h_res_z0_C->Fill(matchtrk_z0->at(it) - tp_z0->at(it));
      else if (fabs(tp_eta->at(it)) < 1.6 && 
               fabs(tp_eta->at(it)) >= 0.8) h_res_z0_I->Fill(matchtrk_z0->at(it) - tp_z0->at(it));
      else if (fabs(tp_eta->at(it)) >= 1.6) h_res_z0_F->Fill(matchtrk_z0->at(it) - tp_z0->at(it));

      if (tp_pt->at(it) < 8.0) {
    		h_res_z0_L->Fill(matchtrk_z0->at(it) - tp_z0->at(it));
    		if (fabs(tp_eta->at(it)) < 1.0) h_res_z0_C_L->Fill(matchtrk_z0->at(it) - tp_z0->at(it));
    		else                            h_res_z0_F_L->Fill(matchtrk_z0->at(it) - tp_z0->at(it));
      }
      else {
    		h_res_z0_H->Fill(matchtrk_z0->at(it) - tp_z0->at(it));
    		if (fabs(tp_eta->at(it)) < 1.0) h_res_z0_C_H->Fill(matchtrk_z0->at(it) - tp_z0->at(it));
    		else                            h_res_z0_F_H->Fill(matchtrk_z0->at(it) - tp_z0->at(it));
      }

      // ----- Fill d0 resolution histograms
      if (matchtrk_d0->at(it) < 999.) {
    		if (fabs(tp_eta->at(it)) < 0.8) h_res_d0_C->Fill(matchtrk_d0->at(it) - tp_d0->at(it));
    		else if (fabs(tp_eta->at(it)) < 1.6 && 
                 fabs(tp_eta->at(it)) >= 0.8) h_res_d0_I->Fill(matchtrk_d0->at(it) - tp_d0->at(it));
    		else if (fabs(tp_eta->at(it)) >= 1.6) h_res_d0_F->Fill(matchtrk_d0->at(it) - tp_d0->at(it));

    		if (tp_pt->at(it) < 8.0) {
  	  		h_res_d0_L->Fill(matchtrk_d0->at(it) - tp_d0->at(it));
  	  		if (fabs(tp_eta->at(it)) < 1.0) h_res_d0_C_L->Fill(matchtrk_d0->at(it) - tp_d0->at(it));
  	  		else                            h_res_d0_F_L->Fill(matchtrk_d0->at(it) - tp_d0->at(it));
    		}
    		else {
  	  		h_res_d0_H->Fill(matchtrk_d0->at(it) - tp_d0->at(it));
  	  		if (fabs(tp_eta->at(it)) < 1.0) h_res_d0_C_H->Fill(matchtrk_d0->at(it) - tp_d0->at(it));
  	  		else                            h_res_d0_F_H->Fill(matchtrk_d0->at(it) - tp_d0->at(it));
		    }

      }

      // ----- Fill resolution vs. pT Histograms
      for (int im=0; im<nRANGE; im++) {
		    if ( (tp_pt->at(it) > (float)im*5.0) && (tp_pt->at(it) < (float)im*5.0+5.0) ) {
  	  		h_resVsPt_pt[im]   ->Fill(matchtrk_pt->at(it)  - tp_pt->at(it));
  	  		h_resVsPt_ptRel[im]->Fill((matchtrk_pt->at(it) - tp_pt->at(it))/tp_pt->at(it));
  	  		h_resVsPt_eta[im]  ->Fill(matchtrk_eta->at(it) - tp_eta->at(it));
  	  		h_resVsPt_phi[im]  ->Fill(matchtrk_phi->at(it) - tp_phi->at(it));
  	  		h_resVsPt_z0[im]   ->Fill(matchtrk_z0->at(it)  - tp_z0->at(it));

        	h_absResVsPt_pt[im]   ->Fill( fabs( matchtrk_pt->at(it)  - tp_pt->at(it) ));
        	h_absResVsPt_ptRel[im]->Fill( fabs( (matchtrk_pt->at(it) - tp_pt->at(it)) )/tp_pt->at(it) );
        	h_absResVsPt_z0[im]   ->Fill( fabs( matchtrk_z0->at(it)  - tp_z0->at(it) ) );
        	h_absResVsPt_phi[im]  ->Fill( fabs( matchtrk_phi->at(it) - tp_phi->at(it) ) );
        	h_absResVsPt_eta[im]  ->Fill( fabs( matchtrk_eta->at(it) - tp_eta->at(it) ) );

  	  		if (fabs(tp_eta->at(it)) < 0.8) {
  	    		h_resVsPt_pt_C[im]   ->Fill(matchtrk_pt->at(it)  - tp_pt->at(it));
  	    		h_resVsPt_ptRel_C[im]->Fill((matchtrk_pt->at(it) - tp_pt->at(it))/tp_pt->at(it));
  	    		h_resVsPt_z0_C[im]   ->Fill(matchtrk_z0->at(it)  - tp_z0->at(it));
  	    		h_resVsPt_phi_C[im]  ->Fill(matchtrk_phi->at(it) - tp_phi->at(it));
  	  		}
  	  		else if (fabs(tp_eta->at(it)) < 1.6 && fabs(tp_eta->at(it)) >= 0.8) {
  	    		h_resVsPt_pt_I[im]   ->Fill(matchtrk_pt->at(it)  - tp_pt->at(it));
  	    		h_resVsPt_ptRel_I[im]->Fill((matchtrk_pt->at(it) - tp_pt->at(it))/tp_pt->at(it));
  	    		h_resVsPt_z0_I[im]   ->Fill(matchtrk_z0->at(it)  - tp_z0->at(it));
  	    		h_resVsPt_phi_I[im]  ->Fill(matchtrk_phi->at(it) - tp_phi->at(it));
  	  		}
  	  		else if (fabs(tp_eta->at(it)) >= 1.6) {
  	    		h_resVsPt_pt_F[im]   ->Fill(matchtrk_pt->at(it)  - tp_pt->at(it));
  	    		h_resVsPt_ptRel_F[im]->Fill((matchtrk_pt->at(it) - tp_pt->at(it))/tp_pt->at(it));
  	    		h_resVsPt_z0_F[im]   ->Fill(matchtrk_z0->at(it)  - tp_z0->at(it));
  	    		h_resVsPt_phi_F[im]  ->Fill(matchtrk_phi->at(it) - tp_phi->at(it));
  	  		}
  	  		if (matchtrk_d0->at(it) < 999) {
  	    		h_resVsPt_d0[im]->Fill(matchtrk_d0->at(it) - tp_d0->at(it));
  	    		h_absResVsPt_d0[im]->Fill( fabs( matchtrk_d0->at(it) - tp_d0->at(it) ) );
  	  		}
		    }
      }

      for (int im=4; im<nRANGE_L+4; im++) {
        if ( (tp_pt->at(it) > (float)im*0.5 ) && (tp_pt->at(it) <= (float)im*0.5+0.5) ) {
          h_absResVsPt_pt_L[im-4]   ->Fill( fabs( matchtrk_pt->at(it)  - tp_pt->at(it) ));
          h_absResVsPt_ptRel_L[im-4]->Fill( fabs( (matchtrk_pt->at(it) - tp_pt->at(it)) )/tp_pt->at(it) );
          h_absResVsPt_z0_L[im-4]   ->Fill( fabs( matchtrk_z0->at(it)  - tp_z0->at(it) ) );
          h_absResVsPt_phi_L[im-4]  ->Fill( fabs( matchtrk_phi->at(it) - tp_phi->at(it) ) );
          h_absResVsPt_eta_L[im-4]  ->Fill( fabs( matchtrk_eta->at(it) - tp_eta->at(it) ) );
        }
      }

      // ----- Fill Resolution vs. eta Histograms
      for (int im=0; im<nETARANGE; im++) {
       	//if ( (fabs(tp_eta->at(it)) > (float)im*0.1) && (fabs(tp_eta->at(it)) < (float)im*0.1+0.1) ) {
        if ( (fabs(tp_eta->at(it)) > (float)im*0.2) && (fabs(tp_eta->at(it)) < (float)im*0.2+0.2) ) {
    	 		h_resVsEta_ptRel[im]->Fill((matchtrk_pt->at(it) - tp_pt->at(it))/tp_pt->at(it));
    	 		h_resVsEta_eta[im]  ->Fill(matchtrk_eta->at(it) - tp_eta->at(it));
    	 		h_resVsEta_phi[im]  ->Fill(matchtrk_phi->at(it) - tp_phi->at(it));
    	 		h_resVsEta_z0[im]   ->Fill(matchtrk_z0->at(it)  - tp_z0->at(it));

         	h_absResVsEta_ptRel[im]->Fill( fabs( (matchtrk_pt->at(it) - tp_pt->at(it)) )/tp_pt->at(it));
         	h_absResVsEta_eta[im]  ->Fill( fabs( matchtrk_eta->at(it) - tp_eta->at(it) ) );
         	h_absResVsEta_phi[im]  ->Fill( fabs( matchtrk_phi->at(it) - tp_phi->at(it) ) );
         	h_absResVsEta_z0[im]   ->Fill( fabs( matchtrk_z0->at(it)  - tp_z0->at(it) ) );

    		 	if (tp_pt->at(it)<8.0) {
  		   		h_resVsEta_ptRel_L[im]->Fill((matchtrk_pt->at(it) - tp_pt->at(it))/tp_pt->at(it));
  		   		h_resVsEta_eta_L[im]->Fill(matchtrk_eta->at(it) - tp_eta->at(it));
  		   		h_resVsEta_z0_L[im]->Fill(matchtrk_z0->at(it) - tp_z0->at(it));
  		   		h_resVsEta_phi_L[im]->Fill(matchtrk_phi->at(it) - tp_phi->at(it));

  		   		h_absResVsEta_ptRel_L[im]->Fill( fabs( (matchtrk_pt->at(it) - tp_pt->at(it)) )/tp_pt->at(it));
  		   		h_absResVsEta_eta_L[im]  ->Fill( fabs( matchtrk_eta->at(it) - tp_eta->at(it) ) );
  		   		h_absResVsEta_phi_L[im]  ->Fill( fabs( matchtrk_phi->at(it) - tp_phi->at(it) ) );
  		   		h_absResVsEta_z0_L[im]   ->Fill( fabs( matchtrk_z0->at(it)  - tp_z0->at(it) ) );

  		   		if (matchtrk_d0->at(it) < 999) {
  		     		h_resVsEta_d0_L[im]->Fill(matchtrk_d0->at(it) - tp_d0->at(it));
  		     		h_absResVsEta_d0_L[im]->Fill( fabs( matchtrk_d0->at(it) - tp_d0->at(it) ) );
  		   		}
	 		    }
		 	    else {
  		   		h_resVsEta_ptRel_H[im]->Fill((matchtrk_pt->at(it) - tp_pt->at(it))/tp_pt->at(it));
  		   		h_resVsEta_eta_H[im]->Fill(matchtrk_eta->at(it) - tp_eta->at(it));
  		   		h_resVsEta_z0_H[im]->Fill(matchtrk_z0->at(it) - tp_z0->at(it));
  		   		h_resVsEta_phi_H[im]->Fill(matchtrk_phi->at(it) - tp_phi->at(it));

  		   		h_absResVsEta_ptRel_H[im]->Fill( fabs( (matchtrk_pt->at(it) - tp_pt->at(it)) )/tp_pt->at(it));
  		   		h_absResVsEta_eta_H[im]  ->Fill( fabs( matchtrk_eta->at(it) - tp_eta->at(it) ) );
  		   		h_absResVsEta_phi_H[im]  ->Fill( fabs( matchtrk_phi->at(it) - tp_phi->at(it) ) );
  		   		h_absResVsEta_z0_H[im]   ->Fill( fabs( matchtrk_z0->at(it)  - tp_z0->at(it) ) );

  			   	if (matchtrk_d0->at(it) < 999) {
  			     	h_resVsEta_d0_H[im]->Fill(matchtrk_d0->at(it) - tp_d0->at(it));
  			     	h_absResVsEta_d0_H[im]->Fill( fabs( matchtrk_d0->at(it) - tp_d0->at(it) ) );
  			   	}
	 		    }
    			if (matchtrk_d0->at(it) < 999) {
    			  h_resVsEta_d0[im]->Fill(matchtrk_d0->at(it) - tp_d0->at(it));
    			  h_absResVsEta_d0[im]->Fill( fabs( matchtrk_d0->at(it) - tp_d0->at(it) ) );
    			}
        }
      }

      // ----- Fill Resolution vs. phi Histograms
      for (int im=0; im<nPHIRANGE; im++) {
        if ( (tp_phi->at(it) > (float)im*0.2-3.2) && (tp_phi->at(it) < (float)im*0.2-3.0) ) {
    	  		h_absResVsPhi_pt[im]   ->Fill( fabs(matchtrk_pt->at(it)  - tp_pt->at(it)) );
    	  		h_absResVsPhi_ptRel[im]->Fill( fabs((matchtrk_pt->at(it) - tp_pt->at(it)) )/tp_pt->at(it));
    		}
      }

    
    } // end of tracking particle loop
    // --------------------------------


    // --------------------------------
    // ----- Fill TP Jet Profile Plot Histograms
    
    // plot # of tracks per jet
    for(int ijet = 0; ijet<njets; ijet++) {

      if(a_tp_per_jet_min[ijet]>0) {
        h_jet_profile_TP ->Fill(jet_njet_pt[ijet],a_tp_per_jet_min[ijet]);
      }
      if(a_match_tp_per_jet_min[ijet]>0) {
        h_jet_profile_MTP->Fill(jet_njet_pt[ijet],a_match_tp_per_jet_min[ijet]);
      }
    }
    // --------------------------------
    


  } // end of event loop
  // ---------------------------------------------------------------------------------------------------------



  // ---------------------------------------------------------------------------------------------------------
  //        * * * * *     E N D    E V E N T S     * * * * *
  // ---------------------------------------------------------------------------------------------------------



  // ---------------------------------------------------------------------------------------------------------
  // Output File for Histograms
  // ---------------------------------------------------------------------------------------------------------

  if (doLooseMatch) type = type+"_loose";
  
  if (TP_select_pdgid != 0) {
    char pdgidtxt[500];
    sprintf(pdgidtxt,"_pdgid%i",TP_select_pdgid);
    type = type+pdgidtxt;
  }

  TString subDIR = "";
  
  TFile* fout = new TFile(subDIR+"plots_"+type+".root","recreate");
   
  


  // ---------------------------------------------------------------------------------------------------------
  // Draw and Save Plots
  // ---------------------------------------------------------------------------------------------------------

  char ctxt[500];
  char ctxtC[500];
  char ctxtI[500];
  char ctxtF[500];
  TCanvas c;

  TString DIR = "TrkPlots/";


  // ---------------------------------------------------------------------------------------------------------
  // Title Plots
  // ---------------------------------------------------------------------------------------------------------

  // To quickly turn off titles, uncomment the gStyle->SetOptTitle(0) line at the end.
  TString htitle;    // This is the title of all histograms
  htitle = type;     // By default, the histogram titles will be the input file name

 
  
  // ---------------------------------------------------------------------------------------------------------
  // Resolution Plots
  // ---------------------------------------------------------------------------------------------------------

  // ----------------------------------------------------------------------------------------------------------------
  // 2D plots  

  // res vs pT
  TH1F* h2_resVsPt_pt   = new TH1F("resVsPt2_pt",   ";Tracking particle p_{T} [GeV]; p_{T} resolution [GeV]", 20,0,100);
  TH1F* h2_resVsPt_pt_C = new TH1F("resVsPt2_pt_C", ";Tracking particle p_{T} [GeV]; p_{T} resolution [GeV]", 20,0,100);
  TH1F* h2_resVsPt_pt_I = new TH1F("resVsPt2_pt_I", ";Tracking particle p_{T} [GeV]; p_{T} resolution [GeV]", 20,0,100);
  TH1F* h2_resVsPt_pt_F = new TH1F("resVsPt2_pt_F", ";Tracking particle p_{T} [GeV]; p_{T} resolution [GeV]", 20,0,100);

  TH1F* h2_resVsPt_ptRel   = new TH1F("resVsPt2_ptRel",     ";Tracking particle p_{T} [GeV]; p_{T} resolution / p_{T}", 20,0,100);
  TH1F* h2_resVsPt_ptRel_C = new TH1F("resVsPt2_ptRel_C", ";Tracking particle p_{T} [GeV]; p_{T} resolution / p_{T}",   20,0,100);
  TH1F* h2_resVsPt_ptRel_I = new TH1F("resVsPt2_ptRel_I", ";Tracking particle p_{T} [GeV]; p_{T} resolution / p_{T}",   20,0,100);
  TH1F* h2_resVsPt_ptRel_F = new TH1F("resVsPt2_ptRel_F", ";Tracking particle p_{T} [GeV]; p_{T} resolution / p_{T}",   20,0,100);

  TH1F* h2_mresVsPt_pt   = new TH1F("mresVsPt2_pt",   ";Tracking particle p_{T} [GeV]; Mean(p_{T} residual) [GeV]", 20,0,100);
  TH1F* h2_mresVsPt_pt_C = new TH1F("mresVsPt2_pt_C", ";Tracking particle p_{T} [GeV]; Mean(p_{T} residual) [GeV]", 20,0,100);
  TH1F* h2_mresVsPt_pt_I = new TH1F("mresVsPt2_pt_I", ";Tracking particle p_{T} [GeV]; Mean(p_{T} residual) [GeV]", 20,0,100);
  TH1F* h2_mresVsPt_pt_F = new TH1F("mresVsPt2_pt_F", ";Tracking particle p_{T} [GeV]; Mean(p_{T} residual) [GeV]", 20,0,100);

  TH1F* h2_resVsPt_z0   = new TH1F("resVsPt2_z0",   ";Tracking particle p_{T} [GeV]; z_{0} resolution [cm]", 20,0,100);
  TH1F* h2_resVsPt_z0_C = new TH1F("resVsPt2_z0_C", ";Tracking particle p_{T} [GeV]; z_{0} resolution [cm]", 20,0,100);
  TH1F* h2_resVsPt_z0_I = new TH1F("resVsPt2_z0_I", ";Tracking particle p_{T} [GeV]; z_{0} resolution [cm]", 20,0,100);
  TH1F* h2_resVsPt_z0_F = new TH1F("resVsPt2_z0_F", ";Tracking particle p_{T} [GeV]; z_{0} resolution [cm]", 20,0,100);

  TH1F* h2_resVsPt_phi   = new TH1F("resVsPt2_phi",   ";Tracking particle p_{T} [GeV]; #phi resolution [rad]", 20,0,100);
  TH1F* h2_resVsPt_phi_C = new TH1F("resVsPt2_phi_C", ";Tracking particle p_{T} [GeV]; #phi resolution [rad]", 20,0,100);
  TH1F* h2_resVsPt_phi_I = new TH1F("resVsPt2_phi_I", ";Tracking particle p_{T} [GeV]; #phi resolution [rad]", 20,0,100);
  TH1F* h2_resVsPt_phi_F = new TH1F("resVsPt2_phi_F", ";Tracking particle p_{T} [GeV]; #phi resolution [rad]", 20,0,100);

  TH1F* h2_resVsPt_eta = new TH1F("resVsPt2_eta",  ";Tracking particle p_{T} [GeV]; #eta resolution",       20,0,100);
  TH1F* h2_resVsPt_d0  = new TH1F("resVsPt2_d0",   ";Tracking particle p_{T} [GeV]; d_{0} resolution [cm]", 20,0,100);

  TH1F* h2_resVsPt_pt_68    = new TH1F("resVsPt2_pt_68",   ";Tracking particle p_{T} [GeV]; p_{T} resolution [GeV]",   20,0,100);
  TH1F* h2_resVsPt_ptRel_68 = new TH1F("resVsPt2_ptRel_68",";Tracking particle p_{T} [GeV]; p_{T} resolution / p_{T}", 20,0,100);
  TH1F* h2_resVsPt_z0_68    = new TH1F("resVsPt2_z0_68",   ";Tracking particle p_{T} [GeV]; z_{0} resolution [cm]",    20,0,100);
  TH1F* h2_resVsPt_phi_68   = new TH1F("resVsPt2_phi_68",  ";Tracking particle p_{T} [GeV]; #phi resolution [rad]",    20,0,100);
  TH1F* h2_resVsPt_eta_68   = new TH1F("resVsPt2_eta_68",  ";Tracking particle p_{T} [GeV]; #eta resolution",          20,0,100);
  TH1F* h2_resVsPt_d0_68    = new TH1F("resVsPt2_d0_68",   ";Tracking particle p_{T} [GeV]; d_{0} resolution [cm]",    20,0,100);

  TH1F* h2_resVsPt_pt_90    = new TH1F("resVsPt2_pt_90",   ";Tracking particle p_{T} [GeV]; p_{T} resolution [GeV]",   20,0,100);
  TH1F* h2_resVsPt_ptRel_90 = new TH1F("resVsPt2_ptRel_90",";Tracking particle p_{T} [GeV]; p_{T} resolution / p_{T}", 20,0,100);
  TH1F* h2_resVsPt_z0_90    = new TH1F("resVsPt2_z0_90",   ";Tracking particle p_{T} [GeV]; z_{0} resolution [cm]",    20,0,100);
  TH1F* h2_resVsPt_phi_90   = new TH1F("resVsPt2_phi_90",  ";Tracking particle p_{T} [GeV]; #phi resolution [rad]",    20,0,100);
  TH1F* h2_resVsPt_eta_90   = new TH1F("resVsPt2_eta_90",  ";Tracking particle p_{T} [GeV]; #eta resolution",          20,0,100);
  TH1F* h2_resVsPt_d0_90    = new TH1F("resVsPt2_d0_90",   ";Tracking particle p_{T} [GeV]; d_{0} resolution [cm]",    20,0,100);

  TH1F* h2_resVsPt_pt_99    = new TH1F("resVsPt2_pt_99",   ";Tracking particle p_{T} [GeV]; p_{T} resolution [GeV]",   20,0,100);
  TH1F* h2_resVsPt_ptRel_99 = new TH1F("resVsPt2_ptRel_99",";Tracking particle p_{T} [GeV]; p_{T} resolution / p_{T}", 20,0,100);
  TH1F* h2_resVsPt_z0_99    = new TH1F("resVsPt2_z0_99",   ";Tracking particle p_{T} [GeV]; z_{0} resolution [cm]",    20,0,100);
  TH1F* h2_resVsPt_phi_99   = new TH1F("resVsPt2_phi_99",  ";Tracking particle p_{T} [GeV]; #phi resolution [rad]",    20,0,100);
  TH1F* h2_resVsPt_eta_99   = new TH1F("resVsPt2_eta_99",  ";Tracking particle p_{T} [GeV]; #eta resolution",          20,0,100);
  TH1F* h2_resVsPt_d0_99    = new TH1F("resVsPt2_d0_99",   ";Tracking particle p_{T} [GeV]; d_{0} resolution [cm]",    20,0,100);

  TH1F* h2_resVsPt_pt_L_68    = new TH1F("resVsPt2_pt_L_68",   ";Tracking particle p_{T} [GeV]; p_{T} resolution [GeV]",   nRANGE_L,2,8);
  TH1F* h2_resVsPt_ptRel_L_68 = new TH1F("resVsPt2_ptRel_L_68",";Tracking particle p_{T} [GeV]; p_{T} resolution / p_{T}", nRANGE_L,2,8);
  TH1F* h2_resVsPt_z0_L_68    = new TH1F("resVsPt2_z0_L_68",   ";Tracking particle p_{T} [GeV]; z_{0} resolution [cm]",    nRANGE_L,2,8);
  TH1F* h2_resVsPt_phi_L_68   = new TH1F("resVsPt2_phi_L_68",  ";Tracking particle p_{T} [GeV]; #phi resolution [rad]",    nRANGE_L,2,8);
  TH1F* h2_resVsPt_eta_L_68   = new TH1F("resVsPt2_eta_L_68",  ";Tracking particle p_{T} [GeV]; #eta resolution",          nRANGE_L,2,8);
  TH1F* h2_resVsPt_d0_L_68    = new TH1F("resVsPt2_d0_L_68",   ";Tracking particle p_{T} [GeV]; d_{0} resolution [cm]",    nRANGE_L,2,8);

  TH1F* h2_resVsPt_pt_L_90    = new TH1F("resVsPt2_pt_L_90",   ";Tracking particle p_{T} [GeV]; p_{T} resolution [GeV]",   nRANGE_L,2,8);
  TH1F* h2_resVsPt_ptRel_L_90 = new TH1F("resVsPt2_ptRel_L_90",";Tracking particle p_{T} [GeV]; p_{T} resolution / p_{T}", nRANGE_L,2,8);
  TH1F* h2_resVsPt_z0_L_90    = new TH1F("resVsPt2_z0_L_90",   ";Tracking particle p_{T} [GeV]; z_{0} resolution [cm]",    nRANGE_L,2,8);
  TH1F* h2_resVsPt_phi_L_90   = new TH1F("resVsPt2_phi_L_90",  ";Tracking particle p_{T} [GeV]; #phi resolution [rad]",    nRANGE_L,2,8);
  TH1F* h2_resVsPt_eta_L_90   = new TH1F("resVsPt2_eta_L_90",  ";Tracking particle p_{T} [GeV]; #eta resolution",          nRANGE_L,2,8);
  TH1F* h2_resVsPt_d0_L_90    = new TH1F("resVsPt2_d0_L_90",   ";Tracking particle p_{T} [GeV]; d_{0} resolution [cm]",    nRANGE_L,2,8);

  TH1F* h2_resVsPt_pt_L_99    = new TH1F("resVsPt2_pt_L_99",   ";Tracking particle p_{T} [GeV]; p_{T} resolution [cm]",    nRANGE_L,2,8);
  TH1F* h2_resVsPt_ptRel_L_99 = new TH1F("resVsPt2_ptRel_L_99",";Tracking particle p_{T} [GeV]; p_{T} resolution / p_{T}", nRANGE_L,2,8);
  TH1F* h2_resVsPt_z0_L_99    = new TH1F("resVsPt2_z0_L_99",   ";Tracking particle p_{T} [GeV]; z_{0} resolution [cm]",    nRANGE_L,2,8);
  TH1F* h2_resVsPt_phi_L_99   = new TH1F("resVsPt2_phi_L_99",  ";Tracking particle p_{T} [GeV]; #phi resolution [rad]",    nRANGE_L,2,8);
  TH1F* h2_resVsPt_eta_L_99   = new TH1F("resVsPt2_eta_L_99",  ";Tracking particle p_{T} [GeV]; #eta resolution",          nRANGE_L,2,8);
  TH1F* h2_resVsPt_d0_L_99    = new TH1F("resVsPt2_d0_L_99",   ";Tracking particle p_{T} [GeV]; d_{0} resolution [cm]",    nRANGE_L,2,8);

  for (int i=0; i<nRANGE; i++) {
    // set bin content and error
    h2_resVsPt_pt  ->SetBinContent(i+1, h_resVsPt_pt[i]  ->GetRMS());
    h2_resVsPt_pt  ->SetBinError(  i+1, h_resVsPt_pt[i]  ->GetRMSError());
    h2_resVsPt_pt_C->SetBinContent(i+1, h_resVsPt_pt_C[i]->GetRMS());
    h2_resVsPt_pt_C->SetBinError(  i+1, h_resVsPt_pt_C[i]->GetRMSError());
    h2_resVsPt_pt_I->SetBinContent(i+1, h_resVsPt_pt_I[i]->GetRMS());
    h2_resVsPt_pt_I->SetBinError(  i+1, h_resVsPt_pt_I[i]->GetRMSError());
    h2_resVsPt_pt_F->SetBinContent(i+1, h_resVsPt_pt_F[i]->GetRMS());
    h2_resVsPt_pt_F->SetBinError(  i+1, h_resVsPt_pt_F[i]->GetRMSError());

    h2_resVsPt_ptRel  ->SetBinContent(i+1, h_resVsPt_ptRel[i]  ->GetRMS());
    h2_resVsPt_ptRel  ->SetBinError(  i+1, h_resVsPt_ptRel[i]  ->GetRMSError());
    h2_resVsPt_ptRel_C->SetBinContent(i+1, h_resVsPt_ptRel_C[i]->GetRMS());
    h2_resVsPt_ptRel_C->SetBinError(  i+1, h_resVsPt_ptRel_C[i]->GetRMSError());
    h2_resVsPt_ptRel_I->SetBinContent(i+1, h_resVsPt_ptRel_I[i]->GetRMS());
    h2_resVsPt_ptRel_I->SetBinError(  i+1, h_resVsPt_ptRel_I[i]->GetRMSError());
    h2_resVsPt_ptRel_F->SetBinContent(i+1, h_resVsPt_ptRel_F[i]->GetRMS());
    h2_resVsPt_ptRel_F->SetBinError(  i+1, h_resVsPt_ptRel_F[i]->GetRMSError());

    h2_mresVsPt_pt  ->SetBinContent(i+1, h_resVsPt_pt[i]  ->GetMean());
    h2_mresVsPt_pt  ->SetBinError(  i+1, h_resVsPt_pt[i]  ->GetMeanError());
    h2_mresVsPt_pt_C->SetBinContent(i+1, h_resVsPt_pt_C[i]->GetMean());
    h2_mresVsPt_pt_C->SetBinError(  i+1, h_resVsPt_pt_C[i]->GetMeanError());
    h2_mresVsPt_pt_I->SetBinContent(i+1, h_resVsPt_pt_I[i]->GetMean());
    h2_mresVsPt_pt_I->SetBinError(  i+1, h_resVsPt_pt_I[i]->GetMeanError());
    h2_mresVsPt_pt_F->SetBinContent(i+1, h_resVsPt_pt_F[i]->GetMean());
    h2_mresVsPt_pt_F->SetBinError(  i+1, h_resVsPt_pt_F[i]->GetMeanError());

    h2_resVsPt_z0  ->SetBinContent(i+1, h_resVsPt_z0[i]  ->GetRMS());
    h2_resVsPt_z0  ->SetBinError(  i+1, h_resVsPt_z0[i]  ->GetRMSError());
    h2_resVsPt_z0_C->SetBinContent(i+1, h_resVsPt_z0_C[i]->GetRMS());
    h2_resVsPt_z0_C->SetBinError(  i+1, h_resVsPt_z0_C[i]->GetRMSError());
    h2_resVsPt_z0_I->SetBinContent(i+1, h_resVsPt_z0_I[i]->GetRMS());
    h2_resVsPt_z0_I->SetBinError(  i+1, h_resVsPt_z0_I[i]->GetRMSError());
    h2_resVsPt_z0_F->SetBinContent(i+1, h_resVsPt_z0_F[i]->GetRMS());
    h2_resVsPt_z0_F->SetBinError(  i+1, h_resVsPt_z0_F[i]->GetRMSError());

    h2_resVsPt_phi  ->SetBinContent(i+1, h_resVsPt_phi[i]  ->GetRMS());
    h2_resVsPt_phi  ->SetBinError(  i+1, h_resVsPt_phi[i]  ->GetRMSError());
    h2_resVsPt_phi_C->SetBinContent(i+1, h_resVsPt_phi_C[i]->GetRMS());
    h2_resVsPt_phi_C->SetBinError(  i+1, h_resVsPt_phi_C[i]->GetRMSError());
    h2_resVsPt_phi_I->SetBinContent(i+1, h_resVsPt_phi_I[i]->GetRMS());
    h2_resVsPt_phi_I->SetBinError(  i+1, h_resVsPt_phi_I[i]->GetRMSError());
    h2_resVsPt_phi_F->SetBinContent(i+1, h_resVsPt_phi_F[i]->GetRMS());
    h2_resVsPt_phi_F->SetBinError(  i+1, h_resVsPt_phi_F[i]->GetRMSError());

    h2_resVsPt_eta->SetBinContent(i+1, h_resVsPt_eta[i]->GetRMS());
    h2_resVsPt_eta->SetBinError(  i+1, h_resVsPt_eta[i]->GetRMSError());

    h2_resVsPt_d0->SetBinContent(i+1, h_resVsPt_d0[i]->GetRMS());
    h2_resVsPt_d0->SetBinError(  i+1, h_resVsPt_d0[i]->GetRMSError());

    h2_resVsPt_pt_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_pt[i], 0.68 ));
    h2_resVsPt_pt_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_pt[i], 0.90 ));
    h2_resVsPt_pt_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_pt[i], 0.99 ));

    h2_resVsPt_ptRel_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_ptRel[i], 0.68 ));
    h2_resVsPt_ptRel_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_ptRel[i], 0.90 ));
    h2_resVsPt_ptRel_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_ptRel[i], 0.99 ));

    h2_resVsPt_eta_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_eta[i], 0.68 ));
    h2_resVsPt_eta_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_eta[i], 0.90 ));
    h2_resVsPt_eta_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_eta[i], 0.99 ));

    h2_resVsPt_z0_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_z0[i], 0.68 ));
    h2_resVsPt_z0_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_z0[i], 0.90 ));
    h2_resVsPt_z0_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_z0[i], 0.99 ));

    h2_resVsPt_phi_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_phi[i], 0.68 ));
    h2_resVsPt_phi_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_phi[i], 0.90 ));
    h2_resVsPt_phi_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_phi[i], 0.99 ));

    h2_resVsPt_d0_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_d0[i], 0.68 ));
    h2_resVsPt_d0_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_d0[i], 0.90 ));
    h2_resVsPt_d0_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_d0[i], 0.99 ));
  }
  
  for (int i=0; i<nRANGE_L; i++) {
    h2_resVsPt_pt_L_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_pt_L[i], 0.68 ));
    h2_resVsPt_pt_L_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_pt_L[i], 0.90 ));
    h2_resVsPt_pt_L_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_pt_L[i], 0.99 ));

    h2_resVsPt_ptRel_L_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_ptRel_L[i], 0.68 ));
    h2_resVsPt_ptRel_L_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_ptRel_L[i], 0.90 ));
    h2_resVsPt_ptRel_L_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_ptRel_L[i], 0.99 ));

    h2_resVsPt_eta_L_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_eta_L[i], 0.68 ));
    h2_resVsPt_eta_L_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_eta_L[i], 0.90 ));
    h2_resVsPt_eta_L_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_eta_L[i], 0.99 ));

    h2_resVsPt_z0_L_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_z0_L[i], 0.68 ));
    h2_resVsPt_z0_L_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_z0_L[i], 0.90 ));
    h2_resVsPt_z0_L_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_z0_L[i], 0.99 ));

    h2_resVsPt_phi_L_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_phi_L[i], 0.68 ));
    h2_resVsPt_phi_L_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_phi_L[i], 0.90 ));
    h2_resVsPt_phi_L_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_phi_L[i], 0.99 ));

    h2_resVsPt_d0_L_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_d0_L[i], 0.68 ));
    h2_resVsPt_d0_L_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_d0_L[i], 0.90 ));
    h2_resVsPt_d0_L_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPt_d0_L[i], 0.99 ));
  }

  // resolution vs. eta histograms
  TH1F* h2_resVsEta_eta   = new TH1F("resVsEta2_eta",   ";Tracking particle |#eta|; #eta resolution", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_eta_L = new TH1F("resVsEta2_eta_L", ";Tracking particle |#eta|; #eta resolution", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_eta_H = new TH1F("resVsEta2_eta_H", ";Tracking particle |#eta|; #eta resolution", nETARANGE,0,2.4);

  TH1F* h2_mresVsEta_eta   = new TH1F("mresVsEta2_eta",   ";Tracking particle |#eta|; Mean(#eta residual)", nETARANGE,0,2.4);
  TH1F* h2_mresVsEta_eta_L = new TH1F("mresVsEta2_eta_L", ";Tracking particle |#eta|; Mean(#eta residual)", nETARANGE,0,2.4);
  TH1F* h2_mresVsEta_eta_H = new TH1F("mresVsEta2_eta_H", ";Tracking particle |#eta|; Mean(#eta residual)", nETARANGE,0,2.4);

  TH1F* h2_resVsEta_z0   = new TH1F("resVsEta_z0",   ";Tracking particle |#eta|; z_{0} resolution [cm]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_z0_L = new TH1F("resVsEta_z0_L", ";Tracking particle |#eta|; z_{0} resolution [cm]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_z0_H = new TH1F("resVsEta_z0_H", ";Tracking particle |#eta|; z_{0} resolution [cm]", nETARANGE,0,2.4);

  TH1F* h2_resVsEta_phi   = new TH1F("resVsEta_phi",   ";Tracking particle |#eta|; #phi resolution [rad]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_phi_L = new TH1F("resVsEta_phi_L", ";Tracking particle |#eta|; #phi resolution [rad]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_phi_H = new TH1F("resVsEta_phi_H", ";Tracking particle |#eta|; #phi resolution [rad]", nETARANGE,0,2.4);

  TH1F* h2_resVsEta_ptRel   = new TH1F("resVsEta_ptRel", ";Tracking particle |#eta|; p_{T} resolution / p_{T}", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_ptRel_L = new TH1F("resVsEta_ptRel_L", ";Tracking particle |#eta|; p_{T} resolution / p_{T}", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_ptRel_H = new TH1F("resVsEta_ptRel_H", ";Tracking particle |#eta|; p_{T} resolution / p_{T}", nETARANGE,0,2.4);

  TH1F* h2_resVsEta_d0  = new TH1F("resVsEta2_d0",  ";Tracking particle |#eta|; d_{0} resolution [cm]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_d0_L  = new TH1F("resVsEta2_d0_L",  ";Tracking particle |#eta|; d_{0} resolution [cm]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_d0_H  = new TH1F("resVsEta2_d0_H",  ";Tracking particle |#eta|; d_{0} resolution [cm]", nETARANGE,0,2.4);


  TH1F* h2_resVsEta_eta_68  = new TH1F("resVsEta_eta_68",   ";Tracking particle |#eta|; #eta resolution", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_eta_90   = new TH1F("resVsEta_eta_90",   ";Tracking particle |#eta|; #eta resolution", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_eta_99   = new TH1F("resVsEta_eta_99",   ";Tracking particle |#eta|; #eta resolution", nETARANGE,0,2.4);

  TH1F* h2_resVsEta_z0_68   = new TH1F("resVsEta_z0_68",   ";Tracking particle |#eta|; z_{0} resolution [cm]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_z0_90   = new TH1F("resVsEta_z0_90",   ";Tracking particle |#eta|; z_{0} resolution [cm]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_z0_99   = new TH1F("resVsEta_z0_99",   ";Tracking particle |#eta|; z_{0} resolution [cm]", nETARANGE,0,2.4);

  TH1F* h2_resVsEta_phi_68   = new TH1F("resVsEta_phi_68",   ";Tracking particle |#eta|; #phi resolution [rad]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_phi_90   = new TH1F("resVsEta_phi_90",   ";Tracking particle |#eta|; #phi resolution [rad]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_phi_99   = new TH1F("resVsEta_phi_99",   ";Tracking particle |#eta|; #phi resolution [rad]", nETARANGE,0,2.4);

  TH1F* h2_resVsEta_ptRel_68   = new TH1F("resVsEta_ptRel_68",   ";Tracking particle |#eta|; p_{T} resolution / p_{T}", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_ptRel_90   = new TH1F("resVsEta_ptRel_90",   ";Tracking particle |#eta|; p_{T} resolution / p_{T}", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_ptRel_99   = new TH1F("resVsEta_ptRel_99",   ";Tracking particle |#eta|; p_{T} resolution / p_{T}", nETARANGE,0,2.4);
  
  TH1F* h2_resVsEta_d0_68   = new TH1F("resVsEta_d0_68",   ";Tracking particle |#eta|; d_{0} resolution [cm]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_d0_90   = new TH1F("resVsEta_d0_90",   ";Tracking particle |#eta|; d_{0} resolution [cm]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_d0_99   = new TH1F("resVsEta_d0_99",   ";Tracking particle |#eta|; d_{0} resolution [cm]", nETARANGE,0,2.4);

  TH1F* h2_resVsEta_eta_L_68   = new TH1F("resVsEta_eta_L_68",  ";Tracking particle |#eta|; #eta resolution", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_z0_L_68    = new TH1F("resVsEta_z0_L_68",   ";Tracking particle |#eta|; z_{0} resolution [cm]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_phi_L_68   = new TH1F("resVsEta_phi_L_68",  ";Tracking particle |#eta|; #phi resolution [rad]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_ptRel_L_68 = new TH1F("resVsEta_ptRel_L_68",";Tracking particle |#eta|; p_{T} resolution / p_{T}", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_d0_L_68    = new TH1F("resVsEta_d0_L_68",   ";Tracking particle |#eta|; d_{0} resolution [cm]", nETARANGE,0,2.4);

  TH1F* h2_resVsEta_eta_L_90   = new TH1F("resVsEta_eta_L_90",  ";Tracking particle |#eta|; #eta resolution", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_z0_L_90    = new TH1F("resVsEta_z0_L_90",   ";Tracking particle |#eta|; z_{0} resolution [cm]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_phi_L_90   = new TH1F("resVsEta_phi_L_90",  ";Tracking particle |#eta|; #phi resolution [rad]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_ptRel_L_90 = new TH1F("resVsEta_ptRel_L_90",";Tracking particle |#eta|; p_{T} resolution / p_{T}", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_d0_L_90    = new TH1F("resVsEta_d0_L_90",   ";Tracking particle |#eta|; d_{0} resolution [cm]", nETARANGE,0,2.4);

  TH1F* h2_resVsEta_eta_L_99   = new TH1F("resVsEta_eta_L_99",  ";Tracking particle |#eta|; #eta resolution", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_z0_L_99    = new TH1F("resVsEta_z0_L_99",   ";Tracking particle |#eta|; z_{0} resolution [cm]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_phi_L_99   = new TH1F("resVsEta_phi_L_99",  ";Tracking particle |#eta|; #phi resolution [rad]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_ptRel_L_99 = new TH1F("resVsEta_ptRel_L_99",";Tracking particle |#eta|; p_{T} resolution / p_{T}", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_d0_L_99    = new TH1F("resVsEta_d0_L_99",   ";Tracking particle |#eta|; d_{0} resolution [cm]", nETARANGE,0,2.4);

  TH1F* h2_resVsEta_eta_H_68   = new TH1F("resVsEta_eta_H_68",  ";Tracking particle |#eta|; #eta resolution", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_z0_H_68    = new TH1F("resVsEta_z0_H_68",   ";Tracking particle |#eta|; z_{0} resolution [cm]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_phi_H_68   = new TH1F("resVsEta_phi_H_68",  ";Tracking particle |#eta|; #phi resolution [rad]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_ptRel_H_68 = new TH1F("resVsEta_ptRel_H_68",";Tracking particle |#eta|; p_{T} resolution / p_{T}", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_d0_H_68    = new TH1F("resVsEta_d0_H_68",   ";Tracking particle |#eta|; d_{0} resolution [cm]", nETARANGE,0,2.4);

  TH1F* h2_resVsEta_eta_H_90   = new TH1F("resVsEta_eta_H_90",  ";Tracking particle |#eta|; #eta resolution", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_z0_H_90    = new TH1F("resVsEta_z0_H_90",   ";Tracking particle |#eta|; z_{0} resolution [cm]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_phi_H_90   = new TH1F("resVsEta_phi_H_90",  ";Tracking particle |#eta|; #phi resolution [rad]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_ptRel_H_90 = new TH1F("resVsEta_ptRel_H_90",";Tracking particle |#eta|; p_{T} resolution / p_{T}", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_d0_H_90    = new TH1F("resVsEta_d0_H_90",   ";Tracking particle |#eta|; d_{0} resolution [cm]", nETARANGE,0,2.4);

  TH1F* h2_resVsEta_eta_H_99   = new TH1F("resVsEta_eta_H_99",  ";Tracking particle |#eta|; #eta resolution", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_z0_H_99    = new TH1F("resVsEta_z0_H_99",   ";Tracking particle |#eta|; z_{0} resolution [cm]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_phi_H_99   = new TH1F("resVsEta_phi_H_99",  ";Tracking particle |#eta|; #phi resolution [rad]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_ptRel_H_99 = new TH1F("resVsEta_ptRel_H_99",";Tracking particle |#eta|; p_{T} resolution / p_{T}", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_d0_H_99    = new TH1F("resVsEta_d0_H_99",   ";Tracking particle |#eta|; d_{0} resolution [cm]", nETARANGE,0,2.4);

  // resolution vs. eta histograms (gaussian fit)
  TH1F* h3_resVsEta_eta_L = new TH1F("resVsEta_eta_L_gaus", ";|#eta|; #sigma(#eta)", nETARANGE,0,2.4);
  TH1F* h3_resVsEta_eta_H = new TH1F("resVsEta_eta_H_gaus", ";|#eta|; #sigma(#eta)", nETARANGE,0,2.4);

  TH1F* h3_resVsEta_z0_L = new TH1F("resVsEta_z0_L_gaus", ";|#eta|; #sigma(z_{0}) [cm]", nETARANGE,0,2.4);
  TH1F* h3_resVsEta_z0_H = new TH1F("resVsEta_z0_H_gaus", ";|#eta|; #sigma(z_{0}) [cm]", nETARANGE,0,2.4);

  TH1F* h3_resVsEta_phi_L = new TH1F("resVsEta_phi_L_gaus", ";|#eta|; #sigma(#phi) [rad]", nETARANGE,0,2.4);
  TH1F* h3_resVsEta_phi_H = new TH1F("resVsEta_phi_H_gaus", ";|#eta|; #sigma(#phi) [rad]", nETARANGE,0,2.4);

  TH1F* h3_resVsEta_ptRel_L = new TH1F("resVsEta_ptRel_L_gaus", ";|#eta|; #sigma(p_{T}) / p_{T}", nETARANGE,0,2.4);
  TH1F* h3_resVsEta_ptRel_H = new TH1F("resVsEta_ptRel_H_gaus", ";|#eta|; #sigma(p_{T}) / p_{T}", nETARANGE,0,2.4);

  TString fitdir = "FitResults/";


  for (int i=0; i<nETARANGE; i++) {
    // set bin content and error
    h2_resVsEta_eta  ->SetBinContent(i+1, h_resVsEta_eta[i]  ->GetRMS());
    h2_resVsEta_eta  ->SetBinError(  i+1, h_resVsEta_eta[i]  ->GetRMSError());
    h2_resVsEta_eta_L->SetBinContent(i+1, h_resVsEta_eta_L[i]->GetRMS());
    h2_resVsEta_eta_L->SetBinError(  i+1, h_resVsEta_eta_L[i]->GetRMSError());
    h2_resVsEta_eta_H->SetBinContent(i+1, h_resVsEta_eta_H[i]->GetRMS());
    h2_resVsEta_eta_H->SetBinError(  i+1, h_resVsEta_eta_H[i]->GetRMSError());

    h2_mresVsEta_eta  ->SetBinContent(i+1, h_resVsEta_eta[i]  ->GetMean());
    h2_mresVsEta_eta  ->SetBinError(  i+1, h_resVsEta_eta[i]  ->GetMeanError());
    h2_mresVsEta_eta_L->SetBinContent(i+1, h_resVsEta_eta_L[i]->GetMean());
    h2_mresVsEta_eta_L->SetBinError(  i+1, h_resVsEta_eta_L[i]->GetMeanError());
    h2_mresVsEta_eta_H->SetBinContent(i+1, h_resVsEta_eta_H[i]->GetMean());
    h2_mresVsEta_eta_H->SetBinError(  i+1, h_resVsEta_eta_H[i]->GetMeanError());

    h2_resVsEta_z0   ->SetBinContent(i+1, h_resVsEta_z0[i]   ->GetRMS());
    h2_resVsEta_z0   ->SetBinError(  i+1, h_resVsEta_z0[i]   ->GetRMSError());
    h2_resVsEta_z0_L ->SetBinContent(i+1, h_resVsEta_z0_L[i] ->GetRMS());
    h2_resVsEta_z0_L ->SetBinError(  i+1, h_resVsEta_z0_L[i] ->GetRMSError());
    h2_resVsEta_z0_H ->SetBinContent(i+1, h_resVsEta_z0_H[i] ->GetRMS());
    h2_resVsEta_z0_H ->SetBinError(  i+1, h_resVsEta_z0_H[i] ->GetRMSError());

    h2_resVsEta_phi  ->SetBinContent(i+1, h_resVsEta_phi[i]  ->GetRMS());
    h2_resVsEta_phi  ->SetBinError(  i+1, h_resVsEta_phi[i]  ->GetRMSError());
    h2_resVsEta_phi_L->SetBinContent(i+1, h_resVsEta_phi_L[i]->GetRMS());
    h2_resVsEta_phi_L->SetBinError(  i+1, h_resVsEta_phi_L[i]->GetRMSError());
    h2_resVsEta_phi_H->SetBinContent(i+1, h_resVsEta_phi_H[i]->GetRMS());
    h2_resVsEta_phi_H->SetBinError(  i+1, h_resVsEta_phi_H[i]->GetRMSError());

    h2_resVsEta_ptRel  ->SetBinContent(i+1, h_resVsEta_ptRel[i]  ->GetRMS());
    h2_resVsEta_ptRel  ->SetBinError(  i+1, h_resVsEta_ptRel[i]  ->GetRMSError());
    h2_resVsEta_ptRel_L->SetBinContent(i+1, h_resVsEta_ptRel_L[i]->GetRMS());
    h2_resVsEta_ptRel_L->SetBinError(  i+1, h_resVsEta_ptRel_L[i]->GetRMSError());
    h2_resVsEta_ptRel_H->SetBinContent(i+1, h_resVsEta_ptRel_H[i]->GetRMS());
    h2_resVsEta_ptRel_H->SetBinError(  i+1, h_resVsEta_ptRel_H[i]->GetRMSError());

    h2_resVsEta_d0  ->SetBinContent(i+1, h_resVsEta_d0[i]   ->GetRMS());
    h2_resVsEta_d0  ->SetBinError(  i+1, h_resVsEta_d0[i]   ->GetRMSError());
    h2_resVsEta_d0_L->SetBinContent(i+1, h_resVsEta_d0_L[i] ->GetRMS());
    h2_resVsEta_d0_L->SetBinError(  i+1, h_resVsEta_d0_L[i] ->GetRMSError());
    h2_resVsEta_d0_H->SetBinContent(i+1, h_resVsEta_d0_H[i] ->GetRMS());
    h2_resVsEta_d0_H->SetBinError(  i+1, h_resVsEta_d0_H[i] ->GetRMSError());

    h2_resVsEta_eta_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_eta[i], 0.68 ));
    h2_resVsEta_eta_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_eta[i], 0.90 ));
    h2_resVsEta_eta_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_eta[i], 0.99 ));

    h2_resVsEta_z0_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_z0[i], 0.68 ));
    h2_resVsEta_z0_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_z0[i], 0.90 ));
    h2_resVsEta_z0_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_z0[i], 0.99 ));

    h2_resVsEta_phi_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_phi[i], 0.68 ));
    h2_resVsEta_phi_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_phi[i], 0.90 ));
    h2_resVsEta_phi_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_phi[i], 0.99 ));

    h2_resVsEta_ptRel_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_ptRel[i], 0.68 ));
    h2_resVsEta_ptRel_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_ptRel[i], 0.90 ));
    h2_resVsEta_ptRel_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_ptRel[i], 0.99 ));

    h2_resVsEta_d0_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_d0[i], 0.68 ));
    h2_resVsEta_d0_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_d0[i], 0.90 ));
    h2_resVsEta_d0_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_d0[i], 0.99 ));


    h2_resVsEta_eta_L_68  ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_eta_L[i],   0.68 ));
    h2_resVsEta_z0_L_68   ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_z0_L[i],    0.68 ));
    h2_resVsEta_phi_L_68  ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_phi_L[i],   0.68 ));
    h2_resVsEta_ptRel_L_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_ptRel_L[i], 0.68 ));
    h2_resVsEta_d0_L_68   ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_d0_L[i],    0.68 ));

    h2_resVsEta_eta_L_90  ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_eta_L[i],   0.90 ));
    h2_resVsEta_z0_L_90   ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_z0_L[i],    0.90 ));
    h2_resVsEta_phi_L_90  ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_phi_L[i],   0.90 ));
    h2_resVsEta_ptRel_L_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_ptRel_L[i], 0.90 ));
    h2_resVsEta_d0_L_90   ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_d0_L[i],    0.90 ));

    h2_resVsEta_eta_L_99  ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_eta_L[i],   0.99 ));
    h2_resVsEta_z0_L_99   ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_z0_L[i],    0.99 ));
    h2_resVsEta_phi_L_99  ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_phi_L[i],   0.99 ));
    h2_resVsEta_ptRel_L_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_ptRel_L[i], 0.99 ));
    h2_resVsEta_d0_L_99   ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_d0_L[i],    0.99 ));

    h2_resVsEta_eta_H_68  ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_eta_H[i],   0.68 ));
    h2_resVsEta_z0_H_68   ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_z0_H[i],    0.68 ));
    h2_resVsEta_phi_H_68  ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_phi_H[i],   0.68 ));
    h2_resVsEta_ptRel_H_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_ptRel_H[i], 0.68 ));
    h2_resVsEta_d0_H_68   ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_d0_H[i],    0.68 ));

    h2_resVsEta_eta_H_90  ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_eta_H[i],   0.90 ));
    h2_resVsEta_z0_H_90   ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_z0_H[i],    0.90 ));
    h2_resVsEta_phi_H_90  ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_phi_H[i],   0.90 ));
    h2_resVsEta_ptRel_H_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_ptRel_H[i], 0.90 ));
    h2_resVsEta_d0_H_90   ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_d0_H[i],    0.90 ));

    h2_resVsEta_eta_H_99  ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_eta_H[i],   0.99 ));
    h2_resVsEta_z0_H_99   ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_z0_H[i],    0.99 ));
    h2_resVsEta_phi_H_99  ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_phi_H[i],   0.99 ));
    h2_resVsEta_ptRel_H_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_ptRel_H[i], 0.99 ));
    h2_resVsEta_d0_H_99   ->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_d0_H[i],    0.99 ));



    // ---------------------------------------------------------------------------------------------------
    // gaussian fit instead
    // ---------------------------------------------------------------------------------------------------
  
    if (doGausFit) {

      TCanvas cfit;
      char text[500];
      
      float sigma = 0;
      float esigma = 0;
      TF1* fit;
      
      float rms = 0;
      float erms = 0;
      
      fit = new TF1("fit", "gaus", -0.01,0.01);
      h_resVsEta_eta_L[i]->Fit("fit","R");
      h_resVsEta_eta_L[i]->SetTitle(htitle);
      sigma  = fit->GetParameter(2);
      esigma = fit->GetParError(2);
      rms = h_resVsEta_eta_L[i]->GetRMS();
      erms = h_resVsEta_eta_L[i]->GetRMSError();
      h3_resVsEta_eta_L->SetBinContent(i+1, sigma);   
      h3_resVsEta_eta_L->SetBinError(i+1, esigma);   
      h_resVsEta_eta_L[i]->Draw();
      sprintf(text,"RMS: %.4f +/- %.4f",rms,erms);
      mySmallText(0.2,0.86,1,text);
      sprintf(text,"Fit: %.4f +/- %.4f",sigma,esigma);
      mySmallText(0.2,0.8,2,text);
      sprintf(text,"p_{T} < 5 GeV");
      mySmallText(0.2,0.7,1,text);
      cfit.SaveAs(fitdir+"resVsEta_eta_L_"+etarange[i]+".png");
      delete fit;
      
      fit = new TF1("fit", "gaus", -0.01,0.01);
      h_resVsEta_eta_H[i]->Fit("fit","R");
      h_resVsEta_eta_H[i]->SetTitle(htitle);
      sigma  = fit->GetParameter(2);
      esigma = fit->GetParError(2);
      rms = h_resVsEta_eta_H[i]->GetRMS();
      erms = h_resVsEta_eta_H[i]->GetRMSError();
      h3_resVsEta_eta_H->SetBinContent(i+1, sigma);   
      h3_resVsEta_eta_H->SetBinError(i+1, esigma);   
      h_resVsEta_eta_H[i]->Draw();
      sprintf(text,"RMS: %.4f +/- %.4f",rms,erms);
      mySmallText(0.2,0.86,1,text);
      sprintf(text,"Fit: %.4f +/- %.4f",sigma,esigma);
      mySmallText(0.2,0.8,2,text);
      sprintf(text,"p_{T} > 15 GeV");
      mySmallText(0.2,0.7,1,text);
      cfit.SaveAs(fitdir+"resVsEta_eta_H_"+etarange[i]+".png");
      delete fit;
      
      fit = new TF1("fit", "gaus", -1,1);
      h_resVsEta_z0_L[i]->Fit("fit","R");
      h_resVsEta_z0_L[i]->SetTitle(htitle);
      sigma  = fit->GetParameter(2);
      esigma = fit->GetParError(2);
      rms = h_resVsEta_z0_L[i]->GetRMS();
      erms = h_resVsEta_z0_L[i]->GetRMSError();
      h3_resVsEta_z0_L->SetBinContent(i+1, sigma);   
      h3_resVsEta_z0_L->SetBinError(i+1, esigma);   
      h_resVsEta_z0_L[i]->Draw();
      sprintf(text,"RMS: %.4f +/- %.4f",rms,erms);
      mySmallText(0.2,0.86,1,text);
      sprintf(text,"Fit: %.4f +/- %.4f",sigma,esigma);
      mySmallText(0.2,0.8,2,text);
      sprintf(text,"p_{T} < 5 GeV");
      mySmallText(0.2,0.7,1,text);
      cfit.SaveAs(fitdir+"resVsEta_z0_L_"+etarange[i]+".png");
      delete fit;
      
      fit = new TF1("fit", "gaus", -1,1);
      h_resVsEta_z0_H[i]->Fit("fit","R");
      h_resVsEta_z0_H[i]->SetTitle(htitle);
      sigma  = fit->GetParameter(2);
      esigma = fit->GetParError(2);
      rms = h_resVsEta_z0_H[i]->GetRMS();
      erms = h_resVsEta_z0_H[i]->GetRMSError();
      h3_resVsEta_z0_H->SetBinContent(i+1, sigma);   
      h3_resVsEta_z0_H->SetBinError(i+1, esigma);   
      h_resVsEta_z0_H[i]->Draw();
      sprintf(text,"RMS: %.4f +/- %.4f",rms,erms);
      mySmallText(0.2,0.86,1,text);
      sprintf(text,"Fit: %.4f +/- %.4f",sigma,esigma);
      mySmallText(0.2,0.8,2,text);
      sprintf(text,"p_{T} > 15 GeV");
      mySmallText(0.2,0.7,1,text);
      cfit.SaveAs(fitdir+"resVsEta_z0_H_"+etarange[i]+".png");
      delete fit;
      
      fit = new TF1("fit", "gaus", -0.005,0.005);
      h_resVsEta_phi_L[i]->Fit("fit","R");
      h_resVsEta_phi_L[i]->SetTitle(htitle);
      sigma  = fit->GetParameter(2);
      esigma = fit->GetParError(2);
      rms = h_resVsEta_phi_L[i]->GetRMS();
      erms = h_resVsEta_phi_L[i]->GetRMSError();
      h3_resVsEta_phi_L->SetBinContent(i+1, sigma);   
      h3_resVsEta_phi_L->SetBinError(i+1, esigma);   
      h_resVsEta_phi_L[i]->Draw();
      sprintf(text,"RMS: %.4f +/- %.4f",rms,erms);
      mySmallText(0.2,0.86,1,text);
      sprintf(text,"Fit: %.4f +/- %.4f",sigma,esigma);
      mySmallText(0.2,0.8,2,text);
      sprintf(text,"p_{T} < 5 GeV");
      mySmallText(0.2,0.7,1,text);
      cfit.SaveAs(fitdir+"resVsEta_phi_L_"+etarange[i]+".png");
      delete fit;
      
      fit = new TF1("fit", "gaus", -0.005,0.005);
      h_resVsEta_phi_H[i]->Fit("fit","R");
      h_resVsEta_phi_H[i]->SetTitle(htitle);
      sigma  = fit->GetParameter(2);
      esigma = fit->GetParError(2);
      rms = h_resVsEta_phi_H[i]->GetRMS();
      erms = h_resVsEta_phi_H[i]->GetRMSError();
      h3_resVsEta_phi_H->SetBinContent(i+1, sigma);   
      h3_resVsEta_phi_H->SetBinError(i+1, esigma);   
      h_resVsEta_phi_H[i]->Draw();
      sprintf(text,"RMS: %.4f +/- %.4f",rms,erms);
      mySmallText(0.2,0.86,1,text);
      sprintf(text,"Fit: %.4f +/- %.4f",sigma,esigma);
      mySmallText(0.2,0.8,2,text);
      sprintf(text,"p_{T} > 15 GeV");
      mySmallText(0.2,0.7,1,text);
      cfit.SaveAs(fitdir+"resVsEta_phi_H_"+etarange[i]+".png");
      delete fit;
      
      
      fit = new TF1("fit", "gaus", -0.5,0.5);
      h_resVsEta_ptRel_L[i]->Fit("fit","R");
      h_resVsEta_ptRel_L[i]->SetTitle(htitle);
      sigma  = fit->GetParameter(2);
      esigma = fit->GetParError(2);
      rms = h_resVsEta_ptRel_L[i]->GetRMS();
      erms = h_resVsEta_ptRel_L[i]->GetRMSError();
      h3_resVsEta_ptRel_L->SetBinContent(i+1, sigma);   
      h3_resVsEta_ptRel_L->SetBinError(i+1, esigma);   
      h_resVsEta_ptRel_L[i]->Draw();
      sprintf(text,"RMS: %.4f +/- %.4f",rms,erms);
      mySmallText(0.2,0.86,1,text);
      sprintf(text,"Fit: %.4f +/- %.4f",sigma,esigma);
      mySmallText(0.2,0.8,2,text);
      sprintf(text,"p_{T} < 5 GeV");
      mySmallText(0.2,0.7,1,text);
      cfit.SaveAs(fitdir+"resVsEta_ptRel_L_"+etarange[i]+".png");
      delete fit;
      
      fit = new TF1("fit", "gaus", -0.5,0.5);
      h_resVsEta_ptRel_H[i]->Fit("fit","R");
      h_resVsEta_ptRel_H[i]->SetTitle(htitle);
      sigma  = fit->GetParameter(2);
      esigma = fit->GetParError(2);
      rms = h_resVsEta_ptRel_H[i]->GetRMS();
      erms = h_resVsEta_ptRel_H[i]->GetRMSError();
      h3_resVsEta_ptRel_H->SetBinContent(i+1, sigma);   
      h3_resVsEta_ptRel_H->SetBinError(i+1, esigma);   
      h_resVsEta_ptRel_H[i]->Draw();
      sprintf(text,"RMS: %.4f +/- %.4f",rms,erms);
      mySmallText(0.2,0.86,1,text);
      sprintf(text,"Fit: %.4f +/- %.4f",sigma,esigma);
      mySmallText(0.2,0.8,2,text);
      sprintf(text,"p_{T} > 15 GeV");
      mySmallText(0.2,0.7,1,text);
      cfit.SaveAs(fitdir+"resVsEta_ptRel_H_"+etarange[i]+".png");
      delete fit;

    }//end doGausFit
    
  }

  // res vs Phi
  TH1F* h2_resVsPhi_pt_68    = new TH1F("resVsPhi2_pt_68",   ";Tracking particle #phi; p_{T} resolution [GeV]",         32,-3.2,3.2);
  TH1F* h2_resVsPhi_pt_90    = new TH1F("resVsPhi2_pt_90",   ";Tracking particle #phi; p_{T} resolution [GeV]",         32,-3.2,3.2);
  TH1F* h2_resVsPhi_pt_99    = new TH1F("resVsPhi2_pt_99",   ";Tracking particle #phi; p_{T} resolution [GeV]",         32,-3.2,3.2);
  TH1F* h2_resVsPhi_ptRel_68 = new TH1F("resVsPhi2_ptRel_68",";Tracking particle #phi; p_{T} resolution / p_{T}", 32,-3.2,3.2);
  TH1F* h2_resVsPhi_ptRel_90 = new TH1F("resVsPhi2_ptRel_90",";Tracking particle #phi; p_{T} resolution / p_{T}", 32,-3.2,3.2);
  TH1F* h2_resVsPhi_ptRel_99 = new TH1F("resVsPhi2_ptRel_99",";Tracking particle #phi; p_{T} resolution / p_{T}", 32,-3.2,3.2);

  for (int i=0; i<nPHIRANGE; i++) {
    h2_resVsPhi_pt_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPhi_pt[i], 0.68 ));
    h2_resVsPhi_pt_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPhi_pt[i], 0.90 ));
    h2_resVsPhi_pt_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPhi_pt[i], 0.99 ));

    h2_resVsPhi_ptRel_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPhi_ptRel[i], 0.68 ));
    h2_resVsPhi_ptRel_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPhi_ptRel[i], 0.90 ));
    h2_resVsPhi_ptRel_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsPhi_ptRel[i], 0.99 ));
  }

  
  // set minimum to zero
  h2_resVsPt_pt  ->SetMinimum(0);
  h2_resVsPt_pt_C->SetMinimum(0);
  h2_resVsPt_pt_I->SetMinimum(0);
  h2_resVsPt_pt_F->SetMinimum(0);

  h2_resVsPt_ptRel  ->SetMinimum(0);
  h2_resVsPt_ptRel_C->SetMinimum(0);
  h2_resVsPt_ptRel_I->SetMinimum(0);
  h2_resVsPt_ptRel_F->SetMinimum(0);

  h2_resVsPt_z0  ->SetMinimum(0);
  h2_resVsPt_z0_C->SetMinimum(0);
  h2_resVsPt_z0_I->SetMinimum(0);
  h2_resVsPt_z0_F->SetMinimum(0);

  h2_resVsPt_phi  ->SetMinimum(0);
  h2_resVsPt_phi_C->SetMinimum(0);
  h2_resVsPt_phi_I->SetMinimum(0);
  h2_resVsPt_phi_F->SetMinimum(0);

  h2_resVsPt_eta->SetMinimum(0);

  h2_resVsEta_eta  ->SetMinimum(0);
  h2_resVsEta_eta_L->SetMinimum(0);
  h2_resVsEta_eta_H->SetMinimum(0);

  h2_resVsEta_z0   ->SetMinimum(0);
  h2_resVsEta_z0_L ->SetMinimum(0);
  h2_resVsEta_z0_H ->SetMinimum(0);

  h2_resVsEta_phi  ->SetMinimum(0);
  h2_resVsEta_phi_L->SetMinimum(0);
  h2_resVsEta_phi_H->SetMinimum(0);

  h2_resVsEta_ptRel  ->SetMinimum(0);
  h2_resVsEta_ptRel_L->SetMinimum(0);
  h2_resVsEta_ptRel_H->SetMinimum(0);

  h2_resVsPt_d0  ->SetMinimum(0);
  h2_resVsEta_d0 ->SetMinimum(0);
  h2_resVsEta_d0_L ->SetMinimum(0);
  h2_resVsEta_d0_H ->SetMinimum(0);


  // plots overlaying 68, 90, 99% confidence levels]
  float max_eta_ptRel = 0.2;
  float max_pt_ptRel  = 0.2;
  float max_pt_pt     = 20;
  float max_z0        = 2.0;
  float max_phi       = 0.01;
  float max_eta       = 0.03;

  if (type.Contains("Electron")) {
    max_pt_ptRel  = 1.0;
    max_eta_ptRel = 1.0;
    max_phi       = 0.1;
  }

  // makeResidualIntervalPlot will save the individual plots to the root file
  makeResidualIntervalPlot( type, DIR, "resVsPt_ptRel", h2_resVsPt_ptRel_68, h2_resVsPt_ptRel_90, h2_resVsPt_ptRel_99, 0, max_pt_ptRel, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsPt_pt", h2_resVsPt_pt_68, h2_resVsPt_pt_90, h2_resVsPt_pt_99, 0, max_pt_pt, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsPt_z0", h2_resVsPt_z0_68, h2_resVsPt_z0_90, h2_resVsPt_z0_99, 0, max_z0, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsPt_phi", h2_resVsPt_phi_68, h2_resVsPt_phi_90, h2_resVsPt_phi_99, 0, max_phi, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsPt_eta", h2_resVsPt_eta_68, h2_resVsPt_eta_90, h2_resVsPt_eta_99, 0, max_eta, htitle );
  //makeResidualIntervalPlot( type, DIR, "resVsPt_d0", h2_resVsPt_d0_68, h2_resVsPt_d0_90, h2_resVsPt_d0_99, 0, 0.02, htitle );
  
  makeResidualIntervalPlot( type, DIR, "resVsPt_L_ptRel", h2_resVsPt_ptRel_L_68, h2_resVsPt_ptRel_L_90, h2_resVsPt_ptRel_L_99, 0, max_pt_ptRel, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsPt_L_pt", h2_resVsPt_pt_L_68, h2_resVsPt_pt_L_90, h2_resVsPt_pt_L_99, 0, max_pt_pt, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsPt_L_z0", h2_resVsPt_z0_L_68, h2_resVsPt_z0_L_90, h2_resVsPt_z0_L_99, 0, max_z0, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsPt_L_phi", h2_resVsPt_phi_L_68, h2_resVsPt_phi_L_90, h2_resVsPt_phi_L_99, 0, max_phi, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsPt_L_eta", h2_resVsPt_eta_L_68, h2_resVsPt_eta_L_90, h2_resVsPt_eta_L_99, 0, max_eta, htitle );
  //makeResidualIntervalPlot( type, DIR, "resVsPt_L_d0", h2_resVsPt_d0_L_68, h2_resVsPt_d0_L_90, h2_resVsPt_d0_L_99, 0, 0.02, htitle );
  
  makeResidualIntervalPlot( type, DIR, "resVsEta_eta", h2_resVsEta_eta_68, h2_resVsEta_eta_90, h2_resVsEta_eta_99, 0, max_eta, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsEta_z0", h2_resVsEta_z0_68, h2_resVsEta_z0_90, h2_resVsEta_z0_99, 0, max_z0, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsEta_phi", h2_resVsEta_phi_68, h2_resVsEta_phi_90, h2_resVsEta_phi_99, 0, max_phi, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsEta_ptRel", h2_resVsEta_ptRel_68, h2_resVsEta_ptRel_90, h2_resVsEta_ptRel_99, 0, max_eta_ptRel, htitle );
  //makeResidualIntervalPlot( type, DIR, "resVsEta_d0", h2_resVsEta_d0_68, h2_resVsEta_d0_90, h2_resVsEta_d0_99, 0, 0.02 );
  
  makeResidualIntervalPlot( type, DIR, "resVsEta_L_eta", h2_resVsEta_eta_L_68, h2_resVsEta_eta_L_90, h2_resVsEta_eta_L_99, 0, max_eta, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsEta_L_z0", h2_resVsEta_z0_L_68, h2_resVsEta_z0_L_90, h2_resVsEta_z0_L_99, 0, max_z0, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsEta_L_phi", h2_resVsEta_phi_L_68, h2_resVsEta_phi_L_90, h2_resVsEta_phi_L_99, 0, max_phi, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsEta_L_ptRel", h2_resVsEta_ptRel_L_68, h2_resVsEta_ptRel_L_90, h2_resVsEta_ptRel_L_99, 0, max_eta_ptRel, htitle );
  //makeResidualIntervalPlot( type, DIR, "resVsEta_L_d0", h2_resVsEta_d0_L_68, h2_resVsEta_d0_L_90, h2_resVsEta_d0_L_99, 0, 0.02, htitle );
  
  makeResidualIntervalPlot( type, DIR, "resVsEta_H_eta", h2_resVsEta_eta_H_68, h2_resVsEta_eta_H_90, h2_resVsEta_eta_H_99, 0, max_eta, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsEta_H_z0", h2_resVsEta_z0_H_68, h2_resVsEta_z0_H_90, h2_resVsEta_z0_H_99, 0, max_z0, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsEta_H_phi", h2_resVsEta_phi_H_68, h2_resVsEta_phi_H_90, h2_resVsEta_phi_H_99, 0, max_phi, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsEta_H_ptRel", h2_resVsEta_ptRel_H_68, h2_resVsEta_ptRel_H_90, h2_resVsEta_ptRel_H_99, 0, max_eta_ptRel, htitle );
  //makeResidualIntervalPlot( type, DIR, "resVsEta_H_d0", h2_resVsEta_d0_H_68, h2_resVsEta_d0_H_90, h2_resVsEta_d0_H_99, 0, 0.02, htitle );
  
  makeResidualIntervalPlot( type, DIR, "resVsPhi_ptRel", h2_resVsPhi_ptRel_68, h2_resVsPhi_ptRel_90, h2_resVsPhi_ptRel_99, 0, 0.5, htitle );
  makeResidualIntervalPlot( type, DIR, "resVsPhi_pt", h2_resVsPhi_pt_68, h2_resVsPhi_pt_90, h2_resVsPhi_pt_99, 0, 20, htitle );
  
  

  // ----------------------------------------------------------------------------------------------------------
  // resoultion vs pt
  // ----------------------------------------------------------------------------------------------------------
  if(doResolutionPlots){
    h2_resVsPt_pt_90->SetTitle(htitle);
    h2_resVsPt_pt_90->SetMinimum(0);
    h2_resVsPt_pt_90->SetMarkerStyle(20);
    h2_resVsPt_pt_90->Draw("p");
    c.SaveAs(DIR+type+"_resVsPt_pt_90.png");

    h2_resVsPt_ptRel_90->SetTitle(htitle);
    h2_resVsPt_ptRel_90->SetMinimum(0);
    h2_resVsPt_ptRel_90->SetMarkerStyle(20);
    h2_resVsPt_ptRel_90->Draw("p");
    c.SaveAs(DIR+type+"_resVsPt_ptRel_90.png");

    h2_resVsPt_z0_90->SetTitle(htitle);
    h2_resVsPt_z0_90->SetMinimum(0);
    h2_resVsPt_z0_90->SetMarkerStyle(20);
    h2_resVsPt_z0_90->Draw("p");
    c.SaveAs(DIR+type+"_resVsPt_z0_90.png");

    h2_resVsPt_phi_90->SetTitle(htitle);
    h2_resVsPt_phi_90->SetMinimum(0);
    h2_resVsPt_phi_90->SetMarkerStyle(20);
    h2_resVsPt_phi_90->Draw("p");
    c.SaveAs(DIR+type+"_resVsPt_phi_90.png");

    h2_resVsPt_eta_90->SetTitle(htitle);
    h2_resVsPt_eta_90->SetMinimum(0);
    h2_resVsPt_eta_90->SetMarkerStyle(20);
    h2_resVsPt_eta_90->Draw("p");
    c.SaveAs(DIR+type+"_resVsPt_eta_90.png");

    /*
    h2_resVsPt_phi_90->SetTitle(htitle);
    h2_resVsPt_phi_90->SetMinimum(0);
    h2_resVsPt_d0_90->SetMarkerStyle(20);
    h2_resVsPt_d0_90->Draw("p");
    c.SaveAs(DIR+type+"_resVsPt_d0_90.pdf");
    c.SaveAs(DIR+type+"_resVsPt_d0_90.png");
    */

    if (doDetailedPlots) {
      h2_resVsPt_eta->Write();

      h2_resVsPt_pt  ->Write();
      h2_resVsPt_pt_C->Write();    
      h2_resVsPt_pt_I->Write();
      h2_resVsPt_pt_F->Write();

      h2_resVsPt_ptRel  ->Write();
      h2_resVsPt_ptRel_C->Write();
      h2_resVsPt_ptRel_I->Write();
      h2_resVsPt_ptRel_F->Write();

      h2_mresVsPt_pt  ->Write();
      h2_mresVsPt_pt_C->Write();
      h2_mresVsPt_pt_I->Write();
      h2_mresVsPt_pt_F->Write();

      h2_resVsPt_d0->Write();

      h2_resVsPt_z0_C->Write();
      h2_resVsPt_z0_I->Write();
      h2_resVsPt_z0_F->Write();

      h2_resVsPt_phi  ->Write();
      h2_resVsPt_phi_C->Write();
      h2_resVsPt_phi_I->Write();
      h2_resVsPt_phi_F->Write();
    }
  
  }

  // ----------------------------------------------------------------------------------------------------------
  // resolution vs eta
  // ----------------------------------------------------------------------------------------------------------
  if(doResolutionPlots){
    h2_resVsEta_eta_90->SetTitle(htitle);
    h2_resVsEta_eta_90->SetMinimum(0);
    h2_resVsEta_eta_90->SetMarkerStyle(20);
    h2_resVsEta_eta_90->Draw("p");
    c.SaveAs(DIR+type+"_resVsEta_eta_90.png");

    h2_resVsEta_eta_68->SetTitle(htitle);
    h2_resVsEta_eta_68->SetMinimum(0);
    h2_resVsEta_eta_68->SetMarkerStyle(20);
    h2_resVsEta_eta_68->Draw("p");
    c.SaveAs(DIR+type+"_resVsEta_eta_68.png");

    if (doDetailedPlots) {
      h2_resVsEta_eta_L_90->SetTitle(htitle);
      h2_resVsEta_eta_L_90->Draw("p");
      sprintf(ctxt,"p_{T} < 8 GeV");
      mySmallText(0.22,0.82,1,ctxt);
      c.SaveAs(DIR+type+"_resVsEta_eta_L_90.png");
      
      h2_resVsEta_eta_H_90->SetTitle(htitle);
      h2_resVsEta_eta_H_90->Draw("p");
      sprintf(ctxt,"p_{T} > 8 GeV");
      mySmallText(0.22,0.82,1,ctxt);
      c.SaveAs(DIR+type+"_resVsEta_eta_H_90.png");
    }

    h2_resVsEta_z0_90->SetTitle(htitle);
    h2_resVsEta_z0_90->SetMinimum(0);
    h2_resVsEta_z0_90->SetMarkerStyle(20);
    h2_resVsEta_z0_90->Draw("p");
    c.SaveAs(DIR+type+"_resVsEta_z0_90.png");

    h2_resVsEta_z0_68->SetTitle(htitle);
    h2_resVsEta_z0_68->SetMinimum(0);
    h2_resVsEta_z0_68->SetMarkerStyle(20);
    h2_resVsEta_z0_68->Draw("p");
    c.SaveAs(DIR+type+"_resVsEta_z0_68.png");

    if (doDetailedPlots) {
      h2_resVsEta_z0_L_90->SetTitle(htitle);
      h2_resVsEta_z0_L_90->Draw();
      sprintf(ctxt,"p_{T} < 8 GeV");
      mySmallText(0.22,0.82,1,ctxt);
      c.SaveAs(DIR+type+"_resVsEta_z0_L_90.png");
      
      h2_resVsEta_z0_H_90->SetTitle(htitle);
      h2_resVsEta_z0_H_90->Draw();
      sprintf(ctxt,"p_{T} > 8 GeV");
      mySmallText(0.22,0.82,1,ctxt);
      c.SaveAs(DIR+type+"_resVsEta_z0_H_90.png");
    }

    /*
    h2_resVsEta_d0_90->SetTitle(htitle);
    h2_resVsEta_d0_90->Draw();
    c.SaveAs(DIR+type+"_resVsEta_d0_90.png");

    h2_resVsEta_d0_L_90->SetTitle(htitle);
    h2_resVsEta_d0_L_90->Draw();
    sprintf(ctxt,"p_{T} < 8 GeV");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs(DIR+type+"_resVsEta_d0_L_90.png");

    h2_resVsEta_d0_H_90->SetTitle(htitle);
    h2_resVsEta_d0_H_90->Draw();
    sprintf(ctxt,"p_{T} > 8 GeV");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs(DIR+type+"_resVsEta_d0_H_90.png");
    */

    h2_resVsEta_phi_90->SetTitle(htitle);
    h2_resVsEta_phi_90->SetMinimum(0);
    h2_resVsEta_phi_90->SetMarkerStyle(20);
    h2_resVsEta_phi_90->Draw("p");
    c.SaveAs(DIR+type+"_resVsEta_phi_90.png");

    h2_resVsEta_phi_68->SetTitle(htitle);
    h2_resVsEta_phi_68->SetMinimum(0);
    h2_resVsEta_phi_68->SetMarkerStyle(20);
    h2_resVsEta_phi_68->Draw("p");
    c.SaveAs(DIR+type+"_resVsEta_phi_68.png");

    if (doDetailedPlots) {
      h2_resVsEta_phi_L_90->SetTitle(htitle);
      h2_resVsEta_phi_L_90->Draw();
      sprintf(ctxt,"p_{T} < 8 GeV");
      mySmallText(0.22,0.82,1,ctxt);
      c.SaveAs(DIR+type+"_resVsEta_phi_L_90.png");
      
      h2_resVsEta_phi_H_90->SetTitle(htitle);
      h2_resVsEta_phi_H_90->Draw();
      sprintf(ctxt,"p_{T} > 8 GeV");
      mySmallText(0.22,0.82,1,ctxt);
      c.SaveAs(DIR+type+"_resVsEta_phi_H_90.png");
    }

    h2_resVsEta_ptRel_90->SetTitle(htitle);
    h2_resVsEta_ptRel_90->SetMinimum(0);
    h2_resVsEta_ptRel_90->SetMarkerStyle(20);
    h2_resVsEta_ptRel_90->Draw("p");
    c.SaveAs(DIR+type+"_resVsEta_ptRel_90.png");

    h2_resVsEta_ptRel_68->SetTitle(htitle);
    h2_resVsEta_ptRel_68->SetMinimum(0);
    h2_resVsEta_ptRel_68->SetMarkerStyle(20);
    h2_resVsEta_ptRel_68->Draw("p");
    c.SaveAs(DIR+type+"_resVsEta_ptRel_68.png");

    if (doDetailedPlots) {
      h2_resVsEta_ptRel_L_90->SetTitle(htitle);
      h2_resVsEta_ptRel_L_90->Draw();
      sprintf(ctxt,"p_{T} < 8 GeV");
      mySmallText(0.22,0.82,1,ctxt);
      c.SaveAs(DIR+type+"_resVsEta_ptRel_L_90.png");
      
      h2_resVsEta_ptRel_H_90->SetTitle(htitle);
      h2_resVsEta_ptRel_H_90->Draw();
      sprintf(ctxt,"p_{T} > 8 GeV");
      mySmallText(0.22,0.82,1,ctxt);
      c.SaveAs(DIR+type+"_resVsEta_ptRel_H_90.png");

      h2_resVsEta_eta  ->Write();
      h2_resVsEta_eta_L->Write();
      h2_resVsEta_eta_H->Write();

      h2_mresVsEta_eta  ->Write();
      h2_mresVsEta_eta_L->Write();
      h2_mresVsEta_eta_H->Write();

      h2_resVsEta_z0  ->Write();
      h2_resVsEta_z0_L->Write();
      h2_resVsEta_z0_H->Write();

      h2_resVsEta_d0  ->Write();
      h2_resVsEta_d0_L->Write();
      h2_resVsEta_d0_H->Write();

      h2_resVsEta_phi  ->Write();
      h2_resVsEta_phi_L->Write();
      h2_resVsEta_phi_H->Write();
      
      h2_resVsEta_ptRel  ->Write();
      h2_resVsEta_ptRel_L->Write();
      h2_resVsEta_ptRel_H->Write();
    }

    if (doGausFit) {
      h3_resVsEta_eta_L  ->Write();
      h3_resVsEta_z0_L   ->Write();
      h3_resVsEta_phi_L  ->Write();
      h3_resVsEta_ptRel_L->Write();
      
      h3_resVsEta_eta_H  ->Write();
      h3_resVsEta_z0_H   ->Write();
      h3_resVsEta_phi_H  ->Write();
      h3_resVsEta_ptRel_H->Write();
    }

    // resolution vs phi
    h2_resVsPhi_pt_90->SetTitle(htitle);
    h2_resVsPhi_pt_90->SetMinimum(0);
    h2_resVsPhi_pt_90->SetMarkerStyle(20);
    h2_resVsPhi_pt_90->Draw("p");
    c.SaveAs(DIR+type+"_resVsPhi_pt_90.png");

    h2_resVsPhi_ptRel_90->SetTitle(htitle);
    h2_resVsPhi_ptRel_90->SetMinimum(0);
    h2_resVsPhi_ptRel_90->SetMarkerStyle(20);
    h2_resVsPhi_ptRel_90->Draw("p");
    c.SaveAs(DIR+type+"_resVsPhi_ptRel_90.png");
  }
  
  

  // ---------------------------------------------------------------------------------------------------------
  // Track Quality Plots
  // ---------------------------------------------------------------------------------------------------------
  if (doDetailedPlots) {
    h_match_trk_nstub   ->Write();
    h_match_trk_nstub_C ->Write();
    h_match_trk_nstub_I ->Write();
    h_match_trk_nstub_F ->Write();
  }

  h_match_trk_chi2->SetTitle(htitle);
  h_match_trk_chi2->Draw();
  h_match_trk_chi2->Write();
    sprintf(ctxt,"|eta| < 2.4");
    mySmallText(0.52,0.82,1,ctxt);
  c.SaveAs(DIR+type+"_match_trk_chi2.png");

  h_match_trk_chi2_dof->SetTitle(htitle);
  h_match_trk_chi2_dof->Draw();
  h_match_trk_chi2_dof->Write();
    sprintf(ctxt,"|eta| < 2.4");
    mySmallText(0.52,0.82,1,ctxt);
  c.SaveAs(DIR+type+"_match_trk_chi2_dof.png");

  if (doDetailedPlots) {
    h_match_trk_chi2_C_L->Write();
    h_match_trk_chi2_I_L->Write();
    h_match_trk_chi2_F_L->Write();
    h_match_trk_chi2_C_H->Write();
    h_match_trk_chi2_I_H->Write();
    h_match_trk_chi2_F_H->Write();

    h_match_trk_chi2_dof_C_L->Write();
    h_match_trk_chi2_dof_I_L->Write();
    h_match_trk_chi2_dof_F_L->Write();
    h_match_trk_chi2_dof_C_H->Write();
    h_match_trk_chi2_dof_I_H->Write();
    h_match_trk_chi2_dof_F_H->Write();
  }



  // ---------------------------------------------------------------------------------------------------------
  // Efficiency Plots
  // ---------------------------------------------------------------------------------------------------------

  // ---------------------------------------------------------------------------------------------------------
  // Rebin Plots

  //-----Tracking Particle Plots-----//
  // General Efficiencies
  h_tp_pt         ->Rebin(2);
  h_tp_pt_L       ->Rebin(2);
  h_tp_pt_H       ->Rebin(2);
  h_match_tp_pt   ->Rebin(2);
  h_match_tp_pt_L ->Rebin(2);
  h_match_tp_pt_H ->Rebin(2);

  h_tp_phi       ->Rebin(2);
  h_match_tp_phi ->Rebin(2);

  // Jet Axis Efficiencies
  h_dpt_tp         ->Rebin(2);
  h_dpt_tp_L       ->Rebin(2);
  h_dpt_tp_H       ->Rebin(2);
  h_dpt_match_tp   ->Rebin(2);
  h_dpt_match_tp_L ->Rebin(2);
  h_dpt_match_tp_H ->Rebin(2);

  h_deta_tp         ->Rebin(2);
  h_deta_tp_L       ->Rebin(2);
  h_deta_tp_H       ->Rebin(2);
  h_deta_match_tp   ->Rebin(2);
  h_deta_match_tp_L ->Rebin(2);
  h_deta_match_tp_H ->Rebin(2);

  //h_eff_tp_jetpt
  h_dpt_jet           ->Rebin(2);
  h_dpt_match_jet     ->Rebin(2);

  h_dpt_jet_C         ->Rebin(5);
  h_dpt_jet_I         ->Rebin(5);
  h_dpt_jet_F         ->Rebin(5);
  h_dpt_match_jet_C   ->Rebin(5);
  h_dpt_match_jet_I   ->Rebin(5);
  h_dpt_match_jet_F   ->Rebin(5);

  h_dpt_jet_H         ->Rebin(3);
  h_dpt_jet_HC        ->Rebin(3);
  h_dpt_jet_HI        ->Rebin(3);
  h_dpt_jet_HF        ->Rebin(3);
  h_dpt_match_jet_H   ->Rebin(3);
  h_dpt_match_jet_HC  ->Rebin(3);
  h_dpt_match_jet_HI  ->Rebin(3);
  h_dpt_match_jet_HF  ->Rebin(3);

  h_dpt_jet_full         ->Rebin(5);
  h_dpt_jet_full_C       ->Rebin(5);
  h_dpt_jet_full_I       ->Rebin(5);
  h_dpt_jet_full_F       ->Rebin(5);
  h_dpt_match_jet_full   ->Rebin(5);
  h_dpt_match_jet_full_C ->Rebin(5);
  h_dpt_match_jet_full_I ->Rebin(5);
  h_dpt_match_jet_full_F ->Rebin(5);


  //-----L1 Track Plots-----//
  // General Purities
  h_trk_pt       ->Rebin(2);
  h_matchtrk_pt  ->Rebin(2);
  h_trk_phi      ->Rebin(2);
  h_matchtrk_phi ->Rebin(2);

  // Jet Axis Purities
  h_dpt_trk       ->Rebin(2);
  h_dpt_matchtrk  ->Rebin(2);

  //h_eff_L1_jetpt
  h_dpt_trk_jet           ->Rebin(2);
  h_dpt_trk_jet_full      ->Rebin(5);
  h_dpt_matchtrk_jet      ->Rebin(2);
  h_dpt_matchtrk_jet_full ->Rebin(5);


  // ---------------------------------------------------------------------------------------------------------
  // Print the pT, eta, phi Event Plots
  
  // ----- Tracking Particles ----- //
  h_tp_p->SetTitle(htitle);
  h_tp_p->Draw();
  h_tp_p->Write();
  c.SaveAs(DIR+type+"_ev_TP_p.png");

  h_tp_pt->SetTitle(htitle);
  h_tp_pt->Draw();
  h_tp_pt->Write();
  c.SaveAs(DIR+type+"_ev_TP_pt.png");

  h_tp_eta->SetTitle(htitle);
  h_tp_eta->Draw();
  h_tp_eta->Write();
  c.SaveAs(DIR+type+"_ev_TP_eta.png");

  h_tp_phi->SetTitle(htitle);
  h_tp_phi->Draw();
  h_tp_phi->Write();
  c.SaveAs(DIR+type+"_ev_TP_phi.png");

  // ----- Matched Tracking Particles ----- //
  h_match_tp_p->SetTitle(htitle);
  h_match_tp_p->Draw();
  h_match_tp_p->Write();
  c.SaveAs(DIR+type+"_ev_TP_match_p.png");

  h_match_tp_pt->SetTitle(htitle);
  h_match_tp_pt->Draw();
  h_match_tp_pt->Write();
  c.SaveAs(DIR+type+"_ev_TP_match_pt.png");

  h_match_tp_eta->SetTitle(htitle);
  h_match_tp_eta->Draw();
  h_match_tp_eta->Write();
  c.SaveAs(DIR+type+"_ev_TP_match_eta.png");

  h_match_tp_phi->SetTitle(htitle);
  h_match_tp_phi->Draw();
  h_match_tp_phi->Write();
  c.SaveAs(DIR+type+"_ev_TP_match_phi.png");

  // ----- L1 Tracks ----- //
  h_trk_p->SetTitle(htitle);
  h_trk_p->Draw();
  h_trk_p->Write();
  c.SaveAs(DIR+type+"_ev_L1_p.png");

  h_trk_pt->SetTitle(htitle);
  h_trk_pt->Draw();
  h_trk_pt->Write();
  c.SaveAs(DIR+type+"_ev_L1_pt.png");

  h_trk_eta->SetTitle(htitle);
  h_trk_eta->Draw();
  h_trk_eta->Write();
  c.SaveAs(DIR+type+"_ev_L1_eta.png");

  h_trk_phi->SetTitle(htitle);
  h_trk_phi->Draw();
  h_trk_phi->Write();
  c.SaveAs(DIR+type+"_ev_L1_phi.png");

  // ----- Matched L1 Tracks ----- //
  h_matchtrk_p->SetTitle(htitle);
  h_matchtrk_p->Draw();
  h_matchtrk_p->Write();
  c.SaveAs(DIR+type+"_ev_L1_match_p.png");

  h_matchtrk_pt->SetTitle(htitle);
  h_matchtrk_pt->Draw();
  h_matchtrk_pt->Write();
  c.SaveAs(DIR+type+"_ev_L1_match_pt.png");

  h_matchtrk_eta->SetTitle(htitle);
  h_matchtrk_eta->Draw();
  h_matchtrk_eta->Write();
  c.SaveAs(DIR+type+"_ev_L1_match_eta.png");

  h_matchtrk_phi->SetTitle(htitle);
  h_matchtrk_phi->Draw();
  h_matchtrk_phi->Write();
  c.SaveAs(DIR+type+"_ev_L1_match_phi.png");

  // ----- Jets ----- //
  h_jet_p->SetTitle(htitle);
  h_jet_p->Draw();
  h_jet_p->Write();
  c.SaveAs(DIR+type+"_ev_jet_p.png");

  h_jet_pt->SetTitle(htitle);
  h_jet_pt->Draw();
  h_jet_pt->Write();
  c.SaveAs(DIR+type+"_ev_jet_pt.png");

  h_jet_eta->SetTitle(htitle);
  h_jet_eta->Draw();
  h_jet_eta->Write();
  c.SaveAs(DIR+type+"_ev_jet_eta.png");

  h_jet_phi->SetTitle(htitle);
  h_jet_phi->Draw();
  h_jet_phi->Write();
  c.SaveAs(DIR+type+"_ev_jet_phi.png");



  // ---------------------------------------------------------------------------------------------------------
  // Efficiency Calculations

  //-----Efficiency - Tracking Particles-----//
  h_match_tp_pt->Sumw2();
  h_tp_pt->Sumw2();
  TH1F* h_eff_TP_pt = (TH1F*) h_match_tp_pt->Clone();
  h_eff_TP_pt->SetName("eff_TP_pt");
  h_eff_TP_pt->GetYaxis()->SetTitle("Efficiency");
  h_eff_TP_pt->Divide(h_match_tp_pt, h_tp_pt, 1.0, 1.0, "B");

    h_match_tp_pt_L->Sumw2();
    h_tp_pt_L->Sumw2();
    TH1F* h_eff_TP_pt_L = (TH1F*) h_match_tp_pt_L->Clone();
    h_eff_TP_pt_L->SetName("eff_TP_pt_L");
    h_eff_TP_pt_L->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_pt_L->Divide(h_match_tp_pt_L, h_tp_pt_L, 1.0, 1.0, "B");

    h_match_tp_pt_LC->Sumw2();
    h_tp_pt_LC->Sumw2();
    TH1F* h_eff_TP_pt_LC = (TH1F*) h_match_tp_pt_LC->Clone();
    h_eff_TP_pt_LC->SetName("eff_TP_pt_LC");
    h_eff_TP_pt_LC->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_pt_LC->Divide(h_match_tp_pt_LC, h_tp_pt_LC, 1.0, 1.0, "B");

    h_match_tp_pt_H->Sumw2();
    h_tp_pt_H->Sumw2();
    TH1F* h_eff_TP_pt_H = (TH1F*) h_match_tp_pt_H->Clone();
    h_eff_TP_pt_H->SetName("eff_TP_pt_H");
    h_eff_TP_pt_H->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_pt_H->Divide(h_match_tp_pt_H, h_tp_pt_H, 1.0, 1.0, "B");

  h_match_tp_eta->Sumw2();
  h_tp_eta->Sumw2();
  TH1F* h_eff_TP_eta = (TH1F*) h_match_tp_eta->Clone();
  h_eff_TP_eta->SetName("eff_TP_eta");
  h_eff_TP_eta->GetYaxis()->SetTitle("Efficiency");
  h_eff_TP_eta->Divide(h_match_tp_eta, h_tp_eta, 1.0, 1.0, "B");

    h_match_tp_eta_L->Sumw2();
    h_tp_eta_L->Sumw2();
    TH1F* h_eff_TP_eta_L = (TH1F*) h_match_tp_eta_L->Clone();
    h_eff_TP_eta_L->SetName("eff_TP_eta_L");
    h_eff_TP_eta_L->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_eta_L->Divide(h_match_tp_eta_L, h_tp_eta_L, 1.0, 1.0, "B");

    h_match_tp_eta_H->Sumw2();
    h_tp_eta_H->Sumw2();
    TH1F* h_eff_TP_eta_H = (TH1F*) h_match_tp_eta_H->Clone();
    h_eff_TP_eta_H->SetName("eff_TP_eta_H");
    h_eff_TP_eta_H->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_eta_H->Divide(h_match_tp_eta_H, h_tp_eta_H, 1.0, 1.0, "B");

  h_match_tp_phi->Sumw2();
  h_tp_phi->Sumw2();
  TH1F* h_eff_TP_phi = (TH1F*) h_match_tp_phi->Clone();
  h_eff_TP_phi->SetName("eff_TP_phi");
  h_eff_TP_phi->GetYaxis()->SetTitle("Efficiency");
  h_eff_TP_phi->Divide(h_match_tp_phi, h_tp_phi, 1.0, 1.0, "B");

  h_match_tp_z0->Sumw2();
  h_tp_z0->Sumw2();
  TH1F* h_eff_TP_z0 = (TH1F*) h_match_tp_z0->Clone();
  h_eff_TP_z0->SetName("eff_TP_z0");
  h_eff_TP_z0->GetYaxis()->SetTitle("Efficiency");
  h_eff_TP_z0->Divide(h_match_tp_z0, h_tp_z0, 1.0, 1.0, "B");

    h_match_tp_z0_L->Sumw2();
    h_tp_z0_L->Sumw2();
    TH1F* h_eff_TP_z0_L = (TH1F*) h_match_tp_z0_L->Clone();
    h_eff_TP_z0_L->SetName("eff_TP_z0_L");
    h_eff_TP_z0_L->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_z0_L->Divide(h_match_tp_z0_L, h_tp_z0_L, 1.0, 1.0, "B");

    h_match_tp_z0_H->Sumw2();
    h_tp_z0_H->Sumw2();
    TH1F* h_eff_TP_z0_H = (TH1F*) h_match_tp_z0_H->Clone();
    h_eff_TP_z0_H->SetName("eff_TP_z0_H");
    h_eff_TP_z0_H->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_z0_H->Divide(h_match_tp_z0_H, h_tp_z0_H, 1.0, 1.0, "B");

  h_match_tp_d0->Sumw2();
  h_tp_d0->Sumw2();
  TH1F* h_eff_TP_d0 = (TH1F*) h_match_tp_d0->Clone();
  h_eff_TP_d0->SetName("eff_TP_d0");
  h_eff_TP_d0->GetYaxis()->SetTitle("Efficiency");
  h_eff_TP_d0->Divide(h_match_tp_d0, h_tp_d0, 1.0, 1.0, "B");

  h_match_tp_absd0->Sumw2();
  h_tp_absd0->Sumw2();
  TH1F* h_eff_TP_absd0 = (TH1F*) h_match_tp_absd0->Clone();
  h_eff_TP_absd0->SetName("eff_TP_absd0");
  h_eff_TP_absd0->GetYaxis()->SetTitle("Efficiency");
  h_eff_TP_absd0->Divide(h_match_tp_absd0, h_tp_absd0, 1.0, 1.0, "B");


  //-----Efficiency - Tracking Particles - Jet Axis-----//

  //h_eff_TP_jetaxis_eta
  h_deta_match_tp->Sumw2();
  h_deta_tp->Sumw2();
  TH1F* h_eff_TP_jetaxis_eta = (TH1F*) h_deta_match_tp->Clone();
  h_eff_TP_jetaxis_eta->SetName("eff_TP_jetaxis_eta");
  h_eff_TP_jetaxis_eta->GetYaxis()->SetTitle("Efficiency");
  h_eff_TP_jetaxis_eta->Divide(h_deta_match_tp, h_deta_tp, 1.0, 1.0, "B");

    h_deta_match_tp_L->Sumw2();
    h_deta_tp_L->Sumw2();
    TH1F* h_eff_TP_jetaxis_eta_L = (TH1F*) h_deta_match_tp_L->Clone();
    h_eff_TP_jetaxis_eta_L->SetName("eff_TP_jetaxis_eta_L");
    h_eff_TP_jetaxis_eta_L->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_jetaxis_eta_L->Divide(h_deta_match_tp_L, h_deta_tp_L, 1.0, 1.0, "B");

    h_deta_match_tp_H->Sumw2();
    h_deta_tp_H->Sumw2();
    TH1F* h_eff_TP_jetaxis_eta_H = (TH1F*) h_deta_match_tp_H->Clone();
    h_eff_TP_jetaxis_eta_H->SetName("eff_TP_jetaxis_eta_H");
    h_eff_TP_jetaxis_eta_H->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_jetaxis_eta_H->Divide(h_deta_match_tp_H, h_deta_tp_H, 1.0, 1.0, "B");

  //h_eff_TP_jetaxis_phi
  h_dphi_match_tp->Sumw2();
  h_dphi_tp->Sumw2();
  TH1F* h_eff_TP_jetaxis_phi = (TH1F*) h_dphi_match_tp->Clone();
  h_eff_TP_jetaxis_phi->SetName("eff_TP_jetaxis_phi");
  h_eff_TP_jetaxis_phi->GetYaxis()->SetTitle("Efficiency");
  h_eff_TP_jetaxis_phi->Divide(h_dphi_match_tp, h_dphi_tp, 1.0, 1.0, "B");

    h_dphi_match_tp_L->Sumw2();
    h_dphi_tp_L->Sumw2();
    TH1F* h_eff_TP_jetaxis_phi_L = (TH1F*) h_dphi_match_tp_L->Clone();
    h_eff_TP_jetaxis_phi_L->SetName("eff_TP_jetaxis_phi_L");
    h_eff_TP_jetaxis_phi_L->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_jetaxis_phi_L->Divide(h_dphi_match_tp_L, h_dphi_tp_L, 1.0, 1.0, "B");

    h_dphi_match_tp_H->Sumw2();
    h_dphi_tp_H->Sumw2();
    TH1F* h_eff_TP_jetaxis_phi_H = (TH1F*) h_dphi_match_tp_H->Clone();
    h_eff_TP_jetaxis_phi_H->SetName("eff_TP_jetaxis_phi_H");
    h_eff_TP_jetaxis_phi_H->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_jetaxis_phi_H->Divide(h_dphi_match_tp_H, h_dphi_tp_H, 1.0, 1.0, "B");

  //h_eff_TP_jetaxis_pt
  h_dpt_match_tp->Sumw2();
  h_dpt_tp->Sumw2();
  TH1F* h_eff_TP_jetaxis_pt = (TH1F*) h_dpt_match_tp->Clone();
  h_eff_TP_jetaxis_pt->SetName("eff_TP_jetaxis_pt");
  h_eff_TP_jetaxis_pt->GetYaxis()->SetTitle("Efficiency");
  h_eff_TP_jetaxis_pt->Divide(h_dpt_match_tp, h_dpt_tp, 1.0, 1.0, "B");

    h_dpt_match_tp_L->Sumw2();
    h_dpt_tp_L->Sumw2();
    TH1F* h_eff_TP_jetaxis_pt_L = (TH1F*) h_dpt_match_tp_L->Clone();
    h_eff_TP_jetaxis_pt_L->SetName("eff_TP_jetaxis_pt_L");
    h_eff_TP_jetaxis_pt_L->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_jetaxis_pt_L->Divide(h_dpt_match_tp_L, h_dpt_tp_L, 1.0, 1.0, "B");

    h_dpt_match_tp_H->Sumw2();
    h_dpt_tp_H->Sumw2();
    TH1F* h_eff_TP_jetaxis_pt_H = (TH1F*) h_dpt_match_tp_H->Clone();
    h_eff_TP_jetaxis_pt_H->SetName("eff_TP_jetaxis_pt_H");
    h_eff_TP_jetaxis_pt_H->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_jetaxis_pt_H->Divide(h_dpt_match_tp_H, h_dpt_tp_H, 1.0, 1.0, "B");

  //h_eff_TP_jetaxis_dR
  h_dR_match_tp->Sumw2();
  h_dR_tp->Sumw2();
  TH1F* h_eff_TP_jetaxis_dR = (TH1F*) h_dR_match_tp->Clone();
  h_eff_TP_jetaxis_dR->SetName("eff_TP_jetaxis_dR");
  h_eff_TP_jetaxis_dR->GetYaxis()->SetTitle("Efficiency");
  h_eff_TP_jetaxis_dR->Divide(h_dR_match_tp, h_dR_tp, 1.0, 1.0, "B");

    h_dR_match_tp_L->Sumw2();
    h_dR_tp_L->Sumw2();
    TH1F* h_eff_TP_jetaxis_dR_L = (TH1F*) h_dR_match_tp_L->Clone();
    h_eff_TP_jetaxis_dR_L->SetName("eff_TP_jetaxis_dR_L");
    h_eff_TP_jetaxis_dR_L->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_jetaxis_dR_L->Divide(h_dR_match_tp_L, h_dR_tp_L, 1.0, 1.0, "B");

    h_dR_match_tp_H->Sumw2();
    h_dR_tp_H->Sumw2();
    TH1F* h_eff_TP_jetaxis_dR_H = (TH1F*) h_dR_match_tp_H->Clone();
    h_eff_TP_jetaxis_dR_H->SetName("eff_TP_jetaxis_dR_H");
    h_eff_TP_jetaxis_dR_H->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_jetaxis_dR_H->Divide(h_dR_match_tp_H, h_dR_tp_H, 1.0, 1.0, "B");

  //h_eff_TP_jetpt
  h_dpt_match_jet->Sumw2();
  h_dpt_jet->Sumw2();
  TH1F* h_eff_TP_jetpt = (TH1F*) h_dpt_match_jet->Clone();
  h_eff_TP_jetpt->SetName("eff_TP_jetpt");
  h_eff_TP_jetpt->GetYaxis()->SetTitle("Efficiency");
  h_eff_TP_jetpt->Divide(h_dpt_match_jet, h_dpt_jet, 1.0, 1.0, "B");

      h_dpt_match_jet_C->Sumw2();
      h_dpt_jet_C->Sumw2();
      TH1F* h_eff_TP_jetpt_C = (TH1F*) h_dpt_match_jet_C->Clone();
      h_eff_TP_jetpt_C->SetName("eff_TP_jetpt_C");
      h_eff_TP_jetpt_C->GetYaxis()->SetTitle("Efficiency");
      h_eff_TP_jetpt_C->Divide(h_dpt_match_jet_C, h_dpt_jet_C, 1.0, 1.0, "B");

      h_dpt_match_jet_I->Sumw2();
      h_dpt_jet_I->Sumw2();
      TH1F* h_eff_TP_jetpt_I = (TH1F*) h_dpt_match_jet_I->Clone();
      h_eff_TP_jetpt_I->SetName("eff_TP_jetpt_I");
      h_eff_TP_jetpt_I->GetYaxis()->SetTitle("Efficiency");
      h_eff_TP_jetpt_I->Divide(h_dpt_match_jet_I, h_dpt_jet_I, 1.0, 1.0, "B");

      h_dpt_match_jet_F->Sumw2();
      h_dpt_jet_F->Sumw2();
      TH1F* h_eff_TP_jetpt_F = (TH1F*) h_dpt_match_jet_F->Clone();
      h_eff_TP_jetpt_F->SetName("eff_TP_jetpt_F");
      h_eff_TP_jetpt_F->GetYaxis()->SetTitle("Efficiency");
      h_eff_TP_jetpt_F->Divide(h_dpt_match_jet_F, h_dpt_jet_F, 1.0, 1.0, "B");

    h_dpt_match_jet_L->Sumw2();
    h_dpt_jet_L->Sumw2();
    TH1F* h_eff_TP_jetpt_L = (TH1F*) h_dpt_match_jet_L->Clone();
    h_eff_TP_jetpt_L->SetName("eff_TP_jetpt_L");
    h_eff_TP_jetpt_L->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_jetpt_L->Divide(h_dpt_match_jet_L, h_dpt_jet_L, 1.0, 1.0, "B");

      h_dpt_match_jet_LC->Sumw2();
      h_dpt_jet_LC->Sumw2();
      TH1F* h_eff_TP_jetpt_LC = (TH1F*) h_dpt_match_jet_LC->Clone();
      h_eff_TP_jetpt_LC->SetName("eff_TP_jetpt_LC");
      h_eff_TP_jetpt_LC->GetYaxis()->SetTitle("Efficiency");
      h_eff_TP_jetpt_LC->Divide(h_dpt_match_jet_LC, h_dpt_jet_LC, 1.0, 1.0, "B");

      h_dpt_match_jet_LI->Sumw2();
      h_dpt_jet_LI->Sumw2();
      TH1F* h_eff_TP_jetpt_LI = (TH1F*) h_dpt_match_jet_LI->Clone();
      h_eff_TP_jetpt_LI->SetName("eff_TP_jetpt_LI");
      h_eff_TP_jetpt_LI->GetYaxis()->SetTitle("Efficiency");
      h_eff_TP_jetpt_LI->Divide(h_dpt_match_jet_LI, h_dpt_jet_LI, 1.0, 1.0, "B");

      h_dpt_match_jet_LF->Sumw2();
      h_dpt_jet_LF->Sumw2();
      TH1F* h_eff_TP_jetpt_LF = (TH1F*) h_dpt_match_jet_LF->Clone();
      h_eff_TP_jetpt_LF->SetName("eff_TP_jetpt_LF");
      h_eff_TP_jetpt_LF->GetYaxis()->SetTitle("Efficiency");
      h_eff_TP_jetpt_LF->Divide(h_dpt_match_jet_LF, h_dpt_jet_LF, 1.0, 1.0, "B");

    h_dpt_match_jet_H->Sumw2();
    h_dpt_jet_H->Sumw2();
    TH1F* h_eff_TP_jetpt_H = (TH1F*) h_dpt_match_jet_H->Clone();
    h_eff_TP_jetpt_H->SetName("eff_TP_jetpt_H");
    h_eff_TP_jetpt_H->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_jetpt_H->Divide(h_dpt_match_jet_H, h_dpt_jet_H, 1.0, 1.0, "B");

      h_dpt_match_jet_HC->Sumw2();
      h_dpt_jet_HC->Sumw2();
      TH1F* h_eff_TP_jetpt_HC = (TH1F*) h_dpt_match_jet_HC->Clone();
      h_eff_TP_jetpt_HC->SetName("eff_TP_jetpt_HC");
      h_eff_TP_jetpt_HC->GetYaxis()->SetTitle("Efficiency");
      h_eff_TP_jetpt_HC->Divide(h_dpt_match_jet_HC, h_dpt_jet_HC, 1.0, 1.0, "B");

      h_dpt_match_jet_HI->Sumw2();
      h_dpt_jet_HI->Sumw2();
      TH1F* h_eff_TP_jetpt_HI = (TH1F*) h_dpt_match_jet_HI->Clone();
      h_eff_TP_jetpt_HI->SetName("eff_TP_jetpt_HI");
      h_eff_TP_jetpt_HI->GetYaxis()->SetTitle("Efficiency");
      h_eff_TP_jetpt_HI->Divide(h_dpt_match_jet_HI, h_dpt_jet_HI, 1.0, 1.0, "B");

      h_dpt_match_jet_HF->Sumw2();
      h_dpt_jet_HF->Sumw2();
      TH1F* h_eff_TP_jetpt_HF = (TH1F*) h_dpt_match_jet_HF->Clone();
      h_eff_TP_jetpt_HF->SetName("eff_TP_jetpt_HF");
      h_eff_TP_jetpt_HF->GetYaxis()->SetTitle("Efficiency");
      h_eff_TP_jetpt_HF->Divide(h_dpt_match_jet_HF, h_dpt_jet_HF, 1.0, 1.0, "B");

    h_dpt_match_jet_full->Sumw2();
    h_dpt_jet_full->Sumw2();
    TH1F* h_eff_TP_jetpt_full = (TH1F*) h_dpt_match_jet_full->Clone();
    h_eff_TP_jetpt_full->SetName("eff_TP_jetpt_full");
    h_eff_TP_jetpt_full->GetYaxis()->SetTitle("Efficiency");
    h_eff_TP_jetpt_full->Divide(h_dpt_match_jet_full, h_dpt_jet_full, 1.0, 1.0, "B");

      h_dpt_match_jet_full_C->Sumw2();
      h_dpt_jet_full_C->Sumw2();
      TH1F* h_eff_TP_jetpt_full_C = (TH1F*) h_dpt_match_jet_full_C->Clone();
      h_eff_TP_jetpt_full_C->SetName("eff_TP_jetpt_full_C");
      h_eff_TP_jetpt_full_C->GetYaxis()->SetTitle("Efficiency");
      h_eff_TP_jetpt_full_C->Divide(h_dpt_match_jet_full_C, h_dpt_jet_full_C, 1.0, 1.0, "B");

      h_dpt_match_jet_full_I->Sumw2();
      h_dpt_jet_full_I->Sumw2();
      TH1F* h_eff_TP_jetpt_full_I = (TH1F*) h_dpt_match_jet_full_I->Clone();
      h_eff_TP_jetpt_full_I->SetName("eff_TP_jetpt_full_I");
      h_eff_TP_jetpt_full_I->GetYaxis()->SetTitle("Efficiency");
      h_eff_TP_jetpt_full_I->Divide(h_dpt_match_jet_full_I, h_dpt_jet_full_I, 1.0, 1.0, "B");

      h_dpt_match_jet_full_F->Sumw2();
      h_dpt_jet_full_F->Sumw2();
      TH1F* h_eff_TP_jetpt_full_F = (TH1F*) h_dpt_match_jet_full_F->Clone();
      h_eff_TP_jetpt_full_F->SetName("eff_TP_jetpt_full_F");
      h_eff_TP_jetpt_full_F->GetYaxis()->SetTitle("Efficiency");
      h_eff_TP_jetpt_full_F->Divide(h_dpt_match_jet_full_F, h_dpt_jet_full_F, 1.0, 1.0, "B");


  //-----Purity - L1 Tracks-----//

  //h_eff_L1_pt
  h_matchtrk_pt->Sumw2();
  h_trk_pt->Sumw2();
  TH1F* h_eff_L1_pt = (TH1F*) h_matchtrk_pt->Clone();
  h_eff_L1_pt->SetName("eff_L1_pt");
  h_eff_L1_pt->GetYaxis()->SetTitle("Purity");
  h_eff_L1_pt->Divide(h_matchtrk_pt, h_trk_pt, 1.0, 1.0, "B");

    h_matchtrk_pt_L->Sumw2();
    h_trk_pt_L->Sumw2();
    TH1F* h_eff_L1_pt_L = (TH1F*) h_matchtrk_pt_L->Clone();
    h_eff_L1_pt_L->SetName("eff_L1_pt_L");
    h_eff_L1_pt_L->GetYaxis()->SetTitle("Purity");
    h_eff_L1_pt_L->Divide(h_matchtrk_pt_L, h_trk_pt_L, 1.0, 1.0, "B");

    h_matchtrk_pt_LC->Sumw2();
    h_trk_pt_LC->Sumw2();
    TH1F* h_eff_L1_pt_LC = (TH1F*) h_matchtrk_pt_LC->Clone();
    h_eff_L1_pt_LC->SetName("eff_L1_pt_LC");
    h_eff_L1_pt_LC->GetYaxis()->SetTitle("Purity");
    h_eff_L1_pt_LC->Divide(h_matchtrk_pt_LC, h_trk_pt_LC, 1.0, 1.0, "B");

    h_matchtrk_pt_H->Sumw2();
    h_trk_pt_H->Sumw2();
    TH1F* h_eff_L1_pt_H = (TH1F*) h_matchtrk_pt_H->Clone();
    h_eff_L1_pt_H->SetName("eff_L1_pt_H");
    h_eff_L1_pt_H->GetYaxis()->SetTitle("Purity");
    h_eff_L1_pt_H->Divide(h_matchtrk_pt_H, h_trk_pt_H, 1.0, 1.0, "B");

  //h_eff_L1_eta
  h_matchtrk_eta->Sumw2();
  h_trk_eta->Sumw2();
  TH1F* h_eff_L1_eta = (TH1F*) h_matchtrk_eta->Clone();
  h_eff_L1_eta->SetName("eff_L1_eta");
  h_eff_L1_eta->GetYaxis()->SetTitle("Purity");
  h_eff_L1_eta->Divide(h_matchtrk_eta, h_trk_eta, 1.0, 1.0, "B");

    h_matchtrk_eta_L->Sumw2();
    h_trk_eta_L->Sumw2();
    TH1F* h_eff_L1_eta_L = (TH1F*) h_matchtrk_eta_L->Clone();
    h_eff_L1_eta_L->SetName("eff_L1_eta_L");
    h_eff_L1_eta_L->GetYaxis()->SetTitle("Purity");
    h_eff_L1_eta_L->Divide(h_matchtrk_eta_L, h_trk_eta_L, 1.0, 1.0, "B");

    h_matchtrk_eta_H->Sumw2();
    h_trk_eta_H->Sumw2();
    TH1F* h_eff_L1_eta_H = (TH1F*) h_matchtrk_eta_H->Clone();
    h_eff_L1_eta_H->SetName("eff_L1_eta_H");
    h_eff_L1_eta_H->GetYaxis()->SetTitle("Purity");
    h_eff_L1_eta_H->Divide(h_matchtrk_eta_H, h_trk_eta_H, 1.0, 1.0, "B");

  //h_eff_L1_phi
  h_matchtrk_phi->Sumw2();
  h_trk_phi->Sumw2();
  TH1F* h_eff_L1_phi = (TH1F*) h_matchtrk_phi->Clone();
  h_eff_L1_phi->SetName("eff_L1_phi");
  h_eff_L1_phi->GetYaxis()->SetTitle("Purity");
  h_eff_L1_phi->Divide(h_matchtrk_phi, h_trk_phi, 1.0, 1.0, "B");


  //-----Purity - L1 Tracks around Jet Axis-----//

  //h_eff_L1_jetaxis_eta
  h_deta_matchtrk->Sumw2();
  h_deta_trk->Sumw2();
  TH1F* h_eff_L1_jetaxis_eta = (TH1F*) h_deta_matchtrk->Clone();
  h_eff_L1_jetaxis_eta->SetName("eff_L1_jetaxis_eta");
  h_eff_L1_jetaxis_eta->GetYaxis()->SetTitle("Purity");
  h_eff_L1_jetaxis_eta->Divide(h_deta_matchtrk, h_deta_trk, 1.0, 1.0, "B");

    h_deta_matchtrk_L->Sumw2();
    h_deta_trk_L->Sumw2();
    TH1F* h_eff_L1_jetaxis_eta_L = (TH1F*) h_deta_matchtrk_L->Clone();
    h_eff_L1_jetaxis_eta_L->SetName("eff_L1_jetaxis_eta_L");
    h_eff_L1_jetaxis_eta_L->GetYaxis()->SetTitle("Purity");
    h_eff_L1_jetaxis_eta_L->Divide(h_deta_matchtrk_L, h_deta_trk_L, 1.0, 1.0, "B");

    h_deta_matchtrk_H->Sumw2();
    h_deta_trk_H->Sumw2();
    TH1F* h_eff_L1_jetaxis_eta_H = (TH1F*) h_deta_matchtrk_H->Clone();
    h_eff_L1_jetaxis_eta_H->SetName("eff_L1_jetaxis_eta_H");
    h_eff_L1_jetaxis_eta_H->GetYaxis()->SetTitle("Purity");
    h_eff_L1_jetaxis_eta_H->Divide(h_deta_matchtrk_H, h_deta_trk_H, 1.0, 1.0, "B");

  //h_eff_L1_jetaxis_phi
  h_dphi_matchtrk->Sumw2();
  h_dphi_trk->Sumw2();
  TH1F* h_eff_L1_jetaxis_phi = (TH1F*) h_dphi_matchtrk->Clone();
  h_eff_L1_jetaxis_phi->SetName("eff_L1_jetaxis_phi");
  h_eff_L1_jetaxis_phi->GetYaxis()->SetTitle("Purity");
  h_eff_L1_jetaxis_phi->Divide(h_dphi_matchtrk, h_dphi_trk, 1.0, 1.0, "B");

    h_dphi_matchtrk_L->Sumw2();
    h_dphi_trk_L->Sumw2();
    TH1F* h_eff_L1_jetaxis_phi_L = (TH1F*) h_dphi_matchtrk_L->Clone();
    h_eff_L1_jetaxis_phi_L->SetName("eff_L1_jetaxis_phi_L");
    h_eff_L1_jetaxis_phi_L->GetYaxis()->SetTitle("Purity");
    h_eff_L1_jetaxis_phi_L->Divide(h_dphi_matchtrk_L, h_dphi_trk_L, 1.0, 1.0, "B");

    h_dphi_matchtrk_H->Sumw2();
    h_dphi_trk_H->Sumw2();
    TH1F* h_eff_L1_jetaxis_phi_H = (TH1F*) h_dphi_matchtrk_H->Clone();
    h_eff_L1_jetaxis_phi_H->SetName("eff_L1_jetaxis_phi_H");
    h_eff_L1_jetaxis_phi_H->GetYaxis()->SetTitle("Purity");
    h_eff_L1_jetaxis_phi_H->Divide(h_dphi_matchtrk_H, h_dphi_trk_H, 1.0, 1.0, "B");

  //h_eff_L1_jetaxis_pt
  h_dpt_matchtrk->Sumw2();
  h_dpt_trk->Sumw2();
  TH1F* h_eff_L1_jetaxis_pt = (TH1F*) h_dpt_matchtrk->Clone();
  h_eff_L1_jetaxis_pt->SetName("eff_L1_jetaxis_pt");
  h_eff_L1_jetaxis_pt->GetYaxis()->SetTitle("Purity");
  h_eff_L1_jetaxis_pt->Divide(h_dpt_matchtrk, h_dpt_trk, 1.0, 1.0, "B");

    h_dpt_matchtrk_L->Sumw2();
    h_dpt_trk_L->Sumw2();
    TH1F* h_eff_L1_jetaxis_pt_L = (TH1F*) h_dpt_matchtrk_L->Clone();
    h_eff_L1_jetaxis_pt_L->SetName("eff_L1_jetaxis_pt_L");
    h_eff_L1_jetaxis_pt_L->GetYaxis()->SetTitle("Purity");
    h_eff_L1_jetaxis_pt_L->Divide(h_dpt_matchtrk_L, h_dpt_trk_L, 1.0, 1.0, "B");

    h_dpt_matchtrk_H->Sumw2();
    h_dpt_trk_H->Sumw2();
    TH1F* h_eff_L1_jetaxis_pt_H = (TH1F*) h_dpt_matchtrk_H->Clone();
    h_eff_L1_jetaxis_pt_H->SetName("eff_L1_jetaxis_pt_H");
    h_eff_L1_jetaxis_pt_H->GetYaxis()->SetTitle("Purity");
    h_eff_L1_jetaxis_pt_H->Divide(h_dpt_matchtrk_H, h_dpt_trk_H, 1.0, 1.0, "B");

  //h_eff_L1_jetaxis_dR
  h_dR_matchtrk->Sumw2();
  h_dR_trk->Sumw2();
  TH1F* h_eff_L1_jetaxis_dR = (TH1F*) h_dR_matchtrk->Clone();
  h_eff_L1_jetaxis_dR->SetName("eff_L1_jetaxis_dR");
  h_eff_L1_jetaxis_dR->GetYaxis()->SetTitle("Purity");
  h_eff_L1_jetaxis_dR->Divide(h_dR_matchtrk, h_dR_trk, 1.0, 1.0, "B");

    h_dR_matchtrk_L->Sumw2();
    h_dR_trk_L->Sumw2();
    TH1F* h_eff_L1_jetaxis_dR_L = (TH1F*) h_dR_matchtrk_L->Clone();
    h_eff_L1_jetaxis_dR_L->SetName("eff_L1_jetaxis_dR_L");
    h_eff_L1_jetaxis_dR_L->GetYaxis()->SetTitle("Purity");
    h_eff_L1_jetaxis_dR_L->Divide(h_dR_matchtrk_L, h_dR_trk_L, 1.0, 1.0, "B");

    h_dR_matchtrk_H->Sumw2();
    h_dR_trk_H->Sumw2();
    TH1F* h_eff_L1_jetaxis_dR_H = (TH1F*) h_dR_matchtrk_H->Clone();
    h_eff_L1_jetaxis_dR_H->SetName("eff_L1_jetaxis_dR_H");
    h_eff_L1_jetaxis_dR_H->GetYaxis()->SetTitle("Purity");
    h_eff_L1_jetaxis_dR_H->Divide(h_dR_matchtrk_H, h_dR_trk_H, 1.0, 1.0, "B");

  //h_eff_L1_jetpt
  h_dpt_matchtrk_jet->Sumw2();
  h_dpt_trk_jet->Sumw2();
  TH1F* h_eff_L1_jetpt = (TH1F*) h_dpt_matchtrk_jet->Clone();
  h_eff_L1_jetpt->SetName("eff_L1_jetpt");
  h_eff_L1_jetpt->GetYaxis()->SetTitle("Purity");
  h_eff_L1_jetpt->Divide(h_dpt_matchtrk_jet, h_dpt_trk_jet, 1.0, 1.0, "B");

    h_dpt_matchtrk_jet_L->Sumw2();
    h_dpt_trk_jet_L->Sumw2();
    TH1F* h_eff_L1_jetpt_L = (TH1F*) h_dpt_matchtrk_jet_L->Clone();
    h_eff_L1_jetpt_L->SetName("eff_L1_jetpt_L");
    h_eff_L1_jetpt_L->GetYaxis()->SetTitle("Purity");
    h_eff_L1_jetpt_L->Divide(h_dpt_matchtrk_jet_L, h_dpt_trk_jet_L, 1.0, 1.0, "B");

    h_dpt_matchtrk_jet_H->Sumw2();
    h_dpt_trk_jet_H->Sumw2();
    TH1F* h_eff_L1_jetpt_H = (TH1F*) h_dpt_matchtrk_jet_H->Clone();
    h_eff_L1_jetpt_H->SetName("eff_L1_jetpt_H");
    h_eff_L1_jetpt_H->GetYaxis()->SetTitle("Purity");
    h_eff_L1_jetpt_H->Divide(h_dpt_matchtrk_jet_H, h_dpt_trk_jet_H, 1.0, 1.0, "B");

    h_dpt_matchtrk_jet_full->Sumw2();
    h_dpt_trk_jet_full->Sumw2();
    TH1F* h_eff_L1_jetpt_full = (TH1F*) h_dpt_matchtrk_jet_full->Clone();
    h_eff_L1_jetpt_full->SetName("eff_L1_jetpt_full");
    h_eff_L1_jetpt_full->GetYaxis()->SetTitle("Purity");
    h_eff_L1_jetpt_full->Divide(h_dpt_matchtrk_jet_full, h_dpt_trk_jet_full, 1.0, 1.0, "B");


  // ---------------------------------------------------------------------------------------------------------
  // Set Axis Range

  // Efficiency - Tracking Particles
  h_eff_TP_pt    ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_pt_L  ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_pt_LC ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_pt_H  ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_eta   ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_eta_L ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_eta_H ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_phi   ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_z0    ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_z0_L  ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_z0_H  ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_d0    ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_absd0 ->SetAxisRange(0,1.1,"Y");
  if (type.Contains("Electron") || type.Contains("Pion")) h_eff_TP_pt->SetAxisRange(0,49,"X");

  // Efficiency - Tracking Particles - Jet Axis
  h_eff_TP_jetaxis_eta   ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetaxis_eta_L ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetaxis_eta_H ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetaxis_phi   ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetaxis_phi_L ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetaxis_phi_H ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetaxis_pt    ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetaxis_pt_L  ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetaxis_pt_H  ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetaxis_dR    ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetaxis_dR_L  ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetaxis_dR_H  ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetpt         ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetpt_C       ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetpt_I       ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetpt_F       ->SetAxisRange(0,1.1,"Y");

  h_eff_TP_jetpt_L       ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetpt_LC      ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetpt_LI      ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetpt_LF      ->SetAxisRange(0,1.1,"Y");

  h_eff_TP_jetpt_H       ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetpt_HC      ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetpt_HI      ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetpt_HF      ->SetAxisRange(0,1.1,"Y");

  h_eff_TP_jetpt_full    ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetpt_full_C  ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetpt_full_I  ->SetAxisRange(0,1.1,"Y");
  h_eff_TP_jetpt_full_F  ->SetAxisRange(0,1.1,"Y");

  // Purity - L1 Tracks
  h_eff_L1_pt    ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_pt_L  ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_pt_LC ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_pt_H  ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_eta   ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_eta_L ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_eta_H ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_phi   ->SetAxisRange(0,1.1,"Y");
  if (type.Contains("Electron") || type.Contains("Pion")) h_eff_L1_pt->SetAxisRange(0,49,"X");

  // Purity - L1 Tracks - Jet Axis
  h_eff_L1_jetaxis_eta   ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_jetaxis_eta_L ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_jetaxis_eta_H ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_jetaxis_phi   ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_jetaxis_phi_L ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_jetaxis_phi_H ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_jetaxis_pt    ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_jetaxis_pt_L  ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_jetaxis_pt_H  ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_jetaxis_dR    ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_jetaxis_dR_L  ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_jetaxis_dR_H  ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_jetpt         ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_jetpt_L       ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_jetpt_H       ->SetAxisRange(0,1.1,"Y");
  h_eff_L1_jetpt_full    ->SetAxisRange(0,1.1,"Y");

  gPad->SetGridx();
  gPad->SetGridy();


  // ---------------------------------------------------------------------------------------------------------
  // Efficiency Plots - Tracking Particles

  //-----eff_TP_pt-----//
  h_eff_TP_pt->SetTitle(htitle);
  h_eff_TP_pt->Draw();
  h_eff_TP_pt->Write();
  c.SaveAs(DIR+type+"_eff_TP_pt.png");

    if (type.Contains("Mu")) {
      h_eff_TP_pt->SetTitle(htitle);
      h_eff_TP_pt->GetYaxis()->SetRangeUser(0.8,1.01);
      c.SaveAs(DIR+type+"_eff_TP_pt_zoom.png");
    }

    if (doDetailedPlots) {
      h_eff_TP_pt_L->SetTitle(htitle);
      h_eff_TP_pt_L->Draw();
      h_eff_TP_pt_L->Write();
      sprintf(ctxt,"p_{T} < 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
      c.SaveAs(DIR+type+"_eff_TP_pt_L.png");

      h_eff_TP_pt_LC->SetTitle(htitle);
      h_eff_TP_pt_LC->Draw();
      h_eff_TP_pt_LC->Write();
      sprintf(ctxt,"p_{T} < 8 GeV, |#eta|<1.0");
      mySmallText(0.45,0.5,1,ctxt);
      c.SaveAs(DIR+type+"_eff_TP_pt_LC.png");

      h_eff_TP_pt_H->SetTitle(htitle);
      h_eff_TP_pt_H->Draw();
      h_eff_TP_pt_H->Write();
      sprintf(ctxt,"p_{T} > 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
      c.SaveAs(DIR+type+"_eff_TP_pt_H.png");
    }

  //-----eff_TP_eta-----//
  h_eff_TP_eta->SetTitle(htitle);
  h_eff_TP_eta->Draw();
  h_eff_TP_eta->Write();
  c.SaveAs(DIR+type+"_eff_TP_eta.png");

    if (type.Contains("Mu")) {
      h_eff_TP_eta->SetTitle(htitle);
      h_eff_TP_eta->GetYaxis()->SetRangeUser(0.8,1.01);
      c.SaveAs(DIR+type+"_eff_TP_eta_zoom.png");
    }

    if (doDetailedPlots) {
      h_eff_TP_eta_L->SetTitle(htitle);
      h_eff_TP_eta_L->Draw();
      h_eff_TP_eta_L->Write();
      sprintf(ctxt,"p_{T} < 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
      c.SaveAs(DIR+type+"_eff_TP_eta_L.png");

      h_eff_TP_eta_H->SetTitle(htitle);
      h_eff_TP_eta_H->Draw();
      h_eff_TP_eta_H->Write();
      sprintf(ctxt,"p_{T} > 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
      c.SaveAs(DIR+type+"_eff_TP_eta_H.png");
    }

  //-----eff_phi-----//
  h_eff_TP_phi->SetTitle(htitle);
  h_eff_TP_phi->Draw();
  h_eff_TP_phi->Write();
  c.SaveAs(DIR+type+"_eff_TP_phi.png");

    if (type.Contains("Mu")) {
      h_eff_TP_phi->SetTitle(htitle);
      h_eff_TP_phi->GetYaxis()->SetRangeUser(0.8,1.01);
      c.SaveAs(DIR+type+"_eff_TP_phi_zoom.png");
    }
  

  //-----eff_z0-----//
  if (doDetailedPlots) {
    h_eff_TP_z0->SetTitle(htitle);
    h_eff_TP_z0->Draw();
    h_eff_TP_z0->Write();
    c.SaveAs(DIR+type+"_eff_TP_z0.png");

      h_eff_TP_z0_L->SetTitle(htitle);
      h_eff_TP_z0_L->Draw();
      h_eff_TP_z0_L->Write();
      c.SaveAs(DIR+type+"_eff_TP_z0_L.png");

      h_eff_TP_z0_H->SetTitle(htitle);
      h_eff_TP_z0_H->Draw();
      h_eff_TP_z0_H->Write();
      c.SaveAs(DIR+type+"_eff_TP_z0_H.png");
  }

  //-----eff_d0-----//
  if (doDetailedPlots) {
    h_eff_TP_d0->SetTitle(htitle);
    h_eff_TP_d0->Draw();
    h_eff_TP_d0->Write();
    c.SaveAs(DIR+type+"_eff_TP_d0.png");

    h_eff_TP_absd0->SetTitle(htitle);
    h_eff_TP_absd0->Draw();
    h_eff_TP_absd0->Write();
    c.SaveAs(DIR+type+"_eff_TP_absd0.png");
  }

  // ---------------------------------------------------------------------------------------------------------
  // Efficiency Plots - Tracking Particles - Jet Axis

  // eff_TP_jetaxis_eta
  h_eff_TP_jetaxis_eta->SetTitle(htitle);
  h_eff_TP_jetaxis_eta->Draw();
  h_eff_TP_jetaxis_eta->Write();
    h_eff_TP_jetaxis_eta->GetXaxis()->SetRangeUser(-0.501,0.501);
  c.SaveAs(DIR+type+"_eff_TP_jetaxis_eta.png");

  if (doDetailedPlots) {
    h_eff_TP_jetaxis_eta_L->SetTitle(htitle);
    h_eff_TP_jetaxis_eta_L->Draw();
    h_eff_TP_jetaxis_eta_L->Write();
      sprintf(ctxt,"p_{T} < 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
      h_eff_TP_jetaxis_eta_L->GetXaxis()->SetRangeUser(-0.501,0.501);
    c.SaveAs(DIR+type+"_eff_TP_jetaxis_eta_L.png");

    h_eff_TP_jetaxis_eta_H->SetTitle(htitle);
    h_eff_TP_jetaxis_eta_H->Draw();
    h_eff_TP_jetaxis_eta_H->Write();
      sprintf(ctxt,"p_{T} > 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
      h_eff_TP_jetaxis_eta_H->GetXaxis()->SetRangeUser(-0.501,0.501);
    c.SaveAs(DIR+type+"_eff_TP_jetaxis_eta_H.png");
  }

  // eff_TP_jetaxis_phi
  h_eff_TP_jetaxis_phi->SetTitle(htitle);
  h_eff_TP_jetaxis_phi->Draw();
  h_eff_TP_jetaxis_phi->Write();
    h_eff_TP_jetaxis_phi->GetXaxis()->SetRangeUser(-0.501,0.501);
  c.SaveAs(DIR+type+"_eff_TP_jetaxis_phi.png");

  if (doDetailedPlots) {
    h_eff_TP_jetaxis_phi_L->SetTitle(htitle);
    h_eff_TP_jetaxis_phi_L->Draw();
    h_eff_TP_jetaxis_phi_L->Write();
      sprintf(ctxt,"p_{T} < 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
      h_eff_TP_jetaxis_phi_L->GetXaxis()->SetRangeUser(-0.501,0.501);
    c.SaveAs(DIR+type+"_eff_TP_jetaxis_phi_L.png");

    h_eff_TP_jetaxis_phi_H->SetTitle(htitle);
    h_eff_TP_jetaxis_phi_H->Draw();
    h_eff_TP_jetaxis_phi_H->Write();
      sprintf(ctxt,"p_{T} > 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
      h_eff_TP_jetaxis_phi_H->GetXaxis()->SetRangeUser(-0.501,0.501);
    c.SaveAs(DIR+type+"_eff_TP_jetaxis_phi_H.png");
  }

  // eff_TP_jetaxis_pt
  h_eff_TP_jetaxis_pt->SetTitle(htitle);
  h_eff_TP_jetaxis_pt->Draw();
  h_eff_TP_jetaxis_pt->Write();
  c.SaveAs(DIR+type+"_eff_TP_jetaxis_pt.png");

  if (doDetailedPlots) {
    h_eff_TP_jetaxis_pt_L->SetTitle(htitle);
    h_eff_TP_jetaxis_pt_L->Draw();
    h_eff_TP_jetaxis_pt_L->Write();
      sprintf(ctxt,"p_{T} < 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
    c.SaveAs(DIR+type+"_eff_TP_jetaxis_pt_L.png");

    h_eff_TP_jetaxis_pt_H->SetTitle(htitle);
    h_eff_TP_jetaxis_pt_H->Draw();
    h_eff_TP_jetaxis_pt_H->Write();
      sprintf(ctxt,"p_{T} > 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
    c.SaveAs(DIR+type+"_eff_TP_jetaxis_pt_H.png");
  }

  // eff_TP_jetaxis_dR
  h_eff_TP_jetaxis_dR->SetTitle(htitle);
  h_eff_TP_jetaxis_dR->Draw();
  h_eff_TP_jetaxis_dR->Write();
  c.SaveAs(DIR+type+"_eff_TP_jetaxis_dR.png");

  if (doDetailedPlots) {
    h_eff_TP_jetaxis_dR_L->SetTitle(htitle);
    h_eff_TP_jetaxis_dR_L->Draw();
    h_eff_TP_jetaxis_dR_L->Write();
      sprintf(ctxt,"p_{T} < 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
    c.SaveAs(DIR+type+"_eff_TP_jetaxis_dR_L.png");

    h_eff_TP_jetaxis_dR_H->SetTitle(htitle);
    h_eff_TP_jetaxis_dR_H->Draw();
    h_eff_TP_jetaxis_dR_H->Write();
      sprintf(ctxt,"p_{T} > 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
    c.SaveAs(DIR+type+"_eff_TP_jetaxis_dR_H.png");
  }

  // eff_TP_jetpt
  h_eff_TP_jetpt->SetTitle(htitle);
  h_eff_TP_jetpt->Draw();
  h_eff_TP_jetpt->Write();
  c.SaveAs(DIR+type+"_eff_TP_jetpt.png");

  if (doDetailedPlots) {
    h_eff_TP_jetpt_C->SetTitle(htitle);
    h_eff_TP_jetpt_C->Draw();
    h_eff_TP_jetpt_C->Write();
      sprintf(ctxtC,"Central |#eta| < 0.8");
      mySmallText(0.6,0.25,1,ctxtC);
    c.SaveAs(DIR+type+"_eff_TP_jetpt_C.png");

    h_eff_TP_jetpt_I->SetTitle(htitle);
    h_eff_TP_jetpt_I->Draw();
    h_eff_TP_jetpt_I->Write();
      sprintf(ctxtI,"Intermediate 0.8 < |#eta| < 1.6");
      mySmallText(0.6,0.25,1,ctxtI);
    c.SaveAs(DIR+type+"_eff_TP_jetpt_I.png");

    h_eff_TP_jetpt_F->SetTitle(htitle);
    h_eff_TP_jetpt_F->Draw();
    h_eff_TP_jetpt_F->Write();
      sprintf(ctxtF,"Forward |#eta| > 1.6");
      mySmallText(0.6,0.25,1,ctxtF);
    c.SaveAs(DIR+type+"_eff_TP_jetpt_F.png");

    h_eff_TP_jetpt_L->SetTitle(htitle);
    h_eff_TP_jetpt_L->Draw();
    h_eff_TP_jetpt_L->Write();
      sprintf(ctxt,"TP p_{T} < 8 GeV");
      mySmallText(0.6,0.3,1,ctxt);
    c.SaveAs(DIR+type+"_eff_TP_jetpt_L.png");

    h_eff_TP_jetpt_LC->SetTitle(htitle);
    h_eff_TP_jetpt_LC->Draw();
    h_eff_TP_jetpt_LC->Write();
      sprintf(ctxt,"TP p_{T} < 8 GeV");
      mySmallText(0.6,0.3,1,ctxt);
      sprintf(ctxtC,"Central |#eta| < 0.8");
      mySmallText(0.6,0.25,1,ctxtC);
    //c.SaveAs(DIR+type+"_eff_TP_jetpt_LC.png");

    h_eff_TP_jetpt_LI->SetTitle(htitle);
    h_eff_TP_jetpt_LI->Draw();
    h_eff_TP_jetpt_LI->Write();
      sprintf(ctxt,"TP p_{T} < 8 GeV");
      mySmallText(0.6,0.3,1,ctxt);
      sprintf(ctxtI,"Intermediate 0.8 < |#eta| < 1.6");
      mySmallText(0.6,0.25,1,ctxtI);
    //c.SaveAs(DIR+type+"_eff_TP_jetpt_LI.png");

    h_eff_TP_jetpt_LF->SetTitle(htitle);
    h_eff_TP_jetpt_LF->Draw();
    h_eff_TP_jetpt_LF->Write();
      sprintf(ctxt,"TP p_{T} < 8 GeV");
      mySmallText(0.6,0.3,1,ctxt);
      sprintf(ctxtF,"Forward |#eta| > 1.6");
      mySmallText(0.6,0.25,1,ctxtF);
    //c.SaveAs(DIR+type+"_eff_TP_jetpt_LF.png");

    h_eff_TP_jetpt_H->SetTitle(htitle);
    h_eff_TP_jetpt_H->Draw();
    h_eff_TP_jetpt_H->Write();
      sprintf(ctxt,"TP p_{T} > 8 GeV");
      mySmallText(0.6,0.3,1,ctxt);
    c.SaveAs(DIR+type+"_eff_TP_jetpt_H.png");

    h_eff_TP_jetpt_HC->SetTitle(htitle);
    h_eff_TP_jetpt_HC->Draw();
    h_eff_TP_jetpt_HC->Write();
      sprintf(ctxt,"TP p_{T} > 8 GeV");
      mySmallText(0.6,0.3,1,ctxt);
      sprintf(ctxtC,"Central |#eta| < 0.8");
      mySmallText(0.6,0.25,1,ctxtC);
    //c.SaveAs(DIR+type+"_eff_TP_jetpt_HC.png");

    h_eff_TP_jetpt_HI->SetTitle(htitle);
    h_eff_TP_jetpt_HI->Draw();
    h_eff_TP_jetpt_HI->Write();
      sprintf(ctxt,"TP p_{T} > 8 GeV");
      mySmallText(0.6,0.3,1,ctxt);
      sprintf(ctxtI,"Intermediate 0.8 < |#eta| < 1.6");
      mySmallText(0.6,0.25,1,ctxtI);
    //c.SaveAs(DIR+type+"_eff_TP_jetpt_HI.png");

    h_eff_TP_jetpt_HF->SetTitle(htitle);
    h_eff_TP_jetpt_HF->Draw();
    h_eff_TP_jetpt_HF->Write();
      sprintf(ctxt,"TP p_{T} > 8 GeV");
      mySmallText(0.6,0.3,1,ctxt);
      sprintf(ctxtF,"Forward |#eta| > 1.6");
      mySmallText(0.6,0.25,1,ctxtF);
    //c.SaveAs(DIR+type+"_eff_TP_jetpt_HF.png");

    h_eff_TP_jetpt_full->SetTitle(htitle);
    h_eff_TP_jetpt_full->Draw();
    h_eff_TP_jetpt_full->Write();
    c.SaveAs(DIR+type+"_eff_TP_jetpt_full.png");

    h_eff_TP_jetpt_full_C->SetTitle(htitle);
    h_eff_TP_jetpt_full_C->Draw();
    h_eff_TP_jetpt_full_C->Write();
      sprintf(ctxtC,"Central |#eta| < 0.8");
      mySmallText(0.6,0.25,1,ctxtC);
    //c.SaveAs(DIR+type+"_eff_TP_jetpt_full_C.png");

    h_eff_TP_jetpt_full_I->SetTitle(htitle);
    h_eff_TP_jetpt_full_I->Draw();
    h_eff_TP_jetpt_full_I->Write();
      sprintf(ctxtI,"Intermediate 0.8 < |#eta| < 1.6");
      mySmallText(0.6,0.25,1,ctxtI);
    //c.SaveAs(DIR+type+"_eff_TP_jetpt_full_I.png");

    h_eff_TP_jetpt_full_F->SetTitle(htitle);
    h_eff_TP_jetpt_full_F->Draw();
    h_eff_TP_jetpt_full_F->Write();
      sprintf(ctxtF,"Forward |#eta| > 1.6");
      mySmallText(0.6,0.25,1,ctxtF);
    //c.SaveAs(DIR+type+"_eff_TP_jetpt_full_F.png");
  }

  // ---------------------------------------------------------------------------------------------------------
  // Purity Plots - L1 Tracks

  //-----eff_L1_pt-----//
  h_eff_L1_pt->SetTitle(htitle);
  h_eff_L1_pt->Draw();
  h_eff_L1_pt->Write();
  c.SaveAs(DIR+type+"_eff_L1_pt.png");

  if (type.Contains("Mu")) {
    h_eff_L1_pt->SetTitle(htitle);
    h_eff_L1_pt->GetYaxis()->SetRangeUser(0.8,1.01);
    c.SaveAs(DIR+type+"_eff_L1_pt_zoom.png");
  }

  if (doDetailedPlots) {
    h_eff_L1_pt_L->SetTitle(htitle);
    h_eff_L1_pt_L->Draw();
    h_eff_L1_pt_L->Write();
      sprintf(ctxt,"p_{T} < 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
    c.SaveAs(DIR+type+"_eff_L1_pt_L.png");

    h_eff_L1_pt_LC->SetTitle(htitle);
    h_eff_L1_pt_LC->Draw();
    h_eff_L1_pt_LC->Write();
      sprintf(ctxt,"p_{T} < 8 GeV, |#eta|<1.0");
      mySmallText(0.45,0.5,1,ctxt);
    c.SaveAs(DIR+type+"_eff_L1_pt_LC.png");

    h_eff_L1_pt_H->SetTitle(htitle);
    h_eff_L1_pt_H->Draw();
    h_eff_L1_pt_H->Write();
      sprintf(ctxt,"p_{T} > 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
    c.SaveAs(DIR+type+"_eff_L1_pt_H.png");
  }
  

  //-----eff_L1_eta-----//
  h_eff_L1_eta->SetTitle(htitle);
  h_eff_L1_eta->Draw();
  h_eff_L1_eta->Write();
  c.SaveAs(DIR+type+"_eff_L1_eta.png");

    if (type.Contains("Mu")) {
      h_eff_L1_eta->SetTitle(htitle);
      h_eff_L1_eta->GetYaxis()->SetRangeUser(0.8,1.01);
      c.SaveAs(DIR+type+"_eff_L1_eta_zoom.png");
    }

    if (doDetailedPlots) {
      h_eff_L1_eta_L->SetTitle(htitle);
      h_eff_L1_eta_L->Draw();
      h_eff_L1_eta_L->Write();
        sprintf(ctxt,"p_{T} < 8 GeV");
        mySmallText(0.45,0.5,1,ctxt);
      c.SaveAs(DIR+type+"_eff_L1_eta_L.png");

      h_eff_L1_eta_H->SetTitle(htitle);
      h_eff_L1_eta_H->Draw();
      h_eff_L1_eta_H->Write();
      sprintf(ctxt,"p_{T} > 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
      c.SaveAs(DIR+type+"_eff_L1_eta_H.png");
    }

  //-----eff_L1_phi-----//
  h_eff_L1_phi->SetTitle(htitle);
  h_eff_L1_phi->Draw();
  h_eff_L1_phi->Write();
  c.SaveAs(DIR+type+"_eff_L1_phi.png");

    if (type.Contains("Mu")) {
      h_eff_L1_phi->SetTitle(htitle);
      h_eff_L1_phi->GetYaxis()->SetRangeUser(0.8,1.01);
      c.SaveAs(DIR+type+"_eff_L1_phi_zoom.png");
    }


  // ---------------------------------------------------------------------------------------------------------
  // Purity Plots - L1 Tracks - Jet Axis

  // eff_L1_jetaxis_eta
  h_eff_L1_jetaxis_eta->SetTitle(htitle);
  h_eff_L1_jetaxis_eta->Draw();
  h_eff_L1_jetaxis_eta->Write();
    h_eff_L1_jetaxis_eta->GetXaxis()->SetRangeUser(-0.501,0.501);
  c.SaveAs(DIR+type+"_eff_L1_jetaxis_eta.png");

  if (doDetailedPlots) {
    h_eff_L1_jetaxis_eta_L->SetTitle(htitle);
    h_eff_L1_jetaxis_eta_L->Draw();
    h_eff_L1_jetaxis_eta_L->Write();
      sprintf(ctxt,"p_{T} < 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
      h_eff_L1_jetaxis_eta_L->GetXaxis()->SetRangeUser(-0.501,0.501);
    c.SaveAs(DIR+type+"_eff_L1_jetaxis_eta_L.png");

    h_eff_L1_jetaxis_eta_H->SetTitle(htitle);
    h_eff_L1_jetaxis_eta_H->Draw();
    h_eff_L1_jetaxis_eta_H->Write();
      sprintf(ctxt,"p_{T} > 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
      h_eff_L1_jetaxis_eta_H->GetXaxis()->SetRangeUser(-0.501,0.501);
    c.SaveAs(DIR+type+"_eff_L1_jetaxis_eta_H.png");
  }

  // eff_L1_jetaxis_phi
  h_eff_L1_jetaxis_phi->SetTitle(htitle);
  h_eff_L1_jetaxis_phi->Draw();
  h_eff_L1_jetaxis_phi->Write();
    h_eff_L1_jetaxis_phi->GetXaxis()->SetRangeUser(-0.501,0.501);
  c.SaveAs(DIR+type+"_eff_L1_jetaxis_phi.png");

  if (doDetailedPlots) {
    h_eff_L1_jetaxis_phi_L->SetTitle(htitle);
    h_eff_L1_jetaxis_phi_L->Draw();
    h_eff_L1_jetaxis_phi_L->Write();
      sprintf(ctxt,"p_{T} < 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
      h_eff_L1_jetaxis_phi_L->GetXaxis()->SetRangeUser(-0.501,0.501);
    c.SaveAs(DIR+type+"_eff_L1_jetaxis_phi_L.png");

    h_eff_L1_jetaxis_phi_H->SetTitle(htitle);
    h_eff_L1_jetaxis_phi_H->Draw();
    h_eff_L1_jetaxis_phi_H->Write();
      sprintf(ctxt,"p_{T} > 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
      h_eff_L1_jetaxis_phi_H->GetXaxis()->SetRangeUser(-0.501,0.501);
    c.SaveAs(DIR+type+"_eff_L1_jetaxis_phi_H.png");
  }

  // eff_L1_jetaxis_pt
  h_eff_L1_jetaxis_pt->SetTitle(htitle);
  h_eff_L1_jetaxis_pt->Draw();
  h_eff_L1_jetaxis_pt->Write();
  c.SaveAs(DIR+type+"_eff_L1_jetaxis_pt.png");

  if (doDetailedPlots) {
    h_eff_L1_jetaxis_pt_L->SetTitle(htitle);
    h_eff_L1_jetaxis_pt_L->Draw();
    h_eff_L1_jetaxis_pt_L->Write();
    sprintf(ctxt,"p_{T} < 8 GeV");
    mySmallText(0.45,0.5,1,ctxt);
    c.SaveAs(DIR+type+"_eff_L1_jetaxis_pt_L.png");

    h_eff_L1_jetaxis_pt_H->SetTitle(htitle);
    h_eff_L1_jetaxis_pt_H->Draw();
    h_eff_L1_jetaxis_pt_H->Write();
    sprintf(ctxt,"p_{T} > 8 GeV");
    mySmallText(0.45,0.5,1,ctxt);
    c.SaveAs(DIR+type+"_eff_L1_jetaxis_pt_H.png");
  }

  // eff_L1_jetaxis_dR
  h_eff_L1_jetaxis_dR->SetTitle(htitle);
  h_eff_L1_jetaxis_dR->Draw();
  h_eff_L1_jetaxis_dR->Write();
  c.SaveAs(DIR+type+"_eff_L1_jetaxis_dR.png");

  if (doDetailedPlots) {
    h_eff_L1_jetaxis_dR_L->SetTitle(htitle);
    h_eff_L1_jetaxis_dR_L->Draw();
    h_eff_L1_jetaxis_dR_L->Write();
      sprintf(ctxt,"p_{T} < 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
    c.SaveAs(DIR+type+"_eff_L1_jetaxis_dR_L.png");

    h_eff_L1_jetaxis_dR_H->SetTitle(htitle);
    h_eff_L1_jetaxis_dR_H->Draw();
    h_eff_L1_jetaxis_dR_H->Write();
      sprintf(ctxt,"p_{T} > 8 GeV");
      mySmallText(0.45,0.5,1,ctxt);
    c.SaveAs(DIR+type+"_eff_L1_jetaxis_dR_H.png");
  }

  // eff_L1_jetpt
  h_eff_L1_jetpt->SetTitle(htitle);
  h_eff_L1_jetpt->Draw();
  h_eff_L1_jetpt->Write();
  c.SaveAs(DIR+type+"_eff_L1_jetpt.png");

  if (doDetailedPlots) {
    h_eff_L1_jetpt_L->SetTitle(htitle);
    h_eff_L1_jetpt_L->Draw();
    h_eff_L1_jetpt_L->Write();
    sprintf(ctxt,"L1 Tracks p_{T} < 8 GeV");
    mySmallText(0.45,0.5,1,ctxt);
    c.SaveAs(DIR+type+"_eff_L1_jetpt_L.png");

    h_eff_L1_jetpt_H->SetTitle(htitle);
    h_eff_L1_jetpt_H->Draw();
    h_eff_L1_jetpt_H->Write();
    sprintf(ctxt,"L1 Tracks p_{T} > 8 GeV");
    mySmallText(0.45,0.5,1,ctxt);
    c.SaveAs(DIR+type+"_eff_L1_jetpt_H.png");

    h_eff_L1_jetpt_full->SetTitle(htitle);
    h_eff_L1_jetpt_full->Draw();
    h_eff_L1_jetpt_full->Write();
    c.SaveAs(DIR+type+"_eff_L1_jetpt_full.png");
  }



  // ---------------------------------------------------------------------------------------------------------
  // Remove Grid Lines
  // ---------------------------------------------------------------------------------------------------------
  gPad->SetGridx(0);
  gPad->SetGridy(0);



  // ---------------------------------------------------------------------------------------------------------
  // Total Track Rates vs pT
  // ---------------------------------------------------------------------------------------------------------
  h_trk_vspt->Scale(1.0/nevt);
  h_tp_vspt->Scale(1.0/nevt);

  h_tp_vspt->SetTitle(htitle);
  h_tp_vspt->GetYaxis()->SetTitle("Tracks / event");
  h_tp_vspt->GetXaxis()->SetTitle("Track p_{T} [GeV]");
  h_tp_vspt->SetLineColor(4);
  h_tp_vspt->SetMarkerColor(4);
  h_tp_vspt->SetLineStyle(2);

  float max = h_tp_vspt->GetMaximum();
  if (h_trk_vspt->GetMaximum() > max) max = h_trk_vspt->GetMaximum();
  h_tp_vspt->SetAxisRange(0,max*1.05,"Y");

  h_tp_vspt ->Draw();
  h_trk_vspt->Draw("same");
  h_tp_vspt ->Draw("same");

  h_trk_vspt->Write();
  h_tp_vspt ->Write();

  char txt[500];
  sprintf(txt,"average # tracks/event = %.1f",h_trk_vspt->GetSum());
  mySmallText(0.5,0.85,1,txt);
  char txt3[500];
  sprintf(txt3,"average # TPs(stubs in #geq 4 layers)/");
  char txt2[500];
  sprintf(txt2,"event = %.1f",h_tp_vspt->GetSum());
  mySmallText(0.5,0.79,4,txt3);
  mySmallText(0.5,0.74,4,txt2);

  c.SaveAs(DIR+type+"_trackrate_vspt.png");



  // ---------------------------------------------------------------------------------------------------------
  // Jet  Plots
  // ---------------------------------------------------------------------------------------------------------



  // ---------------------------------------------------------------------------------------------------------
  // Jet Profile Plots

  // Avoid overwriting
  TProfile* hm_jet_profile_TP  = (TProfile*) h_jet_profile_TP  ->Clone();
  TProfile* hm_jet_profile_MTP = (TProfile*) h_jet_profile_MTP ->Clone();
  TProfile* hm_jet_profile_L1  = (TProfile*) h_jet_profile_L1  ->Clone();
  TProfile* hm_jet_profile_ML1 = (TProfile*) h_jet_profile_ML1 ->Clone();

  if (doDetailedPlots) {
    //h_jet_profile_TP
    h_jet_profile_TP->SetTitle(htitle);
    h_jet_profile_TP->SetMarkerColor(kViolet+3); //deep purple
    h_jet_profile_TP->Draw();
    h_jet_profile_TP->Write();
     //c.SaveAs(DIR+type+"_jet_profile_TP.png");

    //h_jet_profile_MTP
    h_jet_profile_MTP->SetTitle(htitle);
    h_jet_profile_MTP->SetMarkerColor(kAzure-5); //robin's egg blue
    h_jet_profile_MTP->Draw();
    h_jet_profile_MTP->Write();
     //c.SaveAs(DIR+type+"_jet_profile_MTP.png");

    //h_jet_profile_L1
    h_jet_profile_L1->SetTitle(htitle);
    h_jet_profile_L1->SetMarkerColor(kTeal-5); //teal
    h_jet_profile_L1->Draw();
    h_jet_profile_L1->Write();
     //c.SaveAs(DIR+type+"_jet_profile_L1.png");

    //h_jet_profile_L1_matched
    h_jet_profile_ML1->SetTitle(htitle);
    h_jet_profile_ML1->SetMarkerColor(kYellow-6); //mustard
    h_jet_profile_ML1->Draw();
    h_jet_profile_ML1->Write();
     //c.SaveAs(DIR+type+"_jet_profile_ML1.png");
  }


  //-----Jet Profile - Multiplot-----//
  hm_jet_profile_TP->SetTitle(htitle);
  hm_jet_profile_TP->GetYaxis()->SetTitle("Tracks per Jet");
  hm_jet_profile_TP->GetXaxis()->SetTitle("Jet p_{T} [GeV]");

  // Coloring
  hm_jet_profile_TP->SetLineColor(kViolet+3);
  hm_jet_profile_TP->SetMarkerColor(kViolet+3);
  hm_jet_profile_TP->SetMarkerStyle(kFullCircle);

  hm_jet_profile_MTP->SetLineColor(kAzure-5);
  hm_jet_profile_MTP->SetMarkerColor(kAzure-5);
  hm_jet_profile_MTP->SetMarkerStyle(kFullSquare);

  hm_jet_profile_L1->SetLineColor(kTeal-5);
  hm_jet_profile_L1->SetMarkerColor(kTeal-5);
  hm_jet_profile_L1->SetMarkerStyle(kFullTriangleUp);

  hm_jet_profile_ML1->SetLineColor(kYellow-6);
  hm_jet_profile_ML1->SetMarkerColor(kYellow-6);
  hm_jet_profile_ML1->SetMarkerStyle(kFullDiamond);

  // Get maximum for Y
  float profmax = hm_jet_profile_TP->GetMaximum();
  if (hm_jet_profile_MTP ->GetMaximum() > profmax) profmax = hm_jet_profile_MTP ->GetMaximum();
  if (hm_jet_profile_L1  ->GetMaximum() > profmax) profmax = hm_jet_profile_L1  ->GetMaximum();
  if (hm_jet_profile_ML1 ->GetMaximum() > profmax) profmax = hm_jet_profile_ML1 ->GetMaximum();
  hm_jet_profile_TP->SetAxisRange(0,profmax*1.05,"Y");

  // Draw
  hm_jet_profile_TP  ->Draw();
  hm_jet_profile_MTP ->Draw("same");
  hm_jet_profile_L1  ->Draw("same");
  hm_jet_profile_ML1 ->Draw("same");

    // Legend
    float mean_jet_profile_TP  = 0;
    float mean_jet_profile_MTP = 0;
    float mean_jet_profile_L1  = 0;
    float mean_jet_profile_ML1 = 0;

    mean_jet_profile_TP  = hm_jet_profile_TP ->GetMean(2);
    mean_jet_profile_MTP = hm_jet_profile_MTP->GetMean(2);
    mean_jet_profile_L1  = hm_jet_profile_L1 ->GetMean(2);
    mean_jet_profile_ML1 = hm_jet_profile_ML1->GetMean(2);

    TLegend *leg_jet_profile = new TLegend(0.17,0.70,0.35,0.90);
    leg_jet_profile->AddEntry(hm_jet_profile_TP,  Form("TP mean=%.3f",        mean_jet_profile_MTP), "p");
    leg_jet_profile->AddEntry(hm_jet_profile_MTP, Form("Matched TP mean=%.3f",mean_jet_profile_MTP), "p");
    leg_jet_profile->AddEntry(hm_jet_profile_L1,  Form("L1 Track mean=%.3f",  mean_jet_profile_L1),  "p");
    leg_jet_profile->AddEntry(hm_jet_profile_ML1, Form("Matched L1 mean=%.3f",mean_jet_profile_ML1), "p");
    leg_jet_profile->Draw();


    c.SaveAs(DIR+type+"_jet_profile.png");



  // ---------------------------------------------------------------------------------------------------------
  // Jet Sum pT Plots

  if(doDetailedPlots){
    //-----Jet Sum pT vs pT-----//
    h_jet_tp_sumpt_vspt->SetTitle(htitle);
    h_jet_tp_sumpt_vspt->Draw();
    h_jet_tp_sumpt_vspt->Write();
    c.SaveAs(DIR+type+"_jet_tp_sumpt_vspt.png");

    h_jet_trk_sumpt_vspt->SetTitle(htitle);
    h_jet_trk_sumpt_vspt->Draw();
    h_jet_trk_sumpt_vspt->Write();
    c.SaveAs(DIR+type+"_jet_trk_sumpt_vspt.png");

    h_jet_matchtrk_sumpt_vspt->SetTitle(htitle);
    h_jet_matchtrk_sumpt_vspt->Draw();
    h_jet_matchtrk_sumpt_vspt->Write();
    c.SaveAs(DIR+type+"_jet_matchtrk_sumpt_vspt.png");


    //-----Jet Sum pT vs Eta-----//
    h_jet_tp_sumpt_vseta->SetTitle(htitle);
    h_jet_tp_sumpt_vseta->Draw();
    h_jet_tp_sumpt_vseta->Write();
    c.SaveAs(DIR+type+"_jet_tp_sumpt_vseta.png");

    h_jet_trk_sumpt_vseta->SetTitle(htitle);
    h_jet_trk_sumpt_vseta->Draw();
    h_jet_trk_sumpt_vseta->Write();
    c.SaveAs(DIR+type+"_jet_trk_sumpt_vseta.png");

    h_jet_matchtrk_sumpt_vseta->SetTitle(htitle);
    h_jet_matchtrk_sumpt_vseta->Draw();
    h_jet_matchtrk_sumpt_vseta->Write();
    c.SaveAs(DIR+type+"_jet_matchtrk_sumpt_vseta.png");


    //-----Jet Sum pT vs Phi-----//
    h_jet_tp_sumpt_vsphi->SetTitle(htitle);
    h_jet_tp_sumpt_vsphi->Draw();
    h_jet_tp_sumpt_vsphi->Write();
    c.SaveAs(DIR+type+"_jet_tp_sumpt_vsphi.png");

    h_jet_trk_sumpt_vsphi->SetTitle(htitle);
    h_jet_trk_sumpt_vsphi->Draw();
    h_jet_trk_sumpt_vsphi->Write();
    c.SaveAs(DIR+type+"_jet_trk_sumpt_vsphi.png");

    h_jet_matchtrk_sumpt_vsphi->SetTitle(htitle);
    h_jet_matchtrk_sumpt_vsphi->Draw();
    h_jet_matchtrk_sumpt_vsphi->Write();
    c.SaveAs(DIR+type+"_jet_matchtrk_sumpt_vsphi.png");
  }


  // ---------------------------------------------------------------------------------------------------------
  // Fractional Jet Plots

  if(doDetailedPlots){
    //-----Jet Frac pT-----//
    TH1F* h_jet_frac_sumpt_vspt = (TH1F*) h_jet_trk_sumpt_vspt->Clone();
    h_jet_frac_sumpt_vspt->SetName("jet_frac_sumpt_vspt");
    h_jet_frac_sumpt_vspt->GetYaxis()->SetTitle("L1 sum(p_{T}) / TP sum(p_{T})");
    h_jet_frac_sumpt_vspt->Divide(h_jet_trk_sumpt_vspt, h_jet_tp_sumpt_vspt, 1.0, 1.0, "B");
    h_jet_frac_sumpt_vspt->SetTitle(htitle);

    h_jet_frac_sumpt_vspt->Draw();
    h_jet_frac_sumpt_vspt->Write();
    c.SaveAs(DIR+type+"_jet_frac_sumpt_vspt.png");


    //-----Jet Frac Eta-----//
    TH1F* h_jet_frac_sumpt_vseta = (TH1F*) h_jet_trk_sumpt_vseta->Clone();
    h_jet_frac_sumpt_vseta->SetName("jet_frac_sumpt_vseta");
    h_jet_frac_sumpt_vseta->GetYaxis()->SetTitle("L1 sum(p_{T}) / TP sum(p_{T})");
    h_jet_frac_sumpt_vseta->Divide(h_jet_trk_sumpt_vseta, h_jet_tp_sumpt_vseta, 1.0, 1.0, "B");
    h_jet_frac_sumpt_vseta->SetTitle(htitle);

    h_jet_frac_sumpt_vseta->Draw();
    h_jet_frac_sumpt_vseta->Write();
    c.SaveAs(DIR+type+"_jet_frac_sumpt_vseta.png");


    //-----Jet Frac Phi-----//
    TH1F* h_jet_frac_sumpt_vsphi = (TH1F*) h_jet_trk_sumpt_vsphi->Clone();
    h_jet_frac_sumpt_vsphi->SetName("jet_frac_sumpt_vsphi");
    h_jet_frac_sumpt_vsphi->GetYaxis()->SetTitle("L1 sum(p_{T}) / TP sum(p_{T})");
    h_jet_frac_sumpt_vsphi->Divide(h_jet_trk_sumpt_vsphi, h_jet_tp_sumpt_vsphi, 1.0, 1.0, "B");
    h_jet_frac_sumpt_vsphi->SetTitle(htitle);

    h_jet_frac_sumpt_vsphi->Draw();
    h_jet_frac_sumpt_vsphi->Write();
    c.SaveAs(DIR+type+"_jet_frac_sumpt_vsphi.png");


    //-----Jet MatchFrac pT-----//
    TH1F* h_jet_matchfrac_sumpt_vspt = (TH1F*) h_jet_matchtrk_sumpt_vspt->Clone();
    h_jet_matchfrac_sumpt_vspt->SetName("jet_matchfrac_sumpt_vspt");
    h_jet_matchfrac_sumpt_vspt->GetYaxis()->SetTitle("Matched L1 sum(p_{T}) / TP sum(p_{T})");
    h_jet_matchfrac_sumpt_vspt->Divide(h_jet_matchtrk_sumpt_vspt, h_jet_tp_sumpt_vspt, 1.0, 1.0, "B");
    h_jet_matchfrac_sumpt_vspt->SetTitle(htitle);

    h_jet_matchfrac_sumpt_vspt->Draw();
    h_jet_matchfrac_sumpt_vspt->Write();
    c.SaveAs(DIR+type+"_jet_matchfrac_sumpt_vspt.png");


    //-----Jet MatchFrac Eta-----//
    TH1F* h_jet_matchfrac_sumpt_vseta = (TH1F*) h_jet_matchtrk_sumpt_vseta->Clone();
    h_jet_matchfrac_sumpt_vseta->SetName("jet_matchfrac_sumpt_vseta");
    h_jet_matchfrac_sumpt_vseta->GetYaxis()->SetTitle("Matched L1 sum(p_{T}) / TP sum(p_{T})");
    h_jet_matchfrac_sumpt_vseta->Divide(h_jet_matchtrk_sumpt_vseta, h_jet_tp_sumpt_vseta, 1.0, 1.0, "B");
    h_jet_matchfrac_sumpt_vseta->SetTitle(htitle);

    h_jet_matchfrac_sumpt_vseta->Draw();
    h_jet_matchfrac_sumpt_vseta->Write();
    c.SaveAs(DIR+type+"_jet_matchfrac_sumpt_vseta.png");


    //-----Jet MatchFrac Phi-----//
    TH1F* h_jet_matchfrac_sumpt_vsphi = (TH1F*) h_jet_matchtrk_sumpt_vsphi->Clone();
    h_jet_matchfrac_sumpt_vsphi->SetName("jet_matchfrac_sumpt_vsphi");
    h_jet_matchfrac_sumpt_vsphi->GetYaxis()->SetTitle("Matched L1 sum(p_{T}) / TP sum(p_{T})");
    h_jet_matchfrac_sumpt_vsphi->Divide(h_jet_matchtrk_sumpt_vsphi, h_jet_tp_sumpt_vsphi, 1.0, 1.0, "B");
    h_jet_matchfrac_sumpt_vsphi->SetTitle(htitle);

    h_jet_matchfrac_sumpt_vsphi->Draw();
    h_jet_matchfrac_sumpt_vsphi->Write();
      c.SaveAs(DIR+type+"_jet_matchfrac_sumpt_vsphi.png");
  }


  // ---------------------------------------------------------------------------------------------------------
  // NBR Tracks Per Event
  // ---------------------------------------------------------------------------------------------------------

  h_ntrk_pt2->SetTitle(htitle);
  h_ntrk_pt2->Draw();
    sprintf(ctxt,"Mean = %.3f",h_ntrk_pt2->GetMean());
    mySmallText(0.7,0.85,1,ctxt);
  h_ntrk_pt2->Write();
  c.SaveAs(DIR+type+"_ntrk_pt2.png");

  h_ntrk_pt3->SetTitle(htitle);
  h_ntrk_pt3->Draw();
    sprintf(ctxt,"Mean = %.3f",h_ntrk_pt3->GetMean());
    mySmallText(0.7,0.85,1,ctxt);
  h_ntrk_pt3->Write();
  c.SaveAs(DIR+type+"_ntrk_pt3.png");

  h_ntrk_pt10->SetTitle(htitle);
  h_ntrk_pt10->Draw();
    sprintf(ctxt,"Mean = %.3f",h_ntrk_pt10->GetMean());
    mySmallText(0.7,0.85,1,ctxt);
  h_ntrk_pt10->Write();
  c.SaveAs(DIR+type+"_ntrk_pt10.png");


  // ---------------------------------------------------------------------------------------------------------
  // Stub Plots
  // ---------------------------------------------------------------------------------------------------------

  // Avoid overwriting
  TH1F* h_nstub_tp   = (TH1F*) h_tp_nstub   ->Clone();
  TH1F* h_nstub_tp_C = (TH1F*) h_tp_nstub_C ->Clone();
  TH1F* h_nstub_tp_I = (TH1F*) h_tp_nstub_I ->Clone();
  TH1F* h_nstub_tp_F = (TH1F*) h_tp_nstub_F ->Clone();

  TH1F* h_nstub_trk   = (TH1F*) h_trk_nstub   ->Clone();
  TH1F* h_nstub_trk_C = (TH1F*) h_trk_nstub_C ->Clone();
  TH1F* h_nstub_trk_I = (TH1F*) h_trk_nstub_I ->Clone();
  TH1F* h_nstub_trk_F = (TH1F*) h_trk_nstub_F ->Clone();

  TH1F* h_nstub_match_trk   = (TH1F*) h_match_trk_nstub   ->Clone();
  TH1F* h_nstub_match_trk_C = (TH1F*) h_match_trk_nstub_C ->Clone();
  TH1F* h_nstub_match_trk_I = (TH1F*) h_match_trk_nstub_I ->Clone();
  TH1F* h_nstub_match_trk_F = (TH1F*) h_match_trk_nstub_F ->Clone();

  //-----Fill Individual Histograms
  if (doDetailedPlots) {

    // L1 Track nstubs
    h_trk_nstub->SetTitle(htitle);
    h_trk_nstub->SetLineColor(1); //black
    h_trk_nstub->Draw();
    h_trk_nstub->Write();
    //c.SaveAs(DIR+type+"_trk_nstub.png");

      h_trk_nstub_C->SetTitle(htitle);
      h_trk_nstub_C->SetLineColor(1);
      h_trk_nstub_C->Draw();
        sprintf(ctxt,"Central |#eta| < 0.8");
        mySmallText(0.6,0.82,1,ctxt);
      h_trk_nstub_C->Write();
      //c.SaveAs(DIR+type+"_trk_nstub_C.png");

      h_trk_nstub_I->SetTitle(htitle);
      h_trk_nstub_I->SetLineColor(1);
      h_trk_nstub_I->Draw();
        sprintf(ctxt,"Intermediate 0.8 < |#eta| < 1.6");
        mySmallText(0.6,0.82,1,ctxt);
      h_trk_nstub_I->Write();
      //c.SaveAs(DIR+type+"_trk_nstub_I.png");

      h_trk_nstub_F->SetTitle(htitle);
      h_trk_nstub_F->SetLineColor(1);
      h_trk_nstub_F->Draw();
        sprintf(ctxt,"Forward |#eta| > 1.6");
        mySmallText(0.6,0.82,1,ctxt);
      h_trk_nstub_F->Write();
      //c.SaveAs(DIR+type+"_trk_nstub_F.png");

    // TP nstubs
    h_tp_nstub->SetTitle(htitle);
    h_tp_nstub->SetLineColor(4); //blue
    h_tp_nstub->Draw();
    h_tp_nstub->Write();
    //c.SaveAs(DIR+type+"_tp_nstub.png");

      h_tp_nstub_C->SetTitle(htitle);
      h_tp_nstub_C->SetLineColor(4);
      h_tp_nstub_C->Draw();
        sprintf(ctxt,"Central |#eta| < 0.8");
        mySmallText(0.6,0.82,1,ctxt);
      h_tp_nstub_C->Write();
      //c.SaveAs(DIR+type+"_tp_nstub_C.png");

      h_tp_nstub_I->SetTitle(htitle);
      h_tp_nstub_I->SetLineColor(4);
      h_tp_nstub_I->Draw();
        sprintf(ctxt,"Intermediate 0.8 < |#eta| < 1.6");
        mySmallText(0.6,0.82,1,ctxt);
      h_tp_nstub_I->Write();
      //c.SaveAs(DIR+type+"_tp_nstub_I.png");

      h_tp_nstub_F->SetTitle(htitle);
      h_tp_nstub_F->SetLineColor(4);
      h_tp_nstub_F->Draw();
        sprintf(ctxt,"Forward |#eta| > 1.6");
        mySmallText(0.65,0.82,1,ctxt);
      h_tp_nstub_F->Write();
      //c.SaveAs(DIR+type+"_tp_nstub_F.png");

    // Matched Track nstubs
    h_match_trk_nstub->SetTitle(htitle);
    h_match_trk_nstub->SetLineColor(2); //red
    h_match_trk_nstub->Draw();
    h_match_trk_nstub->Write();
    //c.SaveAs(DIR+type+"_match_trk_nstub.png");

      h_match_trk_nstub_C->SetTitle(htitle);
      h_match_trk_nstub_C->SetLineColor(2);
      h_match_trk_nstub_C->Draw();
        sprintf(ctxt,"Central |#eta| < 0.8");
        mySmallText(0.6,0.82,1,ctxt);
      h_match_trk_nstub_C->Write();
      //c.SaveAs(DIR+type+"_match_trk_nstub_C.png");

      h_match_trk_nstub_I->SetTitle(htitle);
      h_match_trk_nstub_I->SetLineColor(2);
      h_match_trk_nstub_I->Draw();
        sprintf(ctxt,"Intermediate 0.8 < |#eta| < 1.6");
        mySmallText(0.6,0.82,1,ctxt);
      h_match_trk_nstub_I->Write();
      //c.SaveAs(DIR+type+"_match_trk_nstub_I.png");

      h_match_trk_nstub_F->SetTitle(htitle);
      h_match_trk_nstub_F->SetLineColor(2);
      h_match_trk_nstub_F->Draw();
        sprintf(ctxt,"Forward |#eta| > 1.6");
        mySmallText(0.6,0.82,1,ctxt);
      h_match_trk_nstub_F->Write();
      //c.SaveAs(DIR+type+"_match_trk_nstub_F.png");

  }

  //-----Stubs Per Event - Multiplot - All eta -----//
  h_nstub_tp->SetTitle(htitle);
  h_nstub_tp->GetYaxis()->SetTitle("Tracks");
  h_nstub_tp->GetXaxis()->SetTitle("Number of Stubs");

  h_nstub_tp->SetLineColor(4); //blue
  h_nstub_tp->SetMarkerColor(4);
  h_nstub_tp->SetLineStyle(1); //solid

  h_nstub_trk->SetLineColor(1); //black
  h_nstub_trk->SetMarkerColor(1);
  h_nstub_trk->SetLineStyle(2); //avg dashes

  h_nstub_match_trk->SetLineColor(2); //red
  h_nstub_match_trk->SetMarkerColor(2);
  h_nstub_match_trk->SetLineStyle(7); //large dashes

  float stubmax = h_nstub_tp->GetMaximum();
  if (h_nstub_trk      ->GetMaximum() > stubmax) stubmax = h_nstub_trk      ->GetMaximum();
  if (h_nstub_match_trk->GetMaximum() > stubmax) stubmax = h_nstub_match_trk->GetMaximum();
  h_nstub_tp->SetAxisRange(0,stubmax*1.05,"Y");

  h_nstub_tp        ->Draw();
  h_nstub_trk       ->Draw("same");
  h_nstub_match_trk ->Draw("same");

  char stubtxt0[500];
  sprintf(stubtxt0,"avg # stubs/track");
  mySmallText(0.6,0.85,1,stubtxt0);

  char stubtxt1[500];
  sprintf(stubtxt1,"TP stubs = %.3f", h_nstub_tp->GetMean());
  mySmallText(0.6,0.75,4,stubtxt1);

  char stubtxt2[500];
  sprintf(stubtxt2,"L1 stubs = %.3f",h_nstub_trk->GetMean());
  mySmallText(0.6,0.80,1,stubtxt2);

  char stubtxt3[500];
  sprintf(stubtxt3,"Matched stubs = %.1f",h_nstub_match_trk->GetMean());
  mySmallText(0.6,0.70,2,stubtxt3);

  
  c.SaveAs(DIR+type+"_nstub.png");
  

  //-----Stubs Per Event - Multiplot - Central eta -----//
  h_nstub_tp_C->SetTitle(htitle);
  h_nstub_tp_C->GetYaxis()->SetTitle("Tracks");
  h_nstub_tp_C->GetXaxis()->SetTitle("Number of Stubs, Central |#eta| < 0.8");

  h_nstub_tp_C->SetLineColor(4); //blue
  h_nstub_tp_C->SetMarkerColor(4);
  h_nstub_tp_C->SetLineStyle(1); //solid

  h_nstub_trk_C->SetLineColor(1); //black
  h_nstub_trk_C->SetMarkerColor(1);
  h_nstub_trk_C->SetLineStyle(2); //avg dashes

  h_nstub_match_trk_C->SetLineColor(2); //red
  h_nstub_match_trk_C->SetMarkerColor(2);
  h_nstub_match_trk_C->SetLineStyle(7); //large dashes

  float stubmax_C = h_nstub_tp_C->GetMaximum();
  if (h_nstub_trk_C      ->GetMaximum() > stubmax_C) stubmax_C = h_nstub_trk_C      ->GetMaximum();
  if (h_nstub_match_trk_C->GetMaximum() > stubmax_C) stubmax_C = h_nstub_match_trk_C->GetMaximum();
  h_nstub_tp_C->SetAxisRange(0,stubmax_C*1.05,"Y");

  h_nstub_tp_C        ->Draw();
  h_nstub_trk_C       ->Draw("same");
  h_nstub_match_trk_C ->Draw("same");

  char stubtxtC0[500];
  sprintf(stubtxtC0,"avg # stubs/track");
  mySmallText(0.6,0.85,1,stubtxtC0);

  char stubtxtC1[500];
  sprintf(stubtxtC1,"TP stubs = %.3f", h_nstub_tp_C->GetMean());
  mySmallText(0.6,0.75,4,stubtxtC1);

  char stubtxtC2[500];
  sprintf(stubtxtC2,"L1 stubs = %.3f",h_nstub_trk_C->GetMean());
  mySmallText(0.6,0.80,1,stubtxtC2);

  char stubtxtC3[500];
  sprintf(stubtxtC3,"Matched stubs = %.1f",h_nstub_match_trk_C->GetMean());
  mySmallText(0.6,0.70,2,stubtxtC3);

  if (doDetailedPlots) {
    c.SaveAs(DIR+type+"_nstub_C.png");
  }

  //-----Stubs Per Event - Multiplot - Intermediate eta -----//
  h_nstub_tp_I->SetTitle(htitle);
  h_nstub_tp_I->GetYaxis()->SetTitle("Tracks");
  h_nstub_tp_I->GetXaxis()->SetTitle("Number of Stubs, Intermediate 0.8 < |#eta| < 1.6");

  h_nstub_tp_I->SetLineColor(4); //blue
  h_nstub_tp_I->SetMarkerColor(4);
  h_nstub_tp_I->SetLineStyle(1); //solid

  h_nstub_trk_I->SetLineColor(1); //black
  h_nstub_trk_I->SetMarkerColor(1);
  h_nstub_trk_I->SetLineStyle(2); //avg dashes

  h_nstub_match_trk_I->SetLineColor(2); //red
  h_nstub_match_trk_I->SetMarkerColor(2);
  h_nstub_match_trk_I->SetLineStyle(7); //large dashes

  float stubmax_I = h_nstub_tp_I->GetMaximum();
  if (h_nstub_trk_I      ->GetMaximum() > stubmax_I) stubmax_I = h_nstub_trk_I      ->GetMaximum();
  if (h_nstub_match_trk_I->GetMaximum() > stubmax_I) stubmax_I = h_nstub_match_trk_I->GetMaximum();
  h_nstub_tp_I->SetAxisRange(0,stubmax_I*1.05,"Y");

  h_nstub_tp_I        ->Draw();
  h_nstub_trk_I       ->Draw("same");
  h_nstub_match_trk_I ->Draw("same");

  char stubtxtI0[500];
  sprintf(stubtxtI0,"avg # stubs/track");
  mySmallText(0.6,0.85,1,stubtxtI0);

  char stubtxtI1[500];
  sprintf(stubtxtI1,"TP stubs = %.3f", h_nstub_tp_I->GetMean());
  mySmallText(0.6,0.75,4,stubtxtI1);

  char stubtxtI2[500];
  sprintf(stubtxtI2,"L1 stubs = %.3f",h_nstub_trk_I->GetMean());
  mySmallText(0.6,0.80,1,stubtxtI2);

  char stubtxtI3[500];
  sprintf(stubtxtI3,"Matched stubs = %.1f",h_nstub_match_trk_I->GetMean());
  mySmallText(0.6,0.70,2,stubtxtI3);

  if (doDetailedPlots) {
    c.SaveAs(DIR+type+"_nstub_I.png");
  }

  //-----Stubs Per Event - Multiplot - Forward eta -----//
  h_nstub_tp_F->SetTitle(htitle);
  h_nstub_tp_F->GetYaxis()->SetTitle("Tracks");
  h_nstub_tp_F->GetXaxis()->SetTitle("Number of Stubs, Forward |#eta| > 1.6");

  h_nstub_tp_F->SetLineColor(4); //blue
  h_nstub_tp_F->SetMarkerColor(4);
  h_nstub_tp_F->SetLineStyle(1); //solid

  h_nstub_trk_F->SetLineColor(1); //black
  h_nstub_trk_F->SetMarkerColor(1);
  h_nstub_trk_F->SetLineStyle(2); //avg dashes

  h_nstub_match_trk_F->SetLineColor(2); //red
  h_nstub_match_trk_F->SetMarkerColor(2);
  h_nstub_match_trk_F->SetLineStyle(7); //large dashes

  float stubmax_F = h_nstub_tp_F->GetMaximum();
  if (h_nstub_trk_F      ->GetMaximum() > stubmax_F) stubmax_F = h_nstub_trk_F      ->GetMaximum();
  if (h_nstub_match_trk_F->GetMaximum() > stubmax_F) stubmax_F = h_nstub_match_trk_F->GetMaximum();
  h_nstub_tp_F->SetAxisRange(0,stubmax_F*1.05,"Y");

  h_nstub_tp_F        ->Draw();
  h_nstub_trk_F       ->Draw("same");
  h_nstub_match_trk_F ->Draw("same");

  char stubtxtF0[500];
  sprintf(stubtxtF0,"avg # stubs/track");
  mySmallText(0.6,0.85,1,stubtxtF0);

  char stubtxtF1[500];
  sprintf(stubtxtF1,"TP stubs = %.3f", h_nstub_tp_F->GetMean());
  mySmallText(0.6,0.75,4,stubtxtF1);

  char stubtxtF2[500];
  sprintf(stubtxtF2,"L1 stubs = %.3f",h_nstub_trk_F->GetMean());
  mySmallText(0.6,0.80,1,stubtxtF2);

  char stubtxtF3[500];
  sprintf(stubtxtF3,"Matched stubs = %.1f",h_nstub_match_trk_F->GetMean());
  mySmallText(0.6,0.70,2,stubtxtF3);

  if (doDetailedPlots) {
    c.SaveAs(DIR+type+"_nstub_F.png");
  }





  // ----------------------------------------------------------------------------------------------------------------
  // More Resolution Plots
  // ----------------------------------------------------------------------------------------------------------------
  float rms = 0;

  // draw and save plots
  h_res_pt->SetTitle(htitle);
  h_res_pt->Draw();
  h_res_pt->Write();
  rms = h_res_pt->GetRMS();
    sprintf(ctxt,"RMS = %.4f",rms);
    mySmallText(0.22,0.82,1,ctxt);
  c.SaveAs(DIR+type+"_res_pt.png");

  h_res_ptRel->SetTitle(htitle);
  h_res_ptRel->Draw();
  h_res_ptRel->Write();
  rms = h_res_ptRel->GetRMS();
    sprintf(ctxt,"RMS = %.4f",rms);
    mySmallText(0.22,0.82,1,ctxt);
  c.SaveAs(DIR+type+"_res_ptRel.png");

  h_res_eta->SetTitle(htitle);
  h_res_eta->Draw();
  h_res_eta->Write();
  rms = h_res_eta->GetRMS();
    sprintf(ctxt,"RMS = %.3e",rms);
    mySmallText(0.22,0.82,1,ctxt);
  c.SaveAs(DIR+type+"_res_eta.png");

  h_res_phi->SetTitle(htitle);
  h_res_phi->Draw();
  h_res_phi->Write();
  rms = h_res_phi->GetRMS();
    sprintf(ctxt,"RMS = %.3e",rms);
    mySmallText(0.22,0.82,1,ctxt);
  c.SaveAs(DIR+type+"_res_phi.png");

  h_res_z0->SetTitle(htitle);
  h_res_z0->Draw();
  h_res_z0->Write();
  rms = h_res_z0->GetRMS();
    sprintf(ctxt,"RMS = %.4f",rms);
    mySmallText(0.22,0.82,1,ctxt);
  c.SaveAs(DIR+type+"_res_z0.png");

  if(doDetailedPlots){
    h_res_z0_C->SetTitle(htitle);
    h_res_z0_C->Draw();
    h_res_z0_C->Write();
    rms = h_res_z0_C->GetRMS();
      sprintf(ctxt,"RMS = %.4f;",rms);
      mySmallText(0.22,0.82,1,ctxt);
      sprintf(ctxt,"|eta| < 0.8");
      mySmallText(0.22,0.76,1,ctxt);
    c.SaveAs(DIR+type+"_res_z0_C.png");

    h_res_z0_I->SetTitle(htitle);
    h_res_z0_I->Draw();
    h_res_z0_I->Write();
    rms = h_res_z0_I->GetRMS();
      sprintf(ctxt,"RMS = %.4f;",rms);
      mySmallText(0.22,0.82,1,ctxt);
      sprintf(ctxt,"0.8 < |eta| < 1.6");
      mySmallText(0.22,0.76,1,ctxt);
    c.SaveAs(DIR+type+"_res_z0_I.png");

    h_res_z0_F->SetTitle(htitle);
    h_res_z0_F->Draw();
    h_res_z0_F->Write();
    rms = h_res_z0_F->GetRMS();
      sprintf(ctxt,"RMS = %.4f;",rms);
      mySmallText(0.22,0.82,1,ctxt);
      sprintf(ctxt,"|eta| > 1.6");
      mySmallText(0.22,0.76,1,ctxt);
    c.SaveAs(DIR+type+"_res_z0_F.png");

    h_res_z0_C_L->SetTitle(htitle);
    h_res_z0_C_L->Draw();
    h_res_z0_C_L->Write();
    rms = h_res_z0_C_L->GetRMS();
      sprintf(ctxt,"RMS = %.4f;",rms);
      mySmallText(0.22,0.82,1,ctxt);
      sprintf(ctxt,"|eta| < 0.8");
      mySmallText(0.22,0.76,1,ctxt);
    c.SaveAs(DIR+type+"_res_z0_C_L.png");

    h_res_z0_I_L->SetTitle(htitle);
    h_res_z0_I_L->Draw();
    h_res_z0_I_L->Write();
    rms = h_res_z0_I_L->GetRMS();
      sprintf(ctxt,"RMS = %.4f;",rms);
      mySmallText(0.22,0.82,1,ctxt);
      sprintf(ctxt,"0.8 < |eta| < 1.6");
      mySmallText(0.22,0.76,1,ctxt);
    c.SaveAs(DIR+type+"_res_z0_I_L.png");

    h_res_z0_F_L->SetTitle(htitle);
    h_res_z0_F_L->Draw();
    h_res_z0_F_L->Write();
    rms = h_res_z0_F_L->GetRMS();
      sprintf(ctxt,"RMS = %.4f;",rms);
      mySmallText(0.22,0.82,1,ctxt);
      sprintf(ctxt,"|eta| > 1.6");
      mySmallText(0.22,0.76,1,ctxt);
    c.SaveAs(DIR+type+"_res_z0_F_L.png");

    h_res_z0_C_H->SetTitle(htitle);
    h_res_z0_C_H->Draw();
    h_res_z0_C_H->Write();
    rms = h_res_z0_C_H->GetRMS();
      sprintf(ctxt,"RMS = %.4f;",rms);
      mySmallText(0.22,0.82,1,ctxt);
      sprintf(ctxt,"|eta| < 0.8");
      mySmallText(0.22,0.76,1,ctxt);
    c.SaveAs(DIR+type+"_res_z0_C_H.png");

    h_res_z0_I_H->SetTitle(htitle);
    h_res_z0_I_H->Draw();
    h_res_z0_I_H->Write();
    rms = h_res_z0_I_H->GetRMS();
      sprintf(ctxt,"RMS = %.4f;",rms);
      mySmallText(0.22,0.82,1,ctxt);
      sprintf(ctxt,"0.8 < |eta| < 1.6");
      mySmallText(0.22,0.76,1,ctxt);
    c.SaveAs(DIR+type+"_res_z0_I_H.png");

    h_res_z0_F_H->SetTitle(htitle);
    h_res_z0_F_H->Draw();
    h_res_z0_F_H->Write();
    rms = h_res_z0_F_H->GetRMS();
      sprintf(ctxt,"RMS = %.4f;",rms);
      mySmallText(0.22,0.82,1,ctxt);
      sprintf(ctxt,"|eta| > 1.6");
      mySmallText(0.22,0.76,1,ctxt);
    c.SaveAs(DIR+type+"_res_z0_F_H.png");

    h_res_z0_L->SetTitle(htitle);
    h_res_z0_L->Draw();
    h_res_z0_L->Write();
    rms = h_res_z0_L->GetRMS();
      sprintf(ctxt,"RMS = %.4f;",rms);
      mySmallText(0.22,0.82,1,ctxt);
      sprintf(ctxt,"p_{T} < 5 GeV");
      mySmallText(0.22,0.76,1,ctxt);
    c.SaveAs(DIR+type+"_res_z0_L.png");

    h_res_z0_H->SetTitle(htitle);
    h_res_z0_H->Draw();
    h_res_z0_H->Write();
    rms = h_res_z0_H->GetRMS();
      sprintf(ctxt,"RMS = %.4f;",rms);
      mySmallText(0.22,0.82,1,ctxt);
      sprintf(ctxt,"p_{T} > 15 GeV");
      mySmallText(0.22,0.76,1,ctxt);
    c.SaveAs(DIR+type+"_res_z0_H.png");
  }

  if (h_res_d0->GetEntries()>0) {
    h_res_d0->SetTitle(htitle);
    h_res_d0->Draw();
    h_res_d0->Write();
    rms = h_res_d0->GetRMS();
      sprintf(ctxt,"RMS = %.4f",rms);
      mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs(DIR+type+"_res_d0.png");

    if(doDetailedPlots){
      h_res_d0_C->SetTitle(htitle);
      h_res_d0_C->Draw();
      h_res_d0_C->Write();
      rms = h_res_d0_C->GetRMS();
        sprintf(ctxt,"RMS = %.4f;",rms);
        mySmallText(0.22,0.82,1,ctxt);
        sprintf(ctxt,"|eta| < 0.8");
        mySmallText(0.22,0.76,1,ctxt);
      c.SaveAs(DIR+type+"_res_d0_C.png");

      h_res_d0_I->SetTitle(htitle);
      h_res_d0_I->Draw();
      h_res_d0_I->Write();
      rms = h_res_d0_I->GetRMS();
        sprintf(ctxt,"RMS = %.4f;",rms);
        mySmallText(0.22,0.82,1,ctxt);
        sprintf(ctxt,"0.8 < |eta| < 1.6");
        mySmallText(0.22,0.76,1,ctxt);
      c.SaveAs(DIR+type+"_res_d0_I.png");

      h_res_d0_F->SetTitle(htitle);
      h_res_d0_F->Draw();
      h_res_d0_F->Write();
      rms = h_res_d0_F->GetRMS();
        sprintf(ctxt,"RMS = %.4f;",rms);
        mySmallText(0.22,0.82,1,ctxt);
        sprintf(ctxt,"|eta| > 1.6");
        mySmallText(0.22,0.76,1,ctxt);
      c.SaveAs(DIR+type+"_res_d0_F.png");

      h_res_d0_C_L->SetTitle(htitle);
      h_res_d0_C_L->Draw();
      h_res_d0_C_L->Write();
      rms = h_res_d0_C_L->GetRMS();
        sprintf(ctxt,"RMS = %.4f;",rms);
        mySmallText(0.22,0.82,1,ctxt);
        sprintf(ctxt,"|eta| < 0.8");
        mySmallText(0.22,0.76,1,ctxt);
      c.SaveAs(DIR+type+"_res_d0_C_L.png");

      h_res_d0_I_L->SetTitle(htitle);
      h_res_d0_I_L->Draw();
      h_res_d0_I_L->Write();
      rms = h_res_d0_I_L->GetRMS();
        sprintf(ctxt,"RMS = %.4f;",rms);
        mySmallText(0.22,0.82,1,ctxt);
        sprintf(ctxt,"0.8 < |eta| < 1.6");
        mySmallText(0.22,0.76,1,ctxt);
      c.SaveAs(DIR+type+"_res_d0_I_L.png");

      h_res_d0_F_L->SetTitle(htitle);
      h_res_d0_F_L->Draw();
      h_res_d0_F_L->Write();
      rms = h_res_d0_F_L->GetRMS();
        sprintf(ctxt,"RMS = %.4f;",rms);
        mySmallText(0.22,0.82,1,ctxt);
        sprintf(ctxt,"|eta| > 1.6");
        mySmallText(0.22,0.76,1,ctxt);
      c.SaveAs(DIR+type+"_res_d0_F_L.png");

      h_res_d0_C_H->SetTitle(htitle);
      h_res_d0_C_H->Draw();
      h_res_d0_C_H->Write();
      rms = h_res_d0_C_H->GetRMS();
        sprintf(ctxt,"RMS = %.4f;",rms);
        mySmallText(0.22,0.82,1,ctxt);
        sprintf(ctxt,"|eta| < 0.8");
        mySmallText(0.22,0.76,1,ctxt);
      c.SaveAs(DIR+type+"_res_d0_C_H.png");

      h_res_d0_I_H->SetTitle(htitle);
      h_res_d0_I_H->Draw();
      h_res_d0_I_H->Write();
      rms = h_res_d0_I_H->GetRMS();
        sprintf(ctxt,"RMS = %.4f;",rms);
        mySmallText(0.22,0.82,1,ctxt);
        sprintf(ctxt,"0.8 < |eta| < 1.6");
        mySmallText(0.22,0.76,1,ctxt);
      c.SaveAs(DIR+type+"_res_d0_I_H.png");

      h_res_d0_F_H->SetTitle(htitle);
      h_res_d0_F_H->Draw();
      h_res_d0_F_H->Write();
      rms = h_res_d0_F_H->GetRMS();
        sprintf(ctxt,"RMS = %.4f;",rms);
        mySmallText(0.22,0.82,1,ctxt);
        sprintf(ctxt,"|eta| > 1.6");
        mySmallText(0.22,0.76,1,ctxt);
      c.SaveAs(DIR+type+"_res_d0_F_H.png");

      h_res_d0_L->SetTitle(htitle);
      h_res_d0_L->Draw();
      h_res_d0_L->Write();
      rms = h_res_d0_L->GetRMS();
        sprintf(ctxt,"RMS = %.4f;",rms);
        mySmallText(0.22,0.82,1,ctxt);
        sprintf(ctxt,"p_{T} < 5 GeV");
        mySmallText(0.22,0.76,1,ctxt);
      c.SaveAs(DIR+type+"_res_d0_L.png");

      h_res_d0_H->SetTitle(htitle);
      h_res_d0_H->Draw();
      h_res_d0_H->Write();
      rms = h_res_d0_H->GetRMS();
        sprintf(ctxt,"RMS = %.4f;",rms);
        mySmallText(0.22,0.82,1,ctxt);
        sprintf(ctxt,"p_{T} > 15 GeV");
        mySmallText(0.22,0.76,1,ctxt);
      c.SaveAs(DIR+type+"_res_d0_H.png");
    }
  }


  fout->Close();



  // ---------------------------------------------------------------------------------------------------------
  // Some Printouts
  // ---------------------------------------------------------------------------------------------------------

  /*
  float k = (float)n_match_eta1p0;
  float N = (float)n_all_eta1p0;
  if (fabs(N)>0) cout << endl << "efficiency for |eta| < 1.0 = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl;
  k = (float)n_match_eta1p75;
  N = (float)n_all_eta1p75;
  if (fabs(N)>0) cout << "efficiency for 1.0 < |eta| < 1.75 = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl;
  k = (float)n_match_eta2p5;
  N = (float)n_all_eta2p5;
  if (fabs(N)>0) cout << "efficiency for 1.75 < |eta| < 2.5 = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl;
  N = (float) n_all_eta1p0 + n_all_eta1p75 + n_all_eta2p5;
  k = (float) n_match_eta1p0 + n_match_eta1p75 + n_match_eta2p5;
  if (fabs(N)>0) cout << "combined efficiency for |eta| < 2.5 = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl << endl;

  // track rates
  cout << "Track Rates for TP and L1 Tracks vs TP pT" << endl;
  cout << "# TP/event (pT > 2.0) =  " << (float)ntp_pt2/nevt << endl;
  cout << "# TP/event (pT > 3.0) =  " << (float)ntp_pt3/nevt << endl;
  cout << "# TP/event (pT > 10.0) = " << (float)ntp_pt10/nevt << endl;

  cout << "# tracks/event (pT > 2.0) =  " << (float)ntrk_pt2/nevt << endl;
  cout << "# tracks/event (pT > 3.0) =  " << (float)ntrk_pt3/nevt << endl;
  cout << "# tracks/event (pT > 10.0) = " << (float)ntrk_pt10/nevt << endl;

  //cout << "# tracks/event (pt > 3.0), 4stubs = " << (float)ntrk_4stub_pt3/nevt << endl;
  //cout << "# tracks/event (pt > 10.0), 4stubs = " << (float)ntrk_4stub_pt10/nevt << endl;
  //cout << "# tracks/event (pt > 3.0), 4stubs+chi2 = " << (float)ntrk_4stubchi2_pt3/nevt << endl;
  //cout << "# tracks/event (pt > 10.0), 4stubs+chi2 = " << (float)ntrk_4stubchi2_pt10/nevt << endl;
  */

}

// ---------------------------------------------------------------------------------------------------------
// Other Functions
// ---------------------------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------------------------
// Plot Style

void SetPlotStyle() {

  // from ATLAS plot style macro

  // use plain black on white colors
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetHistLineColor(1);

  gStyle->SetPalette(1);

  // set the paper & margin sizes
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetLabelFont(42,"xyt");
  gStyle->SetTitleFont(42,"xyt");
  gStyle->SetLabelSize(0.05,"xy");
  gStyle->SetTitleSize(0.05,"xy");
  gStyle->SetLabelSize(0.07,"t");
  gStyle->SetTitleSize(0.07,"t");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]");

  // get rid of error bar caps
  gStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  //gStyle->SetOptTitle(0);   //turns off title
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  // Limit axes values to 4 digits (gets messy otherwise)
  TGaxis::SetMaxDigits(4);

  // Legend Options
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendTextSize(0.04);

}


// ---------------------------------------------------------------------------------------------------------
// Text Overlays on Plots

void mySmallText(Double_t x,Double_t y,Color_t color,char *text) {
  Double_t tsize=0.044;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}


// ---------------------------------------------------------------------------------------------------------
// Resolution Histogram Functions

// Fraction of Entries
double getIntervalContainingFractionOfEntries( TH1* absResidualHistogram, double quantileToCalculate ) {
  
  double totalIntegral = absResidualHistogram->Integral( 0, absResidualHistogram->GetNbinsX() + 1 );

  // Check that the interval is not somewhere in the overflow bin
  double maxAllowedEntriesInOverflow = totalIntegral * ( 1 - quantileToCalculate );
  double nEntriesInOverflow = absResidualHistogram->GetBinContent( absResidualHistogram->GetNbinsX() + 1 );
  if ( nEntriesInOverflow > maxAllowedEntriesInOverflow ) {
        // cout << "WARNING : Cannot compute range corresponding to interval, as it is in the overflow bin" << endl;
        return absResidualHistogram->GetXaxis()->GetXmax()  * 1.2;
  }

  // Calculate quantile for given interval
  double interval[1];
  double quantile[1] = { quantileToCalculate };
  absResidualHistogram->GetQuantiles( 1, interval, quantile);

  return interval[0];
}

// Residual Plot
void makeResidualIntervalPlot( TString type, TString dir, TString variable, TH1F* h_68, TH1F* h_90, TH1F* h_99, double minY, double maxY, TString htitle ) {

  TCanvas c;

  h_68->SetMinimum( minY );
  h_90->SetMinimum( minY );
  h_99->SetMinimum( minY );

  h_68->SetMaximum( maxY );
  h_90->SetMaximum( maxY );
  h_99->SetMaximum( maxY );

  h_68->SetMarkerStyle(20);
  h_90->SetMarkerStyle(26);
  h_99->SetMarkerStyle(24);
  
  h_68->SetTitle(htitle);
  h_68->Draw("P");
  h_68->Write();
  
  h_90->SetTitle(htitle);
  h_90->Draw("P same");
  h_90->Write();
  
  h_99->SetTitle(htitle);
  h_99->Draw("P same");
  h_99->Write();

  TLegend* l = new TLegend(0.65,0.65,0.85,0.85);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetTextSize(0.04);
  l->AddEntry(h_99,"99%","p");
  l->AddEntry(h_90,"90%","p");
  l->AddEntry(h_68,"68%","p");
  l->SetTextFont(42);
  l->Draw();  

  //c.SaveAs(dir+type+"_"+variable+"_interval.png");

  delete l;
}
