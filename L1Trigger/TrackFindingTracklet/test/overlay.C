// -----------------------------------------------------------------------------------------------------------
/*
  HISTOGRAM OVERLAY MACRO
    - Overlay three histograms
    - Use Tracklet L1TrackNtuplePlot output files


    Use by
         root -l
         .L overlay.C++
         overlay("type", file1", "file2", "file3")
      where type is the name of the title of all plots and the prefix of each created png
      where file1, file2, file3 do not include the .root file extension.

    Draw and Save Plot Options
      - Change the output directory (DIR)
      - Mass edit the titles above all plots (htitle)
      - Mass edit the labels in the legends (label1, label2, label3)
      - Mass edit the png filename prefix (overlayType)

    Notes
     - Does not write overlayed plots to root file.
     - d0 plots are only populated if L1Tk_nPar==5 in the NtupleMaker.
     

    In Case of Errors
     - L1 Track plots are not included in the standard plotter. Comment out.
     - Jet eff plots are not included in the standard plotter. Comment out.
     - Chi2 and chi2_dof plots will require write() statements added to the plotter. 
     - Some resolution plots will require write() statements added to the plotter. 
     - ntrk per event plots may require write() statements added to the plotter.
     - Jet sum pT plots are in the standard plotter but may not be written. Check bools and comments.
    
    A supplement to the tracking performance plot script by Louise Skinnari, June 2013
    Created by Rachael Bucci, January 2018 
*/
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

// -----------------------------------------------------------------------------------------------------------
// Main script
// -----------------------------------------------------------------------------------------------------------

void overlay(TString type, TString file1, TString file2, TString file3){
  // type 			        = what type of datasets are being run (e.g. QCD, TTBar, etc.)
  //file1, file2, file3 = file names (excluding .root extension)

  gROOT->SetBatch();
  gErrorIgnoreLevel = kWarning;

  SetPlotStyle();



  // ---------------------------------------------------------------------------------------------------------
  // Histogram Definitions
  // ---------------------------------------------------------------------------------------------------------


  // ---------------------------------------------------------------------------------------------------------
  // File 1

  TH1F* h1_trk_eta = 0;
  TH1F* h1_trk_phi = 0;
  TH1F* h1_trk_pt  = 0; 

  TH1F* h1_tp_eta = 0;
  TH1F* h1_tp_phi = 0;
  TH1F* h1_tp_pt  = 0; 

  TH1F* h1_matchtrk_eta = 0;
  TH1F* h1_matchtrk_phi = 0;
  TH1F* h1_matchtrk_pt  = 0; 

  TH1F* h1_match_tp_eta = 0;
  TH1F* h1_match_tp_phi = 0;
  TH1F* h1_match_tp_pt  = 0; 

  TH1F* h1_eff_L1_eta = 0;
  TH1F* h1_eff_L1_phi = 0;
  TH1F* h1_eff_L1_pt  = 0; 
  TH1F* h1_eff_L1_jetaxis_dR  = 0;
  TH1F* h1_eff_L1_jetaxis_eta = 0;
  TH1F* h1_eff_L1_jetaxis_phi = 0;
  TH1F* h1_eff_L1_jetaxis_pt  = 0;
  TH1F* h1_eff_L1_jetpt = 0;
  
  TH1F* h1_eff_TP_eta = 0;
  TH1F* h1_eff_TP_phi = 0;
  TH1F* h1_eff_TP_pt  = 0;
  TH1F* h1_eff_TP_jetaxis_dR  = 0;
  TH1F* h1_eff_TP_jetaxis_eta = 0;
  TH1F* h1_eff_TP_jetaxis_phi = 0;
  TH1F* h1_eff_TP_jetaxis_pt  = 0;
  TH1F* h1_eff_TP_jetpt   = 0;
  // TH1F* h1_eff_TP_jetpt_C = 0;
  // TH1F* h1_eff_TP_jetpt_I = 0;
  // TH1F* h1_eff_TP_jetpt_F = 0;

  TH1F* h1_res_eta   = 0;
  TH1F* h1_res_phi   = 0;
  TH1F* h1_res_pt    = 0;
  TH1F* h1_res_ptRel = 0;
  TH1F* h1_res_z0    = 0;
  TH1F* h1_res_d0    = 0;

  TH1F* h1_match_trk_chi2     = 0;
  TH1F* h1_match_trk_chi2_dof = 0;

  TH1F* h1_ntrk_pt2  = 0;
  TH1F* h1_ntrk_pt3  = 0;
  TH1F* h1_ntrk_pt10 = 0;

  // TH1F* h1_jet_tp_sumpt_vspt       = 0;
  // TH1F* h1_jet_trk_sumpt_vspt      = 0;
  // TH1F* h1_jet_matchtrk_sumpt_vspt = 0;
  // TH1F* h1_jet_tp_sumpt_vseta       = 0;
  // TH1F* h1_jet_trk_sumpt_vseta      = 0;
  // TH1F* h1_jet_matchtrk_sumpt_vseta = 0;
  // TH1F* h1_jet_tp_sumpt_vsphi       = 0;
  // TH1F* h1_jet_trk_sumpt_vsphi      = 0;
  // TH1F* h1_jet_matchtrk_sumpt_vsphi = 0;
  // TH1F* h1_jet_frac_sumpt_vspt  = 0;
  // TH1F* h1_jet_frac_sumpt_vseta = 0;
  // TH1F* h1_jet_frac_sumpt_vsphi = 0;
  // TH1F* h1_jet_matchfrac_sumpt_vspt  = 0;
  // TH1F* h1_jet_matchfrac_sumpt_vseta = 0;
  // TH1F* h1_jet_matchfrac_sumpt_vsphi = 0;



  // ---------------------------------------------------------------------------------------------------------
  // File 2

  TH1F* h2_trk_eta = 0;
  TH1F* h2_trk_phi = 0;
  TH1F* h2_trk_pt  = 0; 

  TH1F* h2_tp_eta = 0;
  TH1F* h2_tp_phi = 0;
  TH1F* h2_tp_pt  = 0; 

  TH1F* h2_matchtrk_eta = 0;
  TH1F* h2_matchtrk_phi = 0;
  TH1F* h2_matchtrk_pt  = 0; 

  TH1F* h2_match_tp_eta = 0;
  TH1F* h2_match_tp_phi = 0;
  TH1F* h2_match_tp_pt  = 0; 

  TH1F* h2_eff_L1_eta = 0;  
  TH1F* h2_eff_L1_phi = 0;
  TH1F* h2_eff_L1_pt  = 0; 
  TH1F* h2_eff_L1_jetaxis_dR  = 0;
  TH1F* h2_eff_L1_jetaxis_eta = 0;
  TH1F* h2_eff_L1_jetaxis_phi = 0;
  TH1F* h2_eff_L1_jetaxis_pt  = 0;
  TH1F* h2_eff_L1_jetpt = 0;
  
  TH1F* h2_eff_TP_eta = 0;
  TH1F* h2_eff_TP_phi = 0;
  TH1F* h2_eff_TP_pt  = 0;
  TH1F* h2_eff_TP_jetaxis_dR  = 0;
  TH1F* h2_eff_TP_jetaxis_eta = 0;
  TH1F* h2_eff_TP_jetaxis_phi = 0;
  TH1F* h2_eff_TP_jetaxis_pt  = 0;
  TH1F* h2_eff_TP_jetpt   = 0;
  // TH1F* h2_eff_TP_jetpt_C = 0;
  // TH1F* h2_eff_TP_jetpt_I = 0;
  // TH1F* h2_eff_TP_jetpt_F = 0;
  
  TH1F* h2_res_eta   = 0;
  TH1F* h2_res_phi   = 0;
  TH1F* h2_res_pt    = 0;
  TH1F* h2_res_ptRel = 0;
  TH1F* h2_res_z0    = 0;
  TH1F* h2_res_d0    = 0;
  
  TH1F* h2_match_trk_chi2     = 0;
  TH1F* h2_match_trk_chi2_dof = 0;
  
  TH1F* h2_ntrk_pt2  = 0;
  TH1F* h2_ntrk_pt3  = 0;
  TH1F* h2_ntrk_pt10 = 0;

  // TH1F* h2_jet_tp_sumpt_vspt       = 0;
  // TH1F* h2_jet_trk_sumpt_vspt      = 0;
  // TH1F* h2_jet_matchtrk_sumpt_vspt = 0;
  // TH1F* h2_jet_tp_sumpt_vseta       = 0;
  // TH1F* h2_jet_trk_sumpt_vseta      = 0;
  // TH1F* h2_jet_matchtrk_sumpt_vseta = 0;
  // TH1F* h2_jet_tp_sumpt_vsphi       = 0;
  // TH1F* h2_jet_trk_sumpt_vsphi      = 0;
  // TH1F* h2_jet_matchtrk_sumpt_vsphi = 0;
  // TH1F* h2_jet_frac_sumpt_vspt  = 0;
  // TH1F* h2_jet_frac_sumpt_vseta = 0;
  // TH1F* h2_jet_frac_sumpt_vsphi = 0;
  // TH1F* h2_jet_matchfrac_sumpt_vspt = 0;
  // TH1F* h2_jet_matchfrac_sumpt_vseta = 0;
  // TH1F* h2_jet_matchfrac_sumpt_vsphi = 0;


  // ---------------------------------------------------------------------------------------------------------
  // File 3

  TH1F* h3_trk_eta = 0;
  TH1F* h3_trk_phi = 0;
  TH1F* h3_trk_pt  = 0; 

  TH1F* h3_tp_eta = 0;
  TH1F* h3_tp_phi = 0;
  TH1F* h3_tp_pt  = 0; 

  TH1F* h3_matchtrk_eta = 0;
  TH1F* h3_matchtrk_phi = 0;
  TH1F* h3_matchtrk_pt  = 0; 

  TH1F* h3_match_tp_eta = 0;
  TH1F* h3_match_tp_phi = 0;
  TH1F* h3_match_tp_pt  = 0; 

  TH1F* h3_eff_L1_eta = 0;
  TH1F* h3_eff_L1_phi = 0;
  TH1F* h3_eff_L1_pt  = 0; 
  TH1F* h3_eff_L1_jetaxis_dR  = 0;
  TH1F* h3_eff_L1_jetaxis_eta = 0;
  TH1F* h3_eff_L1_jetaxis_phi = 0;
  TH1F* h3_eff_L1_jetaxis_pt  = 0;
  TH1F* h3_eff_L1_jetpt = 0;

  TH1F* h3_eff_TP_eta = 0;
  TH1F* h3_eff_TP_phi = 0;
  TH1F* h3_eff_TP_pt  = 0;
  TH1F* h3_eff_TP_jetaxis_dR  = 0;
  TH1F* h3_eff_TP_jetaxis_eta = 0;
  TH1F* h3_eff_TP_jetaxis_phi = 0;
  TH1F* h3_eff_TP_jetaxis_pt  = 0;
  TH1F* h3_eff_TP_jetpt   = 0;
  // TH1F* h3_eff_TP_jetpt_C = 0;
  // TH1F* h3_eff_TP_jetpt_I = 0;
  // TH1F* h3_eff_TP_jetpt_F = 0;

  TH1F* h3_res_eta   = 0;
  TH1F* h3_res_phi   = 0;
  TH1F* h3_res_pt    = 0;
  TH1F* h3_res_ptRel = 0;
  TH1F* h3_res_z0    = 0;
  TH1F* h3_res_d0    = 0;
  
  TH1F* h3_match_trk_chi2     = 0;
  TH1F* h3_match_trk_chi2_dof = 0;
  
  TH1F* h3_ntrk_pt2  = 0;
  TH1F* h3_ntrk_pt3  = 0;
  TH1F* h3_ntrk_pt10 = 0;

  // TH1F* h3_jet_tp_sumpt_vspt       = 0;
  // TH1F* h3_jet_trk_sumpt_vspt      = 0;
  // TH1F* h3_jet_matchtrk_sumpt_vspt = 0;
  // TH1F* h3_jet_tp_sumpt_vseta       = 0;
  // TH1F* h3_jet_trk_sumpt_vseta      = 0;
  // TH1F* h3_jet_matchtrk_sumpt_vseta = 0;
  // TH1F* h3_jet_tp_sumpt_vsphi       = 0;
  // TH1F* h3_jet_trk_sumpt_vsphi      = 0;
  // TH1F* h3_jet_matchtrk_sumpt_vsphi = 0;
  // TH1F* h3_jet_frac_sumpt_vspt  = 0;
  // TH1F* h3_jet_frac_sumpt_vseta = 0;
  // TH1F* h3_jet_frac_sumpt_vsphi = 0;
  // TH1F* h3_jet_matchfrac_sumpt_vspt  = 0;
  // TH1F* h3_jet_matchfrac_sumpt_vseta = 0;
  // TH1F* h3_jet_matchfrac_sumpt_vsphi = 0;
	


  // ---------------------------------------------------------------------------------------------------------
  // Source Files
  // ---------------------------------------------------------------------------------------------------------
	TString SRCDIR = "~/Tracklet_Emulator/Emulator_932/2018_04_26/";

  TFile f1(SRCDIR+file1+".root");
  TFile f2(SRCDIR+file2+".root");
  TFile f3(SRCDIR+file3+".root");



  // ---------------------------------------------------------------------------------------------------------
  // Get Histograms From Files
  // ---------------------------------------------------------------------------------------------------------
	
  // File 1
  f1.GetObject("trk_eta",             h1_trk_eta);
  f1.GetObject("trk_phi",             h1_trk_phi);
  f1.GetObject("trk_pt",              h1_trk_pt);

  f1.GetObject("tp_eta",              h1_tp_eta);
  f1.GetObject("tp_phi",              h1_tp_phi);
  f1.GetObject("tp_pt",               h1_tp_pt);

  f1.GetObject("matchtrk_eta",        h1_matchtrk_eta);
  f1.GetObject("matchtrk_phi",        h1_matchtrk_phi);
  f1.GetObject("matchtrk_pt",         h1_matchtrk_pt);

  f1.GetObject("match_tp_eta",        h1_match_tp_eta);
  f1.GetObject("match_tp_phi",        h1_match_tp_phi);
  f1.GetObject("match_tp_pt",         h1_match_tp_pt);

  f1.GetObject("eff_L1_eta",          h1_eff_L1_eta);
  f1.GetObject("eff_L1_phi",          h1_eff_L1_phi);
  f1.GetObject("eff_L1_pt",           h1_eff_L1_pt);

  f1.GetObject("eff_L1_jetaxis_dR",   h1_eff_L1_jetaxis_dR);
  f1.GetObject("eff_L1_jetaxis_eta",  h1_eff_L1_jetaxis_eta); 	
  f1.GetObject("eff_L1_jetaxis_phi",  h1_eff_L1_jetaxis_phi);
  f1.GetObject("eff_L1_jetaxis_pt",   h1_eff_L1_jetaxis_pt);
  f1.GetObject("eff_L1_jetpt",        h1_eff_L1_jetpt);

  f1.GetObject("eff_TP_eta",          h1_eff_TP_eta);
  f1.GetObject("eff_TP_phi",          h1_eff_TP_phi);
  f1.GetObject("eff_TP_pt",           h1_eff_TP_pt);

  f1.GetObject("eff_TP_jetaxis_dR",   h1_eff_TP_jetaxis_dR);
  f1.GetObject("eff_TP_jetaxis_eta",  h1_eff_TP_jetaxis_eta);
  f1.GetObject("eff_TP_jetaxis_phi",  h1_eff_TP_jetaxis_phi);
  f1.GetObject("eff_TP_jetaxis_pt",   h1_eff_TP_jetaxis_pt);
  f1.GetObject("eff_TP_jetpt",        h1_eff_TP_jetpt);
  // f1.GetObject("eff_TP_jetpt_C",      h1_eff_TP_jetpt_C);
  // f1.GetObject("eff_TP_jetpt_I",      h1_eff_TP_jetpt_I);
  // f1.GetObject("eff_TP_jetpt_F",      h1_eff_TP_jetpt_F);

  f1.GetObject("res_eta",             h1_res_eta);
  f1.GetObject("res_phi",             h1_res_phi);
  f1.GetObject("res_pt",              h1_res_pt);
  f1.GetObject("res_ptRel",           h1_res_ptRel);
  f1.GetObject("res_z0",              h1_res_z0);
  f1.GetObject("res_d0",              h1_res_d0);

  f1.GetObject("match_trk_chi2",      h1_match_trk_chi2);
  f1.GetObject("match_trk_chi2_dof",  h1_match_trk_chi2_dof);

  f1.GetObject("ntrk_pt2",            h1_ntrk_pt2);
  f1.GetObject("ntrk_pt3",            h1_ntrk_pt3);
  f1.GetObject("ntrk_pt10",           h1_ntrk_pt10);

  // f1.GetObject("jet_tp_sumpt_vspt",           h1_jet_tp_sumpt_vspt);
  // f1.GetObject("jet_trk_sumpt_vspt",          h1_jet_trk_sumpt_vspt);
  // f1.GetObject("jet_matchtrk_sumpt_vspt",     h1_jet_matchtrk_sumpt_vspt);
  // f1.GetObject("jet_tp_sumpt_vseta",          h1_jet_tp_sumpt_vseta);
  // f1.GetObject("jet_trk_sumpt_vseta",         h1_jet_trk_sumpt_vseta);
  // f1.GetObject("jet_matchtrk_sumpt_vseta",    h1_jet_matchtrk_sumpt_vseta);
  // f1.GetObject("jet_tp_sumpt_vsphi",          h1_jet_tp_sumpt_vsphi);
  // f1.GetObject("jet_trk_sumpt_vsphi",         h1_jet_trk_sumpt_vsphi);
  // f1.GetObject("jet_matchtrk_sumpt_vsphi",    h1_jet_matchtrk_sumpt_vsphi);
  // f1.GetObject("jet_frac_sumpt_vspt",         h1_jet_frac_sumpt_vspt);
  // f1.GetObject("jet_frac_sumpt_vseta",        h1_jet_frac_sumpt_vseta);
  // f1.GetObject("jet_frac_sumpt_vsphi",        h1_jet_frac_sumpt_vsphi);
  // f1.GetObject("jet_matchfrac_sumpt_vspt",    h1_jet_matchfrac_sumpt_vspt);
  // f1.GetObject("jet_matchfrac_sumpt_vseta",   h1_jet_matchfrac_sumpt_vseta);
  // f1.GetObject("jet_matchfrac_sumpt_vsphi",   h1_jet_matchfrac_sumpt_vsphi);


  // File 2
  f2.GetObject("trk_eta",             h2_trk_eta);
  f2.GetObject("trk_phi",             h2_trk_phi);
  f2.GetObject("trk_pt",              h2_trk_pt);

  f2.GetObject("tp_eta",              h2_tp_eta);
  f2.GetObject("tp_phi",              h2_tp_phi);
  f2.GetObject("tp_pt",               h2_tp_pt);

  f2.GetObject("matchtrk_eta",        h2_matchtrk_eta);
  f2.GetObject("matchtrk_phi",        h2_matchtrk_phi);
  f2.GetObject("matchtrk_pt",         h2_matchtrk_pt);

  f2.GetObject("match_tp_eta",        h2_match_tp_eta);
  f2.GetObject("match_tp_phi",        h2_match_tp_phi);
  f2.GetObject("match_tp_pt",         h2_match_tp_pt);

  f2.GetObject("eff_L1_eta",          h2_eff_L1_eta);
  f2.GetObject("eff_L1_phi",          h2_eff_L1_phi);
  f2.GetObject("eff_L1_pt",           h2_eff_L1_pt);

  f2.GetObject("eff_L1_jetaxis_dR",   h2_eff_L1_jetaxis_dR);
  f2.GetObject("eff_L1_jetaxis_eta",  h2_eff_L1_jetaxis_eta); 	
  f2.GetObject("eff_L1_jetaxis_phi",  h2_eff_L1_jetaxis_phi);
  f2.GetObject("eff_L1_jetaxis_pt",   h2_eff_L1_jetaxis_pt);
  f2.GetObject("eff_L1_jetpt",        h2_eff_L1_jetpt);

  f2.GetObject("eff_TP_eta",          h2_eff_TP_eta);
  f2.GetObject("eff_TP_phi",          h2_eff_TP_phi);
  f2.GetObject("eff_TP_pt",           h2_eff_TP_pt);

  f2.GetObject("eff_TP_jetaxis_dR",   h2_eff_TP_jetaxis_dR);
  f2.GetObject("eff_TP_jetaxis_eta",  h2_eff_TP_jetaxis_eta);
  f2.GetObject("eff_TP_jetaxis_phi",  h2_eff_TP_jetaxis_phi);
  f2.GetObject("eff_TP_jetaxis_pt",   h2_eff_TP_jetaxis_pt);
  f2.GetObject("eff_TP_jetpt",        h2_eff_TP_jetpt);
  // f2.GetObject("eff_TP_jetpt_C",      h2_eff_TP_jetpt_C);
  // f2.GetObject("eff_TP_jetpt_I",      h2_eff_TP_jetpt_I);
  // f2.GetObject("eff_TP_jetpt_F",      h2_eff_TP_jetpt_F);

  f2.GetObject("res_eta",             h2_res_eta);
  f2.GetObject("res_phi",             h2_res_phi);
  f2.GetObject("res_pt",              h2_res_pt);
  f2.GetObject("res_ptRel",           h2_res_ptRel);
  f2.GetObject("res_z0",              h2_res_z0);
  f2.GetObject("res_d0",              h2_res_d0);

  f2.GetObject("match_trk_chi2",      h2_match_trk_chi2);
  f2.GetObject("match_trk_chi2_dof",  h2_match_trk_chi2_dof);

  f2.GetObject("ntrk_pt2",            h2_ntrk_pt2);
  f2.GetObject("ntrk_pt3",            h2_ntrk_pt3);
  f2.GetObject("ntrk_pt10",           h2_ntrk_pt10);

  // f2.GetObject("jet_tp_sumpt_vspt",           h2_jet_tp_sumpt_vspt);
  // f2.GetObject("jet_trk_sumpt_vspt",          h2_jet_trk_sumpt_vspt);
  // f2.GetObject("jet_matchtrk_sumpt_vspt",     h2_jet_matchtrk_sumpt_vspt);
  // f2.GetObject("jet_tp_sumpt_vseta",          h2_jet_tp_sumpt_vseta);
  // f2.GetObject("jet_trk_sumpt_vseta",         h2_jet_trk_sumpt_vseta);
  // f2.GetObject("jet_matchtrk_sumpt_vseta",    h2_jet_matchtrk_sumpt_vseta);
  // f2.GetObject("jet_tp_sumpt_vsphi",          h2_jet_tp_sumpt_vsphi);
  // f2.GetObject("jet_trk_sumpt_vsphi",         h2_jet_trk_sumpt_vsphi);
  // f2.GetObject("jet_matchtrk_sumpt_vsphi",    h2_jet_matchtrk_sumpt_vsphi);
  // f2.GetObject("jet_frac_sumpt_vspt",         h2_jet_frac_sumpt_vspt);
  // f2.GetObject("jet_frac_sumpt_vseta",        h2_jet_frac_sumpt_vseta);
  // f2.GetObject("jet_frac_sumpt_vsphi",        h2_jet_frac_sumpt_vsphi);
  // f2.GetObject("jet_matchfrac_sumpt_vspt",    h2_jet_matchfrac_sumpt_vspt);
  // f2.GetObject("jet_matchfrac_sumpt_vseta",   h2_jet_matchfrac_sumpt_vseta);
  // f2.GetObject("jet_matchfrac_sumpt_vsphi",   h2_jet_matchfrac_sumpt_vsphi);

  // File 3
  f3.GetObject("trk_eta",             h3_trk_eta);
  f3.GetObject("trk_phi",             h3_trk_phi);
  f3.GetObject("trk_pt",              h3_trk_pt);

  f3.GetObject("tp_eta",              h3_tp_eta);
  f3.GetObject("tp_phi",              h3_tp_phi);
  f3.GetObject("tp_pt",               h3_tp_pt);

  f3.GetObject("matchtrk_eta",        h3_matchtrk_eta);
  f3.GetObject("matchtrk_phi",        h3_matchtrk_phi);
  f3.GetObject("matchtrk_pt",         h3_matchtrk_pt);

  f3.GetObject("match_tp_eta",        h3_match_tp_eta);
  f3.GetObject("match_tp_phi",        h3_match_tp_phi);
  f3.GetObject("match_tp_pt",         h3_match_tp_pt);

  f3.GetObject("eff_L1_eta",          h3_eff_L1_eta);
  f3.GetObject("eff_L1_phi",          h3_eff_L1_phi);
  f3.GetObject("eff_L1_pt",           h3_eff_L1_pt);

  f3.GetObject("eff_L1_jetaxis_dR",   h3_eff_L1_jetaxis_dR);
  f3.GetObject("eff_L1_jetaxis_eta",  h3_eff_L1_jetaxis_eta); 	
  f3.GetObject("eff_L1_jetaxis_phi",  h3_eff_L1_jetaxis_phi);
  f3.GetObject("eff_L1_jetaxis_pt",   h3_eff_L1_jetaxis_pt);
  f3.GetObject("eff_L1_jetpt",        h3_eff_L1_jetpt);

  f3.GetObject("eff_TP_eta",          h3_eff_TP_eta);
  f3.GetObject("eff_TP_phi",          h3_eff_TP_phi);
  f3.GetObject("eff_TP_pt",           h3_eff_TP_pt);

  f3.GetObject("eff_TP_jetaxis_dR",   h3_eff_TP_jetaxis_dR);
  f3.GetObject("eff_TP_jetaxis_eta",  h3_eff_TP_jetaxis_eta);
  f3.GetObject("eff_TP_jetaxis_phi",  h3_eff_TP_jetaxis_phi);
  f3.GetObject("eff_TP_jetaxis_pt",   h3_eff_TP_jetaxis_pt);
  f3.GetObject("eff_TP_jetpt",        h3_eff_TP_jetpt);
  // f3.GetObject("eff_TP_jetpt_C",      h3_eff_TP_jetpt_C);
  // f3.GetObject("eff_TP_jetpt_I",      h3_eff_TP_jetpt_I);
  // f3.GetObject("eff_TP_jetpt_F",      h3_eff_TP_jetpt_F);

  f3.GetObject("res_eta",             h3_res_eta);
  f3.GetObject("res_phi",             h3_res_phi);
  f3.GetObject("res_pt",              h3_res_pt);
  f3.GetObject("res_ptRel",           h3_res_ptRel);
  f3.GetObject("res_z0",              h3_res_z0);
  f3.GetObject("res_d0",              h3_res_d0);

  f3.GetObject("match_trk_chi2",      h3_match_trk_chi2);
  f3.GetObject("match_trk_chi2_dof",  h3_match_trk_chi2_dof);

  f3.GetObject("ntrk_pt2",            h3_ntrk_pt2);
  f3.GetObject("ntrk_pt3",            h3_ntrk_pt3);
  f3.GetObject("ntrk_pt10",           h3_ntrk_pt10);

  // f3.GetObject("jet_tp_sumpt_vspt",           h3_jet_tp_sumpt_vspt);
  // f3.GetObject("jet_trk_sumpt_vspt",          h3_jet_trk_sumpt_vspt);
  // f3.GetObject("jet_matchtrk_sumpt_vspt",     h3_jet_matchtrk_sumpt_vspt);
  // f3.GetObject("jet_tp_sumpt_vseta",          h3_jet_tp_sumpt_vseta);
  // f3.GetObject("jet_trk_sumpt_vseta",         h3_jet_trk_sumpt_vseta);
  // f3.GetObject("jet_matchtrk_sumpt_vseta",    h3_jet_matchtrk_sumpt_vseta);
  // f3.GetObject("jet_tp_sumpt_vsphi",          h3_jet_tp_sumpt_vsphi);
  // f3.GetObject("jet_trk_sumpt_vsphi",         h3_jet_trk_sumpt_vsphi);
  // f3.GetObject("jet_matchtrk_sumpt_vsphi",    h3_jet_matchtrk_sumpt_vsphi);
  // f3.GetObject("jet_frac_sumpt_vspt",         h3_jet_frac_sumpt_vspt);
  // f3.GetObject("jet_frac_sumpt_vseta",        h3_jet_frac_sumpt_vseta);
  // f3.GetObject("jet_frac_sumpt_vsphi",        h3_jet_frac_sumpt_vsphi);
  // f3.GetObject("jet_matchfrac_sumpt_vspt",    h3_jet_matchfrac_sumpt_vspt);
  // f3.GetObject("jet_matchfrac_sumpt_vseta",   h3_jet_matchfrac_sumpt_vseta);
  // f3.GetObject("jet_matchfrac_sumpt_vsphi",   h3_jet_matchfrac_sumpt_vsphi);


  
  // --------------------------------------------------------------------------------------------------------
  // Definitions
  // --------------------------------------------------------------------------------------------------------


  // --------------------------------------------------------------------------------------------------------
  // Resolution Histogram Definitions 

  // Initialize RMS
  float rms_eta1   = 0;
  float rms_eta2   = 0;
  float rms_eta3   = 0;
  float rms_phi1   = 0;
  float rms_phi2   = 0;
  float rms_phi3   = 0;
  float rms_pt1    = 0;
  float rms_pt2    = 0;
  float rms_pt3    = 0;
  float rms_ptRel1 = 0;
  float rms_ptRel2 = 0;
  float rms_ptRel3 = 0;
  float rms_z01    = 0;
  float rms_z02    = 0;
  float rms_z03    = 0;
  float rms_d01    = 0;
  float rms_d02    = 0;
  float rms_d03    = 0;

  // Get RMS
  rms_eta1   = h1_res_eta   ->GetRMS();
  rms_eta2   = h2_res_eta   ->GetRMS();
  rms_eta3   = h3_res_eta   ->GetRMS();
  rms_phi1   = h1_res_phi   ->GetRMS();
  rms_phi2   = h2_res_phi   ->GetRMS();
  rms_phi3   = h3_res_phi   ->GetRMS();
  rms_pt1    = h1_res_pt    ->GetRMS();
  rms_pt2    = h2_res_pt    ->GetRMS();
  rms_pt3    = h3_res_pt    ->GetRMS();
  rms_ptRel1 = h1_res_ptRel ->GetRMS();
  rms_ptRel2 = h2_res_ptRel ->GetRMS();
  rms_ptRel3 = h3_res_ptRel ->GetRMS();
  rms_z01    = h1_res_z0    ->GetRMS();
  rms_z02    = h2_res_z0    ->GetRMS();
  rms_z03    = h3_res_z0    ->GetRMS();
  rms_d01    = h1_res_d0    ->GetRMS();
  rms_d02    = h2_res_d0    ->GetRMS();
  rms_d03    = h3_res_d0    ->GetRMS();


  // --------------------------------------------------------------------------------------------------------
  // Number of Tracks Definitions 

  // Initialize Mean
  float mean_chi2_1      = 0;
  float mean_chi2_2      = 0;
  float mean_chi2_3      = 0;
  float mean_chi2_dof_1  = 0;
  float mean_chi2_dof_2  = 0;
  float mean_chi2_dof_3  = 0;
  float mean_ntrk_pt2_1  = 0;
  float mean_ntrk_pt2_2  = 0;
  float mean_ntrk_pt2_3  = 0;
  float mean_ntrk_pt3_1  = 0;
  float mean_ntrk_pt3_2  = 0;
  float mean_ntrk_pt3_3  = 0;
  float mean_ntrk_pt10_1 = 0;
  float mean_ntrk_pt10_2 = 0;
  float mean_ntrk_pt10_3 = 0;

  // Get Mean
  mean_chi2_1      = h1_match_trk_chi2     ->GetMean();
  mean_chi2_2      = h2_match_trk_chi2     ->GetMean();
  mean_chi2_3      = h3_match_trk_chi2     ->GetMean();
  mean_chi2_dof_1  = h1_match_trk_chi2_dof ->GetMean();
  mean_chi2_dof_2  = h2_match_trk_chi2_dof ->GetMean();
  mean_chi2_dof_3  = h3_match_trk_chi2_dof ->GetMean();
  mean_ntrk_pt2_1  = h1_ntrk_pt2  ->GetMean();
  mean_ntrk_pt2_2  = h2_ntrk_pt2  ->GetMean();
  mean_ntrk_pt2_3  = h3_ntrk_pt2  ->GetMean();
  mean_ntrk_pt3_1  = h1_ntrk_pt3  ->GetMean();
  mean_ntrk_pt3_2  = h2_ntrk_pt3  ->GetMean();
  mean_ntrk_pt3_3  = h3_ntrk_pt3  ->GetMean();
  mean_ntrk_pt10_1 = h1_ntrk_pt10 ->GetMean();
  mean_ntrk_pt10_2 = h2_ntrk_pt10 ->GetMean();
  mean_ntrk_pt10_3 = h3_ntrk_pt10 ->GetMean();


  
  // --------------------------------------------------------------------------------------------------------
  // Draw and Save Plots
  // --------------------------------------------------------------------------------------------------------

  TCanvas c;

  TString DIR = SRCDIR+"overlays/";

  // Histogram Titles
  TString htitle;   // This is the title of all histograms.
  htitle = type;    // The histogram titles are the first parameter of the macro. 
                      // This is set up this way to make it easy to use htitle/type differently 
                      // and to allow changes to individual plot titles
                      // or to set htitle to a file name, etc.

  // Labels for TLegends
  TString label1   = "PU=0";
  TString label2   = "PU=140";
  TString label3   = "PU=200";
  TString overType = ""; //add to file name

	


  // --------------------------------------------------------------------------------------------------------
  // Histograms
  // --------------------------------------------------------------------------------------------------------


  // --------------------------------------------------------------------------------------------------------
  // L1 Plots 

  gPad->SetGridx();
  gPad->SetGridy();

  
  // ----- Event Plots ----- //

  // trk_eta
  h1_trk_eta ->SetLineColor(kBlack);
  h1_trk_eta ->SetLineStyle(1);
  h2_trk_eta ->SetLineColor(kBlue);
  h2_trk_eta ->SetLineStyle(7);
  h3_trk_eta ->SetLineColor(kRed);
  h3_trk_eta ->SetLineStyle(5);
  
  h3_trk_eta ->SetTitle(htitle);
  h3_trk_eta ->Draw();
  h2_trk_eta ->Draw("same");
  h1_trk_eta ->Draw("same");

    TLegend *leg_ev_L1_eta = new TLegend(0.82,0.18,0.94,0.30);
    leg_ev_L1_eta->AddEntry(h1_trk_eta, label1, "l");
    leg_ev_L1_eta->AddEntry(h2_trk_eta, label2, "l");
    leg_ev_L1_eta->AddEntry(h3_trk_eta, label3, "l");
    leg_ev_L1_eta->Draw();

  c.SaveAs(DIR+type+overType+"_ev_L1_eta.png");


  // trk_phi
  h1_trk_phi ->SetLineColor(kBlack);
  h1_trk_phi ->SetLineStyle(1);
  h2_trk_phi ->SetLineColor(kBlue);
  h2_trk_phi ->SetLineStyle(7);
  h3_trk_phi ->SetLineColor(kRed);
  h3_trk_phi ->SetLineStyle(5);
  
  h3_trk_phi ->SetTitle(htitle);
  h3_trk_phi ->Draw();
  h2_trk_phi ->Draw("same");
  h1_trk_phi ->Draw("same");

    TLegend *leg_ev_L1_phi = new TLegend(0.82,0.18,0.94,0.30);
    leg_ev_L1_phi->AddEntry(h1_trk_phi, label1, "l");
    leg_ev_L1_phi->AddEntry(h2_trk_phi, label2, "l");
    leg_ev_L1_phi->AddEntry(h3_trk_phi, label3, "l");
    leg_ev_L1_phi->Draw();

  c.SaveAs(DIR+type+overType+"_ev_L1_phi.png");


  // trk_pt
  h1_trk_pt ->SetLineColor(kBlack);
  h1_trk_pt ->SetLineStyle(1);
  h2_trk_pt ->SetLineColor(kBlue);
  h2_trk_pt ->SetLineStyle(7);
  h3_trk_pt ->SetLineColor(kRed);
  h3_trk_pt ->SetLineStyle(5);
  
  h3_trk_pt ->SetTitle(htitle);
  h3_trk_pt ->Draw();
  h2_trk_pt ->Draw("same");
  h1_trk_pt ->Draw("same");

    TLegend *leg_ev_L1_pt = new TLegend(0.82,0.18,0.94,0.30);
    leg_ev_L1_pt->AddEntry(h1_trk_pt, label1, "l");
    leg_ev_L1_pt->AddEntry(h2_trk_pt, label2, "l");
    leg_ev_L1_pt->AddEntry(h3_trk_pt, label3, "l");
    leg_ev_L1_pt->Draw();

  c.SaveAs(DIR+type+overType+"_ev_L1_pt.png");


  // ----- Purity Plots ----- //

  // eff_L1_eta
  h1_eff_L1_eta ->SetMarkerColor(kBlack);
  h1_eff_L1_eta ->SetMarkerStyle(kCircle);
  h2_eff_L1_eta ->SetMarkerColor(kBlue);
  h2_eff_L1_eta ->SetMarkerStyle(kPlus);
  h3_eff_L1_eta ->SetMarkerColor(kRed);
  h3_eff_L1_eta ->SetMarkerStyle(kMultiply);
  
  h1_eff_L1_eta ->SetTitle(htitle);
  h1_eff_L1_eta ->Draw();
  h2_eff_L1_eta ->Draw("same");
  h3_eff_L1_eta ->Draw("same");

    TLegend *leg_L1_eta = new TLegend(0.82,0.18,0.94,0.30);
    leg_L1_eta->AddEntry(h1_eff_L1_eta, label1, "p");
    leg_L1_eta->AddEntry(h2_eff_L1_eta, label2, "p");
    leg_L1_eta->AddEntry(h3_eff_L1_eta, label3, "p");
    leg_L1_eta->Draw();

  c.SaveAs(DIR+type+overType+"_eff_L1_eta.png");


  // eff_L1_phi
  h1_eff_L1_phi ->SetMarkerColor(kBlack);
  h1_eff_L1_phi ->SetMarkerStyle(kCircle);
  h2_eff_L1_phi ->SetMarkerColor(kBlue);
  h2_eff_L1_phi ->SetMarkerStyle(kPlus);
  h3_eff_L1_phi ->SetMarkerColor(kRed);
  h3_eff_L1_phi ->SetMarkerStyle(kMultiply);
  
  h1_eff_L1_phi ->SetTitle(htitle);
  h1_eff_L1_phi ->Draw();
  h2_eff_L1_phi ->Draw("same");
  h3_eff_L1_phi ->Draw("same");

    TLegend *leg_L1_phi = new TLegend(0.82,0.18,0.94,0.30);
    leg_L1_phi->AddEntry(h1_eff_L1_phi, label1, "p");
    leg_L1_phi->AddEntry(h2_eff_L1_phi, label2, "p");
    leg_L1_phi->AddEntry(h3_eff_L1_phi, label3, "p");
    leg_L1_phi->Draw();
  
  c.SaveAs(DIR+type+overType+"_eff_L1_phi.png");


  // eff_L1_pt
  h1_eff_L1_pt ->SetMarkerColor(kBlack);
  h1_eff_L1_pt ->SetMarkerStyle(kCircle);
  h2_eff_L1_pt ->SetMarkerColor(kBlue);
  h2_eff_L1_pt ->SetMarkerStyle(kPlus);
  h3_eff_L1_pt ->SetMarkerColor(kRed);
  h3_eff_L1_pt ->SetMarkerStyle(kMultiply);
  
  h1_eff_L1_pt ->SetTitle(htitle);
  h1_eff_L1_pt ->Draw();
  h2_eff_L1_pt ->Draw("same");
  h3_eff_L1_pt ->Draw("same");

    TLegend *leg_L1_pt = new TLegend(0.82,0.18,0.94,0.30);
    leg_L1_pt->AddEntry(h1_eff_L1_pt, label1, "p");
    leg_L1_pt->AddEntry(h2_eff_L1_pt, label2, "p");
    leg_L1_pt->AddEntry(h3_eff_L1_pt, label3, "p");
    leg_L1_pt->Draw();
  
  c.SaveAs(DIR+type+overType+"_eff_L1_pt.png");


  // eff_L1_jetaxis_dR
  h1_eff_L1_jetaxis_dR ->SetMarkerColor(kBlack);
  h1_eff_L1_jetaxis_dR ->SetMarkerStyle(kCircle);
  h2_eff_L1_jetaxis_dR ->SetMarkerColor(kBlue);
  h2_eff_L1_jetaxis_dR ->SetMarkerStyle(kPlus);
  h3_eff_L1_jetaxis_dR ->SetMarkerColor(kRed);
  h3_eff_L1_jetaxis_dR ->SetMarkerStyle(kMultiply);
  
  h1_eff_L1_jetaxis_dR ->SetTitle(htitle);
  h1_eff_L1_jetaxis_dR ->Draw();
  h2_eff_L1_jetaxis_dR ->Draw("same");
  h3_eff_L1_jetaxis_dR ->Draw("same");

    TLegend *leg_L1_jetaxis_dR = new TLegend(0.82,0.18,0.94,0.30);
    leg_L1_jetaxis_dR->AddEntry(h1_eff_L1_jetaxis_dR, label1, "p");
    leg_L1_jetaxis_dR->AddEntry(h2_eff_L1_jetaxis_dR, label2, "p");
    leg_L1_jetaxis_dR->AddEntry(h3_eff_L1_jetaxis_dR, label3, "p");
    leg_L1_jetaxis_dR->Draw();
  
  c.SaveAs(DIR+type+overType+"_eff_L1_jetaxis_dR.png");


  // eff_L1_jetaxis_eta
  h1_eff_L1_jetaxis_eta ->SetMarkerColor(kBlack);
  h1_eff_L1_jetaxis_eta ->SetMarkerStyle(kCircle);
  h2_eff_L1_jetaxis_eta ->SetMarkerColor(kBlue);
  h2_eff_L1_jetaxis_eta ->SetMarkerStyle(kPlus);
  h3_eff_L1_jetaxis_eta ->SetMarkerColor(kRed);
  h3_eff_L1_jetaxis_eta ->SetMarkerStyle(kMultiply);
  
  h1_eff_L1_jetaxis_eta ->SetTitle(htitle);
  h1_eff_L1_jetaxis_eta ->Draw();
  h2_eff_L1_jetaxis_eta ->Draw("same");
  h3_eff_L1_jetaxis_eta ->Draw("same");

    TLegend *leg_L1_jetaxis_eta = new TLegend(0.82,0.18,0.94,0.30);
    leg_L1_jetaxis_eta->AddEntry(h1_eff_L1_jetaxis_eta, label1, "p");
    leg_L1_jetaxis_eta->AddEntry(h2_eff_L1_jetaxis_eta, label2, "p");
    leg_L1_jetaxis_eta->AddEntry(h3_eff_L1_jetaxis_eta, label3, "p");
    leg_L1_jetaxis_eta->Draw();
  
  c.SaveAs(DIR+type+overType+"_eff_L1_jetaxis_eta.png");

  
  // eff_L1_jetaxis_phi
  h1_eff_L1_jetaxis_phi ->SetMarkerColor(kBlack);
  h1_eff_L1_jetaxis_phi ->SetMarkerStyle(kCircle);
  h2_eff_L1_jetaxis_phi ->SetMarkerColor(kBlue);
  h2_eff_L1_jetaxis_phi ->SetMarkerStyle(kPlus);
  h3_eff_L1_jetaxis_phi ->SetMarkerColor(kRed);
  h3_eff_L1_jetaxis_phi ->SetMarkerStyle(kMultiply);
  
  h1_eff_L1_jetaxis_phi ->SetTitle(htitle);
  h1_eff_L1_jetaxis_phi ->Draw();
  h2_eff_L1_jetaxis_phi ->Draw("same");
  h3_eff_L1_jetaxis_phi ->Draw("same");

    TLegend *leg_L1_jetaxis_phi = new TLegend(0.82,0.18,0.94,0.30);
    leg_L1_jetaxis_phi->AddEntry(h1_eff_L1_jetaxis_phi, label1, "p");
    leg_L1_jetaxis_phi->AddEntry(h2_eff_L1_jetaxis_phi, label2, "p");
    leg_L1_jetaxis_phi->AddEntry(h3_eff_L1_jetaxis_phi, label3, "p");
    leg_L1_jetaxis_phi->Draw();
  
  c.SaveAs(DIR+type+overType+"_eff_L1_jetaxis_phi.png");


  // eff_L1_jetaxis_pt
  h1_eff_L1_jetaxis_pt ->SetMarkerColor(kBlack);
  h1_eff_L1_jetaxis_pt ->SetMarkerStyle(kCircle);
  h2_eff_L1_jetaxis_pt ->SetMarkerColor(kBlue);
  h2_eff_L1_jetaxis_pt ->SetMarkerStyle(kPlus);
  h3_eff_L1_jetaxis_pt ->SetMarkerColor(kRed);
  h3_eff_L1_jetaxis_pt ->SetMarkerStyle(kMultiply);
  
  h1_eff_L1_jetaxis_pt ->SetTitle(htitle);
  h1_eff_L1_jetaxis_pt ->Draw();
  h2_eff_L1_jetaxis_pt ->Draw("same");
  h3_eff_L1_jetaxis_pt ->Draw("same");

    TLegend *leg_L1_jetaxis_pt = new TLegend(0.82,0.18,0.94,0.30);
    leg_L1_jetaxis_pt->AddEntry(h1_eff_L1_jetaxis_pt, label1, "p");
    leg_L1_jetaxis_pt->AddEntry(h2_eff_L1_jetaxis_pt, label2, "p");
    leg_L1_jetaxis_pt->AddEntry(h3_eff_L1_jetaxis_pt, label3, "p");
    leg_L1_jetaxis_pt->Draw();
  
  c.SaveAs(DIR+type+overType+"_eff_L1_jetaxis_pt.png");


  // eff_L1_jetpt
  h1_eff_L1_jetpt ->SetMarkerColor(kBlack);
  h1_eff_L1_jetpt ->SetMarkerStyle(kCircle);
  h2_eff_L1_jetpt ->SetMarkerColor(kBlue);
  h2_eff_L1_jetpt ->SetMarkerStyle(kPlus);
  h3_eff_L1_jetpt ->SetMarkerColor(kRed);
  h3_eff_L1_jetpt ->SetMarkerStyle(kMultiply);
  
  h1_eff_L1_jetpt ->SetTitle(htitle);
  h1_eff_L1_jetpt ->Draw();
  h2_eff_L1_jetpt ->Draw("same");
  h3_eff_L1_jetpt ->Draw("same");

    TLegend *leg_L1_jetpt = new TLegend(0.82,0.18,0.94,0.30);
    leg_L1_jetpt->AddEntry(h1_eff_L1_jetpt, label1, "p");
    leg_L1_jetpt->AddEntry(h2_eff_L1_jetpt, label2, "p");
    leg_L1_jetpt->AddEntry(h3_eff_L1_jetpt, label3, "p");
    leg_L1_jetpt->Draw();
  
  c.SaveAs(DIR+type+overType+"_eff_L1_jetpt.png");


  // --------------------------------------------------------------------------------------------------------
  // TP Plots

  gPad->SetGridx();
  gPad->SetGridy();


  // ----- Event Plots ----- //

  // ev_TP_eta
  h1_tp_eta ->SetLineColor(kBlack);
  h1_tp_eta ->SetLineStyle(1);
  h2_tp_eta ->SetLineColor(kBlue);
  h2_tp_eta ->SetLineStyle(7);
  h3_tp_eta ->SetLineColor(kRed);
  h3_tp_eta ->SetLineStyle(5);
  
  h3_tp_eta ->SetTitle(htitle);
  h3_tp_eta ->Draw();
  h2_tp_eta ->Draw("same");
  h1_tp_eta ->Draw("same");

    TLegend *leg_ev_TP_eta = new TLegend(0.82,0.18,0.94,0.30);
    leg_ev_TP_eta->AddEntry(h1_tp_eta, label1, "l");
    leg_ev_TP_eta->AddEntry(h2_tp_eta, label2, "l");
    leg_ev_TP_eta->AddEntry(h3_tp_eta, label3, "l");
    leg_ev_TP_eta->Draw();

  c.SaveAs(DIR+type+overType+"_ev_TP_eta.png");


  // ev_TP_phi
  h1_tp_phi ->SetLineColor(kBlack);
  h1_tp_phi ->SetLineStyle(1);
  h2_tp_phi ->SetLineColor(kBlue);
  h2_tp_phi ->SetLineStyle(7);
  h3_tp_phi ->SetLineColor(kRed);
  h3_tp_phi ->SetLineStyle(5);
  
  h3_tp_phi ->SetTitle(htitle);
  h3_tp_phi ->Draw();
  h2_tp_phi ->Draw("same");
  h1_tp_phi ->Draw("same");

    TLegend *leg_ev_TP_phi = new TLegend(0.82,0.18,0.94,0.30);
    leg_ev_TP_phi->AddEntry(h1_tp_phi, label1, "l");
    leg_ev_TP_phi->AddEntry(h2_tp_phi, label2, "l");
    leg_ev_TP_phi->AddEntry(h3_tp_phi, label3, "l");
    leg_ev_TP_phi->Draw();

  c.SaveAs(DIR+type+overType+"_ev_TP_phi.png");


  // ev_TP_pt
  h1_tp_pt ->SetLineColor(kBlack);
  h1_tp_pt ->SetLineStyle(1);
  h2_tp_pt ->SetLineColor(kBlue);
  h2_tp_pt ->SetLineStyle(7);
  h3_tp_pt ->SetLineColor(kRed);
  h3_tp_pt ->SetLineStyle(5);
  
  h3_tp_pt ->SetTitle(htitle);
  h3_tp_pt ->Draw();
  h2_tp_pt ->Draw("same");
  h1_tp_pt ->Draw("same");

    TLegend *leg_ev_TP_pt = new TLegend(0.82,0.18,0.94,0.30);
    leg_ev_TP_pt->AddEntry(h1_tp_pt, label1, "l");
    leg_ev_TP_pt->AddEntry(h2_tp_pt, label2, "l");
    leg_ev_TP_pt->AddEntry(h3_tp_pt, label3, "l");
    leg_ev_TP_pt->Draw();

  c.SaveAs(DIR+type+overType+"_ev_TP_pt.png");



  // ----- Efficiency Plots ----- //

  // eff_TP_eta
  h1_eff_TP_eta ->SetMarkerColor(kBlack);
  h1_eff_TP_eta ->SetMarkerStyle(kCircle);
  h2_eff_TP_eta ->SetMarkerColor(kBlue);
  h2_eff_TP_eta ->SetMarkerStyle(kPlus);
  h3_eff_TP_eta ->SetMarkerColor(kRed);
  h3_eff_TP_eta ->SetMarkerStyle(kMultiply);
  
  h1_eff_TP_eta ->SetTitle(htitle);
  h1_eff_TP_eta ->Draw();
  h2_eff_TP_eta ->Draw("same");
  h3_eff_TP_eta ->Draw("same");

    TLegend *leg_TP_eta = new TLegend(0.82,0.18,0.94,0.30);
    leg_TP_eta->AddEntry(h1_eff_TP_eta, label1, "p");
    leg_TP_eta->AddEntry(h2_eff_TP_eta, label2, "p");
    leg_TP_eta->AddEntry(h3_eff_TP_eta, label3, "p");
    leg_TP_eta->Draw();
  
  c.SaveAs(DIR+type+overType+"_eff_TP_eta.png");


  // eff_TP_phi
  h1_eff_TP_phi ->SetMarkerColor(kBlack);
  h1_eff_TP_phi ->SetMarkerStyle(kCircle);
  h2_eff_TP_phi ->SetMarkerColor(kBlue);
  h2_eff_TP_phi ->SetMarkerStyle(kPlus);
  h3_eff_TP_phi ->SetMarkerColor(kRed);
  h3_eff_TP_phi ->SetMarkerStyle(kMultiply);
  
  h1_eff_TP_phi ->SetTitle(htitle);
  h1_eff_TP_phi ->Draw();
  h2_eff_TP_phi ->Draw("same");
  h3_eff_TP_phi ->Draw("same");

    TLegend *leg_TP_phi = new TLegend(0.82,0.18,0.94,0.30);
    leg_TP_phi->AddEntry(h1_eff_TP_phi, label1, "p");
    leg_TP_phi->AddEntry(h2_eff_TP_phi, label2, "p");
    leg_TP_phi->AddEntry(h3_eff_TP_phi, label3, "p");
    leg_TP_phi->Draw();
  
  c.SaveAs(DIR+type+overType+"_eff_TP_phi.png");


  // eff_TP_pt
  h1_eff_TP_pt ->SetMarkerColor(kBlack);
  h1_eff_TP_pt ->SetMarkerStyle(kCircle);
  h2_eff_TP_pt ->SetMarkerColor(kBlue);
  h2_eff_TP_pt ->SetMarkerStyle(kPlus);
  h3_eff_TP_pt ->SetMarkerColor(kRed);
  h3_eff_TP_pt ->SetMarkerStyle(kMultiply);
  
  h1_eff_TP_pt ->SetTitle(htitle);
  h1_eff_TP_pt ->Draw();
  h2_eff_TP_pt ->Draw("same");
  h3_eff_TP_pt ->Draw("same");

    TLegend *leg_TP_pt = new TLegend(0.82,0.18,0.94,0.30);
    leg_TP_pt->AddEntry(h1_eff_TP_pt, label1, "p");
    leg_TP_pt->AddEntry(h2_eff_TP_pt, label2, "p");
    leg_TP_pt->AddEntry(h3_eff_TP_pt, label3, "p");
    leg_TP_pt->Draw();
  
  c.SaveAs(DIR+type+overType+"_eff_TP_pt.png");


    // eff_TP_jetaxis_dR
    h1_eff_TP_jetaxis_dR ->SetMarkerColor(kBlack);
    h1_eff_TP_jetaxis_dR ->SetMarkerStyle(kCircle);
    h2_eff_TP_jetaxis_dR ->SetMarkerColor(kBlue);
    h2_eff_TP_jetaxis_dR ->SetMarkerStyle(kPlus);
    h3_eff_TP_jetaxis_dR ->SetMarkerColor(kRed);
    h3_eff_TP_jetaxis_dR ->SetMarkerStyle(kMultiply);
    
    h1_eff_TP_jetaxis_dR ->SetTitle(htitle);
    h1_eff_TP_jetaxis_dR ->Draw();
    h2_eff_TP_jetaxis_dR ->Draw("same");
    h3_eff_TP_jetaxis_dR ->Draw("same");

      TLegend *leg_TP_jetaxis_dR = new TLegend(0.82,0.18,0.94,0.30);
      leg_TP_jetaxis_dR->AddEntry(h1_eff_TP_jetaxis_dR, label1, "p");
      leg_TP_jetaxis_dR->AddEntry(h2_eff_TP_jetaxis_dR, label2, "p");
      leg_TP_jetaxis_dR->AddEntry(h3_eff_TP_jetaxis_dR, label3, "p");
      leg_TP_jetaxis_dR->Draw();
    
    c.SaveAs(DIR+type+overType+"_eff_TP_jetaxis_dR.png");


    // eff_TP_jetaxis_eta
    h1_eff_TP_jetaxis_eta ->SetMarkerColor(kBlack);
    h1_eff_TP_jetaxis_eta ->SetMarkerStyle(kCircle);
    h2_eff_TP_jetaxis_eta ->SetMarkerColor(kBlue);
    h2_eff_TP_jetaxis_eta ->SetMarkerStyle(kPlus);
    h3_eff_TP_jetaxis_eta ->SetMarkerColor(kRed);
    h3_eff_TP_jetaxis_eta ->SetMarkerStyle(kMultiply);
    
    h1_eff_TP_jetaxis_eta ->SetTitle(htitle);
    h1_eff_TP_jetaxis_eta ->Draw();
    h2_eff_TP_jetaxis_eta ->Draw("same");
    h3_eff_TP_jetaxis_eta ->Draw("same");

      TLegend *leg_TP_jetaxis_eta = new TLegend(0.82,0.18,0.94,0.30);
      leg_TP_jetaxis_eta->AddEntry(h1_eff_TP_jetaxis_eta, label1, "p");
      leg_TP_jetaxis_eta->AddEntry(h2_eff_TP_jetaxis_eta, label2, "p");
      leg_TP_jetaxis_eta->AddEntry(h3_eff_TP_jetaxis_eta, label3, "p");
      leg_TP_jetaxis_eta->Draw();
    
    c.SaveAs(DIR+type+overType+"_eff_TP_jetaxis_eta.png");

    
    // eff_TP_jetaxis_phi
    h1_eff_TP_jetaxis_phi ->SetMarkerColor(kBlack);
    h1_eff_TP_jetaxis_phi ->SetMarkerStyle(kCircle);
    h2_eff_TP_jetaxis_phi ->SetMarkerColor(kBlue);
    h2_eff_TP_jetaxis_phi ->SetMarkerStyle(kPlus);
    h3_eff_TP_jetaxis_phi ->SetMarkerColor(kRed);
    h3_eff_TP_jetaxis_phi ->SetMarkerStyle(kMultiply);
    
    h1_eff_TP_jetaxis_phi ->SetTitle(htitle);
    h1_eff_TP_jetaxis_phi ->Draw();
    h2_eff_TP_jetaxis_phi ->Draw("same");
    h3_eff_TP_jetaxis_phi ->Draw("same");

      TLegend *leg_TP_jetaxis_phi = new TLegend(0.82,0.18,0.94,0.30);
      leg_TP_jetaxis_phi->AddEntry(h1_eff_TP_jetaxis_phi, label1, "p");
      leg_TP_jetaxis_phi->AddEntry(h2_eff_TP_jetaxis_phi, label2, "p");
      leg_TP_jetaxis_phi->AddEntry(h3_eff_TP_jetaxis_phi, label3, "p");
      leg_TP_jetaxis_phi->Draw();
    
    c.SaveAs(DIR+type+overType+"_eff_TP_jetaxis_phi.png");


    // eff_TP_jetaxis_pt
    h1_eff_TP_jetaxis_pt ->SetMarkerColor(kBlack);
    h1_eff_TP_jetaxis_pt ->SetMarkerStyle(kCircle);
    h2_eff_TP_jetaxis_pt ->SetMarkerColor(kBlue);
    h2_eff_TP_jetaxis_pt ->SetMarkerStyle(kPlus);
    h3_eff_TP_jetaxis_pt ->SetMarkerColor(kRed);
    h3_eff_TP_jetaxis_pt ->SetMarkerStyle(kMultiply);
    
    h1_eff_TP_jetaxis_pt ->SetTitle(htitle);
    h1_eff_TP_jetaxis_pt ->Draw();
    h2_eff_TP_jetaxis_pt ->Draw("same");
    h3_eff_TP_jetaxis_pt ->Draw("same");

      TLegend *leg_TP_jetaxis_pt = new TLegend(0.82,0.18,0.94,0.30);
      leg_TP_jetaxis_pt->AddEntry(h1_eff_TP_jetaxis_pt, label1, "p");
      leg_TP_jetaxis_pt->AddEntry(h2_eff_TP_jetaxis_pt, label2, "p");
      leg_TP_jetaxis_pt->AddEntry(h3_eff_TP_jetaxis_pt, label3, "p");
      leg_TP_jetaxis_pt->Draw();
    
    c.SaveAs(DIR+type+overType+"_eff_TP_jetaxis_pt.png");


    // eff_TP_jetpt
    h1_eff_TP_jetpt ->SetMarkerColor(kBlack);
    h1_eff_TP_jetpt ->SetMarkerStyle(kCircle);
    h2_eff_TP_jetpt ->SetMarkerColor(kBlue);
    h2_eff_TP_jetpt ->SetMarkerStyle(kPlus);
    h3_eff_TP_jetpt ->SetMarkerColor(kRed);
    h3_eff_TP_jetpt ->SetMarkerStyle(kMultiply);
    
    h1_eff_TP_jetpt ->SetTitle(htitle);
    h1_eff_TP_jetpt ->Draw();
    h2_eff_TP_jetpt ->Draw("same");
    h3_eff_TP_jetpt ->Draw("same");

      TLegend *leg_TP_jetpt = new TLegend(0.82,0.18,0.94,0.30);
      leg_TP_jetpt->AddEntry(h1_eff_TP_jetpt, label1, "p");
      leg_TP_jetpt->AddEntry(h2_eff_TP_jetpt, label2, "p");
      leg_TP_jetpt->AddEntry(h3_eff_TP_jetpt, label3, "p");
      leg_TP_jetpt->Draw();
    
    c.SaveAs(DIR+type+overType+"_eff_TP_jetpt.png");

/*
    // eff_TP_jetpt
    h1_eff_TP_jetpt_C ->SetMarkerColor(kBlack);
    h1_eff_TP_jetpt_C ->SetMarkerStyle(kCircle);
    h2_eff_TP_jetpt_C ->SetMarkerColor(kBlue);
    h2_eff_TP_jetpt_C ->SetMarkerStyle(kPlus);
    h3_eff_TP_jetpt_C ->SetMarkerColor(kRed);
    h3_eff_TP_jetpt_C ->SetMarkerStyle(kMultiply);
    
    h1_eff_TP_jetpt_C ->SetTitle(htitle);
    h1_eff_TP_jetpt_C ->Draw();
    h2_eff_TP_jetpt_C ->Draw("same");
    h3_eff_TP_jetpt_C ->Draw("same");

      TLegend *leg_TP_jetpt_C = new TLegend(0.82,0.18,0.94,0.30);
      leg_TP_jetpt_C->AddEntry(h1_eff_TP_jetpt_C, label1, "p");
      leg_TP_jetpt_C->AddEntry(h2_eff_TP_jetpt_C, label2, "p");
      leg_TP_jetpt_C->AddEntry(h3_eff_TP_jetpt_C, label3, "p");
      leg_TP_jetpt_C->Draw();
    
    c.SaveAs(DIR+type+overType+"_eff_TP_jetpt_C.png");


    // eff_TP_jetpt
    h1_eff_TP_jetpt_I ->SetMarkerColor(kBlack);
    h1_eff_TP_jetpt_I ->SetMarkerStyle(kCircle);
    h2_eff_TP_jetpt_I ->SetMarkerColor(kBlue);
    h2_eff_TP_jetpt_I ->SetMarkerStyle(kPlus);
    h3_eff_TP_jetpt_I ->SetMarkerColor(kRed);
    h3_eff_TP_jetpt_I ->SetMarkerStyle(kMultiply);
    
    h1_eff_TP_jetpt_I ->SetTitle(htitle);
    h1_eff_TP_jetpt_I ->Draw();
    h2_eff_TP_jetpt_I ->Draw("same");
    h3_eff_TP_jetpt_I ->Draw("same");

      TLegend *leg_TP_jetpt_I = new TLegend(0.82,0.18,0.94,0.30);
      leg_TP_jetpt_I->AddEntry(h1_eff_TP_jetpt_I, label1, "p");
      leg_TP_jetpt_I->AddEntry(h2_eff_TP_jetpt_I, label2, "p");
      leg_TP_jetpt_I->AddEntry(h3_eff_TP_jetpt_I, label3, "p");
      leg_TP_jetpt_I->Draw();
    
    c.SaveAs(DIR+type+overType+"_eff_TP_jetpt_I.png");


    // eff_TP_jetpt
    h1_eff_TP_jetpt_F ->SetMarkerColor(kBlack);
    h1_eff_TP_jetpt_F ->SetMarkerStyle(kCircle);
    h2_eff_TP_jetpt_F ->SetMarkerColor(kBlue);
    h2_eff_TP_jetpt_F ->SetMarkerStyle(kPlus);
    h3_eff_TP_jetpt_F ->SetMarkerColor(kRed);
    h3_eff_TP_jetpt_F ->SetMarkerStyle(kMultiply);
    
    h1_eff_TP_jetpt_F ->SetTitle(htitle);
    h1_eff_TP_jetpt_F ->Draw();
    h2_eff_TP_jetpt_F ->Draw("same");
    h3_eff_TP_jetpt_F ->Draw("same");

      TLegend *leg_TP_jetpt_F = new TLegend(0.82,0.18,0.94,0.30);
      leg_TP_jetpt_F->AddEntry(h1_eff_TP_jetpt_F, label1, "p");
      leg_TP_jetpt_F->AddEntry(h2_eff_TP_jetpt_F, label2, "p");
      leg_TP_jetpt_F->AddEntry(h3_eff_TP_jetpt_F, label3, "p");
      leg_TP_jetpt_F->Draw();
    
    c.SaveAs(DIR+type+overType+"_eff_TP_jetpt_F.png");
*/

  // --------------------------------------------------------------------------------------------------------
  // Matched Track Events

  // ev_matchtrk_eta
  h1_matchtrk_eta ->SetLineColor(kBlack);
  h1_matchtrk_eta ->SetLineStyle(1);
  h2_matchtrk_eta ->SetLineColor(kBlue);
  h2_matchtrk_eta ->SetLineStyle(7);
  h3_matchtrk_eta ->SetLineColor(kRed);
  h3_matchtrk_eta ->SetLineStyle(5);
  
  h3_matchtrk_eta ->SetTitle(htitle);
  h3_matchtrk_eta ->Draw();
  h2_matchtrk_eta ->Draw("same");
  h1_matchtrk_eta ->Draw("same");

    TLegend *leg_ev_matchtrk_eta = new TLegend(0.82,0.18,0.94,0.30);
    leg_ev_matchtrk_eta->AddEntry(h1_matchtrk_eta, label1, "l");
    leg_ev_matchtrk_eta->AddEntry(h2_matchtrk_eta, label2, "l");
    leg_ev_matchtrk_eta->AddEntry(h3_matchtrk_eta, label3, "l");
    leg_ev_matchtrk_eta->Draw();

  c.SaveAs(DIR+type+overType+"_ev_L1_match_eta.png");


  // ev_matchtrk_phi
  h1_matchtrk_phi ->SetLineColor(kBlack);
  h1_matchtrk_phi ->SetLineStyle(1);
  h2_matchtrk_phi ->SetLineColor(kBlue);
  h2_matchtrk_phi ->SetLineStyle(7);
  h3_matchtrk_phi ->SetLineColor(kRed);
  h3_matchtrk_phi ->SetLineStyle(5);
  
  h3_matchtrk_phi ->SetTitle(htitle);
  h3_matchtrk_phi ->Draw();
  h2_matchtrk_phi ->Draw("same");
  h1_matchtrk_phi ->Draw("same");

    TLegend *leg_ev_matchtrk_phi = new TLegend(0.82,0.18,0.94,0.30);
    leg_ev_matchtrk_phi->AddEntry(h1_matchtrk_phi, label1, "l");
    leg_ev_matchtrk_phi->AddEntry(h2_matchtrk_phi, label2, "l");
    leg_ev_matchtrk_phi->AddEntry(h3_matchtrk_phi, label3, "l");
    leg_ev_matchtrk_phi->Draw();

  c.SaveAs(DIR+type+overType+"_ev_L1_match_phi.png");


  // ev_matchtrk_pt
  h1_matchtrk_pt ->SetLineColor(kBlack);
  h1_matchtrk_pt ->SetLineStyle(1);
  h2_matchtrk_pt ->SetLineColor(kBlue);
  h2_matchtrk_pt ->SetLineStyle(7);
  h3_matchtrk_pt ->SetLineColor(kRed);
  h3_matchtrk_pt ->SetLineStyle(5);
  
  h3_matchtrk_pt ->SetTitle(htitle);
  h3_matchtrk_pt ->Draw();
  h2_matchtrk_pt ->Draw("same");
  h1_matchtrk_pt ->Draw("same");

    TLegend *leg_ev_matchtrk_pt = new TLegend(0.82,0.18,0.94,0.30);
    leg_ev_matchtrk_pt->AddEntry(h1_matchtrk_pt, label1, "l");
    leg_ev_matchtrk_pt->AddEntry(h2_matchtrk_pt, label2, "l");
    leg_ev_matchtrk_pt->AddEntry(h3_matchtrk_pt, label3, "l");
    leg_ev_matchtrk_pt->Draw();

  c.SaveAs(DIR+type+overType+"_ev_L1_match_pt.png");


  // --------------------------------------------------------------------------------------------------------
  // Matched TP Events

  // ev_match_tp_eta
  h1_match_tp_eta ->SetLineColor(kBlack);
  h1_match_tp_eta ->SetLineStyle(1);
  h2_match_tp_eta ->SetLineColor(kBlue);
  h2_match_tp_eta ->SetLineStyle(7);
  h3_match_tp_eta ->SetLineColor(kRed);
  h3_match_tp_eta ->SetLineStyle(5);
  
  h3_match_tp_eta ->SetTitle(htitle);
  h3_match_tp_eta ->Draw();
  h2_match_tp_eta ->Draw("same");
  h1_match_tp_eta ->Draw("same");

    TLegend *leg_ev_match_tp_eta = new TLegend(0.82,0.18,0.94,0.30);
    leg_ev_match_tp_eta->AddEntry(h1_match_tp_eta, label1, "l");
    leg_ev_match_tp_eta->AddEntry(h2_match_tp_eta, label2, "l");
    leg_ev_match_tp_eta->AddEntry(h3_match_tp_eta, label3, "l");
    leg_ev_match_tp_eta->Draw();

  c.SaveAs(DIR+type+overType+"_ev_TP_match_eta.png");


  // ev_match_tp_phi
  h1_match_tp_phi ->SetLineColor(kBlack);
  h1_match_tp_phi ->SetLineStyle(1);
  h2_match_tp_phi ->SetLineColor(kBlue);
  h2_match_tp_phi ->SetLineStyle(7);
  h3_match_tp_phi ->SetLineColor(kRed);
  h3_match_tp_phi ->SetLineStyle(5);
  
  h3_match_tp_phi ->SetTitle(htitle);
  h3_match_tp_phi ->Draw();
  h2_match_tp_phi ->Draw("same");
  h1_match_tp_phi ->Draw("same");

    TLegend *leg_ev_match_tp_phi = new TLegend(0.82,0.18,0.94,0.30);
    leg_ev_match_tp_phi->AddEntry(h1_match_tp_phi, label1, "l");
    leg_ev_match_tp_phi->AddEntry(h2_match_tp_phi, label2, "l");
    leg_ev_match_tp_phi->AddEntry(h3_match_tp_phi, label3, "l");
    leg_ev_match_tp_phi->Draw();

  c.SaveAs(DIR+type+overType+"_ev_TP_match_phi.png");


  // ev_match_tp_pt
  h1_match_tp_pt ->SetLineColor(kBlack);
  h1_match_tp_pt ->SetLineStyle(1);
  h2_match_tp_pt ->SetLineColor(kBlue);
  h2_match_tp_pt ->SetLineStyle(7);
  h3_match_tp_pt ->SetLineColor(kRed);
  h3_match_tp_pt ->SetLineStyle(5);
  
  h3_match_tp_pt ->SetTitle(htitle);
  h3_match_tp_pt ->Draw();
  h2_match_tp_pt ->Draw("same");
  h1_match_tp_pt ->Draw("same");

    TLegend *leg_ev_match_tp_pt = new TLegend(0.82,0.18,0.94,0.30);
    leg_ev_match_tp_pt->AddEntry(h1_match_tp_pt, label1, "l");
    leg_ev_match_tp_pt->AddEntry(h2_match_tp_pt, label2, "l");
    leg_ev_match_tp_pt->AddEntry(h3_match_tp_pt, label3, "l");
    leg_ev_match_tp_pt->Draw();

  c.SaveAs(DIR+type+overType+"_ev_TP_match_pt.png");


  // --------------------------------------------------------------------------------------------------------
  // Resolution Plots

  gPad->SetGridx(0);
  gPad->SetGridy(0);

  // res_eta
  h1_res_eta ->SetLineColor(kBlack);
  h1_res_eta ->SetLineStyle(1);
  h2_res_eta ->SetLineColor(kBlue);
  h2_res_eta ->SetLineStyle(7);
  h3_res_eta ->SetLineColor(kRed);
  h3_res_eta ->SetLineStyle(5);
  
  h1_res_eta ->SetTitle(htitle);
  h1_res_eta ->Draw();
  h2_res_eta ->Draw("same");
  h3_res_eta ->Draw("same");

    TLegend *leg_res_eta = new TLegend(0.73,0.67,0.85,0.87);
    leg_res_eta->AddEntry(h1_res_eta, label1, "l");
    leg_res_eta->AddEntry((TObject*)0, Form("rms=%.3e",rms_eta1), "");
    leg_res_eta->AddEntry(h2_res_eta, label2, "l");
    leg_res_eta->AddEntry((TObject*)0, Form("rms=%.3e",rms_eta2), "");
    leg_res_eta->AddEntry(h3_res_eta, label3, "l");
    leg_res_eta->AddEntry((TObject*)0, Form("rms=%.3e",rms_eta3), "");
    leg_res_eta->Draw();
  
  c.SaveAs(DIR+type+overType+"_res_eta.png");


  // res_phi
  h1_res_phi ->SetLineColor(kBlack);
  h1_res_phi ->SetLineStyle(1);
  h2_res_phi ->SetLineColor(kBlue);
  h2_res_phi ->SetLineStyle(7);
  h3_res_phi ->SetLineColor(kRed);
  h3_res_phi ->SetLineStyle(5);
  
  h1_res_phi ->SetTitle(htitle);
  h1_res_phi ->Draw();
  h2_res_phi ->Draw("same");
  h3_res_phi ->Draw("same");

    TLegend *leg_res_phi = new TLegend(0.73,0.67,0.85,0.87);
    leg_res_phi->AddEntry(h1_res_phi, label1, "l");
    leg_res_phi->AddEntry((TObject*)0, Form("rms=%.3e",rms_phi1), "");
    leg_res_phi->AddEntry(h2_res_phi, label2, "l");
    leg_res_phi->AddEntry((TObject*)0, Form("rms=%.3e",rms_phi2), "");
    leg_res_phi->AddEntry(h3_res_phi, label3, "l");
    leg_res_phi->AddEntry((TObject*)0, Form("rms=%.3e",rms_phi3), "");
    leg_res_phi->Draw();
  
  c.SaveAs(DIR+type+overType+"_res_phi.png");


  // res_pt
  h1_res_pt ->SetLineColor(kBlack);
  h1_res_pt ->SetLineStyle(1);
  h2_res_pt ->SetLineColor(kBlue);
  h2_res_pt ->SetLineStyle(7);
  h3_res_pt ->SetLineColor(kRed);
  h3_res_pt ->SetLineStyle(5);
  
  h1_res_pt ->SetTitle(htitle);
  h1_res_pt ->Draw();
  h2_res_pt ->Draw("same");
  h3_res_pt ->Draw("same");

    TLegend *leg_res_pt = new TLegend(0.73,0.67,0.85,0.87);
    leg_res_pt->AddEntry(h1_res_pt, label1, "l");
    leg_res_pt->AddEntry((TObject*)0, Form("rms=%.4f",rms_pt1), "");
    leg_res_pt->AddEntry(h2_res_pt, label2, "l");
    leg_res_pt->AddEntry((TObject*)0, Form("rms=%.4f",rms_pt2), "");
    leg_res_pt->AddEntry(h3_res_pt, label3, "l");
    leg_res_pt->AddEntry((TObject*)0, Form("rms=%.4f",rms_pt3), "");
    leg_res_pt->Draw();
  
  c.SaveAs(DIR+type+overType+"_res_pt.png");


  // res_ptRel
  h1_res_ptRel ->SetLineColor(kBlack);
  h1_res_ptRel ->SetLineStyle(1);
  h2_res_ptRel ->SetLineColor(kBlue);
  h2_res_ptRel ->SetLineStyle(7);
  h3_res_ptRel ->SetLineColor(kRed);
  h3_res_ptRel ->SetLineStyle(5);
  
  h1_res_ptRel ->SetTitle(htitle);
  h1_res_ptRel ->Draw();
  h2_res_ptRel ->Draw("same");
  h3_res_ptRel ->Draw("same");

    TLegend *leg_res_ptRel = new TLegend(0.73,0.67,0.85,0.87);
    leg_res_ptRel->AddEntry(h1_res_ptRel, label1, "l");
    leg_res_ptRel->AddEntry((TObject*)0, Form("rms=%.4f",rms_ptRel1), "");
    leg_res_ptRel->AddEntry(h2_res_ptRel, label2, "l");
    leg_res_ptRel->AddEntry((TObject*)0, Form("rms=%.4f",rms_ptRel2), "");
    leg_res_ptRel->AddEntry(h3_res_ptRel, label3, "l");
    leg_res_ptRel->AddEntry((TObject*)0, Form("rms=%.4f",rms_ptRel3), "");
    leg_res_ptRel->Draw();
  
  c.SaveAs(DIR+type+overType+"_res_ptRel.png");


  // res_z0
  h1_res_z0 ->SetLineColor(kBlack);
  h1_res_z0 ->SetLineStyle(1);
  h2_res_z0 ->SetLineColor(kBlue);
  h2_res_z0 ->SetLineStyle(7);
  h3_res_z0 ->SetLineColor(kRed);
  h3_res_z0 ->SetLineStyle(5);
  
  h1_res_z0 ->SetTitle(htitle);
  h1_res_z0 ->Draw();
  h2_res_z0 ->Draw("same");
  h3_res_z0 ->Draw("same");

    TLegend *leg_res_z0 = new TLegend(0.73,0.67,0.85,0.87);
    leg_res_z0->AddEntry(h1_res_z0, label1, "l");
    leg_res_z0->AddEntry((TObject*)0, Form("rms=%.4f",rms_z01), "");
    leg_res_z0->AddEntry(h2_res_z0, label2, "l");
    leg_res_z0->AddEntry((TObject*)0, Form("rms=%.4f",rms_z02), "");
    leg_res_z0->AddEntry(h3_res_z0, label3, "l");
    leg_res_z0->AddEntry((TObject*)0, Form("rms=%.4f",rms_z03), "");
    leg_res_z0->Draw();
  
  c.SaveAs(DIR+type+overType+"_res_z0.png");


  // res_d0
  h1_res_d0 ->SetLineColor(kBlack);
  h1_res_d0 ->SetLineStyle(1);
  h2_res_d0 ->SetLineColor(kBlue);
  h2_res_d0 ->SetLineStyle(7);
  h3_res_d0 ->SetLineColor(kRed);
  h3_res_d0 ->SetLineStyle(5);
  
  h1_res_d0 ->SetTitle(htitle);
  h1_res_d0 ->Draw();
  h2_res_d0 ->Draw("same");
  h3_res_d0 ->Draw("same");

    TLegend *leg_res_d0 = new TLegend(0.73,0.67,0.85,0.87);
    leg_res_d0->AddEntry(h1_res_d0, label1, "l");
    leg_res_d0->AddEntry((TObject*)0, Form("rms=%.4f",rms_d01), "");
    leg_res_d0->AddEntry(h2_res_d0, label2, "l");
    leg_res_d0->AddEntry((TObject*)0, Form("rms=%.4f",rms_d02), "");
    leg_res_d0->AddEntry(h3_res_d0, label3, "l");
    leg_res_d0->AddEntry((TObject*)0, Form("rms=%.4f",rms_d03), "");
    leg_res_d0->Draw();
  
  c.SaveAs(DIR+type+overType+"_res_d0.png");



  // --------------------------------------------------------------------------------------------------------
  // Chi2 Distribution Plots

  gPad->SetGridx(0);
  gPad->SetGridy(0);

  // match_trk_chi2
  h1_match_trk_chi2 ->SetLineColor(kBlack);
  h2_match_trk_chi2 ->SetLineStyle(1);
  h2_match_trk_chi2 ->SetLineColor(kBlue);
  h2_match_trk_chi2 ->SetLineStyle(1);
  h3_match_trk_chi2 ->SetLineColor(kRed);
  h3_match_trk_chi2 ->SetLineStyle(1);
  
  h1_match_trk_chi2 ->SetTitle(htitle);
  h1_match_trk_chi2 ->Draw();
  h2_match_trk_chi2 ->Draw("same");
  h3_match_trk_chi2 ->Draw("same");

    TLegend *leg_chi2 = new TLegend(0.73,0.67,0.85,0.87);
    leg_chi2->AddEntry(h1_match_trk_chi2, label1, "l");
    leg_chi2->AddEntry((TObject*)0, Form("mean=%.3f",mean_chi2_1), "");
    leg_chi2->AddEntry(h2_match_trk_chi2, label2, "l");
    leg_chi2->AddEntry((TObject*)0, Form("mean=%.3f",mean_chi2_2), "");
    leg_chi2->AddEntry(h3_match_trk_chi2, label3, "l");
    leg_chi2->AddEntry((TObject*)0, Form("mean=%.3f",mean_chi2_3), "");
    leg_chi2->Draw();
  
  c.SaveAs(DIR+type+overType+"_chi2.png");


  // match_trk_chi2_dof
  h1_match_trk_chi2_dof ->SetLineColor(kBlack);
  h1_match_trk_chi2_dof ->SetLineStyle(1);
  h2_match_trk_chi2_dof ->SetLineColor(kBlue);
  h2_match_trk_chi2_dof ->SetLineStyle(1);
  h3_match_trk_chi2_dof ->SetLineColor(kRed);
  h3_match_trk_chi2_dof ->SetLineStyle(1);
  
  h1_match_trk_chi2_dof ->SetTitle(htitle);
  h1_match_trk_chi2_dof ->Draw();
  h2_match_trk_chi2_dof ->Draw("same");
  h3_match_trk_chi2_dof ->Draw("same");

    TLegend *leg_chi2_dof = new TLegend(0.73,0.67,0.85,0.87);
    leg_chi2_dof->AddEntry(h1_match_trk_chi2_dof, label1, "l");
    leg_chi2_dof->AddEntry((TObject*)0, Form("mean=%.3f",mean_chi2_dof_1), "");
    leg_chi2_dof->AddEntry(h2_match_trk_chi2_dof, label2, "l");
    leg_chi2_dof->AddEntry((TObject*)0, Form("mean=%.3f",mean_chi2_dof_2), "");
    leg_chi2_dof->AddEntry(h3_match_trk_chi2_dof, label3, "l");
    leg_chi2_dof->AddEntry((TObject*)0, Form("mean=%.3f",mean_chi2_dof_3), "");
    leg_chi2_dof->Draw();
  
  c.SaveAs(DIR+type+overType+"_chi2_dof.png");



  // --------------------------------------------------------------------------------------------------------
  // NTrack Plots
  
  gPad->SetGridx(0);
  gPad->SetGridy(0);

  // ntrk_pt2
  h1_ntrk_pt2 ->SetLineColor(kBlack);
  h1_ntrk_pt2 ->SetLineStyle(1);
  h2_ntrk_pt2 ->SetLineColor(kBlue);
  h2_ntrk_pt2 ->SetLineStyle(1);
  h3_ntrk_pt2 ->SetLineColor(kRed);
  h3_ntrk_pt2 ->SetLineStyle(1);
  
  h1_ntrk_pt2 ->SetTitle(htitle);
  h1_ntrk_pt2 ->Draw();
  h2_ntrk_pt2 ->Draw("same");
  h3_ntrk_pt2 ->Draw("same");

    TLegend *leg_ntrk_pt2 = new TLegend(0.73,0.67,0.85,0.87);
    leg_ntrk_pt2->AddEntry(h1_ntrk_pt2, label1, "l");
    leg_ntrk_pt2->AddEntry((TObject*)0, Form("mean=%.3f",mean_ntrk_pt2_1), "");
    leg_ntrk_pt2->AddEntry(h2_ntrk_pt2, label2, "l");
    leg_ntrk_pt2->AddEntry((TObject*)0, Form("mean=%.3f",mean_ntrk_pt2_2), "");
    leg_ntrk_pt2->AddEntry(h3_ntrk_pt2, label3, "l");
    leg_ntrk_pt2->AddEntry((TObject*)0, Form("mean=%.3f",mean_ntrk_pt2_3), "");
    leg_ntrk_pt2->Draw();
  
  c.SaveAs(DIR+type+overType+"_ntrk_pt2.png");


  // ntrk_pt3
  h1_ntrk_pt3 ->SetLineColor(kBlack);
  h1_ntrk_pt3 ->SetLineStyle(1);
  h2_ntrk_pt3 ->SetLineColor(kBlue);
  h2_ntrk_pt3 ->SetLineStyle(1);
  h3_ntrk_pt3 ->SetLineColor(kRed);
  h3_ntrk_pt3 ->SetLineStyle(1);
  
  h1_ntrk_pt3 ->SetTitle(htitle);
  h1_ntrk_pt3 ->Draw();
  h2_ntrk_pt3 ->Draw("same");
  h3_ntrk_pt3 ->Draw("same");

    TLegend *leg_ntrk_pt3 = new TLegend(0.73,0.67,0.85,0.87);
    leg_ntrk_pt3->AddEntry(h1_ntrk_pt2, label1, "l");
    leg_ntrk_pt3->AddEntry((TObject*)0, Form("mean=%.3f",mean_ntrk_pt3_1), "");
    leg_ntrk_pt3->AddEntry(h2_ntrk_pt2, label2, "l");
    leg_ntrk_pt3->AddEntry((TObject*)0, Form("mean=%.3f",mean_ntrk_pt3_2), "");
    leg_ntrk_pt3->AddEntry(h3_ntrk_pt2, label3, "l");
    leg_ntrk_pt3->AddEntry((TObject*)0, Form("mean=%.3f",mean_ntrk_pt3_3), "");
    leg_ntrk_pt3->Draw();
  
  c.SaveAs(DIR+type+overType+"_ntrk_pt3.png");


  // ntrk_pt10
  h1_ntrk_pt10 ->SetLineColor(kBlack);
  h1_ntrk_pt10 ->SetLineStyle(1);
  h2_ntrk_pt10 ->SetLineColor(kBlue);
  h2_ntrk_pt10 ->SetLineStyle(1);
  h3_ntrk_pt10 ->SetLineColor(kRed);
  h3_ntrk_pt10 ->SetLineStyle(1);
  
  h1_ntrk_pt10 ->SetTitle(htitle);
  h1_ntrk_pt10 ->Draw();
  h2_ntrk_pt10 ->Draw("same");
  h3_ntrk_pt10 ->Draw("same");

    TLegend *leg_ntrk_pt10 = new TLegend(0.73,0.67,0.85,0.87);
    leg_ntrk_pt10->AddEntry(h1_ntrk_pt10, label1, "l");
    leg_ntrk_pt10->AddEntry((TObject*)0, Form("mean=%.3f",mean_ntrk_pt10_1), "");
    leg_ntrk_pt10->AddEntry(h2_ntrk_pt10, label2, "l");
    leg_ntrk_pt10->AddEntry((TObject*)0, Form("mean=%.3f",mean_ntrk_pt10_2), "");
    leg_ntrk_pt10->AddEntry(h3_ntrk_pt10, label3, "l");
    leg_ntrk_pt10->AddEntry((TObject*)0, Form("mean=%.3f",mean_ntrk_pt10_3), "");
    leg_ntrk_pt10->Draw();
  
  c.SaveAs(DIR+type+overType+"_ntrk_pt10.png");



  // --------------------------------------------------------------------------------------------------------
  // Jet sum pT Plots
/*
  // jet_tp_sumpt_vspt
  h1_jet_tp_sumpt_vspt ->SetMarkerColor(kBlack);
  h1_jet_tp_sumpt_vspt ->SetMarkerStyle(kCircle);
  h2_jet_tp_sumpt_vspt ->SetMarkerColor(kBlue);
  h2_jet_tp_sumpt_vspt ->SetMarkerStyle(kPlus);
  h3_jet_tp_sumpt_vspt ->SetMarkerColor(kRed);
  h3_jet_tp_sumpt_vspt ->SetMarkerStyle(kMultiply);
  
  h1_jet_tp_sumpt_vspt ->SetTitle(htitle);
  h1_jet_tp_sumpt_vspt ->Draw();
  h2_jet_tp_sumpt_vspt ->Draw("same");
  h3_jet_tp_sumpt_vspt ->Draw("same");

    TLegend *leg_jet_tp_sumpt_vspt = new TLegend(0.82,0.18,0.94,0.30);
    leg_jet_tp_sumpt_vspt->AddEntry(h1_jet_tp_sumpt_vspt, label1, "p");
    leg_jet_tp_sumpt_vspt->AddEntry(h2_jet_tp_sumpt_vspt, label2, "p");
    leg_jet_tp_sumpt_vspt->AddEntry(h3_jet_tp_sumpt_vspt, label3, "p");
    leg_jet_tp_sumpt_vspt->Draw();
  
  c.SaveAs(DIR+type+overType+"_jet_tp_sumpt_vspt.png");

  
  // jet_trk_sumpt_vspt
  h1_jet_trk_sumpt_vspt ->SetMarkerColor(kBlack);
  h1_jet_trk_sumpt_vspt ->SetMarkerStyle(kCircle);
  h2_jet_trk_sumpt_vspt ->SetMarkerColor(kBlue);
  h2_jet_trk_sumpt_vspt ->SetMarkerStyle(kPlus);
  h3_jet_trk_sumpt_vspt ->SetMarkerColor(kRed);
  h3_jet_trk_sumpt_vspt ->SetMarkerStyle(kMultiply);
  
  h1_jet_trk_sumpt_vspt ->SetTitle(htitle);
  h1_jet_trk_sumpt_vspt ->Draw();
  h2_jet_trk_sumpt_vspt ->Draw("same");
  h3_jet_trk_sumpt_vspt ->Draw("same");

    TLegend *leg_jet_trk_sumpt_vspt = new TLegend(0.82,0.18,0.94,0.30);
    leg_jet_trk_sumpt_vspt->AddEntry(h1_jet_trk_sumpt_vspt, label1, "p");
    leg_jet_trk_sumpt_vspt->AddEntry(h2_jet_trk_sumpt_vspt, label2, "p");
    leg_jet_trk_sumpt_vspt->AddEntry(h3_jet_trk_sumpt_vspt, label3, "p");
    leg_jet_trk_sumpt_vspt->Draw();
  
  c.SaveAs(DIR+type+overType+"_jet_trk_sumpt_vspt.png");


  // jet_matchtrk_sumpt_vspt
  h1_jet_matchtrk_sumpt_vspt ->SetMarkerColor(kBlack);
  h1_jet_matchtrk_sumpt_vspt ->SetMarkerStyle(kCircle);
  h2_jet_matchtrk_sumpt_vspt ->SetMarkerColor(kBlue);
  h2_jet_matchtrk_sumpt_vspt ->SetMarkerStyle(kPlus);
  h3_jet_matchtrk_sumpt_vspt ->SetMarkerColor(kRed);
  h3_jet_matchtrk_sumpt_vspt ->SetMarkerStyle(kMultiply);
  
  h1_jet_matchtrk_sumpt_vspt ->SetTitle(htitle);
  h1_jet_matchtrk_sumpt_vspt ->Draw();
  h2_jet_matchtrk_sumpt_vspt ->Draw("same");
  h3_jet_matchtrk_sumpt_vspt ->Draw("same");

    TLegend *leg_jet_matchtrk_sumpt_vspt = new TLegend(0.82,0.18,0.94,0.30);
    leg_jet_matchtrk_sumpt_vspt->AddEntry(h1_jet_matchtrk_sumpt_vspt, label1, "p");
    leg_jet_matchtrk_sumpt_vspt->AddEntry(h2_jet_matchtrk_sumpt_vspt, label2, "p");
    leg_jet_matchtrk_sumpt_vspt->AddEntry(h3_jet_matchtrk_sumpt_vspt, label3, "p");
    leg_jet_matchtrk_sumpt_vspt->Draw();
  
  c.SaveAs(DIR+type+overType+"_jet_matchtrk_sumpt_vspt.png");


  // jet_tp_sumpt_vseta
  h1_jet_tp_sumpt_vseta ->SetMarkerColor(kBlack);
  h1_jet_tp_sumpt_vseta ->SetMarkerStyle(kCircle);
  h2_jet_tp_sumpt_vseta ->SetMarkerColor(kBlue);
  h2_jet_tp_sumpt_vseta ->SetMarkerStyle(kPlus);
  h3_jet_tp_sumpt_vseta ->SetMarkerColor(kRed);
  h3_jet_tp_sumpt_vseta ->SetMarkerStyle(kMultiply);
  
  h1_jet_tp_sumpt_vseta ->SetTitle(htitle);
  h1_jet_tp_sumpt_vseta ->Draw();
  h2_jet_tp_sumpt_vseta ->Draw("same");
  h3_jet_tp_sumpt_vseta ->Draw("same");

    TLegend *leg_jet_tp_sumpt_vseta = new TLegend(0.82,0.18,0.94,0.30);
    leg_jet_tp_sumpt_vseta->AddEntry(h1_jet_tp_sumpt_vseta, label1, "p");
    leg_jet_tp_sumpt_vseta->AddEntry(h2_jet_tp_sumpt_vseta, label2, "p");
    leg_jet_tp_sumpt_vseta->AddEntry(h3_jet_tp_sumpt_vseta, label3, "p");
    leg_jet_tp_sumpt_vseta->Draw();
  
  c.SaveAs(DIR+type+overType+"_jet_tp_sumpt_vseta.png");

  
  // jet_trk_sumpt_vseta
  h1_jet_trk_sumpt_vseta ->SetMarkerColor(kBlack);
  h1_jet_trk_sumpt_vseta ->SetMarkerStyle(kCircle);
  h2_jet_trk_sumpt_vseta ->SetMarkerColor(kBlue);
  h2_jet_trk_sumpt_vseta ->SetMarkerStyle(kPlus);
  h3_jet_trk_sumpt_vseta ->SetMarkerColor(kRed);
  h3_jet_trk_sumpt_vseta ->SetMarkerStyle(kMultiply);
  
  h1_jet_trk_sumpt_vseta ->SetTitle(htitle);
  h1_jet_trk_sumpt_vseta ->Draw();
  h2_jet_trk_sumpt_vseta ->Draw("same");
  h3_jet_trk_sumpt_vseta ->Draw("same");

    TLegend *leg_jet_trk_sumpt_vseta = new TLegend(0.82,0.18,0.94,0.30);
    leg_jet_trk_sumpt_vseta->AddEntry(h1_jet_trk_sumpt_vseta, label1, "p");
    leg_jet_trk_sumpt_vseta->AddEntry(h2_jet_trk_sumpt_vseta, label2, "p");
    leg_jet_trk_sumpt_vseta->AddEntry(h3_jet_trk_sumpt_vseta, label3, "p");
    leg_jet_trk_sumpt_vseta->Draw();
  
  c.SaveAs(DIR+type+overType+"_jet_trk_sumpt_vseta.png");


  // jet_matchtrk_sumpt_vseta
  h1_jet_matchtrk_sumpt_vseta ->SetMarkerColor(kBlack);
  h1_jet_matchtrk_sumpt_vseta ->SetMarkerStyle(kCircle);
  h2_jet_matchtrk_sumpt_vseta ->SetMarkerColor(kBlue);
  h2_jet_matchtrk_sumpt_vseta ->SetMarkerStyle(kPlus);
  h3_jet_matchtrk_sumpt_vseta ->SetMarkerColor(kRed);
  h3_jet_matchtrk_sumpt_vseta ->SetMarkerStyle(kMultiply);
  
  h1_jet_matchtrk_sumpt_vseta ->SetTitle(htitle);
  h1_jet_matchtrk_sumpt_vseta ->Draw();
  h2_jet_matchtrk_sumpt_vseta ->Draw("same");
  h3_jet_matchtrk_sumpt_vseta ->Draw("same");

    TLegend *leg_jet_matchtrk_sumpt_vseta = new TLegend(0.82,0.18,0.94,0.30);
    leg_jet_matchtrk_sumpt_vseta->AddEntry(h1_jet_matchtrk_sumpt_vseta, label1, "p");
    leg_jet_matchtrk_sumpt_vseta->AddEntry(h2_jet_matchtrk_sumpt_vseta, label2, "p");
    leg_jet_matchtrk_sumpt_vseta->AddEntry(h3_jet_matchtrk_sumpt_vseta, label3, "p");
    leg_jet_matchtrk_sumpt_vseta->Draw();
  
  c.SaveAs(DIR+type+overType+"_jet_matchtrk_sumpt_vseta.png");


  // jet_tp_sumpt_vsphi
  h1_jet_tp_sumpt_vsphi ->SetMarkerColor(kBlack);
  h1_jet_tp_sumpt_vsphi ->SetMarkerStyle(kCircle);
  h2_jet_tp_sumpt_vsphi ->SetMarkerColor(kBlue);
  h2_jet_tp_sumpt_vsphi ->SetMarkerStyle(kPlus);
  h3_jet_tp_sumpt_vsphi ->SetMarkerColor(kRed);
  h3_jet_tp_sumpt_vsphi ->SetMarkerStyle(kMultiply);
  
  h1_jet_tp_sumpt_vsphi ->SetTitle(htitle);
  h1_jet_tp_sumpt_vsphi ->Draw();
  h2_jet_tp_sumpt_vsphi ->Draw("same");
  h3_jet_tp_sumpt_vsphi ->Draw("same");

    TLegend *leg_jet_tp_sumpt_vsphi = new TLegend(0.82,0.18,0.94,0.30);
    leg_jet_tp_sumpt_vsphi->AddEntry(h1_jet_tp_sumpt_vsphi, label1, "p");
    leg_jet_tp_sumpt_vsphi->AddEntry(h2_jet_tp_sumpt_vsphi, label2, "p");
    leg_jet_tp_sumpt_vsphi->AddEntry(h3_jet_tp_sumpt_vsphi, label3, "p");
    leg_jet_tp_sumpt_vsphi->Draw();
  
  c.SaveAs(DIR+type+overType+"_jet_tp_sumpt_vsphi.png");

  
  // jet_trk_sumpt_vsphi
  h1_jet_trk_sumpt_vsphi ->SetMarkerColor(kBlack);
  h1_jet_trk_sumpt_vsphi ->SetMarkerStyle(kCircle);
  h2_jet_trk_sumpt_vsphi ->SetMarkerColor(kBlue);
  h2_jet_trk_sumpt_vsphi ->SetMarkerStyle(kPlus);
  h3_jet_trk_sumpt_vsphi ->SetMarkerColor(kRed);
  h3_jet_trk_sumpt_vsphi ->SetMarkerStyle(kMultiply);
  
  h1_jet_trk_sumpt_vsphi ->SetTitle(htitle);
  h1_jet_trk_sumpt_vsphi ->Draw();
  h2_jet_trk_sumpt_vsphi ->Draw("same");
  h3_jet_trk_sumpt_vsphi ->Draw("same");

    TLegend *leg_jet_trk_sumpt_vsphi = new TLegend(0.82,0.18,0.94,0.30);
    leg_jet_trk_sumpt_vsphi->AddEntry(h1_jet_trk_sumpt_vsphi, label1, "p");
    leg_jet_trk_sumpt_vsphi->AddEntry(h2_jet_trk_sumpt_vsphi, label2, "p");
    leg_jet_trk_sumpt_vsphi->AddEntry(h3_jet_trk_sumpt_vsphi, label3, "p");
    leg_jet_trk_sumpt_vsphi->Draw();
  
  c.SaveAs(DIR+type+overType+"_jet_trk_sumpt_vsphi.png");


  // jet_matchtrk_sumpt_vsphi
  h1_jet_matchtrk_sumpt_vsphi ->SetMarkerColor(kBlack);
  h1_jet_matchtrk_sumpt_vsphi ->SetMarkerStyle(kCircle);
  h2_jet_matchtrk_sumpt_vsphi ->SetMarkerColor(kBlue);
  h2_jet_matchtrk_sumpt_vsphi ->SetMarkerStyle(kPlus);
  h3_jet_matchtrk_sumpt_vsphi ->SetMarkerColor(kRed);
  h3_jet_matchtrk_sumpt_vsphi ->SetMarkerStyle(kMultiply);
  
  h1_jet_matchtrk_sumpt_vsphi ->SetTitle(htitle);
  h1_jet_matchtrk_sumpt_vsphi ->Draw();
  h2_jet_matchtrk_sumpt_vsphi ->Draw("same");
  h3_jet_matchtrk_sumpt_vsphi ->Draw("same");

    TLegend *leg_jet_matchtrk_sumpt_vsphi = new TLegend(0.82,0.18,0.94,0.30);
    leg_jet_matchtrk_sumpt_vsphi->AddEntry(h1_jet_matchtrk_sumpt_vsphi, label1, "p");
    leg_jet_matchtrk_sumpt_vsphi->AddEntry(h2_jet_matchtrk_sumpt_vsphi, label2, "p");
    leg_jet_matchtrk_sumpt_vsphi->AddEntry(h3_jet_matchtrk_sumpt_vsphi, label3, "p");
    leg_jet_matchtrk_sumpt_vsphi->Draw();
  
  c.SaveAs(DIR+type+overType+"_jet_matchtrk_sumpt_vsphi.png");


  // jet_frac_sumpt_vspt
  h1_jet_frac_sumpt_vspt ->SetMarkerColor(kBlack);
  h1_jet_frac_sumpt_vspt ->SetMarkerStyle(kCircle);
  h2_jet_frac_sumpt_vspt ->SetMarkerColor(kBlue);
  h2_jet_frac_sumpt_vspt ->SetMarkerStyle(kPlus);
  h3_jet_frac_sumpt_vspt ->SetMarkerColor(kRed);
  h3_jet_frac_sumpt_vspt ->SetMarkerStyle(kMultiply);
  
  h1_jet_frac_sumpt_vspt ->SetTitle(htitle);
  h1_jet_frac_sumpt_vspt ->Draw();
  h2_jet_frac_sumpt_vspt ->Draw("same");
  h3_jet_frac_sumpt_vspt ->Draw("same");

    TLegend *leg_jet_frac_sumpt_vspt = new TLegend(0.82,0.18,0.94,0.30);
    leg_jet_frac_sumpt_vspt->AddEntry(h1_jet_frac_sumpt_vspt, label1, "p");
    leg_jet_frac_sumpt_vspt->AddEntry(h2_jet_frac_sumpt_vspt, label2, "p");
    leg_jet_frac_sumpt_vspt->AddEntry(h3_jet_frac_sumpt_vspt, label3, "p");
    leg_jet_frac_sumpt_vspt->Draw();
  
  c.SaveAs(DIR+type+overType+"_jet_frac_sumpt_vspt.png");


  // jet_frac_sumpt_vseta
  h1_jet_frac_sumpt_vseta ->SetMarkerColor(kBlack);
  h1_jet_frac_sumpt_vseta ->SetMarkerStyle(kCircle);
  h2_jet_frac_sumpt_vseta ->SetMarkerColor(kBlue);
  h2_jet_frac_sumpt_vseta ->SetMarkerStyle(kPlus);
  h3_jet_frac_sumpt_vseta ->SetMarkerColor(kRed);
  h3_jet_frac_sumpt_vseta ->SetMarkerStyle(kMultiply);
  
  h1_jet_frac_sumpt_vseta ->SetTitle(htitle);
  h1_jet_frac_sumpt_vseta ->Draw();
  h2_jet_frac_sumpt_vseta ->Draw("same");
  h3_jet_frac_sumpt_vseta ->Draw("same");

    TLegend *leg_jet_frac_sumpt_vseta = new TLegend(0.82,0.18,0.94,0.30);
    leg_jet_frac_sumpt_vseta->AddEntry(h1_jet_frac_sumpt_vseta, label1, "p");
    leg_jet_frac_sumpt_vseta->AddEntry(h2_jet_frac_sumpt_vseta, label2, "p");
    leg_jet_frac_sumpt_vseta->AddEntry(h3_jet_frac_sumpt_vseta, label3, "p");
    leg_jet_frac_sumpt_vseta->Draw();
  
  c.SaveAs(DIR+type+overType+"_jet_frac_sumpt_vseta.png");


  // jet_frac_sumpt_vsphi
  h1_jet_frac_sumpt_vsphi ->SetMarkerColor(kBlack);
  h1_jet_frac_sumpt_vsphi ->SetMarkerStyle(kCircle);
  h2_jet_frac_sumpt_vsphi ->SetMarkerColor(kBlue);
  h2_jet_frac_sumpt_vsphi ->SetMarkerStyle(kPlus);
  h3_jet_frac_sumpt_vsphi ->SetMarkerColor(kRed);
  h3_jet_frac_sumpt_vsphi ->SetMarkerStyle(kMultiply);
  
  h1_jet_frac_sumpt_vsphi ->SetTitle(htitle);
  h1_jet_frac_sumpt_vsphi ->Draw();
  h2_jet_frac_sumpt_vsphi ->Draw("same");
  h3_jet_frac_sumpt_vsphi ->Draw("same");

    TLegend *leg_jet_frac_sumpt_vsphi = new TLegend(0.82,0.18,0.94,0.30);
    leg_jet_frac_sumpt_vsphi->AddEntry(h1_jet_frac_sumpt_vsphi, label1, "p");
    leg_jet_frac_sumpt_vsphi->AddEntry(h2_jet_frac_sumpt_vsphi, label2, "p");
    leg_jet_frac_sumpt_vsphi->AddEntry(h3_jet_frac_sumpt_vsphi, label3, "p");
    leg_jet_frac_sumpt_vsphi->Draw();
  
  c.SaveAs(DIR+type+overType+"_jet_frac_sumpt_vsphi.png");


  // jet_matchfrac_sumpt_vspt
  h1_jet_matchfrac_sumpt_vspt ->SetMarkerColor(kBlack);
  h1_jet_matchfrac_sumpt_vspt ->SetMarkerStyle(kCircle);
  h2_jet_matchfrac_sumpt_vspt ->SetMarkerColor(kBlue);
  h2_jet_matchfrac_sumpt_vspt ->SetMarkerStyle(kPlus);
  h3_jet_matchfrac_sumpt_vspt ->SetMarkerColor(kRed);
  h3_jet_matchfrac_sumpt_vspt ->SetMarkerStyle(kMultiply);
  
  h1_jet_matchfrac_sumpt_vspt ->SetTitle(htitle);
  h1_jet_matchfrac_sumpt_vspt ->Draw();
  h2_jet_matchfrac_sumpt_vspt ->Draw("same");
  h3_jet_matchfrac_sumpt_vspt ->Draw("same");

    TLegend *leg_jet_matchfrac_sumpt_vspt = new TLegend(0.82,0.18,0.94,0.30);
    leg_jet_matchfrac_sumpt_vspt->AddEntry(h1_jet_matchfrac_sumpt_vspt, label1, "p");
    leg_jet_matchfrac_sumpt_vspt->AddEntry(h2_jet_matchfrac_sumpt_vspt, label2, "p");
    leg_jet_matchfrac_sumpt_vspt->AddEntry(h3_jet_matchfrac_sumpt_vspt, label3, "p");
    leg_jet_matchfrac_sumpt_vspt->Draw();
  
  c.SaveAs(DIR+type+overType+"_jet_matchfrac_sumpt_vspt.png");


  // jet_matchfrac_sumpt_vseta
  h1_jet_matchfrac_sumpt_vseta ->SetMarkerColor(kBlack);
  h1_jet_matchfrac_sumpt_vseta ->SetMarkerStyle(kCircle);
  h2_jet_matchfrac_sumpt_vseta ->SetMarkerColor(kBlue);
  h2_jet_matchfrac_sumpt_vseta ->SetMarkerStyle(kPlus);
  h3_jet_matchfrac_sumpt_vseta ->SetMarkerColor(kRed);
  h3_jet_matchfrac_sumpt_vseta ->SetMarkerStyle(kMultiply);
  
  h1_jet_matchfrac_sumpt_vseta ->SetTitle(htitle);
  h1_jet_matchfrac_sumpt_vseta ->Draw();
  h2_jet_matchfrac_sumpt_vseta ->Draw("same");
  h3_jet_matchfrac_sumpt_vseta ->Draw("same");

    TLegend *leg_jet_matchfrac_sumpt_vseta = new TLegend(0.82,0.18,0.94,0.30);
    leg_jet_matchfrac_sumpt_vseta->AddEntry(h1_jet_matchfrac_sumpt_vseta, label1, "p");
    leg_jet_matchfrac_sumpt_vseta->AddEntry(h2_jet_matchfrac_sumpt_vseta, label2, "p");
    leg_jet_matchfrac_sumpt_vseta->AddEntry(h3_jet_matchfrac_sumpt_vseta, label3, "p");
    leg_jet_matchfrac_sumpt_vseta->Draw();
  
  c.SaveAs(DIR+type+overType+"_jet_matchfrac_sumpt_vseta.png");


  // jet_matchfrac_sumpt_vsphi
  h1_jet_matchfrac_sumpt_vsphi ->SetMarkerColor(kBlack);
  h1_jet_matchfrac_sumpt_vsphi ->SetMarkerStyle(kCircle);
  h2_jet_matchfrac_sumpt_vsphi ->SetMarkerColor(kBlue);
  h2_jet_matchfrac_sumpt_vsphi ->SetMarkerStyle(kPlus);
  h3_jet_matchfrac_sumpt_vsphi ->SetMarkerColor(kRed);
  h3_jet_matchfrac_sumpt_vsphi ->SetMarkerStyle(kMultiply);
  
  h1_jet_matchfrac_sumpt_vsphi ->SetTitle(htitle);
  h1_jet_matchfrac_sumpt_vsphi ->Draw();
  h2_jet_matchfrac_sumpt_vsphi ->Draw("same");
  h3_jet_matchfrac_sumpt_vsphi ->Draw("same");

    TLegend *leg_jet_matchfrac_sumpt_vsphi = new TLegend(0.82,0.18,0.94,0.30);
    leg_jet_matchfrac_sumpt_vsphi->AddEntry(h1_jet_matchfrac_sumpt_vsphi, label1, "p");
    leg_jet_matchfrac_sumpt_vsphi->AddEntry(h2_jet_matchfrac_sumpt_vsphi, label2, "p");
    leg_jet_matchfrac_sumpt_vsphi->AddEntry(h3_jet_matchfrac_sumpt_vsphi, label3, "p");
    leg_jet_matchfrac_sumpt_vsphi->Draw();
  
  c.SaveAs(DIR+type+overType+"_jet_matchfrac_sumpt_vsphi.png");
*/


// ---------------------------------------------------------------------------------------------------------
// close overlay()
// ---------------------------------------------------------------------------------------------------------

} 



// ---------------------------------------------------------------------------------------------------------
// Plot Style
// ---------------------------------------------------------------------------------------------------------

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
