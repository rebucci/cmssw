/*
  Small ROOT macro showing an example of extracted data analysis.

  RECIPE 3: Matching simhits and clusters

  Use:

  root[1]-> .L Cluster2SimHits.C++
  root[2]-> Cluster2SimHits(type,evt) 

  where type is the name of the extracted ROOT file and evt the number of entry to analyze

  Information about this macro may be found in the following page:

  http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto (STEP III)

  Author: viret@in2p3_dot_fr
  Date: 16/05/2013

*/


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


// -----------------------------------------------------------------------------------------------------------
// Main script
// -----------------------------------------------------------------------------------------------------------


void Cluster2SimHits(TString type, int evt) {
  
  gROOT->SetBatch();
  gErrorIgnoreLevel = kWarning;  

  
  // ---------------------------------------------------------------------------------------------------------
  // Retrieve Data

  TString srcDIR = "";

  TChain *MC              = new TChain("MC");  
  TChain *L1TT            = new TChain("L1TrackTrigger");  

  MC  ->Add(srcDIR+type+".root");
  L1TT->Add(srcDIR+type+".root");

  if (MC->GetEntries() == 0) {
    cout << "MC file doesn't exist or is empty, returning..." << endl;
    return;
  }

  if (L1TT->GetEntries() == 0) {
    cout << "L1TT file doesn't exist or is empty, returning..." << endl;
    return;
  }

  // ---------------------------------------------------------------------------------------------------------
  // Definitions for MC Trees

  int    		m_gen_n;    
  int    		m_part_n;    
  int       m_part_nstracks;  
  int    		m_part_nhit;
 
  vector<float>*  m_gen_x;    
  vector<float>*  m_gen_y;    
  vector<float>*  m_gen_z;    
  vector<float>*	m_gen_px;   
  vector<float>*	m_gen_py;   
  vector<float>*  m_gen_pz;   
  vector<int>*    m_gen_pdg;  
  vector<int>*    m_part_pdgId;
  vector<float>* 	m_part_px;   
  vector<float>* 	m_part_py;   
  vector<float>*	m_part_pz;   
  vector<float>*  m_part_eta;  
  vector<float>* 	m_part_phi;  
  vector<float>*  m_part_x;    
  vector<float>*  m_part_y;    
  vector<float>*  m_part_z;    
  vector<int>*    m_st;        
  vector<int>*    m_hits;      
  vector<int>*    m_st_id;  
  vector<float>*  m_hits_x;    
  vector<float>*  m_hits_y;    
  vector<float>*  m_hits_z;    
  vector<float>*  m_hits_e;    
  vector<float>*  m_hits_tof;  
  vector<int>*    m_hits_proc; 
  vector<int>*    m_hits_id;   
  vector<int>*    m_hits_pdgId;
  vector<int>*    m_hits_layer;
  vector<int>*    m_hits_ladder;
  vector<int>*    m_hits_module;

  m_gen_n         =0;    
  m_part_n        =0;    
  m_part_nstracks =0;  
  m_part_nhit     =0;
  m_gen_x    =0;
  m_gen_y    =0;
  m_gen_z    =0;
  m_gen_px   =0;
  m_gen_py   =0;
  m_gen_pz   =0;
  m_gen_pdg     =0;
  m_part_pdgId  =0;
  m_part_px     =0;
  m_part_py     =0;
  m_part_pz     =0;
  m_part_eta    =0;
  m_part_phi    =0;
  m_part_x      =0;
  m_part_y      =0;
  m_part_z      =0;
  m_st          =0;
  m_hits        =0;
  m_st_id       =0;
  m_hits_x      =0;
  m_hits_y      =0;
  m_hits_z      =0;
  m_hits_e      =0;
  m_hits_tof    =0;
  m_hits_proc   =0;
  m_hits_id     =0;
  m_hits_pdgId  =0;
  m_hits_layer  =0;
  m_hits_ladder =0;
  m_hits_module =0;

  MC->SetBranchAddress("gen_n",            &m_gen_n);
  MC->SetBranchAddress("gen_pdg",          &m_gen_pdg);
  MC->SetBranchAddress("gen_px",           &m_gen_px);
  MC->SetBranchAddress("gen_py",           &m_gen_py);
  MC->SetBranchAddress("gen_pz",           &m_gen_pz);
  MC->SetBranchAddress("gen_x",            &m_gen_x);
  MC->SetBranchAddress("gen_y",            &m_gen_y);
  MC->SetBranchAddress("gen_z",            &m_gen_z);
  MC->SetBranchAddress("subpart_n",        &m_part_n);
  MC->SetBranchAddress("subpart_hits",     &m_hits);
  MC->SetBranchAddress("subpart_st",       &m_st);
  MC->SetBranchAddress("subpart_pdgId",    &m_part_pdgId);
  MC->SetBranchAddress("subpart_px",       &m_part_px);
  MC->SetBranchAddress("subpart_py",       &m_part_py);
  MC->SetBranchAddress("subpart_pz",       &m_part_pz);
  MC->SetBranchAddress("subpart_eta",      &m_part_eta);
  MC->SetBranchAddress("subpart_phi",      &m_part_phi);
  MC->SetBranchAddress("subpart_x",        &m_part_x);
  MC->SetBranchAddress("subpart_y",        &m_part_y);
  MC->SetBranchAddress("subpart_z",        &m_part_z);
  MC->SetBranchAddress("subpart_nhit",     &m_part_nhit);
  MC->SetBranchAddress("subpart_hits_x",   &m_hits_x);
  MC->SetBranchAddress("subpart_hits_y",   &m_hits_y);
  MC->SetBranchAddress("subpart_hits_z",   &m_hits_z);
  MC->SetBranchAddress("subpart_hits_e",   &m_hits_e);
  MC->SetBranchAddress("subpart_hits_tof", &m_hits_tof);
  MC->SetBranchAddress("subpart_hits_proc",&m_hits_proc);
  MC->SetBranchAddress("subpart_hits_id",  &m_hits_id);
  MC->SetBranchAddress("subpart_hits_pdgId",  &m_hits_pdgId);
  MC->SetBranchAddress("subpart_hits_layer",  &m_hits_layer);
  MC->SetBranchAddress("subpart_hits_ladder", &m_hits_ladder);
  MC->SetBranchAddress("subpart_hits_module", &m_hits_module);
  MC->SetBranchAddress("subpart_nsimtracks",  &m_part_nstracks);
  MC->SetBranchAddress("subpart_st_ids",      &m_st_id);


  // ---------------------------------------------------------------------------------------------------------
  // Definitions for L1TrackTrigger Trees

  int m_clus_n;

  vector<float>*  m_clus_x;     
  vector<float>*  m_clus_y;     
  vector<float>*  m_clus_z;     
  vector<float>*  m_clus_xmc;   
  vector<float>*  m_clus_ymc;   
  vector<float>*  m_clus_zmc;   
  vector<int>*    m_clus_layer;  
  vector<int>*    m_clus_module; 
  vector<int>*    m_clus_ladder; 
  vector<int>*    m_clus_seg;    
  vector<float>*  m_clus_strip;  
  vector<int>*    m_clus_sat;    
  vector<int>*    m_clus_nstrips;
  vector<int>*    m_clus_matched;
  vector<int>*    m_clus_PS;     
  vector<int>*    m_clus_nrows;  
  vector< vector<int>* >  m_clus_tp;   
  vector< vector<int>* >  m_clus_hits; 
  vector<int>*    m_clus_pid; 

  m_clus_n        =0;
  m_clus_x        =0;     
  m_clus_y        =0;     
  m_clus_z        =0;     
  m_clus_xmc      =0;   
  m_clus_ymc      =0;   
  m_clus_zmc      =0;   
  m_clus_layer    =0;  
  m_clus_module   =0; 
  m_clus_ladder   =0; 
  m_clus_seg      =0;    
  m_clus_strip    =0;  
  m_clus_sat      =0;    
  m_clus_nstrips  =0;
  m_clus_matched  =0;
  m_clus_PS       =0;     
  m_clus_nrows    =0;  
  //m_clus_tp       =0;   
  //m_clus_hits     =0; 
  m_clus_pid      =0; 

  L1TT->SetBranchAddress("CLUS_n",         &m_clus_n);
  L1TT->SetBranchAddress("CLUS_x",         &m_clus_x);
  L1TT->SetBranchAddress("CLUS_y",         &m_clus_y);
  L1TT->SetBranchAddress("CLUS_z",         &m_clus_z);
  L1TT->SetBranchAddress("CLUS_xmc",       &m_clus_xmc);
  L1TT->SetBranchAddress("CLUS_ymc",       &m_clus_ymc);
  L1TT->SetBranchAddress("CLUS_zmc",       &m_clus_zmc);
  L1TT->SetBranchAddress("CLUS_layer",     &m_clus_layer);
  L1TT->SetBranchAddress("CLUS_module",    &m_clus_module);
  L1TT->SetBranchAddress("CLUS_ladder",    &m_clus_ladder);
  L1TT->SetBranchAddress("CLUS_seg",       &m_clus_seg);
  L1TT->SetBranchAddress("CLUS_strip",     &m_clus_strip);
  L1TT->SetBranchAddress("CLUS_nstrip",    &m_clus_nstrips);
  L1TT->SetBranchAddress("CLUS_nsat",      &m_clus_sat);
  L1TT->SetBranchAddress("CLUS_match",     &m_clus_matched);
  L1TT->SetBranchAddress("CLUS_PS",        &m_clus_PS);
  L1TT->SetBranchAddress("CLUS_nrows",     &m_clus_nrows);
  L1TT->SetBranchAddress("CLUS_tp",        &m_clus_tp);
  L1TT->SetBranchAddress("CLUS_hits",      &m_clus_hits);
  L1TT->SetBranchAddress("CLUS_process",   &m_clus_pid);


  // ---------------------------------------------------------------------------------------------------------
  // Data Processing

  int n_entries = MC->GetEntries();
  cout << "Number of Events = " << n_entries << endl;

  if (evt>=n_entries) evt = n_entries-1;
  if (evt<0) evt = 0;

  MC  ->GetEntry(evt);
  L1TT->GetEntry(evt);

  cout << endl;
  cout << "Analyzing event " << evt << endl;
  cout << "where "  << m_clus_n << " cluster(s) where produced" << endl;

  /////////////////////
  // Loop over cluster
  for (int k=0;k<m_clus_n;++k) {
    cout << "-------------------------------------------------"  << endl;
    cout << "Cluster " << k << " properties:"  
         << endl;
    cout << " X/Y/Z (in cm) : " << m_clus_x->at(k) 
         << "/" << m_clus_y->at(k) 
	       << "/" << m_clus_z->at(k) 
         << endl;
    cout << " Width (in strips) : " << m_clus_nstrips->at(k) 
         << endl;

    // If the cluster is matched
    if (m_clus_matched->at(k)!=0) {
      cout << " X/Y/Z (MC) (in cm) : " << m_clus_xmc->at(k) 
	         << "/" << m_clus_ymc->at(k) 
	         << "/" << m_clus_zmc->at(k) 
           << endl;
      cout << " Matched with " << m_clus_matched->at(k) 
            << " simhit(s) (the same simhit could appear twice or more):" 
            << endl;

      for (int l=0;l<m_clus_matched->at(k);++l) {
        cout << "    SimHit " << m_clus_hits.at(k)->at(l) 
             << " from TP "   << m_clus_tp.at(k)->at(l) 
             << " : "         << m_hits_x->at(m_clus_hits.at(k)->at(l)) 
	           << "/"           << m_hits_y->at(m_clus_hits.at(k)->at(l)) 
	           << "/"           << m_hits_z->at(m_clus_hits.at(k)->at(l)) 
             << endl; 
      }
    }
    else {
      cout << "Unmatched" << endl;
    }
  }
}
