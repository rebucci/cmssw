#ifndef EFFICIENCIES_H
#define EFFICIENCIES_H


#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>



#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include <fstream>
#include <string>

using namespace std;

///////////////////////////////////
//
//
// Base class for stub efficiencies calculation (in stubs/module/BX)
//
// This code computes the efficiencies per module or the efficiencies per layer
// both for official and private stub producer
//
// For infos about the efficiency definition, have a look at the following presentation:
//
// https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=263068
//
// Input infos are :
//
// filename : the name and directory of the input ROOT file containing the STUB information
// outfile  : the name of the output ROOT file containing the rates 
// ptype    : the pdg ID of the particle type you want to test 
//
// Info about the code:
//
//  http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto (STEP III)
//
//  Author: viret@in2p3_dot_fr
//  Date: 23/05/2013
//
///////////////////////////////////



class efficiencies
{
 public:

  efficiencies(std::string filename, std::string outfile, int ptype);

  void  get_efficiencies();  // The main method  
  void  initVars();
  void  reset();
  void  initTuple(std::string in,std::string out);

 private:

  TChain *L1TT;      // The trees
  TChain *Pix;
  TChain *MC;      

  TFile *m_outfile;  // The output file
  TTree *m_tree;     // The tree containing the efficiency information

  float pt_val[100];
  float eta_val[50];
  float entries_pt[30][100];
  float entries_eta[30][50];
  float entries_pt_lay[30][100];
  float entries_eta_lay[30][50];

  float digi_pt[30][100]; // Digi efficiencies as function of pT (between 0 and 20 GeV/c)
  float digi_eta[30][50]; // Digi efficiencies as function of eta (between -2.5 and 2.5)

  // Module efficiencies
  float clus_pt[30][100]; // Official cluster efficiencies as function of pT
  float stub_pt[30][100]; // Official stub efficiencies as function of pT
  float clus_eta[30][50]; // Official cluster efficiencies as function of eta
  float stub_eta[30][50]; // Official stub efficiencies as function of eta


  // Layer efficiencies 
  float clus_pt_lay[30][100]; // Official clus efficiencies as function of pT
  float clus_eta_lay[30][50]; // Official clus efficiencies as function of eta
  float stub_pt_lay[30][100]; // Official stub efficiencies as function of pT
  float stub_eta_lay[30][50]; // Official stub efficiencies as function of eta

  // Here are the parameters needed from the data
  // Details on these might be found on
  //
  // https://github.com/sviret/HL_LHC/blob/master/Extractors/RecoExtractor/interface/L1TrackTrigger_analysis.h
  //

  int                 m_type;
  int    	          m_part_n;
  int                 m_pclus;
  int                 m_tkclus;            
  int                 m_tkstub; 

  std::vector<int>    *m_pixclus_row;
  std::vector<int>    *m_pixclus_column;
  std::vector<int>    *m_pixclus_layer;
  std::vector<int>    *m_pixclus_module;
  std::vector<int>    *m_pixclus_bottom;
  std::vector<int>    *m_pixclus_ladder;
  std::vector<int>    *m_pixclus_type;
  std::vector<float>  *m_pixclus_x;
  std::vector<float>  *m_pixclus_y;
  std::vector< std::vector<int> > *m_pixclus_simhitID; 
  std::vector< std::vector<int> > *m_pixclus_evtID; 

  std::vector<int>    *m_part_hits;
  std::vector<int>    *m_part_pdg;
  std::vector<float>  *m_part_px;
  std::vector<float>  *m_part_py;
  std::vector<float>  *m_part_eta;
  std::vector<float>  *m_part_x;
  std::vector<float>  *m_part_y;
  std::vector<int>    *m_part_evtId;//
  std::vector< std::vector<int> > *m_part_stId;// 

  std::vector<int>    *m_tkclus_nstrips;
  std::vector<int>    *m_tkclus_layer; 
  std::vector<int>    *m_tkclus_module;	
  std::vector<int>    *m_tkclus_bottom;								       
  std::vector<int>    *m_tkclus_ladder;				      
  std::vector<int>    *m_tkclus_seg;   
  std::vector<float>  *m_tkclus_strip;
  std::vector<int>    *m_tkstub_layer;
  std::vector<int>    *m_tkstub_type;
  std::vector<int>    *m_tkstub_ladder;
  std::vector<int>    *m_tkstub_clust1;
  std::vector<int>    *m_tkstub_clust2;
  std::vector<int>    *m_tkstub_tp;
  std::vector<int>    *m_tkclus_tp;
  std::vector<int>    *m_tkclus_pdg;
};

#endif

