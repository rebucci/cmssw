/*
  ROOT macro showing an example of extracted data analysis.

  Use:
    root[1]-> .L AnalyzeDigis.C++
    root[2]->  AnalyzeDigis(filename)
  where filename is the name of the extracted ROOT file

  Created by Seb Viret (viret@in2p3_dot_fr), 16/05/2013.
  Modified by Rachael Bucci (rbucci@nd.edu), June 2018.

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

void AnalyzeDigis(TString type) {

  gROOT->SetBatch();
  gErrorIgnoreLevel = kWarning;


  // ---------------------------------------------------------------------------------------------------------
  // Retrieve Data

  TString srcDIR = "";

  TChain *Pix = new TChain("Pixels");   
  Pix->Add(srcDIR+type+".root");

  if (Pix->GetEntries() == 0) {
    cout << "File doesn't exist or is empty, returning..." << endl;
    return;
  }


  // ---------------------------------------------------------------------------------------------------------
  // Definitions for Pixel Trees
  
  int    m_pix_n;
  vector<float>*  m_pix_x;  
  vector<float>*  m_pix_y;  
  vector<float>*  m_pix_z;  
  vector<float>*  m_pix_e;  
  vector<int>*    m_pix_row;   
  vector<int>*    m_pix_column;
  vector<int>*    m_pix_simhit;
  vector<int>*    m_pix_layer; 
  vector<int>*    m_pix_module;
  vector<int>*    m_pix_ladder;

  m_pix_n       = 0;
  m_pix_x       = 0;
  m_pix_y       = 0;
  m_pix_z       = 0;
  m_pix_e       = 0;
  m_pix_row     = 0;
  m_pix_column  = 0;
  m_pix_simhit  = 0;
  m_pix_layer   = 0;
  m_pix_module  = 0;
  m_pix_ladder  = 0;

  Pix->SetBranchAddress("PIX_n",         &m_pix_n);
  Pix->SetBranchAddress("PIX_x",         &m_pix_x);
  Pix->SetBranchAddress("PIX_y",         &m_pix_y);
  Pix->SetBranchAddress("PIX_z",         &m_pix_z);
  Pix->SetBranchAddress("PIX_charge",    &m_pix_e);
  Pix->SetBranchAddress("PIX_row",       &m_pix_row);
  Pix->SetBranchAddress("PIX_column",    &m_pix_column);
  Pix->SetBranchAddress("PIX_simhit",    &m_pix_simhit);
  Pix->SetBranchAddress("PIX_layer",     &m_pix_layer);
  Pix->SetBranchAddress("PIX_module",    &m_pix_module);
  Pix->SetBranchAddress("PIX_ladder",    &m_pix_ladder);



  // ---------------------------------------------------------------------------------------------------------
  // Histogram Definition

  double PI = 4.*atan(1.);
  TH2F* Barrel_hitmap = new TH2F("BHmap",";Digi z [cm]; Digi #phi [rad]", 60,-120.,120.,20,-PI,PI);


  // ---------------------------------------------------------------------------------------------------------
  // Data Processing
  
  int   n_entries = Pix->GetEntries();
  cout << "Number of Events = " << n_entries << endl;

  // Loop over entries
  for (int j=0;j<n_entries;++j){
    if (j%1000==0){
      cout << j << endl;
    }

    Pix->GetEntry(j);

    // Loop over each layer and fill the barrel hitmap
    //int   count_fill = 0; 
    for (int i=0;i<m_pix_n;++i){
      if (m_pix_layer->at(i) == 5){
        //count_fill++;
        Barrel_hitmap->Fill(m_pix_z->at(i),atan2(m_pix_y->at(i),m_pix_x->at(i)));
      }
    }
    //cout << "For entry j = " << j << endl;
    //cout << count_fill << " pixels were added to the barrel hitmap" << endl;
  }
  

  // ---------------------------------------------------------------------------------------------------------
  // Draw and Save Plots

  TString plotDIR = "ExtrPlots/";
  TFile* fout = new TFile(plotDIR+"plotsDigi_"+type+".root","recreate");


  // Title
  TString htitle;        // This is the title of all histograms
  htitle = type;     // By default, the histogram titles will be the input file name

  // Make Color Table
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptStat(0);

  // Canvas
  TCanvas* c = new TCanvas("c","Tracker barrel layer digi map",200,80,1100,640);
  c->Range(-149,-3.9,149.,3.6);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->Range(-149,-3.9,149.,3.6);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetRightMargin(0.066);
  c->SetTopMargin(0.066);
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);

  // Text Labels
  Barrel_hitmap->SetTitle(htitle);
  Barrel_hitmap->GetXaxis()->SetTitle("Digi z (in cm)");
  Barrel_hitmap->GetXaxis()->SetLabelFont(42);
  Barrel_hitmap->GetXaxis()->SetLabelSize(0.035);
  Barrel_hitmap->GetXaxis()->SetTitleSize(0.035);
  Barrel_hitmap->GetXaxis()->SetTitleFont(42);
  Barrel_hitmap->GetYaxis()->SetTitle("Digi #phi (in rad)");
  Barrel_hitmap->GetYaxis()->SetLabelFont(42);
  Barrel_hitmap->GetYaxis()->SetLabelSize(0.035);
  Barrel_hitmap->GetYaxis()->SetTitleSize(0.035);
  Barrel_hitmap->GetYaxis()->SetTitleFont(42);
  Barrel_hitmap->GetZaxis()->SetLabelFont(42);
  Barrel_hitmap->GetZaxis()->SetLabelSize(0.035);
  Barrel_hitmap->GetZaxis()->SetTitleSize(0.035);
  Barrel_hitmap->GetZaxis()->SetTitleFont(42);
  Barrel_hitmap->Draw("col");
  Barrel_hitmap->Write();
  
  c->Update();

  c->SaveAs(plotDIR+type+"_BHmap.png");

}
