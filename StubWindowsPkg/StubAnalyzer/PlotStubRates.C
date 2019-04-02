// -----------------------------------------------------------------------------------------------------------

/*

ROOT macro for stub/cluster rates visualization, to be used with a file created with ./AM_ana -c rates 

Use:
    root[0]-> .L PlotStubRates.C
    root[1]->PlotStubRates("sourcefile",outputname",pu)
where 
    - sourcefile is the name of the root output file produced by AM_ana. It does not include the ".root" suffix.
    - outputname is the name (no type suffix like ".root") of the root file this macro produces. it is also the prefix of all the plot names.
    - both input and output files have customizable preceding paths and prefixes. These options can be changed in the Preliminaries section of this macro.
    - pu is the pileup value of the source file. It is used in the title only.


Options
    - geometry: tilted (default) or flat
    - rate units: in MHz/cm^2 or in stubs/module/BX
    - print each layer: stub, cluster, and ratio plots produced for each of the six barrel layers


Plots: Individual Barrel Layers
    - If PrintEachLayer=true, four canvases (and corresponding .pngs) are printed. 
    - The first canvas contains four histograms: primary, secondary, fake, and total stubs. 
    - The other three canvases have only one histogram each: one for stubs, one for clusters, and one for the ratio of clusters/stubs. 
    - The single histogram stub canvas by default displays primary, secondary, and fake stubs. You can change this by turning on (1) or off (0) the floats pri, sec, and fak in the Preliminaries section of this macro. If you leave all values set to 1, this plot will be the same as the total stub plot in the four-histogram canvas.


*/

// -----------------------------------------------------------------------------------------------------------

#include "TROOT.h"
#include "TStyle.h"
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
#include "TMath.h"
#include "TEfficiency.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TError.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <iomanip>
#include <cmath>
#include <fstream>

void SetPlotStyle();
void SetColorTable();
void fill_histo(TH2F *hist,float val[], double nlad, double dec, int disk, int ring,float max);

void PlotStubRates(TString sourcefile, TString outputname, int pu)
{
  // -----------------------------------------------------------------------------------------------------------
  // Preliminaries
  // -----------------------------------------------------------------------------------------------------------

  gROOT->SetBatch();
  gErrorIgnoreLevel = kError;
  SetPlotStyle();
  SetColorTable();
 
  TString srcDIR  = ""; 
  TString subDIR  = "plots/rates/";
  TString plotDIR = "plots/rates/"; 
  TString plotPRE = "stubRates_";

  // Options
  bool PrintEachLayer = false; // saves each layer plot
  bool PrintEachDisk  = false; // saves each disk plot
  bool PrintTIBRate   = true;
  bool PrintTOBRate   = true;
  bool PrintEndcapRate= true;
  bool RateInHz       = false; // rate is given in MHz/cm^2 (true) or in stubs/module/BX (false)
  bool SavePlotSource = false; // save the macro.C files for the 1D maps

  int  tilt = 1; // Tilted Geometry = 1. Flat geometry = 0.

  // Stub Options
  float pri = 1;    // use primaries = 1
  float sec = 1;    // use secondaries = 1
  float fak = 1;    // use fakes = 1

  // Read in File
  TChain *newtree = new TChain("L1Rates");
  newtree->Add(srcDIR+sourcefile+".root"); 

  // Output file
  TFile* fout = new TFile(subDIR+"plots_"+plotPRE+outputname+".root","recreate");

  
  // -----------------------------------------------------------------------------------------------------------
  // Definitions
  // -----------------------------------------------------------------------------------------------------------


  // -----------------------------------------------------------------------------------------------------------
  // Variables

  // layer map(s)
  float count_sbrp[58000];
  float count_sbrp2[58000]; // needed for TIB and TOB rate
  float count_sbrs[58000];
  float count_sbrf[58000];
  float count_cbr[58000];

  // layer 1D map
  float count_sblr[600];
  float count_cblr[600];

  // disk 1D map
  float count_selr[1500];
  float count_celr[1500];

  // disk map
  float count_serp[142000];
  float count_serp2[142000]; // needed for endcap rate
  float count_sers[142000];
  float count_serf[142000];
  float count_cer[142000];

  float count_e_rate[142000];
  float count_e_clus[142000]; 
  float count_e_ratio[142000];

  // 1D plot colors
  int color[6]  = {24,20,27,34,28,25};


  // -----------------------------------------------------------------------------------------------------------
  // Branch Addresses 
  
  // layers
  newtree->SetBranchAddress("STUB_b_rates_prim",  &count_sbrp);
  newtree->SetBranchAddress("STUB_b_rates_prim2", &count_sbrp2);
  newtree->SetBranchAddress("STUB_b_rates_sec",   &count_sbrs);
  newtree->SetBranchAddress("STUB_b_rates_f",     &count_sbrf);
  newtree->SetBranchAddress("CLUS_b_rates",       &count_cbr);

  newtree->SetBranchAddress("STUB_b_l_rates",     &count_sblr);
  newtree->SetBranchAddress("CLUS_b_l_rates",     &count_cblr);

  // disks
  newtree->SetBranchAddress("STUB_e_l_rates",     &count_selr);
  newtree->SetBranchAddress("CLUS_e_l_rates",     &count_celr);

  newtree->SetBranchAddress("STUB_e_rates_prim",  &count_serp);
  newtree->SetBranchAddress("STUB_e_rates_prim2", &count_serp2);
  newtree->SetBranchAddress("STUB_e_rates_sec",   &count_sers);
  newtree->SetBranchAddress("STUB_e_rates_f",     &count_serf);
  newtree->SetBranchAddress("CLUS_e_rates",       &count_cer);

  // initialize 
  newtree->GetEntry(0);


  // -----------------------------------------------------------------------------------------------------------
  // Geometry

  // active surface of modules (in cm^2)
  double PS_surf = 45.07;   
  double SS_surf = 91.90;

  double n_mod_surf_barrel[6];
  for (int i=0;i<3;++i) n_mod_surf_barrel[i] = PS_surf;
  for (int i=3;i<6;++i) n_mod_surf_barrel[i] = SS_surf;

  double n_mod_surf_endcap[15];
  for (int i=0;i<9;++i) n_mod_surf_endcap[i] = PS_surf;
  for (int i=9;i<15;++i) n_mod_surf_endcap[i] = SS_surf;
 
  float bit_per_s = 1;   // Bit number per stub
  
  // barrel layers
  double n_lad_flat[6]  = {16,24,34,48,62,76};
  double n_mod_flat[6]  = {63,55,54,24,24,24};

  double n_lad_tilt[6]  = {18,26,36,48,60,78};
  double n_mod_tilt[6]  = {31,35,39,24,24,24};

  double n_lad_barrel[6];
  double n_mod_barrel[6];

  if (tilt==1) {
    for (int i=0;i<6;++i) n_lad_barrel[i] = n_lad_tilt[i];
    for (int i=0;i<6;++i) n_mod_barrel[i] = n_mod_tilt[i];
  }
  else {
    for (int i=0;i<6;++i) n_lad_barrel[i] = n_lad_flat[i];
    for (int i=0;i<6;++i) n_mod_barrel[i] = n_mod_flat[i];
  }

  // endcap disks
  double n_lad_endcap[15]   = {20,24,24,28,32,32,36,40,40,44,52,60,64,72,76};
  double n_start_endcap[15] = {0,0,0,0,0,0,1,0,0,0,1,0,0,0,1};

  double n_mod_endcap[15]   = {20,24,24,28,32,32,36,40,40,44,52,60,64,72,76};
  float radii[15]           = {25,30,33,38,41,46,48,53,56,63,71,81,88,98,104};


  // -----------------------------------------------------------------------------------------------------------
  // Text in Plots

  char txt_pu[80];
  sprintf (txt_pu, "#sqrt{s}=14TeV, %d PU", pu);


  
  // ---------------------------------------------------------------------------------------------------------
  //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
  // ---------------------------------------------------------------------------------------------------------



  // -----------------------------------------------------------------------------------------------------------
  // Barrel Layer Maps
  // ---------------------------------------------------------------------------------------------------------

  // Loop over Layers
  for (int layer=5;layer<11;++layer) {
    // Values for EACH layer, not all layers, so they should be reset at the top of the loop.
    
    // from do_layer_maps(). resulting stub plots are four hists (stub pri, sec, fak, tot) on one divided canvas.
    int   idx;
    float b_rate_ALL    = 0;
    float b_maxval_sALL = 0;
    int   b_maxmod_sALL = 0;

    // from do_layer_map(). resulting stub, cluster, ratio plots are one hist per canvas. VAR=variable stub settings. can choose any combination of primary, secondary, and fake stubs for these plots. default setting uses all three.
    float b_rate     = 0;
    float b_ratio    = 0; 
    float b_maxval_s = 0;
    float b_maxval_c = 0; 
    float b_maxval_r = 0;
    int   b_maxmod_s = 0;

    for (int j=(layer-5)*10000;j<(layer-5)*10000+8000;++j) {
      b_rate_ALL = count_sbrp[j]+count_sbrf[j]+count_sbrs[j];
      b_rate     = pri*count_sbrp[j]+fak*count_sbrf[j]+sec*count_sbrs[j];
        
      // ratio, variable stubs (single canvases)
      if (b_rate!=0.) {
        b_ratio = count_cbr[j]/b_rate;
        if (b_ratio>b_maxval_r) b_maxval_r=b_ratio;
      }

      // stub rate for all stubs (divided canvas)
      if (b_rate_ALL>b_maxval_sALL) {
        b_maxval_sALL = b_rate_ALL;
        b_maxmod_sALL = j;
      }

      // stub rate for variable stubs (single canvases)
      if (b_rate>b_maxval_s) {
        b_maxval_s = b_rate;
        b_maxmod_s = j;
      }

      // cluster, variable stubs (single canvases)
      if (count_cbr[j]>b_maxval_c) {
        b_maxval_c = count_cbr[j];
      }
    }
  
    cout << "max module for layer " << layer << " is " << b_maxmod_s << " with " << b_maxval_s << " stubs/mod/BX" << endl;  

    // Histograms
    TH2F *h_layer                = new TH2F("layer",               ";Module Z index; Module #phi index",n_mod_barrel[layer-5],0.,  n_mod_barrel[layer-5],  n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);

    TH2F *h_layer_stub_primary   = new TH2F("layer_stub_primary",  ";Module Z index; Module #phi index",n_mod_barrel[layer-5]+1,0.,n_mod_barrel[layer-5]+1,n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);
    TH2F *h_layer_stub_secondary = new TH2F("layer_stub_secondary",";Module Z index; Module #phi index",n_mod_barrel[layer-5]+1,0.,n_mod_barrel[layer-5]+1,n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);
    TH2F *h_layer_stub_fake      = new TH2F("layer_stub_fake",     ";Module Z index; Module #phi index",n_mod_barrel[layer-5]+1,0.,n_mod_barrel[layer-5]+1,n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);
    TH2F *h_layer_stub_total     = new TH2F("layer_stub_total",    ";Module Z index; Module #phi index",n_mod_barrel[layer-5]+1,0.,n_mod_barrel[layer-5]+1,n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);

    TH2F *h_layer_stub           = new TH2F("layer_stub",      ";Module Z index; Module #phi index",n_mod_barrel[layer-5]+1,0.,n_mod_barrel[layer-5]+1,n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);
    TH2F *h_layer_cluster        = new TH2F("layer_cluster",       ";Module Z index; Module #phi index",n_mod_barrel[layer-5]+1,0.,n_mod_barrel[layer-5]+1,n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);
    TH2F *h_layer_ratio          = new TH2F("layer_ratio",         ";Module Z index; Module #phi index",n_mod_barrel[layer-5]+1,0.,n_mod_barrel[layer-5]+1,n_lad_barrel[layer-5],0.,n_lad_barrel[layer-5]);

    h_layer_stub_primary  ->Fill(n_mod_barrel[layer-5]+0.5,0.5,0.);
    h_layer_stub_primary  ->Fill(n_mod_barrel[layer-5]+0.5,0.5,1.5);
    h_layer_stub_secondary->Fill(n_mod_barrel[layer-5]+0.5,0.5,0.);
    h_layer_stub_secondary->Fill(n_mod_barrel[layer-5]+0.5,0.5,2.5);
    h_layer_stub_fake     ->Fill(n_mod_barrel[layer-5]+0.5,0.5,0.);
    h_layer_stub_fake     ->Fill(n_mod_barrel[layer-5]+0.5,0.5,1.);
    h_layer_stub_total    ->Fill(n_mod_barrel[layer-5]+0.5,0.5,0.);
    h_layer_stub_total    ->Fill(n_mod_barrel[layer-5]+0.5,0.5,4);

    h_layer_stub          ->Fill(n_mod_barrel[layer-5]+0.5,0.5,0.);
    h_layer_stub          ->Fill(n_mod_barrel[layer-5]+0.5,0.5,b_maxval_s);
    h_layer_cluster       ->Fill(n_mod_barrel[layer-5]+0.5,0.5,0.);
    h_layer_cluster       ->Fill(n_mod_barrel[layer-5]+0.5,0.5,b_maxval_c);
    h_layer_ratio         ->Fill(n_mod_barrel[layer-5]+0.5,0.5,0.);
    h_layer_ratio         ->Fill(n_mod_barrel[layer-5]+0.5,0.5,b_maxval_r);

    // Populate Histograms
    for (int j=0;j<n_mod_barrel[layer-5];++j) {
      for (int i=0;i<n_lad_barrel[layer-5];++i)
      {
        idx = 10000*(layer-5) + 100*i + j;
        b_rate_ALL = count_sbrp[idx]+count_sbrf[idx]+count_sbrs[idx];
        b_rate     = pri*count_sbrp[idx]+fak*count_sbrf[idx]+sec*count_sbrs[idx];
  
        h_layer_stub_primary  ->Fill(j+0.5,i+0.5,count_sbrp[idx]);
        h_layer_stub_secondary->Fill(j+0.5,i+0.5,count_sbrs[idx]);
        h_layer_stub_fake     ->Fill(j+0.5,i+0.5,count_sbrf[idx]);
        h_layer_stub_total    ->Fill(j+0.5,i+0.5,b_rate_ALL);

        h_layer_stub          ->Fill(j+0.5,i+0.5,b_rate);
        h_layer_cluster       ->Fill(j+0.5,i+0.5,count_cbr[idx]);
        if (b_rate!=0) h_layer_ratio ->Fill(j+0.5,i+0.5,count_cbr[idx]/b_rate);
      }
    }

    // ----- Layer Map - ALL Stubs ----- //
    TCanvas *c_layer_stubs_all = new TCanvas("c_layer_stubs_all","Layer map - All Stubs",201,77,1470,858);
      c_layer_stubs_all->Range(-5.887851,-1.930603,70.65421,17.37543);
      c_layer_stubs_all->Divide(2,2);
      c_layer_stubs_all->SetGridx();
      c_layer_stubs_all->SetGridy();
        
      // Primaries
      c_layer_stubs_all->cd(1);
        h_layer_stub_primary->GetXaxis()->SetLabelSize(0.02);
        h_layer_stub_primary->GetYaxis()->SetLabelSize(0.02);
        h_layer_stub_primary->GetYaxis()->SetTitleSize(0.03);
        h_layer_stub_primary->GetXaxis()->SetTickLength(1);
        h_layer_stub_primary->GetYaxis()->SetTickLength(0.99);
        h_layer_stub_primary->GetXaxis()->SetNdivisions(n_mod_barrel[layer-5]);
        h_layer_stub_primary->GetYaxis()->SetNdivisions(n_lad_barrel[layer-5]);
        h_layer->Draw("");
        h_layer_stub_primary->Draw("colzsame");
          
        TPaveText *PT_lay_sp = new TPaveText(0.,n_lad_barrel[layer-5],28./55.*n_mod_barrel[layer-5],17./16.*n_lad_barrel[layer-5],"br");
          
        char txt_lay_sp[50];
        sprintf (txt_lay_sp, "Barrel Layer %d primary stub rate (in Stub/Module/BX)",layer);
            
        TLatex TL_lay_sp;
        TL_lay_sp.SetTextSize(0.03);
        TL_lay_sp.DrawLatex(40./55.*n_mod_barrel[layer-5], 16.5/16.*n_lad_barrel[layer-5], "CMS Preliminary Simulation");

        PT_lay_sp->SetFillColor(0);
        PT_lay_sp->SetTextSize(0.03);
        PT_lay_sp->AddText(txt_lay_sp);
        PT_lay_sp->Draw();
        

      // Secondaries
      c_layer_stubs_all->cd(2);
        h_layer_stub_secondary->GetXaxis()->SetLabelSize(0.02);
        h_layer_stub_secondary->GetYaxis()->SetLabelSize(0.02);
        h_layer_stub_secondary->GetYaxis()->SetTitleSize(0.03);
        h_layer_stub_secondary->GetXaxis()->SetTickLength(1);
        h_layer_stub_secondary->GetYaxis()->SetTickLength(0.99);
        h_layer_stub_secondary->GetXaxis()->SetNdivisions(n_mod_barrel[layer-5]);
        h_layer_stub_secondary->GetYaxis()->SetNdivisions(n_lad_barrel[layer-5]);
        h_layer->Draw("");
        h_layer_stub_secondary->Draw("colzsame");
          
        TPaveText *PT_lay_ss = new TPaveText(0.,n_lad_barrel[layer-5],28./55.*n_mod_barrel[layer-5],17./16.*n_lad_barrel[layer-5],"br");

        char txt_lay_ss[50];
        sprintf (txt_lay_ss, "Barrel Layer %d secondary stub rate (in Stub/Module/BX)",layer);
            
        TLatex TL_lay_ss;
        TL_lay_ss.SetTextSize(0.03);
        TL_lay_ss.DrawLatex(40./55.*n_mod_barrel[layer-5], 16.5/16.*n_lad_barrel[layer-5], "CMS Preliminary Simulation");

        PT_lay_ss->SetFillColor(0);
        PT_lay_ss->SetTextSize(0.03);
        PT_lay_ss->AddText(txt_lay_ss);
        PT_lay_ss->Draw();
         

      // Fakes     
      c_layer_stubs_all->cd(3);
        h_layer_stub_fake->GetXaxis()->SetLabelSize(0.02);
        h_layer_stub_fake->GetYaxis()->SetLabelSize(0.02);
        h_layer_stub_fake->GetYaxis()->SetTitleSize(0.03);
        h_layer_stub_fake->GetXaxis()->SetTickLength(1);
        h_layer_stub_fake->GetYaxis()->SetTickLength(0.99);
        h_layer_stub_fake->GetXaxis()->SetNdivisions(n_mod_barrel[layer-5]);
        h_layer_stub_fake->GetYaxis()->SetNdivisions(n_lad_barrel[layer-5]);
        h_layer->Draw("");
        h_layer_stub_fake->Draw("colzsame");
          
        TPaveText *PT_lay_sf = new TPaveText(0.,n_lad_barrel[layer-5],28./55.*n_mod_barrel[layer-5],17./16.*n_lad_barrel[layer-5],"br");

        char txt_lay_sf[50];        
        sprintf (txt_lay_sf, "Barrel Layer %d fake stub rate (in Stub/Module/BX)",layer);
            
        TLatex TL_lay_sf;
        TL_lay_sf.SetTextSize(0.03);
        TL_lay_sf.DrawLatex(40./55.*n_mod_barrel[layer-5], 16.5/16.*n_lad_barrel[layer-5], "CMS Preliminary Simulation");

        PT_lay_sf->SetFillColor(0);
        PT_lay_sf->SetTextSize(0.03);
        PT_lay_sf->AddText(txt_lay_sf);
        PT_lay_sf->Draw();
            

      // Total
      c_layer_stubs_all->cd(4);
        h_layer_stub_total->GetXaxis()->SetLabelSize(0.02);
        h_layer_stub_total->GetYaxis()->SetLabelSize(0.02);
        h_layer_stub_total->GetYaxis()->SetTitleSize(0.03);
        h_layer_stub_total->GetXaxis()->SetTickLength(1);
        h_layer_stub_total->GetYaxis()->SetTickLength(0.99);
        h_layer_stub_total->GetXaxis()->SetNdivisions(n_mod_barrel[layer-5]);
        h_layer_stub_total->GetYaxis()->SetNdivisions(n_lad_barrel[layer-5]);
        h_layer->Draw("");
        h_layer_stub_total->Draw("colzsame");
          
        TPaveText *PT_lay_st = new TPaveText(0.,n_lad_barrel[layer-5],28./55.*n_mod_barrel[layer-5],17./16.*n_lad_barrel[layer-5],"br");

        char txt_lay_st[50];      
        sprintf (txt_lay_st, "Barrel Layer %d total stub rate (in Stub/Module/BX)",layer);
          
        TLatex TL_lay_st;
        TL_lay_st.SetTextSize(0.03);
        TL_lay_st.DrawLatex(40./55.*n_mod_barrel[layer-5], 16.5/16.*n_lad_barrel[layer-5], "CMS Preliminary Simulation");
          
        PT_lay_st->SetFillColor(0);
        PT_lay_st->SetTextSize(0.03);
        PT_lay_st->AddText(txt_lay_st);
        PT_lay_st->Draw();
        

      // Save all four stub plots  
      c_layer_stubs_all->Modified();
      c_layer_stubs_all->Update();
      c_layer_stubs_all->Write();

      if (PrintEachLayer){
        char name_lay_sALL[50];
        sprintf (name_lay_sALL, "_layer_%d_layermap_stubsALL.png", layer);
        c_layer_stubs_all->SaveAs(plotDIR+plotPRE+outputname+name_lay_sALL);
      }

      c_layer_stubs_all->Update();


    // ----- Layer Map Stubs ----- //
    TCanvas *c_layer_stubs = new TCanvas("c_layer_stubs","Layer map - Stubs",201,77,1470,858);
      c_layer_stubs->Range(-5.887851,-1.930603,70.65421,17.37543);
      c_layer_stubs->SetGridx();
      c_layer_stubs->SetGridy();

      h_layer->GetXaxis()->SetLabelSize(0.02);
      h_layer->GetYaxis()->SetLabelSize(0.02);
      h_layer->GetYaxis()->SetTitleSize(0.03);
      h_layer->GetXaxis()->SetTickLength(1);
      h_layer->GetYaxis()->SetTickLength(0.99);
      h_layer->GetXaxis()->SetNdivisions(n_mod_barrel[layer-5]);
      h_layer->GetYaxis()->SetNdivisions(n_lad_barrel[layer-5]);
      h_layer->Draw();
      h_layer_stub->Draw("colzsame");
         
      TPaveText *PT_lay_s = new TPaveText(0.,n_lad_barrel[layer-5],18./55.*n_mod_barrel[layer-5],17./16.*n_lad_barrel[layer-5],"br");

      char txt_lay_s[50];
      sprintf (txt_lay_s, "Barrel Layer %d stub rate (in Stub/Module/BX)",layer);

      TLatex TL_lay;
      TL_lay.SetTextSize(0.02);
      TL_lay.DrawLatex(40./55.*n_mod_barrel[layer-5], 16.5/16.*n_lad_barrel[layer-5], "CMS Preliminary Simulation");  

      PT_lay_s->SetFillColor(0);
      PT_lay_s->SetTextSize(0.02);
      PT_lay_s->AddText(txt_lay_s);
      PT_lay_s->Draw();
        
      c_layer_stubs->Modified();
      c_layer_stubs->Update();
      c_layer_stubs->Write();

      if (PrintEachLayer){
        char name_lay_s[50];
        sprintf (name_lay_s, "_layer_%d_layermap_stubs.png", layer);
        c_layer_stubs->SaveAs(plotDIR+plotPRE+outputname+name_lay_s);
      }

      c_layer_stubs->Update();


    // ----- Layer Map Clusters ----- //
    TCanvas *c_layer_clusters = new TCanvas("c_layer_clusters","Cluster map",201,77,1470,858);
      c_layer_clusters->Range(-5.887851,-1.930603,70.65421,17.37543);

      c_layer_clusters->SetGridx();
      c_layer_clusters->SetGridy();

      h_layer->GetXaxis()->SetLabelSize(0.02);
      h_layer->GetYaxis()->SetLabelSize(0.02);
      h_layer->GetYaxis()->SetTitleSize(0.03);
      h_layer->GetXaxis()->SetTickLength(1);
      h_layer->GetYaxis()->SetTickLength(0.99);
      h_layer->GetXaxis()->SetNdivisions(n_mod_barrel[layer-5]);
      h_layer->GetYaxis()->SetNdivisions(n_lad_barrel[layer-5]);
      h_layer->Draw();
      h_layer_cluster->Draw("colzsame");
         
      TPaveText *PT_lay_c = new TPaveText(0.,n_lad_barrel[layer-5],18./55.*n_mod_barrel[layer-5],17./16.*n_lad_barrel[layer-5],"br");
      
      char txt_lay_c[50];
      sprintf (txt_lay_c, "Barrel Layer %d cluster rate (in Cluster/Module/BX)",layer);

      TL_lay.DrawLatex(40./55.*n_mod_barrel[layer-5], 16.5/16.*n_lad_barrel[layer-5], "CMS Preliminary Simulation");  
      
      PT_lay_c->SetFillColor(0);
      PT_lay_c->SetTextSize(0.02);
      PT_lay_c->AddText(txt_lay_c);
      PT_lay_c->Draw();

      c_layer_clusters->Modified();
      c_layer_clusters->Update();
      c_layer_clusters->Write();

      if (PrintEachLayer) {
        char name_lay_c[50];
        sprintf (name_lay_c, "_layer_%d_layermap_clusters.png", layer);
        c_layer_clusters->SaveAs(plotDIR+plotPRE+outputname+name_lay_c);
      }

      c_layer_clusters->Update();


    // ----- Layer Map Ratio ----- //
    TCanvas *c_layer_ratios = new TCanvas("c_layer_ratios","Cluster over Stub map",201,77,1470,858);
      c_layer_ratios->Range(-5.887851,-1.930603,70.65421,17.37543);

      c_layer_ratios->SetGridx();
      c_layer_ratios->SetGridy();
      
      h_layer->GetXaxis()->SetLabelSize(0.02);
      h_layer->GetYaxis()->SetLabelSize(0.02);
      h_layer->GetYaxis()->SetTitleSize(0.03);
      h_layer->GetXaxis()->SetTickLength(1);
      h_layer->GetYaxis()->SetTickLength(0.99);
      h_layer->GetXaxis()->SetNdivisions(n_mod_barrel[layer-5]);
      h_layer->GetYaxis()->SetNdivisions(n_lad_barrel[layer-5]);
      h_layer->Draw();
      h_layer_ratio->Draw("colzsame");

      TPaveText *PT_lay_r = new TPaveText(0.,n_lad_barrel[layer-5],18./55.*n_mod_barrel[layer-5],17./16.*n_lad_barrel[layer-5],"br");

      char txt_lay_r[50];
      sprintf (txt_lay_r, "Barrel Layer %d cluster/stub ratio",layer);

      TL_lay.DrawLatex(40./55.*n_mod_barrel[layer-5], 16.5/16.*n_lad_barrel[layer-5], "CMS Preliminary Simulation");  

      PT_lay_r->SetFillColor(0);
      PT_lay_r->SetTextSize(0.02);
      PT_lay_r->AddText(txt_lay_r);
      PT_lay_r->Draw();

      c_layer_ratios->Modified();
      c_layer_ratios->Update();
      c_layer_ratios->Write();

      if (PrintEachLayer) {
        char name_lay_r[50];
        sprintf (name_lay_r, "_layer_%d_layermap_ratio.png", layer);
        c_layer_ratios->SaveAs(plotDIR+plotPRE+outputname+name_lay_r);
      }

      c_layer_ratios->Update();
    

  } // close layer loop


  // -----------------------------------------------------------------------------------------------------------
  // 1D Layer Map
  // ---------------------------------------------------------------------------------------------------------

  // Values are for ALL layers and do not get reset at the top of each layer loop.
  float bl_maxval_s = 0;
  float bl_maxval_c = 0;
  float bl_maxval_r = 0;
  float bl_rate_s;
  float bl_rate_c;
  float bl_rate_r;


  // ----- Loop over Layers ----- //

  // First Loop: 
    // Find max values for stub, cluster, and ratio rates for all 6 layers
    // Need these to define the upper bound on Y
  for (int i=0;i<6;++i) { 
    for (int j=i*100;j<(i+1)*100;++j) { 
      bl_rate_s = count_sblr[j]/n_lad_barrel[i];
      bl_rate_c = count_cblr[j]/n_lad_barrel[i];
      bl_rate_r = 0;
      if (bl_rate_s!=0) bl_rate_r = bl_rate_c/bl_rate_s;

      if (RateInHz==true) { 
        bl_rate_s /= n_mod_surf_barrel[i];
        bl_rate_c /= n_mod_surf_barrel[i];
        bl_rate_s *= bit_per_s*40.; // Give the value in MHz/cm2
        bl_rate_c *= bit_per_s*40.; // Give the value in MHz/cm2
      }

      if (bl_rate_s>1000) continue;       

      if (bl_rate_s>bl_maxval_s) bl_maxval_s = bl_rate_s;
      if (bl_rate_c>bl_maxval_c) bl_maxval_c = bl_rate_c;
      if (bl_rate_r>bl_maxval_r) bl_maxval_r = bl_rate_r;
    }
  }

  // Histograms
  TH2F *h_1D_layer_stub    = new TH2F("1D_layer_stub",   ";Module index along Z; Stub rate (number of stubs/module/BX)",1000,-100.,100.,100,0.,1.1*bl_maxval_s);
  TH2F *h_1D_layer_cluster = new TH2F("1D_layer_cluster",";Module index along Z; Cluster rate (number of clusters/module/BX)",1000,-100.,100.,100,0.,1.1*bl_maxval_c);
  TH2F *h_1D_layer_ratio   = new TH2F("1D_layer_ratio",  ";Module index along Z; Cluster over stub ratio",1000,-100.,100.,100,0.,1.1*bl_maxval_r);

  std::vector<TH2F*> h_lay1D_stub_plots; 
  std::vector<TH2F*> h_lay1D_clus_plots; 
  std::vector<TH2F*> h_lay1D_rat_plots; 

  // Second Loop:
    // Define number of bins in X using n_mod_barrel, which is different for each layer
  for (int i=0;i<6;++i) {
    TH2F *h_lay1D_stub = new TH2F("","",20*n_mod_barrel[i],-100.,100.,100,0.,1.1*bl_maxval_s);
    TH2F *h_lay1D_clus = new TH2F("","",20*n_mod_barrel[i],-100.,100.,100,0.,1.1*bl_maxval_c);
    TH2F *h_lay1D_rat = new TH2F("","",20*n_mod_barrel[i],-100.,100.,100,0.,1.1*bl_maxval_r);

    h_lay1D_stub_plots.push_back(h_lay1D_stub);
    h_lay1D_clus_plots.push_back(h_lay1D_clus);
    h_lay1D_rat_plots.push_back(h_lay1D_rat);
  }

  // Populate Histograms
  float bin_w;
  for (int i=0;i<6;++i) {
    bin_w = 200./(20.*n_mod_barrel[i]);

    for (int j=0;j<n_mod_barrel[i];++j) {
      bl_rate_s = count_sblr[100*i+j]/n_lad_barrel[i];
      bl_rate_c = count_cblr[100*i+j]/n_lad_barrel[i];
      if (bl_rate_s!=0) bl_rate_r = bl_rate_c/bl_rate_s;
      if (bl_rate_s==0) continue;
        
      if (RateInHz) { 
        bl_rate_s /= n_mod_surf_barrel[i];
        bl_rate_s *= bit_per_s*40.; // Give the value in MHz/cm2
  
        bl_rate_c /= n_mod_surf_barrel[i];
        bl_rate_c *= bit_per_s*40.; // Give the value in MHz/cm2
      }

      h_lay1D_stub_plots[i]->Fill(-100+20*bin_w*(j+0.5),bl_rate_s);
      h_lay1D_clus_plots[i]->Fill(-100+20*bin_w*(j+0.5),bl_rate_c);
      h_lay1D_rat_plots[i] ->Fill(-100+20*bin_w*(j+0.5),bl_rate_r);
    }
  }

  // ----- 1D Layer Map - Stub Rate ----- //
  TCanvas *c_1D_layer_stubs = new TCanvas("c_1D_layer_stubs","1D Layer map - Stubs",0,0,640,480);
    c_1D_layer_stubs->Range(-1.,-1.,1.,1.);
    c_1D_layer_stubs->SetGridx();
    c_1D_layer_stubs->SetGridy();

    if (RateInHz)  h_1D_layer_stub->GetYaxis()->SetTitle("Stub rate (in MHz/cm^{2})"); 
    h_1D_layer_stub->GetYaxis()->SetTitleOffset(0.83);
    h_1D_layer_stub->Draw();

    for (int i=0;i<6;++i) {
      h_lay1D_stub_plots[i]->SetMarkerStyle(color[i]);
      h_lay1D_stub_plots[i]->SetMarkerSize(1.4);
      h_lay1D_stub_plots[i]->Draw("same");
    }

    TLegend *leg_lay1D_s = new TLegend(0.71,0.69,0.89,0.92);
    leg_lay1D_s->SetTextSize(0.03);
    leg_lay1D_s->SetFillColor(0);
    leg_lay1D_s->AddEntry(h_lay1D_stub_plots[0],"TBPS layer 1","p");
    leg_lay1D_s->AddEntry(h_lay1D_stub_plots[1],"TBPS layer 2","p");
    leg_lay1D_s->AddEntry(h_lay1D_stub_plots[2],"TBPS layer 3","p");
    leg_lay1D_s->AddEntry(h_lay1D_stub_plots[3],"TB2S layer 1","p");
    leg_lay1D_s->AddEntry(h_lay1D_stub_plots[4],"TB2S layer 2","p");
    leg_lay1D_s->AddEntry(h_lay1D_stub_plots[5],"TB2S layer 3","p");
    leg_lay1D_s->Draw();

    TLatex CMS_lay1D_s;
    CMS_lay1D_s.SetTextSize(0.04);
    CMS_lay1D_s.DrawLatex(-80, 1.11*bl_maxval_s, "CMS Phase-2 Simulation");

    TLatex TL_lay1D_s;
    TL_lay1D_s.SetTextSize(0.03);
    TL_lay1D_s.SetTextFont(52);
    TL_lay1D_s.DrawLatex(60, 1.11*bl_maxval_s, txt_pu);

    c_1D_layer_stubs->Modified();
    c_1D_layer_stubs->Update();
    c_1D_layer_stubs->Write();

    char name_lay1D_s[100];
    sprintf (name_lay1D_s, "_1D_layermap_stubs.png");
    if (RateInHz)  sprintf (name_lay1D_s, "_1D_layermap_stubs_inHz.png");
    c_1D_layer_stubs->SaveAs(plotDIR+plotPRE+outputname+name_lay1D_s);

    if (SavePlotSource) {
      char name2_lay1D_s[100];
      sprintf (name2_lay1D_s, "_1D_layermap_stubs.C");
      if (RateInHz)  sprintf (name2_lay1D_s, "_1D_layermap_stubs_inHz.C");
      c_1D_layer_stubs->SaveSource(plotDIR+plotPRE+outputname+name2_lay1D_s);
    }

  // ----- 1D Layer Map - Cluster Rate ----- //
  TCanvas *c_1D_layer_clusters = new TCanvas("c_1D_layer_clusters","1D Layer map - Clusters",0,0,640,480);
    c_1D_layer_clusters->Range(-1.,-1.,1.,1.);
    c_1D_layer_clusters->SetGridx();
    c_1D_layer_clusters->SetGridy();

    if (RateInHz)  h_1D_layer_cluster->GetYaxis()->SetTitle("Cluster rate (in MHz/cm^{2})"); 
    h_1D_layer_cluster->GetYaxis()->SetTitleOffset(0.83);
    h_1D_layer_cluster->Draw();

    for (int i=0;i<6;++i) {
      h_lay1D_clus_plots[i]->SetMarkerStyle(color[i]);
      h_lay1D_clus_plots[i]->SetMarkerSize(1.4);
      h_lay1D_clus_plots[i]->Draw("same");
    }

    TLegend *leg_lay1D_c = new TLegend(0.71,0.69,0.89,0.92);
    leg_lay1D_c->SetTextSize(0.03);
    leg_lay1D_c->SetFillColor(0);
    leg_lay1D_c->AddEntry(h_lay1D_clus_plots[0],"TBPS layer 1","p");
    leg_lay1D_c->AddEntry(h_lay1D_clus_plots[1],"TBPS layer 2","p");
    leg_lay1D_c->AddEntry(h_lay1D_clus_plots[2],"TBPS layer 3","p");
    leg_lay1D_c->AddEntry(h_lay1D_clus_plots[3],"TB2S layer 1","p");
    leg_lay1D_c->AddEntry(h_lay1D_clus_plots[4],"TB2S layer 2","p");
    leg_lay1D_c->AddEntry(h_lay1D_clus_plots[5],"TB2S layer 3","p");
    leg_lay1D_c->Draw();

    TLatex CMS_lay1D_c;
    CMS_lay1D_c.SetTextSize(0.04);
    CMS_lay1D_c.DrawLatex(-80, 1.11*bl_maxval_c, "CMS Phase-2 Simulation");

    TLatex TL_lay1D_c;
    TL_lay1D_c.SetTextSize(0.03);
    TL_lay1D_c.SetTextFont(52);
    TL_lay1D_c.DrawLatex(60, 1.11*bl_maxval_c, txt_pu);

    c_1D_layer_clusters->Modified();
    c_1D_layer_clusters->Update();
    c_1D_layer_clusters->Write();

    char name_lay1D_c[100];
    sprintf (name_lay1D_c, "_1D_layermap_clusters.png");
    if (RateInHz)  sprintf (name_lay1D_c, "_1D_layermap_clusters_inHz.png");
    c_1D_layer_clusters->SaveAs(plotDIR+plotPRE+outputname+name_lay1D_c);

    if (SavePlotSource) {
      char name2_lay1D_c[100];
      sprintf (name2_lay1D_c, "_1D_layermap_clusters.C");
      if (RateInHz)  sprintf (name2_lay1D_c, "_1D_layermap_clusters_inHz.C");
      c_1D_layer_clusters->SaveSource(plotDIR+plotPRE+outputname+name2_lay1D_c);
    }

  // ----- 1D Layer Map - Cluster over Stub Ratio ----- //
  TCanvas *c_1D_layer_ratios = new TCanvas("c_1D_layer_ratios","1D Layer map - Clusters over Stubs",0,0,640,480);
    c_1D_layer_ratios->Range(-1.,-1.,1.,1.);
    c_1D_layer_ratios->SetGridx();
    c_1D_layer_ratios->SetGridy();

    h_1D_layer_ratio->GetYaxis()->SetTitleOffset(0.83);
    h_1D_layer_ratio->Draw();

    for (int i=0;i<6;++i) {
      h_lay1D_rat_plots[i]->SetMarkerStyle(color[i]);
      h_lay1D_rat_plots[i]->SetMarkerSize(1.4);
      h_lay1D_rat_plots[i]->Draw("same");
    }

    TLegend *leg_lay1D_r = new TLegend(0.71,0.69,0.89,0.92);
    leg_lay1D_r->SetTextSize(0.03);
    leg_lay1D_r->SetFillColor(0);
    leg_lay1D_r->AddEntry(h_lay1D_rat_plots[0],"TBPS layer 1","p");
    leg_lay1D_r->AddEntry(h_lay1D_rat_plots[1],"TBPS layer 2","p");
    leg_lay1D_r->AddEntry(h_lay1D_rat_plots[2],"TBPS layer 3","p");
    leg_lay1D_r->AddEntry(h_lay1D_rat_plots[3],"TB2S layer 1","p");
    leg_lay1D_r->AddEntry(h_lay1D_rat_plots[4],"TB2S layer 2","p");
    leg_lay1D_r->AddEntry(h_lay1D_rat_plots[5],"TB2S layer 3","p");
    leg_lay1D_r->Draw();

    TLatex CMS_lay1D_r;
    CMS_lay1D_r.SetTextSize(0.04);
    CMS_lay1D_r.DrawLatex(-80, 1.11*bl_maxval_r, "CMS Phase-2 Simulation");

    TLatex TL_lay1D_r;
    TL_lay1D_r.SetTextSize(0.03);
    TL_lay1D_r.SetTextFont(52);
    TL_lay1D_r.DrawLatex(60, 1.11*bl_maxval_r, txt_pu);

    c_1D_layer_ratios->Modified();
    c_1D_layer_ratios->Update();
    c_1D_layer_ratios->Write();
        
    char name_lay1D_r[100];
    sprintf (name_lay1D_r, "_1D_layermap_ratio.png");
    c_1D_layer_ratios->SaveAs(plotDIR+plotPRE+outputname+name_lay1D_r);

    if (SavePlotSource) {
      char name2_lay1D_r[100];
      sprintf (name2_lay1D_r, "_1D_layermap_ratio.C");
      c_1D_layer_ratios->SaveSource(plotDIR+plotPRE+outputname+name2_lay1D_r);
    }


  // -----------------------------------------------------------------------------------------------------------
  // Endcap Disk Maps
  // ---------------------------------------------------------------------------------------------------------

  // loop over disks
  for (int disk=0;disk<12;++disk) {
    if (disk>=7) { // -Z case
      for (int i=0;i<15;++i) n_start_endcap[i] = 0; // There is no problem in this case
    }

    // Values for EACH disk, not all disks, so they should be reset at the top of the loop.
    int idx;

    // variable stub values
    float e_rate     = 0;
    float e_maxval_s = 0;
    float e_maxval_c = 0;
    float e_maxval_r = 0;

    for (int i=0;i<142000;++i) {
      count_e_rate[i]=0.; 
      count_e_clus[i]=0.;
      count_e_ratio[i]=0.;
    }

     
    for (int i=disk;i<disk+1;++i) { 
      for (int j=10000*i;j<10000*i+2000;++j) {
        e_rate = pri*count_serp[j]+fak*count_serf[j]+sec*count_sers[j];

        // stubs
        count_e_rate[j]=e_rate;
        if (count_e_rate[j] > e_maxval_s) e_maxval_s = count_e_rate[j]; 

        // clusters
        count_e_clus[j]=count_cer[j];
        if (count_e_clus[j] > e_maxval_c) e_maxval_c = count_e_clus[j];

        // ratio
        if (e_rate!=0.){
          count_e_ratio[j] = count_cer[j]/e_rate;
          if (count_e_ratio[j] > e_maxval_r) e_maxval_r = count_e_ratio[j];
        }

      }
    }

    // ----- Histograms----- //

    // Stubs
    TH2F *h_disk_stubs = new TH2F("disk_stubs","zz",30,-15.,15.,30,-15.,15.);
    TH2F *h_ring1_s    = new TH2F("ring1_s","r1",n_lad_endcap[0],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring2_s    = new TH2F("ring2_s","r2",n_lad_endcap[1],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring3_s    = new TH2F("ring3_s","r3",n_lad_endcap[2],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring4_s    = new TH2F("ring4_s","r4",n_lad_endcap[3],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring5_s    = new TH2F("ring5_s","r5",n_lad_endcap[4],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring6_s    = new TH2F("ring6_s","r6",n_lad_endcap[5],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring7_s    = new TH2F("ring7_s","r7",n_lad_endcap[6],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring8_s    = new TH2F("ring8_s","r8",n_lad_endcap[7],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring9_s    = new TH2F("ring9_s","r9",n_lad_endcap[8],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring10_s   = new TH2F("ring10_s","r10",n_lad_endcap[9],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring11_s   = new TH2F("ring11_s","r11",n_lad_endcap[10],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring12_s   = new TH2F("ring12_s","r12",n_lad_endcap[11],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring13_s   = new TH2F("ring13_s","r13",n_lad_endcap[12],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring14_s   = new TH2F("ring14_s","r14",n_lad_endcap[13],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring15_s   = new TH2F("ring15_s","r15",n_lad_endcap[14],0.,8*atan(1.),26,0.,26.);

    fill_histo(h_ring1_s,count_e_rate,n_lad_endcap[0],n_start_endcap[0],disk,0,e_maxval_s);
    fill_histo(h_ring2_s,count_e_rate,n_lad_endcap[1],n_start_endcap[1],disk,1,e_maxval_s);
    fill_histo(h_ring3_s,count_e_rate,n_lad_endcap[2],n_start_endcap[2],disk,2,e_maxval_s);
    fill_histo(h_ring4_s,count_e_rate,n_lad_endcap[3],n_start_endcap[3],disk,3,e_maxval_s);
    fill_histo(h_ring5_s,count_e_rate,n_lad_endcap[4],n_start_endcap[4],disk,4,e_maxval_s);
    fill_histo(h_ring6_s,count_e_rate,n_lad_endcap[5],n_start_endcap[5],disk,5,e_maxval_s);
    fill_histo(h_ring7_s,count_e_rate,n_lad_endcap[6],n_start_endcap[6],disk,6,e_maxval_s);
    fill_histo(h_ring8_s,count_e_rate,n_lad_endcap[7],n_start_endcap[7],disk,7,e_maxval_s);
    fill_histo(h_ring9_s,count_e_rate,n_lad_endcap[8],n_start_endcap[8],disk,8,e_maxval_s);
    fill_histo(h_ring10_s,count_e_rate,n_lad_endcap[9],n_start_endcap[9],disk,9,e_maxval_s);
    fill_histo(h_ring11_s,count_e_rate,n_lad_endcap[10],n_start_endcap[10],disk,10,e_maxval_s);
    fill_histo(h_ring12_s,count_e_rate,n_lad_endcap[11],n_start_endcap[11],disk,11,e_maxval_s);
    fill_histo(h_ring13_s,count_e_rate,n_lad_endcap[12],n_start_endcap[12],disk,12,e_maxval_s);
    fill_histo(h_ring14_s,count_e_rate,n_lad_endcap[13],n_start_endcap[13],disk,13,e_maxval_s);
    fill_histo(h_ring15_s,count_e_rate,n_lad_endcap[14],n_start_endcap[14],disk,14,e_maxval_s);


    // Clusters
    TH2F *h_disk_clusters = new TH2F("disk_clusters","zz",30,-15.,15.,30,-15.,15.);
    TH2F *h_ring1_c       = new TH2F("ring1_c","r1",n_lad_endcap[0],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring2_c       = new TH2F("ring2_c","r2",n_lad_endcap[1],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring3_c       = new TH2F("ring3_c","r3",n_lad_endcap[2],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring4_c       = new TH2F("ring4_c","r4",n_lad_endcap[3],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring5_c       = new TH2F("ring5_c","r5",n_lad_endcap[4],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring6_c       = new TH2F("ring6_c","r6",n_lad_endcap[5],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring7_c       = new TH2F("ring7_c","r7",n_lad_endcap[6],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring8_c       = new TH2F("ring8_c","r8",n_lad_endcap[7],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring9_c       = new TH2F("ring9_c","r9",n_lad_endcap[8],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring10_c      = new TH2F("ring10_c","r10",n_lad_endcap[9],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring11_c      = new TH2F("ring11_c","r11",n_lad_endcap[10],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring12_c      = new TH2F("ring12_c","r12",n_lad_endcap[11],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring13_c      = new TH2F("ring13_c","r13",n_lad_endcap[12],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring14_c      = new TH2F("ring14_c","r14",n_lad_endcap[13],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring15_c      = new TH2F("ring15_c","r15",n_lad_endcap[14],0.,8*atan(1.),26,0.,26.);

    fill_histo(h_ring1_c,count_e_clus,n_lad_endcap[0],n_start_endcap[0],disk,0,e_maxval_c);
    fill_histo(h_ring2_c,count_e_clus,n_lad_endcap[1],n_start_endcap[1],disk,1,e_maxval_c);
    fill_histo(h_ring3_c,count_e_clus,n_lad_endcap[2],n_start_endcap[2],disk,2,e_maxval_c);
    fill_histo(h_ring4_c,count_e_clus,n_lad_endcap[3],n_start_endcap[3],disk,3,e_maxval_c);
    fill_histo(h_ring5_c,count_e_clus,n_lad_endcap[4],n_start_endcap[4],disk,4,e_maxval_c);
    fill_histo(h_ring6_c,count_e_clus,n_lad_endcap[5],n_start_endcap[5],disk,5,e_maxval_c);
    fill_histo(h_ring7_c,count_e_clus,n_lad_endcap[6],n_start_endcap[6],disk,6,e_maxval_c);
    fill_histo(h_ring8_c,count_e_clus,n_lad_endcap[7],n_start_endcap[7],disk,7,e_maxval_c);
    fill_histo(h_ring9_c,count_e_clus,n_lad_endcap[8],n_start_endcap[8],disk,8,e_maxval_c);
    fill_histo(h_ring10_c,count_e_clus,n_lad_endcap[9],n_start_endcap[9],disk,9,e_maxval_c);
    fill_histo(h_ring11_c,count_e_clus,n_lad_endcap[10],n_start_endcap[10],disk,10,e_maxval_c);
    fill_histo(h_ring12_c,count_e_clus,n_lad_endcap[11],n_start_endcap[11],disk,11,e_maxval_c);
    fill_histo(h_ring13_c,count_e_clus,n_lad_endcap[12],n_start_endcap[12],disk,12,e_maxval_c);
    fill_histo(h_ring14_c,count_e_clus,n_lad_endcap[13],n_start_endcap[13],disk,13,e_maxval_c);
    fill_histo(h_ring15_c,count_e_clus,n_lad_endcap[14],n_start_endcap[14],disk,14,e_maxval_c);


    // Ratio
    TH2F *h_disk_ratios = new TH2F("disk_ratios","zz",30,-15.,15.,30,-15.,15.);
    TH2F *h_ring1_r     = new TH2F("ring1_r","r1",n_lad_endcap[0],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring2_r     = new TH2F("ring2_r","r2",n_lad_endcap[1],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring3_r     = new TH2F("ring3_r","r3",n_lad_endcap[2],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring4_r     = new TH2F("ring4_r","r4",n_lad_endcap[3],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring5_r     = new TH2F("ring5_r","r5",n_lad_endcap[4],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring6_r     = new TH2F("ring6_r","r6",n_lad_endcap[5],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring7_r     = new TH2F("ring7_r","r7",n_lad_endcap[6],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring8_r     = new TH2F("ring8_r","r8",n_lad_endcap[7],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring9_r     = new TH2F("ring9_r","r9",n_lad_endcap[8],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring10_r    = new TH2F("ring10_r","r10",n_lad_endcap[9],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring11_r    = new TH2F("ring11_r","r11",n_lad_endcap[10],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring12_r    = new TH2F("ring12_r","r12",n_lad_endcap[11],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring13_r    = new TH2F("ring13_r","r13",n_lad_endcap[12],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring14_r    = new TH2F("ring14_r","r14",n_lad_endcap[13],0.,8*atan(1.),26,0.,26.);
    TH2F *h_ring15_r    = new TH2F("ring15_r","r15",n_lad_endcap[14],0.,8*atan(1.),26,0.,26.);

    fill_histo(h_ring1_r,count_e_ratio,n_lad_endcap[0],n_start_endcap[0],disk,0,e_maxval_r);
    fill_histo(h_ring2_r,count_e_ratio,n_lad_endcap[1],n_start_endcap[1],disk,1,e_maxval_r);
    fill_histo(h_ring3_r,count_e_ratio,n_lad_endcap[2],n_start_endcap[2],disk,2,e_maxval_r);
    fill_histo(h_ring4_r,count_e_ratio,n_lad_endcap[3],n_start_endcap[3],disk,3,e_maxval_r);
    fill_histo(h_ring5_r,count_e_ratio,n_lad_endcap[4],n_start_endcap[4],disk,4,e_maxval_r);
    fill_histo(h_ring6_r,count_e_ratio,n_lad_endcap[5],n_start_endcap[5],disk,5,e_maxval_r);
    fill_histo(h_ring7_r,count_e_ratio,n_lad_endcap[6],n_start_endcap[6],disk,6,e_maxval_r);
    fill_histo(h_ring8_r,count_e_ratio,n_lad_endcap[7],n_start_endcap[7],disk,7,e_maxval_r);
    fill_histo(h_ring9_r,count_e_ratio,n_lad_endcap[8],n_start_endcap[8],disk,8,e_maxval_r);
    fill_histo(h_ring10_r,count_e_ratio,n_lad_endcap[9],n_start_endcap[9],disk,9,e_maxval_r);
    fill_histo(h_ring11_r,count_e_ratio,n_lad_endcap[10],n_start_endcap[10],disk,10,e_maxval_r);
    fill_histo(h_ring12_r,count_e_ratio,n_lad_endcap[11],n_start_endcap[11],disk,11,e_maxval_r);
    fill_histo(h_ring13_r,count_e_ratio,n_lad_endcap[12],n_start_endcap[12],disk,12,e_maxval_r);
    fill_histo(h_ring14_r,count_e_ratio,n_lad_endcap[13],n_start_endcap[13],disk,13,e_maxval_r);
    fill_histo(h_ring15_r,count_e_ratio,n_lad_endcap[14],n_start_endcap[14],disk,14,e_maxval_r);


    // ----- Disk Map Stubs ----- //
    TCanvas *c_disk_stubs = new TCanvas("c_disk_stubs","Disk map - Stubs",166,77,947,880);
      c_disk_stubs->Range(-15.23757,-14.89125,18.52486,15.55968);
      c_disk_stubs->SetLeftMargin(0.03665521);
      c_disk_stubs->SetRightMargin(0.1340206);
      c_disk_stubs->SetTopMargin(0.05121951);
      c_disk_stubs->SetBottomMargin(0.02926829);

      h_disk_stubs->GetXaxis()->SetLabelSize(0.);
      h_disk_stubs->GetXaxis()->SetLabelOffset(999);
      h_disk_stubs->GetYaxis()->SetLabelSize(0.);
      h_disk_stubs->GetYaxis()->SetLabelOffset(999);
      
      h_disk_stubs->Draw("");
      h_ring1_s->Draw("polcolsame");
      h_ring2_s->Draw("polcolsame"); 
      h_ring3_s->Draw("polcolsame");
      h_ring4_s->Draw("polcolsame");
      h_ring5_s->Draw("polcolsame");
      h_ring6_s->Draw("polcolsame");
      h_ring7_s->Draw("polcolsame");
      h_ring8_s->Draw("polcolsame");
      h_ring9_s->Draw("polcolsame");
      h_ring10_s->Draw("polcolsame");
      h_ring11_s->Draw("polcolsame");
      h_ring12_s->Draw("polcolsame");
      h_ring13_s->Draw("polcolsame");
      h_ring14_s->Draw("polcolzsame");
      h_ring15_s->Draw("polcolzsame");
      
      TPaveText *PT_disk_s = new TPaveText(-15.,14.9,3.,16.,"br");

      char txt_disk_s[100];
      if (disk>=7) sprintf(txt_disk_s, "Endcap Disk %d (%sZ side) stub rate (in Stub/Module/BX)", (disk-7)%7+1,"-");
      else         sprintf(txt_disk_s, "Endcap Disk %d (%sZ side) stub rate (in Stub/Module/BX)", disk%7+1,"+");

      PT_disk_s->SetFillColor(0);
      PT_disk_s->SetTextSize(0.02);

      TLatex CMS_disk_s;
      CMS_disk_s.SetTextSize(0.02);
      CMS_disk_s.DrawLatex(9., 15.3, "CMS Preliminary Simulation");     

      TText *ttext_disk_s = PT_disk_s->AddText(txt_disk_s);
      PT_disk_s->Draw();

      c_disk_stubs->Modified();
      c_disk_stubs->Update();
      c_disk_stubs->Write();

      if (PrintEachDisk) {
        char name_disk_s[50];
        sprintf (name_disk_s, "_disk_%d_diskmap_stubs.png", disk+1);
        c_disk_stubs->SaveAs(plotDIR+plotPRE+outputname+name_disk_s);
      }    
      
      c_disk_stubs->Update();


    // ----- Disk Map Clusters ----- //
    TCanvas *c_disk_clusters = new TCanvas("c_disk_clusters","Disk map - Clusters",166,77,947,880);
      c_disk_clusters->Range(-15.23757,-14.89125,18.52486,15.55968);
      c_disk_clusters->SetLeftMargin(0.03665521);
      c_disk_clusters->SetRightMargin(0.1340206);
      c_disk_clusters->SetTopMargin(0.05121951);
      c_disk_clusters->SetBottomMargin(0.02926829);

      h_disk_clusters->GetXaxis()->SetLabelSize(0.);
      h_disk_clusters->GetXaxis()->SetLabelOffset(999);
      h_disk_clusters->GetYaxis()->SetLabelSize(0.);
      h_disk_clusters->GetYaxis()->SetLabelOffset(999);
      
      h_disk_clusters->Draw("");
      h_ring1_c->Draw("polcolsame");
      h_ring2_c->Draw("polcolsame"); 
      h_ring3_c->Draw("polcolsame");
      h_ring4_c->Draw("polcolsame");
      h_ring5_c->Draw("polcolsame");
      h_ring6_c->Draw("polcolsame");
      h_ring7_c->Draw("polcolsame");
      h_ring8_c->Draw("polcolsame");
      h_ring9_c->Draw("polcolsame");
      h_ring10_c->Draw("polcolsame");
      h_ring11_c->Draw("polcolsame");
      h_ring12_c->Draw("polcolsame");
      h_ring13_c->Draw("polcolsame");
      h_ring14_c->Draw("polcolzsame");
      h_ring15_c->Draw("polcolzsame");
      
      TPaveText *PT_disk_c = new TPaveText(-15.,14.9,3.,16.,"br");

      char txt_disk_c[100];
      if (disk>=7) sprintf(txt_disk_c, "Endcap Disk %d (%sZ side) cluster rate (in Cluster/Module/BX)", (disk-7)%7+1,"-");
      else         sprintf(txt_disk_c, "Endcap Disk %d (%sZ side) cluster rate (in Cluster/Module/BX)", disk%7+1,"+");

      PT_disk_c->SetFillColor(0);
      PT_disk_c->SetTextSize(0.02);

      TLatex CMS_disk_c;
      CMS_disk_c.SetTextSize(0.02);
      CMS_disk_c.DrawLatex(9., 15.3, "CMS Preliminary Simulation");     

      TText *ttext_disk_c = PT_disk_c->AddText(txt_disk_c);
      PT_disk_c->Draw();

      c_disk_clusters->Modified();
      c_disk_clusters->Update();
      c_disk_clusters->Write();

      if (PrintEachDisk) {
        char name_disk_c[50];
        sprintf (name_disk_c, "_disk_%d_diskmap_clusters.png", disk+1);
        c_disk_clusters->SaveAs(plotDIR+plotPRE+outputname+name_disk_c);      
      }

      c_disk_clusters->Update();



    // ----- Disk Map Ratio ----- //
    TCanvas *c_disk_ratios = new TCanvas("c_disk_ratios","Disk map - Clusters over Stubs Ratio",166,77,947,880);
      c_disk_ratios->Range(-15.23757,-14.89125,18.52486,15.55968);
      c_disk_ratios->SetLeftMargin(0.03665521);
      c_disk_ratios->SetRightMargin(0.1340206);
      c_disk_ratios->SetTopMargin(0.05121951);
      c_disk_ratios->SetBottomMargin(0.02926829);

      h_disk_ratios->GetXaxis()->SetLabelSize(0.);
      h_disk_ratios->GetXaxis()->SetLabelOffset(999);
      h_disk_ratios->GetYaxis()->SetLabelSize(0.);
      h_disk_ratios->GetYaxis()->SetLabelOffset(999);
      
      h_disk_ratios->Draw("");
      h_ring1_r->Draw("polcolsame");
      h_ring2_r->Draw("polcolsame"); 
      h_ring3_r->Draw("polcolsame");
      h_ring4_r->Draw("polcolsame");
      h_ring5_r->Draw("polcolsame");
      h_ring6_r->Draw("polcolsame");
      h_ring7_r->Draw("polcolsame");
      h_ring8_r->Draw("polcolsame");
      h_ring9_r->Draw("polcolsame");
      h_ring10_r->Draw("polcolsame");
      h_ring11_r->Draw("polcolsame");
      h_ring12_r->Draw("polcolsame");
      h_ring13_r->Draw("polcolsame");
      h_ring14_r->Draw("polcolzsame");
      h_ring15_r->Draw("polcolzsame");

      TPaveText *PT_disk_r = new TPaveText(-15.,14.9,3.,16.,"br");

      char txt_disk_r[100];
      if (disk>=7) sprintf(txt_disk_r, "Endcap Disk %d (%sZ side) cluster/stub ratio", (disk-7)%7+1,"-");
      else         sprintf(txt_disk_r, "Endcap Disk %d (%sZ side) cluster/stub ratio", disk%7+1,"+");

      PT_disk_r->SetFillColor(0);
      PT_disk_r->SetTextSize(0.02);

      TLatex CMS_disk_r;
      CMS_disk_r.SetTextSize(0.02);
      CMS_disk_r.DrawLatex(9., 15.3, "CMS Preliminary Simulation");     

      TText *ttext_disk_r = PT_disk_r->AddText(txt_disk_r);
      PT_disk_r->Draw();

      c_disk_ratios->Modified();
      c_disk_ratios->Update();
      c_disk_ratios->Write();

      if (PrintEachDisk) {
        char name_disk_r[50];
        sprintf (name_disk_r, "_disk_%d_diskmap_ratio.png", disk+1);
        c_disk_ratios->SaveAs(plotDIR+plotPRE+outputname+name_disk_r);      
      }

      c_disk_ratios->Update();
  } // close disk loop


  // -----------------------------------------------------------------------------------------------------------
  // 1D Disk Map
  // ---------------------------------------------------------------------------------------------------------

  // Values are for ALL layers and do not get reset at the top of each layer loop.
  float el_maxval_s = 0;
  float el_maxval_c = 0;
  float el_maxval_r = 0;
  float el_rate_s;
  float el_rate_c;
  float el_rate_r;  

  for (int i=7;i<12;++i) {   
    for (int j=i*100;j<i*100+15;++j) { 
      el_rate_s = count_selr[j]/n_mod_endcap[j-100*i];
      el_rate_c = count_celr[j]/n_mod_endcap[j-100*i];
      el_rate_r = 0;
      if (el_rate_s !=0) el_rate_r = el_rate_c/el_rate_s;

      if (RateInHz) { 
        el_rate_s /= n_mod_surf_endcap[j-100*i];
        el_rate_s *= bit_per_s*40.; // Give the value in MHz/cm2

        el_rate_c /= n_mod_surf_endcap[j-100*i];
        el_rate_c *= bit_per_s*40.; // Give the value in MHz/cm2
      }
      
      if (el_rate_s>el_maxval_s) el_maxval_s = el_rate_s;
      if (el_rate_c>el_maxval_c) el_maxval_c = el_rate_c;
      if (el_rate_r>el_maxval_r) el_maxval_r = el_rate_r;
    }
  }

  TH2F *h_1D_disk_stub    = new TH2F("1D_disk_stub",";r (cm); Stub rate (number of stubs/module/BX)",200,20.,120.,100,0.,1.1*el_maxval_s);
  TH2F *h_1D_disk_cluster = new TH2F("1D_disk_cluster",";r (cm); Cluster rate (number of clusters/module/BX)",200,20.,120.,100,0.,1.1*el_maxval_c);
  TH2F *h_1D_disk_ratio   = new TH2F("1D_disk_ratio",";r (cm); Cluster over stub ratio",200,20.,120.,100,0.,1.1*el_maxval_r); 

  std::vector<TH2F*> h_disk1D_stub_plots; 
  std::vector<TH2F*> h_disk1D_clus_plots; 
  std::vector<TH2F*> h_disk1D_rat_plots; 

  for (int i=0;i<5;++i) {
    TH2F *h_disk1D_stub = new TH2F("","",200,20.,120.,100,0.,1.1*el_maxval_s);
    TH2F *h_disk1D_clus = new TH2F("","",200,20.,120.,100,0.,1.1*el_maxval_c);
    TH2F *h_disk1D_rat  = new TH2F("","",200,20.,120.,100,0.,1.1*el_maxval_r);

    h_disk1D_stub_plots.push_back(h_disk1D_stub);
    h_disk1D_clus_plots.push_back(h_disk1D_clus);
    h_disk1D_rat_plots.push_back(h_disk1D_rat);
  }

  int mv;
  for (int i=7;i<12;++i) {
    for (int j=i*100;j<i*100+15;++j) {
      mv = j%100;
      if (i>=9) mv+=3;
      if (mv>=15) continue;

      el_rate_s = count_selr[j]/n_mod_endcap[mv];
      el_rate_c = count_celr[j]/n_mod_endcap[mv];
      if (el_rate_s!=0) el_rate_r = count_celr[j]/count_selr[j];

      if (RateInHz) { 
        el_rate_s /= n_mod_surf_endcap[mv];
        el_rate_s *= bit_per_s*40.; // Give the value in MHz/cm2

        el_rate_c /= n_mod_surf_endcap[mv];
        el_rate_c *= bit_per_s*40.; // Give the value in MHz/cm2
      }

      if (el_rate_s !=0) h_disk1D_stub_plots[i-7]->Fill(radii[mv],el_rate_s);
      if (el_rate_c !=0) h_disk1D_clus_plots[i-7]->Fill(radii[mv],el_rate_c);
      if (el_rate_r !=0) h_disk1D_rat_plots[i-7] ->Fill(radii[mv],el_rate_r);
    }
  }
  
  // ----- 1D Disk Map - Stub Rate ----- //
  TCanvas *c_1D_disk_stubs = new TCanvas("c_1D_disk_stubs","1D Disk map - Stubs",0,0,640,480);
    c_1D_disk_stubs->Range(-1.,-1.,1.,1.);
    c_1D_disk_stubs->SetGridx();
    c_1D_disk_stubs->SetGridy();

    if (RateInHz)  h_1D_disk_stub->GetYaxis()->SetTitle("Stub rate (in MHz/cm^{2})"); 
    h_1D_disk_stub->GetYaxis()->SetTitleOffset(0.83);
    h_1D_disk_stub->Draw();

    for (int i=0;i<5;++i) {
      h_disk1D_stub_plots[i]->SetMarkerStyle(color[i]);
      h_disk1D_stub_plots[i]->SetMarkerSize(1.4);
      h_disk1D_stub_plots[i]->Draw("same");
      h_disk1D_stub_plots[i]->Write();
    }
    
    TLegend *leg_disk1D_s = new TLegend(0.58,0.63,0.88,0.89);
    leg_disk1D_s->SetTextSize(0.035);
    leg_disk1D_s->SetFillColor(0);
    leg_disk1D_s->AddEntry(h_disk1D_stub_plots[0],"TEDD double-disc 1","p");
    leg_disk1D_s->AddEntry(h_disk1D_stub_plots[1],"TEDD double-disc 2","p");
    leg_disk1D_s->AddEntry(h_disk1D_stub_plots[2],"TEDD double-disc 3","p");
    leg_disk1D_s->AddEntry(h_disk1D_stub_plots[3],"TEDD double-disc 4","p");
    leg_disk1D_s->AddEntry(h_disk1D_stub_plots[4],"TEDD double-disc 5","p");
    leg_disk1D_s->Draw();

    char txt_disk1D[100];
    sprintf (txt_disk1D, "#sqrt{s}=14TeV, %d PU", pu);

    TLatex CMS_disk1D_s;
    CMS_disk1D_s.SetTextSize(0.04);
    CMS_disk1D_s.DrawLatex(30, 1.11*el_maxval_s, "CMS Phase-2 Simulation");

    TLatex TL_disk1D_s;
    TL_disk1D_s.SetTextSize(0.03);
    TL_disk1D_s.SetTextFont(52);
    TL_disk1D_s.DrawLatex(100, 1.11*el_maxval_s, txt_disk1D);

    c_1D_disk_stubs->Modified();
    c_1D_disk_stubs->Update();
    c_1D_disk_stubs->Write();
    
    char name_disk1D_s[100];
    sprintf (name_disk1D_s, "_1D_diskmap_stubs.png");
    if (RateInHz)  sprintf (name_disk1D_s, "_1D_diskmap_stubs_inHz.png");
    c_1D_disk_stubs->SaveAs(plotDIR+plotPRE+outputname+name_disk1D_s);

    if (SavePlotSource) {
      char name2_disk1D_s[100];
      sprintf (name2_disk1D_s, "_1D_diskmap_stubs.C");
      if (RateInHz)  sprintf (name2_disk1D_s, "_1D_diskmap_stubs_inHz.C");
      c_1D_disk_stubs->SaveAs(plotDIR+plotPRE+outputname+name2_disk1D_s);
    }

  // ----- 1D Disk Map - Cluster Rate ----- //
  TCanvas *c_1D_disk_clusters = new TCanvas("c_1D_disk_clusters","1D Disk map - Clusters",0,0,640,480);
    c_1D_disk_clusters->Range(-1.,-1.,1.,1.);
    c_1D_disk_clusters->SetGridx();
    c_1D_disk_clusters->SetGridy();

    if (RateInHz)  h_1D_disk_cluster->GetYaxis()->SetTitle("Cluster rate (in MHz/cm^{2})"); 
    h_1D_disk_cluster->GetYaxis()->SetTitleOffset(0.83);
    h_1D_disk_cluster->Draw();
    
    for (int i=0;i<5;++i) {
      h_disk1D_clus_plots[i]->SetMarkerStyle(color[i]);
      h_disk1D_clus_plots[i]->SetMarkerSize(1.4);
      h_disk1D_clus_plots[i]->Draw("same");
      h_disk1D_clus_plots[i]->Write();
    }
    
    TLegend *leg_disk1D_c = new TLegend(0.58,0.63,0.88,0.89);
    leg_disk1D_c->SetTextSize(0.035);
    leg_disk1D_c->SetFillColor(0);
    leg_disk1D_c->AddEntry(h_disk1D_clus_plots[0],"TEDD double-disc 1","p");
    leg_disk1D_c->AddEntry(h_disk1D_clus_plots[1],"TEDD double-disc 2","p");
    leg_disk1D_c->AddEntry(h_disk1D_clus_plots[2],"TEDD double-disc 3","p");
    leg_disk1D_c->AddEntry(h_disk1D_clus_plots[3],"TEDD double-disc 4","p");
    leg_disk1D_c->AddEntry(h_disk1D_clus_plots[4],"TEDD double-disc 5","p");
    leg_disk1D_c->Draw();

    sprintf (txt_disk1D, "#sqrt{s}=14TeV, %d PU", pu);

    TLatex CMS_disk1D_c;
    CMS_disk1D_c.SetTextSize(0.04);
    CMS_disk1D_c.DrawLatex(30, 1.11*el_maxval_c, "CMS Phase-2 Simulation");

    TLatex TL_disk1D_c;
    TL_disk1D_c.SetTextSize(0.03);
    TL_disk1D_c.SetTextFont(52);
    TL_disk1D_c.DrawLatex(100, 1.11*el_maxval_c, txt_disk1D);

    c_1D_disk_clusters->Modified();
    c_1D_disk_clusters->Update();
    c_1D_disk_clusters->Write();
    
    char name_disk1D_c[100];
    sprintf (name_disk1D_c, "_1D_diskmap_clusters.png");
    if (RateInHz)  sprintf (name_disk1D_c, "_1D_diskmap_clusters_inHz.png");
    c_1D_disk_clusters->SaveAs(plotDIR+plotPRE+outputname+name_disk1D_c);

    if (SavePlotSource) {
      char name2_disk1D_c[100];
      sprintf (name2_disk1D_c, "_1D_diskmap_clusters.C");
      if (RateInHz)  sprintf (name2_disk1D_c, "_1D_diskmap_clusters_inHz.C");
      c_1D_disk_clusters->SaveAs(plotDIR+plotPRE+outputname+name2_disk1D_c);
    }

  // ----- 1D Disk Map - Ratios ----- //
  TCanvas *c_1D_disk_ratios = new TCanvas("c_1D_disk_ratios","1D Disk map - Clusters over Stubs",0,0,640,480);
    c_1D_disk_ratios->Range(-1.,-1.,1.,1.);
    c_1D_disk_ratios->SetGridx();
    c_1D_disk_ratios->SetGridy();

    h_1D_disk_ratio->GetYaxis()->SetTitleOffset(0.83);
    h_1D_disk_ratio->Draw();
    
    for (int i=0;i<5;++i) {
      h_disk1D_rat_plots[i]->SetMarkerStyle(color[i]);
      h_disk1D_rat_plots[i]->SetMarkerSize(1.4);
      h_disk1D_rat_plots[i]->Draw("same");
      h_disk1D_rat_plots[i]->Write();
      h_disk1D_rat_plots[i]->Write();
    }
    
    TLegend *leg_disk1D_r = new TLegend(0.58,0.63,0.88,0.89);
    leg_disk1D_r->SetTextSize(0.035);
    leg_disk1D_r->SetFillColor(0);
    leg_disk1D_r->AddEntry(h_disk1D_rat_plots[0],"TEDD double-disc 1","p");
    leg_disk1D_r->AddEntry(h_disk1D_rat_plots[1],"TEDD double-disc 2","p");
    leg_disk1D_r->AddEntry(h_disk1D_rat_plots[2],"TEDD double-disc 3","p");
    leg_disk1D_r->AddEntry(h_disk1D_rat_plots[3],"TEDD double-disc 4","p");
    leg_disk1D_r->AddEntry(h_disk1D_rat_plots[4],"TEDD double-disc 5","p");
    leg_disk1D_r->Draw();

    sprintf (txt_disk1D, "#sqrt{s}=14TeV, %d PU", pu);

    TLatex CMS_disk1D_r;
    CMS_disk1D_r.SetTextSize(0.04);
    CMS_disk1D_r.DrawLatex(30, 1.11*el_maxval_r, "CMS Phase-2 Simulation");

    TLatex TL_disk1D_r;
    TL_disk1D_r.SetTextSize(0.03);
    TL_disk1D_r.SetTextFont(52);
    TL_disk1D_r.DrawLatex(100, 1.11*el_maxval_r, txt_disk1D);

    c_1D_disk_ratios->Modified();
    c_1D_disk_ratios->Update();
    c_1D_disk_ratios->Write();
    
    char name_disk1D_r[100];
    sprintf (name_disk1D_r, "_1D_diskmap_ratios.png");
    c_1D_disk_ratios->SaveAs(plotDIR+plotPRE+outputname+name_disk1D_r);

    if (SavePlotSource) {
      char name2_disk1D_r[100];
      sprintf (name2_disk1D_r, "_1D_diskmap_ratios.C");
      c_1D_disk_ratios->SaveAs(plotDIR+plotPRE+outputname+name2_disk1D_r);
    }


  // ---------------------------------------------------------------------------------------------------------
  // Printouts

  // -----------------------------------------------------------------------------------------------------------
  // TIB Rate

  if(PrintTIBRate) {
      double TIB_n_lad_barrel[6]  = {18,26,36};
      double TIB_n_mod_barrel[6]  = {31,35,39};

      double TIB_rate   = 0;
      double TIB_rate_p = 0;
      double TIB_rate_p2= 0;
      double TIB_rate_s = 0;
      double TIB_rate_f = 0;
      int TIB_idx;
      
      std::cout.precision(3);

      for (int i=0;i<3;++i) { 
        cout << endl;
        cout << "Layer " << i+1 << endl;
        cout << "        /  TOT  / PRIM>2/ PRIM<2/  SEC  /  FAKE" << endl;
        
        
        for (int j=0;j<TIB_n_mod_barrel[i];++j) {
          cout << "Ring " << j+1 << "  / ";

          TIB_rate = 0;
          TIB_rate_p = 0;
          TIB_rate_p2= 0;
          TIB_rate_s = 0;
          TIB_rate_f = 0;
      
          for (int k=0;k<TIB_n_lad_barrel[i];++k) {
            TIB_idx = 10000*(i) + 100*k + j;

            TIB_rate_p += count_sbrp[TIB_idx]-count_sbrp2[TIB_idx];
            TIB_rate_p2+= count_sbrp2[TIB_idx];
            TIB_rate_s += count_sbrs[TIB_idx];
            TIB_rate_f += count_sbrf[TIB_idx];
          }
     
          float TIB_rate_tot = (TIB_rate_p2+TIB_rate_p+TIB_rate_s+TIB_rate_f);
          float TIB_factor   = 100./TIB_rate_tot;
            
          cout << std::fixed << std::setprecision(1) << std::setw(5) << TIB_rate_tot/TIB_n_lad_barrel[i]
               << " / "<< std::setw(5) << TIB_rate_p2*TIB_factor
               << " / "<< std::setw(5) << TIB_rate_p*TIB_factor
               << " / "<< std::setw(5) << TIB_rate_s*TIB_factor
               << " / "<< std::setw(5) << TIB_rate_f*TIB_factor 
               << endl;
        }
      }
  }
  

  // -----------------------------------------------------------------------------------------------------------
  // TOB Rate
  if(PrintTOBRate) {
      double TOB_n_lad_barrel[6]  = {48,60,78};
      double TOB_n_mod_barrel[6]  = {24,24,24};
      double TOB_n_off_barrel[6]  = {12,12,12};
      
      double TOB_rate   = 0;
      double TOB_rate_p = 0;
      double TOB_rate_p2= 0;
      double TOB_rate_s = 0;
      double TOB_rate_f = 0;
      
      int TOB_idx;
      
      std::cout.precision(3);

      for (int i=0;i<3;++i) { 
        cout << endl;
        cout << "Layer " << i+1 << endl;
        cout << "        /  TOT  / PRIM>2/ PRIM<2/  SEC  /  FAKE" << endl;
        
        
        for (int j=TOB_n_off_barrel[i];j<TOB_n_mod_barrel[i];++j) {
          if (int(j+1-TOB_n_off_barrel[i])<10)
            cout << "Ring " << int(j+1-TOB_n_off_barrel[i]) << "  / ";
          if (int(j+1-TOB_n_off_barrel[i])>=10)
            cout << "Ring " << int(j+1-TOB_n_off_barrel[i]) << " / ";

            TOB_rate   = 0;
            TOB_rate_p = 0;
            TOB_rate_p2= 0;
            TOB_rate_s = 0;
            TOB_rate_f = 0;
      
            for (int k=0;k<TOB_n_lad_barrel[i];++k) {
              TOB_idx = 10000*(i+3) + 100*k + j;

              TOB_rate_p += count_sbrp[TOB_idx]-count_sbrp2[TOB_idx];
              TOB_rate_p2+= count_sbrp2[TOB_idx];
              TOB_rate_s += count_sbrs[TOB_idx];
              TOB_rate_f += count_sbrf[TOB_idx];
            }

            float TOB_rate_tot = (TOB_rate_p2+TOB_rate_p+TOB_rate_s+TOB_rate_f);
            float TOB_factor   = 100./TOB_rate_tot;
            
            cout << std::fixed << std::setprecision(1) << std::setw(5) << TOB_rate_tot/TOB_n_lad_barrel[i]
                  << " / "<< std::setw(5) << TOB_rate_p2*TOB_factor
                  << " / "<< std::setw(5) << TOB_rate_p*TOB_factor
                  << " / "<< std::setw(5) << TOB_rate_s*TOB_factor
                  << " / "<< std::setw(5) << TOB_rate_f*TOB_factor 
                  << endl;

        }
      }
  }


  // -----------------------------------------------------------------------------------------------------------
  // ENDCAP Rate

  if(PrintEndcapRate) {
      double ENDCAP_rate   = 0;
      double ENDCAP_rate_p = 0;
      double ENDCAP_rate_p2= 0;
      double ENDCAP_rate_s = 0;
      double ENDCAP_rate_f = 0;
      int ENDCAP_idx;
      
      std::cout.precision(3);

      for (int i=0;i<5;++i) { 
        cout << endl;
        cout << "DISK " << i+1 << endl;
        cout << "        /  TOT  / PRIM>2/ PRIM<2/  SEC  /  FAKE" << endl;
        
        
        for (int j=0;j<15;++j) {
          if (j<10)  cout << "Ring " << j << "  / ";
          if (j>=10) cout << "Ring " << j << " / ";

          ENDCAP_rate   = 0;
          ENDCAP_rate_p = 0;
          ENDCAP_rate_p2= 0;
          ENDCAP_rate_s = 0;
          ENDCAP_rate_f = 0;
      
          for (int k=0;k<n_lad_endcap[j];++k) {
            ENDCAP_idx = 10000*i + 100*j + k;

            ENDCAP_rate_p += count_serp[ENDCAP_idx]-count_serp2[ENDCAP_idx];
            ENDCAP_rate_p2+= count_serp2[ENDCAP_idx];
            ENDCAP_rate_s += count_sers[ENDCAP_idx];
            ENDCAP_rate_f += count_serf[ENDCAP_idx];
          }
            
            float ENDCAP_rate_tot = (ENDCAP_rate_p2+ENDCAP_rate_p+ENDCAP_rate_s+ENDCAP_rate_f);
            float ENDCAP_factor   = 100./ENDCAP_rate_tot;
            
            cout  << std::fixed << std::setprecision(1) << std::setw(5) << ENDCAP_rate_tot/n_lad_endcap[j]
                  << " / "<< std::setw(5) << ENDCAP_rate_p2*ENDCAP_factor
                  << " / "<< std::setw(5) << ENDCAP_rate_p*ENDCAP_factor
                  << " / "<< std::setw(5) << ENDCAP_rate_s*ENDCAP_factor
                  << " / "<< std::setw(5) << ENDCAP_rate_f*ENDCAP_factor 
                  << endl;
            
        }
      }

  }




} // close RatesPlots



// ---------------------------------------------------------------------------------------------------------
// Other Functions
// ---------------------------------------------------------------------------------------------------------

  // -----------------------------------------------------------------------------------------------------------
  // Fill Histograms
  void fill_histo(TH2F *hist,float val[], double nlad, double dec, int disk, int ring,float max) {
    float scale = 8*atan(1.)/nlad;

    float phi_sec = 0.;

    for (int i=0;i<int(nlad);++i) {
      phi_sec = (i+.5-dec)*scale;

      if (phi_sec<0.) phi_sec+=8*atan(1.); 

      int idx = 10000*disk + 100*ring + i;
      hist->Fill(phi_sec,0.5+ring,val[idx]);
    }
   
    hist->Fill(0.5*scale,25.5,0);
    hist->Fill(1.5*scale,25.5,max);
  }


  // ---------------------------------------------------------------------------------------------------------
  // Plot Style
  void SetPlotStyle() 
  {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasBorderSize(2);

    gStyle->SetFrameBorderMode(0);

    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetTitleFont(42,"xyz");

    gStyle->SetLabelSize(0.035,"xyz");
    gStyle->SetTitleSize(0.035,"xyz");

    gStyle->SetLabelOffset(0.004,"y");

    gStyle->SetPadLeftMargin(0.07692308);
    gStyle->SetPadTopMargin(0.07124352);

    TGaxis::SetMaxDigits(3);

  }

  // ---------------------------------------------------------------------------------------------------------
  // Color Gradient
  void SetColorTable() 
  {
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
