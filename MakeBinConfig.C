#include <iostream>
#include <fstream>

#include "Misc.hh"

using std::string;
using std::vector;

// Script to make simple bin configs.
// 
// Setup:
// label = a string to distinguish between binning schemes 
// signal_def = the cut(s) you apply to define the signal
// sel_def = the cut(s) you apply to apply the selection
// true_variables = the names of the true variables you want to draw
// reco_variables = the names of the reco variables you want to draw
// true_binning = the edges of the bins of the true variables
// reco_binning = the edges of the bins of the reco variables 

void MakeBinConfig(){

  string signal_def = "mc_is_signal";
  string sel_def = "sel_CCNp0pi";

  vector<string> true_variables,reco_variables;
  vector<vector<double>> true_binning,reco_binning;

  // Neutrino energy vs delta pt

  true_variables = {"mc_nu_energy","mc_delta_pT"};
  reco_variables = {"reco_neutrino_energy","delta_pT"};

  true_binning = {
    {0.2,0.4,0.6,0.8,1.0,1.2,1.4},
    {0.0,0.2,0.4,0.6,0.8,1.0}
  };

  reco_binning = {
    {0.2,0.4,0.6,0.8,1.0,1.2,1.4},
    {0.0,0.2,0.4,0.6,0.8,1.0}
  };

  MakeBinSliceConfig("NeutrinoEnergy_DeltaPt",signal_def,sel_def,true_variables,reco_variables,true_binning,reco_binning);

  // Neutrino energy vs delta alphat

  true_variables = {"mc_nu_energy","mc_delta_alphaT"};
  reco_variables = {"reco_neutrino_energy","delta_alphaT"};

  true_binning = {
    {0.2,0.4,0.6,0.8,1.0,1.2,1.4},
    {0.0,0.314,0.628,0.942,1.256,1.57,1.884,2.198,2.512,2.826,3.14}
  };

  reco_binning = {
    {0.2,0.4,0.6,0.8,1.0,1.2,1.4},
    {0.0,0.314,0.628,0.942,1.256,1.57,1.884,2.198,2.512,2.826,3.14}
  };

  MakeBinSliceConfig("NeutrinoEnergy_DeltaAlphat",signal_def,sel_def,true_variables,reco_variables,true_binning,reco_binning);
  
  // Neutrino energy vs delta phit

  true_variables = {"mc_nu_energy","mc_delta_phiT"};
  reco_variables = {"reco_neutrino_energy","delta_phiT"};

  true_binning = {
    {0.2,0.4,0.6,0.8,1.0,1.2,1.4},
    {0.0,0.314,0.628,0.942,1.256,1.57,1.884,2.198,2.512,2.826,3.14}
  };

  reco_binning = {
    {0.2,0.4,0.6,0.8,1.0,1.2,1.4},
    {0.0,0.314,0.628,0.942,1.256,1.57,1.884,2.198,2.512,2.826,3.14}
  };

  MakeBinSliceConfig("NeutrinoEnergy_DeltaPhit",signal_def,sel_def,true_variables,reco_variables,true_binning,reco_binning);

  // Neutrino energy vs delta pL

  true_variables = {"mc_nu_energy","mc_delta_pL"};
  reco_variables = {"reco_neutrino_energy","delta_pL"};

  true_binning = {
    {0.2,0.4,0.6,0.8,1.0,1.2,1.4},
    {-1.0,-0.85,-0.7,-0.55,-0.40,-0.25,-0.10,0.05,0.20,0.35,0.50}
  };

  reco_binning = {
    {0.2,0.4,0.6,0.8,1.0,1.2,1.4},
    {-1.0,-0.85,-0.7,-0.55,-0.40,-0.25,-0.10,0.05,0.20,0.35,0.50}
  };

  MakeBinSliceConfig("NeutrinoEnergy_DeltaPL",signal_def,sel_def,true_variables,reco_variables,true_binning,reco_binning);

}
