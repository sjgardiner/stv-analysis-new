//#define STV_DEBUG

// Standard library includes
#include <iomanip>
#include <iostream>
#include <sstream>

// ROOT includes
#include "TCanvas.h"
#include "TLegend.h"

// STV analysis includes
#include "DAgostiniUnfolder.hh"
#include "FiducialVolume.hh"
#include "MCC9SystematicsCalculator.hh"
#include "SliceBinning.hh"
#include "SliceHistogram.hh"
#include "WienerSVDUnfolder.hh"
#include "FilePropertiesManager.hh"

// Calculates the Mahalanobis distance for the NuWro alt model
// and GENIE multisim variations

void CalcMahalanobisDistance(string inputfile) {

  gSystem->Exec("mkdir -p Plots/");

  auto& fpm = FilePropertiesManager::Instance();
  fpm.load_file_properties( "file_config.txt" );

  // Do the systematics calculations in preparation for unfolding
  auto* mcc9 = new MCC9SystematicsCalculator(
      inputfile,
      "systcalc.conf" );
  auto& syst = *mcc9;

  std::vector<NFT> detvar = {
    NFT::kDetVarMCLYatten, // light-yield attenuation
    NFT::kDetVarMCLYdown, // light-yield down
    NFT::kDetVarMCLYrayl, // light-yield Rayleigh scattering
    NFT::kDetVarMCRecomb2, // light-yield recombination 2
    NFT::kDetVarMCSCE, // space charge effect
    NFT::kDetVarMCWMAngleXZ, // wireMod angle XZ
    NFT::kDetVarMCWMAngleYZ, // wireMod angle YZ
    NFT::kDetVarMCWMdEdx, // wireMod dE/dx
    NFT::kDetVarMCWMX, // wireMod X
    NFT::kDetVarMCWMYZ // wireMod YZ
  };

  std::vector<std::string> reint = {"weight_reint_all","weight_flux_all"};
  
  syst.prepare_mahalanobis_dist_calc(reint,detvar);

  std::unique_ptr<Universe> NuWro = syst.alt_cv_universes_.at(NFT::kAltCVMC)->clone();
  std::pair<int,double> md2 = syst.get_mahalanobis_dist(NuWro);
  std::cout << "NuWro Alt CV: " << md2.first << "  " << md2.second << std::endl;
  std::vector<double> NuWro_X,NuWro_Y;
  NuWro_X = {md2.second/md2.first,md2.second/md2.first};
  NuWro_Y = {0.0,1000.0};
  TGraph* g_NuWro = new TGraph(NuWro_X.size(),&(NuWro_X[0]),&(NuWro_Y[0])); 

  TH1D* h_GENIE_vars = new TH1D("h_GENIE_vars",";Squared MD/dof;Variations",50,0.0,2.0);

  for(size_t i_u=0;i_u<syst.rw_universes_.at("weight_All_UBGenie").size();i_u++){
    std::unique_ptr<Universe> Alt = syst.rw_universes_.at("weight_All_UBGenie").at(i_u)->clone(); 
    md2 = syst.get_mahalanobis_dist(Alt);
    h_GENIE_vars->Fill(md2.second/md2.first);
    std::cout << "Genie Univ " << i_u << ": " << md2.first << "  " << md2.second << std::endl;
  }

 // Draw a plot of the test statistic values
  
 TCanvas* c = new TCanvas("c","c");
 h_GENIE_vars->Draw("HIST");
 g_NuWro->Draw("L same"); 
 
}
