#include <iostream>
#include <fstream>

#include "BinConfigFuncs.hh"

using std::string;
using std::vector;

// Script to make simple bin configs.

void MakeBinConfig(){

  string signal_def = "mc_is_signal";
  string sel_def = "sel_CCNp0pi";
  string testfile = "/uboone/data/users/cthorpe/STV/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";

  vector<double> nu_energy_binning = GetBinning1D(testfile,"",signal_def,sel_def,"mc_nu_energy","reco_neutrino_energy");
  vector<double> delta_pT_binning = GetBinning1D(testfile,"",signal_def,sel_def,"mc_delta_pT","delta_pT");
  vector<double> delta_alphaT_binning = GetBinning1D(testfile,"",signal_def,sel_def,"mc_delta_alphaT","delta_alphaT");
  vector<double> delta_phiT_binning = GetBinning1D(testfile,"",signal_def,sel_def,"mc_delta_phiT","delta_phiT");
  vector<double> delta_pL_binning = GetBinning1D(testfile,"",signal_def,sel_def,"mc_delta_pL","delta_pL");
  vector<double> muon_angle_binning = GetBinning1D(testfile,"",signal_def,sel_def,"mc_p3_mu.CosTheta()","p3_mu.CosTheta()");
  vector<double> muon_mom_binning = GetBinning1D(testfile,"",signal_def,sel_def,"mc_p3_mu.Mag()","p3_mu.Mag()");

  // Various 1D distributions
  MakeBinSliceConfig1D("NeutrinoEnergy",signal_def,sel_def,"mc_nu_energy","reco_neutrino_energy",nu_energy_binning,nu_energy_binning);
  MakeBinSliceConfig1D("DeltaPT",signal_def,sel_def,"mc_delta_pT","delta_pT",delta_pT_binning,delta_pT_binning);
  MakeBinSliceConfig1D("DeltaAlphaT",signal_def,sel_def,"mc_delta_alphaT","delta_alphaT",delta_alphaT_binning,delta_alphaT_binning);
  MakeBinSliceConfig1D("DeltaPhiT",signal_def,sel_def,"mc_delta_phiT","delta_phiT",delta_phiT_binning,delta_phiT_binning);
  MakeBinSliceConfig1D("DeltaPL",signal_def,sel_def,"mc_delta_pL","delta_pL",delta_pL_binning,delta_pL_binning);
  MakeBinSliceConfig1D("MuonAngle",signal_def,sel_def,"mc_p3_mu.CosTheta()","p3_mu.CosTheta()",muon_angle_binning,muon_angle_binning);
  MakeBinSliceConfig1D("MuonMom",signal_def,sel_def,"mc_p3_mu.Mag()","p3_mu.Mag()",muon_mom_binning,muon_mom_binning);

  // Lots of 2D distributions

  // neutrino energy in bins of delta pt
  std::vector<std::pair<double,std::vector<double>>> neutrino_energy_delta_pT_binning = GetBinning2D(testfile,"",signal_def,sel_def,{"mc_nu_energy","mc_delta_pT"},{"reco_neutrino_energy","delta_pT"});
  MakeBinSliceConfig2D("NeutrinoEnergy_DeltaPT",signal_def,sel_def,{"mc_nu_energy","mc_delta_pT"},{"reco_neutrino_energy","delta_pT"},neutrino_energy_delta_pT_binning,neutrino_energy_delta_pT_binning);

  // neutrino energy in bins of delta alphat
  std::vector<std::pair<double,std::vector<double>>> neutrino_energy_delta_alphaT_binning = GetBinning2D(testfile,"",signal_def,sel_def,{"mc_nu_energy","mc_delta_alphaT"},{"reco_neutrino_energy","delta_alphaT"});
  MakeBinSliceConfig2D("NeutrinoEnergy_DeltaAlphaT",signal_def,sel_def,{"mc_nu_energy","mc_delta_alphaT"},{"reco_neutrino_energy","delta_alphaT"},neutrino_energy_delta_alphaT_binning,neutrino_energy_delta_alphaT_binning);

  // neutrino energy in bins of delta phit
  std::vector<std::pair<double,std::vector<double>>> neutrino_energy_delta_phiT_binning = GetBinning2D(testfile,"",signal_def,sel_def,{"mc_nu_energy","mc_delta_phiT"},{"reco_neutrino_energy","delta_phiT"});
  MakeBinSliceConfig2D("NeutrinoEnergy_DeltaPhiT",signal_def,sel_def,{"mc_nu_energy","mc_delta_phiT"},{"reco_neutrino_energy","delta_phiT"},neutrino_energy_delta_phiT_binning,neutrino_energy_delta_phiT_binning);

  // neutrino energy in bins of delta pt
  std::vector<std::pair<double,std::vector<double>>> neutrino_energy_delta_pL_binning = GetBinning2D(testfile,"",signal_def,sel_def,{"mc_nu_energy","mc_delta_pL"},{"reco_neutrino_energy","delta_pL"});
  MakeBinSliceConfig2D("NeutrinoEnergy_DeltaPL",signal_def,sel_def,{"mc_nu_energy","mc_delta_pL"},{"reco_neutrino_energy","delta_pL"},neutrino_energy_delta_pL_binning,neutrino_energy_delta_pL_binning);

  // neutrino energy in bins of muon angle
  std::vector<std::pair<double,std::vector<double>>> neutrino_energy_muon_angle_binning = GetBinning2D(testfile,"",signal_def,sel_def,{"mc_nu_energy","mc_p3_mu.CosTheta()"},{"reco_neutrino_energy","p3_mu.CosTheta()"});
  MakeBinSliceConfig2D("NeutrinoEnergy_MuonAngle",signal_def,sel_def,{"mc_nu_energy","mc_p3_mu.CosTheta()"},{"reco_neutrino_energy","p3_mu.CosTheta()"},neutrino_energy_muon_angle_binning,neutrino_energy_muon_angle_binning);

  // neutrino energy in bins of muon momentum
  std::vector<std::pair<double,std::vector<double>>> neutrino_energy_muon_mom_binning = GetBinning2D(testfile,"",signal_def,sel_def,{"mc_nu_energy","mc_p3_mu.Mag()"},{"reco_neutrino_energy","p3_mu.Mag()"});
  MakeBinSliceConfig2D("NeutrinoEnergy_MuonMom",signal_def,sel_def,{"mc_nu_energy","mc_p3_mu.Mag()"},{"reco_neutrino_energy","p3_mu.Mag()"},neutrino_energy_muon_mom_binning,neutrino_energy_muon_mom_binning);

  // muon angle in bins of delta pT 
  std::vector<std::pair<double,std::vector<double>>> muon_angle_delta_pT_binning = GetBinning2D(testfile,"",signal_def,sel_def,{"mc_p3_mu.CosTheta()","mc_delta_pT"},{"p3_mu.CosTheta()","delta_pT"});
  MakeBinSliceConfig2D("MuonAngle_DeltaPT",signal_def,sel_def,{"mc_p3_mu.CosTheta()","mc_delta_pT"},{"p3_mu.CosTheta()","delta_pT"},muon_angle_delta_pT_binning,muon_angle_delta_pT_binning);

  // muon angle in bins of delta alphaT 
  std::vector<std::pair<double,std::vector<double>>> muon_angle_delta_alphaT_binning = GetBinning2D(testfile,"",signal_def,sel_def,{"mc_p3_mu.CosTheta()","mc_delta_alphaT"},{"p3_mu.CosTheta()","delta_alphaT"});
  MakeBinSliceConfig2D("MuonAngle_DeltaAlphaT",signal_def,sel_def,{"mc_p3_mu.CosTheta()","mc_delta_alphaT"},{"p3_mu.CosTheta()","delta_alphaT"},muon_angle_delta_alphaT_binning,muon_angle_delta_alphaT_binning);

  // muon angle in bins of delta phiT 
  std::vector<std::pair<double,std::vector<double>>> muon_angle_delta_phiT_binning = GetBinning2D(testfile,"",signal_def,sel_def,{"mc_p3_mu.CosTheta()","mc_delta_phiT"},{"p3_mu.CosTheta()","delta_phiT"});
  MakeBinSliceConfig2D("MuonAngle_DeltaPhiT",signal_def,sel_def,{"mc_p3_mu.CosTheta()","mc_delta_phiT"},{"p3_mu.CosTheta()","delta_phiT"},muon_angle_delta_phiT_binning,muon_angle_delta_phiT_binning);

  // muon angle in bins of delta PL 
  std::vector<std::pair<double,std::vector<double>>> muon_angle_delta_pL_binning = GetBinning2D(testfile,"",signal_def,sel_def,{"mc_p3_mu.CosTheta()","mc_delta_pL"},{"p3_mu.CosTheta()","delta_pL"});
  MakeBinSliceConfig2D("MuonAngle_DeltaPL",signal_def,sel_def,{"mc_p3_mu.CosTheta()","mc_delta_pL"},{"p3_mu.CosTheta()","delta_pL"},muon_angle_delta_pL_binning,muon_angle_delta_pL_binning);

  // muon mom in bins of delta pT 
  std::vector<std::pair<double,std::vector<double>>> muon_mom_delta_pT_binning = GetBinning2D(testfile,"",signal_def,sel_def,{"mc_p3_mu.Mag()","mc_delta_pT"},{"p3_mu.Mag()","delta_pT"});
  MakeBinSliceConfig2D("MuonMom_DeltaPT",signal_def,sel_def,{"mc_p3_mu.Mag()","mc_delta_pT"},{"p3_mu.Mag()","delta_pT"},muon_mom_delta_pT_binning,muon_mom_delta_pT_binning);

  // muon mom in bins of delta alphaT 
  std::vector<std::pair<double,std::vector<double>>> muon_mom_delta_alphaT_binning = GetBinning2D(testfile,"",signal_def,sel_def,{"mc_p3_mu.Mag()","mc_delta_alphaT"},{"p3_mu.Mag()","delta_alphaT"});
  MakeBinSliceConfig2D("MuonMom_DeltaAlphaT",signal_def,sel_def,{"mc_p3_mu.Mag()","mc_delta_alphaT"},{"p3_mu.Mag()","delta_alphaT"},muon_mom_delta_alphaT_binning,muon_mom_delta_alphaT_binning);

  // muon mom in bins of delta phiT 
  std::vector<std::pair<double,std::vector<double>>> muon_mom_delta_phiT_binning = GetBinning2D(testfile,"",signal_def,sel_def,{"mc_p3_mu.Mag()","mc_delta_phiT"},{"p3_mu.Mag()","delta_phiT"});
  MakeBinSliceConfig2D("MuonMom_DeltaPhiT",signal_def,sel_def,{"mc_p3_mu.Mag()","mc_delta_phiT"},{"p3_mu.Mag()","delta_phiT"},muon_mom_delta_phiT_binning,muon_mom_delta_phiT_binning);

  // muon mom in bins of delta PL 
  std::vector<std::pair<double,std::vector<double>>> muon_mom_delta_pL_binning = GetBinning2D(testfile,"",signal_def,sel_def,{"mc_p3_mu.Mag()","mc_delta_pL"},{"p3_mu.Mag()","delta_pL"});
  MakeBinSliceConfig2D("MuonMom_DeltaPL",signal_def,sel_def,{"mc_p3_mu.Mag()","mc_delta_pL"},{"p3_mu.Mag()","delta_pL"},muon_mom_delta_pL_binning,muon_mom_delta_pL_binning);

  // muon mom in bins of muon angle 
  std::vector<std::pair<double,std::vector<double>>> muon_mom_muon_angle_binning = GetBinning2D(testfile,"",signal_def,sel_def,{"mc_p3_mu.Mag()","mc_p3_mu.CosTheta()"},{"p3_mu.Mag()","p3_mu.CosTheta()"});
  MakeBinSliceConfig2D("MuonMom_MuonAngle",signal_def,sel_def,{"mc_p3_mu.Mag()","mc_p3_mu.CosTheta()"},{"p3_mu.Mag()","p3_mu.CosTheta()"},muon_mom_muon_angle_binning,muon_mom_muon_angle_binning);

}
