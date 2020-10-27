// Analysis macro for use in the CCNp0pi single transverse variable analysis
// Designed for use with the PeLEE group's "searchingfornues" ntuples
//
// Updated 4 September 2020
// Steven Gardiner <gardiner@fnal.gov>

// Standard library includes
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// ROOT includes
#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"

// Helper function that avoids NaNs when taking square roots of negative
// numbers
double real_sqrt( double x ) {
  if ( x < 0. ) return 0.;
  else return std::sqrt( x );
}

// A few helpful dummy constants
constexpr float BOGUS = 9999.;
constexpr int BOGUS_INT = 9999;
constexpr int BOGUS_INDEX = -1;
constexpr float LOW_FLOAT = -1e30;
constexpr float DEFAULT_WEIGHT = 1.;

// Integer representation of CC versus NC for the ccnc branch
constexpr int CHARGED_CURRENT = 0;
constexpr int NEUTRAL_CURRENT = 1;

// Useful PDG codes
constexpr int ELECTRON_NEUTRINO = 12;
constexpr int MUON = 13;
constexpr int MUON_NEUTRINO = 14;
constexpr int TAU_ANTINEUTRINO = -16;
constexpr int PROTON = 2212;
constexpr int PI_ZERO = 111;
constexpr int PI_PLUS = 211;

// Values of parameters to use in analysis cuts
constexpr float PROTON_CHI2_CUT = 88.;
constexpr unsigned int HITS_Y_CUT = 5u;
constexpr float LEAD_P_MOM_CUT = 0.300; // GeV/c
constexpr float MUON_MOM_CUT = 0.150; // GeV/c
constexpr float CHARGED_PI_MOM_CUT = 0.070; // GeV/c

constexpr float TOPO_SCORE_CUT = 0.1;
constexpr float COSMIC_IP_CUT = 10.; // cm

constexpr float MUON_TRACK_SCORE_CUT = 0.8;
constexpr float MUON_VTX_DISTANCE_CUT = 4.; // cm
constexpr float MUON_LENGTH_CUT = 10.; // cm
constexpr float MUON_PID_CUT = 0.2;

constexpr float TRACK_SCORE_CUT = 0.5;

// Boundaries of the neutrino vertex fiducial volume (cm)
// This is handled the same way for reco (in the NuCCanalyzer)
// and in MC (herein)
double FV_X_MIN =   10.;
double FV_X_MAX =  246.35;

double FV_Y_MIN = -106.5;
double FV_Y_MAX =  106.5;

double FV_Z_MIN =   10.;
double FV_Z_MAX =  968.8;

// Boundaries of the proton containment volume (used in reco only) in cm
double PCV_X_MIN =   10.;
double PCV_X_MAX =  245.;

double PCV_Y_MIN = -100.;
double PCV_Y_MAX =  100.;

double PCV_Z_MIN =   10.;
double PCV_Z_MAX = 1030.;

// Mass values from GENIE v3.0.6
constexpr double TARGET_MASS = 37.215526; // 40Ar, GeV
constexpr double NEUTRON_MASS = 0.93956541; // GeV
constexpr double PROTON_MASS = 0.93827208; // GeV
constexpr double BINDING_ENERGY = 0.0295; // 40Ar, GeV
constexpr double MUON_MASS = 0.10565837; // GeV
constexpr double PI_PLUS_MASS = 0.13957000; // GeV

// Enum used to label event categories of interest for analysis plots
enum EventCategory {

  kUnknown = 0,
  kSignalCCQE = 1,
  kSignalCCMEC = 2,
  kSignalCCRES = 3,
  kSignalOther = 4,
  kOOFV = 5, // True neutrino vertex is outside of fiducial volume
  kNC = 6,
  // True numu CC event with at least one final-state pion above threshold
  kNuMuCCNpi = 7,
  // True nue CC event
  kNuECC = 8,
  kOther = 9
  // True numu CC event with zero protons above threshold
  //kNuMuCC0p,

};

// Class used to hold information from the searchingfornues TTree branches and
// process it for our analysis
class AnalysisEvent {

  public:

    AnalysisEvent() {}

    ~AnalysisEvent() {}

    EventCategory categorize_event();
    void apply_selection();
    void apply_numu_CC_selection();
    void find_muon_candidate();
    void find_lead_p_candidate();
    void compute_observables();
    void compute_mc_truth_observables();

    // Event scores needed for numu CC selection
    float topological_score_ = BOGUS;
    float cosmic_impact_parameter_ = BOGUS;

    // Reco PDG code of the neutrino candidate
    int nu_pdg_ = BOGUS_INT;

    // Reco neutrino vertex coordinates (cm)
    float nu_vx_ = BOGUS;
    float nu_vy_ = BOGUS;
    float nu_vz_ = BOGUS;

    // Reconstructed object counts
    int num_pf_particles_ = BOGUS_INT;
    int num_tracks_ = BOGUS_INT;
    int num_showers_ = BOGUS_INT;

    // PFParticle properties
    std::vector<unsigned int>* pfp_generation_ = nullptr;
    std::vector<float>* pfp_track_score_ = nullptr;
    // Reco PDG code assigned by Pandora
    std::vector<int>* pfp_reco_pdg_ = nullptr;
    // True PDG code found using the backtracker
    std::vector<int>* pfp_true_pdg_ = nullptr;
    // Number of hits on the collection plane
    std::vector<int>* pfp_hitsY_ = nullptr;

    // Shower properties
    std::vector<unsigned long>* shower_pfp_id_ = nullptr;
    std::vector<float>* shower_startx_ = nullptr;
    std::vector<float>* shower_starty_ = nullptr;
    std::vector<float>* shower_startz_ = nullptr;
    std::vector<float>* shower_start_distance_ = nullptr;

    // Track properties
    std::vector<unsigned long>* track_pfp_id_ = nullptr;
    std::vector<float>* track_length_ = nullptr;
    std::vector<float>* track_startx_ = nullptr;
    std::vector<float>* track_starty_ = nullptr;
    std::vector<float>* track_startz_ = nullptr;
    std::vector<float>* track_start_distance_ = nullptr;
    std::vector<float>* track_endx_ = nullptr;
    std::vector<float>* track_endy_ = nullptr;
    std::vector<float>* track_endz_ = nullptr;
    std::vector<float>* track_dirx_ = nullptr;
    std::vector<float>* track_diry_ = nullptr;
    std::vector<float>* track_dirz_ = nullptr;

    // Proton *kinetic* energy using range-based momentum reconstruction
    std::vector<float>* track_kinetic_energy_p_ = nullptr;

    std::vector<float>* track_range_mom_mu_ = nullptr;
    std::vector<float>* track_mcs_mom_mu_ = nullptr;
    std::vector<float>* track_chi2_proton_ = nullptr;
    std::vector<float>* track_llr_pid_score_ = nullptr;

    // True neutrino PDG code
    int mc_nu_pdg_ = BOGUS_INT;

    // True neutrino vertex coordinates (cm)
    float mc_nu_vx_ = BOGUS;
    float mc_nu_vy_ = BOGUS;
    float mc_nu_vz_ = BOGUS;

    // True neutrino 4-momentum
    float mc_nu_energy_ = BOGUS;

    // Whether the event is CC (0) or NC (1)
    int mc_nu_ccnc_ = false;

    // Interaction mode (QE, MEC, etc.)
    int mc_nu_interaction_type_ = BOGUS_INT;

    // Final-state particle PDG codes and energies (post-FSIs)
    std::vector<int>* mc_nu_daughter_pdg_ = nullptr;
    std::vector<float>* mc_nu_daughter_energy_ = nullptr;
    std::vector<float>* mc_nu_daughter_px_ = nullptr;
    std::vector<float>* mc_nu_daughter_py_ = nullptr;
    std::vector<float>* mc_nu_daughter_pz_ = nullptr;

    // GENIE weights
    float spline_weight_ = DEFAULT_WEIGHT;
    float tuned_cv_weight_ = DEFAULT_WEIGHT;

    // Signal definition requirements
    bool is_mc_ = false;
    bool mc_neutrino_is_numu_ = false;
    bool mc_vertex_in_FV_ = false;
    bool mc_pmu_above_threshold_ = false;
    bool mc_has_p_above_threshold_ = false;
    bool mc_no_fs_pi0_ = false;
    bool mc_no_charged_pi_above_threshold_ = false;
    // Intersection of all of these requirements
    bool mc_is_signal_ = false;

    EventCategory category_ = kUnknown;

    // **** Reco selection requirements ****

    // Whether the event passed the upstream numu CC selection (by Wouter)
    bool sel_nu_mu_cc_ = false;
    // False if at least one generation == 2 shower was reconstructed
    bool sel_no_reco_showers_ = false;
    // True if exactly one generation == 2 muon candidate was identified
    bool sel_has_muon_candidate_ = false;
    // Whether the muon candidate has a reco momentum above threshold
    bool sel_muon_above_threshold_ = false;
    // Whether at least one generation == 2 reco track exists that is not the
    // muon candidate
    bool sel_has_p_candidate_ = false;
    // Whether all proton candidates with at least HITS_Y_CUT collection plane
    // hits pass the proton PID cut (chi^2_proton <= PROTON_CHI2_CUT)
    bool sel_passed_proton_pid_cut_ = false;
    // Whether all proton candidates have track end coordinates that lie within
    // the "containment volume"
    bool sel_protons_contained_ = false;
    // Whether the leading proton candidate has at least HITS_Y_CUT collection
    // plane hits
    bool sel_lead_p_passed_hits_cut_ = false;
    // Whether the leading proton candidate has a range-based reco momentum
    // above LEAD_P_MOM_CUT
    bool sel_lead_p_above_mom_cut_ = false;
    // Intersection of all of the above requirements
    bool sel_CCNp0pi_ = false;

    // Muon and leading proton candidate indices (BOGUS_INDEX if not present)
    // in the reco track arrays
    int muon_candidate_idx_ = BOGUS_INDEX;
    int lead_p_candidate_idx_ = BOGUS_INDEX;

    // ** Reconstructed observables **

    // 3-momenta
    TVector3 p3_mu_;
    TVector3 p3_lead_p_;

    // Dummy pointers to assist in setting branch addresses
    TVector3* p3_mu_ptr_ = &p3_mu_;
    TVector3* p3_lead_p_ptr_ = &p3_lead_p_;

    // Reco STVs
    float delta_pT_ = BOGUS;
    float delta_phiT_ = BOGUS;
    float delta_alphaT_ = BOGUS;
    float delta_pL_ = BOGUS;
    float pn_ = BOGUS;

    // ** MC truth observables **
    // These are loaded for signal events whenever we have a GHEP event to use

    // 3-momenta
    TVector3 mc_p3_mu_;
    TVector3 mc_p3_lead_p_;

    // Dummy pointers to assist in setting branch addresses
    TVector3* mc_p3_mu_ptr_ = &mc_p3_mu_;
    TVector3* mc_p3_lead_p_ptr_ = &mc_p3_lead_p_;

    // MC truth STVs
    float mc_delta_pT_ = BOGUS;
    float mc_delta_phiT_ = BOGUS;
    float mc_delta_alphaT_ = BOGUS;
    float mc_delta_pL_ = BOGUS;
    float mc_pn_ = BOGUS;

    bool point_inside_FV( float x, float y, float z ) {
      bool x_inside_FV = ( FV_X_MIN < x ) && ( x < FV_X_MAX );
      bool y_inside_FV = ( FV_Y_MIN < y ) && ( y < FV_Y_MAX );
      bool z_inside_FV = ( FV_Z_MIN < z ) && ( z < FV_Z_MAX );
      return ( x_inside_FV && y_inside_FV && z_inside_FV );
    }

    bool reco_vertex_inside_FV() {
      return point_inside_FV( nu_vx_, nu_vy_, nu_vz_ );
    }

    bool mc_vertex_inside_FV() {
      return point_inside_FV( mc_nu_vx_, mc_nu_vy_, mc_nu_vz_ );
    }

    bool in_proton_containment_vol( float x, float y, float z ) {
      bool x_inside_PCV = ( PCV_X_MIN < x ) && ( x < PCV_X_MAX );
      bool y_inside_PCV = ( PCV_Y_MIN < y ) && ( y < PCV_Y_MAX );
      bool z_inside_PCV = ( PCV_Z_MIN < z ) && ( z < PCV_Z_MAX );
      return ( x_inside_PCV && y_inside_PCV && z_inside_PCV );
    }

};

// Helper function to set branch addresses for reading information
// from the Event TTree
void set_event_branch_addresses(TTree& etree, AnalysisEvent& ev)
{
  // Reco PDG code of primary PFParticle in slice (i.e., the neutrino
  // candidate)
  etree.SetBranchAddress( "slpdg", &ev.nu_pdg_ );

  // Topological score
  etree.SetBranchAddress( "topological_score", &ev.topological_score_ );
  etree.SetBranchAddress( "CosmicIP", &ev.cosmic_impact_parameter_ );
  //etree.SetBranchAddress( "CosmicIPAll3D", &ev.CosmicIPAll3D_ );

  // Reconstructed neutrino vertex position
  etree.SetBranchAddress( "reco_nu_vtx_x", &ev.nu_vx_ );
  etree.SetBranchAddress( "reco_nu_vtx_y", &ev.nu_vy_ );
  etree.SetBranchAddress( "reco_nu_vtx_z", &ev.nu_vz_ );

  // Reconstructed object counts
  etree.SetBranchAddress( "n_pfps", &ev.num_pf_particles_ );
  etree.SetBranchAddress( "n_tracks", &ev.num_tracks_ );
  etree.SetBranchAddress( "n_showers", &ev.num_showers_ );

  // PFParticle properties
  etree.SetBranchAddress( "pfp_generation_v", &ev.pfp_generation_ );
  etree.SetBranchAddress( "trk_score_v", &ev.pfp_track_score_ );
  etree.SetBranchAddress( "pfpdg", &ev.pfp_reco_pdg_ );
  etree.SetBranchAddress( "backtracked_pdg", &ev.pfp_true_pdg_ );
  etree.SetBranchAddress( "pfnplanehits_Y", &ev.pfp_hitsY_ );

  // Shower properties
  etree.SetBranchAddress( "shr_pfp_id_v", &ev.shower_pfp_id_ );
  etree.SetBranchAddress( "shr_start_x_v", &ev.shower_startx_ );
  etree.SetBranchAddress( "shr_start_y_v", &ev.shower_starty_ );
  etree.SetBranchAddress( "shr_start_z_v", &ev.shower_startz_ );
  // Shower start distance from reco neutrino vertex (pre-calculated for
  // convenience)
  etree.SetBranchAddress( "shr_dist_v", &ev.shower_start_distance_ );

  // Track properties
  etree.SetBranchAddress( "trk_pfp_id_v", &ev.track_pfp_id_ );
  etree.SetBranchAddress( "trk_len_v", &ev.track_length_ );
  etree.SetBranchAddress( "trk_start_x_v", &ev.track_startx_ );
  etree.SetBranchAddress( "trk_start_y_v", &ev.track_starty_ );
  etree.SetBranchAddress( "trk_start_z_v", &ev.track_startz_ );
  // Track start distance from reco neutrino vertex (pre-calculated for
  // convenience)
  etree.SetBranchAddress( "trk_distance_v", &ev.track_start_distance_ );
  etree.SetBranchAddress( "trk_end_x_v", &ev.track_endx_ );
  etree.SetBranchAddress( "trk_end_y_v", &ev.track_endy_ );
  etree.SetBranchAddress( "trk_end_z_v", &ev.track_endz_ );
  etree.SetBranchAddress( "trk_dir_x_v", &ev.track_dirx_ );
  etree.SetBranchAddress( "trk_dir_y_v", &ev.track_diry_ );
  etree.SetBranchAddress( "trk_dir_z_v", &ev.track_dirz_ );
  etree.SetBranchAddress( "trk_energy_proton_v", &ev.track_kinetic_energy_p_ );
  etree.SetBranchAddress( "trk_range_muon_mom_v", &ev.track_range_mom_mu_ );
  etree.SetBranchAddress( "trk_mcs_muon_mom_v", &ev.track_mcs_mom_mu_ );
  etree.SetBranchAddress( "trk_pid_chipr_v", &ev.track_chi2_proton_ );
  etree.SetBranchAddress( "trk_llr_pid_score_v", &ev.track_llr_pid_score_ );

  // MC truth information for the neutrino
  etree.SetBranchAddress( "nu_pdg", &ev.mc_nu_pdg_ );
  etree.SetBranchAddress( "true_nu_vtx_x", &ev.mc_nu_vx_ );
  etree.SetBranchAddress( "true_nu_vtx_y", &ev.mc_nu_vy_ );
  etree.SetBranchAddress( "true_nu_vtx_z", &ev.mc_nu_vz_ );
  etree.SetBranchAddress( "nu_e", &ev.mc_nu_energy_ );
  etree.SetBranchAddress( "ccnc", &ev.mc_nu_ccnc_ );
  etree.SetBranchAddress( "interaction", &ev.mc_nu_interaction_type_ );

  // MC truth information for the final-state primary particles
  etree.SetBranchAddress( "mc_pdg", &ev.mc_nu_daughter_pdg_ );
  etree.SetBranchAddress( "mc_E", &ev.mc_nu_daughter_energy_ );
  etree.SetBranchAddress( "mc_px", &ev.mc_nu_daughter_px_ );
  etree.SetBranchAddress( "mc_py", &ev.mc_nu_daughter_py_ );
  etree.SetBranchAddress( "mc_pz", &ev.mc_nu_daughter_pz_ );

  // GENIE weights
  etree.SetBranchAddress( "weightSpline", &ev.spline_weight_ );
  etree.SetBranchAddress( "weightTune", &ev.tuned_cv_weight_ );
}

// Helper function that creates a branch (or just sets a new address) for a
// simple variable in the output TTree
void set_output_branch_address( TTree& out_tree, const std::string& branch_name,
  void* address, bool create = false, const std::string leaf_spec = "")
{
  if ( create ) {
    if ( leaf_spec != "" ) {
      out_tree.Branch( branch_name.c_str(), address, leaf_spec.c_str() );
    }
    else {
      out_tree.Branch( branch_name.c_str(), address );
    }
  }
  else {
    out_tree.SetBranchAddress( branch_name.c_str(), address );
  }
}

// Helper function template that creates a branch (or just sets a new address)
// for a pointer to an object in the output TTree
template <typename T> void set_object_output_branch_address( TTree& out_tree,
  const std::string& branch_name, T*& address, bool create = false )
{
  if ( create ) out_tree.Branch( branch_name.c_str(), &address );
  else out_tree.SetBranchAddress( branch_name.c_str(), &address );
}

// Helper function to set branch addresses for the output TTree
void set_event_output_branch_addresses(TTree& out_tree, AnalysisEvent& ev,
  bool create = false)
{
  // Signal definition flags
  set_output_branch_address( out_tree, "is_mc", &ev.is_mc_, create, "is_mc/O" );

  set_output_branch_address( out_tree, "mc_neutrino_is_numu",
    &ev.mc_neutrino_is_numu_, create, "mc_neutrino_is_numu/O" );

  set_output_branch_address( out_tree, "mc_vertex_in_FV",
    &ev.mc_vertex_in_FV_, create, "mc_vertex_in_FV/O" );

  set_output_branch_address( out_tree, "mc_pmu_above_threshold",
    &ev.mc_pmu_above_threshold_, create, "mc_pmu_above_threshold/O" );

  set_output_branch_address( out_tree, "mc_has_p_above_threshold",
    &ev.mc_has_p_above_threshold_, create, "mc_has_p_above_threshold/O" );

  set_output_branch_address( out_tree, "mc_no_fs_pi0",
    &ev.mc_no_fs_pi0_, create, "mc_no_fs_pi0/O" );

  set_output_branch_address( out_tree, "mc_no_charged_pi_above_threshold",
    &ev.mc_no_charged_pi_above_threshold_, create,
    "mc_no_charged_pi_above_threshold/O" );

  set_output_branch_address( out_tree, "mc_is_signal",
    &ev.mc_is_signal_, create, "mc_is_signal/O" );

  // MC event category
  set_output_branch_address( out_tree, "category",
    &ev.category_, create, "category/I" );

  // Event weights
  set_output_branch_address( out_tree, "spline_weight",
    &ev.spline_weight_, create, "spline_weight/F" );

  set_output_branch_address( out_tree, "tuned_cv_weight",
    &ev.tuned_cv_weight_, create, "tuned_cv_weight/F" );

  // CCNp0pi selection criteria
  set_output_branch_address( out_tree, "sel_nu_mu_cc", &ev.sel_nu_mu_cc_,
    create, "sel_nu_mu_cc/O" );

  set_output_branch_address( out_tree, "sel_no_reco_showers",
    &ev.sel_no_reco_showers_, create, "sel_no_reco_showers/O" );

  set_output_branch_address( out_tree, "sel_has_muon_candidate",
    &ev.sel_has_muon_candidate_, create,
    "sel_has_muon_candidate/O" );

  set_output_branch_address( out_tree, "sel_muon_above_threshold",
    &ev.sel_muon_above_threshold_, create, "sel_muon_above_threshold/O" );

  set_output_branch_address( out_tree, "sel_has_p_candidate",
    &ev.sel_has_p_candidate_, create, "sel_has_p_candidate/O" );

  set_output_branch_address( out_tree, "sel_passed_proton_pid_cut",
    &ev.sel_passed_proton_pid_cut_, create, "sel_passed_proton_pid_cut/O" );

  set_output_branch_address( out_tree, "sel_protons_contained",
    &ev.sel_protons_contained_, create, "sel_protons_contained/O" );

  set_output_branch_address( out_tree, "sel_lead_p_passed_hits_cut",
    &ev.sel_lead_p_passed_hits_cut_, create, "sel_lead_p_passed_hits_cut/O" );

  set_output_branch_address( out_tree, "sel_lead_p_above_mom_cut",
    &ev.sel_lead_p_above_mom_cut_, create, "sel_lead_p_above_mom_cut/O" );

  set_output_branch_address( out_tree, "sel_CCNp0pi",
    &ev.sel_CCNp0pi_, create, "sel_CCNp0pi/O" );

  // Reco 3-momenta (muon, leading proton)
  set_object_output_branch_address< TVector3 >( out_tree,
    "p3_mu", ev.p3_mu_ptr_, create );

  set_object_output_branch_address< TVector3 >( out_tree,
    "p3_lead_p", ev.p3_lead_p_ptr_, create );

  // True 3-momenta (muon, leading proton)
  set_object_output_branch_address< TVector3 >( out_tree,
    "mc_p3_mu", ev.mc_p3_mu_ptr_, create );

  set_object_output_branch_address< TVector3 >( out_tree,
    "mc_p3_lead_p", ev.mc_p3_lead_p_ptr_, create );

  // Reco STVs
  set_output_branch_address( out_tree, "delta_pT",
    &ev.delta_pT_, create, "delta_pT/F" );

  set_output_branch_address( out_tree, "delta_phiT",
    &ev.delta_phiT_, create, "delta_phiT/F" );

  set_output_branch_address( out_tree, "delta_alphaT",
    &ev.delta_alphaT_, create, "delta_alphaT/F" );

  set_output_branch_address( out_tree, "delta_pL",
    &ev.delta_pL_, create, "delta_pL/F" );

  set_output_branch_address( out_tree, "pn",
    &ev.pn_, create, "pn/F" );

  // MC STVs (only filled for signal events)
  set_output_branch_address( out_tree, "mc_delta_pT",
    &ev.mc_delta_pT_, create, "mc_delta_pT/F" );

  set_output_branch_address( out_tree, "mc_delta_phiT",
    &ev.mc_delta_phiT_, create, "mc_delta_phiT/F" );

  set_output_branch_address( out_tree, "mc_delta_alphaT",
    &ev.mc_delta_alphaT_, create, "mc_delta_alphaT/F" );

  set_output_branch_address( out_tree, "mc_delta_pL",
    &ev.mc_delta_pL_, create, "mc_delta_pL/F" );

  set_output_branch_address( out_tree, "mc_pn",
    &ev.mc_pn_, create, "mc_pn/F" );

  // *** Branches copied directly from the input ***

  // Cosmic rejection parameters for numu CC inclusive selection
  set_output_branch_address( out_tree, "topological_score",
    &ev.topological_score_, create, "topological_score/F" );

  set_output_branch_address( out_tree, "CosmicIP",
    &ev.cosmic_impact_parameter_, create, "CosmicIP/F" );

  // Reconstructed neutrino vertex position
  set_output_branch_address( out_tree, "reco_nu_vtx_x",
    &ev.nu_vx_, create, "reco_nu_vtx_x/F" );

  set_output_branch_address( out_tree, "reco_nu_vtx_y",
    &ev.nu_vy_, create, "reco_nu_vtx_y/F" );

  set_output_branch_address( out_tree, "reco_nu_vtx_z",
    &ev.nu_vz_, create, "reco_nu_vtx_z/F" );

  // MC truth information for the neutrino
  set_output_branch_address( out_tree, "mc_nu_pdg", &ev.mc_nu_pdg_,
    create, "mc_nu_pdg/I" );

  set_output_branch_address( out_tree, "mc_nu_vtx_x", &ev.mc_nu_vx_,
    create, "mc_nu_vtx_x/F" );

  set_output_branch_address( out_tree, "mc_nu_vtx_y", &ev.mc_nu_vy_,
    create, "mc_nu_vtx_y/F" );

  set_output_branch_address( out_tree, "mc_nu_vtx_z", &ev.mc_nu_vz_,
    create, "mc_nu_vtx_z/F" );

  set_output_branch_address( out_tree, "mc_nu_energy", &ev.mc_nu_energy_,
    create, "mc_nu_energy/F" );

  set_output_branch_address( out_tree, "mc_ccnc", &ev.mc_nu_ccnc_,
    create, "mc_ccnc/I" );

  set_output_branch_address( out_tree, "mc_interaction",
    &ev.mc_nu_interaction_type_, create, "mc_interaction/I" );

  // PFParticle properties
  set_object_output_branch_address< std::vector<unsigned int> >( out_tree,
    "pfp_generation_v", ev.pfp_generation_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_score_v", ev.pfp_track_score_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfpdg", ev.pfp_reco_pdg_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "backtracked_pdg", ev.pfp_true_pdg_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfnplanehits_Y", ev.pfp_hitsY_, create );

  // Shower properties
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "shr_start_x_v", ev.shower_startx_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "shr_start_y_v", ev.shower_starty_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "shr_start_z_v", ev.shower_startz_, create );

  // Shower start distance from reco neutrino vertex (pre-calculated for
  // convenience)
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "shr_dist_v", ev.shower_start_distance_, create );

  // Track properties
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_len_v", ev.track_length_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_start_x_v", ev.track_startx_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_start_y_v", ev.track_starty_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_start_z_v", ev.track_startz_, create );

  // Track start distance from reco neutrino vertex (pre-calculated for
  // convenience)
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_distance_v", ev.track_start_distance_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_end_x_v", ev.track_endx_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_end_y_v", ev.track_endy_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_end_z_v", ev.track_endz_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_dir_x_v", ev.track_dirx_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_dir_y_v", ev.track_diry_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_dir_z_v", ev.track_dirz_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_energy_proton_v", ev.track_kinetic_energy_p_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_range_muon_mom_v", ev.track_range_mom_mu_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_mcs_muon_mom_v", ev.track_mcs_mom_mu_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_pid_chipr_v", ev.track_chi2_proton_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_llr_pid_score_v", ev.track_llr_pid_score_, create );

  // MC truth information for the final-state primary particles
  set_object_output_branch_address< std::vector<int> >( out_tree, "mc_pdg",
    ev.mc_nu_daughter_pdg_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree, "mc_E",
    ev.mc_nu_daughter_energy_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree, "mc_px",
    ev.mc_nu_daughter_px_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree, "mc_py",
    ev.mc_nu_daughter_py_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree, "mc_pz",
    ev.mc_nu_daughter_pz_, create );

}

void analyze(const std::vector<std::string>& in_file_names,
  const std::string& output_filename)
{
  // Get the TTrees containing the event ntuples and subrun POT information
  // Use TChain objects for simplicity in manipulating multiple files
  TChain events_ch( "nuselection/NeutrinoSelectionFilter" );
  TChain subruns_ch( "nuselection/SubRun" );

  for ( const auto& f_name : in_file_names ) {
    events_ch.Add( f_name.c_str() );
    subruns_ch.Add( f_name.c_str() );
  }

  // OUTPUT TTREE
  // Make an output TTree for plotting (one entry per event)
  TFile* out_file = new TFile( output_filename.c_str(), "recreate" );
  out_file->cd();
  TTree* out_tree = new TTree( "stv_tree", "STV analysis tree" );

  // Get the total POT from the subruns TTree. Save it in the output
  // TFile as a TParameter<float>
  float pot;
  float summed_pot = 0.;
  subruns_ch.SetBranchAddress( "pot", &pot );
  for ( int se = 0; se < subruns_ch.GetEntries(); ++se ) {
    subruns_ch.GetEntry( se );
    summed_pot += pot;
  }

  TParameter<float>* summed_pot_param = new TParameter<float>( "summed_pot",
    summed_pot );

  summed_pot_param->Write();

  // EVENT LOOP
  // TChains can potentially be really big (and spread out over multiple
  // files). When that's the case, calling TChain::GetEntries() can be very
  // slow. I get around this by using a while loop instead of a for loop.
  bool created_output_branches = false;
  long events_entry = 0;
  while ( true ) {

    if ( events_entry % 1000 == 0 ) {
      std::cout << "Processing event #" << events_entry << '\n';
    }

    // Create a new AnalysisEvent object. This will reset all analysis
    // variables for the current event.
    AnalysisEvent cur_event;

    // Set branch addresses for the member variables that will be read
    // directly from the Event TTree.
    set_event_branch_addresses( events_ch, cur_event );

    // Set the output TTree branch addresses, creating the branches if needed
    // (during the first event loop iteration)
    bool create_them = false;
    if ( !created_output_branches ) {
      create_them = true;
      created_output_branches = true;
    }

    set_event_output_branch_addresses( *out_tree, cur_event, create_them );

    // TChain::LoadTree() returns the entry number that should be used with
    // the current TTree object, which (together with the TBranch objects
    // that it owns) doesn't know about the other TTrees in the TChain.
    // If the return value is negative, there was an I/O error, or we've
    // attempted to read past the end of the TChain.
    int local_entry = events_ch.LoadTree( events_entry );

    // If we've reached the end of the TChain (or encountered an I/O error),
    // then terminate the event loop
    if ( local_entry < 0 ) break;

    // Load all of the branches for which we've called
    // TChain::SetBranchAddress() above
    events_ch.GetEntry( events_entry );

    // Apply the CCNp0pi selection criteria and categorize the event.
    cur_event.apply_selection();

    // Compute observables to save to the output TTree
    cur_event.compute_observables();

    // We're done. Save the results and move on to the next event.
    out_tree->Fill();
    ++events_entry;
  }

  out_tree->Write();
}

// Sets the signal definition flags and returns an event category based on MC
// truth information
EventCategory AnalysisEvent::categorize_event() {
  // TODO: switch to using GHEP truth information?
  // At least verify against it if it is available

  // Real data has a bogus true neutrino PDG code that
  // is below the minimum possible value (-16 <--> nutaubar)
  is_mc_ = ( mc_nu_pdg_ >= TAU_ANTINEUTRINO );
  if ( !is_mc_ ) return kUnknown;

  mc_vertex_in_FV_ = mc_vertex_inside_FV();
  mc_neutrino_is_numu_ = ( mc_nu_pdg_ == MUON_NEUTRINO );

  if ( !mc_vertex_in_FV_ ) {
    mc_is_signal_ = false;
    return kOOFV;
  }
  else if ( mc_nu_ccnc_ == NEUTRAL_CURRENT ) {
    mc_is_signal_ = false;
    return kNC;
  }
  else if ( !mc_neutrino_is_numu_ ) {
    mc_is_signal_ = false;
    if ( mc_nu_pdg_ == ELECTRON_NEUTRINO
      && mc_nu_ccnc_ == CHARGED_CURRENT ) return kNuECC;
    else return kOther;
  }

  // Set flags that default to true here
  mc_no_fs_pi0_ = true;
  mc_no_charged_pi_above_threshold_ = true;

  for ( size_t p = 0u; p < mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = mc_nu_daughter_pdg_->at( p );
    float energy = mc_nu_daughter_energy_->at( p );

    if ( pdg == MUON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(MUON_MASS, 2) );
      if ( mom > MUON_MOM_CUT ) {
        mc_pmu_above_threshold_ = true;
      }
    }
    else if ( pdg == PROTON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PROTON_MASS, 2) );
      if ( mom > LEAD_P_MOM_CUT ) {
        mc_has_p_above_threshold_ = true;
      }
    }
    else if ( pdg == PI_ZERO ) {
      mc_no_fs_pi0_ = false;
    }
    else if ( std::abs(pdg) == PI_PLUS ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PI_PLUS_MASS, 2) );
      if ( mom > CHARGED_PI_MOM_CUT ) {
        mc_no_charged_pi_above_threshold_ = false;
      }
    }
  }

  mc_is_signal_ = mc_vertex_in_FV_ && mc_neutrino_is_numu_
    && mc_pmu_above_threshold_ && mc_has_p_above_threshold_
    && mc_no_fs_pi0_ && mc_no_charged_pi_above_threshold_;

  // Sort signal by interaction mode
  if ( mc_is_signal_ ) {
    if ( mc_nu_interaction_type_ == 0 ) return kSignalCCQE; // QE
    else if ( mc_nu_interaction_type_ == 10 ) return kSignalCCMEC; // MEC
    else if ( mc_nu_interaction_type_ == 1 ) return kSignalCCRES; // RES
    //else if ( mc_nu_interaction_type_ == 2 ) // DIS
    //else if ( mc_nu_interaction_type_ == 3 ) // COH
    else return kSignalOther;
  }
  else if ( mc_no_fs_pi0_ || mc_no_charged_pi_above_threshold_ ) {
    return kNuMuCCNpi;
  }
  else return kOther;
}

void AnalysisEvent::apply_numu_CC_selection() {

  bool reco_vertex_ok = this->reco_vertex_inside_FV();
  bool topo_ok = topological_score_ > TOPO_SCORE_CUT;
  bool cip_ok = cosmic_impact_parameter_ > COSMIC_IP_CUT;

  // Apply the containment cut to the starting positions of all
  // reconstructed tracks and showers. Pass this cut by default.
  bool pfp_starts_ok = true;

  // Loop over each PFParticle in the event
  for ( int p = 0; p < num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Get the track score for the current PFParticle
    float tscore = pfp_track_score_->at( p );

    // A PFParticle is considered a track if its score is above the track score
    // cut. Get the track or shower start coordinates as appropriate.
    float x, y, z;
    if ( tscore > TRACK_SCORE_CUT ) {
      x = track_startx_->at( p );
      y = track_starty_->at( p );
      z = track_startz_->at( p );
    }
    else {
      x = shower_startx_->at( p );
      y = shower_starty_->at( p );
      z = shower_startz_->at( p );
    }

    // Verify that the start of the PFParticle lies within the containment
    // volume.
    // TODO: revisit which containment volume to use for PFParticle start
    // positions. See https://stackoverflow.com/a/2488507 for an explanation
    // of the use of &= here. Don't worry, it's type-safe since both operands
    // are bool.
    pfp_starts_ok &= in_proton_containment_vol( x, y, z );
  }

  sel_nu_mu_cc_ = reco_vertex_ok && topo_ok && cip_ok && pfp_starts_ok;
}

// Sets the index of the muon candidate in the track vectors, or BOGUS_INDEX if
// one could not be found. The sel_has_muon_candidate_ flag is also set by this
// function.
void AnalysisEvent::find_muon_candidate() {

  std::vector<int> muon_candidate_indices;
  std::vector<int> muon_pid_scores;

  for ( int p = 0; p < num_pf_particles_; ++p ) {
    float track_score = pfp_track_score_->at( p );
    float start_dist = track_start_distance_->at( p );
    float track_length = track_length_->at( p );
    float pid_score = track_llr_pid_score_->at( p );

    if ( track_score > MUON_TRACK_SCORE_CUT
      && start_dist < MUON_VTX_DISTANCE_CUT
      && track_length > MUON_LENGTH_CUT
      && pid_score > MUON_PID_CUT )
    {
      muon_candidate_indices.push_back( p );
      muon_pid_scores.push_back( pid_score );
    }
  }

  size_t num_candidates = muon_candidate_indices.size();
  if ( num_candidates > 0u ) sel_has_muon_candidate_ = true;

  if ( num_candidates == 1u ) {
    muon_candidate_idx_ = muon_candidate_indices.front();
  }
  else if ( num_candidates > 1u ) {
    // In the case of multiple muon candidates, choose the one with the highest
    // PID score (most muon-like) as the one to use
    float highest_score = LOW_FLOAT;
    int chosen_index = BOGUS_INDEX;
    for ( size_t c = 0; c < num_candidates; ++c ) {
      float score = muon_pid_scores.at( c );
      if ( highest_score < score ) {
        highest_score = score;
        chosen_index = muon_candidate_indices.at( c );
      }
    }
    muon_candidate_idx_ = chosen_index;
  }
  else {
    muon_candidate_idx_ = BOGUS_INDEX;
  }
}

// Sets the analysis cut flags and decides whether the MC truth information
// matches our signal definition
void AnalysisEvent::apply_selection() {

  // If we're working with an MC event, then categorize the event and set the
  // MC signal flags before proceeding with the selection. This function
  // keeps the category as kUnknown for real data.
  category_ = this->categorize_event();

/////////////
///////DEBUG
//  std::cout << "DEBUG: num_pf_particles = " << num_pf_particles_ << '\n';
//  std::cout << "DEBUG: num_tracks = " << num_tracks_ << '\n';
//  std::cout << "DEBUG: num_showers = " << num_showers_ << '\n';
//  for ( int p = 0; p < pfp_track_score_->size(); ++p ) {
//    std::cout << "DEBUG:    PFP #" << p << ": track score = " << pfp_track_score_->at( p ) << '\n';
//    std::cout << "DEBUG:    PFP #" << p << ": generation = " << pfp_generation_->at( p ) << '\n';
//  }
//  for ( int t = 0; t < track_pfp_id_->size(); ++t ) {
//    std::cout << "DEBUG:    Track #" << t << ": track PFP ID = " << track_pfp_id_->at( t ) << ", track length = " << track_length_->at( t ) << '\n';
//    std::cout << "DEBUG:    track start x = " << track_startx_->at( t ) << '\n';
//    std::cout << "DEBUG:    track start y = " << track_starty_->at( t ) << '\n';
//    std::cout << "DEBUG:    track start z = " << track_startz_->at( t ) << '\n';
//  }
//  for ( int s = 0; s < shower_pfp_id_->size(); ++s ) {
//    std::cout << "DEBUG:    Shower #" << s << ": shower PFP ID = "
//      << shower_pfp_id_->at( s ) << ", start distance = "
//      << shower_start_distance_->at( s ) << '\n';
//    std::cout << "DEBUG:    shower start x = " << shower_startx_->at( s ) << '\n';
//    std::cout << "DEBUG:    shower start y = " << shower_starty_->at( s ) << '\n';
//    std::cout << "DEBUG:    shower start z = " << shower_startz_->at( s ) << '\n';
//
//  }
//////////////

  // Set sel_nu_mu_cc_ by applying those criteria
  this->apply_numu_CC_selection();

  this->find_muon_candidate();
  //std::cout << "DEBUG: muon candidate has index " << muon_candidate_idx_ << '\n';

  // Fail the shower cut if any showers were reconstructed
  // NOTE: We could do this quicker like this,
  //   sel_no_reco_showers_ = ( num_showers_ > 0 );
  // but it might be nice to be able to adjust the track score for this cut.
  // Thus, we do it the hard way.
  int reco_shower_count = 0;
  for ( int p = 0; p < num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    float tscore = pfp_track_score_->at( p );
    if ( tscore <= TRACK_SCORE_CUT ) ++reco_shower_count;
  }
  // Check the shower cut
  sel_no_reco_showers_ = ( reco_shower_count == 0 );

  // Set flags that default to true here
  sel_passed_proton_pid_cut_ = true;
  sel_protons_contained_ = true;

  for ( int p = 0; p < num_pf_particles_; ++p ) {

    // Only worry about direct neutrino daughters (PFParticles considered
    // daughters of the reconstructed neutrino)
    unsigned int generation = pfp_generation_->at( p );
    if ( generation != 2 ) continue;

    // Check that we can find a muon candidate in the event. If more than
    // one is found, also fail the cut.
    if ( p == muon_candidate_idx_ ) {
      // Check that the muon candidate is above threshold
      // TODO: revisit muon momentum (use range if contained)
      float muon_mom = track_mcs_mom_mu_->at( p );
      if ( muon_mom >= MUON_MOM_CUT ) sel_muon_above_threshold_ = true;
    }
    else {

      float track_score = pfp_track_score_->at( p );
      if ( track_score <= TRACK_SCORE_CUT ) continue;

      // Bad tracks in the searchingfornues TTree can have
      // bogus track lengths. This skips those.
      float track_length = track_length_->at( p );
      if ( track_length <= 0. ) continue;

      // We found a reco track that is not the muon candidate. All such
      // tracks are considered proton candidates.
      sel_has_p_candidate_ = true;

      int hitsY = pfp_hitsY_->at( p );
      // TODO: revisit proton PID
      float chi2_proton = track_chi2_proton_->at( p );

      // Check whether the current proton candidate fails the proton PID cut
      if ( hitsY >= HITS_Y_CUT && chi2_proton > PROTON_CHI2_CUT ) {
        sel_passed_proton_pid_cut_ = false;
      }

      // Check whether the current proton candidate fails the containment cut
      float endx = track_endx_->at( p );
      float endy = track_endy_->at( p );
      float endz = track_endz_->at( p );
      bool end_contained = this->in_proton_containment_vol( endx, endy, endz );
      if ( !end_contained ) sel_protons_contained_ = false;
    }

  }

  // Don't bother to apply the cuts that involve the leading
  // proton candidate if we don't have one
  if ( !sel_has_p_candidate_ ) {
    sel_CCNp0pi_ = false;
    return;
  }

  // All that remains is to apply the leading proton candidate cuts. We could
  // search for it above, but doing it here makes the code more readable (with
  // likely negligible impact on performance)
  this->find_lead_p_candidate();

  //std::cout << "DEBUG: lead proton has index " << lead_p_candidate_idx_ << '\n';

  // Check whether the leading proton candidate has the required minimum number
  // of collection plane hits
  float lead_p_hitsY = pfp_hitsY_->at( lead_p_candidate_idx_ );
  if ( lead_p_hitsY >= HITS_Y_CUT ) {
    sel_lead_p_passed_hits_cut_ = true;
  }

  // Check the range-based reco momentum for the leading proton candidate
  float lead_p_KE = track_kinetic_energy_p_->at( lead_p_candidate_idx_ );
  float range_mom_lead_p = real_sqrt( lead_p_KE*lead_p_KE
    + 2.*PROTON_MASS*lead_p_KE );
  if ( range_mom_lead_p >= LEAD_P_MOM_CUT ) {
    sel_lead_p_above_mom_cut_ = true;
  }

  // All right, we've applied all selection cuts. Set the flag that indicates
  // whether all were passed (and thus the event is selected as a CCNp0pi
  // candidate)
  sel_CCNp0pi_ = sel_nu_mu_cc_ && sel_no_reco_showers_
    && sel_has_muon_candidate_ && sel_muon_above_threshold_
    && sel_has_p_candidate_ && sel_passed_proton_pid_cut_
    && sel_protons_contained_ && sel_lead_p_passed_hits_cut_
    && sel_lead_p_above_mom_cut_;
}

void AnalysisEvent::find_lead_p_candidate() {
  float lead_p_track_length = LOW_FLOAT;
  size_t lead_p_index = 0u;
  for ( int p = 0; p < num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Skip the muon candidate reco track (this function assumes that it has
    // already been found)
    if ( p == muon_candidate_idx_ ) continue;

    // Skip PFParticles that are shower-like (track scores near 0)
    float track_score = pfp_track_score_->at( p );
    if ( track_score <= TRACK_SCORE_CUT ) continue;

    // All non-muon-candidate reco tracks are considered proton candidates
    float track_length = track_length_->at( p );
    if ( track_length <= 0. ) continue;

    if ( track_length > lead_p_track_length ) {
      lead_p_track_length = track_length;
      lead_p_index = p;
    }
  }

  // If the leading proton track length changed from its initial
  // value, then we found one. Set the index appropriately.
  if ( lead_p_track_length != LOW_FLOAT ) lead_p_candidate_idx_ = lead_p_index;
  // Otherwise, set the index to BOGUS_INDEX
  else lead_p_candidate_idx_ = BOGUS_INDEX;
}

// Helper function for computing STVs (either reco or true)
void compute_stvs( const TVector3& p3mu, const TVector3& p3p, float& delta_pT,
  float& delta_phiT, float& delta_alphaT, float& delta_pL, float& pn )
{
  delta_pT = (p3mu + p3p).Perp();

  delta_phiT = std::acos( (-p3mu.X()*p3p.X() - p3mu.Y()*p3p.Y())
    / (p3mu.XYvector().Mod() * p3p.XYvector().Mod()) );

  TVector2 delta_pT_vec = (p3mu + p3p).XYvector();
  delta_alphaT = std::acos( (-p3mu.X()*delta_pT_vec.X()
    - p3mu.Y()*delta_pT_vec.Y())
    / (p3mu.XYvector().Mod() * delta_pT_vec.Mod()) );

  float Emu = std::sqrt(std::pow(MUON_MASS, 2) + p3mu.Mag2());
  float Ep = std::sqrt(std::pow(PROTON_MASS, 2) + p3p.Mag2());
  float R = TARGET_MASS + p3mu.Z() + p3p.Z() - Emu - Ep;

  // Estimated mass of the final remnant nucleus (CCQE assumption)
  float mf = TARGET_MASS - NEUTRON_MASS + BINDING_ENERGY;
  delta_pL = 0.5*R - (std::pow(mf, 2) + std::pow(delta_pT, 2)) / (2.*R);

  pn = std::sqrt( std::pow(delta_pL, 2) + std::pow(delta_pT, 2) );
}

void AnalysisEvent::compute_observables() {

  // First compute the MC truth observables (if this is a signal MC event)
  this->compute_mc_truth_observables();

  // Abbreviate some of the calculations below by using these handy
  // references to the muon and leading proton 3-momenta
  auto& p3mu = p3_mu_;
  auto& p3p = p3_lead_p_;

  // Set the reco 3-momentum of the muon candidate if we found one
  bool muon = muon_candidate_idx_ != BOGUS_INDEX;
  if ( muon ) {
    float mu_dirx = track_dirx_->at( muon_candidate_idx_ );
    float mu_diry = track_diry_->at( muon_candidate_idx_ );
    float mu_dirz = track_dirz_->at( muon_candidate_idx_ );
    // TODO: revisit use of MCS momentum here
    float mu_mom = track_mcs_mom_mu_->at( muon_candidate_idx_ );
    p3mu = TVector3( mu_dirx, mu_diry, mu_dirz );
    p3mu = p3mu.Unit() * mu_mom;
  }

  // Set the reco 3-momentum of the leading proton candidate if we found one
  bool lead_p = lead_p_candidate_idx_ != BOGUS_INDEX;
  if ( lead_p ) {

    float p_dirx = track_dirx_->at( lead_p_candidate_idx_ );
    float p_diry = track_diry_->at( lead_p_candidate_idx_ );
    float p_dirz = track_dirz_->at( lead_p_candidate_idx_ );
    float KEp = track_kinetic_energy_p_->at( lead_p_candidate_idx_ );
    float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

    p3p = TVector3( p_dirx, p_diry, p_dirz );
    p3p = p3p.Unit() * p_mom;
  }

  // Compute reco STVs if we have both a muon candidate
  // and a leading proton candidate in the event
  if ( muon && lead_p ) {
    compute_stvs( p3mu, p3p, delta_pT_, delta_phiT_,
      delta_alphaT_, delta_pL_, pn_ );
  }
}

void AnalysisEvent::compute_mc_truth_observables() {

  // If this is not an MC event, then just return without doing anything
  if ( !is_mc_ ) return;

  size_t num_mc_daughters = mc_nu_daughter_pdg_->size();

  // Set the true 3-momentum of the final-state muon if there is one
  bool true_muon = ( mc_neutrino_is_numu_ && mc_nu_ccnc_ == CHARGED_CURRENT );
  if ( true_muon ) {
    // Loop over the MC neutrino daughters, find the muon, and get its
    // true 3-momentum. Note that we assume there is only one muon in
    // this loop.
    bool found_muon = false;
    for ( size_t d = 0u; d < num_mc_daughters; ++d ) {
      int pdg = mc_nu_daughter_pdg_->at( d );
      if ( pdg == MUON ) {
        found_muon = true;
        float px = mc_nu_daughter_px_->at( d );
        float py = mc_nu_daughter_py_->at( d );
        float pz = mc_nu_daughter_pz_->at( d );
        mc_p3_mu_ = TVector3( px, py, pz );
        break;
      }
    }

    if ( !found_muon ) {
      std::cout << "WARNING: Missing muon in MC signal event!\n";
      return;
    }
  }

  // Set the true 3-momentum of the leading proton (if there is one)
  float max_mom = LOW_FLOAT;
  for ( int p = 0; p < num_mc_daughters; ++p ) {
    int pdg = mc_nu_daughter_pdg_->at( p );
    if ( pdg == PROTON )
    {
      float px = mc_nu_daughter_px_->at( p );
      float py = mc_nu_daughter_py_->at( p );
      float pz = mc_nu_daughter_pz_->at( p );
      TVector3 temp_p3 = TVector3( px, py, pz );

      float mom = temp_p3.Mag();
      if ( mom > max_mom ) {
        max_mom = mom;
        mc_p3_lead_p_ = temp_p3;
      }
    }
  }

  // If the event contains a leading proton, then set the 3-momentum
  // accordingly
  bool true_lead_p = ( max_mom != LOW_FLOAT );
  if ( !true_lead_p && mc_is_signal_ ) {
    // If it doesn't for a signal event, then something is wrong.
    std::cout << "WARNING: Missing leading proton in MC signal event!\n";
    return;
  }

  // Compute true STVs if the event contains both a muon and a leading
  // proton
  if ( true_muon && true_lead_p ) {
    compute_stvs( mc_p3_mu_, mc_p3_lead_p_, mc_delta_pT_, mc_delta_phiT_,
      mc_delta_alphaT_, mc_delta_pL_, mc_pn_ );
  }
}

void analyzer(const std::string& in_file_name,
 const std::string& output_filename)
{
  std::vector<std::string> in_files = { in_file_name };
  analyze( in_files, output_filename );
}

int main() {
//analyzer("/uboone/data/users/davidc/searchingfornues/v08_00_00_43/0702/run1/prodgenie_CCmuNoPi_overlay_mcc9_v08_00_00_33_all_run1_reco2_reco2.root", "out.root");
analyzer("/uboone/data/users/gardiner/searchingfornues/prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root", "out.root");
//analyzer("/pnfs/uboone/persistent/users/davidc/searchingfornues/v08_00_00_43/0702/run1/prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root", "out.root");
 return 0;

}
