// Analysis macro for use in the CCNp0pi single transverse variable analysis
// Designed for use with the PeLEE group's "searchingfornues" ntuples
//
// Updated 22 April 2023
// Steven Gardiner <gardiner@fnal.gov>

// Standard library includes
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

// ROOT includes
#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"

// STV analysis includes
#include "EventCategory.hh"
#include "FiducialVolume.hh"
#include "TreeUtils.hh"

// Helper function that avoids NaNs when taking square roots of negative
// numbers
double real_sqrt( double x ) {
  if ( x < 0. ) return 0.;
  else return std::sqrt( x );
}

// Helper function that returns true if a given PDG code represents a meson or
// antimeson. Otherwise returns false. Based on points 10, 12, and 13 of the
// Particle Data Group's "Monte Carlo Particle Numbering Scheme"
// (2019 revision).
bool is_meson_or_antimeson( int pdg_code ) {

  // Ignore differences between mesons and antimesons for this test. Mesons
  // will have positive PDG codes, while antimesons will have negative ones.
  int abs_pdg = std::abs( pdg_code );

  // Meson PDG codes have no more than seven digits. Seven-digit
  // codes beginning with "99" are reserved for generator-specific
  // particles
  if ( abs_pdg >= 9900000 ) return false;

  // Mesons have a value of zero for $n_{q1}$, the thousands digit
  int thousands_digit = ( abs_pdg / 1000 ) % 10;
  if ( thousands_digit != 0 ) return false;

  // They also have a nonzero value for $n_{q2}$, the hundreds digit
  int hundreds_digit = ( abs_pdg / 100 ) % 10;
  if ( hundreds_digit == 0 ) return false;

  // Reserved codes for Standard Model parton distribution functions
  if ( abs_pdg >= 901 && abs_pdg <= 930 ) return false;

  // Reggeon and pomeron
  if ( abs_pdg == 110 || abs_pdg == 990 ) return false;

  // Reserved codes for GEANT tracking purposes
  if ( abs_pdg == 998 || abs_pdg == 999 ) return false;

  // Reserved code for generator-specific pseudoparticles
  if ( abs_pdg == 100 ) return false;

  // If we've passed all of the tests above, then the particle is a meson
  return true;
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
constexpr int TAU_NEUTRINO = 16;
constexpr int PROTON = 2212;
constexpr int PI_ZERO = 111;
constexpr int PI_PLUS = 211;

// Values of parameters to use in analysis cuts
constexpr float DEFAULT_PROTON_PID_CUT = 0.2;
constexpr float LEAD_P_MIN_MOM_CUT = 0.250; // GeV/c
constexpr float LEAD_P_MAX_MOM_CUT = 1.; // GeV/c
constexpr float MUON_P_MIN_MOM_CUT = 0.100; // GeV/c
constexpr float MUON_P_MAX_MOM_CUT = 1.200; // GeV/c
constexpr float CHARGED_PI_MOM_CUT = 0.; // GeV/c
constexpr float MUON_MOM_QUALITY_CUT = 0.25; // fractional difference

constexpr float TOPO_SCORE_CUT = 0.1;
constexpr float COSMIC_IP_CUT = 10.; // cm

constexpr float MUON_TRACK_SCORE_CUT = 0.8;
constexpr float MUON_VTX_DISTANCE_CUT = 4.; // cm
constexpr float MUON_LENGTH_CUT = 10.; // cm
constexpr float MUON_PID_CUT = 0.2;

constexpr float TRACK_SCORE_CUT = 0.5;

// Function that defines the track-length-dependent proton PID cut
double proton_pid_cut( double track_length ) {

  double cut = DEFAULT_PROTON_PID_CUT;

  // Piecewise cut removed 27 June 2021
  //// All track length values are in cm
  //if ( track_length >= 0. && track_length <= 10.5 ) {
  //  cut = -0.0034219305*std::pow( track_length, 2 );
  //  cut += 0.018436866*track_length + 0.062718401;
  //}
  //else if ( track_length > 10.5 && track_length <= 33.1776508 ) {
  //  cut = 0.014153245*( track_length - 10.5 ) - 0.12096235;
  //}

  return cut;

}

// Boundaries of the proton containment volume (used in reco only) in cm
double PCV_X_MIN =   10.;
double PCV_X_MAX =  246.35;

double PCV_Y_MIN = -106.5;
double PCV_Y_MAX =  106.5;

double PCV_Z_MIN =   10.;
double PCV_Z_MAX = 1026.8;

// Mass values from GENIE v3.0.6
constexpr double TARGET_MASS = 37.215526; // 40Ar, GeV
constexpr double NEUTRON_MASS = 0.93956541; // GeV
constexpr double PROTON_MASS = 0.93827208; // GeV
constexpr double MUON_MASS = 0.10565837; // GeV
constexpr double PI_PLUS_MASS = 0.13957000; // GeV

// This binding energy value is used in GENIE v3.0.6
//constexpr double BINDING_ENERGY = 0.0295; // 40Ar, GeV

// This value is the shell-occupancy-weighted mean of the $E_{\alpha}$ values
// listed for 40Ar in Table II of arXiv:1609.03530. MINERvA uses an identical
// procedure for 12C to obtain the binding energy value of 27.13 MeV, which is
// adopted in their STV analysis described in arXiv:1910.08658.
constexpr double BINDING_ENERGY = 0.02478; // 40Ar, GeV

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

    // Backtracked purity and completeness of hits (MC only)
    float nu_completeness_from_pfp_ = BOGUS;
    float nu_purity_from_pfp_ = BOGUS;

    // Reco PDG code of the neutrino candidate
    int nu_pdg_ = BOGUS_INT;

    // Number of neutrino slices identified by the SliceID. Allowed values
    // are zero or one.
    int nslice_ = BOGUS_INT;

    // Reco neutrino vertex coordinates (cm). Space charge corrections have
    // been applied for these.
    float nu_vx_ = BOGUS;
    float nu_vy_ = BOGUS;
    float nu_vz_ = BOGUS;

    // Reconstructed object counts
    int num_pf_particles_ = BOGUS_INT;
    int num_tracks_ = BOGUS_INT;
    int num_showers_ = BOGUS_INT;

    // PFParticle properties
    MyPointer< std::vector<unsigned int> > pfp_generation_;
    MyPointer< std::vector<unsigned int> > pfp_trk_daughters_count_;
    MyPointer< std::vector<unsigned int> > pfp_shr_daughters_count_;

    MyPointer< std::vector<float> > pfp_track_score_;

    // Reco PDG code assigned by Pandora
    MyPointer< std::vector<int> > pfp_reco_pdg_;

    // Total number of wire plane hits associated with each PFParticle
    MyPointer< std::vector<int> > pfp_hits_;

    // Number of hits on the three individual planes
    // (Y is the collection plane)
    MyPointer< std::vector<int> > pfp_hitsU_;
    MyPointer< std::vector<int> > pfp_hitsV_;
    MyPointer< std::vector<int> > pfp_hitsY_;

    // True PDG code found using the backtracker
    MyPointer< std::vector<int> > pfp_true_pdg_;

    // True 4-momentum components found using the backtracker
    MyPointer< std::vector<float> > pfp_true_E_;
    MyPointer< std::vector<float> > pfp_true_px_;
    MyPointer< std::vector<float> > pfp_true_py_;
    MyPointer< std::vector<float> > pfp_true_pz_;

    // Shower properties
    MyPointer< std::vector<unsigned long> > shower_pfp_id_;
    MyPointer< std::vector<float> > shower_startx_;
    MyPointer< std::vector<float> > shower_starty_;
    MyPointer< std::vector<float> > shower_startz_;
    MyPointer< std::vector<float> > shower_start_distance_;

    // Track properties
    MyPointer< std::vector<unsigned long> > track_pfp_id_;
    MyPointer< std::vector<float> > track_length_;
    MyPointer< std::vector<float> > track_startx_;
    MyPointer< std::vector<float> > track_starty_;
    MyPointer< std::vector<float> > track_startz_;
    MyPointer< std::vector<float> > track_start_distance_;
    MyPointer< std::vector<float> > track_endx_;
    MyPointer< std::vector<float> > track_endy_;
    MyPointer< std::vector<float> > track_endz_;
    MyPointer< std::vector<float> > track_dirx_;
    MyPointer< std::vector<float> > track_diry_;
    MyPointer< std::vector<float> > track_dirz_;

    // Proton *kinetic* energy using range-based momentum reconstruction
    MyPointer< std::vector<float> > track_kinetic_energy_p_;

    MyPointer< std::vector<float> > track_range_mom_mu_;
    MyPointer< std::vector<float> > track_mcs_mom_mu_;
    MyPointer< std::vector<float> > track_chi2_proton_;

    // Log-likelihood ratio particle ID information

    // Product of muon/proton log-likelihood ratios from all wire three planes
    MyPointer< std::vector<float> > track_llr_pid_;

    // Individual wire plane muon/proton log-likelihood ratios
    MyPointer< std::vector<float> > track_llr_pid_U_;
    MyPointer< std::vector<float> > track_llr_pid_V_;
    MyPointer< std::vector<float> > track_llr_pid_Y_;

    // Rescaled overall PID score (all three planes) that lies
    // on the interval [-1, 1]
    MyPointer< std::vector<float> > track_llr_pid_score_;

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
    MyPointer< std::vector<int> > mc_nu_daughter_pdg_;
    MyPointer< std::vector<float> > mc_nu_daughter_energy_;
    MyPointer< std::vector<float> > mc_nu_daughter_px_;
    MyPointer< std::vector<float> > mc_nu_daughter_py_;
    MyPointer< std::vector<float> > mc_nu_daughter_pz_;

    // General systematic weights
    MyPointer< std::map< std::string, std::vector<double> > > mc_weights_map_;
    // Map of pointers used to set output branch addresses for the elements
    // of the weights map. Hacky, but it works.
    // TODO: revisit this to make something more elegant
    std::map< std::string, std::vector<double>* > mc_weights_ptr_map_;

    // GENIE weights
    float spline_weight_ = DEFAULT_WEIGHT;
    float tuned_cv_weight_ = DEFAULT_WEIGHT;

    // Signal definition requirements
    bool is_mc_ = false;
    bool mc_neutrino_is_numu_ = false;
    bool mc_vertex_in_FV_ = false;
    bool mc_muon_in_mom_range_ = false;
    bool mc_lead_p_in_mom_range_ = false;
    bool mc_no_fs_mesons_ = false;
    // Intersection of all of these requirements
    bool mc_is_signal_ = false;

    // Extra flags for looking specifically at final-state pions
    bool mc_no_fs_pi0_ = false;
    bool mc_no_charged_pi_above_threshold_ = false;

    EventCategory category_ = kUnknown;

    // **** Reco selection requirements ****

    // Whether the event passed the numu CC selection (a subset of the cuts
    // used for the full analysis)
    bool sel_nu_mu_cc_ = false;

    // Whether the reconstructed neutrino vertex lies within the fiducial
    // volume
    bool sel_reco_vertex_in_FV_ = false;
    // Whether the event passed the topological score cut
    bool sel_topo_cut_passed_ = false;
    // Whether the event passed the cosmic impact parameter cut
    bool sel_cosmic_ip_cut_passed_ = false;
    // Whether the start points for all PFParticles lie within the
    // proton containment volume
    bool sel_pfp_starts_in_PCV_ = false;

    // True if a generation == 2 muon candidate was identified
    bool sel_has_muon_candidate_ = false;

    // Whether the end point of the muon candidate track is contained
    // in the "containment volume"
    bool sel_muon_contained_ = false;

    // Whether the muon candidate has MCS- and range-based reco momenta
    // that agree within a given tolerance
    bool sel_muon_quality_ok_ = false;

    // Whether the muon candidate has a reco momentum above threshold
    bool sel_muon_passed_mom_cuts_ = false;

    // False if at least one generation == 2 shower was reconstructed
    bool sel_no_reco_showers_ = false;
    // Whether at least one generation == 2 reco track exists that is not the
    // muon candidate
    bool sel_has_p_candidate_ = false;
    // Whether all proton candidates (i.e., all tracks which are not the muon
    // candidate) pass the proton PID cut or not
    bool sel_passed_proton_pid_cut_ = false;
    // Whether all proton candidates have track end coordinates that lie within
    // the "containment volume"
    bool sel_protons_contained_ = false;
    // Whether the leading proton candidate has a range-based reco momentum
    // above LEAD_P_MIN_MOM_CUT and below LEAD_P_MAX_MOM_CUT
    bool sel_lead_p_passed_mom_cuts_ = false;
    // Intersection of all of the above requirements
    bool sel_CCNp0pi_ = false;

    // Muon and leading proton candidate indices (BOGUS_INDEX if not present)
    // in the reco track arrays
    int muon_candidate_idx_ = BOGUS_INDEX;
    int lead_p_candidate_idx_ = BOGUS_INDEX;

    // ** Reconstructed observables **

    // 3-momenta
    MyPointer< TVector3 > p3_mu_;
    MyPointer< TVector3 > p3_lead_p_;

    // Reconstructed 3-momenta for all proton candidates,
    // ordered from highest to lowest by magnitude
    MyPointer< std::vector<TVector3> > p3_p_vec_;

    // Reco STVs and other variables of interest
    float delta_pT_ = BOGUS;
    float delta_phiT_ = BOGUS;
    float delta_alphaT_ = BOGUS;
    float delta_pL_ = BOGUS;
    float pn_ = BOGUS;
    float delta_pTx_ = BOGUS;
    float delta_pTy_ = BOGUS;
    float theta_mu_p_ = BOGUS;

    // ** MC truth observables **
    // These are loaded for signal events whenever we have MC information
    // to use

    // 3-momenta
    MyPointer< TVector3 > mc_p3_mu_;
    MyPointer< TVector3 > mc_p3_lead_p_;

    // True 3-momenta for all true MC protons, ordered from highest to lowest
    // by magnitude
    MyPointer< std::vector<TVector3> > mc_p3_p_vec_;

    // MC truth STVs and other variables of interest
    float mc_delta_pT_ = BOGUS;
    float mc_delta_phiT_ = BOGUS;
    float mc_delta_alphaT_ = BOGUS;
    float mc_delta_pL_ = BOGUS;
    float mc_pn_ = BOGUS;
    float mc_delta_pTx_ = BOGUS;
    float mc_delta_pTy_ = BOGUS;
    float mc_theta_mu_p_ = BOGUS;

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

  // Number of neutrino slices identified by the SliceID. Allowed values
  // are zero or one.
  etree.SetBranchAddress( "nslice", &ev.nslice_ );

  // Topological score
  etree.SetBranchAddress( "topological_score", &ev.topological_score_ );
  etree.SetBranchAddress( "CosmicIP", &ev.cosmic_impact_parameter_ );
  //etree.SetBranchAddress( "CosmicIPAll3D", &ev.CosmicIPAll3D_ );

  // Reconstructed neutrino vertex position (with corrections for
  // space charge applied)
  etree.SetBranchAddress( "reco_nu_vtx_sce_x", &ev.nu_vx_ );
  etree.SetBranchAddress( "reco_nu_vtx_sce_y", &ev.nu_vy_ );
  etree.SetBranchAddress( "reco_nu_vtx_sce_z", &ev.nu_vz_ );

  // Reconstructed object counts
  etree.SetBranchAddress( "n_pfps", &ev.num_pf_particles_ );
  etree.SetBranchAddress( "n_tracks", &ev.num_tracks_ );
  etree.SetBranchAddress( "n_showers", &ev.num_showers_ );

  // PFParticle properties
  set_object_input_branch_address( etree, "pfp_generation_v",
    ev.pfp_generation_ );

  set_object_input_branch_address( etree, "pfp_trk_daughters_v",
    ev.pfp_trk_daughters_count_ );

  set_object_input_branch_address( etree, "pfp_shr_daughters_v",
    ev.pfp_shr_daughters_count_ );

  set_object_input_branch_address( etree, "trk_score_v", ev.pfp_track_score_ );
  set_object_input_branch_address( etree, "pfpdg", ev.pfp_reco_pdg_ );
  set_object_input_branch_address( etree, "pfnhits", ev.pfp_hits_ );
  set_object_input_branch_address( etree, "pfnplanehits_U", ev.pfp_hitsU_ );
  set_object_input_branch_address( etree, "pfnplanehits_V", ev.pfp_hitsV_ );
  set_object_input_branch_address( etree, "pfnplanehits_Y", ev.pfp_hitsY_ );

  // Backtracked PFParticle properties
  set_object_input_branch_address( etree, "backtracked_pdg", ev.pfp_true_pdg_ );
  set_object_input_branch_address( etree, "backtracked_e", ev.pfp_true_E_ );
  set_object_input_branch_address( etree, "backtracked_px", ev.pfp_true_px_ );
  set_object_input_branch_address( etree, "backtracked_py", ev.pfp_true_py_ );
  set_object_input_branch_address( etree, "backtracked_pz", ev.pfp_true_pz_ );

  // Shower properties
  // These are excluded from some ntuples to ensure blindness for the LEE
  // analyses. We will skip them when not available.
  bool has_shower_branches = ( etree.GetBranch("shr_pfp_id_v") != nullptr );
  if ( has_shower_branches ) {
    set_object_input_branch_address( etree, "shr_pfp_id_v", ev.shower_pfp_id_ );
    set_object_input_branch_address( etree, "shr_start_x_v", ev.shower_startx_ );
    set_object_input_branch_address( etree, "shr_start_y_v", ev.shower_starty_ );
    set_object_input_branch_address( etree, "shr_start_z_v", ev.shower_startz_ );
    // Shower start distance from reco neutrino vertex (pre-calculated for
    // convenience)
    set_object_input_branch_address( etree, "shr_dist_v",
      ev.shower_start_distance_ );
  }
  else {
    // When the shower information is not available, delete the owned vectors
    // to signal that the associated branches should not be written to the
    // output TTree
    ev.shower_pfp_id_.reset( nullptr );
    ev.shower_startx_.reset( nullptr );
    ev.shower_starty_.reset( nullptr );
    ev.shower_startz_.reset( nullptr );
    ev.shower_start_distance_.reset( nullptr );
  }

  // Track properties
  set_object_input_branch_address( etree, "trk_pfp_id_v", ev.track_pfp_id_ );
  set_object_input_branch_address( etree, "trk_len_v", ev.track_length_ );
  set_object_input_branch_address( etree, "trk_sce_start_x_v", ev.track_startx_ );
  set_object_input_branch_address( etree, "trk_sce_start_y_v", ev.track_starty_ );
  set_object_input_branch_address( etree, "trk_sce_start_z_v", ev.track_startz_ );

  // Track start distance from reco neutrino vertex (pre-calculated for
  // convenience)
  set_object_input_branch_address( etree, "trk_distance_v",
    ev.track_start_distance_ );

  set_object_input_branch_address( etree, "trk_sce_end_x_v", ev.track_endx_ );
  set_object_input_branch_address( etree, "trk_sce_end_y_v", ev.track_endy_ );
  set_object_input_branch_address( etree, "trk_sce_end_z_v", ev.track_endz_ );

  set_object_input_branch_address( etree, "trk_dir_x_v", ev.track_dirx_ );
  set_object_input_branch_address( etree, "trk_dir_y_v", ev.track_diry_ );
  set_object_input_branch_address( etree, "trk_dir_z_v", ev.track_dirz_ );

  set_object_input_branch_address( etree, "trk_energy_proton_v",
    ev.track_kinetic_energy_p_ );

  set_object_input_branch_address( etree, "trk_range_muon_mom_v",
    ev.track_range_mom_mu_ );

  set_object_input_branch_address( etree, "trk_mcs_muon_mom_v",
    ev.track_mcs_mom_mu_ );

  // Some ntuples exclude the old proton chi^2 PID score. Only include it
  // in the output if this branch is available.
  bool has_chipr = ( etree.GetBranch("trk_pid_chipr_v") != nullptr );
  if ( has_chipr ) {
    set_object_input_branch_address( etree, "trk_pid_chipr_v",
      ev.track_chi2_proton_ );
  }
  else {
    ev.track_chi2_proton_.reset( nullptr );
  }

  // Log-likelihood-based particle ID information
  set_object_input_branch_address( etree, "trk_llr_pid_v", ev.track_llr_pid_ );

  set_object_input_branch_address( etree, "trk_llr_pid_u_v",
    ev.track_llr_pid_U_ );

  set_object_input_branch_address( etree, "trk_llr_pid_v_v",
    ev.track_llr_pid_V_ );

  set_object_input_branch_address( etree, "trk_llr_pid_y_v",
    ev.track_llr_pid_Y_ );

  set_object_input_branch_address( etree, "trk_llr_pid_score_v",
    ev.track_llr_pid_score_ );

  // MC truth information for the neutrino
  etree.SetBranchAddress( "nu_pdg", &ev.mc_nu_pdg_ );
  etree.SetBranchAddress( "true_nu_vtx_x", &ev.mc_nu_vx_ );
  etree.SetBranchAddress( "true_nu_vtx_y", &ev.mc_nu_vy_ );
  etree.SetBranchAddress( "true_nu_vtx_z", &ev.mc_nu_vz_ );
  etree.SetBranchAddress( "nu_e", &ev.mc_nu_energy_ );
  etree.SetBranchAddress( "ccnc", &ev.mc_nu_ccnc_ );
  etree.SetBranchAddress( "interaction", &ev.mc_nu_interaction_type_ );

  // MC truth information for the final-state primary particles
  set_object_input_branch_address( etree, "mc_pdg", ev.mc_nu_daughter_pdg_ );
  set_object_input_branch_address( etree, "mc_E", ev.mc_nu_daughter_energy_ );
  set_object_input_branch_address( etree, "mc_px", ev.mc_nu_daughter_px_ );
  set_object_input_branch_address( etree, "mc_py", ev.mc_nu_daughter_py_ );
  set_object_input_branch_address( etree, "mc_pz", ev.mc_nu_daughter_pz_ );

  // GENIE and other systematic variation weights
  bool has_genie_mc_weights = ( etree.GetBranch("weightSpline") != nullptr );
  if ( has_genie_mc_weights ) {
    etree.SetBranchAddress( "weightSpline", &ev.spline_weight_ );
    etree.SetBranchAddress( "weightTune", &ev.tuned_cv_weight_ );
  }

  bool has_weight_map = ( etree.GetBranch("weights") != nullptr );
  if ( has_weight_map ) {
    set_object_input_branch_address( etree, "weights", ev.mc_weights_map_ );
  }
  else {
    ev.mc_weights_map_.reset( nullptr );
  }

  // Purity and completeness of the backtracked hits in the neutrino slice
  bool has_pfp_backtracked_purity = ( etree.GetBranch("nu_purity_from_pfp")
    != nullptr );
  if ( has_pfp_backtracked_purity ) {

    etree.SetBranchAddress( "nu_completeness_from_pfp",
      &ev.nu_completeness_from_pfp_ );

    etree.SetBranchAddress( "nu_purity_from_pfp", &ev.nu_purity_from_pfp_ );

  }

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

  set_output_branch_address( out_tree, "mc_muon_in_mom_range",
    &ev.mc_muon_in_mom_range_, create, "mc_muon_in_mom_range/O" );

  set_output_branch_address( out_tree, "mc_lead_p_in_mom_range",
    &ev.mc_lead_p_in_mom_range_, create, "mc_lead_p_in_mom_range/O" );

  set_output_branch_address( out_tree, "mc_no_fs_pi0",
    &ev.mc_no_fs_pi0_, create, "mc_no_fs_pi0/O" );

  set_output_branch_address( out_tree, "mc_no_charged_pi_above_threshold",
    &ev.mc_no_charged_pi_above_threshold_, create,
    "mc_no_charged_pi_above_threshold/O" );

  set_output_branch_address( out_tree, "mc_no_fs_mesons",
    &ev.mc_no_fs_mesons_, create, "mc_no_fs_mesons/O" );

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

  // If MC weights are available, prepare to store them in the output TTree
  if ( ev.mc_weights_map_ ) {

    // Make separate branches for the various sets of systematic variation
    // weights in the map
    for ( auto& pair : *ev.mc_weights_map_ ) {

      // Prepend "weight_" to the name of the vector of weights in the map
      std::string weight_branch_name = "weight_" + pair.first;

      // Store a pointer to the vector of weights (needed to set the branch
      // address properly) in the temporary map of pointers
      ev.mc_weights_ptr_map_[ weight_branch_name ] = &pair.second;

      // Set the branch address for this vector of weights
      set_object_output_branch_address< std::vector<double> >( out_tree,
        weight_branch_name, ev.mc_weights_ptr_map_.at(weight_branch_name),
        create );
    }
  }

  // Backtracked neutrino purity and completeness
  set_output_branch_address( out_tree, "nu_completeness_from_pfp",
    &ev.nu_completeness_from_pfp_, create, "nu_completeness_from_pfp/F" );

  set_output_branch_address( out_tree, "nu_purity_from_pfp",
    &ev.nu_purity_from_pfp_, create, "nu_purity_from_pfp/F" );

  // Number of neutrino slices identified by the SliceID
  set_output_branch_address( out_tree, "nslice", &ev.nslice_, create,
    "nslice/I" );

  // CCNp0pi selection criteria
  set_output_branch_address( out_tree, "sel_nu_mu_cc", &ev.sel_nu_mu_cc_,
    create, "sel_nu_mu_cc/O" );

  set_output_branch_address( out_tree, "sel_reco_vertex_in_FV",
    &ev.sel_reco_vertex_in_FV_, create, "sel_reco_vertex_in_FV/O" );

  set_output_branch_address( out_tree, "sel_topo_cut_passed",
    &ev.sel_topo_cut_passed_, create, "sel_topo_cut_passed/O" );

  set_output_branch_address( out_tree, "sel_cosmic_ip_cut_passed",
    &ev.sel_cosmic_ip_cut_passed_, create, "sel_cosmic_ip_cut_passed/O" );

  set_output_branch_address( out_tree, "sel_pfp_starts_in_PCV",
    &ev.sel_pfp_starts_in_PCV_, create, "sel_pfp_starts_in_PCV/O" );

  set_output_branch_address( out_tree, "sel_no_reco_showers",
    &ev.sel_no_reco_showers_, create, "sel_no_reco_showers/O" );

  set_output_branch_address( out_tree, "sel_has_muon_candidate",
    &ev.sel_has_muon_candidate_, create,
    "sel_has_muon_candidate/O" );

  set_output_branch_address( out_tree, "sel_muon_contained",
    &ev.sel_muon_contained_, create, "sel_muon_contained/O" );

  set_output_branch_address( out_tree, "sel_muon_passed_mom_cuts",
    &ev.sel_muon_passed_mom_cuts_, create, "sel_muon_passed_mom_cuts/O" );

  set_output_branch_address( out_tree, "sel_muon_quality_ok",
    &ev.sel_muon_quality_ok_, create, "sel_muon_quality_ok/O" );

  set_output_branch_address( out_tree, "sel_has_p_candidate",
    &ev.sel_has_p_candidate_, create, "sel_has_p_candidate/O" );

  set_output_branch_address( out_tree, "sel_passed_proton_pid_cut",
    &ev.sel_passed_proton_pid_cut_, create, "sel_passed_proton_pid_cut/O" );

  set_output_branch_address( out_tree, "sel_protons_contained",
    &ev.sel_protons_contained_, create, "sel_protons_contained/O" );

  set_output_branch_address( out_tree, "sel_lead_p_passed_mom_cuts",
    &ev.sel_lead_p_passed_mom_cuts_, create, "sel_lead_p_passed_mom_cuts/O" );

  set_output_branch_address( out_tree, "sel_CCNp0pi",
    &ev.sel_CCNp0pi_, create, "sel_CCNp0pi/O" );

  // Index for the muon candidate in the vectors of PFParticles
  set_output_branch_address( out_tree, "muon_candidate_idx",
    &ev.muon_candidate_idx_, create, "muon_candidate_idx/I" );

  // Index for the leading proton candidate in the vectors of PFParticles
  set_output_branch_address( out_tree, "lead_p_candidate_idx",
    &ev.lead_p_candidate_idx_, create, "lead_p_candidate_idx/I" );

  // Reco 3-momenta (muon, leading proton)
  set_object_output_branch_address< TVector3 >( out_tree,
    "p3_mu", ev.p3_mu_, create );

  set_object_output_branch_address< TVector3 >( out_tree,
    "p3_lead_p", ev.p3_lead_p_, create );

  // Reco 3-momenta (all proton candidates, ordered from highest to lowest
  // magnitude)
  set_object_output_branch_address< std::vector<TVector3> >( out_tree,
    "p3_p_vec", ev.p3_p_vec_, create );

  // True 3-momenta (muon, leading proton)
  set_object_output_branch_address< TVector3 >( out_tree,
    "mc_p3_mu", ev.mc_p3_mu_, create );

  set_object_output_branch_address< TVector3 >( out_tree,
    "mc_p3_lead_p", ev.mc_p3_lead_p_, create );

  // True 3-momenta (all protons, ordered from highest to lowest magnitude)
  set_object_output_branch_address< std::vector<TVector3> >( out_tree,
    "mc_p3_p_vec", ev.mc_p3_p_vec_, create );

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

  set_output_branch_address( out_tree, "delta_pTx",
    &ev.delta_pTx_, create, "delta_pTx/F" );

  set_output_branch_address( out_tree, "delta_pTy",
    &ev.delta_pTy_, create, "delta_pTy/F" );

  set_output_branch_address( out_tree, "theta_mu_p",
    &ev.theta_mu_p_, create, "theta_mu_p/F" );

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

  set_output_branch_address( out_tree, "mc_delta_pTx",
    &ev.mc_delta_pTx_, create, "mc_delta_pTx/F" );

  set_output_branch_address( out_tree, "mc_delta_pTy",
    &ev.mc_delta_pTy_, create, "mc_delta_pTy/F" );

  set_output_branch_address( out_tree, "mc_theta_mu_p",
    &ev.mc_theta_mu_p_, create, "mc_theta_mu_p/F" );

  // *** Branches copied directly from the input ***

  // Cosmic rejection parameters for numu CC inclusive selection
  set_output_branch_address( out_tree, "topological_score",
    &ev.topological_score_, create, "topological_score/F" );

  set_output_branch_address( out_tree, "CosmicIP",
    &ev.cosmic_impact_parameter_, create, "CosmicIP/F" );

  // Reconstructed neutrino vertex position
  set_output_branch_address( out_tree, "reco_nu_vtx_sce_x",
    &ev.nu_vx_, create, "reco_nu_vtx_sce_x/F" );

  set_output_branch_address( out_tree, "reco_nu_vtx_sce_y",
    &ev.nu_vy_, create, "reco_nu_vtx_sce_y/F" );

  set_output_branch_address( out_tree, "reco_nu_vtx_sce_z",
    &ev.nu_vz_, create, "reco_nu_vtx_sce_z/F" );

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

  set_object_output_branch_address< std::vector<unsigned int> >( out_tree,
    "pfp_trk_daughters_v", ev.pfp_trk_daughters_count_, create );

  set_object_output_branch_address< std::vector<unsigned int> >( out_tree,
    "pfp_shr_daughters_v", ev.pfp_shr_daughters_count_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_score_v", ev.pfp_track_score_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfpdg", ev.pfp_reco_pdg_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfnhits", ev.pfp_hits_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfnplanehits_U", ev.pfp_hitsU_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfnplanehits_V", ev.pfp_hitsV_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfnplanehits_Y", ev.pfp_hitsY_, create );

  // Backtracked PFParticle properties
  set_object_output_branch_address< std::vector<int> >( out_tree,
    "backtracked_pdg", ev.pfp_true_pdg_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "backtracked_e", ev.pfp_true_E_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "backtracked_px", ev.pfp_true_px_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "backtracked_py", ev.pfp_true_py_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "backtracked_pz", ev.pfp_true_pz_, create );

  // Shower properties
  // For some ntuples, reconstructed shower information is excluded.
  // In such cases, skip writing these branches to the output TTree.
  if ( ev.shower_startx_ ) {
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
  }

  // Track properties
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_len_v", ev.track_length_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_start_x_v", ev.track_startx_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_start_y_v", ev.track_starty_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_start_z_v", ev.track_startz_, create );

  // Track start distance from reco neutrino vertex (pre-calculated for
  // convenience)
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_distance_v", ev.track_start_distance_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_end_x_v", ev.track_endx_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_end_y_v", ev.track_endy_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_end_z_v", ev.track_endz_, create );

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

  // Some ntuples exclude the old chi^2 proton PID score. Only include it in
  // the output if it is available.
  if ( ev.track_chi2_proton_ ) {
    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_pid_chipr_v", ev.track_chi2_proton_, create );
  }

  // Log-likelihood-based particle ID information
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_llr_pid_v", ev.track_llr_pid_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_llr_pid_u_v", ev.track_llr_pid_U_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_llr_pid_v_v", ev.track_llr_pid_V_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_llr_pid_y_v", ev.track_llr_pid_Y_, create );

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
  // TFile as a TParameter<float>. Real data doesn't have this TTree,
  // so check that it exists first.
  float pot;
  float summed_pot = 0.;
  bool has_pot_branch = ( subruns_ch.GetBranch("pot") != nullptr );
  if ( has_pot_branch ) {
    subruns_ch.SetBranchAddress( "pot", &pot );
    for ( int se = 0; se < subruns_ch.GetEntries(); ++se ) {
      subruns_ch.GetEntry( se );
      summed_pot += pot;
    }
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

    // Set the output TTree branch addresses, creating the branches if needed
    // (during the first event loop iteration)
    bool create_them = false;
    if ( !created_output_branches ) {
      create_them = true;
      created_output_branches = true;
    }

    set_event_output_branch_addresses( *out_tree, cur_event, create_them );

    // Apply the CCNp0pi selection criteria and categorize the event.
    cur_event.apply_selection();

    // Compute observables to save to the output TTree
    cur_event.compute_observables();

    // We're done. Save the results and move on to the next event.
    out_tree->Fill();
    ++events_entry;
  }

  out_tree->Write();
  out_file->Close();
  delete out_file;
}

// Sets the signal definition flags and returns an event category based on MC
// truth information
EventCategory AnalysisEvent::categorize_event() {

  // Real data has a bogus true neutrino PDG code that is not one of the
  // allowed values (±12, ±14, ±16)
  int abs_mc_nu_pdg = std::abs( mc_nu_pdg_ );
  is_mc_ = ( abs_mc_nu_pdg == ELECTRON_NEUTRINO
    || abs_mc_nu_pdg == MUON_NEUTRINO || abs_mc_nu_pdg == TAU_NEUTRINO );
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

  // Set flags to their default values here
  mc_muon_in_mom_range_ = false;
  mc_lead_p_in_mom_range_ = false;
  mc_no_fs_pi0_ = true;
  mc_no_charged_pi_above_threshold_ = true;
  mc_no_fs_mesons_ = true;

  double lead_p_mom = LOW_FLOAT;

  for ( size_t p = 0u; p < mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = mc_nu_daughter_pdg_->at( p );
    float energy = mc_nu_daughter_energy_->at( p );

    // Do the general check for (anti)mesons first before considering
    // any individual PDG codes
    if ( is_meson_or_antimeson(pdg) ) {
      mc_no_fs_mesons_ = false;
    }


    // Check that the muon has a momentum within the allowed range
    if ( pdg == MUON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(MUON_MASS, 2) );
      if ( mom >= MUON_P_MIN_MOM_CUT && mom <= MUON_P_MAX_MOM_CUT ) {
        mc_muon_in_mom_range_ = true;
      }
    }
    else if ( pdg == PROTON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PROTON_MASS, 2) );
      if ( mom > lead_p_mom ) lead_p_mom = mom;
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

  // Check that the leading proton has a momentum within the allowed range
  if ( lead_p_mom >= LEAD_P_MIN_MOM_CUT && lead_p_mom <= LEAD_P_MAX_MOM_CUT ) {
    mc_lead_p_in_mom_range_ = true;
  }

  mc_is_signal_ = mc_vertex_in_FV_ && mc_neutrino_is_numu_
    && mc_muon_in_mom_range_ && mc_lead_p_in_mom_range_
    && mc_no_fs_mesons_;

  // Sort signal by interaction mode
  if ( mc_is_signal_ ) {
    if ( mc_nu_interaction_type_ == 0 ) return kSignalCCQE; // QE
    else if ( mc_nu_interaction_type_ == 10 ) return kSignalCCMEC; // MEC
    else if ( mc_nu_interaction_type_ == 1 ) return kSignalCCRES; // RES
    //else if ( mc_nu_interaction_type_ == 2 ) // DIS
    //else if ( mc_nu_interaction_type_ == 3 ) // COH
    else return kSignalOther;
  }
  else if ( !mc_no_fs_pi0_ || !mc_no_charged_pi_above_threshold_ ) {
    return kNuMuCCNpi;
  }
  else if ( !mc_lead_p_in_mom_range_ ) {
    return kNuMuCC0pi0p;
  }
  else return kNuMuCCOther;
}

void AnalysisEvent::apply_numu_CC_selection() {

  sel_reco_vertex_in_FV_ = this->reco_vertex_inside_FV();
  sel_topo_cut_passed_ = topological_score_ > TOPO_SCORE_CUT;
  sel_cosmic_ip_cut_passed_ = cosmic_impact_parameter_ > COSMIC_IP_CUT;

  // Apply the containment cut to the starting positions of all
  // reconstructed tracks and showers. Pass this cut by default.
  sel_pfp_starts_in_PCV_ = true;

  // Loop over each PFParticle in the event
  for ( int p = 0; p < num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Use the track reconstruction results to get the start point for
    // every PFParticle for the purpose of verifying containment. We could
    // in principle differentiate between tracks and showers here, but
    // (1) we cut out all showers later on in the selection anyway, and
    // (2) the blinded PeLEE data ntuples do not include shower information.
    // We therefore apply the track reconstruction here unconditionally.
    float x = track_startx_->at( p );
    float y = track_starty_->at( p );
    float z = track_startz_->at( p );

    // Verify that the start of the PFParticle lies within the containment
    // volume.
    // TODO: revisit which containment volume to use for PFParticle start
    // positions. See https://stackoverflow.com/a/2488507 for an explanation
    // of the use of &= here. Don't worry, it's type-safe since both operands
    // are bool.
    sel_pfp_starts_in_PCV_ &= in_proton_containment_vol( x, y, z );
  }

  // Sets the sel_has_muon_candidate_ flag as appropriate. The threshold check
  // is handled later.
  this->find_muon_candidate();

  sel_nu_mu_cc_ = sel_reco_vertex_in_FV_ && sel_pfp_starts_in_PCV_
    && sel_has_muon_candidate_ && sel_topo_cut_passed_;
}

// Sets the index of the muon candidate in the track vectors, or BOGUS_INDEX if
// one could not be found. The sel_has_muon_candidate_ flag is also set by this
// function.
void AnalysisEvent::find_muon_candidate() {

  std::vector<int> muon_candidate_indices;
  std::vector<int> muon_pid_scores;

  for ( int p = 0; p < num_pf_particles_; ++p ) {
    // Only direct neutrino daughters (generation == 2) will be considered as
    // possible muon candidates
    unsigned int generation = pfp_generation_->at( p );
    if ( generation != 2u ) continue;

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

  // Set sel_nu_mu_cc_ by applying those criteria
  this->apply_numu_CC_selection();

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

  // Set flags that default to false here
  sel_muon_contained_ = false;

  for ( int p = 0; p < num_pf_particles_; ++p ) {

    // Only worry about direct neutrino daughters (PFParticles considered
    // daughters of the reconstructed neutrino)
    unsigned int generation = pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Check that we can find a muon candidate in the event. If more than
    // one is found, also fail the cut.
    if ( p == muon_candidate_idx_ ) {

      // Check whether the muon candidate is contained. Use the same
      // containment volume as the protons. TODO: revisit this as needed.
      float endx = track_endx_->at( p );
      float endy = track_endy_->at( p );
      float endz = track_endz_->at( p );
      bool end_contained = this->in_proton_containment_vol( endx, endy, endz );

      if ( end_contained ) sel_muon_contained_ = true;

      // Check that the muon candidate is above threshold. Use the best
      // momentum based on whether it was contained or not.

      float muon_mom = LOW_FLOAT;
      float range_muon_mom = track_range_mom_mu_->at( p );
      float mcs_muon_mom = track_mcs_mom_mu_->at( p );

      if ( sel_muon_contained_ ) muon_mom = range_muon_mom;
      else muon_mom = mcs_muon_mom;

      if ( muon_mom >= MUON_P_MIN_MOM_CUT && muon_mom <= MUON_P_MAX_MOM_CUT ) {
        sel_muon_passed_mom_cuts_ = true;
      }

      // Apply muon candidate quality cut by comparing MCS and range-based
      // momentum estimators. Default to failing the cut.
      sel_muon_quality_ok_ = false;

      double frac_diff_range_mcs = std::abs( range_muon_mom - mcs_muon_mom );
      if ( range_muon_mom > 0. ) {
        frac_diff_range_mcs /= range_muon_mom;
        if ( frac_diff_range_mcs < MUON_MOM_QUALITY_CUT ) {
          sel_muon_quality_ok_ = true;
        }
      }

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

      float llr_pid_score = track_llr_pid_score_->at( p );

      // Check whether the current proton candidate fails the proton PID cut
      if ( llr_pid_score > proton_pid_cut(track_length) ) {
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

  // Check the range-based reco momentum for the leading proton candidate
  float lead_p_KE = track_kinetic_energy_p_->at( lead_p_candidate_idx_ );
  float range_mom_lead_p = real_sqrt( lead_p_KE*lead_p_KE
    + 2.*PROTON_MASS*lead_p_KE );
  if ( range_mom_lead_p >= LEAD_P_MIN_MOM_CUT
    && range_mom_lead_p <= LEAD_P_MAX_MOM_CUT )
  {
    sel_lead_p_passed_mom_cuts_ = true;
  }

  // All right, we've applied all selection cuts. Set the flag that indicates
  // whether all were passed (and thus the event is selected as a CCNp0pi
  // candidate)
  sel_CCNp0pi_ = sel_nu_mu_cc_ && sel_no_reco_showers_
    && sel_muon_passed_mom_cuts_ && sel_muon_contained_ && sel_muon_quality_ok_
    && sel_has_p_candidate_ && sel_passed_proton_pid_cut_
    && sel_protons_contained_ && sel_lead_p_passed_mom_cuts_;
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
  float& delta_phiT, float& delta_alphaT, float& delta_pL, float& pn,
  float& delta_pTx, float& delta_pTy )
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

  // Components of the 2D delta_pT vector (see arXiv:1910.08658)

  // We assume that the neutrino travels along the +z direction (also done
  // in the other expressions above)
  TVector3 zUnit( 0., 0., 1. );

  // Defines the x direction for the components of the delta_pT vector
  TVector2 xTUnit = zUnit.Cross( p3mu ).XYvector().Unit();

  delta_pTx = xTUnit.X()*delta_pT_vec.X() + xTUnit.Y()*delta_pT_vec.Y();

  // Defines the y direction for the components of the delta_T vector
  TVector2 yTUnit = ( -p3mu ).XYvector().Unit();

  delta_pTy = yTUnit.X()*delta_pT_vec.X() + yTUnit.Y()*delta_pT_vec.Y();
}

void AnalysisEvent::compute_observables() {

  // First compute the MC truth observables (if this is a signal MC event)
  this->compute_mc_truth_observables();

  // In cases where we failed to find a muon candidate, check whether there are
  // at least two generation == 2 PFParticles. If there are, then compute the
  // usual observables using the longest track as the muon candidate and the
  // second-longest track as the leading proton candidate. This will enable
  // sideband studies of NC backgrounds in the STV phase space.
  if ( !sel_has_muon_candidate_ ) {

    float max_trk_len = LOW_FLOAT;
    int max_trk_idx = BOGUS_INDEX;

    float next_to_max_trk_len = LOW_FLOAT;
    int next_to_max_trk_idx = BOGUS_INDEX;

    for ( int p = 0; p < num_pf_particles_; ++p ) {

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = pfp_generation_->at( p );
      if ( generation != 2u ) continue;

      float trk_len = track_length_->at( p );

      if ( trk_len > next_to_max_trk_len ) {

        next_to_max_trk_len = trk_len;
        next_to_max_trk_idx = p;

        if ( next_to_max_trk_len > max_trk_len ) {

          next_to_max_trk_len = max_trk_len;
          next_to_max_trk_idx = max_trk_idx;

          max_trk_len = trk_len;
          max_trk_idx = p;
        }
      }
    }

    // If we found at least two usable PFParticles, then assign the indices to
    // be used below
    if ( max_trk_idx != BOGUS_INDEX && next_to_max_trk_idx != BOGUS_INDEX ) {
      muon_candidate_idx_ = max_trk_idx;
      lead_p_candidate_idx_ = next_to_max_trk_idx;
    }
  }

  // Abbreviate some of the calculations below by using these handy
  // references to the muon and leading proton 3-momenta
  auto& p3mu = *p3_mu_;
  auto& p3p = *p3_lead_p_;

  // Set the reco 3-momentum of the muon candidate if we found one
  bool muon = muon_candidate_idx_ != BOGUS_INDEX;
  if ( muon ) {
    float mu_dirx = track_dirx_->at( muon_candidate_idx_ );
    float mu_diry = track_diry_->at( muon_candidate_idx_ );
    float mu_dirz = track_dirz_->at( muon_candidate_idx_ );

    // The selection flag indicating whether the muon candidate is contained
    // was already set when the selection was applied. Use it to choose the
    // best momentum estimator to use.
    float muon_mom = LOW_FLOAT;
    if ( sel_muon_contained_ ) {
      muon_mom = track_range_mom_mu_->at( muon_candidate_idx_ );
    }
    else {
      muon_mom = track_mcs_mom_mu_->at( muon_candidate_idx_ );
    }

    p3mu = TVector3( mu_dirx, mu_diry, mu_dirz );
    p3mu = p3mu.Unit() * muon_mom;
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

  // Reset the vector of reconstructed proton candidate 3-momenta
  p3_p_vec_->clear();

  // Set the reco 3-momenta of all proton candidates (i.e., all generation == 2
  // tracks except the muon candidate) assuming we found both a muon candidate
  // and at least one proton candidate.
  if ( muon && lead_p ) {
    for ( int p = 0; p < num_pf_particles_; ++p ) {
      // Skip the muon candidate
      if ( p == muon_candidate_idx_ ) continue;

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = pfp_generation_->at( p );
      if ( generation != 2u ) continue;

      float p_dirx = track_dirx_->at( p );
      float p_diry = track_diry_->at( p );
      float p_dirz = track_dirz_->at( p );
      float KEp = track_kinetic_energy_p_->at( p );
      float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

      TVector3 p3_temp( p_dirx, p_diry, p_dirz );
      p3_temp = p3_temp.Unit() * p_mom;

      p3_p_vec_->push_back( p3_temp );
    }

    // TODO: reduce code duplication by just getting the leading proton
    // 3-momentum from this sorted vector
    // Sort the reco proton 3-momenta in order from highest to lowest magnitude
    std::sort( p3_p_vec_->begin(), p3_p_vec_->end(), [](const TVector3& a,
      const TVector3& b) -> bool { return a.Mag() > b.Mag(); } );
  }

  // Compute reco STVs if we have both a muon candidate
  // and a leading proton candidate in the event
  if ( muon && lead_p ) {
    compute_stvs( p3mu, p3p, delta_pT_, delta_phiT_,
      delta_alphaT_, delta_pL_, pn_, delta_pTx_, delta_pTy_ );

    theta_mu_p_ = std::acos( p3mu.Dot(p3p) / p3mu.Mag() / p3p.Mag() );
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
        *mc_p3_mu_ = TVector3( px, py, pz );
        break;
      }
    }

    if ( !found_muon ) {
      std::cout << "WARNING: Missing muon in MC signal event!\n";
      return;
    }
  }

  // Reset the vector of true MC proton 3-momenta
  mc_p3_p_vec_->clear();

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

      mc_p3_p_vec_->push_back( temp_p3 );

      float mom = temp_p3.Mag();
      if ( mom > max_mom ) {
        max_mom = mom;
        *mc_p3_lead_p_ = temp_p3;
      }
    }
  }

  // TODO: reduce code duplication by just getting the leading proton
  // 3-momentum from this sorted vector
  // Sort the true proton 3-momenta in order from highest to lowest magnitude
  std::sort( mc_p3_p_vec_->begin(), mc_p3_p_vec_->end(), [](const TVector3& a,
    const TVector3& b) -> bool { return a.Mag() > b.Mag(); } );

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
    compute_stvs( *mc_p3_mu_, *mc_p3_lead_p_, mc_delta_pT_, mc_delta_phiT_,
      mc_delta_alphaT_, mc_delta_pL_, mc_pn_, mc_delta_pTx_, mc_delta_pTy_ );

    mc_theta_mu_p_ = std::acos( mc_p3_mu_->Dot(*mc_p3_lead_p_)
      / mc_p3_mu_->Mag() / mc_p3_lead_p_->Mag() );
  }
}

void analyzer(const std::string& in_file_name,
 const std::string& output_filename)
{
  std::vector<std::string> in_files = { in_file_name };
  analyze( in_files, output_filename );
}

int main( int argc, char* argv[] ) {

  if ( argc != 3 ) {
    std::cout << "Usage: analyzer INPUT_PELEE_NTUPLE_FILE OUTPUT_FILE\n";
    return 1;
  }

  std::string input_file_name( argv[1] );
  std::string output_file_name( argv[2] );

  analyzer( input_file_name, output_file_name );

  return 0;
}
