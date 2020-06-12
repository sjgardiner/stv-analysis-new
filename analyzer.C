// Analysis macro for use in the CCNp0pi single transverse variable analysis
// Run with genie -l
//
// Updated 12 June 2020
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
#include "TTree.h"
#include "TVector3.h"

// GENIE includes (v3 headers)
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/GHEP/GHepParticle.h"

// Helper function that avoids NaNs when taking square roots of negative
// numbers
double real_sqrt( double x ) {
  if ( x < 0. ) return 0.;
  else return std::sqrt( x );
}

// A few helpful dummy constants
constexpr float BOGUS = 9999.;
constexpr int BOGUS_INT = 9999;
constexpr float LOW_FLOAT = -1e30;
constexpr float DEFAULT_WEIGHT = 1.;

// Useful PDG codes
constexpr int ELECTRON_NEUTRINO = 12;
constexpr int MUON = 13;
constexpr int MUON_NEUTRINO = 14;
constexpr int PROTON = 2212;
constexpr int PI_ZERO = 111;
constexpr int PI_PLUS = 211;

// Values of parameters to use in analysis cuts
constexpr float PROTON_CHI2_CUT = 88.;
constexpr unsigned int HITS_Y_CUT = 5u;
constexpr float LEAD_P_MOM_CUT = 0.300; // GeV/c
constexpr float MUON_MOM_CUT = 0.150; // GeV/c
constexpr float CHARGED_PI_MOM_CUT = 0.070; // GeV/c

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

// Class to provide temporary storage for processing each entry in the gtree
// TTree from the NuCCanalyzer module
class AnalysisGenieBranches {
  public:

    AnalysisGenieBranches() {}

    ~AnalysisGenieBranches() {
      //if ( gmcrec_ ) delete gmcrec_;
      if ( weights_map_ ) delete weights_map_;
    }

    //genie::NtpMCEventRecord* gmcrec_ = nullptr;
    std::map< std::string, std::vector<double> >* weights_map_ = nullptr;
};

// Class to hold information from each entry in the gtree TTree
// from the NuCCanalyzer module
class AnalysisGenieRecord {

  public:

    AnalysisGenieRecord(const AnalysisGenieBranches& gb)
      { this->init(gb); }

    double cv_tune_weight_;
    double spline_bugfix_weight_;

    // TODO: add flux + GENIE + Geant4 systematic weights here
    // TODO: add other stuff?

  protected:

    void init(const AnalysisGenieBranches& gb);
};

void AnalysisGenieRecord::init(const AnalysisGenieBranches& gb) {
  // These weights are computed for only a single universe, so just take
  // the single vector entry and store it
  cv_tune_weight_ = gb.weights_map_->at( "TunedCentralValue_UBGenie" ).front();
  spline_bugfix_weight_ = gb.weights_map_->at( "splines_general_Spline" ).front();
}

// Class to hold information from each entry in the Daughters TTree
// from the NuCCanalyzer module
class AnalysisDaughter {

  public:

    AnalysisDaughter() {}

    // Collection plane hits for this daughter
    unsigned int hitsY_ = BOGUS_INT;

    // When generation_ == 2, the PFParticle is considered a direct daughter of
    // the reconstructed neutrino
    unsigned int generation_ = BOGUS_INT;

    bool is_shower_ = false; // Whether this PFParticle was labeled as a shower
    bool is_track_ = false; // Whether this PFParticle was labeled as a track
    bool track_is_muon_candidate_ = false; // Is this the muon candidate?

    float track_length_ = BOGUS; // Length of the reco track
    float track_chi2_proton_ = BOGUS; // PID chi^2 score for proton track hypothesis

    // Log-likelihood PID score for proton track hypothesis using all three planes
    // See docDB #23008 and #23348
    float track_3plane_proton_pid_ = BOGUS;

    // Reco momentum magnitude
    float track_mcs_mom_ = BOGUS; // MCS estimate of track momentum (GeV)
    float track_range_mom_p_ = BOGUS; // Range-based momentum for a proton track (GeV)
    float track_range_mom_mu_ = BOGUS; // Range-based momentum for a muon track (GeV)

    // Track reco end coordinates
    float track_endx_ = BOGUS;
    float track_endy_ = BOGUS;
    float track_endz_ = BOGUS;

    // Track reco start direction
    float track_dirx_ = BOGUS;
    float track_diry_ = BOGUS;
    float track_dirz_ = BOGUS;
};


// Class used to hold information from Wouter's TTree branches
// and process it for our analysis
class AnalysisEvent {

  public:

    AnalysisEvent() {
      mc_nu_daughter_pdg_ = new std::vector<int>;
      mc_nu_daughter_energy_ = new std::vector<float>;
      mc_nu_daughter_px_ = new std::vector<float>;
      mc_nu_daughter_py_ = new std::vector<float>;
      mc_nu_daughter_pz_ = new std::vector<float>;
    }

    ~AnalysisEvent() {
      if ( mc_nu_daughter_pdg_ ) delete mc_nu_daughter_pdg_;
      if ( mc_nu_daughter_energy_ ) delete mc_nu_daughter_energy_;
      if ( mc_nu_daughter_px_ ) delete mc_nu_daughter_px_;
      if ( mc_nu_daughter_py_ ) delete mc_nu_daughter_py_;
      if ( mc_nu_daughter_pz_ ) delete mc_nu_daughter_pz_;
    }

    inline void add_daughter( const AnalysisDaughter& ad )
      { daughters_.push_back( ad ); }

    inline void add_genie_record( const AnalysisGenieRecord& gr )
      { genie_records_.push_back( gr ); }

    EventCategory categorize_event();
    void apply_selection();
    void compute_observables();
    void compute_mc_truth_observables();

    const AnalysisDaughter* get_leading_p_candidate() const;
    const AnalysisDaughter* get_muon_candidate() const;

    // Event timestamps
    unsigned int evt_time_sec_;
    unsigned int evt_time_nsec_;

    // Reco neutrino vertex coordinates (cm)
    float nu_vx_ = BOGUS;
    float nu_vy_ = BOGUS;
    float nu_vz_ = BOGUS;

    // Reco PDG code of the neutrino candidate
    int nu_pdg_ = BOGUS_INT;

    // Whether the reco neutrino vertex lies within the fiducial volume
    bool nu_contained_ = false;

    // True neutrino vertex coordinates (cm)
    float mc_nu_vx_ = BOGUS;
    float mc_nu_vy_ = BOGUS;
    float mc_nu_vz_ = BOGUS;

    // Truth information needed for signal definition

    // True neutrino PDG code
    int mc_nu_pdg_ = BOGUS_INT;

    // True neutrino 4-momentum
    float mc_nu_energy_ = BOGUS;
    float mc_nu_px_ = BOGUS;
    float mc_nu_py_ = BOGUS;
    float mc_nu_pz_ = BOGUS;

    // Whether the event is CC (false) or NC (true)
    bool mc_nu_ccnc_ = false;

    // Interaction mode (QE, MEC, etc.)
    int mc_nu_interaction_type_ = BOGUS_INT;

    // Final-state particle PDG codes and energies (post-FSIs)
    std::vector<int>* mc_nu_daughter_pdg_ = nullptr;
    std::vector<float>* mc_nu_daughter_energy_ = nullptr;
    std::vector<float>* mc_nu_daughter_px_ = nullptr;
    std::vector<float>* mc_nu_daughter_py_ = nullptr;
    std::vector<float>* mc_nu_daughter_pz_ = nullptr;

    // This vector stores information retrieved from the Daughters TTree for
    // this event
    std::vector< AnalysisDaughter > daughters_;

    // This vector stores information retrieved from the gtree TTree for
    // this event
    std::vector< AnalysisGenieRecord > genie_records_;

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

    float spline_weight_ = DEFAULT_WEIGHT;
    float tuned_cv_weight_ = DEFAULT_WEIGHT;

    // **** Reco selection requirements ****

    // Whether the event passed the upstream numu CC selection (by Wouter)
    bool sel_nu_mu_cc_ = false;
    // False if at least one generation == 2 shower was reconstructed
    bool sel_no_reco_showers_ = false;
    // True if exactly one generation == 2 muon candidate was identified
    bool sel_has_single_muon_candidate_ = false;
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

    bool mc_vertex_inside_FV() {
      bool x_inside_FV = ( FV_X_MIN < mc_nu_vx_ ) && ( mc_nu_vx_ < FV_X_MAX );
      bool y_inside_FV = ( FV_Y_MIN < mc_nu_vy_ ) && ( mc_nu_vy_ < FV_Y_MAX );
      bool z_inside_FV = ( FV_Z_MIN < mc_nu_vz_ ) && ( mc_nu_vz_ < FV_Z_MAX );
      return ( x_inside_FV && y_inside_FV && z_inside_FV );
    }

    bool in_proton_containment_vol(const AnalysisDaughter& ad) {
      bool x_inside_PCV = ( PCV_X_MIN < ad.track_endx_ )
        && ( ad.track_endx_ < PCV_X_MAX );
      bool y_inside_PCV = ( PCV_Y_MIN < ad.track_endy_ )
        && ( ad.track_endy_ < PCV_Y_MAX );
      bool z_inside_PCV = ( PCV_Z_MIN < ad.track_endz_ )
        && ( ad.track_endz_ < PCV_Z_MAX );
      return ( x_inside_PCV && y_inside_PCV && z_inside_PCV );
    }

};

// Helper function to set branch addresses for reading information
// from the gtree TTree
void set_gtree_branch_addresses(TTree& gtree, AnalysisGenieBranches& gb)
{
  //gtree.SetBranchAddress("gmcrec", &gb.gmcrec_ );
  gtree.SetBranchStatus( "gmcrec", false );
  gtree.SetBranchAddress("weights_map", &gb.weights_map_ );
}

// Helper function to set branch addresses for reading information
// from the Daughters TTree
void set_daughter_branch_addresses(TTree& dtree, AnalysisDaughter& ad)
{
  dtree.SetBranchAddress("hitsY", &ad.hitsY_ );
  dtree.SetBranchAddress("generation", &ad.generation_ );
  dtree.SetBranchAddress("is_track", &ad.is_track_ );
  dtree.SetBranchAddress("is_shower", &ad.is_shower_ );
  dtree.SetBranchAddress("track_length", &ad.track_length_ );
  dtree.SetBranchAddress("track_endx", &ad.track_endx_ );
  dtree.SetBranchAddress("track_endy", &ad.track_endy_ );
  dtree.SetBranchAddress("track_endz", &ad.track_endz_ );
  dtree.SetBranchAddress("track_dirx", &ad.track_dirx_ );
  dtree.SetBranchAddress("track_diry", &ad.track_diry_ );
  dtree.SetBranchAddress("track_dirz", &ad.track_dirz_ );
  dtree.SetBranchAddress("track_is_muon_candidate",
    &ad.track_is_muon_candidate_ );
  dtree.SetBranchAddress("track_range_mom_p", &ad.track_range_mom_p_ );
  dtree.SetBranchAddress("track_range_mom_mu", &ad.track_range_mom_mu_ );
  dtree.SetBranchAddress("track_mcs_mom", &ad.track_mcs_mom_ );
  dtree.SetBranchAddress("track_chi2_proton", &ad.track_chi2_proton_ );
  dtree.SetBranchAddress("track_3plane_proton_pid", &ad.track_3plane_proton_pid_ );
}

// Helper function to set branch addresses for reading information
// from the Event TTree
void set_event_branch_addresses(TTree& etree, AnalysisEvent& ev)
{
  etree.SetBranchAddress("nu_mu_cc_selected", &ev.sel_nu_mu_cc_ );
  etree.SetBranchAddress("evt_time_sec", &ev.evt_time_sec_ );
  etree.SetBranchAddress("evt_time_nsec", &ev.evt_time_nsec_ );
  etree.SetBranchAddress("nu_vx", &ev.nu_vx_ );
  etree.SetBranchAddress("nu_vy", &ev.nu_vy_ );
  etree.SetBranchAddress("nu_vz", &ev.nu_vz_ );
  etree.SetBranchAddress("nu_contained", &ev.nu_contained_ );
  etree.SetBranchAddress("nu_pdg", &ev.nu_pdg_ );

  // Decide whether we are working with data or MC events by checking
  // the Event TTree. If it has a branch named "mc_nu_pdg", then we
  // are using MC events. Otherwise, we're using data. We default
  // to is_mc_ == false in the AnalysisEvent class, so make the
  // switch to true if needed.
  TBranch* temp_mc_branch = etree.GetBranch( "mc_nu_pdg" );
  if ( temp_mc_branch ) ev.is_mc_ = true;
  // If we're working with data, we don't need to set the remaining
  // branch addresses
  else return;

  etree.SetBranchAddress("mc_nu_pdg", &ev.mc_nu_pdg_ );
  etree.SetBranchAddress("mc_nu_vx", &ev.mc_nu_vx_ );
  etree.SetBranchAddress("mc_nu_vy", &ev.mc_nu_vy_ );
  etree.SetBranchAddress("mc_nu_vz", &ev.mc_nu_vz_ );
  etree.SetBranchAddress("mc_nu_energy", &ev.mc_nu_energy_ );
  etree.SetBranchAddress("mc_nu_px", &ev.mc_nu_px_ );
  etree.SetBranchAddress("mc_nu_py", &ev.mc_nu_py_ );
  etree.SetBranchAddress("mc_nu_pz", &ev.mc_nu_pz_ );
  etree.SetBranchAddress("mc_nu_ccnc", &ev.mc_nu_ccnc_ );
  etree.SetBranchAddress("mc_nu_interaction_type",
    &ev.mc_nu_interaction_type_ );
  etree.SetBranchAddress("mc_nu_daughter_pdg", &ev.mc_nu_daughter_pdg_ );
  etree.SetBranchAddress("mc_nu_daughter_energy", &ev.mc_nu_daughter_energy_ );
  etree.SetBranchAddress("mc_nu_daughter_px", &ev.mc_nu_daughter_px_ );
  etree.SetBranchAddress("mc_nu_daughter_py", &ev.mc_nu_daughter_py_ );
  etree.SetBranchAddress("mc_nu_daughter_pz", &ev.mc_nu_daughter_pz_ );

  etree.SetBranchAddress("event_weight", &ev.spline_weight_ );
}

// Helper function that creates a branch (or just sets a new address)
// for a variable in the output TTree
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
  set_output_branch_address( out_tree, "sel_nu_mu_cc", &ev.sel_nu_mu_cc_, create,
    "sel_nu_mu_cc/O" );

  set_output_branch_address( out_tree, "sel_no_reco_showers",
    &ev.sel_no_reco_showers_, create, "sel_no_reco_showers/O" );
  set_output_branch_address( out_tree, "sel_has_single_muon_candidate",
    &ev.sel_has_single_muon_candidate_, create,
    "sel_has_single_muon_candidate/O" );
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
  // Note: these branches contain TVector3 objects
  if ( create ) {
    // Reco
    out_tree.Branch( "p3_mu", &ev.p3_mu_ptr_ );
    out_tree.Branch( "p3_lead_p", &ev.p3_lead_p_ptr_ );

    // True
    out_tree.Branch( "mc_p3_mu", &ev.mc_p3_mu_ptr_ );
    out_tree.Branch( "mc_p3_lead_p", &ev.mc_p3_lead_p_ptr_ );
  }
  else {
    // Reco
    out_tree.SetBranchAddress( "p3_mu", &ev.p3_mu_ptr_ );
    out_tree.SetBranchAddress( "p3_lead_p", &ev.p3_lead_p_ptr_ );

    // True
    out_tree.SetBranchAddress( "mc_p3_mu", &ev.mc_p3_mu_ptr_ );
    out_tree.SetBranchAddress( "mc_p3_lead_p", &ev.mc_p3_lead_p_ptr_ );
  }

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
}

void analyze(const std::vector<std::string>& in_file_names,
  const std::string& output_filename)
{
  // Get the Events, Daughters, subruns, and gtree TTrees
  // Use TChain objects for simplicity in manipulating
  // multiple files
  TChain events_ch("NuCCanalyzer/Event");
  TChain daughters_ch("NuCCanalyzer/Daughters");
  TChain subruns_ch("NuCCanalyzer/subruns");
  TChain gtree_ch("NuCCanalyzer/gtree");

  for ( const auto& f_name : in_file_names ) {
    events_ch.Add( f_name.c_str() );
    daughters_ch.Add( f_name.c_str() );
    subruns_ch.Add( f_name.c_str() );
    gtree_ch.Add( f_name.c_str() );
  }

  // Configure storage for the branch variables that we care
  // about. Start with the Event TTree

  // *** Event TTree ***

  // Reference numbers for the Event TTree
  unsigned int events_event, events_run, events_subrun;

  // Reference indices are not saved for use in the analysis
  events_ch.SetBranchAddress("event", &events_event);
  events_ch.SetBranchAddress("run", &events_run);
  events_ch.SetBranchAddress("subrun", &events_subrun);

  // *** Daughters TTree ***

  // Reference numbers for the Daughters TTree
  unsigned int daughters_event, daughters_run, daughters_subrun;

  // Reference indices are not saved for use in the analysis
  daughters_ch.SetBranchAddress("event", &daughters_event );
  daughters_ch.SetBranchAddress("run", &daughters_run );
  daughters_ch.SetBranchAddress("subrun", &daughters_subrun );

  // *** gtree TTree ***

  // Reference numbers for the gtree TTree
  unsigned int gtree_event, gtree_run, gtree_subrun;

  // Reference indices are not saved for use in the analysis
  gtree_ch.SetBranchAddress("event", &gtree_event );
  gtree_ch.SetBranchAddress("run", &gtree_run );
  gtree_ch.SetBranchAddress("subrun", &gtree_subrun );

  // Use separate entry counters for the different TTrees to keep
  // event, run, and subrun in sync
  int events_entry = 0;
  int daughters_entry = 0;
  int gtree_entry = 0;

  // OUTPUT TTREE
  // Make an output TTree for plotting (one entry per event)
  TFile* out_file = new TFile(output_filename.c_str(), "recreate");
  out_file->cd();
  TTree* out_tree = new TTree("stv_tree", "STV analysis tree");

  // Get the total POT from the subruns TTree. Save it in the output
  // TFile as a TParameter<float>
  float pot;
  float summed_pot = 0.;
  subruns_ch.SetBranchAddress("pot", &pot);
  for ( int se = 0; se < subruns_ch.GetEntries(); ++se ) {
    subruns_ch.GetEntry( se );
    summed_pot += pot;
  }

  TParameter<float>* summed_pot_param = new TParameter<float>("summed_pot",
    summed_pot);

  summed_pot_param->Write();

  // I make the assumption here that the Event and Daughters TTrees
  // start out with their entries aligned. That appears to be the case
  // in Wouter's outputs, but it could cause problems if there are any
  // issues there.

  // EVENT LOOP
  // TChains can potentially be really big (and spread out over multiple
  // files). When that's the case, calling TChain::GetEntries() can be very
  // slow. I get around this by using a while loop instead of a for loop.
  bool created_output_branches = false;
  long event_counter = 0;
  while ( true ) {

    if ( event_counter % 1000 == 0 ) {
      std::cout << "Processing event #" << event_counter << '\n';
    }

    // Create a new AnalysisEvent object. This will reset all analysis
    // variables for the current event.
    AnalysisEvent cur_event;

    // Set branch addresses for the member variables that will be read
    // directly from the Event TTree.
    // NOTE: This function will also set the is_mc_ member flag
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

    // DAUGHTER LOOP
    // TODO: Daughter count doesn't always match TTree entries.
    // Follow up with Wouter about this. We work around this
    // by using a while loop below.

    // Create temporary storage for the Daughters TTree branch variables and
    // set up the branch addresses
    AnalysisDaughter temp_daughter;
    set_daughter_branch_addresses( daughters_ch, temp_daughter );

    while ( daughters_ch.GetEntry( daughters_entry ),
      events_event == daughters_event && events_run == daughters_run
        && events_subrun == daughters_subrun )
    {
      // Store the information about the current daughter for later
      // analysis
      cur_event.add_daughter( temp_daughter );

      // Advance to the next daughter
      ++daughters_entry;

      // daughter_local_entry will be negative if we've
      // reached the end of the Daughter TChain
      int daughter_local_entry = daughters_ch.LoadTree( daughters_entry );

      if ( daughter_local_entry < 0 ) break;
    }

    // Only MC output will have GENIE event records and weights
    if ( cur_event.is_mc_ ) {

      // Create temporary storage for the gtree TTree branch variables and
      // set up the branch addresses
      AnalysisGenieBranches temp_genie_branches;
      set_gtree_branch_addresses( gtree_ch, temp_genie_branches );

      while ( gtree_ch.GetEntry( gtree_entry ),
        events_event == gtree_event && events_run == gtree_run
          && events_subrun == gtree_subrun )
      {
        // Process the information from the current gtree entry and
        // save the results we care about in an AnalysisGenieRecord object
        AnalysisGenieRecord temp_genie_record( temp_genie_branches );

        // Store the extended GENIE event / EventWeight information for later
        // analysis
        cur_event.add_genie_record( temp_genie_record );

        // Avoid memory leaks by manually deleting the owned GENIE EventRecord object
        // (necessary for dumb ROOT reasons)
        //delete temp_genie_branches.gmcrec_->event;

        // Advance to the next GENIE event
        ++gtree_entry;

        // gtree_local_entry will be negative if we've
        // reached the end of the gtree TChain
        int gtree_local_entry = gtree_ch.LoadTree( gtree_entry );

        if ( gtree_local_entry < 0 ) break;
      }

    }

    // We've finished looping over all reco neutrino daughters. Now
    // apply the CCNp0pi selection criteria and categorize the event.
    cur_event.apply_selection();

    // Compute observables to save to the output TTree
    cur_event.compute_observables();

    // We're done. Save the results and move on to the next event.
    out_tree->Fill();
    ++events_entry;
    ++event_counter;
  }

  out_tree->Write();
}

// Sets the signal definition flags and returns an event category based on MC
// truth information
EventCategory AnalysisEvent::categorize_event() {
  // TODO: switch to using GHEP truth information?
  // At least verify against it if it is available

  // This flag will be set before this function is called by
  // set_event_branch_addresses()
  if ( !is_mc_ ) return kUnknown;

  mc_vertex_in_FV_ = mc_vertex_inside_FV();
  mc_neutrino_is_numu_ = ( mc_nu_pdg_ == MUON_NEUTRINO );

  if ( !mc_vertex_in_FV_ ) {
    mc_is_signal_ = false;
    return kOOFV;
  }
  else if ( mc_nu_ccnc_ ) {
    // True is used to represent NC in Wouter's Event TTree
    mc_is_signal_ = false;
    return kNC;
  }
  else if ( !mc_neutrino_is_numu_ ) {
    mc_is_signal_ = false;
    if ( mc_nu_pdg_ == ELECTRON_NEUTRINO && !mc_nu_ccnc_ ) return kNuECC;
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

// Sets the analysis cut flags and decides whether the MC truth information
// matches our signal definition
void AnalysisEvent::apply_selection() {

  // If we're working with an MC event, then categorize the event and set the
  // MC signal flags before proceeding with the selection
  if ( is_mc_ ) {
    category_ = this->categorize_event();
  }

  // The numu CC selection already applied in Wouter's Event TTree,
  // and sel_nu_mu_cc_ will thus already be set appropriately

  // Set flags that default to true here
  sel_no_reco_showers_ = true;
  sel_passed_proton_pid_cut_ = true;
  sel_protons_contained_ = true;

  // Loop over the direct neutrino daughters (generation == 2) and apply the
  // selection cuts
  bool already_found_muon_candidate = false;

  for ( const auto& d : daughters_ ) {

    // Only worry about direct neutrino daughters (PFParticles considered
    // daughters of the reconstructed neutrino)
    if ( d.generation_ != 2 ) continue;

    // If any of the direct neutrino daughters is a shower, then fail the "no
    // showers" cut
    if ( d.is_shower_ ) sel_no_reco_showers_ = false;

    // Skip daughters that are neither a reco track nor a reco shower
    if ( !d.is_track_ ) continue;

    // Check that we can find a muon candidate in the event. If more than
    // one is found, also fail the cut.
    if ( d.track_is_muon_candidate_ ) {
      if ( !already_found_muon_candidate ) {
        sel_has_single_muon_candidate_ = true;
        already_found_muon_candidate = true;

        // Check that the muon candidate is above threshold
        if ( d.track_mcs_mom_ >= MUON_MOM_CUT ) {
          sel_muon_above_threshold_ = true;
        }
      }
      else {
        // More than one muon candidate track found
        sel_has_single_muon_candidate_ = false;
      }
    }
    else {
      // We found a reco track that is not the muon candidate. All such
      // tracks are considered proton candidates.
      sel_has_p_candidate_ = true;

      // Check whether the current proton candidate fails the proton PID cut
      if ( d.hitsY_ >= HITS_Y_CUT && d.track_chi2_proton_ > PROTON_CHI2_CUT ) {
        sel_passed_proton_pid_cut_ = false;
      }

      // Check whether the current proton candidate fails the containment cut
      if ( !this->in_proton_containment_vol(d) ) sel_protons_contained_ = false;
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
  const auto* lead_p_daughter = this->get_leading_p_candidate();

  // TODO: add check for lead_p_daughter == nullptr (shouldn't event happen
  // becaus we already verified above that we have at least one proton
  // candidate in the event)

  // Check whether the leading proton candidate has the required minimum number
  // of collection plane hits
  if ( lead_p_daughter->hitsY_ >= HITS_Y_CUT ) {
    sel_lead_p_passed_hits_cut_ = true;
  }

  // Check the range-based reco momentum for the leading proton candidate
  if ( lead_p_daughter->track_range_mom_p_ >= LEAD_P_MOM_CUT ) {
    sel_lead_p_above_mom_cut_ = true;
  }

  // All right, we've applied all selection cuts. Set the flag that indicates
  // whether all were passed (and thus the event is selected as a CCNp0pi
  // candidate)
  sel_CCNp0pi_ = sel_nu_mu_cc_ && sel_no_reco_showers_
    && sel_has_single_muon_candidate_ && sel_muon_above_threshold_
    && sel_has_p_candidate_ && sel_passed_proton_pid_cut_
    && sel_protons_contained_ && sel_lead_p_passed_hits_cut_
    && sel_lead_p_above_mom_cut_;
}

// Returns a pointer to the leading proton candidate neutrino daughter, or a
// nullptr if one could not be found
const AnalysisDaughter* AnalysisEvent::get_leading_p_candidate() const {
  float lead_p_track_length = LOW_FLOAT;
  size_t lead_p_index = 0u;
  for ( size_t j = 0u; j < daughters_.size(); ++j ) {

    const auto& d = daughters_.at( j );

    // Skip daughters that are not reco tracks
    if ( !d.is_track_ ) continue;

    // Skip the muon candidate reco track
    if ( d.track_is_muon_candidate_ ) continue;

    // All non-muon-candidate reco tracks are considered proton candidates
    if ( d.track_length_ > lead_p_track_length ) {
      lead_p_track_length = d.track_length_;
      lead_p_index = j;
    }
  }

  // If the leading proton track length changed from its initial
  // value, then we found one. Return a pointer to it.
  if ( lead_p_track_length != LOW_FLOAT ) return &daughters_.at( lead_p_index );
  // Otherwise, return a nullptr
  return nullptr;
}

// Returns a pointer to the (first) muon candidate neutrino daughter, or a
// nullptr if one could not be found
const AnalysisDaughter* AnalysisEvent::get_muon_candidate() const {
  for ( const auto& d : daughters_ ) {
    if ( d.is_track_ && d.track_is_muon_candidate_ ) return &d;
  }
  return nullptr;
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
  const auto* muon = this->get_muon_candidate();
  if ( muon ) {
    p3mu = TVector3( muon->track_dirx_, muon->track_diry_, muon->track_dirz_ );
    p3mu = p3mu.Unit() * muon->track_mcs_mom_;
  }

  // Set the reco 3-momentum of the leading proton candidate if we found one
  const auto* lead_p = this->get_leading_p_candidate();
  if ( lead_p ) {

    p3p = TVector3( lead_p->track_dirx_, lead_p->track_diry_,
      lead_p->track_dirz_ );

    p3p = p3p.Unit() * lead_p->track_range_mom_p_;
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
  bool true_muon = ( mc_neutrino_is_numu_ && !mc_nu_ccnc_ );
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
