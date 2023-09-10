#include "ConfigMakerUtils.hh"
#include "HistUtils.hh"
#include "UniverseMaker.hh"
#include "SliceBinning.hh"

// Placeholder value for the block index for bins in which it is irrelevant
constexpr int DUMMY_BLOCK_INDEX = -1;

// Sideband selection cuts
const std::string sel_dirt =
  " !sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
  " && sel_topo_cut_passed && sel_has_muon_candidate"
  " && sel_muon_contained && sel_muon_quality_ok"
  " && sel_muon_passed_mom_cuts && sel_no_reco_showers"
  " && sel_has_p_candidate && sel_protons_contained"
  " && sel_passed_proton_pid_cut && sel_lead_p_passed_mom_cuts";

const std::string sel_NC =
  " sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV && sel_topo_cut_passed"
  " && !sel_has_muon_candidate && muon_candidate_idx != -1"
  " && sel_no_reco_showers && sel_has_p_candidate"
  " && sel_protons_contained && sel_passed_proton_pid_cut"
  " && sel_lead_p_passed_mom_cuts";

const std::string sel_CCNpi =
  " sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
  " && sel_topo_cut_passed && sel_has_muon_candidate"
  " && sel_muon_contained && sel_muon_quality_ok"
  " && sel_muon_passed_mom_cuts && sel_no_reco_showers"
  " && sel_has_p_candidate && sel_protons_contained"
  " && Sum$( pfp_generation_v == 2 && trk_llr_pid_score_v > 0.2 ) > 1"
  " && sel_lead_p_passed_mom_cuts";

const std::string sel_combined =
  "((" + sel_dirt + ") || (" + sel_NC + ") || (" + sel_CCNpi + "))";

struct EdgeDef {

  EdgeDef( std::map< double, std::vector<double> >* edges, bool use_overflow,
    const std::string& var_to_use, const std::string& act_var_name,
    const std::string& oth_var_name) : bin_edges_2d_( edges ),
    needs_overflow_( use_overflow ), branch_var_name_( var_to_use ),
    active_var_name_( act_var_name ), other_var_name_( oth_var_name ) {}

  std::map< double, std::vector<double> >* bin_edges_2d_;
  bool needs_overflow_;
  std::string branch_var_name_;
  std::string active_var_name_;
  std::string other_var_name_;
};

std::vector< EdgeDef > MOM_ANG_2D_EDGES = {
  { &MUON_2D_BIN_EDGES, false, "p3_mu", "cos#theta_{#mu}",
    "p_{#mu}" },
  { &PROTON_2D_BIN_EDGES, false, "p3_lead_p", "cos#theta_{p}",
    "p_{p}" },
};

// delta_pT in delta_alphaT slices
std::map< double, std::vector<double> > pT_in_alphaT_edges = {
  { 0.,  { 0., 0.06, 0.12, 0.18, 0.24, 0.32, 0.4, 0.48 } },
  { 45., { 0., 0.06, 0.12, 0.18, 0.24, 0.32, 0.4, 0.48, 0.55 } },
  { 90., { 0., 0.06, 0.12, 0.18, 0.24, 0.32, 0.4, 0.48, 0.55, 0.63, 0.7 } },
  { 135., { 0., 0.06, 0.12, 0.18, 0.24, 0.32, 0.4, 0.5, 0.6, 0.72, 0.9 } },
  { 180., {} },
};

// delta_alphaT in delta_pT slices
std::map< double, std::vector<double> > alphaT_in_pT_edges = {
  { 0.0, { 0., 25., 60., 95., 120., 145., 165., 180. } },
  { 0.2, { 0., 25., 60., 95., 120., 145., 165., 180. } },
  { 0.3, { 0., 25., 60., 95., 120., 145., 165., 180. } },
  { 0.4, { 0., 25., 60., 95., 120., 145., 165., 180. } },
};

std::vector< double > pT_1D_edges = { 0., 0.06, 0.12, 0.18, 0.24, 0.32,
  0.4, 0.48, 0.55, 0.68, 0.75, 0.9 };

std::vector< double > pTx_1D_edges = { -DBL_MAX, -0.6, -0.45, -0.35, -0.25,
  -0.15, -0.075, 0, 0.075, 0.15, 0.25, 0.35, 0.45, 0.6, DBL_MAX };

// delta_pTx in delta_pTy slices
std::map< double, std::vector<double> > pTx_in_pTy_edges = {
  { -0.15, { -0.6, -0.45, -0.35, -0.25, -0.15, -0.075, 0, 0.075,
    0.15, 0.25, 0.35, 0.45, 0.6 } },
  { 0.15, { -0.6, -0.45, -0.35, -0.25, -0.15, -0.075, 0, 0.075, 0.15, 0.25,
    0.35, 0.45, 0.6 } },
  { DBL_MAX, { -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4 } },
};

// Muon-proton opening angle (deg)
std::vector< double > theta_mu_p_1D_edges = { 0., 30., 40., 50., 60., 70.,
  80., 90., 100., 110., 120., 130., 140., 150., 180. };

// Reconstructed initial nucleon momentum (pn)
std::vector< double > pn_1D_edges = { 0., 0.07, 0.14, 0.21, 0.28, 0.35,
  0.45, 0.54, 0.66, 0.77, 0.9 };

// theta_mu_p in pn slices
std::map< double, std::vector<double> > theta_mu_p_in_pn_edges = {
  { 0., { 0., 60., 70., 80., 90., 100., 110., 120., 130.,
    140., 150., 180. } },
  { 0.21, { 0., 45., 60., 75., 90., 100., 110., 120., 130.,
    140., 150., 180. } },
  { 0.45, { 0., 30., 45., 60., 75., 90., 105., 120., 135.,
    150., 180. } },
};

// cos_theta_mu in 1D
std::vector< double > cos_theta_mu_1D_edges = { -1., -0.85, -0.775, -0.7,
  -0.625, -0.55, -0.475, -0.4, -0.325, -0.25, -0.175, -0.1, -0.025, 0.05,
   0.125, 0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.85,
   0.875, 0.9, 0.925, 0.950, 0.975, 1. };

// cos_theta_p in 1D
std::vector< double > cos_theta_p_1D_edges = { -1., -0.9, -0.75, -0.6,
  -0.45, -0.3, -0.15, 0.0, 0.15, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9,
  0.925, 0.95, 0.975, 1.0 };

// p_p in 1D
std::vector< double > pp_1D_edges = { 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
  0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0 };

// p_mu in 1D
std::vector< double > pmu_1D_edges = { 0.1, 0.175, 0.2, 0.225, 0.25,
  0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.55, 0.6,
  0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2 };

void make_config_all() {

  // Set up an initially empty container to hold the slice definitions. We'll
  // populate it in parallel with defining the bins themselves.
  SliceBinning sb;

  // Set the variables to use when defining phase-space slices
  sb.slice_vars_ = {
    { "p_{#mu}", "p_{\\mu}", "GeV", "\\text{GeV}" },
    { "cos#theta_{#mu}", "\\cos\\theta_{\\mu}", "", "" },
    { "p_{p}", "p_{p}", "GeV", "\\text{GeV}" },
    { "cos#theta_{p}", "\\cos\\theta_{p}", "", "" },
    { "#deltap_{T}", "\\delta p_{T}", "GeV", "\\text{GeV}" },
    { "#delta#alpha_{T}", "\\delta\\alpha_{T}", "deg", "\\text{deg}" },
    { "#deltap_{Tx}", "\\delta p_{T_x}", "GeV", "\\text{GeV}" },
    { "#deltap_{Ty}", "\\delta p_{T_y}", "GeV", "\\text{GeV}" },
    { "#theta_{#mu,p}", "\\theta_{\\mu,p}", "deg", "\\text{deg}" },
    { "p_{n}", "p_{n}", "GeV", "\\text{GeV}" },
    { "cos#theta_{#mu}", "\\cos\\theta_{\\mu}", "", "" },
    { "bin number", "\\text{ bin number}", "", "" }
  };

  std::string selection = "sel_CCNp0pi";
  std::string signal_def = "mc_is_signal";

  // By construction, MC event categories 5-11 contain all beam-correlated
  // backgrounds. This list is therefore comprehensive apart from cosmic
  // overlay stuff which is directly measured in a dedicated sample.
  std::vector< std::string > background_defs = {
    "category == 5", "category == 6", "category == 7", "category == 8",
    "category == 9", "category == 10", "category == 11"
  };

  std::vector< TrueBin > true_bins;
  std::vector< RecoBin > reco_bins;

//// Blocks for 2D momentum/angle measurements for the muon and leading proton

  // Configure kinematic limits for all of the signal bins. Assign them to
  // blocks (for unfolding purposes) using an index which is incremented
  // as we move to each new EdgeDefinition object.
  int block_index = -1;

  for ( const auto& edge_def : MOM_ANG_2D_EDGES ) {

    // Get the indices for the "active" and "other" variables. We will use
    // these to make new slices while also defining the bins
    int act_var_idx = find_slice_var_index( edge_def.active_var_name_,
      sb.slice_vars_ );

    int oth_var_idx = find_slice_var_index( edge_def.other_var_name_,
      sb.slice_vars_ );

    // Increment the block index for the current set of bin edges
    ++block_index;

    // Get the index of the first analysis bin in this block
    int first_block_bin_idx = reco_bins.size();

    // Get an iterator to the last map element. They are sorted numerically, so
    // this will be the upper edge of the last non-overflow bin.
    auto last = edge_def.bin_edges_2d_->cend();
    --last;

    const std::string& var_name = edge_def.branch_var_name_;

    for ( auto iter = edge_def.bin_edges_2d_->cbegin(); iter != last; ++iter )
    {
      // Get an iterator to the map element after the current one. Due to the
      // automatic sorting, this is guaranteed to contain the upper edge of the
      // current momentum bin.
      auto next = iter;
      ++next;

      // Get the current momentum bin limits
      double p_low = iter->first;
      double p_high = next->first;

      // Now iterate over the cosine bins associated with the current momentum
      // bin. Note that we will skip any situations in which the binning is
      // undefined (i.e., because there are less than two bin edges given)
      const auto& cosine_bin_edges = iter->second;

      size_t num_cosine_edges = cosine_bin_edges.size();
      size_t num_cosine_bins = 0u;
      if ( num_cosine_edges >= 2u ) num_cosine_bins = num_cosine_edges - 1u;

      // Before defining each bin, make a new Slice object and set up the
      // corresponding ROOT histogram within it
      auto& cur_slice = add_slice( sb, cosine_bin_edges, act_var_idx,
        oth_var_idx, p_low, p_high );

      for ( size_t b = 0u; b < num_cosine_bins; ++b ) {

        double cos_low = cosine_bin_edges.at( b );
        double cos_high = cosine_bin_edges.at( b + 1u );

        std::stringstream true_ss;
        true_ss << signal_def
          << " && mc_" << var_name << ".Mag() >= " << p_low
          << " && mc_" << var_name << ".Mag() < " << p_high
          << " && mc_" << var_name << ".CosTheta() >= " << cos_low
          << " && mc_" << var_name << ".CosTheta() < " << cos_high;

        std::string true_bin_def = true_ss.str();

        true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

        std::stringstream reco_ss;
        reco_ss << selection
          << " && " << var_name << ".Mag() >= " << p_low
          << " && " << var_name << ".Mag() < " << p_high
          << " && " << var_name << ".CosTheta() >= " << cos_low
          << " && " << var_name << ".CosTheta() < " << cos_high;

        std::string reco_bin_def = reco_ss.str();

        // Here we use a trick: the current analysis bin index is equal
        // to the size of the reco_bins vector before we add the new element.
        size_t ana_bin_idx = reco_bins.size();
        // Here's another trick: the call to operator[]() below will create
        // a new map entry if needed. We then insert the current analysis
        // bin index into the map entry.
        cur_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );

        // Define the new bin and add it to the vector of reco bins
        reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );

        // We don't need an overflow cosine bin because the entire angular
        // range is covered. We'll use a single bin for the momentum overflow
        // if needed.

      } // loop over cosine bins

    } // loop over momentum bins

    if ( edge_def.needs_overflow_ ) {

      // Get the lower limit for the momentum overflow bin
      double p_overflow_min = last->first;

      // Create a dedicated slice for the overflow bin
      auto& cur_slice = add_slice( sb, 1, -1., 1., act_var_idx,
        oth_var_idx, p_overflow_min, p_overflow_min );

      // Use the same tricks as above to set the bin indices for the slice
      size_t ana_bin_idx = reco_bins.size();
      cur_slice.bin_map_[ 1 ].insert( ana_bin_idx );

      // Create the single momentum overflow bin in both true and reco space
      std::stringstream true_ss;
      true_ss << signal_def << " && mc_" + var_name + ".Mag() >= "
        << p_overflow_min;

      std::string true_bin_def = true_ss.str();
      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

      std::stringstream reco_ss;
      reco_ss << selection << " && " + var_name + ".Mag() >= "
        << p_overflow_min;

      std::string reco_bin_def = reco_ss.str();
      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );
    }

    // Save a vector of the momentum bin edges. We need them contiguous
    // in memory to create a new TH1D using them below.
    std::vector< double > mom_bin_edges;
    for ( const auto& edge_pair : *edge_def.bin_edges_2d_ ) {
      double mom = edge_pair.first;
      mom_bin_edges.push_back( mom );
    }

    // Create a slice in which we've integrated over all angular bins
    // to obtain a single-differential cross section in momentum space
    auto& cur_slice = add_slice( sb, mom_bin_edges, oth_var_idx );

    // Pull the angular bins from the previous slices for the current
    // momentum variable. Given how we've generated the slices, they will
    // appear in order and have a one-to-one mapping to bins in the current
    // ROOT histogram.
    int mom_root_bin_idx = 1;
    for ( const auto& temp_slice : sb.slices_ ) {

      const auto& other_vars = temp_slice.other_vars_;

      // Skip existing slices which have no "other" variable defined
      // (this will include the one we're creating itself)
      if ( other_vars.empty() ) continue;
      const auto& other_var_spec = other_vars.front();

      // Skip existing slices which use a different momentum variable
      if ( other_var_spec.var_index_ != oth_var_idx ) continue;

      // Also skip the overflow bin (as signaled by equal lower and upper
      // limits for the "other" variable (i.e., the momentum)
      double plow = other_var_spec.low_bin_edge_;
      double phigh = other_var_spec.high_bin_edge_;
      if ( plow == phigh ) continue;

      for ( const auto& sl_bin_pair : temp_slice.bin_map_ ) {
        const auto& ana_bin_set = sl_bin_pair.second;
        for ( const size_t ana_bin_idx : ana_bin_set ) {
          cur_slice.bin_map_[ mom_root_bin_idx ].insert( ana_bin_idx );
        }
      }

      // Move onward to the next momentum bin
      ++mom_root_bin_idx;
    }

    // Get the index of the final analysis bin in this block
    int last_block_bin_idx = reco_bins.size() - 1;

    // Create a slice in which all bins in the current block are shown as a
    // function of bin number
    int num_block_bins = last_block_bin_idx - first_block_bin_idx + 1;
    int bin_number_var_idx = find_slice_var_index( "bin number",
      sb.slice_vars_ );

    auto& bin_num_slice = add_slice( sb, num_block_bins, first_block_bin_idx,
      last_block_bin_idx + 1, bin_number_var_idx );
    for ( int ab = first_block_bin_idx; ab <= last_block_bin_idx; ++ab ) {
      // The ROOT histogram bins are one-based, so we correct for this here
      bin_num_slice.bin_map_[ ab + 1 - first_block_bin_idx ].insert( ab );
    }

    // For the proton measurement, also make a slice integrated over both
    // momentum and angle
    if ( edge_def.branch_var_name_ != "p3_lead_p" ) continue;
    int pp_var_idx = find_slice_var_index( "p_{p}", sb.slice_vars_ );
    int cosp_var_idx = find_slice_var_index( "cos#theta_{p}",
      sb.slice_vars_ );
    auto& total_slice = add_slice( sb, 1, -1., 1., cosp_var_idx, pp_var_idx,
      PROTON_2D_BIN_EDGES.cbegin()->first,
      PROTON_2D_BIN_EDGES.crbegin()->first );
    for ( int ab = first_block_bin_idx; ab <= last_block_bin_idx; ++ab ) {
      total_slice.bin_map_[ 1 ].insert( ab );
    }

  }

//// 1D delta_pT block

  // Move into the next block
  ++block_index;

  // Get the index for the delta_pT variable definition. We will use
  // it to make a new slice while also defining the 1D bins.
  int pT_var_idx = find_slice_var_index( "#deltap_{T}", sb.slice_vars_ );
  int alphaT_var_idx = find_slice_var_index(
    "#delta#alpha_{T}", sb.slice_vars_ );

  // Get the index of the first analysis bin in this block
  int first_block_bin_idx = reco_bins.size();

  size_t num_pT_1D_edges = pT_1D_edges.size();
  size_t num_pT_1D_bins = 0u;
  if ( num_pT_1D_edges >= 2u ) num_pT_1D_bins = num_pT_1D_edges - 1u;

  // Before defining each bin, make a new Slice object and set up the
  // corresponding ROOT histogram within it
  auto& temp_slice1 = add_slice( sb, pT_1D_edges, pT_var_idx );

  for ( size_t b = 0u; b < num_pT_1D_bins; ++b ) {

    double pT_low = pT_1D_edges.at( b );
    double pT_high = pT_1D_edges.at( b + 1u );

    std::stringstream true_ss;
    true_ss << signal_def
      << " && mc_delta_pT >= " << pT_low << " && mc_delta_pT < " << pT_high;

    std::string true_bin_def = true_ss.str();

    true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

    std::stringstream reco_ss;
    reco_ss << selection
      << " && delta_pT >= " << pT_low << " && delta_pT < " << pT_high;

    std::string reco_bin_def = reco_ss.str();

    // Here we use a trick: the current analysis bin index is equal
    // to the size of the reco_bins vector before we add the new element.
    size_t ana_bin_idx = reco_bins.size();
    // Here's another trick: the call to operator[]() below will create
    // a new map entry if needed. We then insert the current analysis
    // bin index into the map entry.
    temp_slice1.bin_map_[ b + 1 ].insert( ana_bin_idx );

    // Define the new bin and add it to the vector of reco bins
    reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );

  } // loop over 1D delta_pT bins


  // Get the lower limit for the delta_pT overflow bin
  double pT_overflow_min = pT_1D_edges.back();

  // Create a dedicated slice for the overflow bin
  auto& temp_slice2 = add_slice( sb, 1, 0., 180., alphaT_var_idx,
    pT_var_idx, pT_overflow_min, pT_overflow_min );

  // Use the same tricks as above to set the bin indices for the slice
  size_t ana_bin_idx = reco_bins.size();
  temp_slice2.bin_map_[ 1 ].insert( ana_bin_idx );

  // Create the single delta_pT overflow bin in both true and reco space
  std::stringstream true_ss;
  true_ss << signal_def << " && mc_delta_pT >= " << pT_overflow_min;

  std::string true_bin_def = true_ss.str();
  true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

  std::stringstream reco_ss;
  reco_ss << selection << " && delta_pT >= " << pT_overflow_min;

  std::string reco_bin_def = reco_ss.str();
  reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );

  // Get the index of the final analysis bin in this block
  int last_block_bin_idx = reco_bins.size() - 1;

  // Also make a slice integrated over all delta_alphaT values
  auto& total_slice = add_slice( sb, 1, 0., 180., alphaT_var_idx );
  for ( int ab = first_block_bin_idx; ab <= last_block_bin_idx; ++ab ) {
    total_slice.bin_map_[ 1 ].insert( ab );
  }

//// 2D block of delta_pT in delta_alphaT slices

  // Increment the block index for the current set of bin edges
  ++block_index;

  // Get the index of the first analysis bin in this block
  first_block_bin_idx = reco_bins.size();

  // Get an iterator to the last bin edge map element. They are sorted
  // numerically, so this will be the upper edge of the last non-overflow bin.
  auto last = pT_in_alphaT_edges.cend();
  --last;

  for ( auto iter = pT_in_alphaT_edges.cbegin(); iter != last; ++iter )
  {
    // Get an iterator to the map element after the current one. Due to the
    // automatic sorting, this is guaranteed to contain the upper edge of the
    // current delta_alphaT bin
    auto next = iter;
    ++next;

    // Get the current delta_alphaT bin limits
    double aT_low = iter->first;
    double aT_high = next->first;

    // Now iterate over the delta_pT bins associated with the current
    // delta_alphaT bin. Note that we will skip any situations in which the
    // binning is undefined (i.e., because there are less than two bin edges
    // given)
    const auto& pT_bin_edges = iter->second;

    size_t num_pT_edges = pT_bin_edges.size();
    size_t num_pT_bins = 0u;
    if ( num_pT_edges >= 2u ) num_pT_bins = num_pT_edges - 1u;

    // Before defining each bin, make a new Slice object and set up the
    // corresponding ROOT histogram within it
    auto& cur_slice = add_slice( sb, pT_bin_edges, pT_var_idx,
      alphaT_var_idx, aT_low, aT_high );

    for ( size_t b = 0u; b < num_pT_bins; ++b ) {

      double pT_low = pT_bin_edges.at( b );
      double pT_high = pT_bin_edges.at( b + 1u );

      std::stringstream true_ss;
      true_ss << signal_def
        << " && mc_delta_alphaT * 180. / TMath::ACos(-1.) >= " << aT_low
        << " && mc_delta_alphaT * 180. / TMath::ACos(-1.) < " << aT_high
        << " && mc_delta_pT >= " << pT_low
        << " && mc_delta_pT < " << pT_high;

      std::string true_bin_def = true_ss.str();

      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

      std::stringstream reco_ss;
      reco_ss << selection
        << " && delta_alphaT * 180. / TMath::ACos(-1.) >= " << aT_low
        << " && delta_alphaT * 180. / TMath::ACos(-1.) < " << aT_high
        << " && delta_pT >= " << pT_low
        << " && delta_pT < " << pT_high;

      std::string reco_bin_def = reco_ss.str();

      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();
      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry.
      cur_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );

      // Define the new bin and add it to the vector of reco bins
      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );

    } // loop over delta_pT bins

    // Add an overflow bin for delta_pT in the current delta_alphaT slice
    double pT_overflow_min = pT_bin_edges.back();

    std::stringstream true_ss;
    true_ss << signal_def
      << " && mc_delta_alphaT * 180. / TMath::ACos(-1.) >= " << aT_low
      << " && mc_delta_alphaT * 180. / TMath::ACos(-1.) < " << aT_high
      << " && mc_delta_pT >= " << pT_overflow_min;

    std::string true_bin_def = true_ss.str();

    true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

    std::stringstream reco_ss;
    reco_ss << selection
      << " && delta_alphaT * 180. / TMath::ACos(-1.) >= " << aT_low
      << " && delta_alphaT * 180. / TMath::ACos(-1.) < " << aT_high
      << " && delta_pT >= " << pT_overflow_min;

    std::string reco_bin_def = reco_ss.str();

    // Define the new bin and add it to the vector of reco bins
    reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );

  } // loop over delta_alphaT bins

  // Get the index of the final analysis bin in this block
  last_block_bin_idx = reco_bins.size() - 1;

  // Create a slice in which all bins in the current block are shown as a
  // function of bin number
  int num_block_bins = last_block_bin_idx - first_block_bin_idx + 1;
  int bin_number_var_idx = find_slice_var_index( "bin number",
    sb.slice_vars_ );

  auto& block_bin_num_slice = add_slice( sb, num_block_bins,
    first_block_bin_idx, last_block_bin_idx + 1, bin_number_var_idx );

  for ( int block_bin = 0; block_bin < num_block_bins; ++block_bin ) {
    int analysis_bin = first_block_bin_idx + block_bin;
    // The ROOT histogram bins are one-based, so we correct for this here
    block_bin_num_slice.bin_map_[ block_bin + 1 ].insert( analysis_bin );
  }

//// 2D block of delta_alphaT in delta_pT slices

  // Increment the block index for the current set of bin edges
  ++block_index;

  // Create a slice that will be used to store the 1D differential cross
  // section in delta_alphaT (integrated over all of the delta_pT bins)
  auto& at_1d_slice = add_slice( sb, alphaT_in_pT_edges.cbegin()->second,
    alphaT_var_idx );

  // Get the index of the first analysis bin in this block
  first_block_bin_idx = reco_bins.size();

  // Get an iterator to the last bin edge map element. They are sorted
  // numerically, so this will be the upper edge of the last non-overflow bin.
  auto end = alphaT_in_pT_edges.cend();

  for ( auto iter = alphaT_in_pT_edges.cbegin(); iter != end; ++iter )
  {
    // Get an iterator to the map element after the current one. Due to the
    // automatic sorting, this is guaranteed to contain the upper edge of the
    // current delta_pT bin
    auto next = iter;
    ++next;

    // Get the current delta_pT bin limits
    double pT_low = iter->first;
    double pT_high = pT_low;
    if ( next != end ) pT_high = next->first;

    // Now iterate over the delta_alphaT bins associated with the current
    // delta_pT bin. Note that we will skip any situations in which the
    // binning is undefined (i.e., because there are less than two bin edges
    // given)
    const auto& aT_bin_edges = iter->second;

    size_t num_aT_edges = aT_bin_edges.size();
    size_t num_aT_bins = 0u;
    if ( num_aT_edges >= 2u ) num_aT_bins = num_aT_edges - 1u;

    // Before defining each bin, make a new Slice object and set up the
    // corresponding ROOT histogram within it
    auto& cur_slice = add_slice( sb, aT_bin_edges, alphaT_var_idx,
      pT_var_idx, pT_low, pT_high );

    for ( size_t b = 0u; b < num_aT_bins; ++b ) {

      double aT_low = aT_bin_edges.at( b );
      double aT_high = aT_bin_edges.at( b + 1u );

      std::stringstream true_ss;
      true_ss << signal_def << " && mc_delta_pT >= " << pT_low;
      if ( pT_high != pT_low ) {
        true_ss << " && mc_delta_pT < " << pT_high;
      }
      true_ss << " && mc_delta_alphaT * 180. / TMath::ACos(-1.) >= " << aT_low
        << " && mc_delta_alphaT * 180. / TMath::ACos(-1.) < " << aT_high;

      std::string true_bin_def = true_ss.str();

      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

      std::stringstream reco_ss;
      reco_ss << selection << " && delta_pT >= " << pT_low;
      if ( pT_high != pT_low ) {
        reco_ss << " && delta_pT < " << pT_high;
      }
      reco_ss << " && delta_alphaT * 180. / TMath::ACos(-1.) >= " << aT_low
        << " && delta_alphaT * 180. / TMath::ACos(-1.) < " << aT_high;

      std::string reco_bin_def = reco_ss.str();

      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();
      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry.
      cur_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );

      // Also add this bin to the 1D slice of delta_alphaT. Since the same
      // delta_alphaT binning is used throughout every slice in delta_pT,
      // this is convenient to do here.
      at_1d_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );

      // Define the new bin and add it to the vector of reco bins
      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );

    } // loop over delta_pT bins

  } // loop over delta_alphaT bins

  // Get the index of the final analysis bin in this block
  last_block_bin_idx = reco_bins.size() - 1;

  // Create a slice in which all bins in the current block are shown as a
  // function of bin number
  num_block_bins = last_block_bin_idx - first_block_bin_idx + 1;

  auto& block_bin_num_slice2 = add_slice( sb, num_block_bins,
    first_block_bin_idx, last_block_bin_idx + 1, bin_number_var_idx );

  for ( int block_bin = 0; block_bin < num_block_bins; ++block_bin ) {
    int analysis_bin = first_block_bin_idx + block_bin;

    // The ROOT histogram bins are one-based, so we correct for this here
    block_bin_num_slice2.bin_map_[ block_bin + 1 ].insert( analysis_bin );
  }

//// 1D delta_pTx block

  // Increment the block index
  ++block_index;

  // Get the index for the delta_pTx variable definition. We will use
  // it to make a new slice while also defining the 1D bins.
  int pTx_var_idx = find_slice_var_index( "#deltap_{Tx}", sb.slice_vars_ );

  // Get the index of the first analysis bin in this block
  first_block_bin_idx = reco_bins.size();

  size_t num_pTx_1D_edges = pTx_1D_edges.size();
  size_t num_pTx_1D_bins = 0u;
  if ( num_pTx_1D_edges >= 2u ) num_pTx_1D_bins = num_pTx_1D_edges - 1u;

  // Make a copy of the pTx bin edges with the under and overflow parts
  // removed. We will use this to initialize the slice histogram edges
  auto pTx_1D_edges_temp_copy = pTx_1D_edges;
  pTx_1D_edges_temp_copy.pop_back();
  pTx_1D_edges_temp_copy.erase( pTx_1D_edges_temp_copy.begin() );

  // Before defining each bin, make a new Slice object and set up the
  // corresponding ROOT histogram within it
  auto& temp_slice_pTx = add_slice( sb, pTx_1D_edges_temp_copy, pTx_var_idx );

  for ( size_t b = 0u; b < num_pTx_1D_bins; ++b ) {

    double pTx_low = pTx_1D_edges.at( b );
    double pTx_high = pTx_1D_edges.at( b + 1u );

    std::stringstream true_ss;
    true_ss << signal_def;
    if ( pTx_low > -DBL_MAX ) true_ss << " && mc_delta_pTx >= " << pTx_low;
    if ( pTx_high < DBL_MAX ) true_ss << " && mc_delta_pTx < " << pTx_high;

    std::string true_bin_def = true_ss.str();

    true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

    std::stringstream reco_ss;
    reco_ss << selection;
    if ( pTx_low > -DBL_MAX ) reco_ss << " && delta_pTx >= " << pTx_low;
    if ( pTx_high < DBL_MAX ) reco_ss << " && delta_pTx < " << pTx_high;

    std::string reco_bin_def = reco_ss.str();

    // Here we use a trick: the current analysis bin index is equal
    // to the size of the reco_bins vector before we add the new element.
    size_t ana_bin_idx = reco_bins.size();
    // Here's another trick: the call to operator[]() below will create
    // a new map entry if needed. We then insert the current analysis
    // bin index into the map entry.
    // NOTE: we exclude the under/overflow bins here on purpose. This also
    // avoids the need for the usual b + 1 in the [] operator on the bin map
    if ( std::abs(pTx_low) != DBL_MAX && std::abs(pTx_high) != DBL_MAX ) {
      temp_slice_pTx.bin_map_[ b ].insert( ana_bin_idx );
    }

    // Define the new bin and add it to the vector of reco bins
    reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );

  } // loop over 1D delta_pTx bins

//// 2D block delta_pTx in delta_pTy slices
  int pTy_var_idx = find_slice_var_index( "#deltap_{Ty}",
    sb.slice_vars_ );

  // Increment the block index for the current set of bin edges
  ++block_index;

  // Get the index of the first analysis bin in this block
  first_block_bin_idx = reco_bins.size();

  for ( const auto& edge_pair : pTx_in_pTy_edges ) {

    // Get the current delta_pTy bin limits
    // TODO: revisit this super hacky solution
    double pTy_low = edge_pair.first;
    double pTy_high = pTy_low;
    if ( pTy_low == -0.15 ) {
      pTy_low = -DBL_MAX;
    }
    else if ( pTy_low == 0.15 ) {
      pTy_low = -0.15;
    }
    else {
      pTy_low = 0.15;
      pTy_high = DBL_MAX;
    }

    // Now iterate over the delta_pTx bins associated with the current
    // delta_pTy bin. Note that we will skip any situations in which the
    // binning is undefined (i.e., because there are less than two bin edges
    // given)
    const auto& pTx_bin_edges = edge_pair.second;

    int num_pTx_edges = pTx_bin_edges.size();
    int num_pTx_bins = 0;
    if ( num_pTx_edges >= 2 ) num_pTx_bins = num_pTx_edges - 1;

    // Before defining each bin, make a new Slice object and set up the
    // corresponding ROOT histogram within it
    auto& cur_slice = add_slice( sb, pTx_bin_edges, pTx_var_idx,
      pTy_var_idx, pTy_low, pTy_high );

    for ( int b = -1; b < num_pTx_edges; ++b ) {

      double pTx_low = ( b == -1 ? -DBL_MAX : pTx_bin_edges.at(b) );
      double pTx_high = ( b == num_pTx_bins ? DBL_MAX
        : pTx_bin_edges.at(b + 1) );

      std::stringstream true_ss;
      true_ss << signal_def;
      if ( pTx_low != -DBL_MAX ) true_ss << " && mc_delta_pTx >= " << pTx_low;
      if ( pTx_high != DBL_MAX ) true_ss << " && mc_delta_pTx < " << pTx_high;

      if ( pTy_low != -DBL_MAX ) true_ss << " && mc_delta_pTy >= " << pTy_low;
      if ( pTy_high != DBL_MAX ) true_ss << " && mc_delta_pTy < " << pTy_high;

      std::string true_bin_def = true_ss.str();

      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

      std::stringstream reco_ss;
      reco_ss << selection;

      if ( pTx_low != -DBL_MAX ) reco_ss << " && delta_pTx >= " << pTx_low;
      if ( pTx_high != DBL_MAX ) reco_ss << " && delta_pTx < " << pTx_high;

      if ( pTy_low != -DBL_MAX ) reco_ss << " && delta_pTy >= " << pTy_low;
      if ( pTy_high != DBL_MAX ) reco_ss << " && delta_pTy < " << pTy_high;

      std::string reco_bin_def = reco_ss.str();

      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();

      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry. Leave out the under/overflow bins
      // given how the slice is defined.
      if ( b > -1 && b < num_pTx_bins ) {
        cur_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );
      }

      // Define the new bin and add it to the vector of reco bins
      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );

    } // loop over delta_pTx bins

  } // loop over delta_pTy bins

  // Get the index of the final analysis bin in this block
  last_block_bin_idx = reco_bins.size() - 1;

  // Create a slice in which all bins in the current block are shown as a
  // function of bin number
  num_block_bins = last_block_bin_idx - first_block_bin_idx + 1;

  auto& block_bin_num_slice3 = add_slice( sb, num_block_bins,
    first_block_bin_idx, last_block_bin_idx + 1, bin_number_var_idx );

  for ( int block_bin = 0; block_bin < num_block_bins; ++block_bin ) {
    int analysis_bin = first_block_bin_idx + block_bin;

    // The ROOT histogram bins are one-based, so we correct for this here
    block_bin_num_slice3.bin_map_[ block_bin + 1 ].insert( analysis_bin );
  }

//// 1D theta_mu_p block

  // Increment the block index for the current set of bin edges
  ++block_index;

  // Get the index of the first analysis bin in this block
  first_block_bin_idx = reco_bins.size();

  // Get the index for the theta_mu_p variable definition. We will use
  // it to make a new slice while also defining the 1D bins.
  int theta_mu_p_var_idx = find_slice_var_index( "#theta_{#mu,p}",
    sb.slice_vars_ );

  // Use a version of the muon-leading-proton opening angle expressed in degrees
  const std::string mc_theta_mu_p_deg( "mc_theta_mu_p * 180."
    "/ TMath::ACos(-1.)" );

  const std::string reco_theta_mu_p_deg( "theta_mu_p * 180."
    "/ TMath::ACos(-1.)" );

  size_t num_theta_mu_p_1D_edges = theta_mu_p_1D_edges.size();
  size_t num_theta_mu_p_1D_bins = 0u;
  if ( num_theta_mu_p_1D_edges >= 2u ) num_theta_mu_p_1D_bins
    = num_theta_mu_p_1D_edges - 1u;

  // Before defining each bin, make a new Slice object and set up the
  // corresponding ROOT histogram within it
  auto& tmp_slice = add_slice( sb, theta_mu_p_1D_edges, theta_mu_p_var_idx );

  for ( size_t b = 0u; b < num_theta_mu_p_1D_bins; ++b ) {

    double theta_mu_p_low = theta_mu_p_1D_edges.at( b );
    double theta_mu_p_high = theta_mu_p_1D_edges.at( b + 1u );

    std::stringstream true_ss;
    true_ss << signal_def << " && " << mc_theta_mu_p_deg << " >= "
      << theta_mu_p_low << " && " << mc_theta_mu_p_deg << " < "
      << theta_mu_p_high;

    std::string true_bin_def = true_ss.str();

    true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

    std::stringstream reco_ss;
    reco_ss << selection << " && " << reco_theta_mu_p_deg << " >= "
      << theta_mu_p_low << " && " << reco_theta_mu_p_deg << " < "
      << theta_mu_p_high;

    std::string reco_bin_def = reco_ss.str();

    // Here we use a trick: the current analysis bin index is equal
    // to the size of the reco_bins vector before we add the new element.
    size_t ana_bin_idx = reco_bins.size();
    // Here's another trick: the call to operator[]() below will create
    // a new map entry if needed. We then insert the current analysis
    // bin index into the map entry.
    tmp_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );

    // Define the new bin and add it to the vector of reco bins
    reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );

  } // loop over 1D theta_mu_p bins

//// 1D pn block

  // Increment the block index for the current set of bin edges
  ++block_index;

  // Get the index of the first analysis bin in this block
  first_block_bin_idx = reco_bins.size();

  // Get the index for the pn variable definition. We will use
  // it to make a new slice while also defining the 1D bins.
  int pn_var_idx = find_slice_var_index( "p_{n}", sb.slice_vars_ );

  size_t num_pn_1D_edges = pn_1D_edges.size();
  size_t num_pn_1D_bins = 0u;
  if ( num_pn_1D_edges >= 2u ) num_pn_1D_bins = num_pn_1D_edges - 1u;

  // Before defining each bin, make a new Slice object and set up the
  // corresponding ROOT histogram within it
  auto& pn_slice = add_slice( sb, pn_1D_edges, pn_var_idx );

  for ( size_t b = 0u; b < num_pn_1D_edges; ++b ) {

    double pn_low = pn_1D_edges.at( b );
    double pn_high = DBL_MAX;
    if ( b < num_pn_1D_bins ) pn_high = pn_1D_edges.at( b + 1u );

    std::stringstream true_ss;
    true_ss << signal_def << " && mc_pn >= " << pn_low;
    if ( pn_high != DBL_MAX ) true_ss << " && mc_pn < " << pn_high;

    std::string true_bin_def = true_ss.str();

    true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

    std::stringstream reco_ss;
    reco_ss << selection << " && pn >= " << pn_low;
    if ( pn_high != DBL_MAX ) reco_ss << " && pn < " << pn_high;

    std::string reco_bin_def = reco_ss.str();

    // Include all but the overflow bin in a slice of the 1D pn distribution
    if ( pn_high != DBL_MAX ) {
      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();
      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry.
      pn_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );
    }

    // Define the new bin and add it to the vector of reco bins
    reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );

  } // loop over 1D pn bins

//// 2D block of theta_mu_p in pn slices

  // Increment the block index for the current set of bin edges
  ++block_index;

  // Get the index of the first analysis bin in this block
  first_block_bin_idx = reco_bins.size();

  auto begin = theta_mu_p_in_pn_edges.cbegin();
  auto end2 = theta_mu_p_in_pn_edges.cend();

  for ( auto iter = begin; iter != end2; ++iter )
  {
    // Get an iterator to the map element after the current one. Due to the
    // automatic sorting, this will contain the upper edge of the
    // current non-overflow pn bin
    auto next = iter;
    ++next;

    // Get the current delta_alphaT bin limits
    double pn_low = iter->first;
    double pn_high = DBL_MAX;
    if ( next != end2 ) pn_high = next->first;

    // Now iterate over the theta_mu_p bins associated with the current
    // pn bin. Note that we will skip any situations in which the
    // binning is undefined (i.e., because there are less than two bin edges
    // given)
    const auto& theta_mu_p_bin_edges = iter->second;

    size_t num_theta_mu_p_edges = theta_mu_p_bin_edges.size();
    size_t num_theta_mu_p_bins = 0u;
    if ( num_theta_mu_p_edges >= 2u ) {
      num_theta_mu_p_bins = num_theta_mu_p_edges - 1u;
    }

    // Before defining each bin, make a new Slice object and set up the
    // corresponding ROOT histogram within it
    auto& cur_slice = add_slice( sb, theta_mu_p_bin_edges, theta_mu_p_var_idx,
      pn_var_idx, pn_low, pn_high );

    for ( size_t b = 0u; b < num_theta_mu_p_bins; ++b ) {

      double theta_mu_p_low = theta_mu_p_bin_edges.at( b );
      double theta_mu_p_high = theta_mu_p_bin_edges.at( b + 1u );

      std::stringstream true_ss;
      true_ss << signal_def << " && " << mc_theta_mu_p_deg << " >= "
        << theta_mu_p_low << " && " << mc_theta_mu_p_deg << " < "
        << theta_mu_p_high << " && mc_pn >= " << pn_low;

      if ( pn_high != DBL_MAX ) true_ss << " && mc_pn < " << pn_high;

      std::string true_bin_def = true_ss.str();

      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

      std::stringstream reco_ss;
      reco_ss << selection << " && " << reco_theta_mu_p_deg << " >= "
        << theta_mu_p_low << " && " << reco_theta_mu_p_deg << " < "
        << theta_mu_p_high << " && pn >= " << pn_low;

      if ( pn_high != DBL_MAX ) reco_ss << " && pn < " << pn_high;

      std::string reco_bin_def = reco_ss.str();

      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();
      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry.
      cur_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );

      // Define the new bin and add it to the vector of reco bins
      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );

    } // loop over theta_mu_p bins

  } // loop over pn bins

  // Get the index of the final analysis bin in this block
  last_block_bin_idx = reco_bins.size() - 1;

  // Create a slice in which all bins in the current block are shown as a
  // function of bin number
  num_block_bins = last_block_bin_idx - first_block_bin_idx + 1;
  bin_number_var_idx = find_slice_var_index( "bin number",
    sb.slice_vars_ );

  auto& block_bin_num_slice4 = add_slice( sb, num_block_bins,
    first_block_bin_idx, last_block_bin_idx + 1, bin_number_var_idx );

  for ( int block_bin = 0; block_bin < num_block_bins; ++block_bin ) {
    int analysis_bin = first_block_bin_idx + block_bin;
    // The ROOT histogram bins are one-based, so we correct for this here
    block_bin_num_slice4.bin_map_[ block_bin + 1 ].insert( analysis_bin );
  }

//// 1D cos_theta_mu block

  // Increment the block index for the current set of bin edges
  ++block_index;

  // Get the index of the first analysis bin in this block
  first_block_bin_idx = reco_bins.size();

  // Get the index for the cos_theta_mu variable definition. We will use
  // it to make a new slice while also defining the 1D bins.
  int cos_var_idx = find_slice_var_index( "cos#theta_{#mu}",
    sb.slice_vars_ );

  size_t num_cos_1D_edges = cos_theta_mu_1D_edges.size();
  size_t num_cos_1D_bins = 0u;
  if ( num_cos_1D_edges >= 2u ) num_cos_1D_bins
    = num_cos_1D_edges - 1u;

  // Before defining each bin, make a new Slice object and set up the
  // corresponding ROOT histogram within it
  auto& cos_slice = add_slice( sb, cos_theta_mu_1D_edges, cos_var_idx );

  for ( size_t b = 0u; b < num_cos_1D_bins; ++b ) {

    double cos_low = cos_theta_mu_1D_edges.at( b );
    double cos_high = cos_theta_mu_1D_edges.at( b + 1 );

    std::stringstream true_ss;
    true_ss << signal_def << " && mc_p3_mu.CosTheta() >= " << cos_low
     << " && mc_p3_mu.CosTheta() < " << cos_high;

    std::string true_bin_def = true_ss.str();

    true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

    std::stringstream reco_ss;
    reco_ss << selection << " && p3_mu.CosTheta() >= " << cos_low
     << " && p3_mu.CosTheta() < " << cos_high;

    std::string reco_bin_def = reco_ss.str();

    // Here we use a trick: the current analysis bin index is equal
    // to the size of the reco_bins vector before we add the new element.
    size_t ana_bin_idx = reco_bins.size();
    // Here's another trick: the call to operator[]() below will create
    // a new map entry if needed. We then insert the current analysis
    // bin index into the map entry.
    cos_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );

    // Define the new bin and add it to the vector of reco bins
    reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );

  } // loop over 1D cos_theta_mu bins

//// 1D cos_theta_p block

  // Increment the block index for the current set of bin edges
  ++block_index;

  // Get the index of the first analysis bin in this block
  first_block_bin_idx = reco_bins.size();

  // Get the index for the cos_theta_p variable definition. We will use
  // it to make a new slice while also defining the 1D bins.
  int cosp_var_idx = find_slice_var_index( "cos#theta_{p}",
    sb.slice_vars_ );

  size_t num_cosp_1D_edges = cos_theta_p_1D_edges.size();
  size_t num_cosp_1D_bins = 0u;
  if ( num_cosp_1D_edges >= 2u ) num_cosp_1D_bins
    = num_cosp_1D_edges - 1u;

  // Before defining each bin, make a new Slice object and set up the
  // corresponding ROOT histogram within it
  auto& cosp_slice = add_slice( sb, cos_theta_p_1D_edges, cosp_var_idx );

  for ( size_t b = 0u; b < num_cosp_1D_bins; ++b ) {

    double cosp_low = cos_theta_p_1D_edges.at( b );
    double cosp_high = cos_theta_p_1D_edges.at( b + 1 );

    std::stringstream true_ss;
    true_ss << signal_def << " && mc_p3_lead_p.CosTheta() >= " << cosp_low
     << " && mc_p3_lead_p.CosTheta() < " << cosp_high;

    std::string true_bin_def = true_ss.str();

    true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

    std::stringstream reco_ss;
    reco_ss << selection << " && p3_lead_p.CosTheta() >= " << cosp_low
     << " && p3_lead_p.CosTheta() < " << cosp_high;

    std::string reco_bin_def = reco_ss.str();

    // Here we use a trick: the current analysis bin index is equal
    // to the size of the reco_bins vector before we add the new element.
    size_t ana_bin_idx = reco_bins.size();
    // Here's another trick: the call to operator[]() below will create
    // a new map entry if needed. We then insert the current analysis
    // bin index into the map entry.
    cosp_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );

    // Define the new bin and add it to the vector of reco bins
    reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );

  } // loop over 1D cos_theta_p bins

//// 1D p_p block

  // Increment the block index for the current set of bin edges
  ++block_index;

  // Get the index of the first analysis bin in this block
  first_block_bin_idx = reco_bins.size();

  // Get the index for the p_p variable definition. We will use
  // it to make a new slice while also defining the 1D bins.
  int pp_var_idx = find_slice_var_index( "p_{p}", sb.slice_vars_ );

  size_t num_pp_1D_edges = pp_1D_edges.size();
  size_t num_pp_1D_bins = 0u;
  if ( num_pp_1D_edges >= 2u ) num_pp_1D_bins
    = num_pp_1D_edges - 1u;

  // Before defining each bin, make a new Slice object and set up the
  // corresponding ROOT histogram within it
  auto& pp_slice = add_slice( sb, pp_1D_edges, pp_var_idx );

  for ( size_t b = 0u; b < num_pp_1D_bins; ++b ) {

    double pp_low = pp_1D_edges.at( b );
    double pp_high = pp_1D_edges.at( b + 1 );

    std::stringstream true_ss;
    true_ss << signal_def << " && mc_p3_lead_p.Mag() >= " << pp_low
     << " && mc_p3_lead_p.Mag() < " << pp_high;

    std::string true_bin_def = true_ss.str();

    true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

    std::stringstream reco_ss;
    reco_ss << selection << " && p3_lead_p.Mag() >= " << pp_low
     << " && p3_lead_p.Mag() < " << pp_high;

    std::string reco_bin_def = reco_ss.str();

    // Here we use a trick: the current analysis bin index is equal
    // to the size of the reco_bins vector before we add the new element.
    size_t ana_bin_idx = reco_bins.size();
    // Here's another trick: the call to operator[]() below will create
    // a new map entry if needed. We then insert the current analysis
    // bin index into the map entry.
    pp_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );

    // Define the new bin and add it to the vector of reco bins
    reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );

  } // loop over 1D p_p bins


//// 1D p_mu block

  // Increment the block index for the current set of bin edges
  ++block_index;

  // Get the index of the first analysis bin in this block
  first_block_bin_idx = reco_bins.size();

  // Get the index for the p_mu variable definition. We will use
  // it to make a new slice while also defining the 1D bins.
  int pmu_var_idx = find_slice_var_index( "p_{#mu}", sb.slice_vars_ );

  size_t num_pmu_1D_edges = pmu_1D_edges.size();
  size_t num_pmu_1D_bins = 0u;
  if ( num_pmu_1D_edges >= 2u ) num_pmu_1D_bins
    = num_pmu_1D_edges - 1u;

  // Before defining each bin, make a new Slice object and set up the
  // corresponding ROOT histogram within it
  auto& pmu_slice = add_slice( sb, pmu_1D_edges, pmu_var_idx );

  for ( size_t b = 0u; b < num_pmu_1D_bins; ++b ) {

    double pmu_low = pmu_1D_edges.at( b );
    double pmu_high = pmu_1D_edges.at( b + 1 );

    std::stringstream true_ss;
    true_ss << signal_def << " && mc_p3_mu.Mag() >= " << pmu_low
     << " && mc_p3_mu.Mag() < " << pmu_high;

    std::string true_bin_def = true_ss.str();

    true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

    std::stringstream reco_ss;
    reco_ss << selection << " && p3_mu.Mag() >= " << pmu_low
     << " && p3_mu.Mag() < " << pmu_high;

    std::string reco_bin_def = reco_ss.str();

    // Here we use a trick: the current analysis bin index is equal
    // to the size of the reco_bins vector before we add the new element.
    size_t ana_bin_idx = reco_bins.size();
    // Here's another trick: the call to operator[]() below will create
    // a new map entry if needed. We then insert the current analysis
    // bin index into the map entry.
    pmu_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );

    // Define the new bin and add it to the vector of reco bins
    reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );

  } // loop over 1D p_mu bins

//// Slice showing all blocks of ordinary reco bins / signal true bins

  //// Create a slice showing all blocks together as a function of bin number
  //int num_ord_bins = reco_bins.size();

  //auto& ord_bin_num_slice = add_slice( sb, num_ord_bins, 0, num_ord_bins,
  //  bin_number_var_idx );
  //for ( int ab = 0; ab < num_ord_bins; ++ab ) {
  //  // The ROOT histogram bins are one-based, so we correct for this here
  //  ord_bin_num_slice.bin_map_[ ab + 1 ].insert( ab );
  //}

//// Sideband reco bins

  //// Clone each of the ordinary reco bins and apply the sideband selection
  //for ( int ob = 0; ob < num_ord_bins; ++ob ) {
  //  const auto& rbin = reco_bins.at( ob );

  //  // Strip out the initial selection flag ("sel_CCNp0pi") and replace
  //  // it in the selection cuts with the sideband selection
  //  size_t drop_pos = rbin.selection_cuts_.find( "&&" );
  //  std::string sel_new = sel_combined + ' '
  //    + rbin.selection_cuts_.substr( drop_pos );
  //  reco_bins.emplace_back( sel_new, kSidebandRecoBin, DUMMY_BLOCK_INDEX );
  //}

//// Final definitions

  int num_bins = reco_bins.size();

  //// Create a slice showing just the sideband bins as a function of bin number
  //auto& sideband_slice = add_slice( sb, num_ord_bins, num_ord_bins, num_bins,
  //  bin_number_var_idx );
  //for ( int ab = num_ord_bins; ab < num_bins; ++ab ) {
  //  // The ROOT histogram bins are one-based, so we correct for this here
  //  sideband_slice.bin_map_[ ab + 1 ].insert( ab );
  //}

  // Create a slice showing all results together as a function of bin number
  auto& bin_num_slice = add_slice( sb, num_bins, 0, num_bins,
    bin_number_var_idx );
  for ( int ab = 0; ab < num_bins; ++ab ) {
    // The ROOT histogram bins are one-based, so we correct for this here
    bin_num_slice.bin_map_[ ab + 1 ].insert( ab );
  }

  // Add true bins for the background categories of interest. We'll use a
  // dummy block index now since it's only important for signal bins.
  for ( const auto& bdef : background_defs ) {
    true_bins.emplace_back( bdef, kBackgroundTrueBin, DUMMY_BLOCK_INDEX );
  }

  // Dump this information to the output files
  std::ofstream out_file( "myconfig_all.txt" );
  out_file << "ALL\n";
  out_file << "stv_tree\n";
  out_file << true_bins.size() << '\n';
  for ( const auto& tb : true_bins ) out_file << tb << '\n';

  out_file << reco_bins.size() << '\n';
  for ( const auto& rb : reco_bins ) out_file << rb << '\n';

  // Also write a SliceBinning configuration file
  std::ofstream sb_file( "mybins_all.txt" );
  sb_file << sb;
}
