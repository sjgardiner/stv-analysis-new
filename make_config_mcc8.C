#include "ConfigMakerUtils.hh"
#include "UniverseMaker.hh"
#include "SliceBinning.hh"

constexpr int DUMMY_BLOCK_INDEX = -1;

void make_config_mcc8() {

  // Keys are the reco STV ntuple branch names of interest. Values
  // are vectors of bin edges (taken from the MCC8 CCNp0pi paper)
  std::map< std::string, std::vector<double> > mcc8_bin_edge_map = {

    { "p3_mu.CosTheta()", { -1.0, -0.82, -0.66, -0.39, -0.16, 0.05, 0.25, 0.43,
      0.59, 0.73, 0.83, 0.91, 1.0 } },

    { "theta_mu_p", { 0.0, 0.8, 1.2, 1.57, 1.94, 2.34, M_PI } },

    { "p3_lead_p.Mag()", { 0.3, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87,
      0.93, 1.2 } },

    { "p3_mu.Mag()", { 0.1, 0.18, 0.3, 0.48, 0.75, 1.14, 2.5 } },

    { "p3_lead_p.CosTheta()", { -1.0, -0.5, 0.0, 0.27, 0.45, 0.62, 0.76, 0.86,
      0.94, 1.0 } },

  };

  // Used for converting from the variable names given in the bin edge
  // map to the ones used by the SliceBinning object
  std::map< std::string, std::string > var_name_map = {
    { "p3_mu.CosTheta()", "reco cos#theta_{#mu}" },
    { "theta_mu_p", "reco #theta_{#mup}" },
    { "p3_lead_p.Mag()", "reco p_{p}" },
    { "p3_mu.Mag()", "reco p_{#mu}" },
    { "p3_lead_p.CosTheta()", "reco cos#theta_{p}" },
  };

  // Keys are the same reco variable branch expressions as above. Values
  // are bool pairs indicating whether an underflow and overflow bin
  // should be produced.
  std::map< std::string, std::pair<bool,bool> > mcc8_under_overflow_map = {

    // Restricted by valid angular ranges to lie within the bin limits
    { "p3_mu.CosTheta()", { false, false } },

    { "p3_lead_p.CosTheta()", { false, false } },

    { "theta_mu_p", { false, false } },

    // Restricted by the signal definition to lie within the proton momentum
    // bins
    { "p3_lead_p.Mag()", { false, false } },

    // No underflow bin is needed due to the signal definition. No upper limit
    // is imposed, however, so we will make an overflow bin.
    { "p3_mu.Mag()", { false, true } }

  };

  // Set up an initially empty container to hold the slice definitions. We'll
  // populate it in parallel with defining the bins themselves.
  SliceBinning sb;

  // Set the variables to use when defining phase-space slices
  sb.slice_vars_ = {
    { "reco p_{#mu}", "GeV/c", "reco $p_{\\mu}$", "GeV$/c$" },
    { "reco cos#theta_{#mu}", "", "reco $\\cos\\theta_{\\mu}$", "" },
    { "reco p_{p}", "GeV/c", "reco $p_{\\mu}$", "GeV$/c$" },
    { "reco cos#theta_{p}", "", "reco $\\cos\\theta_{\\mu}$", "" },
    { "reco #theta_{#mup}", "", "reco $\\theta_{\\mu p}$", "" },
    { "reco bin number", "", "reco bin number", "" }
  };

  // NOTE: this script assumes that the definitions for the selection flag
  // (sel_CCNp0pi) and signal flag (mc_is_signal) have been suitably changed to
  // match the MCC8 analysis. The user is responsible for using ntuples that
  // have been post-processed consistently. Strictly speaking, this includes
  // the tiny difference between the pionless and mesonless signal definitions
  // (a ~0.06% effect).
  std::string selection = "sel_CCNp0pi";
  std::string signal_def = "mc_is_signal";

  // By construction, MC event categories 5-11 contain all beam-correlated
  // backgrounds. This list is therefore comprehensive apart from cosmic
  // overlay stuff which is directly measured in a dedicated sample.
  // NOTE: We add an extra background bin for events that the MCC9 analysis
  // considers signal and the MCC8 analysis considers background.
  std::vector< std::string > background_defs = {
    "category == 5", "category == 6", "category == 7", "category == 8",
    "category == 9", "category == 10", "category == 11"
  };

  std::vector< TrueBin > true_bins;
  std::vector< RecoBin > reco_bins;

  // Create separate blocks of bins for each kinematic variable using the
  // bin definitions from the MCC8 CCNp0pi analysis
  int block_idx = -1;
  for ( const auto& pair : mcc8_bin_edge_map ) {
    // Start a new block of related bins
    ++block_idx;

    std::string reco_branchexpr = pair.first;
    std::string true_branchexpr = "mc_" + reco_branchexpr;

    // Get the index for the "active" variable in the current block. We will
    // use it below to make a new slice while also defining the bins in the
    // block.
    const std::string& act_var_name = var_name_map.at( reco_branchexpr );
    int act_var_idx = find_slice_var_index( act_var_name, sb.slice_vars_ );

    const auto flag_pair = mcc8_under_overflow_map.at( reco_branchexpr );
    bool needs_underflow_bin = flag_pair.first;
    bool needs_overflow_bin = flag_pair.second;

    // Require at least two bin edges to be present in the input vector.
    // Any variables for which this is not true will be skipped entirely.
    const auto& bin_edges = pair.second;

    size_t num_edges = bin_edges.size();
    size_t num_bins = 0u;
    if ( num_edges >= 2u ) num_bins = num_edges - 1u;
    else continue;

    // Before defining each bin, make a new Slice object and set up the
    // corresponding ROOT histogram within it
    auto& cur_slice = add_slice( sb, bin_edges, act_var_idx );

    // If needed, then create the underflow bin in both true and reco space
    if ( needs_underflow_bin ) {
      double var_underflow_max = bin_edges.front();

      std::stringstream true_ss;
      true_ss << signal_def << " && " << true_branchexpr
        << " < " << var_underflow_max;

      std::string true_bin_def = true_ss.str();
      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_idx );

      std::stringstream reco_ss;
      reco_ss << selection << " && " << reco_branchexpr
        << " < " << var_underflow_max;

      std::string reco_bin_def = reco_ss.str();

      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();
      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry. Note that the ROOT histogram bin
      // indices are one-based, so the underflow bin is always at index zero.
      cur_slice.bin_map_[ 0 ].insert( ana_bin_idx );

      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_idx );
    }

    // Create the ordinary signal bins using the requested edges
    for ( size_t b = 0u; b < num_bins; ++b ) {

      double var_low = bin_edges.at( b );
      double var_high = bin_edges.at( b + 1u );

      std::stringstream true_ss;
      true_ss << signal_def
        << " && " << true_branchexpr << " >= " << var_low
        << " && " << true_branchexpr << " < "  << var_high;

      std::string true_bin_def = true_ss.str();

      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_idx );

      std::stringstream reco_ss;
      reco_ss << selection
        << " && " << reco_branchexpr << " >= " << var_low
        << " && " << reco_branchexpr << " < "  << var_high;

      std::string reco_bin_def = reco_ss.str();

      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();
      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry. Note that the ROOT histogram bin
      // indices are one-based, so we correct for that in the line below.
      cur_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );

      // Add the completed reco bin definition to the vector
      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_idx );

    } // loop over ordinary bins for the current variable

    // If needed, then create the overflow bin in both true and reco space
    if ( needs_overflow_bin ) {
      double var_overflow_min = bin_edges.back();

      std::stringstream true_ss;
      true_ss << signal_def << " && " << true_branchexpr
        << " >= " << var_overflow_min;

      std::string true_bin_def = true_ss.str();
      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_idx );

      std::stringstream reco_ss;
      reco_ss << selection << " && " << reco_branchexpr
        << " >= " << var_overflow_min;

      std::string reco_bin_def = reco_ss.str();

      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();
      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry. Note that the ROOT histogram bin
      // indices are one-based, so the overflow bin has an index equal to
      // the number of bin edges.
      cur_slice.bin_map_[ bin_edges.size() ].insert( ana_bin_idx );

      // Add the completed reco bin definition to the vector
      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_idx );
    }

  } // loop over kinematic variables

  // Add a single set of true bins for the background categories of interest
  for ( const auto& bdef : background_defs ) {
    true_bins.emplace_back( bdef, kBackgroundTrueBin, DUMMY_BLOCK_INDEX );
  }

  // Create a slice showing all blocks together as a function of bin number
  int num_reco_bins = reco_bins.size();
  int bin_number_var_idx = find_slice_var_index( "reco bin number",
    sb.slice_vars_ );

  auto& bin_num_slice = add_slice( sb, num_reco_bins, 0, num_reco_bins,
    bin_number_var_idx );
  for ( int ab = 0; ab < num_reco_bins; ++ab ) {
    // The ROOT histogram bins are one-based, so we correct for this here
    bin_num_slice.bin_map_[ ab + 1 ].insert( ab );
  }

  // Dump this information to the output file
  std::ofstream out_file( "myconfig_mcc8_all.txt" );
  out_file << "mcc8_all" << '\n';
  out_file << "stv_tree\n";
  out_file << true_bins.size() << '\n';
  for ( const auto& tb : true_bins ) out_file << tb << '\n';

  out_file << reco_bins.size() << '\n';
  for ( const auto& rb : reco_bins ) out_file << rb << '\n';

  // Also write a SliceBinning configuration file
  std::ofstream sb_file( "mybins_mcc8_all.txt" );
  sb_file << sb;
}
