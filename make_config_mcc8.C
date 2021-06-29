#include "ResponseMatrixMaker.hh"

void make_config_mcc8() {

  // Keys are the reco STV ntuple branch names of interest. Values
  // are vectors of bin edges (taken from the MCC8 CCNp0pi paper)
  std::map< std::string, std::vector<double> > mcc8_bin_edge_map = {

    { "p3_mu.CosTheta()", { -1.0, -0.82, -0.66, -0.39, -0.16, 0.05, 0.25, 0.43,
      0.59, 0.73, 0.83, 0.91, 1.0 } },

    { "thMuP", { 0.0, 0.8, 1.2, 1.57, 1.94, 2.34, M_PI } },

    { "p3_lead_p.Mag()", { 0.3, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87,
      0.93, 1.2 } },

    { "p3_mu.Mag()", { 0.1, 0.18, 0.3, 0.48, 0.75, 1.14, 2.5 } },

    { "p3_lead_p.CosTheta()", { -1.0, -0.5, 0.0, 0.27, 0.45, 0.62, 0.76, 0.86,
      0.94, 1.0 } },

  // Opening angle
  // "TMath::ACos( (p3_mu.X()*p3_lead_p.X() + "
  // "p3_mu.Y()*p3_lead_p.Y() + p3_mu.Z()*p3_lead_p.Z()) / p3_mu.Mag()"
  //"/ p3_lead_p.Mag() )"

  };

  // Keys are the same reco variable branch expressions as above. Values
  // are bool pairs indicating whether an underflow and overflow bin
  // should be produced.
  std::map< std::string, std::pair<bool,bool> > mcc8_under_overflow_map = {

    // Restricted by valid angular ranges to lie within the bin limits
    { "p3_mu.CosTheta()", { false, false } },

    { "p3_lead_p.CosTheta()", { false, false } },

    { "thMuP", { false, false } },

    // Restricted by the signal definition to lie within the proton momentum
    // bins
    { "p3_lead_p.Mag()", { false, false } },

    // No underflow bin is needed due to the signal definition. No upper limit
    // is imposed, however, so we will make an overflow bin.
    { "p3_mu.Mag()", { false, true } }

  };

  // Variable nicknames to use when naming the output config files
  std::map< std::string, std::string > mcc8_var_name_map = {

    { "p3_mu.CosTheta()", "cth_mu" },

    { "p3_lead_p.CosTheta()", "cth_p" },

    { "thMuP", "th_mu_p" },

    { "p3_lead_p.Mag()", "p_p" },

    { "p3_mu.Mag()", "p_mu" }

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


  // Iterate over each kinematic variable and make a config file for it
  // based on the binning from the MCC8 CCNp0pi analysis
  for ( const auto& pair : mcc8_bin_edge_map ) {

    std::vector< TrueBin > true_bins;
    std::vector< RecoBin > reco_bins;

    std::string reco_branchexpr = pair.first;
    std::string true_branchexpr = "mc_" + reco_branchexpr;

    // I didn't yet make a separate branch for the opening angle between
    // the muon and proton, so fix the branch expressions for that case
    // with this ugly hack
    if ( reco_branchexpr == "thMuP" ) {
      reco_branchexpr = "TMath::ACos( (p3_mu.X()*p3_lead_p.X() + p3_mu.Y()*p3_lead_p.Y() + p3_mu.Z()*p3_lead_p.Z()) / p3_mu.Mag() / p3_lead_p.Mag() )";
      true_branchexpr = "TMath::ACos( (mc_p3_mu.X()*mc_p3_lead_p.X() + mc_p3_mu.Y()*mc_p3_lead_p.Y() + mc_p3_mu.Z()*mc_p3_lead_p.Z()) / mc_p3_mu.Mag() / mc_p3_lead_p.Mag() )";

      // Add to the necessary maps as well so we don't break the lookups below
      mcc8_var_name_map[ reco_branchexpr ] = mcc8_var_name_map.at( "thMuP" );
      mcc8_under_overflow_map[ reco_branchexpr ]
        = mcc8_under_overflow_map.at( "thMuP" );
    }

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

    // If needed, then create the underflow bin in both true and reco space
    if ( needs_underflow_bin ) {
      double var_underflow_max = bin_edges.front();

      std::stringstream true_ss;
      true_ss << signal_def << " && " << true_branchexpr
        << " < " << var_underflow_max;

      std::string true_bin_def = true_ss.str();
      true_bins.emplace_back( true_bin_def, kSignalTrueBin );

      std::stringstream reco_ss;
      reco_ss << selection << " && " << reco_branchexpr
        << " < " << var_underflow_max;

      std::string reco_bin_def = reco_ss.str();
      reco_bins.emplace_back( reco_bin_def );
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

      true_bins.emplace_back( true_bin_def, kSignalTrueBin );

      std::stringstream reco_ss;
      reco_ss << selection
        << " && " << reco_branchexpr << " >= " << var_low
        << " && " << reco_branchexpr << " < "  << var_high;

      std::string reco_bin_def = reco_ss.str();

      reco_bins.emplace_back( reco_bin_def );

    } // loop over ordinary bins for the current variable

    // If needed, then create the overflow bin in both true and reco space
    if ( needs_overflow_bin ) {
      double var_overflow_min = bin_edges.back();

      std::stringstream true_ss;
      true_ss << signal_def << " && " << true_branchexpr
        << " >= " << var_overflow_min;

      std::string true_bin_def = true_ss.str();
      true_bins.emplace_back( true_bin_def, kSignalTrueBin );

      std::stringstream reco_ss;
      reco_ss << selection << " && " << reco_branchexpr
        << " >= " << var_overflow_min;

      std::string reco_bin_def = reco_ss.str();
      reco_bins.emplace_back( reco_bin_def );
    }

    // Add true bins for the background categories of interest
    for ( const auto& bdef : background_defs ) {
      true_bins.emplace_back( bdef, kBackgroundTrueBin );
    }

    std::string var_name = mcc8_var_name_map.at( reco_branchexpr );

    // Dump this information to the output file
    std::ofstream out_file( "myconfig_mcc8_" + var_name + ".txt" );
    out_file << "mcc8_" << var_name << '\n';
    out_file << "stv_tree\n";
    out_file << true_bins.size() << '\n';
    for ( const auto& tb : true_bins ) out_file << tb << '\n';

    out_file << reco_bins.size() << '\n';
    for ( const auto& rb : reco_bins ) out_file << rb << '\n';

  } // loop over kinematic variables

}
