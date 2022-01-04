#include "ResponseMatrixMaker.hh"

// Simple container for holding the configuration information for a single
// kinematic variable
struct VarConfig {

  VarConfig( const std::string& name, const std::vector<double>& edges,
    bool make_under, bool make_over ) : var_name_( name ), bin_edges_( edges ),
    needs_underflow_bin_( make_under ), needs_overflow_bin_( make_over ) {}

  std::string var_name_;
  std::vector<double> bin_edges_;
  bool needs_underflow_bin_;
  bool needs_overflow_bin_;
};

// Keys are the reco STV ntuple branch names of interest. Values are VarConfig
// objects containing the ResponseMatrixMaker configuration settings.
std::map< std::string, VarConfig > mcc9_1D_configs {

  { "p3_mu.CosTheta()", { "cth_mu",
    { -1, -0.925, -0.85, -0.775, -0.7, -0.625, -0.55, -0.475, -0.4, -0.325,
     -0.25, -0.175, -0.1, -0.025, 0.05, 0.125, 0.2, 0.275, 0.35, 0.425, 0.5,
      0.575, 0.65, 0.725, 0.8, 0.85, 0.875, 0.9, 0.925, 0.950, 0.975, 1.00 },
    false, false } },//unmodified

  { "p3_lead_p.Mag()", { "p_p",
    { 0.20, 0.305, 0.365, 0.415, 0.47, 0.525, 0.59, 0.68, 0.79, 0.93, 1.20 },
    false, false } },//modified, down to 10 bins from 14

  { "pn", { "pn", { 0., 0.125, 0.225, 0.325, 0.425, 0.525, 0.65, 0.85 }, false, true } },
  //what's pn?

  { "delta_alphaT", { "delta_alphaT",
    { 0, 0.6, 1.3, 2.0, 2.7, M_PI }, false, false } },
    //modified, down to 5 from 8 bins

    //I haven't seen these broken down in x, y, z components before - how to address?
    //unmodified for now
    //also, what's delta_pL?

  { "delta_pTx", { "delta_pTx",
    { -0.6, -0.45, -0.35, -0.25, -0.15, -0.075, 0, 0.075, 0.15, 0.25, 0.35, 0.45, 0.6 },
    true, true } },

  { "delta_pTy", { "delta_pTy",
    { -0.8, -0.55, -0.39, -0.2125, -0.05, 0.1, 0.225, 0.3375, 0.5 },
    true, true } },

  { "delta_pL", { "delta_pL",
    { -0.8, -0.6, -0.475, -0.35, -0.225, -0.115, -0.0285, 0.0575, 0.145, 0.230,
       0.315, 0.4 },
    true, true } },

  { "delta_phiT", { "delta_phiT",
    { 0., 0.1, 0.3, 0.6, 0.9, 1.3, 1.6, 2.0, 2.7, M_PI },
    false, false } },//updated, 9 bins from 13

  { "delta_pT", { "delta_pT", { 0., 0.1, 0.2, 0.3, 0.4, 0.525, 0.675, 0.9 },
    false, true } },//in progress

  { "p3_lead_p.CosTheta()", { "cth_p",
    { -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8,
      0.85, 0.9, 0.95, 1.0 }, false, false } },

  { "p3_mu.Mag()", { "p_mu", { 0.1, 0.17, 0.24, 0.3, 0.48, 0.75, 1.14, 2.5 },
    false, true } },

  { "thMuP", { "thMuP",
    { 0., 0.52, 0.78, 1.0, 1.15, 1.35, 1.5, 1.65, 1.8, 1.95, 2.1, 2.35, 2.62, M_PI },
    false, false } },

  // Opening angle
  // "TMath::ACos( (p3_mu.X()*p3_lead_p.X() + "
  // "p3_mu.Y()*p3_lead_p.Y() + p3_mu.Z()*p3_lead_p.Z()) / p3_mu.Mag()"
  //"/ p3_lead_p.Mag() )"

};


void make_config_mcc9() {

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
  for ( const auto& pair : mcc9_1D_configs ) {

    std::vector< TrueBin > true_bins;
    std::vector< RecoBin > reco_bins;

    std::string reco_branchexpr = pair.first;
    std::string true_branchexpr = "mc_" + reco_branchexpr;

    const VarConfig& var_config = pair.second;

    bool needs_underflow_bin = var_config.needs_underflow_bin_;
    bool needs_overflow_bin = var_config.needs_overflow_bin_;

    // Require at least two bin edges to be present in the input vector.
    // Any variables for which this is not true will be skipped entirely.
    const auto& bin_edges = var_config.bin_edges_;

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

    std::string var_name = var_config.var_name_;

    // Dump this information to the output file
    std::ofstream out_file( "myconfig_mcc9_" + var_name + ".txt" );
    out_file << "mcc9_" << var_name << '\n';
    out_file << "stv_tree\n";
    out_file << true_bins.size() << '\n';
    for ( const auto& tb : true_bins ) out_file << tb << '\n';

    out_file << reco_bins.size() << '\n';
    for ( const auto& rb : reco_bins ) out_file << rb << '\n';

  } // loop over kinematic variables

}
