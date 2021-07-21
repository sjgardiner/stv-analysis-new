#include "ResponseMatrixMaker.hh"

void make_config_mcc9_2D_muon() {

  // Using floating-point numbers as std::map keys is admittedly evil, but
  // it's safe in this case: all we'll do with this map is iterate over the
  // elements. Keys are muon momentum bin edges, values are muon scattering
  // cosine bin edges.
  std::map< double, std::vector<double> > muon_2D_bin_edges = {

    // No need for an underflow bin: due to the signal definition, all muons
    // with reco momentum below 0.1 GeV/c will be lost

    { 0.1,  { -1., 0., 1. } },

    { 0.17, { -1., -0.2, 0.4, 1. } },

    { 0.21, { -1, -0.2, 0.4, 1. } },

    { 0.24, { -1, -0.1, 0.5, 1. } },

    { 0.27, { -1, -0.1, 0.35, 0.6, 1. } },

    { 0.3,  { -1, -0.4, -0.1, 0.1, 0.35, 0.5, 0.7, 0.85, 1. } },

    { 0.38, { -1, 0, 0.5, 0.65, 0.8, 0.92, 1.00 } },

    { 0.48, { -1, 0.2, 0.5, 0.65, 0.8, 0.875, 0.950, 1.00 } },

    { 0.75, { -1, 0.5, 0.8, 0.875, 0.950, 1.00 } },

    { 1.14, { -1, 0.85, 0.9, 0.950, 1.00 } },

    // Upper edge of the last bin. The script will create an overflow bin above
    // this.
    { 2.5, {} }

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

  // Configure kinematic limits for all of the signal bins

  // Get an iterator to the last map element. They are sorted numerically,
  // so this will be the upper edge of the last non-overflow bin.
  auto last = muon_2D_bin_edges.cend();
  --last;

  auto iter = muon_2D_bin_edges.cbegin();
  for ( auto iter = muon_2D_bin_edges.cbegin(); iter != last; ++iter ) {

    // Get an iterator to the map element after the current one. Due to
    // the automatic sorting, this is guaranteed to contain the upper edge
    // of the current muon momentum bin.
    auto next = iter;
    ++next;

    // Get the current muon momentum bin limits
    double pmu_low = iter->first;
    double pmu_high = next->first;

    // Now iterate over the scattering cosine bins associated with the
    // current momentum bin. Note that we will skip any situations in
    // which the binning is undefined (i.e., because there are less than
    // two bin edges given)
    const auto& cosine_bin_edges = iter->second;

    size_t num_cosine_edges = cosine_bin_edges.size();
    size_t num_cosine_bins = 0u;
    if ( num_cosine_edges >= 2u ) num_cosine_bins = num_cosine_edges - 1u;

    for ( size_t b = 0u; b < num_cosine_bins; ++b ) {

      double cosmu_low = cosine_bin_edges.at( b );
      double cosmu_high = cosine_bin_edges.at( b + 1u );

      std::stringstream true_ss;
      true_ss << signal_def
        << " && mc_p3_mu.Mag() >= " << pmu_low
        << " && mc_p3_mu.Mag() < " << pmu_high
        << " && mc_p3_mu.CosTheta() >= " << cosmu_low
        << " && mc_p3_mu.CosTheta() < " << cosmu_high;

      std::string true_bin_def = true_ss.str();

      true_bins.emplace_back( true_bin_def, kSignalTrueBin );

      std::stringstream reco_ss;
      reco_ss << selection
        << " && p3_mu.Mag() >= " << pmu_low
        << " && p3_mu.Mag() < " << pmu_high
        << " && p3_mu.CosTheta() >= " << cosmu_low
        << " && p3_mu.CosTheta() < " << cosmu_high;

      std::string reco_bin_def = reco_ss.str();

      reco_bins.emplace_back( reco_bin_def );

      // We don't need an overflow cosine bin because the entire angular
      // range is covered. We'll use a single bin for the overflow in pmu.

    } // loop over scattering cosine bins

  } // loop over muon momentum bins


  // Create the single overflow bin for muon momentum in both true and reco
  // space

  double pmu_overflow_min = last->first;

  std::stringstream true_ss;
  true_ss << signal_def << " && mc_p3_mu.Mag() >= " << pmu_overflow_min;

  std::string true_bin_def = true_ss.str();
  true_bins.emplace_back( true_bin_def, kSignalTrueBin );

  std::stringstream reco_ss;
  reco_ss << selection << " && p3_mu.Mag() >= " << pmu_overflow_min;

  std::string reco_bin_def = reco_ss.str();
  reco_bins.emplace_back( reco_bin_def );

  // Add true bins for the background categories of interest
  for ( const auto& bdef : background_defs ) {
    true_bins.emplace_back( bdef, kBackgroundTrueBin );
  }

  // Dump this information to the output file
  std::ofstream out_file( "myconfig_mcc9_2D_muon.txt" );
  out_file << "Muon2D\n";
  out_file << "stv_tree\n";
  out_file << true_bins.size() << '\n';
  for ( const auto& tb : true_bins ) out_file << tb << '\n';

  out_file << reco_bins.size() << '\n';
  for ( const auto& rb : reco_bins ) out_file << rb << '\n';

  // Also write a SliceBinning configuration file
  std::ofstream sb_file( "mybins_mcc9_2D_muon.txt" );
  sb_file << "3\n";
  sb_file << "\"reco p_{#mu}\" \"GeV/c\" \"reco $p_{\\mu}$\" \"GeV$/c$\"\n";
  sb_file << "\"reco cos#theta_{#mu}\" \"\" \"reco $\\cos\\theta_{\\mu}$\""
    " \"\"\n";
  sb_file << "\"reco bin number\" \"\" \"reco bin number\" \"\"\n";
  // Includes a slice for the overflow bin and an extra slice for everything
  // in terms of reco bin number
  size_t num_slices = muon_2D_bin_edges.size() + 1;
  sb_file << num_slices << '\n';

  // Get an iterator to the final entry in the edge map (this is the
  // upper edge of the last bin)
  auto last_edge = muon_2D_bin_edges.cend();
  --last_edge;

  // The reco bins are numbered in the order that their edges appear in the
  // map, so just keep a running counter here to keep track of which reco
  // bin we are on.
  size_t cur_reco_bin = 0u;

  for ( auto iter = muon_2D_bin_edges.cbegin(); iter != last_edge; ++iter ) {
    // Each 1D slice uses the same y-axis units (reco events)
    sb_file << "\"events\"\n";
    const auto& edges = iter->second;
    int num_edges = edges.size();
    int num_bins = num_edges - 1;
    // The muon cosine is the sole "active variable" in each slice
    sb_file << "1 1 " << num_edges;
    for ( const auto& edge : edges ) {
      sb_file << ' ' << edge;
    }
    // The muon momentum is the sole "other variable" in each slice
    double pmu_low = iter->first;
    auto next = iter;
    ++next;
    double pmu_high = next->first;
    sb_file << "\n1 0 " << pmu_low << ' ' << pmu_high << '\n';
    sb_file << num_bins;

    for ( int b = 0; b < num_bins; ++b ) {
      int root_bin_idx = b + 1;
      sb_file << '\n' << cur_reco_bin << " 1 " << root_bin_idx;
      ++cur_reco_bin;
    } // cosine bins
    sb_file << '\n';

  } // pmu slices

  // Handle the overflow bin separately
  sb_file << "\"events\"\n"; // y-axis label
  // Still treat the muon cosine as the single active variable in a single bin
  // spanning the entire range
  sb_file << "1 1 2 -1.0 1.0\n";
  // The muon momentum is the sole "other" variable. This is an overflow bin,
  // which we signal by making the lower and upper edges equal
  double last_pmu_edge_value = last_edge->first;
  sb_file << "1 0 " << last_pmu_edge_value
    << ' ' << last_pmu_edge_value << '\n';
  // A single ResponseMatrixMaker reco bin contributes to the sole ROOT bin
  // in this histogram
  sb_file << "1\n" << cur_reco_bin << " 1 1";

  // Make a final slice with everything expressed in terms of reco bin number
  sb_file << "\"events\"\n"; // y-axis label
  sb_file << "1 2 ";
  // Acount for the zero-based ResponseMatrixMaker bin indices
  size_t num_reco_bins = cur_reco_bin + 1;
  // There is one more edge than the number of bins
  sb_file << num_reco_bins + 1;
  for ( size_t e = 0u; e <= num_reco_bins; ++e ) {
    sb_file << ' ' << e;
  }
  sb_file << '\n';
  // For the "overall slice," there is no other variable apart from reco bin
  // number
  sb_file << "0\n";
  // Loop over each ResponseMatrixMaker bin and assign it to the matching
  // ROOT histogram bin
  sb_file << num_reco_bins << '\n';
  for ( size_t b = 0u; b < num_reco_bins; ++b ) {
    sb_file << b << " 1 " << b + 1 << '\n';
  }
}
