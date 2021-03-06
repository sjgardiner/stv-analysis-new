#include "ResponseMatrixMaker.hh"

void make_config_mcc9_2D_proton() {

  // Using floating-point numbers as std::map keys is admittedly evil, but it's
  // safe in this case: all we'll do with this map is iterate over the
  // elements. Keys are proton momentum bin edges, values are proton cosine bin
  // edges.
  std::map< double, std::vector<double> > proton_2D_bin_edges = {

    // No need for an underflow bin: due to the signal definition, all leading
    // protons with reco momentum below 0.25 GeV/c will be lost
    { 0.250, { -1, -0.5, 0.1, 0.6, 1.0 } },
    { 0.325, { -1, -0.7, -0.4, 0, 0.4, 0.6, 0.8, 1.0 } },
    { 0.4,   { -1, -0.6, -0.2, 0.2, 0.5, 0.65, 0.85, 1.0 } },
    { 0.45,  { -1, -0.4, 0, 0.2, 0.4, 0.55, 0.65, 0.8, 0.92, 1.0 } },
    { 0.5,   { -1, -0.4, 0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 } },
    { 0.550, { -1, 0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 } },
    { 0.6,   { -1, 0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 } },
    { 0.65,  { -1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 } },
    { 0.7,   { -1, 0.4, 0.6, 0.75, 0.82, 0.9, 1.0 } },
    { 0.75,  { -1, 0.4, 0.6, 0.75, 0.87, 1.0 } },
    { 0.8,   { -1, 0.6, 0.73, 0.86, 1.0 } },
    { 0.85,  { -1, 0.6, 0.73, 0.86, 1.0 } },
    { 0.9,   { -1, 0.7, 0.85, 1.0 } },
    { 0.975, { -1, 0.8, 1.0 } },

    // Upper edge of the last bin. We don't need an overflow bin because the
    // signal definition excludes any leading protons with momenta above 1.2
    // GeV/c
    { 1.2, {} }

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
  auto last = proton_2D_bin_edges.cend();
  --last;

  auto iter = proton_2D_bin_edges.cbegin();
  for ( auto iter = proton_2D_bin_edges.cbegin(); iter != last; ++iter ) {

    // Get an iterator to the map element after the current one. Due to
    // the automatic sorting, this is guaranteed to contain the upper edge
    // of the current proton momentum bin.
    auto next = iter;
    ++next;

    // Get the current proton momentum bin limits
    double p_lead_p_low = iter->first;
    double p_lead_p_high = next->first;

    // Now iterate over the scattering cosine bins associated with the
    // current momentum bin. Note that we will skip any situations in
    // which the binning is undefined (i.e., because there are less than
    // two bin edges given)
    const auto& cosine_bin_edges = iter->second;

    size_t num_cosine_edges = cosine_bin_edges.size();
    size_t num_cosine_bins = 0u;
    if ( num_cosine_edges >= 2u ) num_cosine_bins = num_cosine_edges - 1u;

    for ( size_t b = 0u; b < num_cosine_bins; ++b ) {

      double cos_lead_p_low = cosine_bin_edges.at( b );
      double cos_lead_p_high = cosine_bin_edges.at( b + 1u );

      std::stringstream true_ss;
      true_ss << signal_def
        << " && mc_p3_lead_p.Mag() >= " << p_lead_p_low
        << " && mc_p3_lead_p.Mag() < " << p_lead_p_high
        << " && mc_p3_lead_p.CosTheta() >= " << cos_lead_p_low
        << " && mc_p3_lead_p.CosTheta() < " << cos_lead_p_high;

      std::string true_bin_def = true_ss.str();

      true_bins.emplace_back( true_bin_def, kSignalTrueBin );

      std::stringstream reco_ss;
      reco_ss << selection
        << " && p3_lead_p.Mag() >= " << p_lead_p_low
        << " && p3_lead_p.Mag() < " << p_lead_p_high
        << " && p3_lead_p.CosTheta() >= " << cos_lead_p_low
        << " && p3_lead_p.CosTheta() < " << cos_lead_p_high;

      std::string reco_bin_def = reco_ss.str();

      reco_bins.emplace_back( reco_bin_def );

      // We don't need an underflow or overflow cosine bin because the entire
      // angular range is covered in all cases.

    } // loop over cosine bins

  } // loop over leading proton momentum bins


  // Add true bins for the background categories of interest
  for ( const auto& bdef : background_defs ) {
    true_bins.emplace_back( bdef, kBackgroundTrueBin );
  }

  // Dump this information to the output file
  std::ofstream out_file( "myconfig_mcc9_2D_proton.txt" );
  out_file << "Proton2D\n";
  out_file << "stv_tree\n";
  out_file << true_bins.size() << '\n';
  for ( const auto& tb : true_bins ) out_file << tb << '\n';

  out_file << reco_bins.size() << '\n';
  for ( const auto& rb : reco_bins ) out_file << rb << '\n';

  // Also write a SliceBinning configuration file
  std::ofstream sb_file( "mybins_mcc9_2D_proton.txt" );
  sb_file << "3\n";
  sb_file << "\"reco p_{p}\" \"GeV/c\" \"reco $p_{p}$\" \"GeV$/c$\"\n";
  sb_file << "\"reco cos#theta_{p}\" \"\" \"$\\cos\\theta_{p}$\" \"\"\n";
  sb_file << "\"reco bin number\" \"\" \"reco bin number\" \"\"\n";
  // Include an extra slice which shows everything in terms of reco bin number
  size_t num_slices = proton_2D_bin_edges.size();
  sb_file << num_slices << '\n';

  // Get an iterator to the final entry in the edge map (this is the
  // upper edge of the last bin)
  auto last_edge = proton_2D_bin_edges.cend();
  --last_edge;

  // The reco bins are numbered in the order that their edges appear in the
  // map, so just keep a running counter here to keep track of which reco
  // bin we are on.
  size_t cur_reco_bin = 0u;

  for ( auto iter = proton_2D_bin_edges.cbegin(); iter != last_edge; ++iter ) {
    // Each 1D slice uses the same y-axis units (reco events)
    sb_file << "\"events\"\n";
    const auto& edges = iter->second;
    int num_edges = edges.size();
    int num_bins = num_edges - 1;
    // The proton cosine is the sole "active variable" in each slice
    sb_file << "1 1 " << num_edges;
    for ( const auto& edge : edges ) {
      sb_file << ' ' << edge;
    }
    // The proton momentum is the sole "other variable" in each slice
    double pp_low = iter->first;
    auto next = iter;
    ++next;
    double pp_high = next->first;
    sb_file << "\n1 0 " << pp_low << ' ' << pp_high << '\n';
    sb_file << num_bins;

    for ( int b = 0; b < num_bins; ++b ) {
      int root_bin_idx = b + 1;
      sb_file << '\n' << cur_reco_bin << " 1 " << root_bin_idx;
      ++cur_reco_bin;
    } // cosine bins
    sb_file << '\n';

  } // pp slices

  // Make a final slice with everything expressed in terms of reco bin number
  sb_file << "\"events\"\n"; // y-axis label
  sb_file << "1 2 ";
  size_t num_reco_bins = cur_reco_bin;
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
