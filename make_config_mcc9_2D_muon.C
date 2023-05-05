#include "ConfigMakerUtils.hh"
#include "UniverseMaker.hh"

void make_config_mcc9_2D_muon() {

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

  // Also save the bin definitions to an output file that will be used to
  // make a LaTeX table
  std::ofstream tex_bin_table_file( "mybintable_mcc9_muon2D.tex" );
  tex_bin_table_file << "\\documentclass{standalone}\n"
    << "\\usepackage{booktabs}\n"
    << "\\usepackage{siunitx}\n"
    << "\\DeclareSIUnit\\clight{\\text{\\ensuremath{c}}}\n"
    << "\\sisetup{per-mode=symbol}\n"
    << "\\begin{document}\n"
    << "\\begin{tabular}{cSSSScc}\n"
    << "\\toprule\n"
    << "bin number\n"
    << "& {$p_\\mu^\\mathrm{low}$ (\\si{\\GeV\\per\\clight})}"
    << "& {$p_\\mu^\\mathrm{high}$ (\\si{\\GeV\\per\\clight})}"
    << " & {$\\cos\\theta_\\mu^\\mathrm{low}$}\n"
    << "& {$\\cos\\theta_\\mu^\\mathrm{high}$} & efficiency & occupancy \\\\\n"
    << "\\midrule\n";

  // Configure kinematic limits for all of the signal bins

  // The reco bins are numbered in the order that their edges appear in the
  // map, so just keep a running counter here to keep track of which reco
  // bin we are on.
  size_t cur_reco_bin = 0u;

  // Get an iterator to the last map element. They are sorted numerically,
  // so this will be the upper edge of the last non-overflow bin.
  auto last = MUON_2D_BIN_EDGES.cend();
  --last;

  auto iter = MUON_2D_BIN_EDGES.cbegin();
  for ( auto iter = MUON_2D_BIN_EDGES.cbegin(); iter != last; ++iter ) {

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

      true_bins.emplace_back( true_bin_def, kSignalTrueBin, 0 );

      std::stringstream reco_ss;
      reco_ss << selection
        << " && p3_mu.Mag() >= " << pmu_low
        << " && p3_mu.Mag() < " << pmu_high
        << " && p3_mu.CosTheta() >= " << cosmu_low
        << " && p3_mu.CosTheta() < " << cosmu_high;

      std::string reco_bin_def = reco_ss.str();

      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, 0 );

      tex_bin_table_file << cur_reco_bin << " & ";
      ++cur_reco_bin;
      if ( b == 0u ) {
        tex_bin_table_file << pmu_low << " & " << pmu_high << " & ";
      }
      else tex_bin_table_file << " & & ";

      tex_bin_table_file << cosmu_low << " & " << cosmu_high << " &  &"
        << " \\\\";
      // Add extra space at the end of each momentum bin
      if ( b == num_cosine_bins - 1 ) tex_bin_table_file << "[2mm]";
      tex_bin_table_file << '\n';

      // We don't need an overflow cosine bin because the entire angular
      // range is covered. We'll use a single bin for the overflow in pmu.

    } // loop over scattering cosine bins

  } // loop over muon momentum bins

  // No overflow bin is needed due to the muon momentum upper limit in the
  // signal definition

  // Add true bins for the background categories of interest
  for ( const auto& bdef : background_defs ) {
    true_bins.emplace_back( bdef, kBackgroundTrueBin );
  }

  // We're done with all the "ordinary" bins. Add some extra reco bins to use
  // in the sideband control samples.

  // The control samples have significantly lower statistics, so bin in 1D
  // momentum slices only. Use the muon candidate momentum where possible, but
  // for NC switch to the leading proton candidate momentum. This is done
  // because the NC sideband selection excludes events in which a muon candidate
  // has been identified.
  //
  // Keys are selections to use for sidebands, values are the branch names for
  // the reconstructed momentum to use in each case. This is a pretty hacky way
  // to organize the information, but it is simple.
  std::map< std::string, std::string > sideband_selection_to_momentum_map = {
    {  DIRT_SIDEBAND_SELECTION, "p3_mu"     },
    {    NC_SIDEBAND_SELECTION, "p3_lead_p" },
    { CCNPI_SIDEBAND_SELECTION, "p3_mu"     },
  };

  // Loop over the sideband selection definitions. Prepare new reco bin
  // definitions for each in the appropriate 1D reconstructed momentum space.
  for ( const auto& sel_mom_pair : sideband_selection_to_momentum_map ) {
    const auto& side_sel = sel_mom_pair.first;
    const auto& mom_branch = sel_mom_pair.second;
    std::map< double, std::vector<double> >* bin_edge_map = nullptr;
    if ( mom_branch == "p3_mu" ) bin_edge_map = &MUON_2D_BIN_EDGES;
    else if ( mom_branch == "p3_lead_p" ) bin_edge_map = &PROTON_2D_BIN_EDGES;
    else throw std::runtime_error( "Unimplemented sideband momentum!" );

    // Get an iterator to the last map element. They are sorted numerically,
    // so this will be the upper edge of the last non-overflow momentum bin.
    auto last = bin_edge_map->cend();
    --last;

    for ( auto iter = bin_edge_map->cbegin(); iter != last; ++iter ) {

      // Get an iterator to the map element after the current one. Due to the
      // automatic sorting, this is guaranteed to contain the upper edge of the
      // current momentum bin.
      auto next = iter;
      ++next;

      // Get the current momentum bin limits
      double p_low = iter->first;
      double p_high = next->first;

      std::stringstream reco_ss;
      reco_ss << side_sel
        << " && " << mom_branch << ".Mag() >= " << p_low
        << " && " << mom_branch << ".Mag() < " << p_high;

      std::string reco_bin_def = reco_ss.str();

      reco_bins.emplace_back( reco_bin_def, kSidebandRecoBin );
    } // 1D reco momentum bins

    // For sidebands that use bins of muon candidate momentum, create the
    // overflow bin. This isn't needed for the proton momentum due to the
    // upper limit imposed in the signal definition.
    if ( mom_branch != "p3_mu" ) continue;

    double pmu_overflow_min = bin_edge_map->crbegin()->first;

    std::stringstream reco_ss;
    reco_ss << side_sel << " && " << mom_branch << ".Mag() >= "
      << pmu_overflow_min;

    std::string reco_bin_def = reco_ss.str();
    reco_bins.emplace_back( reco_bin_def, kSidebandRecoBin );

  } // sideband selection definitions

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
  // Includes a slice for the overflow bin and three extra slices. One
  // for everything in terms of reco bin number, one integrated over angles,
  // and one showing the sideband control sample results (in terms of reco
  // bin number).
  size_t num_slices = MUON_2D_BIN_EDGES.size() + 3;
  sb_file << num_slices << '\n';

  // Get an iterator to the final entry in the edge map (this is the
  // upper edge of the last bin)
  auto last_edge = MUON_2D_BIN_EDGES.cend();
  --last_edge;

  // The reco bins are numbered in the order that their edges appear in the
  // map, so just keep a running counter here to keep track of which reco
  // bin we are on.
  cur_reco_bin = 0u;

  for ( auto iter = MUON_2D_BIN_EDGES.cbegin(); iter != last_edge; ++iter ) {
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
  // A single UniverseMaker reco bin contributes to the sole ROOT bin
  // in this histogram
  sb_file << "1\n" << cur_reco_bin << " 1 1";

  // Make a final slice with everything expressed in terms of reco bin number
  sb_file << "\"events\"\n"; // y-axis label
  sb_file << "1 2 ";
  // Acount for the zero-based UniverseMaker bin indices
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
  // Loop over each UniverseMaker bin and assign it to the matching
  // ROOT histogram bin
  sb_file << num_reco_bins << '\n';
  for ( size_t b = 0u; b < num_reco_bins; ++b ) {
    sb_file << b << " 1 " << b + 1 << '\n';
  }

  // Make a 1D slice in which the angles have been integrated out. This is a
  // measurement of the leading muon momentum distribution. For this slice,
  // we're still working in terms of reco event counts
  sb_file << "\"events\"\n";
  int num_pmu_edges = MUON_2D_BIN_EDGES.size();
  int num_pmu_bins = num_pmu_edges - 1;
  // The muon momentum is the sole "active variable" in each slice
  sb_file << "1 0 " << num_pmu_edges;
  for ( const auto& pmu_edge_pair : MUON_2D_BIN_EDGES ) {
    const auto pmu_edge = pmu_edge_pair.first;
    sb_file << ' ' << pmu_edge;
  }
  // There is no "other" variable for this slice since we've integrated out
  // the angular information
  sb_file << "\n0\n";
  // Now we're ready to build the 1D muon momentum bins from the 2D ones. We
  // need one entry in the list per reco bin (apart from the overflow bin),
  // although multiple reco bins will contribute to each slice bin in this
  // case.
  sb_file << num_reco_bins - 1;

  // Iterate through the 2D reco bins, noting that they are numbered in the
  // order that their edges appear in the map. In this case, all angular reco
  // bins with the same reco muon momentum should contribute to a particular
  // slice p_mu bin.
  cur_reco_bin = 0u;

  // Keep track of the ROOT slice bin index (one-based) with this counter
  int cur_slice_bin_idx = 1;

  for ( auto iter = MUON_2D_BIN_EDGES.cbegin(); iter != last_edge; ++iter ) {
    const auto& angle_bin_edges = iter->second;
    int num_angle_bins = angle_bin_edges.size() - 1;
    for ( int b = 0; b < num_angle_bins; ++b ) {
      sb_file << '\n' << cur_reco_bin << " 1 " << cur_slice_bin_idx;
      ++cur_reco_bin;
    } // cosine bins

    // Move to the next muon momentum bin in the slice
    ++cur_slice_bin_idx;
  } // pmu slices

  sb_file << '\n';

  // Make a slice containing the sideband results organized by reco bin number
  sb_file << "\"events\"\n"; // y-axis label
  sb_file << "1 2 ";
  // Count the number of sideband reco bins. Also find the index of the
  // first one.
  size_t num_sideband_reco_bins = 0u;
  size_t first_sideband_bin_idx = 0u;
  bool found_first_sideband_bin = false;
  for ( size_t b = 0u; b < reco_bins.size(); ++b ) {
    const auto& rbin = reco_bins.at( b );
    if ( rbin.type_ == kSidebandRecoBin ) {
      ++num_sideband_reco_bins;
      if ( !found_first_sideband_bin ) {
        found_first_sideband_bin = true;
        first_sideband_bin_idx = b;
      }
    }
  }
  // There is one more edge than the number of sideband bins
  sb_file << num_sideband_reco_bins + 1;
  for ( size_t e = 0u; e <= num_sideband_reco_bins; ++e ) {
    sb_file << ' ' << e + first_sideband_bin_idx;
  }
  sb_file << '\n';
  // For the "sideband slice," there is no other variable apart from reco bin
  // number
  sb_file << "0\n";
  // Loop over each UniverseMaker sideband bin and assign it to the
  // matching ROOT histogram bin
  sb_file << num_sideband_reco_bins << '\n';
  for ( size_t b = 0u; b < num_sideband_reco_bins; ++b ) {
    sb_file << b + first_sideband_bin_idx << " 1 " << b + 1 << '\n';
  }

}
