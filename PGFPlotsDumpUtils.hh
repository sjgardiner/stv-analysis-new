#pragma once

// STV analysis includes
#include "SliceBinning.hh"
#include "SliceHistogram.hh"

constexpr double REALLY_BIG_NUMBER = 1e300;

// Utility functions used to prepare input files for use with the LaTeX package
// PGFPlots

const TAxis* get_axis_by_index( int d, const TH1& hist ) {
  const TAxis* axis = nullptr;
  switch ( d ) {
    case 0:
      axis = hist.GetXaxis();
      break;
    case 1:
      axis = hist.GetYaxis();
      break;
    case 2:
      axis = hist.GetZaxis();
      break;
    default:
      throw std::runtime_error( "Unrecognized TAxis index" );
  }
  return axis;
}

void dump_slice_variables( const SliceBinning& sb, size_t slice_idx,
  std::map< std::string, std::string >& pgfplots_params_table )
{
  const auto& slice = sb.slices_.at( slice_idx );
  const auto& vars = sb.slice_vars_;

  int num_active_vars = slice.active_var_indices_.size();
  pgfplots_params_table[ "num_active_vars" ]
    = std::to_string( num_active_vars );

  for ( size_t a = 0u; a < num_active_vars; ++a ) {
    size_t av_idx = slice.active_var_indices_.at( a );
    const SliceVariable& svar = vars.at( av_idx );

    std::string name_for_dump = "{$" + svar.latex_name_ + "$}";
    std::string units_for_dump = "{$" + svar.latex_units_ + "$}";

    std::string col_prefix = "av" + std::to_string( a );
    pgfplots_params_table[ col_prefix + "_name" ] = name_for_dump;
    pgfplots_params_table[ col_prefix + "_units" ] = units_for_dump;
  }


  int num_other_vars = slice.other_vars_.size();
  pgfplots_params_table[ "num_other_vars" ] = std::to_string( num_other_vars );
  for ( size_t b = 0u; b < num_other_vars; ++b ) {
    const auto& ov_spec = slice.other_vars_.at( b );
    size_t ov_idx = ov_spec.var_index_;
    const SliceVariable& svar = vars.at( ov_idx );

    std::string name_for_dump = "{$" + svar.latex_name_ + "$}";
    std::string units_for_dump = "{$" + svar.latex_units_ + "$}";

    double low = ov_spec.low_bin_edge_;
    double high = ov_spec.high_bin_edge_;

    std::string col_prefix = "ov" + std::to_string( b );
    pgfplots_params_table[ col_prefix + "_name" ] = name_for_dump;
    pgfplots_params_table[ col_prefix + "_units" ] = units_for_dump;
    pgfplots_params_table[ col_prefix + "_low" ]
      = low < -REALLY_BIG_NUMBER ? "{inf}" : std::to_string( low );
    pgfplots_params_table[ col_prefix + "_high" ]
      = high > REALLY_BIG_NUMBER ? "{inf}" : std::to_string( high );
  }
}

void write_pgfplots_files( const std::string& out_filename_prefix,
  std::map< std::string, std::vector<double> >& pgfplots_hist_table,
  std::map< std::string, std::string >& pgfplots_params_table )
{
  std::ofstream hist_out_file( out_filename_prefix + "_hist.txt" );

  for ( const auto& pair : pgfplots_hist_table ) {
    const std::string& col_name = pair.first;
    hist_out_file << "  " << col_name;
  }
  hist_out_file << '\n';

  size_t num_rows = pgfplots_hist_table.cbegin()->second.size();
  for ( size_t r = 0u; r < num_rows; ++r ) {
    for ( const auto& pair : pgfplots_hist_table ) {
      hist_out_file << "  " << pair.second.at( r );
    }
    hist_out_file << '\n';
  }

  std::ofstream params_out_file( out_filename_prefix + "_params.txt" );

  for ( const auto& pair : pgfplots_params_table ) {
    const std::string& col_name = pair.first;
    params_out_file << "  " << col_name;
  }
  params_out_file << '\n';
  for ( const auto& pair : pgfplots_params_table ) {
    const auto& param_value = pair.second;
    params_out_file << "  " << param_value;
  }
}

void dump_bin_edges_and_half_widths( const std::string& hist_col_prefix,
  const TH1& hist, int global_bin_idx,
  std::map< std::string, std::vector<double> >& pgf_plots_hist_table )
{
  std::array< int, 3 > temp_axis_bin_indices = { 0, 0, 0 };
  hist.GetBinXYZ( global_bin_idx, temp_axis_bin_indices.at(0),
    temp_axis_bin_indices.at(1), temp_axis_bin_indices.at(2) );

  int dimension = hist.GetDimension();
  for ( int d = 0; d < dimension; ++d ) {
    const TAxis* axis = get_axis_by_index( d, hist );
    int axis_bin_idx = temp_axis_bin_indices.at( d );
    double low_edge = axis->GetBinLowEdge( axis_bin_idx );
    double half_width = 0.5 * axis->GetBinWidth( axis_bin_idx );

    std::string x_col_name = hist_col_prefix;
    if ( !hist_col_prefix.empty() ) x_col_name += '_';
    x_col_name += 'x' + std::to_string( d );
    auto end = pgf_plots_hist_table.end();
    auto iter = pgf_plots_hist_table.find( x_col_name );
    if ( iter == end ) {
      pgf_plots_hist_table[ x_col_name ] = std::vector<double> { low_edge };
    }
    else iter->second.push_back( low_edge );

    std::string x_hw_col_name = x_col_name + "_halfwidth";
    auto end2 = pgf_plots_hist_table.end();
    auto iter2 = pgf_plots_hist_table.find( x_hw_col_name );
    if ( iter2 == end2 ) {
      pgf_plots_hist_table[ x_hw_col_name ]
        = std::vector<double> { half_width };
    }
    else iter2->second.push_back( half_width );
  } // dimension d
}

void dump_slice_plot_limits( const SliceHistogram& slice_bnb,
  const SliceHistogram& slice_mc_plus_ext, const Slice& slice,
  std::map< std::string, std::string >& pgfplots_params_table )
{
  // Get the minimum and maximum coordinates along each of the bin axes
  int dimension = slice_bnb.hist_->GetDimension();
  for ( int d = 0; d < dimension; ++d ) {
    const TAxis* axis = get_axis_by_index( d, *slice_bnb.hist_ );
    // Note the one-based bin indices used by ROOT histograms
    double min = axis->GetBinLowEdge( 1 );
    double max = axis->GetBinLowEdge( axis->GetNbins() + 1 );
    std::string col_prefix = "x" + std::to_string( d );
    pgfplots_params_table[ col_prefix + "_min" ] = std::to_string( min );
    pgfplots_params_table[ col_prefix + "_max" ] = std::to_string( max );
  }

  // TODO: adjust to use a variable y_min if needed
  pgfplots_params_table[ "y_min" ] = std::to_string( 0. );

  // Find a y_max value that allows all bins in the slice to be fully seen in a
  // plot.
  double y_max = std::numeric_limits<double>::lowest();
  std::array< const SliceHistogram*, 2 > temp_slice_hists = { &slice_bnb,
    &slice_mc_plus_ext };

  for ( const auto* sh : temp_slice_hists ) {
    for ( const auto& bin_pair : slice.bin_map_ ) {
      int global_bin_idx = bin_pair.first;
      double y = sh->hist_->GetBinContent( global_bin_idx );
      double yerror = sh->hist_->GetBinError( global_bin_idx );
      double y_up = y + yerror;
      if ( y_up > y_max ) y_max = y_up;
    }
  }

  // Use a margin of an extra 3% on the maximum detected y value
  y_max *= 1.03;
  pgfplots_params_table[ "y_max" ] = std::to_string( y_max );

  // Also find the maximum deviation from unity to use when setting the y-axis
  // range for a data/MC+EXT ratio plot.
  double ratio_max = std::numeric_limits<double>::lowest();
  for ( const auto& bin_pair : slice.bin_map_ ) {
    int global_bin_idx = bin_pair.first;

    double yBNB = slice_bnb.hist_->GetBinContent( global_bin_idx );
    double yBNBerror = slice_bnb.hist_->GetBinError( global_bin_idx );

    double yMC = slice_mc_plus_ext.hist_->GetBinContent( global_bin_idx );
    double yMCerror = slice_mc_plus_ext.hist_->GetBinError( global_bin_idx );

    // Consider several different calculations of the ratio that will appear
    // in the plot. Pick the one that differs the most from unity.
    std::array< double, 4 > abs_diffs;
    abs_diffs[0] = ( yBNB + yBNBerror ) / yMC;
    abs_diffs[1] = ( yBNB - yBNBerror ) / yMC;
    abs_diffs[2] = ( yMC + yMCerror ) / yMC;
    abs_diffs[3] = ( yMC - yMCerror ) / yMC;

    for ( auto& ad : abs_diffs ) {
      ad = std::abs( ad - 1 );
    }

    double max_abs_diff = *std::max_element( abs_diffs.cbegin(),
      abs_diffs.cend() );

    if ( ratio_max < max_abs_diff ) ratio_max = max_abs_diff;
  }

  // Use a margin of an extra 30% on the maximum detected absolute deviation
  // from unity
  ratio_max *= 1.30;
  pgfplots_params_table[ "ratio_max" ] = std::to_string( ratio_max );

}

void dump_slice_histogram( const std::string& hist_col_prefix,
  const SliceHistogram& slice_hist, const Slice& slice,
  std::map< std::string, std::vector<double> >& pgf_plots_hist_table,
  bool include_yerror = true, bool include_x_coords = false )
{
  // Write information about the input SliceHistogram to the input map of
  // PGFPlots table columns
  const auto* hist = slice_hist.hist_.get();
  // TODO: consider adding a check for pre-existing duplicate columns
  // (to avoid accidentally overwriting information)

  std::string bin_col_name = "bin";
  if ( include_x_coords ) {
    pgf_plots_hist_table[ bin_col_name ] = std::vector<double>();
  }

  std::string y_col_name = hist_col_prefix;
  pgf_plots_hist_table[ y_col_name ] = std::vector<double>();

  std::string yerror_col_name;
  if ( include_yerror ) {
    yerror_col_name = hist_col_prefix + "_error";
    pgf_plots_hist_table[ yerror_col_name ] = std::vector<double>();
  }

  std::string y_norm_error_col_name;
  std::string y_shape_error_col_name;
  std::string y_mixed_error_col_name;

  for ( const auto& bin_pair : slice.bin_map_ ) {
    int global_bin_idx = bin_pair.first;
    double y = hist->GetBinContent( global_bin_idx );

    pgf_plots_hist_table.at( y_col_name ).push_back( y );

    if ( include_yerror ) {
      double yerror = hist->GetBinError( global_bin_idx );
      pgf_plots_hist_table.at( yerror_col_name ).push_back( yerror );
    }

    if ( include_x_coords ) {
      pgf_plots_hist_table.at( bin_col_name ).push_back( global_bin_idx );

      dump_bin_edges_and_half_widths( "", *hist, global_bin_idx,
        pgf_plots_hist_table );
    }

  } // slice bins

  // Add a (presumably empty) overflow bin to get certain PGFPlots styles to
  // look right. Cast the TH1* to a TArray* so that we can call the GetSize()
  // member function. See
  // https://root.cern.ch/root/roottalk/roottalk01/0359.html for an
  // explanation.
  const TArray* arr_ptr = dynamic_cast< const TArray* >( hist );
  if ( !arr_ptr ) throw std::runtime_error( "Failed dynamic cast in"
    " dump_slice_histogram()" );

  int last_overflow_global_bin_idx = arr_ptr->GetSize() - 1;

  double overflow_y = hist->GetBinContent( last_overflow_global_bin_idx );
  pgf_plots_hist_table.at( y_col_name ).push_back( overflow_y );

  if ( include_yerror ) {
    double overflow_yerror = hist->GetBinError( last_overflow_global_bin_idx );
    pgf_plots_hist_table.at( yerror_col_name ).push_back( overflow_yerror );
  }

  if ( include_x_coords ) {
    pgf_plots_hist_table.at( bin_col_name ).push_back(
      last_overflow_global_bin_idx );

    dump_bin_edges_and_half_widths( "", *hist, last_overflow_global_bin_idx,
      pgf_plots_hist_table );
  }

}
