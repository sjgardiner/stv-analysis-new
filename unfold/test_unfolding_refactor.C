// Standard library includes
#include <iomanip>
#include <iostream>
#include <sstream>

// ROOT includes
#include "TCanvas.h"
#include "TLegend.h"

// STV analysis includes
#include "../CrossSectionExtractor.hh"
#include "../PGFPlotsDumpUtils.hh"
#include "../SliceBinning.hh"
#include "../SliceHistogram.hh"

constexpr double BIG_DOUBLE = 1e300;

void dump_slice_errors( const std::string& hist_col_prefix,
  const Slice& slice, const std::map< std::string,
  std::unique_ptr<SliceHistogram> >& slice_hist_cov_matrix_map,
  std::map< std::string, std::vector<double> >& pgf_plots_hist_table )
{
  for ( const auto& pair : slice_hist_cov_matrix_map ) {
    std::string err_name = pair.first;
    std::string err_col_name = hist_col_prefix + '_' + err_name + "_error";
    pgf_plots_hist_table[ err_col_name ] = std::vector<double>();
  }

  for ( const auto& bin_pair : slice.bin_map_ ) {
    // TODO: revisit for multi-dimensional slices
    int global_bin_idx = bin_pair.first;

    for ( const auto& err_pair : slice_hist_cov_matrix_map ) {

      std::string err_name = err_pair.first;
      std::string err_col_name = hist_col_prefix + '_' + err_name + "_error";

      const auto* hist = err_pair.second->hist_.get();
      double err = hist->GetBinError( global_bin_idx );

      pgf_plots_hist_table.at( err_col_name ).push_back( err );
    }

  } // slice bins

  // Add a (presumably empty) overflow bin to get certain PGFPlots styles to
  // look right.
  for ( const auto& err_pair : slice_hist_cov_matrix_map ) {
    std::string err_name = err_pair.first;
    std::string err_col_name = hist_col_prefix + '_' + err_name + "_error";

    pgf_plots_hist_table.at( err_col_name ).push_back( 0. );
  }

}

void test_unfolding_refactor() {

  // Use a CrossSectionExtractor object to handle the systematics and unfolding
  auto extr = std::make_unique< CrossSectionExtractor >( "../xsec_config.txt" );
  auto xsec = extr.get_unfolded_events();
  double conv_factor = extr.conversion_factor();

  // Free up the memory used by the CrossSectionExtractor now that we have the
  // results we need
  extr.reset();

  // Plot slices of the unfolded result
  auto* sb_ptr = new SliceBinning( "../mybins_all.txt" );
  auto& sb = *sb_ptr;

  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {

    const auto& slice = sb.slices_.at( sl_idx );

    // Make a histogram showing the unfolded true event counts in the current
    // slice
    SliceHistogram* slice_unf = SliceHistogram::make_slice_histogram(
      *xsec.result_.unfolded_signal_, slice, xsec.result_.cov_matrix_.get() );

    // Temporary copies of the unfolded true event count slices with
    // different covariance matrices
    std::map< std::string, std::unique_ptr<SliceHistogram> > sh_cov_map;
    for ( const auto& uc_pair : unfolded_cov_matrix_map ) {
      const auto& uc_name = uc_pair.first;
      const auto& uc_matrix = uc_pair.second;

      auto& uc_ptr = sh_cov_map[ uc_name ];
      uc_ptr.reset(
        SliceHistogram::make_slice_histogram( *xsec.result_.unfolded_signal_,
          slice, uc_matrix.get() )
      );
    }

    // Also use the GENIE CV model to do the same
    SliceHistogram* slice_cv = SliceHistogram::make_slice_histogram(
      genie_cv_truth_vec, slice, nullptr );

    // If present, also use the truth information from the fake data to do the
    // same
    SliceHistogram* slice_truth = nullptr;
    if ( using_fake_data ) {
      slice_truth = SliceHistogram::make_slice_histogram( fake_data_truth,
        slice, nullptr );
    }

    // Keys are legend labels, values are SliceHistogram objects containing
    // true-space predictions from the corresponding generator models
    auto* slice_gen_map_ptr = new std::map< std::string, SliceHistogram* >();
    auto& slice_gen_map = *slice_gen_map_ptr;

    slice_gen_map[ "unfolded data" ] = slice_unf;
    if ( using_fake_data ) {
      slice_gen_map[ "truth" ] = slice_truth;
    }
    slice_gen_map[ "MicroBooNE Tune" ] = slice_cv;

    for ( const auto& pair : generator_truth_map ) {
      const auto& model_name = pair.first;
      TMatrixD* truth_mat = pair.second;

      SliceHistogram* temp_slice = SliceHistogram::make_slice_histogram(
        *truth_mat, slice, nullptr );

      slice_gen_map[ model_name ] = temp_slice;
    }

    int var_count = 0;
    std::string diff_xsec_denom;
    std::string diff_xsec_units_denom;
    std::string diff_xsec_denom_latex;
    std::string diff_xsec_units_denom_latex;
    double other_var_width = 1.;
    for ( const auto& ov_spec : slice.other_vars_ ) {
      double high = ov_spec.high_bin_edge_;
      double low = ov_spec.low_bin_edge_;
      const auto& var_spec = sb.slice_vars_.at( ov_spec.var_index_ );
      if ( high != low && std::abs(high - low) < BIG_DOUBLE ) {
        ++var_count;
        other_var_width *= ( high - low );
        diff_xsec_denom += 'd' + var_spec.name_;
        diff_xsec_denom_latex += " d" + var_spec.latex_name_;
        const std::string& temp_units = var_spec.units_;
        if ( !temp_units.empty() ) {
          diff_xsec_units_denom += " / " + temp_units;
          diff_xsec_units_denom_latex += " / " + var_spec.latex_units_;
        }
      }
    }

    for ( size_t av_idx : slice.active_var_indices_ ) {
      const auto& var_spec = sb.slice_vars_.at( av_idx );
      const std::string& temp_name = var_spec.name_;
      if ( temp_name != "true bin number" ) {
        var_count += slice.active_var_indices_.size();
        diff_xsec_denom += 'd' + var_spec.name_;
        diff_xsec_denom_latex += " d" + var_spec.latex_name_;

        if ( !var_spec.units_.empty() ) {
          diff_xsec_units_denom += " / " + var_spec.units_;
          diff_xsec_units_denom_latex += " / " + var_spec.latex_units_;
        }
      }
    }

    // NOTE: This currently assumes that each slice is a 1D histogram
    // TODO: revisit as needed
    int num_slice_bins = slice_unf->hist_->GetNbinsX();
    TMatrixD trans_mat( num_slice_bins, num_slice_bins );
    for ( int b = 0; b < num_slice_bins; ++b ) {
      double width = slice_unf->hist_->GetBinWidth( b + 1 );
      width *= other_var_width;
      trans_mat( b, b ) = 1e38 / ( width * integ_flux * num_Ar );
    }

    std::string slice_y_title;
    std::string slice_y_latex_title;
    if ( var_count > 0 ) {
      slice_y_title += "d";
      slice_y_latex_title += "{$d";
      if ( var_count > 1 ) {
        slice_y_title += "^{" + std::to_string( var_count ) + "}";
        slice_y_latex_title += "^{" + std::to_string( var_count ) + "}";
      }
      slice_y_title += "#sigma/" + diff_xsec_denom;
      slice_y_latex_title += "\\sigma / " + diff_xsec_denom_latex;
    }
    else {
      slice_y_title += "#sigma";
      slice_y_latex_title += "\\sigma";
    }
    slice_y_title += " (10^{-38} cm^{2}" + diff_xsec_units_denom + " / Ar)";
    slice_y_latex_title += "\\text{ }(10^{-38}\\text{ cm}^{2}"
      + diff_xsec_units_denom_latex + " / \\mathrm{Ar})$}";

    // Convert all slice histograms from true event counts to differential
    // cross-section units
    for ( auto& pair : slice_gen_map ) {
      auto* slice_h = pair.second;
      slice_h->transform( trans_mat );
      slice_h->hist_->GetYaxis()->SetTitle( slice_y_title.c_str() );
    }

    // Also transform all of the unfolded data slice histograms which have
    // specific covariance matrices
    for ( auto& sh_cov_pair : sh_cov_map ) {
      auto& slice_h = sh_cov_pair.second;
      slice_h->transform( trans_mat );
    }

    // Keys are generator legend labels, values are the results of a chi^2
    // test compared to the unfolded data (or, in the case of the unfolded
    // data, to the fake data truth)
    std::map< std::string, SliceHistogram::Chi2Result > chi2_map;
    std::cout << '\n';
    for ( const auto& pair : slice_gen_map ) {
      const auto& name = pair.first;
      const auto* slice_h = pair.second;

      // Decide what other slice histogram should be compared to this one,
      // then calculate chi^2
      SliceHistogram* other = nullptr;
      // We don't need to compare the unfolded data to itself, so just skip to
      // the next SliceHistogram and leave a dummy Chi2Result object in the map
      if ( name == "unfolded data" ) {
        chi2_map[ name ] = SliceHistogram::Chi2Result();
        continue;
      }
      // Compare all other distributions to the unfolded data
      else {
        other = slice_gen_map.at( "unfolded data" );
      }

      // Store the chi^2 results in the map
      const auto& chi2_result = chi2_map[ name ] = slice_h->get_chi2( *other );

      std::cout << "Slice " << sl_idx << ", " << name << ": \u03C7\u00b2 = "
        << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bin";
      if ( chi2_result.num_bins_ > 1 ) std::cout << 's';
      std::cout << ", p-value = " << chi2_result.p_value_ << '\n';
    }

    TCanvas* c1 = new TCanvas;
    slice_unf->hist_->SetLineColor( kBlack );
    slice_unf->hist_->SetLineWidth( 3 );
    slice_unf->hist_->SetMarkerStyle( kFullCircle );
    slice_unf->hist_->SetMarkerSize( 0.8 );
    slice_unf->hist_->SetStats( false );

    double ymax = -DBL_MAX;
    slice_unf->hist_->Draw( "e" );
    for ( const auto& pair : slice_gen_map ) {
      const auto& name = pair.first;
      const auto* slice_h = pair.second;

      double max = slice_h->hist_->GetMaximum();
      if ( max > ymax ) ymax = max;

      if ( name == "unfolded data" || name == "truth"
        || name == "MicroBooNE Tune" ) continue;

      const auto& file_info = truth_file_map.at( name );
      slice_h->hist_->SetLineColor( file_info.color_ );
      slice_h->hist_->SetLineStyle( file_info.style_ );
      slice_h->hist_->SetLineWidth( 4 );

      slice_h->hist_->Draw( "hist same" );
    }

    slice_cv->hist_->SetStats( false );
    slice_cv->hist_->SetLineColor( kAzure - 7 );
    slice_cv->hist_->SetLineWidth( 5 );
    slice_cv->hist_->SetLineStyle( 5 );
    slice_cv->hist_->Draw( "hist same" );

    if ( using_fake_data ) {
      slice_truth->hist_->SetStats( false );
      slice_truth->hist_->SetLineColor( kOrange );
      slice_truth->hist_->SetLineWidth( 5 );
      slice_truth->hist_->Draw( "hist same" );
    }

    slice_unf->hist_->GetYaxis()->SetRangeUser( 0., ymax*1.07 );
    slice_unf->hist_->Draw( "e same" );

    TLegend* lg = new TLegend( 0.15, 0.6, 0.5, 0.88 );
    for ( const auto& pair : slice_gen_map ) {
      const auto& name = pair.first;
      const auto* slice_h = pair.second;

      std::string label = name;
      std::ostringstream oss;
      const auto& chi2_result = chi2_map.at( name );
      oss << std::setprecision( 3 ) << chi2_result.chi2_ << " / "
        << chi2_result.num_bins_ << " bin";
      if ( chi2_result.num_bins_ > 1 ) oss << 's';

      if ( name != "unfolded data" ) {
        label += ": #chi^{2} = " + oss.str();
      }

      lg->AddEntry( slice_h->hist_.get(), label.c_str(), "l" );
    }

    lg->Draw( "same" );

    // Dump the unfolded results to text files compatible with PGFPlots
    std::map< std::string, std::vector<double> > slice_hist_table;
    std::map< std::string, std::string > slice_params_table;

    dump_slice_variables( sb, sl_idx, slice_params_table );

    for ( const auto& pair : slice_gen_map ) {
      const auto hist_name = samples_to_hist_names.at( pair.first );
      const auto* slice_hist = pair.second;
      bool include_x_coords = ( hist_name == "UnfData" );
      bool include_y_error = include_x_coords;
      dump_slice_histogram( hist_name, *slice_hist, slice, slice_hist_table,
        include_y_error, include_x_coords );
    }

    dump_slice_plot_limits( *slice_unf, *slice_cv, slice, slice_params_table );

    dump_slice_errors( "UnfData", slice, sh_cov_map, slice_hist_table );

    // Dump the chi^2 test results
    for ( const auto& chi2_pair : chi2_map ) {
      const auto hist_name = samples_to_hist_names.at( chi2_pair.first );
      const auto& chi2_result = chi2_pair.second;

      // Comparing the data histogram to itself is trivial, so skip it
      if ( hist_name == "UnfData" ) continue;
      else {
        slice_params_table[ hist_name + "_chi2" ]
          = std::to_string( chi2_result.chi2_ );
        slice_params_table[ hist_name + "_pvalue" ]
          = std::to_string( chi2_result.p_value_ );
      }
    }

    // Dump the total data POT and number of bins in the slice
    slice_params_table[ "bnb_data_pot" ] = std::to_string( total_pot );
    slice_params_table[ "num_bins" ] = std::to_string( num_slice_bins );

    // Dump a LaTeX title for the y-axis
    slice_params_table[ "y_axis_title" ] = slice_y_latex_title;

    // Before moving on to the next slice, dump information about the
    // current one to new pgfplots files that can be used for offline plotting
    std::string output_file_prefix = "dump/pgfplots_slice_";
    // Use at least three digits for numbering the slice output files
    if ( sl_idx < 10 ) output_file_prefix += '0';
    if ( sl_idx < 100 ) output_file_prefix += '0';
    output_file_prefix += std::to_string( sl_idx );

    write_pgfplots_files( output_file_prefix, slice_hist_table,
      slice_params_table );

  } // slices

}

int main() {
  test_unfolding_refactor();
  return 0;
}
