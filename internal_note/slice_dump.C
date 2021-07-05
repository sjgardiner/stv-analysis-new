// Standard library includes
#include <algorithm>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"

// STV analysis includes
#include "../FilePropertiesManager.hh"
#include "../MCC9Unfolder.hh"
#include "../PlotUtils.hh"
#include "../SliceBinning.hh"
#include "../SliceHistogram.hh"

using NFT = NtupleFileType;

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

    std::string name_for_dump = '{' + svar.latex_name_ + '}';
    std::string units_for_dump = '{' + svar.latex_units_ + '}';

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

    std::string name_for_dump = '{' + svar.latex_name_ + '}';
    std::string units_for_dump = '{' + svar.latex_units_ + '}';

    double low = ov_spec.low_bin_edge_;
    double high = ov_spec.high_bin_edge_;

    std::string col_prefix = "ov" + std::to_string( b );
    pgfplots_params_table[ col_prefix + "_name" ] = name_for_dump;
    pgfplots_params_table[ col_prefix + "_units" ] = units_for_dump;
    pgfplots_params_table[ col_prefix + "_low" ] = std::to_string( low );
    pgfplots_params_table[ col_prefix + "_high" ] = std::to_string( high );
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



void slice_dump() {

  auto* syst_ptr = new MCC9Unfolder( "/uboone/data/users/gardiner/"
    "ntuples-stv-MCC9InternalNote/respmat-files/RespMat-mcc9-2D_proton.root",
    "../systcalc.conf" );
  auto& syst = *syst_ptr;

  //// Keys are covariance matrix types, values are CovMatrix objects that
  //// represent the corresponding matrices
  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;

  auto* sb_ptr = new SliceBinning( "../mybins_mcc9_2D_proton.txt" );
  auto& sb = *sb_ptr;

  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {

    const auto& slice = sb.slices_.at( sl_idx );

    // Get access to the relevant histograms owned by the SystematicsCalculator
    // object. These contain the reco bin counts we need to populate the
    // current slice
    TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
    TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();
    TH2D* category_hist = syst.cv_universe().hist_categ_.get();

    // Total MC+EXT prediction in reco bin space. Start by getting EXT.
    TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
      reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
    reco_mc_plus_ext_hist->SetDirectory( nullptr );

    // Add in the CV MC prediction
    reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );

    // We now have all of the reco bin space histograms that we need as input.
    // Use them to make new histograms in slice space.
    SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(
      *reco_bnb_hist, slice, &matrix_map.at("BNBstats") );

    SliceHistogram* slice_ext = SliceHistogram::make_slice_histogram(
      *reco_ext_hist, slice, &matrix_map.at("EXTstats") );

    SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
      *reco_mc_plus_ext_hist, slice, &matrix_map.at("total") );

    // Prepare to dump to a couple of input files suitable for use with the
    // PGFPlots LaTeX package

    // Histogram counts
    std::map< std::string, std::vector<double> > pgfplots_hist_table;
    // Overall slice parameters (data POT, etc.)
    std::map< std::string, std::string > pgfplots_params_table;

    // Store parameters describing the current slice
    pgfplots_params_table[ "bnb_data_pot" ]
      = std::to_string( syst.total_bnb_data_pot_ );

    // Use TH1::Integral() to correctly take event weights into account
    double ext_plus_mc_integral = slice_mc_plus_ext->hist_->Integral();
    double ext_integral = slice_ext->hist_->Integral();
    double ext_fraction = ext_integral / ext_plus_mc_integral;

    pgfplots_params_table[ "ext_fraction" ] = std::to_string( ext_fraction );

    auto chi2_result = slice_bnb->get_chi2( *slice_mc_plus_ext );
    std::cout << "Slice " << sl_idx << ": \u03C7\u00b2 = "
      << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bins,"
      << " p-value = " << chi2_result.p_value_ << '\n';

    pgfplots_params_table[ "chi2" ] = std::to_string( chi2_result.chi2_ );
    pgfplots_params_table[ "pvalue" ] = std::to_string( chi2_result.p_value_ );
    pgfplots_params_table[ "num_bins" ]
      = std::to_string( chi2_result.num_bins_ );
    pgfplots_params_table[ "dof" ] = std::to_string( chi2_result.dof_ );

    dump_slice_plot_limits( *slice_bnb, *slice_mc_plus_ext, slice,
      pgfplots_params_table );

    dump_slice_variables( sb, sl_idx, pgfplots_params_table );

    // Prepare the PGFPlots information for various histograms
    dump_slice_histogram( "BNB", *slice_bnb, slice,
      pgfplots_hist_table, true, true );
    dump_slice_histogram( "MC+EXT", *slice_mc_plus_ext, slice,
      pgfplots_hist_table, true, false );
    dump_slice_histogram( "EXT", *slice_ext, slice,
      pgfplots_hist_table, true, false );

    // Build a stack of categorized central-value MC predictions plus the
    // extBNB contribution in slice space
    const auto& eci = EventCategoryInterpreter::Instance();
    eci.set_ext_histogram_style( slice_ext->hist_.get() );

    THStack* slice_pred_stack = new THStack( "mc+ext", "" );
    slice_pred_stack->Add( slice_ext->hist_.get() ); // extBNB

    const auto& cat_map = eci.label_map();

    // Go in reverse so that signal ends up on top. Note that this index is
    // one-based to match the ROOT histograms
    int cat_bin_index = cat_map.size();
    for ( auto iter = cat_map.crbegin(); iter != cat_map.crend(); ++iter )
    {
      EventCategory cat = iter->first;
      TH1D* temp_mc_hist = category_hist->ProjectionY( "temp_mc_hist",
        cat_bin_index, cat_bin_index );
      temp_mc_hist->SetDirectory( nullptr );

      SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(
        *temp_mc_hist, slice  );

      double category_integral = temp_slice_mc->hist_->Integral();
      double category_fraction = category_integral / ext_plus_mc_integral;
      std::string category_label = "MC" + std::to_string( cat ) + "_fraction";
      pgfplots_params_table[ category_label ]
        = std::to_string( category_fraction );

      eci.set_mc_histogram_style( cat, temp_slice_mc->hist_.get() );

      slice_pred_stack->Add( temp_slice_mc->hist_.get() );

      std::string cat_col_prefix = "MC" + std::to_string( cat );
      dump_slice_histogram( cat_col_prefix, *temp_slice_mc, slice,
        pgfplots_hist_table, false, false );

      --cat_bin_index;
    }

    // We're ready to dump information about the current slice. Write
    // the PGFPlots files now before making a ROOT plot
    std::string output_file_prefix = "pgfplots_slice_";
    if ( sl_idx < 10 ) output_file_prefix += '0';
    output_file_prefix += std::to_string( sl_idx );

    write_pgfplots_files( output_file_prefix, pgfplots_hist_table,
      pgfplots_params_table );

    TCanvas* c1 = new TCanvas;
    slice_bnb->hist_->SetLineColor( kBlack );
    slice_bnb->hist_->SetLineWidth( 3 );
    slice_bnb->hist_->SetMarkerStyle( kFullCircle );
    slice_bnb->hist_->SetMarkerSize( 0.8 );
    slice_bnb->hist_->SetStats( false );
    double ymax = std::max( slice_bnb->hist_->GetMaximum(),
      slice_mc_plus_ext->hist_->GetMaximum() ) * 1.07;
    slice_bnb->hist_->GetYaxis()->SetRangeUser( 0., ymax );

    slice_bnb->hist_->Draw( "e" );

    slice_pred_stack->Draw( "hist same" );

    slice_mc_plus_ext->hist_->SetLineWidth( 3 );
    slice_mc_plus_ext->hist_->Draw( "same hist e" );

    slice_bnb->hist_->Draw( "same e" );

    //std::string out_pdf_name = "plot_slice_";
    //if ( sl_idx < 10 ) out_pdf_name += "0";
    //out_pdf_name += std::to_string( sl_idx ) + ".pdf";
    //c1->SaveAs( out_pdf_name.c_str() );
  }

}

int main() {
  slice_dump();
  return 0;
}
