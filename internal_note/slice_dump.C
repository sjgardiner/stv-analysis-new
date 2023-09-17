// Standard library includes
#include <algorithm>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"

// STV analysis includes
#include "../FilePropertiesManager.hh"
#include "../MCC9SystematicsCalculator.hh"
#include "../PGFPlotsDumpUtils.hh"
#include "../PlotUtils.hh"
#include "../SliceBinning.hh"
#include "../SliceHistogram.hh"

using NFT = NtupleFileType;

void slice_dump() {

  //// Initialize the FilePropertiesManager and tell it to treat the NuWro
  //// MC ntuples as if they were data
  //auto& fpm = FilePropertiesManager::Instance();
  //fpm.load_file_properties( "../nuwro_file_properties.txt" );

  auto* syst_ptr = new MCC9SystematicsCalculator(
    "/uboone/data/users/gardiner/23-sept10-all-universes.root",
    //"../systcalc_unfold_fd.conf" );
    "../systcalc.conf" );
  auto& syst = *syst_ptr;

  // Get access to the relevant histograms owned by the SystematicsCalculator
  // object. These contain the reco bin counts that we need to populate the
  // slices below.
  TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();

  TH2D* category_hist = syst.cv_universe().hist_categ_.get();

  // Total MC+EXT prediction in reco bin space. Start by getting EXT.
  TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
    reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
  reco_mc_plus_ext_hist->SetDirectory( nullptr );

  // Add in the CV MC prediction
  reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );

  // Keys are covariance matrix types, values are CovMatrix objects that
  // represent the corresponding matrices
  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;

  auto* sb_ptr = new SliceBinning( "../mybins_all.txt" );
  auto& sb = *sb_ptr;

  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {

    const auto& slice = sb.slices_.at( sl_idx );

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

    // Get the binning and axis labels for the current slice by cloning the
    // (empty) histogram owned by the Slice object
    TH1* slice_hist = dynamic_cast< TH1* >(
      slice.hist_->Clone("slice_hist") );

    slice_hist->SetDirectory( nullptr );

    // Keys are labels, values are fractional uncertainty histograms
    auto* fr_unc_hists = new std::map< std::string, TH1* >();
    auto& frac_uncertainty_hists = *fr_unc_hists;

    // Show fractional uncertainties computed using these covariance matrices
    // in the ROOT plot. All configured fractional uncertainties will be
    // included in the output pgfplots file regardless of whether they appear
    // in this vector.
    const std::vector< std::string > cov_mat_keys = { "total",
      "detVar_total", "flux", "reint", "xsec_total", "POT", "numTargets",
      "MCstats", "EXTstats", "BNBstats"
    };

    // Loop over the various systematic uncertainties
    int color = 0;
    for ( const auto& pair : matrix_map ) {

      const auto& key = pair.first;
      const auto& cov_matrix = pair.second;

      SliceHistogram* slice_for_syst = SliceHistogram::make_slice_histogram(
        *reco_mc_plus_ext_hist, slice, &cov_matrix );

      // The SliceHistogram object already set the bin errors appropriately
      // based on the slice covariance matrix. Just change the bin contents
      // for the current histogram to be fractional uncertainties. Also set
      // the "uncertainties on the uncertainties" to zero.
      // TODO: revisit this last bit, possibly assign bin errors here
      for ( const auto& bin_pair : slice.bin_map_ ) {
        int global_bin_idx = bin_pair.first;
        double y = slice_for_syst->hist_->GetBinContent( global_bin_idx );
        double err = slice_for_syst->hist_->GetBinError( global_bin_idx );
        double frac = 0.;
        if ( y > 0. ) frac = err / y;
        slice_for_syst->hist_->SetBinContent( global_bin_idx, frac );
        slice_for_syst->hist_->SetBinError( global_bin_idx, 0. );
      }

      // Dump the current fractional errors in new columns for use with
      // pgfplots
      std::string frac_err_column_name( "frac_err_" + key );
      dump_slice_histogram( frac_err_column_name, *slice_for_syst, slice,
        pgfplots_hist_table, false, false );

      // Check whether the current covariance matrix name is present in
      // the vector defined above this loop. If it isn't, don't bother to
      // plot it, and just move on to the next one.
      auto cbegin = cov_mat_keys.cbegin();
      auto cend = cov_mat_keys.cend();
      auto iter = std::find( cbegin, cend, key );
      if ( iter == cend ) continue;

      frac_uncertainty_hists[ key ] = slice_for_syst->hist_.get();

      if ( color <= 9 ) ++color;
      if ( color == 5 ) ++color;
      if ( color >= 10 ) color += 10;

      slice_for_syst->hist_->SetLineColor( color );
      slice_for_syst->hist_->SetLineWidth( 3 );
    }

    TCanvas* c2 = new TCanvas;
    TLegend* lg2 = new TLegend( 0.7, 0.7, 0.9, 0.9 );

    auto* total_frac_err_hist = frac_uncertainty_hists.at( "total" );
    total_frac_err_hist->SetStats( false );
    total_frac_err_hist->GetYaxis()->SetRangeUser( 0.,
      total_frac_err_hist->GetMaximum() * 1.05 );
    total_frac_err_hist->SetLineColor( kBlack );
    total_frac_err_hist->SetLineWidth( 3 );
    total_frac_err_hist->Draw( "hist" );

    lg2->AddEntry( total_frac_err_hist, "total", "l" );

    for ( auto& pair : frac_uncertainty_hists ) {
      const auto& name = pair.first;
      TH1* hist = pair.second;
      // We already plotted the "total" one above
      if ( name == "total" ) continue;

      lg2->AddEntry( hist, name.c_str(), "l" );
      hist->Draw( "same hist" );

      std::cout << name << " frac err in bin #1 = "
        << hist->GetBinContent( 1 )*100. << "%\n";
    }

    lg2->Draw( "same" );

    std::cout << "Total frac error in bin #1 = "
      << total_frac_err_hist->GetBinContent( 1 )*100. << "%\n";

    // Before moving on to the next slice, dump information about the
    // current one to new pgfplots files that can be used for offline plotting
    std::string output_file_prefix = "pgfplots_slice_";
    // Use at least three digits for numbering the slice output files
    if ( sl_idx < 10 ) output_file_prefix += '0';
    if ( sl_idx < 100 ) output_file_prefix += '0';
    output_file_prefix += std::to_string( sl_idx );

    write_pgfplots_files( output_file_prefix, pgfplots_hist_table,
      pgfplots_params_table );

  } // slices

}

int main() {
  slice_dump();
  return 0;
}
