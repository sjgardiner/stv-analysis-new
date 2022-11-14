// Standard library includes
#include <set>

// ROOT includes
#include "TFile.h"

// STV analysis includes
#include "../FilePropertiesManager.hh"
#include "../MCC9SystematicsCalculator.hh"
#include "../PlotUtils.hh"

// ROOT integer code for Arial font
constexpr int FONT_STYLE = 62; // Arial

void dump_pgfplots_smearing_histogram( const std::string& output_table_file,
  const std::string& output_params_file, const TH2D* smear_hist,
  size_t first_bkgd_true_bin_idx )
{
  // Don't do anything if you get a nullptr for the histogram
  if ( !smear_hist ) return;

  std::ofstream out_table_file( output_table_file );
  out_table_file << "xbin  ybin  z\n";
  int num_x_bins = smear_hist->GetXaxis()->GetNbins();
  int num_y_bins = smear_hist->GetYaxis()->GetNbins();
  for ( int xb = 1; xb <= num_x_bins; ++xb ) {
    for ( int yb = 1; yb <= num_y_bins; ++yb ) {
      double z = smear_hist->GetBinContent( xb, yb );
      // Use zero-based bin indices in the dump
      out_table_file << xb - 1 << "  " <<  yb - 1 << "  " << z << '\n';
    }
  }

  // Store the bin counts along each axis in the parameters file as well
  // as the first background bin index along the true (x) axis. This will
  // facilitate easy plotting by PGFPlots.
  std::ofstream out_params_file( output_params_file );
  out_params_file << "numXbins  numYbins  firstBkgdTrueBinIdx\n";
  out_params_file << num_x_bins << "  " << num_y_bins << "  "
    << first_bkgd_true_bin_idx;
}

void smear_matrix() {

  const std::string input_respmat_file_name( "/uboone/data/users/"
    "gardiner/myuniverses-all.root" );

  auto* syst_ptr = new MCC9SystematicsCalculator( input_respmat_file_name,
    "../systcalc.conf" );
  auto& syst = *syst_ptr;

  // Build a set of integer block indices for the reco bins
  const auto& true_bins = syst.true_bins_;
  const auto& reco_bins = syst.reco_bins_;

  std::set< int > reco_blocks;
  size_t num_reco_bins = reco_bins.size();
  for ( const auto& rb : reco_bins ) {
    reco_blocks.insert( rb.block_index_ );
  }

  // Create a new TH2D to store the smearing matrix. Use the 2D event counts
  // (in true and reco bins) from the CV universe to start.
  const auto& cv_univ = syst.cv_universe();
  TH2D* smear_hist = dynamic_cast< TH2D* >( cv_univ.hist_2d_
    ->Clone("smear_hist") );

  // Get the bin index for the first true bin that represents background events
  size_t num_true_bins = true_bins.size();
  size_t first_bkgd_bin_idx = num_true_bins;
  for ( size_t t = 0u; t < num_true_bins; ++t ) {
    const auto& tbin = true_bins.at( t );
    if ( tbin.type_ == TrueBinType::kBackgroundTrueBin ) {
      first_bkgd_bin_idx = t;
      break;
    }
  }

  TH1D* expected_signal_hist = smear_hist->ProjectionY(
    "expected_signal_hist", 1, first_bkgd_bin_idx );
  TCanvas* c = new TCanvas;
  expected_signal_hist->SetLineWidth( 3 );
  expected_signal_hist->SetStats( false );
  expected_signal_hist->GetYaxis()->SetTitle( "expected signal events" );
  expected_signal_hist->Draw( "hist e" );

  // Normalize each block of the smearing matrix elements so that a sum over
  // all reco bins that belong to the same block (including the under/overflow
  // bins) yields a value of one. This means that every selected signal event
  // must end up somewhere in reco space in each block.
  for ( const int cur_block_idx : reco_blocks ) {

    // Loop over all true bins
    for ( size_t bt = 0; bt < num_true_bins; ++bt ) {

      // Skip true bins that do not belong to the current block
      const auto& tb = true_bins.at( bt );
      if ( tb.block_index_ != cur_block_idx ) continue;

      // For the current true (x) bin, compute the sum of all reco (y) bins
      // that belong to the current block
      double y_sum = 0.;
      for ( size_t br = 0; br < num_reco_bins; ++br ) {
        const auto& rb = reco_bins.at( br );
        if ( tb.block_index_ != rb.block_index_ ) continue;
        y_sum += smear_hist->GetBinContent( bt + 1, br + 1 );
      }

      // Normalize each of the reco (y) bins in the current block so that the
      // sum over y is unity.
      for ( size_t br = 0; br < num_reco_bins; ++br ) {

        const auto& rb = reco_bins.at( br );
        if ( tb.block_index_ != rb.block_index_ ) continue;

        // To avoid dividing by zero, set the bin content to zero if the sum of
        // the reco (y) bins is not positive.
        if ( y_sum <= 0. ) {
          smear_hist->SetBinContent( bt + 1, br + 1, 0. );
        }
        else {
          // Otherwise, normalize in the usual way
          double bc = smear_hist->GetBinContent( bt + 1, br + 1 );

          double content = std::max( bc / y_sum, 0. );

          smear_hist->SetBinContent( bt + 1, br + 1, content );
        }
      } // loop over reco (y) bins

    } // loop over true (x) bins

  } // loop over blocks

  // Smearing matrix histogram style options
  smear_hist->GetXaxis()->SetTitleFont( FONT_STYLE);
  smear_hist->GetYaxis()->SetTitleFont( FONT_STYLE );
  smear_hist->GetXaxis()->SetTitleSize( 0.05 );
  smear_hist->GetYaxis()->SetTitleSize( 0.05 );
  smear_hist->GetXaxis()->SetLabelFont( FONT_STYLE );
  smear_hist->GetYaxis()->SetLabelFont( FONT_STYLE );
  smear_hist->GetZaxis()->SetLabelFont( FONT_STYLE );
  smear_hist->GetZaxis()->SetLabelSize( 0.03 );
  smear_hist->GetXaxis()->CenterTitle();
  smear_hist->GetYaxis()->CenterTitle();
  smear_hist->GetXaxis()->SetTitleOffset( 1.2 );
  smear_hist->GetYaxis()->SetTitleOffset( 1.1 );
  smear_hist->SetStats( false );
  smear_hist->SetMarkerSize( 1.8 ); // text size
  smear_hist->SetMarkerColor( kWhite ); // text color

  // Draw the smearing matrix plot
  TCanvas* c_smear = new TCanvas;
  c_smear->SetBottomMargin( 0.15 );
  c_smear->SetLeftMargin( 0.13 );

  smear_hist->Draw( "colz" );

  dump_pgfplots_smearing_histogram( "all_migmat_table.txt",
    "all_migmat_params.txt", smear_hist, first_bkgd_bin_idx );

  // Dump a table of bin occupancies
  std::ofstream occup_table_file( "occupancy_table.txt" );
  for ( size_t bin_idx = 0u; bin_idx < first_bkgd_bin_idx; ++bin_idx ) {
    occup_table_file << bin_idx << " & "
      << smear_hist->GetBinContent( bin_idx + 1, bin_idx + 1 ) << '\n';
  }

}
