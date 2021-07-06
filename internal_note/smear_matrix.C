// ROOT includes
#include "TFile.h"

// STV analysis includes
#include "../FilePropertiesManager.hh"
#include "../PlotUtils.hh"
#include "../SystematicsCalculator.hh"

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

  const std::string input_respmat_file_name(
    "/uboone/data/users/gardiner/respmat_mcc8-cth_p.root" );

  auto* syst_ptr = new SystematicsCalculator( input_respmat_file_name );
  auto& syst = *syst_ptr;

  // Create a new TH2D to store the smearing matrix. Use the 2D event counts
  // (in true and reco bins) from the CV universe to start.
  const auto& cv_univ = syst.cv_universe();
  TH2D* smear_hist = dynamic_cast< TH2D* >( cv_univ.hist_2d_
    ->Clone("smear_hist") );

  // Normalize the smearing matrix elements so that a sum over all reco bins
  // (including the under/overflow bins) yields a value of one. This means that
  // every selected signal event must end up somewhere in reco space.
  int num_bins_x = smear_hist->GetXaxis()->GetNbins();
  int num_bins_y = smear_hist->GetYaxis()->GetNbins();

  // Loop over the true (x) bins. Include the underflow (index zero) and
  // overflow (index num_bins_x + 1) bins.
  for ( int bx = 0; bx <= num_bins_x + 1; ++bx ) {

    // For the current true (x) bin, compute the sum of all reco (y) bins.
    double y_sum = 0.;
    for ( int by = 0; by <= num_bins_y + 1; ++by ) {
      y_sum += smear_hist->GetBinContent( bx, by );
    }

    // Normalize each of the reco (y) bins so that the sum over y is unity.
    for ( int by = 0; by <= num_bins_y + 1; ++by ) {

      // To avoid dividing by zero, set the bin content to zero if the sum of
      // the reco (y) bins is not positive.
      if ( y_sum <= 0. ) {
        smear_hist->SetBinContent( bx, by, 0. );
      }
      else {
        // Otherwise, normalize in the usual way
        double bc = smear_hist->GetBinContent( bx, by );

        double content = std::max( bc / y_sum, 0. );

        smear_hist->SetBinContent( bx, by, content );
      }
    } // loop over reco (y) bins

  } // loop over true (x) bins

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

  // Get the bin index for the first true bin that represents background events
  const auto& true_bins = syst.true_bins_;
  size_t num_true_bins = true_bins.size();
  size_t first_bkgd_bin_idx = num_true_bins;
  for ( size_t t = 0u; t < num_true_bins; ++t ) {
    const auto& tbin = true_bins.at( t );
    if ( tbin.type_ == TrueBinType::kBackgroundTrueBin ) {
      first_bkgd_bin_idx = t;
      break;
    }
  }

  dump_pgfplots_smearing_histogram( "mcc8CthP_respmat_table.txt",
    "mcc8CthP_respmat_params.txt", smear_hist, first_bkgd_bin_idx );

}
