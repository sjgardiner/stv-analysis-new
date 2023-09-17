// Standard library includes
#include <climits>
#include <set>

// ROOT includes
#include "TCanvas.h"
#include "TFile.h"

// STV analysis includes
#include "../FilePropertiesManager.hh"
#include "../MCC9SystematicsCalculator.hh"
#include "../PlotUtils.hh"

// ROOT integer code for Arial font
constexpr int FONT_STYLE = 62; // Arial

void dump_pgfplots_smearing_histogram( const std::string& output_table_stem,
  const TH2D* smear_hist, const std::set< int >& reco_blocks,
  const SystematicsCalculator& syst )
{
  // Don't do anything if you get a nullptr for the histogram
  if ( !smear_hist ) return;

  for ( const int cur_block : reco_blocks ) {

    int num_x_bins_in_block = 0;
    size_t first_x_bin_idx_in_block = UINT_MAX;
    size_t last_x_bin_idx_in_block = 0u;
    for ( size_t tb = 0u; tb < syst.true_bins_.size(); ++tb ) {
      const auto& tbin = syst.true_bins_.at( tb );
      if ( tbin.block_index_ == cur_block ) {
        ++num_x_bins_in_block;
        if ( tb < first_x_bin_idx_in_block ) {
          first_x_bin_idx_in_block = tb;
        }
        if ( tb > last_x_bin_idx_in_block ) {
          last_x_bin_idx_in_block = tb;
        }
      }
    }

    int num_y_bins_in_block = 0;
    size_t first_y_bin_idx_in_block = UINT_MAX;
    size_t last_y_bin_idx_in_block = 0;
    for ( size_t rb = 0u; rb < syst.reco_bins_.size(); ++rb ) {
      const auto& rbin = syst.reco_bins_.at( rb );
      if ( rbin.block_index_ == cur_block ) {
        ++num_y_bins_in_block;
        if ( rb < first_y_bin_idx_in_block ) {
          first_y_bin_idx_in_block = rb;
        }
        if ( rb > last_y_bin_idx_in_block ) {
          last_y_bin_idx_in_block = rb;
        }
      }
    }

    std::string out_file_name = output_table_stem + "_";

    // Use at least three digits for numbering the output files
    if ( cur_block < 10 ) out_file_name += '0';
    if ( cur_block < 100 ) out_file_name += '0';
    out_file_name += std::to_string( cur_block );

    std::ofstream out_table_file( out_file_name + "_hist.txt" );

    out_table_file << "xbin  ybin  z\n";
    int num_x_bins = smear_hist->GetXaxis()->GetNbins();
    int num_y_bins = smear_hist->GetYaxis()->GetNbins();
    for ( int xb = 1; xb <= num_x_bins; ++xb ) {
      for ( int yb = 1; yb <= num_y_bins; ++yb ) {

        int true_block = syst.true_bins_.at( xb - 1 ).block_index_;
        if ( true_block != cur_block ) continue;

        int reco_block = syst.reco_bins_.at( yb - 1 ).block_index_;
        if ( reco_block != cur_block ) continue;

        double z = smear_hist->GetBinContent( xb, yb );

        // Use zero-based bin indices in the dump
        out_table_file << xb - 1 << "  " <<  yb - 1 << "  " << z << '\n';

      } // loop over reco (y) bins
    } // loop over true (x) bins

    // Store the bin counts along each axis in the parameters file. This will
    // facilitate easy plotting by PGFPlots.
    std::ofstream out_params_file( out_file_name + "_params.txt" );
    out_params_file << "numXbins  numYbins  firstXbin  lastXbin";
    out_params_file << "  firstYbin  lastYbin\n";
    out_params_file << num_x_bins_in_block << "  " << num_y_bins_in_block;
    out_params_file << "  " << first_x_bin_idx_in_block << "  "
      << last_x_bin_idx_in_block;
    out_params_file << "  " << first_y_bin_idx_in_block << "  "
      << last_y_bin_idx_in_block;

  } // loop over blocks

}

void smear_matrix() {

  const std::string input_respmat_file_name(
    "/uboone/data/users/gardiner/23-sept10-all-universes.root" );

  auto* syst_ptr = new MCC9SystematicsCalculator( input_respmat_file_name,
    "../systcalc.conf" );
  auto& syst = *syst_ptr;

  // Build a set of integer block indices for the reco bins
  const auto& true_bins = syst.true_bins_;
  const auto& reco_bins = syst.reco_bins_;

  std::set< int > reco_blocks;
  size_t num_true_bins = true_bins.size();
  size_t num_reco_bins = reco_bins.size();
  for ( const auto& rb : reco_bins ) {
    reco_blocks.insert( rb.block_index_ );
  }

  // Create a new TH2D to store the smearing matrix. Use the 2D event counts
  // (in true and reco bins) from the CV universe to start.
  const auto& cv_univ = syst.cv_universe();
  // Use alternate CV universe (currently NuWro) instead
  //const auto& cv_univ = *syst.alt_cv_universes_.at(
  //  NtupleFileType::kAltCVMC );
  TH2D* smear_hist = dynamic_cast< TH2D* >( cv_univ.hist_2d_
    ->Clone("smear_hist") );

  // Make empty 1D histograms with the same reco bin structure as the smearing
  // matrix. These will be populated with the expected selected signal events
  // and the efficiencies below.
  TH1D* expected_sel_signal_hist
    = smear_hist->ProjectionY( "expected_sel_signal_hist" );
  expected_sel_signal_hist->Reset();

  TH1D* efficiency_hist
    = smear_hist->ProjectionY( "efficiency_hist" );
  efficiency_hist->Reset();

  // Calculate the expected number of selected signal events in each reco
  // bin while avoiding double-counting across blocks
  for ( const int cur_block_idx : reco_blocks ) {
    for ( size_t br = 0; br < num_reco_bins; ++br ) {

      // Skip reco bins that do not belong to the current block
      const auto& rb = reco_bins.at( br );
      if ( rb.block_index_ != cur_block_idx ) continue;

      // For the current reco (y) bin, compute the sum of all true (x) bins
      // that belong to the current block
      double x_sum = 0.;
      double x_sum_squared_errors = 0.;
      for ( size_t bt = 0; bt < num_true_bins; ++bt ) {
        const auto& tb = true_bins.at( bt );
        if ( tb.block_index_ != rb.block_index_ ) continue;
        x_sum += smear_hist->GetBinContent( bt + 1, br + 1 );
        double temp_err = smear_hist->GetBinError( bt + 1, br + 1 );
        x_sum_squared_errors += temp_err * temp_err;
      } // loop over true bins

      expected_sel_signal_hist->SetBinContent( br + 1, x_sum );
      double temp_bin_err = std::sqrt( std::max(0., x_sum_squared_errors) );
      expected_sel_signal_hist->SetBinError( br + 1, temp_bin_err );

    } // loop over reco bins

  } // loop over blocks

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

      // Now that we have the total number of reconstructed events for a given
      // true bin index, calculate the efficiency and store it in the 1D
      // histogram created for this purpose
      double num_true_evts = cv_univ.hist_true_->GetBinContent( bt + 1 );
      double efficiency = y_sum / num_true_evts;
      efficiency_hist->SetBinContent( bt + 1, efficiency );

      // See DocDB #32401, Eq. (5.2)
      double eff_stat_err = std::sqrt( std::max(0., efficiency
        * (1. - efficiency) / num_true_evts) );
      efficiency_hist->SetBinError( bt + 1, eff_stat_err );

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

  dump_pgfplots_smearing_histogram( "all_migmat_table", smear_hist,
    reco_blocks, syst );

  // Also dump the expected selected events and efficiencies
  std::ofstream out_1D_file( "all_sig_eff_table.txt" );
  out_1D_file << "bin  expected_sel_sig  ess_stat_error"
    << "  efficiency  eff_stat_error";
  int num_bins_to_dump = expected_sel_signal_hist->GetNbinsX();
  for ( int b = 0; b < num_bins_to_dump; ++b ) {
    double ess = expected_sel_signal_hist->GetBinContent( b + 1 );
    double ess_err = expected_sel_signal_hist->GetBinError( b + 1 );
    double eff = efficiency_hist->GetBinContent( b + 1 );
    double eff_err = efficiency_hist->GetBinError( b + 1 );
    out_1D_file << '\n' << b << "  " << ess << "  " << ess_err << "  "
      << eff << "  " << eff_err;
  }

  TCanvas* c = new TCanvas;
  expected_sel_signal_hist->SetLineWidth( 3 );
  expected_sel_signal_hist->SetStats( false );
  expected_sel_signal_hist->GetYaxis()->SetTitle( "expected signal events" );
  expected_sel_signal_hist->Draw( "hist e" );

  TCanvas* c2 = new TCanvas;
  efficiency_hist->SetLineWidth( 3 );
  efficiency_hist->SetStats( false );
  efficiency_hist->GetYaxis()->SetTitle( "efficiency" );
  efficiency_hist->Draw( "hist e" );
}

int main() {
  smear_matrix();
  return 0;
}
