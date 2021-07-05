#pragma once

// Standard library includes
#include <cmath>
#include <limits>
#include <memory>

// ROOT includes
#include "TH1.h"
#include "TMatrixD.h"
#include "Math/SpecFuncMathCore.h" // Needed for ROOT::Math::inc_gamma_c()

// STV analysis includes
#include "SliceBinning.hh"
#include "SystematicsCalculator.hh"

class SliceHistogram {

  public:

    SliceHistogram() {}

    static SliceHistogram* make_slice_histogram( TH1D& reco_bin_histogram,
      const Slice& slice, const CovMatrix* input_cov_mat = nullptr );

    struct Chi2Result {
      Chi2Result( double chi2, int nbins, int dof, double pval )
        : chi2_( chi2 ), num_bins_( nbins ), dof_( dof ), p_value_( pval ) {}
      double chi2_;
      int num_bins_;
      int dof_;
      double p_value_;
    };

    Chi2Result get_chi2( const SliceHistogram& other );

    std::unique_ptr< TH1 > hist_;
    CovMatrix cmat_;
};

// Creates a new event histogram and an associated covariance matrix for a
// particular slice of phase space. The histogram is filled from the
// appropriate bin(s) of a 1D histogram of reco bin event counts. The mapping
// from reco bin number to the slice histogram bins is described by the input
// Slice object. Bin errors are set according to the reco-bin-space CovMatrix
// object pointed to by the input_cov_mat argument. If it is null, the bin
// errors are set to a default value of zero, and the output CovMatrix object
// owns a nullptr.
SliceHistogram* SliceHistogram::make_slice_histogram( TH1D& reco_bin_histogram,
  const Slice& slice, const CovMatrix* input_cov_mat )
{
  // Get the binning and axis labels for the current slice by cloning the
  // (empty) histogram owned by the Slice object
  TH1* slice_hist = dynamic_cast< TH1* >(
    slice.hist_->Clone("slice_hist") );

  slice_hist->SetDirectory( nullptr );

  // Fill the slice bins based on the input reco bins
  for ( const auto& pair : slice.bin_map_ ) {

    // One-based index for the global TH1 bin number in the slice
    int slice_bin_idx = pair.first;

    const auto& reco_bin_set = pair.second;

    double slice_bin_content = 0.;
    for ( const auto& rb_idx : reco_bin_set ) {
      // The ResponseMatrixMaker reco bin indices are zero-based, so I correct
      // for this here when pulling values from the one-based input ROOT
      // histogram
      slice_bin_content += reco_bin_histogram.GetBinContent( rb_idx + 1 );
    }

    slice_hist->SetBinContent( slice_bin_idx, slice_bin_content );

  } // slice bins

  // If we've been handed a non-null pointer to a CovMatrix object, then
  // we will use it to propagate uncertainties.
  TH2D* covmat_hist = nullptr;
  if ( input_cov_mat ) {

    // Create a new TH2D to hold the covariance matrix elements associated with
    // the slice histogram.
    // NOTE: I assume here that every slice bin is represented in the bin_map.
    // If this isn't the case, the bin counting will be off.
    // TODO: revisit this assumption and perhaps do something better
    int num_slice_bins = slice.bin_map_.size();

    covmat_hist = new TH2D( "covmat_hist", "covariance; slice bin;"
      " slice bin; covariance", num_slice_bins, 0., num_slice_bins,
      num_slice_bins, 0., num_slice_bins );
    covmat_hist->SetDirectory( nullptr );
    covmat_hist->SetStats( false );

    // We're ready. Populate the new covariance matrix using the elements
    // of the one for the reco bin space
    for ( const auto& pair_a : slice.bin_map_ ) {
      // Global slice bin index
      int sb_a = pair_a.first;
      // Set of reco bins that correspond to slice bin sb_a
      const auto& rb_set_a = pair_a.second;
      for ( const auto& pair_b : slice.bin_map_ ) {
        int sb_b = pair_b.first;
        const auto& rb_set_b = pair_b.second;

        double cov = 0.;
        const TH2D* cmat = input_cov_mat->cov_matrix_.get();
        for ( const auto& rb_m : rb_set_a ) {
          for ( const auto& rb_n : rb_set_b ) {
            // The covariance matrix TH2D uses one-based indices even though
            // the ResponseMatrixMaker numbering scheme is zero-based. I
            // correct for this here.
            cov += cmat->GetBinContent( rb_m + 1, rb_n + 1 );
          } // reco bin index m
        } // reco bin index n
        covmat_hist->SetBinContent( sb_a, sb_b, cov );
      } // slice bin index b
    } // slice bin index a


    // We have a finished covariance matrix for the slice. Use it to set
    // the bin errors on the slice histogram.
    for ( const auto& pair : slice.bin_map_ ) {

      int slice_bin_idx = pair.first;
      double bin_variance = covmat_hist->GetBinContent( slice_bin_idx,
        slice_bin_idx );
      double bin_error = std::sqrt( std::max(0., bin_variance) );

      // This works for a multidimensional slice because a global bin index
      // (as returned by TH1::GetBin) is used for slice_bin_idx.
      slice_hist->SetBinError( slice_bin_idx, bin_error );

    } // slice bins

  } // non-null input_cov_mat

  // We're done. Prepare the SliceHistogram object and return it.
  auto* result = new SliceHistogram;
  result->hist_.reset( slice_hist );
  result->cmat_.cov_matrix_.reset( covmat_hist );

  return result;
}

SliceHistogram::Chi2Result SliceHistogram::get_chi2(
  const SliceHistogram& other )
{
  int num_bins = hist_->GetNbinsX();
  if ( other.hist_->GetNbinsX() != num_bins ) {
    throw std::runtime_error( "Incompatible vector sizes in chi^2"
      " calculation" );
  }

  int my_cov_mat_x_bins = cmat_.cov_matrix_->GetNbinsX();
  int my_cov_mat_y_bins = cmat_.cov_matrix_->GetNbinsY();
  int other_cov_mat_x_bins = other.cmat_.cov_matrix_->GetNbinsY();
  int other_cov_mat_y_bins = other.cmat_.cov_matrix_->GetNbinsY();
  if ( my_cov_mat_x_bins != num_bins
    || my_cov_mat_y_bins != num_bins
    || other_cov_mat_x_bins != num_bins
    || other_cov_mat_y_bins != num_bins )
  {
    throw std::runtime_error( "Invalid covariance matrix dimensions"
      " encountered in chi^2 calculation" );
  }

  // The total covariance matrix on the difference between the
  // two histograms is just the sum of each individual SliceHistogram's
  // owned covariance matrix.
  CovMatrix cov_mat;
  cov_mat += cmat_;
  cov_mat += other.cmat_;

  // Get access to a TMatrixDSym object representing the covariance matrix.
  // This supplies a convenient function for handling the inversion.
  auto cov_matrix = cov_mat.get_matrix();

  // Pre-scale before inversion to avoid numerical problems. Here we choose
  // a scaling factor such that the smallest nonzero entry in the original
  // covariance matrix has an absolute value of unity. Note the use of the
  // zero-based element indices for TMatrixDSym.
  constexpr double BIG_DOUBLE = std::numeric_limits<double>::max();
  double min_abs = BIG_DOUBLE;
  for ( int a = 0; a < num_bins; ++a ) {
    for ( int b = 0; b < num_bins; ++b ) {
      double element = cov_matrix->operator()( a, b );
      double abs_el = std::abs( element );
      if ( abs_el > 0. && abs_el < min_abs ) min_abs = abs_el;
    }
  }

  // If all covariance matrix elements are zero, then this scaling won't work
  // and something is wrong. Complain if this is the case.
  if ( min_abs == BIG_DOUBLE ) {
    throw std::runtime_error( "Null covariance matrix encountered in chi^2"
      " calculation" );
  }

  // We're ready. Do the actual pre-scaling here. Keep the first TMatrixDSym
  // and invert a copy. This allows us to check that the inversion was
  // successful below.
  auto inverse_cov_matrix = cov_mat.get_matrix();
  double scaling_factor = 1. / min_abs;
  inverse_cov_matrix->operator*=( scaling_factor );

  // Do the inversion
  inverse_cov_matrix->Invert();

  // Undo the scaling by re-applying it to the inverse matrix
  inverse_cov_matrix->operator*=( scaling_factor );

  // Double-check that we get a unit matrix by multiplying the
  // original by its inverse
  TMatrixD unit_mat( *cov_matrix, TMatrixD::kMult, *inverse_cov_matrix );
  constexpr double INVERSION_TOLERANCE = 1e-4;
  for ( int a = 0; a < num_bins; ++a ) {
    for ( int b = 0; b < num_bins; ++b ) {
      double element = unit_mat( a, b );
      double expected_element = 0.;
      if ( a == b ) expected_element = 1.;
      double abs_diff = std::abs( element - expected_element );
      if ( abs_diff > INVERSION_TOLERANCE ) throw std::runtime_error(
        "Matrix inversion failed in chi^2 calculation\n" );
    }
  }

  // Create a 1D vector containing the difference between the two slice
  // histograms in each bin
  TMatrixD diff_vec( 1, num_bins );
  for ( int a = 0; a < num_bins; ++a ) {
    // Note the one-based bin indices used for ROOT histograms
    double counts = hist_->GetBinContent( a + 1 );
    double other_counts = other.hist_->GetBinContent( a + 1 );
    diff_vec( 0, a ) = counts - other_counts;
  }

  // Multiply diff * covMat^{-1} * diff^{T} to get chi-squared
  TMatrixD temp1( diff_vec, TMatrixD::kMult, *inverse_cov_matrix );
  TMatrixD temp2( temp1, TMatrixD::kMult, diff_vec.T() );

  // We'll now have a 1x1 matrix containing the chi-squared value
  double chi2 = temp2( 0, 0 );

  // Assume that parameter fitting is not done, so that the relevant degrees
  // of freedom for the chi^2 test is just the number of bins minus one
  int dof = num_bins - 1;

  // Calculate a p-value for observing a chi^2 value at least as large as the
  // one actually obtained
  double p_value = ROOT::Math::inc_gamma_c( dof / 2., chi2 / 2. );

  Chi2Result result( chi2, num_bins, dof, p_value );
  return result;

}
