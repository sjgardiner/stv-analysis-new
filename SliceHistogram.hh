#pragma once

// Standard library includes
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

// ROOT includes
#include "TH1.h"
#include "Math/SpecFuncMathCore.h" // Needed for ROOT::Math::inc_gamma_c()

// STV analysis includes
#include "MatrixUtils.hh"
#include "NormShapeCovMatrix.hh"
#include "SliceBinning.hh"
#include "SystematicsCalculator.hh"

class SliceHistogram {

  public:

    SliceHistogram() {}

    static SliceHistogram* make_slice_histogram( TH1D& reco_bin_histogram,
      const Slice& slice, const CovMatrix* input_cov_mat = nullptr );

    static SliceHistogram* make_slice_histogram( TMatrixD& reco_bin_counts,
      const Slice& slice, const TMatrixD* input_cov_mat );

    // TODO: revisit this implementation
    static SliceHistogram* make_slice_efficiency_histogram(
      const TH1D& true_bin_histogram, const TH2D& hist_2d, const Slice& slice );

    // Transform the bin contents by multiplying by the input TMatrixD, which
    // must be a square matrix with a number of columns equal to the number of
    // histogram bins. If present, the owned covariance matrix will also be
    // transformed accordingly.
    void transform( const TMatrixD& mat );

    // Create a column vector with the current histogram bin contents
    TMatrixD get_col_vect() const;

    // Calculates a decomposition of the full covariance matrix into norm,
    // mixed, and shape-only pieces
    void calc_norm_shape_errors();

    struct Chi2Result {

      Chi2Result() {}

      Chi2Result( double chi2, int nbins, int dof, double pval )
        : chi2_( chi2 ), num_bins_( nbins ), dof_( dof ), p_value_( pval ) {}

      double chi2_;
      int num_bins_;
      int dof_;
      double p_value_;
    };

    Chi2Result get_chi2( const SliceHistogram& other,
      const double inversion_tol = DEFAULT_MATRIX_INVERSION_TOLERANCE ) const;

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
      // The UniverseMaker reco bin indices are zero-based, so I correct
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
            // the UniverseMaker numbering scheme is zero-based. I
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

SliceHistogram* SliceHistogram::make_slice_histogram(
  TMatrixD& reco_bin_counts, const Slice& slice,
  const TMatrixD* input_cov_mat )
{
  // TODO: reduce code duplication between this function and the overloaded
  // version that takes an input TH1& and CovMatrix*

  // Check that the reco_bin_counts are given as a column vector
  if ( reco_bin_counts.GetNcols() != 1 ) {
    throw std::runtime_error( "Invalid dimension for bin counts passed"
      "to SliceHistogram::make_slice_histogram()" );
  }

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
      // The UniverseMaker reco bin indices are zero-based like the
      // TMatrixD element indices
      slice_bin_content += reco_bin_counts( rb_idx, 0 );
    }

    slice_hist->SetBinContent( slice_bin_idx, slice_bin_content );

  } // slice bins

  // If we've been handed a non-null pointer to a TMatrixD object representing
  // the covariance matrix, then we will use it to propagate uncertainties.
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
        for ( const auto& rb_m : rb_set_a ) {
          for ( const auto& rb_n : rb_set_b ) {
            // The TMatrixD object uses zero-based indices
            cov += input_cov_mat->operator()( rb_m, rb_n );
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

// TODO: revisit this rough draft. Right now, an assumption is made that the
// true and reco bins are defined in the same way with the same indices. This
// isn't enforced by the UniverseMaker configuration itself, although
// it is currently consistent with what you've done so far.
SliceHistogram* SliceHistogram::make_slice_efficiency_histogram(
  const TH1D& true_bin_histogram, const TH2D& hist_2d, const Slice& slice )
{
  // Get the binning and axis labels for the current slice by cloning the
  // (empty) histogram owned by the Slice object
  TH1* slice_hist = dynamic_cast< TH1* >(
    slice.hist_->Clone("slice_hist") );

  slice_hist->SetDirectory( nullptr );
  slice_hist->GetYaxis()->SetTitle( "efficiency" );
  slice_hist->GetYaxis()->SetRangeUser( 0., 1. );

  // Fill the slice bins based on the input reco bins
  for ( const auto& pair : slice.bin_map_ ) {

    // One-based index for the global TH1 bin number in the slice
    int slice_bin_idx = pair.first;

    const auto& reco_bin_set = pair.second;

    double selected_signal_evts = 0.;
    double all_signal_evts = 0.;
    for ( const auto& rb_idx : reco_bin_set ) {
      // The UniverseMaker reco bin indices are zero-based, so I correct
      // for this here when pulling values from the one-based input ROOT
      // histogram.
      all_signal_evts += true_bin_histogram.GetBinContent( rb_idx + 1 );

      // Include selected signal events in the current true bin that fall into
      // any of the reco bins
      selected_signal_evts += hist_2d.Integral( rb_idx + 1, rb_idx + 1,
        1, hist_2d.GetNbinsY() );
    }

    double bin_efficiency = selected_signal_evts / all_signal_evts;
    // See DocDB #32401, Eq. (5.2)
    double bin_stat_err = std::sqrt( std::max(0., bin_efficiency
      * (1. - bin_efficiency) / all_signal_evts) );
    slice_hist->SetBinContent( slice_bin_idx, bin_efficiency );
    slice_hist->SetBinError( slice_bin_idx, bin_stat_err );

  } // slice bins

  TH2D* covmat_hist = nullptr;

  // We're done. Prepare the SliceHistogram object and return it.
  auto* result = new SliceHistogram;
  result->hist_.reset( slice_hist );
  result->cmat_.cov_matrix_.reset( nullptr );

  return result;
}

SliceHistogram::Chi2Result SliceHistogram::get_chi2(
  const SliceHistogram& other, const double inversion_tol ) const
{
  int num_bins = hist_->GetNbinsX();
  if ( other.hist_->GetNbinsX() != num_bins ) {
    throw std::runtime_error( "Incompatible vector sizes in chi^2"
      " calculation" );
  }

  // If both SliceHistogram objects have a covariance matrix, then
  // check that their dimensions match. If one is missing, it will be assumed
  // to be a null matrix
  if ( cmat_.cov_matrix_ && other.cmat_.cov_matrix_ ) {
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
  }
  else if ( !cmat_.cov_matrix_ && !other.cmat_.cov_matrix_ ) {
    throw std::runtime_error( "Both SliceHistogram objects involved in"
      " a chi^2 calculation have null covariance matrices" );
  }

  // The total covariance matrix on the difference between the
  // two histograms is just the sum of each individual SliceHistogram's
  // owned covariance matrix.
  CovMatrix cov_mat;
  cov_mat += cmat_;
  cov_mat += other.cmat_;

  // Get access to a TMatrixD object representing the covariance matrix.
  auto cov_matrix = cov_mat.get_matrix();

  // Invert the covariance matrix
  auto inverse_cov_matrix = invert_matrix( *cov_matrix, inversion_tol );

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

void SliceHistogram::transform( const TMatrixD& mat ) {

  int dim = hist_->GetDimension();
  if ( dim != 1 ) throw std::runtime_error( "SliceHistogram::transform() is"
    " currently implemented only for 1D histograms." );

  int num_cols = mat.GetNcols();
  int num_bins = hist_->GetNbinsX();
  if ( num_cols != num_bins ) throw std::runtime_error( "Incompatible"
    " transformation matrix passed to SliceHistogram::transform()" );

  int num_rows = mat.GetNrows();
  if ( num_rows != num_cols ) throw std::runtime_error( "Transformations which"
    " change the number of bins are currently unimplemented in"
    " SliceHistogram::transform()" );

  // Create a column vector with the current histogram bin contents
  TMatrixD hist_vec = this->get_col_vect();

  // Apply the transformation matrix to the histogram and store the result in a
  // new column vector
  TMatrixD transformed_hist_vec( mat,
    TMatrixD::EMatrixCreatorsOp2::kMult, hist_vec );

  // Replace the old histogram contents with the new ones
  for ( int b = 0; b < num_bins; ++b ) {
    double val = transformed_hist_vec( b, 0 );
    hist_->SetBinContent( b + 1, val );
  }

  // If the covariance matrix isn't defined, then we're done and can return
  // early. Otherwise, we'll apply a corresponding transformation to the
  // covariance matrix.
  if ( !cmat_.cov_matrix_ ) return;

  // Get the original covariance matrix as a std::unique_ptr< TMatrixD >
  auto orig_cov = cmat_.get_matrix();

  // Take the transpose of the transformation matrix
  TMatrixD tr_mat( TMatrixD::kTransposed, mat );

  // See https://stats.stackexchange.com/q/113700
  TMatrixD transformed_cov = mat * ( *orig_cov ) * tr_mat;

  // Create a new CovMatrix object using the transformed covariance matrix
  CovMatrix transformed_cmat( transformed_cov );

  // Replace the owned CovMatrix object with the new one
  cmat_ = std::move( transformed_cmat );

  // To wrap things up, set the updated histogram bin errors based on the
  // diagonal elements of the covariance matrix
  for ( int b = 0; b < num_bins; ++b ) {
    double variance = cmat_.cov_matrix_->GetBinContent( b + 1, b + 1 );
    double err = std::sqrt( std::max(0., variance) );
    //double err = shape_errors_.at( b );
    hist_->SetBinError( b + 1, err );
  }

}

// Create a column vector with the current histogram bin contents
TMatrixD SliceHistogram::get_col_vect() const {
  int num_bins = hist_->GetNbinsX();
  TMatrixD hist_vec( num_bins, 1 );
  for ( int b = 0; b < num_bins; ++b ) {
    // Note that TH1D bin indices are one based while TMatrixD element indices
    // are zero-based
    double val = hist_->GetBinContent( b + 1 );
    hist_vec( b, 0 ) = val;
  }
  return hist_vec;
}
