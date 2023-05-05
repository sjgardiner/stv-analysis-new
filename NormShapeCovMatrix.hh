#pragma once

// Standard library includes
#include <stdexcept>

// ROOT includes
#include "TMatrixD.h"

// STV analysis includes
#include "UniverseMaker.hh"

// Splits the contributions from an ordinary covariance matrix into
// normalization, shape, and mixed pieces. See
// https://microboone-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=5926. Code
// taken from Afro with minor adjustments and sincere gratitude! ♥
class NormShapeCovMatrix {
  public:
    NormShapeCovMatrix() {}
    NormShapeCovMatrix( const TMatrixD& col_matrix_pred, const TMatrixD& cov_matrix );

    // Covariance matrix components
    TMatrixD mixed_;
    TMatrixD norm_;
    TMatrixD shape_;

    // For convenience, also store the mixed + shape piece in one matrix
    TMatrixD mixed_plus_shape_;
};

NormShapeCovMatrix::NormShapeCovMatrix( const TMatrixD& col_matrix_pred,
  const TMatrixD& cov_matrix )
{
  int n_bins = col_matrix_pred.GetNrows();
  int should_be_one = col_matrix_pred.GetNcols();

  int n_rows = cov_matrix.GetNrows();
  int n_cols = cov_matrix.GetNcols();

  if ( n_bins != n_rows || n_bins != n_cols || should_be_one != 1 ) {
    throw std::runtime_error( "Invalid matrix dimensions encountered in"
      " constructor of NormShapeCovMatrix" );
  }

  shape_.ResizeTo( n_bins, n_bins );
  mixed_.ResizeTo( n_bins, n_bins );
  norm_.ResizeTo( n_bins, n_bins );
  mixed_plus_shape_.ResizeTo( n_bins, n_bins );

  double N_T = 0.;
  for ( int idx = 0; idx < n_bins; ++idx ) {
    N_T += col_matrix_pred( idx, 0 );
  }

  double M_kl = 0;
  for ( int i = 0; i < n_bins; ++i ) {
    for ( int j = 0; j < n_bins; ++j ) {
      M_kl += cov_matrix( i, j );
    }
  }

  for ( int i = 0; i < n_bins; ++i) {
    for ( int j = 0; j < n_bins; ++j ) {
      double N_i = col_matrix_pred( i, 0 );
      double N_j = col_matrix_pred( j, 0 );
      double M_ij = cov_matrix( i, j );
      double M_ik = 0.;
      for ( int k = 0; k < n_bins; ++k ) M_ik += cov_matrix( i, k );
      double M_kj = 0.;
      for( int k = 0; k < n_bins; ++k ) M_kj += cov_matrix( k, j );

      shape_( i, j ) = M_ij - N_j*M_ik/N_T - N_i*M_kj/N_T + N_i*N_j*M_kl/N_T/N_T;
      mixed_( i, j ) = N_j*M_ik/N_T + N_i*M_kj/N_T - 2.*N_i*N_j*M_kl/N_T/N_T;
      norm_( i, j ) = N_i*N_j*M_kl/N_T/N_T;
    }
  }

  mixed_plus_shape_ = shape_ + mixed_;
}

NormShapeCovMatrix make_block_diagonal_norm_shape_covmat(
  const TMatrixD& unfolded_signal, const TMatrixD& unfolded_covmat,
  const std::vector< TrueBin >& true_bins )
{
  // Build a map of block indices to sets of signal true bin indices. This will
  // be used below to extract each individual block.
  std::map< int, std::vector<size_t> > block_map;
  int num_true_signal_bins = 0;
  for ( size_t tb = 0u; tb < true_bins.size(); ++tb ) {
    const auto& tbin = true_bins.at( tb );
    if ( tbin.type_ == TrueBinType::kSignalTrueBin ) {

      ++num_true_signal_bins;

      int b_idx = tbin.block_index_;

      auto& block_bins = block_map[ b_idx ];
      block_bins.push_back( tb );
    }
  }

  // TODO: add sanity checks of the block definitions

  // Create a NormShapeCovMatrix object that we will fill manually with the full
  // results for blockwise diagonal components of the decomposed covariance
  // matrix
  NormShapeCovMatrix bw_ns_cm;
  bw_ns_cm.norm_.ResizeTo( num_true_signal_bins, num_true_signal_bins );
  bw_ns_cm.shape_.ResizeTo( num_true_signal_bins, num_true_signal_bins );
  bw_ns_cm.mixed_.ResizeTo( num_true_signal_bins, num_true_signal_bins );
  bw_ns_cm.mixed_plus_shape_.ResizeTo( num_true_signal_bins,
    num_true_signal_bins );

  bw_ns_cm.norm_.Zero();
  bw_ns_cm.shape_.Zero();
  bw_ns_cm.mixed_.Zero();
  bw_ns_cm.mixed_plus_shape_.Zero();

  // Loop over the blocks. For each block, populate the input matrices and
  // unfold.
  for ( const auto& block_pair : block_map ) {

    int b_idx = block_pair.first;
    const auto& block_bins = block_pair.second;

    std::cout << "Norm/shape decomposition for block " << b_idx << '\n';

    // Get the dimensions of the current block
    int num_block_true_bins = block_bins.size();

    if ( num_block_true_bins < 1 ) throw std::runtime_error( "Block with zero"
      " true bins encountered" );

    // Prepare matrices to store the block contents
    TMatrixD block_signal( num_block_true_bins, 1 );
    TMatrixD block_cm( num_block_true_bins, num_block_true_bins );

    // Populate the matrices for the current block
    for ( int block_tb1 = 0; block_tb1 < num_block_true_bins; ++block_tb1 ) {
      // Convert the current true bin index at the block level to the
      // one at the global level
      int tb1 = block_bins.at( block_tb1 );

      // Copy the prior true bin contents into the block
      block_signal( block_tb1, 0 ) = unfolded_signal( tb1, 0 );

      for ( int block_tb2 = 0; block_tb2 < num_block_true_bins; ++block_tb2 ) {

        // Convert the current true bin index at the block level to the one
        // at the global level
        int tb2 = block_bins.at( block_tb2 );

        // Copy the unfolded covariance matrix element into the block
        block_cm( block_tb1, block_tb2 ) = unfolded_covmat( tb1, tb2 );
      }
    }

    // Decompose the covariance matrix for the current block into normalization,
    // shape, and mixed pieces
    NormShapeCovMatrix block_decomp( block_signal, block_cm );

    // Store the partial results for this block in the appropriate parts of the
    // matrices describing the full measurement
    for ( int block_tb1 = 0; block_tb1 < num_block_true_bins; ++block_tb1 ) {

      // Convert the current true bin index at the block level to the one
      // at the global level
      int tb1 = block_bins.at( block_tb1 );

      for ( int block_tb2 = 0; block_tb2 < num_block_true_bins; ++block_tb2 ) {

        // Convert the current true bin index at the block level to the one
        // at the global level
        int tb2 = block_bins.at( block_tb2 );

        // Copy the decomposed covariance matrix elements from the current block
        bw_ns_cm.norm_( tb1, tb2 ) = block_decomp.norm_( block_tb1, block_tb2 );

        bw_ns_cm.shape_( tb1, tb2 )
          = block_decomp.shape_( block_tb1, block_tb2 );

        bw_ns_cm.mixed_( tb1, tb2 )
          = block_decomp.mixed_( block_tb1, block_tb2 );
      }

    }

  } // block loop

  bw_ns_cm.mixed_plus_shape_ = bw_ns_cm.shape_ + bw_ns_cm.mixed_;

  return bw_ns_cm;
}
