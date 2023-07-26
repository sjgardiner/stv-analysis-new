#pragma once

// Standard library includes
#include <memory>
#include <set>
#include <stdexcept>

// Simple container for the output of Unfolder::unfold()
struct UnfoldedMeasurement {
  UnfoldedMeasurement( TMatrixD* unfolded_signal, TMatrixD* cov_matrix,
    TMatrixD* unfolding_matrix, TMatrixD* err_prop_matrix,
    TMatrixD* add_smear_matrix ) : unfolded_signal_( unfolded_signal ),
    cov_matrix_( cov_matrix ), unfolding_matrix_( unfolding_matrix ),
    err_prop_matrix_( err_prop_matrix ), add_smear_matrix_( add_smear_matrix )
    {}

  std::unique_ptr< TMatrixD > unfolded_signal_;
  std::unique_ptr< TMatrixD > cov_matrix_;
  std::unique_ptr< TMatrixD > unfolding_matrix_;
  std::unique_ptr< TMatrixD > err_prop_matrix_;
  std::unique_ptr< TMatrixD > add_smear_matrix_;
};

// Container for mapping block indices to bin indices in
// Unfolder::blockwise_unfold()
struct BlockBins {
  BlockBins() {}
  std::vector< size_t > true_bin_indices_;
  std::vector< size_t > reco_bin_indices_;
};

// Abstract base class for objects that implement an algorithm for unfolding
// measured background-subtracted event counts from reco space to
// true space, possibly with regularization.
class Unfolder {

  public:

    Unfolder() {}

    // Function that actually implements a specific unfolding algorithm
    virtual UnfoldedMeasurement unfold( const TMatrixD& data_signal,
      const TMatrixD& data_covmat, const TMatrixD& smearcept,
      const TMatrixD& prior_true_signal ) const = 0;

    virtual UnfoldedMeasurement unfold( const TMatrixD& data_signal,
      const TMatrixD& data_covmat, const TMatrixD& smearcept,
      const TMatrixD& prior_true_signal,
      const std::string& block_defs_file ) const final;

    // Helper function that sets up unfolding for multiple blocks of bins,
    // then combines the results
    virtual UnfoldedMeasurement blockwise_unfold(
      const TMatrixD& data_signal, const TMatrixD& data_covmat,
      const TMatrixD& smearcept, const TMatrixD& prior_true_signal,
      std::ifstream& in_block_file ) const final;

  protected:

    // Helper function that does some sanity checks on the dimensions of the
    // input matrices passed to unfold()
    static void check_matrices( const TMatrixD& data_signal,
      const TMatrixD& data_covmat, const TMatrixD& smearcept,
      const TMatrixD& prior_true_signal );
};

UnfoldedMeasurement Unfolder::unfold( const TMatrixD& data_signal,
  const TMatrixD& data_covmat, const TMatrixD& smearcept,
  const TMatrixD& prior_true_signal,
  const std::string& block_defs_file ) const
{
  // Check the signal true bin definitions for the presence of multiple
  // block indices. Store all distinct values in a std::set. We will
  // assume here that the reco bin blocks are defined in a compatible way.
  // TODO: add error handling for bad block configurations
  int num_true_bins, tb_index, block_index;

  std::ifstream in_block_file( block_defs_file );
  in_block_file >> num_true_bins;

  std::set< int > true_blocks;
  for ( size_t tb = 0; tb < num_true_bins; ++tb ) {
    in_block_file >> tb_index >> block_index;
    true_blocks.insert( block_index );
  }

  // If multiple blocks are present, we need to unfold them individually
  // and then combine the results
  size_t num_blocks = true_blocks.size();
  if ( num_blocks > 1u ) {
    return this->blockwise_unfold( data_signal, data_covmat, smearcept,
      prior_true_signal, in_block_file );
  }

  // If there is only one block, we can just handle it directly
  return this->unfold( data_signal, data_covmat, smearcept,
    prior_true_signal );
}

void Unfolder::check_matrices( const TMatrixD& data_signal,
  const TMatrixD& data_covmat, const TMatrixD& smearcept,
  const TMatrixD& prior_true_signal )
{
  // Check the matrix dimensions for sanity
  int num_ordinary_reco_bins = smearcept.GetNrows();
  int num_true_signal_bins = smearcept.GetNcols();

  if ( data_signal.GetNcols() != 1 ) {
    throw std::runtime_error( "The background-subtracted data event counts"
      " must be expressed as a column vector" );
  }

  if ( prior_true_signal.GetNcols() != 1 ) {
    throw std::runtime_error( "The prior true signal event counts must be"
      " expressed as a column vector" );
  }

  if ( data_signal.GetNrows() != num_ordinary_reco_bins ) {
    throw std::runtime_error( "Reco bin mismatch between background-"
      "subtracted data and the smearceptance matrix" );
  }

  if ( data_covmat.GetNrows() != num_ordinary_reco_bins
    || data_covmat.GetNcols() != num_ordinary_reco_bins )
  {
    throw std::runtime_error( "Dimension mismatch between data covariance"
      " matrix and the smearceptance matrix" );
  }

  if ( prior_true_signal.GetNrows() != num_true_signal_bins ) {
    throw std::runtime_error( "Dimension mismatch between prior true signal"
      " event counts and the smearceptance matrix" );
  }

}

UnfoldedMeasurement Unfolder::blockwise_unfold( const TMatrixD& data_signal,
  const TMatrixD& data_covmat, const TMatrixD& smearcept,
  const TMatrixD& prior_true_signal, std::ifstream& in_block_file ) const
{
  // Rewind to the start of the block definition file
  in_block_file.seekg( 0 );

  int block_index;
  size_t num_reco_bins, rb_index, num_true_bins, tb_index;

  in_block_file >> num_true_bins;

  // Build a map of block indices to sets of true bin indices and reco bin
  // indices. This will be used below to extract each individual block from the
  // input matrices.
  std::map< int, BlockBins > block_map;
  for ( int tb = 0; tb < num_true_bins; ++tb ) {
    in_block_file >> tb_index >> block_index;

    auto& block_bins = block_map[ block_index ];
    block_bins.true_bin_indices_.push_back( tb_index );
  }

  in_block_file >> num_reco_bins;

  for ( int rb = 0; rb < num_reco_bins; ++rb ) {
    in_block_file >> rb_index >> block_index;

    auto& block_bins = block_map[ block_index ];
    block_bins.reco_bin_indices_.push_back( rb_index );
  }

  // TODO: add sanity checks of the block definitions

  // Create a single-column TMatrixD with the same number of true bins as the
  // prior. This will be used to combine the unfolded true bin counts from the
  // blocks to produce a final result.
  int num_true_signal_bins = prior_true_signal.GetNrows();
  TMatrixD* unfolded_signal = new TMatrixD( num_true_signal_bins, 1 );
  // Zero out the initial elements, just in case
  unfolded_signal->Zero();

  // Create a TMatrixD to hold the measurement error propagation matrix
  // aggregated across all blocks. This will be used to obtain the full
  // covariance matrix on the unfolded bin counts (including inter-block
  // covariances)
  int num_ordinary_reco_bins = data_signal.GetNrows();
  auto* err_prop = new TMatrixD( num_true_signal_bins, num_ordinary_reco_bins );
  // Zero out the initial elements, just in case
  err_prop->Zero();

  // Create a TMatrixD to hold the full unfolding matrix combined over
  // multiple blocks
  auto* unfold_mat = new TMatrixD( num_true_signal_bins,
    num_ordinary_reco_bins );
  unfold_mat->Zero();

  // Create a TMatrixD to hold the additional smearing matrix used to apply
  // regularization to theoretical predictions
  auto* add_smear = new TMatrixD( num_true_signal_bins, num_true_signal_bins );
  add_smear->Zero();

  // Loop over the blocks. For each block, populate the input matrices and
  // unfold.
  for ( const auto& block_pair : block_map ) {

    int b_idx = block_pair.first;
    const auto& block_bins = block_pair.second;

    std::cout << "Unfolding block " << b_idx << '\n';

    // Get the dimensions of the current block
    int num_block_true_bins = block_bins.true_bin_indices_.size();
    int num_block_reco_bins = block_bins.reco_bin_indices_.size();

    if ( num_block_true_bins < 1 ) throw std::runtime_error( "Block with zero"
      " true bins encountered" );
    if ( num_block_reco_bins < 1 ) throw std::runtime_error( "Block with zero"
      " reco bins encountered" );

    // Prepare matrices to store the block contents
    TMatrixD block_data_signal( num_block_reco_bins, 1 );
    TMatrixD block_data_covmat( num_block_reco_bins, num_block_reco_bins );
    TMatrixD block_smearcept( num_block_reco_bins, num_block_true_bins );
    TMatrixD block_prior_true_signal( num_block_true_bins, 1 );

    // Populate the matrices for the current block
    for ( int block_tb = 0; block_tb < num_block_true_bins; ++block_tb ) {
      // Convert the current true bin index at the block level to the
      // one at the global level
      int tb = block_bins.true_bin_indices_.at( block_tb );

      // Copy the prior true bin contents into the block
      block_prior_true_signal( block_tb, 0 ) = prior_true_signal( tb, 0 );

      for ( int block_rb = 0; block_rb < num_block_reco_bins; ++block_rb ) {
        // Convert the current reco bin index at the block level to the one
        // at the global level
        int rb = block_bins.reco_bin_indices_.at( block_rb );

        // Copy the smearceptance matrix contents into the block
        block_smearcept( block_rb, block_tb ) = smearcept( rb, tb );
      }
    }

    for ( int block_rb1 = 0; block_rb1 < num_block_reco_bins; ++block_rb1 ) {
      // Convert the current reco bin index at the block level to the one
      // at the global level
      int rb1 = block_bins.reco_bin_indices_.at( block_rb1 );

      // Copy the background-subtracted reco bin contents into the block
      block_data_signal( block_rb1, 0 ) = data_signal( rb1, 0 );

      for ( int block_rb2 = 0; block_rb2 < num_block_reco_bins; ++block_rb2 ) {
        // Convert the current reco bin index at the block level to the one
        // at the global level
        int rb2 = block_bins.reco_bin_indices_.at( block_rb2 );

        // Copy the reco-space covariance matrix element into the block
        block_data_covmat( block_rb1, block_rb2 ) = data_covmat( rb1, rb2 );
      }
    }

    // Unfold the measurement for the current block
    auto block_result = this->unfold( block_data_signal, block_data_covmat,
      block_smearcept, block_prior_true_signal );

    // Store the partial results for this block in the appropriate parts of the
    // matrices describing the full measurement
    for ( int block_tb = 0; block_tb < num_block_true_bins; ++block_tb ) {

      // Convert the current true bin index at the block level to the one
      // at the global level
      int tb = block_bins.true_bin_indices_.at( block_tb );

      // Copy the unfolded true bin contents from the current block
      unfolded_signal->operator()( tb, 0 )
        = block_result.unfolded_signal_->operator()( block_tb, 0 );

      for ( int block_rb = 0; block_rb < num_block_reco_bins; ++block_rb ) {

        // Convert the current reco bin index at the block level to the one
        // at the global level
        int rb = block_bins.reco_bin_indices_.at( block_rb );

        // Copy the measurement error propagation matrix element from the
        // current block
        err_prop->operator()( tb, rb ) = block_result.err_prop_matrix_
          ->operator()( block_tb, block_rb );

        // Copy the unfolding matrix element from the current block
        unfold_mat->operator()( tb, rb ) = block_result.unfolding_matrix_
          ->operator()( block_tb, block_rb );
      }

      for ( int block_tb2 = 0; block_tb2 < num_block_true_bins; ++block_tb2 ) {

        // Convert the current true bin index at the block level to the one
        // at the global level
        int tb2 = block_bins.true_bin_indices_.at( block_tb2 );

        // Copy the additional smearing matrix element from the current block
        add_smear->operator()( tb, tb2 ) = block_result.add_smear_matrix_
          ->operator()( block_tb, block_tb2 );
      }

    }

  } // block loop

  // All that remains is to propagate the full covariance matrix on the
  // measurement through the unfolding procedure. Do that transformation
  // using the measurement error propagation matrix built from all of the
  // blocks.
  TMatrixD err_prop_tr( TMatrixD::kTransposed, *err_prop );
  TMatrixD temp_mat( data_covmat, TMatrixD::EMatrixCreatorsOp2::kMult,
    err_prop_tr );

  auto* unfolded_signal_covmat = new TMatrixD( *err_prop,
    TMatrixD::EMatrixCreatorsOp2::kMult, temp_mat );

  UnfoldedMeasurement result( unfolded_signal, unfolded_signal_covmat,
    unfold_mat, err_prop, add_smear );
  return result;
}
