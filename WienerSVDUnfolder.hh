#pragma once

// ROOT includes
#include "TDecompChol.h"
#include "TDecompSVD.h"

// STV analysis includes
#include "MatrixUtils.hh"
#include "Unfolder.hh"

// Implementation of the Wiener-SVD unfolding method
// W. Tang et al., J. Instrum. 12, P10002 (2017)
// https://arxiv.org/abs/1705.03568
class WienerSVDUnfolder : public Unfolder {

  public:

    enum RegularizationMatrixType { kIdentity, kFirstDeriv, kSecondDeriv };

    WienerSVDUnfolder( bool use_wiener_filter = true,
      RegularizationMatrixType type = kIdentity ) : Unfolder(),
      use_filter_( use_wiener_filter ), reg_type_( type ) {}

    // Trick taken from https://stackoverflow.com/a/18100999
    using Unfolder::unfold;

    virtual UnfoldedMeasurement unfold( const TMatrixD& data_signal,
      const TMatrixD& data_covmat, const TMatrixD& smearcept,
      const TMatrixD& prior_true_signal ) const override;

    inline bool use_filter() const { return use_filter_; }
    inline void set_use_filter( bool use_filter ) { use_filter_ = use_filter; }

    inline RegularizationMatrixType get_regularization_type() const
      { return reg_type_; }

    inline void set_regularization_type( const RegularizationMatrixType& type )
      { reg_type_ = type; }

    inline const TMatrixD& additional_smearing_matrix() const
      { return *A_C_; }

  protected:

    // Helper function that sets the contents of the regularization matrix
    // based on the current value of reg_type_
    void set_reg_matrix( TMatrixD& C ) const;

    // Flag indicating whether the Wiener filter should be used. If it is
    // false, then the usual expression will be replaced with an identity
    // matrix
    bool use_filter_ = true;

    // Enum that determines the form to use for the regularization matrix C
    RegularizationMatrixType reg_type_ = kIdentity;

    // Additional smearing matrix from the last call to unfold()
    mutable std::unique_ptr< TMatrixD > A_C_;
};

UnfoldedMeasurement WienerSVDUnfolder::unfold( const TMatrixD& data_signal,
  const TMatrixD& data_covmat, const TMatrixD& smearcept,
  const TMatrixD& prior_true_signal ) const
{
  // Check input matrix dimensions for sanity
  this->check_matrices( data_signal, data_covmat,
    smearcept, prior_true_signal );

  int num_ordinary_reco_bins = smearcept.GetNrows();
  int num_true_signal_bins = smearcept.GetNcols();

  // Sec. 3.1 of the paper mentions that (in their notation) m >= n, i.e.,
  // the number of ordinary reco bins is assumed to be greater than or equal
  // to the number of true signal bins. We could run into trouble below
  // if this assumption is violated, so check it now.
  if ( num_ordinary_reco_bins < num_true_signal_bins ) {
    throw std::runtime_error( "The number of ordinary reco bins must exceed"
      " the number of true signal bins for Wiener-SVD unfolding." );
  }

  // Before doing anything fancy, first invert the input covariance matrix.
  // The utility function used here checks that the inversion was successful
  // and will complain if there's any trouble.
  auto inv_data_covmat = invert_matrix( data_covmat );

  // Perform a Cholesky decomposition of the inverted covariance matrix
  // A into A = U^T * U, where U is an upper-triangular matrix and U^T is
  // its transpose.
  TDecompChol chol( *inv_data_covmat );
  bool cholesky_ok = chol.Decompose();
  if ( !cholesky_ok ) throw std::runtime_error( "Cholesky decomposition failed"
    " during Wiener-SVD unfolding" );

  // The Wiener-SVD paper uses Q as the symbol for the lower-triangular matrix
  // Q = U^T (see Eq. (3.2) and note the opposite convention from ROOT's
  // TDecompChol class). Retrieve it from the Cholesky decomposition object and
  // transpose in place to get the correct definition.
  TMatrixD Q( chol.GetU() );

  // Apply the pre-scaling described in the paper below Eq. (3.3) using Q and
  // the smearceptance matrix. I use the same notation here as in the paper.
  // Note that I refrain from forming M = Q * data_signal (where data_signal is
  // the background-subtracted data) since it is actually more convenient below
  // to form R_tot and use it to transform both data_signal and the covariance
  // matrix.
  //TMatrixD M = Q * data_signal;
  TMatrixD R = Q * smearcept;

  // Also precompute the transpose of R for later convenience
  TMatrixD R_tr( TMatrixD::EMatrixCreatorsOp1::kTransposed, R );

  // Create the regularization matrix C mentioned in Eq. (3.19). Note that it
  // must always be a square matrix with dimension equal to the number of true
  // signal bins.
  TMatrixD C( num_true_signal_bins, num_true_signal_bins );
  this->set_reg_matrix( C );

  // Invert the regularization matrix to obtain C^(-1)
  auto Cinv = invert_matrix( C );

  // Prepare to perform the singular value decomposition (SVD) by multiplying
  // the inverted regularization matrix by the pre-scaled smearceptance matrix
  TMatrixD R_times_Cinv = R * ( *Cinv );

  // Perform a singular value decomposition of R * C^(-1) = U * S * V^T
  // where (switching from ROOT's notation to the notation of the paper)
  // U_C = U, D_C = S, and V_C = V.
  TDecompSVD svd( R_times_Cinv );
  bool svd_ok = svd.Decompose();
  if ( !svd_ok ) throw std::runtime_error( "Singular value decomposition"
    " failed during Wiener-SVD unfolding" );

  // Retrieve the SVD results for use in the unfolding calculation
  const TMatrixD& U_C = svd.GetU();
  const TVectorD& D_C_diag = svd.GetSig(); // diagonal elements of D_C only
  const TMatrixD& V_C = svd.GetV();

  // Precompute the transposes of U_C and V_C as well for later convenience
  TMatrixD U_C_tr( TMatrixD::EMatrixCreatorsOp1::kTransposed, U_C );
  TMatrixD V_C_tr( TMatrixD::EMatrixCreatorsOp1::kTransposed, V_C );

  // We can avoid one numerical matrix inversion by using the trick shown in
  // the Wiener-SVD source code, i.e.,
  // https://github.com/BNLIF/Wiener-SVD-Unfolding/blob/master/src/WienerSVD.C,
  // but not in the paper. To do that, we need to create the transpose of
  // the diagonal matrix D_C from the SVD results
  TMatrixD D_C_tr( num_true_signal_bins, num_ordinary_reco_bins );

  // Initialize all elements to zero, then set just the non-vanishing diagonal
  // elements
  D_C_tr.Zero();
  for ( int t = 0; t < num_true_signal_bins; ++t ) {
    D_C_tr( t, t ) = D_C_diag( t );
  }

  // Build the Wiener filter according to the expression in Eq. (3.24). Note
  // that it is a diagonal matrix.
  // TODO: Revisit whether there is a better way to do this, perhaps with
  // a TMatrixDSparse?
  TMatrixD W_C( num_true_signal_bins, num_true_signal_bins );

  // For simplicity, first calculate a column vector in which the
  // ith element corresponds to the ith value of the numerator in Eq. (3.24)
  // from the paper. We will then use this quantity (which also appears in
  // the denominator of that equation) to populate the Wiener filter matrix.
  TMatrixD numer_vec = V_C_tr * C * prior_true_signal;

  // Square each numerator vector element and multiply by the corresponding
  // squared diagonal element of D_C (see Eq. (3.24) from the paper).
  // NOTE: We don't have to worry about missing D_C_diag vector elements here
  // because the number of reco bins was already checked to equal or exceed the
  // number of true bins. The needed diagonal element of D_C will thus always
  // exist since the number of elements in numer_vec is equal to the number of
  // true bins.
  for ( int e = 0; e < num_true_signal_bins; ++e ) {
    double elem = numer_vec( e, 0 );
    double dC = D_C_diag( e );
    numer_vec( e, 0 ) = dC * dC * elem * elem;
  }

  // We can now fill in the diagonal elements of the Wiener filter matrix W_C
  // according to the expression in Eq. (3.24) from the paper
  for ( int t = 0; t < num_true_signal_bins; ++t ) {
    double numer = numer_vec( t, 0 );
    double denom = numer + 1;
    // Prevent division by zero by setting the numerator to zero
    if ( denom == 0 ) numer = 0.;
    W_C( t, t ) = numer / denom;
  }

  // If the Wiener filter has been disabled, then just replace W_C with
  // an identity matrix
  if ( !use_filter_ ) {
    W_C.UnitMatrix();
  }

  // Make a copy where we divide each diagonal element of W_C by the
  // corresponding diagonal element of (D_C^T * D_C)^(-1)
  TMatrixD W_C_tilde( W_C );
  for ( int t = 0; t < num_true_signal_bins; ++t ) {
    double dC = D_C_diag( t );
    W_C_tilde( t, t ) /= dC * dC;
  }

  // Calculate the additional smearing matrix A_C from Eq. (3.23) in the paper.
  // Matrix multiplication is associative, which is nice because we can chain
  // together a bunch of calls to operator*( const TMatrixD&, const TMatrixD& )
  // below safely.
  A_C_ = std::make_unique< TMatrixD >(
    ( *Cinv ) * V_C * W_C * V_C_tr * C
  );

  // Create the final unfolding matrix R_tot defined in Eq. (3.26) from the
  // paper. Avoid inverting (R^T * R) by using the trick from the Wiener-SVD
  // source code.
  TMatrixD R_tot = ( *Cinv ) * V_C * W_C_tilde * D_C_tr * U_C_tr * Q;

  // Get the unfolded signal event counts as a column vector
  auto* unfolded_signal = new TMatrixD( R_tot,
    TMatrixD::EMatrixCreatorsOp2::kMult, data_signal );

  // Get the covariance matrix on the unfolded signal
  TMatrixD R_tot_tr( TMatrixD::EMatrixCreatorsOp1::kTransposed, R_tot );
  TMatrixD temp_mat = data_covmat * R_tot_tr;

  auto* unfolded_signal_covmat = new TMatrixD( R_tot,
    TMatrixD::EMatrixCreatorsOp2::kMult, temp_mat );

  UnfoldedMeasurement result( unfolded_signal, unfolded_signal_covmat );
  return result;
}

void WienerSVDUnfolder::set_reg_matrix( TMatrixD& C ) const {
  // Zero out any existing matrix contents
  C.Zero();

  int num_rows = C.GetNrows();
  int num_cols = C.GetNcols();
  if ( num_rows != num_cols ) throw std::runtime_error( "Non-square"
    " regularization matrix passed to WienerSVDUnfolder::set_reg_matrix()" );

  // Fill in new contents based on the form requested by the user
  if ( reg_type_ == kIdentity ) {
    C.UnitMatrix();
  }
  else if ( reg_type_ == kFirstDeriv ) {
    for ( int r = 0; r < num_rows; ++r ) {
      C( r, r ) = -1.;
      if ( r < num_rows - 1 ) C( r, r + 1 ) = 1.;
    }
  }
  else if ( reg_type_ == kSecondDeriv ) {
    // Small number added to diagonal to avoid numerical inversion problems
    constexpr double inv_epsilon = 1e-6;
    for ( int r = 0; r < num_rows; ++r ) {
      if ( r < num_rows - 1 ) {
        C( r + 1, r ) = 1.;
        C( r, r + 1 ) = 1.;
      }
      double diag = inv_epsilon - 2.;
      if ( r == 0 || r == num_rows - 1 ) diag += 1.;
      C( r, r ) = diag;
    }
  }
  else throw std::runtime_error( "Unrecognized regularization matrix type"
    " encountered in WienerSVDUnfolder::set_reg_matrix()" );
}
