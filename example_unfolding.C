#include "TFile.h"
#include "TMatrixD.h"

#include "WienerSVDUnfolder.hh"
#include "DAgostiniUnfolder.hh"

using DCC = DAgostiniUnfolder::ConvergenceCriterion;
using RMT = WienerSVDUnfolder::RegularizationMatrixType;

const RMT MY_REGULARIZATION = RMT::kIdentity;
//const RMT MY_REGULARIZATION = RMT::kFirstDeriv;
//const RMT MY_REGULARIZATION = RMT::kSecondDeriv;

void example_unfolding() {

  // Retrieve the pre-calculated unfolding inputs for this example
  TFile* in_tfile = new TFile( "unf_input_dump.root", "read" );

  TMatrixD* data_signal = nullptr;
  TMatrixD* data_covmat = nullptr;
  TMatrixD* smearcept = nullptr;
  TMatrixD* prior_true_signal = nullptr;

  in_tfile->GetObject( "data_signal", data_signal );
  in_tfile->GetObject( "data_covmat", data_covmat );
  in_tfile->GetObject( "smearcept", smearcept );
  in_tfile->GetObject( "prior_true_signal", prior_true_signal );

  // Instantiate an object derived from the Unfolder base class
  WienerSVDUnfolder unfolder( true, MY_REGULARIZATION );

  //DAgostiniUnfolder unfolder( NUM_ITERATIONS );

  //DAgostiniUnfolder unfolder( DCC:FigureOfMerit, 0.025 );

  // Perform the unfolding
  UnfoldedMeasurement result = unfolder.unfold( *data_signal, *data_covmat,
    *smearcept, *prior_true_signal, "blocks.txt" );

  // Output is struct with these data members
  //std::unique_ptr< TMatrixD > unfolded_signal_;
  //std::unique_ptr< TMatrixD > cov_matrix_;
  //std::unique_ptr< TMatrixD > unfolding_matrix_;
  //std::unique_ptr< TMatrixD > err_prop_matrix_;
  //std::unique_ptr< TMatrixD > add_smear_matrix_;

}
