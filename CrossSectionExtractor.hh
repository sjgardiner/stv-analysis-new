#pragma once

// Standard library includes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

// ROOT includes
#include "TCanvas.h"
#include "TLegend.h"

// STV analysis includes
#include "DAgostiniUnfolder.hh"
#include "FiducialVolume.hh"
#include "MatrixUtils.hh"
#include "MCC9SystematicsCalculator.hh"
#include "NormShapeCovMatrix.hh"
#include "PGFPlotsDumpUtils.hh"
#include "SliceBinning.hh"
#include "SliceHistogram.hh"
#include "SystematicsCalculator.hh"
#include "WienerSVDUnfolder.hh"

using WSVD_RMT = WienerSVDUnfolder::RegularizationMatrixType;
using MCC9SystMode = MCC9SystematicsCalculator::SystMode;

constexpr double BIG_DOUBLE = 1e300;
constexpr bool USE_ADD_SMEAR = true;
constexpr bool INCLUDE_BKGD_ONLY_ERRORS = false;
constexpr bool INCLUDE_SIGRESP_ONLY_ERRORS = false;

// Base class for representing predictions of the expected number of events
// in each signal true bin. These can come from multiple sources (e.g.,
// one of the universes or an external calculation).
class PredictedTrueEvents {

  public:

    PredictedTrueEvents( int num_ts_bins, const std::string& name )
      : num_true_signal_bins_( num_ts_bins ), name_( name ),
      pred_( num_ts_bins, 1 ) // Column vector to store predicted event counts
    {
      // The name (used for dumping to pgfplots, etc.) cannot contain
      // a space, tab, or newline character. If it does, then complain.
      for ( size_t c = 0u; c < name.size(); ++c ) {
        char ch = name.at( c );
        if ( ch == ' ' || ch == '\t' || ch == '\n' ) {
          throw std::runtime_error( "Invalid name in constructor of"
            " PredictedTrueEvents" );
        }
      }
    }

    const std::string& name() { return name_; }
    virtual const TMatrixD& prediction() const { return pred_; }
    virtual TMatrixD& get_prediction() { return pred_; }

  protected:

    int num_true_signal_bins_;
    std::string name_;
    TMatrixD pred_;
};

// Predicted true signal event counts from a Universe object
class UniverseTrueEvents : public PredictedTrueEvents {

  public:

    UniverseTrueEvents( int num_ts_bins, const std::string& name,
      const Universe& univ ) : PredictedTrueEvents( num_ts_bins, name )
    {
      // Copy the true event counts from the input Universe object
      for ( int b = 0; b < num_true_signal_bins_; ++b ) {
        double true_evts = univ.hist_true_->GetBinContent( b + 1 );
        pred_( b, 0 ) = true_evts;
      }
    }

};

// Predicted true signal event counts from a ROOT histogram of binwise total
// cross sections (cm^2 / Ar) stored in a file
class FileTrueEvents : public PredictedTrueEvents {
  public:

    FileTrueEvents( int num_ts_bins, const std::string& name,
      const std::string& file_name, const std::string& hist_name,
      double conv_factor ) : PredictedTrueEvents( num_ts_bins, name )
    {
      // Retrieve the raw prediction histogram. This class expects it to be
      // expressed as a total cross section with true bin number along the
      // x-axis
      TFile temp_in_file( file_name.c_str(), "read" );
      TH1D* temp_hist = nullptr;
      temp_in_file.GetObject( hist_name.c_str(), temp_hist );

      if ( !temp_hist ) {
        throw std::runtime_error( "Could not retrieve the histogram \""
          + hist_name + "\" from the file \"" + file_name + '\"' );
      }

      // Convert the content (and error) of each bin to an expected true event
      // count. Do this using the input conversion factor (expected to be
      // calculated as integrated numu flux * number of Ar targets in the
      // fiducial volume)
      int num_bins = temp_hist->GetNbinsX();
      if ( num_bins != num_true_signal_bins_ ) {
        throw std::runtime_error( "Prediction histogram bin count mismatch" );
      }

      for ( int b = 0; b < num_bins; ++b ) {
        double xsec = temp_hist->GetBinContent( b + 1 );
        pred_( b, 0 ) = xsec * conv_factor;
      }

    }

};

// Predicted true signal event counts from an input TH1D
class HistogramTrueEvents : public PredictedTrueEvents {

  public:

    HistogramTrueEvents( int num_ts_bins, const std::string& name,
      const TH1D& hist ) : PredictedTrueEvents( num_ts_bins, name )
    {
      // Check that the input histogram has the correct number of bins
      if ( hist.GetNbinsX() != num_ts_bins ) {
        throw std::runtime_error( "Bin count mismatch in constructor of"
          " HistogramTrueEvents" );
      }

      // Copy the event count values from the input histogram
      for ( int b = 0; b < num_true_signal_bins_; ++b ) {
        double true_evts = hist.GetBinContent( b + 1 );
        pred_( b, 0 ) = true_evts;
      }

    }

};

struct CrossSectionResult {
  CrossSectionResult( UnfoldedMeasurement& meas )
    : result_( std::move(meas) ) {}
  UnfoldedMeasurement result_;
  std::map< std::string, std::unique_ptr<TMatrixD> > unfolded_cov_matrix_map_;
};

class CrossSectionExtractor {

  public:

    CrossSectionExtractor( const std::string& config_file_name );

    // Get the multiplicative factor needed to convert from a total cross
    // section histogram (as in FileTrueEvents) to a true event count histogram
    double conversion_factor() const;

    // Apply background subtraction and the unfolding procedure to obtain
    // a measurement together with full uncertainties
    CrossSectionResult get_unfolded_events();

  protected:

    void prepare_predictions( const std::vector< std::string >& line_vec );

    // Keys are model descriptions for the plot legend (white space allowed),
    // values are PredictedTrueEvents objects
    std::map< std::string, std::unique_ptr< PredictedTrueEvents > > pred_map_;

    // Unfolder object used to correct cross-section measurements for
    // inefficiency and bin migrations
    std::unique_ptr< Unfolder > unfolder_;

    // Systematics calculator object used to compute covariance matrices
    std::unique_ptr< SystematicsCalculator > syst_;
};

CrossSectionExtractor::CrossSectionExtractor(
  const std::string& config_file_name )
{
  std::ifstream config_file( config_file_name );
  if ( !config_file.good() ) {
    throw std::runtime_error( "Could not read CrossSectionExtractor"
      " configuration from the file \"" + config_file_name + '\"' );
  }

  // Strings that will be initialized using the contents of the
  // CrossSectionExtractor configuration file
  std::string syst_config_file_name;
  std::string univ_file_name;

  // Temporary storage for the configuration file lines defining each
  // prediction. We will revisit these once the systematic Universe objects
  // are fully initialized.
  std::vector< std::string > pred_line_vec;

  std::string line;
  while ( std::getline(config_file, line) ) {
    // Skip lines starting with a '#' character
    if ( line.front() == '#' ) continue;

    // Prepare a stringstream to parse the line contents
    std::istringstream iss( line );
    std::string first_word;
    iss >> first_word;

    if ( first_word == "Prediction" ) {
      // Save each line defining a prediction for later
      pred_line_vec.push_back( line );
    }
    else if ( first_word == "FPFile" ) {
      // Read in the non-default setting for the name of the configuration
      // file for the FilePropertiesManager
      std::string file_properties_config;
      iss >> file_properties_config;

      // Reinitialize the FilePropertiesManager using the configuration file
      auto& fpm = FilePropertiesManager::Instance();
      fpm.load_file_properties( file_properties_config );
    }
    else if ( first_word == "SystFile" ) {
      // Get the name of the non-default SystematicsCalculator
      // configuration file
      iss >> syst_config_file_name;
    }
    else if ( first_word == "UnivFile" ) {
      // Get the name of the ROOT file containing the pre-calculated
      // Universe histograms
      iss >> univ_file_name;
    }
    else if ( first_word == "Unfold" ) {
      // Get the string indicating which unfolding method should be used
      std::string unf_type;
      iss >> unf_type;

      // Construct the appropriate Unfolder object
      Unfolder* temp_unfolder = nullptr;
      if ( unf_type == "DAgostini" ) {
        // Determine how to configure the D'Agostini unfolding algorithm
        std::string dagost_mode;
        iss >> dagost_mode;

        if ( dagost_mode == "fm" ) {
          double fig_of_merit;
          iss >> fig_of_merit;
          temp_unfolder = new DAgostiniUnfolder( DAgostiniUnfolder
            ::ConvergenceCriterion::FigureOfMerit, fig_of_merit );
        }
        else if ( dagost_mode == "iter" ) {
          int iterations;
          iss >> iterations;
          temp_unfolder = new DAgostiniUnfolder( iterations );
        }
        else {
          throw std::runtime_error( "Unrecognized D'Agostini unfolding"
            " mode \"" + dagost_mode + '\"' );
        }
      }
      else if ( unf_type == "WienerSVD" ) {
        bool use_filter;
        iss >> use_filter;

        // Default to regularizing using the second derivative. The choice
        // actually doesn't matter unless the Wiener filter is used.
        WSVD_RMT reg_type = WSVD_RMT::kSecondDeriv;

        if ( use_filter ) {
          // Determine the regularization recipe to use with the WSVD
          // method
          std::string reg_mode;
          iss >> reg_mode;

          if ( reg_mode == "identity" ) {
            reg_type = WSVD_RMT::kIdentity;
          }
          else if ( reg_mode == "first-deriv" ) {
            reg_type = WSVD_RMT::kFirstDeriv;
          }
          else if ( reg_mode == "second-deriv" ) {
            reg_type = WSVD_RMT::kSecondDeriv;
          }
          else {
            throw std::runtime_error( "Unrecognized Wiener-SVD"
              " regularization mode \"" + reg_mode + '\"' );
          }
        }

        temp_unfolder = new WienerSVDUnfolder( use_filter, reg_type );
      }
      else {
        throw std::runtime_error( "Unrecognized unfolder type \""
          + unf_type + '\"' );
      }

      // Store the complete Unfolder object for later use
      unfolder_.reset( temp_unfolder );
    }
    else {
      throw std::runtime_error( "Unrecognized CrossSectionExtractor"
        " configuration file command \"" + first_word + '\"' );
    }
  }

  // We've finished parsing the configuration file. Check that we have the
  // required information.
  if ( !unfolder_ ) {
    throw std::runtime_error( "Missing \"Unfold\" command in the"
      " CrossSectionExtractor configuration file" );
  }
  else if ( univ_file_name.empty() ) {
    throw std::runtime_error( "Missing universe file specification in the"
      " CrossSectionExtractor configuration file" );
  }

  // Initialize the owned SystematicsCalculator
  auto* temp_syst = new MCC9SystematicsCalculator( univ_file_name,
    syst_config_file_name );
  syst_.reset( temp_syst );

  // With the SystematicsCalculator in place (including its owned Universe
  // objects), we are ready to initialize the predictions
  this->prepare_predictions( pred_line_vec );
}

CrossSectionResult CrossSectionExtractor::get_unfolded_events() {

  // Evaluate the total and partial covariance matrices in reco space
  auto matrix_map = syst_->get_covariances();

  auto* mcc9 = dynamic_cast< MCC9SystematicsCalculator* >( syst_.get() );

  if ( INCLUDE_BKGD_ONLY_ERRORS ) {

    // Add in covariances calculated by varying only the background MC
    // prediction
    mcc9->set_syst_mode( MCC9SystMode::VaryOnlyBackground );
    auto bkgd_matrix_map = syst_->get_covariances();

    for ( const auto& m_pair : *bkgd_matrix_map ) {
      auto& my_temp_cov_mat = matrix_map->operator[](
        "bkgd_only_" + m_pair.first );
      my_temp_cov_mat += m_pair.second;
    }

    // Restore the default behavior of the SystematicsCalculator
    mcc9->set_syst_mode( MCC9SystMode::ForXSec );
  }

  if ( INCLUDE_SIGRESP_ONLY_ERRORS ) {

    // Add in covariances calculated by varying only the signal response
    // estimated using the CV MC prediction
    mcc9->set_syst_mode( MCC9SystMode::VaryOnlySignalResponse );
    auto sigresp_matrix_map = syst_->get_covariances();

    for ( const auto& m_pair : *sigresp_matrix_map ) {
      auto& my_temp_cov_mat = matrix_map->operator[](
        "sigresp_only_" + m_pair.first );
      my_temp_cov_mat += m_pair.second;
    }

    // Restore the default behavior of the SystematicsCalculator
    mcc9->set_syst_mode( MCC9SystMode::ForXSec );
  }

  // Perform background subtraction and unfolding to get a measurement of event
  // counts in (regularized) true space
  UnfoldedMeasurement result = unfolder_->unfold( *syst_ );
  CrossSectionResult xsec( result );

  if ( USE_ADD_SMEAR ) {

    // Get access to the additional smearing matrix
    const TMatrixD& A_C = *xsec.result_.add_smear_matrix_;

    // Update each of the owned predictions by multiplying them by the
    // additional smearing matrix
    for ( auto& pair : pred_map_ ) {
      //const auto& model_description = pair.first;
      TMatrixD& truth_pred = pair.second->get_prediction();

      TMatrixD ac_temp( A_C, TMatrixD::kMult, truth_pred );
      truth_pred = ac_temp;
    }
  }

  // Propagate all defined covariance matrices through the unfolding procedure
  // using the "error propagation matrix" and its transpose
  const TMatrixD& err_prop = *xsec.result_.err_prop_matrix_;
  TMatrixD err_prop_tr( TMatrixD::kTransposed, err_prop );

  for ( const auto& matrix_pair : *matrix_map ) {
    const std::string& matrix_key = matrix_pair.first;
    auto temp_cov_mat = matrix_pair.second.get_matrix();

    TMatrixD temp_mat( *temp_cov_mat, TMatrixD::EMatrixCreatorsOp2::kMult,
      err_prop_tr );

    xsec.unfolded_cov_matrix_map_[ matrix_key ] = std::make_unique< TMatrixD >(
      err_prop, TMatrixD::EMatrixCreatorsOp2::kMult, temp_mat );
  }

  // Decompose the block-diagonal pieces of the total covariance matrix
  // into normalization, shape, and mixed components (for later plotting
  // purposes)
  NormShapeCovMatrix bd_ns_covmat = make_block_diagonal_norm_shape_covmat(
    *xsec.result_.unfolded_signal_, *xsec.result_.cov_matrix_,
    syst_->true_bins_ );

  // Add the blockwise decomposed matrices into the map
  xsec.unfolded_cov_matrix_map_[ "total_blockwise_norm" ]
    = std::make_unique< TMatrixD >( bd_ns_covmat.norm_ );

  xsec.unfolded_cov_matrix_map_[ "total_blockwise_shape" ]
    = std::make_unique< TMatrixD >( bd_ns_covmat.shape_ );

  xsec.unfolded_cov_matrix_map_[ "total_blockwise_mixed" ]
    = std::make_unique< TMatrixD >( bd_ns_covmat.mixed_ );

  return xsec;
}

double CrossSectionExtractor::conversion_factor() const {
  double total_pot = syst_->total_bnb_data_pot_;
  double integ_flux = integrated_numu_flux_in_FV( total_pot );
  double num_Ar = num_Ar_targets_in_FV();

  double conv_factor = num_Ar * integ_flux;
  return conv_factor;
}

void CrossSectionExtractor::prepare_predictions(
  const std::vector< std::string >& line_vec )
{
  // Just in case, check that we have set up the owned SystematicsCalculator
  if ( !syst_ ) {
    throw std::runtime_error( "Unitializated SystematicsCalculator"
      " encountered in CrossSectionExtractor::prepare_predictions()" );
  }

  // Only the signal true bins really matter for the predictions managed by
  // this function
  int num_bins = syst_->get_num_signal_true_bins();

  double conv_factor = this->conversion_factor();

  for ( const auto& line : line_vec ) {

    // Create a new prediction for each line to be processed
    PredictedTrueEvents* pred = nullptr;

    // Skip the first word on the line (already known to be "Prediction"
    // based on prior processing)
    std::istringstream iss( line );
    std::string word;
    iss >> word;

    std::string name;
    char c;
    std::string description;
    std::string mode;

    iss >> name;

    // Skip to right after the first double quote (")
    while ( iss >> c && c != '\"' );
    // Add the characters that follow until the next double quote
    while ( iss >> c && c != '\"' ) description += c;

    iss >> mode;
    if ( mode == "univ" ) {
      std::string univ_name;
      iss >> univ_name;

      const Universe* univ = nullptr;
      if ( univ_name == "CV" ) {
        univ = &syst_->cv_universe();
      }
      else if ( univ_name == "FakeData" ) {
        univ = syst_->fake_data_universe().get();
        if ( !univ ) {
          throw std::runtime_error( "Fake data prediction requested"
            " when not using fake data" );
        }
      }
      else {
        throw std::runtime_error( "Unrecognized universe specifier" );
      }

      pred = new UniverseTrueEvents( num_bins, name, *univ );
    }
    else if ( mode == "file" ) {
      std::string file_name;
      std::string hist_name;

      iss >> file_name >> hist_name;

      pred = new FileTrueEvents( num_bins, name, file_name,
        hist_name, conv_factor );
    }
    else {
      throw std::runtime_error( "Unrecognized prediction mode \""
        + mode + '\"' );
    }

    // Add the fully-initialized prediction to the owned map. Note that any
    // pre-existing entry with the same description will be replaced.
    pred_map_[ description ].reset( pred );
  }

}
