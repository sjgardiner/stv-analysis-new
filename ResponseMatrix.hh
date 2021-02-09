#pragma once

#include <memory>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TH1D.h"
#include "TH2D.h"

// Enum used to label bin types in true space
enum TrueBinType {

  // Bin contains signal events
  kSignalTrueBin = 0,

  // Bin contains background events
  kBackgroundTrueBin = 1,

  // Placeholder for undefined cases
  kUnknownTrueBin = 2

};

// Defines a histogram bin in true space
struct TrueBin {

  public:

    TrueBin( const std::string& cuts = "",
      TrueBinType bin_type = kSignalTrueBin )
      : signal_cuts_( cuts ), type_( bin_type ) {}

    // Cuts to use in TTree::Draw for filling the bin. An overall event weight
    // used across multiple bins does not need to be included here. It is up to
    // the user to ensure that only true quantities are used for true bin cuts.
    std::string signal_cuts_;

    // Indicates the nature of the events contained in this bin
    TrueBinType type_;

};

inline std::ostream& operator<<( std::ostream& out, const TrueBin& tb ) {
  out << tb.type_ << " \"" << tb.signal_cuts_ << '\"';
  return out;
}

inline std::istream& operator>>( std::istream& in, TrueBin& tb ) {

  int temp_type;
  in >> temp_type;

  tb.type_ = static_cast< TrueBinType >( temp_type );

  // Use two calls to std::getline using a double quote delimiter
  // in order to get the contents of the next double-quoted string
  std::string temp_line;
  std::getline( in, temp_line, '\"' );
  std::getline( in, temp_line, '\"' );
  tb.signal_cuts_ = temp_line;

  return in;
}

// Defines a histogram bin in reco space
struct RecoBin {

  public:

    RecoBin( const std::string& cuts = "" ) : selection_cuts_( cuts ) {}

    // Cuts to use in TTree::Draw for filling the bin. It is up to the user to
    // ensure that only reco quantities are used for reco bin cuts.
    std::string selection_cuts_;

};

inline std::ostream& operator<<( std::ostream& out, const RecoBin& rb ) {
  out << '\"' << rb.selection_cuts_ << '\"';
  return out;
}

inline std::istream& operator>>( std::istream& in, RecoBin& rb ) {
  std::string temp_line;

  // Use two calls to std::getline using a double quote delimiter
  // in order to get the contents of the next double-quoted string
  std::getline( in, temp_line, '\"' );
  std::getline( in, temp_line, '\"' );
  rb.selection_cuts_ = temp_line;

  return in;
}

class ResponseMatrix {

  public:

    ResponseMatrix( const std::string& config_file_name );

    // Updates the overall weighting factor for filling the histograms
    // with MC event counts. Note that this function will call reset()
    // in order to enforce consistency.
    void set_overall_mc_weight( const std::string& weight_expr );

    // Maintains the bin definitions and the definition of the overall MC
    // weight while zeroing out the POT exposure and all histogram contents
    void reset();

    // Fills the owned histograms using the current configuration for binning
    // and MC weighting
    void fill( TTree& in_tree, double mc_pot );

    // Gets the simulated POT exposure used to fill the owned histograms
    inline double pot() { return pot_; }

    // Access the owned histograms
    inline TH1D* true_histogram() { return true_hist_.get(); }
    inline TH2D* smearing_histogram() { return smear_hist_.get(); }

    // Access the bin definitions
    inline const auto& true_bins() { return true_bins_; }
    inline const auto& reco_bins() { return reco_bins_; }

  protected:

    void reset_true_histogram();
    void reset_smearing_histogram();

    // Bin definitions in true space
    std::vector< TrueBin > true_bins_;

    // Bin definitions in reco space
    std::vector< RecoBin > reco_bins_;

    // An overall event weight to use when filling histograms via TTree::Draw.
    // This can be applied to evaluate the response matrix in a particular
    // systematic universe, etc.
    std::string overall_mc_weight_ = "";

    // The simulated total exposure in protons-on-target (POT) used to fill the
    // owned histograms
    double pot_ = 0.;

    // Stores the sum of the generated event weights (and the associated MC
    // statistical uncertainty) in each true bin
    std::unique_ptr< TH1D > true_hist_;

    // Stores the sum of the generated event weights (and the associated MC
    // statistical uncertainty) in each (reco, true) bin pair
    std::unique_ptr< TH2D > smear_hist_;

};

void ResponseMatrix::reset_true_histogram() {

  // Get the current number of bins in true space
  size_t num_true_bins = true_bins_.size();

  // If there aren't any, replace any existing true histogram
  // with a null pointer and return without doing anything else
  if ( num_true_bins == 0u ) {
    true_hist_.reset();
    return;
  }

  // Static counter used to ensure that the new true histogram always has a
  // unique ROOT name
  static int true_hist_count = 0;
  std::string temp_true_hist_name = "RespMat_true_hist"
    + std::to_string( true_hist_count );

  // Otherwise, create a new true histogram using the correct number of bins.
  // Use the bin's index in the true_bins_ vector (zero-based) as the x-axis
  // variable.
  auto temp_true_hist = std::make_unique< TH1D >( temp_true_hist_name.c_str(),
    "; true bin number; sum of event weights", num_true_bins, 0.,
    static_cast<double>(num_true_bins) );

  // Enable automatic tracking of statistical uncertainties (the sum of the
  // squares of the event weights used to fill the histogram)
  temp_true_hist->Sumw2();

  // Take ownership of the temporary histogram
  true_hist_.reset( temp_true_hist.release() );
}

void ResponseMatrix::reset_smearing_histogram() {

  // Get the current number of bins in true space
  size_t num_true_bins = true_bins_.size();

  // Get the current number of bins in reco space
  size_t num_reco_bins = reco_bins_.size();

  // If the number of bins along either axis is zero, then replace any existing
  // smearing histogram with a null pointer and return without doing anything
  // else
  if ( num_true_bins == 0u || num_reco_bins == 0u ) {
    smear_hist_.reset();
    return;
  }

  // Static counter used to ensure that the new smearing histogram always has a
  // unique ROOT name
  static int smear_hist_count = 0;
  std::string temp_smear_hist_name = "RespMat_smear_hist"
    + std::to_string( smear_hist_count );

  // Otherwise, create a new smearing histogram using the correct number of
  // bins along each axis. Use the bin's index in the true_bins_ vector
  // (zero-based) as the x-axis variable, and use the bin's index in the
  // reco_bins_ vector (also zero-based) as the y-axis variable.
  auto temp_smear_hist = std::make_unique< TH2D >( temp_smear_hist_name.c_str(),
    "; true bin number; reco bin number; sum of event weights", num_true_bins,
    0., static_cast<double>(num_true_bins), num_reco_bins, 0.,
    static_cast<double>(num_reco_bins) );

  // Enable automatic tracking of statistical uncertainties (the sum of the
  // squares of the event weights used to fill the histogram)
  temp_smear_hist->Sumw2();

  // Take ownership of the temporary histogram
  smear_hist_.reset( temp_smear_hist.release() );
}

void ResponseMatrix::reset() {
  this->reset_true_histogram();
  this->reset_smearing_histogram();
  pot_ = 0.;
}

void ResponseMatrix::set_overall_mc_weight( const std::string& weight_expr ) {
  overall_mc_weight_ = weight_expr;
  this->reset();
}

// TODO: enable "add-to" filling as well as "overwrite" filling
void ResponseMatrix::fill( TTree& in_tree, double mc_pot ) {

  // Save the simulated POT exposure in case it is needed later for normalizing
  // results
  pot_ = mc_pot;

  // Fill the true bins one by one from the input TTree
  const std::string true_hist_name = true_hist_->GetName();

  // DEBUG
  std::cout << "Filling true bins\n";

  for ( size_t tb = 0u; tb < true_bins_.size(); ++tb ) {
    // The x-axis variable is just the true bin number, so we can
    // fill the appropriate bin via TTree::Draw just by using it
    // as the x value.
    std::string fill_var_expr = std::to_string( tb ) + " >>+ " + true_hist_name;

    // Apply the cuts that define this true bin. Include the factor defined as
    // the overall MC event weight.
    const TrueBin& bin_def = true_bins_.at( tb );
    std::string cuts_expr = '(' + overall_mc_weight_ + ")*("
      + bin_def.signal_cuts_ + ')';

    // We're ready. Do the actual filling using TTree::Draw. Use the "goff"
    // option to improve performance by avoiding graphics function calls.
    in_tree.Draw( fill_var_expr.c_str(), cuts_expr.c_str(), "goff" );

    // DEBUG
    std::cout << "Filled true bin #" << tb << '\n';
  }

  // Fill the smearing bins one by one from the input TTree
  const std::string smear_hist_name = smear_hist_->GetName();

  // DEBUG
  std::cout << "Filling the smearing matrix\n";

  for ( size_t tb = 0u; tb < true_bins_.size(); ++tb ) {
    for ( size_t rb = 0u; rb < reco_bins_.size(); ++rb ) {
      // The y-axis variable is just the reco bin number, so we can fill the
      // appropriate bin via TTree::Draw just by using it as the y value.
      std::string fill_var_expr = std::to_string( rb );

      // Likewise, the x-axis variable is just the true bin number.
      fill_var_expr += ':' + std::to_string( tb ) + " >>+ " + smear_hist_name;

      // Apply the cuts that define this smearing histogram bin. Include the
      // factor defined as the overall MC event weight.
      const TrueBin& true_bin_def = true_bins_.at( tb );
      const RecoBin& reco_bin_def = reco_bins_.at( rb );

      std::string cuts_expr = '(' + overall_mc_weight_ + ")*("
        + true_bin_def.signal_cuts_ + " && "
        + reco_bin_def.selection_cuts_ + ')';

      // We're ready. Do the actual filling using TTree::Draw. Use the "goff"
      // option to improve performance by avoiding graphics function calls.
      in_tree.Draw( fill_var_expr.c_str(), cuts_expr.c_str(), "goff" );

      // DEBUG
      std::cout << "Filled smearing matrix tb = " << tb
        << ", rb = " << rb << '\n';

    } // loop over reco bins
  } // loop over true bins
}

ResponseMatrix::ResponseMatrix( const std::string& config_file_name ) {
  std::ifstream in_file( config_file_name );

  // Load the true bin definitions
  size_t num_true_bins;
  in_file >> num_true_bins;

  for ( size_t tb = 0u; tb < num_true_bins; ++tb ) {
    TrueBin temp_bin;
    in_file >> temp_bin;

    // DEBUG
    std::cout << "tb = " << tb << '\n';
    std::cout << temp_bin << '\n';

    true_bins_.push_back( temp_bin );
  }

  size_t num_reco_bins;
  in_file >> num_reco_bins;
  for ( size_t rb = 0u; rb < num_reco_bins; ++rb ) {
    RecoBin temp_bin;
    in_file >> temp_bin;

    // DEBUG
    std::cout << "rb = " << rb << '\n';
    std::cout << temp_bin << '\n';

    reco_bins_.push_back( temp_bin );
  }
}
