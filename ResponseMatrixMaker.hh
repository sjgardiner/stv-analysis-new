#pragma once

// Standard library includes
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// ROOT includes
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTreeFormula.h"

// STV analysis includes
#include "TreeUtils.hh"
#include "WeightHandler.hh"

// Only load the library in this way if we're using this code from inside the
// ROOT C++ interpreter. We could check for __CINT__ as well, but the specific
// R__LOAD_LIBRARY approach used here only works with ROOT 6.
#ifdef __CLING__
// Pre-load the definition of the TTreeFormula::EvalInstance function from the
// TreePlayer shared library. This approach is based on the trick mentioned
// here: https://tinyurl.com/2s4yuzxm
R__LOAD_LIBRARY(libTreePlayer.so)
#endif

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

    // Cuts to use in TTree::Draw for filling the bin. Any overall event weight
    // included here will be ignored. It is up to the user to ensure that only
    // true quantities are used for true bin cuts.
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

    // Cuts to use in TTree::Draw for filling the bin. Any overall event weight
    // included here will be ignored. It is up to the user to ensure that only
    // reco quantities are used for reco bin cuts.
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

// Provides a set of histograms used to store summed bin counts (with
// associated MC statistical uncertainties) in a given systematic variation
// universe
class Universe {

  public:

    Universe( const std::string& universe_name,
      size_t universe_index, int num_true_bins, int num_reco_bins )
    {
      std::string hist_name_prefix = universe_name + '_'
        + std::to_string( universe_index );

      hist_true_ = std::make_unique< TH1D >(
        (hist_name_prefix + "_true").c_str(), "; true bin number; events",
        num_true_bins, 0., num_true_bins );

      hist_reco_ = std::make_unique< TH1D >(
        (hist_name_prefix + "_reco").c_str(), "; reco bin number; events",
        num_reco_bins, 0., num_reco_bins );

      hist_2d_ = std::make_unique< TH2D >( (hist_name_prefix + "_2d").c_str(),
        "; true bin number; reco bin number; counts", num_true_bins, 0.,
        num_true_bins, num_reco_bins, 0., num_reco_bins );

      // Store summed squares of event weights (for calculations of the MC
      // statistical uncertainty on bin contents)
      hist_true_->Sumw2();
      hist_reco_->Sumw2();
      hist_2d_->Sumw2();
    }

    std::unique_ptr< TH1D > hist_true_;
    std::unique_ptr< TH1D > hist_reco_;
    std::unique_ptr< TH2D > hist_2d_;
};

class ResponseMatrixMaker {

  public:

    // Initialize the ResponseMatrixMaker with true and reco bin definitions
    // stored in a configuration file
    ResponseMatrixMaker( const std::string& config_file_name );

    // Add an ntuple input file to the owned TChain
    void add_input_file( const std::string& input_file_name );

    // Access the bin definitions
    inline const auto& true_bins() const { return true_bins_; }
    inline const auto& reco_bins() const { return reco_bins_; }

    // Access the owned TChain
    inline auto& input_chain() { return input_chain_; }

    // Does the actual calculation of response matrix elements across the
    // various systematic universes
    void build_response_matrices();

  protected:

    // Prepares the TTreeFormula objects needed to test each entry for
    // membership in each bin
    void prepare_formulas();

    // Prepares the Universe objects needed to store summed event weights for
    // each bin in each systematic variation universe
    void prepare_universes( const WeightHandler& wh );

    // Bin definitions in true space
    std::vector< TrueBin > true_bins_;

    // Bin definitions in reco space
    std::vector< RecoBin > reco_bins_;

    // A TChain containing MC event ntuples that will be used to compute
    // response matrix elements
    TChain input_chain_;

    // TTreeFormula objects used to test whether the current TChain entry falls
    // into each true bin
    std::vector< std::unique_ptr<TTreeFormula> > true_bin_formulas_;

    // TTreeFormula objects used to test whether the current TChain entry falls
    // into each reco bin
    std::vector< std::unique_ptr<TTreeFormula> > reco_bin_formulas_;

    // Stores Universe objects used to accumulate event weights
    std::map< std::string, std::vector<Universe> > universes_;
};

ResponseMatrixMaker::ResponseMatrixMaker( const std::string& config_file_name )
{
  std::ifstream in_file( config_file_name );

  // Load the TTree name to use when processing ntuple input files
  std::string ttree_name;
  in_file >> ttree_name;

  // Initialize the owned input TChain with the configured TTree name
  input_chain_.SetName( ttree_name.c_str() );

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

  // Load the reco bin definitions
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

void ResponseMatrixMaker::add_input_file( const std::string& input_file_name )
{
  // Check to make sure that the input file contains the expected ntuple
  TFile temp_file( input_file_name.c_str(), "read" );

  // Temporary storage
  TTree* temp_tree;

  std::string tree_name = input_chain_.GetName();
  temp_file.GetObject( tree_name.c_str(), temp_tree );
  if ( !temp_tree ) throw std::runtime_error( "Missing ntuple TTree "
    + tree_name + " in the input ntuple file " + input_file_name );

  // If we've made it here, then the input file has passed all of the checks.
  // Add it to the input TChain.
  input_chain_.AddFile( input_file_name.c_str() );
}

void ResponseMatrixMaker::prepare_formulas() {
  // Remove any pre-existing TTreeFormula objects from the owned vectors
  true_bin_formulas_.clear();
  reco_bin_formulas_.clear();

  // Create one TTreeFormula for each true bin definition
  for ( size_t tb = 0u; tb < true_bins_.size(); ++tb ) {
    const auto& bin_def = true_bins_.at( tb );
    std::string formula_name = "true_formula_" + std::to_string( tb );

    auto tbf = std::make_unique< TTreeFormula >( formula_name.c_str(),
      bin_def.signal_cuts_.c_str(), &input_chain_ );

    tbf->SetQuickLoad( true );

    true_bin_formulas_.emplace_back( std::move(tbf) );
  }

  // Create one TTreeFormula for each reco bin definition
  for ( size_t rb = 0u; rb < reco_bins_.size(); ++rb ) {
    const auto& bin_def = reco_bins_.at( rb );
    std::string formula_name = "reco_formula_" + std::to_string( rb );

    auto rbf = std::make_unique< TTreeFormula >( formula_name.c_str(),
      bin_def.selection_cuts_.c_str(), &input_chain_ );

    rbf->SetQuickLoad( true );

    reco_bin_formulas_.emplace_back( std::move(rbf) );
  }

}

void ResponseMatrixMaker::build_response_matrices() {

  int num_input_files = input_chain_.GetListOfFiles()->GetEntries();
  if ( num_input_files < 1 ) {
    std::cout << "ERROR: The ResponseMatrixMaker object has not been"
      " initialized with any input files yet.\n";
    return;
  }

  WeightHandler wh;
  wh.set_branch_addresses( input_chain_ );

  this->prepare_formulas();

  bool need_to_prepare_universes = true;

  int treenumber = 0;
  for ( long long entry = 0; entry < input_chain_.GetEntries(); ++entry ) {
    // Load the TTree for the current TChain entry
    input_chain_.LoadTree( entry );

    // If the current entry is in a new TTree, then have all of the
    // TTreeFormula objects make the necessary updates
    if ( treenumber != input_chain_.GetTreeNumber() ) {
      treenumber = input_chain_.GetTreeNumber();
      for ( auto& tbf : true_bin_formulas_ ) tbf->Notify();
      for ( auto& rbf : reco_bin_formulas_ ) rbf->Notify();
    }

    std::vector< size_t > matched_true_bins;
    std::vector< size_t > matched_reco_bins;

    for ( size_t tb = 0u; tb < true_bin_formulas_.size(); ++tb ) {
      auto& tbf = true_bin_formulas_.at( tb );
      if ( tbf->EvalInstance() ) matched_true_bins.push_back( tb );
    }
    for ( size_t rb = 0u; rb < reco_bin_formulas_.size(); ++rb ) {
      auto& rbf = reco_bin_formulas_.at( rb );
      if ( rbf->EvalInstance() ) matched_reco_bins.push_back( rb );
    }

    input_chain_.GetEntry( entry );

    std::cout << "Entry " << entry << '\n';

    if ( need_to_prepare_universes ) {
      this->prepare_universes( wh );
      need_to_prepare_universes = false;
    }

    for ( const auto& pair : wh.weight_map() ) {
      const std::string& wgt_name = pair.first;
      const auto& wgt_vec = pair.second;

      auto& u_vec = universes_.at( wgt_name );

      for ( size_t u = 0u; u < wgt_vec->size(); ++u ) {

        double w = wgt_vec->at( u );
        auto& universe = u_vec.at( u );

        // DEBUG: need to process weight w
        for ( const int& tb : matched_true_bins ) {
          universe.hist_true_->Fill( tb, w );
          for ( const int& rb : matched_reco_bins ) {
            universe.hist_2d_->Fill( tb, rb, w );
          } // reco bins
        } // true bins

        for ( const int& rb : matched_reco_bins ) {
          universe.hist_reco_->Fill( rb, w );
        } // reco bins
      } // universes
    } // weight names
  } // TChain entries

}

void ResponseMatrixMaker::prepare_universes( const WeightHandler& wh ) {
  for ( const auto& pair : wh.weight_map() ) {
    const std::string& weight_name = pair.first;
    size_t num_universes = pair.second->size();

    size_t num_true_bins = true_bins_.size();
    size_t num_reco_bins = reco_bins_.size();

    std::vector< Universe > u_vec;

    for ( size_t u = 0u; u < num_universes; ++u ) {
      u_vec.emplace_back( weight_name, u, num_true_bins, num_reco_bins );
    }

    universes_[ weight_name ] = std::move( u_vec );
  }
}
