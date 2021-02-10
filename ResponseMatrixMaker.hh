#pragma once

// Standard library includes
#include <memory>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// ROOT includes
#include "TH1D.h"
#include "TH2D.h"

// STV analysis includes
#include "TreeUtils.hh"
#include "WeightHandler.hh"

// Define the missing intersection operator for TEntryList objects. Work
// with smart pointers because we're using them elsewhere.
// See tinyurl.com/79q6ztbt for details.
std::unique_ptr< TEntryList > operator*( const std::unique_ptr<TEntryList>& t1,
  const std::unique_ptr<TEntryList>& t2 )
{
  auto intersection = std::make_unique< TEntryList >( *t1 );
  auto temp_copy_t1 = std::make_unique< TEntryList >( *t1 );
  temp_copy_t1->Subtract( t2.get() );
  intersection->Subtract( temp_copy_t1.get() );

  std::string merged_name = t1->GetName();
  merged_name += "_and_";
  merged_name += t2->GetName();

  intersection->SetName( merged_name.c_str() );
  return intersection;
}

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

class ResponseMatrixMaker {

  public:

    // Initialize the ResponseMatrixMaker with true and reco bin definitions
    // stored in a configuration file
    ResponseMatrixMaker( const std::string& config_file_name );

    // Add an ntuple input file to the owned TChain. Update the total
    // simulated POT exposure accordingly.
    // TODO: Consider POT of different classes of input samples appropriately.
    // You can use the same classifications as the FilePropertiesManager class.
    void add_input_file( const std::string& input_file_name );

    // Gets the simulated POT exposure for the MC events used to compute the
    // response matrix elements
    inline double pot() const { return pot_; }

    // Access the bin definitions
    inline const auto& true_bins() const { return true_bins_; }
    inline const auto& reco_bins() const { return reco_bins_; }

    // Populates the owned vectors of TEntryList objects using the input TChain
    // and the current bin configuration
    void build_entry_lists();

    // Does the actual calculation of response matrix elements across the
    // various systematic universes
    void build_response_matrices();

  protected:

    // Bin definitions in true space
    std::vector< TrueBin > true_bins_;

    // Bin definitions in reco space
    std::vector< RecoBin > reco_bins_;

    // The simulated total exposure in protons-on-target (POT) that will used
    // to compute the response matrix elements
    double pot_ = 0.;

    // A TChain containing MC event ntuples that will be used to compute
    // response matrix elements
    TChain input_chain_;

    // Each list in the vector stores the entries in the input TChain which
    // fall into the corresponding true bin
    std::vector< std::unique_ptr<TEntryList> > true_entry_lists_;

    // Each list in the vector stores the entries in the input TChain which
    // fall into the corresponding reco bin
    std::vector< std::unique_ptr<TEntryList> > reco_entry_lists_;

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
  // Check to make sure that the input file contains the expected objects
  TFile temp_file( input_file_name.c_str(), "read" );

  // Temporary storage
  TParameter<float>* summed_pot;
  TTree* temp_tree;

  temp_file.GetObject( "summed_pot", summed_pot );
  if ( !summed_pot ) throw std::runtime_error( "Missing POT parameter"
    " in the input ntuple file " + input_file_name );

  double file_pot = summed_pot->GetVal();
  if ( file_pot <= 0. ) throw std::runtime_error( "Invalid POT parameter"
    " value " + std::to_string(file_pot) + " in the input ntuple file "
    + input_file_name );

  std::string tree_name = input_chain_.GetName();
  temp_file.GetObject( tree_name.c_str(), temp_tree );
  if ( !temp_tree ) throw std::runtime_error( "Missing ntuple TTree "
    + tree_name + " in the input ntuple file " + input_file_name );

  // If we've made it here, then the input file has passed all of the checks.
  // Update the POT exposure tally and add it to the input TChain.
  pot_ += file_pot;
  input_chain_.AddFile( input_file_name.c_str() );
}

void ResponseMatrixMaker::build_entry_lists() {

  // Remove any existing TEntryList objects from the owned vectors. Since we're
  // working with smart pointers to handle storage, pre-existing objects will
  // be automatically deleted.
  true_entry_lists_.clear();
  reco_entry_lists_.clear();

  // Using these counters is a hacky way to ensure that every TEntryList that
  // is automatically generated by this function has a unique ROOT name. Since
  // TTree::Draw relies on this name to create the list, I think it makes sense
  // to play things safe.
  int dummy_tb_counter = 0;
  int dummy_rb_counter = 0;

  // Create one entry list for each true bin
  for ( const auto& tb : true_bins_ ) {

    std::string entry_list_name = "tb_entry_list"
      + std::to_string( dummy_tb_counter );

    std::string temp_var_expr = ">> " + entry_list_name;

    input_chain_.Draw( temp_var_expr.c_str(),
      tb.signal_cuts_.c_str(), "entrylist" );

    auto* temp_el = dynamic_cast<TEntryList*>(
      gDirectory->Get(entry_list_name.c_str()) );

    true_entry_lists_.emplace_back( temp_el );

    ++dummy_tb_counter;
  }

  // Create one entry list for each reco bin
  for ( const auto& rb : reco_bins_ ) {

    std::string entry_list_name = "rb_entry_list"
      + std::to_string( dummy_rb_counter );

    std::string temp_var_expr = ">> " + entry_list_name;

    input_chain_.Draw( temp_var_expr.c_str(),
      rb.selection_cuts_.c_str(), "entrylist" );

    auto* temp_el = dynamic_cast<TEntryList*>(
      gDirectory->Get(entry_list_name.c_str()) );

    reco_entry_lists_.emplace_back( temp_el );

    ++dummy_rb_counter;
  }

  // DEBUG
  for ( const auto& tel : true_entry_lists_ ) {
    std::cout << tel->GetName() << '\n';
    tel->Print();
    std::cout << '\n';
  }

  for ( const auto& rel : reco_entry_lists_ ) {
    std::cout << rel->GetName() << '\n';
    rel->Print();
    std::cout << '\n';
  }

  for ( const auto& tel : true_entry_lists_ ) {
    for ( const auto& rel : reco_entry_lists_ ) {
      // Compute the intersection of the entry lists
      auto inter_tr_el = tel * rel;
      std::cout << '\"' << inter_tr_el->GetName() << "\"\n";
      inter_tr_el->Print();
      std::cout << '\n';
    }
  }

}

void ResponseMatrixMaker::build_response_matrices() {
  // Create temporary storage for the systematic variation event weights
  WeightHandler wh;

  // Sync the branch addresses of the input TChain with the WeightHandler
  // object
  wh.set_branch_addresses( input_chain_ );

  // TESTING CODE
  for ( const auto& tel : true_entry_lists_ ) {
    for ( const auto& rel : reco_entry_lists_ ) {

      // Compute the intersection of the entry lists
      auto inter_tr_el = tel * rel;
      std::cout << '\"' << inter_tr_el->GetName() << "\"\n";
      inter_tr_el->Print();
      std::cout << '\n';

      // Configure the TChain to use it
      input_chain_.SetEntryList( inter_tr_el.get() );

      long long num_entries = inter_tr_el->GetN();

      // Sum the entries in the TChain for the relevant entries
      double sum = 0.;

      for ( long long e = 0; e < num_entries; ++e ) {
        long long entryNumber = input_chain_.GetEntryNumber( e );
        if ( entryNumber < 0 ) break;

        input_chain_.GetEntry( entryNumber );

        sum += wh.bug_fix_weight_map_
          .at("weight_splines_general_Spline")->front();

        std::cout << "  e = " << e << '\n';
      }

      std::cout << "sum = " << sum << '\n';
    }
  }

}
