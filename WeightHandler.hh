#pragma once

// Standard library includes
#include <algorithm>
#include <array>
#include <map>
#include <string>

// ROOT includes
#include "TTree.h"

// STV analysis includes
#include "TreeUtils.hh"

// Class that provides temporary storage for event weights being processed by a
// UniverseMaker object
class WeightHandler {
  public:

    WeightHandler() {};

    // Configure storage and set input branch addresses for processing event
    // weights from the input TTree. If a non-null pointer to a vector of
    // branch names is supplied, use it to decide which branches to include.
    // Otherwise, auto-detect appropriate branches using the TTree itself.
    void set_branch_addresses( TTree& in_tree,
      const std::vector<std::string>* branch_names = nullptr );

    // Add a single branch from the input TTree with the specified name.
    // This function also sets the branch address appropriately.
    void add_branch( TTree& in_tree, const std::string& branch_name,
      bool throw_when_missing = true );

    // Overloaded version that allows for easy configuration of a single input
    // branch
    void set_branch_addresses( TTree& in_tree, const std::string& branch_name );

    // Access the owned map
    inline const auto& weight_map() const { return weight_map_; }
    inline auto& weight_map() { return weight_map_; }

  protected:

    // Keys are branch names in the input TTree, values point to vectors of
    // event weights
    std::map< std::string, MyPointer< std::vector<double> > > weight_map_;
};

void WeightHandler::set_branch_addresses( TTree& in_tree,
  const std::vector<std::string>* branch_names )
{
  // Delete any pre-existing contents of the weight map
  weight_map_.clear();

  // Loop over each of the branches of the input TTree
  auto* lob = in_tree.GetListOfBranches();
  for ( int b = 0; b < lob->GetEntries(); ++b ) {

    auto* branch = dynamic_cast< TBranch* >( lob->At(b) );
    std::string br_name = branch->GetName();

    // If the user specified a vector of branch names when calling this
    // function, then check the current branch against them to decide whether
    // it should be included.
    bool include_branch = false;
    if ( branch_names ) {
      // Include the branch if its name can be found in the user-supplied
      // vector
      auto iter = std::find( branch_names->begin(),
        branch_names->end(), br_name );
      include_branch = ( iter != branch_names->end() );
    }
    else {
      // The user didn't supply a vector of branch names, so resort to
      // automatic detection. Include any branch whose name begins with
      // the string "weight_".
      const std::string wgt_br_prefix = "weight_";
      int compare_result = br_name.compare( 0, wgt_br_prefix.size(),
        wgt_br_prefix );
      include_branch = ( compare_result == 0 );
    }

    // Skip to the next branch name if we don't need to include it
    if ( !include_branch ) continue;

    // Assume that all included branches store a std::vector<double> object.
    // Set the branch address so that the vector can accept input from the
    // TTree.
    weight_map_[ br_name ] = MyPointer< std::vector<double> >();

    auto& wgt_vec = weight_map_.at( br_name );
    set_object_input_branch_address( in_tree, br_name, wgt_vec );

  }

  // TODO: add warning or exception for branches listed in the input vector
  // that could not be found in the TTree?
}

void WeightHandler::set_branch_addresses( TTree& in_tree,
  const std::string& branch_name )
{
  std::vector< std::string > br_names = { branch_name };
  set_branch_addresses( in_tree, &br_names );
}

// Allow the user to manually add a new branch
void WeightHandler::add_branch( TTree& in_tree,
  const std::string& branch_name, bool throw_when_missing )
{
  // If we already have an entry in the map for this branch, just return
  // without doing anything
  auto iter = weight_map_.find( branch_name );
  if ( iter != weight_map_.end() ) return;

  // Don't bother to include the new branch if it doesn't actually exist in the
  // input TTree
  TBranch* br = in_tree.GetBranch( branch_name.c_str() );
  if ( !br ) {
    if ( throw_when_missing ) throw std::runtime_error(
      "Missing TTree branch " + branch_name );
    return;
  }

  // Set up the new branch assuming that it holds a vector of double values
  weight_map_[ branch_name ] = MyPointer< std::vector<double> >();
  auto& wgt_vec = weight_map_.at( branch_name );
  set_object_input_branch_address( in_tree, branch_name, wgt_vec );
}
