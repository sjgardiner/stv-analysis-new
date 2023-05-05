#pragma once

#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>

#include "TH1D.h"

constexpr int DUMMY_BIN_INDEX = 0;

struct SliceVariable {

  SliceVariable() {}

  SliceVariable( const std::string& name, const std::string& latex_name,
    const std::string& units, const std::string& latex_units )
    : name_( name ), latex_name_( latex_name ), units_( units ),
    latex_units_( latex_units ) {}

  std::string name_; // Name to be used when plotting ROOT histograms
  std::string latex_name_; // Form of the name to be used in a LaTeX document
  std::string units_; // Unit specification for ROOT plotting
  std::string latex_units_; // Unit specification (for LaTeX documents)

};


// Reads a string surrounded by double quotes (") from an input stream
std::string get_double_quoted_string( std::istream& in ) {
  std::string temp_str;
  std::getline( in, temp_str, '\"' );
  std::getline( in, temp_str, '\"' );
  return temp_str;
}

std::istream& operator>>( std::istream& in, SliceVariable& svar ) {
  svar.name_ = get_double_quoted_string( in );
  svar.units_ = get_double_quoted_string( in );
  svar.latex_name_ = get_double_quoted_string( in );
  svar.latex_units_ = get_double_quoted_string( in );
  return in;
}

std::ostream& operator<<( std::ostream& out, const SliceVariable& svar )
{
  out << '\"' << svar.name_ << "\" \"" << svar.units_ << "\" \""
    << svar.latex_name_ << "\" \"" << svar.latex_units_ << '\"';
  return out;
}

struct Slice {

  Slice() {}

  // ROOT histogram storing the contents of the slice
  std::unique_ptr< TH1 > hist_;

  // Keys are global bin numbers in hist_ (one-based). Values are the
  // contributing reco bin index (zero-based) as defined in the relevant
  // UniverseMaker configuration
  std::map< int, std::set< size_t > > bin_map_;

  // Indices in the slice_vars_ vector for the "active" SliceVariable
  // definitions used to define axes of the owned TH1
  std::vector< size_t > active_var_indices_;

  // Specification for the "other" relevant SliceVariable values
  struct OtherVariableSpec {

    OtherVariableSpec() {}

    OtherVariableSpec( size_t var_idx, double low, double high )
      : var_index_( var_idx ), low_bin_edge_( low ), high_bin_edge_( high ) {}

    // Index for the relevant SliceVariable definition in the slice_vars_
    // vector
    size_t var_index_;
    // Lower bin edge for the current variable in the current slice
    double low_bin_edge_;
    // Upper bin edge for the current variable in the current slice
    double high_bin_edge_;
  };

  // Specifications for the "other" SliceVariable values used to define
  // the slice
  std::vector< OtherVariableSpec > other_vars_;

};


std::ostream& operator<<( std::ostream& out,
  const Slice::OtherVariableSpec& ovs )
{
  out << ovs.var_index_ << ' ' << ovs.low_bin_edge_ << ' '
    << ovs.high_bin_edge_;

  return out;
}

std::ostream& operator<<( std::ostream& out, const Slice& slice )
{
  // TODO: add handling for multidimensional slices
  std::string final_axis_label = slice.hist_->GetYaxis()->GetTitle();
  size_t num_active_vars = slice.active_var_indices_.size();
  if ( num_active_vars != 1u ) throw std::runtime_error( "Support for"
    " multidimensional slices has not yet been implemented" );

  out << '\"' << final_axis_label << "\"\n";
  out << num_active_vars;

  size_t avar_idx = slice.active_var_indices_.front();
  int num_edges = slice.hist_->GetNbinsX() + 1;
  out << ' ' << avar_idx << ' ' << num_edges;

  // Note that ROOT TH1 bins have one-based indices (bin zero is used for
  // underflow). By using num_edges as the maximum bin index in the loop below,
  // we also get the lower edge of the overflow bin (equivalent to the upper
  // edge of the last regular bin).
  for ( int bin = 1; bin <= num_edges; ++bin ) {
    out << ' ' << slice.hist_->GetBinLowEdge( bin );
  }
  out << '\n';

  size_t num_other_vars = slice.other_vars_.size();
  out << num_other_vars;
  if ( num_other_vars == 0 ) out << '\n';
  else out << ' ';

  for ( const auto& ovar_spec : slice.other_vars_ ) out << ovar_spec << '\n';

  // We need to "invert" the bin map in order for it to have the correct
  // structure for dumping the configuration easily. Here we'll manually
  // build a map from analysis bins to ROOT histogram bins in the slice.
  // TODO: consider other solutions
  std::map< size_t, std::vector< int > > analysis_to_root_bin_map;
  for ( const auto& bin_pair : slice.bin_map_ ) {
    int root_bin_idx = bin_pair.first;
    const auto& analysis_bin_set = bin_pair.second;

    for ( const size_t& analysis_bin_idx : analysis_bin_set ) {
      // Get access to the vector of ROOT histogram bins for the current
      // analysis bin. Note that the call to std::map::operator[] will create
      // a map entry (with an empty vector) if one does not already exist.
      auto& root_bin_vec = analysis_to_root_bin_map[ analysis_bin_idx ];
      root_bin_vec.push_back( root_bin_idx );
    }
  }

  // TODO: Generalize this for multidimensional slices. At present, we
  // assume here that the *global* ROOT bin numbers (as stored in bin_map_)
  // correspond to the active variable bin numbers on the x-axis. This is
  // only reliably true for 1D histograms, so we're currently cheating a bit.
  size_t num_analysis_bins = analysis_to_root_bin_map.size();
  out << num_analysis_bins;

  for ( const auto& ana_bin_pair : analysis_to_root_bin_map ) {
    const size_t ana_bin_idx = ana_bin_pair.first;
    const auto& root_bin_vec = ana_bin_pair.second;
    size_t num_root_bins = root_bin_vec.size();

    out << '\n' << ana_bin_idx << ' ' << num_root_bins;
    for ( const auto& rbin_idx : root_bin_vec ) {
      out << ' ' << rbin_idx;
    }
  }

  return out;
}

// Defines "slices" of a possibly multidimensional phase space to use for
// plotting results calculated in terms of reco/true bin counts or functions
// thereof. These slices are represented by ROOT histograms.
class SliceBinning {

  public:

    // Use this constructor if you want to fill the SliceBinning object
    // programmatically (i.e., not via a pre-existing configuration file)
    SliceBinning() {}

    // Construct the SliceBinning object from a saved configuration file
    SliceBinning( const std::string& config_file_name );

    // Prints the current configuration of the object to a std::ostream.
    // This function can be used to create a new configuration file.
    void print_config( std::ostream& os ) const;

  //protected:

    std::vector< SliceVariable > slice_vars_;

    std::vector< Slice > slices_;
};

SliceBinning::SliceBinning( const std::string& config_file_name ) {
  std::ifstream in_file( config_file_name );

  // First read in the number of slice variables defined for this binning
  // scheme
  int num_variables;
  in_file >> num_variables;

  std::cout << "NUM VARIABLES = " << num_variables << '\n';

  // Get the slice variable specifications from the configuration file
  SliceVariable temp_sv;
  for ( int v = 0; v < num_variables; ++v ) {
    in_file >> temp_sv;
    slice_vars_.push_back ( temp_sv );
  }

  // Read in the number of slices defined for this binning scheme
  int num_slices;
  in_file >> num_slices;
  std::cout << "NUM SLICES = " << num_slices << '\n';

  // Get the slice specifications from the configuration file
  for ( int s = 0; s < num_slices; ++s ) {

    // Get the label for the final axis
    std::string final_axis_label = get_double_quoted_string( in_file );

    // Number of variables for which the current slice provides explicit ROOT
    // histogram bins
    int num_active_variables;
    in_file >> num_active_variables;

    if ( num_active_variables != 1 ) { // || num_active_variables > 2 )
      throw std::runtime_error( "Handling of "
        + std::to_string(num_active_variables) + " active variables for a"
        " slice is unimplemented!" );
    }

    // Keys are indices for the active variable(s) in the slice_vars_ vector.
    // Values are vectors of TH1-style bin edges.
    std::map< size_t, std::vector<double> > edge_map;
    for ( int av = 0; av < num_active_variables; ++av ) {
      size_t var_idx;
      int num_edges;
      in_file >> var_idx >> num_edges;
      int num_bins = num_edges - 1;
      std::cout << "Slice " << s << " with active variable "
        << slice_vars_.at( var_idx ).name_
        << " has " << num_bins << " bins\n";

      // Create a new entry in the map of edges for the current active variable
      // TODO: add error handling for the case where a duplicate active
      // variable index appears in the input
      auto& edge_vec = edge_map[ var_idx ] = std::vector<double>();

      // We have one more edge listed than the number of bins. This allows the
      // upper edge of the last bin to be specified
      for ( int e = 0; e < num_edges; ++e ) {
        double edge;
        in_file >> edge;
        edge_vec.push_back( edge );
      }

      for ( int b = 1; b <= num_bins; ++b ) {
        std::cout << "Bin " << b << ": " << edge_vec.at( b - 1 )
          << ' ' << slice_vars_.at( var_idx ).units_ << " \u2264 "
          << slice_vars_.at( var_idx ).name_ << " < "
          << edge_vec.at( b )
          << ' ' << slice_vars_.at( var_idx ).units_ << '\n';
      }
    }

    // Create an object to represent the current slice
    Slice cur_slice;

    // Read in the specifications for the other (i.e., inactive) variables
    int num_other_variables;
    in_file >> num_other_variables;

    for ( int ov = 0; ov < num_other_variables; ++ov ) {
      Slice::OtherVariableSpec ovs;
      in_file >> ovs.var_index_ >> ovs.low_bin_edge_ >> ovs.high_bin_edge_;

      std::cout << "Slice has other variable "
        << slice_vars_.at( ovs.var_index_ ).name_ << " on ["
        << ovs.low_bin_edge_ << ", " << ovs.high_bin_edge_ << ")\n";

      cur_slice.other_vars_.push_back( ovs );
    }

    // We're ready to build the ROOT histogram owned by the current slice, so
    // do it now before matching ROOT bins with UniverseMaker bins. Start
    // by building the plot title using the bin limits for the "other"
    // variables.
    std::string slice_title = "slice " + std::to_string( slices_.size() );
    if ( !cur_slice.other_vars_.empty() ) {
      bool first_other_var = true;
      for ( const auto& ovs : cur_slice.other_vars_ ) {
        if ( first_other_var ) {
          slice_title += ": ";
          first_other_var = false;
        }
        else slice_title += ", ";

        // Specify the limits for the current variable in the histogram title
        const auto& var_spec = slice_vars_.at( ovs.var_index_ );
        // Use a std::stringstream to easily get reasonable precision on
        // the numerical bin limits
        std::stringstream temp_ss;
        temp_ss << ovs.low_bin_edge_ << ' ' << var_spec.units_ << " #leq "
          << var_spec.name_ << " < " << ovs.high_bin_edge_
          << ' ' << var_spec.units_;

        slice_title += temp_ss.str();
      }
    }

    // Now label the axes appropriately with the active variable names and
    // units. Retrieve the active variable indices in the slice_vars_ vector
    // from the edge_map.
    for ( const auto& pair : edge_map ) {
      size_t var_idx = pair.first;
      const auto& var_spec = slice_vars_.at( var_idx );

      slice_title += "; " + var_spec.name_;
      if ( !var_spec.units_.empty() ) {
        slice_title += " (" + var_spec.units_ + ')';
      }

      // While we're at it, store the active variable indices in the slice
      // for later retrieval after the edge_map goes out of scope
      cur_slice.active_var_indices_.push_back( var_idx );
    }

    // Name the slice histogram
    std::string slice_hist_name = "slice_" + std::to_string( s );

    // Also add the label for the final axis to the title
    slice_title += ';' + final_axis_label;

    if ( num_active_variables == 1 ) {
      // Retrieve the vector of bin edges for the only active variable
      const auto& bin_edges = edge_map.cbegin()->second;
      int num_bins = bin_edges.size() - 1;

      std::unique_ptr< TH1 > temp_hist(
        new TH1D( slice_hist_name.c_str(), slice_title.c_str(), num_bins,
          bin_edges.data() )
      );

      // Move the ready-to-use histogram into the member smart pointer
      cur_slice.hist_.swap( temp_hist );
    }
    //else if ( num_active_variables == 2 ) {
    // TODO: finish
    //}

    // Prevent auto-deletion of the histogram by ROOT by disassociating
    // it from the current directory
    cur_slice.hist_->SetDirectory( nullptr );

    // Also get rid of the default display of the stats box
    cur_slice.hist_->SetStats( false );

    // Now that we have the histogram made, read in the information needed to
    // map from each UniverseMaker reco bin that contributes to this
    // slice to a TH1 global bin number (conveniently converted from x, y, ...
    // bin indices by the new histogram)
    int num_rmm_bins;
    in_file >> num_rmm_bins;

    std::cout << "RMM bin count = " << num_rmm_bins << '\n';

    for ( int cb = 0; cb < num_rmm_bins; ++cb ) {

      // Index of the appropriate reco bin owned by a UniverseMaker
      // object. This index is zero-based.
      size_t rmm_reco_bin_idx;

      // Number of TH1 bins in the current slice to which the
      // current UniverseMaker bin contributes
      size_t num_root_bins;

      in_file >> rmm_reco_bin_idx >> num_root_bins;

      for ( int rtb = 0; rtb < num_root_bins; ++rtb ) {

        // Store the bin indices along each defined axis
        std::vector< int > av_bin_indices;
        for ( int a = 0; a < num_active_variables; ++a ) {
          int idx;
          in_file >> idx;
          av_bin_indices.push_back( idx );
        }

        // Pad with zeros for y and z if needed so that we can use the
        // same function call to get the global bin number for a one-, two-,
        // or three-dimensional histogram
        int pad_entries = 3 - av_bin_indices.size();
        for ( int p = 0; p < pad_entries; ++p ) {
          av_bin_indices.push_back( DUMMY_BIN_INDEX );
        }

        // We're ready. Look up the global bin index in the slice histogram
        // that should be associated with the current UniverseMaker
        // reco bin
        int root_bin_idx = cur_slice.hist_->GetBin(
          av_bin_indices.front(), av_bin_indices.at(1), av_bin_indices.at(2) );

        // Store the indices for both kinds of bins at the appropriate place
        // in the bin map
        auto iter = cur_slice.bin_map_.find( root_bin_idx );
        if ( iter == cur_slice.bin_map_.end() ) {
          cur_slice.bin_map_[ root_bin_idx ]
            = std::set< size_t > { rmm_reco_bin_idx };
        }
        else iter->second.insert( rmm_reco_bin_idx );

        std::cout << "RMM bin " << rmm_reco_bin_idx << " is matched to "
          << root_bin_idx << " in this slice\n";

      } // ROOT bins

    } // UniverseMaker bins

    // Move the completed Slice object into the vector of slices
    slices_.emplace_back( std::move(cur_slice) );

  } // slices

}

void SliceBinning::print_config( std::ostream& out ) const {

  size_t num_variables = slice_vars_.size();
  out << num_variables << '\n';

  for ( const auto& svar : slice_vars_ ) out << svar << '\n';

  size_t num_slices = slices_.size();
  out << num_slices;

  for ( const auto& slice : slices_ ) out << '\n' << slice;
}

std::ostream& operator<<( std::ostream& out, const SliceBinning& sb ) {
  sb.print_config( out );
  return out;
}
