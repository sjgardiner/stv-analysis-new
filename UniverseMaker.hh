#pragma once

// Standard library includes
#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// ROOT includes
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTreeFormula.h"

// STV analysis includes
#include "EventCategory.hh"
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

// Keys used to identify true and reco bin configurations in a response
// matrix file
const std::string TRUE_BIN_SPEC_NAME = "true_bin_spec";
const std::string RECO_BIN_SPEC_NAME = "reco_bin_spec";

// Converts the name of an analysis ntuple file (typically with the full path)
// into a TDirectoryFile name to use as a subfolder of the main output
// TDirectoryFile used for saving response matrices. Since the forward slash
// character '/' cannot be used in a TDirectoryFile name, this function
// replaces all instances of this character by a '+' instead. The technique
// used here is based on https://stackoverflow.com/a/2896627/4081973
std::string ntuple_subfolder_from_file_name( const std::string& file_name ) {
  constexpr char FORWARD_SLASH = '/';
  constexpr char PLUS = '+';

  std::string result = file_name;
  std::replace( result.begin(), result.end(), FORWARD_SLASH, PLUS );
  return result;
}

// Branch names for special event weights
const std::string SPLINE_WEIGHT_NAME = "weight_splines_general_Spline";
const std::string TUNE_WEIGHT_NAME = "weight_TunedCentralValue_UBGenie";

// Special weight name to store the unweighted event counts
const std::string UNWEIGHTED_NAME = "unweighted";

constexpr double MIN_WEIGHT = 0.;
constexpr double MAX_WEIGHT = 30.;

// Event weights that are below MIN_WEIGHT, above MAX_WEIGHT, infinite, or NaN
// are reset to unity by this function. Other weights are returned unaltered.
inline double safe_weight( double w ) {
  if ( std::isfinite(w) && w >= MIN_WEIGHT && w <= MAX_WEIGHT ) return w;
  else return 1.0;
}

// Utility function used to check endings of (trimmed) weight labels based on
// branch names in the weights TTree
bool string_has_end( const std::string& str, const std::string& end ) {
  if ( str.length() >= end.length() ) {
    int comp = str.compare( str.length() - end.length(),
      end.length(), end );
    bool test_result = ( comp == 0 );
    return test_result;
  }
  return false;
}

// Multiplies a given event weight by extra correction factors as appropriate.
// TODO: include the rootino_fix weight as a correction to the central value
void apply_cv_correction_weights( const std::string& wgt_name,
  double& wgt, double spline_weight, double tune_weight )
{
  if ( string_has_end(wgt_name, "UBGenie") ) {
    wgt *= spline_weight;
  }
  else if ( wgt_name == "weight_flux_all"
    || wgt_name == "weight_reint_all"
    || wgt_name == "weight_xsr_scc_Fa3_SCC"
    || wgt_name == "weight_xsr_scc_Fv3_SCC" )
  {
    wgt *= spline_weight * tune_weight;
  }
  else if ( wgt_name == SPLINE_WEIGHT_NAME ) {
    // No extra weight factors needed
    return;
  }
  else throw std::runtime_error( "Unrecognized weight name" );
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

// Enum used to label bin types in reco space
enum RecoBinType {

  // Bin intended for use in the primary measurement of interest
  kOrdinaryRecoBin = 0,

  // Bin intended for use in a sideband constraint
  kSidebandRecoBin = 1,

  // Placeholder for undefined cases
  kUnknownRecoBin = 2

};

// Defines a histogram bin in true space
struct TrueBin {

  public:

    TrueBin( const std::string& cuts = "",
      TrueBinType bin_type = kSignalTrueBin,
      int block = -1 )
      : signal_cuts_( cuts ), type_( bin_type ), block_index_( block ) {}

    // Cuts to use in TTree::Draw for filling the bin. Any overall event weight
    // included here will be ignored. It is up to the user to ensure that only
    // true quantities are used for true bin cuts.
    std::string signal_cuts_;

    // Indicates the nature of the events contained in this bin
    TrueBinType type_;

    // Index used to group true signal bins together for purposes of
    // unfolding. The numerical value is otherwise ignored.
    int block_index_;
};

inline std::ostream& operator<<( std::ostream& out, const TrueBin& tb ) {
  out << tb.type_ << ' ' << tb.block_index_ << " \""
    << tb.signal_cuts_ << '\"';
  return out;
}

inline std::istream& operator>>( std::istream& in, TrueBin& tb ) {

  int temp_type;
  in >> temp_type;

  tb.type_ = static_cast< TrueBinType >( temp_type );

  int block_idx;
  in >> block_idx;

  tb.block_index_ = block_idx;

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

    RecoBin( const std::string& cuts = "",
      RecoBinType bin_type = kOrdinaryRecoBin,
      int block_idx = -1 )
      : type_( bin_type ), selection_cuts_( cuts ),
      block_index_( block_idx ) {}

    // Cuts to use in TTree::Draw for filling the bin. Any overall event weight
    // included here will be ignored. It is up to the user to ensure that only
    // reco quantities are used for reco bin cuts.
    std::string selection_cuts_;

    // The kind of reco bin represented by this object
    RecoBinType type_;

    // Index used to group ordinary reco bins together for purposes of
    // unfolding. The numerical value is otherwise ignored.
    int block_index_;

};

inline std::ostream& operator<<( std::ostream& out, const RecoBin& rb ) {
  out << rb.type_ << ' ' << rb.block_index_ << " \""
    << rb.selection_cuts_ << '\"';
  return out;
}

inline std::istream& operator>>( std::istream& in, RecoBin& rb ) {

  int temp_type;
  in >> temp_type;

  rb.type_ = static_cast< RecoBinType >( temp_type );

  int block_idx;
  in >> block_idx;

  rb.block_index_ = block_idx;

  // Use two calls to std::getline using a double quote delimiter
  // in order to get the contents of the next double-quoted string
  std::string temp_line;
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
      : universe_name_( universe_name ), index_( universe_index )
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

      hist_reco2d_ = std::make_unique< TH2D >(
        (hist_name_prefix + "_reco2d").c_str(),
        "; reco bin number; reco bin number; counts", num_reco_bins, 0.,
        num_reco_bins, num_reco_bins, 0., num_reco_bins );

      // Get the number of defined EventCategory values by checking the number
      // of elements in the "label map" managed by the EventCategoryInterpreter
      // singleton class
      const auto& eci = EventCategoryInterpreter::Instance();
      size_t num_categories = eci.label_map().size();

      hist_categ_ = std::make_unique< TH2D >(
        (hist_name_prefix + "_categ").c_str(),
        "; true event category; reco bin number; counts", num_categories, 0.,
        num_categories, num_reco_bins, 0., num_reco_bins );

      // Store summed squares of event weights (for calculations of the MC
      // statistical uncertainty on bin contents)
      hist_true_->Sumw2();
      hist_reco_->Sumw2();
      hist_2d_->Sumw2();
      hist_categ_->Sumw2();
      hist_reco2d_->Sumw2();
    }

    // Note: the new Universe object takes ownership of the histogram
    // pointers passed to this constructor
    Universe( const std::string& universe_name,
      size_t universe_index, TH1D* hist_true, TH1D* hist_reco, TH2D* hist_2d,
      TH2D* hist_categ, TH2D* hist_reco2d )
      : universe_name_( universe_name ), index_( universe_index ),
      hist_true_( hist_true ), hist_reco_( hist_reco ), hist_2d_( hist_2d ),
      hist_categ_( hist_categ ), hist_reco2d_( hist_reco2d )
    {
      hist_true_->SetDirectory( nullptr );
      hist_reco_->SetDirectory( nullptr );
      hist_2d_->SetDirectory( nullptr );
      hist_categ_->SetDirectory( nullptr );
      hist_reco2d_->SetDirectory( nullptr );
    }

    std::unique_ptr< Universe > clone() const {
      int num_true_bins = hist_2d_->GetXaxis()->GetNbins();
      int num_reco_bins = hist_2d_->GetYaxis()->GetNbins();
      auto result = std::make_unique< Universe >( universe_name_,
        index_, num_true_bins, num_reco_bins );

      result->hist_true_->Add( this->hist_true_.get() );
      result->hist_reco_->Add( this->hist_reco_.get() );
      result->hist_2d_->Add( this->hist_2d_.get() );
      result->hist_categ_->Add( this->hist_categ_.get() );
      result->hist_reco2d_->Add( this->hist_reco2d_.get() );

      return result;
    }

    std::string universe_name_;
    size_t index_;
    std::unique_ptr< TH1D > hist_true_;
    std::unique_ptr< TH1D > hist_reco_;
    std::unique_ptr< TH2D > hist_2d_;
    std::unique_ptr< TH2D > hist_categ_;
    std::unique_ptr< TH2D > hist_reco2d_;
};

class UniverseMaker {

  public:

    // Initialize the UniverseMaker with true and reco bin definitions
    // stored in a configuration file
    UniverseMaker( const std::string& config_file_name );

    // Add an ntuple input file to the owned TChain
    void add_input_file( const std::string& input_file_name );

    // Access the bin definitions
    inline const auto& true_bins() const { return true_bins_; }
    inline const auto& reco_bins() const { return reco_bins_; }

    // Access the owned TChain
    inline auto& input_chain() { return input_chain_; }

    // Does the actual calculation of response matrix elements across the
    // various systematic universes. The optional argument points to a vector
    // of branch names that will be used to retrieve systematic universe
    // weights. If it is omitted, all available ones will be auto-detected and
    // used.
    void build_universes(
      const std::vector<std::string>* universe_branch_names = nullptr );

    // Overloaded version of the function that takes a reference the
    // vector of universe branch names (for convenience). The behavior
    // is the same as the original, but in this case the explicit vector of
    // branch names definitely exists.
    void build_universes(
      const std::vector<std::string>& universe_branch_names );

    // Writes the response matrix histograms to an output ROOT file
    void save_histograms( const std::string& output_file_name,
      const std::string& subdirectory_name, bool update_file = true );

    // Provides read-only access to the map of Universe objects
    const auto& universe_map() const { return universes_; }

    // Returns the name of the TDirectoryFile that will be used to hold the
    // response matrix histograms when they are written to the output ROOT file
    const std::string& dir_name() const { return output_directory_name_; }

  protected:

    // Helper struct that keeps track of bin indices and TTreeFormula weights
    // when filling universe histograms
    struct FormulaMatch {
      // The weight given here is from an evaluation of a TTreeFormula for
      // filling an individual bin (as opposed to an overall event weight which
      // may be given separately)
      FormulaMatch( size_t bin_idx, double wgt ) : bin_index_( bin_idx ),
        weight_( wgt ) {}

      size_t bin_index_;
      double weight_;
    };

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

    // TTreeFormula objects used to test whether the current TChain entry falls
    // into each true EventCategory
    std::vector< std::unique_ptr<TTreeFormula> > category_formulas_;

    // Stores Universe objects used to accumulate event weights
    std::map< std::string, std::vector<Universe> > universes_;

    // Root TDirectoryFile name to use when writing the response matrices to an
    // output ROOT file
    std::string output_directory_name_;
};

UniverseMaker::UniverseMaker( const std::string& config_file_name )
{
  std::ifstream in_file( config_file_name );

  // Load the root TDirectoryFile name to use when writing the response
  // matrices to an output ROOT file
  in_file >> output_directory_name_;

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

void UniverseMaker::add_input_file( const std::string& input_file_name )
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

void UniverseMaker::prepare_formulas() {

  // Remove any pre-existing TTreeFormula objects from the owned vectors
  true_bin_formulas_.clear();
  reco_bin_formulas_.clear();
  category_formulas_.clear();

  // Get the number of defined EventCategory values by checking the number
  // of elements in the "label map" managed by the EventCategoryInterpreter
  // singleton class
  const auto& eci = EventCategoryInterpreter::Instance();
  size_t num_categories = eci.label_map().size();

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

  // Create one TTreeFormula for each true EventCategory
  const auto& category_map = eci.label_map();
  for ( const auto& category_pair : category_map ) {

    EventCategory cur_category = category_pair.first;
    std::string str_category = std::to_string( cur_category );

    std::string category_formula_name = "category_formula_" + str_category;

    std::string category_cuts = "category == " + str_category;

    auto cbf = std::make_unique< TTreeFormula >(
      category_formula_name.c_str(), category_cuts.c_str(), &input_chain_ );

    cbf->SetQuickLoad( true );

    category_formulas_.emplace_back( std::move(cbf) );
  }

}

void UniverseMaker::build_universes(
  const std::vector<std::string>& universe_branch_names )
{
  return this->build_universes( &universe_branch_names );
}

void UniverseMaker::build_universes(
  const std::vector<std::string>* universe_branch_names )
{
  int num_input_files = input_chain_.GetListOfFiles()->GetEntries();
  if ( num_input_files < 1 ) {
    std::cout << "ERROR: The UniverseMaker object has not been"
      " initialized with any input files yet.\n";
    return;
  }

  WeightHandler wh;
  wh.set_branch_addresses( input_chain_, universe_branch_names );

  // Make sure that we always have branches set up for the CV correction
  // weights, i.e., the spline and tune weights. Don't throw an exception if
  // these are missing in the input TTree (we could be working with real data)
  wh.add_branch( input_chain_, SPLINE_WEIGHT_NAME, false );
  wh.add_branch( input_chain_, TUNE_WEIGHT_NAME, false );

  this->prepare_formulas();

  // Set up storage for the "is_mc" boolean flag branch. If we're not working
  // with MC events, then we shouldn't do anything with the true bin counts.
  bool is_mc;
  input_chain_.SetBranchAddress( "is_mc", &is_mc );

  // Get the first TChain entry so that we can know the number of universes
  // used in each vector of weights
  input_chain_.GetEntry( 0 );

  // Now prepare the vectors of Universe objects with the correct sizes
  this->prepare_universes( wh );

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
      for ( auto& cbf : category_formulas_ ) cbf->Notify();
    }

    // Find the reco bin(s) that should be filled for the current event
    std::vector< FormulaMatch > matched_reco_bins;
    for ( size_t rb = 0u; rb < reco_bin_formulas_.size(); ++rb ) {
      auto& rbf = reco_bin_formulas_.at( rb );
      int num_formula_elements = rbf->GetNdata();
      for ( int el = 0; el < num_formula_elements; ++el ) {
        double formula_wgt = rbf->EvalInstance( el );
        if ( formula_wgt ) matched_reco_bins.emplace_back( rb, formula_wgt );
      }
    }

    // Find the EventCategory label(s) that apply to the current event
    std::vector< FormulaMatch > matched_category_indices;
    for ( size_t c = 0u; c < category_formulas_.size(); ++c ) {
      auto& cbf = category_formulas_.at( c );
      int num_formula_elements = cbf->GetNdata();
      for ( int el = 0; el < num_formula_elements; ++el ) {
        double formula_wgt = cbf->EvalInstance( el );
        if ( formula_wgt ) {
          matched_category_indices.emplace_back( c, formula_wgt );
        }
      }
    }

    input_chain_.GetEntry( entry );
    //std::cout << "Entry " << entry << '\n';

    std::vector< FormulaMatch > matched_true_bins;
    double spline_weight = 0.;
    double tune_weight = 0.;

    // If we're working with an MC sample, then find the true bin(s)
    // that should be filled for the current event
    if ( is_mc ) {
      for ( size_t tb = 0u; tb < true_bin_formulas_.size(); ++tb ) {
        auto& tbf = true_bin_formulas_.at( tb );
        int num_formula_elements = tbf->GetNdata();
        for ( int el = 0; el < num_formula_elements; ++el ) {
          double formula_wgt = tbf->EvalInstance( el );
          if ( formula_wgt ) matched_true_bins.emplace_back( tb, formula_wgt );
        }
      } // true bins

      // If we have event weights in the map at all, then get the current
      // event's CV correction weights here for potentially frequent re-use
      // below
      auto& wm = wh.weight_map();
      if ( wm.size() > 0u ) {
        spline_weight = wm.at( SPLINE_WEIGHT_NAME )->front();
        tune_weight = wm.at( TUNE_WEIGHT_NAME )->front();
      }
    } // MC event

    for ( const auto& pair : wh.weight_map() ) {
      const std::string& wgt_name = pair.first;
      const auto& wgt_vec = pair.second;

      auto& u_vec = universes_.at( wgt_name );

      for ( size_t u = 0u; u < wgt_vec->size(); ++u ) {

        // No need to use the slightly slower "at" here since we're directly
        // looping over the weight vector
        double w = wgt_vec->operator[]( u );

        // Multiply by any needed CV correction weights
        apply_cv_correction_weights( wgt_name, w, spline_weight, tune_weight );

        // Deal with NaNs, etc. to make a "safe weight" in all cases
        double safe_wgt = safe_weight( w );

        // Get the universe object that should be filled with the processed
        // event weight
        auto& universe = u_vec.at( u );

        for ( const auto& tb : matched_true_bins ) {
          // TODO: consider including the TTreeFormula weight(s) in the check
          // applied via safe_weight() above
          universe.hist_true_->Fill( tb.bin_index_, tb.weight_ * safe_wgt );
          for ( const auto& rb : matched_reco_bins ) {
            universe.hist_2d_->Fill( tb.bin_index_, rb.bin_index_,
              tb.weight_ * rb.weight_ * safe_wgt );
          } // reco bins
        } // true bins

        for ( const auto& rb : matched_reco_bins ) {
          universe.hist_reco_->Fill( rb.bin_index_, rb.weight_ * safe_wgt );

          for ( const auto& c : matched_category_indices ) {
            universe.hist_categ_->Fill( c.bin_index_, rb.bin_index_,
              c.weight_ * rb.weight_ * safe_wgt );
          }

          for ( const auto& other_rb : matched_reco_bins ) {
            universe.hist_reco2d_->Fill( rb.bin_index_, other_rb.bin_index_,
              rb.weight_ * other_rb.weight_ * safe_wgt );
          }
        } // reco bins
      } // universes
    } // weight names

    // Fill the unweighted histograms now that we're done with the
    // weighted ones. Note that "unweighted" in this context applies to
    // the universe event weights, but that any implicit weights from
    // the TTreeFormula evaluations will still be applied.
    auto& univ = universes_.at( UNWEIGHTED_NAME ).front();
    for ( const auto& tb : matched_true_bins ) {
      univ.hist_true_->Fill( tb.bin_index_, tb.weight_ );
      for ( const auto& rb : matched_reco_bins ) {
        univ.hist_2d_->Fill( tb.bin_index_, rb.bin_index_,
          tb.weight_ * rb.weight_ );
      } // reco bins
    } // true bins

    for ( const auto& rb : matched_reco_bins ) {

      univ.hist_reco_->Fill( rb.bin_index_, rb.weight_ );

      for ( const auto& c : matched_category_indices ) {
        univ.hist_categ_->Fill( c.bin_index_, rb.bin_index_,
          c.weight_ * rb.weight_ );
      }

      for ( const auto& other_rb : matched_reco_bins ) {
        univ.hist_reco2d_->Fill( rb.bin_index_, other_rb.bin_index_,
          rb.weight_ * other_rb.weight_ );
      }

    } // reco bins

  } // TChain entries

  input_chain_.ResetBranchAddresses();
}

void UniverseMaker::prepare_universes( const WeightHandler& wh ) {

  size_t num_true_bins = true_bins_.size();
  size_t num_reco_bins = reco_bins_.size();

  for ( const auto& pair : wh.weight_map() ) {
    const std::string& weight_name = pair.first;
    size_t num_universes = pair.second->size();

    std::vector< Universe > u_vec;

    for ( size_t u = 0u; u < num_universes; ++u ) {
      u_vec.emplace_back( weight_name, u, num_true_bins, num_reco_bins );
    }

    universes_[ weight_name ] = std::move( u_vec );
  }

  // Add the special "unweighted" universe unconditionally
  std::vector< Universe > temp_uvec;
  temp_uvec.emplace_back( UNWEIGHTED_NAME, 0, num_true_bins, num_reco_bins );
  universes_[ UNWEIGHTED_NAME ] = std::move( temp_uvec );

}

void UniverseMaker::save_histograms(
  const std::string& output_file_name,
  const std::string& subdirectory_name,
  bool update_file )
{
  // Decide whether to overwrite the output file or simply update the contents.
  // This difference is only important if the output file already exists before
  // this function is called.
  std::string tfile_option( "recreate" );
  if ( update_file ) {
    tfile_option = "update";
  }

  TFile out_file( output_file_name.c_str(), tfile_option.c_str() );

  // Navigate to the subdirectory within the output ROOT file where the
  // response matrix histograms will be saved. Create new TDirectoryFile
  // objects as needed.
  TDirectoryFile* root_tdir = nullptr;
  TDirectoryFile* sub_tdir = nullptr;

  out_file.GetObject( output_directory_name_.c_str(), root_tdir );
  if ( !root_tdir ) {
    // TODO: add error handling for a forward slash in the root TDirectoryFile
    // name
    root_tdir = new TDirectoryFile( output_directory_name_.c_str(),
      "response matrices", "", &out_file );
  }

  // Save the configuration settings for this class to the main
  // TDirectoryFile before moving on to the appropriate subdirectory. If
  // these settings have already been saved, then double-check that they
  // match the current configuration. In the event of a mismatch, throw an
  // exception to avoid data corruption.
  std::string tree_name, true_bin_spec, reco_bin_spec;
  tree_name = this->input_chain().GetName();

  std::ostringstream oss_true, oss_reco;

  for ( const auto& tbin : true_bins_ ) {
    oss_true << tbin << '\n';
  }

  for ( const auto& rbin : reco_bins_ ) {
    oss_reco << rbin << '\n';
  }

  true_bin_spec = oss_true.str();
  reco_bin_spec = oss_reco.str();

  std::string* saved_tree_name = nullptr;
  std::string* saved_tb_spec = nullptr;
  std::string* saved_rb_spec = nullptr;
  root_tdir->GetObject( "ntuple_name", saved_tree_name );
  root_tdir->GetObject( TRUE_BIN_SPEC_NAME.c_str(), saved_tb_spec );
  root_tdir->GetObject( RECO_BIN_SPEC_NAME.c_str(), saved_rb_spec );

  if ( saved_tree_name ) {
    if ( tree_name != *saved_tree_name ) {
      throw std::runtime_error( "Tree name mismatch: " + tree_name
        + " vs. " + *saved_tree_name );
    }
  }
  else {
    root_tdir->WriteObject( &tree_name, "ntuple_name" );
  }

  if ( saved_tb_spec ) {
    if ( true_bin_spec != *saved_tb_spec ) {
      throw std::runtime_error( "Inconsistent true bin specification!" );
    }
  }
  else {
    root_tdir->WriteObject( &true_bin_spec, TRUE_BIN_SPEC_NAME.c_str() );
  }

  if ( saved_rb_spec ) {
    if ( reco_bin_spec != *saved_rb_spec ) {
      throw std::runtime_error( "Inconsistent reco bin specification!" );
    }
  }
  else {
    root_tdir->WriteObject( &reco_bin_spec, RECO_BIN_SPEC_NAME.c_str() );
  }

  std::string subdir_name = ntuple_subfolder_from_file_name(
    subdirectory_name );

  root_tdir->GetObject( subdir_name.c_str(), sub_tdir );
  if ( !sub_tdir ) {
    sub_tdir = new TDirectoryFile( subdir_name.c_str(), "response matrices",
      "", root_tdir );
  }

  // Now we've found (or created) the TDirectoryFile where the output
  // will be saved. Ensure that it is the active file here before writing
  // out the histograms.
  sub_tdir->cd();

  for ( auto& pair : universes_ ) {
    auto& u_vec = pair.second;
    for ( auto& univ : u_vec ) {
      // Always save the reco histograms
      univ.hist_reco_->Write();
      univ.hist_reco2d_->Write();

      // Save the others if the true histogram was filled at least once
      // (used to infer that we have MC truth information)
      if ( univ.hist_true_->GetEntries() > 0. ) {
        univ.hist_true_->Write();
        univ.hist_2d_->Write();
        univ.hist_categ_->Write();
      }
    } // universes
  } // weight names
}
