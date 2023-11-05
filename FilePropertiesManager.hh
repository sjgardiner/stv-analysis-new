#pragma once

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

// Enum used to label types of analysis ntuple files
enum class NtupleFileType {

  // Data taken with the beam on
  kOnBNB,

  // Data taken with the beam off
  kExtBNB,

  // MC events (standard BNB)
  kNumuMC,

  // MC events (intrinsic-nue-enhanced BNB)
  kIntrinsicNueMC,

  // MC events (dirt BNB)
  kDirtMC,

  // *** DetVar MC ***
  kDetVarMCCV, // central value
  kDetVarMCLYatten, // light-yield attenuation
  kDetVarMCLYdown, // light-yield down
  kDetVarMCLYrayl, // light-yield Rayleigh scattering
  kDetVarMCRecomb2, // light-yield recombination 2
  kDetVarMCSCE, // space charge effect
  kDetVarMCWMAngleXZ, // wireMod angle XZ
  kDetVarMCWMAngleYZ, // wireMod angle YZ
  kDetVarMCWMdEdx, // wireMod dE/dx
  kDetVarMCWMX, // wireMod X
  kDetVarMCWMYZ, // wireMod YZ
  kDetVarMCCVExtra, // alternate CV for small samples

  // An alternate CV MC simulation
  kAltCVMC,

  // Placeholder for invalid values
  kUnknown,
};

// Utility functions for manipulating NtupleFileType values
bool ntuple_type_is_detVar( const NtupleFileType& type ) {
  constexpr std::array< NtupleFileType, 12 > detVar_types = {
    NtupleFileType::kDetVarMCCV, NtupleFileType::kDetVarMCLYatten,
    NtupleFileType::kDetVarMCLYdown, NtupleFileType::kDetVarMCLYrayl,
    NtupleFileType::kDetVarMCRecomb2, NtupleFileType::kDetVarMCSCE,
    NtupleFileType::kDetVarMCWMAngleXZ, NtupleFileType::kDetVarMCWMAngleYZ,
    NtupleFileType::kDetVarMCWMdEdx, NtupleFileType::kDetVarMCWMX,
    NtupleFileType::kDetVarMCWMYZ, NtupleFileType::kDetVarMCCVExtra
  };

  const auto begin = detVar_types.cbegin();
  const auto end = detVar_types.cend();
  const auto iter = std::find( begin, end, type );
  if ( iter != end ) return true;
  return false;
}

bool ntuple_type_is_altCV( const NtupleFileType& type ) {
  if ( type == NtupleFileType::kAltCVMC ) return true;
  return false;
}

bool ntuple_type_is_mc( const NtupleFileType& type ) {
  if ( type != NtupleFileType::kOnBNB && type != NtupleFileType::kExtBNB ) {
    return true;
  }
  return false;
}

bool ntuple_type_is_reweightable_mc( const NtupleFileType& type ) {

  if ( type == NtupleFileType::kNumuMC
    || type == NtupleFileType::kIntrinsicNueMC
    || type == NtupleFileType::kDirtMC )
  {
    return true;
  }

  return false;
}

// Singleton class that keeps track of the various ntuple files to be analyzed
class FilePropertiesManager {

  public:

    // This is a singleton class, so we'll delete the copy constructor
    // the move constructor, and the assignment operators
    FilePropertiesManager( const FilePropertiesManager& ) = delete;
    FilePropertiesManager( FilePropertiesManager&& ) = delete;
    FilePropertiesManager& operator=( const FilePropertiesManager& )
      = delete;
    FilePropertiesManager& operator=( FilePropertiesManager&& )
      = delete;

    // Get a const reference to the singleton instance of the
    // FilePropertiesManager
    inline static FilePropertiesManager& Instance() {

      // Create the FilePropertiesManager object using a static variable.
      // This ensures that the singleton instance is only created once.
      static std::unique_ptr<FilePropertiesManager>
        the_instance( new FilePropertiesManager() );

      // Return a reference to the singleton instance
      return *the_instance;
    }

    // Simple container that stores the number of triggers and the POT exposure
    // represented by a particular data ntuple
    struct TriggersAndPOT {
      TriggersAndPOT() {}
      TriggersAndPOT( int trig, double pot )
        : trigger_count_( trig ), pot_( pot ) {}

      int trigger_count_ = 0;
      double pot_ = 0.;
    };

    inline const std::map< std::string, TriggersAndPOT >& data_norm_map() const
      { return data_norm_map_; }

    inline const std::map< int, std::map<NtupleFileType,
      std::set<std::string> > >& ntuple_file_map() const
      { return ntuple_file_map_; }

    // NOTE: This function assumes that there are no duplicate entries in the
    // nested ntuple file type map. That should always be the case, but beware
    // any odd situations that might break the approach used here.
    inline NtupleFileType get_ntuple_file_type(
      const std::string& file_name ) const
    {
      for ( const auto& run_pair : ntuple_file_map_ ) {
        const auto& type_map = run_pair.second;
        for ( const auto& type_pair : type_map ) {
          const auto& file_set = type_pair.second;
          for ( const auto& name : file_set ) {
            if ( file_name == name ) return type_pair.first;
          }
        }
      }
      throw std::runtime_error( "ntuple file not found" );
      // Bogus value to satisfy the function signature requirements
      return NtupleFileType::kOnBNB;
    }

    // Returns a string representation of an NtupleFileType value, or
    // an empty string if one could not be found.
    std::string ntuple_type_to_string( const NtupleFileType& type ) const {
      std::string result;
      for ( const auto& pair : string_to_file_type_map_ ) {
        if ( pair.second == type ) {
          result = pair.first;
          break;
        }
      }
      return result;
    }

    // Converts a string into an NtupleFileType value. Invalid strings will
    // yield the value NtupleFileType::kUnknown.
    NtupleFileType string_to_ntuple_type( const std::string& str ) const {
      auto end = string_to_file_type_map_.cend();
      auto iter = string_to_file_type_map_.find( str );
      if ( iter != end ) return iter->second;
      return NtupleFileType::kUnknown;
    }

    inline const std::string& analysis_path() const { return analysis_path_; }

    inline void load_file_properties(
      const std::string& input_table_file_name = "" )
    {

      // Clear out any pre-existing contents of the owned maps storing
      // analysis ntuple file properties
      ntuple_file_map_.clear();
      data_norm_map_.clear();

      const char* path = std::getenv( "STV_ANALYSIS_DIR" );
      if ( path == nullptr ) throw std::runtime_error( "The environment"
        " variable STV_ANALYSIS_DIR is not set. Please set it and try again." );

      analysis_path_ = path;

      // If the user didn't manually specify a table of file properties, then
      // use the default one
      std::string in_file_name( input_table_file_name );
      if ( in_file_name.empty() ) {
        in_file_name = analysis_path_ + "/file_properties.txt";
        //std::cout << "Using default file properties which is " << in_file_name << std::endl;
      }

      // Store the name of the configuration file that was used so that (if
      // needed) we can retrieve it later
      config_file_name_ = in_file_name;

      //std::cout << "Loading configuration file " << config_file_name_ << std::endl;

      std::ifstream in_file( in_file_name );

      if ( !in_file ) {
        throw std::runtime_error( "The file properties configuration file \""
          + in_file_name + "\" could not be opened." );
      }

      std::string temp_line;
      while ( std::getline(in_file, temp_line) ) {
        // Ignore lines that begin with the '#' character (this allows for
        // comments in the normalization table file
        if ( temp_line.front() == '#' ) continue;

        // Read in the ntuple file name, the run number, and the file type from
        // the current line of the table file
        std::string file_name;
        int run;
        std::string type_str;
        std::istringstream temp_ss( temp_line );
        temp_ss >> file_name >> run >> type_str;

        // Convert the type string into an enum class value
        NtupleFileType type = string_to_file_type_map_.at( type_str );

        // If there is not an inner map for the current run number, then create
        // one
        if ( !ntuple_file_map_.count(run) ) {
          ntuple_file_map_[ run ] = std::map< NtupleFileType,
            std::set<std::string> >();
        }
        auto& run_map = ntuple_file_map_.at( run );

        // If there is not a set of files for the current ntuple file type,
        // then create one
        if ( !run_map.count(type) ) {
          run_map[ type ] = std::set< std::string >();
        }
        auto& file_set = run_map.at( type );

        // Insert the current file name into the appropriate place in the map
        file_set.insert( file_name );

        // For data files, also read in the trigger count and POT exposure
        // needed for normalization purposes
        if ( type == NtupleFileType::kOnBNB
          || type == NtupleFileType::kExtBNB )
        {
          int trigs;
          double pot;
          temp_ss >> trigs >> pot;

          // Store this information in the normalization map
          data_norm_map_[ file_name ] = TriggersAndPOT( trigs, pot );
        }
      }
    }

    inline const std::string& config_file_name() const
      { return config_file_name_; }

  private:

    inline FilePropertiesManager() {
      this->load_file_properties();
    }

    // Outer keys are run numbers, inner keys are ntuple file types, values are
    // sets of file names
    std::map< int, std::map<NtupleFileType, std::set<std::string> > > ntuple_file_map_;

    // Keys are file names for processed data ntuples, values are objects
    // storing the corresponding number of triggers and POT exposure
    std::map< std::string, TriggersAndPOT > data_norm_map_;

    std::map< std::string, NtupleFileType > string_to_file_type_map_ = {
      { "onBNB", NtupleFileType::kOnBNB },
      { "extBNB", NtupleFileType::kExtBNB },
      { "numuMC", NtupleFileType::kNumuMC },
      { "nueMC", NtupleFileType::kIntrinsicNueMC },
      { "dirtMC", NtupleFileType::kDirtMC },
      { "detVarCV", NtupleFileType::kDetVarMCCV },
      { "detVarLYatten", NtupleFileType::kDetVarMCLYatten },
      { "detVarLYdown", NtupleFileType::kDetVarMCLYdown },
      { "detVarLYrayl", NtupleFileType::kDetVarMCLYrayl },
      { "detVarRecomb2", NtupleFileType::kDetVarMCRecomb2 },
      { "detVarSCE", NtupleFileType::kDetVarMCSCE },
      { "detVarWMAngleXZ", NtupleFileType::kDetVarMCWMAngleXZ },
      { "detVarWMAngleYZ", NtupleFileType::kDetVarMCWMAngleYZ },
      { "detVarWMdEdx", NtupleFileType::kDetVarMCWMdEdx },
      { "detVarWMX", NtupleFileType::kDetVarMCWMX },
      { "detVarWMYZ", NtupleFileType::kDetVarMCWMYZ },
      { "detVarCVExtra", NtupleFileType::kDetVarMCCVExtra },
      { "altCVMC", NtupleFileType::kAltCVMC },
    };

    // Folder that stores the STV analysis configuration files. This is set
    // automatically on construction using the STV_ANALYSIS_DIR environment
    // variable.
    std::string analysis_path_;

    // Name of the file used to configure the singleton class
    std::string config_file_name_;
};
