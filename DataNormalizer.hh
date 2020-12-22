#pragma once

#include <cstdlib>
#include <fstream>
#include <map>
#include <stdexcept>

// Singleton class that keeps track of trigger counts and POT exposure for
// data ntuples
class DataNormalizer {

  public:

    // This is a singleton class, so we'll delete the copy constructor
    // the move constructor, and the assignment operators
    DataNormalizer( const DataNormalizer& ) = delete;
    DataNormalizer( DataNormalizer&& ) = delete;
    DataNormalizer& operator=( const DataNormalizer& )
      = delete;
    DataNormalizer& operator=( DataNormalizer&& )
      = delete;

    // Get a const reference to the singleton instance of the
    // DataNormalizer
    inline static const DataNormalizer& Instance() {

      // Create the DataNormalizer object using a static variable.
      // This ensures that the singleton instance is only created once.
      static std::unique_ptr<DataNormalizer>
        the_instance( new DataNormalizer() );

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

    inline const std::map< std::string, TriggersAndPOT >& norm_map() const
      { return norm_map_; }

  private:

    inline DataNormalizer() {
      this->load_norm_map();
    }

    inline void load_norm_map() {

      const char* path = std::getenv( "STV_ANALYSIS_DIR" );
      if ( path == nullptr ) throw std::runtime_error( "The environment"
        " variable STV_ANALYSIS_DIR is not set. Please set it and try again." );

      std::string in_file_name( path );
      in_file_name += "/data_normalization.txt";
      std::ifstream in_file( in_file_name );

      std::string temp_line;
      while ( std::getline(in_file, temp_line) ) {
        // Ignore lines that begin with the '#' character (this allows for
        // comments in the normalization table file
        if ( temp_line.front() == '#' ) continue;

        // Read in the file basename, the trigger count, and the POT exposure
        // from the current line of the table file
        std::string file_basename;
        int trigs;
        double pot;
        std::istringstream temp_ss( temp_line );
        temp_ss >> file_basename >> trigs >> pot;

        // Store this information in the normalization map
        norm_map_[ file_basename ] = TriggersAndPOT( trigs, pot );
      }
    }

    // Keys are basenames for processed data ntuples, values are objects
    // storing the corresponding number of triggers and POT exposure
    std::map< std::string, TriggersAndPOT > norm_map_;

};
