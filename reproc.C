#include "analyzer.C"

const std::string output_folder = "/uboone/data/users/gardiner/ntuples-stv/";

std::string base_filename( const std::string& filename ) {
  return filename.substr( filename.find_last_of("/") + 1 );
}

void reproc() {

  std::ifstream in_file( "files_to_process.txt" );

  std::string file_to_process;
  while ( std::getline(in_file, file_to_process) ) {
    if ( file_to_process.front() == '#' || file_to_process.empty() ) continue;
    std::cout << "PROCESSING " << file_to_process << '\n';
    std::string out_filename = output_folder + "stv-"
      + base_filename( file_to_process );

    analyzer( file_to_process, out_filename );
  }

}
