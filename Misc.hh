#ifndef _Misc_hh_
#define _Misc_hh_

#include <sstream>

template <typename T> std::string to_string_with_precision(const T a_value, const int n = 6)
{
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return std::move(out).str();
}

////////////////////////////////////////////////////////////////////////////////////////

void MakeBinSliceConfig(string label,string signal_def,string sel_def,vector<string> true_variables,vector<string> reco_variables,vector<vector<double>> true_binning,vector<vector<double>> reco_binning){

  // Check various things

  if(true_variables.size() != reco_variables.size()) 
    throw std::invalid_argument("true_variables.size() != reco_variables.size()");

  if(true_binning.size() != reco_binning.size())
    throw std::invalid_argument("true_binning and reco_binning have different sizes");

  for(size_t i_var=0;i_var<true_binning.size();i_var++){
    if(true_binning.at(i_var).size() != reco_binning.at(i_var).size())
      throw std::invalid_argument("true binning and reco binning different numbers of bins for variable " + std::to_string(i_var));
  }

  // Make bin config file

  std::ofstream output("bin_config_" + label + ".txt"); 
  output << "tutorial_1D_proton" << std::endl;
  output << "stv_tree" << std::endl;

  // count the number of bins
  int true_bins = 0;
  for(size_t i_var=0;i_var<true_variables.size();i_var++){
    true_bins += true_binning.at(i_var).size()-1;
    for(size_t j_var=0;j_var<true_variables.size();j_var++){ 
      if(i_var == j_var) continue;
      true_bins += (true_binning.at(i_var).size()-1)*(true_binning.at(j_var).size()-1);
    }
  }

  output << true_bins + 7 << std::endl;

  int bl_index=0;

  // truth binning for the 1D distributions 
  for(size_t i_var=0;i_var<true_variables.size();i_var++){
    for(size_t i_bin=0;i_bin<true_binning.at(i_var).size()-1;i_bin++){
      output << "0 " << bl_index << " " << "\"" << signal_def << " && " << true_variables.at(i_var) << " >= " << true_binning.at(i_var).at(i_bin) << " && " << true_variables.at(i_var) << " < " << true_binning.at(i_var).at(i_bin+1) << "\"" << std::endl;
    }
    bl_index++;
  } 

  // truth binning for the 2D distributions
  if(true_binning.size() > 1){
    for(size_t i_var=0;i_var<true_variables.size();i_var++){
      for(size_t i_bin=0;i_bin<true_binning.at(i_var).size()-1;i_bin++){
        for(size_t j_var=0;j_var<true_variables.size();j_var++){
          if(i_var == j_var) continue;
          for(size_t j_bin=0;j_bin<true_binning.at(j_var).size()-1;j_bin++){
            output << "0 " << bl_index << " " << "\"" << signal_def << " && " << true_variables.at(i_var) << " >= " << true_binning.at(i_var).at(i_bin) << " && " << true_variables.at(i_var) << " < " << true_binning.at(i_var).at(i_bin+1); 
            output << " && " << true_variables.at(j_var) << " >= " << true_binning.at(j_var).at(j_bin) << " && " << true_variables.at(j_var) << " < " << true_binning.at(j_var).at(j_bin+1) << "\"" << std::endl;
          }
          bl_index++;
        } 
      }
    }
  }

  // Add the background bins
  output << "1 -1 \"category == 5\"" << std::endl;
  output << "1 -1 \"category == 6\"" << std::endl;
  output << "1 -1 \"category == 7\"" << std::endl;
  output << "1 -1 \"category == 8\"" << std::endl;
  output << "1 -1 \"category == 9\"" << std::endl;
  output << "1 -1 \"category == 10\"" << std::endl;
  output << "1 -1 \"category == 11\"" << std::endl;

  // Do the same for the reco bins
  output << true_bins << std::endl;
  bl_index=0;

  // start with the 1D part
  for(size_t i_var=0;i_var<reco_variables.size();i_var++){
    for(size_t i_bin=0;i_bin<reco_binning.at(i_var).size()-1;i_bin++){
      output << "0 " << bl_index << " " << "\"" << signal_def << " && " << reco_variables.at(i_var) << " >= " << reco_binning.at(i_var).at(i_bin) << " && " << reco_variables.at(i_var) << " < " << reco_binning.at(i_var).at(i_bin+1) << "\"" << std::endl;
    }
    bl_index++;
  } 

  // the 2D part
  if(reco_binning.size() > 1){
    for(size_t i_var=0;i_var<reco_variables.size();i_var++){
      for(size_t i_bin=0;i_bin<reco_binning.at(i_var).size()-1;i_bin++){
        for(size_t j_var=0;j_var<reco_variables.size();j_var++){
          if(i_var == j_var) continue;
          for(size_t j_bin=0;j_bin<reco_binning.at(j_var).size()-1;j_bin++){
            output << "0 " << bl_index << " " << "\"" << sel_def << " && " << reco_variables.at(i_var) << " >= " << reco_binning.at(i_var).at(i_bin) << " && " << reco_variables.at(i_var) << " < " << reco_binning.at(i_var).at(i_bin+1); 
            output << " && " << reco_variables.at(j_var) << " >= " << reco_binning.at(j_var).at(j_bin) << " && " << reco_variables.at(j_var) << " < " << reco_binning.at(j_var).at(j_bin+1) << "\"" << std::endl;
          }
          bl_index++;
        } 
      }
    }
  }


  // Make the accompanying slice config

  std::ofstream slice_output("slice_config_" + label + ".txt"); 
  slice_output << bl_index << std::endl; 

  for(size_t i_var=0;i_var<reco_variables.size();i_var++)
    slice_output << "\"" << reco_variables.at(i_var) << "\" \"\" \"" << reco_variables.at(i_var) << "\" \"\"" << std::endl;

  for(size_t i_var=0;i_var<reco_variables.size();i_var++){
    for(size_t i_bin=0;i_bin<reco_binning.at(i_var).size()-1;i_bin++){
      for(size_t j_var=0;j_var<true_variables.size();j_var++){
        if(i_var == j_var) continue;
        slice_output << "\"" << reco_variables.at(i_var) << "_" << i_bin << "_" << reco_variables.at(j_var)  << "\" \"\" \"" << reco_variables.at(i_var) << "_" << i_bin << "_" << reco_variables.at(j_var) << "\" \"\"" << std::endl;
      }
    }
  }

  slice_output << bl_index << std::endl;
  bl_index = 0;
  int bin = 0;

  // the 1D slices
  for(size_t i_var=0;i_var<reco_variables.size();i_var++){
    slice_output << "\"events\"" << std::endl;
    slice_output << 1 << " " << bl_index << " " << reco_binning.at(i_var).size() << " "; 
    for(size_t i_bin=0;i_bin<reco_binning.at(i_var).size();i_bin++)
      slice_output << reco_binning.at(i_var).at(i_bin) << " ";
    slice_output << std::endl;
    slice_output << 0 << std::endl;
    slice_output << reco_binning.at(i_var).size()-1 << std::endl;
    for(size_t i_bin=0;i_bin<reco_binning.at(i_var).size()-1;i_bin++){
      slice_output << bin << " " << 1 << " " << i_bin+1 << std::endl; 
      bin++;
    } 
    bl_index++;
  }

  // the 2D slices
  for(size_t i_var=0;i_var<reco_variables.size();i_var++){
    for(size_t i_bin=0;i_bin<reco_binning.at(i_var).size()-1;i_bin++){
      for(size_t j_var=0;j_var<reco_variables.size();j_var++){
        if(i_var == j_var) continue;
        slice_output << "\"events\"" << std::endl;
        slice_output << 1 << " " << bl_index << " " << reco_binning.at(j_var).size() << " "; 
        for(size_t j_bin=0;j_bin<reco_binning.at(j_var).size();j_bin++)
          slice_output << reco_binning.at(j_var).at(j_bin) << " ";
        slice_output << std::endl;
        slice_output << 0 << std::endl;
        slice_output << reco_binning.at(j_var).size()-1 << std::endl;
        for(size_t j_bin=0;j_bin<reco_binning.at(j_var).size()-1;j_bin++){
          slice_output << bin << " " << 1 << " " << j_bin+1 << std::endl; 
          bin++;
        } 
        bl_index++;
      }
    }
  }

}

#endif
