#ifndef _BinConfigFuncs_hh_
#define _BinConfigFuncs_hh_

#include <sstream>

// Tools for automatically generating new binning configurations

using std::vector;
using std::pair;
using std::map;

const double stat_err_max = 0.032; // maximum acceptable st
const double capture_frac = 0.995; // minimum fraction of dist to contain when setting limits on any variables

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Generates bin/slice config for simple 1D distribution

void MakeBinSliceConfig1D(string label,string signal_def,string sel_def,string true_variable,string reco_variable,vector<double> true_binning,vector<double> reco_binning){

  system("mkdir -p bin_configs/");
  system("mkdir -p slice_configs/");

  if(true_binning.size() != reco_binning.size())
    throw std::invalid_argument("true_binning and reco_binning have different sizes");

  std::ofstream output_1D("bin_configs/bin_config_1D_" + label + ".txt"); 
  output_1D << "tutorial_1D_proton" << std::endl;
  output_1D << "stv_tree" << std::endl;
  output_1D << true_binning.size()-1+7 << std::endl; 

  int global_bin = 0;
  for(size_t i_bin=0;i_bin<true_binning.size()-1;i_bin++){
    output_1D << "0 " << 0 << " " << "\"" << signal_def << " && " << true_variable << " >= " << true_binning.at(i_bin) << " && " << true_variable << " < " << true_binning.at(i_bin+1) << "\"" << std::endl;
    global_bin++;
  }

  output_1D << "1 -1 \"category == 5\"" << std::endl;
  output_1D << "1 -1 \"category == 6\"" << std::endl;
  output_1D << "1 -1 \"category == 7\"" << std::endl;
  output_1D << "1 -1 \"category == 8\"" << std::endl;
  output_1D << "1 -1 \"category == 9\"" << std::endl;
  output_1D << "1 -1 \"category == 10\"" << std::endl;
  output_1D << "1 -1 \"category == 11\"" << std::endl;

  output_1D << true_binning.size()-1 << std::endl;  

  for(size_t i_bin=0;i_bin<reco_binning.size()-1;i_bin++){
    output_1D << "0 " << 0 << " " << "\"" << sel_def << " && " << reco_variable << " >= " << reco_binning.at(i_bin) << " && " << reco_variable << " < " << reco_binning.at(i_bin+1) << "\"" << std::endl;
  }

  std::ofstream slice_output_1D("slice_configs/slice_config_1D_" + label + ".txt"); 

  slice_output_1D << 1 << std::endl; 
  slice_output_1D << "\"" << true_variable << "\" \"\" \"" << true_variable << "\" \"\"" << std::endl;
  slice_output_1D << 1 << std::endl; 

  slice_output_1D << "\"events\"" << std::endl;
  slice_output_1D << 1 << " " << 0 << " " << reco_binning.size() << " "; 
  for(size_t i_bin=0;i_bin<reco_binning.size();i_bin++)
    slice_output_1D << reco_binning.at(i_bin) << " ";
  slice_output_1D << std::endl;
  slice_output_1D << 0 << std::endl;
  slice_output_1D << reco_binning.size()-1 << std::endl;
  for(size_t i_bin=0;i_bin<reco_binning.size()-1;i_bin++){
    slice_output_1D << i_bin << " " << 1 << " " << i_bin+1 << std::endl; 
  } 

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Generates bin/slice config for 2D cross section
// For true_binning/reco_binning, use the output from GetBinning2D below 

void MakeBinSliceConfig2D(string label,string signal_def,string sel_def,vector<string> true_variables,vector<string> reco_variables,vector<pair<double,vector<double>>> true_binning,vector<pair<double,vector<double>>> reco_binning){

  system("mkdir -p bin_configs/");
  system("mkdir -p slice_configs/");

  if(true_variables.size() != 2 || reco_variables.size() != 2)
    throw std::invalid_argument("GetBinning2D: Only supported for two variables at the moment");

  if(true_variables.size() != reco_variables.size()) 
    throw std::invalid_argument("true_variables.size() != reco_variables.size()");

  if(true_binning.size() != reco_binning.size())
    throw std::invalid_argument("true_binning and reco_binning have different sizes");

  // count the total number of bins
  int bins = 0;
  int blocks = true_binning.size()-1;
  for(size_t i_x=0;i_x<true_binning.size()-1;i_x++)
    bins += true_binning.at(i_x).second.size()-1;

  std::ofstream output_2D("bin_configs/bin_config_2D_" + label + ".txt"); 
  output_2D << "tutorial_1D_proton" << std::endl;
  output_2D << "stv_tree" << std::endl;
  output_2D << bins+7 << std::endl; 

  std::ofstream slice_output_2D("slice_configs/slice_config_2D_" + label + ".txt"); 
  slice_output_2D << 2 << std::endl; // number of variables 

  // two variables 
  slice_output_2D << "\"" << true_variables.at(0) << "\" \"\" \"" << true_variables.at(0) << "\" \"\"" << std::endl;
  slice_output_2D << "\"" << true_variables.at(1) << "\" \"\" \"" << true_variables.at(1) << "\" \"\"" << std::endl;
  // number of slices = number of bins in lead variable + 1 for the integral
  slice_output_2D << true_binning.size() << std::endl; // 

  int bl_index = 0;
  int glbl_bin_no = 0;
  for(int i_block=0;i_block<blocks;i_block++){

    double x_limit_low = true_binning.at(i_block).first;
    double x_limit_high = true_binning.at(i_block+1).first;

    slice_output_2D << "\"events\"" << std::endl;
    slice_output_2D << 1 << " " << 1 << " " << true_binning.at(i_block).second.size() << " "; 
    for(size_t i_bin=0;i_bin<true_binning.at(i_block).second.size();i_bin++)
      slice_output_2D << true_binning.at(i_block).second.at(i_bin) << " ";

    slice_output_2D << std::endl;
    slice_output_2D << 1 << std::endl;
    slice_output_2D << "0 " << x_limit_low << " " << x_limit_high << std::endl; 

    slice_output_2D << reco_binning.at(i_block).second.size()-1 << std::endl;

    for(int i_bin=0;i_bin<true_binning.at(i_block).second.size()-1;i_bin++){

      double y_limit_low = true_binning.at(i_block).second.at(i_bin);
      double y_limit_high = true_binning.at(i_block).second.at(i_bin+1);

      output_2D << "0 " << bl_index << " " << "\"" << signal_def << " && " << true_variables.at(0) << " >= " << x_limit_low << " && " << true_variables.at(0) << " < " << x_limit_high; 
      output_2D << " && " << true_variables.at(1) << " >= " << y_limit_low << " && " << true_variables.at(1) << " < " << y_limit_high << "\"" << std::endl;
      slice_output_2D << glbl_bin_no << " " << 1 << " " << i_bin+1 << std::endl;
      glbl_bin_no++;
    }
    bl_index++;
  }

  output_2D << "1 -1 \"category == 5\"" << std::endl;
  output_2D << "1 -1 \"category == 6\"" << std::endl;
  output_2D << "1 -1 \"category == 7\"" << std::endl;
  output_2D << "1 -1 \"category == 8\"" << std::endl;
  output_2D << "1 -1 \"category == 9\"" << std::endl;
  output_2D << "1 -1 \"category == 10\"" << std::endl;
  output_2D << "1 -1 \"category == 11\"" << std::endl;

  output_2D << bins << std::endl; 

  bl_index = 0;
  for(int i_block=0;i_block<blocks;i_block++){
    double x_limit_low = reco_binning.at(i_block).first;
    double x_limit_high = reco_binning.at(i_block+1).first;
    for(int i_bin=0;i_bin<reco_binning.at(i_block).second.size()-1;i_bin++){
      double y_limit_low = reco_binning.at(i_block).second.at(i_bin);
      double y_limit_high = reco_binning.at(i_block).second.at(i_bin+1);
      output_2D << "0 " << bl_index << " " << "\"" << sel_def << " && " << reco_variables.at(0) << " >= " << x_limit_low << " && " << reco_variables.at(0) << " < " << x_limit_high; 
      output_2D << " && " << reco_variables.at(1) << " >= " << y_limit_low << " && " << reco_variables.at(1) << " < " << y_limit_high << "\"" << std::endl;
    }
    bl_index++;
  }

  output_2D.close();

  // 1D histogram by integrating over the 2nd variable

  slice_output_2D << "\"events\"" << std::endl;
  slice_output_2D << 1 << " " << 0 << " " << true_binning.size() << " "; 
  for(int i_block=0;i_block<true_binning.size();i_block++)
    slice_output_2D << reco_binning.at(i_block).first << " ";
   slice_output_2D << std::endl;
   slice_output_2D << 0 << std::endl;
   slice_output_2D << glbl_bin_no << std::endl;
      
   glbl_bin_no = 0;
   for(int i_block=0;i_block<true_binning.size()-1;i_block++){
     for(int i_bin=0;i_bin<true_binning.at(i_block).second.size()-1;i_bin++){
       slice_output_2D << glbl_bin_no << " " << 1 << " " << i_block+1 << std::endl; 
       glbl_bin_no++;
     }
   }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Generates a binning scheme in which each bin has an fractional MC stat error < stat_err
// Histogram holds the original distribution

vector<double> GetBinning1D_hist(TH1F* h_tmp,double stat_err = stat_err_max){

  double tot = h_tmp->Integral();

  int bin_upper_lim = h_tmp->GetNbinsX();
  double frac = (double)h_tmp->GetBinContent(bin_upper_lim)/tot; 
  while(frac < (1.0-capture_frac)/2){
    bin_upper_lim--;  
    frac += (double)h_tmp->GetBinContent(bin_upper_lim)/tot;
  }

  int bin_lower_lim = 1;
  frac = (double)h_tmp->GetBinContent(bin_lower_lim)/tot; 
  while(frac < (1.0-capture_frac)/2){
    bin_lower_lim++;  
    frac += (double)h_tmp->GetBinContent(bin_lower_lim)/tot;
  }

  int bins_set = 0;
  std::vector<int> bin_boundaries;
  bin_boundaries.push_back(bin_lower_lim);
  int bin = bin_lower_lim;
  double current_set_content = h_tmp->GetBinContent(bin);
  while(bin < bin_upper_lim){
    double current_set_content = 0;
    double frac_unc = 1.0;

    while(frac_unc > stat_err && bin < bin_upper_lim){
      bin++;
      current_set_content += h_tmp->GetBinContent(bin);
      if(current_set_content > 0)
        frac_unc = sqrt((double)current_set_content)/current_set_content;
    }

    bin_boundaries.push_back(bin);

  }

  std::vector<double> bin_boundaries_dbl;
  for(size_t i=0;i<bin_boundaries.size();i++)
    bin_boundaries_dbl.push_back(h_tmp->GetBinCenter(bin_boundaries.at(i)));

  return bin_boundaries_dbl;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Generates a  2D binning scheme, uses a test file (in_file) with some mc events to determine
// a sensible scheme

vector<std::pair<double,vector<double>>> GetBinning2D(string in_file,string label,string signal_def,string sel_def,vector<string> true_variables,vector<string> reco_variables){

  if(true_variables.size() != 2 || reco_variables.size() != 2){
    std::cout << "GetBinning2D: Only supported for two variables at the moment" << std::endl;
    return {};
  }

  TFile* f_in = TFile::Open(in_file.c_str());
  TTree* t = static_cast<TTree*>(f_in->Get("stv_tree")); 

  std::string draw_query = true_variables.at(1) + ":" + true_variables.at(0);
  t->Draw(draw_query.c_str(),(signal_def+" && "+sel_def).c_str(),"col");

  // Make a histogram
  TH2F* h_tmp_2 = static_cast<TH2F*>(gPad->GetPrimitive("htemp"));

  // adjust the binning
  double x_low_lim =  h_tmp_2->GetXaxis()->GetBinLowEdge(1);
  double x_high_lim = h_tmp_2->GetXaxis()->GetBinLowEdge(h_tmp_2->GetNbinsX()+1);
  double y_low_lim =  h_tmp_2->GetYaxis()->GetBinLowEdge(1);
  double y_high_lim = h_tmp_2->GetYaxis()->GetBinLowEdge(h_tmp_2->GetNbinsY()+1);
  TH2F* h_tmp = new TH2F("h_tmp","",500,x_low_lim,x_high_lim,500,y_low_lim,y_high_lim);
  t->Draw((draw_query+">>+h_tmp").c_str(),signal_def.c_str(),"col"); 

  double tot = h_tmp->Integral();

  // for each bin in x, perform the 1D binning calc   
  TH1F* h_x_proj = reinterpret_cast<TH1F*>(h_tmp->ProjectionX()); // slightly dangerous but I think it's ok here
  vector<double> x_binning = GetBinning1D_hist(h_x_proj,0.009);

  vector<std::pair<double,vector<double>>> binning;

  // for each x bin, make the y binning
  for(size_t i_b=0;i_b<x_binning.size()-1;i_b++){
    int low_bin = h_tmp->GetXaxis()->FindBin(x_binning.at(i_b));
    int high_bin = h_tmp->GetXaxis()->FindBin(x_binning.at(i_b+1))-1;

    // slightly dangerous but I think it's ok here
    TH1F* h_y_proj = reinterpret_cast<TH1F*>(h_tmp->ProjectionY("h_y_proj",low_bin,high_bin));

    vector<double> y_binning = GetBinning1D_hist(h_y_proj,0.015);
    binning.push_back(std::make_pair(x_binning.at(i_b),y_binning));

    delete h_y_proj;

  }

  binning.push_back(std::make_pair(x_binning.back(),std::vector<double>()));

  f_in->Close();

  return binning;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Generates a 1D binning scheme, uses a test file (in_file) with some mc events to determine
// a sensible scheme

vector<double> GetBinning1D(string in_file,string label,string signal_def,string sel_def,string true_variable,string reco_variable,double stat_err = stat_err_max){

  TFile* f_in = TFile::Open(in_file.c_str());
  TTree* t = static_cast<TTree*>(f_in->Get("stv_tree")); 

  t->Draw(true_variable.c_str(),(signal_def+" && "+sel_def).c_str(),"e1");

  // Make a histogram
  TH1F* h_tmp = static_cast<TH1F*>(gPad->GetPrimitive("htemp"));
  std::vector<double> binning  = GetBinning1D_hist(h_tmp,stat_err);
  f_in->Close();
  return binning;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
