#pragma once

// Standard library includes
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// ROOT includes
#include "TGraph.h"
#include "TH1D.h"

// Helper function that produces the standard MicroBooNE plot legend title
// with the BNB POT displayed
std::string get_legend_title( double bnb_pot ) {
  // Print the data POT exposure to a precision of 3 decimal digits, then
  // separate out the base-ten exponent for the legend header
  std::stringstream temp_ss;
  temp_ss << std::scientific << std::setprecision(3) << bnb_pot;

  std::string pot_str = temp_ss.str();
  size_t e_pos = pot_str.find( 'e' );

  // Digits not including 'e' and the base-ten exponent
  std::string pot_digits_str = pot_str.substr( 0, e_pos );
  // pot_str now contains just the base-ten exponent
  pot_str.erase( 0, e_pos + 1u );
  // If there's a leading '+' in the exponent, erase that too
  if ( pot_str.front() == '+' ) pot_str.erase( 0, 1u );

  std::string legend_title = "MicroBooNE " + pot_digits_str + " #times 10^{"
    + pot_str + "} POT, INTERNAL";

  return legend_title;
}

// Helper function that dumps 1D histogram contents to a map of pgfplotstable
// columns
void dump_1d_histogram( const std::string& hist_col_prefix,
  const TH1D& hist,
  std::map< std::string, std::vector<std::string> >& pgf_plots_hist_table,
  bool include_yerror = true, bool include_x_coords = false )
{
  // TODO: consider adding a check for pre-existing duplicate columns
  // (to avoid accidentally overwriting information)

  const std::string bin_col_name = "bin";
  if ( include_x_coords ) {
    pgf_plots_hist_table[ bin_col_name ] = std::vector< std::string >();
  }

  std::string y_col_name = hist_col_prefix;
  pgf_plots_hist_table[ y_col_name ] = std::vector< std::string >();

  std::string yerror_col_name;
  if ( include_yerror ) {
    yerror_col_name = hist_col_prefix + "_error";
    pgf_plots_hist_table[ yerror_col_name ] = std::vector< std::string >();
  }

  // Include the ordinary bins and also the overflow bin. The latter
  // (presumably empty) will help certain pgfplots plotting styles look good
  int overflow_bin_idx = hist.GetNbinsX() + 1;
  for ( int b = 1; b <= overflow_bin_idx; ++b ) {
    double y = hist.GetBinContent( b );

    pgf_plots_hist_table.at( y_col_name )
      .push_back( std::to_string(y) );

    if ( include_yerror ) {
      double yerror = hist.GetBinError( b );
      pgf_plots_hist_table.at( yerror_col_name )
        .push_back( std::to_string(yerror) );
    }

    if ( include_x_coords ) {
      pgf_plots_hist_table.at( bin_col_name )
        .push_back( std::to_string(b) );

      double low_edge = hist.GetBinLowEdge( b );
      double half_width = 0.5 * hist.GetBinWidth( b );

      std::string x_col_name = hist_col_prefix;
      if ( !hist_col_prefix.empty() ) x_col_name += '_';
      x_col_name += 'x';
      auto end = pgf_plots_hist_table.end();
      auto iter = pgf_plots_hist_table.find( x_col_name );
      if ( iter == end ) {
        pgf_plots_hist_table[ x_col_name ]
          = std::vector< std::string > { std::to_string(low_edge) };
      }
      else iter->second.push_back( std::to_string(low_edge) );

      std::string x_hw_col_name = x_col_name + "_halfwidth";
      auto end2 = pgf_plots_hist_table.end();
      auto iter2 = pgf_plots_hist_table.find( x_hw_col_name );
      if ( iter2 == end2 ) {
        pgf_plots_hist_table[ x_hw_col_name ]
          = std::vector<std::string> { std::to_string(half_width) };
      }
      else iter2->second.push_back( std::to_string(half_width) );
    } // include x coordinates
  } // bin number
}

// Helper function that dumps TGraph contents to a map of pgfplotstable columns
// TODO: add support for TGraphErrors
void dump_tgraph( const std::string& col_prefix, const TGraph& graph,
  std::map< std::string, std::vector<std::string> >& pgf_plots_table )
{
  // TODO: consider adding a check for pre-existing duplicate columns
  // (to avoid accidentally overwriting information)

  std::string x_col_name = col_prefix + "_x";
  pgf_plots_table[ x_col_name ] = std::vector< std::string >();

  std::string y_col_name = col_prefix + "_y";
  pgf_plots_table[ y_col_name ] = std::vector< std::string >();

  int num_points = graph.GetN();
  for ( int p = 0; p < num_points; ++p ) {
    double x = graph.GetX()[ p ];
    double y = graph.GetY()[ p ];

    pgf_plots_table.at( x_col_name ).push_back( std::to_string(x) );
    pgf_plots_table.at( y_col_name ).push_back( std::to_string(y) );
  } // point index
}

void write_pgfplots_file( const std::string& out_filename,
  std::map< std::string, std::vector<std::string> >& pgfplots_table )
{
  std::ofstream out_file( out_filename );

  // Column headings
  for ( const auto& pair : pgfplots_table ) {
    const std::string& col_name = pair.first;
    out_file << "  " << col_name;
  }
  out_file << '\n';

  // Data rows (all columns are assumed to have the same number of rows)
  size_t num_rows = pgfplots_table.cbegin()->second.size();
  for ( size_t r = 0u; r < num_rows; ++r ) {
    for ( const auto& pair : pgfplots_table ) {
      out_file << "  " << pair.second.at( r );
    }
    out_file << '\n';
  }
}
