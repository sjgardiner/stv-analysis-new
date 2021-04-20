#pragma once

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
