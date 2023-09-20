#ifndef _Misc_hh_
#define _Misc_hh_

#include <sstream>

using std::vector;
using std::pair;
using std::map;

template <typename T> std::string to_string_with_precision(const T a_value, const int n = 6)
{
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return std::move(out).str();
}

#endif
