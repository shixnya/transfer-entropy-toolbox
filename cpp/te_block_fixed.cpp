/*=============================================================================
Copyright (c) 2011, The Trustees of Indiana University
All rights reserved.

Authors: Michael Hansen (mihansen@indiana.edu), Shinya Ito

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer in the documentation
     and/or other materials provided with the distribution.

  3. Neither the name of Indiana University nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=============================================================================*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>

#include <boost/lexical_cast.hpp>
#include <boost/limits.hpp>
#include <boost/multi_array.hpp>
#include <boost/program_options.hpp>

#include "transent.hpp"

#ifndef X_ORDER
  #define X_ORDER 1
#endif

#ifndef Y_ORDER
  #define Y_ORDER 1
#endif

// Typedefs
typedef int TimeType;
typedef std::vector<TimeType> TimeSeries;
typedef std::vector< std::vector<TimeType> > TimeSeriesCollection;

typedef boost::multi_array<double, 2> ResultMatrix;
typedef ResultMatrix::index arr_index;

int main(int argc, char *argv[]) {

  namespace opt = boost::program_options;
  std::stringstream desc_stream;
  desc_stream << "Calculates transfer entropy for a block of time series (y -> x)" <<
    " with x-order = " << X_ORDER << ", y-order = " << Y_ORDER;

  opt::options_description desc(desc_stream.str());
  desc.add_options()
    ("help", "Show this help message")
    ("y-delay", opt::value<TimeType>()->default_value(1), "Delay of predictor time series (default 1)")
    ("in-file", opt::value<std::string>(), "Input time series file path")
    ("out-file", opt::value<std::string>(), "Output transfer entropy file path")
    ("col-start", opt::value<arr_index>()->default_value(0), "Column offset of block (default 0)")
    ("cols", opt::value<arr_index>()->default_value(0), "Columns in block (default 0 for remainder)")
    ("row-start", opt::value<arr_index>()->default_value(0), "Row offset of block (default 0)")
    ("rows", opt::value<arr_index>()->default_value(0), "Rows in block (default 0 for remainder)")
    ;

  opt::variables_map opt_vars;
  opt::store(opt::parse_command_line(argc, argv, desc), opt_vars);
  opt::notify(opt_vars);

  if (opt_vars.count("help")) {
    std::cout << desc << std::endl;
    return (0);
  }

  const std::size_t x_order = X_ORDER,
                    y_order = Y_ORDER,
                    y_delay = opt_vars["y-delay"].as<TimeType>();

  const std::size_t num_series = 1 + y_order + x_order;

  assert(x_order > 0);
  assert(y_order > 0);
  assert(y_delay > 0);

  if (num_series > MAX_XY_ORDER) {
    std::cout << "The combined order of x and y cannot exceed " << MAX_XY_ORDER << std::endl;
    return (0);
  }

  // Parse arguments
  if (!opt_vars.count("in-file") || !opt_vars.count("out-file")) {
    std::cout << "Input and output file paths are required" << std::endl;
    return (0);
  }

  std::string in_file_path = opt_vars["in-file"].as<std::string>(),
              out_file_path = opt_vars["out-file"].as<std::string>();

  arr_index col_start = opt_vars["col-start"].as<arr_index>(),
            cols = opt_vars["cols"].as<arr_index>(),
            row_start = opt_vars["row-start"].as<arr_index>(),
            rows = opt_vars["rows"].as<arr_index>();

  // Read in time series block
  TimeSeriesCollection all_series;

  std::ifstream in_file(in_file_path.c_str());
  std::string line;

  getline(in_file, line);
  TimeType duration = boost::lexical_cast<TimeType>(line);

  while (getline(in_file, line)) {

    std::istringstream line_stream(line);
    TimeSeries cur_series;

    // This could be more efficient, but it's fast enough for now
    std::copy(std::istream_iterator<TimeType>(line_stream),
              std::istream_iterator<TimeType>(),
              std::back_inserter(cur_series));

    all_series.push_back(cur_series);
  }

  if (rows == 0) {
    rows = all_series.size();
  }

  if (cols == 0) {
    cols = all_series.size();
  }

  // Calculate TE
  ResultMatrix te_result(boost::extents[rows][cols]);

  transent_ho<TimeSeriesCollection, ResultMatrix, x_order, y_order>
    (all_series, y_delay, duration, te_result,
     row_start, rows, col_start, cols);

  // Write results
  std::ofstream out_file(out_file_path.c_str());

  for (arr_index i = 0; i < rows; ++i) {
    for (arr_index j = 0; j < cols; ++j) {
      out_file << te_result[j][i] << " ";
    }

    out_file << std::endl;
  }

  return (0);
}

