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
#include <vector>
#include <climits>

#include "transent.hpp"

int main(int argc, char **argv) {

  // 1110110100110011101101
  int raw_x[] = { 1, 2, 3, 5, 6, 8, 11, 12, 15, 16, 17, 19, 20, 22 };

  // 1010101011011011011001
  int raw_y[] = { 1, 3, 5, 7, 9, 10, 12, 13, 15, 16, 18, 19, 22 };

  // Time series need to be in proper containers
  std::vector<int> series_x(raw_x, raw_x + 14),
                   series_y(raw_y, raw_y + 13);

  // Add in terminating elements
  series_x.push_back(INT_MAX);
  series_y.push_back(INT_MAX);

  // Package all time series into a single collection
  std::vector< std::vector<int> > all_series;
  all_series.push_back(series_x);
  all_series.push_back(series_y);

  // Calculate TE
  double te_result[2][2];

  // First order (fastest)
  transent_1(all_series,
             1, /* y_delay */
             22, /* duration */
             te_result);

  std::cout << "1st order TE (delay 1)" << std::endl;
  std::cout << "y -> x: " << te_result[0][1] << std::endl;
  std::cout << "x -> y: " << te_result[1][0] << std::endl << std::endl;

  // Higher order (known at compile time -- slower)
  typedef std::vector< std::vector<int> > TimeSeriesCollection;
  typedef double ResultMatrix[2][2];

  transent_ho
    <TimeSeriesCollection, ResultMatrix, 2, 1> /* 2nd order in x */
    (all_series,
     3, /* y_delay */
     22, /* duration */
     te_result);

  std::cout << "2nd order TE in x (delay 3)" << std::endl;
  std::cout << "y -> x: " << te_result[0][1] << std::endl;
  std::cout << "x -> y: " << te_result[1][0] << std::endl << std::endl;

  // Higher order (known at run time -- slowest)
  transent_ho
    (all_series,
     1, /* x_order */
     3, /* y_order */
     2, /* y_delay */
     22, /* duration */
     te_result);

  std::cout << "3rd order TE in y (delay 2)" << std::endl;
  std::cout << "y -> x: " << te_result[0][1] << std::endl;
  std::cout << "x -> y: " << te_result[1][0] << std::endl << std::endl;

  return (0);
}

