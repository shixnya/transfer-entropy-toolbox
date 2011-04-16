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

#include <bitset>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cassert>

#include <boost/mpl/arithmetic.hpp>
#include <boost/static_assert.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/limits.hpp>

#define MAX_XY_ORDER 64

namespace mpl = boost::mpl;
namespace boost { namespace mpl {
  template <std::size_t N, std::size_t Power>
  struct pow
  {
    static const std::size_t value = N * pow<N, Power - 1>::value;
  };

  template <std::size_t N>
  struct pow<N, 0>
  {
    static const std::size_t value = 1;
  };

} }

// Computes the higher-order transfer entropy matrix for all pairs.
// x and y orders must be known at compile time.
template <typename TimeSeriesCollection, typename ResultMatrix,
         std::size_t x_order, std::size_t y_order>
void transent_ho
(const TimeSeriesCollection& all_series,
 const typename TimeSeriesCollection::value_type::value_type y_delay,
 const typename TimeSeriesCollection::value_type::value_type duration,
 ResultMatrix& te_result,
 std::size_t row_start = 0, std::size_t rows = 0,
 std::size_t col_start = 0, std::size_t cols = 0) {

  // Typedefs
  typedef typename TimeSeriesCollection::value_type TimeSeries;
  typedef typename TimeSeries::value_type TimeType;
  typedef typename TimeSeries::const_iterator TimeSeriesIter;

  typedef std::pair<TimeSeriesIter, TimeType> IterShiftPair;
  typedef std::pair<TimeType, std::size_t> TimeIndexPair;

  // Constants
  const std::size_t num_series = 1 + y_order + x_order,
                    num_counts = mpl::pow<2, num_series>::value,
                    num_x = mpl::pow<2, mpl::plus<mpl::int_<x_order>, mpl::int_<1> >::value>::value,
                    num_y = mpl::pow<2, y_order>::value;

  BOOST_STATIC_ASSERT(x_order > 0);
  BOOST_STATIC_ASSERT(y_order > 0);
  assert(y_delay > 0);
  BOOST_STATIC_ASSERT(num_series <= MAX_XY_ORDER);

  if (rows == 0) {
    rows = all_series.size();
  }

  if (cols == 0) {
    cols = all_series.size();
  }

  // Locals
  std::vector<TimeType> counts(num_counts);
  std::bitset<MAX_XY_ORDER> code;
  std::size_t idx = 0;
  double te_final, prob_2, prob_3;

  IterShiftPair ord_iter[num_series];
  TimeType ord_times[num_series];
  TimeSeriesIter ord_end[num_series];

  const std::size_t window = std::max(y_order + y_delay, x_order + 1);
  TimeType cur_time, next_time, end_time = duration - window + 1, shift;

  // Calculate TE
  for (std::size_t i = row_start; i < (rows + row_start); ++i) {
    for (std::size_t j = col_start; j < (cols + col_start); ++j) {

      // NOTE: Time series are assumed to be 1-based, so everything is shifted by 1 too.
      // Order is x^(k+1), y^(l)
      idx = 0;

      // x^(k+1)
      for (std::size_t k = 0; k < (x_order + 1); ++k) {
        shift = (window - 1) - k;
        ord_end[idx] = all_series[i].end();
        ord_iter[idx] = std::make_pair(std::lower_bound(all_series[i].begin(), ord_end[idx], shift + 1), shift);
        ord_times[idx] = *(ord_iter[idx].first) - ord_iter[idx].second;
        ++idx;
      }

      // y^(l)
      for (std::size_t k = 0; k < y_order; ++k) {
        ord_end[idx] = all_series[j].end();
        ord_iter[idx] = std::make_pair(all_series[j].begin(), -k);
        ord_times[idx] = *(ord_iter[idx].first) - ord_iter[idx].second;
        ++idx;
      }

      // Count spikes
      std::fill(counts.begin(), counts.end(), 0);
      cur_time = *(std::min_element(ord_times, ord_times + num_series));

      while (cur_time <= end_time) {

        code.reset();
        next_time = std::numeric_limits<TimeType>::max();

        // Calculate hash code for this time
        for (std::size_t k = 0; k < num_series; ++k) {
          if (ord_times[k] == cur_time) {        
            code[k] = 1;

            // Next spike
            ++(ord_iter[k].first);

            if (ord_iter[k].first == ord_end[k]) {
              ord_times[k] = std::numeric_limits<TimeType>::max();
            }
            else {
              ord_times[k] = *(ord_iter[k].first) - ord_iter[k].second;
            }
          }

          if (ord_times[k] < next_time) {
            next_time = ord_times[k];
          }
        }

        ++(counts[code.to_ulong()]);
        cur_time = next_time;

      } // while spikes left

      // Fill in zero count
      counts[0] = end_time - std::accumulate(counts.begin() + 1, counts.end(), 0);

      // =====================================================================

      // Use counts to calculate TE
      te_final = 0;

      // Order is x^(k), y^(l), x(n+1)
      for (std::size_t k = 0; k < num_counts; ++k) {
        if (counts[k] == 0) {
          continue;
        }

        prob_2 = (double)counts[k] / (double)(counts[k] + counts[k ^ 1]);

        std::size_t c1 = 0, c2 = 0;
        for (std::size_t l = 0; l < num_y; ++l) {
          idx = (k & (num_x - 1)) + (l << (x_order + 1));
          c1 += counts[idx];
          c2 += (counts[idx] + counts[idx ^ 1]);
        }

        prob_3 = (double)c1 / (double)c2;

        te_final += ((double)counts[k] * (log2(prob_2) - log2(prob_3)));
      }

      te_result[i - row_start][j - col_start] = te_final / (double)end_time;

    } // for j

  } // for i

} //transent_ho

// Computes the 1st order transfer entropy matrix for all pairs.
template <typename TimeSeriesCollection, typename ResultMatrix>
void transent_1
(const TimeSeriesCollection& all_series,
 const typename TimeSeriesCollection::value_type::value_type y_delay,
 const typename TimeSeriesCollection::value_type::value_type duration,
 ResultMatrix& te_result,
 std::size_t row_start = 0, std::size_t rows = 0,
 std::size_t col_start = 0, std::size_t cols = 0) {

  return (transent_ho<TimeSeriesCollection, ResultMatrix, 1, 1>
          (all_series, y_delay, duration, te_result,
           row_start, rows, col_start, cols));

} // transent_1


// Computes the higher-order transfer entropy matrix for all pairs.
template <typename TimeSeriesCollection, typename ResultMatrix>
void transent_ho
(const TimeSeriesCollection& all_series,
 const typename std::size_t x_order, std::size_t y_order,
 const typename TimeSeriesCollection::value_type::value_type y_delay,
 const typename TimeSeriesCollection::value_type::value_type duration,
 ResultMatrix& te_result,
 std::size_t row_start = 0, std::size_t rows = 0,
 std::size_t col_start = 0, std::size_t cols = 0) {

  // Typedefs
  typedef typename TimeSeriesCollection::value_type TimeSeries;
  typedef typename TimeSeries::value_type TimeType;
  typedef typename TimeSeries::const_iterator TimeSeriesIter;

  typedef std::pair<TimeSeriesIter, TimeType> IterShiftPair;
  typedef std::pair<TimeType, std::size_t> TimeIndexPair;

  // Constants
  const std::size_t num_series = 1 + y_order + x_order,
                    num_counts = (std::size_t)pow(2, num_series),
                    num_x = (std::size_t)pow(2, x_order + 1),
                    num_y = (std::size_t)pow(2, y_order);

  assert(x_order > 0);
  assert(y_order > 0);
  assert(y_delay > 0);
  assert(num_series <= MAX_XY_ORDER);

  if (rows == 0) {
    rows = all_series.size();
  }

  if (cols == 0) {
    cols = all_series.size();
  }

  // Locals
  std::vector<TimeType> counts(num_counts);
  std::bitset<MAX_XY_ORDER> code;
  std::size_t idx = 0;
  double te_final, prob_2, prob_3;

  IterShiftPair ord_iter[num_series];
  TimeType ord_times[num_series];
  TimeSeriesIter ord_end[num_series];

  const std::size_t window = std::max(y_order + y_delay, x_order + 1);
  TimeType cur_time, next_time, end_time = duration - window + 1, shift;

  // Calculate TE
  for (std::size_t i = row_start; i < (rows + row_start); ++i) {
    for (std::size_t j = col_start; j < (cols + col_start); ++j) {

      // NOTE: Time series are assumed to be 1-based, so everything is shifted by 1 too.
      // Order is x^(k+1), y^(l)
      idx = 0;

      // x^(k+1)
      for (std::size_t k = 0; k < (x_order + 1); ++k) {
        shift = (window - 1) - k;
        ord_end[idx] = all_series[i].end();
        ord_iter[idx] = std::make_pair(std::lower_bound(all_series[i].begin(), ord_end[idx], shift + 1), shift);
        ord_times[idx] = *(ord_iter[idx].first) - ord_iter[idx].second;
        ++idx;
      }

      // y^(l)
      for (std::size_t k = 0; k < y_order; ++k) {
        ord_end[idx] = all_series[j].end();
        ord_iter[idx] = std::make_pair(all_series[j].begin(), -k);
        ord_times[idx] = *(ord_iter[idx].first) - ord_iter[idx].second;
        ++idx;
      }

      // Count spikes
      std::fill(counts.begin(), counts.end(), 0);
      cur_time = *(std::min_element(ord_times, ord_times + num_series));

      while (cur_time <= end_time) {

        code.reset();
        next_time = std::numeric_limits<TimeType>::max();

        // Calculate hash code for this time
        for (std::size_t k = 0; k < num_series; ++k) {
          if (ord_times[k] == cur_time) {        
            code[k] = 1;

            // Next spike
            ++(ord_iter[k].first);

            if (ord_iter[k].first == ord_end[k]) {
              ord_times[k] = std::numeric_limits<TimeType>::max();
            }
            else {
              ord_times[k] = *(ord_iter[k].first) - ord_iter[k].second;
            }
          }

          if (ord_times[k] < next_time) {
            next_time = ord_times[k];
          }
        }

        ++(counts[code.to_ulong()]);
        cur_time = next_time;

      } // while spikes left

      // Fill in zero count
      counts[0] = end_time - std::accumulate(counts.begin() + 1, counts.end(), 0);

      // =====================================================================

      // Use counts to calculate TE
      te_final = 0;

      // Order is x^(k), y^(l), x(n+1)
      for (std::size_t k = 0; k < num_counts; ++k) {
        if (counts[k] == 0) {
          continue;
        }

        prob_2 = (double)counts[k] / (double)(counts[k] + counts[k ^ 1]);

        std::size_t c1 = 0, c2 = 0;
        for (std::size_t l = 0; l < num_y; ++l) {
          idx = (k & (num_x - 1)) + (l << (x_order + 1));
          c1 += counts[idx];
          c2 += (counts[idx] + counts[idx ^ 1]);
        }

        prob_3 = (double)c1 / (double)c2;

        te_final += ((double)counts[k] * (log2(prob_2) - log2(prob_3)));
      }

      te_result[i - row_start][j - col_start] = te_final / (double)end_time;

    } // for j

  } // for i

} // transent_ho

