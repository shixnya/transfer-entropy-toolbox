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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <time.h>
#include <inttypes.h>

typedef int TimeType;

void transent_1
(TimeType **all_series, const size_t series_count,
 const size_t *series_lengths,
 const TimeType y_delay,
 const TimeType duration,
 double *te_result);

void transent_ho
(TimeType **all_series, const size_t series_count,
 const size_t *series_lengths,
 const size_t x_order, const size_t y_order,
 const TimeType y_delay,
 const TimeType duration,
 double *te_result);

// ===========================================================================

int read_int(FILE *fp) {
  unsigned char buffer[4];
  fread(buffer, 1, 4, fp);
  return (((int)buffer[0] << 24) | ((int)buffer[1] << 16) | ((int)buffer[2] << 8) | (int)buffer[3]);
}

// ===========================================================================

int main(int argc, char *argv[]) {

  FILE *fp;
  size_t i, j;

  size_t series_count;
  size_t x_order, y_order;
  TimeType y_delay, duration;
  TimeType **all_series;
  double *te_result;

  if (argc < 4) {
    printf("Usage: transent series_file results_file y_delay [x_order] [y_order]\n");
    return (0);
  }

  // Read in time series
  fp = fopen(argv[1], "r");
  duration = read_int(fp), series_count = read_int(fp);

  all_series = (TimeType**)malloc(sizeof(TimeType*) * series_count);
  size_t* series_lengths = (size_t*)malloc(sizeof(size_t) * series_count);

  for (i = 0; i < series_count; ++i) {
    series_lengths[i] = read_int(fp) + 1; // +1 for terminator
    all_series[i] = (TimeType*)malloc(sizeof(TimeType) * series_lengths[i]);
  }

  for (i = 0; i < series_count; ++i) {
    for (j = 0; j < series_lengths[i] - 1; ++j) {
      all_series[i][j] = read_int(fp);
    }

    all_series[i][j] = duration + 1; // terminator
  }

  fclose(fp);

  /* Extract arguments */
  y_delay = (TimeType)atoi(argv[3]);

  if (argc > 4) {
    x_order = (size_t)atoi(argv[4]);
  }
  else {
    x_order = 1;
  }

  if (argc > 5) {
    y_order = (size_t)atoi(argv[5]);
  }
  else {
    y_order = 1;
  }

  /* Create result matrix */
	te_result = (double*)malloc(sizeof(double) * series_count * series_count);

  /* Do calculation */
  if ((x_order == 1) && (y_order == 1)) {
    transent_1(all_series, series_count,
               series_lengths,
               y_delay, duration,
               te_result);
  } else {
    transent_ho(all_series, series_count,
                series_lengths,
                x_order, y_order,
                y_delay, duration,
                te_result);
  }

  fp = fopen(argv[2], "w");

  for (i = 0; i < series_count; ++i) {
    for (j = 0; j < series_count; ++j) {
      fprintf(fp, "%e ", te_result[(i * series_count) + j]);
    }

    fprintf(fp, "\n");
  }

  fclose(fp);

  // Clean up
  for (i = 0; i < series_count; ++i) {
    free(all_series[i]);
  }

  free(all_series);
  free(series_lengths);
  free(te_result);

  return (0);
}

/* Computes the first-order transfer entropy matrix for all pairs. */
void transent_1
(TimeType **all_series, const size_t series_count,
 const size_t *series_lengths,
 const TimeType y_delay,
 const TimeType duration,
 double *te_result) {

  /* Constants */
  const size_t x_order = 1, y_order = 1,                
               num_series = 3,
               num_counts = 8,
               num_x = 4,
               num_y = 2;

  /* Locals */
  TimeType counts[num_counts];
  uint64_t code;
  size_t k, l, idx, c1, c2;
  double te_final, prob_1, prob_2, prob_3;

  TimeType *ord_iter[num_series];
  TimeType *ord_end[num_series];

  TimeType ord_times[num_series];
  TimeType ord_shift[num_series];

  const size_t window = y_order + y_delay;
  const TimeType end_time = duration - window + 1;
  TimeType cur_time, next_time;

  /* Calculate TE */
  TimeType *i_series, *j_series;
  size_t i_size, j_size;
  size_t i, j;

  for (i = 0; i < series_count; ++i) {
    for (j = 0; j < series_count; ++j) {

      /* Extract series */
      i_size = series_lengths[i];
      i_series = all_series[i];

      j_size = series_lengths[j];
      j_series = all_series[j];

      /* Order is x^(k+1), y^(l) */
      idx = 0;

      /* x^(k+1) */
      for (k = 0; k < (x_order + 1); ++k) {
        ord_iter[idx] = i_series;
        ord_end[idx] = i_series + i_size;
        ord_shift[idx] = (window - 1) - k;

        while (*(ord_iter[idx]) < ord_shift[idx] + 1) {
          ++(ord_iter[idx]);
        }

        ord_times[idx] = *(ord_iter[idx]) - ord_shift[idx];
        ++idx;
      }

      /* y^(l) */
      for (k = 0; k < y_order; ++k) {
        ord_iter[idx] = j_series;
        ord_end[idx] = j_series + j_size;
        ord_shift[idx] = -k;
        ord_times[idx] = *(ord_iter[idx]) - ord_shift[idx];
        ++idx;
      }

      /* Count spikes */
      memset(counts, 0, sizeof(TimeType) * num_counts);

      /* Get minimum next time bin */
      cur_time = ord_times[0];
      for (k = 1; k < num_series; ++k) {
        if (ord_times[k] < cur_time) {
          cur_time = ord_times[k];
        }
      }

      while (cur_time <= end_time) {

        code = 0;
        next_time = end_time + 1;

        /* Calculate hash code for this time bin */
        for (k = 0; k < num_series; ++k) {
          if (ord_times[k] == cur_time) {      
            code |= 1 << k;

            /* Next spike for this neuron */
            ++(ord_iter[k]);

            if (ord_iter[k] == ord_end[k]) {
              ord_times[k] = end_time + 1;
            }
            else {
              ord_times[k] = *(ord_iter[k]) - ord_shift[k];
            }
          }

          /* Find minimum next time bin */
          if (ord_times[k] < next_time) {
            next_time = ord_times[k];
          }
        }

        ++(counts[code]);
        cur_time = next_time;

      } /* while spikes left */

      /* Fill in zero count */
      counts[0] = end_time;
      for (k = 1; k < num_counts; ++k) {
        counts[0] -= counts[k];
      }

      /* ===================================================================== */

      /* Use counts to calculate TE */
      te_final = 0;

      /* Order is x^(k), y^(l), x(n+1) */
      for (k = 0; k < num_counts; ++k) {
        prob_1 = (double)counts[k] / (double)end_time;

        if (prob_1 == 0) {
          continue;
        }

        prob_2 = (double)counts[k] / (double)(counts[k] + counts[k ^ 1]);

        c1 = 0;
        c2 = 0;

        for (l = 0; l < num_y; ++l) {
          idx = (k & (num_x - 1)) + (l << (x_order + 1));
          c1 += counts[idx];
          c2 += (counts[idx] + counts[idx ^ 1]);
        }

        prob_3 = (double)c1 / (double)c2;

        te_final += (prob_1 * log2(prob_2 / prob_3));
      }

      te_result[(j * series_count) + i] = te_final;

    } /* for j */

  } /* for i */
 
} /* transent_1 */


/* Computes the higher-order transfer entropy matrix for all pairs. */
void transent_ho
(TimeType **all_series, const size_t series_count,
 const size_t *series_lengths,
 const size_t x_order, const size_t y_order,
 const TimeType y_delay,
 const TimeType duration,
 double *te_result) {

  /* Constants */
  const size_t num_series = 1 + y_order + x_order,
               num_counts = (size_t)pow(2, num_series),
               num_x = (size_t)pow(2, x_order + 1),
               num_y = (size_t)pow(2, y_order);

  /* Locals */
  TimeType *counts = (TimeType*)malloc(sizeof(TimeType) * num_counts);
  uint64_t code;
  size_t k, l, idx, c1, c2;
  double te_final, prob_1, prob_2, prob_3;

  TimeType *ord_iter[num_series];
  TimeType *ord_end[num_series];

  TimeType ord_times[num_series];
  TimeType ord_shift[num_series];

  const size_t window = (y_order + y_delay) > (x_order + 1) ? (y_order + y_delay) : (x_order + 1);
  const TimeType end_time = duration - window + 1;
  TimeType cur_time, next_time;

  /* Calculate TE */
  TimeType *i_series, *j_series;
  size_t i_size, j_size;
  size_t i, j;

  for (i = 0; i < series_count; ++i) {
    for (j = 0; j < series_count; ++j) {

      /* Extract series */
      i_size = series_lengths[i];
      i_series = all_series[i];

      j_size = series_lengths[j];
      j_series = all_series[j];

      /* Order is x^(k+1), y^(l) */
      idx = 0;

      /* x^(k+1) */
      for (k = 0; k < (x_order + 1); ++k) {
        ord_iter[idx] = i_series;
        ord_end[idx] = i_series + i_size;
        ord_shift[idx] = (window - 1) - k;

        while (*(ord_iter[idx]) < ord_shift[idx] + 1) {
          ++(ord_iter[idx]);
        }

        ord_times[idx] = *(ord_iter[idx]) - ord_shift[idx];
        ++idx;
      }

      /* y^(l) */
      for (k = 0; k < y_order; ++k) {
        ord_iter[idx] = j_series;
        ord_end[idx] = j_series + j_size;
        ord_shift[idx] = -k;
        ord_times[idx] = *(ord_iter[idx]) - ord_shift[idx];
        ++idx;
      }

      /* Count spikes */
      memset(counts, 0, sizeof(TimeType) * num_counts);

      /* Get minimum next time bin */
      cur_time = ord_times[0];
      for (k = 1; k < num_series; ++k) {
        if (ord_times[k] < cur_time) {
          cur_time = ord_times[k];
        }
      }

      while (cur_time <= end_time) {

        code = 0;
        next_time = end_time + 1;

        /* Calculate hash code for this time bin */
        for (k = 0; k < num_series; ++k) {
          if (ord_times[k] == cur_time) {      
            code |= 1 << k;

            /* Next spike for this neuron */
            ++(ord_iter[k]);

            if (ord_iter[k] == ord_end[k]) {
              ord_times[k] = end_time + 1;
            }
            else {
              ord_times[k] = *(ord_iter[k]) - ord_shift[k];
            }
          }

          /* Find minimum next time bin */
          if (ord_times[k] < next_time) {
            next_time = ord_times[k];
          }
        }

        ++(counts[code]);
        cur_time = next_time;

      } /* while spikes left */

      /* Fill in zero count */
      counts[0] = end_time;
      for (k = 1; k < num_counts; ++k) {
        counts[0] -= counts[k];
      }

      /* ===================================================================== */

      /* Use counts to calculate TE */
      te_final = 0;

      /* Order is x^(k), y^(l), x(n+1) */
      for (k = 0; k < num_counts; ++k) {
        prob_1 = (double)counts[k] / (double)end_time;

        if (prob_1 == 0) {
          continue;
        }

        prob_2 = (double)counts[k] / (double)(counts[k] + counts[k ^ 1]);

        c1 = 0;
        c2 = 0;

        for (l = 0; l < num_y; ++l) {
          idx = (k & (num_x - 1)) + (l << (x_order + 1));
          c1 += counts[idx];
          c2 += (counts[idx] + counts[idx ^ 1]);
        }

        prob_3 = (double)c1 / (double)c2;

        te_final += (prob_1 * log2(prob_2 / prob_3));
      }

      te_result[(j * series_count) + i] = te_final;

    } /* for j */

  } /* for i */

  /* Clean up */
  free(counts);
 
} /* transent_ho */

