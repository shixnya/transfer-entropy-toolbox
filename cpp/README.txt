GENERAL INFORMATION
===================
C++ template library for calculating transfer entropy [1] from sparse time
series data.

A vector of N time series is used to calculate an NxN matrix M of real numbers
where M[x][y] is the transfer entropy from time series y to time series x (y ->
x). x is referred to as the "predicted" time series and y is referred to as the
"predictor" time series.

NOTE 1: All time series are assumed to be 1-based. That is, the first possible
time for a "spike" is 1.

[1] Thomas Schreiber. Measuring information transfer. Phys. Rev. Lett.,
    85(2):461–464, Jul 2000.

REQUIREMENTS
============
BOOST 1.42 or higher (may work with previous versions, but not tested).

LIBRARY USAGE
=============
All library functions are available in transent.hpp. Be careful in choosing
which "transent" function to use, as performance between them is significantly
different.

It's recommended to turn on all compiler optimizations for performance (-O3).

The row_start, rows, col_start, and cols parameters for each function can be
used to only compute a portion of the transfer entropy matrix. This is useful if
your calculation will take a long time and you'd like to run several jobs in
parallel.

First Order
-----------

template <typename TimeSeriesCollection, typename ResultMatrix>
void transent_1
(const TimeSeriesCollection& all_series,
 typename TimeSeriesCollection::value_type::value_type y_delay,
 typename TimeSeriesCollection::value_type::value_type duration,
 ResultMatrix& te_result,
 std::size_t row_start = 0, std::size_t rows = 0,
 std::size_t col_start = 0, std::size_t cols = 0)

First order transfer entropy only. If you just need first order, this will be
much faster than the other functions.

[Template Parameters]

TimeSeriesCollection - Vector of all time series containers. Must be indexable
                       and have a size() method (when rows = 0 or cols = 0).

ResultMatrix - Two-dimensional matrix for storing transfer entropy between all
               pairs of time series. Element type must be compatible with double.

[Function Parameters]

y_delay - Number of time bins to shift the predictor time series backwards. Must
          be greater than zero.

duration - Last time bin to take into account.

te_result - 2-D result matrix. Must be at least
            (rows - row_start)x(cols - col_start) in size.

row_start - Offset into all_series for predicted time series (default 0).

rows - Number of predicted time series (default 0 means all).

col_start - Offset into all_series for predictor time series (default 0).

cols - Number of predictor time series (default 0 means all).


Higher Order (compile time)
---------------------------

template <typename TimeSeriesCollection, typename ResultMatrix,
         std::size_t x_order, std::size_t y_order>
void transent_ho
(const TimeSeriesCollection& all_series,
 typename TimeSeriesCollection::value_type::value_type y_delay,
 typename TimeSeriesCollection::value_type::value_type duration,
 ResultMatrix& te_result,
 std::size_t row_start = 0, std::size_t rows = 0,
 std::size_t col_start = 0, std::size_t cols = 0)

Higher order transfer entropy where the orders are known at compile time. If you
only need first order, transent_1 will be MUCH faster.

NOTE: Combined order (x_order + y_order + 1) cannot exceed 32.

[Template Parameters]

TimeSeriesCollection - Vector of all time series containers. Must be indexable
                       and have a size() method (when rows = 0 or cols = 0).

ResultMatrix - Two-dimensional matrix for storing transfer entropy between all
               pairs of time series. Element type must be compatible with double.

x_order - Order of the predicted time series.

y_order - Order of the predictor time series.

[Function Parameters]

y_delay - Number of time bins to shift the predictor time series backwards. Must
          be greater than zero.

duration - Last time bin to take into account.

te_result - 2-D result matrix. Must be at least
            (rows - row_start)x(cols - col_start) in size.

row_start - Offset into all_series for predicted time series (default 0).

rows - Number of predicted time series (default 0 means all).

col_start - Offset into all_series for predictor time series (default 0).

cols - Number of predictor time series (default 0 means all).


Higher Order (run time)
---------------------------

template <typename TimeSeriesCollection, typename ResultMatrix>
void transent_ho
(const TimeSeriesCollection& all_series,
 typename std::size_t x_order, std::size_t y_order,
 typename TimeSeriesCollection::value_type::value_type y_delay,
 typename TimeSeriesCollection::value_type::value_type duration,
 ResultMatrix& te_result,
 std::size_t row_start = 0, std::size_t rows = 0,

Higher order transfer entropy where the orders not known until run time. If you
only need first order, transent_1 will be MUCH faster. If your orders are known
at compile time, the other transent_ho (above) will be faster.

NOTE: Combined order (x_order + y_order + 1) cannot exceed 32.

[Template Parameters]

TimeSeriesCollection - Vector of all time series containers. Must be indexable
                       and have a size() method (when rows = 0 or cols = 0).

ResultMatrix - Two-dimensional matrix for storing transfer entropy between all
               pairs of time series. Element type must be compatible with double.

[Function Parameters]

x_order - Order of the predicted time series.

y_order - Order of the predictor time series.

y_delay - Number of time bins to shift the predictor time series backwards. Must
          be greater than zero.

duration - Last time bin to take into account.

te_result - 2-D result matrix. Must be at least
            (rows - row_start)x(cols - col_start) in size.

row_start - Offset into all_series for predicted time series (default 0).

rows - Number of predicted time series (default 0 means all).

col_start - Offset into all_series for predictor time series (default 0).

cols - Number of predictor time series (default 0 means all).


PROGRAM USAGE
=============
There are three programs included in te_block*.cpp. After compiling them, run
with --help for an explanation of the arguments.

A time series file is ASCII and has the following format: a duration on the
first line and each time series on successive lines. An individual time series
is a sparse list of active time bins in ascending order, separated by spaces.

If we had two time series that went for 10 time bins, one active on the even
bins and one active on the odd bins, the file would look like this:
10
2 4 6 8 10
1 3 5 7 9

Output is written as an ASCII file where each row is a row in the transfer
entropy matrix and columns are separated by spaces. For the previous example,
the output file would be of the form:
0 a
b 0

where "a" is the transfer entropy from the first series to the second and "b" is
the reverse. In this example, both "a" and "b" will be zero.

All programs take row-start, rows, col-start, and cols arguments. These are used
to calculate only a portion of the transfer entropy matrix. This is useful if
you want to split up a long calculation into several parallel jobs.

te_block_1 - Calculates first order transfer entropy for a block of time series.

te_block_fixed - Calculates higher order (fixed at compile time) transfer
                 entropy for a block of time series/

te_block - Calculates higher order transfer entropy for a block of time series.

EXAMPLE
=======
See example.cpp

AUTHORS
=======
Michael Hansen <mihansen@indiana.edu> (Indiana University)
Shinya Ito (Indiana University)

