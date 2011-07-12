GENERAL INFORMATION
===================
MATLAB/C library for calculating first-order and higher-order delayed transfer
entropy [1] from sparse, binary time series data.

[1] Thomas Schreiber. Measuring information transfer. Phys. Rev. Lett.,
    85(2):461–464, Jul 2000.

REQUIREMENTS
============
MATLAB or Octave and the ability to compile MEX files. It's recommended that you
turn on all available compiler optimizations, as they can speed up calculations
significantly.

LIBRARY USAGE
=============
For calculations, all time series data is assumed to be in Another Spike Data
Format (ASDF). This is simply a MATLAB cell array with some extra metadata (see
The ASDF Format section for details).

To produce an ASDF object from a matrix of binary spikes (neurons are rows,
columns are time bins), use the SparseToASDF function:

spikes = [
           1 1 1 0 1 1 0 1 0 0 1 1 0 0 1 1 1 0 1 1 0 1; % neuron 1
           1 0 1 0 1 0 1 0 1 1 0 1 1 0 1 1 0 1 1 0 0 1  % neuron 2
         ];

asdf = SparseToASDF(spikes, 1); % Use a bin size of 1

You can now use the ASDF to calculate transfer entropy by calling ASDFTE. This
function has the following parameters:

  * asdf - ASDF object with N time series
  * j_delay - Number of bins to lag sender (j series) or vector [default 1]
  * i_order - Order of receiver [default 1]
  * j_order - Order of sender [default 1]

ASDFTE returns an NxN matrix T where T(i, j) is the transfer entropy from j->i.

By combining parameters, you can calculate the following kinds of transfer
entropy:

  * Delay 1 first-order - ASDFTE(asdf)
  * Delay d first-order - ASDFTE(asdf, d)
  * Delay d higher-order - ASDFTE(asdf, d, i_order, j_order)
  * Multi-delay (d1-d2) first-order - ASDFTE(asdf, d1:d2)
    * Takes max value, see "Choosing a TE value" below
  * Multi-delay (d1-d2) higher-order - ASDFTE(asdf, d1:d2, i_order, j_order)
    * Takes max value, see "Choosing a TE value" below

Choosing a TE value
-------------------
When computing transfer entropy for real-world data, one does not usually know
in advance the "correct" delays for every pair of connected neurons. It's
necessary, then, to calculate transfer entropy for multiple delay times and use
all of resulting TE matrices to compute the final TE matrix.

By default, ASDFTE will use the max function to choose a final TE value when
multiple delay times are given. The final TE matrix will therefore contain the
maximum TE values for all neuron pairs and all delay times.

ASDFTE(asdf, d1:d2, i_order, j_order) is equivalent to calling
ASDFTE(asdf, d1:d2, i_order, j_order, @max)

Another function, CIReduce, is provided to compute the coincidence index
instead. This is the sum of TE values in a window centered on the peak TE value
divided by the sum of all TE values for a given neuron pair. When using
CIReduce, you must provide a window size as the final argument to ASDFTE:

ASDFTE(asdf, d1:d2, i_order, j_order, @CIReduce, window)

Be careful, as the resulting matrix will now contain coincidence indices for
every neuron pair instead of TE values.

Computation Time
----------------
ASDFTE's calculation time grows differently for its different parameters. In
general, the calculation time grows:

  * Quadradically in the number of neurons (N)
  * Sub-linearly in the lengths of all time series
    * This is heavily dependent on the overlap between series
  * Linearly in the number of delay bins (d1:d2)
  * Exponentially in the total order (i_order + j_order)

Beware especially of the last bullet. Increasing the total order by 1 roughly
doubles the calculation time, so it is easy to run into calculation times
exceeding one's lifespan when going above a total order of 20 or so!

The ASDF Format
---------------
An ASDF object is a MATLAB cell array with two extra cells at the end. You can
create one manually by following the steps in this section. First, you must
create an empty cell array:

asdf = cell();

Each time series inside an ASDF object is an ordered vector of spike times (time
bins where the binary value was 1). If our first of N time series (neuron 1) had
spikes at time bins 1, 3, 6, and 10, we would write:

asdf{1} = [1, 3, 6, 10];

Additional time series are added in the same manner, one after the other:

asdf{2} = [3, 7, 8, 12, 20];
asdf{3} = [5, 11, 24];
...
asdf{N} = [2, 3, 4, 9, 25, 28];

After all of the time series have been added, we must add a cell with the bin
unit. In most cases, this will be 1:

asdf{N + 1} = 1;

Finally, we add the number of neurons and the total duration (number of recorded
time bins):

asdf{N + 2} = [N, 35];

We now have a fully functioning ASDF object!

EXAMPLE
=======
See example.m

AUTHORS
=======
Michael Hansen <mihansen@indiana.edu> (Indiana University)
Shinya Ito <itos@indiana.edu> (Indiana University)

