# Transfer Entropy Toolbox #

A suite of MATLAB/C and C++ tools for computing standard and extended versions of [Thomas Schreiber's transfer entropy](http://prl.aps.org/abstract/PRL/v85/i2/p461_1) on sparse, binary time series.

## What is Transfer Entropy (TE)? ##

From Schreiber, 2000:
> An information theoretic measure is derived that quantifies the statistical coherence between systems evolving in time. The standard time delayed mutual information fails to distinguish information that is actually exchanged from shared information due to common history and input signals. In our new approach, these influences are excluded by appropriate conditioning of transition probabilities. The resulting _transfer entropy_ is able to distinguish effectively driving and responding elements and to detect asymmetry in the interaction of subsystems.

## What Versions of TE are Available in the Toolbox? ##

  * **Delay 1 Transfer Entropy (D1TE)**
    * Standard TE where both the message length and sender delay are a single time bin
  * **Delayed Transfer Entropy (TE)**
    * TE with a message length of one time bin and a variable sender delay
  * **Higher-order Transfer Entropy (HOTE)**
    * TE with a variable message length and sender delay

## How Fast is It? ##

With a standard desktop computer, it took about 5 minutes to compute higher-order TE over 30 different time delays for an hour of data on all the pairs of a 200 neuron network firing at rate of 7 Hz!

In general, the following factors will affect computation time:
  * **Number of delays** - on a single-core machine, doubling the number of delays will double the computation time (linear growth). However, we use the `parfor` looping construct, so multi-core machines will receive a speed-up depending on how many cores are available
  * **Length of the data set** - assuming a fairly consistent firing rate, doubling the duration will double the computation time (linear growth)
  * **Number of neurons** - doubling the number of neurons will quadruple the computation time (quadratic growth)
  * **Message length (order)** - adding one message length will double the computation time (exponential growth)

## How Do I Get Started? ##

Head over to the [Downloads](http://code.google.com/p/transfer-entropy-toolbox/downloads/list) to get the latest version of the TE Toolbox. See the [Documentation](Documentation.md) for details on each MATLAB function and the [Examples](Examples.md) page for a few simple demos. And here is a publication that uses the toolbox for simulated neuronal data. http://dx.plos.org/10.1371/journal.pone.0027431


---


Material in these pages is based upon work supported by the [National Science Foundation](http://www.nsf.gov/div/index.jsp?div=IIS) under [Grant No. 0904912](http://www.nsf.gov/awardsearch/showAward.do?AwardNumber=0904912). Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.