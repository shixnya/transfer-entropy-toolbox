# Getting Started #
Before using the TE Toolbox, you must first compile the mex file (C program) with `gcc` or `lcc`. You can do this within MATLAB by calling
```
mex transent.c
```
If this is your first time compiling a mex file, you may be prompted to choose a compiler.

Next, make sure all of the code files are in your MATLAB path (or the current directory) and you're ready to go! See the [Examples](Examples.md) page for a few simple demos.

# Another Spike Data Format (ASDF) #
In the MATLAB library, we represent spike trains using Another Spike Data Format. This is essentially a cell array where each cell, indexed by neuron number, is an ordered array of time bins in which the neuron was spiking. It is recommended that you create an ASDF object using the `SparseToASDF` function as follows:
```
spikes = [
  1 1 1 0 1 1 0 1 0 0 1 1 0 0 1 1 1 0 1 1 0 1; % neuron 1
  1 0 1 0 1 0 1 0 1 1 0 1 1 0 1 1 0 1 1 0 0 1  % neuron 2
];

asdf = SparseToASDF(spikes, 1); % Use a bin size of 1
```
The last two cells of an ASDF object contain special information:
```
asdf{end-1} % Binning size in miliseconds
asdf{end} % 2-cell array with [# of neurons, total # of time bins]
```
In order to calculate TE correctly, we recommend to use only integers for spike times. To ensure that your data conforms, you can change the binning as so:
```
asdf = ASDFChangeBinning(asdf, 1);
```

# MATLAB/C Library Functions #
| **Function** | **Description** |
|:-------------|:----------------|
| [ASDFChangeBinning](Documentation#ASDFChangeBinning.md) | Change ASDF binning size |
| [ASDFChooseTime](Documentation#ASDFChooseTime.md) | Subsample from ASDF by time |
| [ASDFGetfrate](Documentation#ASDFGetfrate.md) | Get average firing rate for each neuron |
| [ASDFSubsample](Documentation#ASDFSubsample.md) | Subsample neurons from ASDF |
| [ASDFTE](Documentation#ASDFTE.md) | Calculate Transfer Entropy |
| [ASDFToSparse](Documentation#ASDFToSparse.md) | Convert ASDF to a sparse matrix |
| [SparseToASDF](Documentation#SparseToASDF.md) | Convert a sparse matrix to ASDF |


---


## ASDFChangeBinning ##
<a />
**Call Signature:**
> `new_asdf = ASDFChangeBinning(asdf, factor)`

**Parameters:**
> `asdf` - {n\_neu,n\_bin} the old ASDF binned with different bin size

> `factor` - (1,1) the number of bins that you want to consider as one bin in the new ASDF. Must be an integer.

**Returns:**
> `new_asdf` - {n\_neu,n\_bin\_new} the new ASDF binned with different bin size.

**Description:**
> Change binning size. New ASDF will have the size `n_bin_new = floor(n_bin / factor)`. Since it "bins" the spike times, the resulting timing becomes an integer (and unique).

---

## ASDFChooseTime ##
<a />
**Call Signature:**
> `cutasdf = ASDFChoosetime(asdf, startTime, endTime)`

**Parameters:**
> `asdf` - {nNeu+2,1} ASDF to subsample time from

> `startTime` - (scalar) Starting time of new ASDF. If only two arguments are given, this is the time that you wanna choose.

> `endTime` - (scalar) Ending time of new ASDF

**Returns:**
> `cutasdf` - {nNeu+2,1} ASDF with new time and duration

---

## ASDFGetfrate ##
<a />
**Call Signature:**
> `frate = ASDFGetfrate(asdf)`

**Parameters:**
> `asdf` - {nNeu+2,1} ASDF from which to calculate firing rate

**Returns:**
> `frate` - (nNeu,1) Probability of firing at one bin for each neuron

---

## ASDFSubsample ##
<a />
**Call Signature:**
> `subasdf = ASDFSubsample(asdf, subIndex)`

**Parameters:**
> `asdf` - {nNeu+2,1} ASDF to subsample neurons from

> `subIndex` - (subNeu,1) indices of subsampled neurons

**Returns:**
> `subasdf` - {subNeu+2,1} subsampled ASDF

---

## ASDFTE ##
<a />
**Call Signature:**
> `[te_result, ci_result, all_te] = ASDFTE(asdf, j_delay, i_order, j_order, windowsize)`

**Parameters:**
> `asdf` - Time series in Another Spike Data Format (ASDF)

> `j_delay` - Number of bins to lag sender (j series) or a vector [1](default.md)

> `i_order` - Order of receiver [1](default.md)

> `j_order` - Order of sender [1](default.md)

> `windowsize`  - window size used for Coincidence Index (CI) calculation (odd number only)

**Returns:**
> `te_result` - (nNeu, nNeu) NxN matrix where N(i, j) is the peak transfer entropy from j->i over all delays

> `ci_result` - (nNeu, nNeu) NxN matrix where N(i, j) is the CI transfer entropy from j-> over all delays with windowsize

> `all_te` - (nNeu, nNeu, j\_delays) NxNxD matrix where N(i, j, d) is the transfer entropy from j->i at delay number d, where d ranges from 1:length(j\_delay)


---

## ASDFToSparse ##
<a />
**Call Signature:**
> `[raster, binunit] = ASDFToSparse(asdf)`

**Parameters:**
> `asdf` - {nNeu+2, 1} ASDF to convert

**Returns:**
> `raster` - (n\_neu, duration) (Sparse) time raster expressed as a sparse matrix.

> `binunit` - (string) the unit of time in the data (length of a bin in real scale)

**Description:**
> This function converts the time raster of ASDF to sparse matrix.


---

## SparseToASDF ##
<a />
**Call Signature:**
> `asdf = SparseToASDF(raster, binunit)`

**Parameters:**
> `raster` - (n\_neu, duration) (Sparse) time raster expressed as a sparse matrix.

> `binunit` - (scalar, double) the unit of time in the data in ms (length of a bin in real scale)

**Returns:**
> `asdf` - {n\_neu + 2, 1} ASDF version of the data

**Description:**
> This function converts sparse matrix version of the data to ASDF version.