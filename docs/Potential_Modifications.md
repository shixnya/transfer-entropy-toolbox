  * Output full TE delay (without reducing) by ASDFTE.m, or output max TE and CI at the same time. (it is redundant to compute TE delay twice to get both.)
  * ASDFTE Checks the argument of @CIReduce after calculating TE. Shouldn't it check beforehand?
  * When type 'help ASDFTE', I wanna see how to use it.
  * Using Matlab style argument (not very important) (e.g. ASDFTE(asdf, 'JDelay', 1:30, 'ReduceMethod', @max); )
  * We should use the same indentation.