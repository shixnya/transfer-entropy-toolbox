%==============================================================================
% Copyright (c) 2011, The Trustees of Indiana University
% All rights reserved.
% 
% Authors: Michael Hansen (mihansen@indiana.edu), Shinya Ito (itos@indiana.edu)
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
%   1. Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
% 
%   2. Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
% 
%   3. Neither the name of Indiana University nor the names of its contributors
%      may be used to endorse or promote products derived from this software
%      without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%==============================================================================

spikes = [
           1 1 1 0 1 1 0 1 0 0 1 1 0 0 1 1 1 0 1 1 0 1; % neuron 1
           1 0 1 0 1 0 1 0 1 1 0 1 1 0 1 1 0 1 1 0 0 1  % neuron 2
         ];

asdf = SparseToASDF(spikes, 1); % Use a bin size of 1

te_1 = ASDFTE(asdf); % Delay 1
fprintf('Delay 1 1st Order TE\n');
te_1

d = 10;
te_d = ASDFTE(asdf, d); % Delay d
fprintf('\nDelay %d 1st Order TE\n', d);
te_d

ds = 1:10
te_ds = ASDFTE(asdf, ds); % Multiple delays, max TE value
fprintf('\nDelays %d-%d 1st Order TE max\n', 1, length(ds));
te_ds

window = 2;
te_ci = ASDFTE(asdf, ds, 1, 1, @CIReduce, window); % Multiple delays, coincidence index (window 2)
fprintf('\nDelays %d-%d 1st Order Coincidence Index (%d)\n', 1, length(ds), window);
te_ci

i_order = 2;
j_order = 3;
te_ho = ASDFTE(asdf, d, i_order, j_order); % Delay d, higher order
fprintf('\nDelay %d Higher Order TE\n', d);
te_ho

te_ho_ci = ASDFTE(asdf, ds, i_order, j_order, @CIReduce, window); % Multiple delays, higher order, coincidence index (window 2)
fprintf('\nDelays %d-%d Higher Order Coincidence Index (%d)\n', 1, length(ds), window);
te_ho_ci

