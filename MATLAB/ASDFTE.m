% Parameters:
%   asdf        - Time series in Another Spike Data Format (ASDF)
%   j_delay     - Number of bins to lag sender (j series) or a vector [default 1]
%   i_order     - Order of receiver [default 1]
%   j_order     - Order of sender [default 1]
%   reduce      - Function to reduce transfer entropies across multiple delays [default max]
%   reduce_args - Vector of arguments passed to the reduce function [default empty]
%
% Returns:
%   te_result - (nNeu, nNeu) NxN matrix where N(i, j) is the transfer entropy from j->i

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

function te_result = ASDFTE(asdf, j_delay, i_order, j_order, reduce, reduce_args)
	% Set defaults
	if nargin < 2
		        j_delay = 1;
			end

			if nargin < 3
				        i_order = 1;
					end

					if nargin < 4
						        j_order = 1;
							end

							if nargin < 5
								        reduce = @max;
									end

									if length(j_delay) == 1
										        te_result = transent(asdf, j_delay, i_order, j_order); % Single delay
											else

												        % Multiple delays
														        num_delays = length(j_delay);
																        info = asdf{end};
																		        num_neurons = info(1);

																				        % Allocate space for all matrices
																						        all_te = zeros(num_neurons * num_delays, num_neurons);

																								        % Compute TE for delay times
																										        row = 1;
																												        for d = j_delay
																															                all_te(row:(row + num_neurons - 1), :) = transent(asdf, d, i_order, j_order);
																																			                row = row + num_neurons;
																																							        end

																																									        % Reduce to final matrix
																																											        te_result = zeros(num_neurons, num_neurons);

																																													        for i = 1:num_neurons
																																																                rows = i + ((0:(num_delays - 1)) * num_neurons);
																																																				                for j = 1:num_neurons
																																																									                        if nargin < 6
																																																																                                te_result(i, j) = reduce(all_te(rows, j));
																																																																								                        else
																																																																															                                te_result(i, j) = reduce(all_te(rows, j), reduce_args);
																																																																																							                        end
																																																																																													                end
																																																																																																	        end

																																																																																																		end % if multiple delays


