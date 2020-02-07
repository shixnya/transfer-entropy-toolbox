# Simple Demo #
_Sorry about the syntax highlighting. MATLAB does not appear to be supported._
```
%% 1. Prepare some spike trains
% In this example, train2 tends to fire 1 bin after train1.
% Also, train1 tends to fires 2 bins after train2.
disp('Here are some sample trains:')
train1 = [0 1 0 0 1 1 0 1 0 0 0 0 1 0 0 1]
train2 = [0 0 1 0 0 0 1 1 0 0 1 0 0 1 0 0]

%% 2. Combine trains together to make a matrix representation
trains = [train1 ; train2];

%% 3. Make proper data format called ASDF (Another Spiking Data Format)
asdf = SparseToASDF(trains, 1); % use 1 ms bin

disp('Calculating TE for delays 1-3');
[peakTE, CI, TEdelays] = ASDFTE(asdf, 1:3) % using all the pairs in asdf with delay of 1 to 3 bins
% peakTE is a peak value over delay. CI is coincidence index with a default window size of 5 bins.
% TEdelays is a 3-d matrix of (sending train, receiving train, delay).
% TEdelays(1,2,1) and TEdelays(2,1,2) should be high due to the embedded patterns in the sample data.

%% 4. Some Higher Order TE calculations
disp('You can also calculate Higher Order TE with arbitrary delays:');
TE133 = ASDFTE(asdf, 1, 3, 3) % Higher order TE with delay 1, order k=l=3
TE333 = ASDFTE(asdf, 3, 3, 3) % Higher order TE with delay 3, order k=l=3
TE1to312 = ASDFTE(asdf, 1:3, 1, 2) % Higher order TE with delay 1 to 3, order k=1, l=2

%% 5. Let's do the same thing with large spike train.
% One of the simulations from the paper. 100 subsampled neurons with 1.8 million bins. (30 minutes with 1ms bins)
load Izhik_100_0

disp('Calculating TE with 30 different delays for the 30 minute data... it will take a few minutes.')
tic
[peakTE, CI, TEdelays_big] = ASDFTE(asdf, 1:30); % Delays of 1ms to 30ms
toc
% Using one core of a Core 2 duo 2.66GHz, it takes about 135 seconds.

%% 6. Show scatter plots
% removing self connections
peakTE = peakTE - diag(diag(peakTE));
CI = CI - diag(diag(CI));
% plot
from0to10 = find(conmat>0 & conmat<10);
subplot(1,2,1); semilogy(conmat(from0to10), peakTE(from0to10),'r.'); title 'weight vs TEPk'
subplot(1,2,2); plot(conmat(from0to10), CI(from0to10),'r.'); title 'weight vs TECI'

%% 7. Check if the values are consistent with expected values.
cond = 0; % no error
if TEdelays(1,1,1) ~= 0
	cond = 1;
elseif abs(TEdelays(1,2,1) - 0.1953) > 0.001
	cond = 1;
elseif abs(TE333(1,2) - 0.5009) > 0.001
	cond = 1;
end

if cond == 1
	error('Some of your values are not correct. Please check settings!');
end
```