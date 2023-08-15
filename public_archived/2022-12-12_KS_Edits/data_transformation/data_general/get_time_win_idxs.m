function time_idxs = get_time_win_idxs(sub_window,full_window,srate)
%% get_time_win_idxs: get indeces corresponding to a sub-window of a time 
%   range and sample rate
%
% Inputs:
%   sub_window  - sub-window of full_window to find indeces for, in the 
%                 form [t(1) t(2)] (ms)
%   full_window - time window containing sub-window (e.g. segmentation 
%                 window), in the form [t(1) t(2)] (ms)
%   srate       - sample rate of full_window data (Hz)
%
% Outputs:
%   time_idxs   - double containing indeces of each sample in sub_window,
%                 indicating position of sub_window relative to beginning
%                 of full_window, sampled at srate

%%
% get times corresponding to each sample of full time window
sample_times = get_sample_times(full_window,srate); 

% get indeces corresponding to sub_window
time_idxs = get_indexes(sub_window,sample_times,1);

end