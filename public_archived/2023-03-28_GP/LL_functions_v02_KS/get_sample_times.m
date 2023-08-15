function sample_times = get_sample_times(time_window, srate)
%% get_sample_times: get times corresponding to samples of a given time
%   window and sample rate
%
% Inputs:
%   time_window  - time window to find sammple times for, in the form
%                 [t(1) t(2)], in milliseconds
%   srate        - sample rate to give times in (Hz)
%
% Outputs:
%   sample_times - double containing times corresponding to samples in the 
%                  time range given in time_window
%
%%
sample_times = linspace(time_window(1), time_window(2), ...
                        srate*sumabs(time_window)/1000);

end