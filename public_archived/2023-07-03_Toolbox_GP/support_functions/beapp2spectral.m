% Author: Gerardo Parra, Last updated: 2023-03-28
function X = beapp2spectral(C,EEG,FPI)
%% beapp2spectral: convert BEAPP segmented data to spectral events format
% Inputs:
%   - C - struct | set of conversion parameters, with fields:
%       - t_win.samp.i - double array | indeces of data samples to keep
%       - t_win.samp.n - double | number of samples being kept
%   - EEG - cell array | beapp output data (e.g. eeg_w)
%   - FPI - struct | beapp file_proc_info, necesssary fields:
%       - net_vstruct - struct with EEG net information
%       - beapp_nchan - number of channels in data
% Outputs: for each channel i_c in EEG:
%   - X(i_c) - struct with with fields
%       - channel - channel name
%       - data  - matrix with channel data concatenated across conditions
%       - class - vector with numeric classification of trials in channel
 
%% intialize output variables
% create row for each channel in X and trial class outputs
X = struct('channel',{FPI.net_vstruct.labels},'data',NaN,'class',NaN); 
n_cond = size(EEG,3);
% get total number of trials across conds to initialize output matrix
trials.n = cellfun(@(x) size(x,3),EEG,'UniformOutput',true);
trials.total = sum(trials.n); cond.t1_i = [1 trials.n];
%% iterate through channels
for i_ch = 1:FPI.beapp_nchan
    % skip if channel is all NaNs
    if isnan_channel(EEG,i_ch), continue; end
    % initialize channel output matrices
    X(i_ch).data  = NaN(C.t_win.samp.n,trials.total);
    X(i_ch).class = NaN(1,trials.total);
    %% add converted trials to output matrix by condition
    for i_c = 1:n_cond
        % get trial indices for current condition
        cond.t_start = sum(cond.t1_i(1:i_c));
        cond.t_end   = cond.t_start + cond.t1_i(i_c+1) - 1;
        cond.t_range = cond.t_start:cond.t_end; 
        % add channel trials to output matrix
        cond.data = squeeze(EEG{i_c}(i_ch,C.t_win.samp.i,:));
        X(i_ch).data(:,cond.t_range) = squeeze(cond.data);
        % add condition numbers to output class field 
        X(i_ch).class(cond.t_range)  = i_c;
    end
end

end

%% helper functions
% returns whether channel data is NaN
function isnan_ch = isnan_channel(eeg,chan)
    cf_isnan_ch = @(X) all(isnan(squeeze(X(chan,:,1))));
    isnan_ch = cellfun(cf_isnan_ch,eeg,'UniformOutput',true);
end
% returns whether any of selected channels have NaN data
function has_nan_chan = has_nan_channels(eeg,chans)
    % cell function for number of non-NaN channels
    nonNans = @(x) size(rmmissing(x(chans,:,1)),1);
    % get number of good chans in each condition
    n_good_chans = cellfun(nonNans,eeg,'UniformOutput',true);
    % 1 if any condition has fewer good channels than size of input chans
    has_nan_chan = any(n_good_chans<length(chans));
end