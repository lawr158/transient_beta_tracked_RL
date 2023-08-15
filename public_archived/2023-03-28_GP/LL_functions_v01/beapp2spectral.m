% Author: Gerardo Parra, Last updated: 2023-03-28
function [C,X,T] = beapp2spectral(C,channels)
%% beapp2spectral: convert BEAPP segmented data to spectral events format
% Inputs:
%   - C - struct | set of conversion parameters, with fields:
%       - rerun     - logical | toggle re-conversion or load existing data
%       - tag       - char | beapp tag to load data from
%       - path.in   - char | path to directory with task beapp folders
%       - path.out  - char | optional, path to save output data to. default
%                     is: [path.in/[spectral_events C.tag]/inputs]
%       - t_win.to_keep - double array | optional, time window to keep;
%                         relative to stim onset for task-related; relative 
%                         to original epochs for baseline. default keeps
%                         entire original segments
%   - channels - double | index of channel(s) to convert
% Outputs:
%   - C - struct | input parameters updated with batch metadata
%   - for each BEAPP file in path.in i_a:
%       - X{i_a} - [time x trial] matrix of converted data
%       - t_classes{i_a} - numeric classification of trials in X{i_a}

%% load data and iterate through files
% get data directory struct
C.file.dir = dir( fullfile(C.path.in,['segment' C.tag],'*.mat') );
% try HAPPE_V3 folder if directory is emtpy
if isempty(C.file.dir)
    C.file.dir = dir( fullfile(C.path.in,['HAPPE_V3' C.tag],'*.mat') );
end
C.file.n = size(C.file.dir,1);
% initialize output variables
X = []; T = [];
% iterate through files
for i_f = 1:C.file.n
    %% load EEG/metadata for current file
    F = C.file.dir(i_f); F.id = str2double(F.name(1:4));
    load(fullfile(F.folder,F.name),'eeg_w','file_proc_info')
    EEG = file_proc_info; EEG.data = eeg_w; clear eeg_w file_proc_info
    % skip if any condition has less than 9 trials
    if any(get_n_trials(EEG.data)<9), continue; end
    % skip if any of selected channels are NaNs
    if has_nan_channels(EEG.data,channels), continue; end
    %% intialize output variables
    EEG.cond.n = size(EEG.data,3);
    % get total number of trials across conds to initialize output matrix
    EEG.trial.n = cellfun(@(x) size(x,3),EEG.data,'UniformOutput',true);
    EEG.trial.total = sum(EEG.trial.n); cond.t1_i = [1 EEG.trial.n];
    % initialize output matrices
    X(end+1).id = F.id; T(end+1).id = F.id;
    X(end).data  = NaN(C.t_win.samp.n,EEG.trial.total);
    T(end).class = NaN(1,EEG.trial.total);
    %% add converted trials to output matrix by condition
    for i_c = 1:EEG.cond.n
        % get trial indices for current condition
        cond.t_start = sum(cond.t1_i(1:i_c));
        cond.t_end   = cond.t_start + cond.t1_i(i_c+1) - 1;
        cond.t_range = cond.t_start:cond.t_end; 
        % average channels across time window, add trials to output matrix
        cond.data = mean(EEG.data{i_c}(channels,C.t_win.samp.i,:),1);
        X(end).data(:,cond.t_range) = squeeze(cond.data);
        T(end).class(cond.t_range)  = i_c;
    end
end

end

%% helper functions
% returns array with number of trials/epochs in each eeg_w cell
function n_t = get_n_trials(eeg)
    n_t = cellfun(@(x) size(x,3),eeg,'UniformOutput',true);
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