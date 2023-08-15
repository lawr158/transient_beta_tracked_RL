% Author: Gerardo Parra, Last udpated: 2023-03-28
function C = batch_convert_beapp_data(C)
%% BATCH_CONVERT_BEAPP_DATA: convert BEAPP batch to spectral events format
% Inputs:
%   - C - struct | set of conversion parameters, with fields:
%       - tag       - char | beapp tag to load data from
%       - path.in   - char | path to directory with task beapp folders
%       - path.out  - char | optional, path to save output data to. default
%                     is: [path.in/[spectral_events C.tag]/inputs]\
%       - Fs - double | optional, sampling rate of input data. will use
%              beapp_rsamp_srate found in grp_proc_info otherwise
%       - t_win.to_keep - double array | optional, time window to keep;
%                         relative to stim onset for task-related; relative 
%                         to original epochs for baseline. default keeps
%                         entire original segments
%       - rerun     - logical | toggle re-conversion or load existing data
%       - trial_thr - double | optional, trial threshold for converting
%                     data. default 10
% Outputs:
%   - C - struct | input parameters updated with batch metadata
%   - saves a .mat file with variables X, file_proc_info, and C for each 
%     participant to C.path.out, with the following structure for each 
%     channel i_c in that subject's original data:
%       - X(i_c).channel - channel name
%       - X(i_c).data - [time x trial] matrix of converted channel data
%       - X(i_c).class - condition/classification of each trial in
%                        X(i_c).data

%% check optional inputs
% create fields for optional inputs
if ~isfield(C.path,'out'), C.path.out = ''; end
if ~isfield(C,'t_win'), C.t_win.to_keep = []; end
if ~isfield(C,'rerun'), C.rerun = 1; end
if ~isfield(C,'trial_thr'), C.trial_thr = 10; end
% set output folder, create directory if necessary
if isempty(C.path.out)
    C.path.out = fullfile(C.path.in,['spectral_events' C.tag],'inputs');
end
if ~exist(C.path.out,'dir'), mkdir(C.path.out), end

%% get batch metadata
% get original segmentation window from batch beapp output file
if ~isfield(C,'Fs'), C.Fs = []; end
if isempty(C.Fs)
    warning(['No sampling rate was included in parameters. Using ' ...
            'default sampling rate of 1000 Hz'])
    C.Fs = 1000; 
end
[C.t_win.beapp,~,C.Fs] = get_batch_segment_window(C.path.in,C.tag,C.Fs);
C.t_win.beapp = C.t_win.beapp * 1000;
if isempty(C.t_win.to_keep), C.t_win.to_keep = C.t_win.beapp; end
% get indeces of C.window.to_keep, relative to seg_win
C.t_win.samp.i = get_time_win_idxs(C.t_win.to_keep,C.t_win.beapp,C.Fs);
C.t_win.samp.n = length(C.t_win.samp.i);
% get time window string for filename
C.t_win.str = get_range_str(C.t_win.to_keep,'ms');

%% load data and iterate through files
% return if rerun off
if ~C.rerun, return; end
% get data directory struct
C.file.dir = dir( fullfile(C.path.in,['segment' C.tag],'*.mat') );
% try HAPPE_V3 folder if directory is emtpy
if isempty(C.file.dir)
    C.file.dir = dir( fullfile(C.path.in,['HAPPE_V3' C.tag],'*.mat') );
end
C.file.n = size(C.file.dir,1);
% iterate through files
for i_f = 1:C.file.n
    %% load and convert data
    % load EEG/metadata for current file
    F = C.file.dir(i_f); %F.id = str2double(F.name(1:4));
    %%%%%%% KS addition to allow for my particiapnts to be run 6/30/23 %%%%%%%
    id = strrep(F.name,'-','_');
    id = strrep(id,' ','_');
    id = erase(id,'.mat');
    F.id = str2double(id);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(fullfile(F.folder,F.name),'eeg_w','file_proc_info')
    % skip if any condition has less trials than threshold
    if any(get_n_trials(eeg_w)<C.trial_thr), continue; end
    % get conversion to spectral events
    X = beapp2spectral(C,eeg_w,file_proc_info);
    %% save converted data
    save(fullfile(C.path.out,F.name),'C','X','file_proc_info')
    clear X eeg_w file_proc_info
end

end

%% helper functions
% returns array with number of trials/epochs in each EEG cell
function n_t = get_n_trials(eeg)
    n_t = cellfun(@(x) size(x,3),eeg,'UniformOutput',true);
end