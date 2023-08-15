% Author: Gerardo Parra, Last udpated: 2023-03-28
function C = convert_beapp_by_chan(C)
%% convert_beapp_by_chan: convert BEAPP data to spectral events by channel
% Inputs:
%   - C - struct | set of conversion parameters, with fields:
%       - tag       - char | beapp tag to load data from
%       - path.in   - char | path to directory with task beapp folders
%       - path.out  - char | optional, path to save output data to. default
%                     is: [path.in/[spectral_events C.tag]/inputs]
%       - chan.list - double array | optional, channels to convert data for
%                     default: all 10-20 channels
%       - t_win.to_keep - double array | optional, time window to keep;
%                         relative to stim onset for task-related; relative 
%                         to original epochs for baseline. default keeps
%                         entire original segments
%       - rerun     - logical | toggle re-conversion or load existing data
% Outputs:
%   - C - struct | input parameters updated with batch metadata
%   - saves a .mat file with variables X and t_classes to C.path.out, with 
%     the following structure for each participant i_a:
%       - X{i_a} - [time x trial] matrix of converted data
%       - t_classes{i_a} - numeric classification of trials in X{i_a}

%% check optional inputs
% create fields for optional inputs
if ~isfield(C.path,'out'), C.path.out = ''; end
if ~isfield(C,'chan'), C.chan.list = []; end
if ~isfield(C,'t_win'), C.t_win.to_keep = []; end
if ~isfield(C,'rerun'), C.rerun = 1; end
% set output folder, create directory if necessary
if isempty(C.path.out)
    C.path.out = fullfile(C.path.in,['spectral_events' C.tag]);
end
if ~exist(C.path.out,'dir'), mkdir(C.path.out), end
% set channels to process
if isempty(C.chan.list)
    C.chan.list = load_1020_chans;
end
C.chan.n = length(C.chan.list);

%% get batch metadata
% get original segmentation window from batch beapp output file
[C.t_win.beapp,~,C.Fs] = get_batch_segment_window(C.path.in,C.tag);
C.t_win.beapp = C.t_win.beapp * 1000;
if isempty(C.t_win.to_keep), C.t_win.to_keep = C.t_win.beapp; end
% get indeces of C.window.to_keep, relative to seg_win
C.t_win.samp.i = get_time_win_idxs(C.t_win.to_keep,C.t_win.beapp,C.Fs);
C.t_win.samp.n = length(C.t_win.samp.i);
% get time window string for filename
C.t_win.str = get_range_str(C.t_win.to_keep,'ms');
% get string of channels for filename

%% iterate through channels
if ~C.rerun, return; end
for i_c = 1:C.chan.n
    chan = C.chan.list(i_c);
    % convert channel data
    [C,X,T] = beapp2spectral(C,chan);
    % save data
    save(get_converted_fpath(C,chan),'X','T','C')
end

end