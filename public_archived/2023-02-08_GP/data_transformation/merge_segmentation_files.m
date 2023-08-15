%% merge_segmentation_files.m
% For use when same task has been segmented through multiple BEAPP runs and
% different files need to be merged into a single BEAPP compatible file
% The following process was used to generate multiple BEAPP segment files
% that needed to be joined:
%   1. 
% Author: Gerardo Parra
% -------------------------------------------------------------------------

% set data location
data_path = ['X:\Groups\SPA\Levin Lab Internal\Data\EEG\' ...
                        'Participant Data\Processed Data'];
task_dir  = [data_path filesep 'auditory_temporal_habituation'];

% create Matlab dir with data
data_dirs = dir([task_dir filesep 'segment_DI*']);
n_dirs    = size(data_dirs,1);
new_eeg_w{n_dirs,1} = [];
new_fpi = [];

% iterate through directories containing output
for i_d = 1:n_dirs
    stim_struct = data_dirs(i_d);
    stim_path   = [stim_struct.folder filesep stim_struct.name ...
                                                filesep '*.mat'];
    stim_dir    = dir(stim_path);
    
    % load eeg file
    load([stim_dir.folder filesep stim_dir.name],'eeg_w','file_proc_info')
    % create empty file_proc_info struct
    if i_d == 1
        for fn = fieldnames(file_proc_info)'
           new_fpi.(fn{1}) = file_proc_info.(fn{1});
        end
    end
    
    % add EEG data to merged file
    new_eeg_w{i_d,1} = eeg_w{1,1};
    
    % update file_proc_info
    if i_d ~= 1
        new_fpi.beapp_bad_chans{1,1} = union(   ...
            new_fpi.beapp_bad_chans{1,1},       ...
            file_proc_info.beapp_bad_chans{1,1} ...
        );
        new_fpi.epoch{end+1,1} = file_proc_info.epoch{1,1};
        new_fpi.evt_conditions_being_analyzed(end+1,:) = ...
            file_proc_info.evt_conditions_being_analyzed(1,:);
        new_fpi.grp_wide_possible_cond_names_at_segmentation{1,end+1} = ...
          file_proc_info.grp_wide_possible_cond_names_at_segmentation{1,1};
    end
end
clear eeg_w file_proc_info

% save new eeg file
eeg_w = new_eeg_w;
file_proc_info = new_fpi;
out_dir = [task_dir filesep 'segment_seg'];
save([out_dir filesep stim_dir.name],'eeg_w','file_proc_info')
clear