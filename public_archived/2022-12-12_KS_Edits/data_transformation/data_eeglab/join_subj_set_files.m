function [eeg_w,subj_info] = join_subj_set_files(subj_id,tag,data_path)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
echo('off','all');

% get subject file structure/info
subj_info_dir  = dir(fullfile(data_path,[subj_id '*.set']));
subj_info_file = subj_info_dir(1);
subj_info_path = fullfile(subj_info_file.folder,subj_info_file.name);
subj_info = pop_loadset('filename',subj_info_path,'loadmode','info');
clear subj_info_dir subj_info_file subj_info_path

% get data matrices from all condition /set files
subj_dir = dir(fullfile(data_path,[subj_id '*' tag '*.set']));
n_files  = length(subj_dir);
eeg_w    = cell(n_files,1);
for i_f = 1:n_files
    file = subj_dir(i_f);
    eeg  = pop_loadset(fullfile(file.folder,file.name));
    eeg_w{i_f} = eeg.data;
end

end