function path = load_paths()
    % helper paths for below
    path.rcfs = fullfile('X:','Groups','Transient Beta');
    path.eeg  = fullfile(path.rcfs,'01_Data','01_EEG');
    %% these are the important paths
    % path to beapp files
    path.beapp_files = fullfile(path.eeg,'01_BEAPP_Files');
    % BEAPP run tag
    path.beapp_tag = '_run_tag';
    % desired path to output tables to
    path.out_tables = fullfile(path.eeg,'02_Output_tables');
end