function path = load_paths_gp_spa()
    path.rcfs  = fullfile('X:','Groups','SPA');
    path.eeg   = fullfile(path.rcfs,'01_Data_Raw_Summary_Processed','EEG');
    path.tasks = fullfile( path.eeg,'Participant_Data',             ...
                           '03_Processed_Data','01_In_Progress' );
end