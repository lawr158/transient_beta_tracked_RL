clear
%% inputs
% batch information
task = '02_Baseline';
beapp_tag = '_tb_happe';
tbe_batch = 'beta_13to30Hz_0to5000ms_method2';
% features to extract
F_in.channels = [36 104];
F_in.features = {'rate','power','duration','ts_duration','Fspan'};
% filename/id pattern
id_pat = digitsPattern(4);

%% code
% find dataset
P = load_paths_gp_spa; 
P.data = fullfile(P.tasks,task,['spectral_events' beapp_tag],tbe_batch);
D = dir( fullfile(P.data,'*.mat') ); 
% find relevant counts
n_feat = length(F_in.features);
n_file = length(D);
% iterate through files
for i_f = 1:n_file
    % load file event data
    f = D(i_f); load( fullfile(f.folder,f.name),'SE_ts' )
    f.id = extract(f.name(1:6),id_pat); 
    % populate subject row on output table
    T(i_f).id = f.id{1}; T(i_f).filename = f.name;
    % extract features from channels of interest
    for i_c = F_in.channels
        if i_f == 1
            T = extract_channel_features(SE_ts(i_c),F_in,n_feat,T);
        else
            T(i_f) = extract_channel_features(SE_ts(i_c),F_in,n_feat,T(i_f));
        end
    end
end

%% helper functions
function F_out = extract_channel_features(C,F_in,n_feat,F_out)
    ch = ['TBE_' C.channel '_event_'];
    for i_f = 1:n_feat
        feat = F_in.features{i_f};
        F_out.([ch feat '_mean']) = C.mean.(feat);
        F_out.([ch feat '_sd'])   = C.sd.(feat);
    end
end