function batch_create_output_table(C,P,F)
% create output folder
if ~exist(C.path.out,'dir'), mkdir(C.path.out); end
% find dataset
D = dir( fullfile(P.path.out,'*.mat') ); 
% find relevant counts
n_feat = length(F.features);
n_file = length(D);
% iterate through files
for i_f = 1:n_file
    % load file event data
    f = D(i_f); load( fullfile(f.folder,f.name),'SE_ts' )
    f.id = extract(f.name(1:6),F.id_pat); 
    % populate subject row on output table
    T(i_f).id = f.id{1}; T(i_f).filename = f.name;
    % extract features from channels of interest
    for i_c = F.channels
        if i_f == 1
            T = extract_channel_feats(SE_ts(i_c),F,n_feat,T);
        else
            T(i_f) = extract_channel_feats(SE_ts(i_c),F,n_feat,T(i_f));
        end
    end
    clear SE_ts
end
% output table
T_out = struct2table(T); 
out_fname = ['spectral_events_' P.path.fname '.csv'];
writetable( T_out,fullfile(C.path.out_tables,out_fname),                ...
            'WriteRowNames',true )
end

%% helper functions
function F_out = extract_channel_feats(C,F_in,n_feat,F_out)
    ch = ['TBE_' C.channel '_event_'];
    for i_f = 1:n_feat
        feat = F_in.features{i_f};
        try
            F_out.([ch feat '_mean']) = C.mean.(feat);
            F_out.([ch feat '_sd'])   = C.sd.(feat);
        catch
            F_out([ch feat '_mean']) = NaN;
            F_out([ch feat '_sd']) = NaN;
        end
    end
end