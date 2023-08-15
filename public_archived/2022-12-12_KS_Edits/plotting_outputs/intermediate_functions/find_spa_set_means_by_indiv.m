function processed_data = find_spa_set_means_by_indiv(data,features)

% get data information for iterating
groups   = fieldnames(data);
n_groups = length(groups);
n_feats  = length(features{1});

% iterate through groups
for i_g = 1:n_groups
    % get group information
    group   = groups{i_g};
    subjs   = fieldnames(data.(group));
    n_subjs = length(subjs);
    
    % iterate through subjects
    for i_s = 1:n_subjs
        % get subject information
        subj_id = subjs{i_s}(2:5);
        subj_struct_id = ['s' subj_id];
        subj_data = data.(group).(subj_struct_id);
        
        if subj_data.TrialSummary.NumTrials < 10
            data.(group) = rmfield(data.(group),subj_struct_id);
            continue
        end
        
        % iterate through subject features
        for i_f = 1:n_feats
            curr_feat = features{1}{i_f};
            if strcmp(curr_feat,'eventnumber')  
                % get features at trial level
                try
                    trials = subj_data.TrialSummary.TrialSummary;
                catch
                    continue
                end
            else  
                % get features at event level
                try
                    trials = subj_data.Events.Events;
                catch
                    continue
                end
            end
            % get condition/class assignments of trials
            cond_labels = trials.classLabels;
            n_conds     = max(cond_labels);
            
            % iterate through conditions
            feat_mean = zeros(1,n_conds);
            feat_var  = zeros(1,n_conds);
            for i_c = 1:n_conds
                % find trials matching current condition
                matching_trials = cond_labels == i_c;
                matching_trials = logical(sum(matching_trials,2));
                cond_trials     = trials.(curr_feat);
                
                % find mean/var for current condition and save
                cond_mean = nanmean(cond_trials(matching_trials),1);
                cond_var  = nanvar(cond_trials(matching_trials),1);
                if strcmp(features{1}{i_f},'duration')
                    % note: convert s -> ms
                    cond_mean = cond_mean * 1000;
                    cond_var  = cond_var * 1000;
                end
                try
                    feat_mean(i_c) = cond_mean;
                catch
                    break
                end
                feat_var(i_c)  = cond_var;
            end
            
            % save feature means to data structure
            data.(group).(subj_struct_id).([curr_feat '_mean']) = feat_mean;
            data.(group).(subj_struct_id).([curr_feat '_var'])  = feat_var;
        end
    end
end

% output results
processed_data = data;

end