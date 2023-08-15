function [data_mean, data_var] = subj_avg_spec_events( sorted_data,     ...
                                                      feats, conds_to_inc )
%% SUBJ_AVG_SPEC_EVENTS: Finds average of each inputted feature per subject 
% for sorted spec_events structs
%   Inputs:
%     - sorted_data - spec_events struct run through sort_data_by_group_old
%     - feats       - cell of features to extract from struct
%     - conds_to_inc - array of #s of conditions to include, 0 for all
%   Outputs:
%     - data_mean - subject averages of input features, sorted by group
%     - data_var  - subject stdevs of input features, sorted by group
%
%%
% create conds_to_inc if not input
if ~exist('conds_to_inc','var'), conds_to_inc = 0; end

% get number of groups and features
groups  = fieldnames(sorted_data);
n_grps  = length(groups);
n_feats = length(feats);
n_chans = length(sorted_data);
% iterate through channels in data
for i_c = 1:n_chans
    % iterate through groups in data
    for i_g = 1:n_grps
        group    = groups{i_g};
        subjects = fieldnames(sorted_data(i_c).(group));
        n_subj   = length(subjects);
        % iterate through subjects
        for i_s = 1:n_subj
            subj = subjects{i_s};
            subj_data = sorted_data(i_c).(group).(subj);
            % intialize output structs
            if i_s == 1
                data_mean(i_c).(group).id = {};
                data_var(i_c).(group).id  = {};
                for i_f = 1:n_feats
                    data_mean(i_c).(group).(feats{i_f}) = [];
                    data_var(i_c).(group).(feats{i_f}) = [];
                end
            end
            data_mean(i_c).(group).id{end+1,1} = subj(2:end);
            data_var(i_c).(group).id{end+1,1} = subj(2:end);
            % iterate through features
            for i_f = 1:n_feats
                feat = feats{i_f};
                %% add subject features and feature means to arrays
                % trial-level features
                if isequal(feat,'eventnumber')
                    % NaN for trial-less subjects
                    if isempty(subj_data.TrialSummary)
                        data_mean(i_c).(group)(end+1,1) = NaN;
                        data_var(i_c).(group)(end+1,1)  = NaN;
                        continue
                    end
                    % get trial features
                    trial = subj_data.TrialSummary.TrialSummary;
                    trial_feats = trial.(feat);
                    % if input, select specific conditions to average
                    if conds_to_inc > 0
                        trial_matches = trial.classLabels == conds_to_inc;
                        trial_matches = logical(sum(trial_matches,2));
                        trial_feats   = trial_feats(trial_matches);
                    end
                    feat_mean = mean(trial_feats,1,'omitnan');
                    feat_std  = std(trial_feats,1,'omitnan');
                    % save mean & stdev
                    data_mean(i_c).(group).(feat)(end+1,1) = feat_mean;
                    data_var(i_c).(group).(feat)(end+1,1)  = feat_std;
                
                % event-level features
                else
                    % NaN for event-less subjects
                    if isempty(subj_data.Events)
                        data_mean(i_c).(group)(end+1,1) = NaN;
                        data_var(i_c).(group)(end+1,1)  = NaN;
                        continue
                    end
                    % get subject events
                    event = subj_data.Events.Events;
                    event_feats = event.(feat);
                    % if input, select specific conditions to average
                    if conds_to_inc > 0
                        trial_matches = event.classLabels == conds_to_inc;
                        trial_matches = logical(sum(trial_matches,2));
                        event_feats   = event_feats(trial_matches);
                    end
                    if strcmp(feat,'duration')
                        % convert s -> ms
                        event_feats = event_feats * 1000;
                    end
                    feat_mean = mean(event_feats,1,'omitnan');
                    feat_std  = std(event_feats,1,'omitnan');
                    % save mean & stdev
                    data_mean(i_c).(group).(feat)(end+1,1) = feat_mean;
                    data_var(i_c).(group).(feat)(end+1,1)  = feat_std;

                end     % feature type IF
            end     % feature loop
        end     % subject loop
    end     % group loop
end    % channel loop

%% average outputs of multiple channels together
if n_chans > 1
    for i_g = 1:n_grps
        group = groups{i_g};
        chan_mean.(group).id = data_mean(1).(group).id;
        chan_var.(group).id  = data_var(1).(group).id;
        for i_f = 1:n_feats
            feat = feats{i_f};
            % get feature means/sds and average across channels
            group_means = [data_mean(:).(group)];
            group_means = mean([group_means(:).(feat)],2,'omitnan');
            group_vars  = [data_var(:).(group)];
            group_vars  = mean([group_vars(:).(feat)],2,'omitnan');
            % save to new struct
            chan_mean.(group).(feat) = group_means;
            chan_var.(group).(feat)  = group_vars;
        end    % feature loop
    end    % group loop

    % update output struct with channel averages
    data_mean = chan_mean;
    data_var  = chan_var;
end    % number of channels IF

end