function [figs, means, p_vals] = plot_SET_by_indiv_condition(data, groups,    ...
                                    feats, taps_to_inc,evt_per_sec,plot_title,groups_to_plot)
%% plot_SET_scatterplots: Create scatterplots from spectral events data

% inputs
%   - data:      2xn cell, for n groups, containing spec_events structs in
%                first row, other measure to plot in second
%   - groups:    1xn cell containing group names
%   - feats:     2xm cell, row 1 with m fields to extract from spec_events
%                structs; row 2 with m names describing each field
%   - equal_Ns:  0 or 1 to indicate whether to equalize sample sizes by
%                randomly selecting elements from larger datasets
%
% outputs
%   - figs:     matlab figures with scatterplots for each feature
%   - means:    table containing mean and median values of each feature   
%   - p_vals:   table containing significance of group differences

%% get data

% get data information for iterating
groups   = fieldnames(data);
n_groups = length(groups);
n_feats  = length(feats{1});

% initialize data structures
feat_mean_by_subj = cell(n_groups,n_feats);
feats_by_group    = cell(n_groups,n_feats);

% iterate through groups in data
for i_g = 1:n_groups
    % get group information
    group   = groups{i_g};
    subjs   = fieldnames(data.(group));
    n_subj  = length(subjs);
    
    % iterate through features
    for i_f = 1:n_feats
        grp_feats = [];

        % iterate through subjects
        for i_s = 1:n_subj
            % get subject information
            subj_id = subjs{i_s}(2:5);
            subj_struct_id = ['s' subj_id];
            subj_data = data.(group).(subj_struct_id);
            
            % add subject features and feature means to arrays
            if isequal(feats{1}{i_f},'eventnumber')  % trial-level features
                if isempty(subj_data.TrialSummary)
                    continue
                end
                
                trial = subj_data.TrialSummary.TrialSummary;
                % find feature average by condition
                conds   = trial.classLabels; % get cond labels
                n_conds = max(conds);
                subj_mean = zeros(1,n_conds);
                % iterate through conds
                for i_c = 1:n_conds
                    % find trials matching current cond
                    matching_trials = conds == i_c;
                    matching_trials = logical(sum(matching_trials,2));
                    trial_feats     = trial.(feats{1}{i_f});
                    % find mean for current cond and save
                    cond_mean       = nanmean(trial_feats(matching_trials),1);
                    subj_mean(i_c)  = cond_mean;
                end
                
            else  % event-level features
                if isempty(subj_data.Events)
                    continue
                end
                evt       = subj_data.Events.Events;
                % find feature average by condition
                conds     = evt.classLabels;
                if size(conds) < 5
                    continue
                end
                n_conds   = max(conds);
                subj_mean = zeros(1,n_conds);
                for i_c = 1:n_conds
                    matching_trials = conds == i_c;
                    matching_trials = logical(sum(matching_trials,2));
                    evt_feats = evt.(feats{1}{i_f});
                    evt_feats = evt_feats(matching_trials);
                    % find condition mean
                    try     
                        cond_mean = nanmean(evt_feats,1);
                    catch
                        cond_mean = 0;
                    end
                    if isequal(feats{1}{i_f},'duration')
                        % Note: convert s->ms
                        cond_mean = cond_mean * 1000;
                    end
                    try
                        subj_mean(i_c) = cond_mean;
                    catch
                        subj_mean(i_c) = 0;
                    end
                end
            end
            try
                grp_feats = [grp_feats; subj_mean];
            catch
                continue
            end
        end

        % skip if no events occurred
        if isequal(feats{1}{i_f},'eventnumber') && nnz(grp_feats)==0
            close 
            break
        end

        % add group features to cells
        feats_by_group{i_g,i_f} = grp_feats;
    end 
end

%% plot figures
figs = figure('Position',[0 0 500 750]);
sgtitle(plot_title)
for i_f = 1:n_feats
    subplot(n_feats,1,i_f)
    hold on

    % plot feature means
    plot_colors = {'#4444c9','#32a852','#db4214'};
    for i_gt = 1:length(groups_to_plot)
        i_g = groups_to_plot(i_gt);
        grp_feats  = feats_by_group{i_g,i_f};
%         plot(nanmean(grp_feats,1),'o-','Color',plot_colors{i_g})
        errorbar(nanmean(grp_feats,1),var(grp_feats,1),'Color',plot_colors{i_g})
%         n_grp_subs = size(grp_feats,1);
%         for i_gs = 1:n_grp_subs
%             plot(grp_feats(i_gs,:),'o-','Color',plot_colors{i_g})
%         end
    end
    xlabel('Condition','FontSize',5)
    ylabel(feats{2}{i_f})
    legend(groups{groups_to_plot})
    
%     grp_f_t = cellfun(@transpose,grp_feats,'UniformOutput',false);
%     all_feats = [grp_f_t{groups_to_plot}];
%     all_sp2   = [data{2,groups_to_plot}];
%     [r,p] = corrcoef(all_feats,all_sp2);
%     str=sprintf('r=%1.2f,p=%1.2f',r(1,2),p(1,2));
%     T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
%     set(T,'fontsize',10,'verticalalignment','top','horizontalalignment','left');
    
    hold off
end
% shrink legend and move off of plots
% [h,~] = legend(groups{groups_to_plot});
% set(h,'position',[0.03 .1 0.04 0.1]);

end