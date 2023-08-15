function [figs, means, p_vals] = plot_SET_vs_measure(data,groups,feats, ... 
                                    equal_Ns, taps_to_inc, evt_per_sec, ...
                                    plot_title, groups_to_plot, plot_var...
                                    )
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
% get number of groups and features
n_grps  = size(data,2);
n_feats = length(feats{1});

% iterate through groups in data
feat_mean_by_subj = cell(n_grps,n_feats);
feats_by_group    = cell(n_grps,n_feats);
for i_g = 1:n_grps
    specEv_struct = data{1,i_g};
    n_subj = size(specEv_struct,2);
    % iterate through features
    for i_f = 1:n_feats
        grp_feats = [];
        subj_mean = zeros(1,n_subj);
        if isempty(specEv_struct)
            continue
        end
        n_cond = max(specEv_struct(1).TrialSummary.TrialSummary.classLabels);
%         grp_feats = cell(1,n_cond);
%         subj_mean = cell(n_subj,n_cond);
        % iterate through subjects
        for i_s = 1:n_subj
            % add subject features and feature means to arrays
            if isequal(feats{1}{i_f},'eventnumber')  % trial-level features
                if isempty(specEv_struct(i_s).TrialSummary)
                    continue
                end
                
                trial = specEv_struct(i_s).TrialSummary.TrialSummary;
                if taps_to_inc <= 0
                    grp_feats = [grp_feats; trial.(feats{1}{i_f})];
                    if ~plot_var
                        subj_mean(i_s) = nanmean(trial.(feats{1}{i_f}),1);
                    else
                        subj_mean(i_s) = var(trial.(feats{1}{i_f}),1);
                    end
                else
                    % select specific trials
                    matching_trials = trial.classLabels == taps_to_inc;
                    matching_trials = logical(sum(matching_trials,2));
                    trial_feats = trial.(feats{1}{i_f});
                    grp_feats = [grp_feats; trial_feats(matching_trials)];
                    
                    if ~plot_var
                        subj_mean(i_s) = nanmean(trial_feats(matching_trials),1);
                    else
                        subj_mean(i_s) = var(trial_feats(matching_trials),1);
                    end
                end
                
            else  % event-level features
                if isempty(specEv_struct(i_s).Events)
                    continue
                end
                evt = specEv_struct(i_s).Events.Events;
                evt_feats = evt.(feats{1}{i_f});
                if taps_to_inc > 0
                    matching_trials = evt.classLabels == taps_to_inc;
                    matching_trials = logical(sum(matching_trials,2));
                    evt_feats = evt_feats(matching_trials);
                end
                conds = evt.classLabels;
                n_conds = max(conds);
                if isequal(feats{1}{i_f},'duration')
                    % Note: convert s->ms
                    grp_feats = [grp_feats; evt_feats*1000]; 
                    try
                        if ~plot_var
                            subj_mean(i_s) = nanmean(evt_feats,1) * 1000;
                        else
                            subj_mean(i_s) = var(evt_feats,1) * 1000;
                        end
                    catch
                        subj_mean(i_s) = 0;
                    end
                    
                else
                    grp_feats = [grp_feats; evt_feats];
                    try
                        if ~plot_var
                            subj_mean(i_s) = nanmean(evt_feats,1);
                        else
                            subj_mean(i_s) = var(evt_feats,1);
                        end
                    catch
                        subj_mean(i_s) = 0;
                    end
                    
                end
            end
        end

        % skip if no events occurred
        if isequal(feats{1}{i_f},'eventnumber') && nnz(grp_feats)==0
            close 
            break
        end

        % add group features to cells
        feat_mean_by_subj{i_g,i_f} = subj_mean;
        feats_by_group{i_g,i_f}    = grp_feats';
    end 
end

%% plot figures
figs = figure('Position',[0 0 650 250]);
sgtitle(plot_title)
for i_f = 1:n_feats
    subplot(1,n_feats,i_f)
    hold on
    
    % find smallest N to make sample sizes equal, if selected
    if equal_Ns
        Ns = cellfun(@length,feat_mean_by_subj(:,i_f),        ...
                'UniformOutput',false);
        Ns = [Ns{:}];
        min_N = min(Ns);

        % select min_N random entries from subject feature means
        grp_feats = cell(1,n_grps);
        for i_g = 1:n_grps
            grp_means = feat_mean_by_subj{i_g,i_f}';
            perm = randperm(Ns(i_g));
            grp_feats{i_g} = grp_means(perm(1:min_N));
        end
    else
        grp_feats = feat_mean_by_subj(:,i_f);
        grp_feats =cellfun(@transpose,grp_feats,'UniformOutput',0);
    end
    if evt_per_sec ~= 0 && i_f == 1
        grp_feats = cellfun(@(x) x/evt_per_sec,grp_feats,'un',0);
    end

    % plot feature means
    plot_colors = {'#4444c9','#32a852','#db4214'};
    for i_gt = 1:length(groups_to_plot)
        i_g = groups_to_plot(i_gt);
        scatter(grp_feats{i_g},data{2,i_g},'filled','MarkerFaceColor',plot_colors{i_g})
    end
    xlabel(feats{2}{i_f},'FontSize',5)
    ylabel('Score')
    
    grp_f_t = cellfun(@transpose,grp_feats,'UniformOutput',false);
    all_feats = [grp_f_t{groups_to_plot}];
    all_sp2   = [data{2,groups_to_plot}];
    [r,p] = corrcoef(all_feats,all_sp2);
    str=sprintf('r=%1.2f,p=%1.2f',r(1,2),p(1,2));
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    set(T,'fontsize',10,'verticalalignment','top','horizontalalignment','left');
    
    hold off
end
% shrink legend and move off of plots
[h,~] = legend(groups{groups_to_plot});
set(h,'position',[0.03 .1 0.04 0.1]);

end