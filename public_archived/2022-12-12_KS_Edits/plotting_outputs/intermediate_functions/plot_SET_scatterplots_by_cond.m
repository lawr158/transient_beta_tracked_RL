function figs = plot_SET_scatterplots_by_cond(data, groups,         ...
                                    feats, plot_type, equal_Ns,     ...
                                    evt_per_sec, plot_title)
%% plot_SET_scatterplots: Create scatterplots from spectral events data

% inputs
%   - data:      1xn cell, for n groups, containing spec_events structs
%   - groups:    1xn cell containing group names
%   - feats:     2xm cell, row 1 with m fields to extract from spec_events
%                structs; row 2 with m names describing each field
%   - plot_type: 'indiv'|'spread'|'dist'
%                   indiv  - plot individual means using plotSpread
%                   spread - plot group features with plotSpread
%                   dist   - split distribution plots of first 2 groups
%   - equal_Ns:  0 or 1 to indicate whether to equalize sample sizes by
%                randomly selecting elements from larger datasets
%
% outputs
%   - figs:     matlab figures with scatterplots for each feature
%   - means:    table containing mean and median values of each feature   
%   - p_vals:   table containing significance of group differences

%% get data
% get number of groups and features
n_grps  = length(data);
n_feats = length(feats{1});

% iterate through groups in data
feat_mean_by_subj = cell(n_grps,n_feats);
feats_by_group    = cell(n_grps,n_feats);
for i_g = 1:n_grps
    specEv_struct = data{i_g};
    n_subj = size(specEv_struct,2);
    % iterate through features
    for i_f = 1:n_feats
        n_cond = max(specEv_struct(1).TrialSummary.TrialSummary.classLabels);
        grp_feats = cell(1,n_cond);
        subj_mean = cell(n_subj,n_cond);
        % iterate through subjects
        for i_s = 1:n_subj
            % add subject features and feature means to arrays
            if isequal(feats{1}{i_f},'eventnumber')  % trial-level features
                if isempty(specEv_struct(i_s).TrialSummary)
                    continue
                end
                trial = specEv_struct(i_s).TrialSummary.TrialSummary;
                conds = trial.classLabels;
                n_conds = max(conds);
                for i_c = 1:n_cond
                    cond_inds = conds == i_c;
                    cond_feats = trial.(feats{1}{i_f})(cond_inds);
                    grp_feats{i_c} = [grp_feats{i_c}; cond_feats];
                    subj_mean{i_s,i_c} = nanmean(cond_feats,1);
                end
            else  % event-level features
                if isempty(specEv_struct(i_s).Events)
                    continue
                end
                evt = specEv_struct(i_s).Events.Events;
                
                conds = evt.classLabels;
                n_conds = max(conds);
                if isequal(feats{1}{i_f},'duration')
                    % Note: convert s->ms
                    for i_c = 1:n_cond
                        cond_inds = conds == i_c;
                        cond_feats = evt.(feats{1}{i_f})(cond_inds)*1000;
                        grp_feats{i_c} = [grp_feats{i_c}; cond_feats];
                        subj_mean{i_s,i_c} = nanmean(cond_feats,1);
                    end
                else
                    for i_c = 1:n_cond
                        cond_inds = conds == i_c;
                        cond_feats = evt.(feats{1}{i_f})(cond_inds);
                        grp_feats{i_c} = [grp_feats{i_c}; cond_feats];
                        subj_mean{i_s,i_c} = nanmean(cond_feats,1);
                    end
                end
            end
        end

        % skip if no events occurred
        if isequal(feats{1}{i_f},'eventnumber') && min(cellfun(@nnz,grp_feats))==0 % nnz(grp_feats)==0
            close 
            break
        end

        % add group features to cells
        feat_mean_by_subj{i_g,i_f} = subj_mean;
        feats_by_group{i_g,i_f}    = grp_feats';
    end 
end

%% plot figures
figs = {};
for i_f = 1:n_feats
    figs{i_f} = figure('Position',[0 0 800 250]);
    for i_c = 1:n_cond
        ax(i_c) = subplot(1,n_cond,i_c);
        hold on

        switch plot_type
            case 'Subject means'
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
                    all_feats = feat_mean_by_subj(:,i_f);
                    for i_g = 1:n_grps
                        n_grp_feats{i_g} = [all_feats{i_g}{:,i_c}];
                    end
                    n_grp_feats = cellfun(@transpose,n_grp_feats,'UniformOutput',0);
                end
                if evt_per_sec ~= 0 && i_f == 1
                    n_grp_feats = cellfun(@(x) x/evt_per_sec,n_grp_feats,'un',0);
                end

                % plot feature means
                plotSpread( n_grp_feats, 'xNames', groups, 'showMM', 3,...
                            'distributionColors', {'#4444c9','#32a852'} ) %,'#db4214'} )
                title(['Cond. ' int2str(i_c)])

        end
        ylabel(feats{2}{i_f},'FontSize',5)
        hold off
    end
    linkaxes(ax)
    sgtitle(plot_title,'FontSize',10,'FontWeight','bold')
end
% shrink legend and move off of plots
if strcmp(plot_type,'feature dist')
    [h,~] = legend(groups);
    set(h,'position',[0 0 0.07 0.1]);
end

end