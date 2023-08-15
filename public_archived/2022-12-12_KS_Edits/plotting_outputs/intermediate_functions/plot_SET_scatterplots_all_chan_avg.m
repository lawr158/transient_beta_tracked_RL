function [figs, means, p_vals] = plot_SET_scatterplots_all_chan_avg(data, groups,    ...
                                    feats, plot_type, equal_Ns,evt_per_sec)
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
        grp_feats = [];
        subj_mean = zeros(1,n_subj);
        % iterate through subjects
        for i_s = 1:n_subj
            subj_feats = [];
            % add subject features and feature means to arrays
            if isequal(feats{1}{i_f},'eventnumber')  % trial-level features
                trials = specEv_struct(i_s).TrialSummary;
                for i_t = 1:length(trials)
                    trial = trials(i_t).TrialSummary;
                    subj_feats = [subj_feats; trial.(feats{1}{i_f})];
                end
                subj_mean(i_s) = nanmean(subj_feats,1);
%                 subj_mean(i_s) = var(subj_feats,1);
                grp_feats = [grp_feats; subj_feats];
                
            else  % event-level features
                events = specEv_struct(i_s).Events;
                for i_e = 1:length(events)
                    event = events(i_e).Events;
                    subj_feats = [subj_feats; event.(feats{1}{i_f})];
                end
                if isequal(feats{1}{i_f},'duration')
                    % Note: convert s->ms
                    grp_feats = [grp_feats; subj_feats*1000];  
                    subj_mean(i_s) = nanmean(subj_feats,1) * 1000;
%                     subj_mean(i_s) = var(subj_feats,1) * 1000;
                else
                    grp_feats = [grp_feats; subj_feats];  
                    subj_mean(i_s) = nanmean(subj_feats,1);
%                     subj_mean(i_s) = var(subj_feats,1);
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
figs = figure('Position',[0 0 450 200]);
for i_f = 1:n_feats
    subplot(1,n_feats,i_f)
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
                grp_feats = feat_mean_by_subj(:,i_f);
                grp_feats =cellfun(@transpose,grp_feats,'UniformOutput',0);
            end
            if evt_per_sec ~= 0 && i_f == 1
                grp_feats = cellfun(@(x) x/evt_per_sec,grp_feats,'un',0);
            end
            % plot feature means
            plotSpread( grp_feats, 'xNames', groups,                    ...
                    'distributionColors', {'#4444c9','#32a852'} )
%             if i_f==1
%                 title('Beta event number')
%                 ylabel('Events per sec')
%             end
        case 'Group spread'
            % find smallest N to make sample sizes equal, if selected
            if equal_Ns
                Ns = cellfun(@length,feats_by_group(:,i_f),        ...
                        'UniformOutput',false);
                Ns = [Ns{:}];
                min_N = min(Ns);
            
                % select min_N random entries from subject feature means
                grp_feats = cell(1,n_grps);
                for i_g = 1:n_grps
                    grp_means = feats_by_group{i_g,i_f}';
                    perm = randperm(Ns(i_g));
                    grp_feats{i_g} = grp_means(perm(1:min_N));
                end
            else
                grp_feats = feats_by_group(:,i_f);
                grp_feats =cellfun(@transpose,grp_feats,'UniformOutput',0);
            end
            if evt_per_sec ~= 0 && i_f == 1        
                grp_feats = cellfun(@(x) x/evt_per_sec,grp_feats,'un',0);
            end
            
            % plot feature means
            plotSpread( grp_feats, 'xNames', groups,                    ...
                        'distributionColors', {'b','r'} )  
            
        case 'Group distribution'
            % find smallest N to make sample sizes equal, if selected
            if equal_Ns
                Ns = cellfun(@length,feats_by_group(:,i_f),        ...
                        'UniformOutput',false);
                Ns = [Ns{:}];
                min_N = min(Ns);
            
                % select min_N random entries from subject feature means
                grp_feats = cell(1,n_grps);
                for i_g = 1:n_grps
                    grp_means = feats_by_group{i_g,i_f}';
                    perm = randperm(Ns(i_g));
                    grp_feats{i_g} = grp_means(perm(1:min_N));
                end
            else
                grp_feats = feats_by_group(:,i_f);
                grp_feats =cellfun(@transpose,grp_feats,'UniformOutput',0);
            end
            
            if evt_per_sec ~= 0 && i_f == 1
                grp_feats = cellfun(@(x) x/evt_per_sec,grp_feats,'un',0);
            end
            % plot split distribution plots of first 2 groups
            colors = {'b','r','g','o'};
            h_or   = {'right','left'};
            for i_g = 1:n_grps
                h_or_i = mod(i_g,2)+1;
                distributionPlot(grp_feats{i_g},'histOri',h_or{h_or_i}, ...
                    'color',colors{i_g},'widthDiv',[n_grps i_g],        ...
                    'showMM',0,'distWidth',1.5,'histOpt',1)
            end 
            xticklabels('')
    end
    xlabel(feats{2}{i_f},'FontSize',5)
    
    hold off
end
% shrink legend and move off of plots
if strcmp(plot_type,'feature dist')
    [h,~] = legend(groups);
    set(h,'position',[0 0 0.07 0.1]);
end

%% get statistics
% run independent samples t-test on events/trials between groups
for i_f = 1:size(feats_by_group,2)
%     [~,p{1,i_f}] = ttest2(feats_by_group{1,i_f},feats_by_group{2,i_f});
%     [~,p{2,i_f}] = ttest2(feat_mean_by_subj{1,i_f},feat_mean_by_subj{2,i_f});
    x_by_grp = padcat(feats_by_group{1,i_f}',feats_by_group{2,i_f}');
    x_subj_mean = padcat(feat_mean_by_subj{1,i_f}',feat_mean_by_subj{2,i_f}');
    p{1,i_f} = kruskalwallis(x_by_grp);
    p{2,i_f} = kruskalwallis(x_subj_mean);
end

% calculate mean values
var_names = cell(1,size(feats_by_group,2)*3);
for i_g = 1:size(feats_by_group,1)
    for i_f = 1:size(feats_by_group,2)
        means{i_g,3*i_f-2} = nanmean(feats_by_group{i_g,i_f});
        means{i_g,3*i_f-1} = nanmedian(feats_by_group{i_g,i_f});
        means{i_g,3*i_f}   = std(feats_by_group{i_g,i_f});
        var_names{3*i_f-2} = [feats{1}{i_f} ' mean'];
        var_names{3*i_f-1} = [feats{1}{i_f} ' median'];
        var_names{3*i_f}   = [feats{1}{i_f} ' stdev'];
    end
end

means = cell2table(means,'VariableNames',var_names,'RowNames',groups);
p_vals = cell2table(p,'VariableNames',feats{1},'RowNames',{'p, by features','p, by subj means'});

end