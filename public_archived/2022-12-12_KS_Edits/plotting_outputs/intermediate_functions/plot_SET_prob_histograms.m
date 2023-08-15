function [figs, mean_vals, p_vals] = plot_SET_prob_histograms(data_by_grp,groups)
% Event feature probability histograms (see Figure 5 in Shin et al. eLife 2017)
features = {'eventnumber','maximapowerFOM','duration','Fspan'}; %Fields within specEv_struct
feature_names = {'event number','event power (FOM)','event duration (ms)','event F-span (Hz)'}; %Full names describing each field
h_bins{1} = 0:1:15;
h_bins{2} = 5:0.5:20;
h_bins{3} = 50:20:350;
h_bins{4} = 4:1:25;


n_grps = size(groups,2);
fp_by_group = cell(n_grps,numel(features));
fm_by_grp_by_indiv = cell(n_grps,numel(features));
for i_g = 1:n_grps
    specEv_struct = data_by_grp{i_g};
    figure
    numSubj = size(specEv_struct,2);
    for feat_i=1:numel(features)
        feature_agg = [];
        f_indiv     = [];
        for subj_i=1:numSubj
            % Feature-specific considerations
            if isequal(features{feat_i},'eventnumber')
                feature_agg = [feature_agg; specEv_struct(subj_i).TrialSummary.TrialSummary.(features{feat_i})];
                f_indiv     = [f_indiv; nanmean(specEv_struct(subj_i).TrialSummary.TrialSummary.(features{feat_i}),1)];
                
            else
                if isequal(features{feat_i},'duration')
                    feature_agg = [feature_agg; specEv_struct(subj_i).Events.Events.(features{feat_i}) * 1000]; %Note: convert from s->ms
                    f_indiv     = [f_indiv; nanmean(specEv_struct(subj_i).Events.Events.(features{feat_i}),1) * 1000];
                else
                    feature_agg = [feature_agg; specEv_struct(subj_i).Events.Events.(features{feat_i})];
                    f_indiv     = [f_indiv; nanmean(specEv_struct(subj_i).Events.Events.(features{feat_i}),1)];
                end
            end
        end

        % Don't plot if no events occurred
        if isequal(features{feat_i},'eventnumber') && nnz(feature_agg)==0
            close 
            break
        end

        % Calculate probability of aggregate (accross subjects/sessions) and 
        % standardize bins
        [featProb_agg,bins] = histcounts(feature_agg,h_bins{feat_i},'Normalization','probability'); 

        % Correct to show left-side dropoff of histogram if applicable
        if bins(2)-(bins(2)-bins(1))/2>0
            bins = [bins(1)-(bins(2)-bins(1)),bins];
            featProb_agg = histcounts(feature_agg,bins,'Normalization','probability');
        end
        featProb_agg(isnan(featProb_agg)) = 0; %Correct for NaN values resulting from dividing by 0 counts
        fp_by_group{i_g,feat_i} = featProb_agg;
        hp_bins{feat_i} = bins;
        f_by_group{i_g,feat_i} = feature_agg;
        fm_by_grp_by_indiv{i_g,feat_i} = f_indiv;
    end 
end

% plot figures
figs = figure;
for i_f = 1:numel(features)
    subplot(numel(features),1,i_f)
    hold on
    for i_g = 1:n_grps
        histogram('BinEdges',hp_bins{i_f},'BinCounts',fp_by_group{i_g,i_f},'LineWidth',1,'DisplayStyle','stairs')
        xlabel(feature_names{i_f})
        ylabel('probability')
    end
    hold off
    legend(groups)
end

% run independent samples t-test on events/trials between groups
for i_f = 1:size(f_by_group,2)
%     min_length = min(size(f_by_group{1,i_f},1),size(f_by_group{2,i_f},1));
    [~,p{1,i_f}] = ttest2(f_by_group{1,i_f},f_by_group{2,i_f});
    [~,p{2,i_f}] = ttest2(fm_by_grp_by_indiv{1,i_f},fm_by_grp_by_indiv{2,i_f});
end

% calculate mean values
var_names = cell(1,size(f_by_group,2)*3);
for i_g = 1:size(f_by_group,1)
    for i_f = 1:size(f_by_group,2)
        mean_vals{i_g,3*i_f-2} = nanmean(f_by_group{i_g,i_f});
        mean_vals{i_g,3*i_f-1} = nanmedian(f_by_group{i_g,i_f});
        mean_vals{i_g,3*i_f}   = std(f_by_group{i_g,i_f});
        var_names{3*i_f-2} = [features{i_f} ' mean'];
        var_names{3*i_f-1} = [features{i_f} ' median'];
        var_names{3*i_f}   = [features{i_f} ' stdev'];
    end
end

mean_vals = cell2table(mean_vals,'VariableNames',var_names,'RowNames',groups);
p_vals    = cell2table(p,'VariableNames',features,'RowNames',{'p, by features','p, by subj means'});

end