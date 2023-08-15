function [chan_mean,chan_var] = avg_chan_topo(data_mean,data_std,feats)

%% average event number per channel into one number per channel to allow for topoplotting
n_chans = size(data_mean,2);
groups  = fieldnames(data_mean);
n_grps  = length(groups);
n_feats = length(feats); 
%feat = 'eventnumber' ; 


if n_chans > 1
    for i_g = 1:n_grps
        group = groups{i_g};
        %chan_mean.(group).id = data_mean(1).(group).id;
        %chan_var.(group).id  =data_std(1).(group).id;
        for i_f = 1:n_feats
            feat = feats{i_f} ; 
            for i_c = 1:n_chans
            % get feature means/sds and average across channels
            group_means = [data_mean(i_c).(group)];
            group_means = mean([group_means.(feat)],'omitnan');
            group_vars  = [data_std(i_c).(group)];
            group_vars  = mean([group_vars.(feat)],'omitnan');
            % save to new struct
            chan_mean(i_c).(group).(feat) = group_means;
            chan_var(i_c).(group).(feat)  = group_vars;
            end %end chan loop
        end    % feature loop
    end    % group loop

    % update output struct with channel averages
   % s_mean = chan_mean;
   % s_var  = chan_var;
end    % number of channels IF
end