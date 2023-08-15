function [band_pwr_table,mean_pwr] = get_band_power(grp_folder,beapp_tag,avgd,groups_to_plot,group_name,band,chans)

groups  = group_name(groups_to_plot);
n_grps  =length(groups);
n_chans = length(chans);

l = upper(band(1));
band(1) = l ;

if avgd == 1
    excel_file = 'baseline_mean_Pwr_per_Hz.csv';
    x = 'mean_Pwr_Per_Hz';
else
    excel_file = 'baseline_med_Pwr_per_Hz.csv';
    x = 'med_Pwr_Per_Hz';
end

tag = ['out_',beapp_tag];

for i_g = 1:n_grps
    grp = groups{i_g};
    excel_dir = fullfile(grp_folder{i_g},tag,excel_file);
    pwr_table.(grp) = readtable(excel_dir);
    band_pwr_table.(grp)(:,1) = pwr_table.(grp)(:,1);
    for i_c = 1 : n_chans
        ch = num2str(chans(i_c));
        tvar = ['E',ch,'_',band,'_',x];
        temp(:,i_c) = pwr_table.(grp).(tvar);
        mean_pwr.(grp) = mean(temp,2);
        band_pwr_table.(grp)(:,'channel_avg_power') = table(mean_pwr.(grp)) ;
        band_pwr_table.(grp)(:,tvar) = pwr_table.(grp)(:,tvar); 
        clear temp
    end % end chan loop
end %end grp loop
end %end function





