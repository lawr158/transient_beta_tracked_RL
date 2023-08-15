function sorted_data = sort_data_by_group_SET( data,group_map,group_names )
%% sort_data_by_group_SET: Sort spectral events toolbox data into groups

% inputs
%   - data:         1x1 struct, with N fieldnames for N participants
%   - group_map:    1xN double, group_map(i) = assignment for subject i
%                   OR table with id and group columns
%   - group_names:  1xG cell, for G groups containing group names
%
% ouptuts
%   - sorted:       1xM cell, of structs for each group of subjects

%%
if isa(group_map,'double')
    % sort data based on group_map of type double
    n_groups = max(group_map);
    sorted_data   = cell(1,n_groups);
    for i_g = 1:n_groups
        group_idxs = find(group_map == i_g);
        for i_gi = 1:length(group_idxs)
            sorted_data{i_g}(end+1) = data(group_idxs(i_gi));
        end
    end
else
    % sort data based on group_map in table format
    
    ids      = fieldnames(data);
    
    n_subjs  = length(ids);
    n_groups = max(cell2mat(table2cell(group_map(:,2))));%max(group_map.group);
    % create container structure for sorted data
    for i_g = 1:n_groups
        grp_name = group_names{i_g};
        sorted_data.(grp_name) = [];
    end
    for i_s = 1:n_subjs
        % find participant in group assignment table
        try
        subj_id = ids{i_s}(2:9);
        catch
            2;
        end
        subj_struct_id = ['s' subj_id];
        s_idx   = strcmp(group_map.id,subj_id);
        if ~s_idx
            continue  % skip if no matching records
        end
        % get participant group id and name
        grp_id   = cell2mat(group_map.group(s_idx));
        grp_name = group_names{grp_id};
        
        % place subject data in appropriate group
        sorted_data.(grp_name).(subj_struct_id) = data.(subj_struct_id);
    end
end