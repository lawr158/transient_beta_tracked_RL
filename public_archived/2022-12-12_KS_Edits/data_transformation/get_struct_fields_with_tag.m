function matches = get_struct_fields_with_tag( S, tag )
%% GET_STRUCT_FIELDS_WITH_TAG: get struct fields containing input tag
% Inputs:
%  - S   - struct to search for fields
%  - tag - char, tag to search S fieldnames
% Outputs:
%  - matches - cell, array with S fields containing tag

matches = {};
tag_len = length(tag)+2;
all_fields = fieldnames(S);
% iterate through struct fields
for i = 1:length(all_fields)
    field = all_fields{i,1};
    if length(field) <= tag_len, continue; end
    field_tag = field(end-tag_len:end);
    if strcmp(field_tag,tag)
        matches{end+1,1} = field;
    end
end

end