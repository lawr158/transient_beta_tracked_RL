% Author: Gerardo Parra, Last updated: 2023-03-29
function Sa = struct2struct_array(S)
%% STRUCT2STRUCT_ARRAY: converts scalar struct to non-scalar struct array
% Inputs:
%   - S - struct to convert. assumption is that all fields contain data 
%         arrays of the same length n_row
% Outputs:
%   - Sa - non-scalar [n_row x 1] struct array
%%
% get struct fieldnames
fields = fieldnames(S); n_field = length(fields); 
% get number of rows 
n_row = length(S.(fields{1}));
% initialize output struct array
Sa = struct(fields{1},repmat({''},1,n_row));
% iterate through rows in data
for i_r = 1:n_row
    % iterate through struct fields
    for i_f = 1:n_field
        % add row/field data to output array
        field = fields{i_f};
        Sa(i_r).(field) = S.(field)(i_r);
    end
end

end