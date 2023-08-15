function limits = find_struct_limits(data_struct,fields,tag)
    %% find maxima across entries in a struct
    %{
    inputs
        - data_struct:  struct to analyze across fields
        - fields:       fields to analyze in struct
        - tag:          field tags to filter for
    outputs
        - limits: 1x2 double with maxima
    %}
    %%
    limits = [10 -10];
    for i_f = 1:size(fields,1)
        entry_tag = fields{i_f,1};
        if strcmp(tag,entry_tag(length(entry_tag)-length(tag)+1:end))
            entry = data_struct.(fields{i_f,1});
            % compare maxima across entries
            limits(1) = min(limits(1),min(entry));
            limits(2) = max(limits(2),max(entry));
        end
    end
end