function permuted_list = permute_cell_list(a)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n_a = length(a);

permuted_list = cell(1,n_a*(n_a-1));
idx = 1;
for i_a = 1:n_a
    for i_b = 1:n_a
        if i_a == i_b
            continue;
        end
        permuted_list{idx} = [a{i_a} '_x_' a{i_b}];
        idx = idx + 1;
    end
end

end

