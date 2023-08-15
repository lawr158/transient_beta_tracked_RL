function permuted_list = permute_cell_lists(a,b)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n_a = length(a);
n_b = length(b);

permuted_list = cell(1,n_a*n_b);
idx = 1;
for i_a = 1:n_a
    for i_b = 1:n_b
        permuted_list{idx} = [a{i_a} ' ' b{i_b}];
        idx = idx + 1;
    end
end

end

