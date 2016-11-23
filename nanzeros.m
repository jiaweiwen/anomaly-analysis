function [ input_matrix ] = nanzeros(input_matrix)
% Get all the nanzero index
    nanzero_idx = find(input_matrix ~= 0);
    input_matrix = input_matrix(nanzero_idx);
end

