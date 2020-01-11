function [CEA] = CEA_data(CEA_file,num_rows)

    data = zeros(num_rows, 4);  %Cols. = OF, T, MW, gamma respectively.
    [~, data(:, 1), data(:, 2), data(:, 3), ~, data(:, 4)] = ...
    import_cea_data(CEA_file, 2, num_rows + 1);   %Get CEA data pts.
    CEA = zeros(size(unique(data(:, 1)), 1), 4);    %Actual interp. pts matrix.
    CEA(:, 1) = unique(data(:, 1));                 %Get unique OF values.

%Loop to clean up CEA data (average for each OF ratio):
%NOTE: Don't need this loop if you use interp2 with OF and Pcc.
    for i = 1:size(CEA, 1)
        CEA(i, :) = mean(data(data(:, 1) == CEA(i, 1), :));
    end
    
end