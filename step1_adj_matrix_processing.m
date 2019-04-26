% Process the non-windowed adjacency matrices 


% Threshold each non-windowed adjacency matrix at a given density

thresh_bore_full = threshold_proportional(data_bore_full,.3);

thresh_flow_full = threshold_proportional(data_flow_full,.3);

thresh_frust_full = threshold_proportional(data_frust_full,.3);

% Make each thresholded non-windowed adjacency matrix binary

bin_bore_full = weight_conversion(thresh_bore_full,'binarize');

bin_flow_full = weight_conversion(thresh_flow_full,'binarize');

bin_frust_full = weight_conversion(thresh_frust_full,'binarize');

% Sanity check. Is the non-windowed adjacency matrix desnity correct?
% Should match threshold specified above

den_bore_full = density_und(bin_bore_full);

den_flow_full = density_und(bin_flow_full);

den_frust_full = density_und(bin_frust_full);



% Process the windowed adjacency matrices 

% Combine windowed adjacency matrices to cells. This is your multilayer
% adjacency matrix for each condition

data_bore_ml = {data_bore_win1,data_bore_win2,data_bore_win3,data_bore_win4};

data_flow_ml = {data_flow_win1,data_flow_win2,data_flow_win3,data_flow_win4};

data_frust_ml = {data_frust_win1,data_frust_win2,data_frust_win3,data_frust_win4};

% Threshold each multilayer adjacency matrix at a given density

for i = 1:4
    thresh_bore_ml{i} = threshold_proportional(data_bore_ml{1,i},.3);
end

for i = 1:4
    thresh_flow_ml{i} = threshold_proportional(data_flow_ml{1,i},.3);
end

for i = 1:4
    thresh_frust_ml{i} = threshold_proportional(data_frust_ml{1,i},.3);
end

% Make each thresholded multilayer adjacency matrix binary

for i = 1:4
    bin_bore_ml{i} = weight_conversion(thresh_bore_ml{1,i},'binarize');
end

for i = 1:4
    bin_flow_ml{i} = weight_conversion(thresh_flow_ml{1,i},'binarize');
end

for i = 1:4
    bin_frust_ml{i} = weight_conversion(thresh_frust_ml{1,i},'binarize');
end

% Sanity check. Is the multilayer adjacency matrix desnity correct?
% Should match threshold specified above

for i = 1:4
    den_bore_ml{i} = density_und(bin_bore_ml{1,i});
end

for i = 1:4
    den_flow_ml{i} = density_und(bin_flow_ml{1,i});
end

for i = 1:4
    den_frust_ml{i} = density_und(bin_frust_ml{1,i});
end

% Make sanity check results easier to read

den_bore_ml_mat = cell2mat(den_bore_ml);

den_flow_ml_mat = cell2mat(den_flow_ml);

den_frust_ml_mat = cell2mat(den_frust_ml);


% Make thresholded and binarized multilayer adjacency matrix
% multidimensional array

array_bore_ml = cat(3, bin_bore_ml{:});
array_flow_ml = cat(3, bin_flow_ml{:});
array_frust_ml = cat(3, bin_frust_ml{:});
