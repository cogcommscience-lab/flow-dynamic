
% Community detection code from Rubinov & Sporns, 2010 (community_louvain.m)
% Flexibility code from Network Community toolbox (flexibility.m)
% http://commdetect.weebly.com/
% Other code thanks to Sizemore & Bassett, 2018 


% Community detection for thresholded and binarized multilayer adjacency matrix
% Uses the Potts-model Hamiltonian algorithm (for binary networks)

for i = 1:4
    comm_bore_ml{i} = community_louvain(bin_bore_ml{1,i},1,[],'potts');
end

for i = 1:4
    comm_flow_ml{i} = community_louvain(bin_flow_ml{1,i},1,[],'potts');
end

for i = 1:4
    comm_frust_ml{i} = community_louvain(bin_frust_ml{1,i},1,[],'potts');
end

% Convert cells to matrix

comm_bore_ml_mat = cell2mat(comm_bore_ml);

comm_flow_ml_mat = cell2mat(comm_flow_ml);

comm_frust_ml_mat = cell2mat(comm_frust_ml);


% Flexibility for thresholded and binarized multilayer adjacency matrix

flex_bore_ml = flexibility(comm_bore_ml_mat,'temp');

flex_flow_ml = flexibility(comm_flow_ml_mat,'temp');

flex_frust_ml = flexibility(comm_frust_ml_mat,'temp');

