% This code gets the network flexibility for cog control and reward
% networks as identified by the Power 264 atlas

% IMPORTANT: for node_key.mat each column is a network, starting at 1
% (Sensory/somatomotor Hand)and ending at 14 (undefined, -1 in Power Key)


% Make a 1x4 cell with the node key

load('node_key.mat')

key_cell = cell(1,4);

for i = 1:4
    key_cell{1,i} = node_key;
end

% Get the order of the nodes in the Power Atlas

for i = 1:4
    [~,key_cell_sort{1,i}]=sort(key_cell{1,i});
end

% Sort the multilayer adjacency matrices by the Power Atlas

for i = 1:4
    bin_bore_ml_sort{1,i} = bin_bore_ml{1,i}(key_cell_sort{1,i},key_cell_sort{1,i});
end

for i = 1:4
    bin_flow_ml_sort{1,i} = bin_flow_ml{1,i}(key_cell_sort{1,i},key_cell_sort{1,i});
end

for i = 1:4
    bin_frust_ml_sort{1,i} = bin_frust_ml{1,i}(key_cell_sort{1,i},key_cell_sort{1,i});
end

% Sort the multilayer community assignment by the Power Atlas

for i = 1:4
    comm_bore_ml_sort{1,i} = comm_bore_ml{1,i}(key_cell_sort{1,i});
end

for i = 1:4
    comm_flow_ml_sort{1,i} = comm_flow_ml{1,i}(key_cell_sort{1,i});
end

for i = 1:4
    comm_frust_ml_sort{1,i} = comm_frust_ml{1,i}(key_cell_sort{1,i});
end

% Make cells multidimensional array

comm_bore_ml_sort_mat = cat(3, comm_bore_ml_sort{:});
comm_flow_ml_sort_mat = cat(3, comm_flow_ml_sort{:});
comm_frust_ml_sort_mat = cat(3, comm_frust_ml_sort{:});

% Extract FPCN Community Assignment

comm_bore_fpcn_win1 = comm_bore_ml_sort_mat(157:181,1,1);
comm_bore_fpcn_win2 = comm_bore_ml_sort_mat(157:181,1,2);
comm_bore_fpcn_win3 = comm_bore_ml_sort_mat(157:181,1,3);
comm_bore_fpcn_win4 = comm_bore_ml_sort_mat(157:181,1,4);
comm_bore_fpcn_ml_mat = cat(2, comm_bore_fpcn_win1, comm_bore_fpcn_win2, comm_bore_fpcn_win3, comm_bore_fpcn_win4);

comm_flow_fpcn_win1 = comm_flow_ml_sort_mat(157:181,1,1);
comm_flow_fpcn_win2 = comm_flow_ml_sort_mat(157:181,1,2);
comm_flow_fpcn_win3 = comm_flow_ml_sort_mat(157:181,1,3);
comm_flow_fpcn_win4 = comm_flow_ml_sort_mat(157:181,1,4);
comm_flow_fpcn_ml_mat = cat(2, comm_flow_fpcn_win1, comm_flow_fpcn_win2, comm_flow_fpcn_win3, comm_flow_fpcn_win4);

comm_frust_fpcn_win1 = comm_frust_ml_sort_mat(157:181,1,1);
comm_frust_fpcn_win2 = comm_frust_ml_sort_mat(157:181,1,2);
comm_frust_fpcn_win3 = comm_frust_ml_sort_mat(157:181,1,3);
comm_frust_fpcn_win4 = comm_frust_ml_sort_mat(157:181,1,4);
comm_frust_fpcn_ml_mat = cat(2, comm_frust_fpcn_win1, comm_frust_fpcn_win2, comm_frust_fpcn_win3, comm_frust_fpcn_win4);

% Compute Flexibility Cog Control (fronto-parietal control network)

flex_bore_fpcn_ml = flexibility(comm_bore_fpcn_ml_mat,'temp');
flex_flow_fpcn_ml = flexibility(comm_flow_fpcn_ml_mat,'temp');
flex_frust_fpcn_ml = flexibility(comm_frust_fpcn_ml_mat,'temp');

% Extract Subcortical Community Assignment

comm_bore_subcort_win1 = comm_bore_ml_sort_mat(200:212,1,1);
comm_bore_subcort_win2 = comm_bore_ml_sort_mat(200:212,1,2);
comm_bore_subcort_win3 = comm_bore_ml_sort_mat(200:212,1,3);
comm_bore_subcort_win4 = comm_bore_ml_sort_mat(200:212,1,4);
comm_bore_subcort_ml_mat = cat(2, comm_bore_subcort_win1, comm_bore_subcort_win2, comm_bore_subcort_win3, comm_bore_subcort_win4);

comm_flow_subcort_win1 = comm_flow_ml_sort_mat(200:212,1,1);
comm_flow_subcort_win2 = comm_flow_ml_sort_mat(200:212,1,2);
comm_flow_subcort_win3 = comm_flow_ml_sort_mat(200:212,1,3);
comm_flow_subcort_win4 = comm_flow_ml_sort_mat(200:212,1,4);
comm_flow_subcort_ml_mat = cat(2, comm_flow_subcort_win1, comm_flow_subcort_win2, comm_flow_subcort_win3, comm_flow_subcort_win4);

comm_frust_subcort_win1 = comm_frust_ml_sort_mat(200:212,1,1);
comm_frust_subcort_win2 = comm_frust_ml_sort_mat(200:212,1,2);
comm_frust_subcort_win3 = comm_frust_ml_sort_mat(200:212,1,3);
comm_frust_subcort_win4 = comm_frust_ml_sort_mat(200:212,1,4);
comm_frust_subcort_ml_mat = cat(2, comm_frust_subcort_win1, comm_frust_subcort_win2, comm_frust_subcort_win3, comm_frust_subcort_win4);

% Compute Flexibility Cog Control (subcortical network)

flex_bore_subcort_ml = flexibility(comm_bore_subcort_ml_mat,'temp');
flex_flow_subcort_ml = flexibility(comm_flow_subcort_ml_mat,'temp');
flex_frust_subcort_ml = flexibility(comm_frust_subcort_ml_mat,'temp');

% Get multilayer connectivity for fronto-parietal control network subgraph

for i = 1:4
    bin_bore_fpcn_ml{1,i} = bin_bore_ml_sort{1,i}(157:181,157:181);
end

for i = 1:4
    bin_flow_fpcn_ml{1,i} = bin_flow_ml_sort{1,i}(157:181,157:181);
end

for i = 1:4
    bin_frust_fpcn_ml{1,i} = bin_frust_ml_sort{1,i}(157:181,157:181);
end

% Make cells multidimensional array

bin_bore_fpcn_ml_mat = cat(3, bin_bore_fpcn_ml{:});
bin_flow_fpcn_ml_mat = cat(3, bin_flow_fpcn_ml{:});
bin_frust_fpcn_ml_mat = cat(3, bin_frust_fpcn_ml{:});

% Make contact sequence for fronto-parietal control network 

bin_bore_fpcn_ml_ct = arrayToContactSeq(bin_bore_fpcn_ml_mat,0);
bin_flow_fpcn_ml_ct = arrayToContactSeq(bin_flow_fpcn_ml_mat,0);
bin_frust_fpcn_ml_ct = arrayToContactSeq(bin_frust_fpcn_ml_mat,0);

% Get multilayer connectivity for subcortical network subgraph

for i = 1:4
    bin_bore_subcort_ml{1,i} = bin_bore_ml_sort{1,i}(200:212,200:212);
end

for i = 1:4
    bin_flow_subcort_ml{1,i} = bin_flow_ml_sort{1,i}(200:212,200:212);
end

for i = 1:4
    bin_frust_subcort_ml{1,i} = bin_frust_ml_sort{1,i}(200:212,200:212);
end

bin_bore_subcort_full = bin_flow_full(200:212,200:212);

% Make cells multidimensional array

bin_bore_subcort_ml_mat = cat(3, bin_bore_subcort_ml{:});
bin_flow_subcort_ml_mat = cat(3, bin_flow_subcort_ml{:});
bin_frust_subcort_ml_mat = cat(3, bin_frust_subcort_ml{:});

% Make contact sequence for subcortical network 

bin_bore_subcort_ml_ct = arrayToContactSeq(bin_bore_subcort_ml_mat,0);
bin_flow_subcort_ml_ct = arrayToContactSeq(bin_flow_subcort_ml_mat,0);
bin_frust_subcort_ml_ct = arrayToContactSeq(bin_frust_subcort_ml_mat,0);

% Plot, for example, the multilayer network for subcortical subgraphs

figure()
plotDNarc(bin_bore_subcort_ml_ct);
title('Subcortical Multilayer Dynamics During Boredom')

figure()
plotDNarc(bin_flow_subcort_ml_ct);
title('Subcortical Multilayer Dynamics During Flow')

figure()
plotDNarc(bin_frust_subcort_ml_ct);
title('Subcortical Multilayer Dynamics During Frustration')


% Make some ring plots. NB: these look weighted but should look binary
% TODO: try https://github.com/paul-kassebaum-mathworks/circularGraph

color_subcort = [0.2, 0, 0];

figure();
for i = 1:4
    subplot(2, 2, i);
    plotArcNetwork(bin_bore_subcort_ml_mat(:,:,i),color_subcort,color_subcort);
    title(['Boredom Window ',sprintf('%d',i)])
end
    
figure();
for i = 1:4
    subplot(2, 2, i);
    plotArcNetwork(bin_flow_subcort_ml_mat(:,:,i),color_subcort,color_subcort);
    title(['Flow Window ',sprintf('%d',i)])
end

figure();
for i = 1:4
    subplot(2, 2, i);
    plotArcNetwork(bin_frust_subcort_ml_mat(:,:,i),color_subcort,color_subcort);
    title(['Frustration Window ',sprintf('%d',i)])
end

color_fpcn = [1, 1, 0];

figure();
for i = 1:4
    subplot(2, 2, i);
    plotArcNetwork(bin_bore_fpcn_ml_mat(:,:,i),color_fpcn,color_fpcn);
    title(['Boredom Window ',sprintf('%d',i)])
end
    
figure();
for i = 1:4
    subplot(2, 2, i);
    plotArcNetwork(bin_flow_fpcn_ml_mat(:,:,i),color_fpcn,color_fpcn);
    title(['Flow Window ',sprintf('%d',i)])
end

figure();
for i = 1:4
    subplot(2, 2, i);
    plotArcNetwork(bin_frust_fpcn_ml_mat(:,:,i),color_fpcn,color_fpcn);
    title(['Frustration Window ',sprintf('%d',i)])
end