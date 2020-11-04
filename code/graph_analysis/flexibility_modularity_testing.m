% Codes for inferential testing of flexibility and modularity for dynamic
% functional connectivity network of flow state vs. boredom and frustrution
% 
% function used from genLouvain (http://netwiki.amath.unc.edu/GenLouvain/GenLouvain)
% flexibility

clc;
close all;
clear;
load(['/data/time_series/output_01608_gsr_seitzman.mat']);

% load the the adjacency matrices, then store them in a 5-d array
M_bore = cat(4, bore.win_1, bore.win_2, bore.win_3, bore.win_4);
M_flow = cat(4, flow.win_1, flow.win_2, flow.win_3, flow.win_4);
M_frus = cat(4, frus.win_1, frus.win_2, frus.win_3, frus.win_4);
M = cat(5, M_bore, M_flow, M_frus);
M = permute(M, [1 5 4 2 3]);

% construct mean graph from the set of graph for all individuals
M = squeeze(mean(M, 1));
M_weight = squeeze(mean(M_weight, 1));


for i = 1:size(M,1)
    for j = 1:size(M,2)
        % we first set the diag of the matrix to 0
        M(i,j,:,:) = squeeze(M(i,j,:,:)) - diag(diag(squeeze(M(i,j,:,:))));
        % we threthold the matrix at 30%
        M(i,j,:,:) = threshold_proportional( squeeze(M(i,j,:,:)),.3);
        % we binarize the weighted matrix to binary matrix
        M(i,j,:,:) = weight_conversion( squeeze(M(i,j,:,:)),'binarize');
    end
end


%% Define node keys

node_key = readtable('/data/Label.txt', 'ReadVariableNames',false);
node_key = table2array(node_key);
unique(node_key)

Reward = find(strcmp('Reward',node_key));
Auditory = find(strcmp('Auditory',node_key));
CinguloOpercular = find(strcmp('CinguloOpercular',node_key));
DefaultMode = find(strcmp('DefaultMode',node_key));
DorsalAttention = find(strcmp('DorsalAttention',node_key));
FrontoParietal = find(strcmp('FrontoParietal',node_key));
MedialTemporalLobe = find(strcmp('MedialTemporalLobe',node_key));
ParietoMedial = find(strcmp('ParietoMedial',node_key));
Salience = find(strcmp('Salience',node_key));
SomatomotorDorsal = find(strcmp('SomatomotorDorsal',node_key));
SomatomotorLateral = find(strcmp('SomatomotorLateral',node_key));
VentralAttention = find(strcmp('VentralAttention',node_key));
Visual = find(strcmp('Visual',node_key));
unassigned = find(strcmp('unassigned',node_key));
Frontoparietalsubcortical = [FrontoParietal; Reward];



% 1st dim of M should be ses, 2nd should be win, third should be
% nodes
n_ses = size(M,1);
n_win = size(M,2);
n_roi = size(M,3);
%% community detection
gamma = 1;
omega = 1;
n_rep = 1000;

% Empty matrices to store data
modularity_mean = zeros(n_rep, n_ses); % Mean modularity
modules = zeros(n_ses, n_rep, n_roi, n_win); % Module assigment labels    


for ses = 1 : n_ses
%--- define objects ---------------------------------------------------
    A = cell(1, n_win);
    B = spalloc(n_roi * n_win, n_roi * n_win,(n_roi + n_win) * n_roi* n_win);
    twomu = 0;

    %--- null model -------------------------------------------------------
    for win = 1 : n_win
        %--- copy network  --------------
        A{win} = squeeze(M(ses, win, :, :));
        k = sum(A{win});                             % node degree
        twom = sum(k);                               % mean network degree
        twomu = twomu + twom;                        % increment
        indx = [1:n_roi] + (win-1)*n_roi;            % find indices
        B(indx,indx) = A{win} - gamma * [k'*k]/twom;  % fill B matrix
    end
    twomu = twomu + 2*omega* n_roi*(n_win-1);

    B = B + omega/2*spdiags(ones(n_roi*n_win,2),[-n_roi, n_roi], n_roi*n_win, n_roi*n_win);
    B = B + omega*spdiags(ones(n_roi*n_win,2),[-2*n_roi, 2*n_roi], n_roi*n_win, n_roi*n_win);
    %B = B + omega * spdiags(ones(n_roi * n_win,2),[-n_roi, n_roi], n_roi*n_win, n_roi*n_win);

%--- calculate multilayer modules -------------------------------------
    %Qb = 0;
    for rep = 1 : n_rep
        clc;
        [S,Q] = genlouvain(B);
        Q = Q / twomu;
        S = reshape(S, n_roi, n_win);
        
        modularity_mean(rep, ses) = Q;
        modules(ses, rep, :, :) = S;
   end
end

%% flexibility
for i = 1:size(modules,1)
    for j = 1:size(modules,2)
        nodes_flexibility(i,j,:) = flexibility(squeeze(modules(i, j, :, :))');
    end
end

% following the suggestions from Bassett et al., (2011), we compute the flexibility over repetition
% of the modularity maximization process
nodes_flexibility_mean = squeeze(mean(nodes_flexibility, 2));
nodes_mean_flexibility_mean = squeeze(mean(nodes_flexibility_mean, 2));


%% null graph construction

% we first convert the matrices to contact vectors
bin_bore_ml_sort_mat = permute(squeeze(M(1,:,:,:)), [2 3 1]);
bin_flow_ml_sort_mat = permute(squeeze(M(2,:,:,:)), [2 3 1]);
bin_frust_ml_sort_mat = permute(squeeze(M(3,:,:,:)), [2 3 1]);

% Convert empirical adjacency matrix to array
contact_bore_sort = arrayToContactSeq(bin_bore_ml_sort_mat,0);
contact_flow_sort = arrayToContactSeq(bin_flow_ml_sort_mat,0);
contact_frust_sort = arrayToContactSeq(bin_frust_ml_sort_mat,0);

% Make random permuted time graph for each contact sequence

for i=1:1000
    contact_bore_rt_null{i} = randomPermutedTimes(contact_bore_sort);
end

for i=1:1000
    contact_flow_rt_null{i} = randomPermutedTimes(contact_flow_sort);
end

for i=1:1000
    contact_frust_rt_null{i} = randomPermutedTimes(contact_frust_sort);
end

clear bin_bore_ml_sort_mat bin_flow_ml_sort_mat bin_frust_ml_sort_mat;
clear contact_bore_sort contact_flow_sort contact_frust_sort;

%% null graph processing

% Convert all null contact sequences to array

for i = 1:1000
    array_bore_rt_null{i} = networksFromContacts(contact_bore_rt_null{1,i},0);
end

for i = 1:1000
    array_flow_rt_null{i} = networksFromContacts(contact_flow_rt_null{1,i},0);
end

for i = 1:1000
    array_frust_rt_null{i} = networksFromContacts(contact_frust_rt_null{1,i},0);
end

for i = 1:1000
    M_null{i} = permute(cat(4, array_bore_rt_null{i}, array_flow_rt_null{i}, array_frust_rt_null{i}), [4, 3, 1, 2]);
end

clear array_bore_rt_null array_flow_rt_null array_frust_rt_null;

%% Community detection anf flexibility on all null models
clear nodes_flexibility_mean_null nodes_flexibility_null
n_rep = 20; % for computation efficiency the repetition for null models are set to 20

modularity_mean_null = zeros(n_rep, n_ses); % Mean modularity
modules_null = zeros(n_ses, n_rep, n_roi, n_win); % Module assigment labels 

for i = 1:1000
    %--- community detection ---------------------------------------------------
    for ses = 1 : n_ses

        A = cell(1, n_win);
        B = spalloc(n_roi * n_win, n_roi * n_win,(n_roi + n_win) * n_roi* n_win);
        twomu = 0;

        %--- null model -------------------------------------------------------
        for win = 1 : n_win
            %--- copy network  --------------
            A{win} = squeeze(M_null{i}(ses, win, :, :));
            k = sum(A{win});                             % node degree
            twom = sum(k);                               % mean network degree
            twomu = twomu + twom;                        % increment
            indx = [1:n_roi] + (win-1)*n_roi;            % find indices
            B(indx,indx) = A{win} - gamma * [k'*k]/twom;  % fill B matrix
        end
        twomu = twomu + 2*omega* n_roi*(n_win-1);

        B = B + omega/2*spdiags(ones(n_roi*n_win,2),[-n_roi, n_roi], n_roi*n_win, n_roi*n_win);
        B = B + omega*spdiags(ones(n_roi*n_win,2),[-2*n_roi, 2*n_roi], n_roi*n_win, n_roi*n_win);
        %B = B + omega * spdiags(ones(n_roi * n_win,2),[-n_roi, n_roi], n_roi*n_win, n_roi*n_win);

    %--- calculate multilayer modules -------------------------------------
        %Qb = 0;
        for rep = 1 : n_rep
            clc;
            [S,Q] = genlouvain(B);
            Q = Q / twomu;
            S = reshape(S, n_roi, n_win);

            modularity_mean_null(rep, ses) = Q;
            modules_null(ses, rep, :, :) = S;
        end
    end
    
    for j = 1:size(modules_null,1)
        for k = 1:size(modules_null,2)
            nodes_flexibility_null(j, k, :) = flexibility(squeeze(modules_null(j, k, :, :))');
        end
    end
    
    nodes_flexibility_mean_null{i} = squeeze(mean(nodes_flexibility_null, 2));
    nodes_flexibility_mean_null_fron{i} = squeeze(mean(nodes_flexibility_null(:,:,FrontoParietal), 2));
    nodes_flexibility_mean_null_sub{i} = squeeze(mean(nodes_flexibility_null(:,:,Reward), 2));
    nodes_flexibility_mean_null_fron_sub{i} = squeeze(mean(nodes_flexibility_null(:,:,Frontoparietalsubcortical), 2));
    modularity_mean_null_cell{i} = modularity_mean_null;
    
    
end


%%
% The full graph flexibility values are a cell, this makes them a regular
% array
nodes_flexibility_mean_null_array = cat(3,nodes_flexibility_mean_null{:});
nodes_flexibility_mean_null_array_fron = cat(3,nodes_flexibility_mean_null_fron{:});
nodes_flexibility_mean_null_array_sub = cat(3,nodes_flexibility_mean_null_sub{:});
nodes_flexibility_mean_null_array_fron_sub = cat(3,nodes_flexibility_mean_null_fron_sub{:});

full_flexibility_mean_null_array = squeeze(mean(nodes_flexibility_mean_null_array,2));
full_flexibility_mean_null_array_fron = squeeze(mean(nodes_flexibility_mean_null_array_fron,2));
full_flexibility_mean_null_array_sub = squeeze(mean(nodes_flexibility_mean_null_array_sub,2));
full_flexibility_mean_null_array_fron_sub = squeeze(mean(nodes_flexibility_mean_null_array_fron_sub,2));

nodes_modularity_mean_null_array = cat(3,modularity_mean_null_cell{:});
full_modularity_mean_null_array = squeeze(mean(nodes_modularity_mean_null_array,1));

%% inferential testing - preparing
flex_bore_ml = squeeze(mean(nodes_flexibility_mean(1,:), 2))
flex_flow_ml = squeeze(mean(nodes_flexibility_mean(2,:), 2))
flex_frust_ml = squeeze(mean(nodes_flexibility_mean(3,:), 2))

flex_bore_ml_rt_null_t = squeeze(full_flexibility_mean_null_array(1,:));
flex_flow_ml_rt_null_t = squeeze(full_flexibility_mean_null_array(2,:));
flex_frust_ml_rt_null_t = squeeze(full_flexibility_mean_null_array(3,:));

% And lets get some summary stats on each

flex_bore_ml_rt_null_t_mean = mean(flex_bore_ml_rt_null_t);
flex_bore_ml_rt_null_t_sd = std(flex_bore_ml_rt_null_t);
flex_bore_ml_rt_null_t_se = std(flex_bore_ml_rt_null_t)/sqrt(length(flex_bore_ml_rt_null_t));

flex_flow_ml_rt_null_t_mean = mean(flex_flow_ml_rt_null_t);
flex_flow_ml_rt_null_t_sd = std(flex_flow_ml_rt_null_t);
flex_flow_ml_rt_null_t_se = std(flex_flow_ml_rt_null_t)/sqrt(length(flex_flow_ml_rt_null_t));

flex_frust_ml_rt_null_t_mean = mean(flex_frust_ml_rt_null_t);
flex_frust_ml_rt_null_t_sd = std(flex_frust_ml_rt_null_t);
flex_frust_ml_rt_null_t_se = std(flex_frust_ml_rt_null_t)/sqrt(length(flex_frust_ml_rt_null_t));

%% inferential testing - flow less bore
% Dfferences in flexibility flow - bore

flex_ml_rt_null_flow_less_bore = (flex_flow_ml_rt_null_t - flex_bore_ml_rt_null_t)';
mean_flex_ml_rt_null_flow_less_bore = mean(flex_ml_rt_null_flow_less_bore);
vect_ones = ones(1000,1);
vect_mean_flex_ml_rt_null_flow_less_bore = vect_ones * mean_flex_ml_rt_null_flow_less_bore;
diff_mean_flex_ml_rt_null_flow_less_bore = flex_ml_rt_null_flow_less_bore - vect_mean_flex_ml_rt_null_flow_less_bore;
s_diff_mean_flex_ml_rt_null_flow_less_bore = diff_mean_flex_ml_rt_null_flow_less_bore.^2;
ss_diff_mean_flex_ml_rt_null_flow_less_bore = sum(s_diff_mean_flex_ml_rt_null_flow_less_bore);
one_over_m_less_one = 1/999;
mss_diff_mean_flex_ml_rt_null_flow_less_bore = one_over_m_less_one * ss_diff_mean_flex_ml_rt_null_flow_less_bore;
seb_mss_diff_mean_flex_ml_rt_null_flow_less_bore = sqrt(mss_diff_mean_flex_ml_rt_null_flow_less_bore);
seb_mss_diff_mean_flex_ml_rt_null_flow_less_bore_t = seb_mss_diff_mean_flex_ml_rt_null_flow_less_bore.';
num_emp_flex_flow_less_bore = flex_flow_ml - flex_bore_ml;
t_val_flex_flow_less_bore = num_emp_flex_flow_less_bore./seb_mss_diff_mean_flex_ml_rt_null_flow_less_bore_t

%% inferential testing - flow less frust

flex_ml_rt_null_flow_less_frust = (flex_flow_ml_rt_null_t - flex_frust_ml_rt_null_t)';
mean_flex_ml_rt_null_flow_less_frust = mean(flex_ml_rt_null_flow_less_frust);
vect_ones = ones(1000,1);
vect_mean_flex_ml_rt_null_flow_less_frust = vect_ones * mean_flex_ml_rt_null_flow_less_frust;
diff_mean_flex_ml_rt_null_flow_less_frust = flex_ml_rt_null_flow_less_frust - vect_mean_flex_ml_rt_null_flow_less_frust;
s_diff_mean_flex_ml_rt_null_flow_less_frust = diff_mean_flex_ml_rt_null_flow_less_frust.^2;
ss_diff_mean_flex_ml_rt_null_flow_less_frust = sum(s_diff_mean_flex_ml_rt_null_flow_less_frust);
one_over_m_less_one = 1/999;
mss_diff_mean_flex_ml_rt_null_flow_less_frust = one_over_m_less_one * ss_diff_mean_flex_ml_rt_null_flow_less_frust;
seb_mss_diff_mean_flex_ml_rt_null_flow_less_frust = sqrt(mss_diff_mean_flex_ml_rt_null_flow_less_frust);
seb_mss_diff_mean_flex_ml_rt_null_flow_less_frust_t = seb_mss_diff_mean_flex_ml_rt_null_flow_less_frust.';
num_emp_flex_flow_less_frust = flex_flow_ml - flex_frust_ml;
t_val_flex_flow_less_frust = num_emp_flex_flow_less_frust./seb_mss_diff_mean_flex_ml_rt_null_flow_less_frust_t

%% sub graph fron
flex_bore_ml_fron = squeeze(mean(nodes_flexibility_mean(1,FrontoParietal), 2))
flex_flow_ml_fron = squeeze(mean(nodes_flexibility_mean(2,FrontoParietal), 2))
flex_frust_ml_fron = squeeze(mean(nodes_flexibility_mean(3,FrontoParietal), 2))

flex_bore_ml_rt_null_t_fron = squeeze(full_flexibility_mean_null_array_fron(1,:));
flex_flow_ml_rt_null_t_fron = squeeze(full_flexibility_mean_null_array_fron(2,:));
flex_frust_ml_rt_null_t_fron = squeeze(full_flexibility_mean_null_array_fron(3,:));


flex_ml_rt_null_flow_less_bore_fron = (flex_flow_ml_rt_null_t_fron - flex_bore_ml_rt_null_t_fron)';
mean_flex_ml_rt_null_flow_less_bore_fron = mean(flex_ml_rt_null_flow_less_bore_fron);
vect_ones = ones(1000,1);
vect_mean_flex_ml_rt_null_flow_less_bore_fron = vect_ones * mean_flex_ml_rt_null_flow_less_bore_fron;
diff_mean_flex_ml_rt_null_flow_less_bore_fron = flex_ml_rt_null_flow_less_bore_fron - vect_mean_flex_ml_rt_null_flow_less_bore_fron;
s_diff_mean_flex_ml_rt_null_flow_less_bore_fron = diff_mean_flex_ml_rt_null_flow_less_bore_fron.^2;
ss_diff_mean_flex_ml_rt_null_flow_less_bore_fron = sum(s_diff_mean_flex_ml_rt_null_flow_less_bore_fron);
one_over_m_less_one = 1/999;
mss_diff_mean_flex_ml_rt_null_flow_less_bore_fron = one_over_m_less_one * ss_diff_mean_flex_ml_rt_null_flow_less_bore_fron;
seb_mss_diff_mean_flex_ml_rt_null_flow_less_bore_fron = sqrt(mss_diff_mean_flex_ml_rt_null_flow_less_bore_fron);
seb_mss_diff_mean_flex_ml_rt_null_flow_less_bore_t_fron = seb_mss_diff_mean_flex_ml_rt_null_flow_less_bore_fron.';
num_emp_flex_flow_less_bore_fron = flex_flow_ml_fron - flex_bore_ml_fron;
t_val_flex_flow_less_bore_fron = num_emp_flex_flow_less_bore_fron./seb_mss_diff_mean_flex_ml_rt_null_flow_less_bore_t_fron

flex_ml_rt_null_flow_less_frust_fron = (flex_flow_ml_rt_null_t_fron - flex_frust_ml_rt_null_t_fron)';
mean_flex_ml_rt_null_flow_less_frust_fron = mean(flex_ml_rt_null_flow_less_frust_fron);
vect_ones = ones(1000,1);
vect_mean_flex_ml_rt_null_flow_less_frust_fron = vect_ones * mean_flex_ml_rt_null_flow_less_frust_fron;
diff_mean_flex_ml_rt_null_flow_less_frust_fron = flex_ml_rt_null_flow_less_frust_fron - vect_mean_flex_ml_rt_null_flow_less_frust_fron;
s_diff_mean_flex_ml_rt_null_flow_less_frust_fron = diff_mean_flex_ml_rt_null_flow_less_frust_fron.^2;
ss_diff_mean_flex_ml_rt_null_flow_less_frust_fron = sum(s_diff_mean_flex_ml_rt_null_flow_less_frust_fron);
one_over_m_less_one = 1/999;
mss_diff_mean_flex_ml_rt_null_flow_less_frust_fron = one_over_m_less_one * ss_diff_mean_flex_ml_rt_null_flow_less_frust_fron;
seb_mss_diff_mean_flex_ml_rt_null_flow_less_frust_fron = sqrt(mss_diff_mean_flex_ml_rt_null_flow_less_frust_fron);
seb_mss_diff_mean_flex_ml_rt_null_flow_less_frust_t_fron = seb_mss_diff_mean_flex_ml_rt_null_flow_less_frust_fron.';
num_emp_flex_flow_less_frust_fron = flex_flow_ml_fron - flex_frust_ml_fron;
t_val_flex_flow_less_frust_fron = num_emp_flex_flow_less_frust_fron./seb_mss_diff_mean_flex_ml_rt_null_flow_less_frust_t_fron


%% sub graph sub

flex_bore_ml_sub = squeeze(mean(nodes_flexibility_mean(1,Reward), 2))
flex_flow_ml_sub = squeeze(mean(nodes_flexibility_mean(2,Reward), 2))
flex_frust_ml_sub = squeeze(mean(nodes_flexibility_mean(3,Reward), 2))

flex_bore_ml_rt_null_t_sub = squeeze(full_flexibility_mean_null_array_sub(1,:));
flex_flow_ml_rt_null_t_sub = squeeze(full_flexibility_mean_null_array_sub(2,:));
flex_frust_ml_rt_null_t_sub = squeeze(full_flexibility_mean_null_array_sub(3,:));


flex_ml_rt_null_flow_less_bore_sub = (flex_flow_ml_rt_null_t_sub - flex_bore_ml_rt_null_t_sub)';
mean_flex_ml_rt_null_flow_less_bore_sub = mean(flex_ml_rt_null_flow_less_bore_sub);
vect_ones = ones(1000,1);
vect_mean_flex_ml_rt_null_flow_less_bore_sub = vect_ones * mean_flex_ml_rt_null_flow_less_bore_sub;
diff_mean_flex_ml_rt_null_flow_less_bore_sub = flex_ml_rt_null_flow_less_bore_sub - vect_mean_flex_ml_rt_null_flow_less_bore_sub;
s_diff_mean_flex_ml_rt_null_flow_less_bore_sub = diff_mean_flex_ml_rt_null_flow_less_bore_sub.^2;
ss_diff_mean_flex_ml_rt_null_flow_less_bore_sub = sum(s_diff_mean_flex_ml_rt_null_flow_less_bore_sub);
one_over_m_less_one = 1/999;
mss_diff_mean_flex_ml_rt_null_flow_less_bore_sub = one_over_m_less_one * ss_diff_mean_flex_ml_rt_null_flow_less_bore_sub;
seb_mss_diff_mean_flex_ml_rt_null_flow_less_bore_sub = sqrt(mss_diff_mean_flex_ml_rt_null_flow_less_bore_sub);
seb_mss_diff_mean_flex_ml_rt_null_flow_less_bore_t_sub = seb_mss_diff_mean_flex_ml_rt_null_flow_less_bore_sub.';
num_emp_flex_flow_less_bore_sub = flex_flow_ml_sub - flex_bore_ml_sub;
t_val_flex_flow_less_bore_sub = num_emp_flex_flow_less_bore_sub./seb_mss_diff_mean_flex_ml_rt_null_flow_less_bore_t_sub

flex_ml_rt_null_flow_less_frust_sub = (flex_flow_ml_rt_null_t_sub - flex_frust_ml_rt_null_t_sub)';
mean_flex_ml_rt_null_flow_less_frust_sub = mean(flex_ml_rt_null_flow_less_frust_sub);
vect_ones = ones(1000,1);
vect_mean_flex_ml_rt_null_flow_less_frust_sub = vect_ones * mean_flex_ml_rt_null_flow_less_frust_sub;
diff_mean_flex_ml_rt_null_flow_less_frust_sub = flex_ml_rt_null_flow_less_frust_sub - vect_mean_flex_ml_rt_null_flow_less_frust_sub;
s_diff_mean_flex_ml_rt_null_flow_less_frust_sub = diff_mean_flex_ml_rt_null_flow_less_frust_sub.^2;
ss_diff_mean_flex_ml_rt_null_flow_less_frust_sub = sum(s_diff_mean_flex_ml_rt_null_flow_less_frust_sub);
one_over_m_less_one = 1/999;
mss_diff_mean_flex_ml_rt_null_flow_less_frust_sub = one_over_m_less_one * ss_diff_mean_flex_ml_rt_null_flow_less_frust_sub;
seb_mss_diff_mean_flex_ml_rt_null_flow_less_frust_sub = sqrt(mss_diff_mean_flex_ml_rt_null_flow_less_frust_sub);
seb_mss_diff_mean_flex_ml_rt_null_flow_less_frust_t_sub = seb_mss_diff_mean_flex_ml_rt_null_flow_less_frust_sub.';
num_emp_flex_flow_less_frust_sub = flex_flow_ml_sub - flex_frust_ml_sub;
t_val_flex_flow_less_frust_sub = num_emp_flex_flow_less_frust_sub./seb_mss_diff_mean_flex_ml_rt_null_flow_less_frust_t_sub

%% sub graph sub-frontal

flex_bore_ml_fron_sub = squeeze(mean(nodes_flexibility_mean(1,Frontoparietalsubcortical), 2))
flex_flow_ml_fron_sub = squeeze(mean(nodes_flexibility_mean(2,Frontoparietalsubcortical), 2))
flex_frust_ml_fron_sub = squeeze(mean(nodes_flexibility_mean(3,Frontoparietalsubcortical), 2))

flex_bore_ml_rt_null_t_fron_sub = squeeze(full_flexibility_mean_null_array_fron_sub(1,:));
flex_flow_ml_rt_null_t_fron_sub = squeeze(full_flexibility_mean_null_array_fron_sub(2,:));
flex_frust_ml_rt_null_t_fron_sub = squeeze(full_flexibility_mean_null_array_fron_sub(3,:));


flex_ml_rt_null_flow_less_bore_fron_sub = (flex_flow_ml_rt_null_t_fron_sub - flex_bore_ml_rt_null_t_fron_sub)';
mean_flex_ml_rt_null_flow_less_bore_fron_sub = mean(flex_ml_rt_null_flow_less_bore_fron_sub);
vect_ones = ones(1000,1);
vect_mean_flex_ml_rt_null_flow_less_bore_fron_sub = vect_ones * mean_flex_ml_rt_null_flow_less_bore_fron_sub;
diff_mean_flex_ml_rt_null_flow_less_bore_fron_sub = flex_ml_rt_null_flow_less_bore_fron_sub - vect_mean_flex_ml_rt_null_flow_less_bore_fron_sub;
s_diff_mean_flex_ml_rt_null_flow_less_bore_fron_sub = diff_mean_flex_ml_rt_null_flow_less_bore_fron_sub.^2;
ss_diff_mean_flex_ml_rt_null_flow_less_bore_fron_sub = sum(s_diff_mean_flex_ml_rt_null_flow_less_bore_fron_sub);
one_over_m_less_one = 1/999;
mss_diff_mean_flex_ml_rt_null_flow_less_bore_fron_sub = one_over_m_less_one * ss_diff_mean_flex_ml_rt_null_flow_less_bore_fron_sub;
seb_mss_diff_mean_flex_ml_rt_null_flow_less_bore_fron_sub = sqrt(mss_diff_mean_flex_ml_rt_null_flow_less_bore_fron_sub);
seb_mss_diff_mean_flex_ml_rt_null_flow_less_bore_t_fron_sub = seb_mss_diff_mean_flex_ml_rt_null_flow_less_bore_fron_sub.';
num_emp_flex_flow_less_bore_fron_sub = flex_flow_ml_fron_sub - flex_bore_ml_fron_sub;
t_val_flex_flow_less_bore_fron_sub = num_emp_flex_flow_less_bore_fron_sub./seb_mss_diff_mean_flex_ml_rt_null_flow_less_bore_t_fron_sub

flex_ml_rt_null_flow_less_frust_fron_sub = (flex_flow_ml_rt_null_t_fron_sub - flex_frust_ml_rt_null_t_fron_sub)';
mean_flex_ml_rt_null_flow_less_frust_fron_sub = mean(flex_ml_rt_null_flow_less_frust_fron_sub);
vect_ones = ones(1000,1);
vect_mean_flex_ml_rt_null_flow_less_frust_fron_sub = vect_ones * mean_flex_ml_rt_null_flow_less_frust_fron_sub;
diff_mean_flex_ml_rt_null_flow_less_frust_fron_sub = flex_ml_rt_null_flow_less_frust_fron_sub - vect_mean_flex_ml_rt_null_flow_less_frust_fron_sub;
s_diff_mean_flex_ml_rt_null_flow_less_frust_fron_sub = diff_mean_flex_ml_rt_null_flow_less_frust_fron_sub.^2;
ss_diff_mean_flex_ml_rt_null_flow_less_frust_fron_sub = sum(s_diff_mean_flex_ml_rt_null_flow_less_frust_fron_sub);
one_over_m_less_one = 1/999;
mss_diff_mean_flex_ml_rt_null_flow_less_frust_fron_sub = one_over_m_less_one * ss_diff_mean_flex_ml_rt_null_flow_less_frust_fron_sub;
seb_mss_diff_mean_flex_ml_rt_null_flow_less_frust_fron_sub = sqrt(mss_diff_mean_flex_ml_rt_null_flow_less_frust_fron_sub);
seb_mss_diff_mean_flex_ml_rt_null_flow_less_frust_t_fron_sub = seb_mss_diff_mean_flex_ml_rt_null_flow_less_frust_fron_sub.';
num_emp_flex_flow_less_frust_fron_sub = flex_flow_ml_fron_sub - flex_frust_ml_fron_sub;
t_val_flex_flow_less_frust_fron_sub = num_emp_flex_flow_less_frust_fron_sub./seb_mss_diff_mean_flex_ml_rt_null_flow_less_frust_t_fron_sub


%% inferential testing - modularity

modu_bore_ml = squeeze(mean(modularity_mean(:,1), 1))
modu_flow_ml = squeeze(mean(modularity_mean(:,2), 1))
modu_frust_ml = squeeze(mean(modularity_mean(:,3), 1))

modu_bore_ml_rt_null_t = squeeze(full_modularity_mean_null_array(1,:));
modu_flow_ml_rt_null_t = squeeze(full_modularity_mean_null_array(2,:));
modu_frust_ml_rt_null_t = squeeze(full_modularity_mean_null_array(3,:));


modu_ml_rt_null_flow_less_bore = (modu_flow_ml_rt_null_t - modu_bore_ml_rt_null_t)';
mean_modu_ml_rt_null_flow_less_bore = mean(modu_ml_rt_null_flow_less_bore);
vect_ones = ones(1000,1);
vect_mean_modu_ml_rt_null_flow_less_bore = vect_ones * mean_modu_ml_rt_null_flow_less_bore;
diff_mean_modu_ml_rt_null_flow_less_bore = modu_ml_rt_null_flow_less_bore - vect_mean_modu_ml_rt_null_flow_less_bore;
s_diff_mean_modu_ml_rt_null_flow_less_bore = diff_mean_modu_ml_rt_null_flow_less_bore.^2;
ss_diff_mean_modu_ml_rt_null_flow_less_bore = sum(s_diff_mean_modu_ml_rt_null_flow_less_bore);
one_over_m_less_one = 1/999;
mss_diff_mean_modu_ml_rt_null_flow_less_bore = one_over_m_less_one * ss_diff_mean_modu_ml_rt_null_flow_less_bore;
seb_mss_diff_mean_modu_ml_rt_null_flow_less_bore = sqrt(mss_diff_mean_modu_ml_rt_null_flow_less_bore);
seb_mss_diff_mean_modu_ml_rt_null_flow_less_bore_t = seb_mss_diff_mean_modu_ml_rt_null_flow_less_bore.';
num_emp_modu_flow_less_bore = modu_flow_ml - modu_bore_ml;
t_val_modu_flow_less_bore = num_emp_modu_flow_less_bore./seb_mss_diff_mean_modu_ml_rt_null_flow_less_bore_t

modu_ml_rt_null_flow_less_frust = (modu_flow_ml_rt_null_t - modu_frust_ml_rt_null_t)';
mean_modu_ml_rt_null_flow_less_frust = mean(modu_ml_rt_null_flow_less_frust);
vect_ones = ones(1000,1);
vect_mean_modu_ml_rt_null_flow_less_frust = vect_ones * mean_modu_ml_rt_null_flow_less_frust;
diff_mean_modu_ml_rt_null_flow_less_frust = modu_ml_rt_null_flow_less_frust - vect_mean_modu_ml_rt_null_flow_less_frust;
s_diff_mean_modu_ml_rt_null_flow_less_frust = diff_mean_modu_ml_rt_null_flow_less_frust.^2;
ss_diff_mean_modu_ml_rt_null_flow_less_frust = sum(s_diff_mean_modu_ml_rt_null_flow_less_frust);
one_over_m_less_one = 1/999;
mss_diff_mean_modu_ml_rt_null_flow_less_frust = one_over_m_less_one * ss_diff_mean_modu_ml_rt_null_flow_less_frust;
seb_mss_diff_mean_modu_ml_rt_null_flow_less_frust = sqrt(mss_diff_mean_modu_ml_rt_null_flow_less_frust);
seb_mss_diff_mean_modu_ml_rt_null_flow_less_frust_t = seb_mss_diff_mean_modu_ml_rt_null_flow_less_frust.';
num_emp_modu_flow_less_frust = modu_flow_ml - modu_frust_ml;
t_val_modu_flow_less_frust = num_emp_modu_flow_less_frust./seb_mss_diff_mean_modu_ml_rt_null_flow_less_frust_t

%% apply FDR to the p values to correct FWE
t_values = [t_val_modu_flow_less_frust;
    t_val_modu_flow_less_bore;
    t_val_flex_flow_less_bore;
    t_val_flex_flow_less_frust;
    t_val_flex_flow_less_bore_fron;
    t_val_flex_flow_less_frust_fron;
    t_val_flex_flow_less_bore_sub;
    t_val_flex_flow_less_frust_sub;
    t_val_flex_flow_less_bore_fron_sub;
    t_val_flex_flow_less_frust_fron_sub];

% convert the t values to one-tailed t test p values
tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution
tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2;              % 1-tailed t-distribution
p = zeros(size(t_values));

for i = 1:size(t_values)
    p(i) = 1- tdist1T(t_values(i),1000-1);
end

% apply FDR to the raw p values
adjusted_p = fdr_BH(p, .05);



%% Generate the consensus community assignment 
% first we compute the agreement matrix by agreement() in BCT toolbox
agreement_bore = agreement([squeeze(modules(1,:,:,1))' squeeze(modules(1,:,:,2))' ...
    squeeze(modules(1,:,:,3))' squeeze(modules(1,:,:,4))'])./4000;
% second we compute the consensus matrix by consensus_und() in BCT toolbox
consensus_bore = consensus_und(agreement_bore, 0.6, 1000);

% Then we repeate this for flow and frus
agreement_flow = agreement([squeeze(modules(2,:,:,1))' squeeze(modules(2,:,:,2))' ...
    squeeze(modules(2,:,:,3))' squeeze(modules(2,:,:,4))'])./4000;
consensus_flow = consensus_und(agreement_flow, 0.6, 1000);

agreement_frus = agreement([squeeze(modules(3,:,:,1))' squeeze(modules(3,:,:,2))' ...
    squeeze(modules(3,:,:,3))' squeeze(modules(3,:,:,4))'])./4000;
consensus_frus = consensus_und(agreement_frus, 0.6, 1000);

% Then we gather them together and plot the consensus matrix
consensus_matrix = [consensus_bore consensus_flow consensus_frus];
imagesc(consensus_matrix);

