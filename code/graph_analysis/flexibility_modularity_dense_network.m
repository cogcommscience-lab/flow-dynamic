% Codes for inferential testing of flexibility and modularity for dynamic
% functional connectivity network of flow state vs. boredom and frustrution
% 
% function used from genLouvain (http://netwiki.amath.unc.edu/GenLouvain/GenLouvain)
% flexibility

clc;
close all;
clear;
load(['/data/time_series/output_gsr_seitzman.mat']);

M_bore = cat(4, bore.win_1, bore.win_2, bore.win_3, bore.win_4);
M_flow = cat(4, flow.win_1, flow.win_2, flow.win_3, flow.win_4);
M_frus = cat(4, frus.win_1, frus.win_2, frus.win_3, frus.win_4);
M = cat(5, M_bore, M_flow, M_frus);
M = permute(M, [1 5 4 2 3]);


for i = 1:size(M,1)
    for j = 1:size(M,2)
        for k = 1:size(M,3)
            M(i,j,k,:,:) = squeeze(M(i,j,k,:,:)) - diag(diag(squeeze(M(i,j,k,:,:))));
%            M(i,j,k,:,:) = weight_conversion( squeeze(M(i,j,k,:,:)),'normalize');
        end
    end
end

%% Define node keys

node_key = readtable('Label.txt', 'ReadVariableNames',false);
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


% 1st dim of M should be subject, 2nd should be condition, third should be
% window, fourth should be roi

n_sub = size(M,1);
n_ses = size(M,2);
n_win = size(M,3);
n_roi = size(M,4);

%% community detection parameters
gamma = 1;
omega = 1;
n_rep = 100;

%% Empty matrices to store data
modularity_mean = zeros(n_sub, n_rep, n_ses); % Mean modularity
modules = zeros(n_sub, n_ses, n_rep, n_roi, n_win); % Module assigment labels      

%% community detection
for sub = 1 : n_sub
    for ses = 1 : n_ses
    %--- define objects ---------------------------------------------------
        A = cell(1, n_win);
        B = spalloc(n_roi * n_win, n_roi * n_win,(n_roi + n_win) * n_roi* n_win);
        twomu = 0;
    
        %--- null model -------------------------------------------------------
        for win = 1 : n_win
            %--- copy network with positive weights thresholding --------------
            A{win} = squeeze(M(sub, ses, win, :, :) .* (M(sub, ses, win, :, :) > 0));
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
            fprintf('Subject = %i\n',sub);
            [S,Q] = genlouvain(B);
            Q = Q / twomu;
            S = reshape(S, n_roi, n_win);
            
            modularity_mean(sub, rep, ses) = Q;
            modules(sub, ses, rep, :, :) = S;
       end
    end
end



modularity_mean_mean = squeeze(mean(modularity_mean, 2));
figure;
boxplot(modularity_mean_mean);

subjects = [1:35];
t = table(subjects',modularity_mean_mean(:,1),modularity_mean_mean(:,2),modularity_mean_mean(:,3),...
'VariableNames',{'sub','bore','flow','frus'});
Meas = table([1 2 3]','VariableNames',{'Measurements'});
rm = fitrm(t,'bore-frus~sub','WithinDesign',Meas);
ranova(rm)

%% flexibility
for i = 1:size(modules,1)
    for j = 1:size(modules,2)
        for k = 1:size(modules,3)
            nodes_flexibility(i,j,k,:) = flexibility(squeeze(modules(i, j,k, :, :))');
        end
    end
end

%% flexibility inference testing using repeated ANOVA model

% testing the global brain network
nodes_flexibility_mean = squeeze(mean(nodes_flexibility, 3));
nodes_mean_flexibility_mean = squeeze(mean(nodes_flexibility_mean, 3));
subjects = [1:35];
t = table(subjects',nodes_mean_flexibility_mean(:,1),nodes_mean_flexibility_mean(:,2),nodes_mean_flexibility_mean(:,3),...
'VariableNames',{'sub','bore','flow','frus'});
Meas = table([1 2 3]','VariableNames',{'Measurements'});
rm = fitrm(t,'bore-frus~sub','WithinDesign',Meas);
ranova(rm)

% testing the FrontoParietal subnework
flex_bore_ml_fron = squeeze(mean(nodes_flexibility_mean(:,1,FrontoParietal), 3));
flex_flow_ml_fron = squeeze(mean(nodes_flexibility_mean(:,2,FrontoParietal), 3));
flex_frust_ml_fron = squeeze(mean(nodes_flexibility_mean(:,3,FrontoParietal), 3));
flex_fron = [flex_bore_ml_fron flex_flow_ml_fron flex_frust_ml_fron];
boxplot(flex_fron)

subjects = [1:35];
t = table(subjects',flex_fron(:,1),flex_fron(:,2),flex_fron(:,3),...
'VariableNames',{'sub','bore','flow','frus'});
Meas = table([1 2 3]','VariableNames',{'Measurements'});
rm = fitrm(t,'bore-frus~sub','WithinDesign',Meas);
ranova(rm)

% testing the Reward subnework
flex_bore_ml_sub = squeeze(mean(nodes_flexibility_mean(:,1,Reward), 3));
flex_flow_ml_sub = squeeze(mean(nodes_flexibility_mean(:,2,Reward), 3));
flex_frust_ml_sub = squeeze(mean(nodes_flexibility_mean(:,3,Reward), 3));
flex_sub = [flex_bore_ml_sub flex_flow_ml_sub flex_frust_ml_sub];
boxplot(flex_sub)

subjects = [1:35];
t = table(subjects',flex_sub(:,1),flex_sub(:,2),flex_sub(:,3),...
'VariableNames',{'sub','bore','flow','frus'});
Meas = table([1 2 3]','VariableNames',{'Measurements'});
rm = fitrm(t,'bore-frus~sub','WithinDesign',Meas);
ranova(rm)


% testing the FPCN&Reward subnework
flex_bore_ml_fron_sub = squeeze(mean(nodes_flexibility_mean(:,1,Frontoparietalsubcortical), 3));
flex_flow_ml_fron_sub = squeeze(mean(nodes_flexibility_mean(:,2,Frontoparietalsubcortical), 3));
flex_frust_ml_fron_sub = squeeze(mean(nodes_flexibility_mean(:,3,Frontoparietalsubcortical), 3));

flex_fron_sub = [flex_bore_ml_fron_sub flex_flow_ml_fron_sub flex_frust_ml_fron_sub];
boxplot(flex_fron_sub)

subjects = [1:35];
t = table(subjects',flex_fron_sub(:,1),flex_fron_sub(:,2),flex_fron_sub(:,3),...
'VariableNames',{'sub','bore','flow','frus'});
Meas = table([1 2 3]','VariableNames',{'Measurements'});
rm = fitrm(t,'bore-frus~sub','WithinDesign',Meas);
ranova(rm)

%% store the result into a matrix and save
individual_mod_fex = [modularity_mean_mean nodes_mean_flexibility_mean...
    flex_fron flex_sub flex_fron_sub];
% save('individual_mod_flex.mat', 'individual_mod_fex')
