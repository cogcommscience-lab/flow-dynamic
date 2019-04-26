% Process Null Graphs

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

% Community detection on all null models

for i = 1:1000
    for j = 1:4
        comm_bore_ml_rt_null{1,i}(:,:,j) = community_louvain(array_bore_rt_null{1,i}(:,:,j),1,[],'potts');
    end
end

for i = 1:1000
    for j = 1:4
        comm_flow_ml_rt_null{1,i}(:,:,j) = community_louvain(array_flow_rt_null{1,i}(:,:,j),1,[],'potts');
    end
end

for i = 1:1000
    for j = 1:4
        comm_frust_ml_rt_null{1,i}(:,:,j) = community_louvain(array_frust_rt_null{1,i}(:,:,j),1,[],'potts');
    end
end

% Make each multidimensional array in the cell a pxn matrix

for i = 1:1000
    comm_bore_ml_rt_null_mat{1,i} = cat(2, comm_bore_ml_rt_null{1,i}(:,:,1), comm_bore_ml_rt_null{1,i}(:,:,2), comm_bore_ml_rt_null{1,i}(:,:,3),comm_bore_ml_rt_null{1,i}(:,:,4));
end     

for i = 1:1000
    comm_flow_ml_rt_null_mat{1,i} = cat(2, comm_flow_ml_rt_null{1,i}(:,:,1), comm_flow_ml_rt_null{1,i}(:,:,2), comm_flow_ml_rt_null{1,i}(:,:,3),comm_flow_ml_rt_null{1,i}(:,:,4));
end 

for i = 1:1000
    comm_frust_ml_rt_null_mat{1,i} = cat(2, comm_frust_ml_rt_null{1,i}(:,:,1), comm_frust_ml_rt_null{1,i}(:,:,2), comm_frust_ml_rt_null{1,i}(:,:,3),comm_frust_ml_rt_null{1,i}(:,:,4));
end 


% Get flexibility for thresholded and binarized multilayer adjacency matrix

for i = 1:1000
    flex_bore_ml_rt_null{1,i} = flexibility(comm_bore_ml_rt_null_mat{1,i},'temp');
end

for i = 1:1000
    flex_flow_ml_rt_null{1,i} = flexibility(comm_flow_ml_rt_null_mat{1,i},'temp');
end

for i = 1:1000
    flex_frust_ml_rt_null{1,i} = flexibility(comm_frust_ml_rt_null_mat{1,i},'temp');
end

% The full graph flexibility values are a cell, this makes them a regular
% array

flex_bore_ml_rt_null_array = cell2mat(flex_bore_ml_rt_null);
flex_flow_ml_rt_null_array = cell2mat(flex_flow_ml_rt_null);
flex_frust_ml_rt_null_array = cell2mat(flex_frust_ml_rt_null);

% And now lets transpose each matrix

flex_bore_ml_rt_null_t = flex_bore_ml_rt_null_array.';
flex_flow_ml_rt_null_t = flex_flow_ml_rt_null_array.';
flex_frust_ml_rt_null_t = flex_frust_ml_rt_null_array.';

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


% Extract FPCN Community Assignment

for i = 1:1000
    for j = 1:4
        comm_bore_fpcn_ml_rt_null{1,i}(:,:,j) = comm_bore_ml_rt_null{1,i}(157:181,1,j);
    end
end

for i = 1:1000
    for j = 1:4
        comm_flow_fpcn_ml_rt_null{1,i}(:,:,j) = comm_flow_ml_rt_null{1,i}(157:181,1,j);
    end
end

for i = 1:1000
    for j = 1:4
        comm_frust_fpcn_ml_rt_null{1,i}(:,:,j) = comm_frust_ml_rt_null{1,i}(157:181,1,j);
    end
end

% Make each multidimensional array in the cell a pxn matrix

for i = 1:1000
    comm_bore_fpcn_ml_rt_null_mat{1,i} = cat(2, comm_bore_fpcn_ml_rt_null{1,i}(:,:,1), comm_bore_fpcn_ml_rt_null{1,i}(:,:,2), comm_bore_fpcn_ml_rt_null{1,i}(:,:,3), comm_bore_fpcn_ml_rt_null{1,i}(:,:,4));
end

for i = 1:1000
    comm_flow_fpcn_ml_rt_null_mat{1,i} = cat(2, comm_flow_fpcn_ml_rt_null{1,i}(:,:,1), comm_flow_fpcn_ml_rt_null{1,i}(:,:,2), comm_flow_fpcn_ml_rt_null{1,i}(:,:,3), comm_flow_fpcn_ml_rt_null{1,i}(:,:,4));
end

for i = 1:1000
    comm_frust_fpcn_ml_rt_null_mat{1,i} = cat(2, comm_frust_fpcn_ml_rt_null{1,i}(:,:,1), comm_frust_fpcn_ml_rt_null{1,i}(:,:,2), comm_frust_fpcn_ml_rt_null{1,i}(:,:,3), comm_frust_fpcn_ml_rt_null{1,i}(:,:,4));
end

% Get flexibility for thresholded and binarized multilayer adjacency matrix

for i = 1:1000
    flex_bore_fpcn_ml_rt_null{1,i} = flexibility(comm_bore_fpcn_ml_rt_null_mat{1,i},'temp');
end

for i = 1:1000
    flex_flow_fpcn_ml_rt_null{1,i} = flexibility(comm_flow_fpcn_ml_rt_null_mat{1,i},'temp');
end

for i = 1:1000
    flex_frust_fpcn_ml_rt_null{1,i} = flexibility(comm_frust_fpcn_ml_rt_null_mat{1,i},'temp');
end

% The full graph flexibility values are a cell, this makes them a regular
% array

flex_bore_fpcn_ml_rt_null_array = cell2mat(flex_bore_fpcn_ml_rt_null);
flex_flow_fpcn_ml_rt_null_array = cell2mat(flex_flow_fpcn_ml_rt_null);
flex_frust_fpcn_ml_rt_null_array = cell2mat(flex_frust_fpcn_ml_rt_null);

% And now lets transpose each matrix

flex_bore_fpcn_ml_rt_null_t = flex_bore_fpcn_ml_rt_null_array.';
flex_flow_fpcn_ml_rt_null_t = flex_flow_fpcn_ml_rt_null_array.';
flex_frust_fpcn_ml_rt_null_t = flex_frust_fpcn_ml_rt_null_array.';

% And lets get some summary stats on each

flex_bore_fpcn_ml_rt_null_mean = mean(flex_bore_fpcn_ml_rt_null_t);
flex_bore_fpcn_ml_rt_null_sd = std(flex_bore_fpcn_ml_rt_null_t);
flex_bore_fpcn_ml_rt_null_se = std(flex_bore_fpcn_ml_rt_null_t)/sqrt(length(flex_bore_fpcn_ml_rt_null_t));

flex_flow_fpcn_ml_rt_null_mean = mean(flex_flow_fpcn_ml_rt_null_t);
flex_flow_fpcn_ml_rt_null_sd = std(flex_flow_fpcn_ml_rt_null_t);
flex_flow_fpcn_ml_rt_null_se = std(flex_flow_fpcn_ml_rt_null_t)/sqrt(length(flex_flow_fpcn_ml_rt_null_t));

flex_frust_fpcn_ml_rt_null_mean = mean(flex_frust_fpcn_ml_rt_null_t);
flex_frust_fpcn_ml_rt_null_sd = std(flex_frust_fpcn_ml_rt_null_t);
flex_frust_fpcn_ml_rt_null_se = std(flex_frust_fpcn_ml_rt_null_t)/sqrt(length(flex_frust_fpcn_ml_rt_null_t));


% Extract subcort Community Assignment

for i = 1:1000
    for j = 1:4
        comm_bore_subcort_ml_rt_null{1,i}(:,:,j) = comm_bore_ml_rt_null{1,i}(200:212,1,j);
    end
end

for i = 1:1000
    for j = 1:4
        comm_flow_subcort_ml_rt_null{1,i}(:,:,j) = comm_flow_ml_rt_null{1,i}(200:212,1,j);
    end
end

for i = 1:1000
    for j = 1:4
        comm_frust_subcort_ml_rt_null{1,i}(:,:,j) = comm_frust_ml_rt_null{1,i}(200:212,1,j);
    end
end

% Make each multidimensional array in the cell a pxn matrix

for i = 1:1000
    comm_bore_subcort_ml_rt_null_mat{1,i} = cat(2, comm_bore_subcort_ml_rt_null{1,i}(:,:,1), comm_bore_subcort_ml_rt_null{1,i}(:,:,2), comm_bore_subcort_ml_rt_null{1,i}(:,:,3), comm_bore_subcort_ml_rt_null{1,i}(:,:,4));
end

for i = 1:1000
    comm_flow_subcort_ml_rt_null_mat{1,i} = cat(2, comm_flow_subcort_ml_rt_null{1,i}(:,:,1), comm_flow_subcort_ml_rt_null{1,i}(:,:,2), comm_flow_subcort_ml_rt_null{1,i}(:,:,3), comm_flow_subcort_ml_rt_null{1,i}(:,:,4));
end

for i = 1:1000
    comm_frust_subcort_ml_rt_null_mat{1,i} = cat(2, comm_frust_subcort_ml_rt_null{1,i}(:,:,1), comm_frust_subcort_ml_rt_null{1,i}(:,:,2), comm_frust_subcort_ml_rt_null{1,i}(:,:,3), comm_frust_subcort_ml_rt_null{1,i}(:,:,4));
end

% Get flexibility for thresholded and binarized multilayer adjacency matrix

for i = 1:1000
    flex_bore_subcort_ml_rt_null{1,i} = flexibility(comm_bore_subcort_ml_rt_null_mat{1,i},'temp');
end

for i = 1:1000
    flex_flow_subcort_ml_rt_null{1,i} = flexibility(comm_flow_subcort_ml_rt_null_mat{1,i},'temp');
end

for i = 1:1000
    flex_frust_subcort_ml_rt_null{1,i} = flexibility(comm_frust_subcort_ml_rt_null_mat{1,i},'temp');
end

% The full graph flexibility values are a cell, this makes them a regular
% array

flex_bore_subcort_ml_rt_null_array = cell2mat(flex_bore_subcort_ml_rt_null);
flex_flow_subcort_ml_rt_null_array = cell2mat(flex_flow_subcort_ml_rt_null);
flex_frust_subcort_ml_rt_null_array = cell2mat(flex_frust_subcort_ml_rt_null);

% And now lets transpose each matrix

flex_bore_subcort_ml_rt_null_t = flex_bore_subcort_ml_rt_null_array.';
flex_flow_subcort_ml_rt_null_t = flex_flow_subcort_ml_rt_null_array.';
flex_frust_subcort_ml_rt_null_t = flex_frust_subcort_ml_rt_null_array.';

% And lets get some summary stats on each

flex_bore_subcort_ml_rt_null_mean = mean(flex_bore_subcort_ml_rt_null_t);
flex_bore_subcort_ml_rt_null_sd = std(flex_bore_subcort_ml_rt_null_t);
flex_bore_subcort_ml_rt_null_se = std(flex_bore_subcort_ml_rt_null_t)/sqrt(length(flex_bore_subcort_ml_rt_null_t));

flex_flow_subcort_ml_rt_null_mean = mean(flex_flow_subcort_ml_rt_null_t);
flex_flow_subcort_ml_rt_null_sd = std(flex_flow_subcort_ml_rt_null_t);
flex_flow_subcort_ml_rt_null_se = std(flex_flow_subcort_ml_rt_null_t)/sqrt(length(flex_flow_subcort_ml_rt_null_t));

flex_frust_subcort_ml_rt_null_mean = mean(flex_frust_subcort_ml_rt_null_t);
flex_frust_subcort_ml_rt_null_sd = std(flex_frust_subcort_ml_rt_null_t);
flex_frust_subcort_ml_rt_null_se = std(flex_frust_subcort_ml_rt_null_t)/sqrt(length(flex_frust_subcort_ml_rt_null_t));
