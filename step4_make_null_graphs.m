% Make random graphs for null distribution

% Turn Power Atlas sorted empirical adjacency matrix to array

bin_bore_ml_sort_mat = cat(3, bin_bore_ml_sort{:});
bin_flow_ml_sort_mat = cat(3, bin_flow_ml_sort{:});
bin_frust_ml_sort_mat = cat(3, bin_frust_ml_sort{:});

% Convert empirical adjacency matrix to array

contact_bore_sort = arrayToContactSeq(bin_bore_ml_sort_mat,0);
contact_flow_sort = arrayToContactSeq(bin_flow_ml_sort_mat,0);
contact_frust_sort = arrayToContactSeq(bin_flow_ml_sort_mat,0);

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

% Make one random permuted edge graph for each contact sequence with 100
% rewirings - this takes a LONG time

%contact_bore_re_null = randomizedEdges(contact_bore,100);
%contact_flow_re_null = randomizedEdges(contact_flow,100);
%contact_frust_re_null = randomizedEdges(contact_frust,100);