% Load all adjacency (correlation) matrices for each condition and window

% Import group level correlation matrix for bore condition and windows

data_bore_full = dlmread('/home/huskeyadmin/Desktop/merged/bore_full/adjacency_matrix/corel_matrix.txt', ' ');
data_bore_win1 = dlmread('/home/huskeyadmin/Desktop/merged/bore_win1/adjacency_matrix/corel_matrix.txt', ' ');
data_bore_win2 = dlmread('/home/huskeyadmin/Desktop/merged/bore_win2/adjacency_matrix/corel_matrix.txt', ' ');
data_bore_win3 = dlmread('/home/huskeyadmin/Desktop/merged/bore_win3/adjacency_matrix/corel_matrix.txt', ' ');
data_bore_win4 = dlmread('/home/huskeyadmin/Desktop/merged/bore_win4/adjacency_matrix/corel_matrix.txt', ' ');

% Import group level correlation matrix for flow condition and windows

data_flow_full = dlmread('/home/huskeyadmin/Desktop/merged/flow_full/adjacency_matrix/corel_matrix.txt', ' ');
data_flow_win1 = dlmread('/home/huskeyadmin/Desktop/merged/flow_win1/adjacency_matrix/corel_matrix.txt', ' ');
data_flow_win2 = dlmread('/home/huskeyadmin/Desktop/merged/flow_win2/adjacency_matrix/corel_matrix.txt', ' ');
data_flow_win3 = dlmread('/home/huskeyadmin/Desktop/merged/flow_win3/adjacency_matrix/corel_matrix.txt', ' ');
data_flow_win4 = dlmread('/home/huskeyadmin/Desktop/merged/flow_win4/adjacency_matrix/corel_matrix.txt', ' ');

% Import group level correlation matrix for frust condition and windows

data_frust_full = dlmread('/home/huskeyadmin/Desktop/merged/frust_full/adjacency_matrix/corel_matrix.txt', ' ');
data_frust_win1 = dlmread('/home/huskeyadmin/Desktop/merged/frust_win1/adjacency_matrix/corel_matrix.txt', ' ');
data_frust_win2 = dlmread('/home/huskeyadmin/Desktop/merged/frust_win2/adjacency_matrix/corel_matrix.txt', ' ');
data_frust_win3 = dlmread('/home/huskeyadmin/Desktop/merged/frust_win3/adjacency_matrix/corel_matrix.txt', ' ');
data_frust_win4 = dlmread('/home/huskeyadmin/Desktop/merged/frust_win4/adjacency_matrix/corel_matrix.txt', ' ');
