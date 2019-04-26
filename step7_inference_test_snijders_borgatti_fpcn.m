% This code calculates a paired-samples t-test. It is the first formula
% explained in Snijders & Borgatti (1999, Connections) for calculating a
% paired samples t-test.

% Bootstraped SE for paired samples t-test Flow less Bore

% Dfferences in flexibility flow - bore

flex_fpcn_ml_rt_null_flow_less_bore = (flex_flow_fpcn_ml_rt_null_t - flex_bore_fpcn_ml_rt_null_t);

% Get the mean of the differences

mean_flex_fpcn_ml_rt_null_flow_less_bore = mean(flex_fpcn_ml_rt_null_flow_less_bore);

% Make a vector of the mean difference, length of each null vect

vect_ones = ones(1000,1);

vect_mean_flex_fpcn_ml_rt_null_flow_less_bore = vect_ones * mean_flex_fpcn_ml_rt_null_flow_less_bore;

% Subtract each difference by the mean difference

diff_mean_flex_fpcn_ml_rt_null_flow_less_bore = flex_fpcn_ml_rt_null_flow_less_bore - vect_mean_flex_fpcn_ml_rt_null_flow_less_bore;

% Square each element of the vector

s_diff_mean_flex_fpcn_ml_rt_null_flow_less_bore = diff_mean_flex_fpcn_ml_rt_null_flow_less_bore.^2;

% Sum all elements in the vector

ss_diff_mean_flex_fpcn_ml_rt_null_flow_less_bore = sum(s_diff_mean_flex_fpcn_ml_rt_null_flow_less_bore);

% Make 1/(m-1) where m = number of permuations

one_over_m_less_one = 1/999;

% Multiply 1/(m-1) by sum of vector elements

mss_diff_mean_flex_fpcn_ml_rt_null_flow_less_bore = one_over_m_less_one * ss_diff_mean_flex_fpcn_ml_rt_null_flow_less_bore;

% Take root of result and transpose

seb_mss_diff_mean_flex_fpcn_ml_rt_null_flow_less_bore = sqrt(mss_diff_mean_flex_fpcn_ml_rt_null_flow_less_bore);

seb_mss_diff_mean_flex_fpcn_ml_rt_null_flow_less_bore_t = seb_mss_diff_mean_flex_fpcn_ml_rt_null_flow_less_bore.';

% Make numerator in paired t-test flow - bore

num_emp_flex_fpcn_flow_less_bore = flex_flow_fpcn_ml - flex_bore_fpcn_ml;

% Calculate t-value

t_val_flex_fpcn_flow_less_bore = num_emp_flex_fpcn_flow_less_bore./seb_mss_diff_mean_flex_fpcn_ml_rt_null_flow_less_bore_t;



% Bootstraped SE for paired samples t-test Flow less frust

% Dfferences in flexibility flow - frust

flex_fpcn_ml_rt_null_flow_less_frust = (flex_flow_fpcn_ml_rt_null_t - flex_frust_fpcn_ml_rt_null_t);

% Get the mean of the differences

mean_flex_fpcn_ml_rt_null_flow_less_frust = mean(flex_fpcn_ml_rt_null_flow_less_frust);

% Make a vector of the mean difference, length of each null vect

vect_ones = ones(1000,1);

vect_mean_flex_fpcn_ml_rt_null_flow_less_frust = vect_ones * mean_flex_fpcn_ml_rt_null_flow_less_frust;

% Subtract each difference by the mean difference

diff_mean_flex_fpcn_ml_rt_null_flow_less_frust = flex_fpcn_ml_rt_null_flow_less_frust - vect_mean_flex_fpcn_ml_rt_null_flow_less_frust;

% Square each element of the vector

s_diff_mean_flex_fpcn_ml_rt_null_flow_less_frust = diff_mean_flex_fpcn_ml_rt_null_flow_less_frust.^2;

% Sum all elements in the vector

ss_diff_mean_flex_fpcn_ml_rt_null_flow_less_frust = sum(s_diff_mean_flex_fpcn_ml_rt_null_flow_less_frust);

% Make 1/(m-1) where m = number of permuations

one_over_m_less_one = 1/999;

% Multiply 1/(m-1) by sum of vector elements

mss_diff_mean_flex_fpcn_ml_rt_null_flow_less_frust = one_over_m_less_one * ss_diff_mean_flex_fpcn_ml_rt_null_flow_less_frust;

% Take root of result and transpose

seb_mss_diff_mean_flex_fpcn_ml_rt_null_flow_less_frust = sqrt(mss_diff_mean_flex_fpcn_ml_rt_null_flow_less_frust);

seb_mss_diff_mean_flex_fpcn_ml_rt_null_flow_less_frust_t = seb_mss_diff_mean_flex_fpcn_ml_rt_null_flow_less_frust.';

% Make numerator in paired t-test flow - frust

num_emp_flex_fpcn_flow_less_frust = flex_flow_fpcn_ml - flex_frust_fpcn_ml;

% Calculate t-value

t_val_flex_fpcn_flow_less_frust = num_emp_flex_fpcn_flow_less_frust./seb_mss_diff_mean_flex_fpcn_ml_rt_null_flow_less_frust_t;

