function order_parameter = order_parameter(time_series)
% This function is to compute the order_parameter for metastability BOLD signal analysis .
% 
% The inpute time_series is a 3-D array, 1st dimension represent participant index,  2nd
% dimension represent n_TRs, the third dimension represents number of nodes.
% 
% The output order_parameter is a 2-D vector, the first dimension
% represents participants, second dimension represents the TR.
% 
%
%   
hil_time_series = zeros(size(time_series,1), size(time_series,2), size(time_series,3));
theta = zeros(size(time_series,1), size(time_series,2) , size(time_series,3));
real_n = zeros(size(time_series,1), size(time_series,2) , size(time_series,3));
rho = zeros(size(time_series,1), size(time_series,2) , size(time_series,3));
imag_n = zeros(size(time_series,1), size(time_series,2) , size(time_series,3));
complex_n = zeros(size(time_series,1), size(time_series,2) , size(time_series,3));
order_parameter = zeros(size(time_series,1), size(time_series,2));

for i = 1:size(time_series,1)
    hil_time_series(i,:,:) = hilbert(time_series(i,:,:));
    [theta(i,:,:), rho(i,:,:)] = cart2pol(real(hil_time_series(i,:,:)), imag(hil_time_series(i,:,:)));
    [real_n(i,:,:), imag_n(i,:,:)] = pol2cart(theta(i,:,:), 1);
    complex_n(i,:,:) = real_n(i,:,:) + imag_n(i,:,:)*1i;
    order_parameter(i,:) = (abs(sum(squeeze(complex_n(i,:,:))'))/size(complex_n(i,:,:),3))';
end