function [metastability, synchrony] = metastability(order_parameter)
% This function is to compute the metastability and synchrony BOLD signal analysis .
% 
% The inpute data is a 1-D vector of order parameter.
% 
% The output order_parameter is two scalars, the first one is metastability, and second one is synchrony 
% 
%
%   
    metastability = std(order_parameter);
    synchrony = mean(order_parameter);
end