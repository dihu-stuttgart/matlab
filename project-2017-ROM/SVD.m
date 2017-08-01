function [] = SVD()
clc;
%A = csvread('V_m_time.csv');
A = csvread('outVm_time_space_dt_0.01_n_16_O1.csv');
%[U,S,V] = svd(A','econ');
[U,S,V] = svd(A,'econ');
M = V(:,1:5);
csvwrite('Ukt.csv',M);
end