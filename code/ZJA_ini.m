Define_Constants;
% filepath1 = fullfile('data', 'Pseudo_ranges.csv');
% filepath2 = fullfile('data', 'Pseudo_range_rates.csv');
% data1 = table2array(readtable(filepath1));
% data2 = table2array(readtable(filepath2));

data1 = readmatrix('Pseudo_ranges.csv');
data2 = readmatrix('Pseudo_range_rates.csv');

[m,n]=size(data1);
m=m-1;
n=n-1;

sid=data1(1,2:end).';
tid=data1(2:end,1);

r_ej=zeros(3,n);
u_aj=zeros(3,n);

r_aj=zeros(3,n);
r_aj_norm=zeros(1,n);
r_aj_estimated=zeros(3,n);
r_aj_estimated_norm=zeros(n,1);
r_aj_estimated_norm_dot=zeros(n,1);

r_ea=zeros(3,1);

X0 = [  0; 0; 0;0; 0;  0; 0; 0]; 

% Initialise error covariance matrix
P0 =  zeros(8);
P0(1,1) = 1^2;
P0(2,2) = 1^2;
P0(3,3) = 1^2;
P0(4,4) = 1^2;
P0(5,5) = 1^2;
P0(6,6) = 1^2;
P0(7,7) = 1000000^2;
P0(8,8) = 200^2;
P0(8,8) = 100;
P0(7,7) = 100;

% P0 =  zeros(8);
% P0(1,1) = 10^2;
% P0(2,2) = 10^2;
% P0(3,3) = 10^2;
% P0(4,4) = 0.1^2;
% P0(5,5) = 0.1^2;
% P0(6,6) = 0.1^2;
% P0(7,7) = 10^2;
% P0(8,8) = 0.1^2;

