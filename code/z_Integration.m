
ZJA_cw1
clear;
z_DR;
clear

Define_Constants

GNSS_solution = readmatrix('GNSS_solution.csv');
DR_solution = readmatrix('DR_solution.csv');

delta_v = zeros(3,1);
delta_r = zeros(3,1);

x_ini = [delta_v;delta_r];
P_ini = diag([100 100 100 0.01 0.01 0.01]);

t_s = 0.5;

x = x_ini;
P = P_ini;

%transition matrix
I = eye(3);
Z = zeros(3);
phi = [I Z;t_s*I I];

% system noise covariance
s_dr = 0.01;% velocity error 

Q = [s_dr*t_s*I 0.5*s_dr*t_s^2*I;
    0.5*s_dr*t_s^2*I (1/3)*s_dr*t_s^3*I];
%Propagate state & covariance
x_n = phi*x;
P_n = phi*P*phi.'+ Q;

%Calculate measurement matrix
H = [Z -I; -I Z];

%measurement noise covariance matrix
sigma_gr = 10;
sigma_gv = 0.05;
R = [sigma_gr^2*I Z; Z sigma_gv^2*I];

%calculate kalman gain matrix
K = P_n*H.'/(H*P_n*H.'+R);

%measurement innvoation
% delta_z = GNSS_solution()-DR_solution-[Z -I; -I Z]*x_n;
% 
% x_p = x_n+K*delta_z;
% P_p = (eye(6)-K*H)*P_n;
% 
% Correct_solution() = DR_sollution() - x_p;

    


