clc;
clear;

Define_Constants;
ECEF_solution = readmatrix("Workshop2_GNSS_Pos_ECEF.csv");
t_s = 1;
r = [2447019; -5884199; -284783];
v = [184;77;0];
%a Initialise the Kalman filter state
uncer_r = 10;
uncer_x = 5;
p = [uncer_r^2 uncer_r^2 uncer_r^2 uncer_x^2 uncer_x^2 uncer_x^2];
P = diag(p);

%b  Compute the transition matrix using
I = eye(3);
Z = zeros(3);
phi = [I t_s*I; Z I]; %transition_matrix

%c Compute the system noise covariance matrix using
s_a = 5; %acceleration power spectral density
Q = [s_a*t_s^3*I/3 s_a*t_s^2*I/2; s_a*t_s^2*I/2 s_a*t_s*I]; %noise_covariance

record = zeros(181,7);

for i = 1:181

    x = [r;v];
    %d Use the transition matrix to propagate the state estimates
    x_n = phi * x;
    %e Then use this to propagate the error covariance matrix
    P_n = phi * P * phi.' + Q;
    %f Formulate the measurement matrix:
    H = [I Z];
    %g measurement noise covariance matrix
    e_s_d = 2.5;%error standard deviation
    R = diag([e_s_d^2 e_s_d^2 e_s_d^2]);
    %h Kalman gain matrix
    K = P_n*H.'*inv(H*P_n*H.'+ R);
    %i Formulate the measurement innovation vector
    delta_z = ECEF_solution(i,2:4).' - x_n(1:3,1);
    %j Update the state estimates using
    x_p = x_n + K*delta_z;
    %k Update the error covariance matrix using
    P_p = (eye(6)-K*H)*P_n;
    % update state and covariance matrix for next epoch
    r = x_p(1:3,:);
    v = x_p(4:6,:);
    P = P_p;
    %convert cartesian ECEF position
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(r,v);
    L_b = L_b*rad_to_deg;
    lambda_b = lambda_b*rad_to_deg;
    %record solutions for current epoch
    record(i,1)=i-1;
    record(i,2)=L_b;
    record(i,3)=lambda_b;
    record(i,4)= h_b;
    record(i,5:7)=v_eb_n.';

end




