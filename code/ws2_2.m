clc;
clear;
Define_Constants;
%a) Initialise the Kalman filter
[x,P] = Initialise_GNSS_KF;
%b) (Kalman filter Step 1) Compute the transition matrix using
I = eye(3);
Z = zeros(3);
z_31 = zeros(3,1);
z_13 = zeros(1,3);
t_s = 1;
phi = [I t_s*I z_31 z_31;
    Z I z_31 z_31;
    z_13 z_13 1 t_s;
    z_13 z_13 0 1];
%c) (Step 2) Compute the system noise covariance matrix using
s_a = 5;
s_a_phi = 0.01;
s_a_f = 0.04;
Q = [s_a*t_s^3*I/3 s_a*t_s^2*I/2 z_31 z_31;
    s_a*t_s^2*I/2 s_a*t_s*I z_31 z_31;
    z_13 z_13 s_a_phi*t_s+s_a_f*t_s^3/3 s_a_f*t_s^2/2;
    z_13 z_13 s_a_f*t_s^2/2 s_a_f*t_s];

record = zeros(181,7);
for t = 1:181
    %d) (Step 3) Use the transition matrix to propagate the state estimates:
    x_n = phi*x;
    %e) (Step 4) Then use this to propagate the error covariance matrix:
    P_n = phi*P*phi.' + Q;
    %f) Predict the ranges from the approximate user position to each satellite
    Pseudo_ranges = readmatrix('Workshop2_Pseudo_ranges.csv');
    Pseudo_ranges_rate = readmatrix('Workshop2_Pseudo_range_rates.csv');
    time = Pseudo_ranges(t+1,1);
    r_aj_n = zeros(10,1);
    r_dot_aj_n = zeros(10,1);
    H = zeros(20,8);
    r_eb_e = x_n(1:3,:);
    v_eb_e = x_n(4:6,:);
    for i = 1:10
        j = Pseudo_ranges(1,i+1);
        % Calculate satellites position in ECEF
        % Inputs:
        %   time                  Current simulation time(s)
        %   j                     Satellite number
        % Outputs:
        %   sat_r_es_e       ECEF satellite position (m) 3x1 column vector
        %   sat_v_es_e       ECEF satellite velocity (m/s) 3x1 column vector
        [sat_r_es_e,sat_v_es_e] = Satellite_position_and_velocity(time,...
            j);
        r_as = sqrt(sum((sat_r_es_e.' - r_eb_e).^2));
        C_e_I = [1 , omega_ie*r_as/c, 0;
            -omega_ie*r_as/c, 1 , 0;
            0, 0, 1;];
        % consider sagnac effect when calculating predict pesudo ranges
        r_aj = sqrt(transpose(C_e_I*sat_r_es_e.'-r_eb_e)*(C_e_I*sat_r_es_e.'-r_eb_e));
        % Compute the line-of-sight unit vector
        u_aj_e = (C_e_I*sat_r_es_e.'-r_eb_e)/r_aj;
        %
        r_aj_n(i,1) = r_aj;
        r_dot_aj = u_aj_e.'*(C_e_I*(sat_v_es_e.'+Omega_ie*sat_r_es_e.')-(v_eb_e+Omega_ie*r_eb_e));
        r_dot_aj_n(i,:) = r_dot_aj;

        H(i,1:3) = -u_aj_e.';
        H(i,7) = 1;
        H(i+10,4:6) = -u_aj_e.';
        H(i+10,8) = 1;
    end

    %j) (Step 6) Compute the measurement noise covariance matrix
    r_sd = 10;
    rr_sd = 0.05;
    r_sd = r_sd^2*ones(1,10);
    rr_sd = rr_sd^2*ones(1,10);
    R =diag([r_sd rr_sd]);

    %k) (Step 7) Compute the Kalman gain matrix using
    K = P_n*H.'*inv(H*P_n*H.'+R);
    %l) (Step 8) Formulate the measurement innovation vector
    current_p_r = Pseudo_ranges(t+1,2:11).';
    current_p_rr = Pseudo_ranges_rate(t+1,2:11).';
    clock_offset = x_n(7);
    clock_drift = x_n(8);
    delta_z = [current_p_r-r_aj_n-clock_offset*ones(10,1);
        current_p_rr-r_dot_aj_n-clock_drift*ones(10,1)];

    %m) (Step 9) Update the state estimates using
    x_p = x_n + K*delta_z;

    %n) (Step 10) Update the error covariance matrix using
    P_p = (eye(8)-K*H)*P_n;

    % update state and covariance matrix for next epoch
    r = x_p(1:3,:);
    v = x_p(4:6,:);
    P = P_p;
    x = x_p;
    %convert cartesian ECEF position
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(r,v);
    L_b = L_b*rad_to_deg;
    lambda_b = lambda_b*rad_to_deg;
    %record solutions for current epoch
    record(t,1)=t-1;
    record(t,2)=L_b;
    record(t,3)=lambda_b;
    record(t,4)= h_b;
    record(t,5:7)=v_eb_n.';
end




