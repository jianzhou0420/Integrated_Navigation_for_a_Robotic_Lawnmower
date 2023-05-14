

data = readmatrix('Workshop3_GNSS_Pos_Vel_NED.csv');

x = zeros(4,1);
% error covariance
sigma_v = 0.1;% velocity uncertainty
sigma_r = 10;% position uncertainty
L_b = data(1,2)*deg_to_rad;
lambda_b = data(1,3)*deg_to_rad;
h_b = data(1,4);
[R_N,R_E]= Radii_of_curvature(L_b);

P = diag([sigma_v^2 sigma_v^2 ...
    sigma_r^2/(R_N+h_b)^2 sigma_r^2/(R_E+h_b)^2*cos(L_b)^2]);

result = zeros(351,5);
for i = 1:351
    %Kalman filter

    L_b = data(i,2)*deg_to_rad;
    lambda_b = data(i,3)*deg_to_rad;
    h_b = data(i,4);
    [R_N,R_E]= Radii_of_curvature(L_b);

    %transition matrix
    t_s = 0.5;
    phi = eye(4);
    phi(3,1) = t_s/(R_N+h_b);
    phi(4,2) = t_s/((R_E+h_b)*cos(L_b));

    %noise covariance matrix
    s_dr = 0.2;
    Q = zeros(4);
    Q(1,1) = s_dr*t_s;
    Q(2,2) = Q(1,1);
    Q(3,3) = (1/3)*(s_dr*t_s^3)/(R_N+h_b)^2;
    Q(4,4) = (1/3)*(s_dr*t_s^3)/((R_E+h_b)^2*cos(L_b)^2);
    Q(1,3) = (1/2)*(s_dr*t_s^2)/(R_N+h_b);
    Q(3,1) = Q(1,3);
    Q(2,4) = (1/2)*(s_dr*t_s^2)/((R_E+h_b)*cos(L_b));
    Q(4,2) = Q(2,4);

    x_n = phi*x;

    P_n = phi*P*phi.' + Q;

    H = [0 0 -1 0;
        0 0 0 -1;
        -1 0 0 0;
        0 -1 0 0];

    %measurment noise covariance matrix
    sigma_gr = 5;% position measuremnet error standard deviation
    sigma_gv = 0.02;% velocity measuremnet error standard deviation
    h_c = data(i,4);
    L_c = data(i,2)*deg_to_rad;
    lambda_c = data(i,3)*deg_to_rad;
    [R_N,R_E]= Radii_of_curvature(L_c);
    R = diag([sigma_gr^2/(R_N+h_c)^2 sigma_gr^2/((R_E+h_c)^2*cos(L_c)^2) ...
        sigma_gv^2 sigma_gv^2]);


    %Compute the Kalman gain matrix
    K = P_n*H.'*inv(H*P_n*H.'+R);
    %measurement innovation vector
    v_n_G = data(i,5);
    v_e_G = data(i,6);
    DR_record = record(i,2:5).';
    DR_record(1) = DR_record(1)*deg_to_rad;
    DR_record(2) = DR_record(2)*deg_to_rad;

    delta_z = [L_c;lambda_c;v_n_G;v_e_G]-DR_record-H*x_n;

    x_p = x_n + K*delta_z;
    P_p = (eye(4)-K*H)*P_n;

    x_p_record = x_p;
    x_p_record(3) = x_p_record(3)*rad_to_deg;
    x_p_record(4) = x_p_record(4)*rad_to_deg;

    result(i,1) = (i-1)*t_s;
    result(i,2:5) = record(i,2:5).' + H*x_p_record;

    x = x_p;
    P = P_p;

end




