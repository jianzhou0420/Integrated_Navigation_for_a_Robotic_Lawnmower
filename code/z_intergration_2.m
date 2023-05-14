
ZJA_cw1
clear
DR_B
clear
DR_solution = readmatrix('DR_solution.csv');
GNSS_solution = readmatrix('GNSS_solution.csv');
Define_Constants

x = zeros(4,1); % v_n v_e l lamda
% error covariance
sigma_v = 0.1;%0.1;% velocity uncertainty
sigma_r = 10;%10;% position uncertainty
L_b = GNSS_solution(1,2)*deg_to_rad;
lambda_b = GNSS_solution(1,3)*deg_to_rad;
h_b = GNSS_solution(1,4);
[R_N,R_E]= Radii_of_curvature(L_b);

P = diag([sigma_v^2 sigma_v^2 ...
    sigma_r^2/(R_N+h_b)^2 sigma_r^2/(R_E+h_b)^2*cos(L_b)^2]);

result = zeros(851,5);
for i = 1:851
    %Kalman filter

    L_b = GNSS_solution(i,2)*deg_to_rad;
    lambda_b = GNSS_solution(i,3)*deg_to_rad;
    h_b = GNSS_solution(i,4);
    [R_N,R_E]= Radii_of_curvature(L_b);

    %transition matrix
    t_s = 0.5;
    phi = eye(4);
    phi(3,1) = t_s/(R_N+h_b);
    phi(4,2) = t_s/((R_E+h_b)*cos(L_b));

    %noise covariance matrix
    s_dr = 100;%0.001;%0.2;
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
    sigma_gr = 10;%5;% position measuremnet error standard deviation
    sigma_gv = 0.05;%0.02;% velocity measuremnet error standard deviation
    h_c = GNSS_solution(i,4);
    L_c = GNSS_solution(i,2)*deg_to_rad;
    lambda_c = GNSS_solution(i,3)*deg_to_rad;
    [R_N,R_E]= Radii_of_curvature(L_c);
    R = diag([sigma_gr^2/(R_N+h_c)^2 sigma_gr^2/((R_E+h_c)^2*cos(L_c)^2) ...
        sigma_gv^2 sigma_gv^2]);


    %Compute the Kalman gain matrix
    K = P_n*H.'*inv(H*P_n*H.'+R);
    %measurement innovation vector
    v_n_G = GNSS_solution(i,5);
    v_e_G = GNSS_solution(i,6);
    DR_record = DR_solution(i,2:5).';
    DR_record = DR_record*deg_to_rad;

    delta_z = [L_c;lambda_c;v_n_G;v_e_G]-DR_record-H*x_n;

    x_p = x_n + K*delta_z;
    P_p = (eye(4)-K*H)*P_n;

    x_p_record = x_p;
    x_p_record(3) = x_p_record(3)*rad_to_deg;
    x_p_record(4) = x_p_record(4)*rad_to_deg;

    result(i,1) = (i-1)*t_s;
    result(i,2:5) = DR_solution(i,2:5).' + H*x_p_record;

    x = x_p;
    P = P_p;

end
figure;
plot(GNSS_solution(:,3),GNSS_solution(:,2));%
hold on
%graph of longitude and latitude
x =  DR_solution(:,3);
y =  DR_solution(:,2);
% set(gcf,'Color',[0.9 0.9 0.9]);
plot(x,y);
hold on
plot(result(:,3),result(:,2));%
% set(gcf,'Color',[0.9 0.9 0.9]);
title('plot of the latitude and longitude of the lawnmower');
xlabel('longitude');
ylabel('latitude');
legend('GNSS','DR','GNSS')
grid on;



