clear;

GNSS_solution = readmatrix('GNSS_solution.csv');

Define_Constants;
data =  readmatrix('Dead_reckoning.csv');
% seconds fowrward_speed heading_in_degrees

t_s = 0.5;
h = mean(GNSS_solution(:,4));

L_ini = GNSS_solution(1,2)*deg_to_rad;
lambda_ini = GNSS_solution(1,3)*deg_to_rad;
v_forward = (data(1,4)+data(1,5))/2;
v_n_ini = v_forward*cos(data(1,7)*deg_to_rad);
v_e_ini = v_forward*sin(data(1,7)*deg_to_rad);
v_ini = [v_n_ini;v_e_ini];

L_p = L_ini;
lambda_p = lambda_ini;
v_p = v_ini;

record = zeros(851,5);
record(1,1)= 0 ;
record(1,2)=GNSS_solution(1,2);
record(1,3)=GNSS_solution(1,3);
record(1,4:5)=v_ini.';

w = (8*0.5)/4;
psi_p = data(1,7)*deg_to_rad; %previous heading

for i = 1:850

    v_forward = (data(i+1,4)+data(i+1,5))/2;
    
%     psi_c = data(i+1,7)*deg_to_rad; %current heading
    psi_c = w*data(i+1,7)*deg_to_rad + (1-w)*(psi_p+t_s*data(i+1,6));
    v_ned = 0.5*[cos(psi_c)+cos(psi_p);sin(psi_c)+sin(psi_p)]*v_forward;

    [R_N,R_E]= Radii_of_curvature(L_p);

    L_c = L_p+(v_ned(1)*t_s)/(R_N+h);
    lambda_c = lambda_p+(v_ned(2)*t_s)/((R_E+h)*cos(L_c));

    v = 1.7*v_ned-0.7*v_p;

    record(i+1,1)= i*t_s ;
    record(i+1,2)=L_c*rad_to_deg;
    record(i+1,3)=lambda_c*rad_to_deg;
    record(i+1,4:5)=v.';

    L_p = L_c;
    lambda_p = lambda_c;
    v_p = v;
    psi_p = psi_c;
end

writematrix(record,'DR_solution.csv');

