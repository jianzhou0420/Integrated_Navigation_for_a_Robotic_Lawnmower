clear;
Define_Constants;
data =  readmatrix('Workshop3_Speed_Heading.csv');
% seconds fowrward_speed heading_in_degrees

t_s = 0.5;
h = 37.4;

L_ini = 50.4249580*deg_to_rad;
lambda_ini = -3.5957974*deg_to_rad;
v_n_ini = data(1,2)*cos(data(1,3)*deg_to_rad);
v_e_ini = data(1,2)*sin(data(1,3)*deg_to_rad);
v_ini = [v_n_ini;v_e_ini];

L_p = L_ini;
lambda_p = lambda_ini;
v_p = v_ini;

record = zeros(351,5);
record(1,1)= 0 ;
record(1,2)=50.4249580;
record(1,3)=-3.5957974;
record(1,4:5)=v_ini.';

for i = 1:350

    v_forward = data(i+1,2);
    psi_p = data(i,3)*deg_to_rad; %previous heading
    psi_c = data(i+1,3)*deg_to_rad; %current heading
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
end

