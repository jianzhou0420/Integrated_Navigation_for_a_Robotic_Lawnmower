clc;
clear all;
ZJA_ini;

Sa=0.001;
Scphi=0.01;
Scf=0.04;
T=0.5;
F=[eye(3),T*eye(3),zeros(3,1),zeros(3,1);
   zeros(3),eye(3),zeros(3,1),zeros(3,1);
   zeros(1,3),zeros(1,3),1,T;
   zeros(1,3),zeros(1,3),0,1];
Q=[Sa*T^3/3*eye(3),Sa*T^2/2*eye(3),zeros(3,1),zeros(3,1);
    Sa*T^2/2*eye(3),Sa*T*eye(3),zeros(3,1),zeros(3,1);
    zeros(1,3),zeros(1,3), Scphi*T+Scf*T^3/3, Scf*T^2/2;
    zeros(1,3),zeros(1,3), Scf*T^2/2, Scf*T];
X_old=X0;
P_old=P0;
record=zeros(m,7);
X_first=get_EFEC_From_GNSS(data1(1:2,:));

X_old(1:3)=X_first(1:3);
X_old(7)=X_first(4);


for t=1:1:m
    X_estimated=F*X_old;
    P_estimated=F*P_old*F.'+Q;
    
    
    %(f)(g)
    for i=1:1:n
        j=sid(i);
        time=tid(t);
        r_ea_m=X_estimated(1:3);
        v_ea_m=X_estimated(4:6);
    
        [r_ej(:,i),V_ej] = Satellite_position_and_velocity(time,j); % range vector of earth to satellite
        
        r_aj(:,i)=r_ej(:,i)-r_ea_m;
    
        r_aj_norm(1,i)=norm(r_aj(:,i),2);
    
        C=[1,                omega_ie*r_aj_norm(1,i)/c,  0;
               -omega_ie*r_aj_norm(1,i)/c, 1,                0;
               0,                0,                1]; 
    
        r_aj_estimated(:,i)=C*r_ej(:,i)-r_ea_m;
    
        r_aj_estimated_norm(i)=sqrt(r_aj_estimated(:,i).'* r_aj_estimated(:,i)); % predict the ranges from the approximate user position to each satellite.
        
        u_aj(:,i)=r_aj_estimated(:,i)/r_aj_estimated_norm(i); %line of sight unit vector. In other words, unit vector between user and equlivent satellite.
        
        V_ej=V_ej.';
        
        r_aj_estimated_norm_dot(i)=u_aj(:,i).'*(C*(V_ej+Omega_ie*r_ej(:,i))-(v_ea_m+Omega_ie*r_ea_m));
    end
    %(i)
    H=[-u_aj.',zeros(n,3),ones(n,1),zeros(n,1);
        zeros(n,3),-u_aj.',zeros(n,1),ones(n,1)];
    %(j)
    Rk=[10^2*eye(n),zeros(n);
        zeros(n),0.05^2*eye(n)];
    
    K=P_estimated*H.'/(H*P_estimated*H.'+Rk);
    
    this_pesudo_range = data1(t+1,2:end).';
    this_pesudo_range_rate = data2(t+1,2:end).';
    clock_offset=X_estimated(7);
    clock_drift=X_estimated(8);
    
    innovation=[this_pesudo_range - r_aj_estimated_norm - clock_offset;
                this_pesudo_range_rate - r_aj_estimated_norm_dot - clock_drift];
    
    %(m)
    X_new=X_estimated+K*innovation;
    P_new=(eye(8)-K*H)*P_estimated;

    X_old=X_new;
    P_old=P_new;
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(X_new(1:3),X_new(4:6));
    Latitude_estimate=L_b*rad_to_deg;
    Longitude_estimate=lambda_b*rad_to_deg;
    
    record(t,1) = (t-1)*0.5;
    record(t,2:4)=[Latitude_estimate,Longitude_estimate,h_b];
    record(t,5:7)=v_eb_n.';

end
% csvwrite('GNSSresult.csv',record);
writematrix(record,'GNSS_solution.csv');

% TODO Sa未知
% Rk，10，0.05是"reasona but not optimal"，我还在想怎么办


% plot(record(:,2),record(:,3));

