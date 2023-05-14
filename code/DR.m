clc;
clear;
syms x;
Define_Constants;
DR_result = [];
DR_data = readmatrix("Dead_reckoning.csv");
d_t = 0.5;
%length l and width w of the lawnmower
l = 0.5;
d = 0.4;
% heading weight
W_m = (4*0.5)/4; 
%initial height ,latitude and longitude got from GNSS
h   = 38.821189;
l_a = 51.5092543;
l_o = -0.1610451; 
%compute the radii of curvature
[R_N,R_E] = Radii_of_curvature(l_a*deg_to_rad);

% initial statement at time 0
    % smooth heading
  
    DR_result(1,6) = W_m*(DR_data(1,7))*rad_to_deg + (1-W_m)*( DR_data(1,6)*d_t);
    %speed
    % average speed
    v_k = (DR_data(1,4)+DR_data(1,5))/2;
    theta = DR_result(1,6)*d_t*deg_to_rad;
    v_NE = 1/2*[cos(DR_result(1,6)*deg_to_rad);sin(DR_result(1,6)*deg_to_rad)]*v_k;
    DR_result(1,4) = v_NE(1,1);
    DR_result(1,5) = v_NE(2,1);
    % latitude
    DR_result(1,2) =  l_a;
    % longitude
    DR_result(1,3) =  l_o;
    %time
    DR_result(1,1) = 0;
    
 


for i= 2:851
    %time
    DR_result(i,1) = vpa((i-1)*d_t,7);

    %solve the heading
    d_w = 0.003;% growing heading standard deviation of gyro
    W_m = ((1+0.003*d_t)*0.5)/4;
    DR_result(i,6) = W_m*(DR_data(i,7)) + (1-W_m)*(DR_result(i-1,6) + (DR_data(i,6) + DR_data(i-1,6))/2*d_t*rad_to_deg);
    %calculate the latitude , longitude and speed
    v_k = (DR_data(i,4)+DR_data(i,5))/2; 
    v_NE = 1/2*[cos(DR_result(i,6)*deg_to_rad) + cos(DR_result(i,6)*deg_to_rad); 
                sin(DR_result(i,6)*deg_to_rad) + sin(DR_result(i,6)*deg_to_rad)]*v_k;
    DR_result(i,4) = v_NE(1,1);
    DR_result(i,5) = v_NE(2,1);


    [R_N,R_E] = Radii_of_curvature(DR_result(i-1,2));
    % latitude
    DR_result(i,2) = DR_result(i-1,2) + rad_to_deg*v_NE(1,1)*d_t/(R_N +h);
    % longitude
    DR_result(i,3) = DR_result(i-1,3) + rad_to_deg*v_NE(2,1)*d_t/(R_E +h)*cos(DR_result(i,2));

    


end 

writematrix(DR_result,'DR_solution.csv');

%graph of longitude and latitude
% x =  DR_result(:,3);
% y =  DR_result(:,2);
% figure;
% set(gcf,'Color',[0.9 0.9 0.9]);
% plot(x,y);
% xlabel('longitude');
% ylabel('latitude');
% title('plot of the latitude and longitude of the lawnmower');
% 
% 
% %graph of velocity
% x  =  DR_result(:,1);
% y1 =  DR_result(:,4);
% y2 =  DR_result(:,5);
% figure;
% set(gcf,'Color',[0.9 0.9 0.9]);
% plot(x,y1,x,y2,'--');
% legend('line_1','line_2','Location','northeast','Orientation','vertical');
% xlabel('time(s)');
% ylabel('velocity(m/s)');
% title('plot of the velocity');
% 
% 
% %graph of heading
% x  =  DR_result(:,1);
% y =  DR_result(:,6);
% 
% figure;
% set(gcf,'Color',[0.9 0.9 0.9]);
% plot(x,y);
% xlabel('time(s)');
% ylabel('heading(deg)');
% title('plot of the heading');

