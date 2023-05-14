
%losd constants
Define_Constants;
%read pseudo range form csv file 
Pseudo_ranges = readmatrix('Workshop1_Pseudo_ranges');
Pseudo_ranges(:,4) = [];
[m,n] = size(Pseudo_ranges);
% Calculate the phone position in ECEF
% Inputs:
%   L_b           latitude (rad)
%   lambda_b      longitude (rad)
%   h_b           height (m)
%   v_eb_n        velocity of body frame w.r.t. ECEF frame, resolved along
%                 north, east, and down (m/s) 3x1 column vector
%
% Outputs:
%   r_eb_e        Cartesian position of body frame w.r.t. ECEF frame, resolved
%                 along ECEF-frame axes (m) 3x1 column vector
%   v_eb_e        velocity of body frame w.r.t. ECEF frame, resolved along
%                 ECEF-frame axes (m/s) 3x1 column vector

L_b = -33.821075 * deg_to_rad;
lambda_b = 151.188496 * deg_to_rad;
h_b = 120;
v_eb_n = 0;
[r_eb_e,v_eb_e] = pv_NED_to_ECEF(L_b,lambda_b,h_b,v_eb_n);
clock_offset = 0;

% r_eb_e = [0;0;0];
position_record = zeros(11,3);

for t = 1:m-1
    % set variable for following calculation
    time = Pseudo_ranges(t+1,1);
    r_aj_e_n = zeros(n-1,1);
    H = ones(n-1,4);

    % Calculate predict value 
    for i = 2:n
    
        j = Pseudo_ranges(1,i);
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
        r_aj_e_n(i-1,1) = r_aj;
        H(i-1,1:3) = -u_aj_e.';
    end
    %obtaion measured pseudo ranges
    current_p_r = Pseudo_ranges(t+1,2:n).';
    %calculate measurement innovation
    measurment_innovation =current_p_r - r_aj_e_n - ones(n-1,1)*clock_offset;
    x_predict = [r_eb_e;clock_offset];
    %Apply unweight least-squares
    x_update = x_predict + inv(H.'*H)*H.'*measurment_innovation;
    % Convert this Cartesian ECEF position solution to latitude
    % Inputs:
    %   r_eb_e        Cartesian position of body frame w.r.t. ECEF frame, resolved
    %                 along ECEF-frame axes (m) 3x1 column vector
    %   v_eb_e        velocity of body frame w.r.t. ECEF frame, resolved along
    %                 ECEF-frame axes (m/s) 3x1 column vector
    %
    % Outputs:
    %   L_b           latitude (rad)
    %   lambda_b      longitude (rad)
    %   h_b           height (m)
    %   v_eb_n        velocity of body frame w.r.t. ECEF frame, resolved along
    %                 north, east, and down (m/s) 3x1 column vector
    %Convert this Cartesian ECEF position solution to latitude, longitude and height
    [L_b_update,lambda_b_update,h_b_update,v_eb_n_update] =...
        pv_ECEF_to_NED(x_update(1:3,1),0);
    
    
    L_b_update = L_b_update*rad_to_deg;
    lambda_b_update = lambda_b_update*rad_to_deg;
    
    residuals_vector = (H*inv(H.'*H)*H.'-eye(n-1))*measurment_innovation;
    error_standard_deviation = 5;
    residuals_covariance = (eye(n-1)-H*inv(H.'*H)*H.')*error_standard_deviation;
    T = 6;

    for index = 1:n-1
        a = abs(residuals_vector(index,1));
        b = sqrt(residuals_covariance(index,index))*T;
        if  a>b 
            disp(time)
            disp(index)
            disp(a-b)
        end
    end

    position_record(t,1) = L_b_update;
    position_record(t,2) = lambda_b_update;
    position_record(t,3) = h_b_update;

    r_eb_e = x_update(1:3,1);
    clock_offset = x_update(4,1);
end




