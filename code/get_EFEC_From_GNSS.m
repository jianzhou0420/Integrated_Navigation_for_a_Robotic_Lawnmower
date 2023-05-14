function r_ea_all=get_EFEC_From_GNSS(data)
    % input data must be a 2xn matrix with n representing the number of
    % satellite, the first row is satellites' id and first column is the
    % time

    %% Variables definition
    % r_ej: distance vector between earth and jth satellite
    % r_eu: distance vector between earth and user
    % r_aj: distance vector between user and jth satellite
    % r_aj_norm: % absolute distance between user and jth satellite
    
    [m,n]=size(data);
    n=n-1;
    satellites_id=data(1,2:end).';
    time_id=data(2:end,1);
    Pseudo_Range=data(2,2:end).';


    Define_Constants;
    velocity_eb_n=0;
    r_eu=zeros(3,1);
    r_ej=zeros(3,n);
    u_aj=zeros(3,n);
    r_aj=zeros(3,n);
    r_aj_norm=zeros(1,n);
    r_aj_estimate_norm=zeros(n,1);
    X_new=zeros(4,1);
    X_old=zeros(4,1);
    
    for iterations=1:1:100000
        for i=1:1:n
            j=satellites_id(i);
            time=time_id(1);
               
            [r_ej(:,i),~] = Satellite_position_and_velocity(time,j); % range vector of earth to satellite
            r_aj(:,i)=r_ej(:,i)-r_eu;
            r_aj_norm(1,i)=norm(r_aj(:,i),2);
            C=[1,                omega_ie*r_aj_norm(1,i)/c,  0;
                   -omega_ie*r_aj_norm(1,i)/c, 1,                0;
                   0,                0,                1]; 
            media=C*r_ej(:,i)-r_eu;
            r_aj_estimate_norm(i)=sqrt(media.'*media); % predict the ranges from the approximate user position to each satellite.
            media=C*r_ej(:,i)-r_eu;
            u_aj(:,i)=media/r_aj_estimate_norm(i); %line of sight unit vector. In other words, unit vector between user and equlivent satellite.
        end
        
        %% (e)
        
        deltaZ=Pseudo_Range-r_aj_estimate_norm-X_old(4,1);
        H=[-u_aj.',ones(n,1)];
        X_new=X_old+inv(H.'*H)*H.'*deltaZ;
        r_eu=X_new(1:3);
        X_old=X_new;
%         [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(X_new(1:3),0);
%         Latitude_estimate=rad2deg(L_b);
%         Longitude_estimate=rad2deg(lambda_b);
        
    end
    r_ea_all=X_new;
end
