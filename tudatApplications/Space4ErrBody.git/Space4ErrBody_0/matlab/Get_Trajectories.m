function [ evolutions ] = Get_Trajectories(evolutions,prop_path,depvar_path,interp_Ascent_path,interp_Descent_path,DV_mapped_Ascent_path,DV_mapped_Descent_path,v_i,gamma_i,pop_i,lon_i_rad,lat_f_deg,lon_f_deg, startEpoch)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

population =  numel(prop_path)/(numel(pop_i) + 1);

a = linspace(population,numel(prop_path),numel(pop_i) + 1);

% population =  100;
% pop_i = [ 0 1];
% a = linspace(population,200,numel(pop_i) + 1);

start = [1;a(1:end-1)'+1];
%C = [b a'];
%start = b;
finish = a';
%evolutions(numel(pop_i)).evolution(population,1).individual = [];

for k = 1:(numel(pop_i) + 1)
    
    evolutions(k).evolution = k;
    
    pp = 1;
    %     prop_source = prop_path(start(k):finish(k),:);
    depvar_source            = depvar_path(start(k):finish(k),:);
    interp_Ascent_source     = interp_Ascent_path(start(k):finish(k),:);
    interp_Descent_source    = interp_Descent_path(start(k):finish(k),:);
    DV_mapped_Ascent_source  = DV_mapped_Ascent_path(start(k):finish(k),:);
    DV_mapped_Descent_source = DV_mapped_Descent_path(start(k):finish(k),:);
    for p = 1:population
        
        % Open this file and read in the contents
        %         fid = fopen(prop_source{p,:});
        %         prop = dlmread(prop_source{p,:},',');
        %         fclose(fid);
        fid = fopen(depvar_source{p,:});
        depvar = dlmread(depvar_source{p,:},',');
        fclose(fid);
        
        fid = fopen(interp_Ascent_source{p,:});
        interp_Ascent = dlmread(interp_Ascent_source{p,:},',');
        fclose(fid);
        
        fid = fopen(interp_Descent_source{p,:});
        interp_Descent = dlmread(interp_Descent_source{p,:},',');
        fclose(fid);
        
        fid = fopen(DV_mapped_Ascent_source{p,:});
        DV_mapped_Ascent = dlmread(DV_mapped_Ascent_source{p,:},',');
        fclose(fid);
        
                fid = fopen(DV_mapped_Descent_source{p,:});
        DV_mapped_Descent = dlmread(DV_mapped_Descent_source{p,:},',');
        fclose(fid);
        
        % analyzed_simulations(k).tof = data(1,1) - data(end,1);
        %data = [2458419.95833333, -1730261.1974541, 3489942.65996152, 5041452.57112492, 1817.13843876246, 7053.76703722487, 474.632964784125]
        %
        %
        %         t   = prop(:,1);
        %         x_I = prop(:,2);
        %         y_I = prop(:,3);
        %         z_I = prop(:,4);
        %         R_I_vect = [x_I y_I z_I];
        %         %R_I_norm = sqrt(x_I.^2+y_I.^2+z_I.^2);
        %         %     %V_I_norm = sqrt(u_I.^2+v_I.^2+w_I.^2);
        %         tau_I_rad = atan2(y_I,x_I);
        %         %     I = find(tau_I_rad<0);
        %         %     tau_I_rad(I) = 2*pi + tau_I_rad(I);
        %         %     tau_I_rad2 = tau_I_rad;
        %         %     %tau_I_rad = [atan2(data(1,3),data(1,2)) atan2(data(end,3),data(end,2))];
        %         %     tau_I_deg = tau_I_rad*180/pi;
        %         %     tau_I_deg2 = tau_I_deg;
        %         %     delta_I_rad = asin(z_I./R_I_norm);
        %         %     I = find(delta_I_rad<0);
        %         %     delta_I_rad(I) = 2*pi + delta_I_rad(I);
        %         %     delta_I_deg = delta_I_rad*180/pi;
        %         %     delta_I_rad2 = delta_I_rad;
        %         %     delta_I_deg2 = delta_I_deg;
        %
        %
        %         %%%%%%%% This section beloy may be mitigated by transforming with TUDAT
        %         %%%%%%%% first and then exporting/printing.
        %         %%%%%%%%
        %         % Adjust Inertial Frame for plotting. TUDAT does all calculations
        %         % in the Inertial Frame, which is a function of the starting epoch.
        %         % TUDAT uses Earth's rotational ephemeris to identify the
        %         % corresponding state of the initial conditions. This creates a
        %         % basic rotational offset that can be easily reversed by
        %         % calculating the longidutindal difference. This difference is then
        %         % used with a transformation matrix that is then applied to the
        %         % entire trajectory.
        %
        %         % Calculate Transformation Matrix from original Inertial Frame to
        %         % Earth-Fixed Frame
        %
        %         initial_longitude_offset = lon_i_rad - tau_I_rad(1);
        %         st_adjust = [ zeros(numel(t),1) zeros(numel(t),1) -initial_longitude_offset*ones(numel(t),1) ];
        %         C_I_to_EF = DCM(st_adjust); % 3x3xN
        %
        %         % Restructure Inertial Frame position vector to 3x1xN
        %         [r,c] = size(R_I_vect);
        %         pages = numel(t);
        %         R_I_vect = permute(reshape(permute(R_I_vect,[2,1]),[r/pages,c,pages]),[2,1,3]);
        %
        %         % Calculate position vector and related vectors/values in Earth-Fixed Frame
        %         R_EF_vect = mmat(C_I_to_EF,R_I_vect,[1 2]);
        %         R_EF_vect = ipermute(squeeze(R_EF_vect),[2,1]);
        
        % x_EF(:,1) = squeeze(R_EF_vect(1,1,:));
        % y_EF(:,1) = squeeze(R_EF_vect(2,1,:));
        % z_EF(:,1) = squeeze(R_EF_vect(3,1,:));
        % R_EF_norm = sqrt(x_EF.^2+y_EF.^2+z_EF.^2);
        %
        %     %clear R_I_vect;
        %     R_I_vect = [x_I y_I z_I];
        %
        %     % Calcualte adjusted longitude and latitude
        %     delta_I_rad = asin(z_I./R_I_norm);
        %     delta_I_deg = delta_I_rad*180/pi;
        %     tau_I_rad = atan2(y_I,x_I);
        %     %tau_I_rad = [atan2(data(1,3),data(1,2)) atan2(data(end,3),data(end,2))];
        %     tau_I_deg = tau_I_rad*180/pi;
        %
        %
        %
        %     % Extract coordinates of final state from simulation data
        %     lat_f_calc_deg = delta_I_deg(end);
        %     lon_f_calc_deg = tau_I_deg(end);
        %
        %
        %
        %     [lat_f_calc_deg,lon_f_calc_deg] = latlon_orientation(lat_f_calc_deg,lon_f_calc_deg)
        %
        %
        %     % Calculate final state offset from IAD coordinates
        %     dif_lat_deg = lat_IAD_deg - lat_f_calc_deg;
        %     dif_lon_deg = lon_IAD_deg - lon_f_calc_deg_pos;
        %
        %
        % Following the rationale from before, throughout the propagation Earth
        % is also rotating. However, plotting this in a fixed image is not
        % possible. Hence, to plot the trajectory it must then be converted to
        % the Rotating Planetocentric Frame. This transformation would consider
        % the translation of the point mass with respect to the rotating
        % surface of the Earth. Unlike the simple parabolic trajectory that is
        % observed in the Inertial frame, this frame would show some odd
        % trajectories. Anyhow, it is a function of the epoch of the state
        % relative to the starting epoch.
        
        %         % Calculate Transformation Matrix from Inertial Frame to Rotational Frame
        %         omega_E = 7.2921150e-5;
        %         del_t = t-t(1);
        %         st_RI = [ zeros(numel(t),1) zeros(numel(t),1) (omega_E*del_t) ];
        %         C_EF_to_R = DCM(st_RI); % 3x3xN
        %
        %         % Restructure Inertial Frame position vector to 3x1xN
        %         [r,c] = size(R_EF_vect);
        %         pages = numel(t);
        %         R_EF_vect = permute(reshape(permute(R_EF_vect,[2,1]),[r/pages,c,pages]),[2,1,3]);
        %
        %         % Calculate position vector and related vectors/values in Rotational Frame
        %         R_R_vect = mmat(C_EF_to_R,R_EF_vect,[1 2]);
        %         R_R_vect = ipermute(squeeze(R_R_vect),[2,1]);
        %         R_R_norm = sqrt(R_R_vect(:,1).^2 + R_R_vect(:,2).^2 + R_R_vect(:,3).^2 );
        
        %%%%%%%%
        
        %         evolutions(k).trajectories(pp).individual.t = t(:);
        %         evolutions(k).trajectories(pp).individual.x_R = R_R_vect(:,1);
        %         evolutions(k).trajectories(pp).individual.y_R = R_R_vect(:,2);
        %         evolutions(k).trajectories(pp).individual.z_R = R_R_vect(:,3);
        %         evolutions(k).trajectories(pp).individual.R_R_norm = R_R_norm(:);
        
        
        
        
        
        evolutions(k).trajectories(p).individual.t                     = depvar(:,1);
        evolutions(k).trajectories(p).individual.time_vector           = depvar(:,1) - startEpoch;
        evolutions(k).trajectories(p).individual.x_R                   = depvar(:,2);
        evolutions(k).trajectories(p).individual.y_R                   = depvar(:,3);
        evolutions(k).trajectories(p).individual.z_R                   = depvar(:,4);
        evolutions(k).trajectories(p).individual.altitude              = depvar(:,5);
        evolutions(k).trajectories(p).individual.latitude_angle        = rad2deg(depvar(:,6));
        evolutions(k).trajectories(p).individual.longitude_angle       = rad2deg(depvar(:,7));
        %    evolutions(k).trajectories(pp).individual.latitude_angle  = rad2deg(depvar(:,8));
        %    evolutions(k).trajectories(pp).individual.longitude_angle = rad2deg(depvar(:,9));
        evolutions(k).trajectories(p).individual.heading_angle         = rad2deg(depvar(:,8));
        evolutions(k).trajectories(p).individual.flight_path_angle     = rad2deg(depvar(:,9));
        evolutions(k).trajectories(p).individual.angle_of_attack       = rad2deg(depvar(:,10));
        evolutions(k).trajectories(p).individual.angle_of_sideslip     = rad2deg(depvar(:,11));
        evolutions(k).trajectories(p).individual.bank_angle            = rad2deg(depvar(:,12));
        evolutions(k).trajectories(p).individual.height                = depvar(:,13);
        evolutions(k).trajectories(p).individual.mach                  = depvar(:,14);
        evolutions(k).trajectories(p).individual.airspeed              = depvar(:,15);
        evolutions(k).trajectories(p).individual.total_aero_g_load     = depvar(:,16);
        evolutions(k).trajectories(p).individual.acc_aero_x            = depvar(:,17);
        evolutions(k).trajectories(p).individual.acc_aero_y            = depvar(:,18);
        evolutions(k).trajectories(p).individual.acc_aero_z            = depvar(:,19);
        evolutions(k).trajectories(p).individual.acc_grav_x            = depvar(:,20);
        evolutions(k).trajectories(p).individual.acc_grav_y            = depvar(:,21);
        evolutions(k).trajectories(p).individual.acc_grav_z            = depvar(:,22);
        evolutions(k).trajectories(p).individual.acc_thru_x            = depvar(:,23);
        evolutions(k).trajectories(p).individual.acc_thru_y            = depvar(:,24);
        evolutions(k).trajectories(p).individual.acc_thru_z            = depvar(:,25);
        evolutions(k).trajectories(p).individual.acc_x                 = depvar(:,17) + depvar(:,20) + depvar(:,23);
        evolutions(k).trajectories(p).individual.acc_y                 = depvar(:,18) + depvar(:,21) + depvar(:,24);
        evolutions(k).trajectories(p).individual.acc_z                 = depvar(:,19) + depvar(:,22) + depvar(:,25);
        evolutions(k).trajectories(p).individual.acc_aero_M            = sqrt(depvar(:,17).^2 + depvar(:,18).^2 + depvar(:,19).^2);
        evolutions(k).trajectories(p).individual.acc_grav_M            = sqrt(depvar(:,20).^2 + depvar(:,21).^2 + depvar(:,22).^2);
        evolutions(k).trajectories(p).individual.acc_thru_M            = sqrt(depvar(:,23).^2 + depvar(:,24).^2 + depvar(:,25).^2);
        evolutions(k).trajectories(p).individual.dynamic_pressure      = depvar(:,26);
        evolutions(k).trajectories(p).individual.heating_rate          = depvar(:,27);
        %evolutions(k).trajectories(p).individual.stagnation_flux       = depvar(:,27);
        evolutions(k).trajectories(p).individual.mass_rate             = depvar(:,28);
        evolutions(k).trajectories(p).individual.mass                  = depvar(:,29);
        evolutions(k).trajectories(p).individual.E                     = depvar(:,30);
        evolutions(k).trajectories(p).individual.E_hat                 = depvar(:,31);
        evolutions(k).trajectories(p).individual.evaluated_throttle_setting       = depvar(:,32);
        evolutions(k).trajectories(p).individual.evaluated_thrust_elevation_angle = depvar(:,33);
        evolutions(k).trajectories(p).individual.evaluated_thrust_azimuth_angle = depvar(:,34);
        evolutions(k).trajectories(p).individual.evaluated_angle_of_attack  = depvar(:,35);
        evolutions(k).trajectories(p).individual.evaluated_bank_angle    = depvar(:,36);
        evolutions(k).trajectories(p).individual.engine_status           = depvar(:,37);
        evolutions(k).trajectories(p).individual.distance_traveled       = rad2deg(depvar(:,38));
        evolutions(k).trajectories(p).individual.distance_to_go          = rad2deg(depvar(:,39));
        evolutions(k).trajectories(p).individual.heading_to_target       = rad2deg(depvar(:,40));
        evolutions(k).trajectories(p).individual.heading_error           = rad2deg(depvar(:,41));
        evolutions(k).trajectories(p).individual.q_dot_LE                = (depvar(:,42);
        evolutions(k).trajectories(p).individual.body_thrust_x           = depvar(:,43);
        evolutions(k).trajectories(p).individual.body_thrust_y           = depvar(:,44);
        evolutions(k).trajectories(p).individual.body_thrust_z           = depvar(:,45);
        evolutions(k).trajectories(p).individual.bending_moment          = depvar(:,46);
        evolutions(k).trajectories(p).individual.local_gravity_1         = depvar(:,47);
        evolutions(k).trajectories(p).individual.local_gravity_2         = depvar(:,48);
        evolutions(k).trajectories(p).individual.eq_glide_limit          = depvar(:,49);
        evolutions(k).trajectories(p).individual.increment_Cm_bodyflap   = depvar(:,50);
        evolutions(k).trajectories(p).individual.bodyflap_deflection     = depvar(:,51);
        
        
        evolutions(k).trajectories(p).individual.interp_E_mapped_Ascent               = interp_Ascent(:,1);
        evolutions(k).trajectories(p).individual.interp_angle_of_attack_Ascent        = interp_Ascent(:,2);
        evolutions(k).trajectories(p).individual.interp_bank_angle_Ascent             = interp_Ascent(:,3);
        evolutions(k).trajectories(p).individual.interp_thrust_elevation_angle_Ascent = interp_Ascent(:,4);
        evolutions(k).trajectories(p).individual.interp_thrust_azimuth_angle_Ascent   = interp_Ascent(:,5);
        evolutions(k).trajectories(p).individual.interp_throttle_setting_Ascent       = interp_Ascent(:,6);
        
        evolutions(k).trajectories(p).individual.interp_E_mapped_Descent               = interp_Descent(:,1);
        evolutions(k).trajectories(p).individual.interp_angle_of_attack_Descent        = interp_Descent(:,2);
        evolutions(k).trajectories(p).individual.interp_bank_angle_Descent             = interp_Descent(:,3);
        evolutions(k).trajectories(p).individual.interp_thrust_elevation_angle_Descent = interp_Descent(:,4);
        evolutions(k).trajectories(p).individual.interp_thrust_azimuth_angle_Descent   = interp_Descent(:,5);
        evolutions(k).trajectories(p).individual.interp_throttle_setting_Descent       = interp_Descent(:,6);
        
        evolutions(k).trajectories(p).individual.DV_E_mapped_Ascent               = DV_mapped_Ascent(:,1);
        evolutions(k).trajectories(p).individual.DV_angle_of_attack_Ascent        = DV_mapped_Ascent(:,2);
        evolutions(k).trajectories(p).individual.DV_bank_angle_Ascent             = DV_mapped_Ascent(:,3);
        evolutions(k).trajectories(p).individual.DV_thrust_elevation_angle_Ascent = DV_mapped_Ascent(:,4);
        evolutions(k).trajectories(p).individual.DV_thrust_azimuth_angle_Ascent   = DV_mapped_Ascent(:,5);
        evolutions(k).trajectories(p).individual.DV_throttle_setting_Ascent       = DV_mapped_Ascent(:,6);
        evolutions(k).trajectories(p).individual.DV_node_location_Ascent          = DV_mapped_Ascent(:,7);

        
        evolutions(k).trajectories(p).individual.DV_E_mapped_Descent               = DV_mapped_Descent(:,1);
        evolutions(k).trajectories(p).individual.DV_angle_of_attack_Descent        = DV_mapped_Descent(:,2);
        evolutions(k).trajectories(p).individual.DV_bank_angle_Descent             = DV_mapped_Descent(:,3);
        evolutions(k).trajectories(p).individual.DV_thrust_elevation_angle_Descent = DV_mapped_Descent(:,4);
        evolutions(k).trajectories(p).individual.DV_thrust_azimuth_angle_Descent   = DV_mapped_Descent(:,5);
        evolutions(k).trajectories(p).individual.DV_throttle_setting_Descent       = DV_mapped_Descent(:,6);
        evolutions(k).trajectories(p).individual.DV_node_location_Descent          = DV_mapped_Descent(:,7);
        
        
        
        lat_c_rad = depvar(:,6); lon_c_rad = depvar(:,7);
        lat_f_rad = deg2rad(lat_f_deg); lon_f_rad =  deg2rad(lon_f_deg);
        lat_dif = lat_f_rad - lat_c_rad; lon_dif = lon_f_rad - lon_c_rad;
        distance_to_go = rad2deg(acos(sin(lat_c_rad).*sin(lat_f_rad) + cos(lat_c_rad).*cos(lat_f_rad).*cos(lon_dif)));
        chi_req_rad_Y = sin(lon_dif).*cos(lat_f_rad);
        chi_req_rad_X = cos(lat_c_rad)*sin(lat_f_rad) - sin(lat_c_rad).*cos(lat_f_rad).*cos(lon_dif);
        chi_req_deg = rad2deg(atan2(chi_req_rad_Y,chi_req_rad_X));
        chi_err_deg = rad2deg(depvar(:,8)) - chi_req_deg;
        
        %evolutions(k).trajectories(p).individual.E                = 9.80665*depvar(:,13) + 0.5*(depvar(:,15)).^2;
        %evolutions(k).trajectories(p).individual.d_deg            = d_deg;
        %evolutions(k).trajectories(p).individual.heading_error    = chi_err_deg;
        %evolutions(k).trajectories(p).individual.heading_required = chi_req_deg;
        
        evolutions(k).individuals.gamma_i(p)   = evolutions(k).trajectories(p).individual.flight_path_angle(1);
        evolutions(k).individuals.chi_i(p)     = evolutions(k).trajectories(p).individual.heading_angle(1);
        evolutions(k).individuals.lat_f_deg(p) = evolutions(k).trajectories(p).individual.latitude_angle(end);
        evolutions(k).individuals.lon_f_deg(p) = evolutions(k).trajectories(p).individual.longitude_angle(end);
        evolutions(k).individuals.tof(p)       = evolutions(k).trajectories(p).individual.time_vector(end);
        evolutions(k).individuals.interp_E_mapped_Ascent(p)  = evolutions(k).trajectories(p).individual.interp_E_mapped_Ascent(end);
        evolutions(k).individuals.interp_E_mapped_Descent(p)  = evolutions(k).trajectories(p).individual.interp_E_mapped_Descent(end);
        
        evolutions(k).fitness.dif_lat(p)   = evolutions(k).trajectories(p).individual.latitude_angle(end) - lat_f_deg;
        evolutions(k).fitness.dif_lon(p)   = evolutions(k).trajectories(p).individual.longitude_angle(end) - lon_f_deg;
        evolutions(k).fitness.dif_d_deg(p) = evolutions(k).trajectories(p).individual.distance_to_go(end) - 0.75;
        evolutions(k).fitness.dif_h(p)     = evolutions(k).trajectories(p).individual.height(end) - 25000;
        evolutions(k).fitness.tof(p)       = evolutions(k).trajectories(p).individual.time_vector(end);
        
        
        pp = pp + 1;
    end
    [ aaa , idx1 ] = min(evolutions(k).fitness.dif_d_deg);
    [ bbb , idx2 ] = min(evolutions(k).fitness.dif_h);
    [ ccc , idx3 ] = min(evolutions(k).fitness.tof);
    I = [ idx1 idx2 idx3 ];
    criteria = [{'d_deg'} {'h'} {'tof'} ];
    %  I_12 = intersect(idx1,idx2);
    %  I_123 = intersect(I_12,idx3);
    %  I_123 = idx1;
    evolutions(k).max_tof       = max(evolutions(k).individuals.tof);
    evolutions(k).max_interp_E_mapped_Ascent  = max(evolutions(k).individuals.interp_E_mapped_Ascent);
    evolutions(k).max_interp_E_mapped_Descent  = max(evolutions(k).individuals.interp_E_mapped_Descent);
    
    for i = 1:3
        evolutions(k).best(i).criteria  = criteria(i);
        evolutions(k).best(i).index     = I(i);
        % evolutions(k).best(i).v_i       = evolutions(k).individuals.v_i(I(i));
        evolutions(k).best(i).gamma_i   = evolutions(k).individuals.gamma_i(I(i));
        evolutions(k).best(i).chi_i     = evolutions(k).individuals.chi_i(I(i));
        evolutions(k).best(i).dif_d_deg = evolutions(k).fitness.dif_d_deg(I(i));
        evolutions(k).best(i).dif_h     = evolutions(k).fitness.dif_h(I(i));
        evolutions(k).best(i).tof       = evolutions(k).fitness.tof(I(i));
    end
    
end
end

