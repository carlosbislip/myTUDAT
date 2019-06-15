function [evolutions,output] = Load_Data(option,p,Folder_Path_List,pop_path,npop_files)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
disp('LOAD_DATA')



%% Load lists of analyzed simulations

switch option
    
    case 1
        
        % Figure out output name based on what is in folder path
        if contains(convertCharsToStrings(Folder_Path_List{p}),'OUTPUT') == 1
            output = char(strcat('evolutions_',extractAfter(Folder_Path_List{p},'T_'),'.mat'));
        elseif contains(convertCharsToStrings(Folder_Path_List{p}),'VALIDATION') == 1
            output = char(strcat('evolutions_',extractAfter(Folder_Path_List{p},'N_'),'.mat'));
        end
        
        % Open file if it exists. Create variable otherwise.
        if exist(output,'file') == 2
            load(output)
        else
            disp('Run again with option 2 or 3. No .mat file found.');
            
            %             k = npop_files + 1;
            %             evolutions(k).evolution(npop_files)    = nan;
            %             evolutions(k).v_i(npop_files)          = nan;
            %             evolutions(k).gamma_i(npop_files)      = nan;
            %             evolutions(k).chi_i(npop_files)        = nan;
            %             evolutions(k).trajectories(npop_files) = nan;
            %             evolutions(k).fitness(npop_files)      = nan;
            %             evolutions(k).lat(npop_files)          = nan;
            %             evolutions(k).lon(npop_files)          = nan;
            %             evolutions(k).tof(npop_files)          = nan;
        end
        
        idx = [];
        
    case {2, 3}
        
        % Figure out output name based on what is in folder path
        if contains(convertCharsToStrings(Folder_Path_List{p}),'OUTPUT') == 1
            output = char(strcat('evolutions_',extractAfter(Folder_Path_List{p},'T_'),'.mat'));
        elseif contains(convertCharsToStrings(Folder_Path_List{p}),'VALIDATION') == 1
            output = char(strcat('evolutions_',extractAfter(Folder_Path_List{p},'N_'),'.mat'));
        end
        
        % Open file if it exists. Create variable otherwise.
        %if exist(output,'file') == 2
        %    load(output)
       % else
            %if contains(pop_path, 'monteCarlo') == true
                k = npop_files;
            %else
             %   k = npop_files + 1;
            %end
            evolutions(k).evolution                                         = nan;
            
            evolutions(k).population.nodeIntervalAscent1                    = nan;
            evolutions(k).population.nodeIntervalAscent2                    = nan;
            evolutions(k).population.nodeIntervalAscent3                    = nan;
            evolutions(k).population.nodeIntervalAscent4                    = nan;
            evolutions(k).population.nodeIntervalAscent5                    = nan;
            evolutions(k).population.nodeIntervalAscent6                    = nan;
            evolutions(k).population.nodeIntervalAscent7                    = nan;
            evolutions(k).population.nodeIntervalAscent8                    = nan;
            evolutions(k).population.nodeIntervalAscent9                    = nan;
            
            evolutions(k).population.angleOfAttackAscent1                   = nan;
            evolutions(k).population.angleOfAttackAscent2                   = nan;
            evolutions(k).population.angleOfAttackAscent3                   = nan;
            evolutions(k).population.angleOfAttackAscent4                   = nan;
            evolutions(k).population.angleOfAttackAscent5                   = nan;
            evolutions(k).population.angleOfAttackAscent6                   = nan;
            evolutions(k).population.angleOfAttackAscent7                   = nan;
            evolutions(k).population.angleOfAttackAscent8                   = nan;
            evolutions(k).population.angleOfAttackAscent9                   = nan;
            evolutions(k).population.angleOfAttackAscent10                  = nan;
            
            evolutions(k).population.bankAngleAscent1                       = nan;
            evolutions(k).population.bankAngleAscent2                       = nan;
            evolutions(k).population.bankAngleAscent3                       = nan;
            evolutions(k).population.bankAngleAscent4                       = nan;
            evolutions(k).population.bankAngleAscent5                       = nan;
            evolutions(k).population.bankAngleAscent6                       = nan;
            evolutions(k).population.bankAngleAscent7                       = nan;
            evolutions(k).population.bankAngleAscent8                       = nan;
            evolutions(k).population.bankAngleAscent9                       = nan;
            evolutions(k).population.bankAngleAscent10                      = nan;
            
            evolutions(k).population.thrustElevationAngleAscent1            = nan;
            evolutions(k).population.thrustElevationAngleAscent2            = nan;
            evolutions(k).population.thrustElevationAngleAscent3            = nan;
            evolutions(k).population.thrustElevationAngleAscent4            = nan;
            evolutions(k).population.thrustElevationAngleAscent5            = nan;
            evolutions(k).population.thrustElevationAngleAscent6            = nan;
            evolutions(k).population.thrustElevationAngleAscent7            = nan;
            evolutions(k).population.thrustElevationAngleAscent8            = nan;
            evolutions(k).population.thrustElevationAngleAscent9            = nan;
            evolutions(k).population.thrustElevationAngleAscent10           = nan;
            
            evolutions(k).population.thrustAzimuthAngleAscent1              = nan;
            evolutions(k).population.thrustAzimuthAngleAscent2              = nan;
            evolutions(k).population.thrustAzimuthAngleAscent3              = nan;
            evolutions(k).population.thrustAzimuthAngleAscent4              = nan;
            evolutions(k).population.thrustAzimuthAngleAscent5              = nan;
            evolutions(k).population.thrustAzimuthAngleAscent6              = nan;
            evolutions(k).population.thrustAzimuthAngleAscent7              = nan;
            evolutions(k).population.thrustAzimuthAngleAscent8              = nan;
            evolutions(k).population.thrustAzimuthAngleAscent9              = nan;
            evolutions(k).population.thrustAzimuthAngleAscent10             = nan;
            
            evolutions(k).population.throttleSettingAscent1                 = nan;
            evolutions(k).population.throttleSettingAscent2                 = nan;
            evolutions(k).population.throttleSettingAscent3                 = nan;
            evolutions(k).population.throttleSettingAscent4                 = nan;
            evolutions(k).population.throttleSettingAscent5                 = nan;
            evolutions(k).population.throttleSettingAscent6                 = nan;
            evolutions(k).population.throttleSettingAscent7                 = nan;
            evolutions(k).population.throttleSettingAscent8                 = nan;
            evolutions(k).population.throttleSettingAscent9                 = nan;
            evolutions(k).population.throttleSettingAscent10                = nan;
            evolutions(k).population.initialLaunchHeading                   = nan;
            evolutions(k).population.initialVelocity                        = nan;
            evolutions(k).population.maximumVelocity                        = nan;
            evolutions(k).population.maximumHeight                          = nan;
            evolutions(k).population.additionalMass                         = nan;
            evolutions(k).population.terminationDistanceRatio               = nan;
            
            evolutions(k).population.nodeIntervalDescent1                   = nan;
            evolutions(k).population.nodeIntervalDescent2                   = nan;
            evolutions(k).population.nodeIntervalDescent3                   = nan;
            evolutions(k).population.nodeIntervalDescent4                   = nan;
            evolutions(k).population.nodeIntervalDescent5                   = nan;
            evolutions(k).population.nodeIntervalDescent6                   = nan;
            evolutions(k).population.nodeIntervalDescent7                   = nan;
            evolutions(k).population.nodeIntervalDescent8                   = nan;
            evolutions(k).population.nodeIntervalDescent9                   = nan;
            
            evolutions(k).population.angleOfAttackDescent1                  = nan;
            evolutions(k).population.angleOfAttackDescent2                  = nan;
            evolutions(k).population.angleOfAttackDescent3                  = nan;
            evolutions(k).population.angleOfAttackDescent4                  = nan;
            evolutions(k).population.angleOfAttackDescent5                  = nan;
            evolutions(k).population.angleOfAttackDescent6                  = nan;
            evolutions(k).population.angleOfAttackDescent7                  = nan;
            evolutions(k).population.angleOfAttackDescent8                  = nan;
            evolutions(k).population.angleOfAttackDescent9                  = nan;
            evolutions(k).population.angleOfAttackDescent10                 = nan;
            
            evolutions(k).population.bankAngleDescent1                      = nan;
            evolutions(k).population.bankAngleDescent2                      = nan;
            evolutions(k).population.bankAngleDescent3                      = nan;
            evolutions(k).population.bankAngleDescent4                      = nan;
            evolutions(k).population.bankAngleDescent5                      = nan;
            evolutions(k).population.bankAngleDescent6                      = nan;
            evolutions(k).population.bankAngleDescent7                      = nan;
            evolutions(k).population.bankAngleDescent8                      = nan;
            evolutions(k).population.bankAngleDescent9                      = nan;
            evolutions(k).population.bankAngleDescent10                     = nan;
            
            evolutions(k).population.thrustElevationAngleDescent1           = nan;
            evolutions(k).population.thrustElevationAngleDescent2           = nan;
            evolutions(k).population.thrustElevationAngleDescent3           = nan;
            evolutions(k).population.thrustElevationAngleDescent4           = nan;
            evolutions(k).population.thrustElevationAngleDescent5           = nan;
            evolutions(k).population.thrustElevationAngleDescent6           = nan;
            evolutions(k).population.thrustElevationAngleDescent7           = nan;
            evolutions(k).population.thrustElevationAngleDescent8           = nan;
            evolutions(k).population.thrustElevationAngleDescent9           = nan;
            evolutions(k).population.thrustElevationAngleDescent10          = nan;
            
            evolutions(k).population.thrustAzimuthAngleDescent1             = nan;
            evolutions(k).population.thrustAzimuthAngleDescent2             = nan;
            evolutions(k).population.thrustAzimuthAngleDescent3             = nan;
            evolutions(k).population.thrustAzimuthAngleDescent4             = nan;
            evolutions(k).population.thrustAzimuthAngleDescent5             = nan;
            evolutions(k).population.thrustAzimuthAngleDescent6             = nan;
            evolutions(k).population.thrustAzimuthAngleDescent7             = nan;
            evolutions(k).population.thrustAzimuthAngleDescent8             = nan;
            evolutions(k).population.thrustAzimuthAngleDescent9             = nan;
            evolutions(k).population.thrustAzimuthAngleDescent10            = nan;
            
            evolutions(k).population.throttleSettingDescent1                = nan;
            evolutions(k).population.throttleSettingDescent2                = nan;
            evolutions(k).population.throttleSettingDescent3                = nan;
            evolutions(k).population.throttleSettingDescent4                = nan;
            evolutions(k).population.throttleSettingDescent5                = nan;
            evolutions(k).population.throttleSettingDescent6                = nan;
            evolutions(k).population.throttleSettingDescent7                = nan;
            evolutions(k).population.throttleSettingDescent8                = nan;
            evolutions(k).population.throttleSettingDescent9                = nan;
            evolutions(k).population.throttleSettingDescent10               = nan;
            evolutions(k).population.finalVelocity                          = nan;
            evolutions(k).population.skipSuppressionTriggerTime             = nan;
            
            
            evolutions(k).trajectories.individual.t                                = nan;
            evolutions(k).trajectories.individual.t                                = nan;
            evolutions(k).trajectories.individual.x_R                              = nan;
            evolutions(k).trajectories.individual.y_R                              = nan;
            evolutions(k).trajectories.individual.z_R                              = nan;
            evolutions(k).trajectories.individual.lat                              = nan;
            evolutions(k).trajectories.individual.lon                              = nan;
            evolutions(k).trajectories.individual.latitude_angle                   = nan;
            evolutions(k).trajectories.individual.longitude_angle                  = nan;
            evolutions(k).trajectories.individual.heading_angle                    = nan;
            evolutions(k).trajectories.individual.heading_required                 = nan;
            evolutions(k).trajectories.individual.heading_error                    = nan;
            evolutions(k).trajectories.individual.flight_path_angle                = nan;
            evolutions(k).trajectories.individual.angle_of_attack                  = nan;
            evolutions(k).trajectories.individual.angle_of_sideslip                = nan;
            evolutions(k).trajectories.individual.bank_angle                       = nan;
            evolutions(k).trajectories.individual.airspeed                         = nan;
            evolutions(k).trajectories.individual.altitude                         = nan;
            evolutions(k).trajectories.individual.height                           = nan;
            evolutions(k).trajectories.individual.mach                             = nan;
            evolutions(k).trajectories.individual.specificEnergy                   = nan;
            evolutions(k).trajectories.individual.normalizedSpecificEnergy         = nan;
            evolutions(k).trajectories.individual.total_aero_g_load                = nan;
            evolutions(k).trajectories.individual.dynamicPressure                 = nan;
            evolutions(k).trajectories.individual.heat_rate_TUDAT_nose             = nan;
            evolutions(k).trajectories.individual.mass                             = nan;
            
            evolutions(k).trajectories.individual.mass_rate                        = nan;
            evolutions(k).trajectories.individual.evaluated_throttle_setting       = nan;
            evolutions(k).trajectories.individual.evaluated_thrust_elevation_angle = nan;
            evolutions(k).trajectories.individual.evaluated_thrust_azimuth_angle   = nan;
            evolutions(k).trajectories.individual.evaluated_angle_of_attack        = nan;
            evolutions(k).trajectories.individual.evaluated_bank_angle             = nan;
            
            evolutions(k).trajectories.individual.commanded_throttle_setting       = nan;
            evolutions(k).trajectories.individual.commanded_thrust_elevation_angle = nan;
            evolutions(k).trajectories.individual.commanded_thrust_azimuth_angle   = nan;
            evolutions(k).trajectories.individual.commanded_angle_of_attack        = nan;
            evolutions(k).trajectories.individual.commanded_bank_angle             = nan;
            evolutions(k).trajectories.individual.engine_status                    = nan;
            evolutions(k).trajectories.individual.distance_traveled                = nan;
            evolutions(k).trajectories.individual.distance_to_go                   = nan;
            evolutions(k).trajectories.individual.heading_to_target                = nan;
            evolutions(k).trajectories.individual.heading_error                    = nan;
            evolutions(k).trajectories.individual.heat_flux_tauber_leadingedge     = nan;
            evolutions(k).trajectories.individual.body_fixed_thrust_load_x         = nan;
            evolutions(k).trajectories.individual.body_fixed_thrust_load_y         = nan;
            evolutions(k).trajectories.individual.body_fixed_thrust_load_z         = nan;
            evolutions(k).trajectories.individual.bending_moment                   = nan;
            evolutions(k).trajectories.individual.local_gravity_1                  = nan;
            evolutions(k).trajectories.individual.local_gravity_2                  = nan;
            evolutions(k).trajectories.individual.skip_suppression_limit           = nan;
            evolutions(k).trajectories.individual.bodyflap_deflection              = nan;
            evolutions(k).trajectories.individual.increment_Cm_bodyflap            = nan;
            evolutions(k).trajectories.individual.increment_Cm_bodyflap_dif        = nan;
            evolutions(k).trajectories.individual.body_fixed_total_load_x          = nan;
            evolutions(k).trajectories.individual.body_fixed_total_load_y          = nan;
            evolutions(k).trajectories.individual.body_fixed_total_load_z          = nan;
            evolutions(k).trajectories.individual.body_fixed_total_g_load_x        = nan;
            evolutions(k).trajectories.individual.body_fixed_total_g_load_y        = nan;
            evolutions(k).trajectories.individual.body_fixed_total_g_load_z        = nan;
            evolutions(k).trajectories.individual.body_fixed_total_g_load_mag      = nan;
            evolutions(k).trajectories.individual.body_fixed_aero_load_x           = nan;
            evolutions(k).trajectories.individual.body_fixed_aero_load_y           = nan;
            evolutions(k).trajectories.individual.body_fixed_aero_load_z           = nan;
            evolutions(k).trajectories.individual.aero_force_coefficient_C_D       = nan;
            evolutions(k).trajectories.individual.aero_force_coefficient_C_S       = nan;
            evolutions(k).trajectories.individual.aero_force_coefficient_C_L       = nan;
            evolutions(k).trajectories.individual.aero_moment_coefficient_C_l      = nan;
            evolutions(k).trajectories.individual.aero_moment_coefficient_C_m      = nan;
            evolutions(k).trajectories.individual.aero_moment_coefficient_C_n      = nan;
            evolutions(k).trajectories.individual.heat_flux_chapman_nose           = nan;
            evolutions(k).trajectories.individual.localDensity                     = nan;
            evolutions(k).trajectories.individual.passenger_fixed_total_g_load_x   = nan;
            evolutions(k).trajectories.individual.passenger_fixed_total_g_load_y   = nan;
            evolutions(k).trajectories.individual.passenger_fixed_total_g_load_z   = nan;
            evolutions(k).trajectories.individual.currentLiftForce                 = nan;
            evolutions(k).trajectories.individual.headingErrorDeadband             = nan;
            evolutions(k).trajectories.individual.tempBankAngle                    = nan;
            evolutions(k).trajectories.individual.reversalConditional              = nan;
            evolutions(k).trajectories.individual.bank_angle_reversal_trigger      = nan;
            evolutions(k).trajectories.individual.wall_temperature_chapman             = nan;
            evolutions(k).trajectories.individual.wall_temperature_tauber_stagnation   = nan;
            evolutions(k).trajectories.individual.wall_temperature_tauber_flatplate    = nan;
            evolutions(k).trajectories.individual.heat_flux_tauber_stagnation          = nan;
            evolutions(k).trajectories.individual.heat_flux_tauber_flatplate           = nan;
            evolutions(k).trajectories.individual.groundtrack_covered                  = nan;
            evolutions(k).trajectories.individual.groundtrack_difference               = nan;
            evolutions(k).trajectories.individual.time_of_flight                       = nan;
            evolutions(k).trajectories.individual.flight_path_angle_rate               = nan;
            evolutions(k).trajectories.individual.cumulative_distance_travelled        = nan;
            evolutions(k).trajectories.individual.elevon_deflection                    = nan;
            evolutions(k).trajectories.individual.thrustMagnitude                      = nan;
            evolutions(k).trajectories.individual.speedOfSound                         = nan;
            evolutions(k).trajectories.individual.adiabaticWallTemperature             = nan;
            evolutions(k).trajectories.individual.freestreamTemperature                = nan;
            evolutions(k).trajectories.individual.currentDragForce                     = nan;
            evolutions(k).trajectories.individual.estimatedFightPathAngle              = nan;
            evolutions(k).trajectories.individual.aerodyamicFrameAerodynamicLoad_x     = nan;
            evolutions(k).trajectories.individual.aerodyamicFrameAerodynamicLoad_y     = nan;
            evolutions(k).trajectories.individual.aerodyamicFrameAerodynamicLoad_z     = nan;
            evolutions(k).trajectories.individual.aerodyamicFrameTotalLoad_x           = nan;
            evolutions(k).trajectories.individual.aerodyamicFrameTotalLoad_y           = nan;
            evolutions(k).trajectories.individual.aerodyamicFrameTotalLoad_z           = nan;

            
            
            
            evolutions(k).trajectories.individual.acc_aero_x        = nan;
            evolutions(k).trajectories.individual.acc_aero_y        = nan;
            evolutions(k).trajectories.individual.acc_aero_z        = nan;
            evolutions(k).trajectories.individual.acc_grav_x        = nan;
            evolutions(k).trajectories.individual.acc_grav_y        = nan;
            evolutions(k).trajectories.individual.acc_grav_z        = nan;
            evolutions(k).trajectories.individual.acc_thru_x        = nan;
            evolutions(k).trajectories.individual.acc_thru_y        = nan;
            evolutions(k).trajectories.individual.acc_thru_z        = nan;
            evolutions(k).trajectories.individual.acc_x             = nan;
            evolutions(k).trajectories.individual.acc_y             = nan;
            evolutions(k).trajectories.individual.acc_z             = nan;
            evolutions(k).trajectories.individual.acc_aero_M        = nan;
            evolutions(k).trajectories.individual.acc_grav_M        = nan;
            evolutions(k).trajectories.individual.acc_thru_M        = nan;
            evolutions(k).trajectories.individual.interp_E_Ascent                      = nan;
            evolutions(k).trajectories.individual.interp_angle_of_attack_Ascent        = nan;
            evolutions(k).trajectories.individual.interp_bank_angle_Ascent             = nan;
            evolutions(k).trajectories.individual.interp_thrust_elevation_angle_Ascent = nan;
            evolutions(k).trajectories.individual.interp_thrust_azimuth_angle_Ascent   = nan;
            evolutions(k).trajectories.individual.interp_throttle_setting_Ascent       = nan;
            
            evolutions(k).trajectories.individual.DV_E_mapped_Ascent               = nan;
            evolutions(k).trajectories.individual.DV_angle_of_attack_Ascent        = nan;
            evolutions(k).trajectories.individual.DV_bank_angle_Ascent             = nan;
            evolutions(k).trajectories.individual.DV_thrust_elevation_angle_Ascent = nan;
            evolutions(k).trajectories.individual.DV_thrust_azimuth_angle_Ascent   = nan;
            evolutions(k).trajectories.individual.DV_throttle_setting_Ascent       = nan;
            
            evolutions(k).trajectories.individual.interp_E_Descent                      = nan;
            evolutions(k).trajectories.individual.interp_angle_of_attack_Descent        = nan;
            evolutions(k).trajectories.individual.interp_bank_angle_Descent             = nan;
            evolutions(k).trajectories.individual.interp_thrust_elevation_angle_Descent = nan;
            evolutions(k).trajectories.individual.interp_thrust_azimuth_angle_Descent   = nan;
            evolutions(k).trajectories.individual.interp_throttle_setting_Descent = nan;
            
            evolutions(k).trajectories.individual.DV_E_mapped_Descent               = nan;
            evolutions(k).trajectories.individual.DV_angle_of_attack_Descent        = nan;
            evolutions(k).trajectories.individual.DV_bank_angle_Descent             = nan;
            evolutions(k).trajectories.individual.DV_thrust_elevation_angle_Descent = nan;
            evolutions(k).trajectories.individual.DV_thrust_azimuth_angle_Descent   = nan;
            evolutions(k).trajectories.individual.DV_throttle_setting_Descent       = nan;
            
            
            
            evolutions(k).trajectories.individual.headingErrorDeadBand_distance    = nan;
            evolutions(k).trajectories.individual.headingErrorDeadBand_LB          = nan;
            evolutions(k).trajectories.individual.headingErrorDeadBand_UP          = nan;
            
            
            evolutions(k).hardConstraints             = nan;
            evolutions(k).hardConstraints             = nan;
            evolutions(k).max_tof                     = nan;
            evolutions(k).max_interp_E_mapped_Ascent  = nan;
            evolutions(k).max_interp_E_mapped_Descent = nan;
            evolutions(k).population.individual       = nan;
            evolutions(k).individuals.v_i             = nan;
            evolutions(k).individuals.gamma_i         = nan;
            evolutions(k).individuals.chi_i           = nan;
            evolutions(k).fitness.dif_norm            = nan;
            evolutions(k).fitness.dif_lat             = nan;
            evolutions(k).fitness.dif_lon             = nan;
            evolutions(k).fitness.dif_d_deg           = nan;
            evolutions(k).fitness.dif_h               = nan;
            evolutions(k).fitness.tof                 = nan;
            evolutions(k).best(3).criteria            = nan;
            evolutions(k).best(3).v_i                 = nan;
            evolutions(k).best(3).gamma_i             = nan;
            evolutions(k).best(3).chi_i               = nan;
            evolutions(k).best(3).dif_d_deg           = nan;
            evolutions(k).best(3).dif_h               = nan;
            evolutions(k).best(3).tof                 = nan;
            evolutions(k).best(3).max_tof             = nan;
            
        end
        
        %         idx = nan(1,npop_files);
        %         % Anaysis loop
        %         parfor k = 1:npop_files
        %
        %             % Compare each string to the lists of analyzed simulations. There must
        %             % be a way of comparing entire arrays in one go without looping. Then
        %             % just analyzing the remaining indices. This woudl be done before
        %             % entering this loop. Maybe use:    ~any(X)
        %             discriminant = find(strcmp(pop_path(k,:),{analyzed_simulations.path}) == 1);
        %
        %             if isempty(discriminant) == 1
        %                 idx(1,k) = k;
        %             end
        %
        %         end
        %
        %         idx = idx(~isnan(idx));
        
        
end

