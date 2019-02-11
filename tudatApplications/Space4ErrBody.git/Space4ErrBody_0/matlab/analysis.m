
clearvars
close all
clc
dbstop if error

format long
%hold on
%[S_x, S_y, S_z] = sphere(100);
%surf(S_x*r, S_y*r, S_z*r);`
%p = plotearth('Maptype','bluemarble','NEOMap',...
%    '/Users/bislip/Documents/MATLAB/PlotEarth/PlotEarth/BlueMarble.png',...
%    'SampleStep',1,'FullColor',true,'Shape','spherical');

%%
online = 0;

mainpath = '/Users/bislip/Cloud Storage/OneDrive/School/TUDelft/Space Flight/Thesis/tudatBundle.git/tudatApplications/Space4ErrBody.git/Space4ErrBody_0/matlab/';

if online == 1
    mainpath = '..';
end
pop_Location = strcat(mainpath,'SimulationOutput/');
pop_prefix = 'population*';

fit_Location = pop_Location;
fit_prefix = 'fitness*';

Output_Location = pop_Location;
%Folder_prefix = 'HORUS_OUTPUT*';
Folder_prefix = 'OUTPUT*';
%File_prefix = 'HORUSsystemFinalStateGOAL*';
prop_File_prefix = 'HORUSPropHistory*';
depvar_File_prefix = 'HORUSDepVar*';
interp_Ascent_File_prefix = 'interpolators_Ascent*';
interp_Descent_File_prefix = 'interpolators_Descent*';
DV_mapped_Ascent_File_prefix = 'map_DV_mapped_Ascent*';
DV_mapped_Descent_File_prefix = 'map_DV_mapped_Descent*';

% Hardcoded coordinates of origin: AMS
lon_i_deg = 4.76416667;
lon_i_rad = deg2rad(lon_i_deg);
% Hardcoded coordinates of final destination: IAD
lat_f_deg = 38.9444444444444;
lon_f_deg = -77.4558333333;
validation = 0;

% Hardcoded Starting Epoch
startEpoch = 2458419.95833333;

% Coordinates for Validation case: Re-entry towards Kourou
if contains(Folder_prefix,'HORUS_VALIDATION') == 1
    
    lon_i_deg = -106.7;
    lon_i_rad = deg2rad(lon_i_deg);
    lat_f_deg = 5.;
    lon_f_deg = -53;
    validation = 1;
    
end

%% Construct file prefix
[Folder_Path_List,prop_File_Path_List_prefix,...
    depvar_File_Path_List_prefix,interp_Ascent_File_Path_List_prefix,interp_Decent_File_Path_List_prefix,DV_mapped_Ascent_File_Path_List_prefix,DV_mapped_Descent_File_Path_List_prefix,pop_file_path_prefix,...
    fit_file_path_prefix] = Contruct_File_prefix(...
    Output_Location,...
    Folder_prefix,...
    prop_File_prefix,...
    depvar_File_prefix,...
    interp_Ascent_File_prefix,...
    interp_Descent_File_prefix,...
    DV_mapped_Ascent_File_prefix,...
    DV_mapped_Descent_File_prefix,...
    pop_Location,...
    pop_prefix,...
    fit_Location,...
    fit_prefix);

%% Scan directories

% option = 1   Try to plot, indiscriminately
% option = 2   Analyze
% option = 3   Analyze and plot

option = 1;


switch option
    
    case 1
        
        textprogressbar('Loading Data.     ')
        for p = 1:numel(prop_File_Path_List_prefix)
            
            % Prepare what's needed for the analysis
            [evolutions,prop_path,depvar_path,interp_Ascent_path,interp_Descent_path,DV_mapped_Ascent_path,DV_mapped_Descent_path,tof,v_i,gamma_i,chi_i,lat_f,lon_f,pop_path,...
                pop_i,fit_path,fit_i,output] = Path_Prep(option,p,...
                prop_File_Path_List_prefix,depvar_File_Path_List_prefix,interp_Ascent_File_Path_List_prefix,interp_Decent_File_Path_List_prefix,DV_mapped_Ascent_File_Path_List_prefix,...
                DV_mapped_Descent_File_Path_List_prefix,Folder_Path_List,pop_file_path_prefix,fit_file_path_prefix);
            
            compilation(p).set = char(strcat(extractBetween(output,'evolutions','.mat')));
            compilation(p).evolutions = evolutions;
            compilation(p).validation = validation;
            
            % textprogressbar(p*100/numel(File_Path_List_prefix))
            
        end
        
        textprogressbar('    Done. Now plotting/printing.')
        plotsomestuff( compilation , mainpath )
        
    case 2
        
        textprogressbar('Analyzing Simulations.     ')
        for p = 1:numel(prop_File_Path_List_prefix)
            
            % Prepare what's needed for the analysis
            [evolutions,prop_path,depvar_path,interp_Ascent_path,interp_Descent_path,DV_mapped_Ascent_path,DV_mapped_Descent_path,tof,v_i,gamma_i,chi_i,lat_f,lon_f,pop_path,...
                pop_i,fit_path,fit_i,output] = Path_Prep(option,p,...
                prop_File_Path_List_prefix,depvar_File_Path_List_prefix,interp_Ascent_File_Path_List_prefix,interp_Decent_File_Path_List_prefix,DV_mapped_Ascent_File_Path_List_prefix,...
                DV_mapped_Descent_File_Path_List_prefix,Folder_Path_List,pop_file_path_prefix,fit_file_path_prefix);
            
            %             if isempty(idx) == 1
            %
            %                textprogressbar(p*100/numel(File_Path_List_prefix))
            %                continue
            %
            %             end
            
            % Append poplution
            [ evolutions ] = Analyze_Evolution(evolutions,v_i,gamma_i,chi_i,lat_f,lon_f,prop_path,pop_path,pop_i,fit_path,fit_i,lat_f_deg,lon_f_deg);
            
            % Analyze simulations
            [ evolutions ] = Get_Trajectories(evolutions,prop_path,depvar_path,interp_Ascent_path,interp_Descent_path,DV_mapped_Ascent_path,DV_mapped_Descent_path,v_i,gamma_i,pop_i,lon_i_rad,lat_f_deg,lon_f_deg,startEpoch);
            
            save(output,'evolutions');
            %            textprogressbar(p*100/numel(File_Path_List_prefix))
            
        end
        
        % textprogressbar('    Done.')
        
    case 3
        
        %textprogressbar('Analyzing Simulations.     ')
        for p = 1:numel(prop_File_Path_List_prefix)
            
            % Prepare what's needed for the analysis
            [evolutions,prop_path,depvar_path,interp_Ascent_path,interp_Descent_path,DV_mapped_Ascent_path,DV_mapped_Descent_path,tof,v_i,gamma_i,chi_i,lat_f,lon_f,pop_path,...
                pop_i,fit_path,fit_i,output] = Path_Prep(option,p,...
                prop_File_Path_List_prefix,depvar_File_Path_List_prefix,interp_Ascent_File_Path_List_prefix,interp_Decent_File_Path_List_prefix,DV_mapped_Ascent_File_Path_List_prefix,...
                DV_mapped_Descent_File_Path_List_prefix,Folder_Path_List,pop_file_path_prefix,fit_file_path_prefix);
            
            %            if isempty(idx) == 1
            %
            %                compilation(p).set = char(strcat(extractBetween(output,'analyzed_simulations_','.mat')));
            %                compilation(p).analyzed_simulations = analyzed_simulations;
            %
            %                textprogressbar(p*100/numel(File_Path_List_prefix))
            %                continue
            %
            %            end
            
            % Append poplution
            [ evolutions ] = Analyze_Evolution(evolutions,v_i,gamma_i,chi_i,lat_f,lon_f,prop_path,pop_path,pop_i,fit_path,fit_i,lat_f_deg,lon_f_deg);
            
            % Analyze simulations
            [ evolutions ] = Get_Trajectories(evolutions,prop_path,depvar_path,interp_Ascent_path,interp_Descent_path,DV_mapped_Ascent_path,DV_mapped_Descent_path,v_i,gamma_i,pop_i,lon_i_rad,lat_f_deg,lon_f_deg,startEpoch);
            
            save(output, 'evolutions');
            compilation(p).set = char(strcat(extractBetween(output,'evolutions','.mat')));
            compilation(p).evolutions = evolutions;
            compilation(p).validation = validation;
            
            % textprogressbar(p*100/numel(prop_File_Path_List_prefix))
            
        end
        
        textprogressbar('    Done. Now plotting/printing.')
        plotsomestuff( compilation , mainpath )
        
end








'here'