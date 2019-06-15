function [ prop_path,...
    depvar_path,...
    interp_Ascent_path,...
    interp_Descent_path,...
    DV_mapped_Ascent_path,...
    DV_mapped_Descent_path,...
    headingErrorDeadbandBounds_path,...
    alphaMachBounds_path,...
    tof,v_i,gamma_i,chi_i,lat_f,lon_f,...
    pop_path,pop_i,fit_path,fit_i] = ...
    ...
    ...
    Contruct_File_Paths(...
    prop_Output_files,...
    n_prop_Output_files,...
    depvar_Output_files,...
    n_depvar_Output_files,...
    interp_Ascent_Output_files,...
    interp_Descent_Output_files,...
    DV_mapped_Ascent_Output_files,...
    DV_mapped_Descent_Output_files,...
    headingErrorDeadbandBounds_Output_files,...
    alphaMachBounds_Output_files,...
    pop_files,...
    npop_files,...
    fit_files,...
    nfit_files)

%CONSTRUCT_FILE_PATHS Summary of this function goes here
%   Detailed explanation goes here
disp('CONSTRUCT_FILE_PATHS')


prop_path{n_prop_Output_files,1} = [];
depvar_path{n_prop_Output_files,1} = [];
interp_Ascent_path{n_prop_Output_files,1} = [];
interp_Descent_path{n_prop_Output_files,1} = [];
DV_mapped_Ascent_path{n_prop_Output_files,1} = [];
DV_mapped_Descent_path{n_prop_Output_files,1} = [];
headingErrorDeadbandBounds_path{n_prop_Output_files,1} = [];
alphaMachBounds_path{n_prop_Output_files,1} = [];
v_i = nan(n_prop_Output_files,1);
gamma_i = nan(n_prop_Output_files,1);
chi_i = nan(n_prop_Output_files,1);
lat_f = nan(n_prop_Output_files,1);
lon_f = nan(n_prop_Output_files,1);
pop_path{npop_files,1} = [];


[datenum_sorted,I] = sort([prop_Output_files.datenum],'ascend');
for i = 1:n_prop_Output_files
    
    prop_Output_files(i).datenum_sorted = datenum_sorted(i);
    prop_Output_files(i).name_sorted = prop_Output_files(I(i)).name;
    
    %end
    
    %for i = 1:n_prop_Output_files
    
    depvar_Output_files(i).datenum_sorted = datenum_sorted(i);
    depvar_Output_files(i).name_sorted = depvar_Output_files(I(i)).name;
    
    interp_Ascent_Output_files(i).datenum_sorted = datenum_sorted(i);
    interp_Ascent_Output_files(i).name_sorted = interp_Ascent_Output_files(I(i)).name;
    
    interp_Descent_Output_files(i).datenum_sorted = datenum_sorted(i);
    interp_Descent_Output_files(i).name_sorted = interp_Descent_Output_files(I(i)).name;
    
    DV_mapped_Ascent_Output_files(i).datenum_sorted = datenum_sorted(i);
    DV_mapped_Ascent_Output_files(i).name_sorted = DV_mapped_Ascent_Output_files(I(i)).name;
    
    DV_mapped_Descent_Output_files(i).datenum_sorted = datenum_sorted(i);
    DV_mapped_Descent_Output_files(i).name_sorted = DV_mapped_Descent_Output_files(I(i)).name;
    
    headingErrorDeadbandBounds_Output_files(i).datenum_sorted = datenum_sorted(i);
    headingErrorDeadbandBounds_Output_files(i).name_sorted = headingErrorDeadbandBounds_Output_files.name;
    
    alphaMachBounds_Output_files(i).datenum_sorted = datenum_sorted(i);
    alphaMachBounds_Output_files(i).name_sorted = alphaMachBounds_Output_files.name;
    
    
end
for i = 2:n_prop_Output_files
    
    
    headingErrorDeadbandBounds_Output_files(i).name = headingErrorDeadbandBounds_Output_files(i-1).name;
    headingErrorDeadbandBounds_Output_files(i).name_sorted = headingErrorDeadbandBounds_Output_files(i-1).name;
    headingErrorDeadbandBounds_Output_files(i).folder = headingErrorDeadbandBounds_Output_files(i-1).folder;
    headingErrorDeadbandBounds_Output_files(i).date = headingErrorDeadbandBounds_Output_files(i-1).date;
    headingErrorDeadbandBounds_Output_files(i).bytes = headingErrorDeadbandBounds_Output_files(i-1).bytes;
    headingErrorDeadbandBounds_Output_files(i).isdir = headingErrorDeadbandBounds_Output_files(i-1).isdir;
    headingErrorDeadbandBounds_Output_files(i).datenum = headingErrorDeadbandBounds_Output_files(i-1).datenum;
    
    alphaMachBounds_Output_files(i).name = alphaMachBounds_Output_files(i-1).name;
    alphaMachBounds_Output_files(i).name_sorted = alphaMachBounds_Output_files(i-1).name;
    alphaMachBounds_Output_files(i).folder = alphaMachBounds_Output_files(i-1).folder;
    alphaMachBounds_Output_files(i).date = alphaMachBounds_Output_files(i-1).date;
    alphaMachBounds_Output_files(i).bytes = alphaMachBounds_Output_files(i-1).bytes;
    alphaMachBounds_Output_files(i).isdir = alphaMachBounds_Output_files(i-1).isdir;
    alphaMachBounds_Output_files(i).datenum = alphaMachBounds_Output_files(i-1).datenum;
    
    
end


% Loop used to create each file's path and extract the initial conditions
% from the file path/name.
for i = 1:n_prop_Output_files
    
    % Create file path string from known data
    prop_path(i,:) = {strcat(prop_Output_files(i).folder,'/',prop_Output_files(i).name_sorted)};
    depvar_path(i,:) = {strcat(depvar_Output_files(i).folder,'/',depvar_Output_files(i).name_sorted)};
    interp_Ascent_path(i,:) = {strcat(interp_Ascent_Output_files(i).folder,'/',interp_Ascent_Output_files(i).name_sorted)};
    interp_Descent_path(i,:) = {strcat(interp_Descent_Output_files(i).folder,'/',interp_Descent_Output_files(i).name_sorted)};
    DV_mapped_Ascent_path(i,:) = {strcat(DV_mapped_Ascent_Output_files(i).folder,'/',DV_mapped_Ascent_Output_files(i).name_sorted)};
    DV_mapped_Descent_path(i,:) = {strcat(DV_mapped_Descent_Output_files(i).folder,'/',DV_mapped_Descent_Output_files(i).name_sorted)};
    headingErrorDeadbandBounds_path(i,:) = {strcat(headingErrorDeadbandBounds_Output_files(i).folder,'/',headingErrorDeadbandBounds_Output_files(i).name_sorted)};
    alphaMachBounds_path(i,:) = {strcat(alphaMachBounds_Output_files(i).folder,'/',alphaMachBounds_Output_files(i).name_sorted)};
    
    % check this on out to not have to do it with a loop
    %     info{1} = '/data/input/1001_1094.png';
    %     info{2} = '/data/input/1001_1094.png';
    %     info{3} = '/data/input/1209_7856.png';
    %    [~,coord,~] = cellfun(@fileparts,info,'un',0) ;
    %    C = cellfun(@(x) strsplit(x,'_'),coord,'Un',0) ;
    %    C = vertcat(C{:}) ;
    %   coorinates = cellfun(@str2num,C)
    
    % Extract strings of initial conditions from file path/name
    conditions = extractBetween(prop_Output_files(i).name_sorted,'_','.dat');
    conditions_split = strsplit(conditions{:},'_');
    %v_i_string = conditions_split{1};
    %gamma_i_string = conditions_split{2};
    %chi_i_string = conditions_split{3};
    tof_string = conditions_split{1};
    lat_f_string = conditions_split{2};
    lon_f_string = conditions_split{3};
    
    % Create set of Initial velocity and Initial flight-path angle
    % v_i(i) = str2double(v_i_string);
    % gamma_i(i) = str2double(gamma_i_string);
    % chi_i(i) = str2double(chi_i_string);
    tof(i) = str2double(tof_string);
    lat_f(i) = str2double(lat_f_string);
    lon_f(i) = str2double(lon_f_string);
end


% Loop used to create each file's path and extract the evolution from the
% file path/name.
if contains(pop_files(1).name, 'monteCarlo') == true
    % Create file path string from known data
    pop_path(1,:) = {strcat(pop_files(1).folder,'/',pop_files(1).name)};
    
    % Extract string of evolution file path/name
    pop_string = extractBetween(pop_files(1).name,22,'_');
    
    % Convert string to double
    pop_i(1) = str2double(cellstr(pop_string));
    
else
    for i = 1:npop_files
        
        % Create file path string from known data
        pop_path(i,:) = {strcat(pop_files(i).folder,'/',pop_files(i).name)};
        
        % Extract string of evolution file path/name
        pop_string = extractBetween(pop_files(i).name,12,'_');
        
        % Convert string to double
        pop_i(i) = str2double(cellstr(pop_string));
    end
end



% Loop used to create each file's path and extract the evolution from the
% file path/name.

if contains(fit_files(1).name, 'monteCarlo') == true
    
    
    % Create file path string from known data
    fit_path(1,:) = {strcat(fit_files(1).folder,'/',fit_files(1).name)};
    
    % Extract string of evolution file path/name
    fit_string = extractBetween(fit_files(1).name,19,'_');
    
    % Convert string to double
    fit_i(1) = str2double(cellstr(fit_string));
    
    
else
    
    for i = 1:nfit_files
        
        % Create file path string from known data
        fit_path(i,:) = {strcat(fit_files(i).folder,'/',fit_files(i).name)};
        
        % Extract string of evolution file path/name
        fit_string = extractBetween(fit_files(i).name,9,'_');
        
        % Convert string to double
        fit_i(i) = str2double(cellstr(fit_string));
    end
    
    
    
end





end

