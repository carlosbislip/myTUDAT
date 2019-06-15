function [ evolutions ] = ...
    ...
    ...
    Get_Trajectories(...
    evolutions,...
    Folder_Path_List,...
    pop_Location,...
    pop_path,...
    prop_path,...
    depvar_path,...
    interp_Ascent_path,...
    interp_Descent_path,...
    DV_mapped_Ascent_path,...
    DV_mapped_Descent_path,...
    headingErrorDeadbandBounds_path,...
    alphaMachBounds_path,...
    v_i,gamma_i,pop_i,lon_i_rad,lat_f_deg,lon_f_deg,startEpoch)
%GET_TRAJECTORIES Summary of this function goes here
%   Detailed explanation goes here
disp('GET_TRAJECTORIES')


if pop_i > 0
    populationSize = numel(prop_path)/(numel(pop_i) + 1);
    a = linspace(populationSize,numel(prop_path),numel(pop_i) + 1);
    start = [1;a(1:end-1)'+1];
    
    p = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(p)
        parpool(4);
    end
    
    strWindowTitle = convertCharsToStrings('Getting Trajectories');
    nEvolutions = (numel(pop_i) + 1);
    ppm = ParforProgMon(strWindowTitle, nEvolutions);
else
    %    population = numel(prop_path)/(numel(pop_i) + 1);
    
    %      a = linspace(population,numel(prop_path),numel(pop_i) + 1);
    populationSize =  numel(prop_path);
    
    a = populationSize ;
    start = [1;a(end)];
    p = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(p)
        parpool(4);
    end
    nEvolutions = 1;
    strWindowTitle = convertCharsToStrings('Getting Trajectories');
    ppm = ParforProgMon(strWindowTitle, nEvolutions);
    
end

%C = [b a'];
%start = b;
finish = a';
%evolutions(numel(pop_i)).evolution(population,1).individual = [];


p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool(4);
end
Folder_prefix = 'OUTPUT*/';

prop_File_prefix = 'propagationHistory/propagationHistory_';
depvar_File_prefix = 'dependentVariables/dependentVariables_';
interp_Ascent_File_prefix = 'evaluatedInterpolatorsAscent/evaluatedInterpolatorsAscent_';
interp_Ascent_LB_File_prefix = 'evaluatedInterpolatorsAscent_LB/evaluatedInterpolatorsAscent_';
interp_Ascent_UB_File_prefix = 'evaluatedInterpolatorsAscent_UB/evaluatedInterpolatorsAscent_';
interp_Descent_File_prefix = 'evaluatedInterpolatorsDescent/evaluatedInterpolatorsDescent_';
interp_Descent_LB_File_prefix = 'evaluatedInterpolatorsDescent_LB/evaluatedInterpolatorsAscent_';
interp_Descent_UB_File_prefix = 'evaluatedInterpolatorsDescent_UB/evaluatedInterpolatorsAscent_';
%monteCarloPopulation_prefix = 'monteCarloPopulation*';

DV_mapped_Ascent_File_prefix = 'map_DV_mapped_Ascent/map_DV_mapped_Ascent_';
DV_mapped_Descent_File_prefix = 'map_DV_mapped_Descent/map_DV_mapped_Descent_';
headingErrorDeadbandBounds_prefix = 'headingErrorDeadBandBounds';
alphaMachBounds_prefix = 'alphaMachBounds';

headingErrorDeadbandBounds_source = { convertStringsToChars(strcat( Folder_Path_List,"/",headingErrorDeadbandBounds_prefix ) ) };
alphaMachBounds_source   = { convertStringsToChars(strcat( Folder_Path_List,"/",alphaMachBounds_prefix ) ) };

fid = fopen(headingErrorDeadbandBounds_source{:});
headingErrorDeadbandBounds = dlmread(headingErrorDeadbandBounds_source{:},',');
fclose(fid);

fid = fopen(alphaMachBounds_source{:});
alphaMachBounds = dlmread(alphaMachBounds_source{:},',');
fclose(fid);

for k = 1:nEvolutions
    
    evolutions(k).evolution = k;
    
    pp = 1;
    %     prop_source = prop_path(start(k):finish(k),:);
    pop_source               = pop_path(k,:);
    fit_source               = { strrep(pop_source{k,:},'population','fitness') };
    % depvar_source            = depvar_path(start(k):finish(k),:);
    % interp_Ascent_source     = interp_Ascent_path(start(k):finish(k),:);
    % interp_Descent_source    = interp_Descent_path(start(k):finish(k),:);
    % DV_mapped_Ascent_source  = DV_mapped_Ascent_path(start(k):finish(k),:);
    % DV_mapped_Descent_source = DV_mapped_Descent_path(start(k):finish(k),:);
   % headingErrorDeadbandBounds_source = headingErrorDeadbandBounds_path(start(k):finish(k),:);
   % alphaMachBounds_source   = alphaMachBounds_path(start(k):finish(k),:);
    
    
    %monteCarloPopulation_source = { convertStringsToChars(strcat( Folder_Path_List,"/",monteCarloPopulation_prefix,".dat" ) ) };
    
    %fid = fopen(monteCarloPopulation_source{:});
    % monteCarloPopulation = dlmread(monteCarloPopulation_source{:},'\t');
    % fclose(fid);
    % monteCarloPopulation = num2cell(monteCarloPopulation(:,2:end));
    
    % fid = fopen(pop_source{k,:});
    % individuals = dlmread(pop_source{k,:},'\t');
    % fclose(fid);
    
    entirePopulation     = loadPopulationFile( pop_source{k,:} );
    printedIndices       = find([entirePopulation{:,2}] > 0);
    nonPrintedIndices    = find([entirePopulation{:,2}] == 0);
    printedPopulation    = entirePopulation(printedIndices,:);
    nonPrintedPopulation = entirePopulation(nonPrintedIndices,:);
   
    evolutions(k).populationSize           = size(entirePopulation,1);
    evolutions(k).printedIndices           = printedIndices;
    evolutions(k).nonPrintedIndices        = nonPrintedIndices;
    evolutions(k).printedPopulationSize    = size(printedPopulation,1);
    evolutions(k).nonPrintedPopulationSize = size(nonPrintedPopulation,1);
    
    evolutions(k).entirePopulationDV     = loadDecisionVectorFromPopulationCellArray( entirePopulation );
    evolutions(k).printedPopulationDV    = loadDecisionVectorFromPopulationCellArray( printedPopulation );
    evolutions(k).nonPrintedPopulationDV = loadDecisionVectorFromPopulationCellArray( nonPrintedPopulation );
    
    entireFitness     = loadFitnessFile( fit_source{k,:} );
    printedFitness    = entireFitness(printedIndices,:);
    nonPrintedFitness = entireFitness(nonPrintedIndices,:);
    
    evolutions(k).entirePopulationFitness      = loadFitnessVectorFromFitnessCellArray( entireFitness );
    evolutions(k).printedPopulationFitness     = loadFitnessVectorFromFitnessCellArray( printedFitness );
    evolutions(k).nonPrintedPopulationFitnessn = loadFitnessVectorFromFitnessCellArray( nonPrintedFitness );
    
    
    for p = 1:size(printedPopulation,1)
        
        depvar_source            = { convertStringsToChars(strcat( Folder_Path_List,"/",depvar_File_prefix,printedPopulation{p,1},".dat" ) ) };
        interp_Ascent_source     = { convertStringsToChars(strcat( Folder_Path_List,"/",interp_Ascent_File_prefix,printedPopulation{p,1},".dat" ) ) };
        interp_Descent_source    = { convertStringsToChars(strcat( Folder_Path_List,"/",interp_Descent_File_prefix,printedPopulation{p,1},".dat" ) ) };
        DV_mapped_Ascent_source  = { convertStringsToChars(strcat( Folder_Path_List,"/",DV_mapped_Ascent_File_prefix,printedPopulation{p,1},".dat" ) ) };
        DV_mapped_Descent_source = { convertStringsToChars(strcat( Folder_Path_List,"/",DV_mapped_Descent_File_prefix,printedPopulation{p,1},".dat" ) ) };
        
        fid = fopen(depvar_source{:});
        depvar = dlmread(depvar_source{:},',');
        fclose(fid);
        
        fid = fopen(interp_Ascent_source{:});
        interp_Ascent = dlmread(interp_Ascent_source{:},',');
        fclose(fid);
        
        fid = fopen(interp_Descent_source{:});
        interp_Descent = dlmread(interp_Descent_source{:},',');
        fclose(fid);
        
        fid = fopen(DV_mapped_Ascent_source{:});
        DV_mapped_Ascent = dlmread(DV_mapped_Ascent_source{:},',');
        fclose(fid);
        
        fid = fopen(DV_mapped_Descent_source{:});
        DV_mapped_Descent = dlmread(DV_mapped_Descent_source{:},',');
        fclose(fid);
        

      
        evolutions(k).trajectories(p).individual = loadIndividualDataFromDependentVariableFile( depvar, interp_Ascent, interp_Descent, DV_mapped_Ascent, DV_mapped_Descent, headingErrorDeadbandBounds, alphaMachBounds );
        
       
        %evolutions(k).individuals.gamma_i(p)   = evolutions(k).trajectories(p).individual.flight_path_angle(1);
        %evolutions(k).individuals.chi_i(p)     = evolutions(k).trajectories(p).individual.heading_angle(1);
        %evolutions(k).individuals.lat_f_deg(p) = evolutions(k).trajectories(p).individual.latitude_angle(end);
        %evolutions(k).individuals.lon_f_deg(p) = evolutions(k).trajectories(p).individual.longitude_angle(end);
        evolutions(k).individuals.tof(p)       = evolutions(k).trajectories(p).individual.time_vector(end);
        %evolutions(k).individuals.interp_E_mapped_Ascent(p)  = evolutions(k).trajectories(p).individual.interp_E_mapped_Ascent(end);
        %evolutions(k).individuals.interp_E_mapped_Descent(p)  = evolutions(k).trajectories(p).individual.interp_E_mapped_Descent(end);
        
        %evolutions(k).fitness.dif_lat(p)   = evolutions(k).trajectories(p).individual.latitude_angle(end) - lat_f_deg;
        %evolutions(k).fitness.dif_lon(p)   = evolutions(k).trajectories(p).individual.longitude_angle(end) - lon_f_deg;
        %evolutions(k).fitness.dif_d_deg(p) = evolutions(k).trajectories(p).individual.distance_to_go(end) - 0.75;
        %evolutions(k).fitness.dif_h(p)     = evolutions(k).trajectories(p).individual.height(end) - 25000;
        evolutions(k).fitness.tof(p)       = evolutions(k).trajectories(p).individual.time_vector(end);
        
        
        pp = pp + 1;
        % end
    end
    
    %
    %     [ aaa , idx1 ] = min(evolutions(k).fitness.dif_d_deg);
    %     [ bbb , idx2 ] = min(evolutions(k).fitness.dif_h);
    %     [ ccc , idx3 ] = min(evolutions(k).fitness.tof);
    [ ddd , idx4 ] = max(evolutions(k).individuals.tof);
    %I = [ idx1 idx2 idx3 idx4 ];
    I = [ idx4 ];
    % criteria = [ {'d_deg'} {'h'} {'tof'} {'max_tof'} ];
    criteria = [ {'max_tof'} ];
    %  I_12 = intersect(idx1,idx2);
    %  I_123 = intersect(I_12,idx3);
    %  I_123 = idx1;
    evolutions(k).max_tof                     = max(evolutions(k).individuals.tof);
    %evolutions(k).max_interp_E_mapped_Ascent  = max(evolutions(k).individuals.interp_E_mapped_Ascent);
    %evolutions(k).max_interp_E_mapped_Descent = max(evolutions(k).individuals.interp_E_mapped_Descent);
    
    for i = 1:numel(I)
        evolutions(k).best(i).criteria  = criteria(i);
        evolutions(k).best(i).index     = I(i);
        % evolutions(k).best(i).v_i       = evolutions(k).individuals.v_i(I(i));
        %evolutions(k).best(i).gamma_i   = evolutions(k).individuals.gamma_i(I(i));
        %evolutions(k).best(i).chi_i     = evolutions(k).individuals.chi_i(I(i));
        %evolutions(k).best(i).dif_d_deg = evolutions(k).fitness.dif_d_deg(I(i));
        %evolutions(k).best(i).dif_h     = evolutions(k).fitness.dif_h(I(i));
        %evolutions(k).best(i).tof       = evolutions(k).fitness.tof(I(i));
        evolutions(k).best(i).max_tof   = evolutions(k).individuals.tof(I(i));
    end
    %}
    ppm.increment();
    
    %evolutions(k).filteredIndividuals = filteredIndividuals;
    
end



end










