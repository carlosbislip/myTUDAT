function [ evolutions ] = Analyze_Evolution(evolutions,v_i,gamma_i,chi_i,lat_f,lon_f,prop_path,pop_path,pop_i,fit_path,fit_i,lat_f_deg,lon_f_deg)
%UNTITLED Summary of this function goes here
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

for k = 1:(numel(pop_i))
    %      pp = 1;
    %data_source = data_path(start(k):finish(k),:);
 %   evolutions(k).evolution = k;
    
    %if k == 1
   % population =  numel(prop_path)/(numel(pop_i) + 1);
    % population =  1000;
    %
    %     evolutions(k).evolution = [zeros(population,1) v_i(start(k):finish(k)) ...
    %         gamma_i(start(k):finish(k)) chi_i(start(k):finish(k))];
%     fid = fopen(pop_path{k,:});
%     individuals = dlmread(pop_path{k,:});
%     fclose(fid);
%     
%     
%     
%     evolutions(k).individuals.v_i     = individuals(:,1);
%     evolutions(k).individuals.gamma_i = individuals(:,2);
%     evolutions(k).individuals.chi_i   = individuals(:,3);
% 
%     
%     fid = fopen(fit_path{k,:});
%     fitness = dlmread(fit_path{k,:});
%     fclose(fid);
%     
%     evolutions(k).fitness.dif_norm  = fitness(:,1);
%     evolutions(k).fitness.dif_lat   = fitness(:,2);
%     evolutions(k).fitness.dif_lon   = fitness(:,3);
%     evolutions(k).fitness.dif_d_deg = fitness(:,4);
%     evolutions(k).fitness.dif_h     = fitness(:,5);
%     evolutions(k).fitness.tof       = fitness(:,6);
%     
%     %  else
%     
%     %         fid = fopen(convertStringsToChars(pop_path{k-1,:}));
%     %         evolutions(k).evolution = dlmread(convertStringsToChars(pop_path{k-1,:}),'\t');
    %         fclose(fid);
    %         fid = fopen(convertStringsToChars(fit_path{k-1,:}));
    %         evolutions(k).fitness = dlmread(convertStringsToChars(fit_path{k-1,:}),'\t');
    %         fclose(fid);
    %
    %     end
    %     evolutions(k).v_i              = v_i(start(k):finish(k));
    %     evolutions(k).gamma_i          = gamma_i(start(k):finish(k));
    %     evolutions(k).chi_i            = chi_i(start(k):finish(k));
    % %    evolutions(k).aoa_i            = aoa_i(start(k):finish(k));
    %     evolutions(k).lat              = lat_f(start(k):finish(k));
    %     evolutions(k).lon              = lon_f(start(k):finish(k));
    %     evolutions(k).tof              = tof(start(k):finish(k));
    
%     [ aaa , idx1 ] = min(evolutions(k).fitness.dif_d_deg);
%     [ bbb , idx2 ] = min(evolutions(k).fitness.dif_h);
%     [ ccc , idx3 ] = min(evolutions(k).fitness.tof);
%     I = [ idx1 idx2 idx3 ];
%     criteria = [{'d_deg'} {'h'} {'tof'} ];
%   %  I_12 = intersect(idx1,idx2);
%   %  I_123 = intersect(I_12,idx3);
%   %  I_123 = idx1;
% 
%     for i = 1:3
%     evolutions(k).best(i).criteria  = criteria(i);
%     evolutions(k).best(i).index     = I(i);
%     evolutions(k).best(i).v_i       = evolutions(k).individuals.v_i(I(i));
%     evolutions(k).best(i).gamma_i   = evolutions(k).individuals.gamma_i(I(i));
%     evolutions(k).best(i).chi_i     = evolutions(k).individuals.chi_i(I(i));
%     evolutions(k).best(i).dif_d_deg = evolutions(k).fitness.dif_d_deg(I(i));
%     evolutions(k).best(i).dif_h     = evolutions(k).fitness.dif_h(I(i));
%     evolutions(k).best(i).tof       = evolutions(k).fitness.tof(I(i));
%     end
    % end
    %
    %
    % best_range = [min([ evolutions.best_v_i]) ...
    %     min([ evolutions.best_gamma_i]) ...
    %     min([ evolutions.best_chi_i])...
    %     max([ evolutions.best_v_i]) ...
    %     max([ evolutions.best_gamma_i]) ...
    %     max([ evolutions.best_chi_i]) ...
    
    %new_bounds = [0.95*best_range(1,:);1.05*best_range(2,:)];
    
    
    
    
    
    
    
end

