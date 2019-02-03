function [prop_Output_files,n_prop_Output_files,depvar_Output_files,n_depvar_Output_files,interp_Ascent_Output_files,interp_Descent_Output_files,DV_mapped_Ascent_Output_files, DV_mapped_Descent_Output_files,pop_files,npop_files,fit_files,nfit_files] = Read_Folder_Content(p,prop_File_Path_List_prefix,depvar_File_Path_List_prefix,interp_Ascent_File_Path_List_prefix,interp_Decent_File_Path_List_prefix,DV_mapped_Ascent_File_Path_List_prefix,DV_mapped_Descent_File_Path_List_prefix,pop_file_path_prefix,fit_file_path_prefix)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


prop_Output_files = dir(prop_File_Path_List_prefix{p});
n_prop_Output_files = length(prop_Output_files);

depvar_Output_files = dir(depvar_File_Path_List_prefix{p});
n_depvar_Output_files = length(depvar_Output_files);

interp_Ascent_Output_files = dir(interp_Ascent_File_Path_List_prefix{p});
interp_Descent_Output_files = dir(interp_Decent_File_Path_List_prefix{p});

DV_mapped_Ascent_Output_files = dir(DV_mapped_Ascent_File_Path_List_prefix{p});
DV_mapped_Descent_Output_files = dir(DV_mapped_Descent_File_Path_List_prefix{p});


pop_files = dir(pop_file_path_prefix{p});
npop_files = length(pop_files);


fit_files = dir(fit_file_path_prefix{p});
nfit_files = length(fit_files);
end

