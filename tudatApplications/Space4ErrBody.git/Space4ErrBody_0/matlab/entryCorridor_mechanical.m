%% Entry Corridor: Mechanical Load




cleanconfiguration = '/Users/bislip/Cloud Storage/OneDrive/School/TUDelft/Space Flight/Lit Study/Brought to you by Erwin Mooij/Archive/haerobody.mat';
controlsurfaces = '/Users/bislip/Cloud Storage/OneDrive/School/TUDelft/Space Flight/Lit Study/Brought to you by Erwin Mooij/Archive/haerocs.mat';

importfile(cleanconfiguration);
importfile(controlsurfaces);

[ Cm_int_matrix, CL_int_matrix, CD_int_matrix ] = createAeroCoeffInterpolator(alpha_cs, Mach_cs, db_cs, dw_cs, Cm0, dCmb, dCmw, CL0, dCLb, dCLw, CD0, dCDb, dCDw );




%Area = 110;









alpha = deg2rad(40*[1 1 1]);
Ma =[ 20 20 27 ];
bodyflap = deg2rad([ 0 10 20 ]);
elevon = deg2rad([0 0 0]);

%CL = [ interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha(1),alpha(1),0.0,0.0,'linear') ;...

Cm = 10;


%Cm =  interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,40,27,root,0,'linear')
a = -deg2rad(20);
b = deg2rad(30);
alpha = 40;
mach = 27;
p = 0;
root = (a + b)/2;
err = abs(interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,deg2rad(alpha),mach,root,0,'linear'));
while err > 1e-7
    if interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,deg2rad(alpha),mach,a,0,'linear')*interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,deg2rad(alpha),mach,root,0,'linear')<0
        b = root;
    else
        a = root;
    end
    root = (a + b)/2
    err = abs(interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,deg2rad(alpha),mach,root,0,'linear'));
    p = p + 1
end
rad2deg(root)



% 
% 
% 
% %V_q_dyn = nan(numel(compilation(1).evolutions(1).trajectories(1).individual.height),1);
% for i = 1:numel(compilation(1).evolutions(1).trajectories(1).individual.localDensity)
% 
% 
% a = 500;
% b = 8000000;
% f_a =1;
% 
% compilation(1).evolutions(1).trajectories(1).individual.localDensity(i)/
% 
% compilation(1).evolutions(1).trajectories(1).individual.localDensity(i)*(1/2)*a^2 - q_dyn_max;
% f_b = compilation(1).evolutions(1).trajectories(1).individual.localDensity(i)*(1/2)*b^2 - q_dyn_max;
% 
% 
% 
% 
% 
% 
% 
% root = (a + b)/2;
% err = abs(compilation(1).evolutions(1).trajectories(1).individual.localDensity(i)*(1/2)*root^2 - q_dyn_max);
% while err > 1e-7
%     f_a       = compilation(1).evolutions(1).trajectories(1).individual.localDensity(i)*(1/2)*a^2 - q_dyn_max;
%     f_root = compilation(1).evolutions(1).trajectories(1).individual.localDensity(i)*(1/2)*root^2 - q_dyn_max;
%     
%     if f_a*f_root<0
%         b = root;
%     else
%         a = root;
%     end
%     root = (a + b)/2;
%     err = abs(compilation(1).evolutions(1).trajectories(1).individual.localDensity(i)*(1/2)*root^2 - q_dyn_max);
% end
% V_q_dyn(i) = root;
% end

