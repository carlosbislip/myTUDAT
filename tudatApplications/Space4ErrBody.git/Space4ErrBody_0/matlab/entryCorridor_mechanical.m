%% Entry Corridor: Mechanical Load

[ Cm_int_matrix, CL_int_matrix, CD_int_matrix ] = createAeroCoeffInterpolator(  );




Area = 110;









alpha = deg2rad([40 10]);
Ma = a/
        CL = [ interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha(1),alpha(1),0.0,0.0,'spline') ;...





V_q_dyn = nan(numel(compilation(1).evolutions(1).trajectories(1).individual.height),1);
for i = 1:numel(compilation(1).evolutions(1).trajectories(1).individual.localDensity)
    
    
      a = 500;
    b = 8000000;
    f_a = 
    
    compilation(1).evolutions(1).trajectories(1).individual.localDensity(i)/
    
    compilation(1).evolutions(1).trajectories(1).individual.localDensity(i)*(1/2)*a^2 - q_dyn_max;
    f_b = compilation(1).evolutions(1).trajectories(1).individual.localDensity(i)*(1/2)*b^2 - q_dyn_max;

    
    
    

    
    
    root = (a + b)/2;
    err = abs(compilation(1).evolutions(1).trajectories(1).individual.localDensity(i)*(1/2)*root^2 - q_dyn_max);
    while err > 1e-7
        f_a       = compilation(1).evolutions(1).trajectories(1).individual.localDensity(i)*(1/2)*a^2 - q_dyn_max;
        f_root = compilation(1).evolutions(1).trajectories(1).individual.localDensity(i)*(1/2)*root^2 - q_dyn_max;
        
        if f_a*f_root<0
            b = root;
        else
            a = root;
        end
        root = (a + b)/2;
        err = abs(compilation(1).evolutions(1).trajectories(1).individual.localDensity(i)*(1/2)*root^2 - q_dyn_max);
    end
   V_q_dyn(i) = root;
end

