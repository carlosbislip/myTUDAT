%V = linspace(500,8000,1/1e-3);
q_dyn_max = 95000;
rho = 2*q_dyn./V.^2;

compilation(1).evolutions(1).trajectories(1).individual.localDensity;
compilation(1).evolutions(1).trajectories(1).individual.height;

V_q_dyn = nan(numel(compilation(1).evolutions(1).trajectories(1).individual.height),1);
for i = 1:numel(compilation(1).evolutions(1).trajectories(1).individual.localDensity)
    
    
      a = 500;
    b = 8000000;
    f_a = compilation(1).evolutions(1).trajectories(1).individual.localDensity(i)*(1/2)*a^2 - q_dyn_max;
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







V_q_dyn = nan(numel(compilation(1).evolutions(1).trajectories(1).individual.height),1);
for i = 1:numel(compilation(1).evolutions(1).trajectories(1).individual.localDensity)
    
    
      a = 500;
    b = 8000000;
    f_a = compilation(1).evolutions(1).trajectories(1).individual.localDensity(i)*(1/2)*a^2 - q_dyn_max;
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

