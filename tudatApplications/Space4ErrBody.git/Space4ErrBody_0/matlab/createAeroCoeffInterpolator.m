function [ Cm_int_matrix, CL_int_matrix, CD_int_matrix, alpha_cs, Mach_cs, dw_cs, dCmw ] = createAeroCoeffInterpolator(alpha_cs, Mach_cs, db_cs, dw_cs, Cm0, dCmb, dCmw, CL0, dCLb, dCLw, CD0, dCDb, dCDw )



%[Xb,Yb] = meshgrid(Mach_cs,alpha_cs,db_cs);
%%[Xb_grid,Yb_grid,Zb_grid] = ndgrid(Mach_cs,alpha_cs',db_cs);
%XI = [Xb_grid(:) Yb_grid(:) Zb_grid(:)];

% V = dCmb+Cm0;
% V = dCmb;
% ii = 1;
% jj = 1;
% kk = 1;
%close all
alpha    = [0 5 10 15 20 25 30 35 40 45];
mach     = [1.2 1.5 2 3 5 10 20];
bodyflap = [-20 -10 0 10 20 30];
elevon   = [-40 -30 -20 -10 0 10 20 30 40];
Cm_int_0 = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
Cm_int_bodyflap = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
Cm_int_elevon = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
Cm_int_matrix = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));

CL_int_0 = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
CL_int_bodyflap = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
CL_int_elevon = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
CL_int_matrix = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));

CD_int_0 = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
CD_int_bodyflap = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
CD_int_elevon = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
CD_int_matrix = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));

for i = 1:numel(alpha)
    
    xq1 = deg2rad(alpha(i));
    
    for j = 1:numel(mach)
        
        xq2 = mach(j);
        Cm0_int = interpn(alpha_cs',Mach_cs,Cm0,xq1,xq2,'linear');
        CL0_int = interpn(alpha_cs',Mach_cs,CL0,xq1,xq2,'linear');
        CD0_int = interpn(alpha_cs',Mach_cs,CD0,xq1,xq2,'linear');
        
        for k = 1:numel(bodyflap)
            
            xq31 = deg2rad(bodyflap(k));
            dCmb_int = interpn(alpha_cs',Mach_cs,db_cs,dCmb,xq1,xq2,xq31,'linear');
            dCLb_int = interpn(alpha_cs',Mach_cs,db_cs,dCLb,xq1,xq2,xq31,'linear');
            dCDb_int = interpn(alpha_cs',Mach_cs,db_cs,dCDb,xq1,xq2,xq31,'linear');
            
            %Cm_int_0(i,j,k)        = Cm0_int;
            %Cm_int_bodyflap(i,j,k) = dCmb_int;
            %Cm_int_matrix(i,j,k)   = Cm0_int + dCmb_int;
            
            for l = 1:numel(elevon)
                
                xq32 = deg2rad(elevon(l));
                dCme_int = 2*interpn(alpha_cs',Mach_cs,dw_cs,dCmw,xq1,xq2,xq32,'linear');
                dCLe_int = 2*interpn(alpha_cs',Mach_cs,dw_cs,dCLw,xq1,xq2,xq32,'linear');
                dCDe_int = 2*interpn(alpha_cs',Mach_cs,dw_cs,dCDw,xq1,xq2,xq32,'linear');
                
                Cm_int_0(i,j,k,l) = Cm0_int;
                Cm_int_bodyflap(i,j,k,l) = dCmb_int;
                Cm_int_elevon(i,j,k,l) = dCme_int;
                Cm_int_matrix(i,j,k,l) = Cm0_int + dCmb_int + dCme_int;
                
                CL_int_0(i,j,k,l) = CL0_int;
                CL_int_bodyflap(i,j,k,l) = dCLb_int;
                CL_int_elevon(i,j,k,l) = dCLe_int;
                CL_int_matrix(i,j,k,l) = CL0_int + dCLb_int + dCLe_int;
                
                CD_int_0(i,j,k,l) = CD0_int;
                CD_int_bodyflap(i,j,k,l) = dCDb_int;
                CD_int_elevon(i,j,k,l) = dCDe_int;
                CD_int_matrix(i,j,k,l) = CD0_int + dCDb_int + dCDe_int;
                
            end
        end
    end
end





end