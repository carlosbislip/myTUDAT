%% Trim Analysis of HORUS Aerodynamic Database


clearvars
close all
clc
dbstop if error


% cleanconfiguration = '/Users/bislip/Cloud Storage/OneDrive/School/TUDelft/Space Flight/Lit Study/Brought to you by Erwin Mooij/Archive/haerobody.mat';
% controlsurfaces = '/Users/bislip/Cloud Storage/OneDrive/School/TUDelft/Space Flight/Lit Study/Brought to you by Erwin Mooij/Archive/haerocs.mat';
% 
% importfile(cleanconfiguration);
% importfile(controlsurfaces);
% 
% 
% %[Xb,Yb] = meshgrid(Mach_cs,alpha_cs,db_cs);
% %%[Xb_grid,Yb_grid,Zb_grid] = ndgrid(Mach_cs,alpha_cs',db_cs);
% %XI = [Xb_grid(:) Yb_grid(:) Zb_grid(:)];
% 
% % V = dCmb+Cm0;
% % V = dCmb;
% % ii = 1;
% % jj = 1;
% % kk = 1;
% %close all
% alpha    = [0 5 10 15 20 25 30 35 40 45];
% mach     = [1.2 1.5 2 3 5 10 20];
% bodyflap = [-20 -10 0 10 20 30];
% elevon   = [-40 -30 -20 -10 0 10 20 30 40];
% Cm_int_0 = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
% Cm_int_bodyflap = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
% Cm_int_elevon = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
% Cm_int_matrix = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
% 
% CL_int_0 = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
% CL_int_bodyflap = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
% CL_int_elevon = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
% CL_int_matrix = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
% 
% CD_int_0 = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
% CD_int_bodyflap = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
% CD_int_elevon = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
% CD_int_matrix = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
% 
% for i = 1:numel(alpha)
%     
%     xq1 = deg2rad(alpha(i));
%     
%     for j = 1:numel(mach)
%         
%         xq2 = mach(j);
%         Cm0_int = interpn(alpha_cs',Mach_cs,Cm0,xq1,xq2,'spline');
%         CL0_int = interpn(alpha_cs',Mach_cs,CL0,xq1,xq2,'spline');
%         CD0_int = interpn(alpha_cs',Mach_cs,CD0,xq1,xq2,'spline');
%         
%         for k = 1:numel(bodyflap)
%             
%             xq31 = deg2rad(bodyflap(k));
%             dCmb_int = interpn(alpha_cs',Mach_cs,db_cs,dCmb,xq1,xq2,xq31,'spline');
%             dCLb_int = interpn(alpha_cs',Mach_cs,db_cs,dCLb,xq1,xq2,xq31,'spline');
%             dCDb_int = interpn(alpha_cs',Mach_cs,db_cs,dCDb,xq1,xq2,xq31,'spline');
%             
%             %Cm_int_0(i,j,k)        = Cm0_int;
%             %Cm_int_bodyflap(i,j,k) = dCmb_int;
%             %Cm_int_matrix(i,j,k)   = Cm0_int + dCmb_int;
%             
%             for l = 1:numel(elevon)
%                 
%                 xq32 = deg2rad(elevon(l));
%                 dCme_int = 2*interpn(alpha_cs',Mach_cs,dw_cs,dCmw,xq1,xq2,xq32,'spline');
%                 dCLe_int = 2*interpn(alpha_cs',Mach_cs,dw_cs,dCLw,xq1,xq2,xq32,'spline');
%                 dCDe_int = 2*interpn(alpha_cs',Mach_cs,dw_cs,dCDw,xq1,xq2,xq32,'spline');
%                 
%                 Cm_int_0(i,j,k,l) = Cm0_int;
%                 Cm_int_bodyflap(i,j,k,l) = dCmb_int;
%                 Cm_int_elevon(i,j,k,l) = dCme_int;
%                 Cm_int_matrix(i,j,k,l) = Cm0_int + dCmb_int + dCme_int;
%                 
%                 CL_int_0(i,j,k,l) = CL0_int;
%                 CL_int_bodyflap(i,j,k,l) = dCLb_int;
%                 CL_int_elevon(i,j,k,l) = dCLe_int;
%                 CL_int_matrix(i,j,k,l) = CL0_int + dCLb_int + dCLe_int;
%                 
%                 CD_int_0(i,j,k,l) = CD0_int;
%                 CD_int_bodyflap(i,j,k,l) = dCDb_int;
%                 CD_int_elevon(i,j,k,l) = dCDe_int;
%                 CD_int_matrix(i,j,k,l) = CD0_int + dCDb_int + dCDe_int;
%                 
%             end
%         end
%     end
% end

[ Cm_int_matrix, CL_int_matrix, CD_int_matrix ] = createAeroCoeffInterpolator(  );

%machrange = mach;
machrange = linspace(1.2,20,(20-1.2)/1);
%alpharange = deg2rad(alpha);
alpharange = deg2rad(linspace(0,45,(45-0)/1));

%bodyflaprange = deg2rad(bodyflap);
bodyflaprange = deg2rad(linspace(-20,30,(30--20)/1));
%elevonrange = deg2rad(elevon);
elevonrange = deg2rad(linspace(-40,40,(40--40)/1));

%database = nan(numel(machrange)*numel(alpharange),3);
database = [];
result = nan(numel(machrange)*numel(alpharange)*numel(elevonrange),5);
pair1 = [ ];
pair2 = [ ];
pair3 = [ ];
p = 1;


for i = 1:numel(machrange)
    
    mach = machrange(i);

    for j = 1:numel(alpharange)

        alpha = alpharange(j);
        result(p,1) = mach;
        result(p,2) = rad2deg(alpha);
        Cm_int_bounds = [ interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,bodyflaprange(1),0.0,'spline') ;...
            interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,bodyflaprange(end),0.0,'spline')];
        
        if Cm_int_bounds(1)*Cm_int_bounds(2) < 0
            
            a = bodyflaprange(1);
            b = bodyflaprange(end);
            
            
            root_bf = (a + b)/2;
            err = abs(interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,root_bf,0.0,'spline'));
            while err > 1e-7
                f_a       = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,a,0.0,'spline');
                f_root_bf = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,root_bf,0.0,'spline');
                
                if f_a*f_root_bf<0
                    b = root_bf;
                else
                    a = root_bf;
                end
                root_bf = (a + b)/2;
                err = abs(interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,root_bf,0.0,'spline'));
            end
            result(p,3) = rad2deg(root_bf);
            result(p,4) = rad2deg(0);
            result(p,5) = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,root_bf,0.0,'spline');
            
        else
            
            for k = 1:numel(bodyflaprange)
                rad2deg(bodyflaprange(k))
                Cm_int_bounds = [ interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,bodyflaprange(k),elevonrange(1),'spline') ;...
                    interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,bodyflaprange(k),elevonrange(end),'spline')];
                
                if Cm_int_bounds(1)*Cm_int_bounds(2) < 0
                    
                    a = elevonrange(1);
                    b = elevonrange(end);
                    
                    root_el = (a + b)/2;
                    err = abs(interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,bodyflaprange(k),root_el,'spline'));
                    while err > 1e-7
                        f_a       = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,bodyflaprange(k),a,'spline');
                        f_root_el = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,bodyflaprange(k),root_el,'spline');
                        
                        if f_a*f_root_el<0
                            b = root_el;
                        else
                            a = root_el;
                        end
                        root_el = (a + b)/2;
                        err = abs(interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,bodyflaprange(k),root_el,'spline'));
                    end
                    result(p,1) = mach;
                    result(p,2) = rad2deg(alpha);
                    result(p,3) = rad2deg(bodyflaprange(k));
                    result(p,4) = rad2deg(root_el);
                    result(p,5) = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,bodyflaprange(k),root_el,'spline');
                    
                    p = p + 1;
                end
            end
            % p = p -1;
            %{
            Cm_int_bounds = [ interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,bodyflaprange(1),elevonrange(1),'spline') ...
                interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,bodyflaprange(1),elevonrange(end),'spline') ;...
                interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,bodyflaprange(end),elevonrange(1),'spline') ...
                interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,bodyflaprange(end),elevonrange(end),'spline')]
            
            delta = 0.001;
            
            
            % Find local minima generated by body flap deflection angle
            
            % x(n+1) = x(n) - f_prime/f_2prime
            
            x = linspace(bodyflaprange(1),bodyflaprange(end),(bodyflaprange(end) - bodyflaprange(1) / delta));
            y = linspace(elevonrange(1),elevonrange(end),(elevonrange(end) - elevonrange(1) / delta));
            f_eval = (interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x,y,'spline'));
            [X,Y] = meshgrid(rad2deg(x'),rad2deg(y'));
            squeezed = squeeze(f_eval);
            
            y1 = zeros(1,numel(y));
            f_eval1 = (interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x,y1,'spline'));
            [X,Y1] = meshgrid(rad2deg(x'),rad2deg(y1'));
            squeezed1 = squeeze(f_eval1);
            
            
            x1 = zeros(1,numel(x));
            f_eval2 = (interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x1,y,'spline'));
            [X1,Y] = meshgrid(rad2deg(x1'),rad2deg(y'));
            squeezed2 = squeeze(f_eval2);
            
            
            fig_num = 536;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            % title(strcat('Mach through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
            xlim([bodyflap(1) bodyflap(end)])
            ylim([elevon(1) elevon(end)])
            zlim(.1*[-1 1])
            xlabel('bodyflap (deg)') % x-axis label
            ylabel('elevon (deg)') % y-axis label
            zlabel('Cm (-)') % y-axis label
            set(gca,'XTick', bodyflap(1):(bodyflap(end) - bodyflap(1))/10:bodyflap(end));
            set(gca,'YTick', elevon(1):(elevon(end) - elevon(1))/10:elevon(end));
            set(gca,'ZTick', -.1:.01:.1);
            hold on
            mesh(X,Y,squeezed')
            mesh(X,Y,zeros(size(X,1),size(X,2)))
            hold off
            
            
            fig_num = 537;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            % title(strcat('Mach through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
            xlim([bodyflap(1) bodyflap(end)])
            ylim([elevon(1) elevon(end)])
            zlim(.1*[-1 1])
            xlabel('bodyflap (deg)') % x-axis label
            ylabel('elevon (deg)') % y-axis label
            zlabel('Cm (-)') % y-axis label
            set(gca,'XTick', bodyflap(1):(bodyflap(end) - bodyflap(1))/10:bodyflap(end));
            set(gca,'YTick', elevon(1):(elevon(end) - elevon(1))/10:elevon(end));
            set(gca,'ZTick', -.1:.01:.1);
            hold on
            mesh(X,Y1,squeezed1')
            mesh(X,Y,zeros(size(X,1),size(X,2)))
            hold off
            
            
            
            
            fig_num = 538;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            % title(strcat('Mach through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
            xlim([bodyflap(1) bodyflap(end)])
            ylim([elevon(1) elevon(end)])
            zlim(.1*[-1 1])
            xlabel('bodyflap (deg)') % x-axis label
            ylabel('elevon (deg)') % y-axis label
            zlabel('Cm (-)') % y-axis label
            set(gca,'XTick', bodyflap(1):(bodyflap(end) - bodyflap(1))/10:bodyflap(end));
            set(gca,'YTick', elevon(1):(elevon(end) - elevon(1))/10:elevon(end));
            set(gca,'ZTick', -.1:.01:.1);
            hold on
            mesh(X1,Y,squeezed2')
            mesh(X,Y,zeros(size(X,1),size(X,2)))
            
            f_eval = squeeze(interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x,y,'spline'));
            
            plot3(rad2deg(x),rad2deg(y))
            
            plot(rad2deg(x),f_eval)
            [ min, idx ] = min( abs(f_eval) );
            
            
            x = x(idx);
            
            err = abs(interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x,0.0,'spline'));
            while err > 1e-7
                
                f = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x,0,'spline');
                f_p = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x + delta,0.0,'spline');
                f_n = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x - delta,0.0,'spline');
                
                f_prime = (f_p - f_n ) / (2 * delta);
                f_2prime = (f_p - 2 * f + f_n ) / delta^2;
                
                
                x = x - f_prime/f_2prime;
                err = abs(interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x,0.0,'spline'));
                
            end
            
            
            
            
            
            
            
            x = [ 0 ; 0 ];
            x_hist = x';
            err = abs(interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x(1),x(2),'spline'));
            while err > 1e-7
                
                f = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x(1),x(2),'spline');
                
                f_xp = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x(1) + delta,x(2),'spline');
                f_xn = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x(1) - delta,x(2),'spline');
                
                f_x = ( f_xp - f_xn ) / ( 2 * delta);
                
                f_yp = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x(1),x(2) + delta,'spline');
                f_yn = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x(1),x(2) - delta,'spline');
                
                f_pp = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x(1) + delta,x(2) + delta,'spline');
                f_nn = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x(1) - delta,x(2) - delta,'spline');
                
                f_y = ( f_yp - f_yn ) / ( 2 * delta);
                
                J = [ f_x ;f_y ];
                
                f_xx = ( f_xp - 2 * f + f_xn ) / delta^2;
                
                f_yy = ( f_yp - 2 * f + f_yn ) / delta^2;
                
                f_xy = ( f_pp - f_xp - f_yp + 2 * f - f_xn - f_yn + f_nn ) / ( 2 * delta^2 );
                
                H = [ f_xx f_xy ; f_xy f_yy ];
                x_dif = -inv( H ) * J;
                
                x = x + x_dif
                x_hist = [ x_hist ; x'];
                err = abs(interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,x(1),x(2),'spline'));
            end
            
            result(p,3) = rad2deg(x(1));
            result(p,4) = rad2deg(x(2));
            result(p,5) = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,alpha,mach,root_bf,root_el,'spline');
        end
        
        
        
        
        
        
            %}
        end
        
        p = p + 1;
        if j == numel(alpharange)
            'ssss'
        end
    end
end

scatter3(result(:,1),result(:,2),result(:,4))


result(any(isnan(result), 2), :) = [];



x = result(:,1);
y = result(:,2);
z = result(:,3);


tri = delaunay(x,y);
[r,c] = size(tri);
disp(r)

h = trisurf(tri, x, y, z);
axis vis3d

axis off
l = light('Position',[-50 -15 29])
set(gca,'CameraPosition',[208 -50 7687])
lighting phong
shading interp
colorbar EastOutside
