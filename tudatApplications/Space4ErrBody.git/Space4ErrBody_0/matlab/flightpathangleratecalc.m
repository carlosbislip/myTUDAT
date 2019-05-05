%% Flight-path Angle Rate

k = 1;
ii = 10;

trajectory = compilation.evolutions(k).trajectories(ii).individual;

omega_E = 7.2921150e-5;

mass = trajectory.mass;
V = trajectory.airspeed;
r = trajectory.altitude;
c_delta = cos(deg2rad(trajectory.latitude_angle));
s_delta = sin(deg2rad(trajectory.latitude_angle));
c_chi   = cos(deg2rad(trajectory.heading_angle));
s_chi   = sin(deg2rad(trajectory.heading_angle));
c_alpha = cos(deg2rad(trajectory.angle_of_attack));
s_alpha = sin(deg2rad(trajectory.angle_of_attack));
c_sigma = cos(deg2rad(trajectory.bank_angle));
s_sigma = sin(deg2rad(trajectory.bank_angle));

c_gamma = cos(deg2rad(trajectory.flight_path_angle));
s_gamma = sin(deg2rad(trajectory.flight_path_angle));
c_phi   = cos(deg2rad(trajectory.commanded_thrust_azimuth_angle));
s_phi   = sin(deg2rad(trajectory.commanded_thrust_azimuth_angle));
c_eps   = cos(deg2rad(trajectory.commanded_thrust_elevation_angle));
s_eps   = sin(deg2rad(trajectory.commanded_thrust_elevation_angle));

T = trajectory.thrustMagnitude;
L = trajectory.currentLiftForce;
g_d = trajectory.local_gravity_2;
g_n = trajectory.local_gravity_1;


F_gamma1 = L.*c_sigma./(mass.*V) ;
F_gamma2 = T.*(s_alpha.*c_phi.*c_eps + c_alpha.*s_eps).*c_sigma./(mass.*V) ;
F_gamma3 = -mass.*(g_d.*c_gamma)./(mass.*V);
F_gamma4 =  mass.*g_n.*s_gamma.*c_chi./(mass.*V);


F_gamma = F_gamma1 + F_gamma2 + F_gamma3 + F_gamma4;

term1 = 2.*omega_E.*c_delta.*s_chi;
term2 = V.*c_gamma./r;
term3 = (omega_E.*omega_E.*r.*c_delta.*(c_delta.*c_gamma + s_gamma.*s_delta.*c_chi))./V;

gamma_dot = rad2deg(F_gamma + term1 + term2 + term3);

B_names = {'bank_angle','thrust_azimuth' ,'thrust_elevation', 'Velocity' ,'Thrust' ,'Lift' ,'mass', 'F_gamma' ,'term1', 'term2', 'term3', 'gamma_dot'};

B = [trajectory.time_of_flight trajectory.evaluated_bank_angle ,trajectory.commanded_thrust_azimuth_angle, trajectory.commanded_thrust_elevation_angle ,V ,T, L, mass,  F_gamma, term1, term2 ,term3 ,gamma_dot];

F_gamma_table = [ F_gamma1  F_gamma2  F_gamma3  F_gamma4 F_gamma ];

%T = table(B)
