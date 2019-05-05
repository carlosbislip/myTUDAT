

k = 1;
ii = 328;

trajectory = compilation.evolutions(k).trajectories(ii).individual;

omega_E = 7.2921150e-5;

mass = trajectory.mass;
V = trajectory.airspeed;
r = trajectory.altitude;
c_delta = cos(deg2rad(trajectory.latitude_angle));
s_delta = sin(deg2rad(trajectory.latitude_angle));
c_chi   = cos(deg2rad(trajectory.heading_angle));
s_chi   = sin(deg2rad(trajectory.heading_angle));
c_alpha = cos(deg2rad(trajectory.evaluated_angle_of_attack));
s_alpha = sin(deg2rad(trajectory.evaluated_angle_of_attack));

c_gamma = cos(deg2rad(trajectory.flight_path_angle));
s_gamma = sin(deg2rad(trajectory.flight_path_angle));
c_phi   = cos(deg2rad(trajectory.evaluated_thrust_azimuth_angle));
s_phi   = sin(deg2rad(trajectory.evaluated_thrust_azimuth_angle));
c_eps   = cos(deg2rad(trajectory.evaluated_thrust_elevation_angle));
s_eps   = sin(deg2rad(trajectory.evaluated_thrust_elevation_angle));

T = trajectory.thrustMagnitude;
L = trajectory.currentLiftForce;

a = L + T.*(s_alpha.*c_phi.*c_eps + c_alpha.*s_eps);
b = -T.*s_phi.*c_eps;
c = sqrt(a.^2 + b.^2);
phi = atan2(b,a);
g_d = trajectory.local_gravity_2;
g_n = trajectory.local_gravity_1;
argumentNUM = omega_E.*omega_E.*r.*c_delta.*(c_delta.*c_gamma + s_gamma.*s_delta.*c_chi) + V.*V.*c_gamma./r + 2.*omega_E.*V.*c_delta.*s_chi - g_d.*c_gamma + g_n.*s_gamma.*c_chi;
argument = -(mass./c).*argumentNUM;


sigma1 = acos(argument);

sigma1(find(abs(argument) > 1)) = nan;

sigma = rad2deg(sigma1 + phi);


B = [trajectory.evaluated_bank_angle abs(trajectory.commanded_bank_angle) abs(trajectory.skip_suppression_limit) argument phi sigma V T L];

size(sigma)