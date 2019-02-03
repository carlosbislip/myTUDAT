
clc
clear all
close all

      omega_E = 7.2921150e-5;

      delta = deg2rad(30);
      chi = deg2rad(280);
      V = 6000;
      r = 6400000;
      R_E = 6371000;
      J2 = 1082.626523e-6;
      mu = 3.986004418e14;
      gn = -3*J2*(mu/r^2)*((R_E/r)^2)*sin(delta)*cos(delta);
      gamma = linspace(-pi/2,pi/2,50000000);

%       a = (V^2/r)*tan(delta);
%       b = 2*omega_E*V*cos(delta)/tan(chi);
%       c = a + gn + omega_E*omega_E*r*cos(delta)*sin(delta) + 2*omega_E*V*sin(delta);
% 
% root1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
% root2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
% 
% A = -gn - omega_E*omega_E*r*cos(delta)*sin(delta) - 2*omega_E*V*sin(delta);
% 
% LHS = a - A;
%     
%     RHS = (a.*sin(gamma) + 2*omega_E*V*cos(delta)/tan(chi)).*sin(gamma);
% dif = LHS-RHS;
% %plot(dif)
% 
% derp = gn*sin(chi) + omega_E*omega_E*r*cos(delta)*sin(delta)*sin(chi) + (V^2/r)*cos(gamma).*cos(gamma)*tan(delta)*sin(chi) + 2*omega_E*V*(sin(delta)*sin(chi) - cos(delta)*sin(gamma)*cos(chi));
% derp1 = rad2deg(asin(derp));
% 
% 
% 
% plot(rad2deg(gamma),derp)

%LHS = -(V^2/r)*tan(delta)*tan(chi);
%RHS = (tan(chi)*(gn+(omega_E^2)*r*cos(delta)*sin(delta) + 2*omega_E*V*sin(delta)) - 2*omega_E*V*cos(delta)*sin(gamma))./((1-sin(gamma)).^2);
%derp =(sin(chi)./cos(gamma))*(gn+(omega_E^2)*r*cos(delta)*sin(delta)) + 2*omega_E*V*sin(delta) + (V^2/r)*tan(delta)*sin(chi)*cos(gamma) - 2*omega_E*V*cos(delta)*sin(gamma);
derp = gn*sin(chi) + (omega_E^2)*r*cos(delta)*sin(delta)*sin(chi) + (V^2/r)*(cos(gamma).^2)*tan(delta)*sin(chi) + 2*omega_E*V*sin(delta)*cos(gamma) - 2*omega_E*V*cos(delta)*sin(gamma)*cos(chi);
%derp = derp./(V*cos(gamma));

[derp0 idx] = min(abs(derp))

GAMMA1 = rad2deg(gamma(idx))


figure(1)
hold on
ylim(0.5*[-1 1])
plot(rad2deg(gamma),derp)
plot([-100 100],[0 0])
hold off
%plot(gamma,LHS*(ones(numel(gamma),1)))
% 
% 
%LHS = (V*tan(delta)*tan(chi)/(2*omega_E*cos(delta)) + tan(delta)/cos(chi) + (tan(chi)/(2*omega_E*V))*(gn/cos(delta) + (omega_E^2)*r*sin(delta)))/V;
%c = -LHS;
%b = 1/V;
%a = (V*tan(delta)*tan(chi)/(2*omega_E*cos(delta)) + tan(delta)/cos(chi))/V;

% gamma1 = rad2deg(asin((-b + sqrt(b*b - 4*a*c))/(2*a)));
% gamma2 = rad2deg(asin((-b - sqrt(b*b - 4*a*c))/(2*a)));

RHS = 1 + (gn*r)/(V^2*tan(delta)) + (omega_E*r*cos(delta)/V)^2;
phi = atan2(2*omega_E*V*sin(delta),-2*omega_E*V*cos(delta)*cos(chi));
derp333 = (sin(gamma).^2) - ((2*(omega_E*r*cos(delta)/V)*sqrt(1-(cos(delta)*sin(chi))^2))/(tan(delta)*sin(chi)))*sin(gamma + phi) - RHS;
figure(2)
hold on
ylim(.5*[-1 1])
plot([-100 100],[0 0])
[derp00 idx] = min(abs(derp333))

GAMMA_search = rad2deg(gamma(idx))
A = ((2*(omega_E*r*cos(delta)/V)*sqrt(1-(cos(delta)*sin(chi))^2))/(tan(delta)*sin(chi)));
B = RHS;

A= -1.1010609116;
 B= 1.0068134103;
phi= -2.5040055779;
derpAAA = (sin(gamma).^2) - A*sin(gamma + phi) - B;
plot(rad2deg(gamma),derpAAA)

A_p = 2*A;
B_p = 2*B-1;

x1 = 0.1;
x_old1 = 0;
iter1 = 0;
while abs(x_old1-x1) > 1e-7
    x_old1 = x1;
        f = sin(x1)^2 - A*sin(x1 + phi) - B;
        f_p = sin(2*x1) - A*cos(x1 + phi);
%     f = cos(2*x1) + A_p*sin(x1 + phi) + B_p;
%     f_p = -2*sin(2*x1) + A_p*cos(x1 + phi);
    x1 = x1 - f/f_p;
    iter1 = iter1 + 1;
end
GAMMA_NewtRaph = rad2deg(wrapToPi(x1))

x2 = 0.1;
x_old2 = 0;
iter2 = 0;
while abs(x_old2-x2) > 1e-7 
    x_old2 = x2;
    f = sin(x2)^2 - A*sin(x2 + phi) - B;
    f_p = sin(2*x2) - A*cos(x2 + phi);
    f_2p = 2*cos(2*x2) + A*sin(x2 + phi);
%     f = cos(2*x2) + A_p*sin(x2 + phi) + B_p;
%     f_p = -2*sin(2*x2) + A_p*cos(x2 + phi);
%     f_2p = -4*cos(2*x2) - A_p*sin(x2 + phi);
    x2 = x2 - 2*f*f_p/(2*f_p^2 - f*f_2p); %Eq. 2 from  http://sci-hub.tw/10.1002/zamm.19540340110
    iter2 = iter2 + 1;
end
GAMMA_Halley = rad2deg(wrapToPi(x2))

x3 = 0;
x_old3 = .1;
iter3 = 0;
while abs(x_old3-x3) > 1e-7 
    x_old3 = x3;
    f = sin(x3)^2 - A*sin(x3 + phi) - B;
    f_p = sin(2*x3) - A*cos(x3 + phi);
    f_2p = 2*cos(2*x3) + A*sin(x3 + phi);
%         f = cos(2*x3) + A_p*sin(x3 + phi) + B_p;
%     f_p = -2*sin(2*x3) + A_p*cos(x3 + phi);
%     f_2p = -4*cos(2*x3) - A_p*sin(x3 + phi);

    x3 = x3 - (f*f_p)/(f_p^2 - f*f_2p); % http://mathworld.wolfram.com/SchroedersMethod.html      https://www-sciencedirect-com.tudelft.idm.oclc.org/science/article/pii/S0377042709006347#b20
    iter3 = iter3 + 1;
end
GAMMA_Schroder = rad2deg(wrapToPi(x3))

x4 = 0;
x_old4 = .1;
iter4 = 0;
while abs(x_old4-x4) > 1e-7 
    x_old4 = x4;
%     f = sin(x4)^2 - A*sin(x4 + phi) - B;
%     f_p = sin(2*x4) - A*cos(x4 + phi);
%     f_2p = 2*cos(2*x4) + A*sin(x4 + phi);
%     f_3p = -4*sin(2*x4) + A*cos(x4 + phi);
        f = cos(2*x4) + A_p*sin(x4 + phi) + B_p;
    f_p = -2*sin(2*x4) + A_p*cos(x4 + phi);
    f_2p = -4*cos(2*x4) - A_p*sin(x4 + phi);
    f_3p = 8*cos(2*x4) - A_p*cos(x4 + phi);

    x4 = x4 - f*(f_p^2 - .5*f*f_2p)/(f_p^3 - f*f_p*f_2p + (1/6)*f^2*f_3p); %Eq. 4 from  http://sci-hub.tw/10.1002/zamm.19540340110
    iter4 = iter4 + 1;
end
GAMMA_3rd_Order = rad2deg(wrapToPi(x4))


