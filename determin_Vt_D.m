function [Vt, D] = determin_Vt_D(Dp)
%% method I
% d_a = 1.2; % kg m-3, air density
% d_p = 1e3; % kg m-3, particle density
% u_a = 1.8e-5; % kg m-1 s-1
% k = 8.314/6.02e23; % Boltzmann constant
% T = 298; 
% Re = 10.^(-13:0.001:4);
% Ga0 = Drag_coeff(Re).*Re.^2; % loglog(Re,Ga0);
% [Cc,~] = Cunningham(Dp);
% D = k*T.*Cc/3/pi/u_a./Dp;
% Ga1 = 4.*Dp.^3*d_a*d_p*9.8.*Cc/3/u_a^2;
% Vt = zeros(length(Dp),1);
% for i = 1:length(Dp)
%     x = find(Ga1(i) < Ga0,1);
%     r = mean(Re(x-1:x));
%     Vt(i) = r*u_a/d_a/Dp(i);
% end

%% method II
d_a = 1.2; % kg m-3, air density
d_p = 1.4e3; % kg m-3, particle density
u_a = 1.8e-5; % kg m-1 s-1, air viscosity
k = 8.314/6.02e23; % Boltzmann constant
T = 298;
[Cc,~] = Cunningham(Dp);
D = k*T.*Cc/3/pi/u_a./Dp;
Ga = 4.*Dp.^3*d_a*d_p*9.8.*Cc/3/u_a^2;
Re = zeros(length(Dp),1);
for i = 1:length(Dp)
    fRe = @(Re) Drag_coeff(Re).*Re.^2 - Ga(i);
    Re(i) = fzero(fRe,1);
end
Vt = Re*u_a/d_a./Dp;

return

