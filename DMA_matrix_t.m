function [theta_t,lambda] = DMA_matrix_t(volt,flow,cnfg,idx,Frame,don)
% This function returns the matrices of the time and e-mobility of particle
% exiting DMA at time t_e, upscan or downscan.
% However, it only considers the period that voltage changes expotentially,
% no boundary condition is considered.
% It returns both extended and original F space
% -------------------------------------------------------------------------
% tau_t: dimensionless residence matrix
% lambda: dimensionless e-mobility matrix
% -------------------------------------------------------------------------
% volt: electric condition of DMA
% flow: flow condition of DMA
% cnfg: configuration of DMA
% idx: 'u' for upscan and 'd' for downscan, or 's' indicates static mode
% grid_n: grid points in w_i-w_e space, e.g. = 201;
% -------------------------------------------------------------------------
% created: 2017/04/13, YH
% comment: 2017/05/20, YH
% modified: add static mode calculation. set idx = 's'. 2017/05/23, YH
% switch from Omega-frame to FlowFrac-frame. 2017/05/28, YH
% modified: add calculation of all streamlines, 2017/06/03, YH
%% characterize parameters-------------------%
Qa = flow(1); % aerosol inlet flow, m3 s-1
Qc = flow(2); % classified outlet flow, m3 s-1
Qsh = flow(3); % sheath flow, m3 s-1
Qex = flow(4); % excess flow, m3 s-1

L = cnfg(1); % m, length of column
r2 = cnfg(2); % m, outer radius
r1 = cnfg(3); % m, inner radius, from B&W pg 559
gamma = (r1/r2)^2;
Agamma = (-1/2*(1+gamma)*log(gamma)-(1-gamma))^-1;
tau_m = 2*pi*(r2^2-r1^2)*L/(Qsh+Qex+Qa+Qc); % average residence time, s
% Ubar = L/tau_r; % average velocity, m s-1
u_w = @(x) Agamma*((1-gamma)*(x.*log(x)-x)-log(gamma)*(x-x.^2/2));
% dimensionless average velocity = (u_w(we)-u_w(wi))/(we-wi)

t_ramp = volt(1); % ramp up time, s (e.g. = 300)
Vmin = volt(2); % minimum voltage, V (e.g. = 15)
Vmax = volt(3); % maximum voltage, V (e.g. = 9850)
t_s = t_ramp/log(Vmax/Vmin); % s, time constant
tau_s = t_s/tau_m; % dimensionless time constant

%% main part - matrices calculation
omega_i = Frame.omega_i; % from large to small
switch don.ODN
    case 'd'
        omega_e = Frame.omega_e_d; % from large to small
        a = Frame.a;
        grid_n = size(don.gl,2);
    case 'n'
        omega_e = Frame.omega_e; % from large to small
        a = ones(1,length(omega_e));
end
grid_i = length(omega_i);
grid_e = length(omega_e);

theta_t = zeros(grid_i,grid_e);
lambda = zeros(grid_i,grid_e);

if strcmp(idx,'s') % calculate static mode
    wi = repmat(omega_i',1,grid_e); % matrix of omega_i
    we = repmat(omega_e,grid_i,1); % matrix of omega_e
    aa = repmat(a,grid_i,1); % matrix of a
    lambda = tau_s*Agamma*((1-gamma)*(wi.*(log(wi)-1)-...
        we.*(log(we)-1))-(wi-we).*(1-(wi+we)/2)*log(gamma))./aa;
    theta_t = (wi-we)./lambda;
    fprintf('Static mode calculation is done.\n');
else % calculate scan mode
    switch don.ODN
        case 'n'
            for i = 1:grid_i
                for j = 1:grid_e
                    if i == 1 && j == 1
                        ubar = (u_w(omega_e(j))-u_w(omega_i(i)))/...
                            (omega_e(j)-omega_i(i));
                        guess = tau_m/ubar/t_s;
                    elseif j == 1
                        guess = theta_t(i-1,j); % lambda(i-1,j);
                    else
                        guess = theta_t(i,j-1); % lambda(i,j-1);
                    end
                    [theta_t(i,j), lambda(i,j)] = DMA_getTL('t',idx,a(j),...
                        omega_i(i),omega_e(j),tau_s,Agamma,gamma,guess);
                end
                fprintf([idx ', %i (%i) is done.\n'],i,grid_i);
            end
        case 'd'
            lambda(:,grid_n:(2*grid_n-1)) = don.gl;
            theta_t(:,grid_n:(2*grid_n-1)) = don.gt;
            for j = (grid_n-1):-1:1 % 34:-1:1 % 
                for i = 1:grid_i
%                     guess = lambda(i,j+1);
                    guess = theta_t(i,j+1);
                    if isnan(guess) || isinf(guess)
                        guess = 0.01;
                    end
                    [theta_t(i,j), lambda(i,j)] = DMA_getTL('t',idx,a(j),...
                        omega_i(i),omega_e(j),tau_s,Agamma,gamma,guess);
                end
                fprintf([idx ', Part I, %i (1) is done.\n'],j);
            end
            for j = 2*grid_n:(3*grid_n-2)
                for i = 1:grid_i
%                     guess = lambda(i,j-1);
                    guess = theta_t(i,j-1);
                    [theta_t(i,j), lambda(i,j)] = DMA_getTL('t',idx,a(j),...
                        omega_i(i),omega_e(j),tau_s,Agamma,gamma,guess);
                end
                fprintf([idx ', Part II, %i (%i) is done.\n'],j,3*grid_n-2);
            end
    end
end
end % end of FUN DMA_run