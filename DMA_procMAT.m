function DMAinfo = DMA_procMAT(volt,flow,cnfg,idx,Frame,tau_t,lambda,...
    t_e,charge,T,p,don)
% post-process matrices from DMA_matrix
% -------------------------------------------------------------------------
% DMAinfo is a struct: t_i, t_r, Z_i, zeta, Zstar, sigma
% t_i: time matrix at the inlet, s
% t_r: residence time matrix of the trajectories, s
% Z_i: e-mobility matrix corresponding to t_e, m2 V-1 s-1 (C kg-1 s)
% zeta: dimensionless e-mobility matrix
% Zstar: reference e-mobility corresponding to t_e, m2 V-1 s-1 (C kg-1 s)
% sigma: diffusivity caused expanded-distribution
% -------------------------------------------------------------------------
% t_e: the time that particles exit DMA, s
% charge: electron charge
% T = 298; % K
% -------------------------------------------------------------------------
% created: 2017/05/12, YH
% comment: 2017/05/24, YH
% switch from Omega-frame to FlowFrac-frame. 2017/05/28, YH
%% characterize parameters-------------------%
Qa = flow(1); % aerosol inlet flow, m3 s-1
Qc = flow(2); % classified outlet flow, m3 s-1
Qsh = flow(3); % sheath flow, m3 s-1
Qex = flow(4); % excess flow, m3 s-1
beta = (Qa+Qc)/(Qsh+Qex);
delta = (Qc-Qa)/(Qc+Qa);

L = cnfg(1); % m, length of column
r2 = cnfg(2); % m, outer radius
r1 = cnfg(3); % m, inner radius, from B&W pg 559
gamma = (r1/r2)^2;
Agamma = (-1/2*(1+gamma)*log(gamma)-(1-gamma))^-1;
tau_m = 2*pi*(r2^2-r1^2)*L/(Qsh+Qex+Qa+Qc); % average residence time, s
% Ubar = L/tau_r; % average velocity, m s-1
Vzstar = (Qsh+Qex)/4/pi/L/r2; % z-direction reference velocity, m s-1

t_ramp = volt(1); % ramp up time, s (e.g. = 300)
Vmin = volt(2); % minimum voltage, V (e.g. = 15)
Vmax = volt(3); %maximum voltage, V (e.g. = 9850)
t_s = t_ramp/log(Vmax/Vmin); % time constant, s
tau_s = t_s/tau_m; % dimensionless time constant

kB = 1.3806488e-23; % m2 kg s-2 K-1, Boltzmann constant
eC = 1.60217646e-19; % C, electron charge in coulomb

%% transform dimensionless matrices
switch idx
    case 'u' % upscan
        V_e = Vmin*exp(t_e/t_s);
    case 'd' % downscan
        V_e = Vmax*exp(-t_e/t_s);
    case 's' % static mode, same as upscan
        V_e = Vmin*exp(t_e/t_s);
end
Zstar = (Qsh+Qex)/4/pi/L/V_e*log(r2/r1);
% reference e-mobility, m2 V-1 s-1 (C kg-1 s)
tau_e = t_e/tau_m;

% tau_t = -(tau-tau_e)/tau_s;
tau = tau_e - tau_t*tau_s; % intermedia matrix
% tau = t/tau_r;
t_i = tau*tau_m; % time that particle enters DMA, s
t_r = t_e - t_i; % residence time, s
% lambda = zeta*(1-gamma)/(1+beta)*tau_s;
zeta = lambda/tau_s*(1+beta)/(1-gamma); % dimensionless e-mobility
% zeta = Z/Zstar;
Z_i = zeta*Zstar; % e-mobility, m2 V-1 s-1 (C kg-1 s)

%% get dimensionless matrix of sigma
omega_i = Frame.omega_i;
switch don.ODN
    case 'd'
        omega_e = Frame.omega_e_d; % from large to small
    case 'n'
        omega_e = Frame.omega_e; % from large to small
end
grid_i = length(omega_i);
grid_e = length(omega_e);

Dstar = kB*T/charge/eC*Zstar; % reference diffusivity, m2 s-1
Pestar = Vzstar*(r2-r1)/Dstar; % z-direction migration Peclet number
if strcmp(idx,'s')
    Coeff = 4*((1+beta)/(1-gamma))^2*(1-sqrt(gamma))/Pestar;
else
    Coeff = 4*((1+beta)/(1-gamma))^2*(1-sqrt(gamma))/Pestar*lambda;
end

Iz = zeros(size(tau_t));
switch idx
    case 's' % static
        wi = repmat(omega_i',1,grid_e); % matrix of omega_i
        we = repmat(omega_e,grid_i,1); % matrix of omega_e
        Ir = lambda.^2/4/tau_s^2*(r2/L)^2.*(wi-we);
        f_Iz = @(w) Agamma^2*...
            (1/2*w.^2.*((1-gamma)*log(w)-(1-w)*log(gamma)).^2 ...
            -(1/2*w.^2*(1-gamma)+1/3*w.^3*log(gamma)).* ...
            ((1-gamma)*log(w)-(1-w)*log(gamma)) ...
            +1/4*(1-gamma)^2*w.^2 ...
            +5/18*(1-gamma)*log(gamma)*w.^3 ...
            +1/12*log(gamma)^2*w.^4); % analytical sulution
        Iz = f_Iz(wi) - f_Iz(we);
        fprintf('Static Iz is done by analytical solution.\n');
        Stlzbg.Ir = (1/2*(1-gamma)/(1+beta)*r2/L)^2;
        Stlzbg.Iz = (f_Iz(1)-f_Iz(gamma))/(1-gamma);
        Stlzbg.G = 4*(1+beta)^2/(1-gamma)*(Stlzbg.Ir+Stlzbg.Iz);
        Stlzbg.Dc = (1-sqrt(gamma))/Pestar;
        % Stolzenburg form of sigma^2 = G*Dc*zeta,
        % where zeta = Z_i/Zstar
        Stlzbg.beta = beta;
        Stlzbg.delta = delta;
    case 'u' % upscan
        switch don.ODN
            case 'd'
                Ir = lambda.^2/4/tau_s^2*(r2/L)^2*1/2.*(1-exp(-2*tau_t));
                for i = 1:grid_i
                    for j = 1:grid_e
                        f_Iz = @(t) (omega_e(j)+lambda(i,j)*(1-exp(-t))).*...
                            (Agamma*((1-gamma)*log(omega_e(j)+lambda(i,j)*...
                            (1-exp(-t)))...
                            -(1-omega_e(j)-lambda(i,j)*...
                            (1-exp(-t)))*log(gamma))).^2;
                        Iz(i,j) = integral(f_Iz,0,tau_t(i,j));
                    end
                    fprintf('u, calculating Iz %i (%i).\n',i,grid_i);
                end
            case 'n'
                Ir = Iz;
        end
    case 'd' % downscan
        switch don.ODN
            case 'd'
                Ir = lambda.^2/4/tau_s^2*(r2/L)^2*1/2.*(exp(2*tau_t)-1);
                for i = 1:grid_i
                    for j = 1:grid_e
                        f_Iz = @(t) (omega_e(j)-lambda(i,j)*(1-exp(t))).*...
                            (Agamma*((1-gamma)*log(omega_e(j)-lambda(i,j)*...
                            (1-exp(t)))...
                            -(1-omega_e(j)+lambda(i,j)*...
                            (1-exp(t)))*log(gamma))).^2;
                        Iz(i,j) = integral(f_Iz,0,tau_t(i,j));
                    end
                    fprintf('d, calculating Iz %i (%i).\n',i,grid_i);
                end
            case 'n'
                Ir = Iz;
        end
end

sigma2 = Coeff.*(Ir+Iz);

%% write into a struct and tranform matrices so that horizontal is inlet
DMAinfo.t_i = t_i';
DMAinfo.t_r = t_r';
DMAinfo.Z_i = Z_i';
DMAinfo.zeta = zeta';
DMAinfo.sigma = sqrt(sigma2)';
DMAinfo.Zstar = Zstar;
DMAinfo.Dpstar = DMA_getDp(Zstar,charge,T,p); % m, particle diameter
DMAinfo.Dstar = Dstar;
DMAinfo.Pestar = Pestar;
if strcmp(idx,'s')
    DMAinfo.Stlzbg = Stlzbg;
end

end