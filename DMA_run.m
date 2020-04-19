%% Flows, Scan time, Voltages, & DMA can be changed in this section
clear; clc;
% ====================== DMA parameters ====================== %
Qsh = 2.5; % lpm
Qa = 0.25; % lpm
Qc = Qa; % lpm
% Qa is aerosol inlet flow, Qc is classified outlet flow
Qex = Qsh + Qa - Qc; % lpm
% Qsh is sheath flow, Qex is excess flow

% ====================== Scan time here ====================== %
t_ramp = 45; % ramp time, s

% ======================  Voltage here  ====================== %
Vmin = 15; % minimum voltage, V
Vmax = 9850; % maximum voltage, V
V_e = 80; % voltage at the exit time t_e, corresponding to a specific Z+
t_e = t_ramp*log(V_e/Vmin)/log(Vmax/Vmin); % time at V_e

% ====================== TSI column DMA ====================== %
L = 44.369e-2; % m, length of column
r2 = 1.961e-2; % m, outer radius
r1 = 0.937e-2; % m, inner radius, from B&W pg 559

ne = 1; % charge number
T = 298; % K, temperature
p = 745; % mmHg, room pressure

% ====================== Nondimensional ====================== %
gamma = (r1/r2)^2;
beta = (Qa+Qc)/(Qsh+Qex);
tau_m = 2e3*60*pi*(r2^2-r1^2)*L/(Qsh+Qex+Qa+Qc); % mean residence time, s
t_s = t_ramp/log(Vmax/Vmin); % time constant, s
tau_s = t_s/tau_m; % dimensionless time constant

%% compile the above conditions into Frame
volt = [t_ramp; Vmin; Vmax]; % electric condition
flow = [Qa; Qc; Qsh; Qex]*1e-3/60; % flow condition, m3 s-1
cnfg = [L; r2; r1]; % configuration, m
grid_i = 201; % grid point in w_i space
grid_e = 201; % grid point in w_e space

Frame = DMA_frame(flow,cnfg,grid_i,grid_e); % get streamline frame
% linear in flow fraction, nonlinear in omega (need to solve)

%% Obtain TF for non-diffusive case
% This subsection only need to run ONCE, since tau_t and lambda is already
% in general form. This is for non-diffusive case
don.ODN = 'n'; % 'd' diffusion; 'n' non-diffusion
don.gl = []; don.gt = [];

% calculate dimensionless matrices of residence time and e-mobility
idx = 'u'; % upscan
[x.theta_t_up,x.lambda_up] = DMA_matrix_t(volt,flow,cnfg,idx,Frame,don);
x.zeta_up = x.lambda_up/tau_s*(1+beta)/(1-gamma);
x.tr_up = x.theta_t_up*tau_s; % normalized RT, relative to tau_m

idx = 's'; % static
[x.theta_t_st,x.lambda_st] = DMA_matrix_t(volt,flow,cnfg,idx,Frame,don);
x.zeta_st = x.lambda_st/tau_s*(1+beta)/(1-gamma);
x.tr_st = x.theta_t_st*tau_s; % normalized RT, relative to tau_m

% post-process of matrices at V_e
x.DMAinfo_up = DMA_procMAT...
    (volt,flow,cnfg,'u',Frame,x.theta_t_up,x.lambda_up,t_e,ne,T,p,don);

x.DMAinfo_st = DMA_procMAT...
    (volt,flow,cnfg,'s',Frame,x.theta_t_st,x.lambda_st,t_e,ne,T,p,don);

% obtain TF for non-diffusionve cases
ORZ = 'r'; % 'r' relative to Zstar; 'z' absolute e-mobility
ODN = 'n'; % 'n' non-diffusion; 'd' diffusion
[TF.uN,~] = DMA_getTF(flow,Frame,x.DMAinfo_up,ORZ,ODN);
[TF.sN,~] = DMA_getTF(flow,Frame,x.DMAinfo_st,ORZ,ODN);

%% Obtain TF for diffusive case based on non-diffusive info
% This section will take a little longer time, be patient
don.ODN = 'd'; % 'd' diffusion; 'n' non-diffusion
don.gl = x.lambda_up; don.gt = x.theta_t_up;

% calculate residence time and e-mobility in the lower extended space
idx = 'u'; % upscan
[theta_t_up,lambda_up] = DMA_matrix_t(volt,flow,cnfg,idx,Frame,don);
% theta_t_up and lambda_up could have NaN values
y.theta_t_up = theta_t_up;
y.lambda_up = lambda_up;
y.zeta_up = y.lambda_up/tau_s*(1+beta)/(1-gamma);
y.tr_up = y.theta_t_up*tau_s; % normalized RT, relative to tau_m

idx = 's'; % static
[y.theta_t_st,y.lambda_st] = DMA_matrix_t(volt,flow,cnfg,idx,Frame,don);
y.zeta_st = y.lambda_st/tau_s*(1+beta)/(1-gamma);
y.tr_st = y.theta_t_st*tau_s; % normalized RT, relative to tau_m

% post-process of matrices at V_e
DMAinfo_up = DMA_procMAT...
    (volt,flow,cnfg,'u',Frame,y.theta_t_up,y.lambda_up,t_e,ne,T,p,don);

DMAinfo_st = DMA_procMAT...
    (volt,flow,cnfg,'s',Frame,y.theta_t_st,y.lambda_st,t_e,ne,T,p,don);

ODN = 'd'; % 'd' diffusion; 'n' non-diffusion
% calculate TF in the upper extended space
[F_ea,z_a] = DMA_extFe(Frame,DMAinfo_st);
Frame_st = Frame;
Frame_st.F_e = Frame_st.F_e_d;
Frame_st.F_e((2*grid_e-1):end) = F_ea;
y.DMAinfo_st = DMAinfo_st;
y.DMAinfo_st.zeta((2*grid_e-1):end,:) = z_a;
[TF.sD,~] = DMA_getTF(flow,Frame_st,y.DMAinfo_st,ORZ,ODN);

Frame_up = Frame;
Frame_up.F_e = Frame_up.F_e_d;
Frame_up.F_e((2*grid_e-1):end) = (F_ea-min(F_ea))*0.3+min(F_ea);
y.DMAinfo_up = DMAinfo_up;
[TF.uD,~] = DMA_getTF(flow,Frame_up,y.DMAinfo_up,ORZ,ODN);

%% Stolzenburg TF for static DMA
Stlzbg = x.DMAinfo_st.Stlzbg;
G = Stlzbg.G;
Dc = Stlzbg.Dc;
beta = Stlzbg.beta;
delta = Stlzbg.delta;
z_tilde = (0:0.0001:2)';

TF.sS(:,1) = z_tilde;
% non-diffusive form
TF.sS(:,2) = 1/2/beta/(1-delta)...
    *(abs(z_tilde-(1+beta))...
    + abs(z_tilde-(1-beta)) ...
    - abs(z_tilde-(1+beta*delta))...
    - abs(z_tilde-(1-beta*delta)));
% diffusive form
sigma = sqrt(G*Dc*z_tilde);
Int_err = @(x) x.*erf(x)+1/sqrt(pi)*exp(-x.^2);
TF.sS(:,3) = sigma/sqrt(2)/beta/(1-delta)...
    .*(Int_err((z_tilde-(1+beta))/sqrt(2)./sigma)...
     + Int_err((z_tilde-(1-beta))/sqrt(2)./sigma)...
     - Int_err((z_tilde-(1+beta*delta))/sqrt(2)./sigma)...
     - Int_err((z_tilde-(1-beta*delta))/sqrt(2)./sigma));
% clear Stlzbg G Dc delta z_tilde sigma Int_err;

%% Plots
figure(1)
plot(TF.sS(:,1),TF.sS(:,2),'r-','linewidth',2); hold on;
plot(TF.sN(:,1),TF.sN(:,2),'g--','linewidth',2);
plot(TF.sS(:,1),TF.sS(:,3),'m-','linewidth',2);
plot(TF.sD(:,1),TF.sD(:,2),'b--','linewidth',2);
plot(TF.uN(:,1),TF.uN(:,2),'c-','linewidth',2);
plot(TF.uD(:,1),TF.uD(:,2),'k--','linewidth',2); hold off;

legend('SS w/o D','ST w/o D','SS w/ D','ST w/ D','UP w/o D','UP w/ D');

xlabel('\zeta = Z/Z^+');
ylabel('Transfer Function \Omega')
