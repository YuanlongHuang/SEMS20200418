function Frame = DMA_frame(flow,cnfg,grid_i,grid_e)
% This function transforms flow fraction reference to omega reference
% created: 2017/05/28, YH
% modified: add calculation of all streamlines, 2017/06/03, YH

%% characterize parameters-------------------%
Qa = flow(1); % aerosol inlet flow, m3 s-1
Qc = flow(2); % classified outlet flow, m3 s-1
Qsh = flow(3); % sheath flow, m3 s-1
Qex = flow(4); % excess flow, m3 s-1
beta = (Qa+Qc)/(Qsh+Qex);
delta = (Qc-Qa)/(Qc+Qa);

r2 = cnfg(2); % m, outer radius
r1 = cnfg(3); % m, inner radius, from B&W pg 559
gamma = (r1/r2)^2;
Agamma = (-1/2*(1+gamma)*log(gamma)-(1-gamma))^-1;

%% solve for fluid field boundaries -------------------%
Frac = @(w) -Agamma*(log(gamma)/(2*(1-gamma))*(1-w)^2+w*log(w)+(1-w));
w_i = fzero(@(w) Frac(w)-beta*(1-delta)/(1+beta),[gamma 1]);
% down boundary of inlet flow, up limit is 1
w_i_c = fzero(@(w) Frac(w)-beta*(1-delta)/(1+beta)/2,[gamma 1]);
% central boundary of inlet flow
w_e = fzero(@(w) Frac(w)-1+beta*(1+delta)/(1+beta),[gamma 1]);
% up boundary of outlet flow, down limit is gamma
w_e_c = fzero(@(w) Frac(w)-1+beta*(1+delta)/(1+beta)/2,[gamma 1]);
% central boundary of outlet flow

%% transform flow fraction to omega
F_i = linspace(0,beta*(1-delta)/(1+beta),grid_i);
F_e = linspace(1-beta*(1+delta)/(1+beta),1,grid_e);

F_e_ = linspace(beta*(1-delta)/(1+beta),1-beta*(1+delta)/(1+beta),grid_e);
F_e_d = [F_e_(1:end-1), F_e, ones(1,grid_e-1)]; % extended flow fraction

a_ = linspace(0.5,1,grid_e); % extended length correction
% a_ = linspace(1-(r2-r1)/cnfg(1),1,grid_e); % extended length correction
a = [ones(1,2*grid_e-1), fliplr(a_(1:end-1))]; % correction factor

% this variable considers the effect from diffusion in future calculation
omega_i = zeros(1,grid_i);
omega_e = zeros(1,grid_e);
omega_e_d = zeros(1,length(F_e_d));
for ii = 1:grid_i
    omega_i(ii) = fzero(@(w) Frac(w)-F_i(ii),[gamma 1]);
end
for ii = 1:grid_e
    omega_e(ii) = fzero(@(w) Frac(w)-F_e(ii),[gamma 1]);
end
for ii = 1:length(F_e_d)
    if ii < 2*grid_e
        omega_e_d(ii) = fzero(@(w) Frac(w)-F_e_d(ii),[gamma 1]);
    else
        omega_e_d(ii) = gamma;
    end
end
%% write into a struct
Frame.F_i = F_i;
Frame.omega_i = omega_i;
Frame.F_e = F_e;
Frame.omega_e = omega_e;
Frame.F_e_d = F_e_d;
Frame.omega_e_d = omega_e_d;
Frame.a = a;
Frame.w_i = w_i;
Frame.w_e = w_e;
Frame.w_i_c = w_i_c;
Frame.w_e_c = w_e_c;

end