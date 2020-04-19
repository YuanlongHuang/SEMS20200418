function [TF,TF_all] = DMA_getTF(flow,Frame,DMAinfo,ORZ,ODN)
% This function gives transfer function from Z_i and sigma.
% It integrates based on flow fraction frame.
% -------------------------------------------------------------------------
% The idea is that on contour of Z_i:
% First, calculate the possibility of each pair of (f_i,f_e) w.r.t. f_e
% [1-beta*(1+delta)/(1+beta), 1];
% Second, integrate the possibility that isoline spans over f_i [0,
% beta*(1-delta)/(1+beta)].
% ORZ: indicates if calculation is based on dimensionless or absolute zeta
%      'r' (relative) or 'z' (absolute)
% ODN: indicates if diffusion is included, 'd' (diff) or 'n' (non-diff)
% -------------------------------------------------------------------------
% NOTE: the out flow fraction has to be multiplied by a factor of 1+beta to
% keep consistent with sigma, since sigma is scaled by (Qsh+Qex)/4/pi, and
% flow fraction has to time QT/2/pi to get streamline and then scale with
% (Qsh+Qex)/4/pi, which finally gives the factor 1+beta.
% -------------------------------------------------------------------------
% created: 2017/05/25, YH
% comment: 2017/05/29, YH
% working: 2017/05/30, YH
% add diffusion option, 2017/06/01, YH

%% pre-calculate all the probability
Qa = flow(1); % aerosol inlet flow, m3 s-1
Qc = flow(2); % classified outlet flow, m3 s-1
Qsh = flow(3); % sheath flow, m3 s-1
Qex = flow(4); % excess flow, m3 s-1
beta = (Qa+Qc)/(Qsh+Qex);
delta = (Qc-Qa)/(Qc+Qa);

switch ORZ
    case 'r' % relative to Zstar, return dimensionless zeta
        zeta = DMAinfo.zeta;
    case 'z' % absolute e-mobility value, m2 V-1 s-1 (C kg-1 s)
        zeta = DMAinfo.Z_i;
end
sigma = DMAinfo.sigma; % dimensionless expanded-distribution
F_e = Frame.F_e*(1+beta); % outlet flow fraction * (1+beta)
grid_i = length(Frame.F_i);
grid_e = length(F_e);
B_e_lower = 1-beta*delta; % lower boundary of outlet
B_e_upper = 1+beta; % upper boundary of outlet
Pout = 1/2*(erf((B_e_upper-repmat(F_e',1,grid_i))/sqrt(2)./sigma)...
          - erf((B_e_lower-repmat(F_e',1,grid_i))/sqrt(2)./sigma));
% probability at the outlet
Peval = scatteredInterpolant(repelem(Frame.F_i,grid_e)',...
                             repmat(Frame.F_e,1,grid_i)',...
                             Pout(:));
% set the interpolant Peval to evaluate the possibility at given points
range_i = beta*(1-delta)/(1+beta); % streamline range of inlet

%% calculate transfer function
% zeta = DMAinfo.zeta; % dimensionless e-mobility matrix
% zi = linspace(min(min(zeta)),max(max(zeta)),grid_n);
% C = contours(Frame.F_i,Frame.F_e,DMAinfo_up.zeta,[zi(i) zi(i)]);
% Gives contour line without plotting, but may be removed in future.
% Use contourc instead.
grid_nn = grid_e;
C = contourc(Frame.F_i,Frame.F_e,zeta,grid_nn);
TF_all = cell(grid_nn+2,3); % 1. zeta; 2. omega; 3. x,y,P
TF_all{1,1} = min(min(zeta));
TF_all{1,3} = [];
TF_all{end,1} = max(max(zeta));
TF_all{end,3} = [];
TF_all(:,2) = {0};
idx_1 = 1;
i = 2;
while idx_1 <= size(C,2)
    zi = C(1,idx_1);
    if zi == TF_all{i-1,1}
        i = i-1;
    end
    TF_all{i,1} = zi; % get zeta
    
    idx_2 = idx_1 + C(2,idx_1);
    x = C(1,(idx_1+1):idx_2);
    y = C(2,(idx_1+1):idx_2);
    switch ODN
        case 'd'
            P = Peval(x,y);
        case 'n'
            P = ones(size(x));
    end
    xyP = [x; y; P];
    TF_all{i,3} = [TF_all{i,3} xyP]; % get x, y, P
    idx_1 = idx_2 + 1;
    
    i_max = find(x == max(x),1);
    if i_max == length(x) || i_max == 1
        omega = trapz(x,P);
    else
        omega1 = trapz(x(1:i_max),P(1:i_max));
        omega2 = -trapz(x(i_max:end),P(i_max:end));
        omega = omega1 + omega2;
    end
    TF_all{i,2} = TF_all{i,2} + abs(omega)/range_i;
    % get integrated transmission efficiency
    
    i = i + 1;
    
end
TF = cell2mat(TF_all(:,1:2));

end