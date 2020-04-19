function [F_ea,z_a] = DMA_extFe(Frame,DMAinfo)
% this function returns the extended F_e based on (a,zeta,F_i)
% use spline to extrapolation F_e
% created 2020/01/07, yh

grid_i = length(Frame.F_i);
grid_e = length(Frame.F_e);

Fi = repmat(Frame.F_i,grid_e,1);
Fx = Frame.F_e;
zx = DMAinfo.zeta(grid_e:(2*grid_e-1),:); % grid_e x grid_i
za = DMAinfo.zeta((2*grid_e-1):end,:);

Fe = zeros(grid_e,grid_i); % the extended F_e
for i = 1:grid_e
    for j = 1:grid_i % for every F_i
        Fe(i,j) = interp1(zx(:,j),Fx,za(i,j),'linear','extrap');
    end
end

Fe_ = repmat(linspace(min(Fe(1,:)),min(Fe(end,:)),grid_e)',...
    1,grid_i); % new grid
zz = griddata(Fi(:),Fe(:),za(:),Fi(:),Fe_(:)); % new mobility

z_a = reshape(zz,grid_e,grid_i);
F_ea = Fe_(:,1);

end
