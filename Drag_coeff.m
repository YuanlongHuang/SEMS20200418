function Cd = Drag_coeff(Re)

% if Re < 0
%     fprintf('Something is wrong, Reynolds number is negative!\n');
% end

if Re < 0.1
    Cd = 24 ./ Re;
elseif 0.1 <= Re < 0.7
    Cd = 24 ./ Re .* (1+3/16*Re+9/160*Re.^2.*log(2*Re));
elseif 0.7 <= Re < 1000
    Cd = 24 ./ Re .* (1+0.15*Re.^0.687);
else
    Cd = 0.44;
end

return