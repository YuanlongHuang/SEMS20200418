function S = Li_polylog(n,x)
% Note: MATLAB 2014b have built in function polylog(n,x)
% However, function Li_polylog(n,x) runs much faster than polylog(n,x)
% Li_polylog(n,x) is polylogarithm function Li_n(x) of index n at point x
% sum(k=1->inf,x^k/k^n) for |x|<1
% -------------------------------------------------------------------------
% NOTE: this function is used just for |x|<1
% -------------------------------------------------------------------------
% created: 2015/08/05, YH
% comment: 2017/05/20, YH

k = 100; % initially set the series to be 100 terms
i = (1:k)';
v1 = x.^i;
v2 = i.^(-n);
S = dot(v1,v2);
while abs(x.^k/k^n)/S > 1e-12 % compare the relative error
    k = k*10;
    i = (1:k)';
    v1 = x.^i;
    v2 = i.^(-n);
    S = dot(v1,v2);
end

% S = 0;
% for i = 1:k
%     S = S + x.^i/i^n;
% end

return
