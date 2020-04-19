function [tau_f, lambda_f] = ...
    DMA_getTL(TorL,UorD,a,w_i,w_e,tau_s,Agamma,gamma,guess)
% This function is used to find tau and lambda for each (omega_i, omega_e).
% Calculate lambda first, then tau.
% idx = 'u' for upscan and 'd' for downscan
% TorL = 't' calculate tau and then lambda, 'l' the other way
% -------------------------------------------------------------------------
% Either calculating tau_t first or lambda first. Recommend tau_t first.
% -------------------------------------------------------------------------
% created: 2017/04/13, YH
% comment: 2017/05/20, YH
% -------------------------------------------------------------------------
% fixed the problem that 't' and 'l' give different solutions in 'd', 'l'
% has a typo of log(gamma) instead of log(lambda). 2017/05/22, YH
% -------------------------------------------------------------------------
% fixed the problem that calculated tau is too small and lambda is too big,
% it is because the sign +/- in Li_polylog(2,A*exp(-x)) is wrong.
% last modified: 2017/05/23, YH
% tau is theta in the paper, 2019/08/27, YH

flag = [TorL UorD];
fun_x = DMA_TLfun(flag,a,w_i,w_e,tau_s,Agamma,gamma);

tl_f = fzero(fun_x,guess); % either tau or lambda

switch TorL
    case 't' % calculation based on tau_t
        tau_f = tl_f;
        switch UorD
            case 'u'
                % omega_i = omega_e + lambda*(1-exp(-tau_t));
                lambda_f = (w_i-w_e)/(1-exp(-tau_f));
            case 'd' % downscan
                % omega_i = omega_e - lambda*(1-exp(tau_t));
                lambda_f = (w_e-w_i)/(1-exp(tau_f));
        end
        
    case 'l' % calculation based on lambda
        lambda_f = tl_f;
        switch UorD
            case 'u' % upscan
                % omega_i = omega_e + lambda*(1-exp(-tau_t));
                tau_f = log(lambda_f/(-w_i+w_e+lambda_f));
            case 'd' % downscan
                % omega_i = omega_e - lambda*(1-exp(tau_t));
                tau_f = log((w_i-w_e+lambda_f)/lambda_f);
        end
end

end % end of FUN getTauLambda