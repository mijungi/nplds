function Cnew = VBcomputeC_h_NSFR(xyzinpn, r, fromVBEstep, params, opts2)

%% unpack init params

C = params.C;
d = params.d;
h = params.h;

%% numerically compute the mean of C

[p, k] = size(C);

prs0 = C(:); % vectorized C
fun = @(prs) (loglikeliFunc(prs, h, d, fromVBEstep, p, k, r, xyzinpn));
%             loglikeliFunc(prs, h, d, fromVBEstep, p, k, r, xyzinpn)
% prs = fsolve(fun, prs0, opts2);
prs = fminunc(fun, prs0, opts2);

Cnew = reshape(prs, p, k);

% ===========================================================
function [l, dl] = loglikeliFunc(prs, htot, d, fromVBEstep, p, k, r, xyzinpn)

%% (1) reshape C
C = reshape(prs, p, k);

%% 2) compute dlogq(c)/dc
% dlogq(c)/dc = SCtrm - omegaTrm - UpsilonTrm - priorTrm

l = zeros(r,1);
dl = zeros(p*k, r); 

for i=1:r
    
    %% unpack mu/cov of latent variables
    
    if r==1
        mumarg = fromVBEstep.mumarg;
        inv_sigmarg = fromVBEstep.inv_sigmarg;
        suffstat = fromVBEstep.suffstat;
        y = xyzinpn.y;
    else
        mumarg = fromVBEstep{i}.mumarg;
        inv_sigmarg = fromVBEstep{i}.inv_sigmarg;
        suffstat = fromVBEstep{i}.suffstat;
        y = xyzinpn{i}.y;
    end
    
    h = htot(:,i);
    
    T = size(mumarg,2);
    SC = suffstat.SC;
    
    % obj = sum_t( y_t'*(C*(w_t + h) + d) - 1'*sum(exp(g_t)) );
    % where g_t = C*(w_t+h) + 0.5*diag(C*Upsilon_t*C') + d
    
    obj_frstTrm = zeros(T,1);
    obj_scndTrm = zeros(T,1);
    
    % d obj / dC
    % = SC' +  y_t*h' - sum( exp(g_t)* w_t') - sum(diag(exp(g_t))*C*Upsilon_t)
    
    der_scndTrm = zeros(p, k, T);
    der_thrdTrm = zeros(p, k, T);
    der_fourthTrm = zeros(p, k, T);
    
    % highThrsh = 1e6;
    
    for t=1:T
        CUpsilon = C/inv_sigmarg(:,:,t);
        f = diag(CUpsilon*C');
        g = C*(mumarg(:,t)+h) + 0.5*f + d;
        
        tmp = exp(g);
        %     tmp(isinf(tmp)) = highThrsh;
        %     tmp(tmp>highThrsh) = highThrsh;
        
        obj_frstTrm(t) = y(:,t)'*(C*(mumarg(:,t)+h)+d);
        obj_scndTrm(t) = sum(tmp);
        
        der_scndTrm(:,:,t) = y(:,t)*h';
        der_thrdTrm(:,:,t) = tmp*(mumarg(:,t)'+h');
        der_fourthTrm(:,:,t) = diag(tmp)*CUpsilon;
    end
    
    gamma = 0;
    
    obj = sum(obj_frstTrm - obj_scndTrm) - gamma*sum(diag(C*C'));
    l(i) = - obj;
    
    d_obj = SC' + sum(der_scndTrm,3) - sum(der_thrdTrm, 3) - sum(der_fourthTrm, 3) - gamma*C;
    dl(:,i) = - d_obj(:);
    
end

l = sum(l);
dl = sum(dl,2);


