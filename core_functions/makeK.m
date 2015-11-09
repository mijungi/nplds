
function [K, dKdsig2, dKdtau2] = makeK(sig2, epsilon, tau2, r, k, ind_train)

% K = makeK(sig2, epsilon, tau2, r, k);

I_ksqrd = speye(k); % k^2 by k^2 
K = zeros(k*r);
dKdsig2 = zeros(k*r);
dKdtau2 = zeros(k*r);

% complete prior covariance matrix K
for i=1:r
    for j=1:r
        
        ii = ind_train(i);
        jj = ind_train(j);
        
%         [i j ii jj]
        
        if ii==jj
            K((i-1)*k+1:i*k, (j-1)*k+1:j*k) = (sig2+epsilon)*I_ksqrd;
            dKdsig2((i-1)*k+1:i*k, (j-1)*k+1:j*k) = I_ksqrd;
            dKdtau2((i-1)*k+1:i*k, (j-1)*k+1:j*k) = 0*I_ksqrd;
        else
            K((i-1)*k+1:i*k, (j-1)*k+1:j*k) = (sig2)*exp(-0.5*(ii-jj)^2/tau2)*I_ksqrd;
            dKdsig2((i-1)*k+1:i*k, (j-1)*k+1:j*k) = exp(-0.5*(ii-jj)^2/tau2)*I_ksqrd;
            dKdtau2((i-1)*k+1:i*k, (j-1)*k+1:j*k) = (ii-jj)^2/(2*tau2^2)*K((i-1)*k+1:i*k, (j-1)*k+1:j*k);
        end
            
    end
end


%%
% function [K, dKdsig2, dKdtau2] = makeK(sig2, epsilon, tau2, r, k)
% 
% % K = makeK(sig2, epsilon, tau2, r, k);
% 
% I_ksqrd = eye(k); % k^2 by k^2 
% K = zeros(k*r);
% dKdsig2 = zeros(k*r);
% dKdtau2 = zeros(k*r);
% 
% % complete prior covariance matrix K
% for i=1:r
%     for j=1:r
%         if i==j
%             K((i-1)*k+1:i*k, (j-1)*k+1:j*k) = (sig2+epsilon)*I_ksqrd;
%             dKdsig2((i-1)*k+1:i*k, (j-1)*k+1:j*k) = I_ksqrd;
%             dKdtau2((i-1)*k+1:i*k, (j-1)*k+1:j*k) = 0*I_ksqrd;
%         else
%             K((i-1)*k+1:i*k, (j-1)*k+1:j*k) = (sig2)*exp(-0.5*(i-j)^2/tau2)*I_ksqrd;
%             dKdsig2((i-1)*k+1:i*k, (j-1)*k+1:j*k) = exp(-0.5*(i-j)^2/tau2)*I_ksqrd;
%             dKdtau2((i-1)*k+1:i*k, (j-1)*k+1:j*k) = (i-j)^2/(2*tau2^2)*K((i-1)*k+1:i*k, (j-1)*k+1:j*k);
%         end
%     end
% end
% 
