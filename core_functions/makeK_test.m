function Kstar = makeK_test(sig2, epsilon, tau2, rtrain, k2, rtest, ind_test, ind_train)

% K = makeK(sig2, epsilon, tau2, r, k2);

I_ksqrd = eye(k2); % k^2 by k^2 
Kstar = zeros(k2*rtest, k2*rtrain);

% complete prior covariance matrix K
for i=1:rtest
    for j=1:rtrain
        
        ii = ind_test(i);
        jj = ind_train(j);
        
%         [i j ii jj]
        
        if ii==jj
            Kstar((i-1)*k2+1:i*k2, (j-1)*k2+1:j*k2) = (sig2+epsilon)*I_ksqrd;
        else
            Kstar((i-1)*k2+1:i*k2, (j-1)*k2+1:j*k2) = (sig2)*exp(-0.5*(ii-jj)^2/tau2)*I_ksqrd;
        end
            
    end
end