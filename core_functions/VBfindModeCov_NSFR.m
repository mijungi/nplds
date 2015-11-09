function [mufwrd, sigfwrd] = VBfindModeCov_NSFR(mutilde, sigtilde, C, d, h, covh, y, opts)
%                            VBfindModeCov_NSFR(mutilde(:,t), sigtilde(:,:,t), C, d, h, covh, Yn(:,t), opts);

maxexp = 1e2; 

%(1) update mu
a0 = zeros(length(mutilde),1);
fun = @(a)(obj_fwrd_NSFR(a, mutilde, sigtilde, y, C, d, h, covh, maxexp));
mufwrd = fsolve(fun, a0, opts);

%%
%(2) update V

quad = diag(C*covh*C');
g = exp(d + C*(mufwrd+h) + 0.5*quad);
% g = exp(min(d + C*(mufwrd+h) + 0.5*quad, maxexp*ones(p,1)));

p = size(C,1);
k = size(h,1);
cctrp = zeros(k,k,p);
for i=1:p
    cstrp = C(i,:);
    cctrp(:,:,i) = cstrp'*cstrp;
end

frstTrm = bsxfun(@times, cctrp, reshape(g,[1 1 p]));

frstTrm_in_cov = sum(frstTrm,3);

sigfwrd = inv(frstTrm_in_cov + inv(sigtilde)); 

% ===========================================================
function obj = obj_fwrd_NSFR(a, mutilde, sigtilde, y, C, d, h, covh, maxexp)

p = size(C,1);

quad = diag(C*covh*C');
secondTrm_in_mean = exp(min(d + C*(a+h) + 0.5*quad, maxexp*ones(p,1)));

obj = C'*y - C'*secondTrm_in_mean - sigtilde\(a - mutilde);