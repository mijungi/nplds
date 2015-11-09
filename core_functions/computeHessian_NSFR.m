function W = computeHessian_NSFR(mumarg, C, d, h, covh, maxexp)

% W =  computeHessian_NSFR(mumarg(:,1), C, d, h, covh, maxexp);

p = size(C,1);

quad = diag(C*covh*C');
% g = exp(d + C*(mumarg+h) + 0.5*quad);
g = exp(min(d + C*(mumarg+h) + 0.5*quad, maxexp*ones(p,1)));

p = size(C,1);
k = size(h,1);
cctrp = zeros(k,k,p);
for i=1:p
    cstrp = C(i,:);
    cctrp(:,:,i) = cstrp'*cstrp;
end

frstTrm = bsxfun(@times, cctrp, reshape(g,[1 1 p]));

W = sum(frstTrm,3);