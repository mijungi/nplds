function lam_tst = generate_one_step_ahead_prediction_tst_NSFR(fromEstep, fromMstep, xyzinpn, nsamps)
% one-step ahead prediction in the mean firing rate

% note:
% this is for no-input scenarios
% modify this for B = nonempty

%% (1) unpack required quantities

A = fromMstep.A;
C = fromMstep.C;
B = fromMstep.B;
d = fromMstep.d;
h = fromMstep.h;

mu = fromEstep.mumarg;
invcov = fromEstep.inv_sigmarg;

inpn = xyzinpn.inpn;

% if size(xyzinpn.inp_tst,1)==0
%     inpn_tst = xyzinpn.inp_tst;

k = size(A,1);
T = size(mu,2);
p = size(C,1);
pinp = size(B,2);
I = eye(k);

%% (2) generate last

lam_tst = zeros(p, T-1);

for t=1:T-1
    
    dummy = (A/invcov(:,:,t))*A'+I;
    covAnew = 0.5*(dummy + dummy');
    
    if pinp==0
        x_tst = mvnrnd(A*mu(:,t), covAnew, nsamps);
    else
        x_tst = mvnrnd(A*mu(:,t) + B*inpn(:,t+1), covAnew, nsamps);
    end
    
    xtstplush = bsxfun(@plus, x_tst', h);
    lam_tst(:, t) = mean(exp(bsxfun(@plus, C*xtstplush, d)),2);
    
end
    
% else
%     
%     
%     inpn_tst = xyzinpn.inp_tst;
%     
%     k = size(A,1);
%     T = size(mu,2);
%     p = size(C,1);
%     pinp = size(B,2);
%     I = eye(k);
%     
%     %% (2) generate last
%     
%     lam_tst = zeros(p, T);
%     
%     for t=1:T
%         
%         dummy = (A/invcov(:,:,t))*A'+I;
%         covAnew = 0.5*(dummy + dummy');
%         
%         if pinp==0
%             x_tst = mvnrnd(A*mu(:,t), covAnew, nsamps);
%         else
%             if t<T
%                 x_tst = mvnrnd(A*mu(:,t) + B*inpn(:,t+1), covAnew, nsamps);
%             else
%                 x_tst = mvnrnd(A*mu(:,t) + B*inpn_tst, covAnew, nsamps);
%             end
%         end
%         
%         lam_tst(:, t) = mean(exp(bsxfun(@plus, C*x_tst', d)),2);
%         
%     end
%     
% end
% 

