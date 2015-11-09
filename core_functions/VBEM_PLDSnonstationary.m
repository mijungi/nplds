function datastruct = VBEM_PLDSnonstationary(xyzinpn, r, params, Model, varargin)
% everything = VBEM_PLDSnonstationary_A(xyzinpn, r, params);

% May 08, 2014
% wrote by Mijung Park

%% essential quantities

k = size(params.A(:,:,1),1);
p = size(params.C,1);
d = size(params.B,2);

%% initial params

% initparams = generate_params_from_prior(k, p, r);
ind_train = params.ind_train;
initparams = generate_params_from_prior(k, d, p, r, Model, ind_train, params.tau_init);
initparams.d = params.d;
initparams.inpn = params.inpn;

initparams.m_h = params.m_h_init;
initparams.tau2 = params.tau_init;
initparams.sig2 = params.sig_init;

if strcmp(Model, 'PLDS')
    initparams.h = zeros(size(initparams.h));
    initparams.covh = repmat(zeros(k),[1,1,r]);
end


%% VBEM (for now r==1)

nIter = 0;
if ~isfield(params, 'maxIter')
    maxIter = 3;
else
    maxIter = params.maxIter;
end

% to store mumarg, inv_sigmarg, crosscov, A, B, C, and mse values
fromVBEstepcell = cell(maxIter,1);
fromVBMstepcell = cell(maxIter,1);

%%
while (nIter<maxIter)
    
    %% update #Iter
    
    nIter = nIter + 1;
    
    fprintf(['Iteration ' num2str(nIter) '/' num2str(maxIter) ' started ' datestr(clock, 'yy-mm-dd-HH:MM:SS') '\n']);
    
    %Check if iteration has already been done
    if nargin > 4
        %varargin{1} is output folder for intermediate results
        if exist([varargin{1} filesep 'datastruct_' Model '_iter_' num2str(nIter) '.mat'], 'file')
            load([varargin{1} filesep 'datastruct_' Model '_iter_' num2str(nIter) '.mat'], 'datastruct_cur_iter');
            fromVBEstep = datastruct_cur_iter.Estep;
            fromVBMstep = datastruct_cur_iter.Mstep;
            initparams = update_params(fromVBMstep, initparams);
            fromVBEstepcell{nIter} = fromVBEstep;
            fromVBMstepcell{nIter} = fromVBMstep;
            continue;
        end
    end
    
    
    %% (1) VB E-step:
    
    fromVBEstep = runVBEstep(xyzinpn, initparams, r, Model);
    fprintf(['E-step done ' datestr(clock, 'yy-mm-dd-HH:MM:SS') '\n']);
    
    %% (2) VB M-step
    
    fromVBMstep = runVBMstep(xyzinpn, fromVBEstep, initparams, r, ind_train,Model);
    fprintf(['M-step done ' datestr(clock, 'yy-mm-dd-HH:MM:SS') '\n']);
    
    %% store these:
    
    fromVBEstepcell{nIter} = fromVBEstep;
    fromVBMstepcell{nIter} = fromVBMstep;
    
    initparams = update_params(fromVBMstep, initparams);
    
    %Save output at current iter into output folder
    if nargin > 4
        %varargin{1} is output folder for intermediate results
        datastruct_cur_iter.Estep = fromVBEstep;
        datastruct_cur_iter.Mstep = fromVBMstep;
        save([varargin{1} filesep 'datastruct_' Model '_iter_' num2str(nIter) '.mat'], 'datastruct_cur_iter');
    end
    
    
end


%% output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datastruct.Estep = fromVBEstepcell;
datastruct.Mstep = fromVBMstepcell;

% only if varargin is nonempty, we save the results:
if isempty(varargin)~=1
    save([varargin{1} filesep 'datastruct_' Model '_final.mat'], 'datastruct')
end

    function initparams = update_params(fromVBMstep, initparams)
        %% update params (that you want to update)
        
        initparams.h = fromVBMstep.h;
        initparams.covh = fromVBMstep.covh;
        
        initparams.C = fromVBMstep.C;
        
        initparams.A = fromVBMstep.A;
        initparams.AA = fromVBMstep.AA;
        initparams.covA = fromVBMstep.covA;
        initparams.B = fromVBMstep.B;
        initparams.AB = fromVBMstep.AB;
        initparams.covB = fromVBMstep.covB;
        
        initparams.alpha = fromVBMstep.alpha;
        initparams.beta = fromVBMstep.beta;
        initparams.sig2 = fromVBMstep.sig2;
        initparams.tau2 = fromVBMstep.tau2;
    end

end
