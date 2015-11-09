function fromVBEstep = runVBEstep(xyzinpn, initparams, r, Model)


if strcmp(Model, 'NSFR') || strcmp(Model, 'PLDS')  %We do the same for both NSFR and PLDS just in PLDS h will be zero and covh will be eye(k)
    
    opts = optimset('Display','off');
    
    
    if r>1
        
        fromVBEstep = cell(r,1);
        paramseach = cell(r,1);
        
        for i=1:r
            paramseach{i} = initparams;
            paramseach{i}.h = initparams.h(:,i);
            paramseach{i}.covh = initparams.covh(:,:,i);
            if size(initparams.inpn,3)==1
                paramseach{i}.inpn = initparams.inpn;
            else
                paramseach{i}.inpn = initparams.inpn(:,:,i);
            end
        end
    else
        paramseach = initparams;
    end
    % paramseach = initparams;
    
    %%
    
    if r>1
        
%         parpool; 
        for i=1:r
%             [i r]
            %fprintf(['E-step for trial ' num2str(i) ' started ' datestr(clock, 'yy-mm-dd-HH:MM:SS') '\n']);
            Yn = xyzinpn{i}.y;
                        
            %% (1) forward filtering
            
            ff = VBforwardfiltering_NSFR(Yn, paramseach{i}, opts);
            
%             T = size(Yn,2);
%             figure(1); subplot(211); plot(1:T, ff.mufwrd(1,:), 'r', 1:T, xyzinpn{i}.x(1,:), 'k');
%             subplot(212); plot(1:T, ff.mufwrd(2,:), 'r', 1:T, xyzinpn{i}.x(2,:), 'k');
            
            %% (2) backward smoothing
            
            bs = VBbackwardsmoothing_NSFR(Yn, paramseach{i}, opts);
            
%             figure(2); subplot(211); plot(1:T, bs.mubwrd(1,:), 'r', 1:T, xyzinpn{i}.x(1,:), 'k');
%             subplot(212); plot(1:T, bs.mubwrd(2,:), 'r', 1:T, xyzinpn{i}.x(2,:), 'k');

            
            %% (3) marginals / cross-covariance
            
            fromVBEstep{i} = VBcomputingmarginals_NSFR(Yn, ff, bs, paramseach{i});

%             figure(3); 
%             subplot(211); plot(1:T, fromVBEstep{i}.mumarg(1,:), 'r', 1:T, xyzinpn{i}.x(1,:), 'k');
%             subplot(212); plot(1:T, fromVBEstep{i}.mumarg(2,:), 'r', 1:T, xyzinpn{i}.x(2,:), 'k');
            
        end
        
%         myCluster = parcluster('local');
%         delete(myCluster.Jobs);
        
    else
        
        Yn = xyzinpn.y;
        
        %% (1) forward filtering
        
%         ff = VBforwardfiltering_A(Yn, initparams, opts);
        ff = VBforwardfiltering_NSFR(Yn, paramseach, opts);
        
        %% (2) backward smoothing
        
        bs =  VBbackwardsmoothing_NSFR(Yn, paramseach, opts);
        
        %% (3) marginals / cross-covariance
        
        fromVBEstep = VBcomputingmarginals_NSFR(Yn, ff, bs, paramseach);
        
    end
    
%% nonstationary A (without B)
else
   
    opts = optimset('Display','off');
    
    fromVBEstep = cell(r,1);
    
    paramseach = cell(r,1);
    for i=1:r
        paramseach{i} = initparams;
        paramseach{i}.A = initparams.A(:,:,i);
        paramseach{i}.AA = initparams.AA(:,:,i);
        paramseach{i}.covA = initparams.covA(:,:,i);
    end
    % paramseach = initparams;
    
    %%
    
    if r>1
        for i=1:r
            
            Yn = xyzinpn{i}.y;
            
            %     paramseach.A = initparams.A(:,:,i);
            %     paramseach.AA = initparams.AA(:,:,i);
            %     paramseach.covA = initparams.covA(:,:,i);
            
            %% (1) forward filtering
            
            ff = VBforwardfiltering_A(Yn, paramseach{i}, opts);
            
            %% (2) backward smoothing
            
            bs = VBbackwardsmoothing_A(Yn, paramseach{i}, opts);
            
            %% (3) marginals / cross-covariance
            
            fromVBEstep{i} = VBcomputingmarginals_A(Yn, ff, bs, paramseach{i});
            
        end
    else
        Yn = xyzinpn.y;
        
        %     paramseach.A = initparams.A(:,:,i);
        %     paramseach.AA = initparams.AA(:,:,i);
        %     paramseach.covA = initparams.covA(:,:,i);
        
        %% (1) forward filtering
        
        ff = VBforwardfiltering_A(Yn, initparams, opts);
        
        %% (2) backward smoothing
        
        bs = VBbackwardsmoothing_A(Yn, initparams, opts);
        
        %% (3) marginals / cross-covariance
        
        fromVBEstep = VBcomputingmarginals_A(Yn, ff, bs, initparams);
        
    end
    
end