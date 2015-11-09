function xyzinpn = generate_data_PLDS_multiple_recordings(T,params,r,Model)
% generate data from PLDS for several recordings

%% if r==1, then use previous code

xyzinpn = cell(r, 1);

if r==1
    
    xyzinpn{1} = generate_data_PLDS(T, params);
    
else
    
    paramseach = params;
    
    if Model=='NSFR'
        for i=1:r
%             display('here');
            paramseach.h = params.h(:,i);
            xyzinpn{i} = generate_data_PLDS_NSFR(T, paramseach);
            
            % we assume T is same for each recording for now
            xyzinpn{i}.T = T;
        end
        
    else % NSLS
        
        for i=1:r
            paramseach.A = params.A(:,:,i);
            %         paramseach.AA = params.AA(:,:,i);
            %         paramseach.covA = params.covA(:,:,i);
            
            xyzinpn{i} = generate_data_PLDS(T, paramseach);
            
            % we assume T is same for each recording for now
            xyzinpn{i}.T = T;
        end
    end
    
end
