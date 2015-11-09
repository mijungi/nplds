%HARDCODED_DIRECTORY_PATHS
data_dir = '/nfs/data3/gergo/Mijung/Figure4/Output_data_retesting2';
data_file = '/nfs/data3/gergo/Mijung/Figure4/Input_data/alexdata_session2_org0_bin50.mat';
code_dir = '/nfs/nhome/live/gbohner/Git_main/nonstat_plds/nplds_code_NIPS2015/';

addpath functions

Model =  'NSFR'; %FITPARAM - {'NSFR' / 'PLDS'}

on_slurm = 1; % 1 - submit to slurm cluster. 0 - run on local machine (takes around 4 days on i7 2.4 Ghz cpu).

for k = 1:8; %FITPARAM - which dimensionalities to try
    for runnum = 0:9 %FITPARAM - How many runs / dimensionality for crossvalidation
        cd(code_dir);
        output_folder = [data_dir filesep 'k' num2str(k) '_run' num2str(runnum)];
        if ~exist(output_folder, 'dir'), mkdir(output_folder); end
        
        
        
        if on_slurm ==1
          %Create a SLURM script file in the corresponding folder with the
          %correct parameters, then call srun via system
          script_file_path = create_slurm_script( data_file, output_folder, code_dir, k, Model );

          system(['sbatch ' script_file_path]);
        
        else 
          %Call via command line of linux for local runs (within newly initialized Matlab instance) instead of submitting it to the cluster
          matlab_cmd = ['/opt/matlab-R2013b/bin/matlab -nodesktop -nosplash -singleCompThread -logfile ' output_folder filesep 'logfile.txt'];
          func_call = [' -r "cd Figure3_and_4/functions; runAlexsdata_prediction_func ' data_file ' ' output_folder ' ' code_dir ' ' num2str(k) ' ' Model '; exit;" '];
        
          system([matlab_cmd func_call]);
        end

        
    end
end

% exit;

